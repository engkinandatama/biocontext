import asyncio
import time
from typing import Any, Dict, List, Optional

from biocontext.adapters import (
    EnsemblAdapter,
    GeneOntologyAdapter,
    HGNCAdapter,
    LiteratureAdapter,
    MGIAdapter,
    MondoAdapter,
    NCBIAdapter,
    OpenTargetsAdapter,
    ReactomeAdapter,
    UniProtAdapter,
)
from biocontext.base import SQLiteCache
from biocontext.config import ScoringConfig
from biocontext.logging import get_logger
from biocontext.schemas import (
    BatchResolutionSummary,
    DiseaseEntity,
    FunctionalAnnotation,
    GOAnnotation,
    GeneEntity,
    LiteratureContext,
    MatchReason,
    PathwayContext,
    PathwayEntity,
    PublicationEntity,
    ResolutionContext,
    ResolutionResult,
    TargetAssociationContext,
    TargetDiseaseAssociation,
)

logger = get_logger("resolver")


class EntityResolver:
    """Multi-source biological entity resolution engine."""

    def __init__(
        self,
        cache: Optional[SQLiteCache] = None,
        ncbi_api_key: Optional[str] = None,
        email: Optional[str] = None,
        tool: Optional[str] = None
    ):
        self.cache = cache or SQLiteCache()
        self.hgnc = HGNCAdapter(cache=self.cache, email=email)
        self.ncbi = NCBIAdapter(cache=self.cache, api_key=ncbi_api_key, email=email, tool=tool)
        self.uniprot = UniProtAdapter(cache=self.cache, email=email)
        self.ensembl = EnsemblAdapter(cache=self.cache, email=email)
        self.mgi = MGIAdapter(cache=self.cache, email=email)
        self.go = GeneOntologyAdapter(cache=self.cache, email=email)
        self.reactome = ReactomeAdapter(cache=self.cache, email=email)
        self.mondo = MondoAdapter(cache=self.cache, email=email)
        self.opentargets = OpenTargetsAdapter(cache=self.cache, email=email)
        self.literature = LiteratureAdapter(cache=self.cache, email=email)





    def _apply_context_clues(
        self, gene, context: Optional[ResolutionContext], base_confidence: float, reasons: List[MatchReason]
    ) -> float:
        if not context or not gene:
            return base_confidence

        confidence = base_confidence

        # Chromosome disambiguation
        if context.chromosome and gene.location and gene.location.chromosome:
            hint_chr = context.chromosome.lower().replace("chr", "")
            gene_chr = gene.location.chromosome.lower().replace("chr", "")

            # Match prefix e.g. "17" matches "17p13.1" or "17q12"
            if gene_chr.startswith(hint_chr):
                confidence = min(1.0, confidence + ScoringConfig.CHROMOSOME_MATCH_BONUS)
                reasons.append(
                    MatchReason(
                        source="ContextClue",
                        rule="chromosome_match",
                        confidence=confidence,
                        details=f"Candidate chromosome '{gene.location.chromosome}' matches context hint '{context.chromosome}'"
                    )
                )
            else:
                confidence = max(0.1, confidence - ScoringConfig.CHROMOSOME_MISMATCH_PENALTY)
                reasons.append(
                    MatchReason(
                        source="ContextClue",
                        rule="chromosome_mismatch",
                        confidence=confidence,
                        details=f"Candidate chromosome '{gene.location.chromosome}' differs from context hint '{context.chromosome}'"
                    )
                )

        # Locus type / biotype disambiguation
        if context.locus_type and gene.locus_type:
            if context.locus_type.lower() not in gene.locus_type.lower():
                confidence = max(0.1, confidence - ScoringConfig.LOCUS_TYPE_MISMATCH_PENALTY)
                reasons.append(
                    MatchReason(
                        source="ContextClue",
                        rule="locus_type_mismatch",
                        confidence=confidence,
                        details=f"Locus type '{gene.locus_type}' does not match hint '{context.locus_type}'"
                    )
                )

        return round(confidence, 2)

    async def _enrich_entity(self, gene: GeneEntity, query_source: str = "generic") -> GeneEntity:
        """Universally reconcile cross-database IDs and annotations across HGNC, UniProt, and Ensembl."""
        if not gene or not gene.symbol:
            return gene

        # 1. Human (taxon 9606): cross-enrich with HGNC if key IDs or synonyms are missing
        if gene.taxon_id == 9606:
            needs_hgnc = not gene.hgnc_id or not gene.ncbi_gene_id or not gene.ensembl_gene_id or not gene.synonyms
            if needs_hgnc:
                try:
                    hgnc_match = await self.hgnc.fetch_by_symbol(gene.symbol)
                    if hgnc_match:
                        if not gene.hgnc_id and hgnc_match.hgnc_id:
                            gene.hgnc_id = hgnc_match.hgnc_id
                        if not gene.ncbi_gene_id and hgnc_match.ncbi_gene_id:
                            gene.ncbi_gene_id = hgnc_match.ncbi_gene_id
                        if not gene.ensembl_gene_id and hgnc_match.ensembl_gene_id:
                            gene.ensembl_gene_id = hgnc_match.ensembl_gene_id
                        if not gene.name and hgnc_match.name:
                            gene.name = hgnc_match.name
                        if hgnc_match.synonyms:
                            current_syns = set(gene.synonyms or [])
                            for s in hgnc_match.synonyms:
                                if s not in current_syns:
                                    gene.synonyms.append(s)
                                    current_syns.add(s)
                except Exception as e:
                    logger.debug(f"HGNC back-enrichment skipped: {e}")

            # 2. Reconcile UniProt Accession ID if missing
            if not gene.uniprot_ids:
                try:
                    u_res = await self.uniprot.resolve_gene(gene.symbol, taxon_id=9606)
                    if u_res and "accession" in u_res:
                        acc = u_res["accession"]
                        if acc not in gene.uniprot_ids:
                            gene.uniprot_ids.append(acc)
                except Exception as e:
                    logger.debug(f"UniProt back-enrichment skipped: {e}")

            # 3. Reconcile Ensembl model if transcripts or genomic coordinates are missing
            if (not gene.transcripts or not gene.location or not gene.location.start) and gene.ensembl_gene_id:
                try:
                    ens_match = await self.ensembl.fetch_by_id(gene.ensembl_gene_id)
                    if ens_match:
                        if not gene.transcripts and ens_match.transcripts:
                            gene.transcripts = ens_match.transcripts
                        if (not gene.location or not gene.location.start) and ens_match.location:
                            gene.location = ens_match.location
                except Exception as e:
                    logger.debug(f"Ensembl transcript back-enrichment skipped: {e}")

        # 4. Mouse (taxon 10090): cross-enrich MGI and Ensembl IDs
        elif gene.taxon_id == 10090:
            if not gene.mgi_id and gene.ensembl_gene_id:
                try:
                    mgi_id = await self.ensembl.fetch_mgi_id(gene.ensembl_gene_id)
                    if mgi_id:
                        gene.mgi_id = mgi_id
                except Exception as e:
                    logger.debug(f"MGI cross-enrichment skipped: {e}")
            elif not gene.ensembl_gene_id and gene.symbol:
                try:
                    ens_mouse = await self.ensembl.fetch_by_symbol(species="mus_musculus", symbol=gene.symbol)
                    if ens_mouse and ens_mouse.ensembl_gene_id:
                        gene.ensembl_gene_id = ens_mouse.ensembl_gene_id
                        if ens_mouse.transcripts and not gene.transcripts:
                            gene.transcripts = ens_mouse.transcripts
                except Exception as e:
                    logger.debug(f"Mouse Ensembl enrichment skipped: {e}")

        return gene

    async def resolve(
        self, query: str, taxon_id: int = 9606, context: Optional[ResolutionContext] = None
    ) -> ResolutionResult:
        query_clean = query.strip()
        reasons: List[MatchReason] = []

        # Top-level resolution cache (if no context clues supplied)
        cache_key = f"res:{query_clean}:{taxon_id}"
        if not context:
            cached_res = self.cache.get("resolver", cache_key)
            if cached_res:
                return ResolutionResult(**cached_res)

        # 1. Primary Authority for Human: HGNC
        if taxon_id == 9606 and not query_clean.isdigit():
            hgnc_res = await self.hgnc.resolve_gene(query_clean, taxon_id=taxon_id)
            if hgnc_res:
                gene = hgnc_res["entity"]
                rule = hgnc_res["rule"]
                confidence = hgnc_res["confidence"]
                
                reasons.append(
                    MatchReason(
                        source="HGNC",
                        rule=rule,
                        confidence=confidence,
                        details=f"Matched {rule} via HGNC REST API for symbol '{gene.symbol}'"
                    )
                )

                # Cross-reference with UniProt if uniprot_ids missing or to enrich
                await self._enrich_entity(gene, query_source="hgnc")

                # Apply context clues scoring adjustments if hints provided
                final_confidence = self._apply_context_clues(gene, context, confidence, reasons)
                match_status = "exact" if rule == "exact_symbol" else "alias"
                if final_confidence < 0.60:
                    match_status = "ambiguous"

                out = ResolutionResult(
                    query=query_clean,
                    match_status=match_status,
                    confidence_score=final_confidence,
                    resolved_entity=gene,
                    match_reasons=reasons
                )
                if not context:
                    self.cache.set("resolver", cache_key, out.model_dump())
                return out

        # 2. Check direct MGI ID (e.g. MGI:98834 or MGI:xxxx)
        if query_clean.upper().startswith("MGI:"):
            mgi_gene = await self.mgi.fetch_by_mgi_id(query_clean)
            if mgi_gene:
                reasons.append(
                    MatchReason(
                        source="MGI",
                        rule="mgi_id_lookup",
                        confidence=1.0,
                        details=f"Direct MGI ID lookup matched '{query_clean}' to mouse symbol '{mgi_gene.symbol}'"
                    )
                )
                await self._enrich_entity(mgi_gene, query_source="mgi")
                final_confidence = self._apply_context_clues(mgi_gene, context, 1.0, reasons)
                out = ResolutionResult(
                    query=query_clean,
                    match_status="exact",
                    confidence_score=final_confidence,
                    resolved_entity=mgi_gene,
                    match_reasons=reasons
                )
                if not context:
                    self.cache.set("resolver", cache_key, out.model_dump())
                return out

        # 3. Check Ensembl Gene ID (e.g. ENSG... or ENSMUSG...)
        if query_clean.upper().startswith("ENS"):
            ens_gene = await self.ensembl.fetch_by_id(query_clean)
            if ens_gene:
                reasons.append(
                    MatchReason(
                        source="Ensembl",
                        rule="ensembl_id_lookup",
                        confidence=1.0,
                        details=f"Direct Ensembl Gene ID lookup matched '{query_clean}' to symbol '{ens_gene.symbol}'"
                    )
                )
                # Universal reconciliation & back-enrichment (HGNC, UniProt, MGI)
                await self._enrich_entity(ens_gene, query_source="ensembl")

                final_confidence = self._apply_context_clues(ens_gene, context, 1.0, reasons)
                out = ResolutionResult(
                    query=query_clean,
                    match_status="exact",
                    confidence_score=final_confidence,
                    resolved_entity=ens_gene,
                    match_reasons=reasons
                )
                if not context:
                    self.cache.set("resolver", cache_key, out.model_dump())
                return out

        # 4. If query is numeric (Entrez ID) and human (9606), check HGNC first (authoritative and no strict 3 req/s throttle)
        if query_clean.isdigit() and taxon_id == 9606:
            hgnc_entrez_gene = await self.hgnc.fetch_by_entrez_id(query_clean)
            if hgnc_entrez_gene:
                reasons.append(
                    MatchReason(
                        source="HGNC",
                        rule="entrez_id_lookup",
                        confidence=1.0,
                        details=f"Matched Entrez Gene ID {query_clean} via HGNC to symbol '{hgnc_entrez_gene.symbol}'"
                    )
                )
                final_confidence = self._apply_context_clues(hgnc_entrez_gene, context, 1.0, reasons)
                out = ResolutionResult(
                    query=query_clean,
                    match_status="exact",
                    confidence_score=final_confidence,
                    resolved_entity=hgnc_entrez_gene,
                    match_reasons=reasons
                )
                if not context:
                    self.cache.set("resolver", cache_key, out.model_dump())
                return out

        # 5. Check Entrez Gene ID lookup or cross-species / fallback via NCBI
        ncbi_res = await self.ncbi.resolve_gene(query_clean, taxon_id=taxon_id)
        if ncbi_res:
            ncbi_gene = ncbi_res["entity"]
            rule = ncbi_res["rule"]
            confidence = ncbi_res["confidence"]
            reasons.append(
                MatchReason(
                    source="NCBI",
                    rule=rule,
                    confidence=confidence,
                    details=f"Resolved via NCBI Entrez Gene (ID: {ncbi_gene.ncbi_gene_id})"
                )
            )

            # If human, try to enrich with HGNC data (symbol, ensembl, location)
            if taxon_id == 9606 and ncbi_gene.symbol:
                hgnc_gene = await self.hgnc.fetch_by_symbol(ncbi_gene.symbol)
                if hgnc_gene:
                    if not hgnc_gene.ncbi_gene_id:
                        hgnc_gene.ncbi_gene_id = ncbi_gene.ncbi_gene_id
                    out = ResolutionResult(
                        query=query_clean,
                        match_status="exact" if query_clean.isdigit() else "ncbi_matched",
                        confidence_score=confidence,
                        resolved_entity=hgnc_gene,
                        match_reasons=reasons
                    )
                    if not context:
                        self.cache.set("resolver", cache_key, out.model_dump())
                    return out

            # If mouse (10090), check Ensembl to enrich MGI ID and transcripts
            if taxon_id == 10090 and ncbi_gene.symbol:
                ens_mouse = await self.ensembl.fetch_by_symbol(species="mus_musculus", symbol=ncbi_gene.symbol)
                if ens_mouse:
                    if ens_mouse.ensembl_gene_id:
                        ncbi_gene.ensembl_gene_id = ens_mouse.ensembl_gene_id
                        mgi_id = await self.ensembl.fetch_mgi_id(ens_mouse.ensembl_gene_id)
                        if mgi_id:
                            ncbi_gene.mgi_id = mgi_id
                    if ens_mouse.transcripts:
                        ncbi_gene.transcripts = ens_mouse.transcripts

            out = ResolutionResult(
                query=query_clean,
                match_status="exact" if query_clean.isdigit() else "ncbi_matched",
                confidence_score=confidence,
                resolved_entity=ncbi_gene,
                match_reasons=reasons
            )
            if not context:
                self.cache.set("resolver", cache_key, out.model_dump())
            return out

        # 6. Check UniProt directly (accession lookup e.g. P04637)
        uniprot_protein = await self.uniprot.fetch_by_accession(query_clean)
        if uniprot_protein and uniprot_protein.gene_symbol:
            reasons.append(
                MatchReason(
                    source="UniProt",
                    rule="accession_lookup",
                    confidence=0.95,
                    details=f"UniProt accession {uniprot_protein.accession} maps to gene {uniprot_protein.gene_symbol}"
                )
            )
            # Resolve gene symbol back to HGNC
            gene = await self.hgnc.fetch_by_symbol(uniprot_protein.gene_symbol)
            if gene:
                if uniprot_protein.accession not in gene.uniprot_ids:
                    gene.uniprot_ids.append(uniprot_protein.accession)
                await self._enrich_entity(gene, query_source="uniprot")
                out = ResolutionResult(
                    query=query_clean,
                    match_status="exact",
                    confidence_score=0.95,
                    resolved_entity=gene,
                    match_reasons=reasons
                )
                if not context:
                    self.cache.set("resolver", cache_key, out.model_dump())
                return out

        # 5. Fuzzy Match Engine: Check for common typos (e.g. TP54 -> TP53, BRCA -> BRCA1)
        if taxon_id == 9606 and len(query_clean) >= 3 and not query_clean.isdigit():
            fuzzy_matches = await self.hgnc.search_fuzzy(
                query_clean, max_distance=ScoringConfig.FUZZY_MAX_DISTANCE
            )
            if fuzzy_matches:
                best_match = fuzzy_matches[0]
                reasons.append(
                    MatchReason(
                        source="HGNC",
                        rule="fuzzy_levenshtein_match",
                        confidence=ScoringConfig.FUZZY_MATCH_CONFIDENCE,
                        details=f"Approximate typo match resolved query '{query_clean}' to candidate '{best_match.symbol}'"
                    )
                )
                final_confidence = self._apply_context_clues(
                    best_match, context, ScoringConfig.FUZZY_MATCH_CONFIDENCE, reasons
                )
                return ResolutionResult(
                    query=query_clean,
                    match_status="fuzzy",
                    confidence_score=final_confidence,
                    resolved_entity=best_match,
                    alternative_matches=fuzzy_matches[1:],
                    match_reasons=reasons
                )

        return ResolutionResult(
            query=query_clean,
            match_status="unresolved",
            confidence_score=0.0,
            resolved_entity=None,
            match_reasons=[
                MatchReason(
                    source="BioContext",
                    rule="no_match",
                    confidence=0.0,
                    details=f"No matching authoritative entity found for query '{query_clean}'"
                )
            ]
        )

    async def annotate_gene(
        self,
        query: str,
        taxon_id: int = 9606,
        max_terms_per_aspect: int = 25,
        target_aspect: Optional[str] = None
    ) -> Optional[FunctionalAnnotation]:
        """Fetch functional annotations (GO terms) for a gene query across MF, BP, and CC.

        Resolves the gene to an authoritative UniProt accession first, then queries QuickGO.
        """
        clean_q = query.strip()
        uniprot_acc: Optional[str] = None

        # 1. If query is already a UniProt accession (e.g. P04637)
        if len(clean_q) in (6, 10) and (clean_q[0].isalpha() and clean_q[-1].isalnum()):
            protein = await self.uniprot.fetch_by_accession(clean_q)
            if protein:
                uniprot_acc = protein.accession

        # 2. If not an accession, resolve via EntityResolver
        if not uniprot_acc:
            res = await self.resolve(clean_q, taxon_id=taxon_id)
            if res and res.resolved_entity:
                entity = res.resolved_entity
                if entity.uniprot_ids:
                    uniprot_acc = entity.uniprot_ids[0]

        # 3. Fallback: try UniProt adapter resolution directly
        if not uniprot_acc:
            u_res = await self.uniprot.resolve_gene(clean_q, taxon_id=taxon_id)
            if u_res and "entity" in u_res:
                ent = u_res["entity"]
                if hasattr(ent, "uniprot_ids") and ent.uniprot_ids:
                    uniprot_acc = ent.uniprot_ids[0]

        if not uniprot_acc:
            logger.warning("Could not map gene query '%s' to a UniProt accession for GO annotations", clean_q)
            return None

        # Fetch annotations from QuickGO
        aspects_to_fetch = [target_aspect] if target_aspect else ["molecular_function", "biological_process", "cellular_component"]
        all_annots: List[GOAnnotation] = []

        for asp in aspects_to_fetch:
            annots = await self.go.fetch_annotations(
                gene_product_id=uniprot_acc,
                taxon_id=taxon_id,
                aspect=asp,
                limit=max_terms_per_aspect
            )
            all_annots.extend(annots)

        # Categorize by aspect
        mf_list: List[GOAnnotation] = []
        bp_list: List[GOAnnotation] = []
        cc_list: List[GOAnnotation] = []

        for a in all_annots:
            aspect_norm = a.aspect.lower()
            if "molecular_function" in aspect_norm or aspect_norm == "mf":
                mf_list.append(a)
            elif "biological_process" in aspect_norm or aspect_norm == "bp":
                bp_list.append(a)
            elif "cellular_component" in aspect_norm or aspect_norm == "cc":
                cc_list.append(a)

        return FunctionalAnnotation(
            query=clean_q,
            taxon_id=taxon_id,
            uniprot_accession=uniprot_acc,
            molecular_functions=mf_list[:max_terms_per_aspect],
            biological_processes=bp_list[:max_terms_per_aspect],
            cellular_components=cc_list[:max_terms_per_aspect],
            total_annotations=len(mf_list) + len(bp_list) + len(cc_list)
        )

    async def get_pathways(
        self,
        query: str,
        taxon_id: int = 9606,
        species: str = "Homo sapiens",
        limit: int = 25
    ) -> Optional[PathwayContext]:
        """Fetch biological pathways associated with a gene or protein query from Reactome.

        Resolves gene to a UniProt accession first, then queries Reactome Content Service.
        """
        clean_q = query.strip()
        uniprot_acc: Optional[str] = None

        # 1. If query is already a UniProt accession (e.g. P04637)
        if len(clean_q) in (6, 10) and (clean_q[0].isalpha() and clean_q[-1].isalnum()):
            protein = await self.uniprot.fetch_by_accession(clean_q)
            if protein:
                uniprot_acc = protein.accession

        # 2. If not an accession, resolve via EntityResolver
        if not uniprot_acc:
            res = await self.resolve(clean_q, taxon_id=taxon_id)
            if res and res.resolved_entity:
                entity = res.resolved_entity
                if entity.uniprot_ids:
                    uniprot_acc = entity.uniprot_ids[0]

        # 3. Fallback: try UniProt adapter directly
        if not uniprot_acc:
            u_res = await self.uniprot.resolve_gene(clean_q, taxon_id=taxon_id)
            if u_res and "entity" in u_res:
                ent = u_res["entity"]
                if hasattr(ent, "uniprot_ids") and ent.uniprot_ids:
                    uniprot_acc = ent.uniprot_ids[0]

        if not uniprot_acc:
            logger.warning("Could not map gene query '%s' to a UniProt accession for pathway mapping", clean_q)
            return None

        # Fetch pathways from Reactome
        pathways = await self.reactome.fetch_pathways_by_uniprot(
            uniprot_accession=uniprot_acc,
            species=species,
            limit=limit
        )

        return PathwayContext(
            query=clean_q,
            taxon_id=taxon_id,
            uniprot_accession=uniprot_acc,
            source="Reactome",
            pathways=pathways,
            total_pathways=len(pathways)
        )

    async def resolve_batch(
        self,
        queries: List[str],
        taxon_id: int = 9606,
        concurrency: int = 10
    ) -> BatchResolutionSummary:
        """Resolve a batch of gene queries concurrently with semaphore rate limiting.

        Designed to process high-throughput gene lists (e.g. differential expression results).
        """
        import asyncio
        import time

        clean_queries = [q.strip() for q in queries if q and q.strip()]
        if not clean_queries:
            return BatchResolutionSummary(
                total_queries=0,
                resolved_count=0,
                unresolved_count=0,
                success_rate=0.0,
                execution_time_seconds=0.0,
                results=[]
            )

        semaphore = asyncio.Semaphore(max(1, concurrency))
        start_time = time.perf_counter()

        async def _bounded_resolve(query: str) -> ResolutionResult:
            async with semaphore:
                try:
                    return await self.resolve(query=query, taxon_id=taxon_id)
                except Exception as e:
                    logger.error("Error during batch resolution of '%s' | error=%s", query, str(e))
                    return ResolutionResult(
                        query=query,
                        match_status="unresolved",
                        confidence_score=0.0,
                        resolved_entity=None,
                        match_reasons=[
                            MatchReason(
                                source="BioContextBatch",
                                rule="exception_error",
                                confidence=0.0,
                                details=f"Batch resolution raised exception: {str(e)}"
                            )
                        ]
                    )

        tasks = [_bounded_resolve(q) for q in clean_queries]
        results = await asyncio.gather(*tasks)

        elapsed = time.perf_counter() - start_time
        resolved = sum(1 for r in results if r.match_status != "unresolved" and r.resolved_entity is not None)
        unresolved = len(results) - resolved
        success_rate = round(resolved / len(results), 4) if results else 0.0

        return BatchResolutionSummary(
            total_queries=len(results),
            resolved_count=resolved,
            unresolved_count=unresolved,
            success_rate=success_rate,
            execution_time_seconds=round(elapsed, 3),
            results=results
        )

    async def resolve_disease(self, query: str, limit: int = 5) -> List[DiseaseEntity]:
        """Resolve a disease name, synonym, or keyword against MONDO Disease Ontology."""
        clean_q = query.strip()
        if not clean_q:
            return []

        # If query is direct MONDO identifier e.g. MONDO:0018875
        if clean_q.upper().startswith("MONDO:") or clean_q.upper().startswith("MONDO_"):
            entity = await self.mondo.fetch_by_id(clean_q)
            return [entity] if entity else []

        return await self.mondo.search_disease(clean_q, limit=limit)

    async def get_target_diseases(self, gene_query: str, limit: int = 10) -> Optional[TargetAssociationContext]:
        """Retrieve evidence-backed disease associations for a target gene from Open Targets."""
        clean_q = gene_query.strip()
        if not clean_q:
            return None

        ensembl_id = None
        symbol = None

        if clean_q.upper().startswith("ENSG"):
            ensembl_id = clean_q.upper()
        else:
            # Resolve gene entity to acquire authoritative Ensembl Gene ID
            res = await self.resolve(clean_q)
            if res and res.resolved_entity:
                ensembl_id = res.resolved_entity.ensembl_gene_id
                symbol = res.resolved_entity.symbol

        if not ensembl_id:
            logger.warning("Could not map gene query '%s' to an Ensembl Gene ID for Open Targets", clean_q)
            return None

        return await self.opentargets.fetch_target_diseases(
            ensembl_gene_id=ensembl_id,
            symbol=symbol,
            limit=limit
        )

    async def get_supporting_publications(self, query: str, limit: int = 5) -> LiteratureContext:
        """Fetch authoritative supporting scientific publications for a gene, disease, or biomedical query."""
        clean_q = query.strip()
        if not clean_q:
            return LiteratureContext(query=query, total_hits=0, publications=[])

        # If query is a gene, try resolving to approved symbol to enrich query
        resolved_sym = None
        if not clean_q.isdigit() and not clean_q.upper().startswith("PMC") and not clean_q.startswith("10."):
            res = await self.resolve(clean_q)
            if res and res.resolved_entity:
                resolved_sym = res.resolved_entity.symbol

        search_term = resolved_sym or clean_q
        return await self.literature.search_publications(query=search_term, limit=limit)

    async def get_publication_details(self, identifier: str) -> Optional[PublicationEntity]:
        """Fetch detailed publication metadata by PMID, PMCID, or DOI."""
        return await self.literature.fetch_by_id(identifier)



