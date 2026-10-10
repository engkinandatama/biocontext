"""Adapters for authoritative biological databases: HGNC, NCBI, and UniProt."""

import asyncio
import logging
from typing import Any, Dict, List, Optional
import httpx

import os

from biocontext.base import AsyncRateLimiter, BaseBioAdapter, SQLiteCache
from biocontext.config import ClientConfig, RateLimitConfig
from biocontext.logging import get_logger
from biocontext.schemas import (
    DiseaseEntity,
    ExonEntity,
    FunctionalAnnotation,
    GOAnnotation,
    GeneEntity,
    GenomicLocation,
    OrthologEntity,
    PathwayContext,
    PathwayEntity,
    ProteinEntity,
    PublicationEntity,
    LiteratureContext,
    TargetAssociationContext,
    TargetDiseaseAssociation,
    TranscriptEntity,
)


logger = get_logger("adapters")



class HGNCAdapter(BaseBioAdapter):
    """Adapter for HGNC (HUGO Gene Nomenclature Committee) REST API.
    Primary authority for human gene nomenclature.
    """

    BASE_URL = "https://rest.genenames.org"

    def __init__(self, cache: Optional[SQLiteCache] = None, email: Optional[str] = None):
        super().__init__(name="HGNC", cache=cache)
        self.email = ClientConfig.get_email(email)
        self.headers = ClientConfig.get_headers(self.email)


    async def fetch_by_symbol(self, symbol: str) -> Optional[GeneEntity]:
        cache_key = f"symbol:{symbol.upper()}"
        cached = self.cache.get("hgnc", cache_key)
        if cached:
            return GeneEntity(**cached)

        url = f"{self.BASE_URL}/fetch/symbol/{symbol.upper()}"
        async with httpx.AsyncClient(timeout=10.0) as client:
            resp = await client.get(url, headers=self.headers)
            if resp.status_code != 200:
                return None
            data = resp.json()

        docs = data.get("response", {}).get("docs", [])
        if not docs:
            return None

        doc = docs[0]
        uniprot_ids = doc.get("uniprot_ids", [])
        aliases = doc.get("alias_symbol", [])
        prev_symbols = doc.get("prev_symbol", [])
        all_synonyms = list(dict.fromkeys(aliases + prev_symbols))

        location = None
        if "location" in doc:
            location = GenomicLocation(
                chromosome=doc.get("location", ""),
                assembly="GRCh38"
            )

        gene = GeneEntity(
            symbol=doc.get("symbol"),
            name=doc.get("name"),
            taxon_id=9606,
            species="Homo sapiens",
            hgnc_id=doc.get("hgnc_id"),
            ncbi_gene_id=doc.get("entrez_id"),
            ensembl_gene_id=doc.get("ensembl_gene_id"),
            uniprot_ids=uniprot_ids,
            synonyms=all_synonyms,
            locus_type=doc.get("locus_type"),
            location=location,
        )

        self.cache.set("hgnc", cache_key, gene.model_dump())
        return gene

    async def fetch_by_entrez_id(self, entrez_id: str) -> Optional[GeneEntity]:
        cache_key = f"entrez:{entrez_id}"
        cached = self.cache.get("hgnc", cache_key)
        if cached:
            return GeneEntity(**cached)

        url = f"{self.BASE_URL}/fetch/entrez_id/{entrez_id}"
        async with httpx.AsyncClient(timeout=10.0) as client:
            try:
                resp = await client.get(url, headers=self.headers)
                if resp.status_code != 200:
                    return None
                data = resp.json()
            except Exception as e:
                logger.warning(f"HGNC fetch_by_entrez_id failed | error={e}")
                return None

        docs = data.get("response", {}).get("docs", [])
        if not docs:
            return None

        doc = docs[0]
        uniprot_ids = doc.get("uniprot_ids", [])
        aliases = doc.get("alias_symbol", [])
        prev_symbols = doc.get("prev_symbol", [])
        all_synonyms = list(dict.fromkeys(aliases + prev_symbols))

        location = None
        if "location" in doc:
            location = GenomicLocation(
                chromosome=doc.get("location", ""),
                assembly="GRCh38"
            )

        gene = GeneEntity(
            symbol=doc.get("symbol"),
            name=doc.get("name"),
            taxon_id=9606,
            species="Homo sapiens",
            hgnc_id=doc.get("hgnc_id"),
            ncbi_gene_id=str(doc.get("entrez_id")),
            ensembl_gene_id=doc.get("ensembl_gene_id"),
            uniprot_ids=uniprot_ids,
            synonyms=all_synonyms,
            locus_type=doc.get("locus_type"),
            location=location,
        )

        self.cache.set("hgnc", cache_key, gene.model_dump())
        return gene

    async def search_alias(self, alias: str) -> Optional[GeneEntity]:
        cache_key = f"alias:{alias.upper()}"
        cached = self.cache.get("hgnc", cache_key)
        if cached:
            return GeneEntity(**cached)

        # 1. First check prev_symbol (exact historical symbols have higher clinical specificity)
        url_prev = f"{self.BASE_URL}/search/prev_symbol/{alias.upper()}"
        async with httpx.AsyncClient(timeout=10.0) as client:
            try:
                resp_prev = await client.get(url_prev, headers=self.headers)
                if resp_prev.status_code == 200:
                    docs = resp_prev.json().get("response", {}).get("docs", [])
                    if docs:
                        primary_symbol = docs[0].get("symbol")
                        if primary_symbol:
                            return await self.fetch_by_symbol(primary_symbol)
            except Exception as e:
                logger.warning(f"HGNC prev_symbol lookup failed | error={e}")

        # 2. Check alias_symbol (use /fetch endpoint to get full metadata for accurate ranking)
        url = f"{self.BASE_URL}/fetch/alias_symbol/{alias.upper()}"
        async with httpx.AsyncClient(timeout=10.0) as client:
            try:
                resp = await client.get(url, headers=self.headers)
                if resp.status_code != 200:
                    return None
                docs = resp.json().get("response", {}).get("docs", [])
            except Exception as e:
                logger.warning(f"HGNC alias_symbol lookup failed | error={e}")
                return None

        if not docs:
            return None

        # Prioritize:
        # 1. Exact case match in alias_symbol (e.g. "p21" vs "P21")
        # 2. Mention in prev_name / alias_name (e.g. CDKN1A previous name was "cyclin-dependent kinase inhibitor 1A (p21, Cip1)")
        # 3. Protein-coding genes over pseudogenes
        # 4. Authority cross-references and citation prominence
        def alias_rank_key(d: Dict[str, Any]) -> tuple:
            aliases = d.get("alias_symbol", [])
            exact_case = 1 if alias in aliases else 0
            
            # Check prev_name / alias_name mention
            text_context = " ".join(d.get("prev_name", []) + d.get("alias_name", []) + [d.get("name", "")]).lower()
            name_mention = 1 if alias.lower() in text_context else 0

            is_protein_coding = 1 if d.get("locus_group") == "protein-coding gene" or d.get("locus_type") == "gene with protein product" else 0
            pubmed_count = len(d.get("pubmed_id", []))
            has_omim = 1 if d.get("omim_id") else 0
            has_mane = 1 if d.get("mane_select") else 0
            
            return (name_mention, exact_case, is_protein_coding, has_omim + has_mane, pubmed_count)

        docs.sort(key=alias_rank_key, reverse=True)

        primary_symbol = docs[0].get("symbol")
        if primary_symbol:
            resolved_gene = await self.fetch_by_symbol(primary_symbol)
            if resolved_gene:
                self.cache.set("hgnc", cache_key, resolved_gene.model_dump())
            return resolved_gene
        return None

    async def search_fuzzy(self, query: str, max_distance: int = 2) -> List[GeneEntity]:
        """Search HGNC for approximate matches using prefix/wildcard search and Levenshtein distance."""
        query_upper = query.upper().strip()
        if len(query_upper) < 3:
            return []

        # Prefix wildcard search on HGNC (e.g. TP5* or BRC*)
        prefix = query_upper[:max(3, len(query_upper) - 1)]
        url = f"{self.BASE_URL}/search/symbol/{prefix}*"
        candidates: List[str] = []

        async with httpx.AsyncClient(timeout=10.0) as client:
            try:
                resp = await client.get(url, headers=self.headers)
                if resp.status_code == 200:
                    docs = resp.json().get("response", {}).get("docs", [])
                    for d in docs:
                        sym = d.get("symbol")
                        if sym:
                            candidates.append(sym)
            except Exception as e:
                logger.warning(f"HGNC wildcard search failed | query={query} error={e}")

        # Levenshtein distance calculation
        def levenshtein(s1: str, s2: str) -> int:
            if len(s1) < len(s2):
                return levenshtein(s2, s1)
            if len(s2) == 0:
                return len(s1)
            previous_row = range(len(s2) + 1)
            for i, c1 in enumerate(s1):
                current_row = [i + 1]
                for j, c2 in enumerate(s2):
                    insertions = previous_row[j + 1] + 1
                    deletions = current_row[j] + 1
                    substitutions = previous_row[j] + (c1 != c2)
                    current_row.append(min(insertions, deletions, substitutions))
                previous_row = current_row
            return previous_row[-1]

        matched_symbols = []
        for cand in candidates:
            dist = levenshtein(query_upper, cand.upper())
            if 0 < dist <= max_distance:
                matched_symbols.append((dist, cand))

        matched_symbols.sort(key=lambda x: x[0])
        results: List[GeneEntity] = []
        for _, sym in matched_symbols[:5]:
            gene = await self.fetch_by_symbol(sym)
            if gene:
                results.append(gene)

        return results


    async def resolve_gene(self, query: str, taxon_id: int = 9606) -> Optional[Dict[str, Any]]:
        if taxon_id != 9606:
            return None
        gene = await self.fetch_by_symbol(query)
        if gene:
            return {"entity": gene, "rule": "exact_symbol", "confidence": 1.0}
        alias_gene = await self.search_alias(query)
        if alias_gene:
            return {"entity": alias_gene, "rule": "alias_match", "confidence": 0.85}
        return None


class UniProtAdapter(BaseBioAdapter):
    """Adapter for UniProt REST API."""

    BASE_URL = "https://rest.uniprot.org/uniprotkb"

    def __init__(self, cache: Optional[SQLiteCache] = None, email: Optional[str] = None):
        super().__init__(name="UniProt", cache=cache)
        self.email = ClientConfig.get_email(email)
        self.headers = ClientConfig.get_headers(self.email)

    async def fetch_by_accession(self, accession: str) -> Optional[ProteinEntity]:
        cache_key = f"acc:{accession.upper()}"
        cached = self.cache.get("uniprot", cache_key)
        if cached:
            return ProteinEntity(**cached)

        url = f"{self.BASE_URL}/{accession.upper()}.json"
        data = None
        for attempt in range(3):
            try:
                async with httpx.AsyncClient(timeout=15.0) as client:
                    resp = await client.get(url, headers=self.headers)
                    if resp.status_code != 200:
                        return None
                    data = resp.json()
                    break
            except (httpx.TimeoutException, httpx.NetworkError):
                if attempt == 2:
                    return None
                await asyncio.sleep(0.5 * (attempt + 1))

        if not data:
            return None

        primary_acc = data.get("primaryAccession", accession)
        entry_name = data.get("uniProtkbId")
        taxon_id = data.get("organism", {}).get("taxonId", 9606)
        species = data.get("organism", {}).get("scientificName", "Homo sapiens")

        # Protein names
        protein_desc = data.get("proteinDescription", {})
        rec_name = protein_desc.get("recommendedName", {}).get("fullName", {}).get("value")

        # Gene symbol
        genes = data.get("genes", [])
        gene_symbol = None
        if genes:
            gene_symbol = genes[0].get("geneName", {}).get("value")

        # Sequence details
        seq = data.get("sequence", {})
        seq_len = seq.get("length")
        mol_weight = seq.get("molWeight")

        protein = ProteinEntity(
            accession=primary_acc,
            entry_name=entry_name,
            name=rec_name,
            taxon_id=taxon_id,
            species=species,
            gene_symbol=gene_symbol,
            sequence_length=seq_len,
            molecular_weight=mol_weight,
        )

        self.cache.set("uniprot", cache_key, protein.model_dump())
        return protein

    async def resolve_gene(self, query: str, taxon_id: int = 9606) -> Optional[Dict[str, Any]]:
        # UniProt query search for human reviewed entries
        cache_key = f"query:{query.upper()}:{taxon_id}"
        cached = self.cache.get("uniprot", cache_key)
        if cached:
            return cached

        search_url = f"{self.BASE_URL}/search"
        params = {
            "query": f"(gene:{query}) AND (taxonomy_id:{taxon_id}) AND (reviewed:true)",
            "format": "json",
            "size": 1
        }
        async with httpx.AsyncClient(timeout=10.0) as client:
            resp = await client.get(search_url, params=params, headers=self.headers)
            if resp.status_code != 200:
                return None
            results = resp.json().get("results", [])
            if not results:
                return None
            accession = results[0].get("primaryAccession")
            if accession:
                protein = await self.fetch_by_accession(accession)
                if protein:
                    res = {"protein": protein.model_dump(), "accession": accession}
                    self.cache.set("uniprot", cache_key, res)
                    return res
        return None


class NCBIAdapter(BaseBioAdapter):
    """Adapter for NCBI E-utilities (Entrez Gene).
    Authority for cross-species genes and Entrez Gene IDs.
    """

    ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"

    def __init__(
        self,
        cache: Optional[SQLiteCache] = None,
        api_key: Optional[str] = None,
        email: Optional[str] = None,
        tool: Optional[str] = None
    ):
        super().__init__(name="NCBI", cache=cache)
        # Dynamically read from parameter or environment variable
        self.api_key = api_key or os.environ.get("NCBI_API_KEY")
        self.email = ClientConfig.get_email(email)
        self.tool = ClientConfig.get_tool(tool)
        self.headers = ClientConfig.get_headers(self.email)

        # Rate limit: 10 req/s with API key, 3 req/s without key (with safety margin: 2.8 req/s)
        rate = 9.5 if self.api_key else 2.8
        self.rate_limiter = AsyncRateLimiter(requests_per_second=rate)
        if self.api_key:
            logger.info("NCBIAdapter initialized with API key | rate_limit=10_req_sec")
        else:
            logger.info("NCBIAdapter initialized without API key | rate_limit=3_req_sec")

    def _params(self, extra: Dict[str, Any]) -> Dict[str, Any]:
        p = {"retmode": "json", "tool": self.tool}
        if self.email:
            p["email"] = self.email
        if self.api_key:
            p["api_key"] = self.api_key
        p.update(extra)
        return p


    async def fetch_by_id(self, gene_id: str) -> Optional[GeneEntity]:
        cache_key = f"gene_id:{gene_id}"
        cached = self.cache.get("ncbi", cache_key)
        if cached:
            return GeneEntity(**cached)

        params = self._params({"db": "gene", "id": gene_id})
        for attempt in range(3):
            try:
                await self.rate_limiter.acquire()
                async with httpx.AsyncClient(timeout=15.0) as client:
                    resp = await client.get(self.ESUMMARY_URL, params=params, headers=self.headers)
                    if resp.status_code == 200:
                        data = resp.json()
                        break
                    elif resp.status_code == 429:
                        logger.warning("NCBI rate limit 429 encountered, cooling down...")
                        import asyncio
                        await asyncio.sleep(1.0)
            except Exception as e:
                logger.warning(f"NCBI summary attempt {attempt + 1} failed | error={e}")
                if attempt == 2:
                    return None
                import asyncio
                await asyncio.sleep(0.5 * (attempt + 1))
        else:
            return None

        result = data.get("result", {})
        uids = result.get("uids", [])
        if not uids or gene_id not in result:
            return None

        doc = result[gene_id]
        symbol = doc.get("name")
        description = doc.get("description")
        taxon_id = doc.get("taxid", 9606)
        organism = doc.get("organism", {}).get("scientificname", "Homo sapiens")
        aliases_str = doc.get("otheraliases", "")
        aliases = [a.strip() for a in aliases_str.split(",") if a.strip()] if aliases_str else []

        gene = GeneEntity(
            symbol=symbol,
            name=description,
            taxon_id=taxon_id,
            species=organism,
            ncbi_gene_id=str(gene_id),
            synonyms=aliases,
        )

        self.cache.set("ncbi", cache_key, gene.model_dump())
        return gene

    async def resolve_gene(self, query: str, taxon_id: int = 9606) -> Optional[Dict[str, Any]]:
        # If query is purely digits, treat directly as Entrez Gene ID
        if query.isdigit():
            gene = await self.fetch_by_id(query)
            if gene:
                return {"entity": gene, "rule": "ncbi_id_lookup", "confidence": 1.0}

        # Otherwise search by gene symbol & taxid
        cache_key = f"search:{query.upper()}:{taxon_id}"
        cached = self.cache.get("ncbi", cache_key)
        if cached:
            return {"entity": GeneEntity(**cached), "rule": "ncbi_search", "confidence": 0.9}

        term = f"({query}[Gene Name]) AND {taxon_id}[Taxonomy ID]"
        params = self._params({"db": "gene", "term": term, "retmax": 1})
        for attempt in range(3):
            try:
                await self.rate_limiter.acquire()
                async with httpx.AsyncClient(timeout=15.0) as client:
                    resp = await client.get(self.ESEARCH_URL, params=params, headers=self.headers)
                    if resp.status_code == 200:
                        id_list = resp.json().get("esearchresult", {}).get("idlist", [])
                        break
                    elif resp.status_code == 429:
                        logger.warning("NCBI rate limit 429 encountered, cooling down...")
                        import asyncio
                        await asyncio.sleep(1.0)
            except Exception as e:
                logger.warning(f"NCBI search attempt {attempt + 1} failed | error={e}")
                if attempt == 2:
                    return None
                import asyncio
                await asyncio.sleep(0.5 * (attempt + 1))
        else:
            return None

        if not id_list:
            return None

        gene = await self.fetch_by_id(id_list[0])
        if gene:
            self.cache.set("ncbi", cache_key, gene.model_dump())
            return {"entity": gene, "rule": "ncbi_search", "confidence": 0.9}
        return None


class EnsemblAdapter(BaseBioAdapter):
    """Adapter for Ensembl REST API (v15+).
    Primary authority for genomic coordinates, gene models, canonical transcripts, and isoforms.

    Note on Protocol Architecture:
        Full migration to Ensembl Beta GraphQL (beta.ensembl.org/graphql) was evaluated for v0.5.1
        (Issue #4). However, Ensembl's beta GraphQL schema currently lacks cross-species homology/ortholog
        resolvers and remains experimental. BioContext therefore maintains the production REST interface
        with an adaptive multi-tier fallback pipeline (intraservice fallback from /homology/id to
        /homology/symbol and from expand=1 to expand=0) combined with exponential backoff on HTTP 500/503.
    """

    BASE_URL = "https://rest.ensembl.org"

    def __init__(self, cache: Optional[SQLiteCache] = None, email: Optional[str] = None):
        super().__init__(name="Ensembl", cache=cache)
        self.email = ClientConfig.get_email(email)
        self.headers = ClientConfig.get_headers(self.email)
        self.rate_limiter = AsyncRateLimiter(requests_per_second=RateLimitConfig.ENSEMBL_RPS)

    async def _get_json(self, endpoint: str, params: Optional[Dict[str, Any]] = None, max_retries: int = 4) -> Optional[Any]:
        """Perform resilient GET request to Ensembl REST API with rate-limiting and backoff."""
        url = f"{self.BASE_URL}{endpoint}"
        import asyncio
        for attempt in range(max_retries):
            await self.rate_limiter.acquire()
            async with httpx.AsyncClient(timeout=RateLimitConfig.ENSEMBL_TIMEOUT_SEC) as client:
                try:
                    resp = await client.get(url, params=params, headers=self.headers)
                    if resp.status_code == 200:
                        return resp.json()
                    elif resp.status_code == 429:
                        logger.warning(f"Ensembl 429 rate limit hit on attempt {attempt + 1}, cooling down...")
                        await asyncio.sleep(1.0 * (attempt + 1))
                    elif resp.status_code in (500, 502, 503, 504):
                        logger.warning(f"Ensembl server error {resp.status_code} on attempt {attempt + 1} | url={url}")
                        await asyncio.sleep(0.8 * (attempt + 1))
                    elif resp.status_code in (404, 400):
                        return None
                    else:
                        logger.warning(f"Ensembl unexpected status {resp.status_code} on attempt {attempt + 1} | url={url}")
                except Exception as e:
                    logger.warning(f"Ensembl request attempt {attempt + 1} failed | url={url} error={e}")
                    if attempt == max_retries - 1:
                        return None
                    await asyncio.sleep(0.8 * (attempt + 1))
        return None

    def _parse_gene_data(self, data: Dict[str, Any]) -> GeneEntity:
        ensembl_id = data.get("id")
        symbol = data.get("display_name") or ensembl_id
        description = data.get("description")
        biotype = data.get("biotype")
        species = data.get("species", "homo_sapiens").replace("_", " ").capitalize()

        # Coordinates
        seq_region = data.get("seq_region_name", "")
        start = data.get("start")
        end = data.get("end")
        strand_num = data.get("strand")
        strand = "+" if strand_num == 1 else ("-" if strand_num == -1 else None)
        assembly = data.get("assembly_name", "GRCh38")

        location = GenomicLocation(
            chromosome=seq_region,
            start=start,
            end=end,
            strand=strand,
            assembly=assembly
        )

        # Transcripts parsing
        transcripts: List[TranscriptEntity] = []
        raw_transcripts = data.get("Transcript", [])
        for t in raw_transcripts:
            t_id = t.get("id")
            t_name = t.get("display_name")
            is_canon = bool(t.get("is_canonical", 0))
            t_biotype = t.get("biotype")
            t_len = t.get("length")
            t_start = t.get("start")
            t_end = t.get("end")
            
            # Translation / protein ID
            trans_obj = t.get("Translation")
            prot_id = trans_obj.get("id") if trans_obj else None

            # Exons
            exons: List[ExonEntity] = []
            for e in t.get("Exon", []):
                e_id = e.get("id")
                e_start = e.get("start")
                e_end = e.get("end")
                e_strand = "+" if e.get("strand") == 1 else ("-" if e.get("strand") == -1 else None)
                if e_id and e_start and e_end:
                    exons.append(ExonEntity(exon_id=e_id, start=e_start, end=e_end, strand=e_strand))

            transcripts.append(
                TranscriptEntity(
                    transcript_id=t_id,
                    name=t_name,
                    is_canonical=is_canon,
                    biotype=t_biotype,
                    length=t_len,
                    protein_id=prot_id,
                    start=t_start,
                    end=t_end,
                    exons=exons
                )
            )

        # Default taxon_id heuristic
        taxon_id = 9606 if "homo" in species.lower() else (10090 if "mus" in species.lower() else 0)

        # Sanitize description and extract HGNC ID if embedded e.g. [Source:HGNC Symbol;Acc:HGNC:11998]
        hgnc_id = None
        clean_name = description
        if description:
            import re
            m = re.search(r"Acc:(HGNC:\d+)", description)
            if m:
                hgnc_id = m.group(1)
            clean_name = re.sub(r"\s*\[Source:[^\]]+\]", "", description).strip() or description

        return GeneEntity(
            symbol=symbol,
            name=clean_name,
            taxon_id=taxon_id,
            species=species,
            hgnc_id=hgnc_id,
            ensembl_gene_id=ensembl_id,
            locus_type=biotype,
            location=location,
            transcripts=transcripts
        )

    async def fetch_by_id(self, ensembl_id: str, expand: bool = True) -> Optional[GeneEntity]:
        """Fetch gene model and transcript annotations by Ensembl Gene ID (e.g. ENSG00000141510)."""
        cache_key = f"id:{ensembl_id.upper()}:expand={expand}"
        cached = self.cache.get("ensembl", cache_key)
        if cached:
            return GeneEntity(**cached)

        endpoint = f"/lookup/id/{ensembl_id.upper()}"
        params = {"expand": "1" if expand else "0"}
        data = await self._get_json(endpoint, params=params)

        # Fallback: if expand=1 triggers server error 500, fetch compact model and enrich transcripts via symbol
        if not data and expand:
            logger.info(f"Ensembl expand=1 lookup failed for {ensembl_id}; attempting compact lookup with symbol fallback")
            data = await self._get_json(endpoint, params={"expand": "0"})
            if data:
                gene = self._parse_gene_data(data)
                # Attempt transcript enrichment via fetch_by_symbol if symbol is available
                if gene.symbol and gene.symbol != gene.ensembl_gene_id:
                    symbol_gene = await self.fetch_by_symbol(species=data.get("species", "homo_sapiens"), symbol=gene.symbol, expand=True)
                    if symbol_gene and symbol_gene.transcripts:
                        gene.transcripts = symbol_gene.transcripts
                self.cache.set("ensembl", cache_key, gene.model_dump())
                return gene

        if not data:
            return None

        gene = self._parse_gene_data(data)
        self.cache.set("ensembl", cache_key, gene.model_dump())
        return gene

    async def fetch_by_symbol(self, species: str, symbol: str, expand: bool = True) -> Optional[GeneEntity]:
        """Fetch gene model and coordinates by species and approved symbol (e.g. 'homo_sapiens', 'TP53')."""
        species_slug = species.lower().replace(" ", "_")
        cache_key = f"symbol:{species_slug}:{symbol.upper()}:expand={expand}"
        cached = self.cache.get("ensembl", cache_key)
        if cached:
            return GeneEntity(**cached)

        endpoint = f"/lookup/symbol/{species_slug}/{symbol.upper()}"
        params = {"expand": "1" if expand else "0"}
        data = await self._get_json(endpoint, params=params)
        if not data:
            return None

        gene = self._parse_gene_data(data)
        self.cache.set("ensembl", cache_key, gene.model_dump())
        return gene

    async def resolve_gene(self, query: str, taxon_id: int = 9606) -> Optional[Dict[str, Any]]:
        """Resolve a gene query string via Ensembl lookup by ID or symbol."""
        query_clean = query.strip()
        # 1. Direct Ensembl Gene ID lookup (e.g. ENSG... or ENSMUSG...)
        if query_clean.upper().startswith("ENS"):
            gene = await self.fetch_by_id(query_clean)
            if gene:
                return {"entity": gene, "rule": "ensembl_id_lookup", "confidence": 1.0}

        # 2. Symbol lookup
        species = "homo_sapiens" if taxon_id == 9606 else ("mus_musculus" if taxon_id == 10090 else "homo_sapiens")
        gene = await self.fetch_by_symbol(species=species, symbol=query_clean)
        if gene:
            return {"entity": gene, "rule": "ensembl_symbol_lookup", "confidence": 0.95}

        return None

    async def fetch_orthologs(
        self,
        gene_id_or_symbol: str,
        target_species: str = "mus_musculus",
        source_species: str = "homo_sapiens"
    ) -> List[OrthologEntity]:
        """Fetch orthologous genes across species via Ensembl Homology REST API."""
        source_slug = source_species.lower().replace(" ", "_")
        target_slug = target_species.lower().replace(" ", "_")
        query_clean = gene_id_or_symbol.strip()

        # If symbol provided instead of Ensembl Gene ID, resolve symbol first
        gene_id = query_clean
        if not query_clean.upper().startswith("ENS"):
            gene = await self.fetch_by_symbol(species=source_slug, symbol=query_clean, expand=False)
            if not gene or not gene.ensembl_gene_id:
                return []
            gene_id = gene.ensembl_gene_id

        cache_key = f"orthologs:{source_slug}:{gene_id.upper()}:{target_slug}"
        cached = self.cache.get("ensembl", cache_key)
        if cached and "items" in cached:
            return [OrthologEntity(**item) for item in cached["items"]]

        endpoint = f"/homology/id/{source_slug}/{gene_id.upper()}"
        params = {"target_species": target_slug, "type": "orthologues"}
        data = await self._get_json(endpoint, params=params)

        # Fallback: if /homology/id fails due to Ensembl server 500/503, try /homology/symbol
        if not data:
            # If original query was a symbol or if gene symbol can be fetched
            symbol = query_clean if not query_clean.upper().startswith("ENS") else None
            if not symbol:
                gene_obj = await self.fetch_by_id(gene_id, expand=False)
                if gene_obj and gene_obj.symbol and gene_obj.symbol != gene_id:
                    symbol = gene_obj.symbol
            if symbol:
                logger.info(f"Ensembl homology/id failed; attempting fallback via homology/symbol/{source_slug}/{symbol}")
                symbol_ep = f"/homology/symbol/{source_slug}/{symbol.upper()}"
                data = await self._get_json(symbol_ep, params=params)

        if not data:
            return []

        homologies = data.get("data", [{}])[0].get("homologies", [])
        orthologs: List[OrthologEntity] = []

        for h in homologies:
            target = h.get("target", {})
            target_id = target.get("id")
            if not target_id:
                continue

            orthologs.append(
                OrthologEntity(
                    source_gene_id=gene_id,
                    source_species=source_slug,
                    target_gene_id=target_id,
                    target_species=target.get("species", target_slug),
                    orthology_type=h.get("type", "ortholog"),
                    percent_identity=target.get("perc_id"),
                    target_protein_id=target.get("protein_id")
                )
            )

        self.cache.set("ensembl", cache_key, {"items": [o.model_dump() for o in orthologs]})
        return orthologs

    async def fetch_mgi_id(self, ensembl_gene_id: str) -> Optional[str]:
        """Fetch MGI ID for a mouse Ensembl gene via Ensembl xrefs."""
        gene_id = ensembl_gene_id.strip().upper()
        cache_key = f"xref:mgi:{gene_id}"
        cached = self.cache.get("ensembl", cache_key)
        if cached and "mgi_id" in cached:
            return cached["mgi_id"]

        endpoint = f"/xrefs/id/{gene_id}"
        params = {"external_db": "MGI"}
        data = await self._get_json(endpoint, params=params)
        if data:
            for item in data:
                prim_id = item.get("primary_id")
                if prim_id and "MGI:" in prim_id:
                    self.cache.set("ensembl", cache_key, {"mgi_id": prim_id})
                    return prim_id

        return None


class MGIAdapter(BaseBioAdapter):
    """Adapter for Mouse Genome Informatics (MGI) model organism data.
    Authoritative resource for laboratory mouse genetics, markers, and nomenclature.
    Integrates via Alliance of Genome Resources (AGR) API and Ensembl / UniProt cross-references.
    """

    ALLIANCE_API_URL = "https://www.alliancegenome.org/api/gene"

    def __init__(self, cache: Optional[SQLiteCache] = None, email: Optional[str] = None):
        super().__init__(name="MGI", cache=cache)
        self.email = ClientConfig.get_email(email)
        self.headers = ClientConfig.get_headers(self.email)
        self.rate_limiter = AsyncRateLimiter(requests_per_second=RateLimitConfig.MGI_RPS)

    async def fetch_by_mgi_id(self, mgi_id: str) -> Optional[GeneEntity]:
        """Fetch mouse gene details by primary MGI ID (e.g. 'MGI:98834')."""
        clean_id = mgi_id.strip()
        if not clean_id.upper().startswith("MGI:"):
            clean_id = f"MGI:{clean_id}"

        cache_key = f"mgi:{clean_id.upper()}"
        cached = self.cache.get("mgi", cache_key)
        if cached:
            return GeneEntity(**cached)

        url = f"{self.ALLIANCE_API_URL}/{clean_id}"
        data = None
        for attempt in range(3):
            await self.rate_limiter.acquire()
            async with httpx.AsyncClient(timeout=RateLimitConfig.MGI_TIMEOUT_SEC) as client:
                try:
                    resp = await client.get(url, headers=self.headers)
                    if resp.status_code == 200:
                        data = resp.json()
                        break
                    elif resp.status_code in (404, 400):
                        return None
                    elif resp.status_code == 429:
                        import asyncio
                        await asyncio.sleep(1.0)
                except Exception as e:
                    logger.warning(f"MGI AGR lookup attempt {attempt + 1} failed | id={clean_id} error={e}")
                    if attempt == 2:
                        return None
                    import asyncio
                    await asyncio.sleep(0.5 * (attempt + 1))
        else:
            return None

        gene_obj = data.get("gene", {})
        if not gene_obj:
            return None

        # Symbol & Name
        symbol = gene_obj.get("geneSymbol", {}).get("displayText") or clean_id
        full_name = gene_obj.get("geneFullName", {}).get("displayText")

        # Synonyms
        synonyms: List[str] = []
        for syn_item in gene_obj.get("geneSynonyms", []):
            syn_txt = syn_item.get("displayText")
            if syn_txt and syn_txt not in synonyms:
                synonyms.append(syn_txt)

        # Cross-references
        uniprot_ids: List[str] = []
        ncbi_gene_id = None
        for xref in gene_obj.get("crossReferences", []):
            ref_name = xref.get("name") or xref.get("displayName") or ""
            if ref_name.startswith("UniProtKB:"):
                acc = ref_name.replace("UniProtKB:", "").strip()
                if acc and acc not in uniprot_ids:
                    uniprot_ids.append(acc)
            elif ref_name.startswith("NCBI_Gene:") or ref_name.startswith("GeneID:"):
                ncbi_gene_id = ref_name.split(":")[-1].strip()

        # Genomic location
        location = None
        loc_assocs = gene_obj.get("geneGenomicLocationAssociations", [])
        if loc_assocs:
            first_loc = loc_assocs[0]
            start = first_loc.get("start")
            end = first_loc.get("end")
            strand = first_loc.get("strand")
            chr_name = first_loc.get("geneGenomicLocationAssociationObject", {}).get("name", "")
            location = GenomicLocation(
                chromosome=str(chr_name) if chr_name else "unknown",
                start=start,
                end=end,
                strand=strand,
                assembly="GRCm39"
            )

        gene_entity = GeneEntity(
            symbol=symbol,
            name=full_name,
            taxon_id=10090,
            species="Mus musculus",
            mgi_id=clean_id,
            ncbi_gene_id=ncbi_gene_id,
            uniprot_ids=uniprot_ids,
            synonyms=synonyms,
            location=location
        )

        self.cache.set("mgi", cache_key, gene_entity.model_dump())
        return gene_entity

    async def resolve_gene(self, query: str, taxon_id: int = 10090) -> Optional[Dict[str, Any]]:
        """Resolve a mouse gene query via MGI ID."""
        clean_q = query.strip()
        if clean_q.upper().startswith("MGI:"):
            gene = await self.fetch_by_mgi_id(clean_q)
            if gene:
                return {"entity": gene, "rule": "mgi_id_lookup", "confidence": 1.0}

        return None


class GeneOntologyAdapter(BaseBioAdapter):
    """Adapter for EMBL-EBI QuickGO REST API.
    Authoritative resource for Gene Ontology (GO) terms and functional annotations.
    """

    BASE_URL = "https://www.ebi.ac.uk/QuickGO/services"

    def __init__(self, cache: Optional[SQLiteCache] = None, email: Optional[str] = None):
        super().__init__(name="GeneOntology", cache=cache)
        self.email = ClientConfig.get_email(email)
        self.headers = ClientConfig.get_headers(self.email)
        self.rate_limiter = AsyncRateLimiter(requests_per_second=RateLimitConfig.QUICKGO_RPS)

    async def resolve_gene(self, query: str, taxon_id: int = 9606) -> Optional[Dict[str, Any]]:
        """QuickGO does not resolve gene symbols directly; delegating to primary resolver."""
        return None

    async def fetch_term(self, go_id: str) -> Optional[Dict[str, Any]]:
        """Fetch metadata, name, definition, and aspect for a specific GO ID."""
        clean_id = go_id.strip().upper()
        if not clean_id.startswith("GO:"):
            clean_id = f"GO:{clean_id}"

        cache_key = f"term:{clean_id}"
        cached = self.cache.get("go", cache_key)
        if cached:
            return cached

        url = f"{self.BASE_URL}/ontology/go/terms/{clean_id}"
        await self.rate_limiter.acquire()
        try:
            async with httpx.AsyncClient(timeout=RateLimitConfig.QUICKGO_TIMEOUT_SEC) as client:
                resp = await client.get(url, headers=self.headers)
                if resp.status_code != 200:
                    logger.warning("QuickGO term fetch failed | id=%s status=%d", clean_id, resp.status_code)
                    return None
                data = resp.json()
        except Exception as e:
            logger.error("QuickGO term fetch exception | id=%s error=%s", clean_id, str(e))
            return None

        results = data.get("results", [])
        if not results:
            return None

        item = results[0]
        term_data = {
            "go_id": item.get("id"),
            "name": item.get("name"),
            "aspect": item.get("aspect"),
            "is_obsolete": item.get("isObsolete", False),
            "definition": item.get("definition", {}).get("text") if isinstance(item.get("definition"), dict) else None,
            "synonyms": [s.get("name") for s in item.get("synonyms", []) if isinstance(s, dict) and s.get("name")]
        }

        self.cache.set("go", cache_key, term_data)
        return term_data

    async def fetch_annotations(
        self,
        gene_product_id: str,
        taxon_id: int = 9606,
        aspect: Optional[str] = None,
        limit: int = 50
    ) -> List[GOAnnotation]:
        """Fetch GO annotations associated with a UniProt accession or gene product ID."""
        clean_id = gene_product_id.strip()
        # Ensure ID format for QuickGO (e.g. UniProtKB:P04637 or just P04637)
        if clean_id.upper().startswith("UNIPROTKB:"):
            clean_id = clean_id.split(":", 1)[1]

        cache_key = f"annot:{clean_id}:tax{taxon_id}:asp{aspect}:lim{limit}"
        cached = self.cache.get("go", cache_key)
        if cached and isinstance(cached, list):
            return [GOAnnotation(**item) for item in cached if isinstance(item, dict)]

        params: Dict[str, Any] = {
            "geneProductId": clean_id,
            "limit": min(limit, 100),
        }
        if taxon_id:
            params["taxonId"] = taxon_id
        if aspect:
            params["aspect"] = aspect

        url = f"{self.BASE_URL}/annotation/search"
        data = None
        for attempt in range(3):
            await self.rate_limiter.acquire()
            try:
                async with httpx.AsyncClient(timeout=RateLimitConfig.QUICKGO_TIMEOUT_SEC) as client:
                    resp = await client.get(url, params=params, headers=self.headers)
                    if resp.status_code == 200:
                        data = resp.json()
                        break
                    logger.warning("QuickGO annotation search attempt %d failed | gene_product=%s status=%d", attempt + 1, clean_id, resp.status_code)
            except Exception as e:
                logger.warning("QuickGO annotation search attempt %d exception | gene_product=%s error=%s", attempt + 1, clean_id, str(e))
                if attempt == 2:
                    logger.error("QuickGO annotation search failed after 3 attempts | gene_product=%s", clean_id)
                    return []
                import asyncio
                await asyncio.sleep(0.5)

        if not data:
            return []

        annotations: List[GOAnnotation] = []
        raw_results = data.get("results", [])

        # Deduplicate annotations by (go_id, aspect, evidence_code)
        seen_keys = set()
        for res in raw_results:
            go_id = res.get("goId")
            if not go_id:
                continue

            go_aspect = res.get("goAspect") or aspect or "unknown"
            go_name = res.get("goName")
            go_evidence = res.get("goEvidence")
            eco_id = res.get("evidenceCode")
            assigned_by = res.get("assignedBy")

            dedup_key = (go_id, go_aspect, go_evidence)
            if dedup_key in seen_keys:
                continue
            seen_keys.add(dedup_key)

            annotations.append(
                GOAnnotation(
                    go_id=go_id,
                    name=go_name,
                    aspect=go_aspect,
                    evidence_code=go_evidence,
                    eco_id=eco_id,
                    assigned_by=assigned_by,
                    definition=None
                )
            )

        # Batch resolve term names if missing from annotation payload
        missing_names_ids = list({a.go_id for a in annotations if not a.name})
        if missing_names_ids:
            try:
                # Query in batches of up to 50 terms
                term_names_map: Dict[str, str] = {}
                chunk_size = 50
                for i in range(0, len(missing_names_ids), chunk_size):
                    chunk = missing_names_ids[i:i + chunk_size]
                    joined_ids = ",".join(chunk)
                    terms_url = f"{self.BASE_URL}/ontology/go/terms/{joined_ids}"
                    await self.rate_limiter.acquire()
                    async with httpx.AsyncClient(timeout=RateLimitConfig.QUICKGO_TIMEOUT_SEC) as client:
                        terms_resp = await client.get(terms_url, headers=self.headers)
                        if terms_resp.status_code == 200:
                            t_data = terms_resp.json()
                            for t_item in t_data.get("results", []):
                                tid = t_item.get("id")
                                tname = t_item.get("name")
                                if tid and tname:
                                    term_names_map[tid] = tname

                for a in annotations:
                    if not a.name and a.go_id in term_names_map:
                        a.name = term_names_map[a.go_id]
            except Exception as e:
                logger.debug("QuickGO term name batch enrichment skipped | error=%s", str(e))

        # Cache results as list of dicts
        self.cache.set("go", cache_key, [a.model_dump() for a in annotations])
        return annotations


class ReactomeAdapter(BaseBioAdapter):
    """Adapter for Reactome Content Service REST API.
    Authoritative knowledge base for biological pathways, molecular reactions, and processes.
    """

    BASE_URL = "https://reactome.org/ContentService"

    def __init__(self, cache: Optional[SQLiteCache] = None, email: Optional[str] = None):
        super().__init__(name="Reactome", cache=cache)
        self.email = ClientConfig.get_email(email)
        self.headers = ClientConfig.get_headers(self.email)
        self.rate_limiter = AsyncRateLimiter(requests_per_second=RateLimitConfig.REACTOME_RPS)

    async def resolve_gene(self, query: str, taxon_id: int = 9606) -> Optional[Dict[str, Any]]:
        """Reactome does not resolve primary gene symbols; delegating to primary resolver."""
        return None

    async def fetch_pathways_by_uniprot(
        self,
        uniprot_accession: str,
        species: str = "Homo sapiens",
        limit: int = 25
    ) -> List[PathwayEntity]:
        """Fetch biological pathways containing the given UniProt accession."""
        clean_acc = uniprot_accession.strip().upper()
        cache_key = f"uniprot:{clean_acc}:sp:{species}:lim:{limit}"
        cached = self.cache.get("reactome", cache_key)
        if cached and isinstance(cached, list):
            return [PathwayEntity(**item) for item in cached if isinstance(item, dict)]

        url = f"{self.BASE_URL}/data/mapping/UniProt/{clean_acc}/pathways"
        await self.rate_limiter.acquire()
        try:
            async with httpx.AsyncClient(timeout=RateLimitConfig.REACTOME_TIMEOUT_SEC) as client:
                resp = await client.get(url, headers=self.headers)
                if resp.status_code != 200:
                    logger.warning("Reactome pathway mapping failed | acc=%s status=%d", clean_acc, resp.status_code)
                    return []
                raw_data = resp.json()
        except Exception as e:
            logger.error("Reactome pathway mapping exception | acc=%s error=%s", clean_acc, str(e))
            return []

        if not isinstance(raw_data, list):
            return []

        pathways: List[PathwayEntity] = []
        seen_st_ids = set()

        for item in raw_data:
            if not isinstance(item, dict):
                continue
            st_id = item.get("stId")
            if not st_id or st_id in seen_st_ids:
                continue

            # Species filter if provided
            item_species = item.get("speciesName") or "Homo sapiens"
            if species and species.lower() not in item_species.lower():
                continue

            name = item.get("displayName") or (item.get("name", [""])[0] if item.get("name") else st_id)
            is_in_disease = item.get("isInDisease", False)
            url_link = f"https://reactome.org/PathwayBrowser/#/{st_id}"

            seen_st_ids.add(st_id)
            pathways.append(
                PathwayEntity(
                    st_id=st_id,
                    name=name,
                    species=item_species,
                    is_in_disease=is_in_disease,
                    url=url_link,
                    summary=None
                )
            )
            if len(pathways) >= limit:
                break

        self.cache.set("reactome", cache_key, [p.model_dump() for p in pathways])
        return pathways

    async def fetch_pathway_details(self, st_id: str) -> Optional[Dict[str, Any]]:
        """Fetch detailed information for a pathway by its Reactome stable ID."""
        clean_id = st_id.strip().upper()
        cache_key = f"pathway:{clean_id}"
        cached = self.cache.get("reactome", cache_key)
        if cached and isinstance(cached, dict):
            return cached

        url = f"{self.BASE_URL}/data/query/{clean_id}"
        await self.rate_limiter.acquire()
        try:
            async with httpx.AsyncClient(timeout=RateLimitConfig.REACTOME_TIMEOUT_SEC) as client:
                resp = await client.get(url, headers=self.headers)
                if resp.status_code != 200:
                    logger.warning("Reactome pathway query failed | st_id=%s status=%d", clean_id, resp.status_code)
                    return None
                data = resp.json()
        except Exception as e:
            logger.error("Reactome pathway query exception | st_id=%s error=%s", clean_id, str(e))
            return None

        if not isinstance(data, dict):
            return None

        summation_text = None
        summations = data.get("summation", [])
        if summations and isinstance(summations, list):
            first_sum = summations[0]
            if isinstance(first_sum, dict):
                summation_text = first_sum.get("text")

        details = {
            "st_id": data.get("stId", clean_id),
            "name": data.get("displayName") or clean_id,
            "species": data.get("speciesName", "Homo sapiens"),
            "is_in_disease": data.get("isInDisease", False),
            "release_date": data.get("releaseDate"),
            "url": f"https://reactome.org/PathwayBrowser/#/{clean_id}",
            "summary": summation_text,
            "has_diagram": data.get("hasDiagram", False),
            "has_ehld": data.get("hasEHLD", False)
        }

        self.cache.set("reactome", cache_key, details)
        return details


class MondoAdapter(BaseBioAdapter):
    """Adapter for MONDO Disease Ontology via EBI OLS4 API.
    Primary authority for canonical disease identifiers (MONDO:xxxxxxx), preferred names,
    synonyms, and cross-references (OMIM, Orphanet, DOID, UMLS).
    """

    BASE_URL = "https://www.ebi.ac.uk/ols4/api"

    def __init__(self, cache: Optional[SQLiteCache] = None, email: Optional[str] = None):
        super().__init__(name="MONDO", cache=cache)
        self.email = ClientConfig.get_email(email)
        self.headers = ClientConfig.get_headers(self.email)
        self.rate_limiter = AsyncRateLimiter(requests_per_second=RateLimitConfig.MONDO_RPS)

    async def resolve_gene(self, query: str, taxon_id: int = 9606) -> Optional[Dict[str, Any]]:
        """MONDO resolves disease entities rather than genes."""
        return None

    async def fetch_by_id(self, mondo_id: str) -> Optional[DiseaseEntity]:
        """Fetch canonical disease entity by MONDO ID (e.g. 'MONDO:0018875' or 'MONDO_0018875')."""
        clean_id = mondo_id.strip()
        if clean_id.startswith("MONDO:"):
            iri_id = clean_id.replace(":", "_")
        elif clean_id.startswith("MONDO_"):
            iri_id = clean_id
            clean_id = clean_id.replace("_", ":")
        else:
            return None

        cache_key = f"mondo:id:{clean_id.upper()}"
        cached = self.cache.get("mondo", cache_key)
        if cached:
            return DiseaseEntity(**cached)

        iri = f"http://purl.obolibrary.org/obo/{iri_id}"
        url = f"{self.BASE_URL}/ontologies/mondo/terms"
        params = {"iri": iri}

        await self.rate_limiter.acquire()
        try:
            async with httpx.AsyncClient(timeout=RateLimitConfig.MONDO_TIMEOUT_SEC) as client:
                resp = await client.get(url, params=params, headers=self.headers)
                if resp.status_code != 200:
                    return None
                data = resp.json()
        except Exception as e:
            logger.warning("Mondo term lookup exception | id=%s error=%s", clean_id, str(e))
            return None

        terms = data.get("_embedded", {}).get("terms", [])
        if not terms:
            return None

        term = terms[0]
        name = term.get("label") or clean_id
        descriptions = term.get("description", [])
        desc = descriptions[0] if descriptions else None
        synonyms = term.get("synonyms") or []
        
        # Cross references / dbxrefs
        xrefs = []
        annotation = term.get("annotation", {})
        dbxrefs = annotation.get("database_cross_reference", [])
        if isinstance(dbxrefs, list):
            xrefs = [str(x) for x in dbxrefs]

        disease = DiseaseEntity(
            mondo_id=clean_id,
            name=name,
            description=desc,
            synonyms=synonyms,
            cross_references=xrefs
        )
        self.cache.set("mondo", cache_key, disease.model_dump())
        return disease

    async def search_disease(self, query: str, limit: int = 5) -> List[DiseaseEntity]:
        """Search diseases by name, synonym, or keyword within MONDO."""
        clean_query = query.strip()
        if not clean_query:
            return []

        cache_key = f"mondo:search:{clean_query.lower()}:{limit}"
        cached = self.cache.get("mondo", cache_key)
        if cached and "items" in cached:
            return [DiseaseEntity(**item) for item in cached["items"]]

        url = f"{self.BASE_URL}/search"
        params = {
            "q": clean_query,
            "ontology": "mondo",
            "rows": limit,
            "queryFields": "label,synonym"
        }

        await self.rate_limiter.acquire()
        try:
            async with httpx.AsyncClient(timeout=RateLimitConfig.MONDO_TIMEOUT_SEC) as client:
                resp = await client.get(url, params=params, headers=self.headers)
                if resp.status_code != 200:
                    return []
                data = resp.json()
        except Exception as e:
            logger.warning("Mondo search exception | query=%s error=%s", clean_query, str(e))
            return []

        docs = data.get("response", {}).get("docs", [])
        results: List[DiseaseEntity] = []
        for doc in docs:
            short_form = doc.get("short_form", "")
            if not short_form.startswith("MONDO_") and not short_form.startswith("MONDO:"):
                continue
            mondo_id = short_form.replace("_", ":")
            name = doc.get("label") or mondo_id
            descriptions = doc.get("description", [])
            desc = descriptions[0] if descriptions else None
            synonyms = doc.get("synonym", []) or []

            results.append(DiseaseEntity(
                mondo_id=mondo_id,
                name=name,
                description=desc,
                synonyms=synonyms,
                cross_references=[]
            ))

        self.cache.set("mondo", cache_key, {"items": [r.model_dump() for r in results]})
        return results


class OpenTargetsAdapter(BaseBioAdapter):
    """Adapter for Open Targets Platform GraphQL API.
    Primary authority for target-disease association scores, clinical pipeline evidence,
    and genetic target tractability.
    """

    GRAPHQL_URL = "https://api.platform.opentargets.org/api/v4/graphql"

    def __init__(self, cache: Optional[SQLiteCache] = None, email: Optional[str] = None):
        super().__init__(name="OpenTargets", cache=cache)
        self.email = ClientConfig.get_email(email)
        self.headers = ClientConfig.get_headers(self.email)
        self.headers["Content-Type"] = "application/json"
        self.rate_limiter = AsyncRateLimiter(requests_per_second=RateLimitConfig.OPENTARGETS_RPS)

    async def resolve_gene(self, query: str, taxon_id: int = 9606) -> Optional[Dict[str, Any]]:
        """OpenTargets is queried via resolved Ensembl ID rather than as a primary gene symbol resolver."""
        return None

    async def fetch_target_diseases(
        self,
        ensembl_gene_id: str,
        symbol: Optional[str] = None,
        limit: int = 10
    ) -> TargetAssociationContext:
        """Fetch top disease associations for a target gene by Ensembl Gene ID."""
        clean_id = ensembl_gene_id.strip().upper()
        cache_key = f"opentargets:target:{clean_id}:{limit}"
        cached = self.cache.get("opentargets", cache_key)
        if cached:
            return TargetAssociationContext(**cached)

        query = """
        query TargetDiseases($ensemblId: String!, $size: Int!) {
          target(ensemblId: $ensemblId) {
            id
            approvedSymbol
            approvedName
            associatedDiseases(page: {size: $size, index: 0}) {
              count
              rows {
                score
                datatypeScores {
                  id
                  score
                }
                disease {
                  id
                  name
                }
              }
            }
          }
        }
        """

        payload = {
            "query": query,
            "variables": {
                "ensemblId": clean_id,
                "size": limit
            }
        }

        await self.rate_limiter.acquire()
        try:
            async with httpx.AsyncClient(timeout=RateLimitConfig.OPENTARGETS_TIMEOUT_SEC) as client:
                resp = await client.post(self.GRAPHQL_URL, json=payload, headers=self.headers)
                if resp.status_code != 200:
                    logger.warning("OpenTargets GraphQL error | status=%d body=%s", resp.status_code, resp.text[:200])
                    return TargetAssociationContext(query=clean_id, ensembl_gene_id=clean_id, symbol=symbol)
                res_data = resp.json()
        except Exception as e:
            logger.warning("OpenTargets request exception | target=%s error=%s", clean_id, str(e))
            return TargetAssociationContext(query=clean_id, ensembl_gene_id=clean_id, symbol=symbol)

        target_data = res_data.get("data", {}).get("target")
        if not target_data:
            return TargetAssociationContext(query=clean_id, ensembl_gene_id=clean_id, symbol=symbol)

        approved_symbol = target_data.get("approvedSymbol") or symbol
        assoc_block = target_data.get("associatedDiseases", {})
        total_count = assoc_block.get("count", 0)
        rows = assoc_block.get("rows", [])

        associations: List[TargetDiseaseAssociation] = []
        for r in rows:
            disease_info = r.get("disease", {})
            d_id = disease_info.get("id", "")
            d_name = disease_info.get("name", d_id)
            score = float(r.get("score", 0.0))

            dt_scores_dict: Dict[str, float] = {}
            for dt in r.get("datatypeScores", []):
                dt_id = dt.get("id")
                dt_val = dt.get("score")
                if dt_id and dt_val is not None:
                    dt_scores_dict[dt_id] = float(dt_val)

            associations.append(TargetDiseaseAssociation(
                disease_id=d_id.replace("_", ":") if d_id.startswith("MONDO_") else d_id,
                disease_name=d_name,
                score=round(score, 4),
                datatype_scores=dt_scores_dict if dt_scores_dict else None
            ))

        context = TargetAssociationContext(
            query=clean_id,
            ensembl_gene_id=clean_id,
            symbol=approved_symbol,
            total_associations=total_count,
            associations=associations
        )

        self.cache.set("opentargets", cache_key, context.model_dump())
        return context


class LiteratureAdapter(BaseBioAdapter):
    """Adapter for scientific publications via Europe PMC REST API.
    Provides authoritative citation metadata, PMID/PMCID/DOI cross-mapping,
    and open-access abstracts for biological and clinical research.
    """

    BASE_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest"

    def __init__(self, cache: Optional[SQLiteCache] = None, email: Optional[str] = None):
        super().__init__(name="EuropePMC", cache=cache)
        self.email = ClientConfig.get_email(email)
        self.headers = ClientConfig.get_headers(self.email)
        self.rate_limiter = AsyncRateLimiter(requests_per_second=RateLimitConfig.EUROPEPMC_RPS)

    async def resolve_gene(self, query: str, taxon_id: int = 9606) -> Optional[Dict[str, Any]]:
        """Literature adapter queries publications rather than gene symbol models."""
        return None

    def _parse_publication(self, item: Dict[str, Any]) -> PublicationEntity:
        """Parse Europe PMC raw result item into PublicationEntity."""
        pmid = item.get("pmid") or item.get("id")
        pmcid = item.get("pmcid")
        doi = item.get("doi")
        title = item.get("title", "").strip().rstrip(".")
        
        # Author parsing
        author_string = item.get("authorString", "")
        authors = [a.strip() for a in author_string.split(",") if a.strip()] if author_string else []

        journal_title = item.get("journalTitle")
        pub_year = None
        if "pubYear" in item and item["pubYear"]:
            try:
                pub_year = int(item["pubYear"])
            except (ValueError, TypeError):
                pass

        abstract = item.get("abstractText")
        cited_by = None
        if "citedByCount" in item and item["citedByCount"] is not None:
            try:
                cited_by = int(item["citedByCount"])
            except (ValueError, TypeError):
                pass

        url = f"https://europepmc.org/article/MED/{pmid}" if pmid else (
            f"https://doi.org/{doi}" if doi else None
        )

        return PublicationEntity(
            pmid=str(pmid) if pmid else None,
            pmcid=str(pmcid) if pmcid else None,
            doi=str(doi) if doi else None,
            title=title or "Untitled Publication",
            authors=authors,
            journal=journal_title,
            pub_year=pub_year,
            abstract_text=abstract,
            cited_by_count=cited_by,
            url=url
        )

    async def fetch_by_id(self, identifier: str) -> Optional[PublicationEntity]:
        """Fetch publication metadata by PMID, PMCID, or DOI."""
        clean_id = identifier.strip()
        if not clean_id:
            return None

        cache_key = f"europepmc:id:{clean_id.lower()}"
        cached = self.cache.get("europepmc", cache_key)
        if cached:
            return PublicationEntity(**cached)

        # Formulate query based on identifier pattern
        if clean_id.upper().startswith("PMC"):
            query_str = clean_id.upper()
        elif clean_id.startswith("10."):
            query_str = f"DOI:{clean_id}"
        elif clean_id.isdigit():
            query_str = f"EXT_ID:{clean_id} AND SRC:MED"
        else:
            query_str = f"EXT_ID:{clean_id}"


        url = f"{self.BASE_URL}/search"
        params = {
            "query": query_str,
            "format": "json",
            "pageSize": 1,
            "resultType": "core"
        }

        await self.rate_limiter.acquire()
        try:
            async with httpx.AsyncClient(timeout=RateLimitConfig.EUROPEPMC_TIMEOUT_SEC) as client:
                resp = await client.get(url, params=params, headers=self.headers)
                if resp.status_code != 200:
                    logger.warning("Europe PMC id query failed | query=%s status=%d", query_str, resp.status_code)
                    return None
                data = resp.json()
        except Exception as e:
            logger.warning("Europe PMC id query exception | query=%s error=%s", query_str, str(e))
            return None

        results = data.get("resultList", {}).get("result", [])
        if not results:
            return None

        pub = self._parse_publication(results[0])
        self.cache.set("europepmc", cache_key, pub.model_dump())
        return pub

    async def search_publications(self, query: str, limit: int = 5) -> LiteratureContext:
        """Search Europe PMC publications supporting an entity query."""
        clean_query = query.strip()
        if not clean_query:
            return LiteratureContext(query=query, total_hits=0, publications=[])

        cache_key = f"europepmc:search:{clean_query.lower()}:{limit}"
        cached = self.cache.get("europepmc", cache_key)
        if cached:
            return LiteratureContext(**cached)

        # Enrich biological search query (prefer peer-reviewed MEDLINE citations)
        if clean_query.isdigit():
            formatted_query = f"EXT_ID:{clean_query} AND SRC:MED"
        else:
            formatted_query = f"({clean_query}) AND (SRC:MED OR SRC:PMC)"

        url = f"{self.BASE_URL}/search"
        params = {
            "query": formatted_query,
            "format": "json",
            "pageSize": limit,
            "resultType": "core",
            "sort": "CITED desc"
        }

        await self.rate_limiter.acquire()
        try:
            async with httpx.AsyncClient(timeout=RateLimitConfig.EUROPEPMC_TIMEOUT_SEC) as client:
                resp = await client.get(url, params=params, headers=self.headers)
                if resp.status_code != 200:
                    logger.warning("Europe PMC search failed | query=%s status=%d", formatted_query, resp.status_code)
                    return LiteratureContext(query=clean_query, total_hits=0, publications=[])
                data = resp.json()
        except Exception as e:
            logger.warning("Europe PMC search exception | query=%s error=%s", formatted_query, str(e))
            return LiteratureContext(query=clean_query, total_hits=0, publications=[])

        total_count = int(data.get("hitCount", 0))
        raw_items = data.get("resultList", {}).get("result", [])
        publications = [self._parse_publication(item) for item in raw_items]

        context = LiteratureContext(
            query=clean_query,
            source="Europe PMC",
            total_hits=total_count,
            publications=publications
        )

        self.cache.set("europepmc", cache_key, context.model_dump())
        return context







