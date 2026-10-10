"""Tests for BioContext Entity Resolution and Schemas."""

import pytest
from biocontext.resolver import EntityResolver
from biocontext.base import SQLiteCache


@pytest.fixture
def temp_cache(tmp_path):
    db_file = tmp_path / "test_cache.db"
    return SQLiteCache(db_path=str(db_file))


@pytest.fixture
def resolver(temp_cache):
    return EntityResolver(cache=temp_cache)


@pytest.mark.asyncio
async def test_resolve_exact_symbol(resolver):
    result = await resolver.resolve("TP53")
    assert result.match_status == "exact"
    assert result.confidence_score == 1.0
    assert result.resolved_entity is not None
    assert result.resolved_entity.symbol == "TP53"
    assert result.resolved_entity.hgnc_id == "HGNC:11998"
    assert result.resolved_entity.ncbi_gene_id == "7157"
    assert "P04637" in result.resolved_entity.uniprot_ids


@pytest.mark.asyncio
async def test_resolve_alias_symbol(resolver):
    # HER2 is an alias for ERBB2
    result = await resolver.resolve("HER2")
    assert result.match_status == "alias"
    assert result.confidence_score >= 0.8
    assert result.resolved_entity is not None
    assert result.resolved_entity.symbol == "ERBB2"
    assert result.resolved_entity.hgnc_id == "HGNC:3430"


@pytest.mark.asyncio
async def test_uniprot_accession_lookup(resolver):
    # P04637 is UniProt accession for TP53
    result = await resolver.resolve("P04637")
    assert result.match_status == "exact"
    assert result.resolved_entity is not None
    assert result.resolved_entity.symbol == "TP53"


@pytest.mark.asyncio
async def test_ncbi_entrez_id_lookup(resolver):
    # 7157 is Entrez Gene ID for TP53
    result = await resolver.resolve("7157")
    assert result.match_status == "exact"
    assert result.resolved_entity is not None
    assert result.resolved_entity.symbol == "TP53"
    assert result.resolved_entity.ncbi_gene_id == "7157"


@pytest.mark.asyncio
async def test_unresolved_query(resolver):
    result = await resolver.resolve("NONEXISTENT_GENE_123456")
    assert result.match_status == "unresolved"
    assert result.confidence_score == 0.0
    assert result.resolved_entity is None


@pytest.mark.asyncio
async def test_disambiguation_with_matching_chromosome(resolver):
    from biocontext.schemas import ResolutionContext
    # TP53 is on chromosome 17
    ctx = ResolutionContext(chromosome="17")
    result = await resolver.resolve("TP53", context=ctx)
    assert result.match_status == "exact"
    assert result.confidence_score == 1.0
    rules = [r.rule for r in result.match_reasons]
    assert "chromosome_match" in rules


@pytest.mark.asyncio
async def test_disambiguation_with_mismatching_chromosome(resolver):
    from biocontext.schemas import ResolutionContext
    # TP53 is on chr17, but context asserts chromosome 2 (penalty applied)
    ctx = ResolutionContext(chromosome="2")
    result = await resolver.resolve("TP53", context=ctx)
    # 1.0 - 0.40 penalty = 0.60
    assert result.confidence_score == 0.60
    rules = [r.rule for r in result.match_reasons]
    assert "chromosome_mismatch" in rules


@pytest.mark.asyncio
async def test_ensembl_id_resolution(resolver):
    # ENSG00000141510 is Ensembl ID for TP53
    result = await resolver.resolve("ENSG00000141510")
    assert result.match_status == "exact"
    assert result.confidence_score == 1.0
    assert result.resolved_entity is not None
    entity = result.resolved_entity
    assert entity.symbol == "TP53"
    assert entity.ensembl_gene_id == "ENSG00000141510"
    assert entity.hgnc_id == "HGNC:11998"
    assert entity.ncbi_gene_id == "7157"
    assert "P04637" in entity.uniprot_ids
    assert "p53" in entity.synonyms or "LFS1" in entity.synonyms
    assert "[Source:" not in entity.name
    assert entity.location is not None
    assert entity.location.chromosome == "17"
    assert len(entity.transcripts) > 0


@pytest.mark.asyncio
async def test_ensembl_transcripts_retrieval(resolver):
    gene = await resolver.ensembl.fetch_by_symbol(species="homo_sapiens", symbol="TP53")
    assert gene is not None
    assert gene.ensembl_gene_id == "ENSG00000141510"
    assert gene.location.start is not None
    assert gene.location.end is not None
    assert gene.location.strand == "-"
    canonical_list = [t for t in gene.transcripts if t.is_canonical]
    assert len(canonical_list) >= 1
    canon = canonical_list[0]
    assert canon.transcript_id.startswith("ENST")
    assert len(canon.exons) > 0


@pytest.mark.asyncio
async def test_ensembl_ortholog_retrieval(resolver):
    # Test human TP53 (ENSG00000141510) -> mouse ortholog (Trp53 / ENSMUSG00000059552)
    orthologs = await resolver.ensembl.fetch_orthologs(
        gene_id_or_symbol="ENSG00000141510",
        target_species="mus_musculus",
        source_species="homo_sapiens"
    )
    assert len(orthologs) > 0
    mouse_ortholog = orthologs[0]
    assert mouse_ortholog.target_gene_id == "ENSMUSG00000059552"
    assert mouse_ortholog.target_species == "mus_musculus"
    assert mouse_ortholog.orthology_type == "ortholog_one2one"
    assert mouse_ortholog.percent_identity is not None
    assert mouse_ortholog.percent_identity > 70.0


@pytest.mark.asyncio
async def test_fuzzy_resolution_typo(resolver):
    # Test typo "TP54" resolves approximately to "TP53"
    result = await resolver.resolve("TP54")
    assert result.match_status == "fuzzy"
    assert result.confidence_score >= 0.70
    assert result.resolved_entity is not None
    assert result.resolved_entity.symbol == "TP53"
    rules = [r.rule for r in result.match_reasons]
    assert "fuzzy_levenshtein_match" in rules


@pytest.mark.asyncio
async def test_mgi_adapter_direct_lookup(resolver):
    # Test MGI:98834 directly via MGIAdapter
    gene = await resolver.mgi.fetch_by_mgi_id("MGI:98834")
    assert gene is not None
    assert gene.symbol == "Trp53"
    assert gene.taxon_id == 10090
    assert gene.species == "Mus musculus"
    assert gene.mgi_id == "MGI:98834"
    assert "p53" in [s.lower() for s in gene.synonyms]
    assert gene.location is not None
    assert gene.location.chromosome == "11"
    assert gene.location.assembly == "GRCm39"


@pytest.mark.asyncio
async def test_resolve_mgi_id(resolver):
    # Test resolving query starting with MGI:
    res = await resolver.resolve("MGI:98834", taxon_id=10090)
    assert res is not None
    assert res.match_status == "exact"
    assert res.confidence_score == 1.0
    assert res.resolved_entity is not None
    assert res.resolved_entity.symbol == "Trp53"
    assert res.resolved_entity.mgi_id == "MGI:98834"
    assert res.resolved_entity.taxon_id == 10090
