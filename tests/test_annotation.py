from refmasker.annotation import FeatureIndex, detect_format, load_annotation

NCBI_GFF = """##gff-version 3
chr1\tRefSeq\tregion\t1\t30000\t.\t+\t.\tID=chr1:1..30000
chr1\tRefSeq\tgene\t100\t400\t.\t+\t.\tID=gene-A_RS1;Name=dnaA;gene=dnaA;locus_tag=A_RS1
chr1\tProtein Homology\tCDS\t100\t400\t.\t+\t0\tID=cds-1;Parent=gene-A_RS1;product=replication initiator%2C DnaA;locus_tag=A_RS1
chr1\tRefSeq\tgene\t5000\t6000\t.\t-\t.\tID=gene-A_RS2;Name=A_RS2;locus_tag=A_RS2
chr1\tProtein Homology\tCDS\t5000\t6000\t.\t-\t0\tID=cds-2;Parent=gene-A_RS2;product=IS3 family transposase;locus_tag=A_RS2
chr1\tRefSeq\tgene\t7000\t7500\t.\t+\t.\tID=gene-A_RS3;Name=A_RS3;locus_tag=A_RS3;pseudo=true
##FASTA
>chr1
ACGT
"""

ENSEMBL_GFF = """##gff-version 3
NC_1\tena\tgene\t1\t1524\t.\t+\t.\tID=gene:Rv0001;Name=dnaA;biotype=protein_coding;description=Chromosomal replication initiator protein DnaA;gene_id=Rv0001
NC_1\tena\tmRNA\t1\t1524\t.\t+\t.\tID=transcript:Rv0001;Parent=gene:Rv0001
NC_1\tena\tCDS\t1\t1524\t.\t+\t0\tID=CDS:Rv0001;Parent=transcript:Rv0001;protein_id=Rv0001
"""


def test_ncbi_gff(tmp_path):
    p = tmp_path / "a.gff"
    p.write_text(NCBI_GFF)
    assert detect_format(p) == "gff3"
    feats = load_annotation(p)
    cds = [f for f in feats if f.ftype == "CDS"]
    assert cds[0].label == "A_RS1(dnaA)"
    assert cds[0].product == "replication initiator, DnaA"
    assert (cds[0].start, cds[0].end) == (99, 400)
    assert [f.is_mobile for f in cds] == [False, True]
    assert any(f.ftype == "pseudogene" for f in feats)


def test_ensembl_gff_names_come_from_parent_gene(tmp_path):
    p = tmp_path / "a.gff"
    p.write_text(ENSEMBL_GFF)
    (cds,) = load_annotation(p)
    assert cds.locus_tag == "Rv0001"
    assert cds.gene == "dnaA"
    assert "DnaA" in cds.product


def test_genbank(tmp_path):
    from Bio.Seq import Seq
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqIO import write
    from Bio.SeqRecord import SeqRecord

    rec = SeqRecord(Seq("A" * 2000), id="chr1.1", name="chr1", description="test")
    rec.annotations["molecule_type"] = "DNA"
    rec.features = [
        SeqFeature(
            FeatureLocation(10, 400, strand=1),
            type="CDS",
            qualifiers={"locus_tag": ["T1"], "gene": ["abc"], "product": ["ABC protein"]},
        ),
        SeqFeature(
            FeatureLocation(500, 1500, strand=-1),
            type="mobile_element",
            qualifiers={"mobile_element_type": ["insertion sequence:IS256"]},
        ),
    ]
    p = tmp_path / "a.gbk"
    write([rec], p, "genbank")
    assert detect_format(p) == "genbank"
    feats = load_annotation(p)
    assert feats[0].seqid == "chr1.1"
    assert feats[0].label == "T1(abc)"
    assert feats[1].is_mobile


def test_feature_index(tmp_path):
    p = tmp_path / "a.gff"
    p.write_text(NCBI_GFF)
    idx = FeatureIndex(load_annotation(p))
    hits = idx.overlapping("chr1", 350, 5100)
    assert [f.locus_tag for f in hits] == ["A_RS1", "A_RS2"]
    assert idx.overlapping("chr1", 400, 4999) == []
    assert idx.overlapping("other", 0, 10) == []
