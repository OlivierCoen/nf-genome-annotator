import marimo

__generated_with = "0.23.1"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo

    return


@app.cell
def _():
    import polars as pl
    from collections import Counter
    import matplotlib.pyplot as plt
    from pathlib import Path

    return Counter, Path, pl, plt


@app.cell
def _():
    from goatools.obo_parser import GODag
    from goatools.go_enrichment import GOEnrichmentStudy
    from goatools.obo_parser import GODag
    from goatools.mapslim import mapslim

    return GODag, GOEnrichmentStudy, mapslim


@app.cell
def _(GODag):
    # wget http://purl.obolibrary.org/obo/go/go.obo
    # wget http://current.geneontology.org/ontology/subsets/goslim_generic.obo
    #godag = GODag("go-basic.obo")
    godag = GODag("go.obo")
    slim_dag  = GODag("goslim_generic.obo")
    return godag, slim_dag


@app.cell
def _():
    GFF_COLUMNS = [
        'seqname',
        'source',
        'feature',
        'start',
        'end',
        'score',
        'strand',
        'frame',
        'attribute'
    ]
    return (GFF_COLUMNS,)


@app.cell
def _(GFF_COLUMNS, Path, pl):
    def parse_gff3(file: Path) -> pl.LazyFrame:
        return pl.scan_csv(
            file, 
            separator="\t", 
            comment_prefix="#", 
            has_header=False, 
            new_columns=GFF_COLUMNS
        )


    def get_gene_to_go(gff3_lf: pl.LazyFrame):
        attribute_list = (
            gff3_lf
            .filter(pl.col("feature") == "transcript")
            .filter(pl.col("attribute") != "")
            .select(
                pl.col("attribute")
                .str.split(";").list.eval(pl.element().str.strip_chars())
            )
            .collect()
            .to_series().to_list()

        )
        """
        gives:
        [
            [
                "ID=g1226.t1"
                "Parent=g1226"
                "Ontology_term=GO:0003674,GO:0003675"
            ]
        ]
        """
        parsed_attributes = [
            {
                k: v 
                for s in lst if "=" in s 
                for k, v in [s.split("=", 1)]
            }
            for lst in attribute_list
        ]
        """
        gives:
        [
            {
                "ID": "g1226.t1"
                "Parent": "g1226"
                "Ontology_term": "GO:0003674,GO:0003675"
            }
        ]
        """
        return {
            d["Parent"]: set(d["Ontology_term"].split(",")) 
            for d in parsed_attributes 
            if "ID" in d and "Ontology_term" in d
        }

    return get_gene_to_go, parse_gff3


@app.cell
def _():
    gff3 = "/home/olivier/repositories/nf-genome-annotator/data/anasatum.final.gff3"
    gene_ids_file = "/home/olivier/repositories/nf-variant-calling/data/sex_location/distance/strict_with_poolseq/genes_of_interest.txt"
    #gene_ids_file = "/home/olivier/repositories/nf-variant-calling/data/resistance_on_minimum_depth/distance/strict/dominant/genes_of_interest.txt"
    #gene_ids_file = "/home/olivier/repositories/nf-variant-calling/data/resistance_on_minimum_depth/distance/strict/recessive/genes_of_interest.txt"
    return gene_ids_file, gff3


@app.cell
def _(gene_ids_file, get_gene_to_go, gff3, parse_gff3):
    gff3_lf = parse_gff3(gff3)
    gene2go = get_gene_to_go(gff3_lf)
    with open(gene_ids_file) as f:
        study_genes = set(line.strip() for line in f)
    return gene2go, study_genes


@app.cell
def _(GOEnrichmentStudy, gene2go, godag, study_genes):
    population_genes = set(gene2go.keys())

    # Run enrichment (Fisher's exact test + FDR correction)
    goeaobj = GOEnrichmentStudy(
        population_genes,
        gene2go,          # your dict from Step 1
        godag,
        propagate_counts=True,
        alpha=0.05,
        methods=["fdr_bh"]   # Benjamini-Hochberg correction
    )

    results = goeaobj.run_study(study_genes)
    return (results,)


@app.cell
def _(results):
    go2pvalue_uncorrected = {res.GO: res.p_uncorrected for res in results}
    go2pvalue_fdr_bh = {res.GO: res.p_fdr_bh for res in results}
    return go2pvalue_fdr_bh, go2pvalue_uncorrected


@app.cell
def _(go2pvalue_fdr_bh):
    { k: v for k, v in go2pvalue_fdr_bh.items() if v < 1}
    return


@app.cell
def _(gene_ids_file, go2pvalue_uncorrected):
    outfile = gene_ids_file.replace(".txt", ".go2p_value.txt")
    with open(outfile, 'w') as fout:
        for go, pvalue in go2pvalue_uncorrected.items():
            fout.write(f"{go}\t{pvalue}\n")
    return


@app.cell
def _(gene2go, godag, mapslim, slim_dag):
    slim_terms = set(slim_dag.keys())  # all GOslim term IDs

    # Map each gene's GO terms to their slim equivalents
    gene2goslim = {}

    for gene, go_terms in gene2go.items():
        slim_mapped = set()
        for go_term in go_terms:
            if go_term in godag:
                direct_ancestors, all_ancestors = mapslim(go_term, godag, slim_dag)
                slim_mapped.update(direct_ancestors)
        gene2goslim[gene] = slim_mapped
    return (gene2goslim,)


@app.cell
def _(gene2goslim):
    gene2goslim
    return


@app.cell
def _(Counter, gene2goslim, plt, study_genes):
    # Count slim terms across your study genes
    slim_counts = Counter()
    for study_gene in study_genes:
        if study_gene in gene2goslim:
            slim_counts.update(gene2goslim[study_gene])

    # Plot
    labels, values = zip(*slim_counts.most_common(15))
    plt.barh(labels, values)
    plt.xlabel("Number of genes")
    plt.title("GOslim distribution in genes of interest")
    plt.tight_layout()
    plt.savefig("goslim_distribution.png")
    return


@app.cell
def _():
    0
    return


if __name__ == "__main__":
    app.run()
