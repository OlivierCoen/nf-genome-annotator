#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
from pathlib import Path

import logging
import polars as pl
import pandas as pd

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

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

COMMON_ATTRIBUTE_KEYS = [
    'ID',
    'Name',
    'Alias',
    'Parent',
    'Target',
    'Gap',
    'Derives_from',
    'Note',
    'Dbxref',
    'Ontology_term',
    'Is_circular',
]

ATTRIBUTES_TO_PRESERVE = [
    'ID',
    'Name',
    'Parent'
]



# ── GFF3 helpers ──────────────────────────────────────────────────────────────

def parse_gff3(file: Path) -> pl.LazyFrame:
    return pl.scan_csv(
        file, 
        separator="\t", 
        comment_prefix="#", 
        has_header=False, 
        new_columns=GFF_COLUMNS
    )

def parse_attrs(s: str) -> dict:
    attrs = {}
    for field in s.strip().split(";"):
        field = field.strip()
        if not field:
            continue
        if "=" in field:
            k, v = field.split("=", 1)
            attrs[k] = v
        else:
            attrs[field] = ""
    return attrs

def format_attrs(attrs):
    return ";".join(
        f"{k}={v}" if v != "" else k
        for k, v in attrs.items()
    )

def read_gff3_header(path):
    """Return header lines from GFF3"""
    header = []
    with open(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line.startswith("#"):
                header.append(line)
            else:
                break
    return header


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--annot", 
        dest="annot_file",
        required=True, 
        help="Annotation file in GFF3"
    )
    parser.add_argument(
        "--iprscan", 
        dest="iprscan_file",
        required=True, 
        help="InterProScan output GFF3"
    )
    parser.add_argument(
        "--out", 
        dest="outfile",
        required=True, 
        help="Output GFF3"
    )
    return parser.parse_args()


def get_nb_rows(lf: pl.LazyFrame) -> int:
    return lf.select(pl.len()).collect().item(0, 0)


def get_unique_transcript_attribute_keys(annot_lf: pl.LazyFrame) -> list[str]:
    return (
        annot_lf
        .filter(pl.col("feature") == "transcript")
        .filter(pl.col("attribute") != "")
        .select(
            pl.col("attribute")
            .str.split(";")
            .explode()
            .str.split("=")
            .list.to_struct(upper_bound=1)
            .struct[0]
        )
        .unique()
        .collect()
        .to_series()
        .to_list()
    )


def complement_interproscan_attributes(ipr_lf: pl.LazyFrame) -> pl.LazyFrame:
    return (
        ipr_lf
        .with_columns(
            (pl.col("attribute") + ";seqname=" + pl.col("seqname")).alias("completed_attribute")
        )
        .with_columns(
            (pl.col("completed_attribute") + ";source=" + pl.col("source")).alias("completed_attribute")
        )
        .with_columns(
            pl.col("completed_attribute").str.split(";").list.eval(pl.element().str.strip_chars()).alias("attribute_list")
        )
    )


def get_name_and_aliases(name_set: set) -> dict[str, str | list[str]] | None:
    if len(name_set) == 0:
        return None
    # taking the shortest name as the real name, and the rest as aliases
    min_name_length = min(len(name) for name in name_set)
    real_name = next(name for name in name_set if len(name) == min_name_length)
    aliases = [name for name in name_set if len(name) != min_name_length]
    return {"name": real_name, "aliases": aliases}


def parse_attributes_as_dicts(ipr_lf: pl.LazyFrame) -> pd.DataFrame:
    # parsing each list of attributes (one list per original row) into a dictionary
    parsed_attributes = [
        {
            k: v 
            for s in lst if "=" in s 
            for k, v in [s.split("=", 1)]
        }
        for lst in ipr_lf.select("attribute_list").collect().to_series().to_list()
    ]
    # making a dataframe out of it: empty keys will be NaN
    return pd.DataFrame(parsed_attributes)


def unite_multiple_sets(*sets) -> set:
    valid_sets = [s for s in sets if isinstance(s, set)]
    if not valid_sets:
        return set()
    return set().union(*valid_sets)


def get_grouped_set_dataframe(parsed_attr_df: pd.DataFrame, xref_col: str) -> pd.DataFrame:
    s = parsed_attr_df[pd.col(xref_col).notna()].groupby("seqname")[xref_col].apply(set)
    df = pd.DataFrame(s).reset_index()
    df["seqname"] = df["seqname"].astype(str)
    return df


def get_names_and_aliases(parsed_attr_df: pd.DataFrame) -> pd.DataFrame:
    pfam_names = parsed_attr_df[pd.col("source") == "Pfam"].groupby("seqname")["Name"].apply(set)
    seqname_to_names_and_aliases = pfam_names.apply(lambda x: get_name_and_aliases(x))
    seqname_to_names_and_aliases_records = [
        {"seqname": seqname, "Name": d["name"], "Alias": d["aliases"]}
        for seqname, d in seqname_to_names_and_aliases.to_dict().items()
    ]
    seqname_names_aliases_df = pd.DataFrame.from_records(seqname_to_names_and_aliases_records)
    seqname_names_aliases_df["Alias"] = seqname_names_aliases_df["Alias"].apply(','.join)
    return seqname_names_aliases_df


def get_all_dbxrefs(parsed_attr_df: pd.DataFrame):
    # Dbxref (IPR entries, database accessions)
    dbxrefs_df = get_grouped_set_dataframe(parsed_attr_df, "Dbxref")
    # Also store the hit accessions itself (e.g. PF00001, TIGR00001)
    hit_accessions_df = get_grouped_set_dataframe(parsed_attr_df, "hit_accession")
    hit_accession_aliases_df = get_grouped_set_dataframe(parsed_attr_df, "hit_accession_alias")
    
    all_dbxrefs_df = (
        dbxrefs_df
            .merge(hit_accessions_df, how="outer", on="seqname")
            .merge(hit_accession_aliases_df, how="outer", on="seqname")
    )
    all_dbxrefs_df["all"] = all_dbxrefs_df.apply(
        lambda row: unite_multiple_sets(row["Dbxref"], row["hit_accession"], row["hit_accession_alias"]), 
        axis=1
    )
    all_dbxrefs_df["all"] = all_dbxrefs_df["all"].apply(','.join)
    # add dbxrefs column
    return all_dbxrefs_df[["seqname", "all"]].rename(columns={"all": "Dbxref"})


def get_go_terms(parsed_attr_df: pd.DataFrame):
    go_terms_df = get_grouped_set_dataframe(parsed_attr_df, "Ontology_term")
    go_terms_df["Ontology_term"] = go_terms_df["Ontology_term"].apply(','.join)
    return go_terms_df


def merge_new_attributes(transcript_attr_lf: pl.LazyFrame, formated_ipr_attr_df: pl.DataFrame) -> pl.LazyFrame:
    merged_lf = transcript_attr_lf.join(
        formated_ipr_attr_df.lazy(), 
        left_on="ID",
        right_on="seqname", 
        how="left",
        suffix="_ipr",
    )

    attribute_cols = [
        col for col in merged_lf.collect_schema().names()
        if col != "index" and not col.endswith("_ipr")
    ]
    
    # merging original values and InterProScan values
    for col in attribute_cols:
        ipr_col = f'{col}_ipr'
        if ipr_col in attribute_cols:
            if col in ATTRIBUTES_TO_PRESERVE:
                logger.info(f"Merging values for {col}...")
                merged_lf = (
                    merged_lf
                    .with_columns(
                        pl.concat_list([pl.col(col), pl.col(ipr_col)]).alias(col)
                    )
                    .with_columns(
                        pl.col(col).list.join(",")
                    )
                    .drop(ipr_col)
                )
            else:
                logger.info(f"Dropping {ipr_col}...")
                merged_lf = merged_lf.drop(ipr_col)

    # making list of attribute columns in order of COMMON_ATTRIBUTE_KEYS followed by any remaining columns
    final_attribute_cols = [
        col for col in COMMON_ATTRIBUTE_KEYS if col in attribute_cols
    ] + [
        col for col in attribute_cols if col not in COMMON_ATTRIBUTE_KEYS
    ]
    
    # formating new attribute column from individual columns
    merged_lf = merged_lf.with_columns(pl.lit("").alias("attribute"))
    for col in final_attribute_cols:
        merged_lf = (
            merged_lf
            .with_columns(
                pl.when(pl.col(col).is_not_null())
                .then((pl.col("attribute") + f";{col}=" + pl.col(col)))
                .otherwise(pl.col("attribute"))
                .alias("attribute")
            )
            .drop(col)
        )

    # strip leading semicolon from attribute column
    merged_lf = merged_lf.with_columns(
        pl.col("attribute").str.strip_prefix(";").alias("attribute")
    )

    return merged_lf


def add_new_attributes_to_annotation(annot_lf: pl.LazyFrame, merged_lf: pl.LazyFrame) -> pl.LazyFrame:
    return (
        annot_lf
        .join(merged_lf, how="left", on="index", suffix="_new")
        .with_columns(
            pl.when(pl.col("attribute_new").is_not_null())
            .then(pl.col("attribute_new"))
            .otherwise(pl.col("attribute"))
            .alias("attribute")
        )
        .drop(["index", "attribute_new"])
    )
# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    
    args = parse_args()

    annot_header = read_gff3_header(args.annot_file)
    annot_lf = parse_gff3(args.annot_file)
    logger.info(f"Loaded {get_nb_rows(annot_lf)} structural features")
    ipr_lf = parse_gff3(args.iprscan_file)
    logger.info(f"Loaded {get_nb_rows(ipr_lf)} InterProScan features")

    ####################################################################
    # PARSING ANNOTATION
    ####################################################################
    
    annot_lf = annot_lf.with_row_index("index")

    # putting transcript attributes into separate columns
    # eg:
    #┌───────┬──────────┬──────────┬────────────┬───┬───────┬────────────────────┬───────┬────────┐
    #│ index ┆ seqname  ┆ source   ┆ feature    ┆ … ┆ frame ┆ attribute          ┆ ID    ┆ Parent │
    #│ ---   ┆ ---      ┆ ---      ┆ ---        ┆   ┆ ---   ┆ ---                ┆ ---   ┆ ---    │
    #│ u32   ┆ str      ┆ str      ┆ str        ┆   ┆ str   ┆ str                ┆ str   ┆ str    │
    #╞═══════╪══════════╪══════════╪════════════╪═══╪═══════╪════════════════════╪═══════╪════════╡
    #│ 1     ┆ ntLink_1 ┆ AUGUSTUS ┆ transcript ┆ … ┆ .     ┆ ID=g1.t2;Parent=g1 ┆ g1.t2 ┆ g1     │
    #│ 45    ┆ ntLink_1 ┆ AUGUSTUS ┆ transcript ┆ … ┆ .     ┆ ID=g1.t1;Parent=g1 ┆ g1.t1 ┆ g1     │
    #│ 89    ┆ ntLink_1 ┆ AUGUSTUS ┆ transcript ┆ … ┆ .     ┆ ID=g1.t3;Parent=g1 ┆ g1.t3 ┆ g1     │
    #│ 137   ┆ ntLink_1 ┆ AUGUSTUS ┆ transcript ┆ … ┆ .     ┆ ID=g2.t1;Parent=g2 ┆ g2.t1 ┆ g2     │
    #│ 155   ┆ ntLink_1 ┆ AUGUSTUS ┆ transcript ┆ … ┆ .     ┆ ID=g3.t1;Parent=g3 ┆ g3.t1 ┆ g3     │
    #└───────┴──────────┴──────────┴────────────┴───┴───────┴────────────────────┴───────┴────────┘

    unique_attribute_keys = get_unique_transcript_attribute_keys(annot_lf)

    if "ID" not in unique_attribute_keys:
        raise ValueError("ID attribute not found in transcript attributes")
    
    transcript_attr_lf = (
        annot_lf
        .filter(pl.col("feature") == "transcript")
        .with_columns(
            pl.col("attribute").str.extract(rf"{attr}=(.*?)(?:;|$)").alias(attr)
            for attr in unique_attribute_keys
        )
        .select(["index"] + unique_attribute_keys)
    )

    ####################################################################
    # PARSING INTERPROSCAN OUTPUT
    ####################################################################
    
    # removing weird stuff at the end of the file and excluded sources
    EXCLUDED_SOURCES = ["MobiDB-lite"]
    ipr_lf = (
        ipr_lf
        .filter(pl.col("source").is_not_null())
        .filter(~pl.col("source").is_in(EXCLUDED_SOURCES))
    )

    # aranging columns
    TO_REPLACE = {"PROSITE profiles": "PROSITE", "PROSITE patterns": "PROSITE"}
    ipr_lf = ipr_lf.with_columns(pl.col("source").replace(TO_REPLACE))

    ipr_lf = complement_interproscan_attributes(ipr_lf)

    # parsing each list of attributes (one list per original row) into a dictionary
    ipr_attr_df = parse_attributes_as_dicts(ipr_lf)

    # pre-format the hit accession
    ipr_attr_df["hit_accession"] = ipr_attr_df["source"] + ":" + ipr_attr_df["Name"]
    ipr_attr_df["hit_accession_alias"] = ipr_attr_df["source"] + ":" + ipr_attr_df["Alias"]

    EXPECTED_COLUMNS = ["Name", "Ontology_term", "Alias", "Dbxref"]
    for col in EXPECTED_COLUMNS:
        if col not in ipr_attr_df.columns:
            ipr_attr_df[col] = pd.NA

    # Names and aliases
    logger.info("Processing Names and Aliases...")
    seqname_names_aliases_df = get_names_and_aliases(ipr_attr_df)
    
    logger.info("Processing Dbxrefs...")
    all_dbxrefs_df = get_all_dbxrefs(ipr_attr_df)
   
    logger.info("Processing Ontology_term (GO terms)...")
    goterms_df = get_go_terms(ipr_attr_df)

    formated_ipr_attr_df = pl.from_pandas(
        seqname_names_aliases_df
        .merge(all_dbxrefs_df, on="seqname", how="left")
        .merge(goterms_df, on="seqname", how="left")
    )

    del ipr_attr_df
    del seqname_names_aliases_df
    del all_dbxrefs_df
    del goterms_df

    ####################################################################
    # MERGING WITH STRUCTURAL GFF3
    ####################################################################

    logger.info("Merging original attributes with InterProScan attributes...")
    merged_lf = merge_new_attributes(transcript_attr_lf, formated_ipr_attr_df)

    logger.info("Merging original annotation with table containing new attributes...")
    final_annot_lf = add_new_attributes_to_annotation(annot_lf, merged_lf)

    logger.info(f"Writing final annotation to {args.outfile}")
    # write final annotation to output file
    with open(args.outfile, "w") as fout:
        fout.writelines(annot_header)
        final_annot_lf.sink_csv(fout, separator="\t", include_header=False)

    logger.info("Done")


if __name__ == "__main__":
    main()