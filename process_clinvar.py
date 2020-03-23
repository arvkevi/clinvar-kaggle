import csv
import gzip
import re

import pandas as pd


def list_to_dict(l):
    """Convert list to dict."""
    return {k: v for k, v in (x.split("=") for x in l)}


fieldnames = [
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "AF_ESP",
    "AF_EXAC",
    "AF_TGP",
    "CLNDISDB",
    "CLNDISDBINCL",
    "CLNDN",
    "CLNDNINCL",
    "CLNHGVS",
    "CLNSIGINCL",
    "CLNVC",
    "CLNVI",
    "MC",
    "ORIGIN",
    "SSR",
    "CLASS",
    "Allele",
    "Consequence",
    "IMPACT",
    "SYMBOL",
    "Feature_type",
    "Feature",
    "BIOTYPE",
    "EXON",
    "INTRON",
    "cDNA_position",
    "CDS_position",
    "Protein_position",
    "Amino_acids",
    "Codons",
    "DISTANCE",
    "STRAND",
    "BAM_EDIT",
    "SIFT",
    "PolyPhen",
    "MOTIF_NAME",
    "MOTIF_POS",
    "HIGH_INF_POS",
    "MOTIF_SCORE_CHANGE",
    "LoFtool",
    "CADD_PHRED",
    "CADD_RAW",
    "BLOSUM62",
]

cv_columns = {}
with gzip.open("clinvar.vcf.gz", "rt") as f:
    for metaline in f:
        if metaline.startswith("##INFO"):
            colname = re.search("ID=(\w+),", metaline.strip("#\n"))
            coldesc = re.search(".*Description=(.*)>", metaline.strip("#\n"))
            cv_columns[colname.group(1)] = coldesc.group(1).strip('"')

# read clinvar vcf
cv_df = pd.read_csv(
    "clinvar.vcf.gz",
    sep="\t",
    comment="#",
    header=None,
    usecols=[0, 1, 2, 3, 4, 7],
    dtype={0: object},
)

# convert dictionaries to columns
cv_df = pd.concat(
    [
        cv_df.drop([7], axis=1),
        cv_df[7].str.split(";").apply(list_to_dict).apply(pd.Series),
    ],
    axis=1,
)
# rename columns
cv_df.rename(columns={0: "CHROM", 1: "POS", 2: "ID", 3: "REF", 4: "ALT"}, inplace=True)

# drop columns we know we won't need
cv_df = cv_df.drop(columns=["CHROM", "POS", "REF", "ALT"])

# assign classes
cv_df["CLASS"] = 0
cv_df.loc[cv_df["CLNSIGCONF"].notnull(), "CLASS"] = 1

# convert NaN to 0 where allele frequencies are null
cv_df[["AF_ESP", "AF_EXAC", "AF_TGP"]] = cv_df[["AF_ESP", "AF_EXAC", "AF_TGP"]].fillna(
    0
)

# select variants that have beeen submitted by multiple organizations.
cv_df = cv_df.loc[
    cv_df["CLNREVSTAT"].isin(
        [
            "criteria_provided,_multiple_submitters,_no_conflicts",
            "criteria_provided,_conflicting_interpretations",
        ]
    )
]

# Reduce the size of the dataset below
cv_df.drop(columns=["ALLELEID", "RS", "DBVARID"], inplace=True)
# drop columns that would reveal class
cv_df.drop(columns=["CLNSIG", "CLNSIGCONF", "CLNREVSTAT"], inplace=True)
# drop this redundant columns
cv_df.drop(columns=["CLNVCSO", "GENEINFO"], inplace=True)

# dictionary to map ID to clinvar annotations
clinvar_annotations = cv_df.set_index("ID")[
    [col for col in cv_df.columns if col in fieldnames]
].to_dict(orient="index")

# open the output file
outfile = "clinvar_conflicting.csv"
with open(outfile, "w") as fout:
    dw = csv.DictWriter(
        fout, delimiter=",", fieldnames=fieldnames, extrasaction="ignore"
    )
    dw.writeheader()
    # read the VEP-annotated vcf file line-by-line
    with gzip.open("vep/clinvar.annotated.vcf.gz", "rt") as f:
        for line in f:
            if line.startswith("##INFO=<ID=CSQ"):
                m = re.search(r'.*Format: (.*)">', line)
                cols = m.group(1).split("|")
                continue

            if line.startswith("#"):
                continue
            record = line.split("\t")
            (
                chromosome,
                position,
                clinvar_id,
                reference_base,
                alternate_base,
                qual,
                filter_,
                info,
            ) = record
            info_field = info.strip("\n").split(";")

            # to lookup in clivnar_annotaitons
            clinvar_id = int(clinvar_id)

            # only keep the variants that have been evaluated by multiple submitters
            if clinvar_id in clinvar_annotations:
                # initialize a dictionary to hold all the VEP annotation data
                annotation_data = {column: None for column in cols}
                annotation_data.update(clinvar_annotations[clinvar_id])
                # fields directly from the vcf
                annotation_data["CHROM"] = str(chromosome)
                annotation_data["POS"] = position
                annotation_data["REF"] = reference_base
                annotation_data["ALT"] = alternate_base

                for annotations in info_field:
                    column, value = annotations.split("=")

                    if column == "CSQ":
                        for csq_column, csq_value in zip(cols, value.split("|")):
                            annotation_data[csq_column] = csq_value
                        continue

                    annotation_data[column] = value
                dw.writerow(annotation_data)
