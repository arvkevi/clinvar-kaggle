Scripts and data used to prepare a [Kaggle dataset](https://www.kaggle.com/kevinarvai/clinvar-conflicting).

`python process_clinvar.py` will generate the file `clinvar_conflicting.csv`. A version of the file is included in the repo for convenience.

Check out the [notebook](https://github.com/arvkevi/clinvar-kaggle/blob/master/clinvar-conflicting-eda.ipynb) to see some exploratory data analysis.

## Problem Statement

The objective is to predict whether a ClinVar variant will have **conflicting classifications**.

*Conflicting classifications are when two of any of the following three categories are present for one variant, two submissions of one category is not considered conflicting*
This is already handled in the `.csv` file.
1. Likely Benign or Benign
2. VUS
3. Likely Pathogenic or Pathogenic

The `CLASS` feature is the binary representation of whether or not a variant has conflicting classifications where `0` represents consistent classifications and `1` represents conflicting classifications.

Since this problem only relates to variants with multiple classifications, I removed all variants from the original ClinVar vcf which were only had one submission.

## Background

[ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) is a public resource containing annotations about human genetic variations. These *variants* are classified on a spectrum between benign, likely benign, uncertain significance, likely pathogenic, and pathogenic. Variants that have conflicting classifications (defined above) can cause confusion when clinicians or researchers try to interpret whether the variant has an impact on the disease of a given patient.  

I'm exploring ideas for applying machine learning to genomics. I'm hoping this project will encourage others to think about the additional feature engineering that's probably necessary to confidently assess the objective. There could be benefit to identifying *single submission* variants that may yet to have assigned a **conflicting classification**.
