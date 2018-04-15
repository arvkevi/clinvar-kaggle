Scripts and data used to prepare a [Kaggle dataset](https://www.kaggle.com/kevinarvai/clinvar-conflicting).

**Generate dataset using ClinVar .vcf:**  
`python process_clinvar.py` will generate the file `clinvar_conflicting.csv`.

**Generate dataset using ClinVar .vcf w/ VEP annotations:**  
`python process_clinvar.py --vep` will generate a version of the file `clinvar_conflicting.csv` with [vep annotations](https://useast.ensembl.org/Tools/VEP). 
**Note:**  
[vawk](https://github.com/cc2qe/vawk) is required to run the `vep` version.

Check out the [notebook](https://github.com/arvkevi/clinvar-kaggle/blob/master/clinvar-conflicting-eda.ipynb) to see some exploratory data analysis.

## Problem Statement

The objective is to predict whether a ClinVar variant will have **conflicting classifications**.

*Conflicting classifications are when two of any of the following three classification categories are present for one variant, two submissions of one category is not considered conflicting.*

1. Likely Benign or Benign
2. VUS
3. Likely Pathogenic or Pathogenic

The `CLASS` feature in `clinvar_conflicting.csv` is a binary representation of whether or not a variant has conflicting classifications where `0` represents consistent classifications and `1` represents conflicting classifications.

Since this problem only relates to variants with multiple classifications, I removed all variants from the original ClinVar vcf which were only had one submission.

![](https://github.com/arvkevi/clinvar-kaggle/blob/master/clinvar-class-fig.png)

## Background

[ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) is a public resource containing annotations about human genetic variants. These variants are classified on a spectrum between benign, likely benign, uncertain significance, likely pathogenic, and pathogenic. Variants that have conflicting classifications (defined above) can cause confusion when clinicians or researchers try to interpret whether the variant has an impact on the disease of a given patient.  

I'm exploring ideas for applying machine learning to genomics. I'm hoping this project will encourage others to think about the additional feature engineering that's probably necessary to confidently assess the objective. There could be benefit to identifying *single submission* variants that may yet to have assigned a **conflicting classification**.

## VEP annotations

Ensembl's [Variant Effect Predictor (VEP)](http://grch37.ensembl.org/Homo_sapiens/Tools/VEP) was used to annotate the original ClinVar `.vcf`. It provides additional information about variants that can serve as features for the dataset.
#### Step 1:
Download and rename the annotated `.vcf` as `clinvar.annotated.vcf`

#### Step 2:
To make parsing the annotated `.vcf` easier, [vawk](https://github.com/cc2qe/vawk) was used to export the `CSQ` field using the following command:
```vawk --header '{ print $3, I$CSQ }' clinvar.annotated.vcf >clinar.annotated.csq.vcf```

#### Step 3:
Create the new dataset with vep annotations.
```python process_clinvar.py --vep```
