Scripts and data used to prepare for a [Kaggle dataset](). 

The objective is to predict whether ClinVar variants will have **conflicting classifications**.

*of note: Conflicting classifications are when two of any of the following three categories are present for one variant*
1. Likely Benign or Benign
2. VUS
3. Likely Pathogenic or Pathogenic

The `CLASS` feature is the binary representation of whether or not a variant has conflicting classifications where `0` represents consistent classifications and `1` represents conflicting classifications.
