# proteomics_scripts
Data analysis for global proteomics DIA data

Two languages provided, python (.ipynb) and R (.R)

Python version relies on student t-test and use Benjamini-Hochberg to correct for multiple t-test FDR. Imputation is optional but optimized for degrader-related experiment.

### Key Imputation Logic Steps

1. **Control Group Handling**:
    - **If missing values in control group >= 50%**: Discard the protein.
    - **If missing values in control group < 50%**:
        - Impute missing values using a **normal distribution** defined by the **mean** and **standard deviation (std)** of available control group values for that protein.

2. **Treated Group Handling**:
    - **If missing values in treated group < 50%**:
        - Impute missing values using the **mean** and **std** of available treated group values for that protein.
    - **If missing values in treated group >= 50%**:
        - **Case 1**: All treated group values are missing:
            - Check the peptide table to see if the protein has more than 3 peptides.
            - If yes, perform peptide-level statistical significance testing using a **t-test** between treated and control groups for the protein's peptides.
            - If the **median p-value < 0.05** or (all p-values are `NaN` but foldchange contains valid value, meaning n=1 t-test failed), impute treated values using a **uniform distribution** centered around the 1st percentile of the treated group.
        - **Case 2**: Some treated group values are available:
            - If peptide count > 5, impute missing treated values using the **mean** of the treated group for the protein and a **normalized standard deviation (CV)** derived from the treated group.




R version use limma model, a very popular package for analyzing microarray and RNA-seq data. 
LIMMA stands for “linear models for microarray data”. It tends to give more aggressive adj P value than simply BH method.
