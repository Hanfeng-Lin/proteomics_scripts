# proteomics_scripts
Data analysis for global proteomics DIA data

Two languages provided, python (.ipynb) and R (.R)

Python version relies on student t-test and use Benjamini-Hochberg to correct for multiple t-test FDR. Straight and forward.
R version use limma model, a very popular package for analyzing microarray and RNA-seq data. LIMMA stands for “linear models for microarray data”. It tends to give more aggressive adj P value than simply BH method.
