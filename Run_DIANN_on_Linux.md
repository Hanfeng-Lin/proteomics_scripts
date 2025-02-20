### Upload .raw files to server

Use FileZilla to transfer all the .raw data to your directory. (e.g.  $HOME/MS_data/)



### Convert .raw to .mzML format

#### Load ThermoRawFileParser Docker images

Reference: https://github.com/compomics/ThermoRawFileParser?tab=readme-ov-file

##### One single raw file

In the following command, change `/home/user/raw:/data_input` to your .raw file directory `$HOME/MS_data:/data_input`

```bash
docker run -i -t -v /home/user/raw:/data_input quay.io/biocontainers/thermorawfileparser:1.4.5--h05cac1d_1  ThermoRawFileParser.sh -d /data_input

```

##### Multiple raw files

In the following command, change `/home/user/raw:/data_input` to your .raw file directory `$HOME/MS_data:/data_input`

```bash
cd $HOME/MS_data

ls *.raw | parallel -j 20 '
start_time=$(date +%s)
echo "Processing {1} started"
docker run -i -t -v $HOME/20250206_Jurkat_vav1/plate1:/data_input quay.io/biocontainers/thermorawfileparser:1.4.5--h05cac1d_1 ThermoRawFileParser.sh -i /data_input/{1}
end_time=$(date +%s)
elapsed_seconds=$((end_time - start_time))
elapsed_minutes=$(echo "scale=2; $elapsed_seconds / 60" | bc)
echo "Processing {1} finished in ${elapsed_minutes} minutes"
'
```



Alternatively, use local binary file instead of docker (Recommended):

```bash
cd $HOME/MS_data

ls *.raw | parallel -j 12 '
start_time=$(date +%s)
echo "Processing {1} started"
mono /usr/share/ThermoRawFileParser1.4.5/ThermoRawFileParser.exe -i=./{1} -f=2
end_time=$(date +%s)
elapsed_seconds=$((end_time - start_time))
elapsed_minutes=$(echo "scale=2; $elapsed_seconds / 60" | bc)
echo "Processing {1} finished in ${elapsed_minutes} minutes"
'

```



### Run DIA-NN 1.9.2 on Linux

#### For the first time users

​	You need to add `diann-1.9.2` to your system variable`$PATH`  

​	**Open your `.bash_profile`**: You can open the `.bash_profile` file using a text editor. For example, using `nano`:

```bash
nano ~/.bash_profile
```

​	**Add the row to the file**: Once the file is open, add the following line at the end of the file:

```bash
PATH="/usr/share/diann-1.9.2:$PATH"
export DIANN-dir="/usr/share/diann-1.9.2"
```

​	**Save the changes**:

- In `nano`, press `Ctrl + O` to save the file, then press `Enter` to confirm.
- Press `Ctrl + X` to exit the editor.

​	**Apply the changes**: After modifying the `.bash_profile`, you need to apply the changes to the current session. You can do this by running:

```bash
source ~/.bash_profile
```



#### Global proteomics DIANN command example

In the following command, change `--dir` to the directory where `.mzML` files are located in, also `-out` to the output directory.

`--fasta`: Change if you are searching something other than human proteome. Meanwhile, also change `--lib` with a new name for your proteome

`--reanalyse`: enable MBR.

```bash
# DL-based Spectra Prediction from Fasta (in silico digestion)
bash run_diann.sh  --threads 100 --verbose 3  --qvalue 0.01 --matrices --out-lib ./report-lib.parquet --gen-spec-lib --predictor --xic --fasta $DIANN_dir/uniprotkb_proteome_UP000005640_Homo_sapiens_reviewed_2025_01_07.fasta --fasta-search --min-fr-mz 200 --max-fr-mz 2000 --met-excision --min-pep-len 6 --max-pep-len 40 --min-pr-mz 380 --max-pr-mz 980 --min-pr-charge 2 --max-pr-charge 6 --cut K*,R* --missed-cleavages 1 --unimod4 --var-mods 2 --var-mod UniMod:35,15.994915,M --mass-acc 20 --mass-acc-ms1 15 --double-search --relaxed-prot-inf --rt-profiling --pg-level 1 --reanalyse --dir ~/20250206_Jurkat_YDS/mzml/ --out ~/20250206_Jurkat_YDS/mzml/diann.tsv


```

```bash
# use predicted lib, no gen-spec-lib and predictor
bash run_diann.sh --threads 100 --verbose 3 --qvalue 0.01 --matrices --out-lib ./report-lib.parquet --xic --fasta $DIANN_dir/uniprotkb_proteome_UP000005640_Homo_sapiens_reviewed_2025_01_07.fasta --lib $DIANN_dir/report-lib.predicted.speclib --min-fr-mz 200 --max-fr-mz 2000 --met-excision --min-pep-len 6 --max-pep-len 40 --min-pr-mz 380 --max-pr-mz 980 --min-pr-charge 2 --max-pr-charge 6 --cut K*,R* --missed-cleavages 1 --unimod4 --var-mods 2 --var-mod UniMod:35,15.994915,M --mass-acc 20 --mass-acc-ms1 15 --double-search --relaxed-prot-inf --rt-profiling --pg-level 1 --reanalyse --dir ~/20250206_Jurkat_YDS/mzml/ --out ~/20250206_Jurkat_YDS/mzml/diann.tsv

```





### Post-processing (experimental)

```R
library(diann)
library(DIAgui)

df <- diann_load("homo_sapiens_reviewed_highaccuracy_NGT20-12.tsv")
df$Stripped.Sequence <- paste(df$Protein.Group, df$Genes, df$Stripped.Sequence, sep = "_")

peptides.maxlfq <- diann_maxlfq(df[df$Q.Value <= 0.01 & df$PG.Q.Value <= 0.01,], group.header="Stripped.Sequence", id.header = "Precursor.Id", quantity.header = "Precursor.Normalised")

# Matrix to DataFrame
peptides.maxlfq <- as.data.frame(peptides.maxlfq)
# Preserve row names by adding them as a new column
peptides.maxlfq$Uniprot_Gene_Peptide <- rownames(peptides.maxlfq)
Uniprot_Gene_Peptide <- strsplit(peptides.maxlfq$Uniprot_Gene_Peptide, "_")
peptides.maxlfq$Protein.Group <- sapply(Uniprot_Gene_Peptide, `[`, 1)
peptides.maxlfq$Genes <- sapply(Uniprot_Gene_Peptide, `[`, 2)
peptides.maxlfq$peptide <- sapply(Uniprot_Gene_Peptide, `[`, 3)


peptides.maxlfq.imputed <- DIAgui::imputationDIA(peptides.maxlfq, transformation = "log2", method = "MinProb") # log2 transformation is necessary to avoid negative imputation

write.table(peptides.maxlfq.imputed, file = "peptides_maxlfq_imputed.tsv", 
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Then something else to maxLFQ for protein level...
library("iq")

process_wide_format(input_filename = "peptides_maxlfq_imputed.tsv",
                    output_filename = "protein_maxlfq_imputed.tsv",
                    id_column="Protein.Group",
                    quant_columns= 1:8,
                    data_in_log_space = TRUE,
                    annotation_columns = c("Genes"),
                    method = "maxLFQ")

# load protein_maxlfq_imputed.tsv
protein.maxlfq.imputed <- read.table("protein_maxlfq_imputed.tsv", sep="\t", header=TRUE, 
                                     stringsAsFactors=FALSE)
protein.maxlfq.imputed <- as.data.frame(lapply(protein.maxlfq.imputed, function(x) {
    if (is.numeric(x)) {
        return(2^x)
    } else {
        return(x)
    }
}))
peptides.maxlfq.imputed <- as.data.frame(lapply(peptides.maxlfq.imputed, function(x) {
    if (is.numeric(x)) {
        return(2^x)
    } else {
        return(x)
    }
}))
write.table(peptides.maxlfq.imputed, file = "test_maxlfq_imputed.pr_matrix.tsv", 
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(protein.maxlfq.imputed, file = "test_maxlfq_imputed.pg_matrix.tsv", 
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
```

or, use fragpipe_analyst to impute on peptide level.

```R
library("iq")

process_wide_format(input_filename = "NGT20-12_pr_noseq.FragpipeAnalyst_Imputed_matrix.tsv",
                    output_filename = "NGT20-12_pr_noseq_maxLFQ.FragpipeAnalyst_Imputed_matrix.tsv",
                    id_column="ProteinID",
                    quant_columns= 2:9,
                    data_in_log_space = TRUE,
                    annotation_columns = NULL,
                    method = "maxLFQ")

# This leads to the imputation ~1.8 fold log2FC using Perseus-type imputation. No-go.
```

