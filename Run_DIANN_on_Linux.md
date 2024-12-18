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

ls *.raw | parallel -j 12 '
start_time=$(date +%s)
echo "Processing {1} started"
docker run -i -t -v /home/user/raw:/data_input quay.io/biocontainers/thermorawfileparser:1.4.5--h05cac1d_1 ThermoRawFileParser.sh -i /data_input/{1}
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
export DIANN_dir="/usr/share/diann-1.9.2"
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
bash run_diann.sh --dir ~/dir/to/your/mzML/  --lib $DIANN_dir/report-lib.predicted.speclib --threads 100 --verbose 3 --out ~/output/diann.tsv --qvalue 0.01 --matrices --out-lib ./report-lib.parquet --gen-spec-lib --predictor --xic --fasta $DIANN_dir/uniprotkb_proteome_UP000005640_2024_11_22.fasta --fasta-search --min-fr-mz 200 --max-fr-mz 2000 --met-excision --min-pep-len 6 --max-pep-len 40 --min-pr-mz 380 --max-pr-mz 980 --min-pr-charge 2 --max-pr-charge 6 --cut K*,R* --missed-cleavages 1 --unimod4 --var-mods 2 --var-mod UniMod:35,15.994915,M --mass-acc 20 --mass-acc-ms1 15 --double-search --relaxed-prot-inf --rt-profiling --pg-level 1 -reanalyse

```

