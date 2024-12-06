Prepare single-cell BCR/TCR data from the 10x Genomics 5’ solution for TRIBAL analysis.


## create env
```bash
conda env create -f tribal.yml 
conda activate tribal
```


## change R path
To avoid conflict with your local R installation you might have, you would want to run the following command to set R PATH after activated the environment.

```bash
export PATH=$HOME/miniconda3/envs/TRIBAL/bin:$PATH
```


## Parameters
To change parameters used in the file, modify them in `config.yaml`

```bash
	1.	species: Target species for the analysis (only “mouse” or "human" are supported for now).
	2.	data_dir: Path to the main data directory.
	3.	blastdb_dir: Path to the BLAST database directory.
	4.	igblastdb_dir: Path to the IgBLAST database directory.
	5.	germline_dir: Path to the germline sequences directory.
	6.	dandelion_dir: Directory for Dandelion-related files.
	7.	germline_vdj: Path template for V(D)J germline sequences.
	8.	v_germline_file: Path to the V germline FASTA file.
	9.	blastdb_bcr_c: Path to the BCR constant region database file.
	10.	contig_file: Path to the filtered contig file for Dandelion.
	11.	filtered_feature_h5: Name of the HDF5 file with feature data.
	12.	output_columns: List of attributes to include in the output.
	13.	mouse_isotypes: List of mouse antibody isotypes.
	14.	human_isotypes: List of human antibody isotypes.
	15.	tribal_input_table: Name of the TRIBAL input table.
	16.	tribal_root_sequences: Name of the TRIBAL root sequences file.
```

## Run the script 

```bash
python preprocess.py
```