# %%
# Import necessary libraries
import pandas as pd 
import numpy as np
import json
import pickle
import yaml
from scipy.io import mmread
import scipy.io
import gzip

# %%
def consensus(sequences):
    """
    Generate a consensus sequence from a list of sequences.
    Consensus is determined by selecting the most frequent nucleotide at each position.
    """
    if len(sequences) <= 1:  # Return the sequence itself if there's only one sequence
        return sequences[0]
    
    max_len = max(len(seq) for seq in sequences)  # Find the length of the longest sequence
    
    # Initialize a DataFrame to count occurrences of each nucleotide ('A', 'T', 'C', 'G', 'N') at each position
    cons_df = pd.DataFrame(0, index=range(max_len), columns=['A', 'T', 'C', 'G', 'N'])
    
    # Count nucleotide occurrences for each sequence
    for seq in sequences:
        for i, nt in enumerate(seq):
            if nt in cons_df.columns:
                cons_df.loc[i, nt] += 1

    # Determine the most frequent nucleotide at each position
    return ''.join(cons_df.idxmax(axis=1))

# %%
# Define the base directory for the dataset
base = './PARSE_Bcells'

# Load various input files
all_genes = pd.read_csv(f'{base}/all_genes.csv')  # Gene information
metadata = pd.read_csv(f'{base}/cell_metadata.csv')  # Cell metadata
clonotype = pd.read_csv(f'{base}/clonotype_frequency.tsv', sep='\t')  # Clonotype frequency data
annotation = pd.read_csv(f'{base}/bcr_annotation_airr.tsv', sep='\t')  # BCR annotation data
barcodes = pd.read_csv(f'{base}/barcode_report.tsv', sep='\t')  # Barcode data

# %%
barcode_rm_multiplet = barcodes[barcodes['isMultiplet'] == 0]  # Remove multiplets
annotation_rm_multiplet = annotation[annotation['cell_barcode'].isin(barcode_rm_multiplet['Barcode'])]  # Filter annotations for non-multiplets
annotation_full_length = annotation[annotation['full_length'] == 1]  # Filter for full-length annotations

# Filter barcodes with full-length data for heavy and light chains
barcodes_full_length = barcodes[
    (barcodes['IGH_full_length'] == 1) & 
    ((barcodes['IGK_full_length'] == 1) | (barcodes['IGL_full_length'] == 1))
]

# Filter full-length annotations for non-multiplet barcodes
annotation_rm_multiplet_full_length = annotation_rm_multiplet[annotation_rm_multiplet['full_length'] == 1]

# Fill missing values in v_call, d_call, and j_call with 'NA'
annotation_rm_multiplet_full_length.loc[:, ['v_call', 'd_call', 'j_call']] = annotation_rm_multiplet_full_length[['v_call', 'd_call', 'j_call']].fillna('NA')

# %%
# Group annotations by cell barcode and aggregate locus and alleles as lists or concatenated strings
annotation_group = annotation_rm_multiplet_full_length.groupby('cell_barcode').agg(
    {
        'locus': list,
        'v_call': lambda x: ','.join(map(str, sorted(x))),
        'd_call': lambda x: ','.join(map(str, sorted(x))),
        'j_call': lambda x: ','.join(map(str, sorted(x))),
    }
)

# %%
# Select groups with exactly two loci
annotation_group_2_locus = annotation_group[annotation_group['locus'].str.len() == 2]

# Group by V, D, and J alleles to identify clonotypes
clonotype_by_allele = annotation_group_2_locus.groupby(['v_call', 'd_call', 'j_call']).apply(
    lambda x: pd.Series({'count': len(x), 'indices': list(x.index)})
).reset_index()

# Filter for clonotypes with at least two associated cells
clonotype_by_allele_more = clonotype_by_allele[clonotype_by_allele['count'] >= 2]

# %%
import logging

# Configure logging to capture errors
logging.basicConfig(level=logging.ERROR, format='%(asctime)s - %(levelname)s - %(message)s')

# Initialize lists to store output data
records = []
clono_records = []
clonotype_map = {}

# Iterate through each row in the clonotype DataFrame
for i, t in clonotype_by_allele_more.iterrows():
    try:
        consensus_H = []  # Consensus sequences for heavy chains
        consensus_L = []  # Consensus sequences for light chains

        # Define the clonotype for the current row
        clonotype = t['v_call'].replace(',', '_') + ',' + t['d_call'].replace(',', '_') + ',' + t['j_call'].replace(',', '_')
        clonotype_map[i] = clonotype

        # Process each cell associated with the current clonotype
        for cell in t['indices']:
            try:
                # Filter annotation data for the current cell barcode
                temp = annotation_rm_multiplet_full_length[annotation_rm_multiplet_full_length['cell_barcode'] == cell]

                # Extract heavy and light chain sequences and alleles
                try:
                    IGH = temp.loc[temp['locus'] == 'IGH', 'sequence'].values[0]
                    IGL = temp.loc[temp['locus'].isin(['IGL', 'IGK']), 'sequence'].values[0]
                    H_v = temp.loc[temp['locus'] == 'IGH', 'v_call'].values[0]
                    L_v = temp.loc[temp['locus'].isin(['IGL', 'IGK']), 'v_call'].values[0]
                except IndexError as e:
                    logging.error(f"Error extracting sequences or v_call for cell {cell}: {e}")
                    continue

                # Extract heavy chain isotype
                try:
                    IGH_C = barcode_rm_multiplet.loc[barcode_rm_multiplet['Barcode'] == cell, 'IGH_C'].values[0]
                except IndexError as e:
                    logging.error(f"Error extracting IGH_C for cell {cell}: {e}")
                    IGH_C = None  # Default to None if missing

                # Extract germline alignments
                try:
                    H_germ = temp.loc[temp['locus'] == 'IGH', 'germline_alignment'].values[0]
                    L_germ = temp.loc[temp['locus'].isin(['IGL', 'IGK']), 'germline_alignment'].values[0]
                except IndexError as e:
                    logging.error(f"Error extracting germline alignment for cell {cell}: {e}")
                    continue

                # Append germline alignments to consensus lists
                consensus_H.append(H_germ)
                consensus_L.append(L_germ)

                # Create a record for the current cell
                records.append({
                    'cellid': cell,
                    'clonotype': i,
                    'heavy_chain_isotype': IGH_C,
                    'heavy_chain_seq': IGH,
                    'heavy_chain_v_allele': H_v,
                    'light_chain_seq': IGL,
                    'light_chain_v_allele': L_v
                })

            except Exception as cell_error:
                logging.error(f"Error processing cell {cell} in clonotype {i}: {cell_error}")
                continue

        # Create a record for the clonotype's consensus sequences
        try:
            clono_records.append({
                'clonotype': i,
                'heavy_chain_root': consensus(consensus_H),
                'light_chain_root': consensus(consensus_L)
            })
        except Exception as consensus_error:
            logging.error(f"Error creating consensus for clonotype {i}: {consensus_error}")
            continue

    except Exception as row_error:
        logging.error(f"Error processing row {i}: {row_error}")
        continue

# %%
# Convert records to DataFrames and save as CSV files
seq_data = pd.DataFrame.from_records(records)
root_data = pd.DataFrame.from_records(clono_records)

# Save the data to CSV files
seq_data.to_csv(f'./TRIBAL_seq_data_by_alleles.csv', index=False)
root_data.to_csv(f'./TRIBAL_root_data_by_alleles.csv', index=False)

# %%
# Check for invalid entries in the clonotype map
for k, v in clonotype_map.items():
    if not isinstance(v, str):
        print(k, v)

# Save the clonotype map as a JSON file
with open(f'clonotype_by_alleles.json', 'w') as f:
    json.dump(clonotype_map, f)

