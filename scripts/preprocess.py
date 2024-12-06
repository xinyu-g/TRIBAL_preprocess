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
# Filter out multiplets (cells with multiple barcodes)
barcode_rm_multiplet = barcodes[barcodes['isMultiplet'] == 0]

# Filter annotations for barcodes that are not multiplets
annotation_rm_multiplet = annotation[annotation['cell_barcode'].isin(barcode_rm_multiplet['Barcode'])]

# Filter for full-length annotations
annotation_full_length = annotation[annotation['full_length'] == 1]

# Filter barcodes and annotations for full-length and non-multiplets
barcode_rm_multiplet_full_length = barcode_rm_multiplet[barcode_rm_multiplet['Barcode'].isin(annotation_full_length['cell_barcode'].to_list())]
annotation_rm_multiplet_full_length = annotation_rm_multiplet[annotation_rm_multiplet['full_length'] == 1]

# %%
# Group cell barcodes by locus and count the number of loci per barcode
cell_barcode_by_locus_num = annotation_rm_multiplet_full_length.groupby('cell_barcode')['locus'].apply(list).reset_index()
cell_barcode_by_locus_num.loc[:, 'len'] = cell_barcode_by_locus_num['locus'].apply(lambda x: len(x))

# Select cells with exactly two loci
cell_barcode_by_locus_num_more = cell_barcode_by_locus_num[cell_barcode_by_locus_num['len'] == 2]

# Filter barcodes for cells with two loci
barcodes_full_length_more = barcode_rm_multiplet_full_length[barcode_rm_multiplet_full_length['Barcode'].isin(cell_barcode_by_locus_num_more['cell_barcode'].to_list())]

# %%
# Columns to select for further processing
select_columns = [
    'Barcode',
    'IGK_V', 'IGK_D', 'IGK_J', 'IGK_C', 'IGK_cdr3_aa',
    'IGL_V', 'IGL_D', 'IGL_J', 'IGL_C', 'IGL_cdr3_aa',
    'IGH_V', 'IGH_D', 'IGH_J', 'IGH_C', 'IGH_cdr3_aa'
]

# %%
# Select relevant columns and add clonotype information
barcodes_full_length_more = barcodes_full_length_more[select_columns]

# Concatenate CDR3 sequences to define clonotypes
barcodes_full_length_more.loc[:, 'clonotype'] = barcodes_full_length_more[['IGK_cdr3_aa', 'IGL_cdr3_aa', 'IGH_cdr3_aa']].apply(
    lambda row: '_'.join(row.dropna()), axis=1
)

# Create a mapping of clonotype to ID
clonotype['clonotype'] = clonotype[['IGK/L', 'IGH']].apply(
    lambda row: '_'.join(row.dropna()), axis=1
)
clonotype_map = dict(zip(clonotype['clonotype'], clonotype['clonotype_id']))

# Map clonotype IDs to the data
barcodes_full_length_more.loc[:, 'clonotype_id'] = barcodes_full_length_more['clonotype'].map(clonotype_map)

# %%
# Group barcodes by clonotype ID and aggregate as lists
barcodes_group = barcodes_full_length_more.groupby('clonotype_id').agg(
    {
        'Barcode': list
    }
)

# Filter for clonotypes with at least two associated cells
barcodes_group = barcodes_group[barcodes_group['Barcode'].str.len() >= 2]

# Add clonotype information to the grouped DataFrame
barcodes_group.loc[:, 'clonotype'] = barcodes_group.index

# %%
# Initialize lists to store output data
records = []
clono_records = []
repeat_locus = []

# Iterate through each row in the grouped barcodes
for i, t in barcodes_group.iterrows():
    consensus_H = []  # Consensus sequences for heavy chains
    consensus_L = []  # Consensus sequences for light chains
    clonotype = t['clonotype']

    for cell in t['Barcode']:
        try:
            # Filter annotation data for the current cell barcode
            temp = annotation_rm_multiplet_full_length[
                annotation_rm_multiplet_full_length['cell_barcode'] == cell
            ]

            if len(temp['locus'].unique()) < 2:
                print(cell)
                repeat_locus.append(cell)
                continue

            # Extract sequences and related information
            IGH = temp.loc[temp['locus'] == 'IGH', 'sequence'].values[0]
            IGL = temp.loc[temp['locus'].isin(['IGL', 'IGK']), 'sequence'].values[0]
            H_v = temp.loc[temp['locus'] == 'IGH', 'v_call'].values[0]
            L_v = temp.loc[temp['locus'].isin(['IGL', 'IGK']), 'v_call'].values[0]

            IGH_C = barcode_rm_multiplet.loc[
                barcode_rm_multiplet['Barcode'] == cell, 'IGH_C'
            ].values[0]

            H_germ = temp.loc[temp['locus'] == 'IGH', 'germline_alignment'].values[0]
            L_germ = temp.loc[temp['locus'].isin(['IGL', 'IGK']), 'germline_alignment'].values[0]

            # Add germline alignments to consensus lists
            consensus_H.append(H_germ)
            consensus_L.append(L_germ)

            # Add record for the cell
            records.append({
                'cellid': cell,
                'clonotype': clonotype,
                'heavy_chain_isotype': IGH_C,
                'heavy_chain_seq': IGH,
                'heavy_chain_v_allele': H_v,
                'light_chain_seq': IGL,
                'light_chain_v_allele': L_v
            })

        except IndexError as e:
            print(f"Error processing cell {cell}: {e}")
            continue  # Skip this cell if data is missing or invalid

    # Add consensus sequences for the current clonotype
    try:
        clono_records.append({
            'clonotype': clonotype,
            'heavy_chain_root': consensus(consensus_H),
            'light_chain_root': consensus(consensus_L)
        })
    except Exception as e:
        print(f"Error creating consensus for clonotype {clonotype}: {e}")
        continue

# %%
# Convert records to DataFrames and save as CSV files
seq_data = pd.DataFrame.from_records(records)
root_data = pd.DataFrame.from_records(clono_records)

seq_data.to_csv(f'./TRIBAL_seq_data_by_cdrh3.csv', index=False)
root_data.to_csv(f'./TRIBAL_root_data_by_cdrh3.csv', index=False)

# %%
# Save clonotype mapping as a JSON file
with open(f'clonotype_by_cdrh3.json', 'w') as f:
    json.dump(clonotype_map, f)

