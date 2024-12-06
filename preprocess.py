import dandelion as ddl
import os
import yaml
from tqdm import tqdm
import pandas as pd
import numpy as np
import scanpy as sc
import warnings
import matplotlib
matplotlib.use('Agg')
warnings.simplefilter("ignore")
warnings.filterwarnings("ignore")
os.environ["PYTHONWARNINGS"] = "ignore:invalid escape sequence:SyntaxWarning"


def load_config(config_path='config.yaml'):
    """Load configuration file."""
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)
    return config


def consensus(sequences):
    """Generate a consensus sequence from a list of sequences."""
    if len(sequences) <= 1:
        return sequences[0]
    
    max_len = max(len(seq) for seq in sequences)
    
    cons_df = pd.DataFrame(0, index=range(max_len), columns=['A', 'T', 'C', 'G', 'N'])
    
    for seq in sequences:
        for i, nt in enumerate(seq):
            if nt in cons_df.columns:
                cons_df.loc[i, nt] += 1

    return ''.join(cons_df.idxmax(axis=1))


def main():
    # Load the configuration from the YAML file
    config = load_config()

    # Set species and data directory from YAML config
    species = config['species']
    data_dir = config['data_dir']

    # Set environment variables from config
    os.environ["IGDATA"] = config['igblastdb_dir']
    os.environ["GERMLINE"] = config['germline_dir']
    os.environ["BLASTDB"] = config['blastdb_dir']

    blast_db = os.environ["BLASTDB"]
    igblast_db = os.environ["IGDATA"]
    germline_db = os.environ["GERMLINE"]

    # Make BLAST DB
    ddl.utl.makeblastdb(f'{blast_db}/{config["blastdb_bcr_c"].format(species=species)}')

    # Get sample ID names to prefix Dandelion rows
    prefix = os.path.basename(os.path.normpath(data_dir))

    # Format and reannotate contigs                                   
    ddl.pp.format_fastas(f'{data_dir}', prefix=prefix)
    ddl.pp.reannotate_genes(f'{data_dir}', flavour='strict', org=species)

    # Reassign alleles to contigs
    ddl.pp.reassign_alleles(
        data_dir,
        combined_folder=f'{prefix}_comb',
        org=species,
        germline=f'{germline_db}/{config["germline_vdj"].format(species=species)}',
        v_germline=f'{germline_db}/{config["v_germline_file"].format(species=species)}'
    )

    # Reassign Isotypes
    ddl.pp.assign_isotypes(
        data_dir,
        org=species,
        correction_dict={},
        blastdb=f'{blast_db}/{config["blastdb_bcr_c"].format(species=species)}'
    )

    # Mutational Load Analysis
    for s in tqdm([data_dir], desc='Basic mutational load analysis'):
        file_path = f'{s}/{config["contig_file"]}'
        ddl.pp.quantify_mutations(file_path)

    # Load Gene Expression data
    adata = sc.read_10x_h5(f'{data_dir}/{config["filtered_feature_h5"]}')
    adata.obs_names = [f'{prefix}_{j}' for j in adata.obs_names]
    adata.obs_names = [j.split('-')[0] for j in adata.obs_names]
    adata.var_names_make_unique()

    # Remove doublets
    ddl.pp.recipe_scanpy_qc(adata, mito_cutoff=None)
    adata = adata[adata.obs['is_doublet'] == 'False'].copy()

    # Log-Normalize Gene Expression
    sc.pp.filter_genes(adata, min_cells=3)
    adata.raw = adata
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Load BCR data
    bcr = pd.read_csv(f'{data_dir}/{config["contig_file"]}', sep='\t')
    bcr.reset_index(inplace=True, drop=True)

    # Filter for only cells in both BCR and Gene Expression
    vdj, adata = ddl.pp.check_contigs(bcr, adata, library_type='ig')

    # Filter productive BCRs
    filt_vdj = vdj[vdj.metadata['chain_status'] == 'Single pair']
    filt_vdj = filt_vdj[filt_vdj.data['productive'] == 'T']

    # Identify clonotype clusters
    ddl.tl.find_clones(filt_vdj)

    # Get IGH with max expression
    if species == 'mouse':
        isotypes = config['mouse_isotypes']
    else:
        isotypes = config['human_isotypes']
    
    exp_igh = adata[filt_vdj.metadata.index].to_df()[isotypes].idxmax(axis=1)

    filt_vdj.data['cellId'] = [x.split('_contig')[0] for x in filt_vdj.data['sequence_id']]
    # Prepare TRIBE input table
    columns = config['output_columns']
    index = filt_vdj.metadata.index.str.replace('_', '-')
    out_df = pd.DataFrame(index=index, columns=columns)
    out_df.index.name = "cellid"

    for i in index:
        try:
            sub = filt_vdj.data.loc[filt_vdj.data['cellId'] == i]
            hc = sub.loc[sub['v_call'].str.contains('IGH')]
            lc = sub.loc[~sub['v_call'].str.contains('IGH')]
            out_df.loc[i, 'heavy_chain_seq'] = hc['sequence'].values[0]
            out_df.loc[i, 'light_chain_seq'] = lc['sequence'].values[0]
            out_df.loc[i, 'heavy_chain_isotype'] = hc['c_call_10x'].values[0]
            out_df.loc[i, 'light_chain_isotype'] = lc['c_call_10x'].values[0]
            out_df.loc[i, 'heavy_chain_v_allele'] = filt_vdj.metadata.loc[i]['v_call_genotyped_VDJ']
            out_df.loc[i, 'heavy_chain_d_allele'] = filt_vdj.metadata.loc[i]['d_call_VDJ']
            out_df.loc[i, 'heavy_chain_j_allele'] = filt_vdj.metadata.loc[i]['j_call_VDJ']
            out_df.loc[i, 'light_chain_v_allele'] = filt_vdj.metadata.loc[i]['v_call_genotyped_VJ']
            out_df.loc[i, 'light_chain_j_allele'] = filt_vdj.metadata.loc[i]['j_call_VJ']
            out_df.loc[i, 'clonotype'] = filt_vdj.metadata.loc[i]['clone_id']
            out_df.loc[i, 'heavy_chain_isotype_expression'] = exp_igh.loc[i]
        except Exception as e:
            print(f"Error processing cell {i}: {e}")
            break

    # Drop invalid entries and save
    out_df.replace('', np.nan, inplace=True)
    out_df.dropna(inplace=True)
    out_df.to_csv(f'{data_dir}/{config["tribal_input_table"]}')

    # Generate TRIBE clonotype root sequence table
    root_df = pd.DataFrame(index=out_df['clonotype'].unique(), columns=['light_chain_root', 'heavy_chain_root'])
    root_df.index.name = 'clonotype'

    for cid in root_df.index:
        cid_df = filt_vdj.data.loc[filt_vdj.data['clone_id'] == cid]
        heavy_seq = [seq.replace('.', '') for seq in cid_df.loc[cid_df['c_call_10x'].str.contains('IGH', na=False), 'germline_alignment']]
        light_seq = [seq.replace('.', '') for seq in cid_df.loc[~cid_df['c_call_10x'].str.contains('IGH', na=False), 'germline_alignment']]

        heavy_cons = consensus(heavy_seq) if heavy_seq else ''
        light_cons = consensus(light_seq) if light_seq else ''
        
        root_df.loc[cid] = [light_cons, heavy_cons]

    # Save root sequence table
    root_df.to_csv(f'{data_dir}/{config["tribal_root_sequences"]}')


if __name__ == "__main__":
    main()