import anndata
import pandas as pd
import numpy as np
import argparse 
import plotnine as p9
import matplotlib.pyplot as plt
import scanpy as sc
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputs', type=str, nargs='+', help='Input files (10x folders)')
    parser.add_argument('--h5ad', type=str, help='Input file (.h5ad file)')
    parser.add_argument('--experiment_ids', type=str, nargs='+', help='experiment ids for multiple inputs')
    parser.add_argument('--barcodes', type=str, required=True, help='Barcodes file')
    parser.add_argument('--gene', type=str, required=True, help='Barcodes file')
    parser.add_argument('--output', type=str, required=True, help='Output file')
    return parser.parse_args()


def plot_gene(data, barcodes, output, gene):
    """
    On X axis, experiment_id
    On Y axis, gene expression from cp10k layer
    """

    # Get adata.obs as a dataframe
    obs = data.obs.copy()
    obs[f'{gene} (cp10k)'] = data.layers['cp10k'].todense()[:, data.var.gene_symbols == gene]
    
    # Filter data to only those found in barcodes
    print('Filtering data to only those found in barcodes')
    print(f"Before: {obs.shape}")
    obs = obs[obs.index.isin(barcodes)]
    print(f"After: {obs.shape}")

    # Plot
    p = p9.ggplot(obs, p9.aes(x='experiment_id', y=f'{gene} (cp10k)')) + \
        p9.geom_jitter() + \
        p9.theme_bw() + \
        p9.theme(axis_text_x=p9.element_text(angle=90, hjust=1))
 
    p.save(output, width=10, height=10, dpi=300)

    p = p + p9.coord_cartesian(ylim=(0, 100))

    p.save(output.replace('.png', '_zoom.png'), width=10, height=10, dpi=300)


def main():
    # Parse arguments
    args = parse_args()

    # Load data
    data_files = []
    if (args.inputs is None and args.h5ad is None) or (args.inputs is not None and args.h5ad is not None):
        raise ValueError("Must have either 10x folder inputs or .h5ads, not none or both.")

    if args.inputs is not None:
        for f in args.inputs:
            print("Getting file", f)
            data_files.append(sc.read_10x_mtx(Path(f), var_names='gene_ids'))

        data = anndata.AnnData.concatenate(*data_files, batch_key='experiment_id', batch_categories=args.experiment_ids)
    else:

        data = anndata.read_h5ad(Path(args.h5ad))


    # Create the cp10k layer
    print('Creating the cp10k layer')
    data.layers['cp10k'] = data.X.copy()
    sc.pp.normalize_total(data, target_sum=10000, layer='cp10k')

    # Load barcodes
    barcodes = pd.read_csv(args.barcodes, header=None, names=['barcode'])
    barcodes = barcodes['barcode'].tolist()

    # Plot gene
    plot_gene(data, barcodes, args.output, args.gene)


if __name__ == "__main__":
    main()
