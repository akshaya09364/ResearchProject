# scanpy_analysis_starter.py

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import os
import gzip

# Path to the directory containing the .csv.gz files
data_dir = r"C:/Users/Dell/Desktop/Steered Research Project/GSE67835_RAW"
files = [f for f in os.listdir(data_dir) if f.endswith(".csv.gz")]

all_cells = []

for f in files:
    file_path = os.path.join(data_dir, f)
    try:
        with gzip.open(file_path, 'rt') as fh:
            df = pd.read_csv(fh, index_col=0, header=None, sep='\t')

        if df.shape[0] == 0 or df.shape[1] != 1:
            print(f"Skipping empty or malformed file: {f}")
            continue

        df.columns = [f.split("_")[0]]
        all_cells.append(df)

    except Exception as e:
        print(f"Error reading {f}: {e}")
        continue

if not all_cells:
    raise RuntimeError("No valid gene count files were loaded. Please check your input directory.")

combined_df = pd.concat(all_cells, axis=1)
combined_df.to_csv("gene_counts_matrix.csv")

adata = sc.AnnData(combined_df.T)

# Fix BEFORE calling make_unique
adata.var.index.name = None
adata.var_names_make_unique()

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata = adata[:, adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.tl.leiden(adata)

sc.pl.umap(adata, color=['leiden'], save="_clusters.png")
adata.write("processed_data.h5ad")
