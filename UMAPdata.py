import scanpy as sc
import pandas as pd
import json
import os

adata_path = "data/processed_data.h5ad"
if not os.path.exists(adata_path):
    raise FileNotFoundError(f"Cannot find {adata_path} in current directory.")

adata = sc.read(adata_path)

# Ensure PCA and neighbors
if 'X_pca' not in adata.obsm:
    sc.tl.pca(adata)
if 'neighbors' not in adata.uns:
    sc.pp.neighbors(adata)

# Compute UMAP if missing
if 'X_umap' not in adata.obsm:
    sc.tl.umap(adata)

# Build DataFrame for export
umap_df = pd.DataFrame(adata.obsm['X_umap'], columns=['x', 'y'])
umap_df['cell_id'] = adata.obs_names

# Recompute clustering
if 'leiden' not in adata.obs.columns or adata.obs['leiden'].isnull().all():
    print("⚠️ Leiden column missing or null, recomputing...")
    sc.tl.leiden(adata, resolution=0.6)

# Attach leiden cluster (as string)
umap_df['cluster'] = adata.obs['leiden'].astype(str).values

# Save updated AnnData file
adata.write(adata_path)
print("✅ Saved updated AnnData.")

# Export to JSON
output_path = "data/umap_clusters.json"
os.makedirs("data", exist_ok=True)
umap_df.to_json(output_path, orient="records", indent=2)
print(f"✅ Exported: {output_path}")
