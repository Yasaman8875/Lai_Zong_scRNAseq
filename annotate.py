import anndata
import scanpy as sc
from scnym.api import scnym_api, atlas2target

adata = anndata.read_h5ad("results/py_objects/DMSO_seurat.h5ad")
adata.obs = adata.obs.astype('category')
sc.pp.normalize_total(adata, target_sum=1e6) # To get the CPM
sc.pp.log1p(adata) 
adata.raw = adata

joint_adata = atlas2target(
    adata=adata,
    species='human',
    key_added='annotations',
)
scnym_api(
    adata=joint_adata,
    task='train',
    groupby='annotations',
    out_path='./scnym_output',
    config='new_identity_discovery',
)

scnym_api(
    adata=adata,
    task='predict',
    key_added='scNym',
    config='new_identity_discovery',
    trained_model='./scnym_output'
)

sc.pl.umap(
    adata,
    color='scNym',
    show = False, 
    save = 'scNym.png'
)

sc.pl.umap(
  adata,
  color='scNym_confidence',
  show = False,
  save = 'scNym_confidence.png'
)

sc.pl.umap(
    adata, color = "integrated_snn_res.0.3", show = False, 
    size = 50, save = 'clusters.png'
)
