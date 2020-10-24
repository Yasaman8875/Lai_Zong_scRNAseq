import scanpy as sc
import numpy as np
import pandas as pd
import scvelo as scv
from pathlib import Path
import pickle
import os
import cellrank as cr

clusters = 'integrated_snn_res.0.3'

##################
## RNA Velocity ##
##################

## Preparing sample data.

samples = {
  'DMSO' : 'results/py_objects/DMSO_seurat.h5ad',
  'EZH2i' : 'results/py_objects/EZH2i_seurat.h5ad',
  'RACi' : 'results/py_objects/RACi_seurat.h5ad',
  'Combo' : 'results/py_objects/Combo_seurat.h5ad'
}

samples = {x:scv.read(y) for x,y in samples.items()}

## Change the metadata to categorical.

for key in samples.keys():
    samples[key].obs = samples[key].obs.astype('category')

## Preprocess the data.

for key in samples.keys():
    scv.pp.filter_and_normalize(samples[key], min_shared_counts=20, n_top_genes=3000)
    scv.pp.moments(samples[key], n_pcs=30, n_neighbors=30)

## Calculate RNA velocities.

for key in samples.keys():
    scv.tl.velocity(samples[key])
    scv.tl.velocity_graph(samples[key])

## Plot RNA velocity streams.

outdir = 'results/trajectory/velocity/velocity_plots'
if not os.path.exists(outdir):
    os.makedirs(outdir)

scv.settings.figdir = outdir

for key,value in samples.items():
    scv.pl.velocity_embedding_stream(
      value, basis='umap', color=clusters,
      save = '%s.png' % key, title = key, show = False,
      figsize = (10, 10), size = 50, dpi = 300, legend_fontsize = 0
    )

## Plot RNA velocity arrows.

outdir = 'results/trajectory/velocity/velocity_arrows'
if not os.path.exists(outdir):
    os.makedirs(outdir)

scv.settings.figdir = outdir

for key,value in samples.items():
    scv.pl.velocity_embedding(
      value, arrow_length=3, arrow_size=2, dpi=300,
      basis ='umap', color=clusters,
      figsize = (10, 10), size = 50, show = False,
      save = '%s.png' % key, title = key
    )

## Plot velocity speed and coherence.

outdir = 'results/trajectory/velocity/velocity_speed_coherence'
if not os.path.exists(outdir):
    os.makedirs(outdir)

scv.settings.figdir = outdir

metrics = ['velocity_length', 'velocity_confidence']
for key,value in samples.items():
    scv.tl.velocity_confidence(value)
    scv.pl.scatter(
      value, c = metrics, cmap = 'gnuplot', perc=[5, 95],
      size = 50, show = False, dpi = 300, figsize = (10, 10),
      save = '%s.png' % key
    )

## Plot cell connections.

outdir = 'results/trajectory/velocity/velocity_connections'
if not os.path.exists(outdir):
    os.makedirs(outdir)

scv.settings.figdir = outdir

for key,value in samples.items():
    scv.pl.velocity_graph(
      value, threshold = .2, size = 50, show = False, dpi = 300,
      figsize = (10, 10), color = clusters,
      save = '%s.png' % key, title = key
    )

## Plot velocity pseudotime.

outdir = 'results/trajectory/velocity/velocity_pseudotime'
if not os.path.exists(outdir):
    os.makedirs(outdir)

scv.settings.figdir = outdir

for key,value in samples.items():
    scv.tl.velocity_pseudotime(value)
    scv.pl.scatter(
      value, color='velocity_pseudotime', cmap='gnuplot', dpi = 300,
      show = False, figsize = (10, 10), title = key, size = 50,
      save = '%s.png' % key
    )

## PAGA.

outdir = 'results/trajectory/velocity/velocity_paga'
if not os.path.exists(outdir):
    os.makedirs(outdir)

scv.settings.figdir = outdir

for key in samples.keys():
    samples[key].uns['neighbors']['distances'] = samples[key].obsp['distances']
    samples[key].uns['neighbors']['connectivities'] = samples[key].obsp['connectivities']

for key,value in samples.items():
    scv.tl.paga(value, groups = clusters)
    scv.pl.paga(
      value, basis = 'umap', color = clusters,
      dpi = 300, show = False, figsize = (10, 10), title = key, size = 50,
      save = '%s.png' % key
    )

## Export UMAP plots.

outdir = 'results/trajectory/velocity/velocity_umap'
if not os.path.exists(outdir):
    os.makedirs(outdir)

scv.settings.figdir = outdir

for key,value in samples.items():
    scv.pl.umap(
      value, color = clusters, show = False, figsize = (10, 10), title = key,
      size = 50, save = '{}.png'.format(key)
    )
    scv.pl.umap(
      value, show = False, figsize = (10, 10), title = key,
      size = 50, save = '{}_nocolor.png'.format(key)
    )

## Get important genes.

outdir = 'results/trajectory/velocity/velocity_genes'
if not os.path.exists(outdir):
    os.makedirs(outdir)

for value in samples.values():
    scv.tl.rank_velocity_genes(value, groupby = clusters, min_corr=.3)

for key,value in samples.items():
    df = scv.DataFrame(value.uns['rank_velocity_genes']['names'])
    df.to_csv("{}/{}.tsv".format(outdir, key), sep = '\t', header = True, index = False)

## Save the velocities.

with open('results/py_objects/velocities.pickle', 'wb') as handle:
    pickle.dump(samples, handle)

## Genes different between treatments and DMSO.

outdir = 'results/trajectory/velocity/velocity_diff_genes'
if not os.path.exists(outdir):
    os.makedirs(outdir)

comparisons = [
  ('EZH2i', 'DMSO'),
  ('RACi', 'DMSO'),
  ('Combo', 'DMSO')
]
diff_samples = {}
for x,y in  comparisons:
    comp = samples[x].concatenate(samples[y])
    comp.obs['groups'] = comp.obs[['orig.ident', clusters]].astype(str).agg('_'.join, axis = 1)
    diff_samples['{}_vs_{}'.format(x, y)] = comp

for value in diff_samples.values():
    scv.tl.rank_velocity_genes(value, groupby = 'groups', min_corr=.3)

for key,value in diff_samples.items():
    df = scv.DataFrame(value.uns['rank_velocity_genes']['names'])
    df.to_csv("{}/{}.tsv".format(outdir, key), sep = '\t', header = True, index = False)

