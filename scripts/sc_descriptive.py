"""Descriptive replicate-aware summaries for GSE233163 (no diet inference: diet is confounded with run).

Pool = ganglion site x sequencing run (12k/14k loadings and _1/_2 libraries within a run are
merged because GEO/paper describe ganglia pooled per condition; they are not independent animals).
Coverage rule (fixed before inspecting targets): report a pool x cell-type pseudobulk only if it has
>= 20 cells.
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from load_data import ROOT
P = ROOT/'data'/'processed'/'GSE233163'
OUT = ROOT/'results'/'single_cell'; OUT.mkdir(parents=True, exist_ok=True)
MIN_CELLS = 20
GENES = ['Oprl1', 'Pcsk1', 'Pcsk2', 'Cpe', 'Penk', 'Pdyn', 'Pomc', 'Oprm1']

cells = pd.read_csv(P/'target_genes_per_cell.csv.gz', index_col=0)
site = lambda lib: 'SCG' if lib.lower().startswith('scg') else 'stellate'
cells['site'] = cells.library.map(site)
cells['pool'] = cells.site + '_' + cells.run
cells['diet'] = cells.condition.map({'baseline': 'Control', 'obese': 'Hypercaloric'})
pb = pd.read_csv(P/'pseudobulk_library_celltype.csv.gz', index_col=0)
lib_ct = pb.index.to_series().str.split('|', expand=True)
pb['pool'] = lib_ct[0].map(site).values + '_' + lib_ct[0].str.rsplit('_', n=1).str[-1].values
pb['cell_type'] = lib_ct[1].values
genes_all = [c for c in pb.columns if c not in ('pool', 'cell_type')]
pool_pb = pb.groupby(['pool', 'cell_type'])[genes_all].sum()
lib_total = pool_pb.sum(axis=1)
rows = []
for (pool, ct), g in cells.groupby(['pool', 'cell_type']):
    r = dict(pool=pool, site=g.site.iat[0], run=g.run.iat[0], diet=g.diet.iat[0], cell_type=ct, n_cells=len(g),
             n_libraries=g.library.nunique(), libraries=';'.join(sorted(g.library.unique())), pseudobulk_umi=lib_total.loc[(pool, ct)],
             meets_coverage=len(g) >= MIN_CELLS)
    for gene in GENES:
        r[f'{gene}_frac_detected'] = (g[gene] > 0).mean()
        r[f'{gene}_pseudobulk_counts'] = pool_pb.loc[(pool, ct), gene]
        r[f'{gene}_cpm'] = pool_pb.loc[(pool, ct), gene] / lib_total.loc[(pool, ct)] * 1e6
    rows.append(r)
t = pd.DataFrame(rows)
t.to_csv(OUT/'GSE233163_pool_celltype_summary.csv', index=False)
comp = cells.groupby(['pool', 'diet', 'cell_type']).size().unstack(fill_value=0)
(comp.div(comp.sum(axis=1), axis=0)).to_csv(OUT/'GSE233163_celltype_proportions.csv')
print(t[t.cell_type == 'neuron'][['pool', 'diet', 'n_cells', 'n_libraries', 'Oprl1_frac_detected', 'Oprl1_pseudobulk_counts', 'Oprl1_cpm', 'Pcsk1_cpm', 'Cpe_cpm']].round(3).to_string())
print('Pnoc in deposited matrix:', 'Pnoc' in genes_all, '| n genes:', len(genes_all))
