"""Stream the GSE233163 dense cell x gene CSV into library x cell-type pseudobulk sums.

Keeps per-cell UMI values for target genes for detection summaries. Uses deposited
cell annotations and QC as provided; no cells are selected using target-gene values.
"""
from pathlib import Path
import gzip, sys
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from load_data import ALL_TARGETS, ROOT

RAW = ROOT/'data'/'raw'/'GSE233163'
OUT = ROOT/'data'/'processed'/'GSE233163'
OUT.mkdir(parents=True, exist_ok=True)

meta = pd.read_csv(RAW/'GSE233163_cells_metadata.csv.gz', index_col=0)
meta['library'] = meta.ori_sample + '_' + meta.run
with gzip.open(RAW/'GSE233163_cells_counts.csv.gz', 'rt') as f:
    header = [h.strip('"') for h in f.readline().rstrip('\n').split(',')[1:]]
    assert len(header) == len(meta) and set(header) == set(meta.index), 'cell barcodes do not match metadata'
    m = meta.loc[header]
    keys = (m.library + '|' + m.cell_type).values
    groups, gidx = np.unique(keys, return_inverse=True)
    sums, genes, percell, nonint = [], [], {}, 0
    for i, line in enumerate(f):
        name, _, rest = line.partition(',')
        name = name.strip('"')
        v = np.array(rest.split(','), dtype=float)
        if i < 2000: nonint += int(np.any(np.abs(v - np.round(v)) > 1e-9))
        sums.append(np.bincount(gidx, weights=v, minlength=len(groups)))
        genes.append(name)
        if name in ALL_TARGETS:
            percell[name] = v
        if i % 2000 == 0: print(i, name, flush=True)
pb = pd.DataFrame(np.array(sums).T, index=groups, columns=genes)
pb.index.name = 'library|cell_type'
pb.to_parquet(OUT/'pseudobulk_library_celltype.parquet') if False else pb.to_csv(OUT/'pseudobulk_library_celltype.csv.gz')
cells = m[['library', 'tissue', 'run', 'condition', 'cell_type', 'nCount_RNA', 'nFeature_RNA', 'percent.mt']].copy()
for g, v in percell.items():
    cells[g] = v
cells.to_csv(OUT/'target_genes_per_cell.csv.gz')
print('genes', len(genes), 'groups', len(groups), 'non-integer rows among first 2000:', nonint)
