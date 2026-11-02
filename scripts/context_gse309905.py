"""Descriptive Pnoc/Oprl1 detection in HFD-fed vagal snRNA-seq (GSE309905; no control diet).

Putative neuronal nuclei are defined by a marker rule fixed before inspecting targets:
>=1 UMI of Snap25 and >=1 UMI of Phox2b or Prdm12. No diet inference is possible.
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from scipy.io import mmread

sys.path.insert(0, str(Path(__file__).parent))
from load_data import ROOT
X = ROOT/'data'/'raw'/'GSE309905'/'x'
rows = []
for gsm, lab in [('GSM9284683', 'WT1'), ('GSM9284684', 'WT2'), ('GSM9284681', 'KO1'), ('GSM9284682', 'KO2')]:
    f = pd.read_csv(X/f'{gsm}_{lab}features.tsv.gz', sep='\t', header=None)
    m = mmread(X/f'{gsm}_{lab}matrix.mtx.gz').tocsr()  # genes x nuclei
    idx = {g: i for i, g in enumerate(f[1])}
    get = lambda g: np.asarray(m[idx[g]].todense()).ravel()
    umi = np.asarray(m.sum(axis=0)).ravel()
    neuron = (get('Snap25') > 0) & ((get('Phox2b') > 0) | (get('Prdm12') > 0))
    for grp, sel in [('putative neuronal nuclei', neuron), ('other nuclei', ~neuron)]:
        r = dict(gsm=gsm, library=lab, genotype='Pten WT' if lab.startswith('WT') else 'vagal Pten KO', diet='60% HFD (no control diet in series)',
                 nuclei_group=grp, n_nuclei=int(sel.sum()), median_umi=float(np.median(umi[sel])) if sel.any() else np.nan)
        for g in ['Pnoc', 'Oprl1', 'Pcsk1', 'Pcsk2', 'Cpe']:
            v = get(g)[sel]
            r[f'{g}_frac_detected'] = float((v > 0).mean()) if sel.any() else np.nan
            r[f'{g}_umi_per_10k'] = float(v.sum() / umi[sel].sum() * 1e4) if sel.any() else np.nan
        rows.append(r)
out = pd.DataFrame(rows)
(ROOT/'results'/'single_cell').mkdir(parents=True, exist_ok=True)
out.to_csv(ROOT/'results'/'single_cell'/'GSE309905_context_detection.csv', index=False)
print(out.round(4).to_string())
