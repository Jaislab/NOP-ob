"""Load deposited expression matrices into samples x genes tables keyed by manifest columns."""
from pathlib import Path
import gzip
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT/'data'/'raw'
PRIMARY = ['Pnoc', 'Oprl1']
PROCESSING = ['Pcsk1', 'Pcsk2', 'Cpe']
OPIOID = ['Penk', 'Pdyn', 'Pomc', 'Oprm1', 'Oprd1', 'Oprk1']
PANELS = {'primary': PRIMARY, 'secondary_processing': PROCESSING, 'exploratory_opioid': OPIOID}
ALL_TARGETS = PRIMARY + PROCESSING + OPIOID

def manifest(series, included_only=True):
    m = pd.read_csv(ROOT/'metadata'/'sample_manifest.csv')
    m = m[m.series == series]
    return m[m.included] if included_only else m

def ensembl_to_symbol():
    f = RAW/'GSE281091'/'GSE281091_Partek_Bulk_RNA-Seq_Exp_3A_Normalization_Normalized_counts_metabolic_states.txt.gz'
    a = pd.read_csv(f, sep='\t', usecols=['gene_id', 'gene_name'])
    b = pd.read_csv(str(f).replace('3A_Normalization_Normalized_counts_metabolic_states', '3B_Normalization_Normalized_counts_diet_types'), sep='\t', usecols=['gene_id', 'gene_name'])
    return pd.concat([a, b]).drop_duplicates('gene_id').set_index('gene_id')['gene_name']

def per_sample_files(series, xdir):
    m = manifest(series)
    cols = {}
    for _, r in m.iterrows():
        s = pd.read_csv(xdir/r.matrix_column, sep='\t', header=None, comment=None)
        if s.iloc[0, 0] == 'Gene':
            s = s.iloc[1:]
        s = s.set_index(0)[1].astype(float)
        cols[r.gsm] = s
    counts = pd.DataFrame(cols).T
    return counts, m.set_index('gsm')

def load(series):
    """Return (samples x genes matrix, sample metadata indexed by GSM, notes dict)."""
    notes = {}
    if series == 'GSE225224':
        counts, meta = per_sample_files(series, RAW/series/'x')
        notes['ids'] = 'gene symbols (NCBI mm10 annotation)'
    elif series == 'GSE198895':
        counts, meta = per_sample_files(series, RAW/series/'x')
        special = [c for c in counts.columns if str(c).startswith('__')]
        notes['htseq_special_rows_removed'] = {c: int(counts[c].sum()) for c in special}
        counts = counts.drop(columns=special)
        sym = ensembl_to_symbol()
        notes['ensembl_mapped_fraction'] = float(np.mean(counts.columns.isin(sym.index)))
        counts.columns = [sym.get(c, c) for c in counts.columns]
        notes['ids'] = 'Ensembl IDs mapped to symbols via the GSE281091 Ensembl annotation; unmapped IDs kept as Ensembl'
    elif series == 'GSE320302':
        d = pd.read_csv(RAW/series/'GSE320302_rsem_count.txt.gz', sep='\t', index_col=0)
        meta = manifest(series).set_index('gsm')
        counts = d[meta.matrix_column].T
        counts.index = meta.index
        frac = counts.values - np.round(counts.values)
        notes['nonint_fraction'] = float(np.mean(np.abs(frac) > 1e-9))
        notes['ids'] = 'gene symbols (UCSC mm10)'
    elif series == 'GSE281091':
        d = pd.read_csv(RAW/series/'GSE281091_Partek_Bulk_RNA-Seq_Exp_3B_Normalization_Normalized_counts_diet_types.txt.gz', sep='\t')
        meta = manifest(series).set_index('gsm')
        d = d.set_index('gene_name')
        counts = d[meta.matrix_column].T
        counts.index = meta.index
        notes['ids'] = 'Ensembl mm39 gene_name'
    else:
        raise ValueError(series)
    dup = counts.columns[counts.columns.duplicated()].unique()
    notes['duplicated_ids'] = list(map(str, dup))
    if len(dup):
        # Sum duplicated symbols only if none are targets; targets must be unique
        clash = [g for g in dup if g in ALL_TARGETS]
        if clash:
            raise RuntimeError(f'{series}: duplicated target identifiers {clash}')
        counts = counts.T.groupby(level=0).sum().T
    notes['targets_missing'] = [g for g in ALL_TARGETS if g not in counts.columns]
    return counts, meta, notes

if __name__ == '__main__':
    for s in ['GSE225224', 'GSE198895', 'GSE281091', 'GSE320302']:
        c, m, n = load(s)
        lib = c.sum(axis=1)
        print('====', s, c.shape, n)
        print(pd.DataFrame({'group': m.group, 'libsize_M': (lib/1e6).round(2), 'genes_detected': (c > 0).sum(axis=1)}).to_string())
