"""Expression of HFD/fasting counter-regulated candidates in the mouse geniculate ganglion.

Datasets (no diet manipulation; expression only):
- GSE223724: bulk geniculate ganglia, wild-type adults (n = 4; Egr4-/- excluded), raw counts + gene length.
- GSE102443: 96 single geniculate neurons (Fluidigm C1, FPKM).
- GSE135801: Phox2b-expressing geniculate neurons (454 cells, read counts); cells with < 100,000 reads
  excluded (rule fixed before inspecting target genes).
References: nodose (GSE281091, chow fed, n = 6) and DRG (GSE225224, standard diet, n = 11) converted to
TPM-like units with the same length annotation.
Expression calls (bulk, mean TPM): not detected < 0.5; low 0.5-5; moderate 5-50; high >= 50.
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from load_data import ROOT, RAW, load

OUT = ROOT/'results'/'geniculate'; OUT.mkdir(parents=True, exist_ok=True)
TARGETS = ['Acacb', 'Insr', 'Cnr1', 'Pomc', 'Cartpt']
CONTEXT = ['Pnoc', 'Oprl1']
MARKERS = ['Phox2b', 'Prdm12', 'P2rx3', 'Snap25']
GENES = TARGETS + CONTEXT + MARKERS
MIN_READS = 100_000

def tpm(counts, length):
    rpk = counts.div(length / 1e3, axis=1)
    return rpk.div(rpk.sum(axis=1), axis=0) * 1e6

def call(x):
    return 'not detected' if x < 0.5 else 'low' if x < 5 else 'moderate' if x < 50 else 'high'

def bulk():
    x = pd.read_excel(RAW/'GSE223724'/'GSE223724_ILMN_1150_Pierchala_mRNAseq8_Sept2021_processed_data.xlsx', header=[0, 1])
    sym = x[('Gene info', 'GeneSymbol')].astype(str); L = x[('Gene info', 'Length')].astype(float)
    wt = [c for c in x['Raw counts'].columns if '_WT_' in c]
    cnt = x['Raw counts'][wt].astype(float); cnt.index = sym
    L.index = sym
    keep = ~cnt.index.duplicated()
    cnt, L = cnt[keep], L[keep]
    gg = tpm(cnt.T, L)
    length = L
    # pituitary carry-over check (geniculate dissection is at the skull base): Gh/Prl/Cga per sample
    pit = gg[['Gh', 'Prl', 'Cga', 'Pomc']].copy()
    pit['pituitary_contaminated'] = (pit.Gh > 10) | (pit.Prl > 50)
    pit.to_csv(OUT/'geniculate_pituitary_check.csv')
    # nodose: deposited DESeq2-normalized values for fed chow; DRG: raw counts, standard diet, both ages
    nod_raw = pd.read_csv(RAW/'GSE281091'/'GSE281091_Partek_Bulk_RNA-Seq_Exp_3B_Normalization_Normalized_counts_diet_types.txt.gz', sep='\t')
    m = pd.read_csv(ROOT/'metadata'/'sample_manifest.csv')
    fed = m[(m.series == 'GSE281091') & (m.diet_detail == 'Non-fasted_regular chow')].matrix_column.tolist()
    nod = nod_raw.drop_duplicates('gene_name').set_index('gene_name')[fed].T
    drg_c, drg_m, _ = load('GSE225224')
    drg = drg_c.loc[drg_m.diet == 'Control']
    common_n = nod.columns.intersection(length.index); common_d = drg.columns.intersection(length.index)
    ref = {'Geniculate (GSE223724, WT)': gg, 'Nodose (GSE281091, chow fed)': tpm(nod[common_n], length[common_n]),
           'DRG (GSE225224, standard diet)': tpm(drg[common_d], length[common_d])}
    rows, per = [], []
    for name, t in ref.items():
        mean = t.mean()
        expressed = mean[mean > 0]
        for g in GENES:
            v = mean.get(g, np.nan)
            pct = float((expressed < v).mean() * 100) if pd.notna(v) else np.nan
            if name.startswith('Geniculate') and g == 'Pomc':
                clean = pit.index[~pit.pituitary_contaminated]
                rows.append(dict(dataset=name + ' [no pituitary signal]', gene=g, mean_tpm=float(t.loc[clean, g].mean()),
                                 min_tpm=float(t.loc[clean, g].min()), max_tpm=float(t.loc[clean, g].max()), n=len(clean),
                                 percentile_among_detected_genes=float((t.loc[clean].mean()[t.loc[clean].mean() > 0] < t.loc[clean, g].mean()).mean() * 100),
                                 call=call(float(t.loc[clean, g].mean()))))
            rows.append(dict(dataset=name, gene=g, mean_tpm=v, min_tpm=t[g].min() if g in t else np.nan,
                             max_tpm=t[g].max() if g in t else np.nan, n=len(t), percentile_among_detected_genes=pct,
                             call=call(v) if pd.notna(v) else 'absent from matrix'))
            if g in t:
                per += [dict(dataset=name, gene=g, sample=s, tpm=val) for s, val in t[g].items()]
    return pd.DataFrame(rows), pd.DataFrame(per)

def single_cell():
    rows = []
    f = pd.read_csv(RAW/'GSE102443'/'GSE102443_GEO-ID_Dvoryanchikov_2017_Datatable_FPKM.txt.gz', sep='\t')
    cells = [c for c in f.columns if '[[' in c]
    F = f.groupby('Gene')[cells].max().T
    cls = np.where(F['Phox2b'] > 10, 'Phox2b+ (oral/taste)', 'Phox2b- (somatosensory)')
    for grp in ['all neurons', 'Phox2b+ (oral/taste)', 'Phox2b- (somatosensory)']:
        sel = np.ones(len(F), bool) if grp == 'all neurons' else cls == grp
        for g in GENES:
            v = F.loc[sel, g] if g in F else None
            rows.append(dict(dataset='GSE102443 (C1, FPKM)', group=grp, n_cells=int(sel.sum()), gene=g,
                             frac_detected=float((v > 1).mean()) if v is not None else np.nan,
                             median_in_detected=float(v[v > 1].median()) if v is not None and (v > 1).any() else np.nan,
                             unit='FPKM; detected = FPKM > 1'))
    y = pd.read_excel(RAW/'GSE135801'/'x'/'GSM4037432_GG_scRNAseq_Phox2b_expressed454.xlsx').set_index('gene_id')
    y = y.groupby(level=0).sum()
    depth = y.sum(); y = y.loc[:, depth >= MIN_READS]
    cpm = y / y.sum() * 1e6
    for g in GENES:
        v = cpm.loc[g] if g in cpm.index else None
        rows.append(dict(dataset='GSE135801 (Phox2b+ neurons, reads)', group='all neurons (>= 100k reads)', n_cells=int(y.shape[1]), gene=g,
                         frac_detected=float((y.loc[g] > 0).mean()) if v is not None else np.nan,
                         median_in_detected=float(v[v > 0].median()) if v is not None and (v > 0).any() else np.nan,
                         unit='CPM; detected = >= 1 read'))
    return pd.DataFrame(rows)

def main():
    b, per = bulk(); s = single_cell()
    b.to_csv(OUT/'geniculate_bulk_expression.csv', index=False); per.to_csv(OUT/'bulk_per_sample_tpm.csv', index=False)
    s.to_csv(OUT/'geniculate_single_neuron_detection.csv', index=False)
    with pd.option_context('display.width', 250):
        print(b.pivot(index='gene', columns='dataset', values='mean_tpm').reindex(GENES).round(1))
        print(b.pivot(index='gene', columns='dataset', values='percentile_among_detected_genes').reindex(GENES).round(0))
        print(s.pivot_table(index='gene', columns=['dataset', 'group'], values='frac_detected').reindex(GENES).round(2).to_string())
        print(s.groupby(['dataset', 'group']).n_cells.first())

if __name__ == '__main__':
    main()
