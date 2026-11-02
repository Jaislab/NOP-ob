"""Transcriptome-wide diet models for bulk/TRAP series; extract prespecified target tests.

Count series: PyDESeq2 0.5.2 (median-of-ratios size factors, full-transcriptome dispersions,
Wald tests, unshrunk MLE log2 fold changes). GSE281091 (only DESeq2-normalized values deposited):
limma-trend style moderated linear model on log2(normalized + 1), ~ sex + diet.
"""
from pathlib import Path
import json, warnings, sys
import numpy as np
import pandas as pd
from scipy import stats, special
from statsmodels.stats.multitest import multipletests
from statsmodels.nonparametric.smoothers_lowess import lowess
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference

sys.path.insert(0, str(Path(__file__).parent))
from load_data import load, PANELS, ALL_TARGETS, ROOT

OUT = ROOT/'results'
(OUT/'de').mkdir(parents=True, exist_ok=True)
(OUT/'qc').mkdir(parents=True, exist_ok=True)
warnings.filterwarnings('ignore', category=UserWarning)
INF = DefaultInference(n_cpus=4)
CONTRASTS = pd.read_csv(ROOT/'metadata'/'contrast_inventory.csv').set_index('contrast_id')
SEED = 0

def bh(p):
    p = np.asarray(p, float); out = np.full_like(p, np.nan); ok = ~np.isnan(p)
    if ok.any(): out[ok] = multipletests(p[ok], method='fdr_bh')[1]
    return out

# ---------------------------------------------------------------- PyDESeq2
def run_deseq(series, contrast_ids, design='~group', extra=None):
    counts, meta, notes = load(series)
    counts = counts.round().astype(int)          # RSEM expected counts -> integers; others already integer
    counts = counts.loc[:, counts.sum(axis=0) > 0]  # remove wholly unexpressed features only
    meta = meta.copy()
    meta['age'] = meta.group.str.split('_').str[0] if meta.group.str.contains('_').all() else 'na'
    dds = DeseqDataSet(counts=counts, metadata=meta[['group', 'diet', 'age']], design=design,
                       refit_cooks=True, inference=INF, quiet=True)
    dds.deseq2()
    sf = dds.obs['size_factors'].astype(float)
    tg = [g for g in ALL_TARGETS if g in dds.var_names]
    dds.var.loc[tg, ['genewise_dispersions', 'fitted_dispersions', 'MAP_dispersions', 'dispersions']].to_csv(OUT/'qc'/f'{series}_target_dispersions.csv')
    norm = counts.div(sf, axis=0)
    res = {}
    for cid in contrast_ids:
        c = CONTRASTS.loc[cid]
        st = DeseqStats(dds, contrast=['group', c.numerator, c.denominator], inference=INF, quiet=True,
                        cooks_filter=True, independent_filter=True)
        st.summary()
        r = st.results_df.copy()
        r['padj_transcriptome'] = r.pop('padj')
        num = meta.index[meta.group == c.numerator]; den = meta.index[meta.group == c.denominator]
        r['n_control'] = len(den); r['n_hypercaloric'] = len(num)
        r['mean_norm_control'] = norm.loc[den].mean(); r['mean_norm_hypercaloric'] = norm.loc[num].mean()
        r['frac_nonzero_control'] = (counts.loc[den] > 0).mean(); r['frac_nonzero_hypercaloric'] = (counts.loc[num] > 0).mean()
        r['ci_low'] = r.log2FoldChange - 1.959964 * r.lfcSE; r['ci_high'] = r.log2FoldChange + 1.959964 * r.lfcSE
        r['df'] = np.nan
        r = r.rename(columns={'log2FoldChange': 'log2fc', 'lfcSE': 'se', 'pvalue': 'p'})
        r.index.name = 'gene'
        r.to_csv(OUT/'de'/f'{cid}_transcriptome.csv.gz')
        res[cid] = r
    # QC: VST PCA, size factors, Cook's distance for targets
    dds.vst(use_design=False)
    vst = pd.DataFrame(dds.layers['vst_counts'], index=dds.obs_names, columns=dds.var_names)
    qc = qc_tables(series, counts, vst, meta, sf)
    cooks = pd.DataFrame(dds.layers['cooks'], index=dds.obs_names, columns=dds.var_names)
    notes['cooks_max_targets'] = {g: float(cooks[g].max()) for g in ALL_TARGETS if g in cooks}
    notes['n_genes_tested'] = int(counts.shape[1])
    long = norm[[g for g in ALL_TARGETS if g in norm]].assign(gsm=norm.index).melt('gsm', var_name='gene', value_name='normalized')
    long['raw_count'] = [counts.loc[s, g] for s, g in zip(long.gsm, long.gene)]
    long['size_factor'] = long.gsm.map(sf)
    long = long.merge(meta[['group', 'diet', 'title']], left_on='gsm', right_index=True)
    long.insert(0, 'series', series)
    # exploratory interaction (diet x age) where two ages exist
    inter = None
    if meta.age.nunique() == 2 and 'na' not in meta.age.values:
        d2 = DeseqDataSet(counts=counts, metadata=meta[['diet', 'age']], design='~age + diet + age:diet',
                          refit_cooks=True, inference=INF, quiet=True)
        d2.deseq2()
        cols = list(d2.obsm['design_matrix'].columns)
        vec = np.array([1.0 if (':' in col) else 0.0 for col in cols])
        s2 = DeseqStats(d2, contrast=vec, inference=INF, quiet=True); s2.summary()
        inter = s2.results_df.loc[[g for g in PANELS['primary']]].assign(series=series, term=[c for c in cols if ':' in c][0])
    return res, long, qc, notes, inter

def qc_tables(series, counts, logexpr, meta, sf):
    lib = counts.sum(axis=1)
    x = logexpr.loc[:, logexpr.var().sort_values(ascending=False).index[:500]]
    x = x - x.mean()
    u, s, vt = np.linalg.svd(x.values, full_matrices=False)
    var = s**2 / np.sum(s**2)
    pcs = pd.DataFrame(u[:, :3] * s[:3], index=x.index, columns=['PC1', 'PC2', 'PC3'])
    corr = np.corrcoef(logexpr.values)
    np.fill_diagonal(corr, np.nan)
    qc = pd.DataFrame({'series': series, 'group': meta.group, 'diet': meta.diet, 'title': meta.title,
                       'library_size': lib, 'size_factor': sf, 'genes_detected': (counts > 0).sum(axis=1),
                       'median_corr_to_others': np.nanmedian(corr, axis=1)}).join(pcs)
    qc['pc1_var'] = var[0]; qc['pc2_var'] = var[1]
    z = (qc.median_corr_to_others - qc.median_corr_to_others.median()) / (1.4826 * (qc.median_corr_to_others - qc.median_corr_to_others.median()).abs().median())
    qc['corr_robust_z'] = z; qc['flag_low_similarity'] = z < -3.5
    qc.index.name = 'gsm'
    qc.to_csv(OUT/'qc'/f'{series}_sample_qc.csv')
    return qc

# ---------------------------------------------------------------- moderated linear model
def trigamma_inverse(x):
    x = np.asarray(x, float)
    y = 0.5 + 1 / x
    for _ in range(50):
        tri = special.polygamma(1, y)
        dif = tri * (1 - tri / x) / special.polygamma(2, y)
        y = y + dif
        if np.all(-dif / y < 1e-8): break
    return y

def fit_fdist_trend(s2, df, amean):
    z = np.log(s2); e = z - special.digamma(df / 2) + np.log(df / 2)
    trend = lowess(e, amean, frac=0.3, return_sorted=False)
    resid = e - trend
    n = len(e)
    evar = np.mean(resid**2) * n / (n - 1) - np.mean(special.polygamma(1, df / 2))
    if evar > 0:
        d0 = 2 * trigamma_inverse(evar)
        s0 = np.exp(trend + special.digamma(d0 / 2) - np.log(d0 / 2))
    else:
        d0 = np.inf; s0 = np.exp(trend)
    return d0, s0

def run_moderated(series, cid, pseudocount=1.0, save=True):
    vals, meta, notes = load(series)
    vals = vals.loc[:, vals.sum(axis=0) > 0]
    y = np.log2(vals + pseudocount)
    X = pd.DataFrame({'intercept': 1.0, 'female': (meta.sex == 'female').astype(float),
                      'hfd': (meta.diet == 'Hypercaloric').astype(float)}, index=meta.index)
    Xv = X.values; XtXi = np.linalg.inv(Xv.T @ Xv)
    B = XtXi @ Xv.T @ y.values
    resid = y.values - Xv @ B
    dfr = Xv.shape[0] - Xv.shape[1]
    s2 = (resid**2).sum(axis=0) / dfr
    amean = y.values.mean(axis=0)
    ok = s2 > 1e-12
    d0, s0 = fit_fdist_trend(s2[ok], np.full(ok.sum(), dfr, float), amean[ok])
    s0full = np.full(len(s2), np.nan); s0full[ok] = s0
    post = (d0 * s0full + dfr * s2) / (d0 + dfr) if np.isfinite(d0) else s0full
    dft = dfr + d0 if np.isfinite(d0) else np.inf
    j = 2
    se = np.sqrt(post * XtXi[j, j]); beta = B[j]
    t = beta / se
    p = 2 * stats.t.sf(np.abs(t), dft)
    q = stats.t.ppf(0.975, dft)
    num = meta.index[meta.diet == 'Hypercaloric']; den = meta.index[meta.diet == 'Control']
    r = pd.DataFrame({'baseMean': vals.mean(), 'log2fc': beta, 'se': se, 'stat': t, 'p': p,
                      'ci_low': beta - q * se, 'ci_high': beta + q * se}, index=vals.columns)
    r['padj_transcriptome'] = bh(r.p)
    r['n_control'] = len(den); r['n_hypercaloric'] = len(num)
    r['mean_norm_control'] = vals.loc[den].mean(); r['mean_norm_hypercaloric'] = vals.loc[num].mean()
    r['frac_nonzero_control'] = (vals.loc[den] > 0).mean(); r['frac_nonzero_hypercaloric'] = (vals.loc[num] > 0).mean()
    r['df'] = dft
    r.index.name = 'gene'
    notes.update(prior_df=float(d0), residual_df=dfr, pseudocount=pseudocount, n_genes_tested=int(vals.shape[1]))
    if not save:
        return r, notes
    r.to_csv(OUT/'de'/f'{cid}_transcriptome.csv.gz')
    qc = qc_tables(series, vals, y, meta, pd.Series(np.nan, index=meta.index))
    long = vals[[g for g in ALL_TARGETS if g in vals]].assign(gsm=vals.index).melt('gsm', var_name='gene', value_name='normalized')
    long['raw_count'] = np.nan; long['size_factor'] = np.nan
    long = long.merge(meta[['group', 'diet', 'title', 'sex']], left_on='gsm', right_index=True)
    long.insert(0, 'series', series)
    # exploratory sex x diet interaction for primary targets
    Xi = X.assign(inter=X.female * X.hfd).values; XtXii = np.linalg.inv(Xi.T @ Xi)
    inter = []
    for g in PANELS['primary']:
        yy = y[g].values; b = XtXii @ Xi.T @ yy; rr = yy - Xi @ b; d = len(yy) - Xi.shape[1]
        s2g = rr @ rr / d; seg = np.sqrt(s2g * XtXii[3, 3]); tg = b[3] / seg
        inter.append(dict(gene=g, series=series, term='female:HFD (ordinary least squares, unmoderated)', log2FoldChange=b[3], lfcSE=seg, stat=tg, pvalue=2 * stats.t.sf(abs(tg), d)))
    return {cid: r}, long, qc, notes, pd.DataFrame(inter).set_index('gene')

# ---------------------------------------------------------------- main
def main():
    np.random.seed(SEED)
    allres, longs, qcs, notes, inters = {}, [], [], {}, []
    jobs = [('GSE225224', ['DRG_HFD_12wk', 'DRG_HFD_36wk']), ('GSE198895', ['DRG_Nav18_TRAP_WD']),
            ('GSE320302', ['Nodose_earlyHFD_P21', 'Nodose_earlyHFD_P35'])]
    for series, cids in jobs:
        print('fitting', series, flush=True)
        r, l, q, n, i = run_deseq(series, cids)
        allres.update(r); longs.append(l); qcs.append(q); notes[series] = n
        if i is not None: inters.append(i)
    print('fitting GSE281091', flush=True)
    r, l, q, n, i = run_moderated('GSE281091', 'Nodose_HFD_4wk')
    allres.update(r); longs.append(l); qcs.append(q); notes['GSE281091'] = n; inters.append(i)
    # pseudocount sensitivity for GSE281091 targets
    sens, _ = run_moderated('GSE281091', 'Nodose_HFD_4wk', pseudocount=0.5, save=False)
    sens = sens.loc[[g for g in ALL_TARGETS if g in sens.index], ['log2fc', 'se', 'p']].add_suffix('_pseudocount0.5')

    rows = []
    for cid, r in allres.items():
        c = CONTRASTS.loc[cid]
        for panel, genes in PANELS.items():
            for g in genes:
                d = dict(contrast_id=cid, series=c.series, ganglion=c.ganglion, assay=c.assay, exposure=c.exposure,
                         method=c.method, contrast_status=c.status, panel=panel, gene=g)
                if g not in r.index:
                    d.update(result_class='not detected (no counts in any sample)'); rows.append(d); continue
                x = r.loc[g]
                d.update({k: x[k] for k in ['baseMean', 'n_control', 'n_hypercaloric', 'mean_norm_control', 'mean_norm_hypercaloric',
                                             'frac_nonzero_control', 'frac_nonzero_hypercaloric', 'log2fc', 'se', 'ci_low', 'ci_high',
                                             'df', 'stat', 'p', 'padj_transcriptome']})
                low = (x.baseMean < 10) or (min(x.frac_nonzero_control, x.frac_nonzero_hypercaloric) < 0.5)
                d['result_class'] = 'not testable (Cook/filter or no estimate)' if pd.isna(x.p) else ('low-count, unstable' if low else 'testable')
                rows.append(d)
    t = pd.DataFrame(rows)
    for panel in PANELS:
        m = t.panel == panel
        t.loc[m, 'padj_family_BH'] = bh(t.loc[m, 'p'])
        pv = t.loc[m, 'p'].values; ok = ~np.isnan(pv); h = np.full_like(pv, np.nan)
        if ok.any(): h[ok] = multipletests(pv[ok], method='holm')[1]
        t.loc[m, 'padj_family_Holm'] = h
    t['family'] = t.panel.map({'primary': 'Pnoc+Oprl1 x 6 prespecified contrasts (BH)',
                               'secondary_processing': 'Pcsk1/Pcsk2/Cpe x 6 contrasts (BH, separate family)',
                               'exploratory_opioid': 'opioid ligands/receptors x 6 contrasts (BH, separate exploratory family)'})
    t = t.merge(sens, left_on='gene', right_index=True, how='left')
    t.loc[t.series != 'GSE281091', [c for c in t.columns if c.endswith('_pseudocount0.5')]] = np.nan
    t.to_csv(OUT/'target_gene_statistics.csv', index=False)
    pd.concat(longs).to_csv(OUT/'target_gene_sample_values.csv', index=False)
    pd.concat(qcs).to_csv(OUT/'qc'/'all_bulk_sample_qc.csv')
    pd.concat(inters).reset_index().to_csv(OUT/'exploratory_interactions.csv', index=False)
    (OUT/'model_notes.json').write_text(json.dumps(notes, indent=2, default=str) + '\n')
    with pd.option_context('display.width', 250, 'display.max_columns', 30):
        print(t[t.panel == 'primary'][['contrast_id', 'gene', 'baseMean', 'n_control', 'n_hypercaloric', 'log2fc', 'se', 'ci_low', 'ci_high', 'p', 'padj_family_BH', 'padj_family_Holm', 'padj_transcriptome', 'result_class']].round(4).to_string())
        print(pd.concat(inters).round(4))
        print(json.dumps(notes, indent=1, default=str)[:3000])

if __name__ == '__main__':
    main()
