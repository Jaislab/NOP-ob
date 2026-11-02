"""Modules A-C: fasting/refeeding and diet effects on vagal feeding-regulation genes.

Follows metadata/feeding_analysis_rules.md and metadata/feeding_contrast_inventory.csv
(committed before estimation).
"""
from pathlib import Path
from itertools import combinations
import json, sys
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from load_data import ROOT, RAW
from run_bulk_de import fit_fdist_trend, bh

OUT = ROOT/'results'/'feeding'
(OUT/'de').mkdir(parents=True, exist_ok=True)
PANELS = pd.read_csv(ROOT/'metadata'/'feeding_gene_panels.csv')
MIN_MEAN = 5

# ------------------------------------------------------------------ data
def load_combined():
    fa = RAW/'GSE281091'/'GSE281091_Partek_Bulk_RNA-Seq_Exp_3A_Normalization_Normalized_counts_metabolic_states.txt.gz'
    fb = RAW/'GSE281091'/'GSE281091_Partek_Bulk_RNA-Seq_Exp_3B_Normalization_Normalized_counts_diet_types.txt.gz'
    a = pd.read_csv(fa, sep='\t').set_index('gene_id'); b = pd.read_csv(fb, sep='\t').set_index('gene_id')
    sa = [c for c in a if c.startswith('YJB')]; sb = [c for c in b if c.startswith('YJB')]
    shared = sorted(set(sa) & set(sb))
    genes = a.index.intersection(b.index)
    A, B = a.loc[genes, sa], b.loc[genes, sb]
    ok = (A[shared] > 10).all(axis=1) & (B[shared] > 10).all(axis=1)
    factors = {s: float(np.median(B.loc[ok, s] / A.loc[ok, s])) for s in shared}
    f = float(np.median(list(factors.values())))
    A = A * f
    X = pd.concat([B, A[[s for s in sa if s not in shared]]], axis=1)
    sym = b.loc[genes, 'gene_name'].fillna(a.loc[genes, 'gene_name'])
    m = pd.read_csv(ROOT/'metadata'/'sample_manifest.csv')
    m = m[m.series == 'GSE281091'].set_index('matrix_column')
    arm = m.diet_detail.map({'Non-fasted_regular chow': 'NonFasted', 'Fasted_regular chow': 'Fasted',
                             'Feeding_regular chow': 'Refed', 'Non-fasted_high fat diet': 'HFD'})
    meta = pd.DataFrame({'gsm': m.gsm, 'group': arm, 'sex': m.sex, 'title': m.title}).loc[X.columns]
    notes = dict(n_genes_intersection=int(len(genes)), rescale_factor_3A_to_3B=f, per_shared_sample_factors=factors,
                 rescaled_shared_max_abs_log2_dev=float(np.nanmax(np.abs(np.log2((A.loc[ok, shared] + 1) / (B.loc[ok, shared] + 1)).values))))
    return X.T, meta, sym, notes

# ------------------------------------------------------------------ model
def moderated(Y, design, contrasts):
    Xv = design.values; XtXi = np.linalg.inv(Xv.T @ Xv)
    Bc = XtXi @ Xv.T @ Y.values
    res = Y.values - Xv @ Bc
    dfr = Xv.shape[0] - Xv.shape[1]
    s2 = (res ** 2).sum(axis=0) / dfr
    am = Y.values.mean(axis=0)
    ok = s2 > 1e-12
    d0, s0 = fit_fdist_trend(s2[ok], np.full(ok.sum(), dfr, float), am[ok])
    s0f = np.full(len(s2), np.nan); s0f[ok] = s0
    post = (d0 * s0f + dfr * s2) / (d0 + dfr) if np.isfinite(d0) else s0f
    dft = dfr + d0 if np.isfinite(d0) else np.inf
    q = stats.t.ppf(0.975, dft)
    out = {}
    for name, c in contrasts.items():
        c = np.asarray(c, float)
        est = c @ Bc; se = np.sqrt(post * (c @ XtXi @ c))
        t = est / se; p = 2 * stats.t.sf(np.abs(t), dft)
        out[name] = pd.DataFrame({'log2fc': est, 'se': se, 't': t, 'p': p, 'ci_low': est - q * se, 'ci_high': est + q * se,
                                  'df': dft}, index=Y.columns)
    return out, dict(prior_df=float(d0), residual_df=int(dfr))

# ------------------------------------------------------------------ main
def main():
    X, meta, sym, notes = load_combined()
    X = X.loc[:, X.sum(axis=0) > 0]
    Y = np.log2(X + 1)
    g = meta.group
    design = pd.DataFrame({'int': 1.0, 'Fasted': (g == 'Fasted') * 1.0, 'Refed': (g == 'Refed') * 1.0,
                           'HFD': (g == 'HFD') * 1.0, 'female': (meta.sex == 'female') * 1.0}, index=meta.index)
    C = {'Nodose_fasted18h_vs_fed': [0, 1, 0, 0, 0], 'Nodose_refed20min_vs_fasted': [0, -1, 1, 0, 0],
         'Nodose_refed20min_vs_fed': [0, 0, 1, 0, 0], 'Nodose_HFD_4wk_vs_fed': [0, 0, 0, 1, 0]}
    res, mod = moderated(Y, design, C)
    notes.update(mod)
    groups = {'Nodose_fasted18h_vs_fed': ('Fasted', 'NonFasted'), 'Nodose_refed20min_vs_fasted': ('Refed', 'Fasted'),
              'Nodose_refed20min_vs_fed': ('Refed', 'NonFasted'), 'Nodose_HFD_4wk_vs_fed': ('HFD', 'NonFasted')}
    for cid, r in res.items():
        num, den = groups[cid]
        r.insert(0, 'gene', sym.reindex(r.index).values)
        r['mean_norm_num'] = X.loc[g == num].mean().values; r['mean_norm_den'] = X.loc[g == den].mean().values
        r['baseMean'] = X.loc[g.isin([num, den])].mean().values
        r['padj_transcriptome'] = bh(r.p)
        r.index.name = 'gene_id'
        r.to_csv(OUT/'de'/f'{cid}_transcriptome.csv.gz')
        notes[f'{cid}_n_padj05'] = int((r.padj_transcriptome < 0.05).sum())
    # per-sample values for panel genes
    pg = PANELS.gene.unique()
    ids = {gg: sym.index[sym == gg].intersection(X.columns) for gg in pg}
    dup = {gg: len(v) for gg, v in ids.items() if len(v) > 1}
    assert not dup, f'duplicated panel genes {dup}'
    cols = {gg: v[0] for gg, v in ids.items() if len(v) == 1}
    sv = X[list(cols.values())].rename(columns={v: k for k, v in cols.items()})
    sv = sv.join(meta).reset_index(names='matrix_column')
    sv.to_csv(OUT/'GSE281091_panel_sample_values.csv', index=False)

    # ---------------- panel tables (modules A and C)
    rows = []
    def add(cid, family, panel, r, gene, src):
        d = dict(contrast_id=cid, family=family, panel=panel, gene=gene, source=src)
        hit = r[r.gene == gene] if 'gene' in r else r.loc[[gene]] if gene in r.index else r.iloc[0:0]
        if hit.empty:
            d['status'] = 'not in deposited matrix'; rows.append(d); return
        h = hit.iloc[0]
        d.update(log2fc=h.log2fc, se=h.se, ci_low=h.ci_low, ci_high=h.ci_high, p=h.p, baseMean=h.baseMean,
                 padj_transcriptome=h.padj_transcriptome)
        d['status'] = 'low abundance (mean <5): not testable' if h.baseMean < MIN_MEAN else 'testable'
        rows.append(d)
    gp = lambda p: PANELS[PANELS.panel == p].gene.tolist()
    for cid in ['Nodose_fasted18h_vs_fed', 'Nodose_refed20min_vs_fasted']:
        for gene in gp('nop'): add(cid, 'A_nop', 'nop', res[cid], gene, 'GSE281091')
        for gene in gp('vagal_feeding'): add(cid, 'A_panel', 'vagal_feeding', res[cid], gene, 'GSE281091')
    for gene in gp('ieg'): add('Nodose_refed20min_vs_fasted', 'A_ieg', 'ieg', res['Nodose_refed20min_vs_fasted'], gene, 'GSE281091')
    for gene in gp('nop') + gp('vagal_feeding') + gp('ieg'):
        add('Nodose_refed20min_vs_fed', 'secondary', PANELS.set_index('gene').panel.get(gene), res['Nodose_refed20min_vs_fed'], gene, 'GSE281091')
    for gene in gp('vagal_feeding'): add('Nodose_HFD_4wk_vs_fed', 'C_panel', 'vagal_feeding', res['Nodose_HFD_4wk_vs_fed'], gene, 'GSE281091')
    for cid in ['Nodose_earlyHFD_P21', 'Nodose_earlyHFD_P35']:
        r = pd.read_csv(ROOT/'results'/'de'/f'{cid}_transcriptome.csv.gz', index_col=0)
        for gene in gp('vagal_feeding'): add(cid, 'C_panel', 'vagal_feeding', r, gene, 'GSE320302')
    t = pd.DataFrame(rows)
    t['padj_family_BH'] = np.nan
    for fam in ['A_nop', 'A_panel', 'A_ieg', 'C_panel']:
        k = (t.family == fam) & (t.status == 'testable')
        t.loc[k, 'padj_family_BH'] = bh(t.loc[k, 'p'])
    pred = PANELS.set_index('gene')
    def agree(r):
        col = 'predicted_fasting' if r.contrast_id == 'Nodose_fasted18h_vs_fed' else ('predicted_HFD' if r.family == 'C_panel' else None)
        if col is None or r.gene not in pred.index or pd.isna(pred.loc[r.gene, col]) or pd.isna(r.get('log2fc')):
            return ''
        want = pred.loc[r.gene, col]
        return want + (' (agrees)' if (r.log2fc > 0) == (want == 'up') else ' (opposite)')
    t['literature_prediction'] = t.apply(agree, axis=1)
    t.to_csv(OUT/'feeding_panel_statistics.csv', index=False)

    # ---------------- module B: concordance with exact within-sex permutation
    expr_ids = X.columns[X.loc[g.isin(['NonFasted', 'HFD', 'Fasted'])].mean() >= 10]
    Ye = Y[expr_ids]
    def corr_for(hfd_idx, nf_idx, fast_idx):
        l1 = Ye.loc[hfd_idx].mean() - Ye.loc[nf_idx].mean(); l2 = Ye.loc[fast_idx].mean() - Ye.loc[nf_idx].mean()
        return stats.spearmanr(l1, l2)[0]
    pool = meta[g.isin(['NonFasted', 'HFD'])]
    fast = meta.index[g == 'Fasted']
    obs = corr_for(meta.index[g == 'HFD'], meta.index[g == 'NonFasted'], fast)
    m_ids = pool.index[pool.sex == 'male']; f_ids = pool.index[pool.sex == 'female']
    perm = []
    for hm in combinations(m_ids, 3):
        for hf in combinations(f_ids, 3):
            h = list(hm) + list(hf); n = [s for s in pool.index if s not in h]
            perm.append(corr_for(h, n, fast))
    perm = np.array(perm)
    p_perm = float(np.mean(perm >= obs - 1e-12))
    # female-only sensitivity
    fem = meta.sex == 'female'
    obs_f = corr_for(meta.index[(g == 'HFD') & fem], meta.index[(g == 'NonFasted') & fem], meta.index[(g == 'Fasted') & fem])
    perm_f = []
    for hf in combinations(f_ids, 3):
        n = [s for s in f_ids if s not in hf]
        perm_f.append(corr_for(list(hf), n, meta.index[(g == 'Fasted') & fem]))
    perm_f = np.array(perm_f)
    # split-control descriptive estimate (disjoint controls)
    nf = list(meta.index[g == 'NonFasted']); split = []
    for c1 in combinations(nf, 3):
        c2 = [s for s in nf if s not in c1]
        l1 = Ye.loc[meta.index[g == 'HFD']].mean() - Ye.loc[list(c1)].mean(); l2 = Ye.loc[fast].mean() - Ye.loc[c2].mean()
        split.append(stats.spearmanr(l1, l2)[0])
    # also the refed (18 h fasted + 20 min food) arm as a second fasting replicate vs HFD
    obs_refed = corr_for(meta.index[g == 'HFD'], meta.index[g == 'NonFasted'], meta.index[g == 'Refed'])
    B = dict(n_genes=int(len(expr_ids)), spearman_observed=obs, permutation_p_one_sided=p_perm, n_permutations=int(len(perm)),
             perm_null_mean=float(perm.mean()), perm_null_95pct=float(np.quantile(perm, 0.95)),
             female_only_observed=obs_f, female_only_p=float(np.mean(perm_f >= obs_f - 1e-12)), female_only_n_perm=int(len(perm_f)),
             split_control_median=float(np.median(split)), split_control_range=[float(min(split)), float(max(split))],
             spearman_with_refed_arm=obs_refed)
    pd.DataFrame({'perm_spearman': perm}).to_csv(OUT/'moduleB_permutation_null.csv', index=False)
    lf = pd.DataFrame({'gene': sym.reindex(expr_ids).values,
                       'log2fc_HFD_vs_fed': (Ye.loc[g == 'HFD'].mean() - Ye.loc[g == 'NonFasted'].mean()).values,
                       'log2fc_fasted_vs_fed': (Ye.loc[g == 'Fasted'].mean() - Ye.loc[g == 'NonFasted'].mean()).values}, index=expr_ids)
    lf.to_csv(OUT/'moduleB_gene_log2fc.csv.gz')
    notes['moduleB'] = B
    # QC PCA on 24 samples
    v = Ye.loc[:, Ye.var().sort_values().index[-500:]]; v = v - v.mean()
    u, s, _ = np.linalg.svd(v.values, full_matrices=False)
    pd.DataFrame({'PC1': u[:, 0] * s[0], 'PC2': u[:, 1] * s[1], 'pc1_var': s[0] ** 2 / (s ** 2).sum(), 'pc2_var': s[1] ** 2 / (s ** 2).sum()},
                 index=meta.index).join(meta).to_csv(OUT/'GSE281091_24sample_pca.csv')
    (OUT/'feeding_model_notes.json').write_text(json.dumps(notes, indent=2, default=float) + '\n')
    with pd.option_context('display.width', 250, 'display.max_rows', 300):
        print(json.dumps({k: v for k, v in notes.items() if k != 'per_shared_sample_factors'}, indent=1, default=float))
        print(t[t.family != 'secondary'][['contrast_id', 'family', 'gene', 'baseMean', 'log2fc', 'ci_low', 'ci_high', 'p', 'padj_family_BH', 'status', 'literature_prediction']].round(4).to_string())

if __name__ == '__main__' and 'sens' not in sys.argv:
    main()
    length_sensitivity()

# ------------------------------------------------------------------ technical-axis sensitivity (added after QC showed PC1 = gene-length bias)
def length_sensitivity():
    """PC1 in GSE281091 tracks gene length (RNA-integrity/3'-bias signature) and isolates the HFD males.
    Sensitivities for every panel contrast: (i) female-only model; (ii) model adjusted for a per-sample
    length-bias index (Spearman, across expressed genes, of the sample's centred log expression vs log10 gene length)."""
    X, meta, sym, _ = load_combined()
    X = X.loc[:, X.sum(axis=0) > 0]
    Y = np.log2(X + 1)
    ann = pd.read_csv(RAW/'GSE281091'/'GSE281091_Partek_Bulk_RNA-Seq_Exp_3B_Normalization_Normalized_counts_diet_types.txt.gz',
                      sep='\t', usecols=['gene_id', 'Length']).set_index('gene_id').Length
    expr = X.columns[(X.mean() >= 10) & X.columns.isin(ann.index)]
    Z = Y[expr] - Y[expr].mean()
    L = np.log10(ann.loc[expr].values)
    meta = meta.copy()
    meta['length_bias_index'] = [stats.spearmanr(Z.loc[s].values, L)[0] for s in meta.index]
    g = meta.group
    base = {'int': 1.0, 'Fasted': (g == 'Fasted') * 1.0, 'Refed': (g == 'Refed') * 1.0, 'HFD': (g == 'HFD') * 1.0}
    C = {'Nodose_fasted18h_vs_fed': {'Fasted': 1}, 'Nodose_refed20min_vs_fasted': {'Refed': 1, 'Fasted': -1},
         'Nodose_HFD_4wk_vs_fed': {'HFD': 1}}
    out = []
    for label, rows, extra in [('female only', meta.sex == 'female', {}),
                               ('adjusted for length-bias index', meta.sex.notna(), {'lbi': 'length_bias_index'})]:
        mm = meta[rows]; d = pd.DataFrame({k: (v[rows] if hasattr(v, '__len__') else v) for k, v in base.items()}, index=mm.index)
        if label != 'female only':
            d['female'] = (mm.sex == 'female') * 1.0
        for k, col in extra.items():
            d[k] = mm[col]
        cons = {cid: [c.get(col, 0) for col in d.columns] for cid, c in C.items()}
        res, mod = moderated(Y.loc[mm.index], d, cons)
        for cid, r in res.items():
            r = r.assign(gene=sym.reindex(r.index).values)
            n05 = int((bh(r.p) < 0.05).sum())
            tag = 'female_only' if label == 'female only' else 'length_adjusted'
            rr = r.assign(baseMean=X.loc[mm.index, r.index].mean().values, padj_transcriptome=bh(r.p))
            rr.index.name = 'gene_id'
            rr.to_csv(OUT/'de'/f'{cid}_{tag}_transcriptome.csv.gz')
            for gene in PANELS.gene.unique():
                h = r[r.gene == gene]
                if len(h):
                    h = h.iloc[0]
                    out.append(dict(contrast_id=cid, sensitivity=label, gene=gene, log2fc=h.log2fc, ci_low=h.ci_low,
                                    ci_high=h.ci_high, p=h.p, prior_df=mod['prior_df'], residual_df=mod['residual_df'],
                                    n_padj05_transcriptome=n05))
    s = pd.DataFrame(out)
    s.to_csv(OUT/'feeding_panel_sensitivity.csv', index=False)
    meta[['group', 'sex', 'length_bias_index']].to_csv(OUT/'GSE281091_length_bias_index.csv')
    # BH within the prespecified families, applied to each sensitivity model (same testable genes)
    t = pd.read_csv(OUT/'feeding_panel_statistics.csv')
    testable = set(zip(t.loc[t.status == 'testable', 'contrast_id'], t.loc[t.status == 'testable', 'gene']))
    fams = {'A_nop': (['Nodose_fasted18h_vs_fed', 'Nodose_refed20min_vs_fasted'], 'nop'),
            'A_panel': (['Nodose_fasted18h_vs_fed', 'Nodose_refed20min_vs_fasted'], 'vagal_feeding'),
            'A_ieg': (['Nodose_refed20min_vs_fasted'], 'ieg'), 'C_panel_adultHFD': (['Nodose_HFD_4wk_vs_fed'], 'vagal_feeding')}
    q = []
    for sens in s.sensitivity.unique():
        for fam, (cids, pan) in fams.items():
            d = s[(s.sensitivity == sens) & s.contrast_id.isin(cids) & s.gene.isin(PANELS[PANELS.panel == pan].gene)].copy()
            d = d[[(c, gg) in testable for c, gg in zip(d.contrast_id, d.gene)]]
            d['q_within_family'] = bh(d.p); d['family'] = fam; q.append(d)
    pd.concat(q).to_csv(OUT/'feeding_panel_sensitivity_familyq.csv', index=False)
    return s, meta

if __name__ == '__main__' and 'sens' in sys.argv:
    s, m = length_sensitivity()
    print(m.sort_values('length_bias_index')[['group', 'sex', 'length_bias_index']].round(3).to_string())
    print(s.groupby(['contrast_id', 'sensitivity']).n_padj05_transcriptome.first())
