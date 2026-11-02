"""Cross-ganglion HFD consensus of feeding-regulatory genes with fasting counter-regulation.

Implements metadata/hfd_consensus_plan.md (committed before estimation).
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from load_data import ROOT
from run_bulk_de import bh

R = ROOT/'results'; OUT = R/'consensus'; OUT.mkdir(parents=True, exist_ok=True)
GMT = ROOT/'data'/'raw'/'msigdb'/'m5.go.bp.v2024.1.Mm.symbols.gmt'
SETS = ['FEEDING_BEHAVIOR', 'ADULT_FEEDING_BEHAVIOR', 'EATING_BEHAVIOR', 'REGULATION_OF_FEEDING_BEHAVIOR',
        'POSITIVE_REGULATION_OF_FEEDING_BEHAVIOR', 'NEGATIVE_REGULATION_OF_FEEDING_BEHAVIOR', 'REGULATION_OF_EATING_BEHAVIOR',
        'POSITIVE_REGULATION_OF_EATING_BEHAVIOR', 'NEGATIVE_REGULATION_OF_EATING_BEHAVIOR', 'REGULATION_OF_APPETITE',
        'REDUCTION_OF_FOOD_INTAKE_IN_RESPONSE_TO_DIETARY_EXCESS', 'ENERGY_HOMEOSTASIS', 'RESPONSE_TO_LEPTIN',
        'CELLULAR_RESPONSE_TO_LEPTIN_STIMULUS', 'LEPTIN_MEDIATED_SIGNALING_PATHWAY']
MIN_MEAN = 10

def universe():
    gs = {}
    for line in GMT.read_text().splitlines():
        f = line.split('\t')
        if f[0].removeprefix('GOBP_') in SETS:
            gs[f[0].removeprefix('GOBP_')] = set(f[2:])
    assert len(gs) == len(SETS), set(SETS) - set(gs)
    panels = pd.read_csv(ROOT/'metadata'/'feeding_gene_panels.csv')
    cur = set(panels[panels.panel.isin(['vagal_feeding', 'nop'])].gene)
    genes = set().union(*gs.values()) | cur
    src = {g: ';'.join(sorted([k for k, v in gs.items() if g in v] + (['curated_panel'] if g in cur else []))) for g in genes}
    return genes, src

def one_sided(lfc, p2):
    """Convert two-sided P to one-sided P for 'up' and 'down'."""
    up = np.where(lfc > 0, p2 / 2, 1 - p2 / 2); return up, 1 - up

def drg_combined():
    a = pd.read_csv(R/'de'/'DRG_HFD_12wk_transcriptome.csv.gz', index_col=0)
    b = pd.read_csv(R/'de'/'DRG_HFD_36wk_transcriptome.csv.gz', index_col=0)
    g = a.index.intersection(b.index)
    a, b = a.loc[g], b.loc[g]
    w1, w2 = 1 / a.se ** 2, 1 / b.se ** 2
    lfc = (w1 * a.log2fc + w2 * b.log2fc) / (w1 + w2); se = 1 / np.sqrt(w1 + w2)
    z = lfc / se
    Q = w1 * (a.log2fc - lfc) ** 2 + w2 * (b.log2fc - lfc) ** 2
    d = pd.DataFrame({'DRG_lfc': lfc, 'DRG_se': se, 'DRG_p_up': stats.norm.sf(z), 'DRG_p_down': stats.norm.cdf(z),
                      'DRG_het_p': stats.chi2.sf(Q, 1), 'DRG_12wk_lfc': a.log2fc, 'DRG_12wk_p': a.p,
                      'DRG_36wk_lfc': b.log2fc, 'DRG_36wk_p': b.p, 'DRG_baseMean': a.baseMean})
    ok = a.se.notna() & b.se.notna()
    return d[ok]

def nodose(contrast, tag, prefix):
    f = R/'feeding'/'de'/(f'{contrast}_{tag}_transcriptome.csv.gz' if tag else f'{contrast}_transcriptome.csv.gz')
    r = pd.read_csv(f).dropna(subset=['gene'])
    r = r.sort_values('baseMean', ascending=False).drop_duplicates('gene').set_index('gene')
    return pd.DataFrame({f'{prefix}_lfc': r.log2fc, f'{prefix}_p_up': stats.t.sf(r.t, r.df), f'{prefix}_p_down': stats.t.cdf(r.t, r.df),
                         f'{prefix}_p': r.p})

def symp_descriptive():
    pb = pd.read_csv(ROOT/'data'/'processed'/'GSE233163'/'pseudobulk_library_celltype.csv.gz', index_col=0)
    lib = pb.index.str.split('|').str[0]
    site = np.where(lib.str.lower().str.startswith('scg'), 'SCG', 'stellate'); run = lib.str.rsplit('_', n=1).str[-1]
    pool = pd.Series(site, index=pb.index) + '_' + pd.Series(run, index=pb.index)
    P = pb.groupby(pool.values).sum()
    cpm = np.log2(P.div(P.sum(axis=1), axis=0) * 1e6 + 1)
    out = {}
    for s in ['SCG', 'stellate']:
        lean = cpm.loc[[f'{s}_N61', f'{s}_N66']]; ob = cpm.loc[f'{s}_N68']
        out[f'{s}_obese_lfc'] = ob - lean.mean()
        out[f'{s}_lean_run_diff'] = lean.loc[f'{s}_N66'] - lean.loc[f'{s}_N61']
        out[f'{s}_log2cpm'] = cpm.loc[[f'{s}_N61', f'{s}_N66', f'{s}_N68']].mean()
    return pd.DataFrame(out)

def supportive():
    t = pd.read_csv(R/'de'/'DRG_Nav18_TRAP_WD_transcriptome.csv.gz', index_col=0)
    t = t[~t.index.duplicated()]
    e = pd.read_csv(R/'de'/'Nodose_earlyHFD_P21_transcriptome.csv.gz', index_col=0)
    return pd.DataFrame({'TRAP_lfc': t.log2fc, 'TRAP_p': t.p, 'TRAP_baseMean': t.baseMean}).join(
        pd.DataFrame({'P21_lfc': e.log2fc, 'P21_p': e.p, 'P21_baseMean': e.baseMean}), how='outer')

def conj(df, cols_up, cols_down):
    up = df[cols_up].max(axis=1); dn = df[cols_down].max(axis=1)
    p = np.minimum(1, 2 * np.minimum(up, dn))
    return p, np.where(up < dn, 'up', 'down')

def analyse(nod_tag, label, U, src):
    drg = drg_combined()
    nod = nodose('Nodose_HFD_4wk_vs_fed', nod_tag, 'NOD')
    fast = nodose('Nodose_fasted18h_vs_fed', nod_tag, 'FAST')
    nod_mean = pd.read_csv(R/'feeding'/'de'/'Nodose_HFD_4wk_vs_fed_transcriptome.csv.gz').dropna(subset=['gene']) \
        .sort_values('baseMean', ascending=False).drop_duplicates('gene').set_index('gene').baseMean.rename('NOD_baseMean')
    df = drg.join(nod, how='inner').join(fast, how='inner').join(nod_mean, how='left')
    df = df[(df.DRG_baseMean >= MIN_MEAN) & (df.NOD_baseMean >= MIN_MEAN)]
    df['p_conj'], df['conj_dir'] = conj(df, ['DRG_p_up', 'NOD_p_up'], ['DRG_p_down', 'NOD_p_down'])
    df['p_counter'], df['counter_dir'] = conj(df, ['DRG_p_up', 'NOD_p_up', 'FAST_p_down'], ['DRG_p_down', 'NOD_p_up' if False else 'NOD_p_down', 'FAST_p_up'])
    df['same_sign_HFD'] = np.sign(df.DRG_lfc) == np.sign(df.NOD_lfc)
    df['fasting_opposite'] = df.same_sign_HFD & (np.sign(df.FAST_lfc) == -np.sign(df.NOD_lfc))
    tw = df.copy()  # transcriptome-wide context
    tw['q_conj_tw'] = bh(tw.p_conj); tw['q_counter_tw'] = bh(tw.p_counter)
    fu = df[df.index.isin(U)].copy()
    fu['q_conj'] = bh(fu.p_conj); fu['q_counter'] = bh(fu.p_counter)
    fu['gene_sets'] = [src[g] for g in fu.index]
    fu['nodose_model'] = label
    return fu.sort_values('p_conj'), tw

def main():
    U, src = universe()
    pd.DataFrame({'gene': sorted(U), 'sources': [src[g] for g in sorted(U)]}).to_csv(OUT/'feeding_gene_universe.csv', index=False)
    res = {}
    for tag, label in [('length_adjusted', 'length-adjusted (primary)'), (None, 'prespecified 24-sample'), ('female_only', 'female only')]:
        fu, tw = analyse(tag, label, U, src)
        res[label] = fu
        if tag == 'length_adjusted':
            tw.to_csv(OUT/'transcriptome_wide_conjunction.csv.gz')
            ctx = dict(n_tw=len(tw), n_tw_q10=int((tw.q_conj_tw < 0.1).sum()), n_tw_counter_q10=int((tw.q_counter_tw < 0.1).sum()))
    prim = res['length-adjusted (primary)']
    for label in ['prespecified 24-sample', 'female only']:
        o = res[label]
        prim[f'p_conj[{label}]'] = o.p_conj.reindex(prim.index); prim[f'NOD_lfc[{label}]'] = o.NOD_lfc.reindex(prim.index)
        prim[f'p_counter[{label}]'] = o.p_counter.reindex(prim.index)
    prim = prim.join(supportive(), how='left').join(symp_descriptive(), how='left')
    prim.index.name = 'gene'
    prim.to_csv(OUT/'feeding_gene_consensus.csv')
    print('universe', len(U), 'testable', len(prim), ctx)
    cols = ['conj_dir', 'DRG_lfc', 'DRG_12wk_lfc', 'DRG_36wk_lfc', 'NOD_lfc', 'FAST_lfc', 'p_conj', 'q_conj', 'p_counter', 'q_counter',
            'p_conj[female only]', 'p_conj[prespecified 24-sample]', 'TRAP_lfc', 'P21_lfc', 'SCG_obese_lfc', 'SCG_lean_run_diff', 'stellate_obese_lfc', 'stellate_lean_run_diff']
    with pd.option_context('display.width', 300, 'display.max_columns', 30):
        print(prim[cols].head(30).round(3).to_string())
        print(prim.sort_values('p_counter')[cols].head(20).round(3).to_string())

if __name__ == '__main__':
    main()
