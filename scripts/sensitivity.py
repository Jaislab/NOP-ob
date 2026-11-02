"""Robustness checks for the primary Pnoc/Oprl1 tests (reported, never used to select results).

1. Count series: leave-one-biological-unit-out NB GLM for each target with size factors and the
   transcriptome-informed (MAP) dispersion fixed from the full PyDESeq2 fit; plus exact
   Mann-Whitney U on size-factor-normalized counts.
2. GSE281091: sex-stratified OLS on log2(normalized+1); and refit excluding the three samples that
   define the technical PC1 axis (identified from QC on the top-variance transcriptome only).
"""
from pathlib import Path
import sys, warnings
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from load_data import PANELS, ROOT
warnings.filterwarnings('ignore')
R = ROOT/'results'
vals = pd.read_csv(R/'target_gene_sample_values.csv')
inv = pd.read_csv(ROOT/'metadata'/'contrast_inventory.csv').set_index('contrast_id')
rows = []

def nb_fit(y, sf, x, alpha):
    X = sm.add_constant(x.astype(float))
    m = sm.GLM(y, X, family=sm.families.NegativeBinomial(alpha=alpha), offset=np.log(sf)).fit()
    return m.params.iloc[1] / np.log(2), m.bse.iloc[1] / np.log(2), m.pvalues.iloc[1]

for cid, c in inv.iterrows():
    if c.series not in ('GSE225224', 'GSE320302'):
        continue
    disp = pd.read_csv(R/'qc'/f'{c.series}_target_dispersions.csv', index_col=0)
    for g in PANELS['primary']:
        d = vals[(vals.series == c.series) & (vals.gene == g) & vals.group.isin([c.numerator, c.denominator])].reset_index(drop=True)
        x = (d.group == c.numerator).astype(int)
        a = disp.loc[g, 'dispersions']
        lfc, se, p = nb_fit(d.raw_count.values, d.size_factor.values, x, a)
        u = stats.mannwhitneyu(d.normalized[x == 1], d.normalized[x == 0], alternative='two-sided', method='exact' if not d.normalized.duplicated().any() else 'asymptotic')
        loo = []
        for i in range(len(d)):
            k = np.arange(len(d)) != i
            loo.append((d.gsm[i], *nb_fit(d.raw_count.values[k], d.size_factor.values[k], x[k], a)))
        loo = pd.DataFrame(loo, columns=['dropped', 'lfc', 'se', 'p'])
        rows.append(dict(contrast_id=cid, gene=g, check='single-gene NB GLM, fixed MAP dispersion (all samples)', log2fc=lfc, se=se, p=p))
        rows.append(dict(contrast_id=cid, gene=g, check='Mann-Whitney U on normalized counts', log2fc=np.nan, se=np.nan, p=u.pvalue))
        rows.append(dict(contrast_id=cid, gene=g, check='leave-one-out: range of log2FC', log2fc=np.nan, se=np.nan, p=np.nan,
                         loo_lfc_min=loo.lfc.min(), loo_lfc_max=loo.lfc.max(), loo_p_max=loo.p.max(),
                         loo_most_influential=loo.loc[(loo.lfc - lfc).abs().idxmax(), 'dropped']))
# TRAP: targets absent -> nothing to test
# GSE281091
qc = pd.read_csv(R/'qc'/'GSE281091_sample_qc.csv', index_col=0)
tech = qc.sort_values('PC1').index[:3].tolist()
for g in PANELS['primary']:
    d = vals[(vals.series == 'GSE281091') & (vals.gene == g)].set_index('gsm')
    y = np.log2(d.normalized + 1)
    for label, keep in [('female only', d.sex == 'female'), ('male only', d.sex == 'male'),
                        ('excluding technical-PC1 samples ' + ','.join(tech), ~d.index.isin(tech))]:
        dd = d[keep]; yy = y[keep]
        X = pd.DataFrame({'hfd': (dd.diet == 'Hypercaloric').astype(float)})
        if dd.sex.nunique() > 1:
            X['female'] = (dd.sex == 'female').astype(float)
        m = sm.OLS(yy, sm.add_constant(X)).fit()
        rows.append(dict(contrast_id='Nodose_HFD_4wk', gene=g, check='OLS log2(norm+1), ' + label,
                         log2fc=m.params['hfd'], se=m.bse['hfd'], p=m.pvalues['hfd'], n=len(dd)))
out = pd.DataFrame(rows)
out.to_csv(R/'sensitivity_primary.csv', index=False)
with pd.option_context('display.width', 250, 'display.max_columns', 20, 'display.max_colwidth', 70):
    print(out.round(4).to_string())
