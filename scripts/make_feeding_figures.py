"""Figures 4-5 and S4 for the feeding-regulation modules (A-C)."""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import TwoSlopeNorm
from matplotlib.ticker import MaxNLocator

sys.path.insert(0, str(Path(__file__).parent))
from figure_style import configure, save, MM, CONTROL, HYPERCALORIC
from load_data import ROOT

R = ROOT/'results'/'feeding'
FIG = ROOT/'figures'; SRC = FIG/'source_data'
GREY = '#6E6E6E'
FASTED, REFED = '#E69F00', '#009E73'   # Okabe-Ito orange, bluish green
STYLE = {'NonFasted': (CONTROL, 'o', 'Chow, fed'), 'Fasted': (FASTED, 's', 'Chow, fasted 18 h'),
         'Refed': (REFED, 'D', 'Fasted + refed 20 min'), 'HFD': (HYPERCALORIC, '^', 'HFD 4 wk, fed')}
ORDER = ['NonFasted', 'Fasted', 'Refed', 'HFD']

def letter(fig, ax, l, dx=0.05, dy=0.03):
    bb = ax.get_position(); fig.text(bb.x0 - dx, bb.y1 + dy, l, weight='bold', fontsize=8, ha='left', va='bottom')

def per_animal(ax, sv, gene, lbi=None, rng=None):
    rng = rng or np.random.default_rng(5)
    for i, grp in enumerate(ORDER):
        d = sv[sv.group == grp]
        col, mk, _ = STYLE[grp]
        jit = (rng.random(len(d)) - 0.5) * 0.4
        fem = (d.sex == 'female').values
        ax.scatter(i + jit[~fem], d[gene].values[~fem], marker=mk, s=11, facecolors=col, edgecolors=col, linewidths=0.6, zorder=3, clip_on=False)
        ax.scatter(i + jit[fem], d[gene].values[fem], marker=mk, s=11, facecolors='white', edgecolors=col, linewidths=0.7, zorder=3, clip_on=False)
        if lbi is not None:
            low = d.matrix_column.map(lbi).values < -0.2
            if low.any(): ax.scatter(i + jit[low], d[gene].values[low], s=34, facecolors='none', edgecolors='black', linewidths=0.4, zorder=4, clip_on=False)
        ax.plot([i - 0.3, i + 0.3], [np.median(d[gene])] * 2, color='black', lw=0.9, zorder=5)
    ax.set_xticks(range(4)); ax.set_xticklabels(['Fed', 'Fast', 'Refed', 'HFD'])
    ax.set_xlim(-0.6, 3.6); ax.set_ylim(0, sv[gene].max() * 1.12)
    ax.set_title(gene, style='italic', pad=3); ax.yaxis.set_major_locator(MaxNLocator(4))

def legend_states(fig, y=0.0, extra=()):
    h = [Line2D([], [], marker=STYLE[g][1], ls='', mfc=STYLE[g][0], mec=STYLE[g][0], ms=3.6, label=STYLE[g][2]) for g in ORDER]
    h += [Line2D([], [], marker='o', ls='', mfc='white', mec=GREY, ms=3.6, label='Open: female')] + list(extra)
    fig.legend(handles=h, loc='lower center', ncol=3, bbox_to_anchor=(0.5, y), columnspacing=1.2, handletextpad=0.3)

# ------------------------------------------------------------------ Fig 4: fasting
def fig4():
    de = pd.read_csv(R/'de'/'Nodose_fasted18h_vs_fed_transcriptome.csv.gz')
    sv = pd.read_csv(R/'GSE281091_panel_sample_values.csv')
    lbi = pd.read_csv(R/'GSE281091_length_bias_index.csv', index_col=0).length_bias_index
    stats_ = pd.read_csv(R/'feeding_panel_statistics.csv')
    sens = pd.read_csv(R/'feeding_panel_sensitivity.csv')
    fig = plt.figure(figsize=(183 * MM, 130 * MM))
    gs = fig.add_gridspec(2, 4, left=0.07, right=0.985, bottom=0.2, top=0.88, wspace=0.62, hspace=0.75, width_ratios=[1.25, 1, 1, 1])
    # a volcano
    ax = fig.add_subplot(gs[:, 0])
    d = de[de.baseMean >= 10].copy(); d['nl'] = -np.log10(d.p)
    sig = d.padj_transcriptome < 0.05
    ax.scatter(d.log2fc[~sig], d.nl[~sig], s=1.2, c='#C8C8C8', lw=0, rasterized=True)
    ax.scatter(d.log2fc[sig & (d.log2fc > 0)], d.nl[sig & (d.log2fc > 0)], s=1.6, c=FASTED, lw=0, rasterized=True)
    ax.scatter(d.log2fc[sig & (d.log2fc < 0)], d.nl[sig & (d.log2fc < 0)], s=1.6, c=CONTROL, lw=0, rasterized=True)
    lab = ['Fkbp5', 'Hif3a', 'Zbtb16', 'Mt2', 'Plin4', 'Cirbp', 'Atf3', 'Sprr1a', 'Ecel1', 'Ucn', 'Nts', 'Cd74', 'Pnoc', 'Oprl1', 'Htr3a', 'Oxtr']
    for g in lab:
        r = d[d.gene == g]
        if len(r) and abs(r.log2fc.iat[0]) < 5.3:
            r = r.iloc[0]
            off = {'Cirbp': (-4, 7), 'Fkbp5': (4, 3), 'Plin4': (4, 0), 'Mt2': (4, -5)}.get(g, (3, 1))
            ax.annotate(g, (r.log2fc, r.nl), xytext=off, textcoords='offset points', fontsize=5.5, style='italic',
                        color='black' if g not in ('Pnoc', 'Oprl1') else HYPERCALORIC)
            ax.scatter([r.log2fc], [r.nl], s=5, facecolors='none', edgecolors='black', lw=0.4)
    ax.set_xlim(-5.4, 4.4); ax.set_xlabel('log$_2$ fold change\n(fasted 18 h / fed)'); ax.set_ylabel('−log$_{10}$ P')
    ax.set_ylim(0, None)
    ax.text(0.98, 0.02, f'{(sig & (d.log2fc > 0)).sum()} up\n{(sig & (d.log2fc < 0)).sum()} down\n(padj < 0.05)', transform=ax.transAxes,
            ha='right', va='bottom', fontsize=5.5)
    letter(fig, ax, 'a', 0.06)
    # b-c Pnoc/Oprl1 per animal; d-g fasting panel genes
    genes = ['Pnoc', 'Oprl1', 'Htr3a', 'Oxtr', 'Cnr1', 'Cartpt']
    pos = [(0, 1), (0, 2), (0, 3), (1, 1), (1, 2), (1, 3)]
    for k, (g, (i, j)) in enumerate(zip(genes, pos)):
        a = fig.add_subplot(gs[i, j]); per_animal(a, sv, g, lbi=lbi)
        r = stats_[(stats_.contrast_id == 'Nodose_fasted18h_vs_fed') & (stats_.gene == g)].iloc[0]
        s = sens[(sens.contrast_id == 'Nodose_fasted18h_vs_fed') & (sens.gene == g)].set_index('sensitivity')
        a.text(0.5, 1.02, f'fast vs fed {r.log2fc:+.2f}, q = {r.padj_family_BH:.2g}\nlength-adj. P = {s.loc["adjusted for length-bias index", "p"]:.2g}',
               transform=a.transAxes, ha='center', va='bottom', fontsize=5.3, color=GREY)
        a.set_title(g, style='italic', pad=17)
        if j == 1: a.set_ylabel('Normalized counts')
        letter(fig, a, 'bcdefg'[k], 0.045, 0.05)
    ring = Line2D([], [], marker='o', ls='', mfc='none', mec='black', mew=0.4, ms=6, label='Ring: low length-bias index')
    legend_states(fig, extra=[ring])
    save(fig, FIG/'Fig4_fasting_nodose'); plt.close(fig)
    d[['gene', 'baseMean', 'log2fc', 'se', 'p', 'padj_transcriptome']].to_csv(SRC/'Fig4a_fasting_volcano_source_data.csv.gz', index=False)
    sv.to_csv(SRC/'Fig4b-g_GSE281091_per_animal_source_data.csv', index=False)

# ------------------------------------------------------------------ Fig 5: diet vs fasting; satiety panel
def fig5():
    lf = pd.read_csv(R/'moduleB_gene_log2fc.csv.gz', index_col=0)
    null = pd.read_csv(R/'moduleB_permutation_null.csv').perm_spearman
    import json
    B = json.loads((R/'feeding_model_notes.json').read_text())['moduleB']
    st = pd.read_csv(R/'feeding_panel_statistics.csv'); sens = pd.read_csv(R/'feeding_panel_sensitivity.csv')
    sv = pd.read_csv(R/'GSE281091_panel_sample_values.csv')
    lbi = pd.read_csv(R/'GSE281091_length_bias_index.csv', index_col=0).length_bias_index
    fig = plt.figure(figsize=(183 * MM, 158 * MM))
    gs = fig.add_gridspec(2, 3, left=0.075, right=0.985, bottom=0.17, top=0.95, wspace=0.55, hspace=0.42,
                          width_ratios=[1, 1, 1.6], height_ratios=[1, 1.35])
    # a concordance hexbin
    ax = fig.add_subplot(gs[0, 0])
    ax.hexbin(lf.log2fc_fasted_vs_fed, lf.log2fc_HFD_vs_fed, gridsize=45, extent=(-2, 2, -2, 2), bins='log', cmap='Greys', mincnt=1, lw=0.1)
    ax.axhline(0, color=GREY, lw=0.4); ax.axvline(0, color=GREY, lw=0.4)
    ax.set_xlim(-2, 2); ax.set_ylim(-2, 2)
    ax.set_xlabel('log$_2$FC fasted / fed'); ax.set_ylabel('log$_2$FC HFD / fed')
    ax.text(0.03, 0.97, f'Spearman ρ = {B["spearman_observed"]:.2f}\n{B["n_genes"]:,} genes', transform=ax.transAxes, va='top', fontsize=5.8)
    letter(fig, ax, 'a', 0.06)
    # b permutation null
    ax = fig.add_subplot(gs[0, 1])
    ax.hist(null, bins=30, color='#C8C8C8', edgecolor='white', lw=0.3)
    ax.axvline(B['spearman_observed'], color=HYPERCALORIC, lw=1.0)
    ax.set_xlabel('Spearman ρ under relabelling'); ax.set_ylabel('Permutations')
    ax.text(0.97, 0.97, f'observed {B["spearman_observed"]:.2f}\nP = {B["permutation_p_one_sided"]:.2f}\n(400 within-sex\nrelabellings)', transform=ax.transAxes,
            ha='right', va='top', fontsize=5.8)
    letter(fig, ax, 'b', 0.06)
    # c heatmap of panel
    ax = fig.add_subplot(gs[:, 2])
    cols = [('Nodose_fasted18h_vs_fed', None, 'Fasted\nadult'), ('Nodose_HFD_4wk_vs_fed', None, 'HFD\nadult'),
            ('Nodose_HFD_4wk_vs_fed', 'female only', 'HFD\nfemales'), ('Nodose_HFD_4wk_vs_fed', 'adjusted for length-bias index', 'HFD\nlen.-adj.'),
            ('Nodose_earlyHFD_P21', None, 'Early\nP21'), ('Nodose_earlyHFD_P35', None, 'Early\nP35')]
    genes = [g for g in pd.read_csv(ROOT/'metadata'/'feeding_gene_panels.csv').query('panel=="vagal_feeding"').gene]
    M = np.full((len(genes), len(cols)), np.nan); P = np.full_like(M, np.nan); Q = np.full_like(M, np.nan)
    src = []
    for j, (cid, sname, lab) in enumerate(cols):
        for i, g in enumerate(genes):
            if sname is None:
                r = st[(st.contrast_id == cid) & (st.gene == g) & (st.family.isin(['A_panel', 'C_panel']))]
                if len(r) and r.status.iat[0] == 'testable':
                    M[i, j], P[i, j], Q[i, j] = r.log2fc.iat[0], r.p.iat[0], r.padj_family_BH.iat[0]
            else:
                r = sens[(sens.contrast_id == cid) & (sens.gene == g) & (sens.sensitivity == sname)]
                if len(r):
                    M[i, j], P[i, j] = r.log2fc.iat[0], r.p.iat[0]
            src.append(dict(column=lab.replace('\n', ' '), contrast_id=cid, sensitivity=sname or 'prespecified', gene=g,
                            log2fc=M[i, j], p=P[i, j], padj_family_BH=Q[i, j]))
    im = ax.imshow(M, cmap='RdBu_r', norm=TwoSlopeNorm(0, -1.2, 1.2), aspect='auto')
    for i in range(len(genes)):
        for j in range(len(cols)):
            if np.isnan(M[i, j]):
                ax.text(j, i, '×', ha='center', va='center', fontsize=5, color=GREY)
            elif not np.isnan(Q[i, j]) and Q[i, j] < 0.05:
                ax.text(j, i, '●', ha='center', va='center', fontsize=5.5, color='black')
            elif P[i, j] < 0.05:
                ax.text(j, i, '○', ha='center', va='center', fontsize=5.5, color='black')
    ax.set_yticks(range(len(genes))); ax.set_yticklabels(genes, style='italic', fontsize=6)
    ax.set_xticks(range(len(cols))); ax.set_xticklabels([c[2] for c in cols], fontsize=6, rotation=0)
    ax.xaxis.tick_top(); ax.tick_params(length=0)
    for s in ax.spines.values(): s.set_visible(False)
    ax.axvline(1.5, color='white', lw=2); ax.axvline(3.5, color='white', lw=2)
    cb = fig.colorbar(im, ax=ax, orientation='horizontal', fraction=0.03, pad=0.015, aspect=30)
    cb.set_label('log$_2$ fold change'); cb.outline.set_linewidth(0.4); cb.ax.tick_params(labelsize=6)
    bb = ax.get_position(); fig.text(bb.x0 - 0.06, 0.955, 'c', weight='bold', fontsize=8)
    # d-f robust HFD genes per animal
    for k, g in enumerate(['Uts2b', 'Vip', 'Lepr']):
        sub = gs[1, :2].subgridspec(1, 3, wspace=0.75)
        a = fig.add_subplot(sub[k]); per_animal(a, sv, g, lbi=lbi)
        if k == 0: a.set_ylabel('Normalized counts')
        r = st[(st.contrast_id == 'Nodose_HFD_4wk_vs_fed') & (st.gene == g)].iloc[0]
        a.text(0.5, 1.02, f'HFD {r.log2fc:+.2f}, q = {r.padj_family_BH:.2g}', transform=a.transAxes, ha='center', va='bottom', fontsize=5.3, color=GREY)
        a.set_title(g, style='italic', pad=10)
        a.tick_params(axis='x', labelsize=6); plt.setp(a.get_xticklabels(), rotation=45, ha='right', rotation_mode='anchor')
        letter(fig, a, 'def'[k], 0.05, 0.035)
    ring = Line2D([], [], marker='o', ls='', mfc='none', mec='black', mew=0.4, ms=6, label='Ring: low length-bias index')
    fig.legend(handles=[Line2D([], [], marker=STYLE[g][1], ls='', mfc=STYLE[g][0], mec=STYLE[g][0], ms=3.6, label=STYLE[g][2]) for g in ORDER]
               + [Line2D([], [], marker='o', ls='', mfc='white', mec=GREY, ms=3.6, label='Open: female'), ring],
               loc='lower left', ncol=2, bbox_to_anchor=(0.03, 0.0), columnspacing=1.2, handletextpad=0.3)
    cb.ax.text(0.5, -6.3, '● family q < 0.05   ○ nominal P < 0.05\n× not testable / absent', transform=cb.ax.transAxes, ha='center', va='top', fontsize=5.8)
    save(fig, FIG/'Fig5_diet_vs_fasting_satiety_panel'); plt.close(fig)
    pd.DataFrame(src).to_csv(SRC/'Fig5c_panel_heatmap_source_data.csv', index=False)
    lf.to_csv(SRC/'Fig5a_concordance_source_data.csv.gz'); null.to_csv(SRC/'Fig5b_permutation_null_source_data.csv', index=False)

# ------------------------------------------------------------------ Fig S4: technical axis
def figs4():
    lbi = pd.read_csv(R/'GSE281091_length_bias_index.csv', index_col=0)
    pca = pd.read_csv(R/'GSE281091_24sample_pca.csv', index_col=0).join(lbi[['length_bias_index']])
    fig, axes = plt.subplots(1, 2, figsize=(120 * MM, 66 * MM))
    fig.subplots_adjust(left=0.11, right=0.98, bottom=0.42, top=0.9, wspace=0.45)
    ax = axes[0]
    for grp in ORDER:
        d = pca[pca.group == grp]; col, mk, _ = STYLE[grp]
        for fem in [False, True]:
            x = d[(d.sex == 'female') == fem]
            ax.scatter(x.PC1, x.PC2, marker=mk, s=12, facecolors='white' if fem else col, edgecolors=col, linewidths=0.7)
    ax.set_xlabel(f'PC1 ({100 * pca.pc1_var.iat[0]:.0f}%)'); ax.set_ylabel(f'PC2 ({100 * pca.pc2_var.iat[0]:.0f}%)')
    letter(fig, ax, 'a', 0.08)
    ax = axes[1]
    for grp in ORDER:
        d = pca[pca.group == grp]; col, mk, _ = STYLE[grp]
        for fem in [False, True]:
            x = d[(d.sex == 'female') == fem]
            ax.scatter(x.length_bias_index, x.PC1, marker=mk, s=12, facecolors='white' if fem else col, edgecolors=col, linewidths=0.7)
    ax.set_xlabel('Length-bias index\n(ρ, expression vs log gene length)'); ax.set_ylabel('PC1')
    letter(fig, ax, 'b', 0.08)
    legend_states(fig)
    save(fig, FIG/'FigS4_GSE281091_technical_axis'); plt.close(fig)
    pca.to_csv(SRC/'FigS4_source_data.csv')

if __name__ == '__main__':
    configure(); fig4(); fig5(); figs4(); print('done')
