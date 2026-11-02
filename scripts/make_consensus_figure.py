"""Figure 6: cross-ganglion HFD consensus of feeding-regulatory genes and fasting counter-regulation."""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).parent))
from figure_style import configure, save, MM, HYPERCALORIC, CONTROL
from load_data import ROOT

R = ROOT/'results'/'consensus'; FIG = ROOT/'figures'; SRC = FIG/'source_data'
GREY = '#6E6E6E'; FASTED = '#E69F00'

def main():
    configure()
    fu = pd.read_csv(R/'feeding_gene_consensus.csv', index_col=0)
    tw = pd.read_csv(R/'transcriptome_wide_conjunction.csv.gz', index_col=0)
    fig = plt.figure(figsize=(183 * MM, 112 * MM))
    gs = fig.add_gridspec(1, 3, left=0.085, right=0.985, bottom=0.24, top=0.9, wspace=0.42, width_ratios=[1, 1.75, 0.62])
    # a: DRG vs nodose HFD effects for feeding genes, colour by fasting direction relative to nodose HFD
    ax = fig.add_subplot(gs[0])
    ax.scatter(tw.DRG_lfc, tw.NOD_lfc, s=1, c='#E3E3E3', lw=0, rasterized=True, zorder=1)
    opp = np.sign(fu.FAST_lfc) == -np.sign(fu.NOD_lfc)
    ax.scatter(fu.DRG_lfc[~opp], fu.NOD_lfc[~opp], s=9, facecolors='white', edgecolors=GREY, lw=0.5, zorder=2)
    ax.scatter(fu.DRG_lfc[opp], fu.NOD_lfc[opp], s=9, facecolors=FASTED, edgecolors=FASTED, lw=0.5, zorder=3)
    ax.axhline(0, color=GREY, lw=0.4); ax.axvline(0, color=GREY, lw=0.4)
    ax.set_xlim(-2, 2); ax.set_ylim(-1.5, 1.5)
    for g in ['Pparg', 'Ugcg', 'Insr', 'Cnr1', 'Acacb', 'Lepr', 'Uts2b', 'Pnoc']:
        if g in fu.index and abs(fu.loc[g, 'DRG_lfc']) < 2:
            ax.annotate(g, (fu.loc[g, 'DRG_lfc'], fu.loc[g, 'NOD_lfc']), xytext=(2, 2), textcoords='offset points', fontsize=5.3, style='italic')
    ax.set_xlabel('DRG log$_2$FC, HFD/control\n(12 + 36 wk combined)'); ax.set_ylabel('Nodose log$_2$FC, HFD/fed\n(length-adjusted)')
    ax.text(0.02, 0.98, f'{len(fu)} feeding genes\n0 with q < 0.1', transform=ax.transAxes, va='top', fontsize=5.8)
    ax.legend(handles=[Line2D([], [], marker='o', ls='', mfc=FASTED, mec=FASTED, ms=3.3, label='Fasting opposite to nodose HFD'),
                       Line2D([], [], marker='o', ls='', mfc='white', mec=GREY, ms=3.3, label='Fasting same direction'),
                       Line2D([], [], marker='o', ls='', mfc='#E3E3E3', mec='#E3E3E3', ms=3.3, label='All other genes')],
              loc='upper left', bbox_to_anchor=(-0.25, -0.2), fontsize=5.8, handletextpad=0.2, borderpad=0.2)
    fig.text(ax.get_position().x0 - 0.06, 0.93, 'a', weight='bold', fontsize=8)
    # b: heatmap of top candidates (by conjunction P) across all ganglia
    top = fu.head(16).index.tolist()
    for g in ['Lepr', 'Uts2b', 'Pnoc', 'Oprl1']:
        if g in fu.index and g not in top: top.append(g)
    cols = [('DRG_12wk_lfc', 'DRG_12wk_p', 'DRG\n12 wk'), ('DRG_36wk_lfc', 'DRG_36wk_p', 'DRG\n36 wk'),
            ('TRAP_lfc', 'TRAP_p', 'Nav1.8\nTRAP'), ('NOD_lfc', 'NOD_p', 'Nodose\nHFD'), ('P21_lfc', 'P21_p', 'Nodose\nearly P21'),
            ('SCG_obese_lfc', None, 'SCG\nobese*'), ('stellate_obese_lfc', None, 'Stellate\nobese*'), ('FAST_lfc', 'FAST_p', 'Nodose\nfasted')]
    M = np.array([[fu.loc[g, c] if c in fu and pd.notna(fu.loc[g, c]) else np.nan for c, _, _ in cols] for g in top], float)
    ax = fig.add_subplot(gs[1])
    im = ax.imshow(np.clip(M, -1.5, 1.5), cmap='RdBu_r', norm=TwoSlopeNorm(0, -1.5, 1.5), aspect='auto')
    for i, g in enumerate(top):
        for j, (c, pc, _) in enumerate(cols):
            if np.isnan(M[i, j]):
                ax.text(j, i, '×', ha='center', va='center', fontsize=5, color=GREY)
            elif pc and pd.notna(fu.loc[g, pc]) and fu.loc[g, pc] < 0.05:
                ax.text(j, i, '○', ha='center', va='center', fontsize=5.5)
            if c.startswith(('SCG', 'stellate')) and not np.isnan(M[i, j]):
                run = fu.loc[g, c.replace('obese_lfc', 'lean_run_diff')]
                if abs(M[i, j]) <= abs(run):
                    ax.text(j, i, '·', ha='center', va='center', fontsize=7, color=GREY)
    ax.axvline(4.5, color='white', lw=2); ax.axvline(6.5, color='white', lw=2)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels([f'{g}  ' + (('↑' if fu.loc[g, 'conj_dir'] == 'up' else '↓') if fu.loc[g, 'p_conj'] < 0.99 else '–') + f'  P={fu.loc[g, "p_conj"]:.2g}' for g in top], fontsize=5.8)
    for lab in ax.get_yticklabels(): lab.set_fontstyle('italic')
    ax.set_xticks(range(len(cols))); ax.set_xticklabels([c[2] for c in cols], fontsize=5.6); ax.xaxis.tick_top()
    ax.tick_params(length=0)
    for s in ax.spines.values(): s.set_visible(False)
    cb = fig.colorbar(im, ax=ax, orientation='horizontal', fraction=0.035, pad=0.02, aspect=30)
    cb.set_label('log$_2$ fold change (clipped at ±1.5)'); cb.outline.set_linewidth(0.4); cb.ax.tick_params(labelsize=6)
    cb.ax.text(0.5, -4.6, '○ nominal P < 0.05   × not expressed / absent   · within lean between-run difference\n'
               '* sympathetic: obese pool vs lean pools, diet confounded with run (descriptive)', transform=cb.ax.transAxes,
               ha='center', va='top', fontsize=5.5)
    fig.text(ax.get_position().x0 - 0.12, 0.93, 'b', weight='bold', fontsize=8)
    # c: adipocyte-marker check
    ax = fig.add_subplot(gs[2])
    adip = ['Pparg', 'Cidea', 'Cidec', 'Fabp4', 'Plin1', 'Cfd', 'Pck1']
    A = tw.reindex(adip)[['DRG_12wk_lfc', 'DRG_36wk_lfc', 'NOD_lfc']].values
    P = tw.reindex(adip)[['DRG_12wk_p', 'DRG_36wk_p', 'NOD_p']].values
    ax.imshow(np.clip(A, -1.5, 1.5), cmap='RdBu_r', norm=TwoSlopeNorm(0, -1.5, 1.5), aspect='auto')
    for i in range(len(adip)):
        for j in range(3):
            if P[i, j] < 0.05: ax.text(j, i, '○', ha='center', va='center', fontsize=5.5)
    ax.set_yticks(range(len(adip))); ax.set_yticklabels(adip, style='italic', fontsize=6)
    ax.set_xticks(range(3)); ax.set_xticklabels(['DRG\n12 wk', 'DRG\n36 wk', 'Nodose\nHFD'], fontsize=5.6); ax.xaxis.tick_top()
    ax.tick_params(length=0)
    for s in ax.spines.values(): s.set_visible(False)
    ax.set_xlabel('Adipocyte markers', fontsize=6.5)
    fig.text(ax.get_position().x0 - 0.06, 0.93, 'c', weight='bold', fontsize=8)
    save(fig, FIG/'Fig6_cross_ganglion_HFD_consensus'); plt.close(fig)
    fu.to_csv(SRC/'Fig6_consensus_source_data.csv')
    tw.reindex(adip).to_csv(SRC/'Fig6c_adipocyte_markers_source_data.csv')

if __name__ == '__main__':
    main(); print('done')
