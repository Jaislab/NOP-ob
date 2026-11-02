"""Publication figures and per-figure source data from results/ tables.

All values plotted are read from results/*.csv produced by run_bulk_de.py, sensitivity.py,
sc_descriptive.py and context_gse309905.py. Canvas sizes are physical (89/183 mm).
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

sys.path.insert(0, str(Path(__file__).parent))
from figure_style import configure, save, COLORS, MARKERS, MM, CONTROL, HYPERCALORIC
from load_data import ROOT, PANELS

R = ROOT/'results'
FIG = ROOT/'figures'
SRC = FIG/'source_data'
SRC.mkdir(parents=True, exist_ok=True)
GREY = '#6E6E6E'
LIGHT = '#BFBFBF'

ROWS = [  # (contrast_id, ganglion group, label, accession)
    ('DRG_HFD_12wk', 'DRG', 'Adult HFD, 9 wk (age 12 wk)', 'GSE225224'),
    ('DRG_HFD_36wk', 'DRG', 'Adult HFD, 33 wk (age 36 wk)', 'GSE225224'),
    ('DRG_Nav18_TRAP_WD', 'DRG', 'Nav1.8$^+$ TRAP, Western diet 11 wk', 'GSE198895'),
    ('Nodose_HFD_4wk', 'Nodose', 'Adult HFD, 4 wk', 'GSE281091'),
    ('Nodose_earlyHFD_P21', 'Nodose', 'Early-life HFD, at weaning (P21)', 'GSE320302'),
    ('Nodose_earlyHFD_P35', 'Nodose', 'Early-life HFD, 2 wk after weaning (P35)', 'GSE320302'),
    ('Symp_obese_scRNA', 'Sympathetic', 'SCG + stellate scRNA-seq, HFD 12 wk', 'GSE233163'),
]
SHORT = {'DRG_HFD_12wk': 'DRG, HFD\n12 wk of age', 'DRG_HFD_36wk': 'DRG, HFD\n36 wk of age',
         'Nodose_HFD_4wk': 'Nodose, adult\nHFD 4 wk', 'Nodose_earlyHFD_P21': 'Nodose, early-life\nHFD, P21',
         'Nodose_earlyHFD_P35': 'Nodose, early-life\nHFD, P35'}

def fmt_q(q):
    if pd.isna(q):
        return ''
    return f'{q:.2g}' if q >= 0.001 else f'{q:.1e}'

def row_positions():
    y, pos, seps, heads = 0.0, {}, [], []
    last = None
    for cid, grp, _, _ in ROWS:
        if grp != last:
            if last is not None:
                seps.append(y + 0.15); y += 0.35
            heads.append((grp, y)); y += 0.85
            last = grp
        pos[cid] = y; y += 1.0
    return pos, seps, heads, y

def nontestable_text(cid, gene, stats):
    if cid == 'Symp_obese_scRNA':
        if gene == 'Pnoc' or gene in ('Oprd1', 'Oprk1', 'Pdyn'):
            return 'absent from deposited matrix'
        return 'descriptive only (diet confounded with run)'
    r = stats[(stats.contrast_id == cid) & (stats.gene == gene)]
    if len(r) and isinstance(r.result_class.iat[0], str) and r.result_class.iat[0].startswith('not detected'):
        return 'not detected (0 counts, 4/4 libraries)'
    return None

COMPACT = {'not detected (0 counts, 4/4 libraries)': 'not detected', 'absent from deposited matrix': 'not in matrix',
           'descriptive only (diet confounded with run)': 'descriptive only'}

def forest(fig, gs_label, gs_cols, genes, stats, xlims, panel_letters, title_note=None, compact=False):
    pos, seps, heads, ymax = row_positions()
    axl = fig.add_subplot(gs_label)
    axl.set_ylim(ymax - 0.2, -0.4); axl.set_xlim(0, 1); axl.axis('off')
    for grp, y in heads:
        axl.text(0, y, grp, weight='bold', va='center', ha='left')
    for cid, grp, label, acc in ROWS:
        axl.text(0.04, pos[cid], label, va='center', ha='left')
        axl.text(1.0, pos[cid], acc, va='center', ha='right', color=GREY, fontsize=6)
    out = []
    for k, (g, (gsf, gst)) in enumerate(zip(genes, gs_cols)):
        ax = fig.add_subplot(gsf); axt = fig.add_subplot(gst, sharey=ax)
        ax.set_ylim(ymax - 0.2, -0.4)
        ax.axvline(0, color=GREY, lw=0.6, ls=(0, (2, 2)), zorder=0)
        for s in seps:
            ax.axhline(s, color=LIGHT, lw=0.4); axt.axhline(s, color=LIGHT, lw=0.4)
        lo, hi = xlims[k]
        for cid, grp, label, acc in ROWS:
            y = pos[cid]
            r = stats[(stats.contrast_id == cid) & (stats.gene == g)]
            note = nontestable_text(cid, g, stats)
            if note or r.empty or pd.isna(r.log2fc.iat[0]):
                ax.text((lo + hi) / 2, y, (COMPACT.get(note, note) if compact else note) or 'not testable', ha='center', va='center', color=GREY, fontsize=6, style='italic',
                        bbox=dict(fc='white', ec='none', pad=0.6), zorder=5)
                out.append(dict(gene=g, contrast_id=cid, series=acc, status=note or 'not testable'))
                continue
            r = r.iloc[0]
            low = r.result_class == 'low-count, unstable'
            explo = 'exploratory' in str(r.contrast_status)
            cl, ch = max(r.ci_low, lo), min(r.ci_high, hi)
            ax.plot([cl, ch], [y, y], color='black', lw=0.8, solid_capstyle='butt')
            for edge, clipped in [(r.ci_low, r.ci_low < lo), (r.ci_high, r.ci_high > hi)]:
                if clipped:
                    ax.annotate('', xy=(lo if edge < lo else hi, y), xytext=((lo if edge < lo else hi) + (0.06 if edge < lo else -0.06) * (hi - lo), y),
                                arrowprops=dict(arrowstyle='-|>', lw=0.6, color='black', mutation_scale=4))
            if not lo <= r.log2fc <= hi:
                ax.text(hi if r.log2fc > hi else lo, y + 0.42, f'{r.log2fc:.1f}', ha='right' if r.log2fc > hi else 'left', va='center', fontsize=5.5)
            ax.plot(np.clip(r.log2fc, lo, hi), y, marker='s', clip_on=False, ms=3.6, mfc='white' if low else 'black', mec='black', mew=0.7, zorder=3)
            axt.text(0.02, y, f'{int(r.n_control)}/{int(r.n_hypercaloric)}', va='center', ha='left', fontsize=6)
            axt.text(0.98, y, fmt_q(r.padj_family_BH), va='center', ha='right', fontsize=6,
                     weight='bold' if r.padj_family_BH < 0.05 and g in PANELS['primary'] else 'normal')
            out.append(dict(gene=g, contrast_id=cid, series=acc, status=r.result_class, n_control=r.n_control, n_hypercaloric=r.n_hypercaloric,
                            log2fc=r.log2fc, se=r.se, ci_low=r.ci_low, ci_high=r.ci_high, p=r.p, padj_family_BH=r.padj_family_BH,
                            padj_family_Holm=r.padj_family_Holm, padj_transcriptome=r.padj_transcriptome, method=r.method))
        ax.set_xlim(lo, hi)
        ax.spines['left'].set_visible(False); ax.tick_params(axis='y', left=False, labelleft=False)
        ax.set_xlabel('log$_2$ fold change\n(hypercaloric / control)')
        ax.xaxis.set_major_locator(MaxNLocator(5, integer=(hi - lo) > 3))
        if 2 <= hi - lo <= 3:
            ax.set_xticks([-1, 0, 1])
        ax.set_title(g, style='italic', pad=4)
        axt.axis('off'); axt.set_xlim(0, 1)
        axt.text(0.02, -0.4, 'n C/H', ha='left', va='bottom', fontsize=6, color=GREY)
        axt.text(0.98, -0.4, 'q', ha='right', va='bottom', fontsize=6, color=GREY, style='italic')
        fig.text(ax.get_position().x0 - 0.008 if k else axl.get_position().x0, ax.get_position().y1 + 0.035, panel_letters[k],
                 weight='bold', fontsize=8, ha='left' if k == 0 else 'right', va='bottom')
    return pd.DataFrame(out)

# ------------------------------------------------------------------ Figure 1
def figure1(stats):
    fig = plt.figure(figsize=(183 * MM, 80 * MM))
    gs = fig.add_gridspec(1, 5, width_ratios=[2.35, 1.55, 0.42, 1.55, 0.42], left=0.015, right=0.99, bottom=0.26, top=0.86, wspace=0.08)
    src = forest(fig, gs[0], [(gs[1], gs[2]), (gs[3], gs[4])], ['Pnoc', 'Oprl1'], stats, [(-6, 3), (-0.75, 0.75)], ['a', 'b'])
    handles = [Line2D([], [], marker='s', ls='-', color='black', mfc='black', ms=3.6, lw=0.8, label='Estimate, nominal 95% CI'),
               Line2D([], [], marker='s', ls='', color='black', mfc='white', ms=3.6, label='Low abundance (mean normalized count <10)')]
    fig.legend(handles=handles, loc='lower left', bbox_to_anchor=(0.015, 0.0), ncol=2, handlelength=1.6, columnspacing=1.5)
    save(fig, FIG/'Fig1_Pnoc_Oprl1_diet_effects'); plt.close(fig)
    src.to_csv(SRC/'Fig1_forest_source_data.csv', index=False)

# ------------------------------------------------------------------ Figure 2
def strip(ax, d, ycol, sexfill=False, rng=None):
    rng = rng or np.random.default_rng(1)
    for i, diet in enumerate(['Control', 'Hypercaloric']):
        x = d[d.diet == diet]
        jit = (rng.random(len(x)) - 0.5) * 0.36
        filled = np.ones(len(x), bool) if not sexfill else (x.sex == 'male').values
        for f in [True, False]:
            k = filled == f
            if k.any():
                ax.scatter(i + jit[k], x[ycol].values[k], s=11, marker=MARKERS[diet], facecolors=COLORS[diet] if f else 'white',
                           edgecolors=COLORS[diet], linewidths=0.7, zorder=3, clip_on=False)
        med = np.median(x[ycol])
        ax.plot([i - 0.28, i + 0.28], [med, med], color='black', lw=0.9, zorder=4, solid_capstyle='butt')
    ax.set_xlim(-0.6, 1.6); ax.set_xticks([0, 1]); ax.set_xticklabels(['Ctrl', 'HFD'])
    ax.set_ylim(0, d[ycol].max() * 1.12 if d[ycol].max() > 0 else 1)

def figure2(vals, stats):
    cids = ['DRG_HFD_12wk', 'DRG_HFD_36wk', 'Nodose_HFD_4wk', 'Nodose_earlyHFD_P21', 'Nodose_earlyHFD_P35']
    inv = pd.read_csv(ROOT/'metadata'/'contrast_inventory.csv').set_index('contrast_id')
    fig, axes = plt.subplots(2, 5, figsize=(183 * MM, 92 * MM))
    fig.subplots_adjust(left=0.075, right=0.995, bottom=0.135, top=0.88, wspace=0.55, hspace=0.42)
    src = []
    rng = np.random.default_rng(7)
    for i, g in enumerate(['Pnoc', 'Oprl1']):
        for j, cid in enumerate(cids):
            c = inv.loc[cid]; ax = axes[i, j]
            d = vals[(vals.series == c.series) & (vals.gene == g) & vals.group.isin([c.numerator, c.denominator])].copy()
            strip(ax, d, 'normalized', sexfill=(c.series == 'GSE281091'), rng=rng)
            ax.yaxis.set_major_locator(MaxNLocator(4))
            r = stats[(stats.contrast_id == cid) & (stats.gene == g)].iloc[0]
            ax.text(0.5, 1.03, f'q = {fmt_q(r.padj_family_BH)}', transform=ax.transAxes, ha='center', va='bottom', fontsize=6, color=GREY)
            if i == 0:
                ax.set_title(SHORT[cid], fontsize=6.5, pad=13)
            if j == 0:
                ax.set_ylabel(f'$\\it{{{g}}}$ normalized counts')
            d['contrast_id'] = cid; src.append(d)
    fig.text(0.008, 0.9, 'a', weight='bold', fontsize=8); fig.text(0.008, 0.47, 'b', weight='bold', fontsize=8)
    handles = [Line2D([], [], marker='o', ls='', mfc=CONTROL, mec=CONTROL, ms=3.5, label='Control diet'),
               Line2D([], [], marker='^', ls='', mfc=HYPERCALORIC, mec=HYPERCALORIC, ms=3.8, label='Hypercaloric diet'),
               Line2D([], [], color='black', lw=0.9, label='Median'),
               Line2D([], [], marker='^', ls='', mfc='white', mec=GREY, ms=3.8, label='Open: female (GSE281091 only)')]
    fig.legend(handles=handles, loc='lower center', ncol=4, bbox_to_anchor=(0.5, 0.0), columnspacing=1.6, handletextpad=0.4)
    save(fig, FIG/'Fig2_replicate_expression'); plt.close(fig)
    pd.concat(src).to_csv(SRC/'Fig2_replicate_values_source_data.csv', index=False)

# ------------------------------------------------------------------ Figure 3
def figure3():
    t = pd.read_csv(R/'single_cell'/'GSE233163_pool_celltype_summary.csv')
    ctx = pd.read_csv(R/'single_cell'/'GSE309905_context_detection.csv')
    fig = plt.figure(figsize=(183 * MM, 70 * MM))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.55, 0.9, 1.0], left=0.155, right=0.99, bottom=0.25, top=0.88, wspace=0.5)
    # a: Oprl1 detection by cell type, per pool meeting coverage
    ax = fig.add_subplot(gs[0])
    order = (t[t.meets_coverage].groupby('cell_type').Oprl1_frac_detected.mean().sort_values().index.tolist())
    nice = {c: c.replace('_', ' ').replace('schwann cell myelinating', 'myelinating Schwann').replace('schwann cell', 'Schwann')
            .replace('t reg lymphocyte', 'T reg lymphocyte') for c in order}
    tt = t[t.meets_coverage]
    for k, ct in enumerate(order):
        for diet, off in [('Control', -0.15), ('Hypercaloric', 0.15)]:
            x = tt[(tt.cell_type == ct) & (tt.diet == diet)]
            ax.scatter(100 * x.Oprl1_frac_detected, np.full(len(x), k + off), marker=MARKERS[diet], s=10,
                       facecolors=COLORS[diet], edgecolors=COLORS[diet], linewidths=0.5, zorder=3, clip_on=False)
    ax.set_yticks(range(len(order))); ax.set_yticklabels([nice[c] for c in order])
    ax.set_xlabel('Cells with $\\it{Oprl1}$ UMI >0 (%)'); ax.set_xlim(-1.5, 50)
    ax.set_title('Sympathetic ganglia (GSE233163)', fontsize=6.5)
    # b: neuron pseudobulk CPM per pool
    ax2 = fig.add_subplot(gs[1])
    nn = t[(t.cell_type == 'neuron') & t.meets_coverage]
    for k, sitev in enumerate(['SCG', 'stellate']):
        for diet, off in [('Control', -0.17), ('Hypercaloric', 0.17)]:
            x = nn[(nn.site == sitev) & (nn.diet == diet)]
            ax2.scatter(np.full(len(x), k + off), x.Oprl1_cpm, marker=MARKERS[diet], s=13, facecolors=COLORS[diet],
                        edgecolors=COLORS[diet], linewidths=0.5, zorder=3, clip_on=False)
            for _, r in x.iterrows():
                ax2.text(k + off + 0.07, r.Oprl1_cpm, r.run, fontsize=5, va='center', color=GREY)
    ax2.set_xticks([0, 1]); ax2.set_xticklabels(['SCG', 'Stellate']); ax2.set_xlim(-0.55, 1.6)
    ax2.set_ylabel('Neuron pseudobulk $\\it{Oprl1}$ (CPM)'); ax2.set_ylim(-1, 25)
    ax2.set_title('Neurons, per pool', fontsize=6.5)
    # c: GSE309905 detection
    ax3 = fig.add_subplot(gs[2])
    for k, (grp, lab) in enumerate([('putative neuronal nuclei', 'Neuronal'), ('other nuclei', 'Other')]):
        for gi, (gene, off) in enumerate([('Pnoc', -0.17), ('Oprl1', 0.17)]):
            x = ctx[(ctx.nuclei_group == grp) & (ctx.n_nuclei >= 10)]
            ax3.scatter(np.full(len(x), k + off), 100 * x[f'{gene}_frac_detected'], marker='^', s=12,
                        facecolors=HYPERCALORIC if gene == 'Oprl1' else 'white', edgecolors=HYPERCALORIC, linewidths=0.7, zorder=3, clip_on=False)
    ax3.set_xticks([0, 1]); ax3.set_xticklabels(['Neuronal\nnuclei', 'Other\nnuclei']); ax3.set_xlim(-0.55, 1.55)
    ax3.set_ylabel('Nuclei with UMI >0 (%)'); ax3.set_ylim(-0.8, 16)
    ax3.set_title('HFD-fed vagus (GSE309905)', fontsize=6.5)
    for a, l in [(ax, 'a'), (ax2, 'b'), (ax3, 'c')]:
        bb = a.get_position(); fig.text(bb.x0 - (0.14 if l == 'a' else 0.06), bb.y1 + 0.04, l, weight='bold', fontsize=8)
    handles = [Line2D([], [], marker='o', ls='', mfc=CONTROL, mec=CONTROL, ms=3.5, label='Control diet (lean)'),
               Line2D([], [], marker='^', ls='', mfc=HYPERCALORIC, mec=HYPERCALORIC, ms=3.8, label='Hypercaloric diet'),
               Line2D([], [], marker='^', ls='', mfc='white', mec=HYPERCALORIC, ms=3.8, label='c: open, $\\it{Pnoc}$; filled, $\\it{Oprl1}$')]
    fig.legend(handles=handles, loc='lower center', ncol=3, bbox_to_anchor=(0.55, 0.0), columnspacing=1.8, handletextpad=0.4)
    save(fig, FIG/'Fig3_single_cell_context'); plt.close(fig)
    t.to_csv(SRC/'Fig3ab_GSE233163_source_data.csv', index=False)
    ctx.to_csv(SRC/'Fig3c_GSE309905_source_data.csv', index=False)

# ------------------------------------------------------------------ Supplementary QC
def figure_s1():
    qc = pd.read_csv(R/'qc'/'all_bulk_sample_qc.csv')
    series = [('GSE225224', 'DRG (GSE225224)'), ('GSE198895', 'DRG TRAP (GSE198895)'), ('GSE281091', 'Nodose adult (GSE281091)'), ('GSE320302', 'Nodose early life (GSE320302)')]
    fig, axes = plt.subplots(2, 4, figsize=(183 * MM, 88 * MM))
    fig.subplots_adjust(left=0.07, right=0.99, bottom=0.17, top=0.9, wspace=0.55, hspace=0.6)
    for j, (s, title) in enumerate(series):
        q = qc[qc.series == s]
        ax = axes[0, j]
        for diet in ['Control', 'Hypercaloric']:
            x = q[q.diet == diet]
            age = x.group.str.split('_').str[0] if x.group.str.contains('_').all() else None
            filled = np.ones(len(x), bool) if age is None or age.nunique() < 1 else (age == sorted(q.group.str.split('_').str[0].unique())[0]).values
            for f in [True, False]:
                k = filled == f
                ax.scatter(x.PC1[k], x.PC2[k], marker=MARKERS[diet], s=11, facecolors=COLORS[diet] if f else 'white', edgecolors=COLORS[diet], linewidths=0.6)
            fl = x[x.flag_low_similarity]
            ax.scatter(fl.PC1, fl.PC2, s=34, facecolors='none', edgecolors='black', linewidths=0.4)
        ax.set_xlabel(f'PC1 ({100 * q.pc1_var.iat[0]:.0f}%)'); ax.set_ylabel(f'PC2 ({100 * q.pc2_var.iat[0]:.0f}%)')
        ax.set_title(title, fontsize=6.5); ax.xaxis.set_major_locator(MaxNLocator(4)); ax.yaxis.set_major_locator(MaxNLocator(4))
        ax2 = axes[1, j]
        q2 = q.assign(libM=q.library_size / 1e6)
        strip(ax2, q2, 'libM', rng=np.random.default_rng(3))
        ax2.set_ylabel('Library size (M)' if s != 'GSE281091' else 'Sum of normalized counts (M)')
    fig.text(0.005, 0.92, 'a', weight='bold', fontsize=8); fig.text(0.005, 0.46, 'b', weight='bold', fontsize=8)
    handles = [Line2D([], [], marker='o', ls='', mfc=CONTROL, mec=CONTROL, ms=3.5, label='Control'),
               Line2D([], [], marker='^', ls='', mfc=HYPERCALORIC, mec=HYPERCALORIC, ms=3.8, label='Hypercaloric'),
               Line2D([], [], marker='o', ls='', mfc='white', mec=GREY, ms=3.5, label='Open: later age (36 wk / P35)'),
               Line2D([], [], marker='o', ls='', mfc='none', mec='black', ms=6, mew=0.4, label='Ring: low similarity (robust z <−3.5)')]
    fig.legend(handles=handles, loc='lower center', ncol=4, bbox_to_anchor=(0.5, 0.0), columnspacing=1.4, handletextpad=0.3)
    save(fig, FIG/'FigS1_sample_QC'); plt.close(fig)
    qc.to_csv(SRC/'FigS1_QC_source_data.csv', index=False)

def figure_s2(stats):
    fig = plt.figure(figsize=(183 * MM, 80 * MM))
    gs = fig.add_gridspec(1, 7, width_ratios=[3.0, 1.1, 0.42, 1.1, 0.42, 1.1, 0.42], left=0.015, right=0.99, bottom=0.26, top=0.86, wspace=0.08)
    src = forest(fig, gs[0], [(gs[1], gs[2]), (gs[3], gs[4]), (gs[5], gs[6])], PANELS['secondary_processing'], stats,
                 [(-1.5, 1.5)] * 3, ['a', 'b', 'c'], compact=True)
    handles = [Line2D([], [], marker='s', ls='-', color='black', ms=3.6, lw=0.8, label='Estimate, nominal 95% CI (arrow: CI beyond axis)')]
    fig.legend(handles=handles, loc='lower left', bbox_to_anchor=(0.015, 0.0))
    save(fig, FIG/'FigS2_processing_enzymes'); plt.close(fig)
    src.to_csv(SRC/'FigS2_source_data.csv', index=False)

def figure_s3(stats):
    genes = PANELS['exploratory_opioid']
    fig = plt.figure(figsize=(183 * MM, 140 * MM))
    outer = fig.add_gridspec(2, 1, left=0.015, right=0.99, bottom=0.16, top=0.93, hspace=0.75)
    srcs = []
    for i in range(2):
        gs = outer[i].subgridspec(1, 7, width_ratios=[3.0, 1.1, 0.42, 1.1, 0.42, 1.1, 0.42], wspace=0.08)
        g3 = genes[3 * i: 3 * i + 3]
        srcs.append(forest(fig, gs[0], [(gs[1], gs[2]), (gs[3], gs[4]), (gs[5], gs[6])], g3, stats, [(-4, 4)] * 3, list('abcdef'[3 * i: 3 * i + 3]), compact=True))
    handles = [Line2D([], [], marker='s', ls='-', color='black', ms=3.6, lw=0.8, label='Estimate, nominal 95% CI'),
               Line2D([], [], marker='s', ls='', color='black', mfc='white', ms=3.6, label='Low abundance / sparse'),
               Line2D([], [], ls='', label='Number at axis edge: estimate beyond axis')]
    fig.legend(handles=handles, loc='lower left', bbox_to_anchor=(0.015, 0.0), ncol=3)
    save(fig, FIG/'FigS3_opioid_exploratory'); plt.close(fig)
    pd.concat(srcs).to_csv(SRC/'FigS3_source_data.csv', index=False)

if __name__ == '__main__':
    print('font', configure())
    stats = pd.read_csv(R/'target_gene_statistics.csv')
    vals = pd.read_csv(R/'target_gene_sample_values.csv')
    figure1(stats); figure2(vals, stats); figure3(); figure_s1(); figure_s2(stats); figure_s3(stats)
    print('done')
