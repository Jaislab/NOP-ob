"""Figure 7: expression of HFD/fasting counter-regulated candidates in the geniculate ganglion."""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).parent))
from figure_style import configure, save, MM
from load_data import ROOT

R = ROOT/'results'/'geniculate'; FIG = ROOT/'figures'; SRC = FIG/'source_data'
GENES = ['Acacb', 'Insr', 'Cnr1', 'Pomc', 'Cartpt']
DS = [('Geniculate (GSE223724, WT)', 'Geniculate', 'black', 'o'), ('Nodose (GSE281091, chow fed)', 'Nodose', '#7F7F7F', 's'),
      ('DRG (GSE225224, standard diet)', 'DRG', '#BDBDBD', 'D')]

def main():
    configure()
    per = pd.read_csv(R/'bulk_per_sample_tpm.csv'); pit = pd.read_csv(R/'geniculate_pituitary_check.csv', index_col=0)
    sc = pd.read_csv(R/'geniculate_single_neuron_detection.csv')
    fig = plt.figure(figsize=(183 * MM, 72 * MM))
    gs = fig.add_gridspec(1, 2, left=0.075, right=0.985, bottom=0.27, top=0.88, wspace=0.28, width_ratios=[1.25, 1])
    ax = fig.add_subplot(gs[0]); rng = np.random.default_rng(3)
    for i, g in enumerate(GENES):
        for k, (name, lab, col, mk) in enumerate(DS):
            d = per[(per.dataset == name) & (per.gene == g)]
            x = i + (k - 1) * 0.26 + (rng.random(len(d)) - 0.5) * 0.12
            y = np.maximum(d.tpm.values, 0.05)
            cont = d['sample'].map(pit.pituitary_contaminated).fillna(False).values.astype(bool) if (lab == 'Geniculate' and g == 'Pomc') else np.zeros(len(d), bool)
            ax.scatter(x[~cont], y[~cont], s=10, marker=mk, facecolors=col, edgecolors=col, lw=0.5, zorder=3)
            ax.scatter(x[cont], y[cont], s=10, marker=mk, facecolors='white', edgecolors=col, lw=0.6, zorder=3)
            ax.plot([i + (k - 1) * 0.26 - 0.1, i + (k - 1) * 0.26 + 0.1], [np.median(y[~cont]) if (~cont).any() else np.nan] * 2, color='#D55E00', lw=0.9, zorder=4)
    ax.set_yscale('log'); ax.set_ylim(0.04, 3000)
    ax.set_xticks(range(len(GENES))); ax.set_xticklabels(GENES, style='italic')
    for x0 in [0.5, 1.5, 2.5, 3.5]: ax.axvline(x0, color='#E5E5E5', lw=0.5, zorder=0)
    ax.axhspan(0.04, 0.5, color='#F2F2F2', zorder=0, lw=0)
    ax.set_ylabel('Bulk ganglion TPM (log scale)')
    ax.text(-0.45, 0.42, 'not detected', fontsize=5.3, color='#7F7F7F', va='top')
    import matplotlib.transforms as mt
    tr = mt.blended_transform_factory(ax.transData, ax.transAxes)
    ax.text(1, 1.02, 'HFD ↓ / fasting ↑', transform=tr, ha='center', fontsize=6, color='#7F7F7F')
    ax.text(3.5, 1.02, 'HFD ↑ / fasting ↓', transform=tr, ha='center', fontsize=6, color='#7F7F7F')
    ax.plot([-0.4, 2.4], [1.0, 1.0], transform=tr, color='#BDBDBD', lw=0.6, clip_on=False)
    ax.plot([2.6, 4.4], [1.0, 1.0], transform=tr, color='#BDBDBD', lw=0.6, clip_on=False)
    h = [Line2D([], [], marker=mk, ls='', mfc=col, mec=col, ms=3.5, label=lab) for _, lab, col, mk in DS]
    h += [Line2D([], [], marker='o', ls='', mfc='white', mec='black', ms=3.5, label='Geniculate, pituitary carry-over ($\\it{Pomc}$)'),
          Line2D([], [], color='#D55E00', lw=0.9, label='Median (excl. open)')]
    ax.legend(handles=h, loc='upper left', bbox_to_anchor=(0.0, -0.13), ncol=3, fontsize=6, columnspacing=1.0, handletextpad=0.3)
    fig.text(ax.get_position().x0 - 0.06, 0.92, 'a', weight='bold', fontsize=8)
    # b single-neuron detection
    ax = fig.add_subplot(gs[1])
    groups = [('GSE102443 (C1, FPKM)', 'Phox2b+ (oral/taste)', 'C1: Phox2b$^+$ oral/taste (n = 61)', 'black', 'o'),
              ('GSE102443 (C1, FPKM)', 'Phox2b- (somatosensory)', 'C1: Phox2b$^-$ somatosensory (n = 35)', '#9E9E9E', 'o'),
              ('GSE135801 (Phox2b+ neurons, reads)', 'all neurons (>= 100k reads)', 'Phox2b$^+$ neurons, GSE135801 (n = 376)', 'black', '^')]
    w = 0.25
    for k, (ds, grp, lab, col, mk) in enumerate(groups):
        d = sc[(sc.dataset == ds) & (sc.group == grp)].set_index('gene').reindex(GENES)
        x = np.arange(len(GENES)) + (k - 1) * w
        ax.bar(x, 100 * d.frac_detected, width=w * 0.9, color=col if k != 2 else 'white', edgecolor='black', lw=0.5, label=lab,
               hatch=None if k != 2 else '////')
    ax.set_xticks(range(len(GENES))); ax.set_xticklabels(GENES, style='italic')
    ax.set_ylim(0, 105); ax.set_ylabel('Geniculate neurons with expression (%)')
    ax.legend(loc='upper left', bbox_to_anchor=(-0.05, -0.13), fontsize=6, handlelength=1.2)
    fig.text(ax.get_position().x0 - 0.06, 0.92, 'b', weight='bold', fontsize=8)
    save(fig, FIG/'Fig7_geniculate_expression'); plt.close(fig)
    per.to_csv(SRC/'Fig7a_bulk_tpm_source_data.csv', index=False); sc.to_csv(SRC/'Fig7b_single_neuron_source_data.csv', index=False)

if __name__ == '__main__':
    main(); print('done')
