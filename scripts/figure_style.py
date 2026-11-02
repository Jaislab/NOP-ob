"""Shared print figure style; no empirical data or invented example results."""
from pathlib import Path
import matplotlib as mpl
from matplotlib import font_manager

CONTROL = '#0072B2'  # Okabe-Ito blue, RGB 0,114,178
HYPERCALORIC = '#CC79A7'  # Okabe-Ito pink, RGB 204,121,167
COLORS = {'Control': CONTROL, 'Hypercaloric': HYPERCALORIC}
MARKERS = {'Control': 'o', 'Hypercaloric': '^'}
MM = 1 / 25.4

def configure():
    fonts = Path('/usr/share/fonts/opentype/urw-base35')
    for path in fonts.glob('NimbusSans-*.otf'):
        font_manager.fontManager.addfont(str(path))
    resolved = font_manager.findfont('Nimbus Sans', fallback_to_default=False)
    mpl.rcParams.update({
        'font.family':'Nimbus Sans', 'font.size':7, 'axes.labelsize':7,
        'axes.titlesize':8, 'axes.titleweight':'normal', 'xtick.labelsize':6.5,
        'ytick.labelsize':6.5, 'legend.fontsize':6.5, 'axes.linewidth':0.6,
        'lines.linewidth':0.8, 'xtick.major.width':0.6, 'ytick.major.width':0.6,
        'xtick.major.size':2.5, 'ytick.major.size':2.5,
        'axes.spines.top':False, 'axes.spines.right':False,
        'axes.grid':False, 'legend.frameon':False, 'pdf.fonttype':42,
        'ps.fonttype':42, 'svg.fonttype':'none', 'figure.dpi':150,
        'savefig.dpi':600, 'savefig.facecolor':'white',
    })
    return resolved

def save(fig, stem):
    """Keep the specified physical canvas; do not resize with bbox_inches=tight."""
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    for ext in ['pdf','svg','png']:
        fig.savefig(stem.with_suffix('.' + ext), dpi=600)

def specimen():
    """Non-data typography and palette proof, explicitly labelled as such."""
    import matplotlib.pyplot as plt
    font = configure()
    fig, ax = plt.subplots(figsize=(89*MM, 49*MM))
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(0, .94, 'Typography and palette proof', fontsize=8, weight='bold', va='top')
    ax.text(0, .72, 'Nimbus Sans · 7 pt at final print size', va='top')
    for x, label in [(.03, 'Control'), (.53, 'Hypercaloric')]:
        ax.scatter([x], [.46], c=COLORS[label], marker=MARKERS[label], s=28)
        ax.text(x+.055, .46, label, va='center')
        ax.text(x+.055, .29, COLORS[label], va='center', fontsize=6.5)
    ax.text(0, .04, 'Style proof only — no expression data or statistical results.', fontsize=6)
    fig.subplots_adjust(left=.07,right=.98,bottom=.08,top=.98)
    save(fig, Path(__file__).resolve().parents[1] / 'figures' / 'STYLE_PROOF_NOT_DATA')
    print('Verified Nimbus Sans:', font)

if __name__ == '__main__':
    specimen()
