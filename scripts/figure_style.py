"""Shared print figure style (Nimbus Sans, Okabe-Ito blue/pink) and export helpers.

Matplotlib embeds CFF-flavoured OpenType fonts as Type 42 incorrectly (Poppler reports
'Mismatch between font type and embedded font file'). We therefore convert the URW
Nimbus Sans OTF files to TrueType (quadratic) outlines once, cache them under .cache/fonts,
and register only those, so PDFs embed genuine TrueType (Type 42) fonts.
"""
from pathlib import Path
import matplotlib as mpl
from matplotlib import font_manager

ROOT = Path(__file__).resolve().parents[1]
CONTROL = '#0072B2'       # Okabe-Ito blue, RGB 0,114,178
HYPERCALORIC = '#CC79A7'  # Okabe-Ito pink, RGB 204,121,167
COLORS = {'Control': CONTROL, 'Hypercaloric': HYPERCALORIC}
MARKERS = {'Control': 'o', 'Hypercaloric': '^'}
MM = 1 / 25.4
URW = Path('/usr/share/fonts/opentype/urw-base35')
CACHE = ROOT/'.cache'/'fonts'

def _otf_to_ttf(src, dst, max_err=1.0):
    from fontTools.ttLib import TTFont, newTable
    from fontTools.pens.cu2quPen import Cu2QuPen
    from fontTools.pens.ttGlyphPen import TTGlyphPen
    f = TTFont(src)
    order = f.getGlyphOrder()
    gs = f.getGlyphSet()
    f['loca'] = newTable('loca')
    f['glyf'] = glyf = newTable('glyf')
    glyf.glyphOrder = order
    glyf.glyphs = {}
    for name in order:
        pen = TTGlyphPen(gs)
        gs[name].draw(Cu2QuPen(pen, max_err, reverse_direction=True))
        glyf[name] = pen.glyph()
    del f['CFF ']
    if 'VORG' in f:
        del f['VORG']
    glyf.compile(f)
    hmtx = f['hmtx']
    for name, g in glyf.glyphs.items():
        if hasattr(g, 'xMin'):
            hmtx[name] = (hmtx[name][0], g.xMin)
    f['maxp'] = maxp = newTable('maxp')
    maxp.tableVersion = 0x00010000
    for k in ['maxZones', 'maxTwilightPoints', 'maxStorage', 'maxFunctionDefs', 'maxInstructionDefs',
              'maxStackElements', 'maxSizeOfInstructions', 'maxComponentElements']:
        setattr(maxp, k, 1 if k == 'maxZones' else 0)
    maxp.compile(f)
    post = f['post']
    post.formatType = 2.0; post.extraNames = []; post.mapping = {}; post.glyphOrder = order
    f['head'].glyphDataFormat = 0
    f.sfntVersion = '\x00\x01\x00\x00'
    f.save(dst)

def configure():
    CACHE.mkdir(parents=True, exist_ok=True)
    srcs = sorted(URW.glob('NimbusSans-*.otf'))
    if not srcs:
        raise RuntimeError('Nimbus Sans not installed (Debian/Ubuntu package fonts-urw-base35)')
    for src in srcs:
        dst = CACHE/(src.stem + '.ttf')
        if not dst.exists():
            _otf_to_ttf(src, dst)
        font_manager.fontManager.addfont(str(dst))
    fm = font_manager.fontManager
    fm.ttflist = [e for e in fm.ttflist if not (e.name == 'Nimbus Sans' and not e.fname.startswith(str(CACHE)))]
    for attr in ('_findfont_cached',):
        try:
            getattr(fm, attr).cache_clear()
        except AttributeError:
            pass
    font_manager._get_font.cache_clear()
    resolved = font_manager.findfont(font_manager.FontProperties(family='Nimbus Sans'), fallback_to_default=False)
    if not resolved.endswith('.ttf'):
        raise RuntimeError(f'Resolved non-TrueType font {resolved}')
    mpl.rcParams.update({
        'font.family': 'Nimbus Sans', 'font.size': 7, 'axes.labelsize': 7,
        'axes.titlesize': 7, 'axes.titleweight': 'normal', 'xtick.labelsize': 6.5,
        'ytick.labelsize': 6.5, 'legend.fontsize': 6.5, 'axes.linewidth': 0.6,
        'lines.linewidth': 0.8, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
        'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'xtick.minor.width': 0.4, 'ytick.minor.width': 0.4,
        'axes.spines.top': False, 'axes.spines.right': False,
        'axes.grid': False, 'legend.frameon': False, 'pdf.fonttype': 42,
        'ps.fonttype': 42, 'svg.fonttype': 'none', 'figure.dpi': 150,
        'savefig.dpi': 600, 'savefig.facecolor': 'white', 'mathtext.fontset': 'custom',
        'mathtext.rm': 'Nimbus Sans', 'mathtext.it': 'Nimbus Sans:italic', 'mathtext.bf': 'Nimbus Sans:bold', 'mathtext.cal': 'Nimbus Sans', 'mathtext.sf': 'Nimbus Sans', 'mathtext.tt': 'Nimbus Sans',
        'pdf.compression': 6,
    })
    return resolved

def panel_label(ax, letter, x=-0.02, y=1.0, fig=None):
    """Lowercase bold panel letter placed in figure coordinates at the axes' top-left."""
    fig = fig or ax.figure
    bb = ax.get_position()
    fig.text(bb.x0 + x, bb.y1 + 0.01 if y == 1.0 else y, letter, fontsize=8, weight='bold', ha='right', va='bottom')

def save(fig, stem):
    """Keep the specified physical canvas; never resize with bbox_inches='tight'."""
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    meta = {'Creator': 'NOP-ob reanalysis (matplotlib)'}
    fig.savefig(stem.with_suffix('.pdf'), metadata=meta)
    fig.savefig(stem.with_suffix('.svg'))
    fig.savefig(stem.with_suffix('.png'), dpi=600)

if __name__ == '__main__':
    print('Resolved', configure())
