# -*- coding: utf-8 -*-
"""06_model_diagrams.py -- the analysis-step diagram, Fig 2 of the manuscript.

    python 06_model_diagrams.py

  Fig 2  the eight analysis steps

Fig 1, the path diagram of the analysis model, is drawn in PowerPoint
(Fig1.pptx) and exported from there; it is not produced by this script.

The diagram is a schematic and carries no numbers. Until v44 the arrow leaving
step 2 carried the p value of the equality test among the class coefficients;
from v45 the manuscript reports no such test, and the arrow states the
comparison in words.

Output follows the PLOS ONE figure requirements: a TIFF with LZW compression,
RGB at 8 bits per channel and no alpha channel, 600 dpi, 7.5 inches wide, text
in Times New Roman at 8 pt or more at the printed size, and no figure number,
title or legend inside the image. The figure is drawn at its printed size, so
the type sizes below are the printed sizes. Black on white, legible without
colour.
"""
import io
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.transforms import Bbox
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
OUTD = os.path.join(HERE, "outputs", "figures_plos")
os.makedirs(OUTD, exist_ok=True)

mpl.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman"],
    "mathtext.fontset": "custom",
    "mathtext.rm": "Times New Roman",
    "mathtext.it": "Times New Roman:italic",
    "mathtext.bf": "Times New Roman:bold",
    "font.size": 8,
    "pdf.fonttype": 42,         # TrueType in the PDF used for the font check
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
})

GREY = "#4d4d4d"
DPI = 600
WIDTH_IN = 7.5                  # PLOS ONE: 2.63 to 7.5 inches wide
MAX_HEIGHT_IN = 8.75
MIN_PT = 8.0
FONT_FILES = ("TimesNewRoman",)  # embedded PostScript names must start with this


def check_figure(fig, stem):
    """Measure what was drawn: nothing outside the canvas, every text at 8 pt
    or more, and every glyph set in the named typeface (read back from the
    fonts a PDF rendering embeds, so a fallback to another font is caught)."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    w, h = fig.get_size_inches()
    # extent of what is drawn: texts, boxes and arrows (the axes are switched off)
    ext = [a.get_window_extent(r) for ax in fig.axes for a in ax.texts + ax.patches]
    bb = Bbox.union(ext).transformed(fig.dpi_scale_trans.inverted())
    assert abs(w - WIDTH_IN) < 1e-9 and h <= MAX_HEIGHT_IN, (w, h)
    assert bb.x0 >= 0 and bb.y0 >= 0 and bb.x1 <= w + 1e-6 and bb.y1 <= h + 1e-6, \
        "something is drawn outside the canvas: %s" % bb
    texts = [t for t in fig.findobj(mpl.text.Text)
             if t.get_visible() and t.get_text().strip()]
    smallest = min(t.get_fontsize() for t in texts)
    assert smallest >= MIN_PT, smallest
    buf = io.BytesIO()
    fig.savefig(buf, format="pdf")
    fonts = sorted(set(re.findall(rb"/BaseFont /([A-Za-z0-9+,\-]+)", buf.getvalue())))
    fonts = [f.decode() for f in fonts]
    assert fonts and all(f.split("+")[-1].startswith(FONT_FILES) for f in fonts), fonts
    print("  %s: %.2f x %.2f in, content %.2f to %.2f by %.2f to %.2f in, "
          "%d texts, smallest %.2f pt, fonts %s"
          % (stem, w, h, bb.x0, bb.x1, bb.y0, bb.y1, len(texts), smallest, fonts))


def save(fig, stem):
    """TIFF, LZW, RGB without alpha, 600 dpi, at the drawn size."""
    check_figure(fig, stem)
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=DPI)
    plt.close(fig)
    im = Image.open(buf).convert("RGB")
    path = os.path.join(OUTD, "%s.tif" % stem)
    im.save(path, format="TIFF", compression="tiff_lzw", dpi=(DPI, DPI))
    assert im.size[0] == round(WIDTH_IN * DPI), im.size
    print("  wrote %s (%d x %d px)" % (path, im.size[0], im.size[1]))


# ----------------------------------------------- Fig 2, the eight steps
# Type sizes are the printed sizes (the figure is placed at the size drawn).
STEP_PT, PURPOSE_PT, HOW_PT, WHERE_PT, ARROW_PT = 9.84, 8.64, 9.12, 8.28, 8.28
# vertical scale of the drawing, in inches per unit of the 0 to 100 y range
Y_IN = 0.0744
Y_RANGE = (6.5, 97.5)


def flow_diagram():
    """The eight steps run in a fixed order, and the exposure changed from
    the classes to the count at step 2 for a stated reason."""
    steps = [
        ("1", "Combinations of domains",
         "describe which combinations respondents chose",
         "latent class analysis", "Table 2, Fig 3, Table D in S1 Appendix"),
        ("2", "Combination or number",
         "separate the combination from the number",
         "Eq. (3)", "Tables E and F in S1 Appendix"),
        ("3", "Primary model",
         "estimate the shape of the association",
         "Eq. (4)", "Table 3, Fig 4a"),
        ("4", "Commitment intensity",
         "see whether how strongly domains are held explains it",
         "Eq. (5)", "Table 3, Fig 4a"),
        ("5", "Components of well-being",
         "locate which component carries the association",
         "Eq. (6)", "Table G in S1 Appendix"),
        ("6", "Family cancer history",
         "assess whether the association differs by family cancer history",
         "Eq. (7)", "Table 4, Fig 4b"),
        ("7", "Disaster exposure",
         "assess whether the association differs by disaster damage",
         "Eq. (7)", "Table 5, Fig 4c"),
        ("8", "Psychological distress",
         "assess whether the association is separable from distress",
         "Eq. (8)", "Table 3, Fig 4a"),
    ]
    height = (Y_RANGE[1] - Y_RANGE[0]) * Y_IN
    fig = plt.figure(figsize=(WIDTH_IN, height))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 100)
    ax.set_ylim(*Y_RANGE)
    ax.axis("off")

    # Boxes span the width; the arrows run down the middle of the boxes.
    x0, x1 = 1.5, 98.5
    xl, xr, xa = x0 + 3.0, x1 - 3.0, (x0 + x1) / 2
    # The gap after step 2 is wider: the arrow leaving it carries the
    # comparison that moved the exposure, and that label needs room of its own.
    h, gap, extra = 7.7, 10.8, 4.8
    ys, y = [], 96.0 - h / 2
    for i in range(len(steps)):
        if i:
            y -= gap + (extra if i == 2 else 0)
        ys.append(y)
    for i, (num, name, purpose, how, where) in enumerate(steps):
        y = ys[i]
        ax.add_patch(FancyBboxPatch((x0, y - h / 2), x1 - x0, h,
                                    boxstyle="round,pad=0.3,rounding_size=0.8",
                                    fc="white", ec="black", lw=0.9, zorder=2))
        ax.text(xl, y + 2.2, "Step %s   %s" % (num, name), fontsize=STEP_PT,
                ha="left", va="center", zorder=3)
        ax.text(xl, y - 1.6, "to " + purpose, fontsize=PURPOSE_PT, ha="left",
                va="center", color=GREY, zorder=3)
        ax.text(xr, y + 2.2, how, fontsize=HOW_PT, ha="right", va="center",
                zorder=3)
        ax.text(xr, y - 1.6, where, fontsize=WHERE_PT, ha="right", va="center",
                color=GREY, zorder=3)
        if i:
            a, b = ys[i - 1] - h / 2 - 0.3, y + h / 2 + 0.4
            ax.add_patch(FancyArrowPatch((xa, a), (xa, b),
                                         arrowstyle="-|>", mutation_scale=9,
                                         lw=0.9, color="black", zorder=1))
            if i == 2:                   # the arrow that leaves step 2
                ax.text(xa + 2.5, (a + b) / 2,
                        "the three non-Minimal classes differed little from one "
                        "another;\nthe number of "
                        "endorsed domains is the exposure from here on",
                        fontsize=ARROW_PT, ha="left", va="center",
                        color=GREY, linespacing=1.5, zorder=4)
    # no text may run past the box it belongs to, or into the text beside it
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    bx0, bx1 = ax.transData.transform([(x0, 0), (x1, 0)])[:, 0]
    lines = {}
    for t in ax.texts:
        e = t.get_window_extent(r)
        assert bx0 <= e.x0 and e.x1 <= bx1, ("text leaves its box", t.get_text())
        lines.setdefault(round(t.get_position()[1], 3), []).append(e)
    for row in lines.values():
        row.sort(key=lambda e: e.x0)
        for left, right in zip(row, row[1:]):
            assert left.x1 + 0.1 * fig.dpi < right.x0, "texts on one line come within 0.1 in"
    save(fig, "Fig2")


if __name__ == "__main__":
    flow_diagram()
    print("wrote Fig2 (analysis steps) to %s" % OUTD)
