# -*- coding: utf-8 -*-
"""05_figures_python.py -- the two manuscript figures and the figures-only
Word document.

    python 05_figures_python.py          Fig. 3 and Fig. 4
    python 05_figures_python.py fig4     Fig. 4 only

Reads only files already written by the R pipeline. Nothing is re-estimated.

  in   outputs/table2_class_profile.csv      conditional endorsement probabilities
       outputs/table27_endorsement_rates.csv observed rate per domain, from data
       outputs/table22_breadth_ref1.csv      breadth coefficients, reference = 1
                                             (n per level, for the series file)
       outputs/table31_disaster_group.csv    the same, by disaster exposure
       outputs/table37_sd_units.csv          Fig. 4 values, in SD units
                                             (10_sd_units_and_contrasts.R)
       outputs/table38_sd_wb.csv             the divisor, SD of overall
                                             well-being in the primary model

  out  outputs/figures_plos/Fig3.tif, Fig4.tif     the figures
       outputs/figures_plos/Fig3_series.csv        the plotted numbers
       outputs/figures_plos/Fig4_series.csv

Output follows the PLOS ONE figure requirements: a TIFF with LZW compression,
RGB at 8 bits per channel and no alpha channel, 600 dpi, 7.5 inches wide, text
in Arial at 8 pt or more at the printed size, and no figure number, title or
legend inside the image. Each figure is drawn at its printed size, so the type
sizes below are the printed sizes.

Figure contract.
  Fig 3  claim: the four endorsement classes differ in how many domains they
         take in far more than in which ones, and every class tracks the same
         ranking of domains that the sample as a whole shows.
  Fig 4  claim: well-being rises with the number of endorsed value domains and
         stops rising at about four, the rise is untouched by adjustment for
         commitment intensity, it is attenuated by adjustment for concurrent
         distress, and the two family cancer groups are not far apart.
"""
import csv
import io
import os
import re
import sys

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

from docx import Document
from docx.shared import Pt, Inches, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn

HERE = os.path.dirname(os.path.abspath(__file__))
OUTD = os.path.join(HERE, "outputs")
FIGD = os.path.join(OUTD, "figures_plos")
os.makedirs(FIGD, exist_ok=True)
sys.path.insert(0, os.path.join(HERE, "manuscript"))
# figure_legends belongs to the retired document builder. It is imported only
# by figure_document(), which is no longer run: from v1 onward the figures
# document is produced by editing the previous version, not by regenerating it.

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "font.size": 7,
    "axes.titlesize": 8,
    "axes.labelsize": 7.5,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 6.6,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.linewidth": 0.7,
    "xtick.major.width": 0.7,
    "ytick.major.width": 0.7,
    "lines.linewidth": 1.1,
    "legend.frameon": False,
    "figure.dpi": 150,
})

# One restrained palette for both figures. Every series also carries its own
# marker and line style, so the figures survive greyscale printing and the
# common forms of colour blindness, as the submission guidelines require.
BLUE = "#1f5c8b"
ORANGE = "#c46a26"
GREEN = "#2b8069"
PLUM = "#8f3b63"
GREY = "#5c5c5c"
BAND = "#e8eef4"
# x range of the shaded band in panel a: 4, 5 and 6+ endorsed domains
BAND_X = (3.6, 6.45)

# Both panels of Fig. 2 share one vertical scale. Readers compare the gradients
# across the panels, and two scales would make the distress-adjusted series look
# as steep as the primary one. The range covers every interval either panel
# draws, so nothing is clipped.
# The lower bound clears the longest interval any panel draws -- the disaster
# damage series at no endorsed domains, whose lower limit is -0.1185 -- and
# leaves a clear band below it for the "reference" note, so the note sits on
# neither the data nor the axis.
YLIM = (-0.163, 0.212)
YTICKS = [-0.15, -0.10, -0.05, 0.00, 0.05, 0.10, 0.15, 0.20]
# Fig. 4 is drawn in SD units: every coefficient and limit divided by one
# divisor, the model-implied SD of overall well-being in the primary model.
# The axis limits and the note positions above are divided by the same number,
# so the drawing is the coefficient-scale drawing with a relabelled axis.
YTICKS_SD = [-0.75, -0.50, -0.25, 0.00, 0.25, 0.50, 0.75, 1.00]
NOTE_A_Y = -0.112                # panel a note, in coefficient units

# the seven levels of breadth; one endorsed domain is the reference
LEVELS = ["0", "1", "2", "3", "4", "5", "6+"]
XPOS = {lv: i for i, lv in enumerate(LEVELS)}
REF = "1"


def read(name, folder=OUTD):
    with open(os.path.join(folder, name), encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh))


def write_series(name, rows, header):
    with open(os.path.join(FIGD, name), "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)


def draw(ax, xs, b, lo, hi, colour, marker, ls, label, dx=0.0):
    xs = [v + dx for v in xs]
    ax.errorbar(xs, b, yerr=[[bi - l for bi, l in zip(b, lo)],
                             [h - bi for bi, h in zip(b, hi)]],
                fmt="none", ecolor=colour, elinewidth=0.8, capsize=1.6,
                capthick=0.8, zorder=2)
    ax.plot(xs, b, color=colour, linestyle=ls, marker=marker, markersize=3.4,
            markerfacecolor="white", markeredgewidth=0.9, label=label,
            zorder=3, clip_on=False)


def draw_model(ax, rows, model, colour, marker, ls, label, dx=0.0,
               mark_reference=True):
    """One model's coefficients, with the reference level pinned at zero."""
    d = {r["level"]: r for r in rows if r["model"] == model}
    levels = [lv for lv in LEVELS if lv in d or lv == REF]
    xs, b, lo, hi = [], [], [], []
    for lv in levels:
        xs.append(XPOS[lv])
        if lv == REF:
            b.append(0.0)
            lo.append(0.0)
            hi.append(0.0)
        else:
            b.append(float(d[lv]["b"]))
            lo.append(float(d[lv]["ci_low"]))
            hi.append(float(d[lv]["ci_high"]))
    draw(ax, xs, b, lo, hi, colour, marker, ls, label, dx=dx)
    if mark_reference:
        ax.plot([XPOS[REF] + dx], [0.0], marker=marker, markersize=3.4,
                markerfacecolor=colour, markeredgecolor=colour,
                markeredgewidth=0.9, linestyle="none", zorder=4)
    return [(model, lv, d[lv]["n"] if lv in d else "",
             "reference" if lv == REF else d[lv]["b"],
             "" if lv == REF else d[lv]["ci_low"],
             "" if lv == REF else d[lv]["ci_high"],
             "" if lv == REF else d[lv]["p"]) for lv in levels]


def panel_label(ax, text, dx=0.0):
    ax.text(-0.15 + dx, 1.04, text, transform=ax.transAxes,
            fontsize=PANEL_PT, fontweight="bold", va="bottom", ha="left")


def finish_x(ax, note, ha="left"):
    ax.set_xticks(list(XPOS.values()))
    ax.set_xticklabels(LEVELS)
    ax.set_xlim(-0.45, 6.45)
    ax.set_xlabel("Number of endorsed domains")
    ax.axhline(0, color="#9a9a9a", linewidth=0.6, linestyle=(0, (1, 2)),
               zorder=1)
    return ax.text(0.01 if ha == "left" else 0.99, 0.015, note,
                   transform=ax.transAxes, fontsize=NOTE_PT, color=GREY,
                   ha=ha, va="bottom", linespacing=1.3)


# PLOS ONE figure files
DPI = 600
WIDTH_IN = 7.5                  # 2.63 to 7.5 inches wide
MAX_HEIGHT_IN = 8.75
MIN_PT = 8.0
FONT_FILES = ("Arial",)         # embedded PostScript names must start with this


def check_figure(fig, stem):
    """Measure what was drawn: nothing outside the canvas, every text at 8 pt
    or more, and every glyph set in the named typeface (read back from the
    fonts a PDF rendering embeds, so a fallback to another font is caught)."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    w, h = fig.get_size_inches()
    bb = fig.get_tightbbox(r)
    assert abs(w - WIDTH_IN) < 1e-9 and h <= MAX_HEIGHT_IN, (w, h)
    assert bb.x0 >= 0 and bb.y0 >= 0 and bb.x1 <= w + 1e-6 and bb.y1 <= h + 1e-6, \
        "something is drawn outside the canvas: %s" % bb
    texts = [t for t in fig.findobj(mpl.text.Text)
             if t.get_visible() and t.get_text().strip()
             and t.get_window_extent(r).width > 0]
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
    path = os.path.join(FIGD, "%s.tif" % stem)
    im.save(path, format="TIFF", compression="tiff_lzw", dpi=(DPI, DPI))
    assert im.size[0] == round(WIDTH_IN * DPI), im.size
    print("  wrote %s (%d x %d px)" % (path, im.size[0], im.size[1]))


# Every figure is drawn at its printed width, 7.5 inches, and placed at the
# size drawn, so the type sizes below are the printed sizes, all 8 pt or more.
# Constrained layout fits the axes, legends and labels inside that width.
FIG3_SIZE = (WIDTH_IN, 4.9)
FIG3_RC = {"font.size": 8.75,
           "axes.labelsize": 9.375,
           "xtick.labelsize": 8.75,
           "ytick.labelsize": 8.75,
           "legend.fontsize": 8.25}

# Fig 4 carries three panels side by side on one vertical scale.
FIG4_SIZE = (WIDTH_IN, 3.6)
FIG4_RC = {"font.size": 8.5,
           "axes.labelsize": 9.0,
           "xtick.labelsize": 8.5,
           "ytick.labelsize": 8.5,
           "legend.fontsize": 8.0}
NOTE_PT = 8.0
PANEL_PT = 11.0


# ------------------------------------------------------------------- Fig. 3
def figure_endorsement():        # Fig. 3
    prof = read("table2_class_profile.csv")
    rate = read("table27_endorsement_rates.csv")
    domains = [r["domain"] for r in rate]
    labels = [r["label"] for r in rate]
    observed = [float(r["p_endorsed"]) for r in rate]

    styles = [(BLUE, "o", "-"), (ORANGE, "s", "--"), (GREEN, "^", ":"),
              (PLUM, "D", "-.")]
    xs = list(range(len(domains)))
    series = []
    ctx = mpl.rc_context(FIG3_RC)
    ctx.__enter__()
    fig, ax = plt.subplots(figsize=FIG3_SIZE, layout="constrained")
    for r, (colour, marker, ls) in zip(prof, styles):
        y = [float(r["p_" + d]) for d in domains]
        ax.plot(xs, y, color=colour, linestyle=ls, marker=marker,
                markersize=3.6, markerfacecolor="white", markeredgewidth=0.9,
                label="%s %s (n = %s, mean %s domains)"
                      % (r["class"], r["label"], r["n"], r["mean_n_domains"]))
        series += [("%s %s" % (r["class"], r["label"]), lab, "%.3f" % v)
                   for lab, v in zip(labels, y)]
    # the observed rate in the whole sample, counted from the analysis data
    # rather than approximated by a class-share weighted average
    ax.plot(xs, observed, color=GREY, linestyle=(0, (4, 2)), marker="x",
            markersize=4, markeredgewidth=0.9,
            label="Observed rate in the whole sample (N = 2179)")
    series += [("Observed rate, whole sample", lab, "%.3f" % v)
               for lab, v in zip(labels, observed)]
    ax.set_xticks(xs)
    ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.set_ylabel("Probability of endorsement")
    ax.set_ylim(-0.03, 1.05)
    ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    # the legend sits above the panel, in two columns, inside the 7.5-inch width
    fig.legend(*ax.get_legend_handles_labels(), loc="outside upper center",
               ncol=2, handlelength=2.6, labelspacing=0.35, columnspacing=1.6)
    save(fig, "Fig3")
    ctx.__exit__(None, None, None)
    write_series("Fig3_series.csv", series, ["series", "domain", "value"])


# ------------------------------------------------------------------- Fig. 4
def sd_divisor():
    """SD of overall well-being in the primary model, as 10_ wrote it."""
    rows = [r for r in read("table38_sd_wb.csv")
            if r["model"].startswith("Primary model")]
    assert len(rows) == 1
    return float(rows[0]["sd_wb"])


# table37 (table, column) -> the model labels used in table22 and table31
SD_MODELS = {
    ("Table 3", "Primary model"): "Primary model (ref = 1 domain)",
    ("Table 3", "Adjusted for commitment intensity"):
        "Adjusted for commitment intensity (ref = 1 domain)",
    ("Table 3", "Adjusted for psychological distress"):
        "Adjusted for psychological distress (ref = 1 domain)",
    ("Table 4", "No family cancer history"):
        "No family cancer history (ref = 1 domain)",
    ("Table 4", "Family cancer history"):
        "Family cancer history (ref = 1 domain)",
    ("Table 5", "No disaster damage"): "No disaster damage (ref = 1 domain)",
    ("Table 5", "Disaster damage"): "Disaster damage (ref = 1 domain)",
}


def sd_unit_rows():
    """Fig. 4 rows in SD units, at full precision, with n from the b tables."""
    n_of = {(r["model"], r["level"]): r["n"]
            for r in read("table22_breadth_ref1.csv")
            + read("table31_disaster_group.csv")}
    out = []
    for r in read("table37_sd_units.csv"):
        model = SD_MODELS.get((r["table"], r["column"]))
        if model is None:
            continue
        out.append({"model": model, "level": r["row"],
                    "n": n_of[(model, r["row"])], "b": r["sd"],
                    "ci_low": r["sd_low"], "ci_high": r["sd_high"],
                    "p": r["p"]})
    assert len(out) == 6 + 5 + 6 + 6 + 6 + 6 + 6, len(out)
    return out


def check_fig4(fig, axes, note_a, ref_notes, series, ylim=YLIM):
    """Measure what was actually drawn and assert the three layout rules.

    Rendered extents are read back from the figure, not assumed: the checker
    fires on a real violation rather than on an intention.
    """
    YLIM = ylim
    fig.canvas.draw()
    rep = []

    def data_box(t, ax):
        bb = t.get_window_extent(fig.canvas.get_renderer())
        return ax.transData.inverted().transform(bb).tolist()

    # 1. every interval, including its cap, lies inside the vertical scale
    lo = min(float(r[4]) for r in series if r[4] != "")
    hi = max(float(r[5]) for r in series if r[5] != "")
    rep.append("intervals span %.4f to %.4f; axis is %.3f to %.3f"
               % (lo, hi, *YLIM))
    assert YLIM[0] < lo and hi < YLIM[1], "an interval is clipped by the axis"

    # 2. the panel a note sits wholly inside the shaded band and the axes box
    (x0, y0), (x1, y1) = data_box(note_a, axes[0])
    rep.append("panel a note spans x %.2f to %.2f (band %.2f to %.2f), "
               "y %.4f to %.4f" % (x0, x1, *BAND_X, y0, y1))
    assert BAND_X[0] <= x0 and x1 <= BAND_X[1], "the note leaves the band"
    assert axes[0].get_xlim()[0] <= x0 and x1 <= axes[0].get_xlim()[1], \
        "the note leaves the panel"
    assert YLIM[0] <= y0 and y1 <= YLIM[1], "the note leaves the panel"

    # 3. no reference note touches a drawn interval or the frame
    for i, (t, ax) in enumerate(zip(ref_notes, axes)):
        (x0, y0), (x1, y1) = data_box(t, ax)
        near = [float(r[4]) for r in series if r[4] != ""
                and x0 - 0.35 <= XPOS[r[1]] <= x1 + 0.35]
        floor = min(near) if near else hi
        rep.append("panel %s reference note: y %.4f to %.4f, lowest interval "
                   "above it %.4f, axis floor %.3f"
                   % ("abc"[i], y0, y1, floor, YLIM[0]))
        assert y1 < floor, "the reference note runs into the data"
        assert y0 > YLIM[0], "the reference note sits on the axis"
    print("  Fig. 4 layout checks passed:")
    for line in rep:
        print("    " + line)


def figure_coefficients():       # Fig. 4
    # SD units: one divisor for every panel and series, so the three panels
    # still share one vertical scale
    sdwb = sd_divisor()
    rows = sd_unit_rows()
    ylim = (YLIM[0] / sdwb, YLIM[1] / sdwb)
    # The three panels share one vertical scale, so the axis is labelled and
    # ticked once on the left. That frees the width the legends need at 8 pt.
    with mpl.rc_context(FIG4_RC):
        fig, axes = plt.subplots(1, 3, figsize=FIG4_SIZE, sharey=True,
                                 layout="constrained",
                                 gridspec_kw={"width_ratios": [1.0, 1.0, 1.0]})
        series = []

        ax = axes[0]
        ax.axvspan(*BAND_X, color=BAND, zorder=0, linewidth=0)
        # The note sits wholly inside the shaded band, in the empty space below
        # the series. It is broken over four lines so that its width fits
        # between the edges of the band (3.6 to 6.45) rather than running out
        # past the right edge of the panel. `check_fig4` asserts both.
        # From v45 the note carries no test: it describes the estimates, as the
        # Table 3 footnote does with the differences between these levels.
        note_a = ax.text(sum(BAND_X) / 2, NOTE_A_Y / sdwb, "estimates\nchange little\n"
                         "across these\nthree levels",
                         fontsize=NOTE_PT, color=GREY, ha="center",
                         va="bottom", linespacing=1.3)
        # Panel a carries the three columns of Table 3. The three series are
        # offset sideways so that their error bars do not overlap, and are
        # separated by plotting symbol and line type as well as by colour.
        series += draw_model(ax, rows, "Primary model (ref = 1 domain)", BLUE,
                             "o", "-", "Primary model", dx=-0.11)
        series += draw_model(
            ax, rows, "Adjusted for commitment intensity (ref = 1 domain)",
            ORANGE, "s", "--", "Adjusted for commitment intensity", dx=0.0)
        series += draw_model(
            ax, rows, "Adjusted for psychological distress (ref = 1 domain)",
            GREEN, "^", ":", "Adjusted for psychological distress", dx=0.11)
        ax.set_ylabel("Difference in well-being (SD units)")
        ax.set_ylim(*ylim)
        ax.set_yticks(YTICKS_SD)
        ref_notes = [finish_x(ax, "reference: 1 endorsed domain")]
        ax.legend(loc="upper left", bbox_to_anchor=(-0.02, 1.02),
                  handlelength=2.6, labelspacing=0.35)
        panel_label(ax, "a")

        ax = axes[1]
        # Panel b carries Table 4, the primary model fitted separately in the
        # two family cancer groups.
        series += draw_model(ax, rows,
                             "No family cancer history (ref = 1 domain)",
                             GREY, "D", "--", "No family cancer history",
                             dx=-0.06)
        series += draw_model(ax, rows, "Family cancer history (ref = 1 domain)",
                             PLUM, "v", "-.", "Family cancer history", dx=0.06)
        ax.set_ylim(*ylim)
        ax.set_yticks(YTICKS_SD)
        ref_notes.append(finish_x(ax, "reference: 1 endorsed domain"))
        ax.legend(loc="upper left", bbox_to_anchor=(-0.02, 1.02),
                  handlelength=2.6, labelspacing=0.35)
        panel_label(ax, "b", dx=0.11)

        ax = axes[2]
        series += draw_model(ax, rows, "No disaster damage (ref = 1 domain)",
                             GREY, "D", "--", "No disaster damage", dx=-0.06)
        series += draw_model(ax, rows, "Disaster damage (ref = 1 domain)",
                             BLUE, "o", "-", "Disaster damage", dx=0.06)
        ax.set_ylim(*ylim)
        ax.set_yticks(YTICKS_SD)
        ref_notes.append(finish_x(ax, "reference: 1 endorsed domain"))
        ax.legend(loc="upper left", bbox_to_anchor=(-0.02, 1.02),
                  handlelength=2.6, labelspacing=0.35)
        panel_label(ax, "c", dx=0.11)

        check_fig4(fig, axes, note_a, ref_notes, series, ylim=ylim)
        save(fig, "Fig4")
    write_series("Fig4_series.csv", series,
                 ["model", "domains endorsed", "n", "difference in SD units",
                  "ci_low", "ci_high",
                  "p"])


# --------------------------------------------------------------- figure docx
FONT = "Times New Roman"


def figure_document():
    """One figure per page, legend directly under the figure, Normal style."""
    doc = Document()
    st = doc.styles["Normal"]
    st.font.name = FONT
    st.font.size = Pt(12)
    st.font.color.rgb = RGBColor(0, 0, 0)
    st.element.rPr.rFonts.set(qn("w:eastAsia"), FONT)
    st.paragraph_format.space_after = Pt(0)
    for s in doc.sections:
        s.top_margin = s.bottom_margin = Inches(1)
        s.left_margin = s.right_margin = Inches(1)

    def put(run, bold=False):
        run.font.name = FONT
        run.font.size = Pt(12)
        run.bold = bold
        run.font.color.rgb = RGBColor(0, 0, 0)
        run._element.rPr.rFonts.set(qn("w:eastAsia"), FONT)

    for i, (img, title, legend) in enumerate(
            [("Fig1.png", FIG1_TITLE, FIG1_LEGEND),
             ("Fig2.png", FIG2_TITLE, FIG2_LEGEND),
             ("Fig3.png", FIG3_TITLE, FIG3_LEGEND)]):
        if i:
            doc.add_page_break()
        p = doc.add_paragraph()
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        p.add_run().add_picture(os.path.join(FIGD, img), width=Inches(6.4))
        p.paragraph_format.space_after = Pt(10)
        p = doc.add_paragraph()
        put(p.add_run(title), bold=True)
        p.paragraph_format.line_spacing = 1.2
        p.paragraph_format.space_after = Pt(4)
        p = doc.add_paragraph()
        put(p.add_run("\u2733 " + legend))
        p.paragraph_format.line_spacing = 1.2
    out = os.path.join(FIGD, "Figure_value.docx")
    doc.save(out)
    return out


if __name__ == "__main__":
    if sys.argv[1:] != ["fig4"]:
        figure_endorsement()
    figure_coefficients()
    print("wrote %s to %s" % ("Fig4 (coefficients, SD units)" if sys.argv[1:] == ["fig4"]
                              else "Fig3 (endorsement) and Fig4 (coefficients, SD units)", FIGD))
    print("figure_document() is not run: the figures document is carried "
          "forward from the previous version, not regenerated")
