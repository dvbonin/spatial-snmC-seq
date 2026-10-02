#!/usr/bin/env python3
"""
Shared font and export settings for the preprint figures (matplotlib side).

Text stays as text in the vector outputs, so labels can be moved and retyped in
Illustrator rather than arriving as glyph outlines:

  pdf.fonttype = 42   TrueType in the PDF, not Type 3, which Illustrator edits badly
  svg.fonttype = none text as <text> referencing the font by name, not paths

Arial is installed in ~/.fonts, so it is used directly; Nimbus Sans (Helvetica metrics,
system-wide) and DejaVu Sans are fallbacks in case a machine lacks it.

save_figure() writes the PNG for preview and the PDF for editing. SVG is written too
when SNMC_SVG is set -- via the px -> pt patch below, because matplotlib emits
`font: 7px Arial`, a shorthand Illustrator misreads as the wrong size.
"""
import os
import re

import matplotlib

FONT_STACK = ["Arial", "Nimbus Sans", "DejaVu Sans"]

def apply_style():
    """Font and vector-text settings. Call once, before creating any figure."""
    matplotlib.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": FONT_STACK,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })


def _patch_svg(text):
    """`font: [weight ]SIZEpx FAMILY` -> explicit font-size in pt plus font-family."""
    text = re.sub(r"font: (\d+) ([\d.]+)px (.+?)(;|\")",
                  lambda m: f"font-weight: {m.group(1)}; font-size: {m.group(2)}pt; "
                            f"font-family: {m.group(3).strip()}{m.group(4)}", text)
    return re.sub(r"font: ([\d.]+)px (.+?)(;|\")",
                  lambda m: f"font-size: {m.group(1)}pt; "
                            f"font-family: {m.group(2).strip()}{m.group(3)}", text)


def save_figure(fig, png_path, **kwargs):
    """PNG (preview) + PDF (editable). SVG as well when SNMC_SVG is set."""
    os.makedirs(os.path.dirname(png_path) or ".", exist_ok=True)
    fig.savefig(png_path, **kwargs)
    fig.savefig(png_path.replace(".png", ".pdf"), **kwargs)
    if os.environ.get("SNMC_SVG"):
        svg = png_path.replace(".png", ".svg")
        fig.savefig(svg, **kwargs)
        with open(svg, encoding="utf-8") as f:
            patched = _patch_svg(f.read())
        with open(svg, "w", encoding="utf-8") as f:
            f.write(patched)
    return png_path
