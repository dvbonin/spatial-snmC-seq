# Shared font and export settings for the preprint figures (ggplot2 side).
#
# Text stays as text in the vector outputs, so labels can be moved and retyped in
# Illustrator rather than arriving as glyph outlines. That needs specific devices:
#
#   cairo_pdf   goes through fontconfig, takes any system font, embeds a subset with
#               real text. Base pdf() cannot: it only knows its own Type 1 families
#               and errors with "invalid font type" on anything else, Arial included.
#   svglite     writes <text> referencing the font by name. The cairo svg device
#               always converts glyphs to paths, which is what made the earlier
#               SVGs uneditable.
#
# Arial is installed in ~/.fonts. FIGURE_FONT is the single place to change it.

FIGURE_FONT <- "Arial"

# svglite keeps its base fill/stroke in a CSS class block and the per-element colours in
# `style` attributes. Affinity parses neither reliably: it ignores the stylesheet, and it
# lets presentation attributes win over `style` (backwards from the spec). So rewrite each
# shape's effective style as presentation attributes -- class-block defaults overridden by
# that element's own declarations. The `style` attribute is left in place, so renderers
# that follow the spec see exactly what they saw before.
SVG_BASE <- c(fill = "none", stroke = "#000000", "stroke-linecap" = "round",
              "stroke-linejoin" = "round", "stroke-miterlimit" = "10.00")
SVG_SHAPES <- "<(line|polyline|polygon|path|rect|circle)\\b[^>]*>"

# Text is deliberately excluded: it has no presentation attributes of its own, so Affinity
# already falls through to its `style` and renders the right font and size.
svg_attrs_from_style <- function(tag) {
  props <- SVG_BASE
  if (grepl("style='", tag, fixed = TRUE)) {
    decl <- trimws(strsplit(sub(".*style='([^']*)'.*", "\\1", tag), ";")[[1]])
    decl <- decl[nzchar(decl)]
    props[trimws(sub(":.*", "", decl))] <- trimws(sub("^[^:]*:", "", decl))
  }
  sub("^<([a-z]+)", paste0("<\\1 ", paste0(names(props), "='", props, "'", collapse = " ")), tag)
}

inline_svg_defaults <- function(svg_path) {
  txt <- paste(readLines(svg_path, warn = FALSE), collapse = "\n")
  hits <- gregexpr(SVG_SHAPES, txt)
  regmatches(txt, hits) <- list(vapply(regmatches(txt, hits)[[1]],
                                       svg_attrs_from_style, character(1)))
  writeLines(txt, svg_path)
}

# PNG for preview, PDF for editing, SVG as well when SNMC_SVG is set. Pass svg = FALSE
# for a panel whose vector form is too dense to be useful (PubSeq: ~240 MB).
save_figure <- function(plot, png_path, width, height, dpi = 150, svg = TRUE, ...) {
  dir.create(dirname(png_path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(png_path, plot, width = width, height = height, dpi = dpi, ...)
  ggplot2::ggsave(sub("\\.png$", ".pdf", png_path), plot, width = width, height = height,
                  device = grDevices::cairo_pdf, ...)
  if (svg && nzchar(Sys.getenv("SNMC_SVG"))) {
    svg_path <- sub("\\.png$", ".svg", png_path)
    ggplot2::ggsave(svg_path, plot, width = width, height = height,
                    # fix_text_size = FALSE: svglite otherwise pins each string's
                    # textLength, which makes Illustrator squeeze edited labels back to
                    # the original width instead of letting them reflow.
                    device = function(filename, ...) svglite::svglite(filename, fix_text_size = FALSE, ...),
                    ...)
    inline_svg_defaults(svg_path)
  }
  invisible(png_path)
}
