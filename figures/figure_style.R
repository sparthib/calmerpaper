# =============================================================================
# figure_style.R  --  shared style guide for all lrs_cov figures
#
# source() this at the top of every figure script:
#     source(file.path(lrs_cov_dir, "figures/figure_style.R"))
#
# Encodes the auto-enforceable rules from the project figure style guide
# (Illinois "Style Guide for Figures in Science Writing" + the stricter
# publication-grade checklist). Rules that cannot be enforced by a theme
# (units in axis labels, claim-titles being true, stating n per panel,
# <=3 annotations per panel) are a manual checklist in figures/README.md.
#
# What this file provides:
#   theme_lrs()            - shared ggplot theme (sans-serif, role-mapped sizes,
#                            outward ticks, frameless legend). Font sizes map
#                            to ROLE, not space: base / one step down / two
#                            steps down. Minimum 5 pt at default scale.
#   KIT_COLORS             - canonical kit palette (Fig 1 dark PCR ramp wins).
#   PROTOCOL_COLORS        - canonical protocol palette.
#   ISOFORM_COUNT_COLORS   - sequential green ramp for isoform-count bins.
#   TX_LENGTH_COLORS       - categorical tx-length palette.
#   SPIKE_TYPE_COLORS      - spike-in type palette (Fig 3).
#   KIT_LINETYPES / KIT_SHAPES - per-kit linetype + shape so series are
#                            distinguishable WITHOUT relying on color alone.
#   scale_kit_*()          - convenience scales binding color+linetype+shape.
#   scale_fill_log2fc()    - diverging fill centered at semantic zero (0) for
#                            log2 fold-change heatmaps (NOT the data midpoint).
#   italicize_genes()      - helper to italicize gene/species names in labels.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
})

# ---- Font ------------------------------------------------------------------
# Portable sans-serif. Calibri (per the Illinois guide) is not installed on
# JHPCE and would warn-and-fallback, so we use R's generic "sans".
LRS_FONT <- "sans"

# ---- Role-mapped size ladder (>= 5 pt) -------------------------------------
# Four roles: title (bold, one step up from axes), base (axis titles / legend
# title / series identity), one step down (legend text / annotation), two steps
# down (tick labels). Titles use the same sans family as the axes, just bigger
# and bold.
LRS_SIZE_TITLE  <- 14  # plot titles: same sans family as axes, bigger + bold
LRS_SIZE_BASE   <- 12  # axis titles, legend title, series identity
LRS_SIZE_MID    <- 10  # legend text, annotations
LRS_SIZE_TICK   <- 8   # tick labels (>= 5 pt satisfied)
LRS_SIZE_TAG    <- 14  # panel letters (bold)

# ---- Shared theme ----------------------------------------------------------
theme_lrs <- function(base_size = LRS_SIZE_BASE, base_family = LRS_FONT) {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      text            = element_text(family = base_family),
      plot.title      = element_text(size = LRS_SIZE_TITLE, face = "bold",
                                     hjust = 0),
      plot.subtitle   = element_text(size = LRS_SIZE_MID, colour = "grey40"),
      axis.title      = element_text(size = LRS_SIZE_BASE),
      axis.text       = element_text(size = LRS_SIZE_TICK),
      legend.title    = element_text(size = LRS_SIZE_BASE),
      legend.text     = element_text(size = LRS_SIZE_MID),
      strip.text      = element_text(size = LRS_SIZE_BASE, face = "bold"),
      # outward ticks, frameless legend
      axis.ticks      = element_line(colour = "grey30"),
      axis.ticks.length = unit(2.5, "pt"),
      legend.key      = element_blank(),
      legend.background = element_blank(),
      panel.grid.minor = element_blank(),
      plot.tag        = element_text(size = LRS_SIZE_TAG, face = "bold")
    )
}

# ---- Canonical palettes ----------------------------------------------------
# Kit palette: Fig 1 "dark PCR ramp" is canonical (PCR kits stay visually
# dominant). PCB109 keeps the PCS109 shade it had in Fig 2/3.
KIT_COLORS <- c(
  RNA002         = "#7ba3f0",
  RNA004         = "#9bb8f3",
  R2C2           = "#f08a8a",
  captrap_linker = "#ffcc99",
  DCS109         = "#90ee90",
  PCS109         = "#d1b3f3",  # PCR shade (lightest)
  PCS110         = "#a47adf",  # PCR shade (mid)   -- canonical (was #e0ccf5 in Fig2/3)
  PCB114         = "#5e2ba1",  # PCR shade (darkest) -- canonical (was #c299f0 in Fig2/3)
  PCB109         = "#d1b3f3"   # alias shade of PCS109
)

PROTOCOL_COLORS <- c(
  "direct RNA"  = "#2563eb",
  "R2C2"        = "#dc2626",
  "captrap"     = "#059669",
  "direct cDNA" = "#7c2d12",
  "PCR cDNA"    = "#9333ea"
)

ISOFORM_COUNT_COLORS <- c(
  "1"     = "#f7fcf5",
  "2"     = "#e5f5e0",
  "3"     = "#c7e9c0",
  "4"     = "#a1d99b",
  "5"     = "#74c476",
  "6-10"  = "#41ab5d",
  "11-20" = "#238b45",
  "21-30" = "#006d2c",
  "31-40" = "#00441b",
  "41-50" = "#002d12",
  "51-80" = "#001a0a"
)

TX_LENGTH_COLORS <- c(
  "<500"       = "#e41a1c",
  "500-1000"   = "#377eb8",
  "1000-2000"  = "#4daf4a",
  "2000-5000"  = "#984ea3",
  "5000-10000" = "#ff7f00",
  "10000+"     = "#a65628"
)

SPIKE_TYPE_COLORS <- c(
  "SIRV_Short"    = "#e41a1c",
  "SIRV_Long"     = "#377eb8",
  "ERCC"          = "#4daf4a",
  "Sequins_Mix_A" = "#984ea3",
  "Sequins_Mix_B" = "#ff7f00"
)

# ---- Per-kit linetype + shape (redundant encoding, not color-only) ---------
# Strict-guide rule: each series carries a unique line AND marker style so it
# is distinguishable in greyscale / for CVD readers, not by hue alone.
.kit_order <- names(KIT_COLORS)
KIT_LINETYPES <- setNames(
  rep(c("solid", "22", "42", "44", "13", "1343", "73", "2262", "solid"),
      length.out = length(.kit_order)),
  .kit_order)
KIT_SHAPES <- setNames(
  rep(c(16, 17, 15, 18, 3, 4, 8, 7, 1), length.out = length(.kit_order)),
  .kit_order)

# ---- Convenience scales ----------------------------------------------------
scale_color_kit <- function(...) scale_color_manual(values = KIT_COLORS, ...)
scale_fill_kit  <- function(...) scale_fill_manual(values = KIT_COLORS, ...)
scale_linetype_kit <- function(...) scale_linetype_manual(values = KIT_LINETYPES, ...)
scale_shape_kit    <- function(...) scale_shape_manual(values = KIT_SHAPES, ...)

scale_color_protocol <- function(...) scale_color_manual(values = PROTOCOL_COLORS, ...)
scale_fill_protocol  <- function(...) scale_fill_manual(values = PROTOCOL_COLORS, ...)

# Diverging fill for log2 fold-change: centered at semantic zero (0), a
# perceptually reasonable blue-white-red. NEVER center at the data midpoint.
scale_fill_log2fc <- function(limits = NULL, name = expression(log[2]~FC), ...) {
  scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b",
                       midpoint = 0, limits = limits, name = name, ...)
}

# ---- Gene/species italics helper -------------------------------------------
# Wraps a vector of labels so ggplot renders them italic (parse = TRUE contexts)
# or returns an element_text for axis theming. Use where gene names appear.
italicize_genes <- function() element_text(face = "italic")

# ---- Convenience: apply-and-announce ---------------------------------------
message("Loaded figure_style.R: theme_lrs(), KIT_COLORS (dark PCR ramp), ",
        "scale_*_kit(), scale_fill_log2fc(). Font = '", LRS_FONT, "'.")
