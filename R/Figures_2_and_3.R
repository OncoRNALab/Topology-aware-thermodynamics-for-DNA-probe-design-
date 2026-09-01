###############################################################################
# Figure 2 and Figure 3 Generation Script
#
# "Topology-aware thermodynamics for DNA probe design under fixed stringency"
# 
#
# This single script generates both manuscript figures:
#   Figure 2: 6-panel empirical topology proof (A–F)
#   Figure 3: 4-panel allele-specific benchmark (A–D, portrait)
#
# Reads ONLY from supplementary_tables/:
#   Figure 2: S3, S4, S10, S11
#   Figure 3: S13, S14, S16
#
# Requirements: ggplot2, ggprism, patchwork, dplyr, tidyr, readr, svglite
#
# Usage:
#   setwd("/path/to/ECI_02/")
#   source("R/Figures_2_and_3.R")
#
# Output: Fig. 2.pdf/.png/.svg and Fig. 3.pdf/.png/.svg in the working directory
#
###############################################################################

# ---- Libraries ----
library(ggplot2)
library(ggprism)
library(patchwork)
library(dplyr)
library(tidyr)
library(readr)
library(svglite)

# ---- Shared Helpers ----

save_fig <- function(plot, filename, width = 180, height = 120, units = "mm") {
  ggsave(filename = paste0(filename, ".pdf"), plot = plot,
         width = width, height = height, units = units, device = cairo_pdf)
  ggsave(filename = paste0(filename, ".png"), plot = plot,
         width = width, height = height, units = units, dpi = 300)
  svglite::svglite(file = paste0(filename, ".svg"), width = width / 25.4,
                   height = height / 25.4)
  print(plot)
  dev.off()
  message("Saved: ", filename, ".pdf / .png / .svg")
}

###############################################################################
# FIGURE 2: 6-Panel Empirical Topology Proof
###############################################################################

make_figure2_6panel <- function() {

  # Color Palette
  c_navy  <- "#1B3B6F"
  c_coral <- "#E07A5F"
  c_slate <- "#8D99AE"
  c_grey  <- "#D3D3D3"

  # ===========================================================================
  # Panel A: Seringhaus/Gerstein (36-mer)
  # ===========================================================================
  s3 <- read_csv("supplementary_tables/Supplementary table_S3_seringhaus_reconstruction.csv",
                 show_col_types = FALSE) %>%
    filter(gene == "ACT1")

  p2a <- ggplot(s3, aes(x = n_mismatches, y = relative_intensity_loss,
                        color = scheme, linetype = scheme, shape = scheme)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.5) +
    scale_color_manual(
      values = c("centered" = c_navy, "staggered" = c_coral),
      labels = c("Adjacent/centered MM (~7% loss/MM)",
                 "Distributed/staggered MM (~21% loss/MM)")
    ) +
    scale_linetype_manual(
      values = c("centered" = "solid", "staggered" = "dashed"),
      labels = c("Adjacent/centered MM (~7% loss/MM)",
                 "Distributed/staggered MM (~21% loss/MM)")
    ) +
    scale_shape_manual(
      values = c("centered" = 16, "staggered" = 15),
      labels = c("Adjacent/centered MM (~7% loss/MM)",
                 "Distributed/staggered MM (~21% loss/MM)")
    ) +
    scale_x_continuous(breaks = seq(0, 12, by = 2), limits = c(0, 12)) +
    scale_y_continuous(breaks = seq(0, 1.0, by = 0.2), limits = c(0, 1.0)) +
    labs(
      title = "A. 36-mer same-count topology",
      x = "Number of mismatches",
      y = "Relative signal loss from PM",
      color = NULL, linetype = NULL, shape = NULL
    ) +
    theme_prism(base_size = 10) +
    theme(
      aspect.ratio = 0.85,
      legend.position = "inside",
      legend.position.inside = c(0.48, 0.82),
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3)
    )

  # ===========================================================================
  # Panel B: Deng/Zhou (50-mer)
  # ===========================================================================
  s4 <- read_csv("supplementary_tables/Supplementary table_S4_deng_even_random_aggregate.csv",
                 show_col_types = FALSE) %>%
    pivot_longer(
      cols = c(evenly_distributed_relative_signal, randomly_distributed_relative_signal),
      names_to = "distribution", values_to = "signal"
    )

  p2b <- ggplot(s4, aes(x = mismatches, y = signal,
                        color = distribution, linetype = distribution, shape = distribution)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.5) +
    geom_hline(yintercept = 0.1, linetype = "dotted", color = "gray30", linewidth = 0.6) +
    scale_color_manual(
      values = c("randomly_distributed_relative_signal" = c_navy,
                 "evenly_distributed_relative_signal" = c_coral),
      labels = c("Randomly distributed MM", "Evenly distributed MM")
    ) +
    scale_linetype_manual(
      values = c("randomly_distributed_relative_signal" = "solid",
                 "evenly_distributed_relative_signal" = "dashed"),
      labels = c("Randomly distributed MM", "Evenly distributed MM")
    ) +
    scale_shape_manual(
      values = c("randomly_distributed_relative_signal" = 16,
                 "evenly_distributed_relative_signal" = 15),
      labels = c("Randomly distributed MM", "Evenly distributed MM")
    ) +
    scale_x_continuous(breaks = 0:7, limits = c(0, 7.5)) +
    scale_y_continuous(breaks = seq(0, 1.0, by = 0.2), limits = c(0, 1.0)) +
    labs(
      title = "B. 50-mer placement distribution",
      x = "Number of mismatches",
      y = "Relative signal (MM / PM)",
      color = NULL, linetype = NULL, shape = NULL
    ) +
    theme_prism(base_size = 10) +
    theme(
      aspect.ratio = 0.85,
      legend.position = "inside",
      legend.position.inside = c(0.60, 0.82),
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3)
    )

  # ===========================================================================
  # Panel C: Quadratic Sub-box Growth Curves (Nbox vs L)
  # ===========================================================================
  calc_nbox <- function(L, k) {
    pmax(0, (L - k + 1) * (L - k + 2) / 2)
  }

  df_nbox <- expand.grid(L = 1:25, k = c(2, 3, 4, 5)) %>%
    mutate(
      Nbox = calc_nbox(L, k),
      k_label = paste0("k = ", k, " nt")
    )

  p2c <- ggplot(df_nbox, aes(x = L, y = Nbox, color = k_label, linetype = k_label)) +
    geom_line(linewidth = 0.9) +
    scale_color_manual(values = c("k = 2 nt" = c_navy, "k = 3 nt" = c_coral,
                                  "k = 4 nt" = "#2A9D8F", "k = 5 nt" = c_slate)) +
    scale_x_continuous(breaks = seq(0, 25, by = 5)) +
    labs(
      title = "C. Combinatorial sub-box growth",
      x = "Continuous paired island length (L, nt)",
      y = expression("Sub-box pool count, " * N[box] * "(L, k)"),
      color = "Min box threshold", linetype = "Min box threshold"
    ) +
    theme_prism(base_size = 10) +
    theme(
      aspect.ratio = 0.85,
      legend.position = "inside",
      legend.position.inside = c(0.35, 0.75),
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3)
    )

  # ===========================================================================
  # Panel D: Probe Design Rule Schematic
  # ===========================================================================
  seq_len <- 25
  df_clustered <- tibble(
    x = 1:seq_len, y = 2,
    type = factor("Matched Base", levels = c("Matched Base", "Mismatch (Disrupted)")),
    group = "Clustered: L = 25 nt (Nbox = 276)"
  )
  mismatch_pos <- c(5, 9, 13, 17, 21)
  df_distributed <- tibble(
    x = 1:seq_len, y = 1,
    type = factor(if_else(x %in% mismatch_pos, "Mismatch (Disrupted)", "Matched Base"),
                  levels = c("Matched Base", "Mismatch (Disrupted)")),
    group = "Fragmented: 4+3+3+3+4 nt (Nbox = 9)"
  )
  df_all <- bind_rows(df_clustered, df_distributed)

  p2d <- ggplot() +
    geom_tile(data = df_all, aes(x = x, y = y, fill = type),
              color = "white", linewidth = 0.6, height = 0.45) +
    geom_text(data = df_all %>% distinct(y, group),
              aes(x = 13, y = y + 0.38, label = group),
              fontface = "bold", size = 2.8) +
    scale_fill_manual(values = c("Matched Base" = c_grey, "Mismatch (Disrupted)" = c_coral)) +
    scale_y_continuous(limits = c(0.4, 2.7)) +
    scale_x_continuous(limits = c(0, 26)) +
    labs(
      title = "D. Probe-design topology logic",
      x = "Nucleotide position along probe (25-mer)",
      y = NULL
    ) +
    theme_prism(base_size = 10) +
    theme(
      aspect.ratio = 0.85,
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none"
    )

  # ===========================================================================
  # Panel E: Affymetrix Subset Correlations
  # ===========================================================================
  s10 <- read_csv("supplementary_tables/Supplementary table_S10_affymetrix_correlation_summary.csv",
                  show_col_types = FALSE) %>%
    filter(Predictor %in% c("S_ECI", "N_box_k3", "DeltaG37_kcal_mol")) %>%
    mutate(
      Subset = factor(Subset,
                      levels = c("2-4 mismatch subset", "strict 3-mismatch subset"),
                      labels = c("2-4 MM\nsubset (n=8)", "Strict 3-MM\nsubset (n=5)")),
      Predictor = factor(Predictor,
                         levels = c("S_ECI", "N_box_k3", "DeltaG37_kcal_mol"),
                         labels = c("ECI topology score (S_ECI)",
                                    "Box count (N_box, k=3)",
                                    "Scalar \u0394G37"))
    )

  p2e <- ggplot(s10, aes(x = Subset, y = Pearson_r, fill = Predictor)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6, color = "black", linewidth = 0.3) +
    scale_fill_manual(
      values = c(
        "ECI topology score (S_ECI)" = c_navy,
        "Box count (N_box, k=3)"     = c_coral,
        "Scalar \u0394G37"            = c_slate
      )
    ) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
    scale_y_continuous(limits = c(-0.6, 1.0), breaks = seq(-0.6, 1.0, by = 0.2)) +
    labs(
      title = "E. Affymetrix benchmark validation",
      x = NULL,
      y = "Pearson r with log2 intensity",
      fill = NULL
    ) +
    theme_prism(base_size = 10) +
    theme(
      aspect.ratio = 0.85,
      legend.position = "inside",
      legend.position.inside = c(0.35, 0.22),
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3)
    )

  # ===========================================================================
  # Panel F: HPV Diagnostic Probe Audit
  # ===========================================================================
  s11 <- read_csv("supplementary_tables/Supplementary table_S11_HPV_edge_cases.csv",
                  show_col_types = FALSE) %>%
    filter(!is.na(delta_s_eci), !is.na(exploratory_delta_delta_g37))

  p2f <- ggplot(s11, aes(x = exploratory_delta_delta_g37, y = delta_s_eci)) +
    geom_hline(yintercept = 0, color = "gray50", linetype = "dashed", linewidth = 0.4) +
    geom_vline(xintercept = 0, color = "gray50", linetype = "dashed", linewidth = 0.4) +
    geom_point(shape = 16, size = 2.8, color = c_navy) +
    geom_text(aes(label = case), hjust = -0.15, vjust = 0.4, size = 2.5, color = "#222222") +
    scale_x_continuous(limits = c(-6.5, 10.5), breaks = seq(-6, 10, by = 4)) +
    scale_y_continuous(limits = c(-65, 105), breaks = seq(-60, 100, by = 40)) +
    labs(
      title = "F. Clinical HPV probe selection audit",
      x = expression("Thermodynamic margin, " * Delta * Delta * G[37] ~ "(kcal/mol)"),
      y = expression("Topology margin, " * Delta * S[ECI])
    ) +
    theme_prism(base_size = 10) +
    theme(aspect.ratio = 0.85)

  # ===========================================================================
  # Assembly & Export (3 x 2 Grid Layout)
  # ===========================================================================
  fig2_6panel <- (p2a + p2b) / (p2c + p2d) / (p2e + p2f) +
    plot_annotation(
      title = "Figure 2. Empirical mismatch topology proof, physical mechanism, and benchmark audits.",
      theme = theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5))
    )

  return(fig2_6panel)
}

###############################################################################
# FIGURE 3: 4-Panel Allele-Specific Benchmark (Portrait)
###############################################################################

# ---- Color Palettes ----

# Probe Types (Scatter Panels A & B)
probe_type_cols <- c(
  "Mismatched" = "#708090",   # Slate Grey
  "MT_PM"      = "#008080",   # Deep Teal
  "WT_PM"      = "#D97706"    # Burnt Amber / Ochre
)

probe_type_shapes <- c(
  "Mismatched" = 16,  # Solid circle
  "MT_PM"      = 17,  # Solid triangle
  "WT_PM"      = 15   # Solid square
)

probe_type_labels <- c(
  "Mismatched" = "Mismatched probes",
  "MT_PM"      = "MT perfect match",
  "WT_PM"      = "WT perfect match"
)

# Stringency Sweep Conditions (Panel C)
condition_cols <- c(
  "KRAS 600 mM"  = "#E11D48",   # Rose Red
  "KRAS 1000 mM" = "#2563EB",   # Vibrant Blue
  "BRAF 1000 mM" = "#059669"    # Emerald Green
)

# Models (Panel D)
model_cols <- c(
  "S_ECI"     = "#18181B",   # Dark Charcoal / Black
  "scalar_NN" = "#F97316",   # Bright Orange
  "Zbox_NN"   = "#94A3B8"    # Light Slate / Silver
)

model_labels <- c(
  "S_ECI"     = expression(S[ECI]),
  "scalar_NN" = expression(Delta*G[NN]),
  "Zbox_NN"   = expression(Z[box])
)

# ---- Base Theme ----
theme_fig3_clean <- theme_prism(
  base_size = 9,
  base_family = "Liberation Sans"
) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    axis.text = element_text(color = "black", size = 8, face = "bold"),
    axis.title = element_text(color = "black", size = 9, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    legend.key = element_blank(),
    plot.tag = element_text(face = "bold", size = 12)
  )

# ---- Helper: Scatter Plot Sub-Panel ----
make_scatter <- function(df, x_col, x_label, show_y_axis = TRUE, show_x_title = TRUE) {
  fit <- lm(as.formula(paste("log10_measured_ratio ~", x_col)), data = df)
  r2  <- summary(fit)$r.squared

  x_range <- range(df[[x_col]], na.rm = TRUE)
  y_range <- range(df$log10_measured_ratio, na.rm = TRUE)
  r2_x <- x_range[1] + 0.05 * diff(x_range)
  r2_y <- y_range[2] - 0.08 * diff(y_range)

  p <- ggplot(df, aes(x = .data[[x_col]], y = log10_measured_ratio)) +
    geom_smooth(method = "lm", se = FALSE, colour = "#334155", linewidth = 0.6) +
    geom_point(aes(colour = probe_type, shape = probe_type), size = 1.8, alpha = 0.85) +
    scale_colour_manual(values = probe_type_cols, labels = probe_type_labels) +
    scale_shape_manual(values = probe_type_shapes, labels = probe_type_labels) +
    annotate("text", x = r2_x, y = r2_y,
             label = sprintf("R^2 == %.2f", r2),
             parse = TRUE, size = 3.2, hjust = 0, fontface = "bold") +
    labs(x = x_label) +
    theme_fig3_clean +
    theme(aspect.ratio = 1, legend.position = "none")

  if (show_y_axis) {
    p <- p + labs(y = "Relative signal (MT / WT)")
  } else {
    p <- p + theme(axis.title.y = element_blank())
  }

  if (!show_x_title) {
    p <- p + theme(axis.title.x = element_blank())
  }

  return(p)
}

# ---- Assembly Function ----
make_figure3 <- function() {
  s16 <- read_csv("supplementary_tables/Supplementary table_S16_Fig3_per_probe_scatter_data.csv",
                  show_col_types = FALSE) %>%
    mutate(probe_type = factor(probe_type, levels = c("Mismatched", "MT_PM", "WT_PM")))
  s13 <- read_csv("supplementary_tables/Supplementary table_S13_stringency_sweep.csv",
                  show_col_types = FALSE)
  s14 <- read_csv("supplementary_tables/Supplementary table_S14_top_quantile_with_bootstrap_ci.csv",
                  show_col_types = FALSE)

  # Panels A & B: Scatter
  kras <- s16 %>% filter(target == "KRAS")
  braf <- s16 %>% filter(target == "BRAF")

  pA1 <- make_scatter(kras, "log10_S_ECI_ratio", expression(S[ECI]~ratio), show_y_axis = TRUE, show_x_title = FALSE) + labs(tag = "A")
  pA2 <- make_scatter(kras, "log10_NN_ratio", expression(Delta*G[NN]~ratio), show_y_axis = FALSE, show_x_title = FALSE)
  pA3 <- make_scatter(kras, "log10_Zbox_NN_ratio", expression(Z[box]~ratio), show_y_axis = FALSE, show_x_title = FALSE)

  pB1 <- make_scatter(braf, "log10_S_ECI_ratio", expression(S[ECI]~ratio), show_y_axis = TRUE, show_x_title = TRUE) + labs(tag = "B")
  pB2 <- make_scatter(braf, "log10_NN_ratio", expression(Delta*G[NN]~ratio), show_y_axis = FALSE, show_x_title = TRUE)
  pB3 <- make_scatter(braf, "log10_Zbox_NN_ratio", expression(Z[box]~ratio), show_y_axis = FALSE, show_x_title = TRUE)

  # Panel C: Stringency Sweep
  sweep <- s13 %>%
    mutate(condition_short = case_when(
      target == "KRAS" & salt == 600  ~ "KRAS 600 mM",
      target == "KRAS" & salt == 1000 ~ "KRAS 1000 mM",
      target == "BRAF" & salt == 1000 ~ "BRAF 1000 mM",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(condition_short)) %>%
    mutate(condition_short = factor(condition_short, levels = c("KRAS 600 mM", "KRAS 1000 mM", "BRAF 1000 mM")))

  pC <- ggplot(sweep, aes(x = temperature, y = partial_r_Zbox_beyond_NN, colour = condition_short)) +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey40", linewidth = 0.5) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.0) +
    scale_colour_manual(values = condition_cols) +
    labs(
      tag = "C",
      x = expression(Temperature~(degree*C)),
      y = expression(partial~r~(Z[box]~beyond~scalar~NN))
    ) +
    theme_fig3_clean +
    theme(
      aspect.ratio = 1,
      legend.position = c(0.65, 0.25),
      legend.background = element_rect(fill = alpha("white", 0.8), color = NA)
    )

  # Panel D: Top-Quantile
  topq <- s14 %>%
    filter(quantile %in% c("top_5", "top_10", "top_25")) %>%
    mutate(
      model = factor(model, levels = c("S_ECI", "scalar_NN", "Zbox_NN")),
      quantile = factor(quantile, levels = c("top_5", "top_10", "top_25"),
                        labels = c("Top 5%", "Top 10%", "Top 25%")),
      target = factor(target, levels = c("KRAS", "BRAF"))
    )

  pD <- ggplot(topq, aes(x = quantile, y = identification_rate_pct, fill = model)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6, color = "black", linewidth = 0.4) +
    facet_wrap(~ target) +
    scale_fill_manual(values = model_cols, labels = model_labels) +
    labs(
      tag = "D",
      x = "Quantile threshold",
      y = "Identification rate (%)"
    ) +
    theme_fig3_clean +
    theme(
      aspect.ratio = 1,
      legend.position = "bottom",
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 9)
    )

  # Assembly
  fig3 <- ((pA1 | pA2 | pA3) /
             (pB1 | pB2 | pB3) /
             (pC  | pD)) +
    plot_layout(heights = c(1, 1, 1.1))

  return(fig3)
}

###############################################################################
# Execute: Generate Both Figures
###############################################################################

fig2 <- make_figure2_6panel()
save_fig(fig2, "Fig. 2", width = 241, height = 318, units = "mm")

fig3 <- make_figure3()
save_fig(fig3, "Fig. 3", width = 180, height = 230, units = "mm")
