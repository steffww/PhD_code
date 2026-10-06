# =============================================================================
# enbi_analysis.R
# Ecological Network Balance Index (ENBI) using the current mt_all_449 cohort.
#
# Adapted from:
#   calculate_enbi_ibd_dep_purePHQ9.R
#   calculate_enbi_ibd_anx_pureGAD7.R
#   ENBI_*_sciencestyle.R
#
# Design (symptom cutoff >=5):
#   Healthy only, symptom only, IBD only, and IBD + symptom.
#
# Default computation:
#   B = 500 overlapping 80% subsamples per group; FlashWeave sensitive=true,
#   heterogeneous=false; ENBI = rho - median(rho in the healthy-only group).
#
# Important interpretation boundary:
#   bootstrap/subsample rho values are repeated, overlapping resamples from the
#   same participants. Their Mann-Whitney p values are descriptive and do not
#   represent inference from independent cohorts or independent participants.
#
# Run:
#   Rscript enbi_analysis.R
# Smoke test:
#   ENBI_B=5 ENBI_FORCE=1 Rscript enbi_analysis.R
# =============================================================================

suppressPackageStartupMessages({
  library(microeco)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

source("R/00_setup.R")

B <- as.integer(Sys.getenv("ENBI_B", "500"))
FRACTION <- as.numeric(Sys.getenv("ENBI_FRACTION", "0.8"))
SEED <- as.integer(Sys.getenv("ENBI_SEED", "1"))
FORCE <- identical(Sys.getenv("ENBI_FORCE", "0"), "1")
CORES <- max(1L, as.integer(Sys.getenv("ENBI_CORES", "4")))

JULIA <- Sys.getenv("ENBI_JULIA", unname(Sys.which("julia")))
JULIA_WORKER <- "scripts/run_flashweave_bootstrap.jl"
INPUT_RDS <- "microbiome/mt_all_449.rds"

OUT_DIR <- file.path("results", "enbi_analysis")
WORK_DIR <- file.path(OUT_DIR, "intermediate")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(WORK_DIR, recursive = TRUE, showWarnings = FALSE)

stopifnot(
  B > 0L,
  FRACTION > 0,
  FRACTION <= 1,
  nzchar(JULIA),
  file.exists(JULIA),
  file.exists(JULIA_WORKER),
  file.exists(INPUT_RDS)
)

save_figure <- function(plot, stem, width, height, dpi = 300) {
  ggsave(
    file.path(OUT_DIR, paste0(stem, ".pdf")),
    plot, width = width, height = height, device = cairo_pdf, bg = "white"
  )
  ggsave(
    file.path(OUT_DIR, paste0(stem, ".png")),
    plot, width = width, height = height, dpi = dpi, bg = "white"
  )
  message("Saved: ", file.path(OUT_DIR, stem), " (.pdf + .png)")
}

cliffs_delta <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (!length(x) || !length(y)) return(NA_real_)
  (sum(outer(x, y, ">")) - sum(outer(x, y, "<"))) /
    (length(x) * length(y))
}

compare_distributions <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  tibble::tibble(
    n_boot_x = length(x),
    n_boot_y = length(y),
    p_value = wilcox.test(x, y, exact = FALSE)$p.value,
    cliffs_delta = cliffs_delta(x, y)
  )
}

format_p <- function(p) {
  ifelse(p < 0.001, formatC(p, format = "e", digits = 1),
         formatC(p, format = "f", digits = 3))
}

mt <- readRDS(INPUT_RDS)
species <- as.matrix(mt$otu_table)
metadata <- mt$sample_table

stopifnot(
  ncol(species) == 449L,
  nrow(metadata) == 449L,
  identical(colnames(species), rownames(metadata)),
  all(colSums(species) > 0)
)

# The current microtable contains retained species-level relative abundance in
# percentage-like units. Renormalize every sample to a sum of 1, matching the
# referenced publication pipeline.
species <- sweep(species, 2, colSums(species), "/")

input_manifest <- tibble::tibble(
  input = INPUT_RDS,
  md5 = unname(tools::md5sum(INPUT_RDS)),
  n_species = nrow(species),
  n_samples = ncol(species),
  B = B,
  fraction = FRACTION,
  seed = SEED,
  cores = CORES,
  julia = JULIA,
  julia_worker = JULIA_WORKER
)
data.table::fwrite(
  input_manifest,
  file.path(OUT_DIR, "source_manifest.csv")
)

make_outcome_config <- function(
    id, score_var, score_label, symptom_name, symptom_abbrev,
    healthy_group, symptom_only_group, ibd_only_group, combined_group) {
  list(
    id = id,
    score_var = score_var,
    score_label = score_label,
    symptom_name = symptom_name,
    symptom_abbrev = symptom_abbrev,
    healthy_group = healthy_group,
    symptom_only_group = symptom_only_group,
    ibd_only_group = ibd_only_group,
    combined_group = combined_group,
    severity_groups = c(symptom_only_group, ibd_only_group, combined_group)
  )
}

configs <- list(
  make_outcome_config(
    id = "phq9",
    score_var = "PHQ9_sum",
    score_label = "PHQ-9",
    symptom_name = "depression",
    symptom_abbrev = "Dep.",
    healthy_group = "HC_pure",
    symptom_only_group = "dep_no_IBD",
    ibd_only_group = "IBD_ND",
    combined_group = "IBD_D"
  ),
  make_outcome_config(
    id = "gad7",
    score_var = "GAD7_sum",
    score_label = "GAD-7",
    symptom_name = "anxiety",
    symptom_abbrev = "Anx.",
    healthy_group = "HC_pure_anx",
    symptom_only_group = "anx_no_IBD",
    ibd_only_group = "IBD_NA",
    combined_group = "IBD_A"
  )
)

build_group_membership <- function(config) {
  score <- metadata[[config$score_var]]
  is_ibd <- !is.na(metadata$IBD) & metadata$IBD == "IBD"
  has_score <- !is.na(score)

  groups <- list(
    !is_ibd & score < 5 & has_score,
    !is_ibd & score >= 5 & has_score,
    is_ibd & score < 5 & has_score,
    is_ibd & score >= 5 & has_score
  )
  names(groups) <- c(
    config$healthy_group,
    config$symptom_only_group,
    config$ibd_only_group,
    config$combined_group
  )
  groups
}

run_flashweave_group <- function(config, group_name, membership, group_index) {
  sample_ids <- colnames(species)[membership]
  input_path <- file.path(
    WORK_DIR,
    paste0(config$id, "_", group_name, "_input.tsv")
  )
  output_path <- file.path(
    WORK_DIR,
    paste0(config$id, "_", group_name, "_rho_B", B, ".tsv")
  )

  table_for_julia <- t(species[, sample_ids, drop = FALSE])
  write.table(
    table_for_julia,
    input_path,
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )

  valid_cached_output <- FALSE
  if (file.exists(output_path) && !FORCE) {
    cached <- tryCatch(read.delim(output_path), error = function(e) NULL)
    valid_cached_output <- !is.null(cached) && nrow(cached) == B &&
      all(c("bootstrap", "rho", "n_edges", "n_samples") %in% names(cached))
  }

  if (!valid_cached_output) {
    message(
      "FlashWeave: ", config$id, " / ", group_name,
      " (participants=", length(sample_ids), ", B=", B, ")"
    )
    status <- system2(
      JULIA,
      args = c(
        JULIA_WORKER,
        input_path,
        output_path,
        B,
        FRACTION,
        SEED + 100L * match(config$id, c("phq9", "gad7")) + group_index
      )
    )
    if (status != 0L || !file.exists(output_path)) {
      stop("FlashWeave failed for ", config$id, " / ", group_name)
    }
  } else {
    message("Reusing cached result: ", output_path)
  }

  read.delim(output_path) %>%
    mutate(
      outcome = config$id,
      group = group_name,
      participant_n = length(sample_ids)
    )
}

build_plot_and_statistics <- function(config, rho_df, group_sizes) {
  baseline_median <- median(
    rho_df$rho[rho_df$group == config$healthy_group],
    na.rm = TRUE
  )
  plot_df <- rho_df %>%
    mutate(enbi = rho - baseline_median)

  healthy <- plot_df$enbi[plot_df$group == config$healthy_group]
  any_symptom_or_ibd <- plot_df$enbi[
    plot_df$group %in% config$severity_groups
  ]

  group_1 <- config$severity_groups[[1]]
  group_2 <- config$severity_groups[[2]]
  group_3 <- config$severity_groups[[3]]

  comparisons <- bind_rows(
    compare_distributions(any_symptom_or_ibd, healthy) %>%
      mutate(
        comparison = "Any vs Healthy",
        group_x = "Any",
        group_y = "Healthy",
        participant_n_x = sum(unlist(group_sizes[config$severity_groups])),
        participant_n_y = group_sizes[[config$healthy_group]]
      ),
    compare_distributions(
      plot_df$enbi[plot_df$group == group_2],
      plot_df$enbi[plot_df$group == group_1]
    ) %>%
      mutate(
        comparison = paste(group_2, "vs", group_1),
        group_x = group_2,
        group_y = group_1,
        participant_n_x = group_sizes[[group_2]],
        participant_n_y = group_sizes[[group_1]]
      ),
    compare_distributions(
      plot_df$enbi[plot_df$group == group_3],
      plot_df$enbi[plot_df$group == group_2]
    ) %>%
      mutate(
        comparison = paste(group_3, "vs", group_2),
        group_x = group_3,
        group_y = group_2,
        participant_n_x = group_sizes[[group_3]],
        participant_n_y = group_sizes[[group_2]]
      )
  ) %>%
    mutate(
      outcome = config$id,
      q_value = p.adjust(p_value, method = "BH"),
      .before = 1
    )

  left_df <- bind_rows(
    tibble::tibble(group = "Healthy", enbi = healthy),
    tibble::tibble(group = "Any", enbi = any_symptom_or_ibd)
  ) %>%
    mutate(group = factor(group, levels = c("Healthy", "Any")))

  right_df <- plot_df %>%
    filter(group %in% config$severity_groups) %>%
    mutate(group = factor(group, levels = config$severity_groups))

  finite_enbi <- plot_df$enbi[is.finite(plot_df$enbi)]
  y_min <- min(finite_enbi)
  y_max <- max(finite_enbi)
  y_range <- y_max - y_min
  if (y_range == 0) y_range <- 1

  healthy_n <- group_sizes[[config$healthy_group]]
  any_n <- sum(unlist(group_sizes[config$severity_groups]))
  left_labels <- c(
    Healthy = paste0("Healthy\n(n=", healthy_n, ")"),
    Any = paste0("Any\n(n=", any_n, ")")
  )
  right_labels <- c(
    stats::setNames(
      paste0(config$symptom_abbrev, " alone\n(n=", group_sizes[[group_1]], ")"),
      group_1
    ),
    stats::setNames(
      paste0("IBD alone\n(n=", group_sizes[[group_2]], ")"),
      group_2
    ),
    stats::setNames(
      paste0("IBD + ", config$symptom_abbrev, "\n(n=", group_sizes[[group_3]], ")"),
      group_3
    )
  )

  left_stat <- comparisons %>% filter(comparison == "Any vs Healthy")
  adjacent_stats <- comparisons %>% filter(comparison != "Any vs Healthy")

  significance_data <- tibble::tibble(
    xmin = c(1, 2),
    xmax = c(2, 3),
    y_position = y_max + c(0.14, 0.26) * y_range,
    annotations = paste0(
      "p = ", format_p(adjacent_stats$p_value),
      "\nCliff's delta = ", sprintf("%.2f", adjacent_stats$cliffs_delta)
    )
  )

  theme_science <- theme_classic(base_size = 13) +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      axis.line = element_line(colour = "black", linewidth = 0.8),
      axis.ticks = element_line(colour = "black", linewidth = 0.8),
      axis.text = element_text(colour = "black", size = 10),
      axis.title.y = element_text(colour = "black", size = 13, face = "bold"),
      axis.title.x = element_blank(),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.position = "none",
      panel.grid = element_blank()
    )

  colour_healthy <- "#1B6E8C"
  colour_any <- "#D6423F"
  severity_colours <- grDevices::colorRampPalette(c("#F6B8B4", "#7A1620"))(3)
  names(severity_colours) <- config$severity_groups

  p_left <- ggplot(left_df, aes(group, enbi, fill = group)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_boxplot(
      width = 0.52, linewidth = 0.8, colour = "black", outlier.shape = NA
    ) +
    annotate(
      "text",
      x = 0.55,
      y = y_max + 0.04 * y_range,
      label = paste0(
        "p = ", format_p(left_stat$p_value),
        "\nCliff's delta = ", sprintf("%.2f", left_stat$cliffs_delta)
      ),
      hjust = 0,
      vjust = 1,
      size = 3.4
    ) +
    scale_fill_manual(values = c(Healthy = colour_healthy, Any = colour_any)) +
    scale_x_discrete(labels = left_labels) +
    coord_cartesian(
      ylim = c(y_min - 0.05 * y_range, y_max + 0.14 * y_range),
      clip = "off"
    ) +
    labs(y = "ENBI") +
    theme_science

  p_right <- ggplot(right_df, aes(group, enbi, fill = group)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_boxplot(
      width = 0.52, linewidth = 0.8, colour = "black", outlier.shape = NA
    ) +
    geom_segment(
      data = significance_data,
      inherit.aes = FALSE,
      aes(x = xmin, xend = xmax, y = y_position, yend = y_position),
      linewidth = 0.55,
      colour = "black"
    ) +
    geom_segment(
      data = significance_data,
      inherit.aes = FALSE,
      aes(
        x = xmin, xend = xmin,
        y = y_position, yend = y_position - 0.025 * y_range
      ),
      linewidth = 0.55,
      colour = "black"
    ) +
    geom_segment(
      data = significance_data,
      inherit.aes = FALSE,
      aes(
        x = xmax, xend = xmax,
        y = y_position, yend = y_position - 0.025 * y_range
      ),
      linewidth = 0.55,
      colour = "black"
    ) +
    geom_text(
      data = significance_data,
      inherit.aes = FALSE,
      aes(
        x = (xmin + xmax) / 2,
        y = y_position + 0.01 * y_range,
        label = annotations
      ),
      size = 3.1,
      vjust = 0
    ) +
    scale_fill_manual(values = severity_colours) +
    scale_x_discrete(labels = right_labels) +
    coord_cartesian(
      ylim = c(y_min - 0.05 * y_range, y_max + 0.38 * y_range),
      clip = "off"
    ) +
    labs(y = NULL) +
    theme_science +
    theme(
      axis.line.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )

  panel <- (p_left | p_right) +
    plot_layout(widths = c(2, 3)) +
    plot_annotation(
      title = paste0(
        "Inferred network balance: IBD x ", config$symptom_name,
        " (", config$score_label, ">=5)"
      ),
      subtitle = paste0(
        "B=", B, ", 80% overlapping subsamples; ",
        "p values describe bootstrap distributions"
      ),
      theme = theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 9, hjust = 0.5, colour = "grey35")
      )
    )

  rho_summary <- plot_df %>%
    group_by(outcome, group, participant_n) %>%
    summarise(
      n_boot = sum(is.finite(rho)),
      median_rho = median(rho, na.rm = TRUE),
      median_enbi = median(enbi, na.rm = TRUE),
      mean_edges = mean(n_edges, na.rm = TRUE),
      .groups = "drop"
    )

  list(
    plot = panel,
    plot_data = plot_df,
    comparisons = comparisons,
    rho_summary = rho_summary,
    baseline_median = baseline_median
  )
}

analysis_results <- vector("list", length(configs))
names(analysis_results) <- vapply(configs, `[[`, character(1), "id")

for (config in configs) {
  membership <- build_group_membership(config)
  group_sizes <- lapply(membership, sum)

  message("\n", toupper(config$id), " group sizes:")
  for (group_name in names(group_sizes)) {
    message("  ", group_name, ": ", group_sizes[[group_name]])
  }

  group_names <- names(membership)
  rho_list <- parallel::mclapply(
    X = seq_along(group_names),
    FUN = function(group_index) {
      group_name <- group_names[[group_index]]
      run_flashweave_group(
        config = config,
        group_name = group_name,
        membership = membership[[group_name]],
        group_index = group_index
      )
    },
    mc.cores = min(CORES, length(group_names)),
    mc.preschedule = FALSE
  )

  rho_df <- bind_rows(rho_list)
  if (any(!is.finite(rho_df$rho))) {
    warning(config$id, ": non-finite rho values detected; excluded from tests")
  }

  result <- build_plot_and_statistics(
    config = config,
    rho_df = rho_df,
    group_sizes = group_sizes
  )
  analysis_results[[config$id]] <- result

  data.table::fwrite(
    result$plot_data,
    file.path(OUT_DIR, paste0("enbi_bootstrap_", config$id, ".csv")),
    na = "NA"
  )
  data.table::fwrite(
    result$rho_summary,
    file.path(OUT_DIR, paste0("enbi_group_summary_", config$id, ".csv")),
    na = "NA"
  )
  data.table::fwrite(
    result$comparisons,
    file.path(OUT_DIR, paste0("enbi_comparisons_", config$id, ".csv")),
    na = "NA"
  )

  save_figure(
    result$plot,
    paste0("enbi_panel_f_", config$id),
    width = 9,
    height = 5.5
  )
}

p_combined <- patchwork::wrap_elements(full = analysis_results$phq9$plot) /
  patchwork::wrap_elements(full = analysis_results$gad7$plot) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "A")

save_figure(
  p_combined,
  "enbi_panel_f_combined",
  width = 9.5,
  height = 11
)

all_comparisons <- bind_rows(
  analysis_results$phq9$comparisons,
  analysis_results$gad7$comparisons
) %>%
  mutate(q_value_across_both_outcomes = p.adjust(p_value, method = "BH"))

data.table::fwrite(
  all_comparisons,
  file.path(OUT_DIR, "enbi_comparisons_all.csv"),
  na = "NA"
)

cat("\nENBI analysis complete.\n")
cat("Input: ", INPUT_RDS, "\n", sep = "")
cat("Output: ", normalizePath(OUT_DIR), "\n", sep = "")
