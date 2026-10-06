library(data.table)
library(dplyr)
library(microeco)

source("R/00_setup.R")

if (!requireNamespace("ppcor", quietly = TRUE)) {
  stop("Package 'ppcor' is required. Install it with install.packages('ppcor').")
}

# =============================================================================
# 1. Cohort and prevalence
# =============================================================================

mt_all <- readRDS("microbiome/mt_all_449.rds")

# baseline_full contains the hematology-derived inflammatory indices that are
# not present in baseline_core. Restrict it to participants retained in mt_all
# so that this script continues to use the same 449-person analysis cohort.
meta_all <- fread("metadata/baseline_full_450.csv") %>%
  filter(participant_id %in% mt_all$sample_table$participant_id)

stopifnot(
  nrow(meta_all) == nrow(mt_all$sample_table),
  !anyDuplicated(meta_all$participant_id)
)

cohort_ibd = filter(meta_all, IBD == 'IBD')
mt_ibd = readRDS('./microbiome/mt_ibd_283.rds')
mt_ibd$sample_table = filter(mt_ibd$sample_table, participant_id %in% cohort_ibd$participant_id)
mt_ibd$tidy_dataset()
mt_ibd


saveRDS(meta_all,file = './metadata/cohort_all_449.rds')
saveRDS(cohort_ibd,file = './metadata/cohort_ibd_282.rds')
saveRDS(mt_ibd,file = './microbiome/mt_ibd_282.rds')

cat(sprintf("PHQ-9 >= 10: %d/%d (%.1f%%)\n",
            sum(meta_all$PHQ9_sum >= 10, na.rm = TRUE),
            sum(!is.na(meta_all$PHQ9_sum)),
            100 * mean(meta_all$PHQ9_sum >= 10, na.rm = TRUE)))

cat(sprintf("GAD-7 >= 10: %d/%d (%.1f%%)\n",
            sum(meta_all$GAD7_sum >= 10, na.rm = TRUE),
            sum(!is.na(meta_all$GAD7_sum)),
            100 * mean(meta_all$GAD7_sum >= 10, na.rm = TRUE)))


# =============================================================================
# 2. IBDQ domains, inflammation and disease severity vs emotional symptoms
#    Pairwise partial Spearman correlation, adjusted for
#    Age, Sex, BMI, site and hPDI
# =============================================================================

# Disease-specific measures are meaningful only among participants with IBD.
meta_ibd <- meta_all %>%
  filter(IBD == "IBD") %>%
  mutate(
    Sex_num = as.integer(factor(Sex)),
    site_num = as.integer(factor(site)),
    across(all_of(c("Age", "BMI", "hPDI")), as.numeric)
  )

# Per request, IBDQ total is recalculated as the sum of the four domain scores.
# rowSums(..., na.rm = FALSE) deliberately returns NA unless all four domains
# are available; no domain-level imputation is performed.
ibdq_domain_vars <- c(
  "IBDQ_bowel", "IBDQ_systemic", "IBDQ_emotional", "IBDQ_social"
)
meta_ibd$IBDQ_total_domains <- rowSums(
  as.data.frame(dplyr::select(meta_ibd, all_of(ibdq_domain_vars))),
  na.rm = FALSE
)

# QA: where the item-recomputed total is available, it must equal the sum of
# the four domains. The original locked IBDQ_total is retained in the metadata
# but is intentionally not used below.
ibdq_check <- !is.na(meta_ibd$IBDQ_total_domains) &
  !is.na(meta_ibd$IBDQ_item_total_recomputed)
stopifnot(all(
  meta_ibd$IBDQ_total_domains[ibdq_check] ==
    meta_ibd$IBDQ_item_total_recomputed[ibdq_check]
))

score_specs <- tibble::tribble(
  ~score,      ~score_label,
  "PHQ9_sum", "PHQ-9",
  "GAD7_sum", "GAD-7"
)

# population controls the clinically applicable subset for each metric:
# CDAI and SES-CD are evaluated in CD only; Mayo indices in UC only.
metric_specs <- tibble::tribble(
  ~metric,                ~metric_label,             ~metric_group,     ~population,
  "IBDQ_total_domains",   "IBDQ total (domain sum)", "IBDQ",            "All IBD",
  "IBDQ_bowel",           "IBDQ bowel",              "IBDQ",            "All IBD",
  "IBDQ_systemic",        "IBDQ systemic",           "IBDQ",            "All IBD",
  "IBDQ_emotional",       "IBDQ emotional",          "IBDQ",            "All IBD",
  "IBDQ_social",          "IBDQ social",             "IBDQ",            "All IBD",
  "hsCRP",                "hs-CRP",                  "Inflammation",     "All IBD",
  "Calprotectin",         "Fecal calprotectin",      "Inflammation",     "All IBD",
  "WBC",                  "White blood cells",       "Inflammation",     "All IBD",
  "Fibrinogen",           "Fibrinogen",              "Inflammation",     "All IBD",
  "SII",                  "SII",                     "Inflammation",     "All IBD",
  "NLR",                  "NLR",                     "Inflammation",     "All IBD",
  "PLR",                  "PLR",                     "Inflammation",     "All IBD",
  "LMR",                  "LMR",                     "Inflammation",     "All IBD",
  "CDAI_recalc",          "CDAI",                    "Disease severity", "CD",
  "SES_CD",               "SES-CD",                  "Disease severity", "CD",
  "Mayo_total_recalc",    "Mayo total",              "Disease severity", "UC",
  "Mayo_partial_recalc",  "Partial Mayo",            "Disease severity", "UC"
)

covariates <- c("Age", "Sex_num", "BMI", "site_num", "hPDI")

partial_spearman <- function(data, score, metric, population,
                             covariates, min_n = 10L) {
  analysis_data <- data
  if (population == "CD") {
    analysis_data <- analysis_data %>% filter(IBD_subtype == "CD")
  } else if (population == "UC") {
    analysis_data <- analysis_data %>% filter(IBD_subtype == "UC")
  }

  required_vars <- c(score, metric, covariates)
  analysis_data <- analysis_data %>%
    dplyr::select(all_of(required_vars)) %>%
    filter(if_all(everything(), ~ !is.na(.x))) %>%
    as.data.frame()

  n_complete <- nrow(analysis_data)

  # A constant exposure/outcome cannot be correlated. Constant covariates are
  # dropped only within the relevant pairwise-complete subset and are recorded.
  if (n_complete == 0L ||
      dplyr::n_distinct(analysis_data[[score]]) < 2L ||
      dplyr::n_distinct(analysis_data[[metric]]) < 2L) {
    return(tibble::tibble(
      rho = NA_real_, p_value = NA_real_, n = n_complete,
      covariates_used = NA_character_, covariates_dropped = NA_character_,
      status = "not_estimable"
    ))
  }

  varying_covariates <- covariates[
    vapply(analysis_data[covariates], dplyr::n_distinct, integer(1)) > 1L
  ]
  dropped_covariates <- setdiff(covariates, varying_covariates)

  # Five adjustment variables require at least eight observations
  # mathematically; min_n = 10 adds a small stability floor. Low-n rows remain
  # in the exported table and are shown as non-estimable cells in the heatmap.
  if (n_complete < min_n ||
      n_complete <= length(varying_covariates) + 2L) {
    return(tibble::tibble(
      rho = NA_real_, p_value = NA_real_, n = n_complete,
      covariates_used = paste(varying_covariates, collapse = ";"),
      covariates_dropped = paste(dropped_covariates, collapse = ";"),
      status = "insufficient_n"
    ))
  }

  test <- tryCatch(
    ppcor::pcor.test(
      x = analysis_data[[score]],
      y = analysis_data[[metric]],
      z = analysis_data[, varying_covariates, drop = FALSE],
      method = "spearman"
    ),
    error = function(e) e
  )

  if (inherits(test, "error")) {
    return(tibble::tibble(
      rho = NA_real_, p_value = NA_real_, n = n_complete,
      covariates_used = paste(varying_covariates, collapse = ";"),
      covariates_dropped = paste(dropped_covariates, collapse = ";"),
      status = paste0("error: ", conditionMessage(test))
    ))
  }

  tibble::tibble(
    rho = unname(test$estimate),
    p_value = test$p.value,
    n = n_complete,
    covariates_used = paste(varying_covariates, collapse = ";"),
    covariates_dropped = paste(dropped_covariates, collapse = ";"),
    status = "ok"
  )
}

partial_cor_results <- tidyr::crossing(
  metric_specs,
  score_specs
) %>%
  purrr::pmap_dfr(function(metric, metric_label, metric_group, population,
                           score, score_label) {
    partial_spearman(
      data = meta_ibd,
      score = score,
      metric = metric,
      population = population,
      covariates = covariates
    ) %>%
      mutate(
        metric = metric,
        metric_label = metric_label,
        metric_group = metric_group,
        population = population,
        score = score,
        score_label = score_label,
        .before = 1
      )
  })

# BH correction is applied once across all estimable cells displayed in the
# combined heatmap. Both nominal p and BH-FDR q are retained in the CSV.
partial_cor_results$q_value <- p.adjust(
  partial_cor_results$p_value,
  method = "BH"
)

partial_cor_results <- partial_cor_results %>%
  mutate(
    fdr_symbol = case_when(
      q_value < 0.001 ~ "***",
      q_value < 0.01  ~ "**",
      q_value < 0.05  ~ "*",
      TRUE ~ ""
    ),
    cell_label = if_else(
      status == "ok",
      sprintf("%.2f%s\nn=%d", rho, fdr_symbol, n),
      sprintf("NE\nn=%d", n)
    )
  ) %>%
  arrange(match(metric, metric_specs$metric), match(score, score_specs$score))

results_dir <- file.path("results", "cohort_analysis")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(
  partial_cor_results %>% dplyr::select(-cell_label, -fdr_symbol),
  file.path(results_dir, "partial_spearman_results.csv"),
  na = "NA"
)

group_palettes <- list(
  # Match the reference figure's group-specific green/orange/purple families.
  # Within every group, |rho| = 0 is white and |rho| = 1 is the dark endpoint.
  "IBDQ" = "#16876C",
  "Inflammation" = "#D95F0E",
  "Disease severity" = "#5E3C99"
)

make_partial_cor_panel <- function(result_data, group_name, show_y = FALSE) {
  metric_levels <- metric_specs %>%
    filter(metric_group == group_name,
           metric %in% unique(result_data$metric)) %>%
    pull(metric_label)

  plot_data <- result_data %>%
    filter(metric_group == group_name) %>%
    mutate(
      metric_label = factor(metric_label, levels = metric_levels),
      # Reverse the factor so PHQ-9 is displayed above GAD-7.
      score_label = factor(score_label, levels = rev(score_specs$score_label)),
      abs_rho = abs(rho),
      # q significance takes precedence. # denotes nominal p significance only
      # when the corresponding BH-FDR q value is not below 0.05.
      significance = case_when(
        !is.na(q_value) & q_value < 0.05 ~ "*",
        !is.na(p_value) & p_value < 0.05 ~ "#",
        TRUE ~ ""
      )
    )

  group_colour <- group_palettes[[group_name]]

  p <- ggplot(plot_data, aes(x = metric_label, y = score_label, fill = abs_rho)) +
    # Dark borders separate every individual cell, including adjacent columns.
    geom_tile(colour = "grey35", linewidth = 0.45) +
    geom_text(
      aes(label = significance),
      colour = "black",
      size = 7,
      fontface = "bold",
      vjust = 0.65
    ) +
    scale_fill_gradient(
      low = "white",
      high = group_colour,
      limits = c(0, 1),
      na.value = "grey90",
      guide = "none"
    ) +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 12)) +
    labs(x = NULL, y = NULL, title = group_name) +
    theme_pub +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "grey25", fill = NA, linewidth = 0.6),
      axis.text.x = element_text(
        angle = 55, hjust = 1, vjust = 1, size = 9, colour = "black"
      ),
      axis.text.y = element_text(size = 11, colour = "black"),
      axis.ticks = element_blank(),
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      legend.position = "none",
      plot.margin = margin(4, 4, 4, 4)
    )

  if (!show_y) {
    p <- p + theme(axis.text.y = element_blank())
  }

  p
}

# The tiles use category-specific hues, while this single neutral legend shows
# their shared intensity scale. Its vertical bar runs from dark (|rho| = 1) at
# the top to light (|rho| = 0) at the bottom, as in the reference figure.
make_shared_intensity_legend <- function() {
  legend_data <- tibble::tibble(
    x = 1,
    y = seq(0, 1, length.out = 100),
    intensity = seq(0, 1, length.out = 100)
  )

  ggplot(legend_data, aes(x = x, y = y, fill = intensity)) +
    geom_tile(show.legend = TRUE) +
    scale_fill_gradient(
      low = "white",
      high = "black",
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.50, 0.75, 1),
      name = "Association strength\n|partial Spearman rho|"
    ) +
    guides(
      fill = guide_colorbar(
        direction = "vertical",
        title.position = "top",
        title.hjust = 0.5,
        barwidth = grid::unit(0.38, "cm"),
        barheight = grid::unit(1.6, "cm")
      )
    ) +
    coord_cartesian(xlim = c(2, 3)) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      plot.margin = margin(0, 0, 0, 0)
    )
}

make_partial_cor_heatmap <- function(result_data) {
  group_order <- c("IBDQ", "Inflammation", "Disease severity")
  groups_present <- group_order[group_order %in% unique(result_data$metric_group)]
  group_widths <- vapply(
    groups_present,
    function(g) dplyr::n_distinct(result_data$metric[result_data$metric_group == g]),
    numeric(1)
  )

  panels <- lapply(seq_along(groups_present), function(i) {
    make_partial_cor_panel(
      result_data = result_data,
      group_name = groups_present[[i]],
      show_y = i == 1L
    )
  })

  panels[[length(panels) + 1L]] <- make_shared_intensity_legend()

  patchwork::wrap_plots(
    panels,
    nrow = 1,
    widths = c(group_widths, 1.8)
  )
}

p_heatmap_all <- make_partial_cor_heatmap(
  partial_cor_results
)
save_plot(
  p_heatmap_all,
  name = "partial_spearman_all_heatmap",
  subdir = "cohort_analysis",
  width = 14,
  height = 4.8
)

p_heatmap_ibdq <- make_partial_cor_heatmap(
  partial_cor_results %>% filter(metric_group == "IBDQ")
)
save_plot(
  p_heatmap_ibdq,
  name = "partial_spearman_ibdq_heatmap",
  subdir = "cohort_analysis",
  width = 8,
  height = 4.2
)

print(partial_cor_results %>%
        dplyr::select(metric_group, metric_label, population, score_label,
                      rho, p_value, q_value, n, status),
      n = Inf)
print(p_heatmap_all)


# =============================================================================
# 3. Group-meeting descriptive figures
#    Integrated with the current eligible mt_all-matched cohort above
# =============================================================================

suppressPackageStartupMessages({
  library(ggforce)
  library(cowplot)
  library(ggpubr)
  library(PResiduals)
  library(Formula)
  library(rms)
})

# Use the device-independent sans family. On macOS this renders as a
# Helvetica/Arial-like font and avoids PostScript font-database warnings when
# cowplot calculates text bounds before the Cairo PDF device is opened.
figure_font <- "sans"

out_dir <- "results/group_meeting/add"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =========================
# 3a. 数据预处理
# =========================
# Reuse the current, already validated inputs. This avoids silently falling
# back to stale/missing intermediate RDS files from another project tree.
samp_tab <- as.data.frame(meta_all)  # 449 participants in mt_all
cohort <- as.data.frame(meta_ibd)    # 282 IBD participants in mt_all

stopifnot(
  nrow(samp_tab) == nrow(mt_all$sample_table),
  all(cohort$IBD == "IBD"),
  !anyDuplicated(samp_tab$participant_id),
  !anyDuplicated(cohort$participant_id)
)

plot_dat <- samp_tab %>%
  mutate(
    IBD_group = case_when(
      IBD == "nonIBD" ~ "Healthy",
      IBD == "IBD" ~ "IBD",
      TRUE ~ NA_character_
    ),
    
    Sex_group = case_when(
      tolower(as.character(Sex)) %in% c("male", "m", "男") ~ "Male",
      tolower(as.character(Sex)) %in% c("female", "f", "女") ~ "Female",
      TRUE ~ as.character(Sex)
    ),
    Sex_group = factor(Sex_group, levels = c("Female", "Male")),
    
    Age_group = case_when(
      Age < 30 ~ "<30",
      Age >= 30 & Age < 45 ~ "30-44",
      Age >= 45 & Age < 60 ~ "45-59",
      Age >= 60 ~ ">=60",
      TRUE ~ NA_character_
    ),
    Age_group = factor(
      Age_group,
      levels = c("<30", "30-44", "45-59", ">=60")
    ),
    BMI_group = case_when(
      BMI < 18.5 ~ "Underweight",
      BMI >= 18.5 & BMI < 24 ~ "Normal",
      BMI >= 24 & BMI < 28 ~ "Overweight",
      BMI >= 28 ~ "Obese",
      TRUE ~ NA_character_
    ),
    BMI_group = factor(
      BMI_group,
      levels = c("Underweight", "Normal", "Overweight", "Obese")
    ),
    Depression_group = case_when(
      PHQ9_sum >= 10 ~ "Yes",
      PHQ9_sum < 10 ~ "No",
      TRUE ~ NA_character_
    ),
    
    Anxiety_group = case_when(
      GAD7_sum >= 10 ~ "Yes",
      GAD7_sum < 10 ~ "No",
      TRUE ~ NA_character_
    ),

    Depression_tri = case_when(
      PHQ9_sum >= 10 ~ "Moderate-Severe",
      PHQ9_sum >= 5  ~ "Mild",
      PHQ9_sum < 5   ~ "Minimal",
      TRUE ~ NA_character_
    ),
    Depression_tri = factor(
      Depression_tri,
      levels = c("Minimal", "Mild", "Moderate-Severe")
    ),

    Anxiety_tri = case_when(
      GAD7_sum >= 10 ~ "Moderate-Severe",
      GAD7_sum >= 5  ~ "Mild",
      GAD7_sum < 5   ~ "Minimal",
      TRUE ~ NA_character_
    ),
    Anxiety_tri = factor(
      Anxiety_tri,
      levels = c("Minimal", "Mild", "Moderate-Severe")
    )
  )

plot_donut_by_ibd <- function(data, var, title = NULL, palette = NULL) {
  
  var <- rlang::ensym(var)
  
  df <- data %>%
    filter(!is.na(IBD_group), !is.na(!!var)) %>%
    count(IBD_group, !!var, name = "n") %>%
    group_by(IBD_group) %>%
    mutate(
      total_n = sum(n),
      ratio = n / total_n,
      ymax = cumsum(ratio),
      ymin = lag(ymax, default = 0),
      start = 2 * pi * ymin,
      end = 2 * pi * ymax
    ) %>%
    ungroup()
  
  center_df <- df %>%
    distinct(IBD_group, total_n) %>%
    mutate(
      label = paste0(IBD_group, "\n(n=", total_n, ")")
    )
  
  p <- ggplot(df) +
    geom_arc_bar(
      aes(
        x0 = 0,
        y0 = 0,
        r0 = 0.45,
        r = 1,
        start = start,
        end = end,
        fill = !!var
      ),
      color = "grey40",
      linewidth = 0.3
    ) +
    geom_text(
      data = center_df,
      aes(x = 0, y = 0, label = label),
      size = 3.2,
      lineheight = 0.9,
      family = figure_font
    ) +
    facet_wrap(~ IBD_group, nrow = 1) +
    coord_fixed() +
    labs(
      title = title,
      fill = NULL
    ) +
    theme_void(base_family = figure_font) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        size = 13,
        face = "bold"
      ),
      strip.text = element_blank(),
      legend.position = "right",
      legend.text = element_text(size = 9),
      legend.key.width = unit(0.5, "cm"),
      legend.key.height = unit(0.28, "cm"),
      legend.spacing.y = unit(0.1, "cm"),
      plot.margin = margin(5, 5, 5, 5)
    )
  
  if (!is.null(palette)) {
    p <- p + scale_fill_manual(values = palette)
  }
  
  return(p)
}

pal_sex <- c(
  "Female" = "#F4C7C3",
  "Male" = "#9ED0F6"
)

pal_age <- c(
  "<30" = "#FFE3B3",
  "30-44" = "#FDBF6F",
  "45-59" = "#E69F00",
  ">=60" = "#BDBDBD"
)

pal_bmi <- c(
  "Underweight" = "#9E9E9E",
  "Normal" = "#7AC6A4",
  "Overweight" = "#B7E1B0",
  "Obese" = "#7AA33D"
)

pal_dep <- c(
  "No" = "#BDBDBD",
  "Yes" = "#F4A3A3"
)

pal_anx <- c(
  "No" = "#BDBDBD",
  "Yes" = "#8FC9F2"
)

pal_dep_tri <- c(
  "Minimal" = "#F8DDD4",
  "Mild" = "#F4A3A3",
  "Moderate-Severe" = "#E76F51"
)

pal_anx_tri <- c(
  "Minimal" = "#D0E8F9",
  "Mild" = "#8FC9F2",
  "Moderate-Severe" = "#4A90C2"
)

p_sex <- plot_donut_by_ibd(
  plot_dat,
  Sex_group,
  title = "Gender",
  palette = pal_sex
)

p_age <- plot_donut_by_ibd(
  plot_dat,
  Age_group,
  title = "Age",
  palette = pal_age
)

p_bmi <- plot_donut_by_ibd(
  plot_dat,
  BMI_group,
  title = "BMI",
  palette = pal_bmi
)

p_dep <- plot_donut_by_ibd(
  plot_dat,
  Depression_group,
  title = "Depression",
  palette = pal_dep
)

p_anx <- plot_donut_by_ibd(
  plot_dat,
  Anxiety_group,
  title = "Anxiety",
  palette = pal_anx
)

p_dep_tri <- plot_donut_by_ibd(
  plot_dat,
  Depression_tri,
  title = "Depression (3-level)",
  palette = pal_dep_tri
)

p_anx_tri <- plot_donut_by_ibd(
  plot_dat,
  Anxiety_tri,
  title = "Anxiety (3-level)",
  palette = pal_anx_tri
)


plot_score_compare <- function(data, score_var, ylab, 
                               ibd_fill = "#A7D3F2",
                               healthy_fill = "#BDBDBD",
                               test_method = c("wilcox", "t.test"),
                               y_limit = NULL) {
  
  test_method <- match.arg(test_method)
  score_var <- rlang::ensym(score_var)
  score_name <- rlang::as_name(score_var)
  
  df <- data %>%
    filter(!is.na(IBD_group), !is.na(!!score_var)) %>%
    mutate(
      IBD_group = factor(IBD_group, levels = c("Healthy", "IBD"))
    ) %>%
    dplyr::select(IBD_group, !!score_var)
  
  # summary: mean ± SE
  sum_df <- df %>%
    group_by(IBD_group) %>%
    summarise(
      n = n(),
      mean = mean(!!score_var, na.rm = TRUE),
      sd = sd(!!score_var, na.rm = TRUE),
      se = sd / sqrt(n),
      .groups = "drop"
    )
  
  # statistical test
  if (test_method == "wilcox") {
    pval <- wilcox.test(df[[score_name]] ~ df$IBD_group)$p.value
  } else {
    pval <- t.test(df[[score_name]] ~ df$IBD_group)$p.value
  }
  
  sig_lab <- dplyr::case_when(
    pval < 0.001 ~ "***",
    pval < 0.01  ~ "**",
    pval < 0.05  ~ "*",
    TRUE         ~ "ns"
  )
  
  y_max_data <- max(df[[score_name]], na.rm = TRUE)
  y_bar <- max(sum_df$mean + sum_df$se, na.rm = TRUE)
  y_sig <- max(y_max_data, y_bar) * 1.10
  
  p <- ggplot(sum_df, aes(x = IBD_group, y = mean, fill = IBD_group)) +
    geom_col(
      width = 0.62,
      color = NA
    ) +
    geom_errorbar(
      aes(ymin = mean - se, ymax = mean + se),
      width = 0.15,
      linewidth = 0.5
    ) +
    geom_jitter(
      data = df,
      aes(x = IBD_group, y = !!score_var),
      width = 0.15,
      height = 0,
      shape = 21,
      size = 2,
      stroke = 0.4,
      fill = "white",
      color = "grey25",
      inherit.aes = FALSE
    ) +
    annotate(
      "text",
      x = 1.5,
      y = y_sig,
      label = sig_lab,
      size = 5,
      fontface = "bold",
      family = figure_font
    ) +
    scale_fill_manual(
      values = c(
        "Healthy" = healthy_fill,
        "IBD" = ibd_fill
      )
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.12))
    ) +
    labs(
      x = NULL,
      y = ylab
    ) +
    theme_classic(base_family = figure_font) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 10, color = "black"),
      axis.text.y = element_text(size = 9, color = "black"),
      axis.title.y = element_text(size = 10, color = "black"),
      axis.line = element_line(linewidth = 0.6, color = "black"),
      axis.ticks = element_line(linewidth = 0.6, color = "black"),
      plot.margin = margin(8, 8, 8, 8)
    )

  if (!is.null(y_limit)) {
    p <- p + coord_cartesian(ylim = y_limit)
  }

  return(p)
}

plot_score_compare_violin <- function(data, score_var, ylab,
                                      fill_color = NULL,
                                      ibd_fill = NULL,
                                      healthy_fill = "#BDBDBD",
                                      test_method = c("wilcox", "t.test"),
                                      y_limit = NULL) {
  
  test_method <- match.arg(test_method)
  score_var <- rlang::ensym(score_var)
  score_name <- rlang::as_name(score_var)
  
  if (is.null(ibd_fill)) {
    ibd_fill <- ifelse(is.null(fill_color), "#A7D3F2", fill_color)
  }
  
  df <- data %>%
    filter(!is.na(IBD_group), !is.na(!!score_var)) %>%
    mutate(
      IBD_group = factor(IBD_group, levels = c("Healthy", "IBD"))
    ) %>%
    dplyr::select(IBD_group, !!score_var)
  
  # statistical test
  if (test_method == "wilcox") {
    pval <- wilcox.test(df[[score_name]] ~ df$IBD_group)$p.value
  } else {
    pval <- t.test(df[[score_name]] ~ df$IBD_group)$p.value
  }
  
  sig_lab <- dplyr::case_when(
    pval < 0.001 ~ "***",
    pval < 0.01  ~ "**",
    pval < 0.05  ~ "*",
    TRUE         ~ "ns"
  )
  
  y_max <- max(df[[score_name]], na.rm = TRUE)
  y_sig <- y_max * 1.10
  
  p <- ggplot(df, aes(x = IBD_group, y = !!score_var, fill = IBD_group)) +
    geom_violin(
      width = 0.75,
      alpha = 0.75,
      color = "grey30",
      linewidth = 0.4,
      trim = FALSE
    ) +
    geom_boxplot(
      width = 0.18,
      outlier.shape = NA,
      color = "grey20",
      fill = "white",
      linewidth = 0.45
    )  +
    annotate(
      "text",
      x = 1.5,
      y = y_sig,
      label = sig_lab,
      size = 5,
      fontface = "bold",
      family = figure_font
    ) +
    scale_fill_manual(
      values = c(
        "Healthy" = healthy_fill,
        "IBD" = ibd_fill
      )
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.15))
    ) +
    labs(
      x = NULL,
      y = ylab
    ) +
    theme_classic(base_family = figure_font) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 10, color = "black"),
      axis.text.y = element_text(size = 9, color = "black"),
      axis.title.y = element_text(size = 10, color = "black"),
      axis.line = element_line(linewidth = 0.6, color = "black"),
      axis.ticks = element_line(linewidth = 0.6, color = "black"),
      plot.margin = margin(8, 8, 8, 8)
    )
  
  if (!is.null(y_limit)) {
    p <- p + coord_cartesian(ylim = y_limit)
  }
  
  return(p)
}

p_dep_score <- plot_score_compare(
  plot_dat,
  PHQ9_sum,
  ylab = "Depression scores",
  ibd_fill = "#F4A3A3",
  healthy_fill = "#BDBDBD",
  test_method = "wilcox",
  y_limit = c(0, 30)
)

p_anx_score <- plot_score_compare(
  plot_dat,
  GAD7_sum,
  ylab = "Anxiety scores",
  ibd_fill = "#8FC9F2",
  healthy_fill = "#BDBDBD",
  test_method = "wilcox",
  y_limit = c(0, 30)
)

p_dep_score <- plot_score_compare_violin(
  plot_dat,
  PHQ9_sum,
  ylab = "Depression scores",
  fill_color = "#F4A3A3",
  test_method = "wilcox",
  y_limit = c(0, 30)
)

p_anx_score <- plot_score_compare_violin(
  plot_dat,
  GAD7_sum,
  ylab = "Anxiety scores",
  fill_color = "#8FC9F2",
  test_method = "wilcox",
  y_limit = c(0, 30)
)


p_top <- p_sex + p_age + p_bmi + 
  plot_layout(widths = c(1, 1, 1))

p_bottom <- p_dep + p_anx + 
  plot_layout(widths = c(1, 1))

p_all <- p_top / p_bottom +
  plot_annotation(
    title = "Cohort characteristics",
    theme = theme(
      plot.title = element_text(
        hjust = 0.5,
        size = 16,
        face = "bold",
        family = figure_font
      )
    )
  )

p_all
ggsave(file.path(out_dir, "cohort_characteristics_v1.pdf"),
       p_all, width = 14, height = 8, device = cairo_pdf)

p_top <- p_sex + p_age + p_bmi +
  plot_layout(widths = c(1, 1, 1))

p_bottom <- p_dep + p_dep_score + p_anx + p_anx_score +
  plot_layout(widths = c(1.15, 0.85, 1.15, 0.85))

p_all <- p_top / p_bottom +
  plot_annotation(
    title = "Cohort characteristics",
    theme = theme(
      plot.title = element_text(
        hjust = 0.5,
        size = 16,
        face = "bold",
        family = figure_font
      )
    )
  )

p_all
ggsave(file.path(out_dir, "cohort_characteristics_v2.pdf"),
       p_all, width = 16, height = 8, device = cairo_pdf)

p_top <- p_sex + p_age + p_bmi +
  plot_layout(widths = c(1, 1, 1))

p_bottom <- p_dep_tri + p_dep_score + p_anx_tri + p_anx_score +
  plot_layout(widths = c(1.15, 0.85, 1.15, 0.85))

p_all <- p_top / p_bottom +
  plot_annotation(
    title = "Cohort characteristics",
    theme = theme(
      plot.title = element_text(
        hjust = 0.5,
        size = 16,
        face = "bold",
        family = figure_font
      )
    )
  )

p_all
ggsave(file.path(out_dir, "cohort_characteristics_v3.pdf"),
       p_all, width = 16, height = 8, device = cairo_pdf)

### Density figure ####

# -----------------------------
# 0. 可随时修改的主题颜色变量
# -----------------------------
plot_theme_cols <- list(
  phq_fill   = "#E69F00",
  phq_line   = "#B36B00",
  
  gad_fill   = "#56B4E9",
  gad_line   = "#0072B2",
  
  cutoff     = "#D55E00",
  text       = "#222222",
  title      = "#111111",
  axis       = "#333333",
  grid       = "#EAEAEA",
  border     = "#333333",
  
  bg         = "white",
  panel_bg   = "white"
)

# -----------------------------
# 1. 统一主题
# -----------------------------
my_density_theme <- theme_classic(base_size = 11, base_family = figure_font) +
  theme(
    plot.background = element_rect(
      fill = plot_theme_cols$bg,
      color = NA
    ),
    panel.background = element_rect(
      fill = plot_theme_cols$panel_bg,
      color = NA
    ),
    panel.grid.major.y = element_line(
      color = plot_theme_cols$grid,
      linewidth = 0.3
    ),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),

    axis.line = element_blank(),
    axis.ticks = element_line(
      color = plot_theme_cols$axis,
      linewidth = 0.4
    ),
    axis.text = element_text(
      color = plot_theme_cols$axis,
      size = 9
    ),
    axis.title = element_text(
      color = plot_theme_cols$text,
      size = 10
    ),
    plot.title = element_text(
      color = plot_theme_cols$title,
      size = 11,
      hjust = 0
    ),
    legend.text = element_text(
      color = plot_theme_cols$text,
      size = 9
    ),
    legend.background = element_blank(),
    legend.key = element_blank(),
    panel.border = element_rect(
      color = plot_theme_cols$border,
      fill = NA,
      linewidth = 0.4
    )
  )

# -----------------------------
# 2. 整理成长格式数据
# -----------------------------
density_data <- cohort %>%
  dplyr::select(PHQ9_sum, GAD7_sum, IBD_subtype) %>%
  pivot_longer(
    cols = c(PHQ9_sum, GAD7_sum),
    names_to = "scale",
    values_to = "score"
  ) %>%
  filter(!is.na(score)) %>%
  mutate(
    scale = recode(
      scale,
      PHQ9_sum = "PHQ-9",
      GAD7_sum = "GAD-7"
    ),
    IBD_subtype = as.character(IBD_subtype)
  )

# -----------------------------
# 3. PHQ-9 / GAD-7 颜色
# -----------------------------
scale_cols <- c(
  "PHQ-9" = plot_theme_cols$phq_line,
  "GAD-7" = plot_theme_cols$gad_line
)

scale_fills <- c(
  "PHQ-9" = plot_theme_cols$phq_fill,
  "GAD-7" = plot_theme_cols$gad_fill
)

# -----------------------------
# 4. 主图
# -----------------------------
main_plot <- ggplot(
  density_data,
  aes(x = score, fill = scale, color = scale)
) +
  geom_density(
    alpha = 0.35,
    linewidth = 0.8,
    adjust = 1
  ) +
  geom_vline(
    xintercept = 10,
    linetype = "dashed",
    color = plot_theme_cols$cutoff,
    linewidth = 0.6
  ) +
  annotate(
    "text",
    x = 10.5,
    y = Inf,
    label = "Cutoff = 10",
    vjust = 2,
    hjust = 0,
    color = plot_theme_cols$cutoff,
    size = 3.5,
    family = figure_font
  ) +
  scale_fill_manual(values = scale_fills) +
  scale_color_manual(values = scale_cols) +
  scale_x_continuous(
    limits = c(0, 27),
    breaks = seq(0, 27, 3),
    expand = c(0, 0)
  ) +
  labs(
    x = "Score",
    y = "Density",
    fill = NULL,
    color = NULL
  ) +
  my_density_theme +
  theme(
    legend.position = c(0.18, 0.85)
  ) +
  scale_y_continuous(
    limits = c(0,0.20),
    expand = c(0,0)
  )

# -----------------------------
# 5. CD 子图
# -----------------------------
cd_plot <- ggplot(
  density_data %>% filter(IBD_subtype == "CD"),
  aes(x = score, fill = scale, color = scale)
) +
  geom_density(
    alpha = 0.35,
    linewidth = 0.6,
    adjust = 1
  ) +
  geom_vline(
    xintercept = 10,
    linetype = "dashed",
    color = plot_theme_cols$cutoff,
    linewidth = 0.4
  ) +
  scale_fill_manual(values = scale_fills) +
  scale_color_manual(values = scale_cols) +
  scale_x_continuous(
    limits = c(0, 27),
    breaks = seq(0, 27, 9),
    expand = c(0, 0)
  ) +
  labs(
    title = "CD",
    x = "",
    y = ""
  ) +
  my_density_theme +
  theme(
    legend.position = "none",
    plot.margin = margin(2, 2, 8, 8),
    plot.title = element_text(
      color = plot_theme_cols$title,
      size = 10,
      hjust = 0
    )
  ) + scale_y_continuous(
    limits = c(0,0.20),
    expand = c(0, 0))

# -----------------------------
# 6. UC 子图
# -----------------------------
uc_plot <- ggplot(
  density_data %>% filter(IBD_subtype == "UC"),
  aes(x = score, fill = scale, color = scale)
) +
  geom_density(
    alpha = 0.35,
    linewidth = 0.6,
    adjust = 1
  ) +
  geom_vline(
    xintercept = 10,
    linetype = "dashed",
    color = plot_theme_cols$cutoff,
    linewidth = 0.4
  ) +
  scale_fill_manual(values = scale_fills) +
  scale_color_manual(values = scale_cols) +
  scale_x_continuous(
    limits = c(0, 27),
    breaks = seq(0, 27, 9),
    expand = c(0, 0)
  ) +
  labs(
    title = "UC",
    x = "",
    y = ""
  ) +
  my_density_theme +
  theme(
    legend.position = "none",
    plot.margin = margin(2, 2, 8, 8),
    plot.title = element_text(
      color = plot_theme_cols$title,
      size = 10,
      hjust = 0
    )
  ) + scale_y_continuous(
    limits = c(0,0.20),
    expand = c(0, 0)
  )

# -----------------------------
# 7. 合并图
# -----------------------------
final_plot <- ggdraw() +
  draw_plot(main_plot) +
  
  draw_plot(
    cd_plot,
    x = 0.60,
    y = 0.58,
    width = 0.35,
    height = 0.30
  ) +
  
  draw_plot(
    uc_plot,
    x = 0.60,
    y = 0.24,
    width = 0.35,
    height = 0.30
  ) +
  
  draw_label(
    label = paste0("n = ", format(nrow(cohort), big.mark = ",")),
    x = 0.90,
    y = 0.92,
    hjust = 1,
    size = 9,
    color = plot_theme_cols$text,
    fontfamily = figure_font
  ) +

  draw_label(
    label = paste0("n = ", format(sum(cohort$IBD_subtype == "CD", na.rm = TRUE), big.mark = ",")),
    x = 0.92,
    y = 0.83,
    hjust = 1,
    size = 8,
    color = plot_theme_cols$text,
    fontfamily = figure_font
  ) +

  draw_label(
    label = paste0("n = ", format(sum(cohort$IBD_subtype == "UC", na.rm = TRUE), big.mark = ",")),
    x = 0.92,
    y = 0.49,
    hjust = 1,
    size = 8,
    color = plot_theme_cols$text,
    fontfamily = figure_font
  )


final_plot
ggsave(file.path(out_dir, "phq_gad_density.pdf"),
       final_plot, width = 8, height = 6, device = cairo_pdf)

###饼图####
# -----------------------------
# 1. 可修改参数
# -----------------------------
pie_cols <- list(
  high_phq = "#E76F51",
  mid_phq  = "#F09E83",
  low_phq  = "#F4D6C6",
  high_gad = "#2A9D8F",
  mid_gad  = "#7EC8BE",
  low_gad  = "#D7F0EA",
  border   = "white",
  text     = "#222222",
  title    = "#111111"
)

explode_value <- 0.18
pie_radius <- 2
label_size <- 5
title_size <- 12


# -----------------------------
# 2. 构造饼图数据函数（二分类，保留备用）
# -----------------------------
make_pie_data <- function(data, score_var, cutoff = 10, high_label = ">=10") {

  score_vec <- data[[score_var]]

  n_total <- sum(!is.na(score_vec))
  n_high  <- sum(score_vec >= cutoff, na.rm = TRUE)
  n_low   <- n_total - n_high

  pie_data <- data.frame(
    group = c(
      paste0(high_label),
      paste0("<10")
    ),
    n = c(n_high, n_low)
  ) %>%
    mutate(
      percentage = 100 * n / sum(n),
      label = paste0(group, "\n", n, "/", sum(n), "\n", sprintf("%.1f%%", percentage)),
      focus = ifelse(group == high_label, explode_value, 0),
      end_angle = 2 * pi * cumsum(percentage) / 100,
      start_angle = lag(end_angle, default = 0),
      mid_angle = 0.5 * (start_angle + end_angle),
      label_x = 1.05 * sin(mid_angle),
      label_y = 1.05 * cos(mid_angle)
    )

  return(pie_data)
}

# -----------------------------
# 2b. 构造三分类饼图数据函数
# -----------------------------
make_pie_data_tri <- function(data, score_var, prefix = "PHQ-9") {

  score_vec <- data[[score_var]]
  score_vec <- score_vec[!is.na(score_vec)]
  n_total <- length(score_vec)

  labels <- c(
    paste0(prefix, " < 5"),
    paste0(prefix, " 5-9"),
    paste0(prefix, " >= 10")
  )

  n_low  <- sum(score_vec < 5)
  n_mid  <- sum(score_vec >= 5 & score_vec < 10)
  n_high <- sum(score_vec >= 10)

  pie_data <- data.frame(
    group = factor(labels, levels = labels),
    n = c(n_low, n_mid, n_high)
  ) %>%
    mutate(
      percentage = 100 * n / sum(n),
      label = paste0(
        group, "\n", n, "/", sum(n), "\n",
        sprintf("%.1f%%", percentage)
      ),
      focus = case_when(
        group == labels[3] ~ explode_value,
        group == labels[2] ~ explode_value * 0.5,
        TRUE ~ 0
      ),
      end_angle = 2 * pi * cumsum(percentage) / 100,
      start_angle = lag(end_angle, default = 0),
      mid_angle = 0.5 * (start_angle + end_angle),
      label_x = 1.05 * sin(mid_angle),
      label_y = 1.05 * cos(mid_angle)
    )

  return(pie_data)
}


# -----------------------------
# 3. 生成 PHQ-9 和 GAD-7 三分类数据
# -----------------------------
phq9_pie_data <- make_pie_data_tri(
  data = cohort,
  score_var = "PHQ9_sum",
  prefix = "PHQ-9"
)

gad7_pie_data <- make_pie_data_tri(
  data = cohort,
  score_var = "GAD7_sum",
  prefix = "GAD-7"
)


# -----------------------------
# 4. 饼图函数
# -----------------------------
plot_exploded_pie <- function(pie_data, title, fill_values) {
  
  ggplot() +
    geom_arc_bar(
      data = pie_data,
      stat = "pie",
      aes(
        x0 = 0,
        y0 = 0,
        r0 = 0,
        r = pie_radius,
        amount = n,
        fill = group,
        explode = focus
      ),
      color = pie_cols$border,
      linewidth = 0.8,
      show.legend = FALSE
    ) +
    geom_text(
      data = pie_data,
      aes(
        x = label_x,
        y = label_y,
        label = label
      ),
      size = label_size,
      color = pie_cols$text,
      lineheight = 0.9,
      family = figure_font
    ) +
    scale_fill_manual(values = fill_values) +
    coord_fixed() +
    labs(title = title) +
    theme_no_axes() +
    theme(
      text = element_text(family = figure_font),
      plot.title = element_text(
        color = pie_cols$title,
        size = title_size,
        hjust = 0.5,
        face = "bold"
      ),
      plot.margin = margin(5, 5, 5, 5)
    )
}


# -----------------------------
# 5. 分别绘制两个饼图
# -----------------------------
p_phq9_pie <- plot_exploded_pie(
  pie_data = phq9_pie_data,
  title = "PHQ-9 Depression Symptoms",
  fill_values = c(
    "PHQ-9 < 5"  = pie_cols$low_phq,
    "PHQ-9 5-9"  = pie_cols$mid_phq,
    "PHQ-9 >= 10" = pie_cols$high_phq
  )
)

p_gad7_pie <- plot_exploded_pie(
  pie_data = gad7_pie_data,
  title = "GAD-7 Anxiety Symptoms",
  fill_values = c(
    "GAD-7 < 5"  = pie_cols$low_gad,
    "GAD-7 5-9"  = pie_cols$mid_gad,
    "GAD-7 >= 10" = pie_cols$high_gad
  )
)


# -----------------------------
# 6. 合并两个饼图
# -----------------------------
final_pie_plot <- p_phq9_pie + p_gad7_pie +
  plot_layout(ncol = 2)

final_pie_plot

p <- final_pie_plot + plot_layout(widths = c(1,1))
p
ggsave(file.path(out_dir, "phq_gad_pie.pdf"),
       p, width = 12, height = 6, device = cairo_pdf)
# library(export)
# graph2ppt(file="Pie chart.pptx", width=7, height=5) #导出为PPT，继续调整细节
#我调整的细节包括：饼图与图例之间的间距，基因名称斜体、饼图中的标签和百分比的位置等

###cohort ibd：infla_index ~ phq9; gad7 :spearman correlation####
bio_markers <- c(
  "Calprotectin",
  "SII",
  "NLR",
  "PLR",
  "LMR",
  "hsCRP"
)

psych_scores <- c(
  "PHQ9_sum" = "PHQ-9",
  "GAD7_sum" = "GAD-7"
)

# 每个 biomarker 一套色系：
# PHQ-9 = 深色
# GAD-7 = 浅色
pal_marker_score <- c(
  "Calprotectin_PHQ-9" = "#D95F5F",
  "Calprotectin_GAD-7" = "#F4A3A3",
  
  "SII_PHQ-9" = "#E6862D",
  "SII_GAD-7" = "#FDBF6F",
  
  "NLR_PHQ-9" = "#C69214",
  "NLR_GAD-7" = "#FFE3B3",
  
  "PLR_PHQ-9" = "#4F9A63",
  "PLR_GAD-7" = "#B7E1B0",
  
  "LMR_PHQ-9" = "#2F9B85",
  "LMR_GAD-7" = "#7AC6A4",
  
  "hsCRP_PHQ-9" = "#4A90C2",
  "hsCRP_GAD-7" = "#8FC9F2"
)

cor_dat <- cohort %>%
  dplyr::select(all_of(c(names(psych_scores), bio_markers))) %>%
  pivot_longer(
    cols = all_of(bio_markers),
    names_to = "marker",
    values_to = "marker_value"
  ) %>%
  pivot_longer(
    cols = all_of(names(psych_scores)),
    names_to = "score_type",
    values_to = "score_value"
  ) %>%
  mutate(
    score_type = recode(score_type, !!!psych_scores),
    marker = factor(marker, levels = bio_markers),
    score_type = factor(score_type, levels = c("PHQ-9", "GAD-7")),
    color_group = paste(marker, score_type, sep = "_"),
    marker_value_log = log10(marker_value + 1)
  ) %>%
  filter(
    !is.na(marker_value),
    !is.na(score_value),
    is.finite(marker_value_log)
  )

cor_res <- cor_dat %>%
  group_by(marker, score_type) %>%
  summarise(
    n = n(),
    rho = cor(
      score_value,
      marker_value,
      method = "spearman",
      use = "complete.obs"
    ),
    p = cor.test(
      score_value,
      marker_value,
      method = "spearman",
      exact = FALSE
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    q = p.adjust(p, method = "BH"),
    color_group = paste(marker, score_type, sep = "_"),
    label = sprintf("%s: rho = %.2f, p = %.2e", score_type, rho, p)
  )

fwrite(
  as.data.frame(cor_res),
  file.path(out_dir, "inflammation_spearman_results.csv"),
  na = "NA"
)

label_dat <- cor_res %>%
  mutate(
    label_x = -Inf,
    label_y = Inf,
    hjust_value = -0.05,
    vjust_value = case_when(
      score_type == "PHQ-9" ~ 1.15,
      score_type == "GAD-7" ~ 2.45
    )
  )

p_cor_marker <- ggplot(
  cor_dat,
  aes(
    x = score_value,
    y = marker_value_log,
    color = color_group,
    fill = color_group
  )
) +
  geom_point(
    alpha = 0.45,
    size = 1.8
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    alpha = 0.18,
    linewidth = 0.9
  ) +
  geom_text(
    data = label_dat,
    aes(
      x = label_x,
      y = label_y,
      label = label,
      color = color_group,
      vjust = vjust_value
    ),
    inherit.aes = FALSE,
    hjust = -0.05,
    size = 3.2,
    fontface = "bold",
    family = figure_font,
    show.legend = FALSE
  ) +
  facet_wrap(
    ~ marker,
    scales = "free_y",
    ncol = 3
  ) +
  scale_color_manual(values = pal_marker_score) +
  scale_fill_manual(values = pal_marker_score) +
  labs(
    x = "Psychological symptom score",
    y = "log10(Biomarker + 1)",
    title = "Correlation between inflammatory markers and psychological symptoms",
    subtitle = "PHQ-9 and GAD-7 are shown with different shades within each biomarker-specific color palette"
  ) +
  theme_classic(base_size = 13, base_family = figure_font) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11, color = "grey35"),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "none",
    panel.spacing = unit(1.0, "lines")
  )

p_cor_marker
ggsave(file.path(out_dir, "inflammation_correlation.pdf"),
       p_cor_marker, width = 14, height = 9, device = cairo_pdf)







#### Partial Spearman correlation: PHQ-9/GAD-7 ~ inflammation indices ####
# Adjusted for Age, Sex, BMI, hPDI, and site.

cohort <- cohort %>%
  mutate(
    Sex = factor(Sex)
  )

covars <- c("Age", "Sex", "BMI", "hPDI", "site")

run_partial_spearman <- function(data, marker, score_var, score_label) {
  
  dat <- data %>%
    dplyr::select(all_of(c(marker, score_var, covars))) %>%
    filter(complete.cases(.))
  
  fml <- as.formula(
    paste0(marker, " | ", score_var, " ~ ", paste(covars, collapse = " + "))
  )
  
  fit <- partial_Spearman(
    formula = fml,
    data = dat,
    fit.x = "orm",
    fit.y = "orm",
    fisher = TRUE,
    conf.int = 0.95
  )
  
  tibble(
    marker = marker,
    score_var = score_var,
    score_type = score_label,
    n = fit$data.points,
    rho = as.numeric(fit$TS$TB$ts),
    stderr = sqrt(as.numeric(fit$TS$TB$var)),
    p = as.numeric(fit$TS$TB$pval),
    lower = as.numeric(fit$TS$TB$lower),
    upper = as.numeric(fit$TS$TB$upper)
  )
}

partial_cor_res <- expand.grid(
  marker = bio_markers,
  score_var = names(psych_scores),
  stringsAsFactors = FALSE
) %>%
  as_tibble() %>%
  mutate(
    score_type = recode(score_var, !!!psych_scores)
  ) %>%
  pmap_dfr(function(marker, score_var, score_type) {
    run_partial_spearman(
      data = cohort,
      marker = marker,
      score_var = score_var,
      score_label = score_type
    )
  }) %>%
  mutate(
    q = p.adjust(p, method = "BH"),
    marker = factor(marker, levels = bio_markers),
    score_type = factor(score_type, levels = c("PHQ-9", "GAD-7")),
    color_group = paste(marker, score_type, sep = "_"),
    label = sprintf(
      "%s: partial rho = %.2f, p = %.2e",
      score_type, rho, p
    )
  )

fwrite(
  as.data.frame(partial_cor_res),
  file.path(out_dir, "inflammation_partial_spearman_results.csv"),
  na = "NA"
)

partial_cor_res

cor_dat <- cohort %>%
  dplyr::select(all_of(c(names(psych_scores), bio_markers))) %>%
  pivot_longer(
    cols = all_of(bio_markers),
    names_to = "marker",
    values_to = "marker_value"
  ) %>%
  pivot_longer(
    cols = all_of(names(psych_scores)),
    names_to = "score_type",
    values_to = "score_value"
  ) %>%
  mutate(
    score_type = recode(score_type, !!!psych_scores),
    marker = factor(marker, levels = bio_markers),
    score_type = factor(score_type, levels = c("PHQ-9", "GAD-7")),
    color_group = paste(marker, score_type, sep = "_")
  ) %>%
  filter(
    !is.na(marker_value),
    !is.na(score_value),
    is.finite(marker_value),
    is.finite(score_value)
  )

label_dat <- partial_cor_res %>%
  mutate(
    marker = factor(marker, levels = bio_markers),
    score_type = factor(score_type, levels = c("PHQ-9", "GAD-7")),
    color_group = paste(marker, score_type, sep = "_"),
    label = sprintf(
      "%s: rho = %.2f, p = %.2e",
      score_type, rho, p
    ),
    label_x = -Inf,
    label_y = Inf,
    vjust_value = case_when(
      score_type == "PHQ-9" ~ 1.15,
      score_type == "GAD-7" ~ 2.45
    )
  )

p_partial_cor_marker <- ggplot(
  cor_dat,
  aes(
    x = score_value,
    y = marker_value,
    color = color_group,
    fill = color_group,
    shape = score_type
  )
) +
  geom_point(
    alpha = 0.55,
    size = 2.0
  ) +
  geom_smooth(
    aes(
      group = color_group,
      linetype = score_type
    ),
    method = "lm",
    se = TRUE,
    alpha = 0.18,
    linewidth = 0.9,
    show.legend = FALSE
  ) +
  geom_text(
    data = label_dat,
    aes(
      x = label_x,
      y = label_y,
      label = label,
      color = color_group,
      vjust = vjust_value
    ),
    inherit.aes = FALSE,
    hjust = -0.05,
    size = 3.2,
    fontface = "bold",
    family = figure_font,
    show.legend = FALSE
  ) +
  ggh4x::facet_wrap2(
    ~ marker,
    scales = "free_y",
    ncol = 3,
    axes = "x",
    remove_labels = "none"
  ) +
  scale_color_manual(
    values = pal_marker_score,
    guide = "none"
  ) +
  scale_fill_manual(
    values = pal_marker_score,
    guide = "none"
  ) +
  scale_shape_manual(
    values = c(
      "PHQ-9" = 16,
      "GAD-7" = 17
    ),
    name = NULL
  ) +
  scale_linetype_manual(
    values = c(
      "PHQ-9" = "dashed",
      "GAD-7" = "solid"
    ),
    guide = "none"
  ) +
  guides(
    shape = guide_legend(
      override.aes = list(
        size = 3,
        alpha = 1,
        color = "grey30"
      )
    )
  ) +
  labs(
    x = "Psychological symptom score",
    y = "Biomarker level"
  ) +
  theme_classic(base_size = 13, base_family = figure_font) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11, color = "grey35"),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.key.width = unit(1.0, "cm"),
    legend.key.height = unit(0.4, "cm"),
    panel.spacing = unit(1.0, "lines")
  )

p_partial_cor_marker
ggsave(file.path(out_dir, "inflammation_partial_correlation.pdf"),
       p_partial_cor_marker, width = 14, height = 9, device = cairo_pdf)
