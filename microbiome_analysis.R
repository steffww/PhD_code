# 03_abundance_analysis.R — 丰度展示
# =============================================================================
library(microeco)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(showtext)

font_add("Arial",
         regular    = "/System/Library/Fonts/Supplemental/Arial.ttf",
         bold       = "/System/Library/Fonts/Supplemental/Arial Bold.ttf",
         italic     = "/System/Library/Fonts/Supplemental/Arial Italic.ttf",
         bolditalic = "/System/Library/Fonts/Supplemental/Arial Bold Italic.ttf")
showtext_auto()
showtext_opts(dpi = 300)

source("R/00_setup.R")

# =============================================================================
# 2. 构建 microtable 对象（统一 QC）
# =============================================================================
samp_all  = readRDS('./metadata/cohort_all_449.rds')
mt_all = readRDS("./microbiome/mt_all_449.rds")

# --- 2a. 全部样本（IBD + HC），用于 IBD vs HC 对比 ---

mt_all$cal_abund(select_cols = 1:6, rel = TRUE)
cat("\nmt_all:", nrow(mt_all$otu_table), "taxa x", ncol(mt_all$otu_table), "samples\n")

# --- 2b. 仅 IBD 样本，用于 dep/anx 内部对比 ---

samp_ibd  = readRDS("./metadata/cohort_ibd_282.rds")
mt_ibd = readRDS("./microbiome/mt_ibd_282.rds")

mt_ibd$cal_abund(select_cols = 1:6, rel = TRUE)
cat("mt_ibd:", nrow(mt_ibd$otu_table), "taxa x", ncol(mt_ibd$otu_table), "samples\n")


# 导入后，把新建的分组列同步到 microtable 的 sample_table
mt_ibd$sample_table$dep_group <- factor(
  mt_ibd$sample_table$PHQ9_bi10,
  levels = c(0, 1),
  labels = c("Non-depressed", "Depressed")
)

mt_ibd$sample_table$anx_group <- factor(
  mt_ibd$sample_table$GAD7_bi10,
  levels = c(0, 1),
  labels = c("Non-anxious", "Anxious")
)

# --- 三分组：Minimal (0-4) / Mild-Moderate (5-9) / Moderate-Severe (>=10) ---
mt_ibd$sample_table$dep_3grp <- cut(
  mt_ibd$sample_table$PHQ9_sum,
  breaks = c(-Inf, 4, 9, Inf),
  labels = c("Minimal(0-4)", "Mild(5-9)", "Moderate-Severe(>=10)")
)

mt_ibd$sample_table$anx_3grp <- cut(
  mt_ibd$sample_table$GAD7_sum,
  breaks = c(-Inf, 4, 9, Inf),
  labels = c("Minimal(0-4)", "Mild(5-9)", "Moderate-Severe(>=10)")
)

mt_all$sample_table$Disease <- ifelse(mt_all$sample_table$IBD == "IBD", "IBD", "HC")
mt_all$sample_table$Disease <- factor(mt_all$sample_table$Disease, levels = c("HC", "IBD"))

# 验证新建列在 microtable 里存在
cat("\nmt_ibd dep_group:\n"); print(table(mt_ibd$sample_table$dep_group, useNA = "always"))
cat("mt_ibd anx_group:\n"); print(table(mt_ibd$sample_table$anx_group, useNA = "always"))
cat("mt_ibd dep_3grp:\n"); print(table(mt_ibd$sample_table$dep_3grp, useNA = "always"))
cat("mt_ibd anx_3grp:\n"); print(table(mt_ibd$sample_table$anx_3grp, useNA = "always"))
cat("mt_all Disease:\n"); print(table(mt_all$sample_table$Disease, useNA = "always"))

# =========================
# 主题
# =========================
theme_half_open <- function() {
  theme_classic(base_size = 18) +
    theme(
      panel.background  = element_rect(fill = "white", color = NA),
      text = element_text(family = "Arial"),
      plot.background   = element_rect(fill = "white", color = NA),
      panel.border      = element_blank(),
      axis.line.x       = element_line(color = "black", linewidth = 1.0),
      axis.line.y       = element_line(color = "black", linewidth = 1.0),
      axis.ticks        = element_line(color = "black", linewidth = 1.0),
      axis.ticks.length = unit(0.2, "cm"),
      axis.text         = element_text(color = "black", size = 16),
      axis.text.x       = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.title        = element_text(color = "black", size = 10),
      plot.title        = element_text(size = 10, face = "bold", hjust = 0.5),
      legend.title      = element_text(size = 14, face = "bold"),
      legend.text       = element_text(size = 12),
      legend.background = element_rect(fill = "white", color = NA),
      panel.grid        = element_blank()
    )
}

# # 物种堆积图配色（15色，高区分度）
# pal_taxa <- c(
#   "#D95F5F", "#E6862D", "#C69214", "#4F9A63", "#2F9B85",
#   "#4A90C2", "#6A5ACD", "#E57373", "#FDBF6F", "#B7E1B0",
#   "#8FC9F2", "#DDA0DD", "#F4A3A3", "#7AC6A4", "#FFE3B3",
#   "grey70"  # 最后一个留给 "Others"
# )


pal_taxa <- c(
  "#F4A3A3",  # 淡红
  "#FDBF6F",  # 淡橙
  "#FFE3B3",  # 淡黄
  "#B7E1B0",  # 淡绿
  "#7AC6A4",  # 薄荷绿
  "#8FC9F2",  # 淡蓝
  "#B8A9D4",  # 淡紫
  "#F2C1D1",  # 淡粉
  "#C5D8A4",  # 黄绿
  "#A8D8DC",  # 青色
  "#D4C5A9",  # 米色
  "#C9B1FF",  # 浅薰衣草
  "#FFD6AA",  # 浅杏
  "#A3C4BC",  # 灰绿
  "#E0BBE4",  # 浅兰花
  "grey80"    # Others
)


# =============================================================================
# 3. 堆积图函数
# =============================================================================
plot_stacked_bar <- function(mt_obj, group_var, tax_level = "Phylum",
                             top_n = 10, title_text = "") {
  
  t_abund <- trans_abund$new(
    dataset   = mt_obj,
    taxrank   = tax_level,
    ntaxa     = top_n,
    groupmean = group_var
  )
  
  plot_data <- t_abund$data_abund
  plot_data$Taxonomy <- factor(plot_data$Taxonomy, levels = rev(unique(plot_data$Taxonomy)))
  
  p <- ggplot(plot_data, aes(x = Sample, y = Abundance, fill = Taxonomy)) +
    geom_bar(stat = "identity", position = "stack", width = 0.7) +
    scale_y_continuous(expand = c(0, 0), labels = function(x) paste0(x, "%")) +
    labs(x = NULL, y = "Relative Abundance (%)",
         title = title_text, fill = tax_level) +
    scale_fill_manual(values = pal_taxa) +
    theme_half_open() +
    theme(legend.position = "right")
  
  return(p)
}

# =============================================================================
# 4. Phylum 层面堆积图
# =============================================================================

p_phylum_disease <- plot_stacked_bar(
  mt_all, group_var = "Disease", tax_level = "Phylum",
  top_n = 10, title_text = "Phylum — IBD vs HC"
)

p_phylum_dep <- plot_stacked_bar(
  mt_ibd, group_var = "dep_group", tax_level = "Phylum",
  top_n = 10, title_text = "Phylum — Depression"
)

p_phylum_anx <- plot_stacked_bar(
  mt_ibd, group_var = "anx_group", tax_level = "Phylum",
  top_n = 10, title_text = "Phylum — Anxiety"
)

p_phylum_all <- p_phylum_disease + p_phylum_dep + p_phylum_anx
print(p_phylum_all)

ggsave("./results/abundance_phylum_stacked.pdf", p_phylum_all,
       width = 20, height = 8, device = cairo_pdf)

# --- Phylum 三分组 ---
p_phylum_dep3 <- plot_stacked_bar(
  mt_ibd, group_var = "dep_3grp", tax_level = "Phylum",
  top_n = 10, title_text = "Phylum — Depression (3-group)"
)

p_phylum_anx3 <- plot_stacked_bar(
  mt_ibd, group_var = "anx_3grp", tax_level = "Phylum",
  top_n = 10, title_text = "Phylum — Anxiety (3-group)"
)

p_phylum_3grp <- p_phylum_dep3 + p_phylum_anx3
print(p_phylum_3grp)

ggsave("./results/abundance_phylum_stacked_3grp.pdf", p_phylum_3grp,
       width = 14, height = 8, device = cairo_pdf)

# =============================================================================
# 5. Genus 层面堆积图
# =============================================================================

# 方案一：Nature 风格（饱和度适中，经典）
pal_nature <- c(
  "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
  "#8491B4", "#91D1C2", "#DC9041", "#7E6148", "#B09C85",
  "grey80"
)

# 方案二：Lancet 风格（沉稳专业）
pal_lancet <- c(
  "#00468B", "#ED0000", "#42B540", "#0099B4", "#925E9F",
  "#FDAF91", "#AD002A", "#ADB6B6", "#1B1919", "#E6A024",
  "grey80"
)

# 方案三：JAMA 风格（柔和清爽）
pal_jama <- c(
  "#374E55", "#DF8F44", "#00A1D5", "#B24745", "#79AF97",
  "#6A6599", "#80796B", "#E6A024", "#5A9EC4", "#CD534C",
  "grey80"
)

# 方案四：ggsci NPG 风格（Nature Publishing Group）
pal_npg <- c(
  "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
  "#8491B4", "#91D1C2", "#DC9041", "#7E6148", "#B09C85",
  "grey80"
)


plot_stacked_bar <- function(mt_obj, group_var, tax_level = "Phylum",
                             top_n = 8, title_text = "") {
  
  t_abund <- trans_abund$new(
    dataset   = mt_obj,
    taxrank   = tax_level,
    ntaxa     = top_n,
    groupmean = group_var
  )
  
  plot_data <- t_abund$data_abund
  
  # 把非 top_n 的全部归为 Others
  top_taxa <- plot_data %>%
    group_by(Taxonomy) %>%
    summarise(mean_abund = mean(Abundance)) %>%
    arrange(desc(mean_abund)) %>%
    head(top_n) %>%
    pull(Taxonomy)
  
  plot_data$Taxonomy <- ifelse(plot_data$Taxonomy %in% top_taxa,
                               plot_data$Taxonomy, "Others")
  
  # 重新汇总 Others
  plot_data <- plot_data %>%
    group_by(Sample, Taxonomy) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop")
  
  # 排序：top taxa 按丰度排，Others 放最底
  taxa_order <- c(rev(top_taxa), "Others")
  plot_data$Taxonomy <- factor(plot_data$Taxonomy, levels = taxa_order)
  
  # 配色
  taxa_colors <- setNames(pal_taxa[1:top_n], rev(top_taxa))
  taxa_colors["Others"] <- "grey80"
  
  p <- ggplot(plot_data, aes(x = Sample, y = Abundance, fill = Taxonomy)) +
    geom_bar(stat = "identity", position = "stack", width = 0.7) +
    scale_fill_manual(values = taxa_colors) +
    scale_y_continuous(expand = c(0, 0), labels = function(x) paste0(x, "%")) +
    labs(x = NULL, y = "Relative Abundance (%)",
         title = title_text, fill = tax_level) +
    theme_half_open() +
    theme(legend.position = "right")
  
  return(p)
}

# =============================================================================
# 6. Family 层面堆积图
# =============================================================================

# ---- 并集：三张图出现过的所有 family（按大致丰度排序）----
all_taxa <- c(
  "Lachnospiraceae",
  "Bacteroidaceae",
  "Oscillospiraceae",
  "Enterobacteriaceae",
  "Bifidobacteriaceae",
  "Streptococcaceae",
  "Enterococcaceae",
  "Prevotellaceae",        # 仅 panel 1
  "Erysipelotrichaceae",   # 仅 panel 2
  "Tannerellaceae"         # 仅 panel 3
)


# 图例顺序：高丰度在上（rev），Others 垫底
legend_order <- c(rev(all_taxa), "Others")

# legend 顺序不变
legend_order <- c(rev(all_taxa), "Others")
#shared_colors <- pal_family[legend_order]

shared_colors <- c(
  "Tannerellaceae"       = pal_taxa[1],
  "Erysipelotrichaceae"  = pal_taxa[2],
  "Prevotellaceae"       = pal_taxa[3],
  "Enterococcaceae"      = pal_taxa[4],
  "Streptococcaceae"     = pal_taxa[5],
  "Bifidobacteriaceae"   = pal_taxa[6],
  "Enterobacteriaceae"   = pal_taxa[7],
  "Oscillospiraceae"     = pal_taxa[8],
  "Bacteroidaceae"       = pal_taxa[9],
  "Lachnospiraceae"      = pal_taxa[10],
  "Others"               = "grey80"
)

plot_stacked_bar <- function(mt_obj, group_var, tax_level = "Phylum",
                             top_n = 8, title_text = "",
                             taxa_colors = NULL, keep_taxa = NULL) {
  
  t_abund <- trans_abund$new(
    dataset   = mt_obj,
    taxrank   = tax_level,
    ntaxa     = top_n,
    groupmean = group_var
  )
  
  plot_data <- t_abund$data_abund
  
  if (!is.null(keep_taxa)) {
    plot_data$Taxonomy <- ifelse(plot_data$Taxonomy %in% keep_taxa,
                                 plot_data$Taxonomy, "Others")
  } else {
    top_taxa <- plot_data %>%
      filter(Taxonomy != "Others") %>%
      group_by(Taxonomy) %>%
      summarise(mean_abund = mean(Abundance)) %>%
      arrange(desc(mean_abund)) %>%
      head(top_n) %>%
      pull(Taxonomy)
    plot_data$Taxonomy <- ifelse(plot_data$Taxonomy %in% top_taxa,
                                 plot_data$Taxonomy, "Others")
  }
  
  plot_data <- plot_data %>%
    group_by(Sample, Taxonomy) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop") %>%
    tidyr::complete(Sample, Taxonomy, fill = list(Abundance = 0))  # ★ 补零
  
  if (!is.null(taxa_colors)) {
    plot_data$Taxonomy <- factor(plot_data$Taxonomy, levels = names(taxa_colors))
    fill_scale <- scale_fill_manual(values = taxa_colors, drop = FALSE)
  } else {
    taxa_order <- c(rev(top_taxa), "Others")
    plot_data$Taxonomy <- factor(plot_data$Taxonomy, levels = taxa_order)
    local_colors <- setNames(pal_taxa[seq_along(top_taxa)], rev(top_taxa))
    local_colors["Others"] <- "grey80"
    fill_scale <- scale_fill_manual(values = local_colors)
  }
  
  p <- ggplot(plot_data, aes(x = Sample, y = Abundance, fill = Taxonomy)) +
    geom_bar(stat = "identity", position = "stack", width = 0.7) +
    fill_scale +
    scale_y_continuous(expand = c(0, 0), labels = function(x) paste0(x, "%")) +
    labs(x = NULL, y = "Relative Abundance (%)",
         title = title_text, fill = tax_level) +
    theme_half_open() +
    theme(legend.position = "right")
  
  return(p)
}
# ---- 调用：传入 keep_taxa = all_taxa ----
p1 <- plot_stacked_bar(mt_all, "Disease",   "Family", 8,
                       "Family — IBD vs HC",
                       taxa_colors = shared_colors, keep_taxa = all_taxa)
p2 <- plot_stacked_bar(mt_ibd, "dep_group", "Family", 8,
                       "Family — Depression",
                       taxa_colors = shared_colors, keep_taxa = all_taxa)
p3 <- plot_stacked_bar(mt_ibd, "anx_group", "Family", 8,
                       "Family — Anxiety",
                       taxa_colors = shared_colors, keep_taxa = all_taxa)

p_family_all <- (p1 + p2 + p3)  &
  theme(legend.position = "right")

ggsave("./results/abundance_family_stacked.pdf", p_family_all,
       width = 20, height = 8, device = cairo_pdf)

# --- Family 三分组 ---
p4 <- plot_stacked_bar(mt_ibd, "dep_3grp", "Family", 8,
                       "Family — Depression (3-group)",
                       taxa_colors = shared_colors, keep_taxa = all_taxa)
p5 <- plot_stacked_bar(mt_ibd, "anx_3grp", "Family", 8,
                       "Family — Anxiety (3-group)",
                       taxa_colors = shared_colors, keep_taxa = all_taxa)

p_family_3grp <- (p1 + p4 + p5) &
  theme(legend.position = "right")

ggsave("./results/abundance_family_stacked_3grp.pdf", p_family_3grp,
       width = 20, height = 8, device = cairo_pdf)

# =============================================================================
# 6. Genus层面堆积图
# =============================================================================


p_genus_disease <- plot_stacked_bar(
  mt_all, group_var = "Disease", tax_level = "Genus",
  top_n = 8, title_text = "Genus — IBD vs HC"
)

p_genus_dep <- plot_stacked_bar(
  mt_ibd, group_var = "dep_group", tax_level = "Genus",
  top_n = 8, title_text = "Genus — Depression"
)

p_genus_anx <- plot_stacked_bar(
  mt_ibd, group_var = "anx_group", tax_level = "Genus",
  top_n = 8, title_text = "Genus — Anxiety"
)

p_genus_all <- p_genus_disease + p_genus_dep + p_genus_anx 

print(p_genus_all)


# 直接用默认 pdf 设备保存
ggsave("./results/abundance_genus_stacked.pdf", p_genus_all,
       width = 20, height = 8, device = cairo_pdf)

# --- Genus 三分组 ---
p_genus_dep3 <- plot_stacked_bar(
  mt_ibd, group_var = "dep_3grp", tax_level = "Genus",
  top_n = 8, title_text = "Genus — Depression (3-group)"
)

p_genus_anx3 <- plot_stacked_bar(
  mt_ibd, group_var = "anx_3grp", tax_level = "Genus",
  top_n = 8, title_text = "Genus — Anxiety (3-group)"
)

p_genus_3grp <- p_genus_disease + p_genus_dep3 + p_genus_anx3
print(p_genus_3grp)

ggsave("./results/abundance_genus_stacked_3grp.pdf", p_genus_3grp,
       width = 16, height = 8, device = cairo_pdf)


# =============================================================================
# 7. 单样本堆积图（展示每个样本的组成，按组排列）
# =============================================================================

plot_sample_bar <- function(mt_obj, group_var, tax_level = "Phylum",
                            top_n = 10, title_text = "") {
  
  t_abund <- trans_abund$new(
    dataset   = mt_obj,
    taxrank   = tax_level,
    ntaxa     = top_n,
    groupmean = NULL
  )
  
  plot_data <- t_abund$data_abund
  
  # 归类 Others
  top_taxa <- plot_data %>%
    group_by(Taxonomy) %>%
    summarise(mean_abund = mean(Abundance)) %>%
    arrange(desc(mean_abund)) %>%
    head(top_n) %>%
    pull(Taxonomy)
  
  plot_data$Taxonomy <- ifelse(plot_data$Taxonomy %in% top_taxa,
                               plot_data$Taxonomy, "Others")
  
  plot_data <- plot_data %>%
    group_by(Sample, Taxonomy) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop")
  
  taxa_order <- c(rev(top_taxa), "Others")
  plot_data$Taxonomy <- factor(plot_data$Taxonomy, levels = taxa_order)
  
  # 配色
  taxa_colors <- setNames(pal_taxa[1:length(top_taxa)], rev(top_taxa))
  taxa_colors["Others"] <- "grey80"
  
  # 按组排序样本
  sample_order <- mt_obj$sample_table %>%
    arrange(.data[[group_var]]) %>%
    rownames()
  plot_data$Sample <- factor(plot_data$Sample, levels = sample_order)
  
  # facet 需要 group 信息，从 sample_table 合并进来
  group_info <- mt_obj$sample_table[, group_var, drop = FALSE]
  group_info$Sample <- rownames(group_info)
  plot_data <- left_join(plot_data, group_info, by = "Sample")
  
  p <- ggplot(plot_data, aes(x = Sample, y = Abundance, fill = Taxonomy)) +
    geom_bar(stat = "identity", position = "stack", width = 1) +
    scale_fill_manual(values = taxa_colors) +
    scale_y_continuous(expand = c(0, 0), labels = function(x) paste0(x, "%")) +
    facet_grid(as.formula(paste("~", group_var)), scales = "free_x", space = "free_x") +
    labs(x = NULL, y = "Relative Abundance (%)",
         title = title_text, fill = tax_level) +
    theme_half_open() +
    theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      strip.text   = element_text(size = 16, face = "bold"),
      legend.position = "right"
    )
  
  return(p)
}
# IBD vs HC 单样本
p_sample_disease <- plot_sample_bar(
  mt_all, group_var = "Disease", tax_level = "Phylum",
  top_n = 10, title_text = "Individual Composition — IBD vs HC"
)
print(p_sample_disease)
ggsave("abundance_phylum_individual_disease.pdf", p_sample_disease,
       width = 20, height = 7, device = cairo_pdf)

# Depression 单样本
p_sample_dep <- plot_sample_bar(
  mt_ibd, group_var = "dep_group", tax_level = "Phylum",
  top_n = 10, title_text = "Individual Composition — Depression"
)
ggsave("abundance_phylum_individual_dep.pdf", p_sample_dep,
       width = 18, height = 7, device = cairo_pdf)

# Anxiety 单样本
p_sample_anx <- plot_sample_bar(
  mt_ibd, group_var = "anx_group", tax_level = "Phylum",
  top_n = 10, title_text = "Individual Composition — Anxiety"
)
ggsave("abundance_phylum_individual_anx.pdf", p_sample_anx,
       width = 18, height = 7, device = cairo_pdf)

# --- 单样本堆积图 三分组 ---
p_sample_dep3 <- plot_sample_bar(
  mt_ibd, group_var = "dep_3grp", tax_level = "Phylum",
  top_n = 10, title_text = "Individual Composition — Depression (3-group)"
)
ggsave("abundance_phylum_individual_dep_3grp.pdf", p_sample_dep3,
       width = 20, height = 7, device = cairo_pdf)

p_sample_anx3 <- plot_sample_bar(
  mt_ibd, group_var = "anx_3grp", tax_level = "Phylum",
  top_n = 10, title_text = "Individual Composition — Anxiety (3-group)"
)
ggsave("abundance_phylum_individual_anx_3grp.pdf", p_sample_anx3,
       width = 20, height = 7, device = cairo_pdf)

cat("\n=== 03_abundance_analysis.R 完成 ===\n")



####nested####
# 安装 ggnested（如果没有）
library(ggnested)

# =========================
# 嵌套堆积图函数
# =========================
plot_nested_bar <- function(mt_obj, group_var, top_n = 20, fix_nsub = 3,
                            title_text = "", use_groupmean = TRUE) {
  
  t_abund <- trans_abund$new(
    dataset    = mt_obj,
    taxrank    = "Genus",
    ntaxa      = top_n,
    show       = 0,
    high_level = "Phylum",
    high_level_fix_nsub = fix_nsub,
    groupmean  = if (use_groupmean) group_var else NULL
  )
  
  p <- t_abund$plot_bar(
    ggnested        = TRUE,
    high_level_add_other = TRUE,
    xtext_angle     = 45
  )
  
  p <- p +
    labs(x = NULL, y = "Relative Abundance (%)", title = title_text) +
    theme_half_open() +
    theme(
      legend.position = "right",
      legend.text     = element_text(size = 10),
      legend.title    = element_text(size = 12, face = "bold")
    )
  
  return(p)
}

# =========================
# 三组对比：组均值嵌套堆积图
# =========================

# IBD vs HC
p_nested_disease <- plot_nested_bar(
  mt_all, group_var = "Disease",
  top_n = 20, fix_nsub = 3,
  title_text = "IBD vs HC"
)

# Depression
p_nested_dep <- plot_nested_bar(
  mt_ibd, group_var = "dep_group",
  top_n = 20, fix_nsub = 3,
  title_text = "Depression"
)

# Anxiety
p_nested_anx <- plot_nested_bar(
  mt_ibd, group_var = "anx_group",
  top_n = 20, fix_nsub = 3,
  title_text = "Anxiety"
)

print(p_nested_disease)
print(p_nested_dep)
print(p_nested_anx)

ggsave("nested_bar_disease.pdf", p_nested_disease,
       width = 14, height = 8, device = cairo_pdf)
ggsave("nested_bar_depression.pdf", p_nested_dep,
       width = 14, height = 8, device = cairo_pdf)


plot_nested_bar <- function(mt_obj, group_var, top_n = 30, fix_nsub = 5,
                            title_text = "", use_groupmean = TRUE) {
  
  t_abund <- trans_abund$new(
    dataset    = mt_obj,
    taxrank    = "Genus",
    ntaxa      = top_n,
    show       = 0,
    high_level = "Phylum",
    high_level_fix_nsub = fix_nsub,
    groupmean  = if (use_groupmean) group_var else NULL
  )
  
  p <- t_abund$plot_bar(
    ggnested        = TRUE,
    high_level_add_other = TRUE,
    xtext_angle     = 45
  )
  
  p <- p +
    labs(x = NULL, y = "Relative Abundance (%)", title = title_text) +
    theme_half_open() +
    theme(
      legend.position = "right",
      legend.text     = element_text(size = 10),
      legend.title    = element_text(size = 12, face = "bold")
    )
  
  return(p)
}

# IBD vs HC
p_nested_disease <- plot_nested_bar(
  mt_all, group_var = "Disease",
  top_n = 20, fix_nsub = 3,
  title_text = "IBD vs HC"
)

# Depression
p_nested_dep <- plot_nested_bar(
  mt_ibd, group_var = "dep_group",
  top_n = 20, fix_nsub = 3,
  title_text = "Depression"
)

# Anxiety
p_nested_anx <- plot_nested_bar(
  mt_ibd, group_var = "anx_group",
  top_n = 20, fix_nsub = 3,
  title_text = "Anxiety"
)

# Depression 三分组
p_nested_dep3 <- plot_nested_bar(
  mt_ibd, group_var = "dep_3grp",
  top_n = 20, fix_nsub = 3,
  title_text = "Depression (3-group)"
)

# Anxiety 三分组
p_nested_anx3 <- plot_nested_bar(
  mt_ibd, group_var = "anx_3grp",
  top_n = 20, fix_nsub = 3,
  title_text = "Anxiety (3-group)"
)

ggsave("nested_bar_depression_3grp.pdf", p_nested_dep3,
       width = 16, height = 8, device = cairo_pdf)
ggsave("nested_bar_anxiety_3grp.pdf", p_nested_anx3,
       width = 16, height = 8, device = cairo_pdf)


# 待整理

# ==============================================================================
# 03_alpha_diversity_plot.R — Alpha 多样性箱线图 (Shannon + Simpson)
#
# Part A: PHQ9/GAD7 二分组 (≥10 vs <10)
# Part B: PHQ9/GAD7 三分组 (0-4 / 5-9 / ≥10)
#
# 输出:
#   results/02_diversity/alpha_diversity_dep_anx.pdf/png        (二分组)
#   results/02_diversity/alpha_diversity_dep_anx_3grp.pdf/png   (三分组)
# ==============================================================================

source("R/00_setup.R")

library(ggpubr)
library(microeco)

# =============================================================================
# 1. 数据导入与 Alpha 多样性计算
# =============================================================================

mt_ibd <- readRDS("./microbiome/mt_ibd_282.rds")

mt1 <- clone(mt_ibd)
mt1$cal_abund(select_cols = 1:7, rel = FALSE)
mt1$tidy_dataset()

# PHQ9 二分组 (Depression)
t1 <- trans_alpha$new(dataset = mt1, group = "PHQ9_bi10")
t1$cal_diff(method = "wilcox")

# GAD7 二分组 (Anxiety)
t2 <- trans_alpha$new(dataset = mt1, group = "GAD7_bi10")
t2$cal_diff(method = "wilcox")

cat("=== Alpha diff (Depression) ===\n")
print(head(t1$res_diff))
cat("\n=== Alpha diff (Anxiety) ===\n")
print(head(t2$res_diff))

# =============================================================================
# 2. 数据整理：合并 Depression + Anxiety，统一列名
# =============================================================================

# Depression
df_dep <- subset(t1$data_alpha, Measure %in% c("Shannon", "Simpson"))
df_dep$PHQ9_bi10 <- factor(df_dep$PHQ9_bi10,
                           levels = c("0", "1"),
                           labels = c("Non-depressed", "Depressed"))
df_dep$Group <- df_dep$PHQ9_bi10
df_dep$Group_type <- "Depression"

# Anxiety
df_anx <- subset(t2$data_alpha, Measure %in% c("Shannon", "Simpson"))
df_anx$GAD7_bi10 <- factor(df_anx$GAD7_bi10,
                           levels = c("0", "1"),
                           labels = c("Non-anxious", "Anxious"))
df_anx$Group <- df_anx$GAD7_bi10
df_anx$Group_type <- "Anxiety"
pacman::p_unload(MASS)
df_all <- bind_rows(
  df_dep %>% select(Measure, Value, Group, Group_type),
  df_anx %>% select(Measure, Value, Group, Group_type)
)
df_all$Group_type <- factor(df_all$Group_type, levels = c("Depression", "Anxiety"))
df_all$Measure    <- factor(df_all$Measure, levels = c("Shannon", "Simpson"))

df_shannon <- subset(df_all, Measure == "Shannon")
df_simpson <- subset(df_all, Measure == "Simpson")

# =============================================================================
# 3. 绘图主题（半开放轴线）
# =============================================================================

theme_half_open <- function() {
  theme_classic(base_size = 18) +
    theme(
      panel.background  = element_rect(fill = "white", color = NA),
      plot.background   = element_rect(fill = "white", color = NA),
      panel.border      = element_blank(),
      axis.line.x       = element_line(color = "black", linewidth = 1.0),
      axis.line.y       = element_line(color = "black", linewidth = 1.0),
      axis.ticks        = element_line(color = "black", linewidth = 1.0),
      axis.ticks.length = unit(0.2, "cm"),
      axis.text         = element_text(color = "black", size = 16),
      axis.text.x       = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.title        = element_text(color = "black", size = 20, face = "bold"),
      plot.title        = element_text(size = 20, face = "bold", hjust = 0.5),
      legend.position   = "none",
      strip.background  = element_blank(),
      strip.text        = element_text(size = 16, face = "bold", color = "black"),
      panel.grid        = element_blank()
    )
}

# =============================================================================
# 4. 配色与统计设置
# =============================================================================

fill_cols <- c(
  "Non-depressed" = "#2F9B85",
  "Depressed"     = "#D95F5F",
  "Non-anxious"   = "#6FA8DC",
  "Anxious"       = "#E6862D"
)

comp_dep <- list(c("Non-depressed", "Depressed"))
comp_anx <- list(c("Non-anxious", "Anxious"))

symnum_args <- list(
  cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
  symbols   = c("***", "**", "*", "#", "ns")
)

# =============================================================================
# 5. 构建单个面板的辅助函数
# =============================================================================

make_alpha_panel <- function(df, y_label, title_label) {
  ggplot(df, aes(x = Group, y = Value, fill = Group)) +
    geom_violin(trim = FALSE, color = NA, alpha = 0.28, width = 0.9) +
    geom_boxplot(width = 0.22, linewidth = 1.2, color = "black",
                 outlier.shape = NA, alpha = 1) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 2,
                 fill = "white", color = "black", stroke = 1.4) +
    stat_compare_means(
      comparisons = comp_dep,
      data = subset(df, Group_type == "Depression"),
      label = "p.signif", size = 5, method = "wilcox.test",
      symnum.args = symnum_args
    ) +
    stat_compare_means(
      comparisons = comp_anx,
      data = subset(df, Group_type == "Anxiety"),
      label = "p.signif", size = 5, method = "wilcox.test",
      symnum.args = symnum_args
    ) +
    facet_wrap(~ Group_type, scales = "free_x") +
    scale_fill_manual(values = fill_cols) +
    labs(x = NULL, y = y_label, title = title_label) +
    theme_half_open()
}

# =============================================================================
# 6. 绘图与拼接
# =============================================================================

p_shannon <- make_alpha_panel(df_shannon, "Shannon Index", "Shannon")
p_simpson <- make_alpha_panel(df_simpson, "Simpson Index", "Simpson")

p_final <- p_shannon + p_simpson +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "A")

print(p_final)

# =============================================================================
# 7. 保存
# =============================================================================

save_plot(p_final, "alpha_diversity_dep_anx", "02_diversity", width = 14, height = 6)

# =============================================================================
# ==================== Part B: 三分组 (0-4 / 5-9 / ≥10) ======================
# =============================================================================

# =============================================================================
# 8. 创建三分组变量并计算 Alpha 多样性
# =============================================================================

mt2 <- clone(mt_ibd)

mt2$sample_table$phq9_3cat <- cut(
  mt2$sample_table$PHQ9_sum,
  breaks = c(-1, 4, 9, Inf),
  labels = c("None (0-4)", "Mild (5-9)", "Mod-Severe (>=10)")
)

mt2$sample_table$gad7_3cat <- cut(
  mt2$sample_table$GAD7_sum,
  breaks = c(-1, 4, 9, Inf),
  labels = c("None (0-4)", "Mild (5-9)", "Mod-Severe (>=10)")
)

mt2$cal_abund(select_cols = 1:7, rel = FALSE)
mt2$tidy_dataset()

# PHQ9 三分组
t3 <- trans_alpha$new(dataset = mt2, group = "phq9_3cat")
t3$cal_diff(method = "KW")

# GAD7 三分组
t4 <- trans_alpha$new(dataset = mt2, group = "gad7_3cat")
t4$cal_diff(method = "KW")

cat("\n=== Alpha diff — PHQ9 三分组 (KW) ===\n")
print(head(t3$res_diff))
cat("\n=== Alpha diff — GAD7 三分组 (KW) ===\n")
print(head(t4$res_diff))

# =============================================================================
# 9. 三分组数据整理
# =============================================================================

grp_levels <- c("None (0-4)", "Mild (5-9)", "Mod-Severe (>=10)")

# Depression 三分组
df_dep3 <- subset(t3$data_alpha, Measure %in% c("Shannon", "Simpson"))
df_dep3$Group <- factor(df_dep3$phq9_3cat, levels = grp_levels)
df_dep3$Group_type <- "PHQ-9 (Depression)"

# Anxiety 三分组
df_anx3 <- subset(t4$data_alpha, Measure %in% c("Shannon", "Simpson"))
df_anx3$Group <- factor(df_anx3$gad7_3cat, levels = grp_levels)
df_anx3$Group_type <- "GAD-7 (Anxiety)"

df_all3 <- bind_rows(
  df_dep3 %>% select(Measure, Value, Group, Group_type),
  df_anx3 %>% select(Measure, Value, Group, Group_type)
)
df_all3$Group_type <- factor(df_all3$Group_type,
                             levels = c("PHQ-9 (Depression)", "GAD-7 (Anxiety)"))
df_all3$Measure <- factor(df_all3$Measure, levels = c("Shannon", "Simpson"))

df_shannon3 <- subset(df_all3, Measure == "Shannon")
df_simpson3 <- subset(df_all3, Measure == "Simpson")

# =============================================================================
# 10. 三分组配色与统计设置
# =============================================================================

fill_cols3 <- c(
  "None (0-4)"        = "#2F9B85",
  "Mild (5-9)"        = "#F5C242",
  "Mod-Severe (>=10)" = "#D95F5F"
)

comp_3grp <- list(
  c("None (0-4)", "Mild (5-9)"),
  c("Mild (5-9)", "Mod-Severe (>=10)"),
  c("None (0-4)", "Mod-Severe (>=10)")
)

# =============================================================================
# 11. 三分组绘图辅助函数
# =============================================================================

make_alpha_panel_3grp <- function(df, y_label, title_label) {
  ggplot(df, aes(x = Group, y = Value, fill = Group)) +
    geom_violin(trim = FALSE, color = NA, alpha = 0.28, width = 0.9) +
    geom_boxplot(width = 0.22, linewidth = 1.2, color = "black",
                 outlier.shape = NA, alpha = 1) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 2,
                 fill = "white", color = "black", stroke = 1.4) +
    stat_compare_means(
      comparisons = comp_3grp,
      label = "p.signif", size = 4, method = "wilcox.test",
      symnum.args = symnum_args,
      step.increase = 0.08
    ) +
    facet_wrap(~ Group_type, scales = "free_x") +
    scale_fill_manual(values = fill_cols3) +
    labs(x = NULL, y = y_label, title = title_label) +
    theme_half_open()
}

# =============================================================================
# 12. 三分组绘图与拼接
# =============================================================================

p_shannon3 <- make_alpha_panel_3grp(df_shannon3, "Shannon Index", "Shannon")
p_simpson3 <- make_alpha_panel_3grp(df_simpson3, "Simpson Index", "Simpson")

p_final3 <- p_shannon3 / p_simpson3 +
  plot_annotation(tag_levels = "A")

print(p_final3)

# =============================================================================
# 13. 保存三分组图
# =============================================================================

save_plot(p_final3, "alpha_diversity_dep_anx_3grp", "02_diversity", width = 8, height = 12)

cat("\n=== 03_alpha_diversity_plot.R 完成 ===\n")




library(ggplot2)
library(vegan)
library(dplyr)
library(patchwork)
library(ggpubr)

perm_dep <- load_intermediate("permanova_results_dep.rds")
perm_anx <- load_intermediate("permanova_results_anx.rds")

dm <- load_intermediate("distance_matrices.rds")
bc_complete <- dm$bc_complete
ait_dist    <- dm$ait_dist
cohort      <- dm$cohort

# 从 tibble 中提取指定距离+变量的 R² 和 P 值
extract_perm <- function(tbl, dist, var) {
  row <- tbl %>% filter(Distance == dist, Variable == var)
  list(r2 = row$R2, pval = row$P)
}


# =========================
# 半封闭主题
# =========================
theme_half_open <- function() {
  theme_classic(base_size = 18) +
    theme(
      panel.background  = element_rect(fill = "white", color = NA),
      plot.background   = element_rect(fill = "white", color = NA),
      panel.border      = element_blank(),
      axis.line.x       = element_line(color = "black", linewidth = 1.0),
      axis.line.y       = element_line(color = "black", linewidth = 1.0),
      axis.ticks        = element_line(color = "black", linewidth = 1.0),
      axis.ticks.length = unit(0.2, "cm"),
      axis.text         = element_text(color = "black", size = 16),
      axis.title        = element_text(color = "black", size = 20, face = "bold"),
      plot.title        = element_text(size = 20, face = "bold", hjust = 0.5),
      plot.subtitle     = element_text(size = 14, hjust = 0.5, color = "grey30"),
      legend.title      = element_text(size = 16, face = "bold"),
      legend.text       = element_text(size = 15),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key        = element_rect(fill = "white", color = NA),
      panel.grid        = element_blank()
    )
}

# =========================
# 配色（与之前统一）
# =========================
colors_dep <- c("Non-depressed" = "#2F9B85", "Depressed" = "#D95F5F")
colors_anx <- c("Non-anxious" = "#6FA8DC", "Anxious" = "#E6862D")

# =========================
# PCoA 绘图函数
# =========================
plot_pcoa <- function(dist_obj, sample_data, group_var, group_labels,
                      fill_colors, title_text, r2, pval) {

  pcoa_res <- cmdscale(dist_obj, k = 2, eig = TRUE)
  eig_pct <- round(100 * pcoa_res$eig[1:2] / sum(pcoa_res$eig[pcoa_res$eig > 0]), 1)

  df <- tibble(
    PC1   = pcoa_res$points[, 1],
    PC2   = pcoa_res$points[, 2],
    Group = factor(sample_data[[group_var]],
                   levels = names(group_labels),
                   labels = group_labels)
  )

  p_label <- ifelse(pval < 0.001, "P < 0.001", sprintf("P = %.3f", pval))

  ggplot(df, aes(x = PC1, y = PC2, color = Group, fill = Group)) +
    geom_point(size = 2.5, alpha = 0.8) +
    stat_ellipse(level = 0.95, geom = "polygon", alpha = 0.1,
                 linewidth = 1, linetype = "solid") +
    scale_color_manual(values = fill_colors) +
    scale_fill_manual(values = fill_colors) +
    labs(
      x = sprintf("PCoA1 (%.1f%%)", eig_pct[1]),
      y = sprintf("PCoA2 (%.1f%%)", eig_pct[2]),
      title = title_text,
      subtitle = sprintf("PERMANOVA: R² = %.4f, %s", r2, p_label)
    ) +
    theme_half_open() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.0),
      axis.line = element_blank(),
      legend.position = "right"
    )
}

# =========================
# Beta dispersion 箱线图（distance to centroid）
# =========================
plot_beta_boxplot <- function(dist_obj, sample_data, group_var, group_labels,
                              fill_colors, title_text) {

  groups <- factor(sample_data[[group_var]],
                   levels = names(group_labels),
                   labels = group_labels)
  names(groups) <- rownames(sample_data)

  shared <- intersect(names(groups), labels(dist_obj))
  dist_sub <- as.dist(as.matrix(dist_obj)[shared, shared])
  groups   <- groups[shared]

  bd <- betadisper(dist_sub, groups)

  df <- data.frame(
    Group    = bd$group,
    Distance = bd$distances
  )

  comps <- combn(levels(df$Group), 2, simplify = FALSE)

  ggplot(df, aes(x = Group, y = Distance, fill = Group)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.5, width = 0.6) +
    stat_compare_means(comparisons = comps, method = "wilcox.test",
                       label = "p.signif", size = 4, step.increase = 0.08) +
    scale_fill_manual(values = fill_colors) +
    labs(x = NULL, y = "Distance to centroid", title = title_text) +
    theme_half_open() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.0),
      axis.line = element_blank(),
      legend.position = "none",
      axis.text.x = element_text(size = 13, angle = 45, hjust = 1, vjust = 1)
    )
}

# =========================
# PERMANOVA 结果提取
# =========================
pm_dep_bc  <- extract_perm(perm_dep, "Bray-Curtis", "PHQ9_sum")
pm_dep_ait <- extract_perm(perm_dep, "Aitchison", "PHQ9_sum")
pm_anx_bc  <- extract_perm(perm_anx, "Bray-Curtis", "GAD7_sum")
pm_anx_ait <- extract_perm(perm_anx, "Aitchison", "GAD7_sum")

# =========================
# 两组标签 & 配色
# =========================
dep_labels <- c("0" = "Non-depressed", "1" = "Depressed")
anx_labels <- c("0" = "Non-anxious", "1" = "Anxious")

# =========================
# 三组标签 & 配色 (severity)
# PHQ-9: 0-4 None, 5-9 Mild, >=10 Moderate+
# GAD-7: 0-4 None, 5-9 Mild, >=10 Moderate+
# =========================
cohort$dep_sev3 <- cut(cohort$PHQ9_sum,
                       breaks = c(-Inf, 4, 9, Inf),
                       labels = c("None", "Mild", "Moderate+"))
cohort$anx_sev3 <- cut(cohort$GAD7_sum,
                       breaks = c(-Inf, 4, 9, Inf),
                       labels = c("None", "Mild", "Moderate+"))

dep3_labels <- c("None" = "None", "Mild" = "Mild", "Moderate+" = "Moderate+")
anx3_labels <- c("None" = "None", "Mild" = "Mild", "Moderate+" = "Moderate+")

colors_dep3 <- c("None" = "#2F9B85", "Mild" = "#F5A623", "Moderate+" = "#D95F5F")
colors_anx3 <- c("None" = "#6FA8DC", "Mild" = "#F5A623", "Moderate+" = "#E6862D")

# =====================================================================
# A. 两组: PCoA + Within/Between boxplot
# =====================================================================

# --- Depression (Bray-Curtis) ---
p_pcoa_dep_bc <- plot_pcoa(bc_complete, cohort, "PHQ9_bi10", dep_labels,
                           colors_dep, "Bray-Curtis — Depression",
                           pm_dep_bc$r2, pm_dep_bc$pval)
p_box_dep_bc  <- plot_beta_boxplot(bc_complete, cohort, "PHQ9_bi10", dep_labels,
                                       colors_dep, "Beta Dispersion")

# --- Depression (Aitchison) ---
p_pcoa_dep_ait <- plot_pcoa(ait_dist, cohort, "PHQ9_bi10", dep_labels,
                            colors_dep, "Aitchison — Depression",
                            pm_dep_ait$r2, pm_dep_ait$pval)
p_box_dep_ait  <- plot_beta_boxplot(ait_dist, cohort, "PHQ9_bi10", dep_labels,
                                        colors_dep, "Beta Dispersion")

# --- Anxiety (Bray-Curtis) ---
p_pcoa_anx_bc <- plot_pcoa(bc_complete, cohort, "GAD7_bi10", anx_labels,
                           colors_anx, "Bray-Curtis — Anxiety",
                           pm_anx_bc$r2, pm_anx_bc$pval)
p_box_anx_bc  <- plot_beta_boxplot(bc_complete, cohort, "GAD7_bi10", anx_labels,
                                       colors_anx, "Beta Dispersion")

# --- Anxiety (Aitchison) ---
p_pcoa_anx_ait <- plot_pcoa(ait_dist, cohort, "GAD7_bi10", anx_labels,
                            colors_anx, "Aitchison — Anxiety",
                            pm_anx_ait$r2, pm_anx_ait$pval)
p_box_anx_ait  <- plot_beta_boxplot(ait_dist, cohort, "GAD7_bi10", anx_labels,
                                        colors_anx, "Beta Dispersion")

# 拼接: PCoA(宽) + boxplot(窄)
p_2grp <- (p_pcoa_dep_bc | p_box_dep_bc | p_pcoa_dep_ait | p_box_dep_ait) /
          (p_pcoa_anx_bc | p_box_anx_bc | p_pcoa_anx_ait | p_box_anx_ait) +
  plot_layout(widths = rep(c(3, 1.5), 2)) +
  plot_annotation(tag_levels = "A",
                  title = "PCoA & Within/Between Distance — Binary Groups",
                  theme = theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5)))

print(p_2grp)
ggsave("results/group_meeting/add/beta_pcoa_2group.pdf", p_2grp,
       width = 28, height = 14)

# =====================================================================
# B. 三组 (severity): PCoA + Within/Between boxplot
# =====================================================================

# --- Depression severity (Bray-Curtis) ---
p_pcoa_dep3_bc <- plot_pcoa(bc_complete, cohort, "dep_sev3", dep3_labels,
                            colors_dep3, "Bray-Curtis — Depression Severity",
                            pm_dep_bc$r2, pm_dep_bc$pval)
p_box_dep3_bc  <- plot_beta_boxplot(bc_complete, cohort, "dep_sev3", dep3_labels,
                                        colors_dep3, "Beta Dispersion")

# --- Depression severity (Aitchison) ---
p_pcoa_dep3_ait <- plot_pcoa(ait_dist, cohort, "dep_sev3", dep3_labels,
                             colors_dep3, "Aitchison — Depression Severity",
                             pm_dep_ait$r2, pm_dep_ait$pval)
p_box_dep3_ait  <- plot_beta_boxplot(ait_dist, cohort, "dep_sev3", dep3_labels,
                                         colors_dep3, "Beta Dispersion")

# --- Anxiety severity (Bray-Curtis) ---
p_pcoa_anx3_bc <- plot_pcoa(bc_complete, cohort, "anx_sev3", anx3_labels,
                            colors_anx3, "Bray-Curtis — Anxiety Severity",
                            pm_anx_bc$r2, pm_anx_bc$pval)
p_box_anx3_bc  <- plot_beta_boxplot(bc_complete, cohort, "anx_sev3", anx3_labels,
                                        colors_anx3, "Beta Dispersion")

# --- Anxiety severity (Aitchison) ---
p_pcoa_anx3_ait <- plot_pcoa(ait_dist, cohort, "anx_sev3", anx3_labels,
                             colors_anx3, "Aitchison — Anxiety Severity",
                             pm_anx_ait$r2, pm_anx_ait$pval)
p_box_anx3_ait  <- plot_beta_boxplot(ait_dist, cohort, "anx_sev3", anx3_labels,
                                         colors_anx3, "Beta Dispersion")

p_3grp <- (p_pcoa_dep3_bc | p_box_dep3_bc | p_pcoa_dep3_ait | p_box_dep3_ait) /
          (p_pcoa_anx3_bc | p_box_anx3_bc | p_pcoa_anx3_ait | p_box_anx3_ait) +
  plot_layout(widths = rep(c(3, 1.5), 2)) +
  plot_annotation(tag_levels = "A",
                  title = "PCoA & Within/Between Distance — Severity Groups",
                  theme = theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5)))

print(p_3grp)
ggsave("results/group_meeting/add/beta_pcoa_3group_severity.pdf", p_3grp,
       width = 28, height = 14)

###R2####
library(tidyverse)

# --- 加载数据 ---
permanova_dep <- load_intermediate("permanova_results_dep.rds") %>%
  mutate(Model = "Depression (PHQ-9)") %>%
  filter(Variable != "GAD7_sum")

permanova_anx <- load_intermediate("permanova_results_anx.rds") %>%
  mutate(Model = "Anxiety (GAD-7)") %>%
  filter(Variable != "PHQ9_sum")

df_plot <- bind_rows(permanova_dep, permanova_anx) %>%
  mutate(
    Variable_label = case_when(
      Variable == "PHQ9_sum" ~ "PHQ-9",
      Variable == "GAD7_sum" ~ "GAD-7",
      Variable == "site"     ~ "Study site",
      Variable == "Age"      ~ "Age",
      Variable == "Sex"      ~ "Sex",
      Variable == "BMI"      ~ "BMI",
      Variable == "hPDI"     ~ "hPDI",
      TRUE ~ Variable
    ),
    Sig = case_when(
      P < 0.001 ~ "***",
      P < 0.01  ~ "**",
      P < 0.05  ~ "*",
      TRUE      ~ ""
    ),
    R2_pct = R2 * 100
  )

# Distance 顺序：Bray-Curtis 在左
df_plot$Distance <- factor(df_plot$Distance, levels = c("Bray-Curtis", "Aitchison"))
# Model 顺序
df_plot$Model <- factor(df_plot$Model, levels = c("Anxiety (GAD-7)", "Depression (PHQ-9)"))

# --- 按 Bray-Curtis 的 R2 排序 y 轴（每个 Model 独立排序）---
# 取 Bray-Curtis 的 R2 作为排序依据
bc_order <- df_plot %>%
  filter(Distance == "Bray-Curtis") %>%
  dplyr::select(Model, Variable_label, R2_pct)

# 创建排序用的 interaction variable，保证每个 facet 行独立
df_plot <- df_plot %>%
  left_join(bc_order %>% rename(R2_bc = R2_pct),
            by = c("Model", "Variable_label")) %>%
  mutate(facet_var = paste(Model, Variable_label, sep = "___"))

# 按 Model + Bray-Curtis R2 排序
ordered_levels <- df_plot %>%
  distinct(Model, Variable_label, facet_var, R2_bc) %>%
  arrange(Model, R2_bc) %>%
  pull(facet_var)

df_plot$facet_var <- factor(df_plot$facet_var, levels = ordered_levels)

# --- Nature 配色 ---
pal_nature <- c(
  "GAD-7"       = "#E64B35",
  "PHQ-9"       = "#E64B35",
  "Study site"  = "#8491B4",
  "Age"         = "#4DBBD5",
  "Sex"         = "#00A087",
  "BMI"         = "#3C5488",
  "hPDI"        = "#F39B7F"
)

# --- 画图 ---
p <- ggplot(df_plot, aes(x = R2_pct, y = facet_var, fill = Variable_label)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(aes(label = Sig), hjust = -0.3, size = 4, fontface = "bold") +
  facet_grid(Model ~ Distance, scales = "free") +
  scale_fill_manual(values = pal_nature) +
  scale_y_discrete(labels = function(x) sub("^.*___", "", x)) +   # 去掉 Model 前缀，只显示变量名
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.15)),
    labels = function(x) paste0(x, "%")
  ) +
  labs(
    x = expression("Explained variance (" * R^2 * ", %)"),
    y = NULL,
    title = "Proportion of variance in microbial composition\nexplained by host characteristics",
    subtitle = "PERMANOVA with marginal effects (9,999 permutations)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title         = element_text(face = "bold", size = 14, hjust = 0),
    plot.subtitle      = element_text(size = 11, color = "grey40", hjust = 0),
    strip.background   = element_rect(fill = "grey95", color = "grey70"),
    strip.text         = element_text(face = "bold", size = 11),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.y        = element_text(size = 11),
    plot.margin        = margin(10, 15, 10, 10)
  )

print(p)


# ======================================================================
# ============  microeco 官方接口重做 beta 多样性分析  ==================
# ======================================================================
library(microeco)
library(ggplot2)
library(patchwork)

# --- 数据导入 ---
mt <- load_intermediate("mt_ibd.rds")
mt$tidy_dataset()

# --- 配色 ---
colors_dep2 <- c("Non-depressed" = "#2F9B85", "Depressed" = "#D95F5F")
colors_anx2 <- c("Non-anxious" = "#6FA8DC", "Anxious" = "#E6862D")
colors_dep3 <- c("None" = "#2F9B85", "Mild" = "#F5A623", "Moderate+" = "#D95F5F")
colors_anx3 <- c("None" = "#6FA8DC", "Mild" = "#F5A623", "Moderate+" = "#E6862D")

# --- 创建分组变量 ---
# 两组 (binary)
mt$sample_table$depression <- ifelse(mt$sample_table$PHQ9_bi10 == 1,
                                     "Depressed", "Non-depressed")
mt$sample_table$anxiety    <- ifelse(mt$sample_table$GAD7_bi10 == 1,
                                     "Anxious", "Non-anxious")
# 三组 (severity)
mt$sample_table$dep_sev3 <- cut(mt$sample_table$PHQ9_sum,
                                breaks = c(-Inf, 4, 9, Inf),
                                labels = c("None", "Mild", "Moderate+"))
mt$sample_table$anx_sev3 <- cut(mt$sample_table$GAD7_sum,
                                breaks = c(-Inf, 4, 9, Inf),
                                labels = c("None", "Mild", "Moderate+"))

# --- 计算 Bray-Curtis 距离 ---
mt$cal_betadiv(unifrac = FALSE)

# =====================================================================
# 辅助函数: 对指定分组跑 PCoA + PERMANOVA + group distance boxplot
# =====================================================================
run_beta_microeco <- function(mt_obj, group_col, colors, title_prefix) {

  tb <- trans_beta$new(dataset = mt_obj, group = group_col, measure = "bray")

  # --- PERMANOVA ---
  tb$cal_manova(manova_all = TRUE, p_adjust_method = "fdr")
  cat("\n===", title_prefix, "PERMANOVA ===\n")
  print(tb$res_manova)

  perm_row <- tb$res_manova[group_col, ]
  r2_val   <- perm_row$R2
  p_val    <- perm_row$`Pr(>F)`
  p_text   <- if (p_val < 0.001) "P < 0.001" else sprintf("P = %.3f", p_val)
  anno_lab <- sprintf("PERMANOVA: R² = %.4f, %s", r2_val, p_text)

  # --- PCoA + PERMANOVA 注释 ---
  tb$cal_ordination(ordination = "PCoA")
  p_pcoa <- tb$plot_ordination(plot_color = group_col, plot_shape = group_col,
                               plot_type = c("point", "ellipse")) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    annotate("text", x = Inf, y = Inf, label = anno_lab,
             hjust = 1.05, vjust = 1.5, size = 4.2, fontface = "italic") +
    labs(title = paste(title_prefix, "— PCoA (Bray-Curtis)")) +
    theme_classic(base_size = 14) +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.0),
      axis.line = element_blank(),
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      legend.position = "right"
    )

  # --- Beta dispersion ---
  tb$cal_betadisper()
  cat("\n===", title_prefix, "Betadisper ===\n")
  print(tb$res_betadisper)

  # --- Group distance boxplot + 组间显著性 ---
  tb$cal_group_distance(within_group = TRUE)

  grp_lvls <- levels(factor(mt_obj$sample_table[[group_col]]))
  comps    <- combn(grp_lvls, 2, simplify = FALSE)

  p_dist <- tb$plot_group_distance(distance_pair_stat = FALSE) +
    stat_compare_means(comparisons = comps, method = "wilcox.test",
                       label = "p.signif", size = 4.5, step.increase = 0.08) +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    labs(title = paste(title_prefix, "— Bray-Curtis Distance"),
         y = "Bray-Curtis distance") +
    theme_classic(base_size = 14) +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.0),
      axis.line = element_blank(),
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
      legend.position = "none"
    )

  list(pcoa = p_pcoa, dist = p_dist, tb = tb)
}

# =====================================================================
# 1. 两组比较
# =====================================================================
res_dep2 <- run_beta_microeco(mt, "depression", colors_dep2, "Depression")
res_anx2 <- run_beta_microeco(mt, "anxiety",    colors_anx2, "Anxiety")

p_2grp_me <- (res_dep2$pcoa | res_dep2$dist) /
             (res_anx2$pcoa | res_anx2$dist) +
  plot_layout(widths = c(3, 1.5)) +
  plot_annotation(tag_levels = "A")

print(p_2grp_me)
ggsave("results/group_meeting/add/beta_microeco_2group.pdf", p_2grp_me,
       width = 16, height = 14)

# =====================================================================
# 2. 三组比较 (severity)
# =====================================================================
res_dep3 <- run_beta_microeco(mt, "dep_sev3", colors_dep3, "Depression Severity")
res_anx3 <- run_beta_microeco(mt, "anx_sev3", colors_anx3, "Anxiety Severity")

p_3grp_me <- (res_dep3$pcoa | res_dep3$dist) /
             (res_anx3$pcoa | res_anx3$dist) +
  plot_layout(widths = c(3, 1.5)) +
  plot_annotation(tag_levels = "A")

print(p_3grp_me)
ggsave("results/group_meeting/add/beta_microeco_3group_severity.pdf", p_3grp_me,
       width = 16, height = 14)
