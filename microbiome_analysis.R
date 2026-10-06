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
