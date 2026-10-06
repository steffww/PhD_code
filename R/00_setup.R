# ==============================================================================
# 00_setup.R — 全局设置：色卡、R包、辅助函数
# ==============================================================================

# --- R packages ---------------------------------------------------------------

suppressPackageStartupMessages({
  # 数据处理
  library(tidyverse)
  library(openxlsx)

  # 微生物组分析 (核心)
  library(vegan)          # PERMANOVA, alpha diversity
  library(microeco)       # 微生物组容器
  # Bioconductor 包: 按需加载 (需 BiocManager::install)
  # library(Maaslin2)     # 差异丰度 — 在 04_differential_abundance.R 中加载
  # library(ANCOMBC)      # 差异丰度 — 在 04_differential_abundance.R 中加载
  library(phyloseq)     # 微生物组容器 — 按需加载

  # 统计
  library(lme4)           # 混合效应模型 (备用)
  library(mediation)      # 中介分析
  library(lavaan)         # SEM

  # 可视化
  library(ggplot2)
  library(patchwork)      # 拼图
  library(ggrepel)        # 标签避让
  library(pheatmap)       # 热图
  library(forestploter)   # 森林图
  library(RColorBrewer)
  library(scales)
})

# 安全加载 Bioconductor 包
load_bioc <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Package '", pkg, "' not installed. Install with: BiocManager::install('", pkg, "')")
    return(FALSE)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  TRUE
}

# --- 色卡定义 (来自 色卡.png) ------------------------------------------------

pal_main <- c(
  dark_green  = "#2e3d28",
  orange      = "#ea7f2d",
  red         = "#d62f2f",
  lime        = "#87ae41",
  purple      = "#694898",
  lavender    = "#9c89bb",
  pink        = "#f09393",
  yellow      = "#cec73b",
  brown       = "#7d3326",
  blue        = "#286a9e"
)

# 常用分组配色
col_cd_uc  <- c(CD = "#ea7f2d", UC = "#286a9e", HC = "#87ae41")
col_dep    <- c(No = "#87ae41", Yes = "#d62f2f")         # depression Yes/No
col_site   <- c(Zhejiang = "#286a9e", Hunan = "#ea7f2d")
col_sex    <- c(male = "#286a9e", female = "#d62f2f")

# 连续色阶 (PHQ9 / GAD7)
col_gradient_low  <- "#87ae41"
col_gradient_high <- "#d62f2f"

# --- 全局绘图主题 ------------------------------------------------------------

theme_pub <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", colour = NA),
    legend.background = element_rect(fill = NA),
    plot.title = element_text(face = "bold", size = 13),
    axis.title = element_text(size = 11),
    legend.title = element_text(size = 10),
    legend.text  = element_text(size = 9)
  )
theme_set(theme_pub)

# --- 路径常量 -----------------------------------------------------------------

dir_data    <- "data"
dir_clean   <- "data/cleaned"
dir_inter   <- "intermediate"
dir_results <- "results"

# --- 辅助函数 -----------------------------------------------------------------

# 安全地加载中间数据 (如果存在就读，不存在就返回 NULL)
load_intermediate <- function(fname) {
  fpath <- file.path(dir_inter, fname)
  if (file.exists(fpath)) readRDS(fpath) else NULL
}

# 保存中间数据
save_intermediate <- function(obj, fname) {
  saveRDS(obj, file.path(dir_inter, fname))
  message("Saved: ", file.path(dir_inter, fname))
}

# 保存图片 (同时输出 PDF 和 PNG)
save_plot <- function(p, name, subdir, width = 8, height = 6, dpi = 300) {
  dir.create(file.path(dir_results, subdir), showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(dir_results, subdir, paste0(name, ".pdf")),
         p, width = width, height = height)
  ggsave(file.path(dir_results, subdir, paste0(name, ".png")),
         p, width = width, height = height, dpi = dpi)
  message("Saved: ", file.path(dir_results, subdir, name), " (.pdf + .png)")
}

# 标准化 "From" 到二分类 site
encode_site <- function(from_vec) {
  ifelse(from_vec %in% c("SYF", "XS"), "Zhejiang", "Hunan")
}

# 标准化 Region (合并大小写)
clean_region <- function(region_vec) {
  tolower(region_vec)
}

# BH 校正后提取显著项
extract_sig <- function(df, p_col = "pval", q_col = "qval", q_thresh = 0.1) {
  df[[q_col]] <- p.adjust(df[[p_col]], method = "BH")
  df[df[[q_col]] < q_thresh & !is.na(df[[q_col]]), ]
}

message("00_setup.R loaded successfully.")
