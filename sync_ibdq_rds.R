#!/usr/bin/env Rscript
# Synchronize the new baseline IBDQ domains into microeco sample_table and
# faecal RDS meta_enriched. Model-only meta and numerical matrices are unchanged.
arg <- grep("^--file=", commandArgs(FALSE), value=TRUE)
script <- sub("^--file=", "", arg)
if (!file.exists(script)) script <- gsub("~\\+~", " ", script)
root <- dirname(normalizePath(script))
meta <- read.csv(file.path(root,"metadata","baseline_full_450.csv"),
                 check.names=FALSE)
stopifnot(nrow(meta)==450L,!anyDuplicated(meta$sample_id))
cols <- c("IBDQ_bowel","IBDQ_systemic","IBDQ_emotional","IBDQ_social",
          "IBDQ_items_valid_n","IBDQ_item_total_recomputed",
          "IBDQ_total_locked_discordant","IBDQ_domain_scoring_status")
stopifnot(all(cols%in%names(meta)))
rownames(meta) <- meta$sample_id
for (name in c("mt_ibd_283.rds","mt_all_449.rds")) {
  path <- file.path(root,"microbiome",name)
  mt <- readRDS(path)
  ids <- rownames(mt$sample_table)
  stopifnot(inherits(mt,"microtable"),all(ids%in%rownames(meta)))
  for (col in cols) mt$sample_table[[col]] <- meta[ids,col]
  stopifnot(identical(rownames(mt$sample_table),colnames(mt$otu_table)))
  saveRDS(mt,path)
}
for (panel in c("full_panel","HMDB_only")) {
  path <- file.path(root,"metabolomics",panel,"input_eligible.rds")
  x <- readRDS(path)
  ids <- rownames(x$meta_enriched)
  stopifnot(length(ids)==299L,all(ids%in%rownames(meta)))
  for (col in cols) x$meta_enriched[[col]] <- meta[ids,col]
  stopifnot(identical(rownames(x$meta_enriched),colnames(x$mat)))
  saveRDS(x,path)
}
cat("IBDQ_RDS_SYNC_PASS: two microeco and two metabolomics objects\n")
