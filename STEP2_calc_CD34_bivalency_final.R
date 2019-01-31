library(peakrefine)
library(BiocFileCache)
library(seqsetvis)
library(GenomicRanges)
library(data.table)
setwd("~/R/waldron/")
source("setup_files.R")

load("CD34_consensus_v1.save")

k4_bams = bams[grepl("CD34", bams) & grepl("K4", bams)]
k4_bws = bws[grepl("CD34", bws) & grepl("K4", bws)]
k27_bws = k27_bws[grepl("CD", k27_bws)]
k27_bams = bams[grepl("CD34", bams) & grepl("K4", bams)]

subsetByOverlaps(k4_consenus, k27_consenus)

bw_dt = ssvFetchBigwig(c(k4_bws, k27_bws), resize(k4_consenus, 1e4, fix = "center"), 
                       win_size = 50, win_method = "summary", return_data.table = TRUE)
cap = 20
bw_dt[, ycap := y]
bw_dt[ycap > cap, ycap := cap]

ssvSignalHeatmap(bw_dt[grepl("K27", sample)], fill_ = "ycap")

qgr = resize(subsetByOverlaps(k4_consenus, k27_consenus), 1e4, fix = "center")
names(qgr) = paste0("region_", seq_along(qgr))
bw_dt = ssvFetchBigwig(c(k4_bws, k27_bws), resize(subsetByOverlaps(k4_consenus, k27_consenus), 1e4, fix = "center"), 
                       win_size = 50, win_method = "summary", return_data.table = TRUE)
cap = 20
bw_dt[, ycap := y]
bw_dt[ycap > cap, ycap := cap]

biv_clust = ssvSignalClustering(bw_dt[grepl("K27", sample)], fill_ = "ycap", max_rows = Inf)
valid = unique(biv_clust[, .(id, cluster_id)])[cluster_id != 4]$id
qgr[valid]
bw_dt = ssvFetchBigwig(c(k4_bws, k27_bws), qgr[valid], 
                       win_size = 50, win_method = "summary", return_data.table = TRUE)
cap = 20
bw_dt[, ycap := y]
bw_dt[ycap > cap, ycap := cap]
biv_clust = ssvSignalClustering(bw_dt, fill_ = "ycap", max_rows = Inf)
ssvSignalHeatmap(biv_clust, fill_ = "ycap")
biv_clust_assign = unique(biv_clust[, .(cluster_id, id)])
valid = biv_clust_assign[!cluster_id %in% 5:6]$id


ssvSignalHeatmap(biv_clust[id %in% valid], fill_ = "ycap")
biv_gr = qgr[valid]
save(biv_gr, file = "CD34_bivalent.save")
