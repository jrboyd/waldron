library(peakrefine)
library(BiocFileCache)
library(seqsetvis)
library(GenomicRanges)
library(data.table)
setwd("~/R/waldron/")
source("setup_files.R")

ssvFetchPeak = function(grs, qgr){
    qgr = resize(qgr, 5000, fix = "center")
    gr_dt = lapply(grs, function(x){
        vgr = viewGRangesWinSample_dt(GRanges(coverage(x)), 
                                      qgr = qgr, window_size = 50)
        # mdt = merge(vgr, k4_clust_assignment)
        # ggplot(mdt[id %in% sample(unique(id), 100)], 
        #        aes(x = x, y = id, fill = y == 1)) + 
        #     geom_raster() + 
        #     scale_fill_manual(values = c("TRUE" = "black", "FALSE" = "white"))
        # mdt
        vgr
    })
    rbindlist(gr_dt, use.names = TRUE, idcol = "sample") 
}

k4_bams = bams[grepl("CD34", bams) & grepl("K4", bams)]
k4_bws = bws[grepl("CD34", bws) & grepl("K4", bws)]
k27_bams = bams[grepl("CD34", bams) & grepl("K4", bams)]
k4_peaks = k4_peaks[grepl("CD34", names(k4_peaks)) & grepl("K4", names(k4_peaks))]
k27_peaks = k27_peaks[grepl("CD34", names(k27_peaks)) & grepl("K27", names(k27_peaks))]

k4_peaks = easyLoad_narrowPeak(k4_peaks)
k27_peaks = easyLoad_narrowPeak(k27_peaks)

k27_olaps = ssvOverlapIntervalSets(k27_peaks)

k4_peaks.focused = lapply(k4_peaks, function(x){
    start(x) = start(x) + x$relSummit
    end(x) = start(x)
    resize(x, width = 300, fix = "center")
})
k4_olaps = ssvOverlapIntervalSets(k4_peaks.focused)

ssvFeatureBinaryHeatmap(k4_olaps)
ssvFeatureBinaryHeatmap(k27_olaps)

hist(width(k4_olaps), breaks = 200)
k4_dt = ssvFetchBigwig(k4_bws, resize(k4_olaps, 20000, fix = "center"), 
                       return_data.table = TRUE, win_method = "summary", win_size = 50)
k4_dt.raw = copy(k4_dt)
k4_dt = centerAtMax(k4_dt, view_size = .1)

k4_clust = ssvSignalClustering(k4_dt, max_rows = Inf)
k4_clust_assignment = unique(k4_clust[, .(id, cluster_id)])
clust_memb = unique(k4_clust[, .(id, cluster_id)])[, .(memb_id = list(id)), by = .(cluster_id)][order(cluster_id)]
pcall = t(sapply(clust_memb$memb_id, function(x){
    clust_olap = k4_olaps[names(k4_olaps) %in% x]
    colSums(as.data.frame(mcols(clust_olap))) / length(x)
}))
round(pcall, 2)
cap = 50
k4_clust[, ycap := y]
k4_clust[ycap > cap, ycap := cap]

ssvSignalHeatmap(k4_clust[id %in% sample(unique(id), 300)], fill_ = "ycap") + labs(title = "FE")


valid_id = k4_clust_assignment[cluster_id %in% 3:6]$id
k4_consenus = k4_clust[id %in% valid_id] %>% GRanges %>% reduce

k4_peak_dt = ssvFetchPeak(k4_peaks, k4_olaps)

k4_peak_olap_dt = lapply(k4_peaks, function(x){
    vgr = viewGRangesWinSample_dt(GRanges(coverage(x)), 
                                  qgr = resize(k4_olaps, 5000, fix = "center"), window_size = 50)
    mdt = merge(vgr, k4_clust_assignment)
    # ggplot(mdt[id %in% sample(unique(id), 100)], 
    #        aes(x = x, y = id, fill = y == 1)) + 
    #     geom_raster() + 
    #     scale_fill_manual(values = c("TRUE" = "black", "FALSE" = "white"))
    mdt
})
k4_peak_olap_dt = rbindlist(k4_peak_olap_dt, use.names = TRUE, idcol = "sample")
k4_peak_olap_dt$id = factor(k4_peak_olap_dt$id, levels = levels(k4_clust$id))

ssvSignalHeatmap(k4_peak_olap_dt[id %in% sample(unique(id), 500)]) + labs(title = "peak call")

k4_frac_dt = k4_peak_olap_dt[, .(frac = sum(y) / .N), by = .(sample, cluster_id, x)]
ggplot(k4_frac_dt, aes(x = x, y = frac, color = sample)) + facet_wrap("cluster_id") + geom_path() +
    labs(title = "positional peak call fraction by cluster")

# ggplot(k4_peak_olap_dt[id %in% sample(unique(id), 500)], 
#        aes(x = x, fill = y == 1, y = id)) + 
#     facet_wrap("sample") + 
#     geom_raster() + theme(panel.background = element_blank()) +
#     scale_fill_manual(values = c("TRUE" = "black", "FALSE" = "white"))

k4_avg_dt = k4_clust[, .(y = mean(y)),  by = .(sample, cluster_id, x)]
ggplot(k4_avg_dt, aes(x = x, y = y, color = sample)) + facet_wrap("cluster_id") + geom_path() +
    labs(title = "mean FE by cluster")



valid_k4 = k4_clust[cluster_id != 6]$id %>% unique %>% as.character()

ssvSignalHeatmap(k4_clust[id %in% valid_k4][id %in% sample(unique(id), 100)]) + labs(title = "FE")

reduced_k4 = reduce(resize(k4_olaps[valid_k4], 2000, fix = "center"))
hist(width(reduced_k4))


k27_peaks.focused = lapply(k27_peaks, function(x){
    start(x) = start(x) + x$relSummit
    end(x) = start(x)
    resize(x, width = 600, fix = "center")
})
k27_olaps = ssvOverlapIntervalSets(k27_peaks.focused)
hist(width(k27_olaps), breaks = 200)
k27_bws = k27_bws[grepl("CD", k27_bws)]
# k27_dt = ssvFetchBigwig(k27_bws, resize(k27_olaps, 5000, fix = "center"), 
#                         return_data.table = TRUE, win_size = 50)
k27_dt = ssvFetchBigwig(k27_bws, resize(k27_olaps, 50000, fix = "center"), 
                        return_data.table = TRUE, win_size = 50, 
                        win_method = "summary", 
                        summary_FUN = weighted.mean)


k27_clust = ssvSignalClustering(k27_dt, max_rows = Inf)
k27_clust_assignment = unique(k27_clust[, .(id, cluster_id)])

ssvSignalHeatmap(k27_clust[id %in% sample(unique(id), 500)])

# k4_of_27_dt = ssvFetchBigwig(k4_bws, resize(k27_olaps, 5000, fix = "center"), return_data.table = TRUE)
k4_of_27_dt = ssvFetchBigwig(k4_bws, resize(k27_olaps, 50000, fix = "center"), 
                        return_data.table = TRUE, win_size = 50, 
                        win_method = "summary", 
                        summary_FUN = weighted.mean)

k4_of_27_dt = merge(k4_of_27_dt, k27_clust_assignment)
k4_of_27_dt$id = factor(k4_of_27_dt$id, levels = levels(k27_clust$id))
ssvSignalHeatmap(k4_of_27_dt[id %in% sample(unique(id), 500)])

ssvSignalHeatmap(k4_of_27_dt[cluster_id == 1, .(x, y, id, sample)], nclust = 3)
ssvSignalHeatmap(k4_of_27_dt[cluster_id == 2, .(x, y, id, sample)], nclust = 3)
ssvSignalHeatmap(k4_of_27_dt[cluster_id == 3, .(x, y, id, sample)], nclust = 3)
ssvSignalHeatmap(k4_of_27_dt[cluster_id == 4, .(x, y, id, sample)], nclust = 3)
ssvSignalHeatmap(k4_of_27_dt[cluster_id == 5, .(x, y, id, sample)], nclust = 3)
ssvSignalHeatmap(k4_of_27_dt[cluster_id == 6, .(x, y, id, sample)], nclust = 3)

k27peaks_of_27_dt = ssvFetchPeak(k27_peaks, k27_olaps)
k27peaks_of_27_dt = merge(k27peaks_of_27_dt, k27_clust_assignment)
ssvSignalHeatmap(k27peaks_of_27_dt[id %in% sample(unique(id), 500)])

k4peaks_of_27_dt = ssvFetchPeak(k4_peaks, k27_olaps)
k4peaks_of_27_dt = merge(k4peaks_of_27_dt, k27_clust_assignment)
ssvSignalHeatmap(k4peaks_of_27_dt[id %in% sample(unique(id), 500)])

# ssvSignalHeatmap(k4peaks_of_27_dt[cluster_id == 1, .(x, y, id, sample)], nclust = 3)
# ssvSignalHeatmap(k4peaks_of_27_dt[cluster_id == 2, .(x, y, id, sample)], nclust = 3)
# ssvSignalHeatmap(k4peaks_of_27_dt[cluster_id == 3, .(x, y, id, sample)], nclust = 3)
# ssvSignalHeatmap(k4peaks_of_27_dt[cluster_id == 4, .(x, y, id, sample)], nclust = 3)
# ssvSignalHeatmap(k4peaks_of_27_dt[cluster_id == 5, .(x, y, id, sample)], nclust = 3)
# ssvSignalHeatmap(k4peaks_of_27_dt[cluster_id == 6, .(x, y, id, sample)], nclust = 3)

k4_of_27_frac_dt = k4peaks_of_27_dt[, .(y = sum(y) / .N), by = .(sample, x, cluster_id)]
k27_of_27_frac_dt = k27peaks_of_27_dt[, .(y = sum(y) / .N), by = .(sample, x, cluster_id)]
k_of_27_frac_dt = rbind(k4_of_27_frac_dt, k27_of_27_frac_dt)
k_of_27_frac_dt[, c("cell", "mark") := tstrsplit(sample, "_")]

ggplot(k4_of_27_frac_dt, aes(x = x, y = y, color = sample)) + 
    facet_grid("cluster_id~.") + geom_path() +
    labs(title = "k4 positional peak call fraction by k27 cluster")

ggplot(k_of_27_frac_dt, aes(x = x, y = y, color = mark)) + 
    facet_grid("cluster_id~cell") + geom_path() +
    labs(title = "k positional peak call fraction by k27 cluster")

k27_avg_dt = k27_clust[, .(y = mean(y)),  by = .(sample, cluster_id, x)]
ggplot(k27_avg_dt, aes(x = x, y = y, color = sample)) + facet_wrap("cluster_id") + geom_path() +
    labs(title = "mean FE by cluster")

# k27_peak_olap_dt = lapply(k27_peaks, function(x){
#     vgr = viewGRangesWinSample_dt(GRanges(coverage(x)), 
#                                   qgr = resize(k27_olaps, 5000, fix = "center"), window_size = 50)
#     mdt = merge(vgr, k27_clust_assignment)
#     # ggplot(mdt[id %in% sample(unique(id), 100)], 
#     #        aes(x = x, y = id, fill = y == 1)) + 
#     #     geom_raster() + 
#     #     scale_fill_manual(values = c("TRUE" = "black", "FALSE" = "white"))
#     mdt
# })
# k27_peak_olap_dt = rbindlist(k27_peak_olap_dt, use.names = TRUE, idcol = "sample")
# ssvFetchPeak(k27)
# k27_peak_olap_dt$id = factor(k27_peak_olap_dt$id, levels = levels(k27_clust$id))

ssvSignalHeatmap(k27peaks_of_27_dt[id %in% sample(unique(id), 500)]) + labs(title = "peak call")

k27_frac_dt = k27peaks_of_27_dt[, .(frac = sum(y) / .N), by = .(sample, cluster_id, x)]
ggplot(k27_frac_dt, aes(x = x, y = frac, color = sample)) + facet_wrap("cluster_id") + geom_path() +
    labs(title = "k27 positional peak call fraction by cluster")

# ggplot(k27_peak_olap_dt[id %in% sample(unique(id), 500)], 
#        aes(x = x, fill = y == 1, y = id)) + 
#     facet_wrap("sample") + 
#     geom_raster() + theme(panel.background = element_blank()) +
#     scale_fill_manual(values = c("TRUE" = "black", "FALSE" = "white"))


