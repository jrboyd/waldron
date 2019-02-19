library(magrittr)
library(GenomicRanges)
library(data.table)
library(seqsetvis)
library(BiocFileCache)

bfc = BiocFileCache()
bfcif = peakrefine::bfcif

wd = "/slipstream/galaxy/uploads/working/qc_framework/output_waldron_bivalency/"
bams = dir(wd, pattern = "ed$", full.names = TRUE) %>% dir(., pattern = ".bam$", full.names = TRUE)
names(bams) = basename(bams) %>% sub("_pool.+", "", .) %>% sub("ME3", "me3", .)
names(bams)

bws = dir(wd, pattern = "ed$", full.names = TRUE) %>% dir(., pattern = "FE.bw$", full.names = TRUE)
names(bws) = basename(bws) %>% sub("_pool.+", "", .) %>% sub("ME3", "me3", .)
names(bws)

k4_peaks = dir(wd, pattern = "ed$", full.names = TRUE) %>% dir(., pattern = "H3K4.+peaks.narrowPeak$", full.names = TRUE)
names(k4_peaks) = basename(k4_peaks) %>% sub("_pooled.+", "", .) %>% sub("ME3", "me3", .)
k27_peaks = dir(wd, pattern = "ed$", full.names = TRUE) %>% dir(., pattern = "H3K27.+peaks.narrowPeak$", full.names = TRUE)
names(k27_peaks) = basename(k27_peaks) %>% sub("_pooled.+", "", .) %>% sub("ME3", "me3", .)


k4_peaks = easyLoad_narrowPeak(k4_peaks)
# k27_peaks = easyLoad_broadPeak(k27_peaks)
k27_peaks = easyLoad_narrowPeak(k27_peaks)
k27_peaks = lapply(k27_peaks, function(x)subset(x, pValue > 10))

biv_peaks =lapply(seq_along(k4_peaks), function(i){
    olaps = ssvOverlapIntervalSets(list("K4" = k4_peaks[[i]], "K27" = k27_peaks[[i]]))
    olaps = subset(olaps, K4 & K27)
    mcols(olaps) = NULL
    olaps
})
names(biv_peaks) = sub("_.+", "", names(k4_peaks))

k4_olaps = ssvOverlapIntervalSets(k4_peaks)
k27_olaps = ssvOverlapIntervalSets(k27_peaks)
biv_olaps = ssvOverlapIntervalSets(biv_peaks)

ssvFeatureBars(k4_olaps)
ssvFeatureBars(k27_olaps)
ssvFeatureBars(biv_olaps)

### check k4
hist(width(k4_olaps))
k4_olaps = resize(k4_olaps, 1500, fix = "center")

k4_bw = bws[grepl("K4", names(bws))]
stopifnot(all(names(k4_peaks) == names(k4_bw)))

k4_dt = lapply(names(k4_peaks), function(nam){
    qgr = k4_peaks[[nam]]
    qgr = resize(qgr, 1500, fix = "center")
    FUN = function()
        ssvFetchBigwig(k4_bw[nam], qgr, return_data.table = TRUE, summary_FUN = function(x, w)max(x), win_method = "summary")
    bfcif(bfc, digest::digest(list(qgr, k4_bw[nam])), FUN = FUN)
})
k4_dt = rbindlist(k4_dt)
# k4_dt = ssvFetchBigwig(bws[grepl("K4", names(bws))], k4_olaps, return_data.table = TRUE)
k4_dt[, c("cell", "mark") := tstrsplit(sample, "_")]
k4_max_dt = k4_dt[, .(maxFE = max(y)), by = .(cell, mark, id)]
ggplot(k4_max_dt, aes(x = maxFE, group = cell)) + 
    geom_histogram(bins = 200) + 
    labs(title = "H3K4me3")  + 
    facet_wrap("cell", scales = "free_y") + 
    coord_cartesian(xlim = c(0,150))

### check k27
k27_bw = bws[grepl("K27", names(bws))]
stopifnot(all(names(k27_peaks) == names(k27_bw)))

k27_dt = lapply(names(k27_peaks), function(nam){
    qgr = k27_peaks[[nam]]
    qgr = resize(qgr, 1500, fix = "center")
    
    FUN = function()
        ssvFetchBigwig(k27_bw[nam], qgr, return_data.table = TRUE, summary_FUN = function(x, w)max(x), win_method = "summary")
    bfcif(bfc, digest::digest(list(qgr, k4_bw[nam])), FUN = FUN)
})
k27_dt = rbindlist(k27_dt)
# k27_dt = ssvFetchBigwig(bws[grepl("k27", names(bws))], k27_olaps, return_data.table = TRUE)
k27_dt[, c("cell", "mark") := tstrsplit(sample, "_")]
k27_max_dt = k27_dt[, .(maxFE = max(y)), by = .(cell, mark, id)]
ggplot(k27_max_dt, aes(x = maxFE, group = cell)) + 
    geom_histogram(bins = 200) + 
    labs(title = "H3K27me3") +
    facet_wrap("cell", scales = "free_y") + 
    coord_cartesian(xlim = c(0,150))

###
# 
# bw_dt = ssvFetchBigwig(bws, subset(biv_olaps, `CD34-01562` == TRUE), return_data.table = TRUE)
# bw_dt[, c("cell", "mark") := tstrsplit(sample, "_")]
# 
# bw_dt[, max(y), by = .(cell, mark, id)]
# 
# theme_set(theme_classic())
# 
# key = 1:8 + 0
# biv_olaps[key]
# 
# ggplot(bw_dt[id %in% unique(bw_dt$id)[key]], aes(x = x, y = y, color = mark)) + 
#     geom_path(size= 1.2) + 
#     facet_grid("id~cell") + 
#     scale_color_manual(values = c("H3K27me3" = "darkred", "H3K4me3" = "darkgreen"))


