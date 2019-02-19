library(BiocFileCache)
library(seqsetvis)
library(GenomicRanges)
setwd("~/R/waldron/")
source("setup_files.R")

peak_files = k4_peaks[grepl("^CD", names(k4_peaks))]
bam_files = k4_bams[grepl("^CD", names(k4_bams))]

bam27_files = k27_bams[grepl("^CD", names(k4_bams))]

peak_grs = easyLoad_narrowPeak(peak_files)
peak_grs = lapply(peak_grs, function(x){
    start(x) = start(x) + x$relSummit
    end(x) = start(x)
    resize(x, 300, fix = "center")
})
olaps = ssvOverlapIntervalSets(peak_grs)
ssvFeatureBars(olaps)
ssvFeatureEuler(olaps)
ssvFeatureBinaryHeatmap(olaps, raster_approximation = TRUE)

options(mc.cores = 8)

q_olaps = resize(sample(olaps, 2000), 2500, fix = "center")

bam_dt = ssvFetchBam(bam_files, q_olaps, target_strand = "+",
                     return_data.table = TRUE, fragLens = 50, max_dupes = 1)
bam_dt$y_minus = ssvFetchBam(bam_files, q_olaps, target_strand = "-",
                             return_data.table = TRUE, fragLens = 50, max_dupes = 1)$y
bam_dt[, diff := y - y_minus]

cap = 10
bam_dt[, diff_cap := diff]
bam_dt[diff_cap > cap, diff_cap := cap]
bam_dt[diff_cap < -cap, diff_cap := -cap]

clust_dt = ssvSignalClustering(bam_dt, max_rows = 2000, fill_ = "diff_cap")

ssvSignalHeatmap(clust_dt, fill_ = "diff_cap")

agg_dt = clust_dt[, .(agg = mean(diff)) , by = .(cluster_id, x)][order(cluster_id)]
ggplot(agg_dt, aes(x = x, y = agg)) + 
    geom_path() + facet_grid("cluster_id~.") + 
    theme(panel.background = element_blank()) +
    theme(panel.grid.major = element_line(color = "gray"), panel.grid.minor = element_blank())

###
bam27_dt = ssvFetchBam(bam27_files, q_olaps, target_strand = "+",
                     return_data.table = TRUE, fragLens = 50, max_dupes = 1)
bam27_dt$y_minus = ssvFetchBam(bam27_files, q_olaps, target_strand = "-",
                             return_data.table = TRUE, fragLens = 50, max_dupes = 1)$y
bam27_dt[, diff := y - y_minus]

cap = 10
bam27_dt[, diff_cap := diff]
bam27_dt[diff_cap > cap, diff_cap := cap]
bam27_dt[diff_cap < -cap, diff_cap := -cap]

cluster_assignments = unique(clust_dt[, .(id, cluster_id)])
# clust27_dt = ssvSignalClustering(bam27_dt, max_rows = 2000, fill_ = "diff_cap")

clust27_dt = merge(bam27_dt, cluster_assignments)

ssvSignalHeatmap(clust27_dt, fill_ = "diff_cap")

agg27_dt = clust27_dt[, .(agg = mean(diff)) , by = .(cluster_id, x, sample)][order(cluster_id)]
unique(clust27_dt[, .(id, cluster_id)])[, .(.N), cluster_id]
ggplot(agg27_dt, aes(x = x, y = agg)) + 
    geom_path() + facet_grid("cluster_id~sample") + 
    theme(panel.background = element_blank()) +
    theme(panel.grid.major = element_line(color = "gray"), panel.grid.minor = element_blank()) +
    scale_y_continuous(limits = c(-6,6))
