library(peakrefine)
library(BiocFileCache)
library(seqsetvis)
library(GenomicRanges)
library(data.table)
setwd("~/R/waldron/")
source("setup_files.R")

k4_peaks = easyLoad_narrowPeak("/slipstream/galaxy/uploads/working/qc_framework/output_bivalency_redo_patients_H7/H7_H3K4ME3_pooled/H7_H3K4ME3_pooled_peaks.narrowPeak")
k27_peaks = easyLoad_narrowPeak("/slipstream/galaxy/uploads/working/qc_framework/output_bivalency_redo_patients_H7/H7_H3K27ME3_pooled/H7_H3K27ME3_pooled_peaks.narrowPeak")

peaks = append(k4_peaks, k27_peaks)
names(peaks) = sub("_pool.+", "", names(peaks))

biv_olaps = ssvOverlapIntervalSets(peaks)


ssvFeatureBinaryHeatmap(biv_olaps)

bws = c("/slipstream/galaxy/uploads/working/qc_framework/output_bivalency_redo_patients_H7/H7_H3K4ME3_pooled/H7_H3K4ME3_pooled_FE.bw",
        "/slipstream/galaxy/uploads/working/qc_framework/output_bivalency_redo_patients_H7/H7_H3K27ME3_pooled/H7_H3K27ME3_pooled_FE.bw")
names(bws) = basename(bws) %>% sub("_pool.+", "", .)

bws = c(bws, k4_bws[1], k27_bws[1])

bams = c("/slipstream/galaxy/uploads/working/qc_framework/output_bivalency_redo_patients_H7/H7_H3K4ME3_pooled/H7_H3K4ME3_pooled.bam",
         "/slipstream/galaxy/uploads/working/qc_framework/output_bivalency_redo_patients_H7/H7_H3K27ME3_pooled/H7_H3K27ME3_pooled.bam")
names(bams) = basename(bams) %>% sub("_pool.+", "", .) %>% sub("ME3", "me3", .)

hist(width(biv_olaps), breaks = 200)
# qgr = resize(subset(biv_olaps, H7_H3K4ME3 & H7_H3K27ME3), 20000, fix = "center")
qgr = resize(subsetByOverlaps(peaks$H7_H3K4ME3, peaks$H7_H3K27ME3), 20000, fix = "center")
biv_dt = ssvFetchBigwig(bws, qgr, 
                        return_data.table = TRUE, win_size = 50, 
                        win_method = "summary", summary_FUN = weighted.mean)
biv_dt.raw = copy(biv_dt)
# biv_dt = biv_dt.raw

cent_dt = centerAtMax(biv_dt[grepl("K4", sample)], by_ = c("id"), view_size = c(-.1, .1), replace_x = FALSE, trim_to_valid = FALSE)
biv_dt = merge(biv_dt, cent_dt[, .(id, x, x_summitPosition)])#, by = c("id", "x"))
biv_dt[, x := x_summitPosition]
biv_dt$x_summitPosition = NULL

xcounts = biv_dt[, .N, x]
xcounts = xcounts[N == max(N)]
biv_dt = biv_dt[x %in% xcounts$x]

# biv_dt.trim = centerAtMax(biv_dt, by_ = c("id"), view_size = c(-.1, .1))
# biv_dt.trim
# biv_dt = biv_dt.trim
biv_dt[, ycap := y]
cap = 50
biv_dt[ycap > cap, ycap := cap]



biv_clust_cd34 = ssvSignalClustering(biv_dt[grepl("CD34", sample)], fill_ = "ycap", max_rows = Inf)
biv_clust_assigned = unique(biv_clust_cd34[, .(id, cluster_id)])

biv_clust = merge(biv_dt, biv_clust_assigned)
biv_clust$id = factor(biv_clust$id, levels = levels(biv_clust_cd34$id))

biv_clust[, ycap := y]
cap = 150
biv_clust[ycap > cap, ycap := cap]

ssvSignalHeatmap(biv_clust[biv_clust$id %in% sample(unique(biv_clust$id), 300) ], fill_ = "ycap")

biv_agg = biv_clust[, .(y = mean(y)) , by = .(sample, cluster_id, x)]
biv_agg[grepl("K27", sample), y := y * 2]
biv_agg.sm = applySpline(biv_agg, n = 5, by_ = c("sample", "cluster_id"))
biv_agg.sm[, c("cell", "mark") := tstrsplit(sample, "_")]
ggplot(biv_agg.sm, aes(x = x, y = y, color = mark)) + facet_grid(cluster_id~cell) + geom_path()

# ref_gr = rtracklayer::import.gff("~/gencode.v28.annotation.gtf.gz", feature.type = "gene", format = "gtf")
ref_gr = rtracklayer::import.gff("~/gencode.v28.annotation.gtf.gz", feature.type = "transcript", format = "gtf")
ssvFetchPeak = function(grs, qgr){
    qgr = resize(qgr, 25000, fix = "center")
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
ref_gr = subset(ref_gr, gene_type == "protein_coding")
ref_hit = ssvFetchPeak(list("genes" = ref_gr, 
                            "tssPlus" = promoters(subset(ref_gr, strand == "+"), 1000, 1000), 
                            "tssMinus" = promoters(subset(ref_gr, strand == "-"), 1000, 1000), 
                            "tss" = promoters(ref_gr, 1000, 1000)), qgr)
ref_hit = merge(ref_hit, biv_clust_assigned)
ref_hit$id = factor(ref_hit$id, levels = levels(biv_clust$id))
# ref_hit[y > 1, y := 1]
ssvSignalHeatmap(ref_hit[id %in% sample(unique(id), 500)])
ref_frac = ref_hit[, .(y = sum(y > 0) / .N), by = .(x, cluster_id, sample)]
ggplot(ref_frac, aes(x = x, y = y)) + facet_grid("cluster_id~sample") + geom_path()


        