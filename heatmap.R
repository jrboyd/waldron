library(BiocFileCache)
source("setup_files.R")
library(seqsetvis)

k4_peaks_hESC = "/slipstream/galaxy/uploads/working/qc_framework/output_bivalency_redo_patients_H7/H7_H3K4ME3_pooled/H7_H3K4ME3_pooled_peaks.narrowPeak"
k27_peaks_hESC = "/slipstream/galaxy/uploads/working/qc_framework/output_bivalency_redo_patients_H7/H7_H3K27ME3_pooled/H7_H3K27ME3_pooled_peaks.narrowPeak"
k4_bw_hESC = "/slipstream/galaxy/uploads/working/qc_framework/output_bivalency_redo_patients_H7/H7_H3K4ME3_pooled/H7_H3K4ME3_pooled_FE.bw"
k27_bw_hESC = "/slipstream/galaxy/uploads/working/qc_framework/output_bivalency_redo_patients_H7/H7_H3K27ME3_pooled/H7_H3K27ME3_pooled_FE.bw"
names(k4_peaks_hESC) = "H7_H3K4me3"
names(k27_peaks_hESC) = "H7_H3K27me3"
names(k4_bw_hESC) = "H7_H3K4me3"
names(k27_bw_hESC) = "H7_H3K27me3"

k4_peaks = c(k4_peaks_hESC, k4_peaks)


bws = c(k4_bw_hESC, k4_bws, k27_bw_hESC, k27_bws)

k4_gr.raw = easyLoad_narrowPeak((c(k4_peaks_hESC, k27_peaks_hESC)))
k4_gr = lapply(k4_gr.raw, function(x)subset(x, qValue > 10))
k4_olap = ssvOverlapIntervalSets(k4_gr)


rs = apply(mcols(k4_olap), 1, sum)
rs = rep(2, length(k4_olap))
names(rs) = NULL

bfcif = peakrefine::bfcif
bfc = BiocFileCache()

k4_olap = resize(k4_olap, 8000, fix = "center")


ssvFeatureBinaryHeatmap(k4_olap)

rname = digest::digest(list(bws, k4_olap[rs > 1]))
bw_dt = bfcif(bfc, rname, function(){
    ssvFetchBigwig(bws, k4_olap[rs > 1], return_data.table = TRUE)
})

vs = 1000
rname2 = digest::digest(list(bws, k4_olap[rs > 1], vs, "centered"))
bw_dt = bfcif(bfc, rname2, function(){
    centerAtMax(bw_dt, view_size = vs)
})


biv = names(subset(k4_olap, H7_H3K4me3 & H7_H3K27me3))

clust_dt = ssvSignalClustering(bw_dt[id %in% biv], nclust = 8, max_rows = Inf, max_cols = Inf)
clust_dt[, c("cell", "mark") := tstrsplit(sample, "_")]
clust_dt$sample = factor(clust_dt$sample)
clust_dt$sample = factor(clust_dt$sample, levels = levels(clust_dt$sample)[c(13:14, 1:12, 15:28)])
clust_dt$sample = factor(clust_dt$sample, levels = levels(clust_dt$sample)[c(c(1:14)*2, c(1:14)*2-1)])

# tmp = clust_dt[, quantile(y, .99), by = .(cell, mark, id)]
# ggplot(tmp, aes(x = V1)) + facet_grid("cell~mark", scales = "free_x") + geom_histogram()

clust_dt[, ycap := y / quantile(y, .999), by = .(cell, mark)]
clust_dt[ycap > 1, ycap := 1]

id_o = levels(clust_dt$id)
id_k = id_o[seq_along(id_o) %% 8 == 1]

ssvSignalHeatmap(clust_dt[mark == "H3K27me3" & id %in% id_k], fill_ = "ycap", facet_ = "cell")
ggsave("hmap_k27.pdf", width = 12, height = 8)
ssvSignalHeatmap(clust_dt[mark == "H3K4me3" & id %in% id_k], fill_ = "ycap", facet_ = "cell")
ggsave("hmap_k4.pdf", width = 12, height = 8)

ssvSignalHeatmap(clust_dt[mark == "H3K27me3" & id %in% id_k], fill_ = "y", facet_ = "cell")
ggsave("hmap_k27_raw.pdf", width = 12, height = 8)
ssvSignalHeatmap(clust_dt[mark == "H3K4me3" & id %in% id_k], fill_ = "y", facet_ = "cell")
ggsave("hmap_k4_raw.pdf", width = 12, height = 8)

agg_dt = clust_dt[, .(y = mean(y)), by = .(x, sample, cluster_id)]
agg_dt[, c("cell", "mark") := tstrsplit(sample, "_")]
agg_dt[, ycap := y / quantile(y, .999), by = .(cell, mark)]
agg_dt[ycap > 1, ycap := 1]

ggplot(agg_dt, aes(x = x, y = ycap, color = mark)) + geom_path() + facet_grid("cluster_id~cell")
ggsave("hmap_profile.pdf", width = 18, height = 8)


agg_dt[, ycap := y / quantile(y, .999), by = .(mark)]
agg_dt[ycap > 1, ycap := 1]

ggplot(agg_dt, aes(x = x, y = ycap, color = mark)) + geom_path() + facet_grid("cluster_id~cell")
ggsave("hmap_profile_raw.pdf", width = 18, height = 8)
# unique(clust_dt[, .(cluster_id, id)])[, .N, by = .(cluster_id)][order(cluster_id)]

cluster_assignment = unique(clust_dt[, .(id, cluster_id)])
cluster_assignment = cbind(cluster_assignment, as.data.table(k4_olap[cluster_assignment$id]))

ref_gr = rtracklayer::import.gff("~/gencode.v28.annotation.gtf.gz", feature.type = "transcript", format = "gtf")
library(GenomicRanges)
max_dist = 5e3
dists = distanceToNearest(ref_gr, resize(GRanges(cluster_assignment), 2000, fix = "center"))
dists = subset(dists, distance <= max_dist)
anno_dt = cbind(as.data.table(ref_gr[queryHits(dists)])[, .(gene_name, gene_id, transcript_id)], 
                cluster_assignment[subjectHits(dists)], distance = mcols(dists)[[1]])
anno_dt = anno_dt[, .(list(unique(gene_name))), by = .(cluster_id)]
anno_lists = anno_dt$V1
names(anno_lists) = paste0("cluster_", anno_dt$cluster_id)
anno_lists = anno_lists[order(anno_dt$cluster_id)]
names(anno_lists)
lengths(anno_lists)

bg_genes = unique(unlist(anno_lists))
gene_lists = anno_lists

library(GOfuncR)

if(file.exists("go_res.save")){
    load("go_res.save")
}else{
    options(mc.cores = min(length(gene_lists), 8))
    go_res = parallel::mclapply(gene_lists, function(gl){
        message(gl[1])
        input_hyper = rbind(
            data.frame(gene_ids = gl, is_candidate=1),
            data.frame(gene_ids = setdiff(bg_genes, gl), is_candidate=0)
        )
        suppressWarnings(go_enrich(input_hyper, n_randset=1000, silent = TRUE))
    })
    save(go_res, file = "go_res.save")
}



sig_bp_over = lapply(go_res, function(res){
    stats = res[[1]]
    bp_dt = as.data.table(stats)
    bp_dt[ontology == "biological_process", ][FWER_overrep < .05]
})

all_sig_bp_over = unique(unlist(lapply(sig_bp_over, function(x){
    x$node_id
})))

#grab node_id for plotting for each cell
sig_bp_toplot = lapply(go_res, function(res){
    stats = res[[1]]
    bp_dt = as.data.table(stats)
    bp_dt[node_id %in% all_sig_bp_over]
})
sig_bp_toplot = rbindlist(sig_bp_toplot, use.names = TRUE, idcol = "sample")
# no idea why only 3 decimals for FWER
sig_bp_toplot[FWER_overrep == 0, FWER_overrep := .001]
sig_bp_toplot[, lgFWER_overrep := -log10(FWER_overrep)]

clust = ssvSignalClustering(sig_bp_toplot, 
                            nclust = 10,
                            row_ = "node_name", 
                            column_ = "sample", 
                            fill = "lgFWER_overrep")


clust$sample_label = clust$sample

pgo = ggplot(clust) + 
    geom_raster(aes(x = sample_label, y = node_name, fill = lgFWER_overrep)) +
    labs(fill = "-log10 FWER", title = "Biological Process GO", y = "", x= "") +
    scale_x_discrete(position = "top") + 
    scale_y_discrete(position = "right") +
    scale_fill_gradientn(colors = safeBrew(8, "Reds")) +
    guides(fill = guide_colorbar(title.position = "top", label.position = "bottom")) +
    theme(
        axis.text.x.top = element_text(
            angle = 90, 
            # color = col[c(1:2, 2, 3:7)], 
            hjust = 0, vjust = .5), 
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0),
        legend.position = "bottom",
        axis.line = element_blank() 
    ) 
pgo
