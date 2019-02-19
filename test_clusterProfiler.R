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

bams = c("/slipstream/galaxy/uploads/working/qc_framework/output_bivalency_redo_patients_H7/H7_H3K4ME3_pooled/H7_H3K4ME3_pooled.bam",
         "/slipstream/galaxy/uploads/working/qc_framework/output_bivalency_redo_patients_H7/H7_H3K27ME3_pooled/H7_H3K27ME3_pooled.bam")
names(bams) = basename(bams) %>% sub("_pool.+", "", .)

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
biv_dt[ycap > 30, ycap := 30]



biv_clust = ssvSignalClustering(biv_dt, fill_ = "ycap", max_rows = Inf)
biv_clust_assigned = unique(biv_clust[, .(id, cluster_id)])

ssvSignalHeatmap(biv_clust[biv_clust$id %in% sample(unique(biv_clust$id), 300)], fill_ = "ycap")

biv_agg = biv_clust[, .(y = mean(y)) , by = .(sample, cluster_id, x)]
biv_agg[sample == "H7_H3K27ME3", y := y * 2]
biv_agg.sm = applySpline(biv_agg, n = 5, by_ = c("sample", "cluster_id"))
ggplot(biv_agg.sm, aes(x = x, y = y, color = sample)) + facet_grid(cluster_id~.) + geom_path()

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

## GO heatmap

my_annotate = function(gr, 
                       ref_gr = rtracklayer::import.gff("~/gencode.v28.annotation.gtf.gz", 
                                                        feature.type = "transcript", format = "gtf"), 
                       gr_size = 2000){
    
    
    max_dist = 5e3
    dists = distanceToNearest(ref_gr, resize(gr, gr_size, fix = "center"))
    dists = subset(dists, distance <= max_dist)
    anno_dt = cbind(as.data.table(ref_gr[queryHits(dists)])[, .(gene_name, gene_id, transcript_id)], 
                    as.data.table(gr[subjectHits(dists)]), distance = mcols(dists)[[1]])
    anno_dt
}

peak_dt = as.data.table(peaks$H7_H3K4ME3)
mdt = merge(peak_dt, biv_clust_assigned, by.x = "name", by.y = "id")
anno_dt = my_annotate(GRanges(mdt))

anno_dt = anno_dt[, .(list(unique(gene_name))), by = .(cluster_id)]
anno_lists = anno_dt$V1
names(anno_lists) = paste0("cluster_", anno_dt$cluster_id)
anno_lists = anno_lists[order(anno_dt$cluster_id)]
names(anno_lists)
lengths(anno_lists)

# each cluster vs total genes in heatmap
bg_genes = unique(unlist(anno_lists))

# each cluster vs all marked by k4me3
bg_dt = my_annotate(peaks$H7_H3K4ME3)
# bg_dt = my_annotate(k4_gr$`CD34-01517_H3K4me3`)
bg_genes = unique(bg_dt$gene_name)

gene_lists = anno_lists

BiocManager::install("clusterProfiler", version = "3.8")
library(clusterProfiler)
data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]
gene.df <- bitr(gene, fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Hs.eg.db)

gene.univ <- names(geneList)
gene.df.univ <- bitr(gene.univ, fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Hs.eg.db)

head(gene.df)
ggo <- groupGO(gene     = gene,
               OrgDb    = org.Hs.eg.db,
               ont      = "CC",
               level    = 3,
               readable = TRUE)

head(ggo)

ego3 <- enrichGO(gene         = gene.df$SYMBOL,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)
head(summary(ego3))
dotplot(ego3)


ego <- enrichGO(gene          = gene,
                universe      = names(geneList),
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)

ego4 <- enrichGO(gene          = gene.df$SYMBOL,
                keyType       = "SYMBOL",
                universe      = gene.df.univ$SYMBOL,
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)
head(ego4)

sego4 = simplify(ego4)

barplot(sego4, drop=TRUE, showCategory=12)
dotplot(sego4)
emapplot(ego4)

kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               nPerm        = 1000,
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
head(kk2)

gseaplot(kk2, geneSetID = "hsa04145")

data(gcSample)
lapply(gcSample, head)
gcSample.sym = lapply(gcSample, function(x){
    bitr(x, fromType = "ENTREZID",
                         toType = c("ENSEMBL", "SYMBOL"),
                         OrgDb = org.Hs.eg.db)$SYMBOL
})
ck <- compareCluster(geneCluster = gcSample.sym, 
                     universe      = gene.df.univ$SYMBOL,
                     fun = "enrichGO", 
                     OrgDb = org.Hs.eg.db,                  
                     keyType       = 'SYMBOL',
                     ont           = "CC",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05)
head(as.data.frame(ck))

dotplot(ck)
sck = simplify(ck)
dotplot(sck)

library(GOfuncR)
library(BiocFileCache)
bfc = BiocFileCache()
bfcif = peakrefine::bfcif
# if(file.exists("go_res.save")){
#     load("go_res.save")
# }else{
# options(mc.cores = min(length(gene_lists), 8))
go_res = lapply(gene_lists, function(gl){
    # go_res = parallel::mclapply(gene_lists, function(gl){
    bfcif(bfc, digest::digest(list(gl, bg_genes, "v1")), function(){
        message(length(gl), " / ", length(bg_genes) )
        input_hyper = rbind(
            data.frame(gene_ids = gl, is_candidate=1),
            data.frame(gene_ids = setdiff(bg_genes, gl), is_candidate=0)
        )
        suppressWarnings(go_enrich(input_hyper, n_randset=1000, silent = TRUE))
    })
})
# save(go_res, file = "go_res.save")
# }

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
# ssvSignalHeatmap(ref_hit[cluster_id == 1, .(x, y = ifelse(y > 0, 1, 0), sample, id)], max_cols = Inf)

        