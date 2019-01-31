library(peakrefine)
library(BiocFileCache)
library(seqsetvis)
library(GenomicRanges)
library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)
library(cowplot)
setwd("~/R/waldron/")
source("setup_files.R")

load("CD34_bivalent.save")
k4_bw_hESC = c("H7_H3K4me3" = "/slipstream/galaxy/uploads/working/qc_framework/output_bivalency_redo_patients_H7/H7_H3K4ME3_pooled/H7_H3K4ME3_pooled_FE.bw")
k27_bw_hESC = c("H7_H3K27me3" = "/slipstream/galaxy/uploads/working/qc_framework/output_bivalency_redo_patients_H7/H7_H3K27ME3_pooled/H7_H3K27ME3_pooled_FE.bw")
k4_bws  = c(k4_bws, k4_bw_hESC)
k27_bws = c(k27_bws, k27_bw_hESC)

qgr = resize(biv_gr, 50000, fix = "center")
bw_dt = ssvFetchBigwig(c(k4_bws, k27_bws), qgr, 
                       return_data.table = TRUE, win_method = "summary", win_size = 100)

cap = 10
bw_dt[, ycap := y]
bw_dt[ycap > cap, ycap := cap]
bw_dt[, sample := sub("_", "\n", sample)]
bw_clust = ssvSignalClustering(bw_dt[grepl("K27", sample)])
# ssvSignalHeatmap(bw_clust, fill_ = "ycap")
# ggsave("plot4.pdf", width = 15, height = 6)

bw_clust = ssvSignalClustering(bw_dt)
# ssvSignalHeatmap(bw_clust, fill_ = "ycap")
# ggsave("plot5.pdf", width = 30, height = 6)

ref_gr = rtracklayer::import.gff("~/gencode.v28.annotation.gtf.gz", feature.type = "transcript", format = "gtf")
ref_gr = subset(ref_gr, gene_type == "protein_coding")

all_gr = append(easyLoad_narrowPeak(k27_peaks), easyLoad_narrowPeak(k4_peaks)) %>% ssvOverlapIntervalSets %>% reduce
all_genes = my_annotate(all_gr, ref_gr)$gene_name %>% unique

qhits = list("genes" = ref_gr, 
             "tssPlus" = promoters(subset(ref_gr, strand == "+"), 1000, 1000), 
             "tssMinus" = promoters(subset(ref_gr, strand == "-"), 1000, 1000), 
             "tss" = promoters(ref_gr, 1000, 1000))
ref_hit = ssvFetchPeak(qhits, qgr)


bw_dt[, c("cell", "mark") := tstrsplit(sample, "\n")]
td = bw_dt$cell %>% unique

library(BiocFileCache)
bfcif = peakrefine::bfcif
bfc = BiocFileCache()

nclusts = rep(4, length(td))
names(nclusts) = td
nclusts[c("Loucy", "DND-41", "DOHH2", "Kasumi1", "Loucy", "Nalm6", "OCI-LY1", "PDX2", "SU-DHL-6")] = 3
nclusts[c("U937")] = 2

odir = "biv_varClust_wKegg"
dir.create(odir)

bw_dt[, x := x *50000]

pdf(file.path(odir, "All_combined_plots.pdf"), width = 14, height = 14)
for(cl in td){
    message(cl, " ", which(cl == td), "/", length(td))
    bw_clust = ssvSignalClustering(bw_dt[cell == cl], fill_ = "ycap", nclust = nclusts[cl])
    clust_assign = unique(bw_clust[, .(id, cluster_id)])
    
    p_heat = ssvSignalHeatmap(bw_clust, fill_ = "ycap", facet_ = "mark") + 
        labs(x = "kbp", y = "", fill = "FE") + 
        scale_x_continuous(labels = function(x)x/1000)
    # plot(p_heat + labs(title = cl))
    
    bw_agg = bw_clust[, .(y = mean(ycap)), by = .(cell, mark, cluster_id, x)]
    bw_agg.sm = applySpline(bw_agg, n = 10, by_ = c("cell", "mark", "cluster_id"))
    # bw_agg.sm[mark == "H3K27me3", y := y * 3]
    
    p_agg = ggplot(bw_agg.sm, aes(x = x, y = y, color = mark)) + 
        geom_path(size = 2.5) + facet_grid("cluster_id~.") + 
        scale_color_manual(values = c("H3K4me3" = "forestgreen", "H3K27me3" = "firebrick1")) +
        labs(x = "kbp", y = "mean FE") +
        scale_x_continuous(labels = function(x)x/1000) +
        scale_y_continuous(breaks = c(0,5,10), labels = round)
    # plot(p_agg + labs(title = cl))
    
    mdt = merge(ref_hit, clust_assign)
    mdt_frac = mdt[, .(y = sum(y > 0) /  .N), by = .(cluster_id, x, sample)]
    # ssvSignalHeatmap(mdt_frac)
    p_ref = ggplot(mdt_frac, aes(x = x, y = y)) + geom_path() + facet_grid("cluster_id~sample")
    # plot(p_ref + labs(title = cl))
    
    # go_res = my_clusterProfiler(qgr, clust_assign)
    # p_go = go_res[[2]]
    # plot(p_go + labs(title = cl, subtitle = "vs CD34 bivalent genes"))
    
    go_res2 = my_clusterProfiler(qgr, clust_assign, bg_genes = all_genes)
    p_go2 = go_res2[[2]]
    # plot(p_go2 + labs(title = cl, subtitle = "vs all genes"))
    
    kegg_res = my_clusterProfiler_kegg(qgr, clust_assign, bg_genes = all_genes)
    p_kegg = kegg_res[[2]]
    
    pg = plot_grid(p_heat, p_agg)

    pg = ggdraw() + 
        draw_plot(p_heat + labs(title = cl), x = 0, y = .6, width = .5, height = .4) +
        draw_plot(p_agg, x = .5, y = .6, width = .5, height = .4) +
        draw_plot(p_kegg + labs(title = "KEGG enrichment"), x = 0, y = .3, width = 1, height = .3) +
        draw_plot(p_go2 + labs(title = "GO BP enrichment"), x = 0, y = 0, width = 1, height = .3)
    ggsave(filename = file.path(odir, paste0("combined_plots_", cl, ".pdf")), plot = pg, width = 14, height = 14)
    # pg = cowplot::plot_grid(p_go, p_go2, ncol = 1)
    plot(pg)
    
    # cdat = kegg_res[[1]]@compareClusterResult
    # cdat = as.data.table(cdat)
    # 
    # gs_key = unique(cdat[, .(ID, Description)])
    # setkey(gs_key, ID)
    # 
    # cmat = cdat[, .(uniprot = tstrsplit(geneID, "/")), by = .(ID,Cluster)]
    # cmat$uniprot = unlist(cmat$uniprot)
    # ctab =  dcast(cmat[Cluster == "c3"], "ID~uniprot")#, value.var = length)
    # dmat = as.matrix(ctab[,-1])
    # dmat = ifelse(is.na(dmat), 0, 1)
    # rownames(dmat) = gs_key[.(ctab$ID)]$Description
    # 
    # colnames(dmat) = bitr(colnames(dmat), fromType = "UNIPROT", toType = "SYMBOL", OrgDb = org.Hs.eg.db)$SYMBOL
    # 
    # library(DT)
    # datatable(t(dmat)) %>% formatStyle(
    #     T,
    #     backgroundColor = styleEqual(c(0, 1), c('white', 'blue')), 
    #     color = "00000000"
    #     
    # )
    # cdat[2, list(tstrsplit(geneID, "/"))]
}
dev.off()
