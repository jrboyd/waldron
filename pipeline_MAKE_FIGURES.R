#MAKE_FIGURES
# Figure1: venn of A) H7 H3K4me3/27me3 B) CD34+ H3K4me3/K27me3 C) H7 bv/CD34+ bv
# Figure2: Current figure of CD34+ bv-marked genes; snapshot of a region
# Figure3: comparison to several (not necessary all?) cell lines via clustering signal.
# Figure4: Kasumi and another leukemia lineage (PDX2?) maintain most but not all bv domains. BV domains lost + retained pathways.
# Figure5: prognostication of genes
source('setup_files.R')
source("functions.R")
tx_gr =rtracklayer::import.gff("~/gencode.v28.annotation.gtf.gz", format = "gtf", feature.type  = "transcript")
tx_gr = subset(tx_gr, transcript_support_level %in% 1:2)
pr_gr = promoters(tx_gr, 1, 1)
pr_gr = split(pr_gr, pr_gr$gene_name)
pr_gr = resize(unlist(reduce(resize(pr_gr, 200, fix = "center"))), 1, fix = "center")

# Figure1: venn of A) H7 H3K4me3/27me3 B) CD34+ H3K4me3/K27me3 C) H7 bv/CD34+ bv
np_grs = easyLoad_narrowPeak(c(k4_peaks["H7_H3K4me3"], k4_peaks[grepl("CD34", k4_peaks)],
                               k27_peaks["H7_H3K27me3"], k27_peaks[grepl("CD34", k4_peaks)]))
names(np_grs)
olaps_h7 = ssvOverlapIntervalSets(np_grs[grepl("H7", names(np_grs))])                    
p_vennH7 = ssvFeatureVenn(olaps_h7, circle_colors = c("red1", "red4"), fill_alpha = .5)
biv_h7 = subset(olaps_h7, H7_H3K4me3 & H7_H3K27me3)

#CD34 k4
olaps_cd34_k4 = ssvOverlapIntervalSets(np_grs[grepl("CD34", names(np_grs)) & grepl("K4", names(np_grs))])                    
ssvFeatureEuler(olaps_cd34_k4)
ssvFeatureBinaryHeatmap(olaps_cd34_k4)
is_k4_consensus = rowSums(as.data.frame(mcols(olaps_cd34_k4))) > 2
sum(is_k4_consensus)
olaps_cd34_k4 = olaps_cd34_k4[is_k4_consensus]

#CD34 k27
olaps_cd34_k27 = ssvOverlapIntervalSets(np_grs[grepl("CD34", names(np_grs)) & grepl("K27", names(np_grs))])                    
ssvFeatureEuler(olaps_cd34_k27)
ssvFeatureBinaryHeatmap(olaps_cd34_k27)
is_k4_consensus = rowSums(as.data.frame(mcols(olaps_cd34_k27))) > 2
sum(is_k4_consensus)
olaps_cd34_k27 = olaps_cd34_k27[is_k4_consensus]
ssvFeatureEuler(olaps_cd34_k27)
ssvFeatureBinaryHeatmap(olaps_cd34_k27)

#CD34 bivalent
olaps_cd34 = ssvOverlapIntervalSets(list("CD34_H3K4me3" = olaps_cd34_k4, "CD34_H3K27me3" = olaps_cd34_k27))
p_vennCD34 = ssvFeatureVenn(olaps_cd34,  circle_colors = c("blue1", "blue4"), fill_alpha = .5)
biv_cd34 = subset(olaps_cd34, CD34_H3K4me3 & CD34_H3K27me3)

#H7 and CD34 bivalent comparison
olaps_biv = ssvOverlapIntervalSets(list("H7_bivalent" = biv_h7, "CD34_bivalent" = biv_cd34))
p_vennBiv = ssvFeatureVenn(olaps_biv, circle_colors = c("red2", "blue2"), fill_alpha = .5)
grps_biv = ssvFactorizeMembTable(olaps_biv)

cowplot::plot_grid(p_vennH7, p_vennCD34, p_vennBiv)
ggsave("fig1_venns.pdf", width = 8, height = 4)

fig1_bws = c(k4_bws[grepl("H7", k4_bws) | grepl("CD34", k4_bws)],
             k27_bws[grepl("H7", k27_bws) | grepl("CD34", k27_bws)])
fig1_bws_dt = data.table(files = fig1_bws, sample = names(fig1_bws))
fig1_bws_dt[, c("cell", "mark") := tstrsplit(sample, "_")]
fig1_bws_dt[, c("cell") := tstrsplit(cell, "-", keep = 1)]
hist(width(olaps_biv))
olaps_biv = resize(olaps_biv, 15000, fix = "center")
fig1_dt = bfcif(bfc, "fig1_dt_v1", function(){
    ssvFetchBigwig(fig1_bws_dt, olaps_biv, 
                   win_size = 50, win_method = "summary", 
                   return_data.table = TRUE, n_cores = 16)
})

# fig1_dt[, cell := tstrsplit(cell, "-", keep = 1)]
fig1_dt = fig1_dt[, .(y = mean(y)), .(cell, mark, x, id)]

cap = 30
fig1_dt[, ycap := y]
fig1_dt[ycap > cap, ycap := cap]
fig1_dt = merge(fig1_dt, grps_biv)
fig1_dt[, sample := paste(cell, mark)]

#Heatmaps of CD34 and H7 bivalency
pdf("fig1_venn_clustering.pdf")
theme_set(theme_classic())
clust_bivH7 = ssvSignalClustering(
    fig1_dt[group == "H7_bivalent"], 
    fill_ = "ycap", nclust = 2, max_rows = Inf) 
ssvSignalHeatmap(clust_bivH7[id %in% sampleCap(id)], 
                 fill_ = "ycap") + 
    labs(title = "bivalent in H7, not CD34")
ggplot(clust_bivH7[, .(y = mean(ycap)), .(cell, mark, x, cluster_id)], 
       aes(x = x, y = y, color = mark)) + 
    geom_path(size = 2) + 
    scale_color_manual(values = 
                           c("H3K27me3" = "red", "H3K4me3" = "forestgreen")) +
    facet_grid(cluster_id~cell) + 
    labs(title = "bivalent in H7, not CD34")

clust_bivCD34 = ssvSignalClustering(
    fig1_dt[group == "CD34_bivalent"], 
    fill_ = "ycap", nclust = 2, max_rows = Inf) 
ssvSignalHeatmap(clust_bivCD34[id %in% sampleCap(id)], 
                 fill_ = "ycap") + 
    labs(title = "bivalent in CD34, not H7")
ggplot(clust_bivCD34[, .(y = mean(ycap)), .(cell, mark, x, cluster_id)], 
       aes(x = x, y = y, color = mark)) + 
    geom_path(size = 2) + 
    scale_color_manual(values = 
                           c("H3K27me3" = "red", "H3K4me3" = "forestgreen")) +
    facet_grid(cluster_id~cell)+ 
    labs(title = "bivalent in CD34, not H7")

clust_bivBoth = ssvSignalClustering(
    fig1_dt[group == "H7_bivalent & CD34_bivalent"], 
    fill_ = "ycap", nclust = 2, max_rows = Inf) 
ssvSignalHeatmap(clust_bivBoth[id %in% sampleCap(id)], 
                 fill_ = "ycap") + 
    labs(title = "bivalent in both H7 and CD34")
ggplot(clust_bivBoth[, .(y = mean(ycap)), .(cell, mark, x, cluster_id)], 
       aes(x = x, y = y, color = mark)) + 
    geom_path(size = 2) + 
    scale_color_manual(values = 
                           c("H3K27me3" = "red", "H3K4me3" = "forestgreen")) +
    facet_grid(cluster_id~cell) + 
    labs(title = "bivalent in both H7 and CD34")
dev.off()

# Figure2: Current figure of CD34+ bv-marked genes; snapshot of a region
#view centered around any tss marked by k4 or k27
pr_gr_cd34 = subsetByOverlaps(pr_gr, olaps_cd34)
pr_gr_cd34$id = paste0("tss_", seq_along(pr_gr_cd34))
pr_gr_cd34$gene_name = names(pr_gr_cd34)
names(pr_gr_cd34) = NULL
fig2_dt = bfcif(bfc, "fig2_dt_v1", function(){
    ssvFetchBigwig(fig1_bws_dt[cell == "CD34"], resize(pr_gr_cd34, 5000, fix = "center"), 
                   win_size = 50, win_method = "summary", 
                   return_data.table = TRUE, n_cores = 16)
})

fig2_dt = fig2_dt[, .(y = mean(y)), .(id, x, cell, mark)]
fig2_dt[, sample := paste(cell, mark)]
fig2_dt[, ycap := y]
fig2_dt[ycap > cap, ycap := cap]
clust_cd34 = ssvSignalClustering(
    fig2_dt, 
    fill_ = "ycap", 
    nclust = 14, facet_ = "mark",
    max_rows = Inf) 
p_fig2 = ssvSignalHeatmap(clust_cd34[id %in% sampleCap(id)], facet_ = "mark",
                          fill_ = "ycap") + 
    labs(title = "CD34 k4 and k27 at TSSes")
ggsave("fig2_CD34_heatmap.pdf", p_fig2, width = 6, height = 4)

unique(clust_cd34[, .(id, cluster_id)])[, .N, .(cluster_id)][order(cluster_id)]

mdt = merge(unique(clust_cd34[, .(id, cluster_id)]), as.data.table(pr_gr_cd34)[, .(id, gene_name)])
cluster_genes = split(mdt$gene_name, mdt$cluster_id)

fig2_clust_res = my_clusterProfiler_fromGenes(cluster_genes)#, force_overwrite = TRUE)

sego = simplify(fig2_clust_res[[1]])
p = dotplot(sego, showCategory = Inf)
ggsave("fig2_cluster_go.pdf", p, width = 20, height = 60, limitsize = FALSE)

fig2_kegg_res = my_clusterProfiler_kegg_fromGenes(cluster_genes)#$, force_overwrite = TRUE)
p = dotplot(fig2_kegg_res[[1]], showCategory = 150)
ggsave("fig2_cluster_kegg.pdf", p, width = 20, height = 25, limitsize = FALSE)


view_size = 15e3
n_points = 16
theme_set(theme_classic())
set.seed(0)
qgr = pr_gr_cd34
qgr = resize(qgr, view_size, fix = "center")
qbw = c(k4_bws, k27_bws)
qdt = data.table(qbw = qbw)
qdt[, c("cell", "mark") := tstrsplit(basename(qbw), "_", keep = 1:2)]
qdt[, mark := sub("ME3", "me3", mark)]
qdt[grepl("CD34", cell), cell := "CD34"]
# stopifnot(length(unique(table(qdt$cell))) == 1)
stopifnot(length(unique(table(qdt$mark))) == 1)

message("fetch tsne input")
source("functions_tsne.R")
options(mc.cores = 32)
tsne_input = fetch_tsne_mat(qdt, qgr, 
                            qwin = 50, 
                            qmet = "summary", 
                            cap_value = 30, 
                            high_on_right = FALSE)

message("run tsne")
tsne_res = run_tsne(tsne_input$tsne_mat, perplexity = 100)

# message("run tsne")
# tsne_res = run_tsne(tsne_input$tsne_mat, perplexity = 2000)
tp = sample(unique(tsne_res$id), min(500, length(unique(tsne_res$id))))
tsne_res.tp = tsne_res[id %in% tp]



p_basic = ggplot() + 
    annotate("point", x = tsne_res.tp$tx, y = tsne_res.tp$ty, color = "lightgray")
p_basic
ggplot() + 
    annotate("point", x = tsne_res.tp$tx, y = tsne_res.tp$ty, color = "lightgray") +
    geom_point(data = tsne_res.tp[grepl("CD34", cell)], aes(x = tx, y = ty, color = cell)) 

message("make images")
img_res = make_tsne_img(
    bw_dt = tsne_input$bw_dt,
    tdt = tsne_res, #force_rewrite = TRUE, 
    n_points = n_points
)

ggplot(tsne_res, aes(x = tx, y = ty, color = cell)) + geom_density2d() + facet_wrap("cell")
# p_profiles = ggplot(img_res$images_dt, aes(x = tx, y = ty, image = png_file)) + geom_image()

p_density = plot_tsne_img(img_res$images_dt, n_points = n_points, 
                          N_ceiling = NULL, N_floor = -50, min_size = 0,
                          show_plot = FALSE)$plot
spc = 1/n_points/2*.9
p_profiles = ggplot(img_res$images_dt, aes(xmin = tx - spc, xmax = tx + spc, 
                                           ymin = ty - spc, ymax = ty + spc, 
                                           image = png_file)) + geom_image.rect()
ggsave("fig3_profiles.pdf", 
       cowplot::plot_grid(p_profiles, p_density), 
       width = 12, height = 6)
plot_tsne_img_byCell(img_res$images_dt, tsne_dt = img_res$tsne_dt[grepl("CD34", cell)], n_points = n_points, min_size = .05)

cell_a = "H7"
cell_b = "CD34"

p = plot_tsne_img_byCell(img_res$images_dt, 
                         tsne_dt = img_res$tsne_dt[cell %in% c(cell_a, cell_b)], 
                         N_floor = 0, 
                         # N_ceiling = 300, 
                         n_points = n_points, min_size = 0)
ggsave("fig3_sideBySide.pdf", p$plot, width = 8, height = 4)

tsne_res$btx = mybin(tsne_res$tx, n_points = n_points)
tsne_res$bty = mybin(tsne_res$ty, n_points = n_points)

plot_velocity_arrows = function(tsne_res, cell_a, cell_b,
                                p = NULL,
                                id_to_plot = NULL,
                                max_plotted = 500,
                                min_delta = .1,
                                angle.min = 0, angle.max = 360){
    v_dt = calc_delta(tsne_res, cell_a, cell_b, n_points)$velocity_dt
    if(is.null(id_to_plot)){
        v_dt.tp = copy(v_dt)
    }else{
        v_dt.tp = v_dt[id %in% id_to_plot]    
    }
    v_dt.tp = v_dt.tp[id %in% sampleCap(v_dt.tp$id, max_plotted)]
    
    v_dt.tp[, angle := xy2deg(x1 = tx_cell_a, x2 = tx_cell_b, y1 = ty_cell_a, y2 = ty_cell_b)]
    if(angle.min > angle.max){
        v_dt.tp[, foreground := angle <= angle.min & angle >= angle.max]
    }else{
        v_dt.tp[, foreground := angle >= angle.min & angle <= angle.max]
    }

    bins = 36
    v_dt.tp[, angle_bin := ceiling((angle)/(360/bins))]
    b_dt = v_dt.tp[, .N, angle_bin]
    p_key = ggplot(b_dt, aes(x = angle_bin, y = N, fill = angle_bin)) + 
        geom_bar(width = 1, stat = "identity") + coord_polar() +
        scale_fill_gradientn(colours = c("orange", "red", "purple", "blue", 
                                          "green", "orange"), limits = c(0, 360)/(360/bins), 
                             breaks = 0:4*90/(360/bins), labels = function(x)x*360/bins) + 
        scale_x_continuous(labels = function(x)x*360/bins, breaks = 1:4*90/(360/bins))
    p_key
    bg = v_dt.tp[foreground == FALSE]
    # v_dt.tp[, grp1 := tx_cell_a > tx_cell_b]
    # v_dt.tp[, grp2 := ty_cell_a > ty_cell_b]
    if(is.null(p)) p = ggplot()
    p_arrows = p + 
        labs(title = paste("from", cell_a, "to", cell_b), subtitle = "color mapped to angle") +
        annotate("segment", x = bg$tx_cell_a, xend = bg$tx_cell_b, y = bg$ty_cell_a, yend = bg$ty_cell_b, color = "lightgray") +
        geom_segment(data = v_dt.tp[foreground == TRUE], aes(x = tx_cell_a, xend = tx_cell_b, 
                                         y = ty_cell_a, yend = ty_cell_b, 
                                         color = angle), 
                     arrow = arrow(length = unit(0.1,"cm"))) + 
        
        scale_color_gradientn(colours = c("orange", "red", "purple", "blue", 
                                          "green", "orange"), limits = c(0, 360), breaks = 0:4*90) 
    list(p_arrows, p_key)
}

plot_velocity_arrows(tsne_res, "H7", "CD34", id_to_plot = tp, angle.min = 0, angle.max = 360)
plot_velocity_arrows(tsne_res, "CD34", "H7")
plot_velocity_arrows(tsne_res, "CD34", "Kasumi1")

delta_res = calc_delta(tsne_res, cell_a, cell_b, n_points)

v_dt = delta_res$velocity_dt
v_dt.tp = v_dt[id %in% tp]

v_dt.tp[, angle := xy2deg(x1 = tx_cell_a, x2 = tx_cell_b, y1 = ty_cell_a, y2 = ty_cell_b)]
v_dt.tp[, grp1 := tx_cell_a > tx_cell_b]
v_dt.tp[, grp2 := ty_cell_a > ty_cell_b]

p_arrows = ggplot(v_dt.tp, aes(x = tx_cell_a, xend = tx_cell_b, 
                               y = ty_cell_a, yend = ty_cell_b, 
                               color = angle)) + 
    labs(title = paste("from", cell_a, "to", cell_b), subtitle = "color mapped to angle") +
    geom_segment(arrow = arrow(length = unit(0.1,"cm"))) + 
    scale_color_gradientn(colours = c("orange", "red", "purple", "blue", 
                                      "green", "orange"), limits = c(0, 360), breaks = 0:4*90) #+ facet_wrap("grp1~grp2")
ggsave("fig3_arrows.pdf")

av_dt = delta_res$agg_velocity_dt

p_velocity = ggplot() + 
    annotate("point", x = tsne_res.tp$tx, y = tsne_res.tp$ty, color = "lightgray") +
    geom_point(data = tsne_res.tp[cell %in% c(cell_a, cell_b)], aes(x = tx, y = ty, color = cell)) +
    geom_segment(data = av_dt[N > 6], aes(x = tx_cell_a, xend = tx_cell_b, y = ty_cell_a, yend = ty_cell_b, size = N), arrow = arrow()) + 
    coord_cartesian(xlim = range(tsne_res$tx), ylim = range(tsne_res$ty)) +
    labs(x = "x", y = "y") +
    scale_size_continuous(range = c(.5, 2), breaks = range(av_dt$N)) + theme_classic()

pg = cowplot::plot_grid(p_basic + labs(title = "t-sne: 2 ChIP-seq mark, 14 cell lines, 386 sites"), 
                        p_density + labs(title = "t-sne: profile frequency"), 
                        p_velocity + labs(title = "t-sne: changes between two cell lines"), nrow = 1, rel_widths = c(1,1, 1.3))
ggsave("fig3_changes.pdf", pg, width = 14, height = 4)



pr_gr = pr_gr_cd34
pr_gr = resize(pr_gr, 2000, fix = "center")
# pr_gr = reduce(pr_gr)
hist(width(pr_gr))
pr_dt = ssvFetchGRanges(list(prom2k = pr_gr), tsne_input$query_gr, 
                        return_data.table = TRUE, 
                        target_strand = "*",
                        win_size = 50, 
                        win_method = "summary")
pr_dt = pr_dt[order(x)]
pr_dt$x = round(pr_dt$x, 3)
# pr_dt = pr_dt[, .(y = sum(y)), by = .(id, x, sample)]
pr_dt[y > 1, y := 1]
p_pr_global = ggplot(pr_dt[, .(y = mean(y)), by = .(x)], aes(x = x, y = y)) + geom_path()
# p_pr_global
# ssvSignalHeatmap(pr_dt)

# pr_img_res = make_tsne_img(
#     bw_dt = pr_dt, apply_norm = FALSE,
#     tdt = tsne_res, #force_rewrite = TRUE, 
#     n_points = n_points, line_colors = c("signal" = "black")
# )


pdf("fig3_tss_plot.pdf", width = 8, height = 4)
make_tss_plot("H7")
make_tss_plot("CD34")
dev.off()
#


# plot(pr_img_res$tsne_dt[cell == "H7"]$bx,
# img_res$tsne_dt[cell == "H7"]$bx)
