#MAKE_FIGURES
# Figure1: venn of A) H7 H3K4me3/27me3 B) CD34+ H3K4me3/K27me3 C) H7 bv/CD34+ bv
# Figure2: Current figure of CD34+ bv-marked genes; snapshot of a region
# Figure3: comparison to several (not necessary all?) cell lines via clustering signal.
# Figure4: Kasumi and another leukemia lineage (PDX2?) maintain most but not all bv domains. BV domains lost + retained pathways.
# Figure5: prognostication of genes
source('setup_files.R')
source("functions.R")
source("functions_tsne.R")


tx_gr = rtracklayer::import.gff("~/gencode.v28.annotation.gtf.gz", format = "gtf", feature.type  = "transcript")
tx_gr = subset(tx_gr, transcript_support_level %in% 1:2)
pr_gr = promoters(tx_gr, 1, 1)
pr_gr = split(pr_gr, pr_gr$gene_name)
# pr_gr = resize(unlist(reduce(resize(pr_gr, 200, fix = "center"))), 1, fix = "center")
pr_gr = resize(unlist(reduce(resize(pr_gr, 1000, fix = "center"))), 1, fix = "center")


# view_size = 15e3
view_size = 6e3
# n_points = 20
n_points = 16
theme_set(theme_classic())
set.seed(0)

np_grs = easyLoad_narrowPeak(c(k4_peaks["H7_H3K4me3"], k4_peaks[grepl("CD34", k4_peaks)],
                               k27_peaks["H7_H3K27me3"], k27_peaks[grepl("CD34", k4_peaks)]))
#CD34 k4
olaps_cd34_k4 = ssvOverlapIntervalSets(np_grs[grepl("CD34", names(np_grs)) & grepl("K4", names(np_grs))])                    
is_k4_consensus = rowSums(as.data.frame(mcols(olaps_cd34_k4))) > 2
olaps_cd34_k4 = olaps_cd34_k4[is_k4_consensus]

#CD34 k27
olaps_cd34_k27 = ssvOverlapIntervalSets(np_grs[grepl("CD34", names(np_grs)) & grepl("K27", names(np_grs))])                    
is_k4_consensus = rowSums(as.data.frame(mcols(olaps_cd34_k27))) > 2
olaps_cd34_k27 = olaps_cd34_k27[is_k4_consensus]

olaps_cd34 = ssvOverlapIntervalSets(list("CD34_H3K4me3" = olaps_cd34_k4, "CD34_H3K27me3" = olaps_cd34_k27))
pr_gr_cd34 = subsetByOverlaps(pr_gr, olaps_cd34)
pr_gr_cd34$id = paste0("tss_", seq_along(pr_gr_cd34))
pr_gr_cd34$gene_name = names(pr_gr_cd34)
names(pr_gr_cd34) = NULL

#prep tsne args
qgr = pr_gr_cd34
qgr = resize(qgr, view_size, fix = "center")
qbw = c(k4_bws, k27_bws)
qdt = data.table(qbw = qbw)
qdt[, c("cell", "mark") := tstrsplit(basename(qbw), "_", keep = 1:2)]
qdt[, mark := sub("ME3", "me3", mark)]
qdt[grepl("CD34", cell), cell := "CD34"]
qdt$norm_factor = 1
qdt[mark == "H3K4me3"]$norm_factor = 1/5
qdt
# qdt[cell == "mm1s"]$qbw = rev(qdt[cell == "mm1s"]$qbw)

# stopifnot(length(unique(table(qdt$cell))) == 1)
stopifnot(length(unique(table(qdt$mark))) == 1)

#retrieve tidy profiles and cast to wide matrix for tsne
message("fetch tsne input")
options(mc.cores = 32)
tsne_input = fetch_tsne_mat(qdt, qgr, 
                            qwin = 50, 
                            qmet = "summary", 
                            cap_value = 30, 
                            high_on_right = FALSE,
                            force_overwrite = FALSE)

#run tsne
message("run tsne")
tsne_res = run_tsne(tsne_input$tsne_mat, perplexity = 100, force_overwrite = FALSE)
tsne_res$cell = factor(tsne_res$cell)
tsne_res$cell = factor(tsne_res$cell, levels = levels(tsne_res$cell)[c(4, 1, 5, 2:3, 6:12)])

#preset sampled id to minimize points plotted
tp = sample(unique(tsne_res$id), min(500, length(unique(tsne_res$id))))
tsne_res.tp = tsne_res[id %in% tp]


#configure for bam retrieval
bam_qdt = data.table(filepath = bams, sample = names(bams))
bam_qdt[, c("cell", "mark") := tstrsplit(sample, "_")]
bam_qdt[, cell := sub("CD34.+", "CD34", cell)]
# bam_qdt[cell == "mm1s" & mark != "input"]$filepath = rev(bam_qdt[cell == "mm1s" & mark != "input"]$filepath)
# bam_qdt = bam_qdt[mark == "input"]

bams_hESC = c("H7_H3K4me3" = "/slipstream/galaxy/uploads/working/qc_framework/output_bivalency_redo_patients_H7/H7_H3K4ME3_pooled/H7_H3K4ME3_pooled.bam",
              "H7_H3K27me3" = "/slipstream/galaxy/uploads/working/qc_framework/output_bivalency_redo_patients_H7/H7_H3K27ME3_pooled/H7_H3K27ME3_pooled.bam",
              "H7_input" = "/slipstream/galaxy/uploads/working/qc_framework/output_bivalency_redo_patients_H7/H7_input_pooled/H7_input_pooled.bam")
bam_qdt2 = data.table(filepath = bams_hESC, sample = names(bams_hESC))
bam_qdt2[, c("cell", "mark") := tstrsplit(sample, "_")]

bam_pile_dt = bfcif(bfc, #force_overwrite = TRUE,
                    digest::digest(list("input_bams", 
                                        bam_qdt, tsne_input$query_gr)), 
                    function(){
                        ssvFetchBam(bam_qdt,
                                    tsne_input$query_gr,
                                    return_data.table = TRUE,
                                    target_strand = "*",
                                    fragLens = NA,
                                    win_size = 50,
                                    win_method = "summary", 
                                    n_cores = 20)
                    })

bam_pile_dt2 = bfcif(bfc, #force_overwrite = TRUE,
                     digest::digest(list("input_bams", 
                                         bam_qdt2, tsne_input$query_gr)), 
                     function(){
                         ssvFetchBam(bam_qdt2,
                                     tsne_input$query_gr,
                                     return_data.table = TRUE,
                                     target_strand = "*",
                                     fragLens = NA,
                                     win_size = 50,
                                     win_method = "summary", 
                                     n_cores = 20)
                     })

# bam_stranded_dt = bfcif(bfc, #force_overwrite = TRUE,
#                     digest::digest(list("input_bams", 
#                                         bam_qdt, tsne_input$query_gr)), 
#                     function(){
#                         ssvFetchBam(bam_qdt,
#                                     tsne_input$query_gr,
#                                     return_data.table = TRUE,
#                                     target_strand = "both",
#                                     fragLens = NA,
#                                     win_size = 50,
#                                     win_method = "summary", 
#                                     n_cores = 20)
#                     })
# 
# bam_stranded_dt2 = bfcif(bfc, #force_overwrite = TRUE,
#                      digest::digest(list("input_bams", 
#                                          bam_qdt2, tsne_input$query_gr)), 
#                      function(){
#                          ssvFetchBam(bam_qdt2,
#                                      tsne_input$query_gr,
#                                      return_data.table = TRUE,
#                                      target_strand = "both",
#                                      fragLens = NA,
#                                      win_size = 50,
#                                      win_method = "summary", 
#                                      n_cores = 20)
#                      })

#retrieve bam sizes
bam_qdt$count = unlist(
    bfcif(bfc, 
          digest::digest(list("bam counts",
                              bam_qdt$filepath)),
          function(){
              parallel::mclapply(bam_qdt$filepath, function(x){
                  Rsamtools::countBam(x)$records
              })
          })
)

bam_qdt2$count = unlist(
    bfcif(bfc, 
          digest::digest(list("bam counts",
                              bam_qdt2$filepath)),
          function(){
              parallel::mclapply(bam_qdt2$filepath, function(x){
                  Rsamtools::countBam(x)$records
              })
          })
)

#agg CD34
bam_pile_dt = bam_pile_dt[, .(y = sum(y)), .(id, x, cell, mark)]
#depth normalize
bam_pile_dt = merge(bam_pile_dt, bam_qdt[, .(count = sum(count)), .(cell, mark)])
bam_pile_dt = bam_pile_dt[, .(id, x, y = y / count * 1e6, cell, mark)]
bam_pile_dt2 = merge(bam_pile_dt2, bam_qdt2[, .(count = sum(count)), .(cell, mark)])
bam_pile_dt2 = bam_pile_dt2[, .(id, x, y = y / count * 1e6, cell, mark)]
#combine
bam_pile_dt = rbind(bam_pile_dt, bam_pile_dt2)
remove(bam_pile_dt2)
#trim outliers
qcap = quantile(bam_pile_dt[, max(y), by = .(id, cell)]$V1, .95)
bam_pile_dt[y > qcap, y := qcap]
#view mean profiles
mean_dt = bam_pile_dt[, .(y = mean(y)), .(cell, mark, x)]
ggplot(mean_dt, 
       aes(x = x, y = y, color = mark)) + 
    geom_path() + facet_grid("mark~cell", scales = "free_y")

tmp1 = mean_dt[mark != "H3K27me3"]
tmp2 = mean_dt[mark != "H3K4me3"]
tmp1$group = "H3K4me3"
tmp2$group = "H3K27me3"
mean2_dt = rbind(tmp1, tmp2)
remove(tmp1)
remove(tmp2)
mean2_dt$group = factor(mean2_dt$group, levels = c("H3K4me3", "H3K27me3"))
ggplot(mean2_dt, 
       aes(x = x, y = y, color = mark)) + 
    geom_path() + facet_grid("group~cell", scales = "free_y") +
    labs(title = "rpm normalized pileups at tsses", 
         subtitle = paste0("view is ", view_size / 1e3, "kb window centered on tss"), 
         x = paste0(view_size / 1e3, "kb window"), y = "rpm") +
    theme(axis.line.x = element_blank(), 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank())
ggsave("qc_cell_pileups.pdf", width = 11, height = 6)

piles_img_res = make_tsne_img(profiles_dt = bam_pile_dt, 
                              position_dt = tsne_res, #force_rewrite = TRUE, 
                              apply_norm = FALSE, 
                              
                              ylim = c(0,.4),  
                              # xrng = zoom_x,
                              # yrng = zoom_y,
                              n_points = n_points, 
                              line_colors = c("input" = "blue", "H3K4me3" = "forestgreen", "H3K27me3" = "red")
)

p_basic = make_img_plots(img_results = list(piles_img_res), 
               qcell = NULL, 
               min_size = 0, 
               N_ceiling = NULL,
               as_facet = FALSE)

qcells = c("CD34", "Kasumi1", "mm1s", "H7")
ver = "wH7"
# rects = list(
#     c(-.32, -.095, -.5, -.37),
#     c(-.4, -.25, -.25, -.13),
#     c(-.05, .25, -.07, .2),
#     c(-.08, .08, -.5, -.42),
#     c(-.1, .02, -.2, -.15),
#     c(.12, .3, -.27, -.15),
#     c(-.3, -.05, .05, .45),
#     c(.3, .5, .0, .25)
# )
p_alpha = .5
cell_tsne_res = tsne_res[cell %in% qcells]
cell_tsne_res$cell = factor(cell_tsne_res$cell, levels = qcells)
p_full = ggplot() + geom_point(data = cell_tsne_res, 
                               aes(x = tx, y = ty, color = cell),
                               alpha = p_alpha, shape = 16, size = .3) + 
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 2, shape = 16))) + 
    labs(x = "", y = "") +
    guides(color = "none")

rects = list(
    c(-.4, .1, -.5, -.13),
    c(-.5, -.17, -.13, .25),
    c(-.17, 0, -.13, .05),
    c(-.17, .17, .05, .5),
    c(.17, .5, -.05, .5),
    c(0, .5, -.5, -.05),
    c(0, .17, -.05, .05)
)

xr = sapply(rects, function(x)min(x[1:2]))
yr = sapply(rects, function(x)min(x[3:4]))
rects = rects[order(-yr)]
p_ann = annotate_rects(p_full, rects)
annotate_rects(p_basic[[1]], rects)
p_ann

png(paste0("regions_master.", ver, ".png"), width = 10, height = 4.5, units = "in", res = 300)
pg = cowplot::plot_grid(p_ann, p_full + facet_wrap("cell", nrow = 2), nrow = 1)
print(pg)
dev.off()

profile_dt = bam_pile_dt

for(i in seq_along(rects)){
    message(i)
    zoom_xi = rects[[i]][1:2]
    zoom_yi = rects[[i]][3:4]
    zoom_np = 12
    p_points = ggplot() + geom_point(data = cell_tsne_res, 
                                     aes(x = tx, y = ty),#, color = cell),
                                     alpha = .5, shape = 16, size = 1) + 
        guides(color = guide_legend(override.aes = list(alpha = 1, size = 2, shape = 16))) + 
        coord_cartesian(zoom_xi, zoom_yi) + 
        facet_wrap("cell", drop = FALSE) +
        labs(x = "", y = "", title = paste("t-sne region", i))
    pile_img_res = make_tsne_img(
        profiles_dt = profile_dt,#[cell %in% qcells], 
        position_dt = tsne_res[cell %in% qcells],#[cell %in% qcells], #force_rewrite = TRUE, 
        apply_norm = FALSE, 
        ylim = c(0, .8),  
        xrng = zoom_xi,
        yrng = zoom_yi,
        n_points = zoom_np, 
        line_colors = c("input" = "blue", 
                        "H3K4me3" = "forestgreen", 
                        "H3K27me3" = "red"),
        facet_by = "cell"
    )
    pile_img_res$images_dt$cell = factor(pile_img_res$images_dt$cell, levels = qcells)
    pile_img_res$tsne_dt$cell = factor(pile_img_res$tsne_dt$cell, levels = qcells)
    p_images = make_img_plots_facet(img_results = pile_img_res, 
                                    xrng = zoom_xi, 
                                    yrng = zoom_yi,
                                    qcell = NULL, 
                                    min_size = .5, 
                                    N_ceiling = quantile((pile_img_res$images_dt$N), .5)*2)+
        labs(x = "", y = "", title = paste("t-sne region", i))
    
    pg = cowplot::plot_grid(p_points, p_images, nrow = 1)
    # print(pg)
    png(paste0("regions_", i, ".", zoom_np, "np.", ver, ".png"), width = 10, height = 4.5, units = "in", res = 300)
    print(pg)
    dev.off()
}

cowplot::plot_grid(plotlist = plot_velocity_arrows(tsne_res, "CD34", "Nalm6"))
cowplot::plot_grid(plotlist = plot_velocity_arrows(tsne_res, "CD34", "mm1s"))
cowplot::plot_grid(plotlist = plot_velocity_arrows(tsne_res, "H7", "CD34"))
cowplot::plot_grid(plotlist = plot_velocity_arrows(tsne_res, "H7", "mm1s"))

plot_velocity_arrows_selected(tsne_res, qcells[1:3], "RUNX1")
plot_profiles_selected(tsne_input$bw_dt, qcells[1:3], "RUNX1")
plot_profiles_selected(bam_pile_dt, qcells[1:3], "RUNX1")

plot_velocity_arrows_selected(tsne_res, qcells[1:3], "VIM")
plot_profiles_selected(tsne_input$bw_dt, qcells[1:3], "VIM")


# sdt = pile_img_res$summary_profiles_dt
# cap = .05
# sdt[, ycap := y]
# sdt[y > cap, ycap := cap]
# sdt[y < -cap, ycap := -cap]
# hist(sdt$ynorm)
# library(GGally)
# glyph_df = GGally::glyphs(sdt, x_major = "bx", x_minor = "x", y_major = "by", y_minor = "ycap")
# ggplot(glyph_df, aes(gx, gy, group = paste(gid, mark), color = mark)) + 
#     geom_path() + 
#     scale_color_manual(values =  c("input" = "blue", 
#                                    "H3K4me3" = "forestgreen", 
#                                    "H3K27me3" = "red")) + 
#     facet_wrap("cell")

