library(BiocFileCache)
library(ggimage)
library(Rtsne)
library(magrittr)
setwd("~/R/waldron/")
source("setup_files.R")
source("functions_tsne.R")
source("geom_image.rect.R")
bfcif = peakrefine::bfcif
bfc = BiocFileCache()
options("mc.cores" = 36)
load("CD34_bivalent.save")
load("CD34_consensus_v1.save")

ref_gr = rtracklayer::import.gff("~/gencode.v28.annotation.gtf.gz", format = "gtf", feature.type = "transcript")
trans_gr = subset(ref_gr, transcript_support_level %in% 1:2)
promoters(trans_gr, 1, 1)

view_size = 15e3
n_points = 16
theme_set(theme_classic())
peaks = easyLoad_narrowPeak(c(k4_peaks, k27_peaks))
set.seed(0)

# qgr = biv_gr
# qgr = sample(peaks[["CD34-01562_H3K4me3"]], 15000)
# qgr = peaks[["CD34-01562_H3K4me3"]]

# qgr = reduce(c(k4_consenus, subsetByOverlaps(k27_consenus, k4_consenus, invert = TRUE)))
# qgr = reduce(c(peaks[["H7_H3K4me3"]], subsetByOverlaps(peaks[["H7_H3K27me3"]], peaks[["H7_H3K4me3"]], invert = TRUE)))
qgr = promoters(trans_gr, 1, 1) %>% reduce %>% resize(., view_size, fix = "center")
# qgr = sample(qgr, 15000)

qgr = resize(qgr, view_size, fix = "center")
qbw = c(k4_bws, k27_bws)
qdt = data.table(qbw = qbw)
qdt[, c("cell", "mark") := tstrsplit(basename(qbw), "_", keep = 1:2)]
qdt[, mark := sub("ME3", "me3", mark)]

stopifnot(length(unique(table(qdt$cell))) == 1)
stopifnot(length(unique(table(qdt$mark))) == 1)

message("fetch tsne input")
tsne_input = fetch_tsne_mat(qdt, qgr, 
                            qwin = 50, 
                            qmet = "summary", 
                            cap_value = 30, 
                            high_on_right = FALSE)
message("run tsne")
tsne_res = run_tsne(tsne_input$tsne_mat, perplexity = 100)

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

p_profiles = ggplot(img_res$images_dt, aes(x = tx, y = ty, image = png_file)) + geom_image()

ggplot(img_res$images_dt, aes(xmin = tx - .05, xmax = tx + .05, ymin = ty - .05, ymax = ty + .05, image = png_file)) + geom_image.rect()

p_density = plot_tsne_img(img_res$images_dt, n_points = n_points, 
                          N_ceiling = NULL, N_floor = 50, 
                          show_plot = TRUE)$plot
p_density + theme(panel.background = element_rect(fill= "lightblue"))
cowplot::plot_grid(p_profiles, p_density)
# plot_tsne_img_byCell(img_res$images_dt, tsne_dt = img_res$tsne_dt, N_ceiling = 30, n_points = n_points, min_size = .05)
plot_tsne_img_byCell(img_res$images_dt, tsne_dt = img_res$tsne_dt[grepl("CD34", cell)], n_points = n_points, min_size = .05)

cell_a = "H7"
cell_b = "CD34-01517"

p = plot_tsne_img_byCell(img_res$images_dt, 
                         tsne_dt = img_res$tsne_dt[cell %in% c(cell_a, cell_b)], 
                         N_floor = 0, 
                         # N_ceiling = 300, 
                         n_points = n_points, min_size = 0)

ggsave("tmp_sideBySide.pdf", p$plot, width = 8, height = 4)

delta_res = calc_delta(tsne_res, cell_a, cell_b, n_points)

v_dt = delta_res$velocity_dt
v_dt.tp = v_dt[id %in% tp]
xy2deg = function(x1, y1, x2, y2){
    x = x2 - x1
    y = y2 - y1
    deg = atan(y/x) * 180 / pi + 180
    deg[x < 0] = deg[x < 0] + 180
    deg[deg > 360] = deg[deg > 360] - 360
    deg
}
v_dt.tp[, angle := xy2deg(x1 = tx_cell_a, x2 = tx_cell_b, y1 = ty_cell_a, y2 = ty_cell_b)]
v_dt.tp[, grp1 := tx_cell_a > tx_cell_b]
v_dt.tp[, grp2 := ty_cell_a > ty_cell_b]
p_arrows = ggplot(v_dt.tp, aes(x = tx_cell_a, xend = tx_cell_b, 
                    y = ty_cell_a, yend = ty_cell_b, 
                    color = angle)) + 
    labs(title = "color mapped to angle") +
    geom_segment(arrow = arrow(length = unit(0.1,"cm"))) + 
    scale_color_gradientn(colours = c("red", "purple", "blue", 
                                      "green", "yellow", "orange"), limits = c(0, 360), breaks = 1:4*90) #+ facet_wrap("grp1~grp2")
ggsave("tmp_arrows.pdf")

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
ggsave("tmp_changes.pdf", pg, width = 14, height = 4)



pr_gr = promoters(ref_gr, 1000, 1000)
pr_gr = reduce(pr_gr)
hist(width(pr_gr))
pr_dt = ssvFetchGRanges(list(prom2k = pr_gr), tsne_input$query_gr, 
                        return_data.table = TRUE, 
                        target_strand = "*",
                        win_size = 50, 
                        win_method = "summary")
pr_dt = pr_dt[order(x)]
pr_dt$x = round(pr_dt$x, 3)
# pr_dt = pr_dt[, .(y = sum(y)), by = .(id, x, sample)]
# pr_dt[y > 1, y := 1]
p_pr_global = ggplot(pr_dt[, .(y = mean(y)), by = .(x)], aes(x = x, y = y)) + geom_path()
p_pr_global
ssvSignalHeatmap(pr_dt)

# pr_img_res = make_tsne_img(
#     bw_dt = pr_dt, apply_norm = FALSE,
#     tdt = tsne_res, #force_rewrite = TRUE, 
#     n_points = n_points, line_colors = c("signal" = "black")
# )
make_tss_plot = function(qcell){
    xrng = range(tsne_res$tx)
        yrng = range(tsne_res$ty)
    
    pr_img_res = make_tsne_img(
        bw_dt = pr_dt, apply_norm = FALSE, 
        tdt = tsne_res[cell == qcell], #force_rewrite = TRUE, 
        xrng = xrng,
        yrng = yrng,
        n_points = n_points, line_colors = c("signal" = "black")
    )
    p_pr_density = plot_tsne_img_byCell(pr_img_res$images_dt, 
                                        pr_img_res$tsne_dt[cell == qcell],
                                        n_points = n_points, N_ceiling = NULL)$plot +
        coord_cartesian(xlim = xrng, ylim = yrng)
    p_h7_density = plot_tsne_img_byCell(img_res$images_dt, 
                                        img_res$tsne_dt[cell == qcell], 
                                        n_points = n_points, N_ceiling = NULL)$plot +
        coord_cartesian(xlim = xrng, ylim = yrng)
    
    pg = cowplot::plot_grid(p_h7_density + labs(title = paste(qcell, "k4me3+k27me3")), 
                       p_pr_density + labs(title = paste(qcell, "tss frequency")))
    ggsave(paste0("tmp_", qcell, "_tss.pdf"), plot = pg, width = 8, height = 4)
    
    # head(pr_img_res$tsne_dt[cell == cell, .N, by = .(bx, by)][order(N, decreasing = TRUE)])
    # head(img_res$tsne_dt[cell == cell, .N, by = .(bx, by)][order(N, decreasing = TRUE)])
}

make_tss_plot("H7")
make_tss_plot("CD34-01562")

#


# plot(pr_img_res$tsne_dt[cell == "H7"]$bx,
# img_res$tsne_dt[cell == "H7"]$bx)
