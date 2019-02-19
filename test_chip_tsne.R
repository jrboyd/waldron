library(peakrefine)
library(BiocFileCache)
library(seqsetvis)
library(GenomicRanges)
library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)
library(cowplot)
library(Rtsne)
setwd("~/R/waldron/")
source("setup_files.R")
library(ggimage)
library(BiocFileCache)
bfcif = peakrefine::bfcif
bfc = BiocFileCache()

load("CD34_bivalent.save")
k4_bw_hESC = c("H7_H3K4me3" = "/slipstream/galaxy/uploads/working/qc_framework/output_bivalency_redo_patients_H7/H7_H3K4ME3_pooled/H7_H3K4ME3_pooled_FE.bw")
k27_bw_hESC = c("H7_H3K27me3" = "/slipstream/galaxy/uploads/working/qc_framework/output_bivalency_redo_patients_H7/H7_H3K27ME3_pooled/H7_H3K27ME3_pooled_FE.bw")
k4_bws  = c(k4_bws, k4_bw_hESC)
k27_bws = c(k27_bws, k27_bw_hESC)

view_size = 3e3

# function(qgr, view_size, bws)

qgr = resize(biv_gr, view_size, fix = "center")
qbw = c(k4_bws, k27_bws)
qdt = data.table(qbw = qbw)
qdt[, c("cell", "mark") := tstrsplit(basename(qbw), "_", keep = 1:2)]
qdt[, mark := sub("ME3", "me3", mark)]
stopifnot(length(unique(table(qdt$cell))) == 1)
stopifnot(length(unique(table(qdt$mark))) == 1)

qwin = 50
qmet = "summary"
rname = digest::digest(list(qgr, qdt, qwin, qmet))

# bw_dt = bfcif(bfc, rname, function(){
#     ssvFetchBigwig(qbw, qgr, 
#                    return_data.table = TRUE, win_method = qmet, win_size = qwin)
# })
bw_dt$x = round(bw_dt$x, 3)

bw_dt = bw_dt[, .(id, x, y, sample)]
bw_dt[, c("cell", "mark") := tstrsplit(sample, "_")]
bw_dt$sample = NULL

#if capped
cap = 20
bw_dt[y > cap, y := cap]

#if flip (inappropriate if qgr is inherently stranded, say TSSes)
balance_dt = bw_dt[, .(right_sum = sum(y[x > 0]), left_sum = sum(y[x < 0])), by = .(cell, id)]
balance_dt = balance_dt[, .(needs_flip = left_sum > right_sum, cell, id)]
bw_dt = merge(bw_dt, balance_dt)
bw_dt[needs_flip == TRUE, x := -x]
bw_dt$needs_flip = NULL

k4_dt = bw_dt[mark == "H3K4me3"]
k4_mat = dcast(k4_dt, id+cell~x, value.var = "y")
k27_dt = bw_dt[mark == "H3K27me3"]
k27_mat = dcast(k27_dt, id+cell~x, value.var = "y")

rn = paste(k4_mat$cell, k4_mat$id)
rn2 = paste(k27_mat$cell, k27_mat$id)
stopifnot(all(rn == rn2))
mat = cbind(as.matrix(k4_mat[, -1:-2]), as.matrix(k27_mat[, -1:-2])*5)



perp = 100
tsne_rname = digest::digest(list(qgr, qbw, qwin, qmet, perp, head(mat)))
res_tsne = bfcif(bfc, tsne_rname, function(){
    Rtsne(mat, num_threads = 20, perplexity = 100)
})

plot(res_tsne$Y, pch = 16, col = rgb(0,0,0,.2))

npoints = 20
min2agg = 5

rname2 = digest::digest(list(qgr, qbw, qwin, qmet, res_tsne))
# unique(mdt[, .(bx, by)])
odir = file.path("tsne_components", rname2)
dir.create(odir, recursive = TRUE, showWarnings = FALSE)

if(TRUE){ #aggregation method
    tdt = as.data.table(res_tsne$Y)
    colnames(tdt) = c("tx", "ty")
    
    norm1 = function(x){
        (x - min(x)) / (max(x) - min(x))
    }
    mybin = function(x, npoints){
        floor(norm1(x) * (npoints-.00001))+1
    }
    tdt[, bx := mybin(tx, npoints)]
    tdt[, by := mybin(ty, npoints)]
    tdt$rn = rn
    tdt[, c("cell", "id") := tstrsplit(rn, " ")]
    tdt$cell = factor(tdt$cell)
    tdt$cell = factor(tdt$cell, levels = levels(tdt$cell)[c(7, 1:4, 8, 5:6, 9:14)])
    
    tdt[, N := .N, .(bx, by)]
    tdt = tdt[N >= min2agg]
    tdt$N = NULL
    
    mdt = merge(bw_dt, tdt[, .(bx, by, cell, id)])
    mdt = mdt[, .(y = mean(y)), .(bx, by, x, mark)]
    mdt[, plot_id := paste(bx, by, sep = "_")]
    
    
    
    img_dt = unique(mdt[, .(bx, by, plot_id)])
    img_dt[, png_file := file.path(odir, paste0(plot_id, ".png"))]
    
    xrng = range(res_tsne$Y[,1])
    yrng = range(res_tsne$Y[,2])
    xspc = diff(xrng)/npoints/2
    yspc = diff(yrng)/npoints/2
    xs = seq(min(xrng)+xspc, max(xrng)-xspc, diff(xrng)/(npoints))
    ys = seq(min(yrng)+yspc, max(yrng)-yspc, diff(yrng)/(npoints))
    
    plot(expand.grid(xs, ys), xlim = xrng, ylim = yrng)
    rect(min(xrng), min(yrng), max(xrng), max(yrng), col = rgb(0,0,1,.1))
    
    img_dt[, tx := xs[bx]]
    img_dt[, ty := ys[by]]
    
    mdt[, ynorm := y / quantile(y, .95), by = .(mark)]
    mdt[ynorm > 1, ynorm := 1]
    # for(i in seq_len(nrow(img_dt))){
    options(mc.cores = 20)
    if(any(!file.exists(img_dt$png_file))){
        hidden = parallel::mclapply(which(!file.exists(img_dt$png_file)), function(i){
            fpath = img_dt$png_file[i]
            p_id = img_dt$plot_id[i]
            pdt = mdt[plot_id == p_id]
            # pdt[, ysm := seqsetvis:::movingAverage(y, n = 8), by = .(mark)]
            pdt[, ysm := seqsetvis:::movingAverage(ynorm, n = 8), by = .(mark)]
            pdt = applySpline(pdt, n = 10, by_ = "mark", y_ = "ysm")
            p = ggplot(pdt, aes(x = x, y = ysm, color = mark)) + 
                geom_path(size = .6) +
                scale_color_manual(values = c("H3K4me3" = "forestgreen", "H3K27me3" = "firebrick1")) +
                theme_void() + guides(color = "none") + 
                coord_cartesian(ylim = c(0, 1))
            p
            ggsave(fpath, p, width = 2, height = 2, units = "cm")
            p
        })
    }
    
    p_describe = ggplot(img_dt, aes(x = tx, y = ty, image = png_file)) + geom_image()
    

    img_dt = merge(img_dt, tdt[, .N, .(bx, by)])
    img_dt[, img_size := (N - min(N))]
    img_dt[, img_size := img_size / max(img_size)]
    img_dt[, img_size := img_size / npoints]
    p_dens = ggplot(img_dt, aes(x = tx, y = ty, image = png_file)) + geom_image(size = img_dt$img_size) +
        # annotate("tile", x = img_dt$tx, y = img_dt$ty, size = img_dt$img_size, color = "black", fill = NA)
        annotate("rect", 
                 xmin = img_dt$tx - xspc * img_dt$img_size * npoints,
                 xmax = img_dt$tx + xspc * img_dt$img_size * npoints,
                 ymin = img_dt$ty - yspc * img_dt$img_size * npoints,
                 ymax = img_dt$ty + yspc * img_dt$img_size * npoints, 
                 # size = img_dt$img_size, 
                 color = "black", fill = NA)
    # p_describe
    # ggsave("tsne_summary.pdf", plot = p_describe, width = 8, height = 8)
    
   
    cols = safeBrew(length(unique(tdt$cell)))
    
    p1 = ggplot(tdt, aes(x = tx, y = ty)) + geom_point() + 
        labs(title = "t-sne of 381 CD34+ bivalent sites in all cell lines")
    
    p2 = p_describe + 
        labs(title = "t-sne averaged signals")
    
    
    
    p3 = ggplot() + geom_image(data = img_dt, aes(x = tx, y = ty, image = png_file)) + 
        geom_point(data = tdt[cell == "CD34-01536" | cell == "H7"], aes(x = tx, y = ty, color = cell)) + 
        scale_color_discrete(drop = FALSE) +
        labs(title = "CD34-01536 and H7 highlighted")
    
    p4 = ggplot() + geom_image(data = img_dt, aes(x = tx, y = ty, image = png_file)) + 
        geom_point(data = tdt[cell == "Nalm6" | cell == "Kasumi1"], aes(x = tx, y = ty, color = cell)) + 
        scale_color_discrete(drop = FALSE) + 
        labs(title = "Nalm6 and Kasumi1 highlighted")
    
    # ggplot() + geom_image(data = img_dt, aes(x = tx, y = ty, image = png_file)) + 
    #     geom_point(data = tdt[cell == "H7"], aes(x = tx, y = ty, color = cell)) + 
    #     scale_color_discrete(drop = FALSE)
    
    
    
    
    
    p5 = ggplot(tdt, aes(x = tx, y = ty, color = cell)) + geom_point(size = 1) + facet_wrap("cell") + 
        theme(panel.grid.major = element_line(color = "darkgray")) + 
        labs("t-sne - every cell line", subtitle = "scatterplot")
    p6 = ggplot(tdt, aes(x = tx, y = ty, color = cell)) + geom_density2d(bins = 5) + 
        facet_wrap("cell") + 
        theme(panel.grid.major = element_line(color = "darkgray"))  + 
        labs("t-sne - every cell line", subtitle = "density contours")
    p7 = ggplot(tdt, aes(x = tx, y = ty)) + 
        stat_density_2d(aes(fill = stat(level)), geom = "polygon", bins = 50) + 
        facet_wrap("cell") + 
        theme(panel.grid.major = element_line(color = "darkgray")) + 
        scale_fill_viridis_c() + 
        labs("t-sne - every cell line", subtitle = "density fill")
    p8 = ggplot(tdt, aes(x = tx, y = ty, color = cell)) + 
        stat_density_2d(geom = "point", 
                        aes(size = stat(density)), 
                        n = 20, contour = FALSE) + 
        facet_wrap("cell") + 
        scale_size(range = c(.02, 2), limits = c(.00005, NA)) + 
        labs("t-sne - every cell line", subtitle = "density weighted points")
    
    pdf("tsne_report3.pdf")
    print(p1)
    print(p2)
    print(p_dens)
    print(p3)
    print(p4)
    print(p5)
    print(p6)
    print(p7)
    print(p8)
    dev.off()
}
# }
if(FALSE){ #sampling method
    
    
    tdt$png_file = file.path(odir, paste0(paste(mdt$bx, mdt$by, sep = "_"), ".png"))
    
    
    ex = sample(seq_len(nrow(res_tsne$Y)), 500)
    
    xrng = range(res_tsne$Y[,1])
    yrng = range(res_tsne$Y[,2])
    
    xs = seq(min(xrng), max(xrng), diff(xrng)/npoints)
    ys = seq(min(yrng), max(yrng), diff(yrng)/npoints)
    xy_dt = as.data.table(base::expand.grid(xs, ys))
    rdt
    
    # res_nn = RANN::nn2(xy_dt, res_tsne$Y, k = 1)
    
    rdt = data.table(tx = res_tsne$Y[ex,1], 
                     ty = res_tsne$Y[ex,2], 
                     pid = rn[ex])
    rdt[, c("cell", "id") := tstrsplit(pid, " ")]
    
    
    
    rdt$png_file = file.path(odir, paste0(paste(rdt$id, rdt$cell, sep = "_"), ".png"))
    
    if(!all(file.exists(rdt$png_file))){
        hidden = pbapply::pbsapply(which(!file.exists(rdt$png_file)), function(i){
            pcell = rdt$cell[i]
            pid = rdt$id[i]
            pfile = rdt$png_file[i]
            pdt = bw_dt[id == pid & cell == pcell]
            pdt[, ysm := seqsetvis:::movingAverage(y, n = 8), by = .(mark)]
            pdt = applySpline(pdt, n = 10, by_ = "mark", y_ = "ysm")
            p = ggplot(pdt, aes(x = x, y = ysm, color = mark)) + 
                geom_path(size = .6) +
                scale_color_manual(values = c("H3K4me3" = "forestgreen", "H3K27me3" = "firebrick1")) +
                theme_void() + guides(color = "none")
            p
            ggsave(pfile, p, width = 2, height = 2, units = "cm")
            
        })
    }
    
    
    ggplot(rdt, aes(x = tx, y = ty, image = png_file)) + geom_image()
}
