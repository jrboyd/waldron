library(Rtsne)
library(ggimage)

annotate_rects = function(p, rects){
    anns = matrix(unlist(rects), ncol = 4, byrow = TRUE)
    # anns = anns[order(anns[,3], decreasing = TRUE),]
    p_alpha = .1
    
    # anns = data.frame(xmin = c(-.33, -.4, -.05), xmax = c(-.08, -.25, .25), ymin = c(-.5, -.25, -.07), ymax = c(-.35, -.13, .2))
    # p
    
    for(i in seq_len(nrow(anns))){
        p = p + 
            annotate("rect", 
                     xmin = anns[i, 1],
                     xmax = anns[i, 2],
                     ymin = anns[i, 3],
                     ymax = anns[i, 4], color = "black", fill = NA) +
            annotate("label", 
                     # x = (anns[i, 1] + anns[i, 2])/2,
                     # y = (anns[i, 3] + anns[i, 4])/2,
                     x = anns[i, 1],
                     y = anns[i, 4],
                     hjust = 1,
                     vjust = 1,
                     label = i)
        
        
    }
    p
}

norm1 = function(x, xrng = range(x)){
    stopifnot(length(xrng) == 2)
    (x - min(xrng)) / (max(xrng) - min(xrng))
}
mybin = function(x, n_points, xrng = range(x)){
    stopifnot(length(xrng) == 2)
    floor(norm1(x, xrng) * (n_points-.00001))+1
}

mybin_centers = function(x, n_points, xrng = range(x)){
    stopifnot(length(xrng) == 2)
    # xrng = range(x)
    xspc = diff(xrng)/n_points/2
    xs = seq(min(xrng)+xspc, max(xrng)-xspc, diff(xrng)/(n_points))
    xs
}


fetch_tsne_mat = function(qdt, 
                          qgr, 
                          qwin = 50,
                          qmet = "summary",
                          cap_value = 20,
                          high_on_right = TRUE,
                          bfc = BiocFileCache::BiocFileCache(),
                          n_cores = getOption("mc.cores", 1),
                          rname = digest::digest(list(qgr, qdt[, 1:3], qwin, qmet, cap_value, high_on_right)),
                          force_overwrite = FALSE){
    bw_dt = bfcif(bfc, rname, function(){
        ssvFetchBigwig(qdt[, 1:3], qgr, 
                       return_data.table = TRUE, 
                       win_method = qmet, 
                       win_size = qwin, n_cores = n_cores)
    }, force_overwrite = force_overwrite)   
    bw_dt$sample = NULL
    bw_dt = bw_dt[, .(y = mean(y)), .(cell, id, mark, x)]
    if(!all(qdt$norm_factor == 1)){
        norm_dt = unique(qdt[, .(cell, mark, norm_factor)])
        bw_dt = merge(bw_dt, norm_dt)
        bw_dt[, y := y * norm_factor]
        bw_dt$norm_factor = NULL
    }
    
    bw_dt[y > cap_value, y := cap_value]
    
    
    
    if(high_on_right){
        balance_dt = bw_dt[, .(right_sum = sum(y[x > 0]), left_sum = sum(y[x < 0])), by = .(cell, id)]
        balance_dt = balance_dt[, .(needs_flip = left_sum > right_sum, cell, id)]
        most_flipped = balance_dt[, .(fraction_flipped = sum(needs_flip) / .N), by = .(id)]
        most_flipped[, flip_strand := fraction_flipped > .5]
        strand(qgr) = "+"
        strand(qgr)[most_flipped$flip_strand] = "-"
        bw_dt = merge(bw_dt, balance_dt)
        remove(balance_dt)
        bw_dt[needs_flip == TRUE, x := -x]
        bw_dt$needs_flip = NULL
    }
    
    bw_dt$x = round(bw_dt$x, 3)
    
    for(m in unique(qdt$mark)){
        if(!exists("...tsne_mat")){
            dt = dcast( bw_dt[mark == m], id+cell~x, value.var = "y")
            ...tsne_mat = as.matrix(dt[, -1:-2])
            rn = paste(dt$id, dt$cell)
        }else{
            dt = dcast( bw_dt[mark == m], id+cell~x, value.var = "y")
            stopifnot(all(paste(dt$id, dt$cell) == rn))
            ...tsne_mat = cbind(...tsne_mat, as.matrix(dt[, -1:-2]))
            
            
        }
    }
    rownames(...tsne_mat) = rn
    return(list(bw_dt = bw_dt, tsne_mat = ...tsne_mat, query_gr = qgr))
}

run_tsne = function(tsne_mat, perplexity = 100,
                    n_cores = getOption("mc.cores", 1),
                    high_topright = TRUE,
                    norm1 = TRUE,
                    rname = digest::digest(list(
                        tsne_mat[sample(1:nrow(tsne_mat), 20),
                                 sample(1:ncol(tsne_mat), 20)],
                        perplexity
                    )),
                    force_overwrite = FALSE){
    set.seed(0)
    res_tsne = bfcif(bfc, rname, force_overwrite = force_overwrite,
                     FUN = function(){
                         Rtsne(tsne_mat, num_threads = n_cores, perplexity = perplexity, check_duplicates = FALSE)
                     })
    
    tdt = as.data.table(res_tsne$Y)
    colnames(tdt) = c("tx", "ty")
    tdt$rn = rownames(tsne_mat)
    tdt[, c("id", "cell") := tstrsplit(rn, " ", keep = 1:2)]
    
    if(norm1){
        tdt$tx = norm1(tdt$tx)-.5
        tdt$ty = norm1(tdt$ty)-.5
    }
    
    
    if(high_topright){
        rs = rowSums(tsne_mat)
        tdt$rs = rs[tdt$rn]
        x_cutoff = mean(range(tdt$tx))
        x_flip = sum(tdt[tx > x_cutoff]$rs) < sum(tdt[tx < x_cutoff]$rs)
        if(x_flip){
            tdt[, tx := max(tx) - tx + min(tx)]
        }
        y_cutoff = mean(range(tdt$ty))
        y_flip = sum(tdt[ty > y_cutoff]$rs) < sum(tdt[ty < y_cutoff]$rs)
        if(y_flip){
            tdt[, ty := max(ty) - ty + min(ty)]
        }
        tdt$rs = NULL
    }
    
    tdt
}

make_tsne_img = function(profiles_dt, position_dt, n_points, 
                         xrng = range(position_dt$tx), 
                         yrng = range(position_dt$ty),
                         rname = digest::digest(list(
                             profiles_dt, position_dt, 
                             n_points, apply_norm, 
                             ylim, line_colors, facet_by)), 
                         odir = file.path("tsne_images/", rname), 
                         force_rewrite = FALSE,
                         apply_norm = TRUE,
                         ylim = c(0, 1),
                         facet_by = NULL,
                         # view_rect = list(),
                         ma_size = 2,
                         n_cores = getOption("mc.cores", 1),
                         line_colors = c("H3K4me3" = "forestgreen", 
                                         "H3K27me3" = "firebrick1")){
    # stopifnot(is.list(view_rect))
    # if(is.null(view_rect$xmin))
    position_dt = copy(position_dt)
    position_dt = position_dt[tx >= min(xrng) & tx <= max(xrng) & ty >= min(yrng) & ty <= max(yrng)]
    #use positional info from position_dt to bin points
    position_dt[, bx := mybin(tx, n_points, xrng = xrng)]
    position_dt[, by := mybin(ty, n_points, xrng = yrng)]
    #merge binning info to profiles
    mdt = merge(profiles_dt, position_dt[, .(bx, by, cell, id)], allow.cartesian=TRUE, by = intersect(colnames(profiles_dt), c("cell", "id")))#, by = c("cell", "id"))
    if(is.null(mdt$mark)) mdt$mark = "signal"
    
    if(is.null(facet_by)){
        mdt = mdt[, .(y = mean(y)), .(bx, by, x, mark)]
    }else{
        mdt = mdt[, .(y = mean(y)), .(bx, by, x, mark, get(facet_by))]
        colnames(mdt)[colnames(mdt) == "get"] = facet_by
    }
    #each combination of bx and by is a unique plot_id
    mdt[, plot_id := paste(bx, by, sep = "_")]
    dir.create(odir, recursive = TRUE, showWarnings = FALSE)
    
    if(is.null(facet_by)){
        img_dt = unique(mdt[, .(bx, by, plot_id)])
        img_dt[, png_file := file.path(odir, paste0(plot_id, ".png"))]
    }else{
        img_dt = unique(mdt[, .(bx, by, plot_id, get(facet_by))])
        colnames(img_dt)[4] = facet_by
        img_dt[, png_file := file.path(odir, paste0(get(facet_by), "_", plot_id, ".png"))]
    }
    
    xs = mybin_centers(tsne_res$tx, n_points, xrng = xrng)
    ys = mybin_centers(tsne_res$ty, n_points, xrng = yrng)
    
    # plot(expand.grid(xs, ys), xlim = xrng, ylim = yrng)
    # rect(min(xrng), min(yrng), max(xrng), max(yrng), col = rgb(0,0,1,.1))
    
    img_dt[, tx := xs[bx]]
    img_dt[, ty := ys[by]]
    
    
    if(apply_norm){
        mdt[, ynorm := y / quantile(y, .95), by = .(mark)]
        mdt[ynorm > 1, ynorm := 1]    
    }else{
        mdt[, ynorm := y]
    }
    
    
    if(force_rewrite){
        file.remove(img_dt$png_file[file.exists(img_dt$png_file)])
    }
    
    if(is.null(facet_by)){
        img_dt = merge(img_dt, position_dt[, .N, .(bx, by)])
    }else{
        tmp_dt = position_dt[, .N, .(bx, by, get(facet_by))]
        colnames(tmp_dt)[colnames(tmp_dt) == "get"] = facet_by
        img_dt = merge(img_dt, tmp_dt)
    }
    
    if(any(!file.exists(img_dt$png_file))){
        plot_info = lapply(which(!file.exists(img_dt$png_file)), function(i){
            # figure out how not to copy global env
            # hidden = parallel::mclapply(which(!file.exists(img_dt$png_file)), function(i){}, mc.cores = n_cores)
            fpath = img_dt$png_file[i]
            p_id = img_dt$plot_id[i]
            if(is.null(facet_by)){
                pdt = mdt[plot_id == p_id]    
            }else{
                pdt = mdt[plot_id == p_id & get(facet_by) == img_dt[[facet_by]][i]]
            }
            
            # pdt[, ysm := seqsetvis:::movingAverage(y, n = 8), by = .(mark)]
            pdt[, ysm := seqsetvis:::movingAverage(ynorm, n = ma_size), by = .(mark)]
            pdt = applySpline(pdt, n = 10, by_ = "mark", y_ = "ysm")
            list(pdt, fpath)
        })
        # hidden = lapply(plot_info, function(x){
        # figure out how not to copy global env
        hidden = parallel::mclapply(plot_info, function(x){
            fpath = x[[2]]
            pdt = x[[1]]
            # pdt[, ysm := seqsetvis:::movingAverage(y, n = 8), by = .(mark)]
            p = ggplot(pdt, aes(x = x, y = ysm, ymin = 0, ymax = ysm, color = mark, fill = mark)) + 
                geom_ribbon(alpha = .3) +
                geom_path(size = .6, alpha = 1) +
                scale_color_manual(values = line_colors) +
                scale_fill_manual(values = line_colors) +
                theme_void() + guides(color = "none", fill = 
                                          'none') + 
                coord_cartesian(ylim = ylim, xlim = c(-.5, .5), expand = FALSE)
            ggsave(fpath, p, width = 2, height = 2, units = "cm")
            # p
            NULL
            # })
        }, mc.cores = n_cores)
    }
    
    if(is.factor(position_dt$cell)){
        if(!is.null(img_dt$cell)){
            img_dt$cell = factor(img_dt$cell, levels = levels(position_dt$cell))
        }
        if(!is.null(mdt$cell)){
            mdt$cell = factor(mdt$cell, levels = levels(position_dt$cell))
        }
    }
    
    return(list(images_dt = img_dt, summary_profiles_dt = mdt, tsne_dt = position_dt, n_points = n_points, xrng = xrng, yrng = yrng))
}

prep_tsne_img = function(simg_dt, 
                         n_points,
                         xrng,
                         yrng,
                         N_floor = 0,
                         N_ceiling = NULL,
                         min_size = .3
){
    if(is.null(N_ceiling)){
        N_ceiling = max(simg_dt$N)
    }
    simg_dt[, img_size := N]
    simg_dt[img_size > N_ceiling, img_size := N_ceiling]
    simg_dt[img_size < N_floor, img_size := N_floor]
    
    simg_dt[, img_size := img_size - N_floor]
    simg_dt[, img_size := img_size / N_ceiling]
    simg_dt[, img_size := img_size]
    simg_dt = simg_dt[img_size >= min_size]
    
    xspc = diff(xrng)/n_points/2
    yspc = diff(yrng)/n_points/2
    
    simg_dt[, xmin := tx - xspc * img_size]
    simg_dt[, xmax := tx + xspc * img_size]
    simg_dt[, ymin := ty - yspc * img_size]
    simg_dt[, ymax := ty + yspc * img_size]
    simg_dt
}

plot_tsne_img = function(images_dt,
                         n_points,
                         xrng = c(-.5, .5),
                         yrng = c(-.5, .5),
                         N_floor = 0,
                         N_ceiling = NULL,
                         min_size = .3, 
                         show_plot = TRUE
){
    simg_dt = copy(images_dt) 
    simg_dt = prep_tsne_img(simg_dt, 
                            n_points = n_points,
                            xrng = xrng,
                            yrng = yrng,
                            N_floor = N_floor,
                            N_ceiling = N_ceiling,
                            min_size = min_size
    )
    
    p = ggplot() + 
        geom_image.rect(data = simg_dt, 
                        aes(xmin = xmin, xmax = xmax, 
                            ymin = ymin, ymax = ymax, 
                            image = png_file)) +#, color = rgb(0,0,1,.2)) + 
        geom_rect(data = simg_dt, 
                  aes(xmin = xmin, xmax = xmax, 
                      ymin = ymin, ymax = ymax),
                  fill = NA, color = "black")
    if(show_plot) print(p)
    invisible(list(plot = p, plot_data = simg_dt))
}

plot_tsne_img_byCell = function(images_dt,
                                tsne_dt,
                                n_points,
                                xrng = c(-.5, .5),
                                yrng = c(-.5, .5),
                                N_floor = 0,
                                N_ceiling = NULL,
                                min_size = .3, 
                                show_plot = TRUE
){
    simg_dt = merge(images_dt[, .(bx, by, plot_id, png_file, tx, ty)], 
                    tsne_dt[, .(N = .N), .(cell, bx, by)])
    simg_dt = prep_tsne_img(simg_dt, 
                            n_points = n_points,
                            xrng = xrng,
                            yrng = yrng,
                            N_floor = N_floor,
                            N_ceiling = N_ceiling,
                            min_size = min_size
    )
    p = ggplot(simg_dt, aes(xmin = xmin, xmax = xmax, 
                            ymin = ymin, ymax = ymax, 
                            image = png_file)) + 
        geom_image.rect() +#, color = rgb(0,0,1,.2)) + 
        geom_rect(fill = NA, color = "black") +
        facet_wrap("cell", drop = FALSE)
    if(show_plot) print(p)
    invisible(list(plot = p, plot_data = simg_dt))
}

calc_delta = function(tsne_res, cell_a, cell_b, n_points){
    v_dt = dcast(tsne_res[cell %in% c(cell_a, cell_b)], "id~cell", value.var = c("tx", "ty"))
    colnames(v_dt) = sub(cell_a, "cell_a", colnames(v_dt))
    colnames(v_dt) = sub(cell_b, "cell_b", colnames(v_dt))
    v_dt$bx_cell_a = mybin(v_dt$tx_cell_a, n_points = n_points, xrng = range(tsne_res$tx))
    xs = mybin_centers(v_dt$tx_cell_a, n_points = n_points, xrng = range(tsne_res$tx))
    v_dt$btx_cell_a = xs[v_dt$bx_cell_a]
    
    v_dt$by_cell_a = mybin(v_dt$ty_cell_a, n_points = n_points, xrng = range(tsne_res$ty))
    ys = mybin_centers(v_dt$ty_cell_a, n_points = n_points, xrng = range(tsne_res$ty))
    v_dt$bty_cell_a = ys[v_dt$by_cell_a]
    
    av_dt = v_dt[, .(tx_cell_b = mean(tx_cell_b), ty_cell_b = mean(ty_cell_b), N = .N), .(bx_cell_a, by_cell_a)]
    av_dt$tx_cell_a = xs[av_dt$bx_cell_a]
    av_dt$ty_cell_a = ys[av_dt$by_cell_a]
    return(list(velocity_dt = v_dt, agg_velocity_dt = av_dt))
}


bfcif = function(bfc, rname, FUN, force_overwrite = FALSE){
    # is rname in cache?
    if(nrow(BiocFileCache::bfcquery(bfc, query = rname, field = "rname")) == 0){
        cache_path = BiocFileCache::bfcnew(bfc, rname = rname)
        
    }else{
        cache_path = BiocFileCache::bfcrpath(bfc, rname)
    }
    # does cached file exist?
    if(file.exists(cache_path) && !force_overwrite){
        message("loading cached results.")
        load(BiocFileCache::bfcrpath(bfc, rname))
    }else{
        message("running function and caching results...")
        res = FUN()
        save(res, file = cache_path)
    }
    # return either new results or cached results
    res
}

make_tss_plot = function(qcell, as_facet = TRUE){
    xrng = range(tsne_res$tx)
    yrng = range(tsne_res$ty)
    
    pr_img_res = make_tsne_img(
        bw_dt = pr_dt, apply_norm = FALSE, 
        tdt = tsne_res[cell %in% qcell], #force_rewrite = TRUE, 
        xrng = xrng,
        yrng = yrng,
        n_points = n_points, line_colors = c("signal" = "black")
    )
    if(as_facet){
        p_pr_density = plot_tsne_img_byCell(pr_img_res$images_dt, 
                                            pr_img_res$tsne_dt[cell %in% qcell],
                                            n_points = n_points, N_ceiling = NULL)$plot +
            coord_cartesian(xlim = xrng, ylim = yrng)
        p_h7_density = plot_tsne_img_byCell(img_res$images_dt, 
                                            img_res$tsne_dt[cell %in% qcell], 
                                            n_points = n_points, N_ceiling = NULL)$plot +
            coord_cartesian(xlim = xrng, ylim = yrng)
    }else{
        p_pr_density = plot_tsne_img(pr_img_res$images_dt, 
                                     # pr_img_res$tsne_dt[cell %in% qcell],
                                     n_points = n_points, N_ceiling = NULL)$plot +
            coord_cartesian(xlim = xrng, ylim = yrng)
        p_h7_density = plot_tsne_img(img_res$images_dt, 
                                     # img_res$tsne_dt[cell %in% qcell], 
                                     n_points = n_points, N_ceiling = NULL)$plot +
            coord_cartesian(xlim = xrng, ylim = yrng)
    }
    
    
    pg = cowplot::plot_grid(p_h7_density + labs(title = paste(paste(qcell, collapse = ", "), ": k4me3+k27me3")), 
                            p_pr_density + labs(title = paste(paste(qcell, collapse = ", "), ": tss frequency")))
    # ggsave(paste0("tmp_", qcell, "_tss.pdf"), plot = pg, width = 8, height = 4)
    pg
    # head(pr_img_res$tsne_dt[cell == cell, .N, by = .(bx, by)][order(N, decreasing = TRUE)])
    # head(img_res$tsne_dt[cell == cell, .N, by = .(bx, by)][order(N, decreasing = TRUE)])
}

make_img_plots_facet = function(img_results, qcell = NULL, 
                                xrng = c(-.5, .5), 
                                yrng = c(-.5, .5), 
                                N_floor = 0,
                                N_ceiling = NULL,
                                min_size = .3,
                                facet_by = "cell"){
    return_list = TRUE
    if(all(c("images_dt", "summary_profiles_dt", "tsne_dt") %in% names(img_results))){
        img_results = list(img_results)
        return_list = FALSE
    }
    stopifnot(is.list(img_results))
    # stopifnot(is.list(img_results$img_res))
    
    if(is.null(qcell)) 
        qcell = 
        as.character(unique(
            img_results[[1]]$tsne_dt$cell
        ))
    
    plots = lapply(img_results, function(x){
        img_dt = copy(x$images_dt)
        # img_dt$N = NULL
        # tdt = x$tsne_dt[cell %in% qcell, .(.N), .(bx, by)]
        # img_dt = merge(img_dt, tdt, by = c("bx", "by"))
        plot_tsne_img(img_dt, 
                      n_points = x$n_points, 
                      N_floor = N_floor,
                      N_ceiling = N_ceiling, 
                      min_size = min_size,
                      show_plot = FALSE, 
                      xrng = xrng, 
                      yrng = yrng)$plot +
            coord_cartesian(xlim = xrng, ylim = yrng) +
            facet_wrap(facet_by, drop = FALSE)
    })
    
    
    # pg = cowplot::plot_grid(plotlist = plots, nrow = length(plots))
    # pg
    if(return_list){
        plots    
    }else{
        plots[[1]]
    }
    
}

make_img_plots = function(img_results, qcell = NULL, 
                          xrng = c(-.5, .5), 
                          yrng = c(-.5, .5), 
                          N_floor = 0,
                          N_ceiling = NULL,
                          min_size = .3,
                          as_facet = TRUE){
    return_list = TRUE
    if(all(c("images_dt", "summary_profiles_dt", "tsne_dt") %in% names(img_results))){
        img_results = list(img_results)
        return_list = FALSE
    }
    stopifnot(is.list(img_results))
    
    # img_results = lapply(img_results, function(x){
    #     if(is.null(x$cell)){
    #         x$cell = factor("cell")
    #     }
    #     x
    # })
    if(is.null(qcell)) 
        qcell = 
            levels(img_results[[1]]$tsne_dt$cell)
        
    
    if(as_facet){
        plots = lapply(img_results, function(x){
            pdt = x$tsne_dt[cell %in% qcell]
            pdt$cell = factor(pdt$cell, levels = qcell)
            plot_tsne_img_byCell(x$images_dt, 
                                 pdt,
                                 n_points = x$n_points, 
                                 N_floor = N_floor,
                                 N_ceiling = N_ceiling, 
                                 min_size = min_size,
                                 show_plot = FALSE, 
                                 xrng = xrng, yrng = 
                                     yrng)$plot +
                coord_cartesian(xlim = xrng, ylim = yrng)
        })
    }else{
        plots = lapply(img_results, function(x){
            img_dt = copy(x$images_dt)
            img_dt$N = NULL
            tdt = x$tsne_dt[cell %in% qcell, .(.N), .(bx, by)]
            img_dt = merge(img_dt, tdt)
            plot_tsne_img(img_dt, 
                          n_points = x$n_points, 
                          N_floor = N_floor,
                          N_ceiling = N_ceiling, 
                          min_size = min_size,
                          show_plot = FALSE, 
                          xrng = xrng, 
                          yrng = yrng)$plot +
                coord_cartesian(xlim = xrng, ylim = yrng)
        })
    }
    
    
    # pg = cowplot::plot_grid(plotlist = plots, nrow = length(plots))
    # pg
    if(return_list){
        plots    
    }else{
        plots[[1]]
    }
}

ang_cap = function(angle){
    # angle[angle > 360]
    angle - floor(angle / 360)*360
}

xy2deg = function(x1, y1, x2, y2){
    x = x2 - x1
    y = y2 - y1
    deg = atan(y/x) * 180 / pi + 180
    deg[x < 0] = deg[x < 0] + 180
    # deg[deg > 360] = deg[deg > 360] - 360
    deg = 360 - (deg + 90)
    ang_cap(deg)
}

xy2dist = function(x1, y1, x2, y2){
    x = x2 - x1
    y = y2 - y1
    (x^2 + y^2)^.5
}

plot_velocity_arrows_binned = function(tsne_res, 
                                       cell_a, cell_b,
                                       n_points,
                                       p = NULL,
                                       min_N = 0,
                                       id_to_plot = NULL
){
    if(is.null(id_to_plot)){
        tsne_res.tp = copy(tsne_res)
    }else{
        tsne_res.tp = tsne_res[id %in% id_to_plot]    
    }
    
    av_dt = calc_delta(tsne_res.tp, cell_a, cell_b, n_points)$agg_velocity_dt
    
    p_velocity = ggplot() + 
        # geom_density2d(data = tsne_res, aes(x = tsne_res$tx, y = tsne_res$ty, color = "lightgray") +
        # geom_density2d(data = tsne_res, aes(x = tx, y = ty), color = "lightgray") +
        stat_density_2d(data = tsne_res, aes(x = tx, y = ty, fill = stat(level)), geom = "polygon", bins = 7) +
        # geom_point(data = tsne_res.tp[cell %in% c(cell_a, cell_b)][id %in% sampleCap(id, 500)], aes(x = tx, y = ty, color = cell)) +
        geom_point(data = tsne_res.tp[cell %in% c(cell_a, cell_b)], aes(x = tx, y = ty, color = cell)) +
        geom_segment(data = av_dt[N >= min_N], 
                     aes(x = tx_cell_a, xend = tx_cell_b, 
                         y = ty_cell_a, yend = ty_cell_b, 
                         size = N), arrow = arrow()) + 
        coord_cartesian(xlim = range(tsne_res$tx), ylim = range(tsne_res$ty)) +
        scale_fill_gradient2(low = "gray", high = "black") +
        labs(x = "x", y = "y", fill = "density") +
        scale_size_continuous(range = c(.5, 2), 
                              breaks = range(av_dt$N)) + 
        theme_classic()
    p_velocity
}


plot_velocity_arrows = function(tsne_res, cell_a, cell_b,
                                p = NULL,
                                id_to_plot = NULL,
                                max_plotted = 500,
                                delta.min = 0,
                                delta.max = Inf,
                                angle.min = 0, 
                                angle.max = 360){
    v_dt = calc_delta(tsne_res, cell_a, cell_b, n_points)$velocity_dt
    v_dt[, distance := xy2dist(x1 = tx_cell_a, x2 = tx_cell_b, y1 = ty_cell_a, y2 = ty_cell_b)]
    v_dt = v_dt[distance >= delta.min & distance <= delta.max]
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
    bg = v_dt.tp[foreground == FALSE & distance >= delta.min]
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

plot_velocity_arrows_selected = function(tsne_res, qgr, qcells, tss_ids,
                                         grp_var = c("id", "grp")[1],
                                         line_type = c("curve", "spline", "straight")[2],
                                         label_type = c("text", "label", "none")[2]){
    # qcells = c("H7", "CD34", "Kasumi1", "mm1s", "Nalm6")
    # tss_ids = subset(qgr, gene_name == "RUNX1")$id[1]
    stopifnot(qcells %in% unique(tsne_res$cell))
    if(!tss_ids %in% tsne_res$id){
        tmp = unlist(strsplit(tss_ids, " "))
        if(length(tmp) > 1){
            tss_ids = tmp[1]
            tmp = as.numeric(tmp[-1])
            tss_ids = subset(qgr, gene_name == tss_ids)$id[tmp]
        }else{
            tss_ids = subset(qgr, gene_name == tss_ids)$id
        }
        
    }
    names(qgr) = qgr$id
    message(paste(as.character(qgr[tss_ids]), collapse = "\n"))
    stopifnot(tss_ids %in% tsne_res$id)
    
    lines_dt = tsne_res[cell %in% qcells & id %in% tss_ids]
    
    lines_dt$cell = factor(lines_dt$cell, levels = qcells)
    lines_dt = lines_dt[order(cell)][order(id)][]
    lines_dt[, cell_o := seq(.N), by = .(id)]
    # lines_dt
    lines_dt$id = as.character(qgr[lines_dt$id])
    
    
    # lines_dt$tx = scales::rescale(c(1,1, 2,2,1,1:5), to = c(-.5,.5))
    # lines_dt$ty = scales::rescale(c(0,1, 1.5,2.5,3,1:5), to = c(-.5,.5))
    # 
    # 
    
    # lines_dt[, .(tx = spline(x = pid, y = tx, n = n*(length(qcells)-1)), 
    # ty = spline(x = pid, y = ty, n = n*(length(qcells)-1))), by = id]
    p =     ggplot() + 
        geom_point(data = tsne_res[sample(seq(nrow(tsne_res)), 5000),],
                   aes(x = tx, y = ty), color = "gray") + 
        labs(title = paste(qcells, collapse = ", ")) +
        theme_classic() +
        scale_color_brewer(palette = "Dark2")
    switch(line_type, 
           curve = {
               plot_dt = merge(lines_dt[seq_along(qcells)[-length(qcells)],.(tx, ty, id, cell_o)],
                               lines_dt[seq_along(qcells)[-1], .(tx_end = tx, ty_end = ty, id, cell_o = cell_o -1)])
               switch(grp_var, 
                      grp = {
                          p = p +
                              geom_curve(data = plot_dt,
                                         aes(x = tx, y = ty, xend = tx_end, yend = ty_end, color = id),
                                         size = 1, arrow = arrow())
                      },
                      id = {
                          p = p +
                              geom_curve(data = plot_dt[cell_o < max(cell_o)],
                                         aes(x = tx, y = ty, xend = tx_end, yend = ty_end, color = id),
                                         size = 1) +
                              geom_curve(data = plot_dt[cell_o == max(cell_o)],
                                         aes(x = tx, y = ty, xend = tx_end, yend = ty_end, color = id),
                                         size = 1, arrow = arrow())
                      })
           },
           spline = {
               n = 20
               sp_y = lines_dt[, spline(x = cell_o, y = ty, 
                                        n = n*(length(qcells)-1)), by = id][
                                            , .(pid = seq(.N), ty = y), by = .(id)]
               sp_x = lines_dt[, spline(x = cell_o, y = tx, 
                                        n = n*(length(qcells)-1)), by = id][
                                            , .(pid = seq(.N), tx = y), by = .(id)]
               sp_dt = merge(sp_x, sp_y, by = c("id", "pid"))
               ceiling(sp_dt$pid/n)
               
               sp_dt[, grp := ceiling(pid / n)]
               sp_dt[, grp_o := seq(.N), by = .(grp, id)]
               start_dt = merge(lines_dt[cell_o < length(qcells), .(tx, ty, grp = cell_o, id)], 
                                unique(sp_dt[, .(id, grp)]))[, grp_o := 0]
               end_dt = merge(lines_dt[cell_o > 1 & cell_o < length(qcells), .(tx, ty, grp = cell_o-1, id)], 
                              unique(sp_dt[, .(id, grp = grp)]))[, grp_o := n+1]
               plot_dt = rbind(
                   sp_dt[, .(grp, id, tx, ty, grp_o)],
                   start_dt, 
                   end_dt)[order(grp_o)][order(id)][order(grp)]
               switch(grp_var, 
                      grp = {
                          p = p +
                              geom_path(data = plot_dt,
                                        aes(x = tx, y = ty, color = id, group = paste(grp,id)),
                                        arrow = arrow(),
                                        size = 1.2, alpha = 1,
                                        show.legend = FALSE)
                      },
                      id = {
                          p = p +
                              geom_path(data = plot_dt,
                                        aes(x = tx, y = ty, color = id, group = id),
                                        arrow = arrow(),
                                        size = 1.2, alpha = 1,
                                        show.legend = FALSE)
                      })
               
           }, 
           straight = {
               switch(grp_var, 
                      grp = {
                          plot_dt = merge(lines_dt[seq_along(qcells)[-length(qcells)],.(tx, ty, id, cell_o)],
                                          lines_dt[seq_along(qcells)[-1], .(tx_end = tx, ty_end = ty, id, cell_o = cell_o -1)])
                          p = p +
                              geom_segment(data = plot_dt,
                                           aes(x = tx, y = ty, xend = tx_end, yend = ty_end, color = id),
                                           size = 1, arrow = arrow())
                      },
                      id = {
                          plot_dt = lines_dt
                          p = p + geom_path(data = plot_dt, aes(x = tx, y = ty), arrow = arrow())
                      })
               
           })
    p = p + geom_point(data = lines_dt, 
                       aes(x = tx, y = ty, color = id),
                       size = 3, shape = 21, fill = "white")
    switch(label_type,
           text = {
               p = p + ggrepel::geom_text_repel(data = lines_dt,
                                                aes(x = tx, y = ty, color = id, label = cell),
                                                show.legend = FALSE)
           },
           label = {
               p = p + ggrepel::geom_label_repel(data = lines_dt,
                                                 aes(x = tx, y = ty, color = id, label = cell),
                                                 fill = "white", show.legend = FALSE)
           }, 
           none = {
               p = p
           })
    p
}

plot_profiles_selected = function(data_dt, qgr, qcells, tss_ids){
    # qcells = c("H7", "CD34", "Kasumi1", "mm1s", "Nalm6")
    # tss_ids = subset(qgr, gene_name == "RUNX1")$id[1]
    stopifnot(qcells %in% unique(data_dt$cell))
    if(!tss_ids %in% data_dt$id){
        tmp = unlist(strsplit(tss_ids, " "))
        if(length(tmp) > 1){
            tss_ids = tmp[1]
            tmp = as.numeric(tmp[-1])
            tss_ids = subset(qgr, gene_name == tss_ids)$id[tmp]
        }else{
            tss_ids = subset(qgr, gene_name == tss_ids)$id
        }
        
    }
    names(qgr) = qgr$id
    message(paste(as.character(qgr[tss_ids]), collapse = "\n"))
    plot_dt = data_dt[id %in% tss_ids & cell %in% qcells]
    plot_dt$cell = factor(plot_dt$cell, levels = qcells)
    plot_dt$id = factor(plot_dt$id, levels = tss_ids)
    
    p = ggplot(plot_dt, aes(x = x, ymin = 0, ymax = y, y = y, color = mark, fill = mark)) + 
        facet_grid("cell~id", switch = "y") + 
        geom_ribbon(alpha = .5) +
        geom_path(show.legend = FALSE) + 
        theme_classic() + 
        theme(strip.background = element_blank(), strip.placement = "outside",
              strip.text.y = element_text(angle = 180)) +
        labs(x = "", y = "") +
        scale_color_manual(values = c("H3K27me3" = "firebrick", "H3K4me3" = "forestgreen")) +
        scale_fill_manual(values = c("H3K27me3" = "firebrick", "H3K4me3" = "forestgreen"))
    p
}

library(magick)
library(grid)
ggname = ggimage:::ggname
color_image = ggimage:::color_image
# color_image = function (img, color, alpha = NULL) 
# {
#     if (is.null(color)) 
#         return(img)
#     if (length(color) > 1) {
#         stop("color should be a vector of length 1")
#     }
#     bitmap <- img[[1]]
#     col <- col2rgb(color)
#     bitmap[1, , ] <- as.raw(col[1])
#     bitmap[2, , ] <- as.raw(col[2])
#     bitmap[3, , ] <- as.raw(col[3])
#     if (!is.null(alpha) && alpha != 1) 
#         browser()
#         if(dim(bitmap)[1] == 3){
#             bitmap[,4] = as.raw(255)
#         }
#         bitmap[4, , ] <- as.raw(as.integer(bitmap[4, , ]) * alpha)
#     image_read(bitmap)
# }
##' geom layer for visualizing image files
##'
##'
##' @title geom_image.rect
##' @param mapping aes mapping
##' @param data data
##' @param stat stat
##' @param position position
##' @param inherit.aes logical, whether inherit aes from ggplot()
##' @param na.rm logical, whether remove NA values
##' @param by one of 'width' or 'height'
##' @param nudge_x horizontal adjustment to nudge image
##' @param ... additional parameters
##' @return geom layer
##' @importFrom ggplot2 layer
##' @export
##' @examples
##' library("ggplot2")
##' library("ggimage")
##' set.seed(2017-02-21)
##' d <- data.frame(x = rnorm(10),
##'                 y = rnorm(10),
##'                 image = sample(c("https://www.r-project.org/logo/Rlogo.png",
##'                                 "https://jeroenooms.github.io/images/frink.png"),
##'                               size=10, replace = TRUE)
##'                )
##' ggplot(d, aes(x, y)) + geom_image(aes(image=image))
##' @author guangchuang yu
geom_image.rect <- function(mapping=NULL, data=NULL, stat="identity",
                            position="identity", inherit.aes=TRUE,
                            na.rm=FALSE, 
                            # by="width", 
                            # nudge_x = 0, 
                            ...) {
    
    # by <- match.arg(by, c("width", "height"))
    
    layer(
        data=data,
        mapping=mapping,
        geom=GeomImage.rect,
        stat=stat,
        position=position,
        show.legend=NA,
        inherit.aes=inherit.aes,
        params = list(
            na.rm = na.rm,
            # by = by,
            # nudge_x = nudge_x,
            ##angle = angle,
            ...),
        check.aes = FALSE
    )
}


##' @importFrom ggplot2 ggproto
##' @importFrom ggplot2 Geom
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 draw_key_blank
##' @importFrom grid gTree
##' @importFrom grid gList
GeomImage.rect <- ggproto("GeomImage.rect", Geom,
                          setup_data = function(data, params) {
                              if (is.null(data$subset))
                                  return(data)
                              data[which(data$subset),]
                          },
                          
                          default_aes = aes(image=system.file("extdata/Rlogo.png", package="ggimage"), 
                                            #size=0.05, 
                                            colour = NULL, #angle = 0, 
                                            alpha=1),
                          
                          draw_panel = function(data, panel_params, coord, by, na.rm=FALSE,
                                                .fun = NULL, height, image_fun = NULL,
                                                # hjust=0.5, 
                                                # nudge_x = 0, nudge_y = 0, 
                                                asp=1) {
                              # data$x <- data$x + nudge_x
                              # data$y <- data$y + nudge_y
                              data <- coord$transform(data, panel_params)
                              
                              if (!is.null(.fun) && is.function(.fun))
                                  data$image <- .fun(data$image)
                              
                              groups <- split(data, factor(data$image))
                              imgs <- names(groups)
                              grobs <- lapply(seq_along(groups), function(i) {
                                  d <- groups[[i]]
                                  imageGrob.rect(d$xmin, d$xmax, d$ymin, d$ymax, imgs[i], #by, 
                                                 # hjust,
                                                 d$colour, d$alpha, image_fun, #d$angle, 
                                                 asp)
                              })
                              grobs <- do.call("c", grobs)
                              class(grobs) <- "gList"
                              
                              ggname("geom_image.rect",
                                     gTree(children = grobs))
                          },
                          non_missing_aes = c(#"size", 
                              "image"),
                          required_aes = c("xmin", "xmax", "ymin", "ymax"),
                          draw_key = draw_key_image ## draw_key_blank ## need to write the `draw_key_image` function.
)



##' @importFrom magick image_read
##' @importFrom magick image_read_svg
##' @importFrom magick image_read_pdf
##' @importFrom magick image_transparent
##' @importFrom magick image_rotate
##' @importFrom grid rasterGrob
##' @importFrom grid viewport
##' @importFrom grDevices rgb
##' @importFrom grDevices col2rgb
##' @importFrom methods is
##' @importFrom tools file_ext
imageGrob.rect <- function(xmin, xmax, ymin, ymax, img, #by, hjust, 
                           colour, alpha, image_fun, #angle, 
                           asp=1) {
    if (!is(img, "magick-image")) {
        if (tools::file_ext(img) == "svg") {
            img <- image_read_svg(img)
        } else if (tools::file_ext(img) == "pdf") {
            img <- image_read_pdf(img)
        } else {
            img <- image_read(img)
        }
        asp <- getAR2(img)/asp
    }
    
    unit <- "native"
    width = xmax - xmin
    height = ymax - ymin
    # if (any(size == Inf)) {
    #     x <- 0.5
    #     y <- 0.5
    #     width <- 1
    #     height <- 1
    #     unit <- "npc"
    # } else if (by == "width") {
    #     width <- size
    #     height <- size/asp
    # } else {
    #     width <- size * asp
    #     height <- size
    # }
    # 
    # if (hjust == 0 || hjust == "left") {
    #     x <- x + width/2
    # } else if (hjust == 1 || hjust == "right") {
    #     x <- x - width/2
    # }
    
    if (!is.null(image_fun)) {
        img <- image_fun(img)
    }
    
    if (is.null(colour)) {
        grobs <- list()
        grobs[[1]] <- rasterGrob(x = xmin,
                                 y = ymin,
                                 just = c(0,0),
                                 image = img,
                                 default.units = unit,
                                 height = height,
                                 width = width,
                                 interpolate = FALSE)
    } else {
        cimg <- lapply(seq_along(colour), function(i) {
            color_image(img, colour[i], alpha[i])
        })
        
        grobs <- lapply(seq_along(xmin), function(i) {
            img <- cimg[[i]]
            # if (angle[i] != 0) {
            #     img <- image_rotate(img, angle[i])
            #     img <- image_transparent(img, "white")
            # }
            rasterGrob(x = xmin[i],
                       y = ymin[i],
                       just = c(0,0),
                       image = img,
                       default.units = unit,
                       height = height[i],
                       width = width[i],
                       interpolate = FALSE
                       ## gp = gpar(rot = angle[i])
                       ## vp = viewport(angle=angle[i])
            )
        })
    }
    return(grobs)
}


##' @importFrom magick image_info
getAR2 <- function(magick_image) {
    info <- image_info(magick_image)
    info$width/info$height
}


compute_just <- getFromNamespace("compute_just", "ggplot2")


