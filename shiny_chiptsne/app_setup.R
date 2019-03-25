if(!exists("LOADED")){
    LOADED = TRUE
    library(shiny)
    library(data.table)
    library(seqsetvis)
    options(mc.cores = 16)
    
    source("functions.R")
    source("functions_tsne.R")
    load("dataset1.save")
    
    n_points = 16
    piles_img_res = make_tsne_img(profiles_dt = tsne_input$bw_dt, 
                                  position_dt = tsne_res, #force_rewrite = TRUE, 
                                  apply_norm = FALSE, 
                                  ylim = c(0,10),  
                                  # xrng = zoom_x,
                                  # yrng = zoom_y,
                                  n_points = n_points, 
                                  line_colors = c(
                                      "H3K4me3" = "forestgreen", 
                                      "H3K27me3" = "firebrick")
    )
    
    p_basic = make_img_plots(img_results = list(piles_img_res), 
                             qcell = NULL, 
                             min_size = 0, 
                             N_ceiling = NULL,
                             as_facet = FALSE)[[1]]
    
    
    message(getwd())
    
    tsne_res
    tsne_input$tsne_mat = NULL
    
    UI_CELLS = unique(tsne_input$bw_dt$cell)
    UI_CELLS = UI_CELLS[order(UI_CELLS != "Kasumi1")]
    UI_CELLS = UI_CELLS[order(UI_CELLS != "mm1s")]
    UI_CELLS = UI_CELLS[order(UI_CELLS != "CD34")]
    UI_CELLS = UI_CELLS[order(UI_CELLS != "H7")]
    UI_MARKS = unique(tsne_input$bw_dt$mark)
    UI_MARKS = UI_MARKS[order(UI_MARKS != "H3K4me3")]
    n_tp = 5000
    set.seed(1)
    UI_TP = sample(unique(tsne_res$id), n_tp / length(UI_CELLS))
    tsne_tp = tsne_res[id %in% UI_TP]
    
    GLOBAL_VIEW_POINTS = "points"
    GLOBAL_VIEW_PROFILES_FAST = "profiles (fast)"
    GLOBAL_VIEW_PROFILES_SLOW = "profiles (slow)"
    GLOBAL_VIEW_DENSITY = "density"
    
    tsne_res$bx = mybin(tsne_res$tx, n_points = 16)
    tsne_res$by = mybin(tsne_res$ty, n_points = 16)
    
    szt = tsne_res[, .N, .(bx, by)]
    szt[, frac := N / max(N)]
    szt = szt[frac > .2]
    mdt = merge(tsne_input$bw_dt, tsne_res[, .(id, cell, bx, by)], by = c("id", "cell"))
    mdt = mdt[, .(y = mean(y)), by = .(bx, by, x, mark)]
    mdt = merge(mdt, szt)
    mdt$gx = mybin_centers(tsne_res$tx, n_points = 16)[mdt$bx]
    mdt$gy = mybin_centers(tsne_res$ty, n_points = 16)[mdt$by]
    
    # mdt = piles_img_res$summary_profiles_dt
    
    glyph_df = GGally::glyphs(mdt, x_major = "gx", x_minor = "x", y_major = "gy", y_minor = "y")
    ggplot(glyph_df, aes(gx, gy, group = paste(gid, mark), color = mark)) + 
        geom_path() +
        scale_color_manual(values =  c("input" = "blue",
                                       "H3K4me3" = "forestgreen",
                                       "H3K27me3" = "red")) 
    
    qgr = tsne_input$query_gr
    names(qgr) = qgr$id
}