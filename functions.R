ssvFetchPeak = function(grs, qgr){
    qgr = resize(qgr, 50000, fix = "center")
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

symbol2uniprot = function(x){
    bres = bitr(x, fromType = "SYMBOL",
                toType = c("UNIPROT"),
                OrgDb = org.Hs.eg.db)
    bres[!duplicated(bres$SYMBOL),]$UNIPROT
}

my_clusterProfiler_kegg = function(qgr, clust_assign, bg_genes = NULL){
    peak_dt = as.data.table(qgr)
    peak_dt$name = names(qgr)
    mdt = merge(peak_dt, clust_assign, by.x = "name", by.y = "id")
    anno_dt = my_annotate(GRanges(mdt), ref_gr)
    
    anno_dt = anno_dt[, .(list(unique(gene_name))), by = .(cluster_id)]
    anno_lists = anno_dt$V1
    names(anno_lists) = paste0("c", anno_dt$cluster_id)
    anno_lists = anno_lists[order(anno_dt$cluster_id)]
    names(anno_lists)
    lengths(anno_lists)
    # each cluster vs total genes in heatmap
    if(is.null(bg_genes)){
        bg_genes = unique(unlist(anno_lists))    
    }
    
    # # each cluster vs all marked by k4me3
    # bg_dt = my_annotate(peaks$H7_H3K4ME3)
    # # bg_dt = my_annotate(k4_gr$`CD34-01517_H3K4me3`)
    # bg_genes = unique(bg_dt$gene_name)
    
    gene_lists = lapply(anno_lists, symbol2uniprot)
    
    rname_go_dat = digest::digest(list(gene_lists, bg_genes, "kegg"))
    rname_go_plot = digest::digest(list(gene_lists, bg_genes, "kegg", "plot"))
    
    res = bfcif(bfc, rname_go_plot, 
                force_overwrite = TRUE, 
                function(){
                    message("calc KEGG res")
                    tryCatch(
                        expr = {
                            ck = bfcif(bfc, rname_go_dat, 
                                       force_overwrite = TRUE, function(){
                                           message("calc compareCluster")
                                           compareCluster(geneCluster = gene_lists,
                                                          universe      = symbol2uniprot(bg_genes),
                                                          fun = "enrichKEGG",
                                                          # OrgDb = org.Hs.eg.db,
                                                          keyType       = 'uniprot',
                                                          # ont           = "BP",
                                                          pAdjustMethod = "BH",
                                                          pvalueCutoff  = 0.05,
                                                          qvalueCutoff  = 0.1)
                                       })
                            # ck = my_clusterProfiler(qgr, clust_assign)
                            p = ck %>% dotplot    
                            list(ck, p)
                        }, error = {
                            function(e){
                                ck = NULL
                                p = ggplot() + annotate("text", x = 0, y = 0, label ="no KEGG results")
                                list(ck, p)
                            }
                        })
                })
    res
}

my_clusterProfiler = function(qgr, clust_assign, bg_genes = NULL){
    peak_dt = as.data.table(qgr)
    peak_dt$name = names(qgr)
    mdt = merge(peak_dt, clust_assign, by.x = "name", by.y = "id")
    anno_dt = my_annotate(GRanges(mdt), ref_gr)
    
    anno_dt = anno_dt[, .(list(unique(gene_name))), by = .(cluster_id)]
    anno_lists = anno_dt$V1
    names(anno_lists) = paste0("c", anno_dt$cluster_id)
    anno_lists = anno_lists[order(anno_dt$cluster_id)]
    names(anno_lists)
    lengths(anno_lists)
    # each cluster vs total genes in heatmap
    if(is.null(bg_genes)){
        bg_genes = unique(unlist(anno_lists))    
    }
    
    # # each cluster vs all marked by k4me3
    # bg_dt = my_annotate(peaks$H7_H3K4ME3)
    # # bg_dt = my_annotate(k4_gr$`CD34-01517_H3K4me3`)
    # bg_genes = unique(bg_dt$gene_name)
    
    gene_lists = anno_lists
    
    rname_go_dat = digest::digest(list(gene_lists, bg_genes, "BP"))
    rname_go_plot = digest::digest(list(gene_lists, bg_genes, "BP", "plot"))
    
    res = bfcif(bfc, rname_go_plot, function(){
        message("calc go res")
        tryCatch(
            expr = {
                ck = bfcif(bfc, rname_go_dat, function(){
                    message("calc compareCluster")
                    compareCluster(geneCluster = gene_lists, 
                                   universe      = bg_genes,
                                   fun = "enrichGO", 
                                   OrgDb = org.Hs.eg.db,                  
                                   keyType       = 'SYMBOL',
                                   ont           = "BP",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = 0.01,
                                   qvalueCutoff  = 0.05)
                })
                # ck = my_clusterProfiler(qgr, clust_assign)
                p = ck %>% simplify %>% dotplot    
                list(ck, p)
            }, error = {
                function(e){
                    ck = NULL
                    p = ggplot() + annotate("text", x = 0, y = 0, label ="no GO results")
                    list(ck, p)
                }
            })
    })
    res
}