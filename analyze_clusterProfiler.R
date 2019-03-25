setwd("~/R/RNAseq_JonRamsey_melanoma/")
source("~/R/waldron/functions.R")
source("analyze_DEseq2_heatmap.R")

cm_ensg = hres@cluster_members
lengths(cm_ensg)
cm_sym = lapply(cm_ensg, function(x){
    unique(ref[x]$gene_name)
})

library(BiocFileCache)
bfc = BiocFileCache()

cp_kegg = my_clusterProfiler_kegg_fromGenes(cm_sym)
cp_go = my_clusterProfiler_fromGenes(cm_sym)
