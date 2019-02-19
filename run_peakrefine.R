library(peakrefine)
library(BiocFileCache)
library(seqsetvis)
library(GenomicRanges)
library(data.table)
setwd("~/R/waldron/")
source("setup_files.R")

options(mc.cores = 32)

# k4_bams
# k27_bams
# k4_peaks
# k27_peaks

stopifnot(names(k4_bams) == names(k4_peaks))
stopifnot(names(k27_bams) == names(k27_peaks))

bams = c(k4_bams, k27_bams)
peaks = append(easyLoad_narrowPeak(k4_peaks), easyLoad_narrowPeak(k27_peaks))

stopifnot(names(bams) == names(peaks))


f3 = "scc_dt_0-300.save"
if(file.exists(f3)){
    load(f3)
}else{
    scc_dt = lapply(seq_along(bams), function(i){
        message(names(bams)[i])
        peakrefine::calcSCCMetrics(bams[i], qgr = peaks[[i]], 
                                   frag_sizes = 0:60*5, fetch_size = 1800)
    })
    save(scc_dt, file = f3)
}
scc_dt1 = scc_dt
# 
# f1 = "scc_dt_100-300.save"
# if(file.exists(f1)){
#     load(f1)
# }else{
#     scc_dt = lapply(seq_along(bams), function(i){
#         message(names(bams)[i])
#         peakrefine::calcSCCMetrics(bams[i], qgr = peaks[[i]], frag_sizes = 100:300)
#     })
#     save(scc_dt, file = f1)
# }
# scc_dt1 = scc_dt
# 
# 
f2 = "scc_dt_301-600.save"
if(file.exists(f2)){
    load(f2)
}else{
    scc_dt = lapply(seq_along(bams), function(i){
        message(names(bams)[i])
        peakrefine::calcSCCMetrics(bams[i], qgr = peaks[[i]], frag_sizes = 300+1:60*5)
    })

    save(scc_dt, file = f2)
}

scc_dt = append(scc_dt1, scc_dt)
remove(scc_dt1)
# 
scc_dt = lapply(scc_dt, function(x)x$full_correlation_results) %>% rbindlist
# 
# scc_dt = scc_dt[shift %% 5 == 0 | shift < 100]
# 
scc_dt = scc_dt[order(shift)][order(id)]
save(scc_dt, file = "scc_dt_combined2.save")

# scc_dt = rbindlist(lapply(scc_dt, function(x)x$full_correlation_results))
# # scc_dt[order(shift)][order(id)]
# tmp = scc_dt[shift == 150, .N, by = .(cell, mark)]
# dcast(tmp, "cell~mark")
# 
# pdf("scc_k4.pdf")
# for(cl in unique(scc_dt$cell)){
#     for(m in "H3K4me3"){#unique(scc_dt$mark)){
#         message(cl, " ", m)
#         dt = scc_dt[cell == cl & mark == m]
#         p = ggplot(dt[id %in% sample(unique(id), 9)], 
#                    aes(x = shift, y = correlation, group = id)) + 
#             geom_path() + 
#             facet_wrap("id") + 
#             labs(paste(cl, m))        
#         print(p)
#     }    
# }
# dev.off()
# 
# pdf("scc_k27.pdf")
# for(cl in unique(scc_dt$cell)){
#     for(m in "H3K27me3"){#unique(scc_dt$mark)){
#         message(cl, " ", m)
#         dt = scc_dt[cell == cl & mark == m]
#         p = ggplot(dt[id %in% sample(unique(id), min(9, lenght(unique(id))))], 
#                    aes(x = shift, y = correlation, group = id)) + 
#             geom_path() + 
#             facet_wrap("id") + 
#             labs(paste(cl, m))        
#         print(p)
#     }    
# }
# dev.off()










