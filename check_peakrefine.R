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

load("scc_dt_combined2.save")
eshift_dt = scc_dt[, .(shift = shift[correlation == max(correlation)], 
                       correlation = max(correlation)), 
                   .(cell, mark, pid)]

eshift_dt[, quantile(correlation), .(cell, mark)]

# scc_dt[, c("cell", "mark", "pid") := tstrsplit(id, "_", keep = c(1:2, 5))]
# scc_dt = scc_dt[order(shift)]

# scc_dt[order(shift)][order(id)]
# tmp = scc_dt[shift == 150, .N, by = .(cell, mark)]
# dcast(tmp, "cell~mark")

bad_id = c("DOHH2_H3K27me3_pooled_peak_283",
           "DOHH2_H3K27me3_pooled_peak_343",
           "DOHH2_H3K27me3_pooled_peak_606")

good_id = c("DOHH2_H3K27me3_pooled_peak_2104",
            "DOHH2_H3K27me3_pooled_peak_2340",
            "DOHH2_H3K27me3_pooled_peak_2558")
q_dt = scc_dt[id %in% c(bad_id, good_id)]
q_dt[shift < 100]
ggplot(q_dt, aes(x = shift, y = correlation)) + facet_wrap("id") + geom_path()

gr = peaks$DOHH2_H3K27me3
gr = subset(gr, name %in% c(bad_id, good_id))       
gr$group = "bad"
names(gr) = gr$name
gr[good_id]$group = "good"

bam_file = bams["DOHH2_H3K27me3"]
bam_dt = ssvFetchBam(bam_file, resize(gr, 40000, fix = "center"), 
                     target_strand = "both", fragLens = NA, 
                     win_size = 50, return_data.table = TRUE)
bam_dt[, id := tstrsplit(id, "_", keep = 5)]
ggplot(bam_dt, aes(x = x, y = y, color = strand)) + geom_path() + facet_wrap("id", scales = "free_y")

diff_dt = dcast(bam_dt, "id+x+sample~strand", value.var = "y")
diff_dt[, ydiff := `+` - `-`]

ggplot(diff_dt[abs(x) < 5000], aes(x = x, ymax = ydiff, ymin = 0)) + geom_ribbon() + facet_wrap("id")

#shifted

bam_dt = ssvFetchBam(bam_file, GenomicRanges::shift(resize(gr, 40000, fix = "center"), 1e5), 
                     target_strand = "both", fragLens = NA, win_size = 50, return_data.table = TRUE)
bam_dt[, id := tstrsplit(id, "_", keep = 5)]
ggplot(bam_dt, aes(x = x, y = y, color = strand)) + geom_path() + facet_wrap("id", scales = "free_y")

diff_dt2 = dcast(bam_dt, "id+x+sample~strand", value.var = "y")
diff_dt2[, ydiff := `+` - `-`]

ggplot(diff_dt2, aes(x = x, ymax = ydiff, ymin = 0)) + geom_ribbon() + facet_wrap("id")

diff_dt$group = "peak"
diff_dt2$group = "ctrl"
diff_dt_c = rbind(diff_dt, diff_dt2)

diff_dt_c[, ysum := `+` + `-`]

ggplot(diff_dt_c, aes(x = x, ymax = ydiff, ymin = 0, fill = group)) + geom_ribbon() + facet_wrap("id~group")
ggplot(diff_dt_c, aes(x = x, ymax = ysum, ymin = 0, fill = group)) + geom_ribbon() + facet_wrap("id~group")

#extended
bam_dt = ssvFetchBam(bam_file, resize(gr, 40000, fix = "center"), 
                     target_strand = "both", fragLens = 250, 
                     win_size = 50, return_data.table = TRUE)
bam_dt[, id := tstrsplit(id, "_", keep = 5)]
ggplot(bam_dt, aes(x = x, y = y, color = strand)) + geom_path() + facet_wrap("id", scales = "free_y")

diff_dt = dcast(bam_dt, "id+x+sample~strand", value.var = "y")
diff_dt[, ydiff := `+` - `-`]

bam_dt = ssvFetchBam(bam_file, GenomicRanges::shift(resize(gr, 40000, fix = "center"), 1e5), 
                     target_strand = "both", fragLens = 250, win_size = 50, return_data.table = TRUE)
bam_dt[, id := tstrsplit(id, "_", keep = 5)]
ggplot(bam_dt, aes(x = x, y = y, color = strand)) + geom_path() + facet_wrap("id", scales = "free_y")

diff_dt2 = dcast(bam_dt, "id+x+sample~strand", value.var = "y")
diff_dt2[, ydiff := `+` - `-`]

diff_dt$group = "peak"
diff_dt2$group = "ctrl"
diff_dt_c = rbind(diff_dt, diff_dt2)
diff_dt_c[, ymin := min(`+`, `-`), .(id, x, sample, ydiff, group)]

ggplot(diff_dt_c, aes(x = x, ymax = ydiff, ymin = 0, fill = group)) + geom_ribbon() + facet_wrap("id~group")
ggplot(diff_dt_c, aes(x = x, ymax = ysum, ymin = 0, fill = group)) + geom_ribbon() + facet_wrap("id~group")
ggplot(diff_dt_c, aes(x = x, ymax = ymin, ymin = 0, fill = group)) + geom_ribbon() + facet_wrap("id~group")
