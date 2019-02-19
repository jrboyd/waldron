library(BiocFileCache)
library(seqsetvis)
library(GenomicRanges)
setwd("~/R/waldron/")
source("setup_files.R")

input_bams = bams[grepl("input", bams)]
input_bam_hESC = "/slipstream/galaxy/uploads/working/qc_framework/output_bivalency_redo_patients_H7/H7_input_pooled/H7_input_pooled.bam"
k4_bam_hESC = "/slipstream/galaxy/uploads/working/qc_framework/output_bivalency_redo_patients_H7/H7_H3K4ME3_pooled/H7_H3K4ME3_pooled.bam"
k27_bam_hESC = "/slipstream/galaxy/uploads/working/qc_framework/output_bivalency_redo_patients_H7/H7_H3K27ME3_pooled/H7_H3K27ME3_pooled.bam"
names(input_bam_hESC) = "H7_input"
names(k4_bam_hESC) = "H7_H3K4me3"
names(k27_bam_hESC) = "H7_H3K27me3"


qgr = GRanges("chr9:133984242-134000152")
qgr = GRanges("chr9:133990000-133995000")

qgr = GRanges("chr9:134064533-134072226")
qgr = resize(GRanges("chr9:133271076-133279300"), 5000, fix = "center")

qgr = GRanges("chr6:45416255-45426709")
qgr = GRanges("chr6:45420882-45424735")

qgr = GRanges("chr6:45230109-45615508")
qgr = GRanges("chr6:44259629-44317234")
qgr = GRanges("chr6:44262165-44277280")

bams = c(k4_bam_hESC, k4_bams, 
         k27_bam_hESC, k27_bams, 
         input_bam_hESC, input_bams)

bam_dt = ssvFetchBam(bams, qgr, 
                     fragLens = NA, 
                     win_size = 10,
                     target_strand = "+", 
                     return_data.table = TRUE, 
                     max_dupes = 1)
bam_dt[["y_minus"]] =
    ssvFetchBam(bams, qgr, 
                fragLens = NA, 
                win_size = 10,
                target_strand = "-", 
                return_data.table = TRUE, 
                max_dupes = 1)[["y"]]

bam_dt[, diff := y - y_minus]

bam_dt[, c("cell", "mark") := tstrsplit(sample, "_")]
# bam_dt[, z := (y - mean(y)) / sd(y), by = .(cell, mark)]
bam_dt[, z := (diff - mean(diff)) / sd(diff), by = .(cell, mark)]

ggplot(bam_dt, aes(x = (start + end) / 2, y = z, color = strand)) + 
    facet_grid("cell~mark") +#, scales = "free_y") + 
    geom_path() + 
    scale_y_continuous(limits = c(-10, 10)) + 
    theme_classic()
ggsave("tmp.pdf", width = 40, height = 15)
