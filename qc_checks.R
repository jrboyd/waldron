library(ChIPQC)
library(BiocFileCache)
source("setup_files.R")
bfc = BiocFileCache()
bfcif = peakrefine::bfcif
samples = data.table(SampleID = names(c(k4_bams.reps, k27_bams.reps)),
                     bamReads = c(k4_bams.reps, k27_bams.reps), 
                     Peaks = c(k4_peaks.reps, k27_peaks.reps.broad))
samples[, c("Tissue", "Factor", "Replicate") := tstrsplit(SampleID, "_")]
samples$Replicate = as.numeric(sub("R", "", samples$Replicate))
samples[, Tissue := gsub("-", "", Tissue)]
samples[, SampleID := gsub("-", "", SampleID)]

cache_qc = function(samples, bfc){
    FUN = function()
        ChIPQC(as.data.frame(samples), annotaiton="hg38", chromosomes = NULL)
    bfcif(bfc, digest::digest(samples), FUN)
}

k27_exp = cache_qc(samples[Factor == "H3K27me3"], bfc)
k4_exp = cache_qc(samples[Factor == "H3K4me3"], bfc)

ChIPQCreport(k27_exp, reportFolder = "ChIPQC_k27")
ChIPQCreport(k4_exp, reportFolder = "ChIPQC_k4")

plotCoverageHist(k4_exp, facetBy=c("Tissue"))

plotCC(k27_exp, facetBy=c("Tissue"))
plotCC(k4_exp, facetBy=c("Tissue"))
