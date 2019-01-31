library(magrittr)
library(GenomicRanges)
library(data.table)
library(seqsetvis)

### pooled 
wd = "/slipstream/galaxy/uploads/working/qc_framework/output_waldron_bivalency/"
bams = dir(wd, pattern = "ed$", full.names = TRUE) %>% dir(., pattern = ".bam$", full.names = TRUE)
names(bams) = basename(bams) %>% sub("_pool.+", "", .) %>% sub("ME3", "me3", .)
names(bams)

k4_bams = bams[grepl("K4", names(bams))]
k27_bams = bams[grepl("K27", names(bams))]

k4_bws = dir(wd, pattern = "K4.+ed$", full.names = TRUE) %>% dir(., pattern = "FE.bw$", full.names = TRUE)
names(k4_bws) = basename(k4_bws) %>% sub("_pool.+", "", .) %>% sub("ME3", "me3", .)
k27_bws = dir(wd, pattern = "K27.+ed$", full.names = TRUE) %>% dir(., pattern = "FE.bw$", full.names = TRUE)
names(k27_bws) = basename(k27_bws) %>% sub("_pool.+", "", .) %>% sub("ME3", "me3", .)
bws = c(k4_bws, k27_bws)

names(bws)

k4_peaks = dir(wd, pattern = "ed$", full.names = TRUE) %>% dir(., pattern = "H3K4.+peaks.narrowPeak$", full.names = TRUE)
names(k4_peaks) = basename(k4_peaks) %>% sub("_pooled.+", "", .) %>% sub("ME3", "me3", .)
k27_peaks.narrow = dir(wd, pattern = "ed$", full.names = TRUE) %>% dir(., pattern = "H3K27.+peaks.narrowPeak$", full.names = TRUE)
names(k27_peaks.narrow) = basename(k27_peaks.narrow) %>% sub("_pooled.+", "", .) %>% sub("ME3", "me3", .)

k27_peaks.broad = dir(wd, pattern = "ed$", full.names = TRUE) %>% dir(., pattern = "H3K27.+peaks.broadPeak$", full.names = TRUE)
names(k27_peaks.broad) = basename(k27_peaks.broad) %>% sub("_pooled.+", "", .) %>% sub("ME3", "me3", .)

k27_peaks = k27_peaks.narrow

stopifnot(names(k4_bams) == names(k4_peaks))
stopifnot(names(k27_bams) == names(k27_peaks))

### reps
bams.reps = dir(wd, pattern = "R[0-9]$", full.names = TRUE) %>% dir(., pattern = ".bam$", full.names = TRUE)
names(bams.reps) = basename(bams.reps) %>% sub(".bam", "", .) %>% sub("ME3", "me3", .)
names(bams.reps)

k4_bams.reps = bams.reps[grepl("K4", names(bams.reps))]
k27_bams.reps = bams.reps[grepl("K27", names(bams.reps))]

bws.reps = dir(wd, pattern = "R[0-9]$", full.names = TRUE) %>% dir(., pattern = "FE.bw$", full.names = TRUE)
names(bws.reps) = basename(bws.reps) %>% sub("_FE.bw", "", .) %>% sub("ME3", "me3", .)
names(bws.reps)

k4_peaks.reps = dir(wd, pattern = "R[0-9]$", full.names = TRUE) %>% dir(., pattern = "H3K4.+peaks.narrowPeak$", full.names = TRUE)
k4_peaks.reps = k4_peaks.reps[!grepl("loose", k4_peaks.reps)]
names(k4_peaks.reps) = basename(k4_peaks.reps) %>% sub("_peaks.+", "", .) %>% sub("ME3", "me3", .)

k27_peaks.reps.narrow = dir(wd, pattern = "R[0-9]$", full.names = TRUE) %>% dir(., pattern = "H3K27.+peaks.narrowPeak$", full.names = TRUE)
k27_peaks.reps.narrow = k27_peaks.reps.narrow[!grepl("loose", k27_peaks.reps.narrow)]
names(k27_peaks.reps.narrow) = basename(k27_peaks.reps.narrow) %>% sub("_peaks.narrowPeak", "", .) %>% sub("ME3", "me3", .)

k27_peaks.reps.broad = dir(wd, pattern = "R[0-9]$", full.names = TRUE) %>% dir(., pattern = "H3K27.+peaks.broadPeak$", full.names = TRUE)
k27_peaks.reps.broad = k27_peaks.reps.broad[!grepl("loose", k27_peaks.reps.broad)]
names(k27_peaks.reps.broad) = basename(k27_peaks.reps.broad) %>% sub("_peaks.broadPeak", "", .) %>% sub("ME3", "me3", .)

stopifnot(names(k4_bams.reps) == names(k4_peaks.reps))
stopifnot(names(k27_bams.reps) == names(k27_peaks.reps.narrow))
stopifnot(names(k27_bams.reps) == names(k27_peaks.reps.broad))
