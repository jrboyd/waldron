setClass(Class = "ssvTsne", 
         slots = c(
             bw_dt = "data.table",
             qgr = "GRanges",
             view_size = "integer",
             bfc = "BiocFileCache",
             rname = "character"
             
         ),
         
         validity = function(object){
             errors <- character()
             # mat_cnames = c("i", "j", "val")
             # if (length(intersect(colnames(object@hic_2d), mat_cnames)) != length(mat_cnames)){
             #     msg <- "colnames of hic_2d must be c(i, j, val)"
             #     errors <- c(errors, msg)
             # }
             # reg_cnames = c("seqnames", "start", "end", "index")
             # if (length(intersect(colnames(object@hic_1d), reg_cnames)) != length(reg_cnames)){
             #     msg <- "colnames of hic_1d must be c(seqnames, start, end, index)"
             #     errors <- c(errors, msg)
             # }
             if (length(errors) == 0) TRUE else errors
         }
)

setMethod("initialize", "ssvTsne", function(.Object, matrix_file, regions_file, parameters) {
    # print(matrix_file)
    if(missing(matrix_file) & missing(regions_file) & missing(parameters)){
        return(.Object)
    }
    .Object@matrix_file = matrix_file
    .Object@regions_file = regions_file
    .Object@parameters = parameters
    
    .Object
}

ssvTsne = function(matrix_file, regions_file = NULL, parameters = NULL){
    if(is.null(regions_file)){
        #fix ice first
        mfile = matrix_file
        mfile = sub("_iced.matrix$", ".matrix", mfile)
        mfile = sub("/iced/", "/raw/", mfile)
        regions_file = sub("\\.matrix$", "_abs.bed", mfile)
        if(!file.exists(regions_file)){
            stop("could not auto determine valid regions_file, please supply.")
        }
    }
    if(is.null(parameters)){
        tst = read.table(regions_file, nrows = 2)
        bin_size = tst[2, 2] - tst[1, 2]
        parameters = HiC_parameters(bin_size = bin_size)
    }
    new("ssvTsne",
        matrix_file = matrix_file,
        regions_file = regions_file,
        parameters = parameters)
}

setMethod("show", "ssvTsne",
          function(object) {
              nr_mat = nrow(object@hic_2d)
              nr_reg = nrow(object@hic_1d)
              covered = nr_mat / ((nr_reg^2 - nr_reg) / 2)
              print(paste("size is", format(object.size(object), units = "GB")))
              print(paste0(round(covered*100, 2), "% of bins have signal"))
          }
)