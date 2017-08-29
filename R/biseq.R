

#' Read in files from biseq meth caller
#' 
#' Parses the x/y format of methylation calls, 
#' splitting them into individual columns: "methylCount" column for
#' number of methylated reads for site and "coverage" column for total
#' number of reads covering that site. Input files should have the following
#' columns: "chr", "start", "end", "meth", "rate", "strand".
#' 
#' 
#' This can run into memory problems if there are too many files...
#' because of the way parallel lacks long vector support. The solution is
#' to just use a single core; or to pass mc.preschedule = FALSE; This
#' makes it so that each file is processed as a separate job. Much better.
#' 
#' @param files a list of filenames (use parseInputArg if necessary)
#' @param contrastList Generally not needed for MIRA. 
#' A list of named character vectors, 
#' each with length equal to the number of items in files. 
#' These will translate into column names in the final table.
#' @param sampleNames   a vector of length length(files), name for each file. 
#' @param cores number of processors.
#' @param returnAsList Whether to return the output as a list 
#' or as one big data.table.
#' @return Data from each input file joined together into one big data.table.
#' If returnAsList = TRUE, then input from each file will be 
#' in its own data.table in a list.
#' @examples
#' shortBSDTFile = system.file("extdata", "shortRRBS.bed", package = "MIRA") 
#' shortBSDT = BSreadBiSeq(shortBSDTFile)
#'
#' @export
BSreadBiSeq = function(files, contrastList = NULL, 
                       sampleNames = tools::file_path_sans_ext(basename(files)), 
                       cores = 4, returnAsList = FALSE) {
    
    cores = min(length(files), cores); # not more cores than files!
    setLapplyAlias(cores);
    if (!is.null(contrastList)) {
        if (any(sapply(contrastList, length) != length(files))) {
            stop("contrastList must be a list, 
                 with each value having the same number of elements as files.");
        }
    }
    message("Reading ", length(files), " files..");
    freadList = lapplyAlias(files, fread, mc.preschedule = FALSE);
    colNames = names(contrastList)
    message("File reading finished (",
        length(files),
        " files). Parsing Biseq format...",
        appendLF = FALSE);
    # TODO: This parsing takes awhile, and could be done in parallel.
    freadListParsed = lapplyAlias(freadList, parseBiseq, mc.preschedule = FALSE)

    message("Parsing complete, building final tables and cleaning up...")
    numberOfFiles = length(freadListParsed);
    for (i in 1:numberOfFiles) {
        if (numberOfFiles > 1) {
            message(i, ": ", sampleNames[i], "; ", appendLF = FALSE)
        }
        if (numberOfFiles > 1 && i == numberOfFiles) {
            message("", appendLF = TRUE)
        }
        DT = freadListParsed[[i]]; # convenience alias.
        if (!is.null(contrastList)) {
            DT[, get("colNames") := as.list(sapply(contrastList, "[[", i))]
        }
        if (!is.null(sampleNames)) {
            DT[, sampleName := sampleNames[i]]
        }
        freadListParsed[[i]] = DT
    }

    # filteredList = do.call(rbind, freadListParsed)
    # gc(); # rbind call is memory inefficient; this helps.
    # rbindlist supposedly does the same thing as do.call(rbind, list) but 
    # faster
    # default (returnAsList = FALSE) is to return as 
    # one combined data.table/data.frame
    if (!returnAsList) { 
        filteredList = rbindlist(freadListParsed)
    }else{
        filteredList = freadListParsed
    }

    return(filteredList);
}

# Takes a data.table from BSreadBiSeq and parses the strange x/y format
# of methylation calls, splitting them into individual columns
# @param DT data.table to parse
# "chr", "start", "end", "meth", "rate", "strand" columns expected
# in that order.
# @return data.table with separate methylated and unmethylated columns.
# Specific col names are set
parseBiseq = function(DT) {
    message(".", appendLF = FALSE);
    setnames(DT, paste0("V", 1:6), 
             c("chr", "start", "end", "meth", "rate", "strand"))
    DT[, meth := gsub("'", "", meth)]
    # split the '12/12' format of meth calls
    ll = unlist(strsplit(DT$meth, "/", fixed = TRUE))
    idx = seq(1, length(ll), by = 2)
    DT[, `:=` (methylCount = as.integer(ll[idx]), 
              coverage = as.integer(ll[idx + 1]))]
    DT[, start := as.integer(start + 1)] # re-index
    DT[, c("rate", "end", "meth" ) := NULL] # remove unnecessary columns
    DT[, strand := NULL]
    DT = DT[, list(methylCount = sum(methylCount), coverage = sum(coverage)), 
          by = list(chr, start)] # smash measurements
    setcolorder(DT, c("chr", "start", "methylCount", "coverage"));
    DT = DT[ !grep("_", chr), ]; # clean Chrs
    return(DT)
}
