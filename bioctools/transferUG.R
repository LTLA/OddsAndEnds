args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
    stop("need source and destination paths")
}
first <- args[1]
second <- args[2]

copied <- file.copy(first, second, overwrite=TRUE)
if (!copied) {
    stop("failed to copy to destination")
}

while (1) {
    out <- tools::compactPDF(second, gs_quality="ebook")
    if (nrow(out)==0 || out$old[1]==out$new[1]) { break }
}
