files <- unlist(snakemake@input)
d <- purrr::map_dfc(files, read.table,header=T,stringsAsFactors=F)
write.table(d, snakemake@output[[1]], sep="\t", quote=F, row.names=F)