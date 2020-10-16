suppressPackageStartupMessages(library(tidyverse))

files <- snakemake@input[["cells"]]
cat("^^^ Processing", length(files), "RNA files\n")
head(files)

fx <- function(x){
	info = file.info(x)
	if(file.exists(x) & info$size!=0){
  		read.table(
    		x, sep = "\t", colClasses = c("NULL","integer"), 
    		col.names = c("NULL", sub("(.*)\\.htseq\\.out", "\\1", basename(x)))
  		)
	} else {
  		cat("Missing",x,"\n")
	}
}

genes <- read.table(files[1], sep = "\t", colClasses = c("character", "NULL"), col.names = c("HGN", "NULL"))
head(genes)

d <- map_dfc(files, fx)
dim(d)

d <- bind_cols(genes, d)
dim(d)

write_tsv(d, snakemake@output[["counts"]])

# Write out basic RNA QC too
ercc.idx <- grepl("^ERCC-",d$HGN)
htseq.idx <- grepl("^__",d$HGN)
ercc.frac <- colSums(d[ercc.idx,-1])/colSums(d[!htseq.idx,-1])
filtered <- as.matrix(d[!htseq.idx & !ercc.idx,-1])

rna_qc <- enframe(ercc.frac,name="rna_library_id",value="ercc_frac") %>% 
  mutate(rna_counts=colSums(filtered),
         rna_features=colSums(apply(filtered, 2, function(x) x > 0)))

write_tsv(rna_qc, snakemake@output[["qc"]])
