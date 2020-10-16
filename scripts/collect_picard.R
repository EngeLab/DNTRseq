suppressPackageStartupMessages(library(tidyverse))
cellid_to_plateid <- function(x){
    sub("([A-Z]{3}[0-9]{5}).*","\\1", x)
}

# Get dup metrics
dup_files <- unlist(snakemake@input)
cat("Length:",length(dup_files))
names(dup_files) <- sub("(.*)\\-picard-mark_duplicates.txt", "\\1", basename(dup_files))
d.dup <- map_dfr(dup_files, read_tsv, skip=6, n_max=1, .id="dna_library_id", col_types="cnnnnnnnnnn") %>%
    mutate(plate_id=cellid_to_plateid(dna_library_id))

dim(d.dup)
head(d.dup)
write_tsv(d.dup, snakemake@output[[1]])
