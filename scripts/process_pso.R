suppressPackageStartupMessages(library(tidyverse))
library(parallel)
source("scripts/functions.R") # includes modified PSO algorithm

file=snakemake@input[[1]]
outfile=snakemake@output[[1]]
threads=snakemake@threads[[1]]

cat("Loading", file,"\n")
cat("Outfile:",outfile,"\n")
cat("Threads:",threads,"\n")
dnaobj <- read_rds(file)

# Helper functions
dna_to_cellid <- function(x){
    sub("([A-Z]{3}[0-9]{5})(D?)_([A-Z][0-9]{2}).*","\\1.\\3", x)
}
rna_to_cellid <- function(x){
    sub("([A-Z]{3}[0-9]{5})(R?).([A-Z][0-9]{2}).*","\\1.\\3", x)
}
cellid_to_plateid <- function(x){
    sub("([A-Z]{3}[0-9]{5}).*","\\1", x)
}
# Calculates distance between actual value and predicted ploidy. Factor is the scaling factor (to be optimized), and segs is per-bin values before scaling
euclid.dist <- function(factor, segs) {
    sum(((segs*factor)-round(segs*factor))^2)
}
# Helper function to transform raw log values into predicted ploidy using the current tables, etc.
name.to.ploidy <- function(cellname, bigtable) {
    cat("^Processing",cellname,"..")
    start_time <- Sys.time()
    o <- bigtable$cell_id == cellname # Select correct segments for cell
    segs <- rep(2^bigtable[o,"seg.median"], bigtable[o, 'num.mark']) # Make vector of values per-bin instead of segments (and de-log2)
    res <- ploidy.by.bin(segs)
    end_time <- Sys.time()
    cat("done in",end_time - start_time,"s\n")
    return(res)
}
# Actual function to transform raw vals to ploidy
ploidy.by.bin <- function(cell.segs) {
    res <- psoptim1(NA, fn=euclid.dist, segs=cell.segs, lower=1.5, upper=3, 
                    Xinit = matrix(seq(1.5, 3, length.out=100), nrow=1),
                    control =  list(maxit = 1000, s = 100))
    vals <- round(cell.segs * res$par)
    resid <- (cell.segs*res$par-round(cell.segs*res$par))^2 # Residuals
    par <- res$par
    return(list(vals,resid,par))
}

# Filter objects
d <- dnaobj %>% 
    unnest(segstats) %>%
    rename(dna_library_id=ID) %>% 
    mutate(cell_id=dna_to_cellid(dna_library_id),
           plate_id=cellid_to_plateid(cell_id))
    #filter(dna_to_cellid(dna_library_id) %in% top10) # DEBUG FILTER

distinct(d, cell_id, .keep_all=T) %>% select(plate_id) %>% table
samples <- unique(d$cell_id)
dnaids <- unique(d$dna_library_id)
plates <- d %>% distinct(cell_id, .keep_all = T) %>% select(plate_id) %>% unlist(use.names = F)

# Keep LRR object together for later plotting, subsetting for same samples as above
dnalrr <- dnaobj$lrr[[1]][,dnaids]
colnames(dnalrr) <- dna_to_cellid(colnames(dnalrr)) # Convert DNA library IDs to cell_ids to keep consistency

dnaraw <- dnaobj$filtered[[1]][,dnaids]
colnames(dnaraw) <- dna_to_cellid(colnames(dnaraw)) # Convert DNA library IDs to cell_ids to keep consistency

dnanorm <- dnaobj$normal.gc[[1]][,dnaids]
colnames(dnanorm) <- dna_to_cellid(colnames(dnanorm)) # Convert DNA library IDs to cell_ids to keep consistency

# Run pso optimization for each sample separately, to get scaling factor. Then scale and round to nearest integer!
cat("[",outfile,"] ^^Running PSO optimization\n")
dnaseg <- dplyr::select(d, chrom:ncol(d))
psoptim.list <- mclapply(samples, name.to.ploidy, bigtable=as.data.frame(dnaseg), mc.cores=threads, mc.silent=F)
psoptim <- as.matrix(as.data.frame(lapply(psoptim.list, `[[`, 1), col.names=samples)) # PSO matrix
psoptim.resid <- as.matrix(as.data.frame(lapply(psoptim.list, `[[`, 2), col.names=samples)) # Residuals
psoptim.par <- as.matrix(as.data.frame(lapply(psoptim.list, `[[`, 3), col.names=samples)) # Scaling parameter

save(list=c("psoptim","psoptim.resid","psoptim.par","dnaseg","dnalrr","dnaraw", "dnanorm", "plates","samples"),file=outfile)