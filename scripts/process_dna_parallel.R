# Process DNA counts in parallel for many binsizes and parameters in segmentation and other processing
suppressMessages(library(tidyverse))
library(DNAcopy)

# Collect parameters
input <- snakemake@input[[1]]
output <- snakemake@output[[1]]
binsize <- as.integer(snakemake@wildcards["binsize"])
readlen <- as.integer(snakemake@params["readlength"])
mincounts <- as.integer(snakemake@params[["min_count"]])
minwidth <- list(as.integer(snakemake@params[["cbs_minwidth"]]))
alpha <- list(as.numeric(snakemake@params[["cbs_alpha"]]))
undosplits <- list(snakemake@params[["cbs_undosplits"]])
sdprune <- list(as.integer(snakemake@params[["cbs_sdprune"]]))
gc_file <- snakemake@params[["gc_file"]]
bin_file <- snakemake@params[["bin_file"]]
bnd_file <- snakemake@params[["bnd_file"]]
removebadbins <- snakemake@params[["rm_outliers"]]
threads <- snakemake@threads[[1]]

run <- basename(input)
cat("[",run,"] Loading files\n")
data <- tibble(file=input,
               binsize=binsize,
               readlength=readlen,
               ref="hg38",
               raw=map(file, ~suppressMessages(read_tsv(.))), 
               bins=map_int(raw, nrow), 
               cells=map_int(raw, ncol),
               gk.gc=map(gc_file, ~read.table(., header=FALSE)),
               gk.loc=map(bin_file, ~read.table(., header=TRUE, sep="\t", as.is=TRUE)),
               gk.bnd=map(bnd_file, ~read.table(., header=F, sep="\t", as.is=TRUE)),
               gk.pos=map2(gk.bnd, bins, function(.x, .y){ cbind(c(1,.x[,2]), c(.x[,2], .y)) }))

any_empty_cells <- map(data$raw, function(.x){ 
    empty_bool <- lapply(.x, function(x){ all(is.na(x)) }) %>% unlist 
    empty_bool[empty_bool]
}) %>% unlist %>% names
if(length(any_empty_cells)){
    cat("[",run,"] WARNING! One or more samples had NULL/empty count matrices (will be filtered out):",paste(any_empty_cells,collapse=","),"\n")
}

any_na_cells <- map(data$raw, function(.x){ 
    na_bool <- lapply(.x, function(x){ any(is.na(x)) }) %>% unlist 
    na_bool[na_bool]
}) %>% unlist %>% names
if(length(any_na_cells)){
    cat("[",run,"] WARNING! One or more samples had NA entries (will be filtered out):",paste(any_na_cells,collapse=","),"\n")
}

cat("[",run,"] Filter by min count", mincounts, "\n")
# check_empty[check_empty]
cat("[",run,"] Before filtering: n=", data$cells[[1]], "\n")
data <- data %>% 
    mutate(
        filtered=map(raw, function(.x){
            keep <- colSums(.x,na.rm=T) > mincounts & ! names(.x) %in% any_empty_cells & ! names(.x) %in% any_na_cells
            .x[,keep]
        })
    )
cat("[",run,"] After filtering: n=", ncol(data$filtered[[1]]), "\n")

if(removebadbins==1){
    cat("[",run,"] Calculate outlier bins\n")
    data <- data %>% 
        mutate(
            bb=map(filtered, function(.x){
            pr <- .x %>% 
            sweep(2, colSums(.x), "/") %>% 
            rowMeans()
            qt <- quantile(pr, probs=0.995)
            which(pr > qt)
            }))
    # Filter out bad bins
    cat("[",run,"] Filter out bad bins\n")
    filter_rows <- function(df,filter_vector){df[-filter_vector,]}
    data <- data %>% 
        mutate(
            filtered = map2(filtered, bb, ~filter_rows(.x, .y)),
            bins_filt = map2_int(bins, bb, function(.x, .y){ .x - length(.y) }), # Recalc bin count after filtering
            gk.gc = map2(gk.gc, bb, ~filter_rows(.x, .y) %>% unlist %>% unname),
            gk.loc = map2(gk.loc, bb, ~filter_rows(.x, .y))
        )
}

# Add stats
cat("[",run,"] Calculate overall statistics\n")
data <- data %>% 
    mutate(
        stats=map(filtered, function(.x){
            tibble(
                cell=colnames(.x),
                cell_id=sub("([A-Z]{3}[0-9]{5})(D?)_([A-Z][0-9]{2}).*","\\1.\\3", cell),
                reads=sapply(.x, sum),
                bins=dim(.x)[1],
                mean=sapply(.x, function(x){ round(mean(x), digits=2) }),
                var=sapply(.x, function(x){ round(sd(x), digits=2) }),
                disp=round(var/mean, digits=2),
                min=sapply(.x, min),
                `25th`=sapply(.x, function(x){ quantile(x, c(0.25))[[1]] }),
                median=sapply(.x, median),
                `75th`=sapply(.x, function(x){ quantile(x, c(0.75))[[1]] }),
                max=sapply(.x, max)
            )})
    )

### GC normalize and calculate log2 ratio
cat("[",run,"] Normalizing and log2 transformation\n")
timestart <- Sys.time()
lowess.gc = function(x, y) {
    low = lowess(x, log(y), f=0.05);
    z = approx(low$x, low$y, x)
    return(exp(log(y) - z$y))
}
data <- data %>% ungroup %>% 
    mutate(
        normal=map(filtered, function(.x){ sweep(.x+1, 2, colMeans(.x+1), '/') }),
        normal.gc=map2(filtered, gk.gc, ~apply(.x, 2, function(x){ lowess.gc(.y[[1]], x+1/mean(x+1)) })),
        lrr=map(normal.gc, ~log2(.x))
    )
cat("[",run,"] Finished normalization\n")
diff(c(timestart, Sys.time()))

### Segmentation with CBS
data <- data %>% 
    mutate(
        cna.obj=map2(lrr, gk.loc, function(.x,.y){ CNA(genomdat=.x, sampleid=colnames(.x), chrom=ordered(.y[,"CHR"], levels=paste0("chr", c(1:22, "X", "Y"))), maploc=as.numeric(.y[,"END"]), data.type="logratio") }),
        cna.smooth=map(cna.obj, ~smooth.CNA(.x, smooth.region=10, outlier.SD.scale=4, smooth.SD.scale=2, trim=0.025))
    )
# NOTE1! CBS algo applies its own trimming - on top of the raw filtering we do, when calling segments
# NOTE2! Have to make the chromosomes an ordered variable to get correct (natural) order of chr naming
# Collect parameters in tibble, expanded to even lenghts (becomes list of lists)
seg_params <- expand.grid(minwidth=minwidth, alpha=alpha, undosplits=undosplits, sdprune=sdprune) %>% mutate(param_id=row_number())
cat("[",run,"] Segmentation parameters:\n")
print(as.data.frame(seg_params))
# Nested calls (eg for loop within for loop) where:
# 1. Outer loop is over `cna.smooth` list 
# 2. Inner loop is over the parameter tibble
# To expose the outer variable to the innermost function call (segment) .x is assigned to .y
# Run time ~10 mins for 250k/500k bins with 2 parameter sets
timestart <- Sys.time()
cat("[",run,"] Starting segmentation\n")
data <- data %>% ungroup %>% 
    mutate(res=map(cna.smooth, ~pmap(seg_params, function(.y, minwidth, alpha, undosplits, sdprune, ...){
                       segment(.y, verbose=0, min.width=minwidth, alpha=alpha, undo.splits=undosplits, undo.SD=sdprune)
                   },.y=.x)),
           segparams=list(seg_params)) # To keep with bindata when saving
cat("[",run,"] Finished segmentation\n")
diff(c(timestart, Sys.time()))

# Collect results (still working on two levels)
# data$segs_list[[binsize_id]][[param_id]][[cell_id]]
# data$segs_list2[[binsize_id]][[cell_id]] - all param_ids in same df for ~faceting
seg_params_j <- seg_params %>% mutate(param_id = as.character(param_id)) # To allow joing in segmentation parameters
data <- data %>% 
    mutate(
        segs=map(res, function(.x){
            map_dfr(.x, function(.y){ cbind(.y$output, .y$segRows) }, .id="param_id")}),
        segs_list=map(res, function(.x){
            map(.x, function(.y){ d=cbind(.y$output, .y$segRows); split(d, d$ID) })}),
        segs_list2=map(res, function(.x){
            d=map_dfr(.x, function(.y){ cbind(.y$output, .y$segRows) }, .id="param_id"); d=left_join(d, seg_params_j, by="param_id"); split(d, d$ID) }),
        segstats=map(res, function(.x){ map_dfr(.x, ~segments.summary(.x), .id="param_id") })
    )

### Save data
cat("[",run,"] Saving data\n")
saveRDS(object=data, file=output)
