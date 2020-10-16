# Plot specific cells using DNAcopy internals
options(scipen=999)
library(tidyverse)
library(scales)
library(cowplot)
library(gridGraphics)
# For QC plots
library(ggExtra)
library(ggrepel)
library(gtools)

# Setup
source(snakemake@params[["r_global"]])
print(load(snakemake@input[["pso"]]))

# Metadata
dna_stats <- tibble(
    cell_id=colnames(dnaraw),
    bin_counts=sapply(dnaraw, sum),
    bins=dim(dnaraw)[1],
    bin_mean=sapply(dnaraw, function(x){ round(mean(x), digits=2) }),
    bin_var=sapply(dnaraw, function(x){ round(sd(x), digits=2) }),
    bin_disp=round(bin_var/bin_mean, digits=2),
    bin_0=sapply(dnaraw, min),
    bin_25=sapply(dnaraw, function(x){ quantile(x, c(0.25))[[1]] }),
    bin_median=sapply(dnaraw, median),
    bin_75=sapply(dnaraw, function(x){ quantile(x, c(0.75))[[1]] }),
    bin_max=sapply(dnaraw, max)
)
    
meta <- read_tsv(snakemake@input[["meta"]]) %>% 
    left_join(dna_stats, by="cell_id")
bins <- read_tsv(snakemake@params[["bin_file_bed"]],col_names = c("chr","start","end")) %>% 
    mutate(bin=seq_along(chr), end_cum=cumsum((end-start)+1))
chr_bounds <- bins %>% 
    group_by(chr) %>% summarize(min=min(bin),max=max(bin),chrlen_bp=max(end))
chr_bounds <- chr_bounds[gtools::mixedorder(chr_bounds$chr),] %>% 
    mutate(mid=round(min+(max-min)/2,0),end_bp=cumsum(chrlen_bp), start_bp=end_bp-chrlen_bp, mid_bp=round((chrlen_bp/2)+start_bp,0))

# Meta ordered by PSO and other DNA matrices
m.sort <- meta[match(colnames(psoptim), as.character(meta$cell_id)),]
m.sort$pso_par <- unlist(psoptim.par[1,]) # Scaling factor from PSO

# Residuals of bin-level vs integer copy number
lrr_delog <- lapply(colnames(psoptim.par), function(x) 2^dnalrr[,x]*psoptim.par[1,x]) # de-log and scale by PSO factor
lrr.m <- matrix(unlist(lrr_delog, use.names=T), ncol = ncol(psoptim.par), byrow = TRUE) # make matrix
colnames(lrr.m) <- colnames(psoptim.par)
lrr_residuals <- (lrr.m - psoptim)^2 # Calculate squared residuals

# DNA damage
pso.med <- apply(psoptim[,m.sort$source=='Control'], 1, median) # Median profile of control cells
psoptim.rel <- apply(psoptim, 2, function(x) {x-pso.med}) # Calculate relative difference 
save("psoptim.rel",file=snakemake@output[["pso_relative_ctrl"]])
damage.burden <- colSums(psoptim.rel != 0)/nrow(psoptim.rel) # Calculate damage as fraction of bins different from control
damage.burden.df <- damage.burden %>% enframe("cell_id","dna_damage") 
write_tsv(damage.burden.df, snakemake@output[["damage"]])
#plot(sort(damage.burden))

# Make damage rank
damage.rank.unsorted <- enframe(sort(damage.burden,decreasing=T),name="cell_id",value="damage") %>% rownames_to_column("rank") %>% mutate(rank=as.integer(rank), damage_bin=ntile(damage,100))
damage.rank <- damage.rank.unsorted[match(m.sort$cell_id, damage.rank.unsorted$cell_id),]

# Probability of loss or gain by condition
all.cond <-  c("Control","X-ray 1Gy", "X-ray 3Gy", "X-ray 5Gy", "ETO 1uM", "ETO 3uM", "ETO 10uM", "ETO 30uM", "ETO 100uM")
p.loss.df <- as.data.frame(lapply(all.cond, function(i) {
    apply(psoptim.rel[,m.sort$source == i & m.sort$sort_gate=='G1'], 1, function(x) { sum(x<0)/length(x)})
}))
names(p.loss.df) <- all.cond
write_tsv(p.loss.df, snakemake@output[["p_loss"]])
p.gain.df <- as.data.frame(lapply(all.cond, function(i) {
    apply(psoptim.rel[,m.sort$source == i & m.sort$sort_gate=='G1'], 1, function(x) { sum(x>0)/length(x)})
}))
names(p.gain.df) <- all.cond
write_tsv(p.gain.df, snakemake@output[["p_gain"]])

# Integer plots
cn_high <- rep("red4",95)
names(cn_high) <- seq.int(6,100)
cn_colors <- c("0"="darkblue","1"="lightblue","2"="darkgrey","3"="red1","4"="red2","5"="red3",cn_high)

# Plot functions
plot_cell <- function(cell.idx,ymax=T,coord_bp=T,chr_lines=T,chr_labels=T,residuals=T,residuals_level="segment",plot=T,scale_factor="pso",show_lrr_segs=T){
    if(!exists("psoptim") | !exists("dnalrr") | !exists("psoptim.resid")) stop("Missing required data: psoptim, psoptim.resid and dnalrr")
    if(!is.numeric(cell.idx)) cell.idx <- which(m.sort$cell_id==cell.idx); if(length(cell.idx)==0) stop("No match for cell id in metadata")
    
    if(coord_bp) {
        coord_x = "end_cum"
        coord_x_chr = "end_bp"
        coord_x_chrLab = "mid_bp"
    } else { 
        coord_x = "bin" 
        coord_x_chr = "max"
        coord_x_chrLab = "mid"
    }
    
    cell.label <-  m.sort$cell_id[cell.idx]
    if(scale_factor=="pso"){ 
        rescaled=FALSE
        scale_factor <- m.sort$pso_par[cell.idx] 
        pso_integers <- psoptim[,cell.idx]
        pso_residuals <- psoptim.resid[,cell.idx]
        pso_residuals_lrr <- ((2^dnalrr[,cell.idx]*scale_factor)-psoptim[,cell.idx])^2
    } else if(is.numeric(scale_factor)){ 
        rescaled=TRUE
        cat("^^Rescaling",cell.label,"with factor",scale_factor," (original",m.sort$pso_par[cell.idx],")\n")
        cell.idx.seg <- which(dnaseg$cell_id==cell.label)
        cell_seg_bins <- rep(2^unlist(dnaseg[cell.idx.seg,"seg.median"],use.names=F), unlist(dnaseg[cell.idx.seg, 'num.mark'],use.names=F))
        pso_integers <- round(cell_seg_bins * scale_factor)
        pso_residuals <- (cell_seg_bins*scale_factor-round(cell_seg_bins*scale_factor))^2 # Residuals
        cat("^^Residual sum",sum(pso_residuals)," vs original", sum(psoptim.resid[,cell.idx]),"\n")
        pso_residuals_lrr <- ((2^dnalrr[,cell.idx]*scale_factor)-psoptim[,cell.idx])^2
    } else { stop("Scale factor",scale_factor,"not recoginzed\n") }
    cell.label_ext <- paste0(m.sort$source[cell.idx],", ",round(m.sort$dna_read_pairs[cell.idx]/1e6,digits=2),"M, ",
                             round(m.sort$dna_dups_pct[cell.idx]*100,2),"% dups, ",round(m.sort$bin_mean[cell.idx],0)," reads/bin, ",
                             m.sort$bin_disp[cell.idx]," disp, ",m.sort$bin_var[cell.idx]," var, sf ",round(scale_factor,3),", ", round(damage.burden[cell.idx],3)," damage (",damage.rank$rank[cell.idx],"/",nrow(damage.rank),")")
    pso <- pso_integers %>% enframe(name="bin",value="cn") %>% left_join(bins, by="bin")
    lrr <- dnalrr[,cell.idx] %>% enframe(name="bin",value="log2") %>% left_join(bins, by="bin")
    ptop <- ggplot(pso, aes(!!ensym(coord_x), cn, color=as.character(cn))) + 
        geom_hline(yintercept=2, linetype="dashed", color="grey") +
        geom_point(data=lrr,aes(y=2^log2*scale_factor),color="black",size=0.5,alpha=0.2,pch=19,show.legend=F) +
        geom_point(size=0.8, alpha=0.5,pch=19,show.legend=F) +
        scale_color_manual(values=cn_colors) +
        scale_y_continuous(labels=comma_format(accuracy=1), breaks=pretty_breaks(5)) +
        scale_x_continuous(expand=c(0,0)) +
        labs(y="Copy number") +
        annotate("text",x=-Inf,y=Inf,label=cell.label,hjust=-0.05,vjust=1.5,color="grey20") +
        annotate("text",x=Inf,y=Inf,label=cell.label_ext,hjust=1.01,vjust=1.5,color="grey20") +
        theme_dntr() +
        theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank(), plot.margin = margin(5,5,0,5,unit="pt"))
    if(show_lrr_segs){
        seg <- dnaseg[dnaseg$cell_id==cell.label,] %>% 
            group_by(chrom) %>% mutate(max=max(loc.end)) %>% group_by(chrom,max) %>% nest() %>% ungroup() %>% 
            mutate(sum=cumsum(max),prev_sum=lag(sum,default=0)) %>% unnest(data) %>% 
            mutate(new_start=prev_sum+loc.start,new_end=prev_sum+loc.end)
        ptop <- ptop + geom_segment(data=seg, aes(x=new_start, xend=new_end,y=2^seg.median*scale_factor,yend=2^seg.median*scale_factor),color="orange",size=2,alpha=0.5)
    }
    if(!is.null(ymax) && is.numeric(ymax)) {
        ptop <- ptop + coord_cartesian(ylim=c(0,ymax))
    } else if(!is.null(ymax) && ymax==TRUE) {
        ymax <- max(pso$cn) + 1
        ptop <- ptop + coord_cartesian(ylim=c(0,ymax))
    }
    if(chr_lines){
        if(!exists("chr_bounds")) stop("Missing chr/bin coordinates")
        ptop <- ptop + geom_vline(data=chr_bounds,aes(xintercept=!!ensym(coord_x_chr)),color="grey80", alpha=0.5)
    }
    if(chr_labels){
        if(!exists("chr_bounds")) stop("Missing chr/bin coordinates")
        ptop <- ptop + geom_text(data=chr_bounds,aes(x=!!ensym(coord_x_chrLab),y=-Inf,label=sub("chr","",chr)),vjust = -0.5,hjust="center",inherit.aes = F)
    }
    
    pso.res <- pso_residuals %>% enframe(name="bin",value="res") %>% left_join(bins, by="bin")
    pso.res.lrr <- pso_residuals_lrr %>% enframe(name="bin",value="res") %>% left_join(bins, by="bin")
    pbot <- pso.res %>% 
        ggplot(aes(!!ensym(coord_x), res)) + 
        geom_line(size=0.8, alpha=1) +
        labs(y="r seg") + theme_dntr() +
        theme(plot.margin=margin(0,0,0,0,unit="pt"), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_blank()) +
        scale_y_continuous(breaks=pretty_breaks(2)) +
        scale_x_continuous(expand=c(0,0))
    pbot2 <- pso.res.lrr %>% 
        ggplot(aes(!!ensym(coord_x), res)) + 
        geom_line(size=0.8, alpha=1) +
        labs(y="r bin") + theme_dntr() +
        theme(plot.margin=margin(0,0,0,0,unit="pt"), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_blank()) +
        scale_y_continuous(breaks=pretty_breaks(2)) +
        scale_x_continuous(expand=c(0,0))
    if(chr_lines){
        if(!exists("chr_bounds")) stop("Missing chr/bin coordinates")
        pbot <- pbot + geom_vline(data=chr_bounds,aes(xintercept=!!ensym(coord_x_chr)),color="grey80", alpha=0.5)
    }
    if(coord_bp){
        pbot <- pbot + scale_x_continuous(expand=c(0,0), labels=comma_format(accuracy = 100,scale=1e-6, suffix="M", big.mark=" "), breaks=pretty_breaks(10))
        pbot2 <- pbot2 + scale_x_continuous(expand=c(0,0), labels=comma_format(accuracy = 100,scale=1e-6, suffix="M", big.mark=" "), breaks=pretty_breaks(10))
    }
    # if(residuals==F & coord_bp==T){
    #     ptop <- ptop + scale_x_continuous(labels=comma_format(accuracy = 100,scale=1e-6, suffix="M", big.mark=" "), breaks=pretty_breaks(10)) + theme(axis.text.x=element_text(),axis.ticks.x=element_line())
    # }
    
    if(residuals) p <- plot_grid(ptop,pbot,pbot2,ncol=1,align="v",rel_heights = c(1,0.3,0.3)) else p <- ptop
    if(plot) plot(p) else return(p)
}
plot_raw_gg <- function(cell,ylim=c(-3,2),plot=F,coord_bp=T,show_lrr_seg=F){
    if(coord_bp) {
        coord_x = "end_cum"
        coord_x_chr = "end_bp"
        coord_x_chrLab = "mid_bp"
    } else { 
        coord_x = "bin" 
        coord_x_chr = "max"
        coord_x_chrLab = "mid"
    }
    raw_median <- median(unlist(dnaraw[,cell]))
    p1 <- dnaraw[,cell] %>% unlist(use.names=F) %>% enframe(name="bin",value="reads") %>% 
        ggplot(aes(bin,reads)) + geom_hline(yintercept=raw_median,linetype="dashed",color="grey20") +
        geom_point(size=0.5,alpha=0.5) + labs(subtitle="Raw reads",caption=cell) + theme_dntr() + theme(axis.title = element_blank())
    
    pso <- psoptim[,cell] %>% enframe(name="bin",value="cn") %>% left_join(bins, by="bin")
    lrr <- dnalrr[,cell] %>% enframe(name="bin",value="lrr") %>% left_join(bins, by="bin")
    
    # Compress outliers at ylim borders
    lrr_outl <- lrr %>% filter(lrr < min(ylim) | lrr > max(ylim)) %>% 
        mutate(lrr_pch=case_when(lrr < min(ylim) ~ "bot", lrr > max(ylim) ~ "top"),
               lrr=case_when(lrr < min(ylim) ~ min(ylim), lrr > max(ylim) ~ max(ylim))) 
    p2 <- dnalrr[,cell] %>% enframe(name="bin",value="lrr") %>% 
        filter(lrr > min(ylim), lrr < max(ylim)) %>% 
        ggplot(aes(bin,lrr)) + geom_point(size=0.5,alpha=0.5) + 
        geom_point(data=lrr_outl,aes(shape=lrr_pch,color=lrr_pch), alpha=0.5, show.legend = F,size=0.5) +
        geom_hline(yintercept=c(-1,0,0.5,1),linetype="dashed",color="grey20") + 
        scale_shape_manual(values=c("top"=2,"bot"=6))+
        scale_color_manual(values=c("top"="red","bot"="blue","NA"="black")) +
        coord_cartesian(ylim=ylim) +
        theme_dntr() + labs(subtitle="LRR",caption=cell) + theme(axis.title = element_blank())
    if(show_lrr_seg){
        seg <- dnaseg[dnaseg$cell_id==cell,] %>% 
            group_by(chrom) %>% mutate(max=max(loc.end)) %>% group_by(chrom,max) %>% nest() %>% 
            mutate(sum=cumsum(max),prev_sum=lag(sum,default=0)) %>% unnest() %>% 
            mutate(new_start=prev_sum+loc.start,new_end=prev_sum+loc.end)
        p2 <- p2 + geom_segment(data=seg, aes(x=new_start, xend=new_end,y=2^seg.median*scale_factor,yend=2^seg.median*scale_factor),
                                color="orange",size=2,alpha=0.5)
    }
    p3 <- psoptim[,cell] %>% enframe(name="bin",value="cn_int") %>% 
        ggplot(aes(bin,cn_int)) + geom_point(size=0.5,alpha=0.5) +
        geom_hline(yintercept=2,linetype="dashed",color="grey20") +
        theme_dntr() + labs(subtitle="PSO integers",caption=cell) + theme(axis.title = element_blank())
    return(plot_grid(p1,p2,p3,nrow=1))
}

# Print PDFs
pdf(snakemake@output[["fig_top20"]],width=12,height=5)
for(cell in names(sort(damage.burden,decreasing = T))[1:10]){
    p1 <- plot_cell(cell,ymax=T,plot=F)
    p2 <- plot_raw_gg(cell,plot=F)
    print(plot_grid(p1,p2,ncol=1,rel_heights=c(1.5,1),axis = "lr"))
}
dev.off()
pdf(snakemake@output[["fig_bot20"]],width=12,height=5)
for(cell in names(sort(damage.burden, decreasing=T))[1:10]){
    p1 <- plot_cell(cell,ymax=T,plot=F)
    p2 <- plot_raw_gg(cell,plot=F)
    print(plot_grid(p1,p2,ncol=1,rel_heights=c(1.5,1),axis = "lr"))
}
dev.off()

### DNA QC and plotting
sum_lrr_residuals <- colSums(lrr_residuals) # Sum per cell
robust_lrr_residuals <- apply(lrr_residuals, 2, function(x) { sum(sort(x, decreasing=T)[-1:-10]) })
median_robust_lrr_residuals <- apply(lrr_residuals, 2, function(x) { median(sort(x, decreasing=T)[-1:-10]) })
#
bin0 <- colSums(psoptim == 0)
raw0 <- colSums(dnaraw %>% as.matrix == 0)
raw0_frac <- raw0/nrow(bins)
totRes <- colSums(psoptim.resid)
#
nonzero_mean <- apply(psoptim, 2, function(x) { mean(x[x!=0]) })
zero <- apply(psoptim, 2, function(x) {mean(x==0)})
nonzero_bins <- psoptim!=0
nonzero_median_lrr_residuals <- c()
for(i in 1:ncol(lrr_residuals)){
    nonzero_median_lrr_residuals[i] <- median(sort(lrr_residuals[i,][nonzero_bins[i,]],decreasing=T)[-1:-10])
}

qc_ext <- tibble(cell_id=colnames(psoptim),
                 nonzero_mean=nonzero_mean,
                 zero_frac=zero,
                 bin0=bin0,
                 raw0_frac=raw0_frac,
                 sum_res_sq_seg=totRes,
                 sum_res_sq_lrr=sum_lrr_residuals,
                 sum_res_robust=robust_lrr_residuals,
                 med_res_robust=median_robust_lrr_residuals,
                 nonzero_med_res=nonzero_median_lrr_residuals)

### Add new DNA filters to metadata
meta_qc <- meta %>% 
    left_join(qc_ext,by="cell_id") %>% 
    mutate(dnaQC_reads=dna_read_pairs > snakemake@params[["min_count"]],
           dnaQC_dups=dna_dups_pct < snakemake@params[["max_dup_frac"]],
           dnaQC_residual=med_res_robust < quantile(med_res_robust,0.995,names=F,na.rm=T),
           dna_pass_qc=dnaQC_reads & dnaQC_dups & dnaQC_residual)

## Make cross-table of DNA qc exclusions
meta_qc %>% select(starts_with("dnaQC"),dna_pass_qc,sort_count) %>% 
    mutate(wells=TRUE,
           sorted=sort_count>0) %>% 
    select(-sort_count) %>% 
    gather(var,val) %>% 
    table(useNA = "always")

qc_cols <- c("TRUE"=rgb(0,0,0,0.5),"FALSE"=rgb(1,0,0,0.5))
d <- select(meta_qc, dna_read_pairs, nonzero_mean, zero_frac, sum_res_robust, med_res_robust)
pdf(snakemake@output[["fig_crossplot"]])
pairs(d,pch=21,cex=0.4,lower.panel=NULL,col=qc_cols[as.character(meta_qc$dna_pass_qc)],log="xy")
dev.off()

p_reads_vs_residuals <- meta_qc %>% 
    ggplot(aes(dna_read_pairs, med_res_robust, color=dna_pass_qc)) + geom_point(size=0.8,pch=19) +
    geom_point(data=filter(meta_qc, dnaQC_dups==FALSE), size=1.8,pch=21, color="blue") +
    scale_y_log10() +
    scale_x_log10() +
    geom_hline(yintercept=quantile(median_robust_lrr_residuals,0.995), color="grey", linetype="dotted") +
    geom_vline(xintercept=20000, color="grey", linetype="dotted") +
    scale_color_manual(values=c("TRUE"="black","FALSE"="red")) +
    theme_dntr()
ggsave(snakemake@output[["fig_filterScatter"]], p_reads_vs_residuals, width=6,height=4)

# Plot by worst robust residuals
pdf(snakemake@output[["fig_badResiduals"]],width=12,height=6)
for(cell in names(sort(median_robust_lrr_residuals,decreasing=T)[1:25])){
    p1 <- plot_cell(cell,ymax=T,plot=F)
    p2 <- plot_raw_gg(cell,plot=F)
    print(plot_grid(p1,p2,ncol=1,rel_heights=c(1.5,1),axis = "lr"))
}
dev.off()

# Plot damage boxplot
dmg.df.filt <- damage.burden.df %>% 
    left_join(meta_qc, by="cell_id") %>% 
    filter(dna_pass_qc==TRUE) %>% 
    mutate(source=factor(source, levels=c("Control",
                                          "ETO 1uM", "ETO 3uM", "ETO 10uM", "ETO 30uM", "ETO 100uM",
                                          "X-ray 1Gy", "X-ray 3Gy", "X-ray 5Gy")))
p_dmg <- dmg.df.filt %>% 
    ggplot(aes(source, dna_damage, color=source)) + geom_boxplot(show.legend=F) + 
    facet_grid(~timepoint, space="free_x", scales = "free_x") +
    theme_dntr() +
    theme(axis.title.x = element_blank(), axis.text.x=element_text(angle=90,vjust=0.5, hjust=1)) +
    labs(y="DNA damage fraction", caption=paste0("n=",nrow(dmg.df.filt))) +
    scale_color_manual(values=source_cols)
ggsave(snakemake@output[["fig_damageBoxplot"]], p_dmg, width=4.5,height=4.5,useDingbats=F)

# Plot mapped read pairs vs duplicate fraction for supplementary
tab <- table(meta_qc[meta_qc$dna_pass_qc,"source"])[names(source_cols)]
labs <- paste0(names(tab)," (n=",tab,")")

p_mappedvsdup_pass <- meta_qc %>% 
    filter(dna_pass_qc) %>% 
    mutate(mapped_readpairs=((dna_read_pairs*2+dna_unpaired_reads)-dna_unmapped)/2,
           source_broad=case_when(grepl("X-ray",source)~"X-ray",
                                  grepl("ETO",source)~"ETO",
                                  TRUE ~ "Control"),
           source_broad=factor(source_broad,levels=c("Control","ETO","X-ray")),
           source=factor(source,levels=names(source_cols))) %>% 
    ggplot(aes(mapped_readpairs, dna_dups_pct, color=source)) +
    geom_point(pch=19, alpha=0.5, size=0.8) +
    theme(panel.border = element_rect(fill = NA, colour = "black", size=1),
          rect = element_blank(), axis.line=element_blank(),
          legend.position="right", legend.justification = "top", legend.title=element_blank(),
          legend.direction = "vertical", legend.margin=margin(0,0,0,0,unit="pt"), 
          axis.text=element_text(size=8)) +
    scale_x_continuous(breaks=pretty_breaks(3), labels=unit_format(unit="M",scale=1e-6,sep="",accuracy=.1)) +
    scale_color_manual(labels=labs, values=source_cols) +
    labs(y="Duplicate fraction",x="Mapped read-pairs") +
    facet_wrap(~source_broad, nrow=1,scales="free_x") +
    guides(color=guide_legend(override.aes=list(size=2,title.position="top",title.hjust=0,alpha=1)))
ggsave(snakemake@output[["figext_mappedVSdup"]], p_mappedvsdup_pass, width=7, height=3)    

# Plot example HCT116 cells at different read depths
p70 <- plot_cell("VZA01712.G02",plot=F,residuals=F)
p70ext <- plot_raw_gg("VZA01712.G02")
p70g <- plot_grid(p70,p70ext,nrow=2)

p200 <- plot_cell("VZA01702.J20",plot=F,residuals=F)
p200ext <- plot_raw_gg("VZA01702.J20")
p200g <- plot_grid(p200,p200ext,nrow=2)

p1m <- plot_cell("VZA01701.J10",plot=F,residuals=F)
p1mext <- plot_raw_gg("VZA01701.J10",coord_bp=T)
p1mg <- plot_grid(p1m,p1mext,nrow=2)
pdf("fig/dna_example_profiles_70k200k1m.pdf",width=12,height=6)
plot(p70g)
plot(p200g)
plot(p1mg)
dev.off()



# Print out metadata table, with extended DNA QC filters
write.table(meta_qc, snakemake@output[["metadata_filtered"]], sep="\t", quote=F, row.names=F)

# Filter psoptim files and save
#meta <- read.table(snakemake@input[["meta"]], sep="\t", header=T, stringsAsFactors=F)
#meta.o <- meta[match(colnames(psoptim), as.character(meta$cell_id)),] # Order as PSO matrices
cells.f <- meta_qc$cell_id[meta_qc$dna_pass_qc] # Get vector of DNA libraries that pass QC

psoptim <- psoptim[,cells.f]
psoptim.resid <- psoptim.resid[,cells.f]
psoptim.par <- psoptim.par[,cells.f]
dnaseg <- dnaseg[dnaseg$cell_id %in% cells.f,]
dnaraw <- dnaraw[,cells.f]
dnanorm <- dnanorm[,cells.f]
dnalrr <- dnalrr[,cells.f]

cat("Saving filtered PSO/DNA sets\n")
save(file=snakemake@output[["pso_filtered"]],
     list=c("psoptim","psoptim.resid","psoptim.par","dnaseg","dnaraw","dnanorm","dnalrr"))

# Debug
# Warning! This will load the already filtered pso set, which does not work with plot_cell (by index). To plot, load debug file and rerun file from top.
if(as.logical(snakemake@params[["debug"]])) cat("^Saving Rda for debug\n"); save(list = ls(all.names = TRUE), file="debug/dna_processing.Rda")