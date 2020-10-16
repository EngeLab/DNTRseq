# Prepare ploidy data. Run on server
library(tidyverse)
library(cowplot)
library(scales)
options(scipen=999)

fxfile <- snakemake@params[["fxfile"]]
files <- snakemake@input[["files"]]
metafile <- snakemake@input[["meta"]]
outfile <- snakemake@output[["rda"]]

# Debug
fxfile <- "scripts/functions.R"
f <- readLines("meta/ploidy_cells_list.tsv")
files <- paste0("/wrk/data/human_genomic/aligned/",f,"/",f,".fragtrimcov.tsv")
metafile <- "meta/metadata.tsv"
outfile <- "out/ploidy.Rda"
rna <- "out/rna.Rda"

source(fxfile)
meta <- read_tsv(metafile)
print(load(rna))

dat <- map_dfr(files, read_tsv, col_names=c("path","avg_len","total_frags","total_bp","chr","cov","bp_covered","bp_space","cov_frac"),col_types="cniicinnn") %>% 
    mutate(file=basename(path)) %>%
    separate(file, c("dna_library_id",NA, NA), "\\.", remove=T) %>% # Ex: HCA00102D_A10.frag_trim_dedup_655000.cov
    mutate(n_frags=sub(".*_([0-9]+).cov","\\1",path))
# OR 
#dat <- read_rds("out/ploidy_200528.Rds")
#dat <- read_rds("out/ploidy_lnscale_200531.Rds")

dat.m <- dat %>% left_join(meta,by="dna_library_id")
# Add other DNA QC measures?

dat_ext <- dat.m %>% 
    filter(sort_count==1, cov <= 5, dna_read_pairs>80e3, dna_dups_pct < 0.15) %>% 
    group_by(source, path) %>% 
    mutate(cov_bin=cov) %>% 
    #mutate(cov_bin=ifelse(cov < 5, cov, 5)) %>% # Bin any coverage at >=5x
    ungroup() %>%
    mutate(n_frags=as.numeric(n_frags),
           mapped_bp=n_frags*avg_len,
           mapped_frac=mapped_bp/bp_space) %>% 
    group_by(source, path, cov_bin) %>% # Calculate new coverages at 0,1,2,3,4,>=5x coverage
    mutate(bp_covered_bin=sum(bp_covered),
           cov_frac_bin=bp_covered_bin/bp_space) %>% 
    ungroup %>% 
    mutate(subset=factor(paste(source, sort_gate)))

# Colors
colors <- Polychrome::kelly.colors(length(levels(dat_ext$subset))+1)[-1] # Skip white
names(colors) <- levels(dat_ext$subset)
gate_cols <- c("G1"="#03a9f4","non-G1"="#ffa000")
phase_cols <- c("G1"="grey","S"="#FF5722","G2M"="darkblue")
source_cols <- c("Control"="#9E9E9E", 
                 "ETO 1uM"="#FFEB3B", "ETO 3uM"="#FFC107", "ETO 10uM"="#FF9800", "ETO 30uM"="#FF5722", "ETO 100uM"="#d50000", 
                 "X-ray 1Gy"="#03A9F4", "X-ray 3Gy"="#2196F3", "X-ray 5Gy"="#3F51B5")

# Fit model
df_fit_model <- dat_ext %>%
    mutate(mapped_frac_sq=mapped_frac^2) %>% 
    group_by(subset, dna_library_id, cov_bin) %>% # Create 1 model per cell and coverage bin (eg cell A at 1x, cell A at 2x, cell B at 1x, cell B at 2x etc)
    nest() %>% 
    mutate(model=map(data, ~lm(cov_frac_bin ~ mapped_frac + mapped_frac_sq, data=.x)))
df_fit <- df_fit_model %>% 
    mutate(fit=map(model, ~tibble(x=seq(0,0.2,0.001),y=predict(.x, list(mapped_frac=seq(0,0.2,0.001), mapped_frac_sq=seq(0,0.2,0.001)^2))))) %>% 
    unnest(fit)

# SVD fit 0,1,2x coverage fraction
# df_svd <- df_fit %>% filter(x==0.2) %>% spread(cov_bin, y)
# #df_svd$filter <- !rowSums(is.na(select(df_svd, "0","1","2","3")))>0
# df_svd$filter <- !rowSums(is.na(select(df_svd, "1","2")))>0
# df_svd_filt$gate <- ifelse(grepl("non-G1",df_svd_filt$subset),"non-G1","G1")
# df_svd_filt$source <- lapply(strsplit(as.character(df_svd_filt$subset), split=" "), function(x){ ifelse(length(x)>2,paste0(x[1:2],sep=""),x[1]) }) %>% unlist
# svd_res <- svd(as.matrix(df_svd_filt[,colnames(df_svd_filt) %in% 1:2]))

# res_df <- as.data.frame(svd_res$u[,1]) # Select first SVD component
# colnames(res_df) <- "svd"
# res_df2 <- res_df %>% mutate(dna_library_id=df_svd_filt$dna_library_id)

save(list = c("dat_ext", "df_fit"), file=outfile)

# SVD plot
p0 <- res_df2 %>% 
    left_join(meta, by="dna_library_id") %>% 
    mutate(subset=factor(paste(source, sort_gate))) %>% 
    #gather(coverage, svd_1, -dna_library_id, -sort_gate, -source) %>% 
    ggplot(aes(subset,svd,color=sort_gate)) + 
    ggbeeswarm::geom_quasirandom(width = 0.15,pch=19,alpha=.8, size=0.6) +
    scale_color_manual(values=gate_colors) +
    facet_wrap(~source,scales="free",nrow=1) +
    theme(strip.text=element_text(hjust=0,size=12),
          panel.border = element_rect(fill = NA, colour = "black", size=1),
          rect = element_blank(), axis.line=element_blank(),
          legend.position="right", legend.justification = "top", legend.title = element_blank(),
          legend.direction = "vertical", legend.margin=margin(0,0,0,0,unit="pt"),
          title=element_text(size=10), axis.text=element_text(size=10),
          axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
          plot.margin=margin(0,0,0,0,unit="pt"))

# Orig plot
# Remove (not delete) outliers top/bot 0.5%
filt <- df_fit %>% filter(x==0.2 & cov_bin==2)
y <- filt$y
y.m <- median(y)
y.sd <- sqrt(sum((y[y > y.m] - y.m)^2)/(sum(y  > y.m)-1))
bot <- qnorm(p=0.05, mean=y.m, sd=y.sd)
top <- qnorm(p=0.995, mean=y.m, sd=y.sd)

library(ggridges)

# Ridgeplot of distributions
ridge1 <- filt %>% 
    left_join(meta_wclass, by="dna_library_id") %>% 
    mutate(source=factor(source, levels=c("Control",
                                          "ETO 1uM", "ETO 3uM", "ETO 10uM", "ETO 30uM", "ETO 100uM",
                                          "X-ray 1Gy", "X-ray 3Gy", "X-ray 5Gy"))) %>% 
    filter(source %in% c("Control") | grepl("X-ray",source)) %>%
    #filter(source %in% c("Control")) %>% 
    ggplot(aes(y,rna_ccscore,fill=sort_gate)) + 
    geom_density_ridges(aes(y=rna_ccscore),alpha=0.5) +
    scale_fill_manual(values=gate_cols) +
    facet_wrap(~source) +
    theme_dntr() +
    labs(x="Genome fraction covered at >2x (modelled 0.2x coverage)",
         y="Seurat cell cycle score")

ridge2 <- filt %>% 
    left_join(meta_wclass, by="dna_library_id") %>% 
    mutate(source=factor(source, levels=c("Control",
                                          "ETO 1uM", "ETO 3uM", "ETO 10uM", "ETO 30uM", "ETO 100uM",
                                          "X-ray 1Gy", "X-ray 3Gy", "X-ray 5Gy"))) %>% 
    filter(sort_gate=="G1") %>%
    ggplot(aes(y,fill=source)) + 
    geom_density_ridges(aes(y=source)) +
    scale_fill_manual(values=source_cols) +
    theme_dntr() +
    labs(x="Genome fraction covered at >2x (modelled 0.2x coverage)",
         y="Seurat cell cycle score")

meta_wclass %>% 
    mutate(source=factor(source, levels=c("Control",
                                          "ETO 1uM", "ETO 3uM", "ETO 10uM", "ETO 30uM", "ETO 100uM",
                                          "X-ray 1Gy", "X-ray 3Gy", "X-ray 5Gy"))) %>% 
    ggplot(aes(dna_damage,fill=source)) +
    geom_density_ridges(aes(y=source)) +
    scale_fill_manual(values=source_cols) +
    theme_dntr() +
    labs(x="DNA damage fraction")
    
    

p1 <- filt %>% 
    left_join(meta, by="dna_library_id") %>%
    mutate(source=factor(source, levels=c("Control",
                                          "ETO 1uM", "ETO 3uM", "ETO 10uM", "ETO 30uM", "ETO 100uM",
                                          "X-ray 1Gy", "X-ray 3Gy", "X-ray 5Gy"))) %>% 
    ggplot(aes(subset,y,color=sort_gate)) + 
    ggbeeswarm::geom_quasirandom(width = 0.15,pch=19,alpha=.8, size=0.6) +
    scale_color_manual(values=gate_cols) +
    facet_wrap(~source,scales="free_x",nrow=1) +
    coord_cartesian(ylim=c(bot,top)) +
    theme(strip.text=element_text(hjust=0,size=12),
          panel.border = element_rect(fill = NA, colour = "black", size=1),
          rect = element_blank(), axis.line=element_blank(),
          legend.position="right", legend.justification = "top",
          legend.direction = "vertical", legend.margin=margin(0,0,0,0,unit="pt"),
          axis.text=element_text(size=10),
          axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
          plot.margin=margin(0,0,0,0,unit="pt")) +
    labs(subitle="Covered frac of genome at >=2x depth",x="",y="Covered fraction of genome")

# Add plots grouped by Seurat CC score. FACS Hoechst etc
# Load extra data
print(load(paste0("out/dnaobj-psoptim-500000.Rda")))

m.sort <- meta[match(colnames(psoptim), as.character(meta$cell_id)),]
m.sort$pso_par <- unlist(psoptim.par[1,]) # Scaling factor from PSO
#m.sort$pso_par_bins <- cut_interval(m.sort$pso_par, 3)
m.sort$pso_par_bins <- cut(m.sort$pso_par, breaks = c(1.5,1.8,2.2,3))
pso.med <- apply(psoptim[,m.sort$source=='Control'], 1, median)
psoptim.rel <- apply(psoptim, 2, function(x) {x-pso.med})
m.sort$dna_damage_frc <- colSums(psoptim.rel != 0)/nrow(psoptim.rel)

filt %>% 
    left_join(m.sort, by="dna_library_id") %>% 
    ggplot(aes(reorder(seq_along(cell_id),-y),y,color=factor(pso_par_bins))) + geom_point()

p2 <- filt %>% 
    left_join(m.sort, by="dna_library_id") %>% 
    ggplot(aes(subset,y,color=factor(pso_par_bins))) + 
    ggbeeswarm::geom_quasirandom(width = 0.15,pch=19,alpha=.8, size=0.6) +
    #scale_color_manual(values=gate_colors) +
    scale_color_brewer(palette="Set1") +
    facet_wrap(~source,scales="free_x",nrow=1) +
    #coord_cartesian(ylim=c(bot,top)) +
    theme(strip.text=element_text(hjust=0,size=12),
          panel.border = element_rect(fill = NA, colour = "black", size=1),
          rect = element_blank(), axis.line=element_blank(),
          legend.position="right", legend.justification = "top",
          legend.direction = "vertical", legend.margin=margin(0,0,0,0,unit="pt"),
          axis.text=element_text(size=10),
          axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
          plot.margin=margin(0,0,0,0,unit="pt")) +
    labs(subitle="Covered frac of genome at >=2x depth",x="",y="Covered fraction of genome") +
    guides(color=guide_legend(override.aes = list(size=2,alpha=1)))
    
p3 <- res_df2 %>% 
    left_join(m.sort, by="dna_library_id") %>% 
    ggplot(aes(subset,y,color=factor(pso_par_bins))) + 
    ggbeeswarm::geom_quasirandom(width = 0.15,pch=19,alpha=.8, size=0.6) +
    #scale_color_manual(values=gate_colors) +
    scale_color_brewer(palette="Set1") +
    facet_wrap(~source,scales="free_x",nrow=1) +
    #coord_cartesian(ylim=c(bot,top)) +
    theme(strip.text=element_text(hjust=0,size=12),
          panel.border = element_rect(fill = NA, colour = "black", size=1),
          rect = element_blank(), axis.line=element_blank(),
          legend.position="right", legend.justification = "top",
          legend.direction = "vertical", legend.margin=margin(0,0,0,0,unit="pt"),
          axis.text=element_text(size=10),
          axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
          plot.margin=margin(0,0,0,0,unit="pt")) +
    labs(subitle="Covered frac of genome at >=2x depth",x="",y="Covered fraction of genome")
    

# > d[rowSums(is.na(as.matrix(d[,5:7])))==1,]
# # A tibble: 6 x 9
# subset      dna_library_id     x   `0`   `1`    `2`   `3`   `4`   `5`
# <fct>       <chr>          <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl>
#     1 ETO 10uM G1 VZA01707D_B11    0.2 0.769 0.214 0.0179    NA    NA    NA
# 2 ETO 10uM G1 VZA01707D_E11    0.2 0.872 0.114 0.0136    NA    NA    NA
# 3 ETO 10uM G1 VZA01707D_F05    0.2 0.845 0.143 0.0152    NA    NA    NA
# 4 ETO 1uM G1  VZA01706D_E05    0.2 0.802 0.177 0.0209    NA    NA    NA
# 5 ETO 1uM G1  VZA01706D_H09    0.2 0.825 0.159 0.0159    NA    NA    NA
# 6 ETO 1uM G1  VZA01706D_N05    0.2 0.875 0.110 0.0141    NA    NA    NA

# Hoechst data
facs_files <- list.files("meta/facs_hoechst/", full.names = T)
names(facs_files) <- sub(".*vasilios_(.*)_Live singlets.csv","\\1", facs_files)
facs_df <- map_dfr(facs_files, read.csv, na.strings=c("", NA), stringsAsFactors=F, .id="id") %>% 
    mutate(hoechst=case_when(!is.na(Hoechst.33342.A.Compensated) ~ Hoechst.33342.A.Compensated,
                             !is.na(Hoechst.33342.A) ~ Hoechst.33342.A.Compensated),
           pi=case_when(!is.na(PI.A.Compensated) ~ PI.A.Compensated,
                        !is.na(PI.A) ~ PI.A),
           plate_id=sub(".*(VZA[0-9]{5}).*","\\1",id),
           well_id=ifelse(!is.na(Index), paste0(str_sub(Index,1,1), str_pad(str_sub(Index, 2,-1),2,"left",0)), NA), # Pad well numbers (H1 --> H01)
           cell_id = ifelse(!is.na(Index), paste0(plate_id, ".", well_id), NA))

facs_df %>% 
    ggplot(aes(hoechst, plate_id)) + geom_density_ridges() + theme_dntr() +
    coord_cartesian(xlim=c(0.2e6,1e6)) 

meta_facs <- facs_df %>% 
    left_join(meta_wclass, by="cell_id")

# PLOT HOECHST VS ploidy score 0.2X


### Plotting
plot_lines <- function(dat){
    dat <- dat %>% mutate(cov_bin=paste0(cov_bin,"x"))
    dat1 <- dat %>% filter(grepl("Non",subset))
    dat2 <- dat %>% filter(grepl("doublet",subset))
    dat3 <- dat %>% filter(grepl("S",subset))
    dat4 <- dat %>% filter(grepl("G2",subset))
    dat4 %>% 
        ggplot(aes(x=mapped_frac, y=cov_frac_bin)) +
        geom_line(aes(group=dna_library_id, color=subset), alpha=0.5) +
        geom_line(dat=dat1, aes(group=dna_library_id, color=subset), alpha=0.5) +
        geom_line(dat=dat2, aes(group=dna_library_id, color=subset), alpha=0.5) +
        geom_line(dat=dat3, aes(group=dna_library_id, color=subset), alpha=0.5) +
        scale_y_continuous(breaks=pretty_breaks(4)) + 
        facet_wrap(~cov_bin, scales="free_y",nrow=1) +
        theme(strip.text=element_text(hjust=0,size=12),
              panel.border = element_rect(fill = NA, colour = "black", size=1),
              rect = element_blank(), axis.line=element_blank(),
              legend.position="right", legend.justification = "top", legend.title = element_blank(),
              legend.direction = "vertical", legend.margin=margin(0,0,0,0,unit="pt"),
              title=element_text(hjust=0,size=10),
              plot.margin=margin(0,0,0,0,unit="pt")) +
        #labs(x=paste0("Mapped bp\n(as fraction of genome space; ",round(unique(dat$bp_space)/1e9,2),"Gbp)"), y="Fraction covered") +
        labs(x=paste0("Mapped fraction"), y="Fraction covered") +
        #geom_label_repel(data=labels_df, aes(label=sub("HCA00102D_","",cell)), nudge_x = 0.05) +
        scale_color_manual(values=colors) +
        guides(colour = guide_legend(override.aes = list(size=2, alpha=1, fill=NA)))
}
plot_lines2 <- function(dat, line="dashed", ao=0.5, am=0.5,rows=2,show_avg_frac=T,legend=T){
    dat <- dat %>% mutate(cov_bin=paste0(cov_bin,"x"))
    dat1 <- dat %>% filter(grepl("S",subset))
    dat2 <- dat %>% filter(grepl("Non",subset))
    dat3 <- dat %>% filter(grepl("doublet",subset))
    dat4 <- dat %>% filter(grepl("G2",subset))
    # Calc median mapped fraction
    avg_mapped_frac <- distinct(dat, dna_library_id, .keep_all=T) %>% pull(max_mapped) %>% mean
    p <- ggplot() +
        geom_line(data=dat1[dat1$obs=="real", ], aes(x,y,group=dna_library_id, color=subset), alpha=ao) +
        geom_line(data=dat1[dat1$obs=="model",], aes(x,y,group=dna_library_id, color=subset), linetype=line, alpha=am) +
        geom_line(data=dat2[dat2$obs=="real", ], aes(x,y,group=dna_library_id, color=subset), alpha=ao) +
        geom_line(data=dat2[dat2$obs=="model",], aes(x,y,group=dna_library_id, color=subset), linetype=line, alpha=am) +
        geom_line(data=dat3[dat3$obs=="real", ], aes(x,y,group=dna_library_id, color=subset), alpha=ao) +
        geom_line(data=dat3[dat3$obs=="model",], aes(x,y,group=dna_library_id, color=subset), linetype=line, alpha=am) +
        geom_line(data=dat4[dat4$obs=="real", ], aes(x,y,group=dna_library_id, color=subset), alpha=ao) +
        geom_line(data=dat4[dat4$obs=="model",], aes(x,y,group=dna_library_id, color=subset), linetype=line, alpha=am) +
        scale_y_continuous(breaks=pretty_breaks(3)) + 
        facet_wrap(~cov_bin, scales="free_y",nrow=rows) +
        theme(strip.text=element_text(hjust=0,size=12),
              panel.border = element_rect(fill = NA, colour = "black", size=1),
              rect = element_blank(), axis.line=element_blank(),
              legend.position="right", legend.justification = "top", legend.title = element_blank(),
              legend.direction = "vertical", legend.margin=margin(0,0,0,0,unit="pt"),
              title=element_text(size=10), axis.text=element_text(size=10),
              plot.margin=margin(0,0,0,0,unit="pt")) +
        labs(x=paste0("Mapped reads (as fraction of genome space)"), y="Fraction covered") +
        #geom_label_repel(data=labels_df, aes(label=sub("HCA00102D_","",cell)), nudge_x = 0.05) +
        scale_color_manual(values=colors) +
        guides(colour = guide_legend(override.aes = list(size=2, alpha=1, fill=NA)))
    if(show_avg_frac) p <- p + geom_vline(xintercept=avg_mapped_frac, color="grey70",linetype="dotted",size=0.5)
    if(!legend) p <- p + theme(legend.position="none")
    return(p)
}
plot_boxes <- function(dat){
    dat %>% 
        mutate(cov_bin=paste0(cov_bin,"x")) %>% 
        ggplot(aes(subset, y, color=subset)) +
        geom_boxplot(position="dodge", outlier.shape = NA) +
        #geom_point(position=position_jitterdodge(jitter.width = 0.2),pch=21) +
        ggbeeswarm::geom_quasirandom(width = 0.25, alpha=0.7) +
        facet_wrap(~cov_bin,scale="free_y",nrow=1) +
        scale_color_manual(values=colors) +
        scale_y_continuous(breaks=pretty_breaks(4)) + 
        theme(strip.text=element_text(hjust=0,size=12),
              panel.border = element_rect(fill = NA, colour = "black", size=1),
              rect = element_blank(), axis.line=element_blank(),
              legend.position="right", legend.justification = "top", legend.title = element_blank(),
              legend.direction = "vertical", legend.margin=margin(0,0,0,0,unit="pt"),
              title=element_text(size=10),axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
              plot.margin=margin(0,0,0,0,unit="pt")) +
        labs(x="",y="Cov frac at 0.2x")
}
plot_crosses <- function(dat, ptsize=0.6, rows=2){
    dat <- dat %>% mutate(cov_bin=paste0(cov_bin,"x"))
    median <- dat %>% group_by(subset, cov_bin) %>% summarize(median=median(y))
    error <- dat %>% group_by(subset, cov_bin) %>% summarize(median=median(y), sd=sd(y), mad=mad(y)) %>% 
        #mutate(sd_low=median-sd, sd_high=median+sd)
        mutate(sd_low=median-mad, sd_high=median+mad)
    dat %>% 
        ggplot(aes(subset, y, color=subset)) +
        ggbeeswarm::geom_quasirandom(width = 0.15,pch=19,alpha=.8, size=ptsize) +
        geom_errorbar(data=error, aes(ymin=sd_low, ymax=sd_high, x=subset), inherit.aes=F, color="grey30", size=0.7, width=0.01, alpha=.5) +
        geom_point(data=median, aes(y=median, x=subset), shape=95, color="grey30", size=5, inherit.aes = F, alpha=.7) +
        facet_wrap(~cov_bin,scale="free_y",nrow=rows) +
        scale_color_manual(values=colors) +
        scale_y_continuous(breaks=pretty_breaks(3)) + 
        theme(strip.text=element_text(hjust=0,size=12),
              panel.border = element_rect(fill = NA, colour = "black", size=1),
              rect = element_blank(), axis.line=element_blank(),
              legend.position="right", legend.justification = "top", legend.title = element_blank(),
              legend.direction = "vertical", legend.margin=margin(0,0,0,0,unit="pt"),
              title=element_text(size=10), axis.text=element_text(size=10),
              axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
              plot.margin=margin(0,0,0,0,unit="pt")) +
        labs(x="",y="Cov frac at 0.2x") +
        guides(colour = guide_legend(override.aes = list(size=3)))
}

# p1 <- dat_ext %>% filter(patient_id=="HCT116") %>% plot_lines
# p2 <- df_fit %>% 
#     filter(grepl("HCT116",subset) & x==0.2 & cov_bin %in% c(0:2)) %>% 
#     plot_crosses(rows=1) + theme(axis.text.x=element_blank())

pX <- df_fit %>% 
    filter(grepl("X-ray",subset) & x==0.2 & cov_bin %in% c(1:3) & y>0) %>% 
    plot_crosses(rows=1) + theme(axis.text.x=element_blank())
pE <- df_fit %>% 
    filter(grepl("ETO",subset) & x==0.2 & cov_bin %in% c(1:3) & y>0) %>% 
    plot_crosses(rows=1) + theme(axis.text.x=element_blank())
pC <- df_fit %>% 
    filter(grepl("Control",subset) & x==0.2 & cov_bin %in% c(1:3) & y>0) %>% 
    plot_crosses(rows=1) + theme(axis.text.x=element_blank())

plot_grid(pX,pE,pC,ncol=1,align="hv")