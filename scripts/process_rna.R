suppressPackageStartupMessages(library(tidyverse))
library(Seurat)
library(cowplot)
theme_set(theme_cowplot())
library(RColorBrewer)

# RNA settings
countsfile <- snakemake@input[["counts"]]
metafile <- snakemake@input[["meta"]]
dnaclass_file <- snakemake@input[["dnaclass"]]
dnadamage_file <- snakemake@input[["dnadamage"]]
source(snakemake@params[["fxfile"]])

mincount <- as.integer(snakemake@params[["mincount"]])
meanexpr_cutoff <- as.numeric(snakemake@params[["meanexpr_cutoff"]])
disp_cutoff <- as.numeric(snakemake@params[["disp_cutoff"]])
minquant <- as.numeric(snakemake@params[["qcut"]])
genequant <- as.character(snakemake@params[["qgene"]])
perpl_rna <- as.integer(snakemake@params[["tsne_perpl"]])
print(paste("File:",countsfile))
print(paste("Mininum counts/cell:",mincount))
print(paste("Quantile cutoff:",minquant,"using",genequant))

# Meta
meta_dna_class <- read_tsv(dnaclass_file)
meta_dna_damage <- read_tsv(dnadamage_file)
meta <- read_tsv(metafile) %>% 
    left_join(meta_dna_class, by="cell_id") %>% 
    left_join(meta_dna_damage, by="cell_id")

# Load RNA
counts <- data.table::fread(countsfile, sep="\t", stringsAsFactors=F, header=T, data.table=F)
dim(counts)
rownames(counts) <- counts[[1]]
counts[[1]] <- NULL
colnames(counts) <- rna_to_cellid(colnames(counts)) # !! Remove *R suffix
singlets <- meta$cell_id[meta$sort_count==1]
include_cells <- colnames(counts) %in% singlets
ercc.idx <- grepl("^ERCC-", rownames(counts)) # Select row indices and not ERCC names 
counts.f <- counts[!grepl("^__|_",rownames(counts)) & !ercc.idx,include_cells] # Throws HTseq rows, ERCC and underscored genes (2)
counts.ercc <- counts[ercc.idx,include_cells]
ercc.frac <- colSums(counts.ercc)/(colSums(counts.f)+colSums(counts.ercc))
counts.f[1:5,1:5]
meta_seurat <- meta %>% filter(cell_id %in% colnames(counts.f)) %>% mutate(ercc_frac=ercc.frac[match(cell_id, names(ercc.frac))]) %>% as.data.frame
rownames(meta_seurat) <- meta_seurat$cell_id
colnames(meta_seurat)

dat_list <- load_seurat_rna(counts.f, mincount, minquant, genequant, meanexpr_cutoff, disp_cutoff, assay="RNA", verbose = T)
rna.countfilter <- dat_list$countfilter
rna.genefilter <- dat_list$genefilter
dat <- dat_list$seurat
dim(dat)
dat <- AddMetaData(dat, meta_seurat)
dat <- NormalizeData(dat, normalization.method = "LogNormalize", scale.factor = 1e4, verbose=F)
dat <- FindVariableFeatures(dat, selection.method="mean.var.plot", mean.cutoff=c(0.3, +Inf), dispersion.cutoff = c(0.5, Inf), verbose=F)
dat <- ScaleData(dat, verbose=F)
dat <- RunPCA(dat, npcs=50, verbose=F)
dat <- FindNeighbors(dat, dims=1:10, verbose=F)
dat <- FindClusters(dat, resolution = 1, algorithm=1, n.iter=1000, n.start=100, verbose=F)
dat <- RunTSNE(dat, n.components=2, dims=1:10,perplexity=perpl_rna, verbose=F)
dat <- RunUMAP(dat, n.components=2, dims=1:10,verbose=F)

# Cell cycle scoring
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
dat <- CellCycleScoring(dat, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

# Save classes
rna_classes <- unlist(dat$RNA_snn_res.1)
rna_classes_df <- rna_classes %>% enframe("cell_id", "rna_class") %>% 
    mutate(rna_ccscore=dat$Phase)
meta_wclass <- left_join(meta, rna_classes_df, by="cell_id")

# Colors
class_cols <- Polychrome::kelly.colors(12)[-1] # Skip white
names(class_cols) <- 0:10
meta_wclass <- meta_wclass %>% mutate_at(vars(contains("class")), as.character)
source_cols <- c("Control"="#9E9E9E", 
                 "ETO 1uM"="#FFEB3B", "ETO 3uM"="#FFC107", "ETO 10uM"="#FF9800", "ETO 30uM"="#FF5722", "ETO 100uM"="#d50000", 
                 "X-ray 1Gy"="#03A9F4", "X-ray 3Gy"="#2196F3", "X-ray 5Gy"="#3F51B5","HCT-SS2"="#EEEE00","HCT-DNTR"="#000000")
gate_cols <- c("G1"="#4CAF50", "non-G1"="#FF5722","HCT Singlets"="#CCCCCC","HCT S Phase"="#FF0000","HCT G2 Phase"="#FFF000")
phase_cols <- c("G1"="grey","S"="#FF5722","G2M"="darkblue")

# Embeddings
rna_tsne_df <- dat@reductions$tsne@cell.embeddings %>% as.data.frame %>% rownames_to_column("cell_id")
rna_umap_df <- dat@reductions$umap@cell.embeddings %>% as.data.frame %>% rownames_to_column("cell_id")

plot_tsne_rna <- function(df, expr_obj=NULL, meta_df=meta_wclass, marker, colorset=T, title="RNA", show_legend=F, show_title=T, show_axis=T){
    df <- df %>% 
        left_join(meta_df, by="cell_id")
    if(!is.null(expr_obj)){
        expr <- FetchData(object=expr_obj, vars=marker) %>% rownames_to_column("cell_id") %>% 
            mutate(cell_id=rna_to_cellid(cell_id))
        df <- df %>% 
            left_join(expr, by="cell_id")
    }
    p <- df %>% 
        ggplot(aes(tSNE_1,tSNE_2)) + 
        geom_point(color="black",size=1.8) +
        geom_point(aes(color=!!ensym(marker)), size=1.5, alpha=0.85) +
        labs(x="RNA tSNE-1", y="RNA tSNE-2", title=title, subtitle = paste0("n=",nrow(df),", perplexity ",perpl_rna,", ",marker)) +
        theme(axis.title = element_text(hjust=0), axis.ticks=element_blank(), axis.text=element_blank(), 
              legend.justification = "right", legend.position="bottom", legend.title=element_text(hjust=0, size=10), legend.box.margin = margin(-0.5, unit="cm"),
              plot.title = element_text(hjust=0), panel.border = element_rect(fill = NA, colour = "black", size=1),
              rect = element_rect(fill = "transparent", colour = "black", size = 1, linetype = "solid"), axis.line = element_blank())
    if(is.character(colorset)){
        p <- p + scale_color_manual(values=colorset, na.value="#d3d3d3") +
            theme(legend.text=element_text(size=8)) +
            guides(colour = guide_legend(override.aes = list(size=2),title.position="top",title.hjust = 0))
    } else if(colorset & length(unique(df[[marker]])) < 10) {
        p <- p + scale_color_brewer(palette="Set1", na.value="#d3d3d3") +
            theme(legend.text=element_text(size=8), legend.title=element_blank(),
                  legend.position=c(0,0),legend.justification = c(0,0)) +
            guides(color = guide_legend(override.aes = list(size=2, alpha=1, fill=NA)))
    } else if(grepl("^chr", marker) && is.numeric(df[[marker]])){
        p <- p + scale_color_gradient2(low="blue", mid="lightgrey", high="red", midpoint=0, na.value="#d3d3d3") +
            theme(legend.direction="horizontal", legend.position=c(0,0),legend.justification=c(0,0)) +
            guides(color = guide_colorbar(barwidth=0, barheight=0, default.unit = "lines", title.position = "left", label = F))
    } else if(grepl("FACS$", marker) && is.numeric(df[[marker]])){
        p <- p + scale_color_gradient(low="lightgrey", high="darkblue", na.value="#ffffff") +
            theme(legend.direction="horizontal", legend.position=c(0,0),legend.justification=c(0,0)) +
            guides(color = guide_colorbar(barwidth=0, barheight=0, default.unit = "lines", title.position = "left", label = F))
    } else if(!is.null(expr_obj) && is.numeric(expr[[2]])){
        p <- p + 
            scale_color_gradient(low="lightgrey", high="red", na.value="#d3d3d3") +
            theme(legend.direction="horizontal", legend.position=c(0,0),legend.justification=c(0,0)) +
            guides(color = guide_colorbar(barwidth=0, barheight=0, default.unit = "lines", title.position = "left", label = F))
    } else {
        n <- length(unique(df[[marker]])) + 1
        cols <- c(Polychrome::kelly.colors(n)[2:n]) # Skip white
        names(cols) <- unique(df[[marker]])
        p <- p + scale_color_manual(values=cols, na.value="#d3d3d3") +
            theme(legend.text=element_text(size=8)) +
            guides(colour = guide_legend(override.aes = list(size=1.2)))
    }
    if(show_legend=="text") {
        p <- p + annotate("text",x=-Inf,y=-Inf,label=marker,hjust=-0.1,vjust=-1,size=3) + theme(legend.position="none")
    } else if(show_legend==FALSE) {
        p <- p + theme(legend.position="none")
    }
    if(!show_title) p <- p + theme(plot.title=element_blank(), plot.subtitle=element_blank())
    if(!show_axis) p <- p + theme(axis.title=element_blank())
    return(p)
}
plot_umap_rna <- function(df, expr_obj=NULL, meta_df=meta_wclass, marker, colorset=T, title="RNA", show_legend=F, show_title=T, show_axis=T){
    df <- df %>% 
        left_join(meta_df, by="cell_id")
    if(!is.null(expr_obj)){
        expr <- FetchData(object=expr_obj, vars=marker) %>% rownames_to_column("cell_id") %>% 
            mutate(cell_id=rna_to_cellid(cell_id))
        df <- df %>% 
            left_join(expr, by="cell_id")
    }
    p <- df %>% 
        ggplot(aes(UMAP_1,UMAP_2)) + 
        geom_point(color="black",size=1.8) +
        geom_point(aes(color=!!ensym(marker)), size=1.5, alpha=0.85) +
        labs(x="RNA UMAP-1", y="RNA UMAP-2", title=title, subtitle = paste0("n=",nrow(df),", perplexity ",perpl_rna,", ",marker)) +
        theme(axis.title = element_text(hjust=0), axis.ticks=element_blank(), axis.text=element_blank(), 
              legend.justification = "right", legend.position="bottom", legend.title=element_text(hjust=0, size=10), legend.box.margin = margin(-0.5, unit="cm"),
              plot.title = element_text(hjust=0), panel.border = element_rect(fill = NA, colour = "black", size=1),
              rect = element_rect(fill = "transparent", colour = "black", size = 1, linetype = "solid"), axis.line = element_blank())
    if(is.character(colorset)){
        p <- p + scale_color_manual(values=colorset, na.value="#d3d3d3") +
            theme(legend.text=element_text(size=8)) +
            guides(colour = guide_legend(override.aes = list(size=2),title.position="top",title.hjust = 0))
    } else if(colorset & length(unique(df[[marker]])) < 10) {
        p <- p + scale_color_brewer(palette="Set1", na.value="#d3d3d3") +
            theme(legend.text=element_text(size=8), legend.title=element_blank(),
                  legend.position=c(0,0),legend.justification = c(0,0)) +
            guides(color = guide_legend(override.aes = list(size=2, alpha=1, fill=NA)))
    } else if(grepl("^chr", marker) && is.numeric(df[[marker]])){
        p <- p + scale_color_gradient2(low="blue", mid="lightgrey", high="red", midpoint=0, na.value="#d3d3d3") +
            theme(legend.direction="horizontal", legend.position=c(0,0),legend.justification=c(0,0)) +
            guides(color = guide_colorbar(barwidth=0, barheight=0, default.unit = "lines", title.position = "left", label = F))
    } else if(grepl("FACS$", marker) && is.numeric(df[[marker]])){
        p <- p + scale_color_gradient(low="lightgrey", high="darkblue", na.value="#ffffff") +
            theme(legend.direction="horizontal", legend.position=c(0,0),legend.justification=c(0,0)) +
            guides(color = guide_colorbar(barwidth=0, barheight=0, default.unit = "lines", title.position = "left", label = F))
    } else if(!is.null(expr_obj) && is.numeric(expr[[2]])){
        p <- p + 
            scale_color_gradient(low="lightgrey", high="red", na.value="#d3d3d3") +
            theme(legend.direction="horizontal", legend.position=c(0,0),legend.justification=c(0,0)) +
            guides(color = guide_colorbar(barwidth=0, barheight=0, default.unit = "lines", title.position = "left", label = F))
    } else if(is.numeric(df[[marker]])){
        p <- p + 
            scale_color_gradient(low="white", high="red4", na.value="#d3d3d3") +
            theme(legend.direction="horizontal", legend.position=c(0,0),legend.justification=c(0,0)) +
            guides(color = guide_colorbar(barwidth=0, barheight=0, default.unit = "lines", title.position = "left", label = F))
    } else {
        n <- length(unique(df[[marker]])) + 1
        cols <- c(Polychrome::kelly.colors(n)[2:n]) # Skip white
        names(cols) <- unique(df[[marker]])
        p <- p + scale_color_manual(values=cols, na.value="#d3d3d3") +
            theme(legend.text=element_text(size=8)) +
            guides(colour = guide_legend(override.aes = list(size=1.2)))
    }
    if(show_legend=="text") {
        p <- p + annotate("text",x=-Inf,y=-Inf,label=marker,hjust=-0.1,vjust=-1,size=3) + theme(legend.position="none")
    } else if(show_legend==FALSE) {
        p <- p + theme(legend.position="none")
    }
    if(!show_title) p <- p + theme(plot.title=element_blank(), plot.subtitle=element_blank())
    if(!show_axis) p <- p + theme(axis.title=element_blank())
    return(p)
}

pdf(snakemake@output[["fig_boxplotByPlate"]])
par(mfrow=c(1,2))
boxplot(split(log10(dat$nCount_RNA),dat$source),main="RNA counts",ylab="log10(counts)", las=3)
boxplot(split(dat$nFeature_RNA,dat$source),main="RNA features",ylab="n genes",las=3)
dev.off()
pdf(snakemake@output[["fig_dupsvsERCC"]])
plot(meta_seurat$dna_dups_pct, meta_seurat$ercc_frac,ylab="ERCC fraction",xlab="DNA duplicate fraction",main=paste0("n=",ncol(dat)))
dev.off()

pr <- plot_umap_rna(rna_umap_df, expr_obj=NULL, meta_df = meta_wclass, marker="rna_class", colorset=class_cols, title="", show_legend=T, show_title=F, show_axis=T)
pd <- plot_umap_rna(rna_umap_df, expr_obj=NULL, meta_df = meta_wclass, marker="dna_class", colorset=class_cols, title="", show_legend=T, show_title=F, show_axis=T)
ps <- plot_umap_rna(rna_umap_df, expr_obj=NULL, meta_df = meta_wclass, marker="source", colorset=source_cols, title="", show_legend=T, show_title=F, show_axis=T)
pg <- plot_umap_rna(rna_umap_df, expr_obj=NULL, meta_df = meta_wclass, marker="sort_gate", colorset=gate_cols, title="", show_legend=T, show_title=F, show_axis=T)
pdam <- plot_umap_rna(rna_umap_df, expr_obj=NULL, meta_df = meta_wclass, marker="dna_damage", colorset=T, title="", show_legend=T, show_title=F, show_axis=T)
pplate <- plot_umap_rna(rna_umap_df, expr_obj=NULL, meta_df = meta_wclass, marker="plate_id", colorset=T, title="", show_legend=T, show_title=F, show_axis=T)
pcc <- plot_umap_rna(rna_umap_df, expr_obj=NULL, meta_df = meta_wclass, marker="rna_ccscore", colorset=phase_cols, title="", show_legend=T, show_title=F, show_axis=T)
pgrid <- plot_grid(pr,pd,ps,pg,pdam,pcc,ncol=2,align="hv",axis="lrbt")
ggsave(snakemake@output[["fig_umapGrid"]], pgrid, width=12, height=16, useDingbats=F)

pr <- plot_tsne_rna(rna_tsne_df, expr_obj=NULL, meta_df = meta_wclass, marker="rna_class", colorset=class_cols, title="", show_legend=T, show_title=F, show_axis=T)
pd <- plot_tsne_rna(rna_tsne_df, expr_obj=NULL, meta_df = meta_wclass, marker="dna_class", colorset=class_cols, title="", show_legend=T, show_title=F, show_axis=T)
ps <- plot_tsne_rna(rna_tsne_df, expr_obj=NULL, meta_df = meta_wclass, marker="source", colorset=source_cols, title="", show_legend=T, show_title=F, show_axis=T)
pg <- plot_tsne_rna(rna_tsne_df, expr_obj=NULL, meta_df = meta_wclass, marker="sort_gate", colorset=gate_cols, title="", show_legend=T, show_title=F, show_axis=T)
pgrid <- plot_grid(pr,pd,ps,pg,ncol=2,align="hv",axis="lrbt")
ggsave(snakemake@output[["fig_tsneGrid"]], pgrid, width=12, height=12, useDingbats=F)

# Gene markers
mark <- FindAllMarkers(dat)
top5 <- mark %>% group_by(cluster) %>% top_n(5, avg_logFC)

p1 <- DoHeatmap(dat,features = top5$gene)
p2 <- DoHeatmap(dat,features = top5$gene, group.by = "source")
p3 <- DoHeatmap(dat,features = top5$gene, group.by = "sort_gate")
ggsave(snakemake@output[["fig_seuratHeatmaps"]], plot_grid(p1,p2,p3, ncol=1), width=10, height=16, useDingbats=F)

ph1<- DoHeatmap(dat,features = top5$gene)
phb <- DimPlot(dat, reduction="umap", label = T)
plot_grid(ph1,phb)
p_counts <- FeaturePlot(dat, reduction="umap","nCount_RNA")
p_genes <- FeaturePlot(dat, reduction="umap", "nFeature_RNA")
p_ercc <- FeaturePlot(dat, reduction="umap", "ercc_frac")
pgrid <- plot_grid(p_counts, p_genes, p_ercc,nrow=1)
ggsave(snakemake@output[["fig_seuratHeatmapDimplot"]], pgrid, width=10, height=5, useDingbats=F)

# Take cells in [subset] -- find markers for dna_class 0
# c0 <- FindMarkers(dat, ident.1 = "0", group.by="dna_class")
# VlnPlot(dat, features = c("S100A6","MKI67"), group.by = "source", cols = source_cols)
# VlnPlot(dat, features = c("TOP2A","MKI67","PLK1"), group.by = "sort_gate", cols = gate_cols)
# VlnPlot(dat, features = c("S100A6","MKI67"), group.by = "dna_class", cols = class_cols)

p1 <- plot_umap_rna(rna_umap_df, expr_obj=dat, meta_df = meta_wclass, marker="CDC20", colorset=T, title="", show_legend=T, show_title=F, show_axis=F)
p2 <- plot_umap_rna(rna_umap_df, expr_obj=dat, meta_df = meta_wclass, marker="TOP2A", colorset=T, title="", show_legend=T, show_title=F, show_axis=F)
p3 <- plot_umap_rna(rna_umap_df, expr_obj=dat, meta_df = meta_wclass, marker="MDM2", colorset=T, title="", show_legend=T, show_title=F, show_axis=F)
p4 <- plot_umap_rna(rna_umap_df, expr_obj=dat, meta_df = meta_wclass, marker="MKI67", colorset=T, title="", show_legend=T, show_title=F, show_axis=F)
p5 <- plot_umap_rna(rna_umap_df, expr_obj=dat, meta_df = meta_wclass, marker="CDKN1A", colorset=T, title="", show_legend=T, show_title=F, show_axis=F)
p6 <- plot_umap_rna(rna_umap_df, expr_obj=dat, meta_df = meta_wclass, marker="SAT1", colorset=T, title="", show_legend=T, show_title=F, show_axis=F)
p7 <- plot_umap_rna(rna_umap_df, expr_obj=dat, meta_df = meta_wclass, marker="FOS", colorset=T, title="", show_legend=T, show_title=F, show_axis=F)
p8 <- plot_umap_rna(rna_umap_df, expr_obj=dat, meta_df = meta_wclass, marker="HIST1H1E", colorset=T, title="", show_legend=T, show_title=F, show_axis=F)
p9 <- plot_umap_rna(rna_umap_df, expr_obj=dat, meta_df = meta_wclass, marker="HIST1H4C", colorset=T, title="", show_legend=T, show_title=F, show_axis=F)
pgrid <- plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9)
ggsave(snakemake@output[["fig_markerGrid"]], pgrid, width=10, height=10, useDingbats=F)

save(list=c("dat","meta","meta_wclass","rna.countfilter","rna.genefilter"),file=snakemake@output[["rda"]])
if(as.logical(snakemake@params[["debug"]])) cat("^Saving Rda for debug\n"); save(list = ls(all.names = TRUE), file="debug/rna_processing.Rda")