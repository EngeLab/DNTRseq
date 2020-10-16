# DNTR-seq processing workflow (Vasilios Zachariadis 2020)
# Run with: 
# conda env create -f dntr.yml
# conda activate dntr
# snakemake

import pandas as pd
import os
import glob
#------------------------------------------ Set variables
configfile: "config.yaml"
cells_dna = pd.read_table(config["cells_dna"], sep="\t", header=0).set_index("cell", drop=False)
cells_rna = pd.read_table(config["cells_rna"], sep="\t", header=0).set_index("cell", drop=False)

dnaout = config["path"]["dnaout"]
rnaout = config["path"]["rnaout"]
binsizes = config["dna"]["binsizes"]
#------------------------------------------ Main
rule all: 
	input: 
		# DNA
		expand(os.path.join("out","bincounts-{binsize}.tsv"), binsize=binsizes),
		expand(os.path.join("out","dnaobj-{binsize}.Rds"), binsize=binsizes),
		expand(os.path.join("out","dnaobj-psoptim-{binsize}.Rda"), binsize=binsizes),
		expand(os.path.join("out","dna-filtered-{binsize}.Rda"), binsize=binsizes),
		
		# RNA
		"out/rna_counts.tsv.gz",
		"out/rna.Rda",

		# Metadata
		"meta/metadata.tsv",
		expand(os.path.join("meta","metadata_dnafilters_{binsize}.tsv"),binsize=binsizes),
		
#------------------------------------------ DNA setup
# Create necessary genome files and variable bins using Ginkgo's setup
# Scripts are here: https://github.com/robertaboukhalil/ginkgo/tree/master/genomes/scripts
# scripts/build_genome_files_bybin.sh [resources/genome_files] [resources/genome_files/ginkgo_scripts] [hg38] [hg38-iGenome-ERCC.fa] [hg38_RefSeq.bed] [read_length] [binsize]
# Pre-made files (bins, genes, boundaries etc) are provided

#------------------------------------------ DNA
def get_dna_fq(wildcards):
	'''
	Returns paths coded into cells_dna.tsv under "fq1" and "fq2" for the respective cell
	Uses "cell" key from wildcards to access "cells" table, and pick out "fq1" + "fq2"
	'''
	return sorted(cells_dna.loc[wildcards.cell, ['fq1', 'fq2']])

rule fastq_screen_dna: 
	'''
	Quick sample and map against multiple genomes as QC. 
	Using untrimmed fastqs for now.
	'''
	input: get_dna_fq
	output: 
		dir=directory(os.path.join(dnaout,"{cell}/fastq_screen"))
	threads: 2
	params: 
		conf_file=config["ref"]["fastq_screen_config"]
	shell:
		'''
		fastq_screen \
		--threads {threads} \
		--conf {params.conf_file} \
		--aligner 'bwa' \
		--outdir {output.dir} {input[0]}
		'''

rule cutadapt_pe:
	'''
	Trim adapter sequences and poor quality bases
	'''
	input: get_dna_fq
	output: 
		temp(os.path.join(dnaout,"{cell}/trimmed/{cell}_val_1.fq.gz")),
		temp(os.path.join(dnaout,"{cell}/trimmed/{cell}_val_2.fq.gz"))
	shell:
		'''
		outdir=$(dirname {output[0]})
		mkdir -vp $outdir
		trim_galore --paired --nextera -q 20 -e 0.1 --basename {wildcards.cell} -o $outdir {input}
		'''

rule align_bwa: 
	'''
	Align paired DNA reads with bwa
	'''
	input: 
		fq1=os.path.join(dnaout,"{cell}/trimmed/{cell}_val_1.fq.gz"),
		fq2=os.path.join(dnaout,"{cell}/trimmed/{cell}_val_2.fq.gz")
	params:
		rg="'@RG\\tID:{cell}\\tPL:ILLUMINA\\tPU:{cell}\\tLB:{cell}\\tSM:{cell}'",
		ref=config["ref"]["fasta"]
	output: 
		temp(os.path.join(dnaout,"{cell}/{cell}.sorted.bam"))
	threads: 4
	shell: 
		'''
		bwa mem -t {threads} -M -k 25 -R {params.rg} {params.ref} {input.fq1} {input.fq2} | \
		samtools sort -@{threads} -O BAM -o {output} -
		'''

rule dedup:
	input: "{cell}.sorted.bam"
	output: 
		bam="{cell}.dedup.bam",
		qc="{cell}-picard-mark_duplicates.txt"
	params:
		tmp=config["path"]["temp"]
	threads: 2
	shell: 
		'''
		picard MarkDuplicates \
		MAX_RECORDS_IN_RAM=5000000 \
		TMP_DIR={params.tmp} \
		INPUT={input} \
		OUTPUT={output.bam} \
		REMOVE_DUPLICATES=true \
		VALIDATION_STRINGENCY=LENIENT \
		METRICS_FILE={output.qc}
		'''
		
rule index_bam:
	input: "{cell}.dedup.bam"
	output: "{cell}.dedup.bam.bai"
	shell:
		"samtools index {input}"

rule create_bed:
	input: "{cell}.dedup.bam"
	output: "{cell}.bed.gz"
	params: 
		blacklist=config["ref"]["blacklist"],
		min_mapq=config["dna"]["min_mapq"]
	shell:
		'''
		bedtools bamtobed -i {input} | \
		awk -F'\t' -v mapq={params.min_mapq} '$1~/^chr([0-9]+|[X,Y])$/ && $5>=mapq{{print}}' | \
		bedtools intersect -a stdin -b {params.blacklist} -v | gzip -c > {output}
		'''

rule bincount:
	'''
	Bincounting by bedtools intersect, over variable sized bins (bed-file)
	'''
	input: "{cell}.bed.gz"
	output: "{cell}-bincounts_{binsize}.tsv"
	params: 
		bedbin=os.path.join(config["path"]["dna_bin_resources"],"variable_{binsize}_37_bwa.bed"),
		cellid=lambda wildcards: os.path.basename(wildcards.cell) # Retrieve cellid from path for header
	shell: 
		'''
		echo {params.cellid} > {output} # Set header with sample name
		bedtools intersect -sorted -nonamecheck -F 0.5 -c -a {params.bedbin} -b {input} | cut -f4 >> {output}
		'''

rule concat_bincounts:
	input: lambda wildcards: expand(os.path.join(dnaout,"{cell}/{cell}-bincounts_{binsize}.tsv"), cell=cells_dna.index, binsize=wildcards.binsize)
	output: "out/bincounts-{binsize}.tsv"
	script: "scripts/collect_bincounts.R"

rule segment_dna:
	'''
	Segment DNA read counts by variable bins using Bioconductor/DNAcopy (circular binary segmentation)
	'''
	input: "out/bincounts-{binsize}.tsv"
	output: "out/dnaobj-{binsize}.Rds"
	threads: 2
	params:
		bin_file=config["path"]["dna_bin_resources"] + "/variable_{binsize}_" + str(config["dna"]["bin_readlength"]) + "_bwa",
		gc_file=config["path"]["dna_bin_resources"] + "/GC_variable_{binsize}_" + str(config["dna"]["bin_readlength"]) + "_bwa",
		bnd_file=config["path"]["dna_bin_resources"] + "/bounds_variable_{binsize}_" + str(config["dna"]["bin_readlength"]) + "_bwa",
		readlength=config["dna"]["bin_readlength"],
		min_count=config["dna"]["min_count"],
		cbs_minwidth=config["dna"]["cbs_minwidth"],
		cbs_alpha=config["dna"]["cbs_alpha"],
		cbs_undosplits=config["dna"]["cbs_undosplits"],
		cbs_sdprune=config["dna"]["cbs_sdprune"],
		rm_outliers=config["dna"]["cbs_remove_outlier_bins"]
	script: 
		"scripts/process_dna_parallel.R"

rule integer_dna:
	'''
	Find copy number integers by PSO
	'''
	input: "out/dnaobj-{binsize}.Rds"
	output: "out/dnaobj-psoptim-{binsize}.Rda"
	threads: 7
	script: "scripts/process_pso.R"

rule process_dna_qc: 
	'''
	DNA filtering and QC plots
	'''
	input: 
		pso="out/dnaobj-psoptim-{binsize}.Rda",
		meta="meta/metadata.tsv"
	output: 
		p_loss="out/p_loss_{binsize}.tsv",
		p_gain="out/p_gain_{binsize}.tsv",
		pso_relative_ctrl="out/pso_relative_control_{binsize}.Rda",
		damage="out/damage_burden_{binsize}.tsv",
		fig_top20="fig/dna_top20byDamage_{binsize}.pdf",
		fig_bot20="fig/dna_bot20byDamage_{binsize}.pdf",
		fig_crossplot="fig/dna_qc_crossplot_{binsize}.pdf",
		fig_badResiduals="fig/dna_bot25byResiduals_{binsize}.pdf",
		fig_filterScatter="fig/dna_qc_filterScatter_{binsize}.pdf",
		fig_damageBoxplot="fig/dna_qc_damageBySource_{binsize}.pdf",
		figext_mappedVSdup="fig/dnaQC_mapped_vs_dup_passQC_{binsize}.pdf",
		metadata_filtered="meta/metadata_dnafilters_{binsize}.tsv",
		pso_filtered="out/dna-filtered-{binsize}.Rda"
	wildcard_constraints: 
		binsize=500000
	params: 
		bin_file_bed=config["path"]["dna_bin_resources"] + "/variable_{binsize}_" + str(config["dna"]["bin_readlength"]) + "_bwa.bed",
		r_global=config["path"]["r_global"],
		min_count=config["dna"]["min_count"],
		max_dup_frac=config["dna"]["max_dup_frac"],
		debug=1
	script: "scripts/process_dna_qc.R"
		
rule ploidy_create_bed:
	'''
	Create BED-files with 5bp trimming for fragment overlap analysis
	''' 
	input: 
		bam="{cell}.dedup.bam"
	output: "{cell}.frag_trim_dedup.bed"
	params:
		trim_bp=5,
		chr_auto=config["ref"]["chr_autosomal"],
		min_mapq=config["dna"]["min_mapq"],
		max_insert=config["dna"]["max_insert"]
	threads: 1
	shell: 
		'''
		samtools sort -n -m 2G {input.bam} | samtools view -L {params.chr_auto} -bf 0x2 - | bedtools bamtobed -bedpe -i stdin | \
		awk -F'\t' -v mapq={params.min_mapq} -v ins={params.max_insert} '$1==$4 && ($6-$2)<ins && ($6-$2)>30 && $8>=mapq{{OFS="\t"; print $1,$2,$6,$7,$8}}' | \
		awk -F'\t' -v bp={params.trim_bp} '{{if(($2+bp)>= ($3-bp)){{ next }};OFS="\t";print $1,($2+bp),($3-bp),$4,$5}}' | \
		awk -F'\t' '!_[$2]++' | \
		sort -S 2G -k1,1 -k2,2n > {output}
		'''

rule ploidy_run:
	'''
	Select cells with at least 200k read pairs
	Random sample at 20 fixed intervals up to 200k read pairs
	Calculate breadth of coverage at 0,1,2,3x depth
	'''
	input: "{cell}.frag_trim_dedup.bed"
	output: "{cell}.fragtrimcov.tsv"
	params: 
		chr_auto=config["ref"]["chr_autosomal"]
	shell: "scripts/run_ploidy.sh {input} {params.chr_auto} {output}"

rule ploidy_score: 
	input: 
		files=expand(os.path.join(dnaout,"{cell}/{cell}.fragtrimcov.tsv"), cell=pd.read_csv("meta/ploidy_cells_list.tsv",header=None)[0].tolist()),
		meta="meta/metadata.tsv"
	output: 
		rda="out/ploidy.Rda",
		tsv="out/ploidy_scores.tsv"
	params: 
		fxfile=config["path"]["r_global"]
	script: "scripts/ploidy_score.R"

rule picard_collect: 
	input: expand(os.path.join(dnaout,"{cell}/{cell}-picard-mark_duplicates.txt"), cell=cells_dna.index)
	output: "out/qc_dupstats.tsv"
	script: "scripts/collect_picard.R"

#------------------------------------------ RNA
def get_rna_fq1(wildcards):
	return sorted(cells_rna.loc[wildcards.cell, ['fq1']])

def get_rna_fq2(wildcards):
	return sorted(cells_rna.loc[wildcards.cell, ['fq2']])

rule rna_unzip: 
	input: 
		fq1=get_rna_fq1,
		fq2=get_rna_fq2
	output:
		fq1=temp(config["path"]["temp"] + "/{cell}.1.fq"),
		fq2=temp(config["path"]["temp"] + "/{cell}.2.fq")
	shell:
		'''
		gzip -d -c {input.fq1} > {output.fq1}
		gzip -d -c {input.fq2} > {output.fq2}
		'''

rule rna_prinseq1:
	'''
	Trim first 10bp, low complexity reads and for base quality
	'''
	input:
		fq1=config["path"]["temp"] + "/{cell}.1.fq",
		fq2=config["path"]["temp"] + "/{cell}.2.fq"
	output:
		good1=temp(os.path.join(rnaout,"{cell}/pre/{cell}.prinseqRun1_1.fastq")),
		good2=temp(os.path.join(rnaout,"{cell}/pre/{cell}.prinseqRun1_2.fastq"))
	log: os.path.join(rnaout,"{cell}/log/{cell}.prinseqRun1.log")
	shell: 
		'''
		mkdir -vp $(dirname {output[0]})
		f1={output.good1}
		pfx=${{f1%_1.fastq}} # Create output prefix for prinseq files
		prinseq-lite.pl -fastq {input.fq1} -fastq2 {input.fq2} \
		-out_format 3 -trim_left 3 -min_len 15 -trim_qual_right 25 \
		-lc_method entropy -lc_threshold 65 \
		-out_good $pfx -out_bad null 2> {log}
		'''

rule rna_prinseq2: 
	'''
	Run prinseq again to remove orphans of <30bp
	'''
	input:
		fq1=os.path.join(rnaout,"{cell}/pre/{cell}.prinseqRun1_1.fastq"),
		fq2=os.path.join(rnaout,"{cell}/pre/{cell}.prinseqRun1_2.fastq")
	output:
		good1=temp(os.path.join(rnaout,"{cell}/pre/{cell}.prinseqRun2_1.fastq")),
		good2=temp(os.path.join(rnaout,"{cell}/pre/{cell}.prinseqRun2_2.fastq"))
	log: os.path.join(rnaout,"{cell}/log/{cell}.prinseqRun2.log")
	shell:
		'''
		f1={output.good1}
		pfx=${{f1%_1.fastq}} # Create output prefix for prinseq files
		prinseq-lite.pl -fastq {input.fq1} -fastq2 {input.fq2} \
		-out_format 3 -min_len 15 -out_good $pfx -out_bad null 2> {log}
		'''

rule rna_cutadapt: 
	'''
	Remove Nextera adapter sequence
	'''
	input: os.path.join(rnaout,"{cell}/pre/{cell}.prinseqRun2_1.fastq"), os.path.join(rnaout,"{cell}/pre/{cell}.prinseqRun2_2.fastq")
	output: os.path.join(rnaout,"{cell}/pre/{cell}_val_1.fq"), os.path.join(rnaout,"{cell}/pre/{cell}_val_2.fq")
	log: os.path.join(rnaout,"{cell}/log/{cell}.cutadapt.log")
	shell: 
		'''
		outdir=$(dirname {output[0]})
		fqcdir="$outdir/fastqc_output"
		mkdir -vp $outdir $fqcdir
		trim_galore --paired -a 'CTGTCTCTTATACACATCT' --length 15 --stringency 1 \
		--fastqc_args "-j java --outdir=$fqcdir --extract" \
		--basename {wildcards.cell} -o $outdir {input} 2> {log}
		'''

rule rna_cutadapt_zip: 
# TODO: Add wildcard exceptions to exclude DNA, can then revert to simple paths (without parent path)
	input: 
		fq1=os.path.join(rnaout,"{cell}/pre/{cell}_val_1.fq"),
		fq2=os.path.join(rnaout,"{cell}/pre/{cell}_val_2.fq")
	output: 
		fq1=os.path.join(rnaout,"{cell}/pre/{cell}_val_1.fq.gz"),
		fq2=os.path.join(rnaout,"{cell}/pre/{cell}_val_2.fq.gz")
	shell:
		'''
		gzip {input.fq1}
		gzip {input.fq2}
		'''

rule rna_align:
	'''
	Align preprocessed fastqs with STAR
	'''
	input: os.path.join(rnaout,"{cell}/pre/{cell}_val_1.fq.gz"), os.path.join(rnaout,"{cell}/pre/{cell}_val_2.fq.gz")
	output: temp(os.path.join(rnaout,"{cell}/align/Aligned.sortedByCoord.out.bam"))
	params: ref_dir=config["ref"]["rna_ref_dir"]
	threads: 3
	shell: 
		'''
		outdir=$(dirname {output})
		mkdir -vp $outdir
		echo $outdir
		STAR \
		--genomeDir {params.ref_dir} \
		--readFilesIn {input} \
		--readFilesCommand zcat \
		--outFileNamePrefix "$outdir"/ \
		--outFilterType BySJout \
		--outFilterMultimapNmax 20 \
		--alignSJoverhangMin 8 \
		--alignSJDBoverhangMin 1 \
		--outFilterMismatchNmax 999 \
		--outFilterMismatchNoverLmax 0.04 \
		--alignIntronMin 20 \
		--alignIntronMax 1000000 \
		--alignMatesGapMax 1000000 \
		--outSAMstrandField intronMotif \
		--runThreadN {threads} \
		--outSAMtype BAM SortedByCoordinate \
		--outSAMattributes NH HI NM MD \
		--genomeLoad LoadAndKeep \
		--limitBAMsortRAM 12000000000
		'''
rule rna_dedup: 
	input: os.path.join(rnaout,"{cell}/align/Aligned.sortedByCoord.out.bam")
	output: 
		bam=os.path.join(rnaout,"{cell}/{cell}.sorted.bam"),
		qc=os.path.join(rnaout,"{cell}/log/{cell}-picard-mark_duplicates.txt")
	params: 
		tmp=config["path"]["temp"]
	threads: 2
	log: os.path.join(rnaout,"{cell}/log/{cell}.picard.log")
	shell: 
		'''
		picard MarkDuplicates \
		MAX_RECORDS_IN_RAM=5000000 \
		TMP_DIR={params.tmp} \
		INPUT={input} \
		OUTPUT={output.bam} \
		REMOVE_DUPLICATES=true \
		VALIDATION_STRINGENCY=LENIENT \
		METRICS_FILE={output.qc} 2> {log}
		'''

rule rna_unload_genome: 
	'''
	When all RNAs have been aligned - unload STAR genome to free memory
	Use dedup bam as input to allow deletion of original files when not needed
	'''
	input: expand(os.path.join(rnaout,"{cell}/{cell}.sorted.bam"), cell=cells_rna.index)
	params: ref_dir=config["ref"]["rna_ref_dir"]
	shell: 
		'''
		STAR --genomeDir {params.ref_dir} --genomeLoad Remove
		'''

rule rna_index_bam:
	input: "{cell}.sorted.bam"
	output: "{cell}.sorted.bam.bai"
	shell: "samtools index {input}"

rule rna_sum_expr: 
	input: 
		bam="{cell}.sorted.bam",
		bai="{cell}.sorted.bam.bai"
	output: "{cell}.htseq.out"
	params: 
		genes=config["ref"]["rna_gtf"]
	shell: 
		'''
		samtools view {input.bam} | htseq-count -m intersection-nonempty -s no - {params.genes} > {output}
		'''

rule rna_collect:
	input: 
		cells=expand(os.path.join(rnaout,"{cell}","{cell}.htseq.out"),cell=cells_rna.index)
	output: 
	    counts="out/rna_counts.tsv.gz",
	    qc="out/qc_rnastats.tsv"
	script: 
		"scripts/collect_rna.R"


rule rna_process:
	'''
	Take countsfile + metadata and preprocess RNA with Seurat. Output Seurat Rda object and QC report. 
	'''
	input: 
		counts="out/rna_counts.tsv.gz",
		meta="meta/metadata.tsv",
		dnaclass="out/hct116_dna_lvclasses_500000.tsv",
		dnadamage="out/damage_burden_500000.tsv"
	params: 
		mincount=config["rna"]["min_count"],
		qcut=config["rna"]["quantile_cut"],
		qgene=config["rna"]["quantile_gene"],
		tsne_perpl=config["rna"]["tsne_perplexity"],
		fxfile=config["path"]["r_global"],
		debug=1
	output: 
		rda="out/rna.Rda",
		fig_boxplotByPlate="fig/rnaQC_boxplotByPlate.pdf",
		fig_dupsvsERCC="fig/rnaQC_dups_vs_ercc.pdf",
		fig_umapGrid="fig/rnaQC_umap_grids.pdf",
		fig_tsneGrid="fig/rnaQC_tsne_grids.pdf",
		fig_seuratHeatmaps="fig/rnaQC_seuratHeatmaps.pdf",
		fig_seuratHeatmapDimplot="fig/rnaQC_seuratHeatmapDimplot.pdf",
		fig_markerGrid="fig/rnaQC_markerGeneGrid.pdf"
	script: "scripts/process_rna.R"

#------------------------------------------ Metadata
rule filter_ploidy_cells: 
	'''
	Make list of DNA libraries with >X read pairs to use for ploidy analysis
	'''
	input: "meta/metadata.tsv"
	output: "meta/ploidy_cells_list.tsv"
	params: 
		min_reads=config["dna"]["ploidy_min_reads"]
	run:
		meta=pd.read_table("meta/metadata.tsv", sep="\t", header=0) 
		filt=meta['dna_library_id'][(meta.dna_read_pairs > params["min_reads"])]
		filt.to_csv(output[0], sep="\t", index=False, header=False)