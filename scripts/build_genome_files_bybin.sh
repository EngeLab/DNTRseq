#!/bin/bash
set -e

DIR=$1
SCRIPTS=$2
BUILD=$3
MASTERREF=$4
REFGENES=$5
READLEN=$6
BINSIZE=$7

[[ -z $DIR ]] && echo "DIR not set, reverting to default" && DIR=/wrk/resources/dntr/resources/genome_files
[[ -z $SCRIPTS ]] && echo "SCRIPTS not set, reverting to default" && SCRIPTS=$DIR/ginkgo_scripts
[[ -z $BUILD ]] && BUILD=hg38
[[ -z $MASTERREF ]] && echo "MASTERREF not set, reverting to default" && MASTERREF=/wrk/resources/genomes/hg38-iGenome-ERCC/hg38-iGenome-ERCC.fa
[[ -z $REFGENES ]] && echo "REFGENES not set, reverting to default" && REFGENES=/wrk/resources/genomes/hg38-iGenome-ERCC/hg38_RefSeq.bed
[[ -z $READLEN ]] && READLEN=37
[[ -z $BINSIZE ]] && BINSIZE=100000

# Prerequisites: ginkgo .cpp files, compiled
# Bed-file with genes
# gnu parallel

if [[ $HOSTNAME = *.uppmax.uu.se ]]; then 
	uppmax=T
	echo "^^^ Runnning on Uppmax"
	module load bioinfo-tools gnuparallel bwa BEDTools
	MASTERREF=/proj/sens2018557/resources/hg38-iGenome-ERCC/hg38-iGenome-ERCC.fa
	GENES=/proj/sens2018557/resources/hg38-iGenome-ERCC/gencode.v28.annotation.bed
	DIR=/proj/sens2018557/resources/dntr
fi

# Create softlinks to current dir (to allow indexing etc)
REF=$DIR/$(basename $MASTERREF)
GENES=$DIR/$(basename $REFGENES)
[[ ! -f $REF ]] && ln -s $MASTERREF $REF
[[ ! -f $GENES ]] && ln -s $REFGENES $GENES

cd $DIR

echo $DIR
echo $SCRIPTS

# Split fasta file
if [[ ! -s chrlist ]]
then
	mkdir -vp split
	cd split
	csplit -qz $REF /\>chr.*/ {*}
	for a in x*; do contig=$(head -1 $a | awk '{gsub(">","");print $1}'); mv -v $a ${contig}.fa; done
	cd ..
	ls split | egrep chr[0-9,X,Y]\{1,2\}.fa | awk -F".fa" '{print $1}' | sort -V > chrlist
	echo ">> Selected contigs:"
	cat chrlist
else 
	echo ">>> Split fasta files found"
fi

# Index with bwa ~2h
if [[ ! -s ${REF}.bwt || ! -s ${REF}.sa ]]
then 
	nohup ./index_bwa.sh $REF > index-bwa.out &
else 
	echo ">>> bwa index found"
fi

# Simulate reads in parallel
# TODO: gzip if simReads is made to output to stdout instead of hardcoding to file
mkdir -vp simreads
if [ ! -f frag_${READLEN}_done ] #TODO: Actually check each chr_frag file
then
echo -e "\n 1. Simulating $READLEN bp Reads (binsize: $BINSIZE)"
parallel --joblog simread.log --resume-failed --progress \
$SCRIPTS/simReads split/{1}.fa simreads/{1}_${READLEN}_frags ${READLEN} :::: chrlist
touch frag_${READLEN}_done
else 
	echo ">>> Simulated reads for ${READLEN}bp found"
fi

# Generate accessory chromosome files
if [ ! -s lengths ]
then
	echo -e "\n 2. Generating Essential Chromosome Files (binsize: $BINSIZE)"
	echo -e "  [Computing chromomsome lengths]"
	while read CHROM; do 
	echo $CHROM $((`grep -v ">" split/${CHROM}.fa | wc -c`-`wc -l < split/${CHROM}.fa`+1))
	done < chrlist > lengths
else
	echo ">>> Chromosome length file (lengths) found"
fi

if [ ! -s centromeres ] #34s
then
  echo "  [Computing centromere positions] (binsize: $BINSIZE)"
  for i in `cat chrlist`; do 
    ${SCRIPTS}/findCentromeres split/${i}.fa cenout${i}
    if [ `wc -l < cenout${i}` -eq 0 ]; then
      echo -e "${i}\t0\t0"
    else
      awk '{print $2-$1"\t"$1"\t"$2}' cenout${i} | sort -rg | head -1 | awk -v chr=$i '{print chr"\t"$2"\t"$3}'
    fi
  done > centromeres
  rm -f cenout*
else
	echo ">>> Centromere position file (centromeres) found"
fi

echo -e "\n 3. Generating fixed-length interval files (binsize: $BINSIZE)"
if [ ! -s fixed_${BINSIZE} ]
then
	echo ">>> Creating fixed bins for $BINSIZE (binsize: $BINSIZE)"
	${SCRIPTS}/fixed lengths fixed_${BINSIZE} ${BINSIZE}
else 
	echo "^^^ Fixed bins for $BINSIZE already exist (binsize: $BINSIZE)"
fi
echo ">>> Fixed bounds (binsize: $BINSIZE)"
[[ ! -s bounds_fixed_${BINSIZE} ]] && ${SCRIPTS}/bounds fixed_${BINSIZE} bounds_fixed_${BINSIZE} || echo "^^^ Fixed bounds for $BINSIZE already exist"
echo ">>> Fixed bins GC content (binsize: $BINSIZE)"
if [ ! -s GC_fixed_${BINSIZE} ]
then
	mkdir -vp split
	cd split
	echo ">>> Calculating GC tables for fixed bins (binsize: $BINSIZE)"
	${SCRIPTS}/GC ../fixed_${BINSIZE} ../GC_fixed_${BINSIZE} lengths
	cd ..
else 
	echo "^^^ GC table for fixed bins for $BINSIZE already exist"
fi
echo ">>> Fixed bins genes (binsize: $BINSIZE)"
[[ ! -s genes_fixed_${BINSIZE} ]] && ${SCRIPTS}/match_genes_to_bins fixed_${BINSIZE} genes genes_fixed_${BINSIZE} || echo "^^^ Genelist for fixed $BINSIZE already exist"

echo -e "\n 4. Mapping sampled reads back to reference (binsize: $BINSIZE)"
# On uppmax
if [[ $uppmax == T ]]; then
	echo ">>> Distributing mapping jobs on slurm"
	mkdir -vp bwa_mapped
	running=0
	for chr in `cat chrlist`; do
		[[ -s bwa_mapped/${chr}_${READLEN}_bwa_done ]] && echo "^^^ $chr at ${READLEN}bp already mapped" && continue
		jobrunning=$(sacct -u vasilios -o 'JobID,JobName,Elapsed,State' | grep RUNNING | grep -w "map-${chr}")
		[[ ! -z $jobrunning ]] && echo -e "^^^ $chr at ${READLEN}bp is currently running:\n$jobrunning" && running=$((running+1)) && continue
		echo "sbatch -o ${chr}_${READLEN}_bwa.out map_bwa.sh $chr $READLEN"
		sbatch -o bwa_mapped/${chr}_${READLEN}_bwa.out -J map-${chr} map-bwa.sh $DIR $REF $chr $READLEN
		running=$((running+1))
	done
	[[ $running > 0 ]] && echo "Wait for $running jobs to finish!" && exit 0
else 
	echo ">>> Running @ $HOSTNAME (binsize: $BINSIZE)"
	mkdir -vp bwa_mapped
	for chr in `cat chrlist`; do
		[[ -s bwa_mapped/${chr}_${READLEN}_bwa_done ]] && echo "^^^ $chr at ${READLEN}bp already mapped" && continue
		bwa mem -t 2 -M -k 16 ${REF} simreads/${chr}_${READLEN}_frags | \
		awk -v chr=$chr '$12==\"NM:i:0\" && $3==chr && $5>=30{ print $3\"\t\"$4 }' | \
		sort -k2,2g > bwa_mapped/${chr}_${READLEN}_bwa_done
	done
fi

echo -e "\n 5. Bin mapped reads (binsize: $BINSIZE)"
# ex: ${SCRIPTS}/bin results_${binsize}_${READLEN}_bwa_${chr} centromeres ${binsize} ${READLEN} ${chrlen} bwa_mapped/${chr}_${READLEN}_bwa_done
# chrlist and lengths are linked (by ::::+)
# colsep allows use of only the 2nd column in lengths file (which then becomes argument 5 instead of 4)
while read CHROM; do
if [ ! -s results_${BINSIZE}_${READLEN}_bwa_${CHROM} ]; then
echo ">>> Binning reads for $CHROM (binsize: $BINSIZE)"
CHRLEN=$(awk -v c=$CHROM '$1==c{print $2}' lengths)
cmd="${SCRIPTS}/bin results_${BINSIZE}_${READLEN}_bwa_${CHROM} centromeres $BINSIZE $READLEN $CHRLEN bwa_mapped/${CHROM}_${READLEN}_bwa_done"
echo "Command: $cmd"
$cmd
else 
echo "^^^ Binned reads for $CHROM found"
fi
done < chrlist

### Create variable bins
echo ">>> Variable bins (binsize: $BINSIZE)"
# Run if file is empty/doesnt exist or has only 1 line
if [[ ! -s variable_${BINSIZE}_${READLEN}_bwa || $(wc -l <variable_${BINSIZE}_${READLEN}_bwa) -le 1 ]]
then 
	awk 'BEGIN{print "CHR\tEND"}{print}' results_${BINSIZE}_${READLEN}_bwa_* | sort -V > variable_${BINSIZE}_${READLEN}_bwa
fi

echo ">>> Variable bin bounds (binsize: $BINSIZE)"
if [[ ! -s bounds_variable_${BINSIZE}_${READLEN}_bwa ]]
then 
	${SCRIPTS}/bounds variable_${BINSIZE}_${READLEN}_bwa bounds_variable_${BINSIZE}_${READLEN}_bwa
fi

echo ">>> Variable bin GC content (binsize: $BINSIZE)"
if [ ! -s GC_variable_${BINSIZE}_${READLEN}_bwa ]
then
	mkdir -vp split
	cd split
	${SCRIPTS}/GC ../variable_${BINSIZE}_${READLEN}_bwa ../GC_variable_${BINSIZE}_${READLEN}_bwa ../lengths
	cd ..
fi

echo ">>> Variable gene listings (binsize: $BINSIZE)"
if [[ ! -s genes_variable_${BINSIZE}_${READLEN}_bwa ]]
then
	${SCRIPTS}/match_genes_to_bins variable_${BINSIZE}_${READLEN}_bwa genes genes_variable_${BINSIZE}_${READLEN}_bwa
fi

echo -e "\n\n\n\t\t >>> Finished creating genome files for ${BINSIZE}bp @ ${READLEN}bp \n\n\n"
touch build_genome_files_done_${BINSIZE}
