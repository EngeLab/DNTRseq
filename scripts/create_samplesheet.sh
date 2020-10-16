#!/bin/bash -l
dnain="input/dna"
rnain="input/rna"

samplesheet_dna="cells_dna.tsv"
samplesheet_rna="cells_rna.tsv"

print_usage() {
  printf "Usage: ./create_samplesheet.sh [-d dna filter] [-r rna filter]"
  printf "\n\t-d Required: one or more DNA filter strings. Ex 'VZA001 VZD002 VZA017' or HCA00102D_A01" 
  printf "\n\t-r Required: one or more RNA filter strings. Ex 'VZA001 VZD002 VZA017' or HCA00102R_A01" 
}

while getopts 'd:r:o:u:' flag; do
  case "${flag}" in
    #i) i_flag=TRUE ;;
    #f) filter=(${OPTARG}) ;;
    d) dna=(${OPTARG}) ;; # Puts them into array
    r) rna=(${OPTARG}) ;;
    o) samplesheet_dna=(${OPTARG}) ;; # Output file
    u) samplesheet_rna=(${OPTARG}) ;; # Output file RNA
    *) print_usage
       exit 1 ;;
  esac
done

#[[ -z $filter ]] && echo "No prefix/filter defined." && print_usage && exit 1;
[[ -z $dna && -z $rna ]] && echo "No prefix/filter provided. Provide either DNA and/or RNA filter" && print_usage && exit 1;

if [[ ! -z $dna ]]; then 
  echo "^^Collecting DNA files with prefix: $dna"
  # Find files and create array
  for f in "${dna[@]}"; do
   files_dna+=($(find ${dnain} -maxdepth 2 -type f -regextype posix-extended -regex ".*(${f}).*R1_001.fastq.gz" | sort))
  done
  # Print sample sheet DNA
  echo "^Creating samplesheet ($samplesheet_dna) with filter: ${dna[@]}"
  echo "^Number of DNA samples: ${#files_dna[@]}"
  echo -e "cell\tfq1\tfq2" > $samplesheet_dna
  for fq1 in ${files_dna[@]}; do
      fq2=${fq1/R1/R2}
      d=$(dirname $fq1)
      file=$(basename $fq1)
      pfx=${file%*_S*_R1*.fastq.gz}
      echo -e "$pfx\t$fq1\t$fq2" >> $samplesheet_dna
  done
fi

if [[ ! -z $rna ]]; then 
  echo "^^Collecting RNA files with prefix: ${rna[@]}"
  for f in "${rna[@]}"; do
   files_rna+=($(find ${rnain} -maxdepth 2 -type f -regextype posix-extended -regex ".*(${f}).*R1_001.fastq.gz" | sort))
  done 
  echo "^Creating samplesheet ($samplesheet_rna) with filter: ${rna[@]}"
  echo "^Number of DNA samples: ${#files_rna[@]}"
  echo -e "cell\tfq1\tfq2" > $samplesheet_rna
  for fq1 in ${files_rna[@]}; do
      fq2=${fq1/R1/R2}
      d=$(dirname $fq1)
      file=$(basename $fq1)
      pfx=${file%*_S*_R1*.fastq.gz}
      echo -e "$pfx\t$fq1\t$fq2" >> $samplesheet_rna
  done
fi

exit 0;
