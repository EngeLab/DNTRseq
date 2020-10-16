#!/bin/bash
set -e

bed=$1
ref=$2
out=$3

outdir=$(dirname $bed)/fragcov
mkdir -vp $outdir

pfx=$(basename $bed .bed)

echo "[${pfx}] Calculating read stats"
avg_frag=$(awk '{sum+=($3-$2)}END{print (sum/NR)}' $bed)
tot_frags=$(wc -l $bed | cut -d' ' -f1)
totbp=$(awk -v avg=$avg_frag -v tot=$tot_frags 'BEGIN{printf "%.0f", (avg*tot)}')
echo "[${pfx}] avg=$avg_frag; tot_frags=$tot_frags; totbp=$totbp"
echo "[$pfx] Total nr of fragments: $tot_frags (mean ${avg_frag}bp)"

# Loop over 10 equal parts up to 200k read pairs = fragments
for k in {1..10}; do 
	nfrag=$(printf "%.0f" $(bc -l <<< "scale=2; 200000*$k/10"))
	[[ $nfrag -gt $tot_frags ]] && echo "[$pfx] ^^No more fragments to sample. Stopping" && break
	echo "[$pfx] Running n=${nfrag}"
	[[ -s ${outdir}/${pfx}_${nfrag}.cov ]] && echo "[$pfx] ^^Already run ${nfrag}" && continue
	shuf -n $nfrag $bed | bedtools coverage -hist -a $ref -b stdin | grep ^all > ${outdir}/${pfx}_${nfrag}.cov
done

# # Loop over linear sequence of 20 random samples 
# for k in {1..20}; do 
# 	nfrag=$(printf "%.0f" $(bc -l <<< "scale=2; $tot_frags*$k/20"))
# 	[[ $nfrag -gt $tot_frags ]] && echo "[$pfx] ^^No more fragments to sample. Stopping" && break
# 	echo "[$pfx] Running n=${nfrag}"
# 	[[ -f ${outdir}/${pfx}_${nfrag}.cov ]] && echo "[$pfx] ^^Already run ${nfrag}" && continue
# 	shuf -n $nfrag $bed | bedtools coverage -hist -a $ref -b stdin | grep ^all > ${outdir}/${pfx}_${nfrag}.cov
# done

# if [[ -s ${outdir}/${pfx}_${tot_frags}.cov ]]; then 
# 	echo "[$pfx] ^^Already run full file"
# else 
#     echo "[$pfx] Running full file (n=${tot_frags})"
# 	bedtools coverage -hist -a $ref -b $bed | grep ^all > ${outdir}/${pfx}_${tot_frags}.cov
# fi

echo "[$pfx] Collect fragment stats and combine coverage files"
awk -v avg=$avg_frag -v tot=$tot_frags -v bp=$totbp '{OFS="\t"; print FILENAME,avg,tot,bp,$0}' ${outdir}/*.cov > $out

echo "[$pfx] Finished @ $(date)"
exit 0;
