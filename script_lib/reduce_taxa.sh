#!/bin/bash
# TLDR: reduce problemeatic taxa in gene sets to maximize Single Copy Orthologs
# This little script is meant to filter many multi-fasta alignments in $alignemt_dir
#  for taxa with a sequence in $num_genes OGs
#  then subset genes for single copy orthologs
#
# Run from intended output directory
# Usage: reduce_taxa.sh FULL_PATH_TO_ALIGNMENT_DIR min_num_of_gens_per_taxa
#

alignment_dir=$1
num_genes=$2

# get taxa with sequences in >= $num_genes OGs
cat $alignment_dir/*.fa | \
	grep -e ">" | \
	sort | \
	uniq -c | \
	sed 's/^ *//g;s/>//g' | \
	sort -rnk1,1 | \
	awk -v num_genes=${num_genes} '{if ($1 > num_genes) print $2}' | \
	sort \
	> taxa_with_${num_genes}.tmp

export num_taxa=$(cat taxa_with_${num_genes}.tmp | wc -l)

:> OG_SCO_${num_taxa}_taxa_${num_genes}_genes
for I in $alignment_dir/*.fa
do
	# get list of taxa for each Alignment
	cat $I | grep -e ">" | sed 's/>//g' | sort | uniq \
		> ${I%.*}.taxa
	# compare alignment list of taxa with taxa with sequences in >= $num_genes OGs
	if [ $(comm -12 taxa_with_${num_genes}.tmp ${I%.*}.taxa | wc -l) -eq $num_taxa ]
	then
		base=$(basename ${I%.*.*.*.*})
		echo $base >> ./OG_SCO_${num_taxa}_taxa_${num_genes}_genes
	fi
done
cat ./OG_SCO_${num_taxa}_taxa_${num_genes}_genes | wc -l
