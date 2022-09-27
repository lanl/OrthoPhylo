
alignment_dir=~/scratch/orthophylo/Brucella.full.6.03.2022/phylo_current/AlignmentsTrans.trm.nm/
full_taxa_list=full_taxa_list

cat $alignment_dir/*.fa | \
        grep -e ">" | \
        sort | \
        uniq | \
	sed 's/>//g' \
	> $full_taxa_list

:> taxa_missing_from_alignments
for I in $(ls $alignment_dir/*.taxa)
do
	comm -13 $I $full_taxa_list
done > taxa_missing_from_alignments

cat taxa_missing_from_alignments | sort | uniq -c | \
sed 's/^ *//g' | sort -rnk1,1 > taxa_missing_from_alignments.counts
