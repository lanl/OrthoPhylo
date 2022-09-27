#!/bin/bash
# --------------------------------------------------------------
### PART 1: Requests resources to run your job.
# --------------------------------------------------------------
### Optional. Set the job name
#SBATCH --job-name=JOBNAME.checkm
### Optional. Set the output filename.
### SLURM reads %x as the job name and %j as the job ID
#SBATCH --output=slurm_out/%x-%j.out
### REQUIRED. Specify the PI group for this job
#SBATCH --account=USER
### REQUIRED. Set the partition for your job.
#SBATCH --partition=PARTITION
### REQUIRED. Set the number of cores that will be used for this job.
#SBATCH --nodes=1
#SBATCH --ntasks=30
### REQUIRED. Set the memory required for this job.
#SBATCH --mem=200gb
### REQUIRED. Specify the time required for this job, hhh:mm:ss
#SBATCH --time=200:00:00

conda activate gather_genomes

wd=~/PATH_TO_WORKING_DIRECTORY
input=$wd/INPUT/
output=$wd/OUTPUT/
checkM_type="genome"
suffix="fna"


threads=30
# Set these if you hav a good idea of the expected values
# if not, set to "default" and they will be set as the average+-(3*stddev)
# of course stddev is a bit weird because the values are almost surely not normal..
# NOTE: N50 is checkM's definition (length of contig, not number)
MIN_LEN="500000"
MAX_LEN="5500000"
MIN_N50="1500"
MIN_GC="0.20"
MAX_GC="0.750"
MAX_dup="0.03"
MIN_completeness="97"
MAX_contam="2.0"

MAX_dup_default="0.01"
MIN_completness_default="98"
MAX_contam_default="0.2"


get_stats_with_checkM () {
	echo "############################################################"
        echo "###### Run CheckM to get completeness, contamination  ######"
        echo "####### and general assembly  metrics for filtering  #######"
        echo "############################################################" 
	# checkM has the option to write stdout to file...
	#   dont know if this would make aggregating its output
	#   stdout is really messy
	max_genomes=200
	threads=$threads
	checkM_input=$1
        checkM_dir=$2
	checkM_type=$3
	suffix=$4
	if [[ $checkM_type == "protien" ]]
	then
		checkM_args="-t $threads -g -x $suffix"
	elif [[ $checkM_type == "genome" ]]
	then
		checkM_args="-t $threads -x $suffix"
	else
		echo "Unknown checkM input type" && exit
	fi
        mkdir $checkM_dir
        cd $checkM_dir || exit
        # split assemblies into different directories
        J=0
	K=0
	for I in $(ls $checkM_input/*.$suffix)
        do
          	if [ $((J % max_genomes)) -eq 0 ]
                then
                    	K=$((K+1))
                        mkdir $checkM_dir/${checkM_type}_${K}
                        mkdir $checkM_dir/${checkM_type}_${K}_out
                        cd $checkM_dir/${checkM_type}_${K}
                fi
                # make simlinks for assembly/proteome subsets
                ln -s $I ./
                J=$((J+1))
        done

        # run checkM on each proteome subset
        cd $checkM_dir
        J=1
	while [ $J -le $K ]
        do
          	in=$checkM_dir/${checkM_type}_${J}
                out=$checkM_dir/${checkM_type}_${J}_out
                checkm lineage_wf ${checkM_args} $in $out
                J=$((J+1))
        done

        # Aggregate checkM output
	echo "acc lineage #markerGenes #genomes_based_on missing 1copy 2copy 3copy 4copy 5+copy duplication_ratio completeness contamination GC GC_std Genome-size #scaffs scaff_N50" \
		> $wd/assemblies_all.stats.txt
	cat $checkM_dir/${checkM_type}_*_out/storage/bin_stats_ext.tsv | \
		sed 's/{,//g;s/,//g' | \
		awk '{print $1,$4,$10,$14,$16,$18,$20,$22,$24,$26,($26*5+$24*4+$22*3+$20*2)/$10,$28,$30,$32,$35,$38,$45,$57}' \
		>> $wd/assemblies_all.stats.txt
}


filter_asm_by_stats () {
	echo "###############################################"
        echo "#### filter out assemblies with low N50s, ####"
        echo "#### low total length and extra length  ######"
	echo "##############################################"
	cd $wd
	MIN_LEN=$1
	MAX_LEN=$2
	MIN_N50=$3
	MIN_GC=$4
        MAX_GC=$5
	MAX_dup=$6
	MIN_completeness=$7
	MAX_contam=$8
	if [ $MAX_contam == '' ]
        then
            	echo "!!!!! Not enough args given...needs 8 !!!!!!!!"
                echo "please look at function filter_asm_by_stats for details"
                echo "Set variables....
                MIN_LEN=$1
                MAX_LEN=$2
                MIN_N50=$3
                MIN_GC=$4
                MAX_GC=$5
                MAX_dup=$6
                MIN_completeness=$7
                MAX_contam=$8
                "
                exit
        fi
	echo "MIN_LEN=$MIN_LEN"
        echo "MAX_LEN=$MAX_LEN"
        echo "MIN_N50=$MIN_N50"
        echo "MIN_GC=$MIN_GC"
        echo "MAX_GC=$MAX_GC"
	echo "MAX_dup=$MAX_dup"
        echo "MIN_completeness=$MIN_completeness"
        echo "MAX_contam=$MAX_contam"

	#create a file with averages and stdDevs for columns
	#   in checkM_stats_aggregated
	#   scaf_bp  scaf_N50 gc_avg dup complete contam
	cat assemblies_all.stats.txt | \
	awk '{print $16,$18,$14,$11,$12,$13}' |\
	awk '{for(i=1;i<=NF;i++) {sum[i] += $i; sumsq[i] += ($i)^2}}
        	END {for (i=1;i<=NF;i++) {
         	printf "%f %f ", sum[i]/NR, sqrt((sumsq[i]-sum[i]^2/NR)/NR)}
         }' > assemblies_all.stats.ave_stddev

	# set filtering defaults if not set at start of script
	#   This is too complicated for bash....should have writen in python.
	while read scaf_bp scaf_bp_stddev scaf_N50 scaf_N50_stddev gc_avg gc_avg_stddev dup_ave dup_stddev complete_ave complete_stddev contam_ave contam_ave_stddev
	do
		# if [ "default" in variables ]
		if [ $MIN_LEN == "default" ]
		then
                        MIN_LEN=$(bc -l <<< "$scaf_bp-($scaf_bp_stddev*3)")
		fi
		if [ $MAX_LEN == "default" ]
                then
                        MAX_LEN=$(bc -l <<< "$scaf_bp+($scaf_bp_stddev*3)")
		fi
		if [ $MIN_N50 == "default" ]
		then
			MIN_N50=$(bc -l <<< "$ctg_N50+($ctg_N50_stddev*3)")
		fi
		if [ $MIN_GC == "default" ]
                then
			MIN_GC=$(bc -l <<< "$gc_avg-($gc_avg_stddev*3)")
                fi
                if [ $MAX_GC == "default" ]
                then
			MAX_GC=$(bc -l <<< "$gc_avg+($gc_avg_stddev*3)")
                fi
		if [ $MAX_dup == "default" ]
                then
                    	MAX_dup=$MAX_dup_default
                fi
		if [ $MIN_completeness == "default" ]
                then
                    	MAX_complete=$MAX_complete_default
                fi
		if [ $MAX_contam == "default" ]
                then
                    	MIN_contam=$MIN_contam_default
                fi
	done <<< $(cat assemblies_all.stats.ave_stddev)


	# Filter assemblies
	#   This is also too complicated for bash...really slow
	echo "Filtering assemblies based on:"
        echo "MIN_LEN=$MIN_LEN"
        echo "MAX_LEN=$MAX_LEN"
        echo "MIN_N50=$MIN_N50"
        echo "MIN_GC=$MIN_GC"
        echo "MAX_GC=$MAX_GC"
	echo "MAX_dup=$MAX_dup"
        echo "MIN_completeness=$MIN_completeness"
        echo "MAX_contam=$MAX_contam"
	echo "Paths to assemblies being filtered out are found in assemblies_to_remove.stats"
	# Emply output file
	:> assemblies_to_remove.stats
	cat assemblies_all.stats.txt |\
	tail -n +2 |\
	awk '{print $1,$16,$18,$14,$11,$12,$13}' |\
	while read acc scaf_bp ctg_N50 gc_avg dup complete contam
	do
		if (( $(echo "$scaf_bp < $MIN_LEN" |bc -l) )) || \
                   (( $(echo "$scaf_bp > $MAX_LEN" |bc -l) )) || \
                   (( $(echo "$ctg_N50 < $MIN_N50" |bc -l) )) || \
                   (( $(echo "$gc_avg < $MIN_GC" |bc -l) )) || \
                   (( $(echo "$gc_avg > $MAX_GC" |bc -l) )) || \
		   (( $(echo "$dup > $MAX_dup" |bc -l) ))  || \
		   (( $(echo "$complete < $MIN_completeness" |bc -l) )) || \
		   (( $(echo "$contam > $MAX_contam" |bc -l) ))
		then
			#echo "$scaf_bp $ctg_N50 $gc_avg $dup $complete $contam"
			echo ${acc} >> assemblies_to_remove.stats
		fi
	done
}

#get_stats_with_checkM $input $output $checkM_type $suffix
filter_asm_by_stats $MIN_LEN $MAX_LEN $MIN_N50 $MIN_GC $MAX_GC $MAX_dup $MIN_completeness $MAX_contam
