# OrthoPhylo
### Developed at Los Alamos National Labs (LANL - C22064)
#### Written by Earl Middlbrook with input from Robab Katani at Penn State and Jeanne Fair at LANL.
#### The software is available through a GLPv3 open source licence. 
####
#### This software is designed to generate phylogenetic trees from bacterial genome assemblies. While many methods use whole genome alignments to generate informative sites to base tress on, OrthoPhylo annotates bacterial genes, identifies orthologous sequences, aligns related proteins to inform transcript alignments, then builds species trees with two methods. The first is a conventional gene concatenation and ML tree estemation method. The second attempts to reconcile gene trees with a unified species tree using quartets (ASTRAL). Both methods allow filtering of gene lists on number of species represented, length, and gappiness in order to tune signal-to-noise ratio for tree estimation. 
#### The main advantages of this software pipeline are three fold: 1) It extends the evolutionary distance input species can represent (over whole genome alignment and k-mer methods) while maintaining phylogenetic resolution, 2) this software is designed to by very user friendly, requiring mostly "conda" install-able dependencies and just a single input directory of genomes for which to reconstruct the evolutionary history and 3) this pipleline is amenable to estimating trees for 1000s of bacterial genome asselblies. To handle large numbers of genomes, the pipline calculates a diversity-representing subset of genomes to run orthofinder on, then expands the OrthoGroups to all samples with HMM searches.

### More detailed discription: 
#### This software is largely a wrapper for open source bioinformatic programs (not packaged within). Genomes are annotated with Prodigal.  If a large number of genomes are to be analyzed, fastANI is used to calculate the pairwise Average Nucleotide Identity (ANI) and, with a custom script, is used to subset the genomes to a minimal number which represent the diversity of the full set. Protein sequences from this subset (or full set for small numbers of genomes) are used as input for Orthofinder to identify orthologous gene families. For the subset method, a HMM search, HMMER, is used to generalize the resulting orthogroup to the full data set. From there, orthogroup proteins are realigned with MAFFT, and these aignments are used to "codon" align the transcript sequences (generated by earlier annotation). These orthogroup transcript alignments are then trimmed (with TRIMAL) and filtered for strict single copy orthologs (SCO_strict) or SCOs found in at least X% of the input genomes (with X being tunable). Next, transcript alignmants are concatenated to generate super-maticies and used as input for species tree generation with either RAxML or fastTREE. Aditionally, gene trees are generated from the individual transcript alignments, which are used in gene tree to species tree estimation with ASTRAL.  

### Cleaning input genomes
#### The script orthophylo/utils/gather_filter_asms.XX.sh was writen to streamiline aquiring all available assemblies for a specific taxon. It takes NCBI's taxID as and input (BrucellaTaxID=234). After the genomes are DL'd there are several simple genome filtering steps: length,N50,GC, etc. There are also a few more soficticated filtering methods: CheckM is used to assess completeness and contamination, while dedupe from BBmap is used to reduce duplicated contigs. 
#### It is absolutely imparative to clean up assemblies to get the best results from OrthoPhylo. For instance, when looking at the output for all Brucella accessions, checkM output shows many assemblies have duplicated "marker genes", if these are from falsely duplicated contigs in the asm, they will lead to removal of the ortholog group from both the strict SCO and relaxed SCO gene sets within the OrthoPhylo workflow. Furthermore, dedupe.sh (from bbmap) identifies many Duplicated or Contained contigs, some of them >100kb long. Again, if any of the duplicated contigs contain what would otherwise be SCOs, the SCOs will be removed. In essence, having duplicates poisons the analysis by severely reducing the number of Orthologs for downstream analysis. This problem is exacerbated by including 1) less curated assemblies (GenBank) and 2) the total number of assemblies fed to the analysis pipeline.

### What this software does not do: 
#### Generate trees that are ready for publication without parameter tuning or manual inspection. Reconstructing trees from whole genomes requires many many steps, all of which have parameters that will differ based on the input sequences. Some importent outputs too look at: input genomes quality (checkM output), assembly subset used for Ortholog model generation (assembly shortlist), number of strict/relaxed single copy orthologs (drops quickly with additional assemblies), phylogenetic signal for transcript alignments, missing data in alignments (per gene and concatenated alignments)...to name a few.

#### This workflow has been tested on a CentOS8 machine, but should be pretty portable. Testing is starting on diverse *nix systems. There will likely be some errors due to the wrapper having moderate bash complexity. 
## Dependencies
+ prodigal 
+ orthofinder 
+ bbmap
+ fastTree
+ hmmer
+ pal2nal
+ prodigal 
+ ete3 
+ raxml*
+ trimal
+ parallel
+ catfasta2phyml
+ ASTRAL
+ fastANI
+ R
+ Alignment_Assessment

*depending on tree method used
## Install
### Conda installable dependencies
If you need to install conda...
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```
#### follow on screen instructions

Set up auto initialize. If you don't, it is up to you to figure out how to make the script happy with your decision

### create conda environment and install dependencies
might take some time....
```
conda create -n orthophylo -c bioconda -c conda-forge \
git prodigal orthofinder \
bbmap fasttree hmmer pal2nal \
prodigal ete3 raxml trimal \
parallel # add bash if you are on macOS
```
to allow changes to take effect. Alternatively, restart shell session
```
source ~/.bashrc
```
to remove parallel's citation reminder (BUT DONT FORGET TO CITE!)
```
conda activate orthophylo
parallel --citation
```
### install R
#### if you are brave, you can install R through conda...it tends to break things, so be warned.
```
conda install -c conda-forge r-essentials
```

#### OR 
#### Go to https://cran.r-project.org/mirrors.html and pick an appropriate mirror. Then chose the linux distro you are using, download package files, and follow install instructions.
#### Test R install
```
Rscript --help
```

### Other dependencies
This reflects how I like ot organize my machine, pick what works for you.
```
Path_to_gits=~/gits
mkdir ~/$Path_to_gits
cd ~/$Path_to_gits/
git clone https://github.com/smirarab/ASTRAL.git
cd ASTRAL
unzip Astral.5.7.8.zip #change to curren version if needed

cd ~/$Path_to_gits/
git clone https://github.com/nylander/catfasta2phyml.git
git clone https://github.com/dportik/Alignment_Assessment.git
cd Alignment_Assessment/
# convert script to python3
2to3 -w Alignment_Assessment_v2.py 
mv Alignment_Assessment_v2.py Alignment_Assessment_v2.py3

cd ~/Downloads
wget https://github.com/ParBLiSS/FastANI/releases/download/v1.33/fastANI-Linux64-v1.33.zip
unzip fastANI-Linux64-v1.33.zip
# For macOS !!!!! NOT SOLVED !!!!!!
# git clone https://github.com/ParBLiSS/FastANI.git
# and follow install instructions. This has been a huge pain on my machine...not solved
mkdir ~/apps/
mv fastANI ~/apps/ # or anywhere else you would like to put it. Change control_file.required to reflect path
```

### Clone the Orthophylo repo, set up control files and run short* test
Requires a username and HTTP key or collaborator status

Cloning takes a minute because of the test files
```
cd ~/$Path_to_gits/
git clone https://github.com/eamiddlebrook/OrthoPhylo.git
cd ortho_phylo1

# to test conda activate orthophylo (within scirpt)
conda deactivate
``` 
## Test install
### run Chloroplast test locally
#### !!! You must change TESTER/control_file.required_chloroplast to reflect the number of core you would like to use
#### Tested on RHEL 8.5 machine with Intel Core i7-8700 CPU and 16gb ram (~8 minute runtime using 3 cores)
```
nano TESTER/control_file.required_chloroplast
bash orthophylo.1.25.control.sh TESTER_chloroplast
```
### run Chloroplast test with slurm scheduler
#### !!! You must change TESTER/control_file.required_chloroplast to reflect the number of core you would like to use
####  change partition, group, ntasks(cores) and mem accordingly in the orthophylo.XXX.control.slurmstart
```
nano TESTER/control_file.required_chloroplast
nano orthophylo.1.25.control.slurmstart
sbatch orthophylo.1.25.control.slurmstart TESTER_chloroplast
```
### !! This got messed up.
#### Should create 4 *.tree files in orthophylo/TESTER//Workflow_test.chloroplast$(date +%m-%d-%Y)/FINAL_SPECIES_TREES/
#### If these are not empty, things are looking really good
#### If they are identical to TESTER/REFERENCE_TESTER_TREES/*.tree, orthophylo/TESTER_chloroplast.pass will be created
#### and "###### It looks like the install was successful ######" should be sent to stdout at the end of the run 



## run bigger test locally
####   this script takes about 20 hr to complete with 20 cores
####   Most of this is concat tree building
####   will make an artificial set of truncated genomes later
#### !!! You must change TESTER/control_file.required to reflect the number of core you would like to use
```
bash orthophylo.1.25.control.sh TESTER
```
## run bigger test with slurm scheduler
####  open orthophylo*.slurmstart and change partition, group, ntasks(cores)  and mem accordingly
```
nano orthophylo.1.25.control.slurmstart
sbatch orthophylo.1.25.control.slurmstart TESTER
```
#### if the test was successful, there should be a file at "orthophylo/TESTER.pass"
#### and "###### It looks like the install was successful ######" should be sent to stdout 
#### towards the end of the run



## known errors:
#### If you have special characters in you fasta file names
Filenames with special characters will likely make this workflow fail. As I encountered these characters, i will add fixes accordingly in orthophylo.XX.sh, function SET_UP_DIR_STRUCTURE.
#### If the combination of fasta file name and contigs within are very long, 
...mafft will truncate the sequence names generated by prodigal. This will make the PAL2NAL module fail on many orthogroups. It manifests as a lot of cat commands failing during trimming because they cant find the transcript file PAL2NAL is suposed to spit out. This can be seen in the log file (in slurm_out). It's normally caused by very redundant sequence identifiers in the contig names. Can be fixed with something like this:
```
sed -i 's/REDUNDANT_STUFF/_/g' < SEQUENCE.fasta
```
#### If genomes are less than ~75% identity... 
FastANI doesnt calculate an ANI value for very divergent genomes. In this case my script ANI_picking.py assigns an ANI of 50 (hard coded) for these comparisons. This should not pose a problem if the number of clades with 75% divergence internally are ~= ANI_shortlist number. The goal of the ANI stuff is to find representetives of the evolutionaly breadth of the dataset, so missing representetives from some clades will likelly not affect the final result (there is a lot of filtering later on in the workflow).
#### If you get an error like: 
``` 
/cm/local/apps/slurm/var/spool/job49381/slurm_script: line 511: J % percent: division by 0 (error token is "percent") 
```
It means that an upstream process failed and the directory used to enumerate a loop is empty

I will attempt to make errors easier to track...

#### I have gotten an error with gather_filter_genomes when pulling extra genomes via wget:
```
dyld: Symbol not found: _libiconv_open
```
This appeared to be a problem with wget, which was resolved by updating through conda.

## Future modules:
+ DONE: add assembly filtering, now very incomplete assemblies will be used, so the number of strict single copy orthologs could be drastically reduced (maybe to zero). Filtering done with gather_filter_asm
+ DONE: Quantify phylogenetic info per gene (https://github.com/dportik/Alignment_Assessment.git)
+ Test evolutionalty models of genes with ete3. Then test GO enrichment
+ Look at tree wide paralog numbers. Do GO analysis...
+ Identify HGT from Transcript alignments. HGTs specific to any group?
+ allow "protected assemblies" when filtering. Perhaps your favorite assembly is crap, but you really want it in the tree. A couple of bad assemblies shouldnt reduce the number of SCOs that much