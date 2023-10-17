# Commands used in evalauting UCE sequence data

## Count the read data

for i in *_R_1*.fastq.gz;
do echo $i; gunzip -c $i | wc -l | awk '{print $1/4}';
done

## Clean the read data

source activate phyluce

### we need to create a configuration file that we will use to inform the program about adapters

illumiprocessor \
    --input raw-fastq/ \
    --output clean-fastq \
    --config illumiprocessor.conf \
    --cores 4
## Quality control

# move to the directory holding our cleaned reads
cd clean-fastq/

# run this script against all directories of reads

for i in *;
do
    phyluce_assembly_get_fastq_lengths --input $i/split-adapter-quality-trimmed/
     --csv;
done

#Assemble the data
#create the file assembly.conf with the input below

[samples]
alligator_mississippiensis:/work/LAS/phylo-lab/kevinq/git_repositories/Fig_wasp_phylogeography/etc/uce-tutorial/clean-fastq/alligator_mississippiensis/split-adapter-quality-trimmed/
anolis_carolinensis:/work/LAS/phylo-lab/kevinq/git_repositories/Fig_wasp_phylogeography/etc/uce-tutorial/clean-fastq/anolis_carolinensis/split-adapter-quality-trimmed/
gallus_gallus:/work/LAS/phylo-lab/kevinq/git_repositories/Fig_wasp_phylogeography/etc/uce-tutorial/clean-fastq/gallus_gallus/split-adapter-quality-trimmed/
mus_musculus:/work/LAS/phylo-lab/kevinq/git_repositories/Fig_wasp_phylogeography/etc/uce-tutorial/clean-fastq/mus_musculus/split-adapter-quality-trimmed/

#make sure we are at te top-levl of our uce tutorial directory
cd uce-tutorial


module unuse /opt/rit/spack-modules/lmod/linux-rhel7-x86_64/Core
module use /opt/rit/modules

module load java/1.7.0_55
module load bowtie/1.1.2

# run the assembly
phyluce_assembly_assemblo_trinity \
    --conf assembly.conf \
    --output trinity-assemblies \
    --clean \
    --cores 12

##Finding UCE Loci

#getting UCE probe set (I will get this from Jordan)
wget https://raw.githubusercontent.com/faircloth-lab/uce-probe-sets/master/uce-5k-probe-set/uce-5k-probes.fasta

phyluce_assembly_match_contigs_to_probes \
    --contigs trinity-assemblies/contigs \
    --probes uce-5k-probes.fasta \
    --output uce-search-results


 #Extracting UCE Loci

#create a conf file , which specifies the taxa we want
#see taxon-set.conf

#create an output directory for this taxon set - this just keeps
#things neat
cd uce-tutorial
mkdir -p taxon-sets/all

# create the data matrix configuration file (imcomplete)
phyluce_assembly_get_match_counts \
    --locus-db uce-search-results/probe.matches.sqlite \
    --taxon-list-config taxon-set.conf \
    --taxon-group 'all' \
    --incomplete-matrix \
    --output taxon-sets/all/all-taxa-incomplete.conf


 ## complete matrix
 phyluce_assembly_get_match_counts \
    --locus-db uce-search-results/probe.matches.sqlite  \
    --taxon-list-config taxon-set.conf \
    --taxon-group 'all' \
    --output taxon-sets/all/all-taxa-complete.conf

##extra FASTA data that correspond to the loci in all-taxa incomplete.conf

# change to the taxon-sets/all directory
cd taxon-sets/all

# make a log directory to hold our log files - this keeps things neat
mkdir log

# get FASTA data for taxa in our taxon set (using incomplete matrix)
phyluce_assembly_get_fastas_from_match_counts \
    --contigs ../../trinity-assemblies/contigs \
    --locus-db ../../uce-search-results/probe.matches.sqlite \
    --match-count-output all-taxa-incomplete.conf \
    --output all-taxa-incomplete.fasta \
    --incomplete-matrix all-taxa-incomplete.incomplete \
    --log-path log


#using complete matrix

phyluce_assembly_get_fastas_from_match_counts \
    --contigs ../../trinity-assemblies/contigs \
    --locus-db ../../uce-search-results/probe.matches.sqlite \
    --match-count-output all-taxa-complete.conf\
    --output all-taxa-complete.fasta

# explode the monolithic FASTA by taxon (you can also do by locus)
phyluce_assembly_explode_get_fastas_file \
    --input all-taxa-incomplete.fasta \
    --output exploded-fastas \
    --by-taxon


# get summary stats on the FASTAS
for i in exploded-fastas/*.fasta;
do
    phyluce_assembly_get_fasta_lengths --input $i --csv;
done

# samples,contigs,total bp,mean length,95 CI length,min length,max length,median legnth,contigs >1kb
alligator-mississippiensis.unaligned.fasta,4069,2694866,662.291963627,3.29432860587,224,2579,675.0,165
anolis-carolinensis.unaligned.fasta,685,388032,566.470072993,9.18273050108,224,1039,537.0,7
gallus-gallus.unaligned.fasta,3883,2772837,714.096574813,3.80175813481,224,1594,729.0,442
mus-musculus.unaligned.fasta,730,517085,708.335616438,10.3709178277,224,1150,809.5,100


## Aligning UCE loci
#when taxa are "closely" related (<30-50 MYA, perhaps), I think that edge-trimming alignments is reasonable
# when the taxa span a wider range of divergence time (>50 MYA), you mat want to think about interanl
#trimmin

#runn this the taxon-sets/all directory
cd uce-tutorial/taxon-sets/all

# align the data (incomplete matrix)
#did not work for me
phyluce_align_seqcap_align \
    --fasta all-taxa-incomplete.fasta \
    --output mafft-nexus-edge-trimmed \
    --taxa 4 \
    --aligner mafft \
    --cores 12 \
    --incomplete-matrix \
    --output-format fasta \
    --log-path log
#summary stats for alignment

phyluce_align_get_align_summary_data \
    --alignments mafft-nexus-edge-trimmed \
    --cores 12 \
    --log-path log
#aligning data internal trimming
# make sure we are in the correct directory
cd uce-tutorial/taxon-sets/all

# align the data - turn off trimming and output FASTA
phyluce_align_seqcap_align \
    --fasta all-taxa-incomplete.fasta \
    --output mafft-nexus-internal-trimmed \
    --taxa 4 \
    --aligner mafft \
    --cores 12 \
    --incomplete-matrix \
    --output-format fasta \
    --no-trim \
    --log-path log

# run gblocks trimming on the alignments
    phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed \
        --alignments mafft-nexus-internal-trimmed \
        --output mafft-nexus-internal-trimmed-gblocks \
        --cores 12 \
        --log log

#summmary stats
phyluce_align_get_align_summary_data \
    --alignments mafft-nexus-internal-trimmed-gblocks \
    --cores 12 \
    --log-path log

#Alignment cleaning
# make sure we are in the correct directory
cd uce-tutorial/taxon-sets/all

# align the data - turn off trimming and output FASTA
phyluce_align_remove_locus_name_from_nexus_lines \
    --alignments mafft-nexus-internal-trimmed-gblocks \
    --output mafft-nexus-internal-trimmed-gblocks-clean \
    --cores 12 \
    --log-path log


##Final data matrices

# make sure we are in the correct directory
cd uce-tutorial/taxon-sets/all

# the integer following --taxa is the number of TOTAL taxa
# and I use "75p" to denote the 75% complete matrix
#change --percet to change amount of missing data
phyluce_align_get_only_loci_with_min_taxa \
    --alignments mafft-nexus-internal-trimmed-gblocks-clean \
    --taxa 4 \
    --percent 0.75 \
    --output mafft-nexus-internal-trimmed-gblocks-clean-75p \
    --cores 12 \
    --log-path log
# 90% complete matrix
phyluce_align_get_only_loci_with_min_taxa \
    --alignments mafft-nexus-internal-trimmed-gblocks-clean \
    --taxa 4 \
    --percent 0.90 \
    --output mafft-nexus-internal-trimmed-gblocks-clean-90p \
    --cores 12 \
    --log-path log


#Preparing data for RAxML and ExaML

# make sure we are in the correct directory
cd uce-tutorial/taxon-sets/all

# build the concatenated data matrix
phyluce_align_format_nexus_files_for_raxml \
    --alignments mafft-nexus-internal-trimmed-gblocks-clean-75p \
    --output mafft-nexus-internal-trimmed-gblocks-clean-75p-raxml \
    --charsets \
    --log-path log

# build the concatenated data matrix with 90p matrix
phyluce_align_format_nexus_files_for_raxml \
    --alignments mafft-nexus-internal-trimmed-gblocks-clean-90p \
    --output mafft-nexus-internal-trimmed-gblocks-clean-90p-raxml \
    --charsets \
    --log-path log
