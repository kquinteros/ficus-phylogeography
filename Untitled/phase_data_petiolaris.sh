#!/bin/bash

#Submit this script with: sbatch thefilename

#SBATCH -t 12:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 16   # total number of processor cores in this job
#SBATCH -J "phase_petiolaris"   # job name
#SBATCH --mail-user=jsatler@iastate.edu   # email address
#SBATCH --mail-type=END



# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
source activate phyluce167


####################################
###          Petiolaris          ###
####################################

#################################
### Phasing loci into alleles ###
#################################

#current working directory
cd /ptmp/LAS/phylo-lab/jsatler/phyluce/taxon-sets-trinity/petiolaris

#explode the edge trimmed alignments
phyluce_align_explode_alignments \
    --alignments mafft-nexus-edge-trimmed \
    --input-format fasta \
    --output mafft-nexus-edge-trimmed-exploded \
    --by-taxon
    
#get summary stats on the FASTAS
for i in mafft-nexus-edge-trimmed-exploded/*.fasta;
do
    phyluce_assembly_get_fasta_lengths --input $i --csv;
done

#map reads to contigs
phyluce_snp_bwa_multiple_align \
    --config phasing_petiolaris.conf \
    --output multialign-bams \
    --cores 16 \
    --log-path log \
    --mem

#phase alleles
phyluce_snp_phase_uces \
    --config phasing_petiolaris.conf \
    --bams multialign-bams \
    --output multialign-bams-phased-reads

#realign phased loci, but do not trim
phyluce_align_seqcap_align \
    --fasta multialign-bams-phased-reads/fastas/joined_allele_sequences_all_samples.fasta \
    --output phased-mafft-nexus-edge-trimmed \
    --taxa 52 \
    --aligner mafft \
    --no-trim \
    --ambiguous \
    --cores 16 \
    --incomplete-matrix \
    --output-format fasta \
    --log-path log

#summary stats for the phased alignments
phyluce_align_get_align_summary_data \
    --alignments phased-mafft-nexus-edge-trimmed \
    --cores 16 \
    --input-format fasta \
    --log-path log

# run gblocks with internal trimming on the phased alignments
phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed \
    --alignments phased-mafft-nexus-edge-trimmed \
    --output phased-mafft-nexus-edge_internal-trimmed-gblocks \
    --cores 16 \
    --log log

#summary stats for the phased gblock-trimmed alignments
phyluce_align_get_align_summary_data \
    --alignments phased-mafft-nexus-edge_internal-trimmed-gblocks \
    --cores 16 \
    --log-path log

# align the data - turn off trimming and output FASTA
phyluce_align_remove_locus_name_from_nexus_lines \
    --alignments phased-mafft-nexus-edge_internal-trimmed-gblocks \
    --output phased-mafft-nexus-edge_internal-trimmed-gblocks-clean \
    --cores 16 \
    --log-path log


####################################################
### Output final matrices with % of missing data ###
####################################################

# the integer following --taxa is the number of TOTAL taxa
# and I use "100p" to denote the 100% complete matrix
phyluce_align_get_only_loci_with_min_taxa \
    --alignments phased-mafft-nexus-edge_internal-trimmed-gblocks-clean \
    --taxa 52 \
    --percent 1.0 \
    --output phased-mafft-nexus-edge_internal-trimmed-gblocks-clean-100p \
    --cores 16 \
    --log-path log


# the integer following --taxa is the number of TOTAL taxa
# and I use "90p" to denote the 90% complete matrix
phyluce_align_get_only_loci_with_min_taxa \
    --alignments phased-mafft-nexus-edge_internal-trimmed-gblocks-clean \
    --taxa 52 \
    --percent 0.90 \
    --output phased-mafft-nexus-edge_internal-trimmed-gblocks-clean-90p \
    --cores 16 \
    --log-path log

# the integer following --taxa is the number of TOTAL taxa
# and I use "70p" to denote the 70% complete matrix
phyluce_align_get_only_loci_with_min_taxa \
    --alignments phased-mafft-nexus-edge_internal-trimmed-gblocks-clean \
    --taxa 52 \
    --percent 0.70 \
    --output phased-mafft-nexus-edge_internal-trimmed-gblocks-clean-70p \
    --cores 16 \
    --log-path log

# 50% complete matrix
phyluce_align_get_only_loci_with_min_taxa \
    --alignments phased-mafft-nexus-edge_internal-trimmed-gblocks-clean \
    --taxa 52 \
    --percent 0.50 \
    --output phased-mafft-nexus-edge_internal-trimmed-gblocks-clean-50p \
    --cores 16 \
    --log-path log


#################################################
### Build input files for RAxML/svdq analysis ###
#################################################

# build the concatenated data matrix for raxml
phyluce_align_format_nexus_files_for_raxml \
    --alignments phased-mafft-nexus-edge_internal-trimmed-gblocks-clean-100p \
    --output phased-mafft-nexus-edge_internal-trimmed-gblocks-clean-100p-raxml \
    --charsets \
    --log-path log

# build the concatenated data matrix for svdquartets
phyluce_align_format_nexus_files_for_raxml \
    --alignments phased-mafft-nexus-edge_internal-trimmed-gblocks-clean-100p \
    --output phased-mafft-nexus-edge_internal-trimmed-gblocks-clean-100p-svdq \
    --charsets \
    --nexus \
    --log-path log

# build the concatenated data matrix for raxml
phyluce_align_format_nexus_files_for_raxml \
    --alignments phased-mafft-nexus-edge_internal-trimmed-gblocks-clean-90p \
    --output phased-mafft-nexus-edge_internal-trimmed-gblocks-clean-90p-raxml \
    --charsets \
    --log-path log

# build the concatenated data matrix for svdquartets
phyluce_align_format_nexus_files_for_raxml \
    --alignments phased-mafft-nexus-edge_internal-trimmed-gblocks-clean-90p \
    --output phased-mafft-nexus-edge_internal-trimmed-gblocks-clean-90p-svdq \
    --charsets \
    --nexus \
    --log-path log

# build the concatenated data matrix for raxml
phyluce_align_format_nexus_files_for_raxml \
    --alignments phased-mafft-nexus-edge_internal-trimmed-gblocks-clean-70p \
    --output phased-mafft-nexus-edge_internal-trimmed-gblocks-clean-70p-raxml \
    --charsets \
    --log-path log

# build the concatenated data matrix for svdquartets
phyluce_align_format_nexus_files_for_raxml \
    --alignments phased-mafft-nexus-edge_internal-trimmed-gblocks-clean-70p \
    --output phased-mafft-nexus-edge_internal-trimmed-gblocks-clean-70p-svdq \
    --charsets \
    --nexus \
    --log-path log

# build the concatenated data matrix for raxml
phyluce_align_format_nexus_files_for_raxml \
    --alignments phased-mafft-nexus-edge_internal-trimmed-gblocks-clean-50p \
    --output phased-mafft-nexus-edge_internal-trimmed-gblocks-clean-50p-raxml \
    --charsets \
    --log-path log

# build the concatenated data matrix for svdquartets
phyluce_align_format_nexus_files_for_raxml \
    --alignments phased-mafft-nexus-edge_internal-trimmed-gblocks-clean-50p \
    --output phased-mafft-nexus-edge_internal-trimmed-gblocks-clean-50p-svdq \
    --charsets \
    --nexus \
    --log-path log


