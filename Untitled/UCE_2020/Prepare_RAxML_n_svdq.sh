
#set working directory

################################################################################
                        #Preparing data for RAxML
################################################################################

# build the concatenated data matrix for raxml
phyluce_align_format_nexus_files_for_raxml \
    --alignments mafft-nexus-edge_internal-trimmed-gblocks-clean-50p \
    --output mafft-nexus-edge_internal-trimmed-gblocks-clean-50p-raxml \
    --charsets \
    --log-path log

# build the concatenated data matrix for raxml
phyluce_align_format_nexus_files_for_raxml \
    --alignments mafft-nexus-edge_internal-trimmed-gblocks-clean-70p \
    --output mafft-nexus-edge_internal-trimmed-gblocks-clean-70p-raxml \
    --charsets \
    --log-path log

# build the concatenated data matrix for raxml
phyluce_align_format_nexus_files_for_raxml \
    --alignments mafft-nexus-edge_internal-trimmed-gblocks-clean-90p \
    --output mafft-nexus-edge_internal-trimmed-gblocks-clean-90p-raxml \
    --charsets \
    --log-path log

# build the concatenated data matrix for raxml
phyluce_align_format_nexus_files_for_raxml \
    --alignments mafft-nexus-edge_internal-trimmed-gblocks-clean-100p \
    --output mafft-nexus-edge_internal-trimmed-gblocks-clean-100p-raxml \
    --charsets \
    --log-path log

################################################################################
                        #Preparing data for svdq
################################################################################

# build the concatenated data matrix for svdquartets
phyluce_align_format_nexus_files_for_raxml \
    --alignments mafft-nexus-edge_internal-trimmed-gblocks-clean-50p \
    --output mafft-nexus-edge_internal-trimmed-gblocks-clean-50p-svdq \
    --charsets \
    --nexus \
    --log-path log

# build the concatenated data matrix for svdquartets
phyluce_align_format_nexus_files_for_raxml \
    --alignments mafft-nexus-edge_internal-trimmed-gblocks-clean-70p \
    --output mafft-nexus-edge_internal-trimmed-gblocks-clean-70p-svdq \
    --charsets \
    --nexus \
    --log-path log

# build the concatenated data matrix for svdquartets
phyluce_align_format_nexus_files_for_raxml \
    --alignments mafft-nexus-edge_internal-trimmed-gblocks-clean-90p \
    --output mafft-nexus-edge_internal-trimmed-gblocks-clean-90p-svdq \
    --charsets \
    --nexus \
    --log-path log

# build the concatenated data matrix for svdquartets
phyluce_align_format_nexus_files_for_raxml \
    --alignments mafft-nexus-edge_internal-trimmed-gblocks-clean-100p \
    --output mafft-nexus-edge_internal-trimmed-gblocks-clean-100p-svdq \
    --charsets \
    --nexus \
    --log-path log
