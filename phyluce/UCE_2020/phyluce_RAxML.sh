
################################################################################
                                  RAxML
################################################################################

# make sure we are in the correct directory
cd uce-tutorial/taxon-sets/all/mafft-nexus-internal-trimmed-gblocks-clean-75p-raxml

# get two random numbers
for i in 1 2; do echo $RANDOM; done

# that output the following
#21869
#1046

# run the search for the "best" ML tree
raxmlHPC-PTHREADS-SSE3 \
    -f a \
    -m GTRGAMMA \
    -N 100 \
    -x 12345 \
    -p 21869 \
    -n results_KQ_pet \
    -s phased-mafft-nexus-post-edge_internal-trimmed-gblocks-clean-70p-raxml/phased-mafft-nexus-post-edge-internal-trimmed-gblocks-clean-70p.phylip \
    -T 36
