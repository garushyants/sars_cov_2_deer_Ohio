### This repository contains scripts and data for mutational analysis presented in "SARS-CoV-2 spatial dispersal and fast evolution in free-ranging white-tailed deer in Ohio" by McBride et al. 2023


#### *scripts*
Contains R, python and shell scripts  to reproduce the analysis.\
General scripts are provided at the directory root.\
More specific sripts are grouped by folders depending on the type of analysis.

**root-to-tip_regression** contains *root_to_tip_regression.py* that allows to calculate number of mutations from root to tip by type and gene, assuming that phylogenetic tree as well states at each node are known. *root_to_tip_regression_global.R* allows to visuzlize obtained counts and perform some statistical analysis. *convert_dates_to_float.R* is a technical script to convert dates to floats if needed.
\
**mutational_contexts** allows to reproduce analysis of mutational contexts provided in the paper. The main script is *reconstruct_mutational_contexts.R*\
\
Scripts in **lineages_frequencies** allows to manipulate with outbreak.info data for mutations of interest.

*subset_tree_and_alignment.py* removes all sequences with missing characters for particular gene from alignment and phylogenetic tree. This script is required to run HyPhy.\
*vizualize_clusters.py* allows to vizualize clusters and mutations happening at each branch and requires the set of clusters, vcf and phylognetic tree.\
*parse_treetime.py* generates vcf files from treetime trees with reconstracted ancestral states.

#### *data* 
 Contains phylogenetic trees, information about clusters and singletons, and vcf files with reconstructed mutations.

*data/Delta/Delta.clusters.snpeff.corrected.vcf* - corrected vcf file with mutations in Delta clusters\
*data/Delta/Delta.singletons.spneff.vcf* - corrected vcf file with mutations in Delta singletons\
*data/Delta/Delta.clusters.final.csv* contains information about Delta\
*data/Alpha/Alpha.clusters.snpeff.corrected.vcf* - corrected vcf file with mutations in Alpha clusters\
*data/Alpha/Alpha.clusters* contains information about Alpha clusters\
*data/Alpha/global_dataset* contains files for root-to-tip regression analysis\
\
\
\
Folder *figures* contains graphs generated with the provided scripts.

Folder *supplementary_tables* contains Table S6 describing all mutations in clusters and Table S8 describing recurrent mutations in deer clusters.
