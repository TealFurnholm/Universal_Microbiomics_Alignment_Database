# UMRAD: Universal Microbiomics Reference and Alignment Database
This script compiles the data from databases 1-3 below, plus functional information from UniProt and other public databases, to annotate the UniRef100 set of deduplicated proteins. Collectively these annotations enable the unification of various meta-omics studies: metabolomes, metaproteomes, metatranscriptomes, and metagenomes. This is #4 of a series of pipelines that comprise this Universal Reference for various 'omics data analysis.
    <br>1. Universal Taxonomy Database: [found here](https://github.com/TealFurnholm/Universal-Taxonomy-Database)
    <br>2. Universal Compounds Database: [found here](https://github.com/TealFurnholm/Universal_Biological_Compounds_Database)
    <br>3. Universal Reactions Database: [found here](https://github.com/TealFurnholm/Universal_Biological_Reactions_Database)
    <br>4. Universal Protein Alignment Database: *this repository*
    <br>5. Universal ncRNA Alignment Database: [found here](https://github.com/TealFurnholm/Fix_RNACentral_Taxonomy)


### Goal
Ideally there should be a Universal Reference Database that:
- can be used on any NGS/MS data type
- is comprised of unified compounds, proteins, reactions and pathway reference data
- has as complete as possible functional annotations and cross-linked ids 
- has all kingdoms of life eukaryotes, bacteria, archaea and MoNA (mobile nucleic acids: viruses, phage, plamids, IS, constructs)
- and thus relies on a well curated and standardized taxonomy database

### UMRAD usage:
1. With this data you can align any NGS reads or proteomics data using [**Diamond**](https://github.com/bbuchfink/diamond) alignment software, and get simultaneous whole community phylogenetic and functional outputs.
2. A centralized repo allows cross-study comparisons, since they would have the same compound, function, pathway, protein, and taxon IDs.
3. By containing all sequenced organisms, it can be used on either environmental and/or host-associated data sets.

- For primary **metagenome** analysis: https://github.com/TealFurnholm/Strain-Level_Metagenome_Analysis
- For primary **metatranscriptome** analysis: https://github.com/TealFurnholm/Strain-Level_Metatranscriptome_Analysis
- For primary **proteomics** analysis: _TBD_
- To integrate NGS and **metabolomics**: _TBD_ - use compounds and reactions databases to link the two data sets in the meantime
- For secondary **pathway/functional** analysis: https://github.com/TealFurnholm/Meta-omics_Functional_Analysis


Get Started [HERE](https://github.com/TealFurnholm/Universal_Microbiomics_Alignment_Database/wiki) 
