# Universal Microbiomics Alignment Database
Ideally there should be a Universal Reference Database that:
- can be used on any NGS/MS data type: (meta)-genomics, transcriptomics, metabolomics, proteomics
- and thus is comprised of unified compounds, proteins, and pathway reference data
- has as complete as possible functional annotations and cross-linked ids 
- has all kingdoms of life eukaryotes, bacteria, archaea and MoNA (mobile nucleic acids: viruses, phage, plamids, IS, constructs)
- and thus relies on a well curated and standardized taxonomy database

This URD enables:
1. cross-study comparisons, since they would have the same compound, function, pathway, protein, and organism IDs.
2. integrate different data types: compounds, DNA, RNA, and proteins
3. be used on either environmental and/or host-associated data sets.

### Purpose
The herein described process is designed to create the Protein alignment reference portion of the URD. I've broken them up for my own sanity's sake - there's a lot going on under the hood!!<br>
The 3 other Universal Reference Database creation scripts can be found here:
1. Taxonomy: https://github.com/TealFurnholm/Universal-Taxonomy-Database
2. Non-Coding RNAs: https://github.com/TealFurnholm/Fix_RNACentral_Taxonomy
3. Metabolites and Enzyme Reactions: https://github.com/TealFurnholm/Universal_Biological_Functions_Database

With this data you can align any NGS reads or proteomics data using [**Diamond**](https://github.com/bbuchfink/diamond) alignment software, and get simultaneous taxonomic and functional outputs.
- For primary **metagenome** analysis: https://github.com/TealFurnholm/Strain-Level_Metagenome_Analysis
- For primary **metatranscriptome** analysis: https://github.com/TealFurnholm/Strain-Level_Metatranscriptome_Analysis
- For primary **proteomics** analysis: _TBD_
- To integrate NGS and **metabolomics**: _TBD_ - use compounds and reactions databases to link the two data sets in the meantime
- For secondary **pathway/functional** analysis: https://github.com/TealFurnholm/Meta-omics_Functional_Analysis



# STEP 1. Download UniProt and UniRef
The UniProt-UniRef databases are quite nice as a base. Unfortunately they do not directly link to reactions, compounds, and pathways in a helpful way. Using their Kegg genes, Rhea IDs, and BioCyc monomer IDs, all proteins can be linked with both their compounds and reactions through the 
### General Annotations:
wget -O uniprot-all_6-16-21.tab.gz 'https://www.uniprot.org/uniprot/?query=*&format=tab&force=true&columns=id,protein%20names,length,lineage-id,lineage(GENUS),lineage(SPECIES),organism,feature(SIGNAL),feature(TRANSMEMBRANE),database(TCDB),database(eggNOG),database(Pfam),database(TIGRFAMs),go-id,database(InterPro),ec,database(BioCyc),feature(DNA%20BINDING),feature(METAL%20BINDING),comment(SUBCELLULAR%20LOCATION),database(KEGG),rhea-id&compress=yes'
Sorry, you have to do this MANUALLY.<br>
Go to UniProt, do an empty search, select edit columns, search/select the following IN THIS ORDER!:
 -  Entry
 -  Protein names
 -  Length    
 -  Taxonomic lineage IDs
 -  Taxonomic lineage (GENUS)    
 -  Taxonomic lineage (SPECIES)
 -  Organism
 -  Signal peptide
 -  Transmembrane
 -  TCDB
 -  eggNOG
 -  Pfam
 -  TIGRFAMs
 -  Gene ontology IDs
 -  InterPro
 -  EC number
 -  BioCyc
 -  DNA binding
 -  Metal binding
 -  Subcellular location [CC]
 -  KEGG
 -  Rhea ID
<p></p>
Now save, check box to select all proteins, hit download, change to tab delimited, compressed, and download (be patient!). Move the file to the server where you will be working.

### Create Compounds and Reactions Databases: 
Run this process https://github.com/TealFurnholm/Universal_Biological_Functions_Database/wiki <br>
The scripts outputs
 -  **all_compound_info_[month]_[year].txt:** crosslinks compound ids from PubChem, Chebi, Kegg, Inchi, HMDB, and BioCyc as well as providing compound name, formula, charge, and molecular weight
 -  **all_reaction_info_[month]_[year].txt:** has all reactant compound ids and product ids, reaction direction, and the Kegg, Rhea, and BioCyc reaction ids
 -  






















