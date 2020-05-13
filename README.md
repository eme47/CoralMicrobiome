# Characterization of a thermally tolerant _Orbicella faveolata_ reef in Abaco, The Bahamas

**Authors**: Katherine E Parker, Jeremy O Ward, Erin M Eggleston, Evan Fedorov, John Everett Parkinson, Craig P Dahlgren, Ross Cunning
**Journal**: _Coral Reefs_
**Link**: <https://doi.org/10.1007/s00338-020-01948-0>

This work analyses the abiotic and biotic factors contributing to higher thermal tolerance documented at Mermaid Reef in comparison to the neighboring less heat tolerant Sandy Cay Reef in Abaco, Bahamas. Temperature, depth, light intensity, and coral composition were analyzed, as well as host genotype diversity, Symbiodiniaceae composition, and bacterial communities in O. faveolata colonies at each site. This repository contains all analyses of environmental variables and Symbiodiniaceae composition. Metadata for coral samples, raw qPCR data for Symbiodiniaceae quantification, and coral reef benthic community composition and colony size from AGRRA surveys at both sites is available at <https://github.com/kparker96/Sandy_Mermaid>.

**Repository contents:**

**Data:**
* Raw 16S rRNA amplicon data can be downloaded from the NCBI SRA repository via accession number PRJNA591009 <https://www.ncbi.nlm.nih.gov/bioproject/PRJNA591009>

**Code:**
* /coral_microbiome_FsOnly.R: Code for amplicon sequencing data QA/QC, dada2 ASV generation and taxonomic assignment, phyloseq analysis and removal of chloroplast and mitochondrial ASVs, visualization of taxonomic data, and diversity estimates using breakaway and divnet.
