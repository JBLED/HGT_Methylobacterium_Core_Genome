# HGT_Methylobacterium_Core_Genome
Data and code for "Life on leaves left a lasting mark on the Methylobacterium core genome"

## 1-Assemblies-Methylo-India
## Assembly curation of newly sequenced Methylobacterium genomes
This workflow was used to assess the quality of draft genome assemblies and identify contaminant scaffolds in newly sequenced Methylobacterium isolates included in the study "Life on leaves left a lasting mark on the Methylobacterium core genome."

### A) Assembly statistics
Calculation of assembly size, scaffold number, maximum scaffold size, N50, GC content, and average sequencing coverage for each selected assembly.

### B) Automated scaffold filtering
Removal of short scaffolds (<500 bp).
Removal of low-coverage scaffolds when scaffold length was <2 kb and coverage was <10×.
Recalculation of assembly statistics after filtering.

### C) Manual curation of contaminated assemblies
Identification of assemblies containing genomic fragments with distinct sequencing coverage profiles.
Partitioning of scaffolds into putative genomes based primarily on coverage distributions.
Comparison of GC content among scaffold groups.
Generation of cleaned assemblies and summary statistics for each putative genome.

### Input files
- List-final-assemblies.txt (List of selected assemblies to evaluate)
- Genome assemblies in FASTA format (not provided here; clean assemblies could be retrieved under NBCI BioProject PRJNA1548628) - Draft genome assemblies generated from short-read, long-read, or hybrid assembly approaches.

### Output files
- Sum-Final-assemblies.txt (Assembly statistics calculated from the original assemblies.)
- Sum-Final-assemblies-Megahit-filtered.txt (Assembly statistics after automated scaffold filtering.)
- Sum-Filtering-conta.txt (Statistics for contaminated assemblies after separation into putative genomes.)
- Megahit-Size-Cov-GC.pdf (Global overview of scaffold size, coverage, and GC content distributions across Megahit assemblies.)
- Filtering-Coverage-Isolate_*.pdf and Isolate_*-Megahit-SR*.pdf (Diagnostic plots used during manual inspection of contaminated assemblies.)

### Notes
This workflow is intended to document the assembly curation process leading to the genome dataset used for downstream comparative genomic analyses. It does not reproduce the testing and benchmarking of alternative assembly methods performed during pipeline development. Only the assembly quality assessment, scaffold filtering, and contamination screening steps retained for the final dataset are included.

## 2-Core-genome 
## Definition and characterization of the Methylobacterium core genome
This workflow was used to define the Methylobacterium core genome and generate the sequence, functional, and evolutionary datasets used throughout the study Life on leaves left a lasting mark on the Methylobacterium core genome.

### A) Genome and annotation statistics
Genome assemblies and RAST annotations are summarized to calculate genome size; number of contigs; N50; GC content; numbers of protein-coding genes, RNA genes, repeats, and hypothetical proteins.

### B) Core genome definition
The core genome is defined from RAST annotations across all Methylobacterium genomes. Candidate core genes are selected according to: 1) occurrence in at least 80% of genomes within each major phylogenetic group; 2) average copy number ≤ 1.5; 3) conserved sequence length across genomes. Genes remaining present in at least 70% of genomes after sequence-length quality filtering are retained in the final core genome.

### C) Core gene extraction
For each retained core gene, nucleotide sequences are extracted from all genomes. One FASTA file is generated per genome;
one FASTA file is generated per core gene.

### D) Core gene alignment
RNA genes are aligned directly using ClustalW. Protein-coding genes are aligned through amino-acid-guided alignments before being converted back to codon-preserving nucleotide alignments. Sites with excessive missing data are removed to generate the final alignments used in downstream analyses.

### E) Average nucleotide identity (ANI)
Pairwise average nucleotide identity is calculated using FastANI. ANI values are used to characterize genomic similarity among isolates, delimit species-level clusters and generate species assignment files for downstream phylogenetic analyses.

### F) Functional annotation of the core genome
Representative protein sequences of core genes are exported for annotation with KEGG Automatic Annotation Server (KAAS). The resulting annotations are used to characterize the functional composition of the core genome.

### G) Evolutionary statistics
For each protein-coding core gene, the workflow calculates: median dN/dS ratio; alignment size; number of polymorphic sites; number of informative sites. These variables are used in subsequent analyses of gene evolution and phylogenetic discordance.

### Input files
- data.txt (Metadata table describing the genomes included in the analysis)
- RAST_files/ (RAST annotation files used to identify and characterize core genes)
- fastANI-out-raw.tsv (Raw pairwise ANI values generated with FastANI)
- KEGG.txt (KEGG/KAAS annotations assigned to representative core genes)

**not provided here:** genome assemblies (FASTA genome assemblies used for genome statistics, ANI calculations, and core-gene extraction)

### Main output files
- Stat-genomes.txt (Summary statistics for each genome and its annotation content)
- Annotation-occurence-raw.txt (Occurrence matrix of all RAST annotations across genomes)
- Methylobacterium-core-gene-occ.txt (Presence/absence and copy-number matrix of the final core genes)
- Core-gene-code.txt (Reference table linking core-gene identifiers to their functional annotations)
- core-genes/*_alignedNF.fas (Final filtered alignments of individual core genes; intermediate files are not provided)
- Summary-per-core-gene-after-alignment.txt (Summary statistics for all core-gene alignments)
- Methylobacterium-ListSpecies-Sept-2025.txt (Species assignment of each genome based on ANI)
- Classif-Methylo-ANI-sept2025.txt (group assignment of each genome based on different ANI thresholds)
- SEQaa-CoreMeth-CBMB20.fa (Representative protein sequences submitted for KEGG annotation)
- MAP-core-genes-CBMB20.txt (position of core genes on the CBMB20 reference chromosome)
- CoreGenes_nPOL_dNdS.txt (Per-gene evolutionary statistics, including dN/dS ratios and sequence polymorphism metrics)

### Notes
This workflow defines the core genome using annotation-based orthology inferred from RAST annotations. Core genes are subsequently aligned, functionally annotated, and characterized using ANI, dN/dS, and sequence polymorphism metrics, generating the datasets used for downstream comparative and phylogenomic analyses.
