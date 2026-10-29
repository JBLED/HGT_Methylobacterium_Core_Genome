# HGT_Methylobacterium_Core_Genome
Data and code for "Life on leaves left a lasting mark on the Methylobacterium core genome"

# 1-Assemblies-Methylo-India

## Assembly curation of newly sequenced Methylobacterium genomes
This workflow was used to assess the quality of draft genome assemblies and identify contaminant scaffolds in newly sequenced Methylobacterium isolates included in the study "Life on leaves left a lasting mark on the Methylobacterium core genome."

The workflow consists of three main steps:

## A) Assembly statistics
Calculation of assembly size, scaffold number, maximum scaffold size, N50, GC content, and average sequencing coverage for each selected assembly.

## B) Automated scaffold filtering
Removal of short scaffolds (<500 bp).
Removal of low-coverage scaffolds when scaffold length was <2 kb and coverage was <10×.
Recalculation of assembly statistics after filtering.

## C) Manual curation of contaminated assemblies
Identification of assemblies containing genomic fragments with distinct sequencing coverage profiles.
Partitioning of scaffolds into putative genomes based primarily on coverage distributions.
Comparison of GC content among scaffold groups.
Generation of cleaned assemblies and summary statistics for each putative genome.

## Input files
- List-final-assemblies.txt (List of selected assemblies to evaluate)
- Genome assemblies in FASTA format (not provided here; clean assemblies could be retrieved under NBCI BioProject PRJNA1548628) - Draft genome assemblies generated from short-read, long-read, or hybrid assembly approaches.

## Output files
- Sum-Final-assemblies.txt (Assembly statistics calculated from the original assemblies.)
- Sum-Final-assemblies-Megahit-filtered.txt (Assembly statistics after automated scaffold filtering.)
- Sum-Filtering-conta.txt (Statistics for contaminated assemblies after separation into putative genomes.)
- Megahit-Size-Cov-GC.pdf (Global overview of scaffold size, coverage, and GC content distributions across Megahit assemblies.)
- Filtering-Coverage-Isolate_*.pdf and Isolate_*-Megahit-SR*.pdf (Diagnostic plots used during manual inspection of contaminated assemblies.)

## Notes
This workflow is intended to document the assembly curation process leading to the genome dataset used for downstream comparative genomic analyses. It does not reproduce the testing and benchmarking of alternative assembly methods performed during pipeline development. Only the assembly quality assessment, scaffold filtering, and contamination screening steps retained for the final dataset are included.
