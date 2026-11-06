# HGT_Methylobacterium_Core_Genome
Data and code for **"Life on leaves left a lasting mark on the Methylobacterium core genome"**

**Preprint available here: [https://doi.org/10.64898/2026.11.02.763778](https://www.biorxiv.org/content/10.64898/2026.11.02.763778v1)**

**not provided here:** Large datasets are available through Zenodo (**doi [10.5281/zenodo.23768560](https://doi.org/10.5281/zenodo.23768560)**)
- genome assemblies (FASTA genome assemblies used for genome statistics, ANI calculations, and core-gene extraction)
- RAST_files (RAST annotation files used to identify and characterize core genes)

## 1-Assemblies-Methylo-India - *Assembly curation of newly sequenced Methylobacterium genomes*
This workflow was used to assess the quality of draft genome assemblies and identify contaminant scaffolds in newly sequenced Methylobacterium isolates included in the study "Life on leaves left a lasting mark on the Methylobacterium core genome."

#### A) Assembly statistics
Calculation of assembly size, scaffold number, maximum scaffold size, N50, GC content, and average sequencing coverage for each selected assembly.
#### B) Automated scaffold filtering
Removal of short scaffolds (<500 bp).
Removal of low-coverage scaffolds when scaffold length was <2 kb and coverage was <10×.
Recalculation of assembly statistics after filtering.
#### C) Manual curation of contaminated assemblies
Identification of assemblies containing genomic fragments with distinct sequencing coverage profiles.
Partitioning of scaffolds into putative genomes based primarily on coverage distributions.
Comparison of GC content among scaffold groups.
Generation of cleaned assemblies and summary statistics for each putative genome.

#### Input files
- List-final-assemblies.txt (List of selected assemblies **used fore Table S2**)
- Genome assemblies in FASTA format (not provided here; clean assemblies could be retrieved under NBCI BioProject PRJNA1548628) - Draft genome assemblies generated from short-read, long-read, or hybrid assembly approaches.

#### Output files
- Sum-Final-assemblies.txt (Assembly statistics calculated from the original assemblies)
- Sum-Final-assemblies-Megahit-filtered.txt (Assembly statistics after automated scaffold filtering)
- Sum-Filtering-conta.txt (Statistics for contaminated assemblies after separation into putative genomes)
- Megahit-Size-Cov-GC.pdf (Global overview of scaffold size, coverage, and GC content distributions across Megahit assemblies)
- Filtering-Coverage-Isolate_xxx.pdf and Isolate_xxx-Megahit-SRxx.pdf (Diagnostic plots used during manual inspection of contaminated assemblies **used for Figure S1**)

#### Notes
This workflow is intended to document the assembly curation process leading to the genome dataset used for downstream comparative genomic analyses. It does not reproduce the testing and benchmarking of alternative assembly methods performed during pipeline development. Only the assembly quality assessment, scaffold filtering, and contamination screening steps retained for the final dataset are included.

## 2-Core-genome - *Definition and characterization of the Methylobacterium core genome*
This workflow was used to define the Methylobacterium core genome and generate the sequence, functional, and evolutionary datasets used throughout the study Life on leaves left a lasting mark on the Methylobacterium core genome.

#### A) Genome and annotation statistics
Genome assemblies and RAST annotations are summarized to calculate genome size; number of contigs; N50; GC content; numbers of protein-coding genes, RNA genes, repeats, and hypothetical proteins.
#### B) Core genome definition
The core genome is defined from RAST annotations across all Methylobacterium genomes. Candidate core genes are selected according to: 1) occurrence in at least 80% of genomes within each major phylogenetic group; 2) average copy number ≤ 1.5; 3) conserved sequence length across genomes. Genes remaining present in at least 70% of genomes after sequence-length quality filtering are retained in the final core genome.
#### C) Core gene extraction
For each retained core gene, nucleotide sequences are extracted from all genomes. One FASTA file is generated per genome;
one FASTA file is generated per core gene.
#### D) Core gene alignment
RNA genes are aligned directly using ClustalW. Protein-coding genes are aligned through amino-acid-guided alignments before being converted back to codon-preserving nucleotide alignments. Sites with excessive missing data are removed to generate the final alignments used in downstream analyses.
#### E) Average nucleotide identity (ANI)
Pairwise average nucleotide identity is calculated using FastANI. ANI values are used to characterize genomic similarity among isolates, delimit species-level clusters and generate species assignment files for downstream phylogenetic analyses.
#### F) Functional annotation of the core genome
Representative protein sequences of core genes are exported for annotation with KEGG Automatic Annotation Server (KAAS). The resulting annotations are used to characterize the functional composition of the core genome.
#### G) Evolutionary statistics
For each protein-coding core gene, the workflow calculates: median dN/dS ratio; alignment size; number of polymorphic sites; number of informative sites. These variables are used in subsequent analyses of gene evolution and phylogenetic discordance.

#### Input files
- data.txt (Metadata table describing the genomes included in the analysis; **used for Table S1**)
- fastANI-out-raw.tsv (Raw pairwise ANI values generated with FastANI)
- KEGG.txt (KEGG/KAAS annotations assigned to representative core genes)

#### Main output files
- Stat-genomes.txt (Summary statistics for each genome and its annotation content; **used for Table S1**)
- Annotation-occurence-raw.txt (Occurrence matrix of all RAST annotations across genomes)
- Methylobacterium-core-gene-occ.txt (Presence/absence and copy-number matrix of the final core genes; **used for Table S1**)
- Core-gene-code.txt (Reference table linking core-gene identifiers to their functional annotations; **used for Table S3**)
- core-genes/*_alignedNF.fas (Final filtered alignments of individual core genes; intermediate files are not provided)
- Summary-per-core-gene-after-alignment.txt (Summary statistics for all core-gene alignments; **used for Table S3**)
- Methylobacterium-ListSpecies-Sept-2025.txt (Species assignment of each genome based on ANI)
- Classif-Methylo-ANI-sept2025.txt (group assignment of each genome based on different ANI thresholds)
- SEQaa-CoreMeth-CBMB20.fa (Representative protein sequences submitted for KEGG annotation)
- MAP-core-genes-CBMB20.txt (position of core genes on the CBMB20 reference chromosome)
- CoreGenes_nPOL_dNdS.txt (Per-gene evolutionary statistics, including dN/dS ratios and sequence polymorphism metrics; **used for Table S3**)
- CoreGenome_N50.pdf **used for Figure S2**
- DistANI.pdf **used for Figure S3**

#### Notes
This workflow defines the core genome using annotation-based orthology inferred from RAST annotations. Core genes are subsequently aligned, functionally annotated, and characterized using ANI, dN/dS, and sequence polymorphism metrics, generating the datasets used for downstream comparative and phylogenomic analyses.

## 3-Phylogeny - *Phylogenomic of Methylobacterium core genome* 
This workflow reconstructs Methylobacterium phylogenies from individual core genes and quantifies phylogenetic discordance among gene trees.

#### A) Core-gene phylogenies
Maximum-likelihood phylogenies are inferred independently for each core gene using RAxML. Gene trees are filtered according to bootstrap support and reformatted for downstream species-tree reconstruction.
#### B) Species-tree reconstruction
Collections of core-gene phylogenies are combined using ASTRAL/ASTER to infer lineage and species trees while accounting for gene-tree discordance. Separate analyses are performed using all core genes or only protein-encoding genes.
#### C) Quantification of gene-tree discordance
For each protein-coding core gene, Robinson-Foulds (RF) distances are calculated relative to the consensus lineage tree. RF values are then modeled against gene-specific properties, including: number of polymorphic sites; dN/dS ratio; number of taxa represented; median node support. Residual RF distances are used as corrected measures of phylogenetic discordance.

#### Main input and output files
(list excluding files from the previous workflows)
- RAxML_bipartitionsBranchLabels.CG_***-GTRGAMMA-1000 (raw Maximum-likelihood gene trees generated with RAxML)
- Methylobacterium-Multi-coregene-sept-2025_mBS=10.tre (multiple core-gene trees combined in one file used as input for ASTER)
- Methylobacterium-Multi-coregene-PEG-sept-2025_mBS=10.tre (all protein-encoding core-gene trees combined in one file used as input for ASTER)
- Methylobacterium-Multi-coregene-PEG-SPE-sept-2025_mBS=10.tre (all protein-encoding core gene trees combined and labels simplified at the species level to determine the species tree)
- Methylobacterium-AsterLinTree-sept-2025_mBS=10.tre (Consensus lineage tree inferred from all core genes)
- Methylobacterium-AsterLinTree-PEG-sept-2025_mBS=10.tre (Consensus lineage tree inferred from protein-coding genes only - **used for Figure S4**)
- Methylobacterium-AsterSpeTree-PEG-SPE-sept-2025_mBS=10.tre (Consensus species tree inferred from protein-coding genes - **used for Figure 1A**)
- RFdist-Correction.pdf (Diagnostic plots used to model RF distances as a function of gene evolutionary properties - **used for Figure S7**)
- ANOVA-RFdist-Correction.txt (Statistical model used to correct RF distances **used for Table S7**)
- GeneTreeStat.txt (Summary table containing RF distances, node support, dN/dS, polymorphism statistics, and corrected RF residuals for each core gene; **used for Table S3**)

## 4 - CoreGeneColocalization - *Core gene colocalization, phylogenetic discordance, and functional coherence*
This workflow quantifies the relationship between physical distance among core genes on the CBMB20 reference genome, similarity of evolutionary histories. Pairwise comparisons among all protein-coding core genes are performed using gene positions in the reference genome, gene-tree discordance (RF distances), and KEGG annotations. Breakpoint analyses are then used to identify the genomic scale at which gene colocalization is associated with shared evolutionary histories and functions.

**Input files are from the previous workflows**

#### Output files
- Summary-KEGG-RF-Pos-per-PEG.txt (Summary table containing, for each protein-coding core gene, its annotation, genomic position, RF residual value, and KEGG assignment; **used for Table S3**)
- pDistanceRF.txt (Pairwise comparison table containing physical distance, RF distance, corrected RF residual distance, and KEGG overlap for every pair of protein-coding core genes)
- PairwiseDistRF-KEGG-per-window.txt (Statistics calculated within discrete genomic distance windows, including mean RF residual distance and proportion of gene pairs sharing at least one KEGG pathway)
- PairwiseDistRF-KEGG-per-window-cumulated.txt (Same statistics as above, calculated using cumulative genomic distance windows)
- PairwiseDistRF-KEGG.pdf (Diagnostic plots describing relationships among genomic distance, RF distance, and KEGG similarity)
- PairwiseDistRF-KEGG+Segmented-raw.pdf (Segmented regression analysis performed on the complete set of pairwise gene comparisons)
- PairwiseDistRF-KEGG+Segmented.pdf (Final segmented regression analysis identifying genomic breakpoints associated with increased phylogenetic discordance and reduced functional similarity among core genes. This figure is used to define the genomic scale of conserved gene clusters in subsequent analyses - **Used for Figures S5 A and B** )

## 5 - KEGG-enrichment - *KEGG pathway enrichment analysis*
This workflow identifies KEGG pathways enriched in genes with elevated phylogenetic discordance (high residual RF values, RES). For each KEGG pathway, the mean RES value is compared against a null distribution generated by random permutation of RES values among core genes while preserving KEGG pathway assignments. Statistical significance is evaluated using permutation tests and corrected for multiple testing.

#### A) Permutation analysis
For each KEGG pathway: calculate the observed mean RES value; randomly permute RES values among core genes; recompute mean RES values for each permutation; estimate enrichment p-values and adjusted p-values.
#### B) Visualization
Generate the final enrichment figure from permutation results.

**Input files are from the previous workflows**

#### Output files
- Perm-RF-Summary-Per-KEGG.txt (Permutation-based enrichment statistics for all KEGG pathways, including observed mean RES values, permutation counts, p-values, and Bonferroni-adjusted p-values - **Used for Table S5**)
- KEGG-enrichment.pdf (Graphical summary of KEGG enrichment results showing the relationship between mean RES values and statistical significance for each KEGG pathway. Point size is proportional to the number of genes assigned to each pathway. Significantly enriched pathways are highlighted - **used for Figure S8**)

## 6 - ConservedClusterOfGenes - *Conserved gene clusters and genome organization*
This workflow identifies pairs of core genes that remain colocalized across Methylobacterium genomes and quantifies the conservation of gene neighborhoods among species and lineages. Cluster conservation (ACC) is compared with genome similarity (ANI) and phylogenetic discordance to investigate the relationship between genome organization and gene evolutionary history.

#### A) Calculation of Average Cluster Conservation (ACC) and ANI among species
Core-gene pairs located within 25 kb are identified in each genome. Pairwise cluster conservation among genomes is calculated and summarized at the species level. Average ANI values among species are calculated from FastANI results.
#### B) – Conserved gene clusters across lineages
Core-gene pairs are aggregated across species and lineages to identify clusters that remain conserved at different phylogenetic depths. Permutation tests are used to identify clusters enriched in genes with elevated residual RF values.
#### C) – Calculation of lineage-level colocalization scores
For each pair of core genes, PEGs colocalization conservation across all lineages (CCL) is calculated. These values are used for network reconstruction and downstream analyses.
#### D) – Gene colocalization and phylogenetic discordance
CCL is compared with pairwise phylogenetic dissimilarity (pRF residual distance) to determine whether genes that remain colocated across lineages also exhibit more similar evolutionary histories.
#### E) – Species and lineage summary statistics
ANI and ACC values are summarized at the species and lineage levels. The relationship between genome similarity and gene-cluster conservation is quantified.

#### Input files (Non-listed input files are generated by previous workflows)
- Taxonomy.txt (Species and lineage assignments for all genomes; **Used for Table S4**)

#### Outputs
- Average-cluster-cons-per-species.txt and Average-cluster-per-species_yellow.pdf - Heatmap showing average conservation of gene clusters among species (**used for Figure 1B**)
- Average-ANI-per-species.txt and Average-ANI-per-species_GREY.pdf - Heatmap showing average ANI among species (**used for Figure 1B**)
- clustering-genes-lineages.txt - Classification of core genes according to ACC
- clustering-genes-lineages-pVal.txt - Statistical support for ACC across lineages (**used for Table S8**)
- AverageCLU-per-LIN.txt - Cluster conservation statistics summarized at the lineage level
- Pairwise_RF-vs_colocalization-raw.pdf - Relationship between pairwise RF distances and gene colocalization using raw values (**used for Figure S5C**)
- Pairwise_RF-vs_colocalization-bins.pdf - Relationship between pairwise RF distances and gene colocalization after binning gene pairs by conservation level (**used for Figure S5C**)
- Stat-per-species.txt - Summary statistics of ANI and cluster conservation at the species level (**used for Table S4**)
- Stat-per-lineage.txt - Summary statistics of ANI and cluster conservation at the lineage level (**used for Table S4**)
- Mean_ANI_ACC_among_lineage.txt - Pairwise ANI and cluster conservation values among lineages (**used for Table S4**)
- ANI_vs_ACC.pdf - Relationship between genome similarity (ANI) and gene-cluster conservation (**used for Figure S6**)
- Pairs-CG-minLin=1.txt - Network file listing recurrently colocalized core-gene pairs - Visualization of the conserved core-gene network (Pairs-CG-minLin=1.txt_3.pdf in cytoscape **used for Figure 3**).

## 7 - EvolutionaryNetwork - *Evolutionary network and lineage-level phylogenetic discordance*
This workflow quantifies phylogenetic discordance among core-gene trees at the lineage level and reconstructs an evolutionary network describing alternative phylogenetic histories across the Methylobacterium core genome. Gene-tree topologies are summarized as lineage partitions, their frequencies are quantified across the core genome, and enrichment of evolutionary network edges in KEGG functions is assessed using permutation tests.

#### A) Monophyly of groups, lineages, and species across core-gene trees
For each protein-coding core gene: evaluate monophyly of groups A, B, C and D; evaluate monophyly of all lineages; calculate the proportion of monophyletic species; generate annotated gene-tree visualizations; generate species and lineage tree visualizations.
#### B) – Reconstruction of the evolutionary network
Gene trees are simplified at the lineage level and decomposed into monophyletic partitions. Consensus and alternative partitions are identified and their frequencies are quantified across the core genome. Network files are generated for pre-visualization in Cytoscape.
#### C) Statistics on partition categories
Partition frequencies are summarized for different categories of evolutionary events, including consensus partitions, lineage-specific partitions, and alternative partitions. Additional ANOVA, Tukey and proportion tests are performed to compare partition categories.
#### D) KEGG enrichment of evolutionary network partitions
KEGG pathways associated with core genes are mapped onto lineage partitions. Permutation tests are performed to identify partitions enriched in specific functional categories.
#### E) – Generation of annotated Cytoscape networks
Significant KEGG enrichments are integrated into the evolutionary network and new Cytoscape files are generated for visualization and downstream analyses.

#### Input files (Non-listed input files are generated by previous workflows.)
- nodes.txt : Node definitions used to build the evolutionary network.
- consensus.txt : Consensus lineage partitions used to classify network edges.

#### Output files
- Species_tree.pdf : Species tree colored by lineage.
- Lineage_tree.pdf : Lineage tree colored by lineage (**used for Figure S4A**).
- PlotTrees/CG_xxx-tree.pdf: Annotated phylogeny for each core gene.
- Monophyly+BranchLength-per-PEG.txt:  Monophyly statistics and branch-length summaries for all core genes.
- Other-Partitions_ALT_Sept2026.txt: Partitions not represented in the final network.
- Stat-network.txt: Summary statistics for partition categories (**used for Figure 2C, Table S6**).
- Statistics-Networks.pdf: Distribution of partition frequencies across network categories, Tukey test and proportion tests.
- PartKEGG_K=100000.txt: Combined permutation results across all replicates.
- Sum_PartKEGG_K=100000.txt: KEGG pathways significantly enriched in network partitions (**used for Figure 2B, Table S6**).
- Partitions-for-Cyto_Sept2026-New.txt: All partition (**Used for Table S6**)
- Partitions-for-Cyto_Sept2026-New2.txt: All partition (formatted for Cytoscape)
- Partitions-for-Cyto_Sept2026-ALT_New.txt: Annotated network containing only alternative partitions - input for Cytoscape (**Used for Figure 2A**)
- Partitions-for-Cyto_Sept2026-CONS_New.txt: Annotated network containing only consensus partitions  (backbone of the network) - input for Cytoscape (**Used for Figure 2A**)
- Merged Network.pdf: output from cytoscape, combining above partition files (**Used for Figure 2A**)

## 8 - BAPS - *Hierarchical population structure and phylogenetic incongruence*
This workflow uses rhierBAPS to identify population structure within each core-gene alignment and compare BAPS-defined populations with established Methylobacterium lineages. The resulting classifications are used to estimate the proportion of phylogenetically incongruent populations across genes, visualize population structure across the core genome, and evaluate the relationship between BAPS incongruence and HGT likelihood (RES).

#### A) Selection of the number of BAPS populations
Core-gene alignments are analyzed across a range of predefined population numbers (n.pops = 2–60) to evaluate the number of BAPS lineages recovered for each gene. Additional statistics quantify the proportion of conserved BAPS lineages and species across different values of n.pops.
#### B) Classification of strains using rhierBAPS
Each core-gene alignment is analyzed with rhierBAPS using the number of Methylobacterium lineages as the target number of populations (nLIN = 24). BAPS populations are compared with species-tree lineages and classified according to their phylogenetic consistency.
#### C) Visualization of BAPS classifications across the core genome
BAPS classifications are summarized across all core genes and displayed in a circular representation showing genome position, HGT likelihood (RES), and lineage assignments. Special emphasis is given to core-gene clusters of interest (KDP, BCHL1 and BCHL2).
#### D) Distribution of BAPS classifications among species
For each species, the proportion of core genes assigned to different BAPS categories (conserved, within-group incongruence, among-group incongruence, low-resolution, missing) is calculated and visualized.
#### E) Distribution of BAPS classifications among Methylobacterium lineages
BAPS classifications are summarized at the lineage level to quantify variation in phylogenetic incongruence among lineages.
#### F) Relationship between BAPS incongruence and HGT likelihood
For each core gene, the proportion of phylogenetically incongruent BAPS populations is compared with its HGT likelihood (RES) using linear and logarithmic regression models.

#### Input files (Non-listed input files are generated by previous workflows.)
- CG_FastGear.txt: Positions and annotations of KDP and BCHL core-gene clusters in the reference genome

#### Output files
- CG_xxx-rhierbaps_ALL.txt: BAPS classifications for each tested value of n.pops.
- Npop_Combined.pdf: Number of recovered BAPS lineages as a function of the tested number of populations (**used for Figure S10A**).
- Npop_Combined_VATstat.pdf: Proportion of conserved BAPS lineages and species across different values of n.pops (**used for Figure S10B**).
- CG_xxx-rhierbaps.txt: rhierBAPS classification for each core gene using nLIN = 24.
- SUM-CG_xxx.txt: Summary statistics and BAPS classification categories for each core gene.
- Colors-rhierbaps-species.txt: Species-level color matrix summarizing BAPS classifications across all core genes.
- Circular_BAPS_nLIN=24.pdf: Circular overview of BAPS classifications, lineage assignments, core-gene positions, and RES values (**used for Figure 4B-C**).
- proportion_species_BAPS_nLIN=24.pdf: Distribution of BAPS classification categories among species (**used for Figure 4D**).
- BAPSlin_mean-per-METHlin.pdf: Average distribution of BAPS classification categories among Methylobacterium lineages (**used for Figure S10C**).
- BAPSlin_mean-per-METHlin.txt: Numerical summary corresponding to the previous figure.
- HGT_rhierbaps_vs_RF_linear.pdf: Relationship between the proportion of incongruent BAPS lineages and HGT likelihood (RES) (**used for Figure 4A**).
- HGT_rhierbaps_vs_RF_log.pdf: Alternative logarithmic model of the same relationship.

## 9 - BCHL_Hyphomicrobiales - *Evolution and organization of bacteriochlorophyll and carotenoid biosynthesis gene clusters*
This workflow investigates the evolution, phylogeny, and genomic organization of bacteriochlorophyll and carotenoid biosynthesis (BCHL) genes across Methylobacterium and closely related Hyphomicrobiales. BCHL gene occurrence, cluster architecture, gene phylogenies, and genomic context are compared across species, lineages, and genera to identify patterns consistent with horizontal gene transfer and cluster rearrangement.

**Genome dataset used in this section are from a previous study and can be retrieved here: https://github.com/JBLED/Hyphomicrobiales-phyllosphere**

#### A) BCHL gene occurrence across Hyphomicrobiales
Occurrences of the 32 BCHL core genes are quantified across Methylobacterium, related Hyphomicrobiales genomes, and MAGs. Gene prevalence is summarized at the genus level and compared with HGT likelihood (RES values).
#### B) Identification of BCHL and other outlier gene clusters
Core genes associated with BCHL biosynthesis and other highly discordant functional clusters are identified from core-gene annotations and RF residual values. These gene sets are used in subsequent phylogenetic and cluster-structure analyses.
#### C) BCHL cluster architecture
BCHL gene copy number and genomic organization are characterized across genomes. Genes located on the same scaffold and within 25 kb are grouped into subclusters following the criteria used in the conserved-cluster analyses. Cluster architecture is summarized for both Methylobacterium and related Hyphomicrobiales.
#### D) BCHL phylogenies
All copies of BCHL-associated genes are extracted from genome annotations, aligned, filtered, and individually analyzed phylogenetically using RAxML. Because BCHL copy-number variation is common outside Methylobacterium, all detected copies are retained for analysis.
#### E) BCHL phylogeny visualization
Individual BCHL gene trees are reformatted and visualized to compare gene histories across Methylobacterium lineages and related Hyphomicrobiales taxa.
#### F) Comparative BCHL cluster organization
Representative BCHL clusters are extracted from selected genomes and MAGs, re-annotated, and prepared for comparative gene-cluster alignment with clinker. Cluster architecture is summarized across Methylobacterium lineages and Hyphomicrobiales genera.
#### G) Combined Hyphomicrobiales phylogeny
A consensus phylogeny combining Methylobacterium lineages and related Hyphomicrobiales families is generated to place BCHL cluster evolution in a broader evolutionary context.

#### Input files (Non-listed input files are generated by previous workflows.)
- Hyphomicrobiales-Genus-Order.txt: Reference order of Hyphomicrobiales genera.
- Lichenihabitantaceae-AsterLinTree-nov-2025_mBS=50.tre: Consensus Hyphomicrobiales lineage tree.

#### Output files
- HGTlikeliness-distribution.pdf: Distribution of HGT likelihood values among BCHL genes.
- BCHLhypOcc.pdf: Heatmap showing BCHL gene occurrence across genera.
- BCHLhypOcc.txt: BCHL gene occurrence matrix (**used for used to generate Figure S14 top panels**).
- Outliers-for-Phylo-Hyphomicrobiales.txt: List of BCHL and other high-discordance gene clusters used in downstream analyses.
- BCHL-gene-occurence-per-genome.txt: BCHL gene copy numbers across genomes (**used for Table S9**).
- BCHL-arch-Hyphomicr.txt: BCHL cluster architecture summary across genomes.
- BchClusters_4&8-per-spe.pdf (**used for Figure S11A**)
- BchClusters_4&8-per-lin.pdf (**used for Figure S11B**)
- Summary-per-outlier-core-gene-Hyphomycrobiales.txt: Alignment and sequence statistics for BCHL-associated genes.
- CG-Hypho_*.fas: Final alignments submitted to RAxML.
- RAxML_tree-BCHL1-CG_799_reformated: bchH phylogeny (**used for Figure 5D**).
- RAxML_tree-BCHL1-CG_223_reformated: bchG phylogeny (**used for Figure 5E**).
- RAxML_tree-BCHL2-CG_743_reformated: pufC phylogeny (**used for Figure S14A**).
- RAxML_tree-BCHL2-CG_745_reformated: pufL phylogeny (**used for Figure S14B**).
- RAxML_tree-BCHL2-CG_746_reformated: pufM phylogeny (**used for Figure S14C**).
- Plot-Phylo-Outliers-Hyphomicrobiales.pdf: Overview of BCHL and outlier-gene phylogenies.
- BCHL-cluster-sequences-Hyphomicrobiales.txt: Coordinates and sequences of extracted BCHL subclusters.
- BCHL-cluster-sequences-Hyphomicrobiales.fas: FASTA file submitted for re-annotation.
- For_Clinker/*.gbk: GenBank files prepared for clinker analyses.
- Output-W678X869S241X23.pdf: Clinker comparison of BCHL clusters across Methylobacterium lineages (**used for Figure 5C, Figure S12B**).
- Output-N467B133V970G77.pdf: Clinker comparison of BCHL clusters in Hyphomicrobiales (**used for Figure 5C, Figure S13**).
- HYPHOMICROBIALES+METHYLOBACERIUM-Manually-Simplified-Rooted.tre: Combined Hyphomicrobiales–Methylobacterium lineage phylogeny.
- Hyphomicrobiales+Methylob-rooted_FigTree.pdf: Final rooted phylogeny used to display BCHL evolution (**used for Figure 5A**).
- Methylobacterium-AsterSpeTree-PEG-SPE-sept-2025_mBS=10-FIGTREE-Lineages.pdf: Species tree annotated with lineage assignments (**used for Figure S12A**).

## Files used to produce figures and tables
#### - Figure 1
- A) Methylobacterium-AsterSpeTree-PEG-SPE-sept-2025.pdf (Methylobacterium-AsterSpeTree-PEG-SPE-sept-2025_mBS=10.tre)
- B) Average-ANI-per-species.pdf & Average-cluster-cons-per-species.pdf
#### - Figure 2
- A) Merged Network.pdf (need to identifiy text files used as input in cytoscape)
- B&C) Sum_PartKEGG_K=1e+05.txt & Stat-network.txt
#### - Figure 3
- Pairs-CG-minLin=1.txt_3.pdf (Pairs-CG-minLin=1.txt) 
#### - Figure 4
- A) HGT_rhierbaps_vs_RF_linear.pdf;
- B) Circular_BAPS_nLIN=24.pdf (externe ring); C) (intern ring)
- D) proportion_species_BAPS_nLIN=24.pdf
#### - Figure 5
- A) Hyphomicrobiales+Methylob-rooted_FigTree.pdf (Hyphomicrobiales+Methylob-rooted.tre combining Lichenihabitantaceae-AsterLinTree-nov-2025_mBS=50.tre rooted on Rhabdaerophylum and Methylobacterium-Species-Manually-Simplified-Rooted.tre)
- B) BCHL-cluster-occurence.xlsx
- C) Modified from Clinker raw outputs: Output-W678X869S241X23.pdf (Methylobacterium) +  and Output-N467B133V970G77.pdf (Hyphomicrobiales)
- D) bchH: RAxML_tree-BCHL1-CG_799_reformated_FigTree_Collapsed.pdf
- E) bchG: RAxML_tree-BCHL1-CG_223_reformated_FigTree_Collapsed.pdf
#### - Figure S1
- Isolate_xxx-Megahit-SRxx.pdf
#### - Figure S2
- CoreGenome_N50.pdf
#### - Figure S3
- DistANI.pdf
#### - Figure S4
- Lineage_tree.pdf (Methylobacterium-AsterLinTree-PEG-sept-2025_mBS=10.tre)
#### - Figure S5
- A-B: PairwiseDistRF-KEGG+Segmented.pdf (CBMB20 genome)
- C: Pairwise_RF-vs_colocalization-raw ; Pairwise_RF-vs_colocalization-bins.pdf
#### - Figure S6
- ANI_vs_ACC.pdf
#### - Figure S7
- RFdist-Correction.pdf
#### - Figure S8
- KEGG-enrichment.pdf
#### - Figure S9
- BCHL_diagram.png (no script needed)
#### - Figure S10
- A) Npop_Combined.pdf (BAPS lineage reach a plateau)
- B) Npop_Combined_VATstat.pdf (max proportion of lineages arround n.pops =24)
- C) BAPSlin_mean-per-METHlin.pdf
#### - Figure S11
- A) BchClusters_4&8-per-spe.pdf
- B) BchClusters_4&8-per-lin.pdf
#### - Figure S12
- A) Lineage tree : Methylobacterium-AsterSpeTree-PEG-SPE-sept-2025_mBS=10-FIGTREE-Lineages.pdf
- B) Clinker raw output: Output-W678X869S241X23.pdf
#### - Figure S13
- Output-N467B133V970G77
#### - Figure S14
- A) RAxML_tree-BCHL2-CG_743_reformated_figtree_pufC.pdf
- B) RAxML_tree-BCHL2-CG_745_reformated_figtree_pufL.pdf
- C) RAxML_tree-BCHL2-CG_746_reformated_figtree_pufM.pdf
- and BCHLhypOcc.xlsx for A,B and C top panels
#### - Table S1
- Data.txt
- Stat-genomes.txt
- Methylobacterium-core-gene-occ.txt
#### - Table S2
- List-final-assemblies.txt (and metadata for NCBI submissions)
#### - Table S3
- Core-gene-code.txt
- Summary-per-core-gene-after-alignment.txt
- CoreGenes_nPOL_dNdS.txt
- GeneTreeStat.txt
- Summary-KEGG-RF-Pos-per-PEG.txt
#### - Table S4
- Stat-per-species.txt
- Stat-per-lineage.txt
- Taxonomy.txt
- Mean_ANI_ACC_among_lineage.txt
#### - Table S5
- Perm-RF-Summary-Per-KEGG.txt
#### - Table S6
- Partitions-for-Cyto_Sept2026-New2.txt
- Stat-network.txt (and missing one file with tukey and proportion test)
- Sum_PartKEGG_K=1e+05.txt
#### - Table S7
- ANOVA-RFdist-Correction.txt
#### - Table S8
- clustering-genes-lineages-pVal.txt
#### - Table S9
- BCHL-gene-occurence-per-genome.txt
