# HybTree-Pipeline
Here we introduce a pipeline that is designed to construct phylogenetic tree directly by Whole Genome Sequencing reads. Hybtree works by skipping the genome assembly and annotation step, and directly generating high-accuracy species-level trees using reference gene assemblies from across the genome.
# Description of the Pipeline
The pipeline involves pre-processing sequencing reads to ensure high-quality data for analysis. Raw sequencing reads are analyzed using FastQC and poor quality bases are trimmed using Trimmomatic v0.36 to ensure only high-quality reads are utilized. Hybpiper performs reference-based assembly of pre-processed reads, using specific aligners for nucleotide or amino acid sequences. The most suitable reference sequence for each gene and sample is identified based on alignment scores. Reads are sorted into directories corresponding to each gene in the reference set, and SPAdes reconstructs the reference gene sequence from the sorted reads. Extracted coding sequences are translated into amino acids using biopython, amino acid sequences are aligned using mafft v7.505 and then converted back to nucleotide sequences. The aligned nucleotide sequences are partitioned by genes and concatenated into a phylogenetic supermatrix using the catsequences tool. IQ-TREE v1.6.12 constructs a phylogenetic tree from the concatenated alignment, using the ModelFinder feature to identify the optimal substitution model for the data. The tree is then bootstrapped using UFBoot2 to determine the support for the tree's branches. Tree files can be then visualized in ETE toolkit or iTOL.

# Dependencies
1. Snakemake :
https://snakemake.readthedocs.io/en/stable/getting_started/installation.html

2. FastQC :
https://github.com/s-andrews/FastQC

3. Trimmomatic : 
http://www.usadellab.org/cms/?page=trimmomatic

4. Hybpiper :
https://github.com/mossmatters/HybPiper/wiki/Installation

5. Biopython :
https://biopython.org/wiki/Download

6. MAFFT :
https://mafft.cbrc.jp/alignment/software/

7. IQ-TREE :
http://www.iqtree.org/

# Pipeline Inputs
1. High Throughput DNA sequencing reads
Raw sequencing reads in FASTQ format. The cleaning and removal of contaminant sequences will be carried out by Hybtree pipeline.

2. Reference genes
Along with the reads, users also require a set of reference genes that cab be either nucleotide or amino acid sequences in FASTA format. These reference sequences can be obtained from OrthoDB, OMA browser or custom set of genes expected to be present across all the samples used for constructing phylogenetic tree. The reference sequence file can contain multiple sequences for each locus using sequence IDs. Each sequence ID should include the organism name and the gene identifier, separated by a hyphen.

example of a reference genes file :
"""
>Arabidopsis-atpB
MATNRAIRLSSMIFILGAFLGAGAVVGSAQGKFVEKIARQPELQGKIRGTLTNSLAFLMG
>Zea-atpB
MASTKAIRLSAIVFLGAFLGAGAVVGAAQGKFVEKIARQPELQGKIRGTLTNSLAFLMG
>Arabidopsis-rbcS
MASTKEILKAYPETMRRFLVGYQGCTEEYSSLRDGRDLIGNDTGAASTERVLAKYGPRPL
>Zea-rbcS
MASTKEVLKAYPETMRRFLVGYQGCTEEYSSLRDGRDLVGNDTGAASTERVLAKYGPRPL
"""
