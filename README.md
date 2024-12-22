
Identification of tumorigenic drivers from nanopore DNA sequencing with adaptive sampling
=========================================================================================

Calls primary drivers with a focus on pediatric acute leukemia, including:

*  Translocations/fusion genes (ex. *ETV6*::*RUNX1*, *BCR*::*ABL1*)
*  Gross karyotype abnormalities (aneuploidy)
*  Focal copy-number variation (ex. *RUNX1*, *CDKN2A*, *ERG*)
*  Single-nucleotide variation and insertions/deletions (ex. *FLT3* internal tandem duplication [*FLT3*-ITD], *TPMT* mutation)


Installation
------------

Tested with Python 3.11

    git clone https://github.com/jwanglab/long-read-leukemia-dna
    cd long-read-leukemia-dna
    pip install -r requirements.txt
    chmod +x lrl


Quick start
-----------

To run the full pipeline, including fusion calling, karyotyping, SNV calling in relevant enriched genes, and targeted *FLT*-ITD detection on an existing, sorted BAM file:

    ./lrl <sorted.bam> <targets.bed> <reference.fasta> <repeats.bed> <annotation.gff> <output_prefix>

Input:
*  **sorted.bam** - BAM file including all reads (enriched and others), aligned to the appropriate reference and sorted
*  **reference.fasta** - Reference FASTA file, ex. [GRCh38](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GCF_000001405.26_GRCh38_genomic.fna.gz)
*  **targets.bed** - BED file including relevant genes that should match the file used for adaptive sampling except *WITHOUT* the 50Kbp padding
*  **repeats.bed** - BED file including repetitive or poorly mappable elements, ex. RepeatMasker, cen/sat/tel regions, genomicSuperDups, etc.
*  **annotation.gff** - GFF annotation from NCBI, ex. for [GRCh38](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GCF_000001405.26_GRCh38_genomic.gff.gz)
*  **output\_prefix** - output path for coverage profile, enriched read alignments

The output folder will contain the logs and results of the component analyses, as described below:


Fusion calling and karyotyping
------------------------------

**fusion\_karyotype.py** requires an optionally sorted BAM file, BED file with enriched targets (including putative translocation genes), and a BED file of annotated repetitive elements. The output log file contains sequencing and enrichment statistics including throughput and relative enrichment from adaptive sampling, sequencing depth across target genes and chromosomes, and detected translocations among the target genes including read support and genomic breakpoint(s). A figure is also produced showing relative genome-wide sequencing depth for karyotype assessment. If the input BAM file is sorted and the script is not run with the "--full" flag, it will only assess alignments in target regions, making sequencing statistics and karyotype unreliable, however translocation detection will be performed much faster.

    usage: Call karyotype and gene fusions/translocations from whole-genome nanopore seq BAM file [-h] [--read_lim READ_LIM] [--outdir OUTDIR] [--full] [--genes GENES [GENES ...]] bam bed repeats cov
    
    positional arguments:
      bam                   BAM of reads to hg38 (mm2 -x map-ont)
      bed                   BED file including target genes
      repeats               RepeatMasker BED file
      cov                   Output file for coverage figure
    
    options:
      -h, --help            show this help message and exit
      --read_lim READ_LIM   Maximum reads to process (default is all)
      --outdir OUTDIR       Output directory for bigger results
      --full                Assess full file, including off-target reads, in sorted BAM
      --genes GENES [GENES ...]
                        Evaluate only the given genes (must be in the enrichment set)

If "--outdir" is specific, raw per-nucleotide coverage profiles will be generated and stored as exported Numpy arrays in the specific directory. Also a PAF file will be created containing the alignments exclusively within the target regions. This file can be used to interrogate and visualize translocations using our fusion visualization interface:
[https://github.com/jwanglab/jwanglab.github.io/tree/master/fusion](https://github.com/jwanglab/jwanglab.github.io/tree/master/fusion), with a publicly available instance at [https://jwanglab.org/fusion/](https://jwanglab.org/fusion/)


SNV detection in enriched targets
---------------------------------

**adaptive\_SNVs.py** reports high-MAF (>0.3) SNVs in enriched genes from nanopore WGS with adaptive sampling. It is *not* a replacement for a comprehensive small variant detection pipeline but it designed to be fast, lightweight, and targeted to relevant genes in pediatric acute leukemia.

    usage: Identify coding variation in nanopore adaptive sampling BAM file [-h] [--genes GENES [GENES ...]] [--phase] bam ref gff

    positional arguments:
      bam                   BAM of reads to hg38 (mm2 -x map-ont)
      ref                   Reference FASTA file
      gff                   Matching reference annotation GFF

    options:
      -h, --help            show this help message and exit
      --genes GENES [GENES ...]
                            Gene name
      --phase               Perform variant phasing where possible


Utilities
---------

**focal\_depth.py** plots sequencing depth by total nucleotides or number of reads over a given region and is often useful for closer assessment of focal copy-number variation. Required annotation GFF file can be obtained from NCBI, ex. [GRCh38](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GCF_000001405.26_GRCh38_genomic.gff.gz).

    usage: Draw coverage from whole-genome nanopore seq BAM file [-h] [--binsize BINSIZE] [--ymax YMAX] bam gff region fig

    positional arguments:
      bam                BAM reads to ref (mm2 -x map-ont)
      gff                Gene annotation (GFF) matching the appropriate reference, it should include gene, mRNA, and CDS boundaries
      region             {chrom[:start-end] | gene[+margin]}
      fig                Filename for plot/figure

    options:
      -h, --help         show this help message and exit
      --binsize BINSIZE  Bin width if binning read counts
      --ymax YMAX        Y axis maximum

For example, to assess sequencing depth within ERG enrichment window (&plusmn; 50Kbp):

    python focal_depth.py <sample.sorted.bam> GCF_000001405.26_GRCh38_genomic.gff ERG+50000 ERG_50Kbp_margin.png


**flt3\_insilico\_pcr.py** identifies reads containing commonly used *FLT3* primer sites for capillary electrophoresis-based fragment size analysis and plots the distribution of these WGS-based "amplicon" sizes to enable identification of *FLT3*-ITD and allelic ratio.

    usage: Predict FLT3-ITD and allelic ratio from targeted nanopore sequencing [-h] fastq out

    positional arguments:
      fastq       Reads in FASTQ format
      out         Output image file for pseudo-capillary electrophoresis sizing of FLT3 amplicons

    options:
      -h, --help  show this help message and exit


Citation
--------

Julie Geyer, Kofi Opoku, John Lin, Lori Ramkissoon, Charles Mullighan, Nickhill Bhakta, Thomas B. Alexander, Jeremy R. Wang.
*Real-time genomic characterization of pediatric acute leukemia using adaptive sampling*.
[https://www.biorxiv.org/content/10.1101/2024.10.11.617690v1](https://www.biorxiv.org/content/10.1101/2024.10.11.617690v1) 
