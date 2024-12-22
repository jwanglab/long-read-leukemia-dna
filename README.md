
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


Quick start
-----------

    python3 callmemaybe.py *sample.bam* *targets.bed* *repeats.bed* *output_prefix*

Input:
*  **sample.bam** - BAM file including all reads (enriched and others), aligned to the appropriate reference
*  **targets.bed** - BED file including relevant genes that should match the file used for adaptive sampling except *WITHOUT* the 50Kbp padding
*  **repeats.bed** - BED file including repetitive or poorly mappable elements, ex. RepeatMasker, cen/sat/tel regions, genomicSuperDups, etc.
*  **output\_prefix** - output path for coverage profile, enriched read alignments


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
