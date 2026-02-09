
Identification of tumorigenic drivers from nanopore DNA sequencing with adaptive sampling
=========================================================================================

Calls primary drivers with a focus on pediatric acute leukemia, including:

*  Translocations/fusion genes (ex. *ETV6*::*RUNX1*, *BCR*::*ABL1*)
*  Gross karyotype abnormalities (aneuploidy)
*  Focal copy-number variation (ex. *RUNX1*, *CDKN2A*, *ERG*)
*  Single-nucleotide variation and insertions/deletions (ex. *FLT3* internal tandem duplication [*FLT3*-ITD], *TPMT* mutation)


Installation
------------

Tested with Python 3.11 and Rust 1.87.0 (other modern versions of Rust will probably work)

    git clone https://github.com/jwanglab/long-read-leukemia-dna
    cd long-read-leukemia-dna
    pip install -r requirements.txt
    cd maf_rust
    cargo build --release
    cd ..
    chmod +x lrl


Quick start
-----------

To run the full pipeline, including fusion calling, karyotyping, SNV calling in relevant enriched genes, and targeted *FLT*-ITD detection on an existing, sorted BAM file:

    ./lrl <sorted.bam> <target_genes.bed> <enriched_regions.bed> <reference.fasta> <repeats.bed> <annotation.gff> <output_prefix>

Input:
*  **sorted.bam** - BAM file including all reads (enriched and others), aligned to the appropriate reference and sorted
*  **reference.fasta** - Reference FASTA file, ex. [T2T-CHM13v2.0](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz)
*  **target_genes.bed** - BED file including relevant genes that should match the file used for adaptive sampling except *WITHOUT* the 50Kbp padding
*  **enriched_regions.bed** - BED file used for adaptive sampling
*  **repeats.bed** - BED file including repetitive or poorly mappable elements, ex. RepeatMasker, cen/sat/tel regions, genomicSuperDups, etc.
*  **annotation.gff** - GFF annotation from NCBI, ex. for [T2T-CHM13v2.0](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz)
*  **output\_prefix** - output path for coverage profile, enriched read alignments

The output folder will contain the logs and results of the component analyses, as described below:


Fusion calling
--------------

**fusions.py** requires an optionally sorted BAM file, BED file with enriched targets (including putative translocation genes), and a BED file of annotated repetitive elements. The output log file contains sequencing and enrichment statistics including throughput and relative enrichment from adaptive sampling, sequencing depth across target genes and chromosomes, and detected translocations among the target genes including read support and genomic breakpoint(s). If the input BAM file is sorted and the script is not run with the "--full" flag, it will only assess alignments in target regions, making sequencing statistics unreliable, however translocation detection will be performed much faster. Similarly, if the BAM file includes only on-target reads, the genome-wide coverage and enrichment will be incorrect.

    usage: Call gene fusions/translocations from whole-genome nanopore seq BAM file [-h] [--outdir OUTDIR] [--full] [--genes GENES [GENES ...]] [--one_sided ONE_SIDED] [--verbose] bed repeats bams [bams ...]

    positional arguments:
      bed                   BED file including target genes
      repeats               RepeatMasker BED file
      bams                  BAM file(s) of aligned reads (minimap2 -cx map-ont)

    options:
      -h, --help            show this help message and exit
      --outdir OUTDIR       Output directory for bigger results
      --full                Assess full file, including off-target reads, in sorted BAM
      --genes GENES [GENES ...]
                            Evaluate only the given genes (must be in the enrichment set)
      --one_sided ONE_SIDED
                            Given a file with a list of genes, searches for 'one-sided' fusions with these genes - will be slower
      --verbose             Print a huge amount of extra information


If "--outdir" is specified a PAF file will be created containing the alignments exclusively within the target regions. This file can be used to interrogate and visualize translocations using our fusion visualization interface:
[https://github.com/jwanglab/jwanglab.github.io/tree/master/fusion](https://github.com/jwanglab/jwanglab.github.io/tree/master/fusion), with a publicly available instance at [https://jwanglab.org/fusion/](https://jwanglab.org/fusion/)


Karyotyping
-----------

**lrdk_with_maf.py** requires a BAM file, BED file of annotated repetitive elements, BED file with enriched regions, and minor allele frequency data from 'nano_maf'.

    usage: Digital karyotyping (i.e. chromosome arm-level copy-number variation) from BAM or binned depth file [-h] [--repeats REPEATS] [--adaptive ADAPTIVE] [--maf MAF] [--outdir OUTDIR] [--read_lim READ_LIM] [--json] [-v] input [input ...]

    positional arguments:
      input                BAM/SAM of reads aligned to T2T or TSV of binned sequencing depth

    options:
      -h, --help           show this help message and exit
      --repeats REPEATS    RepeatMasker BED file
      --adaptive ADAPTIVE  Adaptive sequencing BED file
      --maf MAF            MAF file
      --outdir OUTDIR      Output directory/prefix
      --read_lim READ_LIM  Maximum reads to process (default is all)
      --json               Output JSON format (default: TSV)
      -v                   Verbose

In the specified output directory, a binned coverage profile (\<OUTDIR\>/coverage\_1Mbp\_bins.tsv) is produced, and a figure (\<OUTDIR\>/karytype.png) are produced showing relative genome-wide sequencing depth and minor allele frequency distribution for karyotype assessment. Output includes arm-level copy number estimates for ~metacentric chromosomes and whole-chromosome copy number for the others and an ISCN-compliant karyotype profile string.

Minor allele frequencies are calculated using the included nano_maf utility (Rust) and includes a filtered set of highly variable sites derived from 1000 Genomes data by default, on the T2T-CHM13v2.0 reference. See ./maf_rust/target/release/nano_maf or ./lrl for its usage.


SNV detection in enriched targets
---------------------------------

**targeted\_SNVs.py** reports high-MAF (>0.3) SNVs in enriched genes from nanopore WGS with adaptive sampling. It is *not* a replacement for a comprehensive small variant detection pipeline but it designed to be fast, lightweight, and targeted to relevant genes in pediatric acute leukemia.

    usage: Identify coding variation in nanopore adaptive sampling BAM file [-h] [--genes GENES [GENES ...]] [--maf MAF] [--phase] bam ref gff

    positional arguments:
      bam                   BAM of reads to T2T (mm2 -x map-ont)
      ref                   Reference FASTA file
      gff                   Matching reference annotation GFF

    options:
      -h, --help            show this help message and exit
      --genes GENES [GENES ...]
                            Gene name
      --maf MAF             Minimum MAF to report (default: 0.3)
      --phase               Perform variant phasing where possible


Internal tandem duplications
----------------------------

**itd.py** identifies intra-read insertions consistent with characteristic *FLT3* and *UBTF* internal tandem duplications (ITDs).

    usage: Predict FLT3 and UBTF ITD and allelic ratio from targeted nanopore sequencing [-h] bams [bams ...]

    positional arguments:
      bams        BAM of aligned reads (minimap2 -cx map-ont)

    options:
      -h, --help  show this help message and exit


Focal copy-number variation
---------------------------

**call_CNVs.py** calls relevant focal CNVs based on sequencing depth in enriched regions, local (~1Mbp scale) read depth, and/or detection of large intergenic deletions. It focuses on relevant variation in *CDKN2A*, *CDKN2B*, *ERG*, *IKZF1*, *PAX5*, and *RUNX1*, but is extensible to arbitrary regions.

    usage: Call focal CNVs based on regional depth, specific deletions, duplications [-h] bam bed covg fusion

    positional arguments:
      bam         BAM file, aligned to T2T
      bed         BED file of adaptive target genes
      covg        Binned read depth file (coverage_1Mbp_bins.tsv, from lrdk_with_maf.py)
      fusion      JSON output from fusion caller (includes target depths)

    options:
      -h, --help  show this help message and exit


Utilities
---------

**focal\_depth.py** plots sequencing depth by total nucleotides or number of reads over a given region and is often useful for closer assessment of focal copy-number variation. Required annotation GFF file can be obtained from NCBI, ex. [T2T-CHM13v2.0](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz).

    usage: Draw coverage from whole-genome nanopore seq BAM file [-h] [--binsize BINSIZE] [--ymax YMAX] [--anno ANNO [ANNO ...]] bam gff region fig

    positional arguments:
      bam                   BAM reads to ref (mm2 -x map-ont)
      gff                   Gene annotation (GFF) matching the appropriate reference, it should include gene, mRNA, and CDS boundaries
      region                {chrom[:start-end] | gene[+margin]}
      fig                   Filename for plot/figure

    options:
      -h, --help            show this help message and exit
      --binsize BINSIZE     Bin width if binning read counts
      --ymax YMAX           Y axis maximum
      --anno ANNO [ANNO ...]
                            Draw region annotation

For example, to assess sequencing depth within ERG enrichment window (&plusmn; 50Kbp):

    python focal_depth.py <sample.sorted.bam> GCF_000001405.26_GRCh38_genomic.gff ERG+50000 ERG_50Kbp_margin.png


Citation
--------

Julie Geyer, Kofi B. Opoku, John Lin, Lori Ramkissoon, Charles Mullighan, Nickhill Bhakta, Thomas B. Alexander, Jeremy R. Wang.
**Real-time genomic characterization of pediatric acute leukemia using adaptive sampling**. *Leukemia* 39, 1069-1077 (2025).
[https://www.nature.com/articles/s41375-025-02565-y](https://www.nature.com/articles/s41375-025-02565-y) 
