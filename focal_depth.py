from collections import defaultdict
import argparse
import sys
import pysam
import numpy
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

hg38_chroms = {
  "NC_000001.11": ("1", 248956422),
  "NC_000002.12": ("2", 242193529),
  "NC_000003.12": ("3", 198295559),
  "NC_000004.12": ("4", 190214555),
  "NC_000005.10": ("5", 181538259),
  "NC_000006.12": ("6", 170805979),
  "NC_000007.14": ("7", 159345973),
  "NC_000008.11": ("8", 145138636),
  "NC_000009.12": ("9", 138394717),
  "NC_000010.11": ("10", 133797422),
  "NC_000011.10": ("11", 135086622),
  "NC_000012.12": ("12", 133275309),
  "NC_000013.11": ("13", 114364328),
  "NC_000014.9": ("14", 107043718),
  "NC_000015.10": ("15", 101991189),
  "NC_000016.10": ("16", 90338345),
  "NC_000017.11": ("17", 83257441),
  "NC_000018.10": ("18", 80373285),
  "NC_000019.10": ("19", 58617616),
  "NC_000020.11": ("20", 64444167),
  "NC_000021.9": ("21", 46709983),
  "NC_000022.11": ("22", 50818468),
  "NC_000023.11": ("X", 156040895),
  "NC_000024.10": ("Y", 57227415),
  "NC_012920.1": ("MT", 16569)
}
t2t_chroms = {
  "NC_060925.1": ("1", 248387328),
  "NC_060926.1": ("2", 242696752),
  "NC_060927.1": ("3", 201105948),
  "NC_060928.1": ("4", 193574945),
  "NC_060929.1": ("5", 182045439),
  "NC_060930.1": ("6", 172126628),
  "NC_060931.1": ("7", 160567428),
  "NC_060932.1": ("8", 146259331),
  "NC_060933.1": ("9", 150617247),
  "NC_060934.1": ("10", 134758134),
  "NC_060935.1": ("11", 135127769),
  "NC_060936.1": ("12", 133324548),
  "NC_060937.1": ("13", 113566686),
  "NC_060938.1": ("14", 101161492),
  "NC_060939.1": ("15", 99753195),
  "NC_060940.1": ("16", 96330374),
  "NC_060941.1": ("17", 84276897),
  "NC_060942.1": ("18", 80542538),
  "NC_060943.1": ("19", 61707364),
  "NC_060944.1": ("20", 66210255),
  "NC_060945.1": ("21", 45090682),
  "NC_060946.1": ("22", 51324926),
  "NC_060947.1": ("X", 154259566),
  "NC_060948.1": ("Y", 62460029),
}
chrom_sets = [
    [str(i+1) for i in range(22)] + ["X", "Y"],
    ["chr"+str(i+1) for i in range(22)] + ["chrX", "chrY"],
    [c for c in hg38_chroms],
    [c for c in t2t_chroms]
]
chr_idx = {}
for chroms in chrom_sets:
    for c in range(len(chroms)):
        chr_idx[chroms[c]] = c

def get_gene_annotation(gene, gff):
    gene_id = None
    gene_name = None
    chrom = None
    gene_strand = None
    mrna_id = None
    mrna_id = None
    cds = []
    for line in open(gff):
        g = line.strip().split('\t')
        if g[0][0] == '#':
            continue
        if gene_id is None and g[2] == "gene":
            fields = {a[:a.index('=')]:a[a.index('=')+1:] for a in g[8].split(';')}
            if (("ID" in fields and gene.lower() == fields["ID"].lower()) or
               ("Name" in fields and gene.lower() == fields["Name"].lower()) or
               ("gene_synonym" in fields and gene.lower() in [s.lower() for s in fields["gene_synonym"].split(',')])):
                gene_id = fields["ID"]
                gene_name = fields["Name"]
            if gene_id is not None:
                gene_strand = g[6]
        elif gene_id is not None and mrna_id is None and g[2] == "mRNA":
            fields = {a[:a.index('=')]:a[a.index('=')+1:] for a in g[8].split(';')}
            if fields["Parent"] == gene_id:
                mrna_id = fields["ID"]
                mrna_name = fields["Name"]
                chrom = g[0]
        elif mrna_id is not None:
            if g[2] == "gene":
                break
            elif g[2] == "CDS":
                fields = {a[:a.index('=')]:a[a.index('=')+1:] for a in g[8].split(';')}
                if fields["Parent"] == mrna_id:
                    cds.append((int(g[3]), int(g[4]))) # 1-indexed start and end, both inclusive
    return (gene_id, gene_name, mrna_id, mrna_name, chrom, gene_strand, cds)

def main(af, chrom, st, en, binsize, outfig, ymax=None):
    c = chrom

    # convert chrom to appropriate accession/name
    found = False
    for r in af.references:
        for chroms in chrom_sets:
            if r in chroms:
                chrom = chroms[chrom_sets[0].index(chrom)]
                found = True
                break
        if found:
            break
    if not found:
        sys.stderr.write("ERROR: {chrom} could not be converted to BAM references using hg38 or T2T names: {af.references}\n")
        sys.exit(1)

    alns = []
    fusions = defaultdict(list)
    ct = 0
    covg = numpy.zeros((en-st+1) if binsize == 0 else (en//binsize-st//binsize+1), dtype='u4')
    plt.figure(figsize=(10,8), dpi=300)
    plt.rc('font', size=14)
    plt.rc('axes', titlesize=24)
    plt.rc('axes', labelsize=14)
    plt.rc('ytick', labelsize=14)
    plt.rc('xtick', labelsize=14)
    for a in af.fetch(chrom, st, en):
        if a.is_secondary: # use only primary alignments
            continue
        if binsize == 0:
          covg[max(a.reference_start,st)-st:min(a.reference_end,en)-st+1] += 1
        else:
          mid = a.reference_start + (a.reference_end-a.reference_start)//2
          if mid < en and mid > st:
            covg[mid//binsize - st//binsize] += 1

    if binsize == 0:
      plt.scatter(range(st,en+1), covg, s=3)
    else:
      #plt.plot(range(st//binsize * binsize, (en//binsize+1)*binsize, binsize), covg)
      plt.scatter(range(st//binsize * binsize, (en//binsize+1)*binsize, binsize), covg, s=3)
    plt.xlim((st,en))
    plt.ylim((0,plt.ylim()[1] if ymax is None else ymax))
    val = binsize
    unit = 1
    while val >= 1000:
        val /= 1000
        unit *= 1000
    if int(val) * unit == binsize:
        val = int(val)
    unit = {1: "bp", 1000: "Kbp", 1000000: "Mbp", 1000000000: "Gbp"}[unit]
    plt.ylabel("Nucleotide sequencing depth" if binsize == 0 else f"Reads / {val} {unit}")
    plt.xlabel(f"Chromosome {c}")
    ticks, labels = plt.xticks()
    diff = ticks[1] - ticks[0]
    div = 1
    while diff/div > 1000:
        div *= 1000
    unit = {1: "bp", 1000: "Kbp", 1000000: "Mbp", 1000000000: "Gbp"}[div]
    plt.xticks(ticks, [f"{(t/div) if t//div * div != t else int(t//div)} {unit}" for t in ticks], rotation=45)
    plt.tight_layout()
    plt.savefig(outfig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Draw coverage from whole-genome nanopore seq BAM file")
    parser.add_argument("bam", help="BAM reads to ref (mm2 -x map-ont)")
    parser.add_argument("gff", help="Gene annotation (GFF) matching the appropriate reference, it should include gene, mRNA, and CDS boundaries")
    parser.add_argument("region", help="{chrom[:start-end] | gene[+margin]}")
    parser.add_argument("fig", help="Filename for plot/figure")
    parser.add_argument("--binsize", help="Bin width if binning read counts", type=int, default=0)
    parser.add_argument("--ymax", help="Y axis maximum", type=int, default=None)
    args = parser.parse_args()

    # parse region string
    r = args.region
    gene = None
    chrom = None
    start = None
    end = None
    if ':' in r:
        chrom = r[:r.index(':')]
        start, end = [int(a) for a in r[r.index(':')+1:].split('-')]
    elif '+' in r:
        gene, margin = r.split('+')
        margin = int(margin)
    else:
        for chroms in chrom_sets:
            if r in chroms:
                chrom = r
                break
        else:
            gene = r
            margin = 0

    if gene is not None:
        (gene_id, gene_name, mrna_id, mrna_name, chrom, gene_strand, cds) = get_gene_annotation(gene, args.gff)
        if gene_id is None:
            sys.stderr.write(f"'{gene}' is not a recognized chromosome or gene\n")
            sys.exit(1)
        start = min(c[0] for c in cds) - margin
        end = max(c[1] for c in cds) + margin
        sys.stderr.write(f"Plotting {gene} on {chrom}:{start}-{end}...\n")

    if chrom is not None:
        if start is None:
            start = 0
        if chrom in hg38_chroms:
            if end is None:
                end = hg38_chroms[chrom][1]
            chrom = hg38_chroms[chrom][0] # convert to 1..22,X,Y
        elif chrom in t2t_chroms:
            if end is None:
                end = t2t_chroms[chrom][1]
            chrom = t2t_chroms[chrom][0] # convert to 1..22,X,Y
        else:
            for chroms in chrom_sets:
                if chrom in chroms:
                    chrom = chrom_sets[0][chroms.index(chrom)] # convert to 1..22,X,Y
                    if end is None:
                        end = hg38_chroms[chrom_sets[2][chroms.index(chrom)]][1]
                    break

    main(af, chrom, start, end, args.binsize, args.fig, args.ymax)
