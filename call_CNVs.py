import sys
import json
import argparse
import pysam
from collections import defaultdict

min_anchor = 500
gene_del_margin = 50000
breakpoint_margin = 50

class Reference:
    def __init__(self, chrom_list):
        self.chroms = chrom_list
        self.chrom_names = ["chr"+str(i+1) for i in range(22)] + ["chrX", "chrY"]
        self.chr_idx = {self.chroms[c]:c for c in range(len(self.chroms))}
        for c in range(len(self.chrom_names)):
            self.chr_idx[self.chrom_names[c]] = c

t2t = Reference(["NC_060925.1", "NC_060926.1", "NC_060927.1", "NC_060928.1", "NC_060929.1", "NC_060930.1", "NC_060931.1", "NC_060932.1", "NC_060933.1", "NC_060934.1", "NC_060935.1", "NC_060936.1", "NC_060937.1", "NC_060938.1", "NC_060939.1", "NC_060940.1", "NC_060941.1", "NC_060942.1", "NC_060943.1", "NC_060944.1", "NC_060945.1", "NC_060946.1", "NC_060947.1", "NC_060948.1"])

def main(bam_file, bed_genes, covg_tsv, fusion_json):

    # parse BED file
    data = [line.strip().split('\t') for line in open(bed_genes)]
    targets = {d[3]: (t2t.chr_idx[d[0]], int(d[1]), int(d[2])) for d in data}

    # read focal on-target sequencing depth
    f = json.loads(open(fusion_json).read())

    # read regional read depth bins
    # chromosome	bin_start	bin_end	n_reads
    # chr1	0	1000000	2183
    covg = [line.strip().split('\t') for line in open(covg_tsv)][1:]

    result = {}
    for g in ["RUNX1", "CDKN2A", "CDKN2B", "IKZF1"]:
        (focal, local) = depth(g, targets, f, covg)
        result[g] = {"focal": focal, "local": local}

    chrom, ERG_dels = internal_del(bam_file, "ERG", targets)
    result["ERG"] = {}
    result["ERG"]["deletions"] = []
    for d in ERG_dels:
        sys.stderr.write(f"ERG deletion {d[0]} - {d[1]} ({d[1]-d[0]} nt): {ERG_dels[d]} reads\n")
        if ERG_dels[d] >= 3:
            result["ERG"]["deletions"].append({"chrom": t2t.chrom_names[chrom], "start": d[0], "end": d[1], "reads": ERG_dels[d]})

    chrom, IKZF1_dels = internal_del(bam_file, "IKZF1", targets)
    result["IKZF1"] = {}
    result["IKZF1"]["deletions"] = []
    for d in IKZF1_dels:
        sys.stderr.write(f"IKZF1 deletion {d[0]} - {d[1]}: {IKZF1_dels[d]} reads\n")
        if IKZF1_dels[d] >= 3:
            result["IKZF1"]["deletions"].append({"chrom": t2t.chrom_names[chrom], "start": d[0], "end": d[1], "reads": IKZF1_dels[d]})

    chrom, PAX5_dels = internal_del(bam_file, "PAX5", targets)
    result["PAX5"] = {}
    result["PAX5"]["deletions"] = []
    for d in PAX5_dels:
        sys.stderr.write(f"PAX5 deletion {d[0]} - {d[1]}: {PAX5_dels[d]} reads\n")
        if PAX5_dels[d] >= 3:
            result["PAX5"]["deletions"].append({"chrom": t2t.chrom_names[chrom], "start": d[0], "end": d[1], "reads": PAX5_dels[d]})

    print(json.dumps(result))

# returns CN by (focal sequencing depth, local 1Mbp bins)
def depth(gene, targets, f, covg, exclude=["RUNX1", "CDKN2A"]):
    focal_depth = f["qc"]["gene_coverage"][gene]
    all_gene_depths = sorted([f["qc"]["gene_coverage"][g] for g in f["qc"]["gene_coverage"] if g not in exclude and "DUX4" not in g])
    median_other = all_gene_depths[len(all_gene_depths)//2]
    sys.stderr.write(f"Est. {gene} focal CN (assuming broadly diploid): {focal_depth:.2f}/{median_other:.2f} ({focal_depth/median_other*2:.2f}x)\n")

    bins = [int(c[3]) for c in covg if t2t.chr_idx[c[0]] == targets[gene][0] and int(c[1]) < targets[gene][2] and int(c[2]) > targets[gene][1]] # should be at most 2 bins unless the gene is >>1Mbp
    bin_depth = sum(bins)/len(bins)
    all_gene_bins = [[int(c[3]) for c in covg if t2t.chr_idx[c[0]] == targets[g][0] and int(c[1]) < targets[g][2] and int(c[2]) > targets[g][1]] for g in targets if g not in exclude and "DUX4" not in g]
    # some genes fall in bins that are filtered out for being broadly repetitive (BCL9, TCRA, IG*, PAR1 genes, etc) - ignore them
    sorted_bins = sorted([sum(b)/len(b) for b in all_gene_bins if len(b) > 0])
    median_bin = sorted_bins[len(sorted_bins)//2]
    sys.stderr.write(f"Est. {gene} bin CN (assuming broadly diploid): {bin_depth}/{median_bin} ({bin_depth/median_bin*2:.2f}x)\n")

    return (focal_depth/median_other*2, bin_depth/median_bin*2)

def internal_del(bam_file, gene, targets):
    try:
        af = pysam.AlignmentFile(bam_file, mode='rb', ignore_truncation=True)
    except FileNotFoundError as e:
        sys.stderr.write(f"Error: File not found: {bam_file}\n")
        sys.exit(1)

    hits = defaultdict(list)

    chrom, st, en = targets[gene]
    if t2t.chroms[chrom] in af.references:
        for a in af.fetch(t2t.chroms[chrom], st-gene_del_margin, en+gene_del_margin):
            _, qs, qe, ts, te = fix_sam_coords(a)
            hits[a.query_name].append((a, (qs, qe, ts, te)))
    elif t2t.chrom_names[chrom] in af.references:
        for a in af.fetch(t2t.chrom_names[chrom], st-gene_del_margin, en+gene_del_margin):
            _, qs, qe, ts, te = fix_sam_coords(a)
            hits[a.query_name].append((a, (qs, qe, ts, te)))
    else:
        sys.stderr.write("ERROR: unknown reference IDs in bam file '{bam_file}'")
        sys.exit(1)

    # filter duplex from hits list
    to_remove = []
    for r in hits:
        if ';' in r: # a duplex read
            to_remove.extend(r.split(';'))
    for r in to_remove:
        if r in hits:
            del hits[r]

    dels = {}
    for read in hits:
        for i in range(len(hits[read])):
            h0, (qs0, qe0, ts0, te0) = hits[read][i]
            for j in range(i+1, len(hits[read])):
                h1, (qs1, qe1, ts1, te1) = hits[read][j]

                # file must have been sorted, and alignments in order, so h0 must be the left and h1 the right side of a putative deletion
                '''
                if te0 - ts1 > -1000: # they are too close in the target (ERG deletions are typically ~47Kbp)
                    continue
                '''

                if h0.is_reverse != h1.is_reverse or qs0 < qs1 == h0.is_reverse: # not oriented properly
                    continue

                # the read orientation is arbitrary however
                if min(qe0, qe1) - max(qs0, qs1) > 100: # they significantly overlap in the read
                    continue

                if te0 - ts0 > min_anchor and te1 - ts1 > min_anchor:
                    for d in dels:
                        if abs(d[0] - te0) < breakpoint_margin and abs(d[1] - ts1) < breakpoint_margin: # merge to existing
                            dels[d] += 1
                            break
                    else:
                        dels[(te0, ts1)] = 1
    return chrom, dels

# PAF and BAM files store query coordinates differently - BAM have to be converted back to the original strand coordinates before comparison
def fix_sam_coords(aln):
    qlen = aln.infer_read_length()
    # query_alignment_start IGNORES hard-clipped bases, so we have to add them back to the beginning by checking the cigar string
    hard_clipped = 0
    for i in range(len(aln.cigarstring)):
        if aln.cigarstring[i] in 'MIDNSHP=X':
            if aln.cigarstring[i] == 'H':
                hard_clipped = int(aln.cigarstring[0:i])
            break
    qs = aln.query_alignment_start + hard_clipped # 0-indexed
    qe = aln.query_alignment_end + hard_clipped
    ts = aln.reference_start
    te = aln.reference_end
    if aln.is_reverse:
        tmp = qlen - qe
        qe = qlen - qs
        qs = tmp
    return (aln, qs, qe, ts, te)


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Call focal CNVs based on regional depth, specific deletions, duplications")
    parser.add_argument("bam", help="BAM file, aligned to T2T")
    parser.add_argument("bed", help="BED file of adaptive target genes")
    parser.add_argument("covg", help="Binned read depth file (from lrdk)")
    parser.add_argument("fusion", help="JSON output from fusion caller (includes target depths)")
    args = parser.parse_args()
    main(args.bam, args.bed, args.covg, args.fusion)
