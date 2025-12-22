import sys
import json
import argparse
import pysam
import numpy
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

def main(bam_file, bed_genes, covg_tsv, fusion_json,karyo_json):

    # read focal on-target sequencing depth
    f = json.loads(open(fusion_json).read())

    # read karyotype info
    karyo = json.loads(open(karyo_json).read())
    # get 2n
    cn2s = []
    # get segments that are 2n
    cn2_segments = []
    for s in karyo:
        if not isinstance(karyo[s], str):
            if 'n' in karyo[s]:
                if karyo[s]['n']=='2':
                    cn2s.append(float(karyo[s]['x']))
                    cn2_segments.append(s)
    bin_cov_2n = numpy.mean(cn2s)

    sys.stderr.write(f"Diploid segment count: {len(cn2_segments)}\n")

    # read blast ratio
    if 'blast_ratio' in karyo:
        blast_ratio = float(karyo['blast_ratio'])
        # as a sanity check, cap at 1
        if blast_ratio > 1:
            blast_ratio = 1
    else:
        blast_ratio = None

    # parse BED file
    data = [line.strip().split('\t') for line in open(bed_genes)]
    targets = {d[3]: (t2t.chr_idx[d[0]], int(d[1]), int(d[2])) for d in data}
    diploid_gene = []
    for d in data:
        seg = get_segment([t2t.chrom_names[t2t.chr_idx[d[0]]], int(d[1]), int(d[2])])
        # only keep genes in diploid segments
        if seg is not None and seg in cn2_segments:
            diploid_gene.append(d[3])

    # read regional read depth bins
    # chromosome	bin_start	bin_end	n_reads
    # chr1	0	1000000	2183
    covg = [line.strip().split('\t') for line in open(covg_tsv)][1:]

    result = {"RUNX1":{}, "CDKN2A":{}, "CDKN2B":{}, "ERG":{}, "RUNX1":{}, "PAX5":{}, "KMT2A":{}}
    for g in ["RUNX1", "CDKN2A", "CDKN2B", "IKZF1"]:
        (focal, local) = depth(g, targets, f, diploid_gene, covg, bin_cov_2n, blast_ratio)
        result[g] = {"focal": focal, "local": local}

    chrom, ERG_dels, ERG_dups = internal_del_dup(bam_file, "ERG", targets)
    result["ERG"]["deletions"] = []
    for d in ERG_dels:
        sys.stderr.write(f"ERG deletion {d[0]} - {d[1]} ({d[1]-d[0]} nt): {ERG_dels[d]} reads\n")
        if ERG_dels[d] >= 3:
            result["ERG"]["deletions"].append({"chrom": t2t.chrom_names[chrom], "start": d[0], "end": d[1], "reads": ERG_dels[d]})

    chrom, IKZF1_dels, IKZF1_dups = internal_del_dup(bam_file, "IKZF1", targets)
    result["IKZF1"]["deletions"] = []
    for d in IKZF1_dels:
        sys.stderr.write(f"IKZF1 deletion {d[0]} - {d[1]}: {IKZF1_dels[d]} reads\n")
        if IKZF1_dels[d] >= 3:
            result["IKZF1"]["deletions"].append({"chrom": t2t.chrom_names[chrom], "start": d[0], "end": d[1], "reads": IKZF1_dels[d]})

    chrom, PAX5_dels, PAX5_dups = internal_del_dup(bam_file, "PAX5", targets)
    result["PAX5"]["deletions"] = []
    result["PAX5"]["duplications"] = []
    for d in PAX5_dels:
        sys.stderr.write(f"PAX5 deletion {d[0]} - {d[1]}: {PAX5_dels[d]} reads\n")
        if PAX5_dels[d] >= 3:
            result["PAX5"]["deletions"].append({"chrom": t2t.chrom_names[chrom], "start": d[0], "end": d[1], "reads": PAX5_dels[d]})
    for d in PAX5_dups:
        sys.stderr.write(f"PAX5 duplication {d[0]} - {d[1]}: {PAX5_dups[d]} reads\n")
        if PAX5_dups[d] >= 3:
            result["PAX5"]["duplications"].append({"chrom": t2t.chrom_names[chrom], "start": d[0], "end": d[1], "reads": PAX5_dups[d]})

    chrom, KMT2A_dels, KMT2A_dups = internal_del_dup(bam_file, "KMT2A", targets)
    result["KMT2A"]["duplications"] = []
    for d in KMT2A_dups:
        sys.stderr.write(f"KMT2A duplication {d[0]} - {d[1]}: {KMT2A_dups[d]} reads\n")
        if KMT2A_dups[d] >= 3:
            result["KMT2A"]["duplications"].append({"chrom": t2t.chrom_names[chrom], "start": d[0], "end": d[1], "reads": KMT2A_dups[d]})

    print(json.dumps(result))

# returns CN by (focal sequencing depth, local 1Mbp bins)
def depth(gene, targets, f, diploid_gene, covg, covg_2n, blast_ratio, exclude=["RUNX1", "CDKN2A", "CDKN2B", "IKZF1"]):
    focal_depth = f["qc"]["gene_coverage"][gene]
    # calculate focal depth from diploid genes only
    all_gene_depths = sorted([f["qc"]["gene_coverage"][g] for g in f["qc"]["gene_coverage"] if g in diploid_gene and g not in exclude and "DUX4" not in g])
    median_other = all_gene_depths[len(all_gene_depths)//2]
    if blast_ratio is not None:
        # caculate expected 3n cov
        focal_expected_3n = blast_ratio*3*median_other/2 + (1-blast_ratio)*median_other
        # calculate delta
        focal_delta = focal_expected_3n - median_other
        focal_cn = (focal_depth-median_other)/focal_delta +2 # should work for both high (3+) or low (1) CN
        sys.stderr.write(f"Est. {gene} focal CN (with estimated blast ratio {blast_ratio}): {focal_depth:.2f}/{median_other:.2f} ({focal_cn:.2f}x)\n")
        # write 0 if cn is negative
        if focal_cn < 0:
            focal_cn = 0
        #sys.stderr.write(f"For testing purpose, original focal CN (assuming broadly diploid): {focal_depth:.2f}/{median_other:.2f} ({focal_depth/median_other*2:.2f}x)\n")
    else:
        focal_cn = focal_depth/median_other*2
        sys.stderr.write(f"Est. {gene} focal CN: {focal_depth:.2f}/{median_other:.2f} ({focal_depth/median_other*2:.2f}x)\n")

    bins = [int(c[3]) for c in covg if t2t.chr_idx[c[0]] == targets[gene][0] and int(c[1]) < targets[gene][2] and int(c[2]) > targets[gene][1]] # should be at most 2 bins unless the gene is >>1Mbp
    bin_depth = sum(bins)/len(bins)
    # switch to using covg_2n from karyo
    if blast_ratio is not None:
        # caculate expected 3n cov
        expected_3n = blast_ratio*3*covg_2n/2 + (1-blast_ratio)*covg_2n
        # calculate delta
        delta = expected_3n - covg_2n
        cn = (bin_depth-covg_2n)/delta +2 # should work for both high (3+) or low (1) CN
        sys.stderr.write(f"Est. {gene} bin CN (with estimated blast ratio {blast_ratio}): {bin_depth:.2f}/{covg_2n:.2f} ({cn:.2f}x)\n")
        #sys.stderr.write(f"For testing purpose, original (assuming broadly diploid): {bin_depth}/{median_bin} ({bin_depth/median_bin*2:.2f}x)\n")
        #sys.stderr.write(f"Est. {gene} bin CN: {bin_depth}/{covg_2n} ({bin_depth/covg_2n*2:.2f}x)\n")
    else:
        cn = bin_depth/covg_2n*2
        sys.stderr.write(f"Est. {gene} bin CN: {bin_depth:.2f}/{covg_2n:.2f} ({cn:.2f}x)\n")
        #sys.stderr.write(f"For testing purpose, original (assuming broadly diploid): {bin_depth}/{median_bin} ({bin_depth/median_bin*2:.2f}x)\n")

    return (focal_cn, cn)


# calls intra-genic deletions with split reads (not very small dels), and *tandem* duplications
def internal_del_dup(bam_file, gene, targets):
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
    dups = {}
    for read in hits:
        for i in range(len(hits[read])):
            h0, (qs0, qe0, ts0, te0) = hits[read][i]
            for j in range(i+1, len(hits[read])):
                h1, (qs1, qe1, ts1, te1) = hits[read][j]

                # file must have been sorted, and alignments in order, so h0 must be the left and h1 the right side of a putative deletion
                # the read orientation is arbitrary however

                q_overlap = min(qe0, qe1) - max(qs0, qs1)
                if q_overlap > 100: # they significantly overlap in the read
                    continue

                if h0.is_reverse != h1.is_reverse: # implies an inversion of some kind, we don't care (at the moment)
                    continue

                if te0 - ts0 > min_anchor and te1 - ts1 > min_anchor:

                    # dups:
                    # 0 and 1 are always in order in the reference; if they overlap on the reference and are *out of order* in the read, they are a potential tandem dup
                    # TODO: what's a reasonable minimum size?? (how big would not be encoded as a simple intra-read insertion?)

                    # tandem duplications:
                    #         -------->
                    #        | aln 1  |
                    #  -------===--------
                    #    -------->
                    #   | aln 0  |
                    #  -------===--------
                    if (not h0.is_reverse and ((qe0 < qe1 and te0 - ts1 > 1000) or (h0.is_reverse and (qe1 < qe0 and te0 - ts1 > 1000)))):
                        for d in dups:
                            if abs(d[0] - ts1) < breakpoint_margin and abs(d[1] - te0) < breakpoint_margin: # merge to existing
                                dups[d] += 1
                                break
                        else:
                            dups[(ts1, te0)] = 1

                    if (h0.is_reverse and ((qe0 < qe1 and te1 - ts0 > 1000) or (not h0.is_reverse and (qe1 < qe0 and te1 - ts0 > 1000)))):
                        for d in dups:
                            if abs(d[0] - ts0) < breakpoint_margin and abs(d[1] - te1) < breakpoint_margin: # merge to existing
                                dups[d] += 1
                                break
                        else:
                            dups[(ts0, te1)] = 1

                    # dels:
                    # order in reference matches order in query, and there is a nonzero gap in the reference
                    elif ((not h0.is_reverse and qe0 < qe1) or (h0.is_reverse and qe0 > qe1)) and te0 - ts1 < 0:
                        # proper orientation for deletions or (non-dup/non-tandem) insertions:
                        #       -------->
                        #  aln / 0 /\ 1 \
                        #  --------__--------

                        for d in dels:
                            if abs(d[0] - te0) < breakpoint_margin and abs(d[1] - ts1) < breakpoint_margin: # merge to existing
                                dels[d] += 1
                                break
                        else:
                            dels[(te0, ts1)] = 1

                    # not a rational thing we care about at the moment
                    else:
                        pass
    return chrom, dels, dups

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

# use this to find segment of chrom & position
def get_segment(b):
    # these are T2T/CHM13v2 coordinates
    centromeres = {
      "chr1": (121405145, 144008994),
      "chr2": (87403406, 98345164),
      "chr3": (90339372, 96498140),
      "chr4": (49050159, 57435402),
      "chr5": (46161727, 51082874),
      "chr6": (57216843, 62665581),
      "chr7": (55900922, 66709007),
      "chr8": (43506737, 48466602),
      "chr9": (40119875, 81269165),
      "chr10": (38527422, 43443011),
      "chr11": (48819609, 55216869),
      "chr12": (30997656, 38541929),
      "chr16": (32625033, 52263389),
      "chr17": (20527834, 28139193),
      "chr18": (11350790, 21125454),
      "chr19": (19876180, 30396484),
      "chr20": (25835813, 37114981)
    }
    par1 = {
      "chrX": (0, 2394410),
      "chrY": (0, 2458320)
    }
    # PAR2 is only ~300Kbp, so we won't consider it separately, but it will be excluded from X/Y
    par2 = {
      "chrX": (153925834, 154259566),
      "chrY": (62122809, 62460029)
    }
    #b: (chromosome, bin_start, bin_end, ...)
    if b[0] in centromeres and b[2] < centromeres[b[0]][0]:
        name = b[0][3:]+"p"
    elif b[0] in centromeres and b[1] > centromeres[b[0]][1]:
        name = b[0][3:]+"q"
    elif b[0] in centromeres and b[1] < centromeres[b[0]][1] and b[2] > centromeres[b[0]][0]: # overlaps the centromere, skip it
        return None
    elif b[0] in par1 and b[1] < par1[b[0]][1] and b[2] > par1[b[0]][0]: # if any part of it overlaps PAR1
        if b[1] < par1[b[0]][0] or b[2] > par1[b[0]][1]: # if it does not completely overlap, ignore it
            return None
        else: # completely contained
            name = "PAR1"
    elif b[0] in par2 and b[2] > par2[b[0]][0] and b[1] < par2[b[0]][1]: # overlaps PAR2, skip it
        return None
    else:
        name = b[0][3:]
    return name



if __name__ == "__main__":
    parser = argparse.ArgumentParser("Call focal CNVs based on regional depth, specific deletions, duplications")
    parser.add_argument("bam", help="BAM file, aligned to T2T")
    parser.add_argument("bed", help="BED file of adaptive target genes")
    parser.add_argument("covg", help="Binned read depth file (from lrdk)")
    parser.add_argument("fusion", help="JSON output from fusion caller (includes target depths)")
    parser.add_argument("karyo", help="JSON output from lrdk")
    args = parser.parse_args()
    main(args.bam, args.bed, args.covg, args.fusion,args.karyo)
