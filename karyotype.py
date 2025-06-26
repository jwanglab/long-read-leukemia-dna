from collections import defaultdict
import sys
import argparse
import numpy
import numpy.ma as ma
import pysam
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from scipy import stats

class Reference:
    def __init__(self, chrom_list):
        self.chroms = chrom_list
        self.chrom_names = ["chr"+str(i+1) for i in range(22)] + ["chrX", "chrY"]
        self.chr_idx = {self.chroms[c]:c for c in range(len(self.chroms))}
        for c in range(len(self.chrom_names)):
            self.chr_idx[self.chrom_names[c]] = c
        self.chr_len = [None for c in self.chroms]
        self.chr_covg = [None for c in self.chroms]

T2T_accessions = ["NC_060925.1", "NC_060926.1", "NC_060927.1", "NC_060928.1", "NC_060929.1", "NC_060930.1", "NC_060931.1", "NC_060932.1", "NC_060933.1", "NC_060934.1", "NC_060935.1", "NC_060936.1", "NC_060937.1", "NC_060938.1", "NC_060939.1", "NC_060940.1", "NC_060941.1", "NC_060942.1", "NC_060943.1", "NC_060944.1", "NC_060945.1", "NC_060946.1", "NC_060947.1", "NC_060948.1"]

bin_size = 1000000


# covert bins from list of chrom,st,en,ct to segment(chr p/q arm, PAR, etc) with array of cts
def parse_tsv(tsv_file):
    #chromosome	bin_start	bin_end	n_reads
    data = [line.strip().split('\t') for line in open(tsv_file)]
    bins = [(d[0], int(d[1]), int(d[2]), int(d[3])) for d in data if d[0] != "chromosome"]
    new_bins = defaultdict(list)
    for b in bins:
        name = get_segment(b)
        if name is None:
            continue
        new_bins[name].append(int(b[3]))
    for chrom in new_bins:
        new_bins[chrom] = numpy.array(new_bins[chrom], dtype='u4')
    return new_bins


# returns the segment (p/q arm, chromosome, PAR) or None if it does not cleanly fit in one
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


# bins: {chrom: numpy array of u8, ...}
def find_levels(bins):
    a = numpy.concatenate([bins[c] for c in bins])
    binsize = int(numpy.percentile(a, 99.5) // 20)
    hist = numpy.histogram(a, bins=range(0, a.max(), binsize))

    # find peaks
    #peak1 = hist[1][numpy.argmax(hist[0])]
    maxima = []
    for i in range(1, hist[0].shape[0]-1):
        if hist[0][i] > hist[0][i-1] and hist[0][i] > hist[0][i+1]: # local maxima
            maxima.append((hist[1][i], hist[0][i]))
    maxima.sort(key=lambda a: -a[1])
    # fine-tune maxima peaks: median of pts under peak +- 1.5 bins, total area as num points within new peak +- 1
    for m in range(len(maxima)):
        new_peak = numpy.median(a[numpy.logical_and(a > maxima[m][0]-binsize*1.5, a < maxima[m][0]+binsize*1.5)])
        maxima[m] = (new_peak, numpy.count_nonzero(a[numpy.logical_and(a > new_peak-binsize, a < new_peak+binsize)]))
    return [m for m in maxima if m[1] > maxima[0][1]*0.01] # peaks have to have >1% maximum


def predict_karyo_v2(levels, medians, chrom_names, verbose):
    if verbose:
        sys.stderr.write(f"{medians} {levels} {len(medians)}\n")

    if len(levels) < 2: # close to 46XX, but assess each chromosome/arm independently assuming the peak is 2X
        cn2 = levels[0][0]
        cn1 = cn2 // 2
        cn3 = cn2 + cn1
    else:
        l1 = min(levels[0][0], levels[1][0])
        l2 = max(levels[0][0], levels[1][0])

        cn1 = None
        cn2 = None
        cn3 = None
        if l2 < l1*1.5:
            cn2 = l1
            cn3 = l2
        else:
            cn1 = l1
            cn2 = l2

    # estimate missing CNs
    if cn3 is None:
        cn3 = cn2 + (cn2 - cn1)
        if verbose:
            sys.stderr.write(f"CN1 is {cn1}, CN2 is {cn2}, est. CN3 at {cn3}\n")
    if cn1 is None:
        cn1 = cn2 - (cn3 - cn2)
        # this estimates the blast 1X, if the blast count isn't 100%, this will be different from the germline 1X that might be present in ex. X and Y
        if verbose:
            sys.stderr.write(f"CN2 is {cn2}, CN3 is {cn3}, est. CN1 at {cn1}\n")
    delta = cn3 - cn2

    cns = []
    cn_depths = defaultdict(list)
    for i in range(len(medians)):
        if medians[i] > cn2 - delta/2 and medians[i] < cn2 + delta/2:
            if verbose:
                sys.stderr.write(f"{chrom_names[i]} {medians[i]} x2\n")
            cns.append(2)
        elif medians[i] > cn1 - delta/2 and medians[i] < cn1 + delta/2:
            if verbose:
                sys.stderr.write(f"{chrom_names[i]} {medians[i]} x1\n")
            cns.append(1)
        elif medians[i] > cn3 - delta/2 and medians[i] < cn3 + delta/2:
            if verbose:
                sys.stderr.write(f"{chrom_names[i]} {medians[i]} x3\n")
            cns.append(3)
        elif medians[i] > cn3 + delta/2:
            est_cn = (medians[i] - cn2) / delta + 2
            if verbose:
                sys.stderr.write(f"{chrom_names[i]} {medians[i]} est CN {est_cn:.2f}x\n")
            cns.append(est_cn)
        elif medians[i] < delta / 2:
            if verbose:
                sys.stderr.write(f"{chrom_names[i]} {medians[i]} x0\n")
            cns.append(0)
        else: # median must be somewhere in between one of these bins
            if verbose:
                sys.stderr.write(f"{chrom_names[i]} {medians[i]} ???\n")
            cns.append(-1)
        cn_depths[cns[-1]].append(medians[i])

    # fix a specific kind of lower-blast hypodiploid where the single-copy loss is called diploid, the original diploid is called triploid, and the germline haploid (X,Y probably) are called haploid
    # this is only possible when it's XY, otherwise nothing comes up germline haploid
    if cns[-1] == 1 and cns[:-2].count(1) == 0:
        cn_1 = sum(cn_depths[1]) / len(cn_depths[1])
        cn_2 = sum(cn_depths[2]) / len(cn_depths[2])
        cn_3 = (sum(cn_depths[3]) / len(cn_depths[3])) if len(cn_depths[3]) > 0 else 0
        # if it actually looks like triploid is actually ~2x haploid,
        if cn_3 / 2 > cn_1 * 0.9 and cn_3 / 2 < cn_1 * 1.1:
            sys.stderr.write("---------- Adjusting to low-pentrance hypodiploid! ----------\n")
            # make "diploid" haploid, and "triploid" diploid
            cns = [({2:1, 3:2}[a] if a in [2,3] else a) for a in cns]

    # make some heuristic adjustments if there are way more 3x than expected and no 1x
    # ...suggesting a lower-blast aneuploidy, so levels should be shifted down
    if cns.count(1) == 0:
        if cns.count(3) > 25:
            cns = [c-1 if c != 0 else 0 for c in cns]

    return {chrom_names[i]:int(round(cns[i])) for i in range(len(chrom_names))}


def digital_karyotype(bins, verbose):
    chrom_names = [c for c in bins] # keys include chromosome arms, whole chromosomes for those that are acrocentric, and PAR1
    medians = [None for i in range(len(chrom_names))]
    edges = [] # contains tuple pairs representing significantly different coverage

    levels = find_levels(bins)
    for i in range(len(chrom_names)):
        a = bins[chrom_names[i]]
        medians[i] = numpy.median(a)
    cns = predict_karyo_v2(levels, medians, chrom_names, verbose)

    return cns


def from_bam(bam_file, repeatmasker_file, ref, read_lim=None, outdir=None, verbose=False):

    # ------------- Open BAM and detect reference ----------------
    sys.stderr.write("Opening BAM file...\n")
    # !!!! BAM files do not follow our normal assumptions from PAF files that all alignments of a reads are together - this is routinely violated in a BAM generated by dorado!
    # ignore apparently truncated BAM files so that it can be run on in-progress runs
    try:
        af = pysam.AlignmentFile(bam_file, mode='rb', ignore_truncation=True)
    except FileNotFoundError as e:
        sys.stderr.write(f"Error: File not found: {bam_file}\n")
        sys.exit(1)

    # check BAM file for the appropriate reference to use
    if af.references[0] not in ref.chroms and not af.references[0].startswith("chr"):
        sys.stderr.write("ERROR: BAM references do not match T2T accessions or 'chr' names, check that reads are aligned to T2T/CHM13v2\n")
        sys.exit(1)

    # ------------ Build repeat/adaptive sampling mask -------------
    mask = [[] for i in range(len(ref.chroms))]

    sys.stderr.write("Reading repeats...\n")
    # collect repeats from repeatmasker file and Ns
    with open(repeatmasker_file) as rf:
        for line in rf:
            parts = line.strip().split('\t')
            chrom, st, en = parts[0], int(parts[1]), int(parts[2])
            if chrom not in ref.chrom_names:
                continue
            chrom = ref.chrom_names.index(chrom)
            mask[chrom].append((st,en))

    # sort masked features by start position
    for i in range(len(mask)):
        mask[i].sort(key = lambda a: a[0])

    # ------------ Read BAM and count reads/coverage -------------
    last = None
    skipped_refs = set()

    bitmask = [numpy.zeros(af.lengths[af.references.index(ref.chroms[c] if ref.chroms[c] in af.references else ref.chrom_names[c])], dtype='u1') for c in range(len(mask))]
    for c in range(len(mask)):
        for a in mask[c]:
            bitmask[c][a[0]:a[1]] = 1

    sys.stderr.write("Processing reads...\n")
    n_reads = 0
    for a in af.fetch(until_eof=True):
        if read_lim is not None and n_reads > read_lim:
            break
        if a.query_name != last:
            last = a.query_name

        if a.reference_name not in ref.chr_idx:
            skipped_refs.add(a.reference_name)
            continue
        c = ref.chr_idx[a.reference_name]
        if ref.chr_covg[c] is None:
            l = af.lengths[af.references.index(a.reference_name)]
            ref.chr_len[c] = l
            ref.chr_covg[c] = numpy.zeros(l+1, dtype='u4')
            if verbose:
                sys.stderr.write(f"{c}\t{a.reference_name}\t{l}\n")
        ref.chr_covg[c][int(a.reference_start)] += 1
        n_reads += 1

    sys.stderr.write(f"WARNING: {len(list(skipped_refs))} references with alignments not found and skipped: [{', '.join(str(r) for r in skipped_refs)}]\n")

    # ------------ binning and masking, writing coverage files -------------
    if outdir is not None:
        sys.stderr.write(f"Saving chromosome coverages to {outdir}/\n")
    bins = [None for i in range(len(ref.chroms))]
    for c in range(len(ref.chroms)):
        if ref.chr_covg[c] is None:
            sys.stderr.write(f"WARNING: chrom '{ref.chr_covg[c]}' coverage is missing - maybe there were no assigned reads\n")
            continue
        # mask repeats, Ns, enriched genes
        ref.chr_covg[c] = ma.masked_array(ref.chr_covg[c])
        for region in mask[c]:
            ref.chr_covg[c][region[0]:region[1]] = ma.masked
        bins[c] = numpy.zeros((ref.chr_covg[c].shape[0]//bin_size + 1, 4), dtype='u4') # total_reads, start, end, [unused]
        bin_mask = numpy.zeros((ref.chr_covg[c].shape[0]//bin_size + 1), dtype='bool')
        for i in range(bins[c].shape[0]):
            st = i * bin_size
            en = st + bin_size
            bin_tot = ref.chr_covg[c][st:en].sum()
            bin_ct = ref.chr_covg[c][st:en].count() # count unmasked values
            # more than 25% of sites in this bin must be unmasked to include the bin at all
            if bin_ct > bin_size * 0.25:
                bin_mask[i] = True
                bin_tot *= bin_size/bin_ct # normalize the total to adjust for masked sites
                bins[c][i,0] = bin_tot
                bins[c][i,1] = st
                bins[c][i,2] = en
        bins[c] = bins[c][bin_mask,:]
        if verbose:
            sys.stderr.write(f"{ref.chroms[c]}: {ref.chr_len[c]} nt, masked reads / {bin_size}: {numpy.median(bins[c][:,0]):.2f}\n")

    if outdir is not None:
        fout = open(f"{outdir}/coverage_1Mbp_bins.tsv", 'w')
        fout.write(f"chromosome\tbin_start\tbin_end\tn_reads\n")
        for c in range(len(ref.chr_covg)):
            if ref.chr_covg[c] is not None:
                #ref.chr_covg[c].tofile(f"{outdir}/{ref.chrom_names[c]}_covg.npy")
                for i in range(bins[c].shape[0]):
                    fout.write(f"{ref.chrom_names[c]}\t{bins[c][i,1]}\t{bins[c][i,2]}\t{bins[c][i,0]}\n")
        fout.close()

    new_bins = defaultdict(list)
    for c in range(len(bins)):
        for i in range(bins[c].shape[0]):
            name = get_segment((ref.chrom_names[c], bins[c][i,1], bins[c][i,2]))
            if name is None:
                continue
            new_bins[name].append(bins[c][i,0])
    for chrom in new_bins:
        new_bins[chrom] = numpy.array(new_bins[chrom], dtype='u4')
    return new_bins


def ISCN_string(cns):
    s = "seq"
    cts = defaultdict(list)
    for c in cns:
        cts[cns[c]].append(c)
    # merge p and q if they have the same copy #
    for n in cts:
        for aut in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 16, 17, 18, 19, 20]: # ~metacentric, we count arms separately
            if f"{aut}p" in cts[n] and f"{aut}q" in cts[n]:
                cts[n].insert(cts[n].index(f"{aut}p"), f"{aut}") # this preserves a logical order
                cts[n].remove(f"{aut}p")
                cts[n].remove(f"{aut}q")
    if len(cts[1]) > 0:
        s += f"({', '.join(c for c in cts[1])})x1 "
    # check if autosomes are ALL diploid
    aut_diploid = len(set(cts[2]) & set([f"{c}" for c in range(1,23)])) == 22
    if aut_diploid or 'X' in cts[2] or 'Y' in cts[2]:
        dips = []
        if aut_diploid:
            dips.append('1-22')
        if 'X' in cts[2]:
            dips.append('X')
        if 'Y' in cts[2]:
            dips.append('Y')
        s += f"({', '.join(dips)})x2 "
    for n in range(3, max(n for n in cts)+1):
        if len(cts[n]) > 0:
            s += f"({', '.join(c for c in cts[n])})x{n} "
    return s


# bins should be in the format {'1p': [numpy array], '1q':..., 'PAR1': ...}
def plot_karyo(bins, outfile):
    plt.figure(figsize=(20,4), dpi=100)
    chr_locs = []
    x = 0
    for c in bins:
        med = numpy.median(bins[c])
        plt.scatter([x+i for i in range(bins[c].shape[0]) if bins[c][i] != 2**16-1], [b for b in bins[c] if b != 2**16-1], color="blue", marker='.')
        plt.plot([x, x+bins[c].shape[0]-1],[med,med], color="red", linestyle="solid", linewidth=3)
        chr_locs.append(x+bins[c].shape[0]/2)
        x += bins[c].shape[0]
        plt.plot([x+5, x+5], [0, 2**16-1], color="black")
        x += 10
    ymax = numpy.percentile(numpy.concatenate([bins[c] for c in bins]), 99.9)

    # side histogram
    a = numpy.concatenate([bins[c] for c in bins])
    binsize = int(numpy.percentile(a, 99.5) // 20)
    hist = numpy.histogram(a, bins=range(0, a.max(), binsize))
    plt.plot(hist[0], [hist[1][i]+(hist[1][i+1]-hist[1][i])//2 for i in range(hist[0].shape[0])], color="orange")

    plt.xticks(chr_locs, [c for c in bins], rotation='vertical')
    plt.ylim((0,ymax))
    plt.xlim((0,x-10))
    plt.xlabel("Chromosome")
    plt.ylabel("Reads per Mbp")
    plt.title("Karyotype")
    plt.tight_layout()
    plt.savefig(outfile)


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Digital karyotyping from BAM or binned depth file")
    parser.add_argument("input", help="BAM/SAM of reads aligned to T2T or TSV of binned sequencing depth", nargs='+')
    parser.add_argument("--repeats", help="RepeatMasker BED file")
    parser.add_argument("--outdir", help="Output directory/prefix")
    parser.add_argument("--read_lim", help="Maximum reads to process (default is all)", type=int)
    parser.add_argument("-v", help="Verbose", action="store_true", default=False)
    args = parser.parse_args()

    t2t = Reference(T2T_accessions) # no covg data herein, just for reference
    for a in args.input:
        if a.endswith(".tsv"):
            bins = parse_tsv(a)
        elif a.endswith(".bam"):
            if args.repeats is None:
                sys.stderr.write("ERROR: --repeats required if using a BAM file\n")
                sys.exit(1)
            # each get its own copy of the reference with empty covg lists
            bins = from_bam(a, args.repeats, Reference(T2T_accessions), args.read_lim, args.outdir, args.v)
        else:
            sys.stderr.write("ERROR: unrecognized extension '{a[a.rindex('.'):]}'\n")
            sys.exit(1)
        if args.outdir is not None:
            plot_karyo(bins, f"{args.outdir}/karyotype.png")

        cns = digital_karyotype(bins, args.v)
        for c in cns:
            sys.stderr.write(f"{c}: {cns[c]}\n")
        print("sample\t" + "\t".join(c for c in cns) + "\tISCN")
        print(f"{a}\t" + "\t".join(str(cns[c]) for c in cns) + f"\t{ISCN_string(cns)}")
        sys.stderr.write(f"{ISCN_string(cns)}\n")

