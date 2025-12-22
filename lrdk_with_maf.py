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
import math
from scipy.signal import find_peaks

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

def find_levels(bins, trim_percentile, binsize):
    a = numpy.concatenate([bins[c] for c in bins])
    # use mid point instead of left point
    counts, bin_edges=numpy.histogram(a, bins=range(0, int(max(numpy.percentile(bins[c], trim_percentile) for c in bins))+2*binsize, int(binsize)))
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    # find peaks; adjust for peaks that are very close together (separated by >2 binwidth)
    local_max = find_peaks(counts,distance=3)

    maxima = []
    # go through the indices and append bin_center and count
    for i in local_max[0]:
        maxima.append((bin_centers[i], counts[i]))
    # sort by counts (find the largest peaks)
    maxima.sort(key=lambda a: -a[1])
    # fine-tune maxima peaks: median of pts under peak +- 1.5 bins, total area as num points within new peak +- 1
    for m in range(len(maxima)):
        new_peak = numpy.median(a[numpy.logical_and(a > maxima[m][0]-binsize*1.5, a < maxima[m][0]+binsize*1.5)])
        maxima[m] = (new_peak, numpy.count_nonzero(a[numpy.logical_and(a > new_peak-binsize, a < new_peak+binsize)]))
    levels = [m for m in maxima if m[1] > maxima[0][1]*0.05] # peaks have to have >5% maximum

    # filter peaks
    medians = {c: numpy.median(values) for c, values in bins.items()}
    # find segments for each level
    lvls_seg = defaultdict(list)
    top_levels = numpy.array([i[0] for i in levels])
    # assigning a level to each segment
    for m in medians:
        arr = top_levels/medians[m] # reverse so denominator is the same
        idx = numpy.argmin(numpy.abs(arr - 1))
        if numpy.abs(arr[idx]-1)<0.3: # sanity check, cannot be too far away
            lvls_seg[top_levels[idx]].append(m)
    # remove level if it doesnt have any segment
    new_levels = []
    for l in levels:
        if len(lvls_seg[l[0]])>0:
            new_levels.append(l)

    # filter out peaks that are too close (by 9.5%)
    just_levels = numpy.array([l[0] for l in new_levels])
    tolerance = 0.095 # consider "small fraction gain", value that works for test samples
    keep_ind = set()
    discard_ind = set()
    for i in range(len(just_levels)):
        lower = just_levels[i] * (1 - tolerance)
        upper = just_levels[i] * (1 + tolerance)
        indices = numpy.where((just_levels >= lower) & (just_levels <= upper))[0]
        keep_ind.add(indices[0]) # only keep the first one
        if len(indices)>1:
            for i in indices[1:]:
                discard_ind.add(i) # make sure not to include these
    res = []
    keep_ind = keep_ind-discard_ind
    for i in sorted(keep_ind):
        res.append(new_levels[i])
    return res

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
        if l2 < l1*1.51: # one real sample was 1.5003, so we need a little bit of wiggle room, but they are usually <1.5x
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

    # Estimate blast ratio. Use copy numbers after adjustment.
    cn_depths_auto = defaultdict(list)
    for i in range(len(medians)-2): # exclude X and Y
        cn_depths_auto[cns[i]].append(medians[i])

    blast_ratio = estimate_blast_ratio(cn_depths_auto)
    sys.stderr.write(f"blast_ratio in predict_karyo_v2: {blast_ratio}\n")

    return {chrom_names[i]:int(round(cns[i])) for i in range(len(chrom_names))},blast_ratio


def digital_karyotype(bins, trim_percentile, binsize, verbose):
    chrom_names = [c for c in bins] # keys include chromosome arms, whole chromosomes for those that are acrocentric, and PAR1
    medians = [None for i in range(len(chrom_names))]
    edges = [] # contains tuple pairs representing significantly different coverage

    levels = find_levels(bins, trim_percentile, binsize)
    for i in range(len(chrom_names)):
        a = bins[chrom_names[i]]
        medians[i] = numpy.median(a)
    cns, br = predict_karyo_v2(levels, medians, chrom_names, verbose)

    return cns, br


def from_bam(bam_file, repeatmasker_file, ref, read_lim=None, outdir=None, verbose=False):

    # ------------- Open BAM and detect reference ----------------
    sys.stderr.write("Opening BAM file...\n")
    # !!!! BAM files do not follow our normal assumptions from PAF files that all alignments of a reads are together - this is routinely violated in a BAM generated by dorado!
    # ignore apparently truncated BAM files so that it can be run on in-progress runs
    try:
        if bam_file.endswith(".bam"):
            af = pysam.AlignmentFile(bam_file, mode='rb', ignore_truncation=True)
        elif bam_file.endswith(".sam"):
            af = pysam.AlignmentFile(bam_file, mode='r', ignore_truncation=True)
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


def karyotype_string(cns):
    cts = defaultdict(list)
    for c in cns:
        cts[cns[c]].append(c)
    mods = []
    chrom_ct = 0
    for aut in range(1,23):
        if f'{aut}' in cts[2] or f'{aut}p' in cts[2] or f'{aut}q' in cts[2]:
            ct = 2
        elif f'{aut}' in cts[1] or f'{aut}p' in cts[1] or f'{aut}q' in cts[1]:
            ct = 1
        else:
            ct = cns[f'{aut}'] if f'{aut}' in cns else min((cns[f'{aut}p'], cns[f'{aut}q']))
        chrom_ct += ct
    baseline = 1 if chrom_ct <= 30 else 2
    for aut in range(1,23):
        if f'{aut}' in cns or cns[f'{aut}p'] == cns[f'{aut}q']:
            for i in range(baseline - cns[f'{aut}' if f'{aut}' in cns else f'{aut}p']):
                mods.append(f"-{aut}")
            for i in range(cns[f'{aut}' if f'{aut}' in cns else f'{aut}p'] - baseline):
                mods.append(f"+{aut}")
        else:
            for i in range(baseline - cns[f'{aut}p']):
                mods.append(f"-{aut}p")
            for i in range(cns[f'{aut}p'] - baseline):
                mods.append(f"+{aut}p")
            for i in range(baseline - cns[f'{aut}q']):
                mods.append(f"-{aut}q")
            for i in range(cns[f'{aut}q'] - baseline):
                mods.append(f"+{aut}q")
    chrom_ct += cns['X'] + cns['Y']
    xy = "XY" if cns['Y'] >= 1 else ("XX" if cns['X'] >= 2 else ("X" if cns['X'] == 1 else "??"))
    for i in range(xy.count('X') - cns['X']):
        mods.append("-X")
    for i in range(cns['X'] - xy.count('X')):
        mods.append("+X")
    for i in range(xy.count('Y') - cns['Y']):
        mods.append("-Y")
    for i in range(cns['Y'] - xy.count('Y')):
        mods.append("+Y")
    return f"{chrom_ct}{xy}; " +("(1n) " if baseline == 1 else "") + ', '.join(mods)


def ISCN_string(cns):
    s = "seq "
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
    groups = []
    if len(cts[1]) > 0:
        groups.append(f"({', '.join(c for c in cts[1])})x1 ")
    # check if autosomes are ALL diploid
    aut_diploid = len(set(cts[2]) & set([f"{c}" for c in range(1,23)])) == 22
    if aut_diploid or 'X' in cts[2] or 'Y' in cts[2]:
        dips = []
        if 'X' in cts[2]:
            dips.append('X')
        if 'Y' in cts[2]:
            dips.append('Y')
        if aut_diploid:
            dips.append('1-22')
        groups.append(f"({', '.join(dips)})x2")
    for n in range(3, max(n for n in cts)+1):
        if len(cts[n]) > 0:
            groups.append(f"({', '.join(c for c in cts[n])})x{n}")
    return s + ", ".join(groups)

# 20251118 MAF code
def find_bin_width2(bins, alpha=0.12, trim_percentile=95, deviation=35):
    """find bin width function that takes range into consideration"""
    a = numpy.concatenate([bins[c] for c in bins])
    iqr = numpy.percentile(a, 50+deviation) - numpy.percentile(a, 50-deviation)
    n = len(a)
    maxi = numpy.percentile(a, trim_percentile)
    binsize = ((1-alpha)*iqr + alpha*maxi/2)*2*(n**(-1/3))
    return binsize

def within_segment_spread(bins):
    """Calculate spread within segment coverage. If spread is too high, raise a warning flag."""
    seg_iqr = []
    for c in bins:
        seg_iqr.append(numpy.percentile(bins[c], 75)-numpy.percentile(bins[c], 25))
    # remove outlier and y; get mean
    spread = numpy.mean((numpy.percentile(seg_iqr[:-1], 75)-numpy.percentile(seg_iqr[:-1], 25)))
    # get medians
    med = numpy.median(numpy.concatenate([bins[c] for c in bins]))
    return spread/med

# hard code the plot width if we don't plan to change that
seg_plot_width = {'1p': 115, '1q': 99, '2p': 86, '2q': 139, '3p': 85, '3q': 99, '4p': 44, '4q': 132, '5p': 41, '5q': 125,
                  '6p': 53, '6q': 106, '7p': 55, '7q': 88, '8p': 37, '8q': 94, '9p': 39, '9q': 68, '10p': 35, '10q': 86,
                  '11p': 46, '11q': 75, '12p': 30, '12q': 91, '13': 92, '14': 82, '15': 75, '16p': 24, '16q': 39, '17p': 18,
                  '17q': 50, '18p': 11, '18q': 55, '19p': 19, '19q': 26, '20p': 21, '20q': 28, '21': 29, '22': 29, 'X': 121, 'Y': 6}

def get_segment_pos(chrom, pos):
    """Get segment position for MAF use"""
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
    new_pos = ''
    if chrom in centromeres and pos < centromeres[chrom][0]:
        name = chrom[3:]+"p"
        new_pos = pos
    elif chrom in centromeres and pos > centromeres[chrom][1]:
        name = chrom[3:]+"q"
        new_pos = pos-centromeres[chrom][1]
    # within centromere regions
    elif chrom in centromeres and pos < centromeres[chrom][1] and pos > centromeres[chrom][0]:
        name = None
    # don't worry about PAR; XY stay XY
    elif chrom in par1:
        name = chrom[3:]
        new_pos = pos
    else:
        name = chrom[3:]
        new_pos = pos

    if name:
        return name,new_pos
    else:
        return name


def find_bin_width_maf(array, deviation=25):
    """find bin width with MAF data"""
    # default deviation 25 means 75 - 25 percentile
    iqr = numpy.percentile(array, 50+deviation) - numpy.percentile(array, 50-deviation)
    n = len(array)
    binsize = 2*iqr*(n**(-1/3))
    return binsize


def predict_karyo_maf(levels, lvls_maf_peak, medians, chrom_names, verbose):
    if verbose:
        sys.stderr.write(f"{medians} {levels} {len(medians)}\n")
    # if no maf info, return None
    if lvls_maf_peak is None:
        return None

    blast_ratio = None

    if len(levels) < 2: # close to 46XX, but assess each chromosome/arm independently assuming the peak is 2X
        cn2 = levels[0][0]
        cn1 = cn2 // 2
        cn3 = cn2 + cn1
        if lvls_maf_peak[cn2]>0.4:
            if verbose:
                sys.stderr.write(f"only 1 coverage peak (2n), MAF checks out: {lvls_maf_peak[cn2]}\n")
        else:
            if verbose:
                sys.stderr.write(f"only 1 coverage peak (2n), MAF is too low: {lvls_maf_peak[cn2]}\n")
    else:
        cn1 = None
        cn2 = None
        cn3 = None
        # mostly care about the top 2 peaks; consider where the peaks are too
        # this is the lower level
        l1 = min(levels[0][0], levels[1][0])
        # this is the higher level
        l2 = max(levels[0][0], levels[1][0])

        # deal with different cases
        # if MAF is available for both levels
        if l1 in lvls_maf_peak and l2 in lvls_maf_peak:
            # l1 is 2n, l2 is 3n
            if lvls_maf_peak[l1]>0.4: #set lower level to 2n
                if l2/l1 >=1.6: # coverage ratio sanity check
                    cn1 = l1
                    cn2 = l2
                    if verbose:
                        sys.stderr.write(f"1n maf high, but ratio>=1.6. cn2/cn1 = {cn2/cn1}\n")
                else:
                    cn2 = l1
                    cn3 = l2
                    if verbose:
                        sys.stderr.write(f"cn3/cn2 = {cn3/cn2}\n")
            # l1 is 1n, l2 is 2n
            elif lvls_maf_peak[l2]>0.4 and lvls_maf_peak[l2]>lvls_maf_peak[l1]:
                cn1 = l1
                cn2 = l2
                if verbose:
                    sys.stderr.write(f'cn2/cn1 = {cn2/cn1}\n')
            else:
                if verbose:
                    sys.stderr.write(f"I dont know whats going on. lvls_maf_peak: {lvls_maf_peak}\n")
                return None

        # after more stringent maf quantity filtering, we should not run into the following cases
        elif l1 not in lvls_maf_peak and l2 in lvls_maf_peak:
            # 2n with 1n; 1n no MAF data
            if lvls_maf_peak[l2]>0.4:
                cn1 = l1
                cn2 = l2
                if verbose:
                    sys.stderr.write(f"only 2n MAF peak available; 1n no MAF. cn2/cn1 = {cn2/cn1}\n")
            else:
                cn2 = l1
                cn3 = l2
                if verbose:
                    sys.stderr.write(f"only 3n MAF peak available, 2n no MAF. cn3/cn2 = {cn3/cn2}\n")
        # 2n with 3n; 3n no MAF data
        elif l2 not in lvls_maf_peak and l1 in lvls_maf_peak:
            if lvls_maf_peak[l1]>0.4:
                cn2 = l1
                cn3 = l2
                if verbose:
                    sys.stderr.write(f"only 2n MAF peak available, 3n no MAF. cn3/cn2 = {cn3/cn2}\n")
            else:
                if verbose:
                    sys.stderr.write(f"I dont know whats going on. lvls_maf_peak: {lvls_maf_peak}\n")
                return None
        else:
            if verbose:
                sys.stderr.write('no MAF info for first two levels\n')
            return None


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

    if verbose:
        sys.stderr.write(f"cn1, cn2, cn3: {cn1} {cn2} {cn3}\n")

    cns = []
    cn_depths = defaultdict(list)
    cn_depths_auto = defaultdict(list)
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
        if i < len(medians)-2: # autosomes only
            cn_depths_auto[cns[-1]].append(medians[i])

    # fix autosomes first
    # parse cns, change remaining -1 to 1 (autosomes)
    for i in range(len(cns[:-2])):
        if cns[i] == -1:
            cns[i] = 1
            cn_depths[1].append(medians[i])
            cn_depths_auto[1].append(medians[i])
            sys.stderr.write(f"changing {chrom_names[i]} from -1 to 1\n")

    # calculate blast ratio
    blast_ratio = estimate_blast_ratio(cn_depths_auto)

    # adjust for low blast XY?
    # we assume they cannot have zero copy of chrom except for Y; make them 1 copy
    if cns[-1]==-1: # y is -1
        cn_2 = sum(cn_depths[2]) / len(cn_depths[2])
        cn_1 = cn_2/2
        if  medians[-1] > cn_1 * 0.9 and  medians[-1] < cn_1 * 1.1: # is Y actually around 1n?
            y_cn1 = medians[-1]
            cns[-1] = 1
            sys.stderr.write(f"y ~ cn2/2;set y to 1. ")
            if cns[-2]==-1 or cns[-2]==1: # only try to adjust X if X equals -1, 1
                if blast_ratio:
                    estimated_x_cn2 = 2*y_cn1*blast_ratio + (1-blast_ratio)*y_cn1
                    if medians[-2] > estimated_x_cn2 * 0.9 and  medians[-2] < estimated_x_cn2 * 1.1:
                        cns[-2] = 2 # adjusting X to 2n
                        sys.stderr.write(f"set X to 2n, with blast percentage\n")
                    else:
                        cns[-2] = 1
                        sys.stderr.write(f"set X to 1n, with blast percentage\n")
                else:
                    # there was only 2n available
                    if (cn_2-medians[-2]) > (medians[-2]-cn_1): # if X cov closer to 1n, then 1n
                        cns[-2] = 1
                        sys.stderr.write(f"set X to 1n, without blast percentage\n")
                    else:
                        cns[-2] = 2
                        sys.stderr.write(f"set X to 2n, without blast percentage\n")
        else:
            # we don't have a Y 1n reference, give up
            sys.stderr.write(f"no y 1n reference, set XY to 1.\n")
            cns[-1] = 1
            cns[-2] = 1


    return [{chrom_names[i]:int(round(cns[i])) for i in range(len(chrom_names))}, blast_ratio]


# package blast ratio code in a function
def estimate_blast_ratio(cn_medians_auto):
    has2 = len(cn_medians_auto[2])>0
    has1 = len(cn_medians_auto[1])>0
    has3 = len(cn_medians_auto[3])>0
    sys.stderr.write(f"{len(cn_medians_auto[2])} segments are 2n\n")
    sys.stderr.write(f"{len(cn_medians_auto[1])} segments are 1n\n")
    sys.stderr.write(f"{len(cn_medians_auto[3])} segments are 3n\n")
    if has2:
        CN2 = sum(cn_medians_auto[2]) / len(cn_medians_auto[2])
    if has3>0:
        CN3 = sum(cn_medians_auto[3]) / len(cn_medians_auto[3])
    if has1:
        CN1 = sum(cn_medians_auto[1]) / len(cn_medians_auto[1])
    # if they are all available
    if has1 and has2 and has3:
        if len(cn_medians_auto[1]) >= len(cn_medians_auto[3]):
            # calculate blast ratio from 1n and 2n
            blast_ratio = (1-(CN1/CN2))*2
            sys.stderr.write(f"blast_ratio from 1n 2n: {blast_ratio}\n")
        else:
            # calculate blast ratio from 2n and 3n
            blast_ratio = (CN3/CN2*2)-2
            sys.stderr.write(f"blast_ratio from 2n 3n: {blast_ratio}\n")
    elif has1 and has2:
        # calculate blast ratio from 1n and 2n
        blast_ratio = (1-(CN1/CN2))*2
        sys.stderr.write(f"blast_ratio from 1n 2n: {blast_ratio}\n")
    elif has2 and has3:
        # calculate blast ratio from 2n and 3n
        blast_ratio = (CN3/CN2*2)-2
        sys.stderr.write(f"blast_ratio from 2n 3n: {blast_ratio}\n")
    else:
        blast_ratio = None
        sys.stderr.write("There is only autosomal diploid, we cannot estimate blast ratio\n")

    return blast_ratio


# read enriched region information for plotting MAF data
seg_bases = {'1p': 20725, '1q': 9866, '2p': 6068, '2q': 15246, '3p': 6912, '3q': 14854, '4p': 981, '4q': 14323, '5p': 2773,
             '5q': 17629, '6p': 3578, '6q': 20454, '7p': 11499, '7q': 24100, '8p': 2515, '8q': 12487, '9p': 22318, '9q': 10544,
             '10p': 2854, '10q': 8013, '11p': 7732, '11q': 7528, '12p': 9862, '12q': 5815, '13': 7174, '14': 25099, '15': 9105,
             '16p': 7366, '16q': 5599, '17p': 3163, '17q': 13396, '18q': 9210, '19p': 12326, '19q': 5317, '20p': 4169, '20q': 3038,
             '21': 7285, '22': 15080, 'X': 10688}

# for chrom conversion
t2t_chrom = {}
chr_names = ["chr"+str(i+1) for i in range(22)] + ["chrX", "chrY"]
for i in range(24):
    t2t_chrom[T2T_accessions[i]] = chr_names[i]

# read AS file to get segment information
def read_AS(file):
    seg_region_start = defaultdict(list)
    seg_region_len = defaultdict(list)
    with open(file,'r') as pf:
        for line in pf:
            fs = line.split('\t')
            start = int(fs[1])
            chr = t2t_chrom[fs[0]]
            region_len = int(fs[2])-int(fs[1])
            # check if segment in centromeres
            if get_segment_pos(chr, start):
                seg, pos = get_segment_pos(chr, start)
                seg_region_start[seg].append(pos)
                seg_region_len[seg].append(region_len)
    # each region cumsum; start with zero
    seg_region_cumsum = {}
    for c in seg_region_len:
        seg_region_cumsum[c] = numpy.concatenate(([0],numpy.cumsum(seg_region_len[c])))
    # seg total region length
    seg_total_len = defaultdict(int)
    for c in seg_region_len:
        seg_total_len[c] = numpy.sum(seg_region_len[c])
    seg_conversion_factor = {}
    for s in seg_plot_width:
        if s!='18p': # 18p has no enriched region
            seg_conversion_factor[s] = math.ceil(seg_total_len[s]/seg_plot_width[s])
    return seg_region_start,seg_region_len,seg_region_cumsum,seg_total_len,seg_conversion_factor


def find_levels_maf(maf_data, levels, bins):
    """Given coverage levels, find MAF data associated with each level."""
    seg_maf = defaultdict(list)
    for s in maf_data:
        for m in maf_data[s]:
            seg_maf[s].append(m)
    # find medians
    medians = {}
    for c in bins:
        medians[c]=numpy.median(bins[c])
    # find segments for each level; only look at the first three levels
    if len(levels)>3:
        levels = levels[:3]
    lvls_seg = defaultdict(list)
    top_levels = numpy.array([i[0] for i in levels])
    # assigning a level to each segment
    for m in medians:
        arr = medians[m]/top_levels
        #print(m, arr)
        idx = numpy.argmin(numpy.abs(arr - 1))
        if numpy.abs(arr[idx]-1)<0.3: # sanity check, cannot be too far away
            lvls_seg[top_levels[idx]].append(m)
    # make sure each level has some seg
    levels_seg = {}
    for l in top_levels:
        if len(lvls_seg[l])>1:
            levels_seg[l] = lvls_seg[l]
        else:
            levels_seg[l] = [] # set as empty list for now
    # combine maf for each level
    levels_maf = defaultdict(list)
    for l in levels_seg:
        for s in levels_seg[l]:
            levels_maf[l] = levels_maf[l]+seg_maf[s]
    # check if there are enough maf data points. If too few data points, don't use this level
    for l in levels_maf:
        if len(levels_maf[l])<50:
            levels_maf[l] = []
    return levels_maf


# for each group of MAF data, use this function to find peak MAF
def find_maf_peak(levels_maf):
    """find peak in MAF data (per segment)"""
    levels_maf_peak = {}
    for l in levels_maf:
        maf_arr = levels_maf[l]
        # in here we are filtering out low MAF levels
        if len(maf_arr)>=10:
            binsize=find_bin_width_maf(maf_arr,25)
            # adjust binsize if too big; can't find peak with big bin size
            if binsize >0.05:
                binsize = 0.05
            counts, bin_edges = numpy.histogram(maf_arr, bins=numpy.arange(0, 0.6, binsize))
            bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
            # find peaks
            # adjust for peaks that are very close together (separated by >2 binwidth)
            local_max = find_peaks(counts,distance=2)
            maxima = []
            # go through the indices and append bin_center and count
            for i in local_max[0]:
                maxima.append((bin_centers[i], counts[i]))
            # need to filter out large peak in the bottom if there is a peak around 0.5/0.33
            # get the peak with greatest MAF, at least half of max peak
            maxima.sort(key=lambda a: -a[1])
            max_peak_size = maxima[0][1]
            maxima.sort(key=lambda a: -a[0])
            for p in maxima:
                if p[1]>0.5*max_peak_size:
                    levels_maf_peak[l] = p[0]
                    break

    return levels_maf_peak

def read_maf_file(maf_file,threshold=0.1):
    # read maf file, get maf data for each segment; also calculate callable percentages
    seg_maf_pos = defaultdict(list)
    seg_sites = defaultdict(int)
    with open(maf_file, 'r') as pf:
        for line in pf:
            if line[:7]=='WARNING':
                continue
            else:
                fs = line.split('\t')
                major = int(fs[2])
                minor = int(fs[3])
                maf = minor/(major+minor)
                loc=int(fs[1])
                chr=fs[0]
                # some maf files use accession numbers instead
                if chr[:3]!='chr':
                    chr = t2t_chrom[chr]
                new_name = get_segment_pos(chr, loc)
                if new_name:
                    new_pos = new_name[1]
                    seg = new_name[0]
                    seg_sites[seg]+=1
                    if maf>threshold:
                        # transform location on chr to location in region
                        bool = numpy.array(seg_region_start[seg])>new_pos
                        if numpy.any(bool):
                            index = numpy.where(numpy.array(seg_region_start[seg])>new_pos)[0][0] -1
                        else:
                            # deal with bases in the last region of a chr
                            index = len(seg_region_start[seg])-1
                        # new location equals
                        seg_loc = seg_region_cumsum[seg][index]+ new_pos - seg_region_start[seg][index]
                        #print(fs, new_loc)
                        seg_maf_pos[seg].append([seg_loc, maf])
    seg_ratio = numpy.array([seg_sites[s]/seg_bases[s] for s in seg_sites])
    return seg_maf_pos, seg_ratio


def cov_and_maf2_file(outfile,bins,seg_maf_pos,binsize,trim_percentile=95,threshold=0.1):
    # Create two subplots that share the same x-axis
    fig, (cov, maf_plot) = plt.subplots(2, 1, sharex=True, figsize=(20, 10))
    # plot cov
    medians = {}
    chr_locs = []
    x = 0
    vertical_line = []
    for c in bins:
        med = numpy.median(bins[c])
        medians[c]=med
        cov.scatter([x+i for i in range(bins[c].shape[0]) if bins[c][i] != 2**16-1], [b for b in bins[c] if b != 2**16-1], color="blue", marker='.')
        cov.plot([x, x+bins[c].shape[0]-1],[med,med], color="red", linestyle="solid", linewidth=3)
        chr_locs.append(x+bins[c].shape[0]/2)
        x += bins[c].shape[0]
        cov.plot([x+5, x+5], [0, 2**16-1], color="black") # adjust this to see if it matches my maf plot
        vertical_line.append(x+5)
        x += 10
    ymax = max(numpy.percentile(bins[c], trim_percentile) for c in bins)

    # side histogram
    a = numpy.concatenate([bins[c] for c in bins])
    hist = numpy.histogram(a, bins=range(0, a.max(), binsize))
    cov.plot(hist[0], [hist[1][i]+(hist[1][i+1]-hist[1][i])//2 for i in range(hist[0].shape[0])], color="orange")

    cov.set_xticks(chr_locs, [c for c in bins], rotation='vertical')
    cov.set_ylim((0,ymax))
    cov.set_xlim((0,x)) # adjust for extra Y space
    cov.set_xlabel("Chromosome")
    cov.set_ylabel("Reads per Mbp")
    cov.set_title("Karyotype")
    cov.tick_params(labelbottom=True)

    # adjusted for segment plot width
    plot_start_pos = numpy.concatenate(([0],numpy.cumsum([seg_plot_width[c] for c in seg_plot_width])))
    # build seg index
    seg_ind = {c: i for i, c in enumerate(seg_plot_width)}

    x_values = []
    for c in seg_maf_pos:
        ind = seg_ind[c]
        for i in range(len(seg_maf_pos[c])):
            # adjust to plot width
            seg_pos_adj = seg_maf_pos[c][i][0]/seg_conversion_factor[c]
            # adjust for segment position
            x_values.append(plot_start_pos[ind]+seg_pos_adj+10*ind) # add gap space between segments
    # where to put seg labels
    seg_locs = []
    for i in seg_plot_width:
        seg_locs.append(plot_start_pos[seg_ind[i]] + seg_plot_width[i]/2 + seg_ind[i]*10)
    seg_locs = numpy.array(seg_locs)

    # getting MAF
    y_values = []
    for c in seg_maf_pos:
        for i in range(len(seg_maf_pos[c])):
            y_values.append(seg_maf_pos[c][i][1])

    maf_plot.scatter(x_values, y_values, color="blue", marker='.')
    # vertical lines
    for i in vertical_line:
        maf_plot.plot([i, i], [0, 1], color="black")

    maf_plot.set_xticks(seg_locs, [c for c in seg_plot_width], rotation='vertical')
    maf_plot.set_ylim((threshold,0.6))
    maf_plot.set_xlabel("Chromosome")
    maf_plot.set_ylabel("MAF")

    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Digital karyotyping from BAM or binned depth file")
    parser.add_argument("input", help="BAM/SAM of reads aligned to T2T or TSV of binned sequencing depth", nargs='+')
    parser.add_argument("--repeats", help="RepeatMasker BED file")
    parser.add_argument("--adaptive", help="Adaptive sequencing BED file")
    parser.add_argument("--maf", help="MAF file")
    parser.add_argument("--outdir", help="Output directory/prefix")
    parser.add_argument("--read_lim", help="Maximum reads to process (default is all)", type=int)
    parser.add_argument("--json", help="Output JSON format (default: TSV)", action="store_true", default=False)
    parser.add_argument("-v", help="Verbose", action="store_true", default=False)
    args = parser.parse_args()

    warnings = []
    used_maf = True
    Blast_ratio = None

    t2t = Reference(T2T_accessions) # no covg data herein, just for reference
    for a in args.input:
        if a.endswith(".tsv"):
            bins = parse_tsv(a)
        elif a.endswith(".bam") or a.endswith(".sam"):
            if args.repeats is None:
                sys.stderr.write("ERROR: --repeats required if using a BAM file\n")
                sys.exit(1)
            # each get its own copy of the reference with empty covg lists
            bins = from_bam(a, args.repeats, Reference(T2T_accessions), args.read_lim, args.outdir, args.v)
        else:
            sys.stderr.write("ERROR: unrecognized extension '{a[a.rindex('.'):]}'\n")
            sys.exit(1)
        # check AS file
        if args.adaptive is None:
                sys.stderr.write("ERROR: --adaptive required\n")
                sys.exit(1)
        else:
            # read AS file
            seg_region_start,seg_region_len,seg_region_cumsum,seg_total_len,seg_conversion_factor = read_AS(args.adaptive)
        # check MAF file
        if args.maf is None:
                sys.stderr.write("ERROR: --maf required\n")
                sys.exit(1)

        trim_percentile = 95
        binsize = round(find_bin_width2(bins,trim_percentile=trim_percentile,alpha=0.12))

        # check if sample has too much variation in coverage
        if within_segment_spread(bins)>0.08:
            warning = "Too much variation in coverage to call confidently."
            sys.stderr.write(warning+'\n')
            warnings.append(warning)

        levels=find_levels(bins,trim_percentile,binsize)
        medians = [numpy.median(bins[c]) for c in bins]
        # read maf data
        seg_maf_pos, seg_ratios = read_maf_file(args.maf)
        if numpy.sum(seg_ratios>0.01)<20: # a more strigent maf filter. at least 20 seg with >1% variable sites callable (enough coverage)
            warning = 'Not enough MAF data for karyotyping!'
            sys.stderr.write(warning+'\n')
            warnings.append(warning)
            cns, Blast_ratio = digital_karyotype(bins, trim_percentile, binsize, args.v)
            used_maf = False
        else:
            levels_maf = find_levels_maf(seg_maf_pos, levels, bins)
            levels_maf_peaks = find_maf_peak(levels_maf)
            chrom_names = [c for c in bins]
            cns_results = predict_karyo_maf(levels, levels_maf_peaks, medians, chrom_names, args.v)
            if cns_results is None:
                warning = 'Not enough minor allele frequency data for confident automated karyotype estimation.'
                sys.stderr.write(warning+'\n')
                warnings.append(warning)
                cns, Blast_ratio = digital_karyotype(bins, trim_percentile, binsize, args.v)
                used_maf = False
            else:
                cns = cns_results[0]
                Blast_ratio = cns_results[1]
        # plot
        if args.outdir is not None:
            cov_and_maf2_file(f"{args.outdir}/karyotype.png",bins,seg_maf_pos,binsize,trim_percentile)

        for c in cns:
            sys.stderr.write(f"{c}: {cns[c]}\n")
        if args.json:
            import json
            res = {f"{c}":{'n':f"{cns[c]}", 'x':f"{numpy.median(bins[c])}"} for c in cns}
            res["ISCN_karyotype"] = f"{ISCN_string(cns)}"
            res["karyotype_string"] = f"{karyotype_string(cns)}"
            res["reads_aligned"] = f"{sum([numpy.sum(bins[c]) for c in bins])}"
            if len(warnings)>0:
                res['warning'] = '\n'.join([i for i in warnings])
            if Blast_ratio is not None:
                res['blast_ratio'] = f"{Blast_ratio:.4f}"
            print(json.dumps(res))

        else:
            if len(warnings)>0:
                for i in warnings:
                    print(i)
            print("sample\t" + "\t".join(c for c in cns) + "\tISCN")
            print(f"{a}\t" + "\t".join(str(cns[c]) for c in cns) + f"\t{ISCN_string(cns)}")
        sys.stderr.write(f"{ISCN_string(cns)}\n")

