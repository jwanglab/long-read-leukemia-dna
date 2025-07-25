import json
import argparse
import sys
import edlib
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from collections import defaultdict
import pysam


# https://www.nature.com/articles/s41375-023-01906-z#MOESM2
primer_fw = "CCAAGAAGCCAGCCCAGGAAGG"
primer_rv = "GCCTCTCGGGCCTTGTACTTGG"

try:
  from string import maketrans
  COMPL = maketrans("ATGC","TACG")
except: # python 3
  COMPL = str.maketrans("ATGC","TACG")

def rc(seq):
  return seq.translate(COMPL)[::-1]

def main(bam_files, out):
    result = {
        "reads": {},
        "length_histogram": {}
    }
    hist = defaultdict(int)
    for bam_file in bam_files:
        try:
            af = pysam.AlignmentFile(bam_file, mode=('rb' if bam_file.endswith(".bam") else 'r'), ignore_truncation=True)
        except FileNotFoundError as e:
            sys.stderr.write(f"Error: File not found: {bam_file}\n")
            sys.exit(1)
        for a in af.fetch(until_eof=True):
            s = a.query_sequence
            if s is None or len(s) < 20:
                continue
            pos = []
            for primer in (primer_fw, primer_rv):
                res = edlib.align(primer, s, "HW")
                ed = res["editDistance"]
                if ed < len(primer)*0.1:
                    pos.append(res["locations"][0])
                else:
                    primer = rc(primer)
                    res = edlib.align(primer, s, "HW")
                    ed = res["editDistance"]
                    if ed < len(primer)*0.1:
                        pos.append(res["locations"][0])
            if len(pos) > 1:
                pos.sort(key=lambda a: a[1])
                l = pos[1][1] - pos[0][1]
                hist[l//10] += 1
                result["reads"][a.query_name] = l
    mx = max(k for k in hist)+5
    hist_array = [0]*(mx+1)
    for k in hist:
        hist_array[k] = hist[k]
    plt.plot([a*10 for a in range(len(hist_array))], hist_array)
    result["length_histogram"] = {a*10:hist_array[a] for a in range(min(len(hist_array), 80))}
    plt.ylabel("# reads")
    plt.xlabel("$\it{UBTF}$ F/R amplicon size (nt)")
    plt.xlim((400,800))
    plt.savefig(out)
    sys.stdout.write(json.dumps(result)+"\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser("Predict UBRF-TD and allelic ratio from targeted nanopore sequencing")
    parser.add_argument("out", help="Output image file for pseudo-capillary electrophoresis sizing of UBTF amplicons")
    parser.add_argument("bams", help="BAM of aligned reads (minimap2 -cx map-ont)", nargs='+')
    args = parser.parse_args()
    main(args.bams, args.out)
