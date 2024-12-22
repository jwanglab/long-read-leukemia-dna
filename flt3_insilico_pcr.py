import argparse
import sys
import edlib
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from collections import defaultdict


# https://ashpublications.org/blood/article/93/9/3074/266503/Prognostic-Implication-of-FLT3-and-N-RAS-Gene
primer_11f = "GCAATTTAGGTATGAAAGCCAGC"
primer_12r = "CTTTCAGCATTTTGACGGCAACC"

try:
  from string import maketrans
  COMPL = maketrans("ATGC","TACG")
except: # python 3
  COMPL = str.maketrans("ATGC","TACG")

def rc(seq):
  return seq.translate(COMPL)[::-1]

def main(fq, out):
    n = 0
    hist = defaultdict(int)
    for line in open(fq):
        if n%4 == 0:
            name = line[1:-1]
        if n%4 == 1:
            s = line.strip()
            pos = []
            for primer in (primer_11f, primer_12r):
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
                l = abs(pos[0][1] - pos[1][1])
                hist[l//10] += 1
                print(name, l)
        n += 1
    mx = max(k for k in hist)+5
    hist_array = [0]*(mx+1)
    for k in hist:
        hist_array[k] = hist[k]
    plt.plot([a*10 for a in range(len(hist_array))], hist_array)
    print(hist_array)
    plt.ylabel("# reads")
    plt.xlabel("$\it{FLT3}$ 11F/12R amplicon size (nt)")
    plt.savefig(out)

if __name__ == "__main__":
    parser = argparse.ArgumentParser("Predict FLT3-ITD and allelic ratio from targeted nanopore sequencing")
    parser.add_argument("fastq", help="Reads in FASTQ format")
    parser.add_argument("out", help="Output image file for pseudo-capillary electrophoresis sizing of FLT3 amplicons")
    args = parser.parse_args()
    main(args.fastq, args.out)
