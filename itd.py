import json
import argparse
import sys
import pysam

flt3_hg38 = ("NC_000013.11", 28003274-50000, 28100576+50000)
flt3_t2t = ("NC_060937.1", 27225623, 27323127)
flt3_t2t_exon14= ("NC_060937.1", 27256456, 27256588)
flt3_t2t_exon15 = ("NC_060937.1", 27256261, 27256365)
flt3_t2t_focal = (flt3_t2t[0], flt3_t2t_exon15[1], flt3_t2t_exon14[2])
flt3_chr_t2t = ("chr13", flt3_t2t_focal[1], flt3_t2t_focal[2])

ubtf_t2t_exon13 = ("NC_060941.1", 45063895, 45064050)
ubtf_chr_t2t = ("chr17", ubtf_t2t_exon13[1], ubtf_t2t_exon13[2])

def main(bam_files):
    ins_map = [("FLT3", {}, flt3_t2t_focal, flt3_chr_t2t, [0]*(flt3_t2t_focal[2]-flt3_t2t_focal[1])), ("UBTF", {}, ubtf_t2t_exon13, ubtf_chr_t2t, [0]*(ubtf_t2t_exon13[2]-ubtf_t2t_exon13[1]))]

    for bam_file in bam_files:
        try:
            af = pysam.AlignmentFile(bam_file, mode=('rb' if bam_file.endswith(".bam") else 'r'), ignore_truncation=True)
        except FileNotFoundError as e:
            sys.stderr.write(f"Error: File not found: {bam_file}\n")
            sys.exit(1)
        for name, insertions, t2t_region, chr_region, covg in ins_map:
            try:
                ch, st, en = t2t_region
                region = af.fetch(ch, st, en)
            except Exception:
                ch, st, en = chr_region
                region = af.fetch(ch, st, en)
            for a in region:
                # walk CIGAR, looking for insertions
                q = a.query_alignment_start
                t = a.reference_start
                val = ""
                for c in a.cigarstring:
                    if c in "MIDNSHP=X": # op
                        val = int(val)
                        if c in "M=X":
                            q += val
                            if t <= en and t+val >= st:
                                for i in range(max(st, t), min(t+val, en)):
                                    covg[i-st] += 1;
                            t += val
                        elif c in "IS":
                            if t >= st and t < en:
                                if t not in insertions:
                                    insertions[t] = {}
                                if val not in insertions[t]:
                                    insertions[t][val] = 0
                                insertions[t][val] += 1
                            q += val
                        elif c in "DN":
                            if t <= en and t+val >= st:
                                for i in range(max(st, t), min(t+val, en)):
                                    covg[i-st] += 1;
                            t += val
                        else: # H (hard clip), do nothing
                            pass
                        val = ""
                    else:
                        val += c

    result = {}
    for name, insertions, t2t_region, _, covg in ins_map:
        ch, st, en = t2t_region
        result[name] = []
        for t in sorted(t for t in insertions):
            cts = [(a, insertions[t][a]) for a in insertions[t] if a >= 3 and insertions[t][a] > 1]
            if len(cts) == 0: # all <3nt or 1 read, skip it
                continue
            for l in sorted(l for l in insertions[t]):
                if l >= 3 and insertions[t][l] > 1:
                    #print(f"pos {t} insertion of {l} nt, {insertions[t][l]} times")
                    result[name].append({
                            "position": t,
                            "length": l,
                            "exact": insertions[t][l],
                            "merged": insertions[t][l],
                            "coverage": covg[t-st]
                    })
        # merge insertions within <ins size+10%> by length and distance, they may the same event
        for t in sorted(t for t in insertions):
            for l in sorted(l for l in insertions[t]):
                for ins in result[name]:
                    if abs(t-ins["position"]) < ins["length"] * 1.1 and l*0.9 < ins["length"] < l*1.1 and (t!=ins["position"] or l!=ins["length"]):
                        ins["merged"] += insertions[t][l]

    print(json.dumps(result))


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Predict FLT3 and UBTF ITD and allelic ratio from targeted nanopore sequencing")
    parser.add_argument("bams", help="BAM of aligned reads (minimap2 -cx map-ont)", nargs='+')
    args = parser.parse_args()
    main(args.bams)
