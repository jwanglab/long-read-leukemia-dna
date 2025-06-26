import sys
import argparse
from collections import defaultdict
import pysam
import pyfastx
import numpy as np

class Reference:
    def __init__(self, chrom_list):
        self.chroms = chrom_list
        self.chrom_names = ["chr"+str(i+1) for i in range(22)] + ["chrX", "chrY"]
        self.chr_idx = {self.chroms[c]:c for c in range(len(self.chroms))}
        for c in range(len(self.chrom_names)):
            self.chr_idx[self.chrom_names[c]] = c
        self.chr_len = [None for c in self.chroms]
        self.chr_covg = [None for c in self.chroms]

hg38 = Reference(["NC_000001.11", "NC_000002.12", "NC_000003.12", "NC_000004.12", "NC_000005.10", "NC_000006.12", "NC_000007.14", "NC_000008.11", "NC_000009.12", "NC_000010.11", "NC_000011.10", "NC_000012.12", "NC_000013.11", "NC_000014.9", "NC_000015.10", "NC_000016.10", "NC_000017.11", "NC_000018.10", "NC_000019.10", "NC_000020.11", "NC_000021.9", "NC_000022.11", "NC_000023.11", "NC_000024.10"])
t2t = Reference(["NC_060925.1", "NC_060926.1", "NC_060927.1", "NC_060928.1", "NC_060929.1", "NC_060930.1", "NC_060931.1", "NC_060932.1", "NC_060933.1", "NC_060934.1", "NC_060935.1", "NC_060936.1", "NC_060937.1", "NC_060938.1", "NC_060939.1", "NC_060940.1", "NC_060941.1", "NC_060942.1", "NC_060943.1", "NC_060944.1", "NC_060945.1", "NC_060946.1", "NC_060947.1", "NC_060948.1"])

use_transcripts = {
    "TPMT": "rna-NM_001346817.1",
}

# this defaults to the first mRNA/CDS listed after the gene in the GFF file
#   it is presumed to be the canonical transcript, but in looks like sometimes it isn't
#   so for the time being, we are manually overwriting those when we find them (ex. TP53)
def get_gene_annotation(gene, gff):
    gene_id = None
    gene_name = None
    chrom = None
    gene_strand = None
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
            # requires specific mRNA/transcripts for specified genes
            if gene.upper() in use_transcripts and fields["ID"] != use_transcripts[gene.upper()]:
                continue
            if fields["Parent"] == gene_id:
                if gene_name == "TP53" and fields["ID"] != "rna-NM_001126112.3": # this is the canonical transcript but is listed second
                    continue
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

def translate(nt):
    aa = ""
    skipped = 0
    for i in range(0, len(nt), 3):
        if 'N' in nt[i:i+3]:
            aa += ' '
            skipped += 1
            continue
        if '-' in nt[i:i+3]:
            sys.stderr.write(f"DELETION in codon {i/3}: '{nt[i:i+3]}'\n")
            aa += '-'
            continue
        if len(nt) < i+3:
            sys.stderr.write(f"WARNING: not an even number of codons, ignoring '{nt[i:]}' at the end\n")
            break
        aa += codontab[nt[i:i+3]]
    sys.stderr.write(f"skipped {skipped} codons with Ns in them.\n")
    return aa

codontab = {
    'TCA': 'S',
    'TCC': 'S',
    'TCG': 'S',
    'TCT': 'S',
    'TTC': 'F',
    'TTT': 'F',
    'TTA': 'L',
    'TTG': 'L',
    'TAC': 'Y',
    'TAT': 'Y',
    'TAA': '*',
    'TAG': '*',
    'TGC': 'C',
    'TGT': 'C',
    'TGA': '*',
    'TGG': 'W',
    'CTA': 'L',
    'CTC': 'L',
    'CTG': 'L',
    'CTT': 'L',
    'CCA': 'P',
    'CCC': 'P',
    'CCG': 'P',
    'CCT': 'P',
    'CAC': 'H',
    'CAT': 'H',
    'CAA': 'Q',
    'CAG': 'Q',
    'CGA': 'R',
    'CGC': 'R',
    'CGG': 'R',
    'CGT': 'R',
    'ATA': 'I',
    'ATC': 'I',
    'ATT': 'I',
    'ATG': 'M',
    'ACA': 'T',
    'ACC': 'T',
    'ACG': 'T',
    'ACT': 'T',
    'AAC': 'N',
    'AAT': 'N',
    'AAA': 'K',
    'AAG': 'K',
    'AGC': 'S',
    'AGT': 'S',
    'AGA': 'R',
    'AGG': 'R',
    'GTA': 'V',
    'GTC': 'V',
    'GTG': 'V',
    'GTT': 'V',
    'GCA': 'A',
    'GCC': 'A',
    'GCG': 'A',
    'GCT': 'A',
    'GAC': 'D',
    'GAT': 'D',
    'GAA': 'E',
    'GAG': 'E',
    'GGA': 'G',
    'GGC': 'G',
    'GGG': 'G',
    'GGT': 'G'
}

def main(bam_file, gene_names, refseq, gff, phase, maf_threshold):
    refseq = pyfastx.Fasta(refseq)
    af = pysam.AlignmentFile(bam_file, mode='rb')
    if gene_names is not None and len(gene_names) > 0:
        for gene_name in gene_names:
            sys.stderr.write("\n")
            sys.stderr.write(f"------ {gene_name} ------\n")
            depth, rpm, muts = compare(bam_file, gene_name, refseq, gff, af, phase, maf_threshold)
            sys.stdout.write(f"{gene_name}\t{depth}\t" + '\t'.join(muts) + "\n")
    else:
        for gene_name in ["PAX5", "IKZF1", "CDKN2A", "CDKN2B", "CRLF2", "ETV6", "EBF1", "ERG", "RUNX1", "TPMT", "NUDT15", "FLT3", "NPM1"]:
            sys.stderr.write("\n")
            sys.stderr.write(f"------ {gene_name} ------\n")
            depth, rpm, muts = compare(bam_file, gene_name, refseq, gff, af, phase, maf_threshold)
            sys.stdout.write(f"{gene_name}\t{depth}\t" + '\t'.join(muts) + "\n")

def compare(bam_file, gene_name, refseq, gff, af, do_phase, maf_threshold, verbose=True):
    gene_id, gene_name, mrna_id, mrna_name, chrom, gene_strand, cds = get_gene_annotation(gene_name, gff)
    cds_len = 0
    if verbose:
        sys.stderr.write(f"{gene_name} {mrna_name}\n")
        sys.stderr.write(f"{chrom} {gene_strand}\n")
        sys.stderr.write("CDS:\n")
    for c in cds:
        if verbose:
            sys.stderr.write(f"{c}\n")
        cds_len += c[1]-c[0]+1
    cds.sort(key = lambda a:a[0])
    if verbose:
        sys.stderr.write(f"CDS length: {cds_len}\n")

    bam_chrom = chrom
    if chrom not in af.references:
        for ref in [hg38, t2t]:
            if chrom in ref.chroms and ref.chrom_names[ref.chroms.index(chrom)] in af.references:
                bam_chrom = ref.chrom_names[ref.chroms.index(chrom)]
                break

    reads = {} # dict of colinear CDS per read - for phasing
    covg = np.zeros((cds_len, 7)) # (covg, A, C, G, T, N, del)
    insertions = [defaultdict(int) for c in range(cds_len+1)]
    ref_cds = ""
    offset = 0
    simplex_accs = []
    duplex_accs = []
    n_reads = 0
    read_starts = 0
    for alignment in af.fetch(bam_chrom, cds[0][0]-1, cds[-1][1]):
        if alignment.reference_start >= cds[0][0]-1 and alignment.reference_start <= cds[-1][1]:
            read_starts += 1
    for c in cds:
        st = c[0]-1 # pysam pileup() uses 0-based indexing
        en = c[1]
        for alignment in af.fetch(bam_chrom, st, en):
            if not alignment.query_name in reads:
                n_reads += 1
                reads[alignment.query_name] = [' ']*cds_len
                ps = alignment.get_cigar_stats()[0]
                # pileup stats: MIDNSHP=XB + NM (edits)
                acc = (sum(ps[:10]) - ps[10]) / sum(ps[:10])
                if ';' in alignment.query_name:
                    duplex_accs.append((ps[10], sum(ps[:10]), ps[1]+ps[2])) # NM (edits), total alignment length, indels
                else:
                    simplex_accs.append((ps[10], sum(ps[:10]), ps[1]+ps[2]))
            if alignment.query_sequence is None:
                sys.stderr.write(f"read '{alignment.query_name}' has no sequence\n")
                continue
            ins = ""
            last_t = None
            for (q,t) in alignment.get_aligned_pairs():
                if (t is not None and t >= en) or (t is None and last_t is not None and last_t+len(ins)+1 >= en):
                    break;
                if t is None:
                    if last_t is not None and last_t >= st:
                        ins += " ACGTN-"[(" ACGTN" if gene_strand == '+' else " TGCAN").index(alignment.query_sequence[q])]
                else:
                    if t >= st:
                        cds_pos = offset + t - st
                        if gene_strand == '-': #alignment.is_forward != (gene_strand == '+'):
                            cds_pos = cds_len - cds_pos - 1
                        # complement here if necessary
                        #al = 6 if q is None else (" ACGTN" if alignment.is_forward == (gene_strand == '+') else " TGCAN").index(alignment.query_sequence[q])
                        al = 6 if q is None else (" ACGTN" if gene_strand == '+' else " TGCAN").index(alignment.query_sequence[q])
                        reads[alignment.query_name][cds_pos] = " ACGTN-"[al]
                        covg[cds_pos,0] += 1 # total coverage
                        covg[cds_pos,al] += 1
                        if ins != "":
                            insertions[cds_pos][ins] += 1
                            ins = ""
                    last_t = t
        ref_cds += refseq.fetch(chrom, (c[0], c[1])).upper()
        offset += (en-st)
    if gene_strand == '-':
        ref_cds = ''.join([{'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}[c] for c in ref_cds[::-1]])

    muts = []
    sys.stderr.write(f"{n_reads} reads aligned\n")
    if n_reads == 0:
        return 0, 0, []
    sys.stderr.write("\n")
    sys.stderr.write("Insertions (>30% AF):\n")
    for i in range(len(insertions)):
        for ins in insertions[i]:
            if insertions[i][ins] > covg[i,0] * 0.1:
                sys.stderr.write(f"  {i+1}: {ins} x{insertions[i][ins]} ({insertions[i][ins]/covg[i,0]:.2f} AF)\n")
            if insertions[i][ins] >= covg[i,0] * 0.3:
                muts.append(f"{i+1}: {ins} x{insertions[i][ins]} ({insertions[i][ins]/covg[i,0]:.2f} AF)\n")
    sys.stderr.write("\n")
    sys.stderr.write("Average accuracy across all reads covering this gene (including introns):\n")
    if len(simplex_accs) > 0:
        simplex_mismatches = sum(a[0] for a in simplex_accs)
        simplex_alns = sum(a[1] for a in simplex_accs)
        simplex_indels = sum(a[2] for a in simplex_accs)
        sys.stderr.write(f"  Simplex: {(simplex_alns-simplex_mismatches)/simplex_alns*100:.2f} % (excluding indels: {(simplex_alns-simplex_mismatches)/(simplex_alns-simplex_indels)*100:.2f} %)\n")
    if len(duplex_accs) > 0:
        duplex_mismatches = sum(a[0] for a in duplex_accs)
        duplex_alns = sum(a[1] for a in duplex_accs)
        duplex_indels = sum(a[2] for a in duplex_accs)
        sys.stderr.write(f"  Duplex: {(duplex_alns-duplex_mismatches)/duplex_alns*100:.2f} % (excluding indels: {(duplex_alns-duplex_mismatches)/(duplex_alns-duplex_indels)*100:.2f} %)\n")
    cons = []
    alt = []
    depth = 0
    hets = []
    homs = []
    low = 0
    for row in range(covg.shape[0]):
        als = sorted([(i,covg[row,i]) for i in range(1,covg.shape[1])], key=lambda a:a[1])
        locus_depth = sum(a[1] for a in als)
        depth += locus_depth

        if als[-1][1] < 5:
            low += 1
            if verbose:
                sys.stderr.write(f"low coverage at {row} {covg[row,:]}\n")
            cons.append(('N', 0))
            alt.append(('N', 0))
        else:
            if als[-2][1] < locus_depth * 0.1 and als[-1][1] > 1 and ' ACGTN-'[als[-1][0]] != ref_cds[row].upper(): # ~homozygous if the next highest is low
                sys.stderr.write(f"HOM {ref_cds[row].upper()} -> {' ACGTN-'[als[-1][0]]} at {row+1}\n")
                sys.stderr.write(f"  {covg[row]}\n")
                homs.append(row)
                alt.append((" ACGTN-"[als[-1][0]], als[-1][0]/locus_depth))
            elif als[-1][1] >= locus_depth * maf_threshold and als[-2][1] >= locus_depth * maf_threshold and als[-2][1] > 1: # ~heterozygous
                sys.stderr.write(f"HET {ref_cds[row].upper()} -> {' ACGTN-'[als[-1][0]]}/{' ACGTN-'[als[-2][0]]} at {row+1}\n")
                sys.stderr.write(f"  {covg[row]}\n")
                hets.append(row)
                if ' ACGTN-'[als[-1][0]] != ref_cds[row].upper():
                    alt.append((" ACGTN-"[als[-1][0]], als[-1][1]/locus_depth))
                else:
                    alt.append((" ACGTN-"[als[-2][0]], als[-2][1]/locus_depth))
            else:
                alt.append((" ACGTN-"[als[-1][0]], als[-1][1]/locus_depth))

            cons.append((" ACGTN-"[als[-1][0]], als[-1][1]/locus_depth))

    sys.stderr.write(f"{low} sites with coverage < 5\n")

    # phasing
    if do_phase:
        for i in range(len(hets)):
            if sum(covg[i]) < 10:
                continue
            for j in range(i+1, len(hets)):
                if sum(covg[j]) < 10:
                    continue
                phase = defaultdict(int)
                spanning = 0
                sys.stderr.write(f"Phasing {hets[i]+1} : {hets[j]+1}\n")
                for r in reads:
                    if reads[r][hets[i]] != ' ' and reads[r][hets[j]] != ' ':
                        spanning += 1
                        phase[f"{reads[r][hets[i]]}/{reads[r][hets[j]]}"] += 1
                sys.stderr.write(f"  spanning reads: {spanning}\n")
                for p in phase:
                    sys.stderr.write(f"  {p}: {phase[p]} ({phase[p]/spanning*100:.2f}%)\n")
    depth /= covg.shape[0]
    sys.stderr.write(f"Average depth: {depth}\n")
    rpm = read_starts / covg.shape[0] * 1000000 # reads per million nt
    if len(cons) != len(ref_cds):
        sys.stderr.write(f"WARNING: reference ({len(ref_cds)}) and consensus ({len(cons)}) CDS are not the same length!\n")
    if len(cons)%3 != 0:
        sys.stderr.write(f"WARNING: Coding sequence is not an even number of codons: {len(cons)}\n")

    sys.stderr.write("Nucleotide substitutions:\n")
    for i in range(len(cons)):
        if i >= len(ref_cds):
            sys.stderr.write(f"  template truncated at {i}: {cons[i]}\n")
        elif cons[i][0] != ref_cds[i] and cons[i][1] > 0:
            sys.stderr.write(f"  {ref_cds[i]}{i+1}{cons[i][0]} ({cons[i][1]:.2f} AF)\n")
        elif alt[i][0] != ref_cds[i] and alt[i][1] > 0:
            sys.stderr.write(f"  HET {ref_cds[i]}{i+1}{alt[i][0]} ({alt[i][1]:.2f} AF)\n")

    trans_af = [min([c[1] for c in cons[i:i+3]]) for i in range(0, len(cons), 3)]
    cons = ''.join([a[0] for a in cons])
    trans = translate(cons)

    alt_trans_af = [min([c[1] for c in alt[i:i+3]]) for i in range(0, len(alt), 3)]
    alt = ''.join([a[0] for a in alt])
    alt_trans = translate(alt)

    if trans[-1] != "*":
        sys.stderr.write(f"WARNING: AA sequence does not end in a stop codon: ...{trans[-10:]}\n")

    sys.stderr.write("Amino acid substitutions:\n")
    ref_trans = translate(ref_cds)
    for i in range(len(trans)):
        if i >= len(ref_trans):
            sys.stderr.write(f"  template truncated at {i}: {trans[i]}\n")
        elif trans[i] != ref_trans[i] and trans[i] != ' ':
            sys.stderr.write(f"  {ref_trans[i]}{i+1}{trans[i]} ({trans_af[i]:.2f} AF)\n")
            if trans_af[i] >= maf_threshold:
                muts.append(f"{ref_trans[i]}{i+1}{trans[i]} ({trans_af[i]:.2f} AF)")
        elif alt_trans[i] != ref_trans[i] and alt_trans[i] != ' ':
            sys.stderr.write(f"  HET {ref_trans[i]}{i+1}{alt_trans[i]} ({alt_trans_af[i]:.2f} AF)\n")
            if alt_trans_af[i] >= maf_threshold:
                muts.append(f"{ref_trans[i]}{i+1}{alt_trans[i]} ({alt_trans_af[i]:.2f} AF)")

    return depth, rpm, muts


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Identify coding variation in nanopore adaptive sampling BAM file")
    parser.add_argument("bam", help="BAM of reads to hg38 (mm2 -x map-ont)")
    parser.add_argument("ref", help="Reference FASTA file")
    parser.add_argument("gff", help="Matching reference annotation GFF")
    parser.add_argument("--genes", help="Gene name", nargs='+')
    parser.add_argument("--maf", help="Minimum MAF to report (default: 0.3)", type=float, default=0.3)
    parser.add_argument("--phase", help="Perform variant phasing where possible", action="store_true", default=False)
    args = parser.parse_args()
    main(args.bam, args.genes, args.ref, args.gff, args.phase, args.maf)
