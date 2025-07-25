import sys
import json
import argparse
import datetime
from docx import Document
from docx.shared import Inches, RGBColor

class Reference:
    def __init__(self, chrom_list):
        self.chroms = chrom_list
        self.chrom_names = ["chr"+str(i+1) for i in range(22)] + ["chrX", "chrY"]
        self.chr_idx = {self.chroms[c]:c for c in range(len(self.chroms))}
        for c in range(len(self.chrom_names)):
            self.chr_idx[self.chrom_names[c]] = c
        self.chr_len = [None for c in self.chroms]
        self.chr_covg = [None for c in self.chroms]

t2t = Reference(["NC_060925.1", "NC_060926.1", "NC_060927.1", "NC_060928.1", "NC_060929.1", "NC_060930.1", "NC_060931.1", "NC_060932.1", "NC_060933.1", "NC_060934.1", "NC_060935.1", "NC_060936.1", "NC_060937.1", "NC_060938.1", "NC_060939.1", "NC_060940.1", "NC_060941.1", "NC_060942.1", "NC_060943.1", "NC_060944.1", "NC_060945.1", "NC_060946.1", "NC_060947.1", "NC_060948.1"])

def add_text(p, txts, add_newline=True):
    if type(p) == type(Document()):
        p = p.add_paragraph()
    for t in range(len(txts)):
        bold = False
        italic = False
        color = None # default
        if type(txts[t]) == type((0,)):
            txt = txts[t][0]
            bold = txts[t][1]
            if len(txts[t]) >= 3:
                italic = txts[t][2]
            if len(txts[t]) >= 4: # color included as a hex string ex. #FF00AA
                color = txts[t][3]
        else:
            txt = txts[t]
        r = p.add_run(txt)
        r.bold = bold
        r.italic = italic
        if color is not None:
            r.font.color.rgb = RGBColor(
                    int.from_bytes(bytes.fromhex(color[1:3]), byteorder=sys.byteorder),
                    int.from_bytes(bytes.fromhex(color[3:5]), byteorder=sys.byteorder),
                    int.from_bytes(bytes.fromhex(color[5:7]), byteorder=sys.byteorder)
            )
    if add_newline:
        r.add_break()

def bold(cells):
    for c in cells:
        for p in c.paragraphs:
            for r in p.runs:
                r.font.bold = True

def main(run_dir, outfile):
    document = Document()
    sections = document.sections
    margin = 0.5 # inches
    for section in sections:
        section.top_margin = Inches(margin)
        section.bottom_margin = Inches(margin)
        section.left_margin = Inches(margin)
        section.right_margin = Inches(margin)

    document.add_heading('Adaptive Whole Genome Sequencing Report', 0)

    p = document.add_paragraph()
    add_text(p, [("FOR RESEARCH USE ONLY", True, False, "#FF0000")])

    qf = json.loads(open(f"{run_dir}/fusions.json").read())
    k = json.loads(open(f"{run_dir}/karyotype.json").read())

    now = datetime.datetime.now(datetime.timezone.utc)
    add_text(document, [("Date/time of report: ", True), f"{now.strftime('%a, %d %b %Y %X UTC')}"])

    document.add_heading("Specimen details", level=4)
    p = document.add_paragraph()
    add_text(p, [("Specimen ID: ", True)])
    add_text(p, [("Date of specimen collection: ", True)])
    add_text(p, [("Specimen source: ", True)], False)

    if "meta" in qf and "run_start_time" in qf["meta"]:
        document.add_heading("Sequencing details", level=4)
        p = document.add_paragraph()
        add_text(p, [("Date/time of sequencing start: ", True), qf["meta"]["run_start_time"]])
        add_text(p, [("Sequencing run ID: ", True), qf["meta"]["run_id"]])
        add_text(p, [("Basecalling model: ", True), qf["meta"]["basecall_model"]])
        add_text(p, [("Library ID: ", True), qf["meta"]["library_id"]])
        add_text(p, [("Sequencer ID: ", True), qf["meta"]["sequencer_id"]])
        add_text(p, [("Flow cell ID: ", True), qf["meta"]["flow_cell_id"]], False)

    document.add_heading("Results", level=1)
    document.add_heading("Quality control metrics:", level=3)

    document.add_paragraph("Blast percentage reported by hematopathology: ", style="List Bullet")
    document.add_paragraph(f"Genome-wide reads aligned: {int(float(k['reads_aligned'])):,}", style="List Bullet")
    document.add_paragraph(f"Mean on-target read length: {int(qf['qc']['nt_on_target'] / qf['qc']['reads_on_target'])} nt", style="List Bullet")
    document.add_paragraph(f"Mean coverage over target regions: {qf['qc']['nt_on_target']/qf['qc']['target_regions_nt']:.2f}x", style="List Bullet")
    #document.add_paragraph(f"Relative enrichment: {qf['qc']['nt_on_target']/qf['qc']['target_regions_nt']/(qf['qc']['nt_aligned']/3.1e9):.2f}x", style="List Bullet")


    # ----------------------- Translocations ----------------------------
    document.add_heading('Structural Variants', level=1)
    # qf["fusions"] = {... "CBFA2T3-GLIS2": {"supporting_reads": 19, "duplicate_reads": 1, "repetitive_reads": 2} ...}

    document.add_heading('Fusions (by genomic translocation or deletion)', level=3)
    table = document.add_table(rows=1, cols=6)
    hdr_cells = table.rows[0].cells
    hdr_cells[0].text = 'Gene 1'
    hdr_cells[1].text = 'Gene 2'
    hdr_cells[2].text = 'Deduplicated reads'
    hdr_cells[3].text = 'Repetitive reads'
    hdr_cells[4].text = 'Supporting reads'
    hdr_cells[5].text = 'Breakpoint(s)(?)'
    bold(hdr_cells)
    for fus in qf["fusions"]:
        if qf["fusions"][fus]["supporting_reads"]-qf["fusions"][fus]["duplicate_reads"] < 3:
            continue
        row_cells = table.add_row().cells
        row_cells[0].text = fus[:fus.index('-')]
        row_cells[1].text = fus[fus.index('-')+1:]
        row_cells[2].text = f'{qf["fusions"][fus]["supporting_reads"]-qf["fusions"][fus]["duplicate_reads"]}'
        row_cells[3].text = f'{qf["fusions"][fus]["repetitive_reads"]}'
        row_cells[4].text = f'{qf["fusions"][fus]["supporting_reads"]-qf["fusions"][fus]["duplicate_reads"]-qf["fusions"][fus]["repetitive_reads"]}'
        if qf["fusions"][fus]["supporting_reads"]-qf["fusions"][fus]["duplicate_reads"]-qf["fusions"][fus]["repetitive_reads"] >= 3:
            bold([row_cells[4]])
    table.style = "Table Grid"
    document.add_paragraph("")


    # ----------------------- Karyotype ----------------------------

    document.add_heading('Karyotype', level=1)
    p = document.add_paragraph()
    add_text(p, ["ISCN karyotype: ", (f"{k['ISCN_karyotype']}", True)])
    document.add_picture(f"{run_dir}/karyotype.png", width=Inches(7))

    table = document.add_table(rows=1, cols=3)
    hdr_cells = table.rows[0].cells
    hdr_cells[0].text = 'Chromosome'
    #hdr_cells[1].text = 'Accession ID'
    hdr_cells[1].text = 'Copy number'
    hdr_cells[2].text = 'Median depth (reads/Mbp)'
    bold(hdr_cells)
    for chrom in k:
        if chrom in ["ISCN_karyotype", "reads_aligned"]:
            continue
        row_cells = table.add_row().cells
        #row_cells[0].text = t2t.chrom_names[t2t.chroms.index(chrom)]
        row_cells[0].text = chrom
        row_cells[1].text = f'{k[chrom]["n"]}'
        row_cells[2].text = f'{k[chrom]["x"]}'
    table.style = "Table Grid"

    p = document.add_paragraph()
    add_text(p, [""])
    add_text(p, [("Test description", True)])
    add_text(p, [f"Relative copy number across the genome at the chromosome level ('digital karyotype') is inferred based on relative sequencing depth by assessing genome-wide and chromosome-level depth of coverage. To avoid the potentially confounding effect of adaptive sampling on relative sequencing depth assessment, coarse-scale depth is constructed as a function of reads per million base pairs (Mbp), where each read contributes a count of one to the bin in which the center of the read aligns. A pairwise Kolmogorov-Smirnov test is applied with a threshold of p < 1e-9 and a minimum median divergence of 20% to determine which chromosomes exhibit significantly divergent copy number."], False)


    # ----------------------- FLT3-ITD ----------------------------

    document.add_heading('Internal/partial tandem duplications', level=1)

    document.add_heading('FLT3 interal tandem duplication (ITD)', level=3)
    size0 = 329 - len("GCAATTTAGGTATGAAAGCCAGC") # based on the reference
    ct0 = 1
    size1 = 0
    ct1 = 0
    document.add_paragraph(f"Normal region size (bp): {size0}")
    document.add_paragraph(f"Predicted abnormal region size (bp): {size1}")
    document.add_paragraph(f"Ratio of abnormal to normal reads (allelic ratio/AR): {ct1/ct0}")
    document.add_picture(f"{run_dir}/flt3_itd.png", width=Inches(4))

    p = document.add_paragraph()
    add_text(p, [""])
    add_text(p, [("Test description", True)])
    add_text(p, [f"We recapitulate commonly used capillary electrophoresis methods to characterize FLT3-ITD (Kiyoi et al., 1999). Read segments bounded by FLT3 11F (GCAATTTAGGTATGAAAGCCAGC) and 12R (CTTTCAGCATTTTGACGGCAACC) primers are identified and the inter-primer fragment size is defined as the last aligning nucleotide of one primer to the last aligning nucleotide of the other. The distribution of fragment sizes shown above is used to define normal/abnormal fragment populations and the respective allelic ratio between them."], False)


    # ----------------------- UBTF-ITD ----------------------------

    document.add_heading('UBTF interal tandem duplication (ITD)', level=3)
    size0 = 617 - len("CCAAGAAGCCAGCCCAGGAAGG") # based on the reference
    ct0 = 1
    size1 = 0
    ct1 = 0
    document.add_paragraph(f"Normal region size (bp): {size0}")
    document.add_paragraph(f"Predicted abnormal region size (bp): {size1}")
    document.add_paragraph(f"Ratio of abnormal to normal reads (allelic ratio/AR): {ct1/ct0}")
    document.add_picture(f"{run_dir}/ubtf_itd.png", width=Inches(4))
    p = document.add_paragraph()
    add_text(p, [""])
    add_text(p, [("Test description", True)])
    add_text(p, [f"We recapitulate capillary electrophoresis method to characterize UBTF internal tandem duplication (Duployez et al., 2023). Read segments bounded by UBTF_Forward (CCAAGAAGCCAGCCCAGGAAGG) and UBTF_Reverse (GCCTCTCGGGCCTTGTACTTGG) primers are identified and the inter-primer fragment size is defined as the last aligning nucleotide of one primer to the last aligning nucleotide of the other. The distribution of fragment sizes shown above is used to define normal/abnormal fragment populations and the respective allelic ratio between them."], False)


    # ----------------------- KMT2A-ITD ----------------------------

    #document.add_heading('KMT2A partial tandem duplication (PTD)', level=3)
    #document.add_paragraph("TBD")


    # ----------------------- SNVs / genotypes ----------------------------

    gens = json.loads(open(f"{run_dir}/genotypes.json").read())
    # coverage, mutations, phase, genotype
    document.add_heading('SNV Genotypes', level=3)
    for g in gens:
        p = document.add_paragraph()
        add_text(p, [""])
        add_text(p, [(f"{g}", True)])
        add_text(p, [f"Sequencing depth: {gens[g]['coverage']:.2f}x"])
        add_text(p, [f"Mutations:"])
        for m in gens[g]['mutations']:
            add_text(p, [f"  {m}: {gens[g]['mutations'][m]} AF"])
        if len(gens[g]['mutations']) == 0:
            add_text(p, ["(none detected)"])
        add_text(p, [f"Genotype: ", (f"{gens[g]['genotype']}", True)], False)


    '''

Small deletions / insertions.
Gene
Gain (number of copies) / Loss (hemizygous or homozygous)
Deletion insertion coordinates

Single nucleotide variants:
Gene, mutation (protein and genomic), VAF

Internal tandem duplication / partial tandem duplication:
Gene, length, allelic ratio.
Test Explanation: (see validation plan for list of reportable fusions)

Karyotype:
Report down to the chromosomal arm level for the karyotype
Fusion detection:
Methods
Reportable alteration list
ETV::RUNX1
TCF3:PBX1
KMT2Ar::***
MEF2D
ZNF384
Small deletions / insertions:
Methods
Reportable alteration list
iAMP21
IKZF1, PAX5, CDKN2A, CDKN2B, Xp22.33/Yp11.31 (PAR1 region), ERG
RUNX1
ETV6, RB1, BTG1, EBF1
SNV:
Methods
Reportable alteration list:
B-ALL
PAX5 (P80R), IKZF1 (N159Y)
Pharmacogenomic
TPMT, NUDT15
Tandem Duplications:
Methods
Reportable alteration list
FLT3, KMT2A, UBTF




Test Limitations:

Targeted:

Depth:

Low VAF:

Low blast percentage:
    '''

    document.save(outfile)

if __name__ == "__main__":
    parser = argparse.ArgumentParser("Make WGS report (docx) from analysis logs/directory")
    parser.add_argument("dir", help="Analysis output directory")
    parser.add_argument("out", help="Output log file (.docx)")
    args = parser.parse_args()
    main(args.dir, args.out)
