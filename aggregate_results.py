import sys
import json
import argparse
import datetime

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


def main(run_dirs, version):

    aggregate_result = []

    for run_dir in run_dirs:
        name = run_dir[run_dir.rindex('/')+1:run_dir.index('.')]

        # columns in tabular output
        result = {
            "name": name,
            "fusions": "",
            "one_sided_fusions": ""
        }

        try:
            qf = json.loads(open(f"{run_dir}/fusions.json").read())
            k = json.loads(open(f"{run_dir}/karyotype.json").read())
        except Exception as e:
            sys.stderr.write(f"Error reading {run_dir}/fusions.json -- skipping this sample\n")
            continue

        blacklist = [("CYP2C19", "CYP2C9"), ("BCR", "IGL"), ("IGK", "IGL"), ("IGK", "KMT2C")]
        filtered_fusions = []
        for fus in qf["fusions"]:
            g0 = fus["gene1"]
            g1 = fus["gene2"]
            if "name" in g0 and "name" in g1 and ((g0["name"], g1["name"]) in blacklist or (g1["name"], g0["name"]) in blacklist):
                continue
            if fus["supporting_reads"] < 3:
                continue
            for b in fus["breakpoints"]:
                if b["n_reads"] > 2 or "DUX4" in g0["name"] or "DUX4" in g1["name"]: # DUX4 is excused from the breakpoint filter
                    filtered_fusions.append(fus)
                    break

        two_sided = [f for f in filtered_fusions if len(f["gene1"]["name"]) > 0 and len(f["gene1"]["name"]) > 0]
        one_sided = [f for f in filtered_fusions if len(f["gene1"]["name"]) == 0 or len(f["gene1"]["name"]) == 0]

        if len(two_sided) > 0:
            for fus in two_sided:
                g0 = fus["gene1"]
                g1 = fus["gene2"]

                g0_txt = f"{g0['name'] if 'name' in g0 else ''} ({g0['chr']} ~{g0['pos']/100:.2f}Mbp)"
                g1_txt = f"{g1['name'] if 'name' in g1 else ''} ({g1['chr']} ~{g1['pos']/100:.2f}Mbp)"
                n_reads = f'{fus["supporting_reads"]}'
                breakpoints = fus["breakpoints"]
                breakpoint_text = '; '.join([f'{breakpoints[i]["gene0_chr"]}:{breakpoints[i]["gene0_pos"]}{breakpoints[i]["gene0_dir"]} -- {breakpoints[i]["gene1_chr"]}:{breakpoints[i]["gene1_pos"]}{breakpoints[i]["gene1_dir"]} ({breakpoints[i]["n_reads"]} reads)' for i in range(len(breakpoints))]);

                result["fusions"] += f"{g0_txt} :: {g1_txt} ({n_reads} reads); "


        if len(one_sided) > 0:
            for fus in one_sided:
                g0 = fus["gene1"]
                g1 = fus["gene2"]

                g0_txt = f"{g0['name'] if 'name' in g0 else ''} ({g0['chr']} ~{g0['pos']/100:.2f}Mbp)"
                g1_txt = f"{g1['name'] if 'name' in g1 else ''} ({g1['chr']} ~{g1['pos']/100:.2f}Mbp)"
                n_reads = f'{fus["supporting_reads"]}'
                breakpoints = fus["breakpoints"]
                breakpoint_text = '; '.join([f'{breakpoints[i]["gene0_chr"]}:{breakpoints[i]["gene0_pos"]}{breakpoints[i]["gene0_dir"]} -- {breakpoints[i]["gene1_chr"]}:{breakpoints[i]["gene1_pos"]}{breakpoints[i]["gene1_dir"]} ({breakpoints[i]["n_reads"]} reads)' for i in range(len(breakpoints))]);

                result["one_sided_fusions"] += f"{g0_txt} :: {g1_txt} ({n_reads} reads); "


        # ----------------------- Karyotype ----------------------------

        result["warnings"] = ';'.join(k['warning'].split('\n') if 'warning' in k else [])

        result["karyotype"] = f"{k['karyotype_string']}"
        result["ISCN"] = f"{k['ISCN_karyotype']}"


        # ----------------------- CNVs ----------------------------

        cnvs = json.loads(open(f"{run_dir}/cnv.json").read())

        result["RUNX1"] = f"{int(round(cnvs['RUNX1']['local']))}x ({cnvs['RUNX1']['local']:.1f}x)"

        result["CDKN2A"] = f"{int(round(cnvs['CDKN2A']['focal']))}x ({cnvs['CDKN2A']['focal']:.1f}x)"
        result["CDKN2B"] = f"{int(round(cnvs['CDKN2B']['focal']))}x ({cnvs['CDKN2B']['focal']:.1f}x)"

        if "ERG" in cnvs and len(cnvs["ERG"]["deletions"]) > 0:
            result["ERG"] = ';'.join(f"{d['chrom']}: {d['start']} - {d['end']} ({d['end']-d['start']} nt): {d['reads']} supporting reads" for d in cnvs["ERG"]["deletions"])
        else:
            result["ERG"] = ""

        if "IKZF1" in cnvs and len(cnvs["IKZF1"]["deletions"]) > 0:
            result["IKZF1_deletions"] = '; '.join(f"{d['chrom']}: {d['start']} - {d['end']} ({d['end']-d['start']} nt): {d['reads']} supporting reads" for d in cnvs["IKZF1"]["deletions"])
        else:
            result["IKZF1_deletions"] = f"{int(round(cnvs['IKZF1']['focal']))}x  ({cnvs['IKZF1']['focal']:.1f}x)"

        if "PAX5" in cnvs and len(cnvs["PAX5"]["deletions"]) > 0:
            result["PAX5_deletions"] = '; '.join(f"{d['chrom']}: {d['start']} - {d['end']} ({d['end']-d['start']} nt): {d['reads']} supporting reads" for d in cnvs["PAX5"]["deletions"])
        else:
            result["PAX5_deletions"] = ""

        if "PAX5" in cnvs and len(cnvs["PAX5"]["duplications"]) > 0:
            result["PAX5_ITD"] = '; '.join(f"{d['chrom']}: {d['start']} - {d['end']} ({d['end']-d['start']} nt): {d['reads']} supporting reads" for d in cnvs["PAX5"]["duplications"])
        else:
            result["PAX5_ITD"] = ""

        if "KMT2A" in cnvs and len(cnvs["KMT2A"]["duplications"]) > 0:
           result["KMT2A_ITD"] = '; '.join(f"{d['chrom']}: {d['start']} - {d['end']} ({d['end']-d['start']} nt): {d['reads']} supporting reads" for d in cnvs["KMT2A"]["duplications"])
        else:
            result["KMT2A_ITD"] = ""


        # ----------------------- ITD ----------------------------

        itds = json.loads(open(f"{run_dir}/itd.json").read())

        result["FLT3_ITD"] = '; '.join(f"{ins['length']} nt insertion at position {ins['position']}: {ins['merged']} reads (of {ins['coverage']}); {ins['merged']/(ins['coverage']-ins['merged']):.2f} AR" for ins in itds["FLT3"])
        result["UBTF_ITD"] = '; '.join(f"{ins['length']} nt insertion at position {ins['position']}: {ins['merged']} reads (of {ins['coverage']}); {ins['merged']/(ins['coverage']-ins['merged']):.2f} AR" for ins in itds["UBTF"])


        # ----------------------- SNVs / genotypes ----------------------------

        gens = json.loads(open(f"{run_dir}/genotypes.json").read())
        for g in ["PAX5", "IKZF1", "TPMT", "NUDT15"]:
            result[f"{g}_depth"] = f"{gens[g]['coverage']:.2f}"
            result[f"{g}_mutations"] = '; '.join(f"{m}: {gens[g]['mutations'][m]} AF" for m in gens[g]['mutations'])
            result[f"{g}_genotype"] = f"{gens[g]['genotype']}"

        result["reads_aligned"] = f"{int(float(k['reads_aligned']))}"
        result["on_target_read_length"] = f"{int(qf['qc']['nt_on_target'] / qf['qc']['reads_on_target'])}"
        result["on_target_depth"] = f"{qf['qc']['nt_on_target']/qf['qc']['target_regions_nt']:.2f}"

        if "meta" in qf and "run_start_time" in qf["meta"]:
            result["seq_date"] = qf["meta"]["run_start_time"]
            result["seq_ID"] = qf["meta"]["run_id"]
            result["basecalling_model"] = qf["meta"]["basecall_model"]
            result["library_ID"] = qf["meta"]["library_id"]
            result["sequencer_ID"] = qf["meta"]["sequencer_id"]
            result["flow_cell_ID"] = qf["meta"]["flow_cell_id"]

        aggregate_result.append(result)

    cols = list(set(k for r in aggregate_result for k in r))
    cols = ["name", "karyotype", "ISCN", "warnings", "fusions", "one_sided_fusions", "on_target_depth", "reads_aligned", "on_target_read_length",
 "PAX5_mutations", "PAX5_genotype", "PAX5_depth", "IKZF1_mutations", "IKZF1_genotype", "IKZF1_depth", "NUDT15_mutations", "NUDT15_genotype", "NUDT15_depth", "TPMT_mutations", "TPMT_genotype", "TPMT_depth", "FLT3_ITD", "UBTF_ITD", "KMT2A_ITD", "RUNX1", "CDKN2A", "CDKN2B", "ERG", "PAX5_deletions", "PAX5_ITD", "IKZF1_deletions", "seq_date", "library_ID", "basecalling_model", "sequencer_ID", "flow_cell_ID", "seq_ID"]
    print('\t'.join(cols))
    for r in aggregate_result:
        print('\t'.join(r[c] if c in r else "MISSING" for c in cols))



if __name__ == "__main__":
    parser = argparse.ArgumentParser("Make aggregate report (TSV) from analysis logs/directory")
    parser.add_argument("version", help="Analysis version (string)")
    parser.add_argument("dirs", help="Analysis output directories", nargs='+')
    args = parser.parse_args()
    main(args.dirs, args.version)
