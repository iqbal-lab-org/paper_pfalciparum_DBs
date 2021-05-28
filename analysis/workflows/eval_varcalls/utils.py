from typing import List
from csv import DictReader
from pathlib import Path

GMTOOLS_COMMIT = get_gmtools_commit(config["container_gramtools"])
gram_adju = f"gram_adju_{GMTOOLS_COMMIT}"
gram_jointgeno = f"gram_jointgeno_{GMTOOLS_COMMIT}"


def ev_get_tool_vcf(wildcards):
    varcall_tool = map(
        lambda target: wildcards.tool.startswith(target),
        ["cortex", "minospb", "gram_adju"],
    )
    if wildcards.tool == "baseline":
        return f'{config["input_data"]}/template.vcf.gz'
    elif wildcards.tool == "myo_7_pf_genes":
        return f'{config["input_data"]}/barry_lab/{wildcards.tool}.vcf.gz'
    elif wildcards.tool == "paolo_pvgv":
        return f'{config["input_data"]}/barry_lab/{wildcards.tool}.vcf.gz'
    elif any(varcall_tool):
        tool_path = f'{config["varcall_dir"]}/{wildcards.tool}/{wildcards.dataset_name}/{wildcards.sample_name}/'
        if wildcards.tool == "cortex":
            tool_path += "cortex.vcf.gz"
        else:
            tool_path += "final.vcf.gz"
        return tool_path
    elif wildcards.tool.startswith("gram_jointgeno"):
        elems = wildcards.tool.split("__")
        gram_version = elems[0]
        dataset_name = elems[1]
        prg_min_match_len = elems[2]
        prg_kmer_size = elems[3]
        return f'{config["jointgeno_dir"]}/{gram_version}/{dataset_name}/{wildcards.gene_list_name}_{prg_min_match_len}_{prg_kmer_size}/{wildcards.sample_name}/final.vcf.gz'
    else:
        raise ValueError(f"Unsupported toolname {wildcards.tool}")


def ev_get_expected_alignments(wildcards):
    """
    For mapping-to-assembly-based call validation
    """
    pacb_ilmn_records = load_pacb_ilmn_pf(config["pacb_ilmn_pf_tsv"])
    pacb_ilmn_snames = [rec.sample_name for rec in pacb_ilmn_records]
    if wildcards.gene_list_name.startswith("pf6"):
        # baseline: runs mapping on an empty vcf, giving a baseline to compare tools to
        tools = [
            "baseline",
            "minospb",
            gram_adju,
            "cortex",
            f"{gram_jointgeno}__pacb_ilmn_pf@pf6_analysis_set__7__13",
            f"{gram_jointgeno}__pacb_ilmn_pf@pf6_analysis_set_1500__7__13",
            f"{gram_jointgeno}__pacb_ilmn_pf@pf6_analysis_set_1500__12__13",
            f"{gram_jointgeno}__pacb_ilmn_pf@pf6_analysis_set_3000__7__13",
            f"{gram_jointgeno}__pacb_ilmn_pf@pf6_analysis_set_3000__12__13",
        ]
    else:
        raise ValueError(f"Unsupported dataset name: {wildcards.dataset_name}")
    return expand(
        f"{output_bowtie2_pf}/{wildcards.gene_list_name}/{{tool}}/{{sample_name}}.sam",
        sample_name=pacb_ilmn_snames,
        tool=tools,
    )


def ev_get_expected_stats(wildcards):
    """
    For induced_ref mapping based call validation
    """
    if wildcards.dataset_name.startswith("pf6"):
        tools = [
            "baseline",
            "cortex",
            "myo_7_pf_genes",
            f"{gram_jointgeno}__pf6_analysis_set__7__13",
            # f"{gram_jointgeno}__pf6_analysis_set_3000__7__13",
            # f"{gram_jointgeno}__pf6_analysis_set_3000__12__13",
            # f"{gram_jointgeno}__pf6_analysis_set_1500__7__13",
            # f"{gram_jointgeno}__pf6_analysis_set_1500__12__13",
        ]
        samples = record_to_sample_names(load_pf6(config["pf6_validation_tsv"]))
    elif wildcards.dataset_name.startswith("pvgv"):
        tools = [f"cortex", "paolo_pvgv", f"{gram_jointgeno}__pvgv__7__13"]
        samples = record_to_sample_names(load_pvgv(config["pvgv_validation_tsv"]))
    else:
        raise ValueError(f"Unsupported dataset name: {wildcards.dataset_name}")
    return expand(
        f"{output_ir}/{wildcards.dataset_name}/{wildcards.gene_list}/{{tool}}/{{sample_name}}_ir_stats.tsv",
        tool=tools,
        sample_name=samples,
    )
