from typing import List
from csv import DictReader
from pathlib import Path

GMTOOLS_COMMIT = cu_get_gmtools_commit(config["container_gramtools"])
gram_adju = f"gram_adju_{GMTOOLS_COMMIT}"
gram_jointgeno = f"gram_jointgeno_{GMTOOLS_COMMIT}"

def ev_get_tool_vcf(wildcards):
    varcall_tool = map(
        lambda target: wildcards.tool.startswith(target),
        ["cortex", "octopus", "gram_adju"],
    )
    malariagen_vcf = wildcards.tool in {"pf6","pf7"}
    if wildcards.tool == "baseline":
        return f'{config["input_data"]}/template.vcf.gz'
    elif wildcards.tool == "myo_7_pf_genes":
        return f'{config["input_data"]}/barry_lab/{wildcards.tool}.vcf.gz'
    elif wildcards.tool == "paolo_pvgv":
        return f'{config["input_data"]}/barry_lab/{wildcards.tool}.vcf.gz'
    elif malariagen_vcf:
        return f'{config["dl_output_dir"]}/vcfs/{wildcards.tool}/combined_{wildcards.gene_list_name}_filterPASS.vcf.gz'
    elif any(varcall_tool):
        tool_path = f'{config["varcall_dir"]}/{wildcards.tool}/{wildcards.dataset_name}/{wildcards.sample_name}/'
        if wildcards.tool.startswith("gram"):
            tool_path += "final.vcf.gz"
        else:
            tool_path += f"{wildcards.tool}.vcf.gz"
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
    pacb_ilmn_snames = cu_record_to_sample_names(cu_load_pacb_ilmn_pf(config["pacb_ilmn_pf_tsv"]))
    if wildcards.gene_list_name.startswith("pf6"):
        # baseline: runs mapping on an empty vcf, giving a baseline to compare tools to
        tools = [
            "baseline",
            "octopus",
            gram_adju,
            "cortex",
            "pf7",
            f"{gram_jointgeno}__pacb_ilmn_pf@pf6_analysis_set__7__13",
            f"{gram_jointgeno}__pacb_ilmn_pf@pf6_analysis_set_fws95__7__13",
        ]
    else:
        raise ValueError(f"Unsupported dataset name: {wildcards.dataset_name}")
    return expand(
        f"{output_bowtie2_pf}/{wildcards.gene_list_name}/{{tool}}/{{sample_name}}.sam",
        sample_name=pacb_ilmn_snames,
        tool=tools,
    )

def ev_get_expected_varifier_outputs(gene_list_name):
    pacb_ilmn_snames = cu_record_to_sample_names(cu_load_pacb_ilmn_pf(config["pacb_ilmn_pf_tsv"]))
    if gene_list_name.startswith("pf6"):
        tools = [
            "octopus",
            gram_adju,
            "cortex",
            "pf7",
            f"{gram_jointgeno}__pacb_ilmn_pf@pf6_analysis_set_fws95__7__13",
        ]
    else:
        raise ValueError(f"Unsupported dataset name: {dataset_name}")
    return expand(
        f"{output_varifier_pf}/{gene_list_name}/{{tool}}/{{sample_name}}/summary_stats.json",
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
            "octopus",
            "pf6",
            gram_adju,
            f"{gram_jointgeno}__pf6_analysis_set_fws95__7__13",
            #f"{gram_jointgeno}__pf6_analysis_set__12__13",
        ]
        samples = cu_record_to_sample_names(cu_load_pf6(config["pf6_validation_tsv"]))
    elif wildcards.dataset_name.startswith("pvgv"):
        tools = [f"cortex", "paolo_pvgv", f"{gram_jointgeno}__pvgv__7__13"]
        samples = cu_record_to_sample_names(cu_load_pvgv(config["pvgv_validation_tsv"]))
    elif wildcards.dataset_name.startswith("pacb_ilmn_pf"):
        tools = [
            "baseline",
            "cortex",
            "octopus",
            "pf7",
            gram_adju,
            f"{gram_jointgeno}__pacb_ilmn_pf@pf6_analysis_set_fws95__7__13",
            ]
        samples = cu_record_to_sample_names(cu_load_pacb_ilmn_pf(config["pacb_ilmn_pf_tsv"]))
    else:
        raise ValueError(f"Unsupported dataset name: {wildcards.dataset_name}")
    return expand(
        f"{output_ir}/{wildcards.dataset_name}/{wildcards.gene_list}/{{tool}}/{{sample_name}}_ir_stats.tsv",
        tool=tools,
        sample_name=samples,
    )
