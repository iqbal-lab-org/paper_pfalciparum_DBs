from typing import List
from csv import DictReader
from pathlib import Path


def get_tool_vcf(wildcards):
    varcall_tool = map(
        lambda target: wildcards.tool.startswith(target),
        ["cortex", "minospb", "gram_adju"],
    )
    if wildcards.tool == "baseline":
        return f'{config["input_data"]}/template.vcf.gz'
    elif wildcards.tool == "myo_7_pf_genes":
        return f'{config["input_data"]}/barry_lab/{wildcards.tool}.vcf.gz'
    elif any(varcall_tool):
        tool_path = f'{config["varcall_dir"]}/{wildcards.tool}/{wildcards.dataset_name}/{wildcards.sample_name}/'
        if wildcards.tool == "cortex":
            tool_path += "cortex.vcf.gz"
        else:
            tool_path += "final.vcf.gz"
        return tool_path
    elif wildcards.tool.startswith("gram_jointgeno"):
        return f'{config["jointgeno_dir"]}/{wildcards.tool}/pf6_analysis_set/{wildcards.gene_list_name}_7_13/{wildcards.sample_name}/final.vcf.gz'
    else:
        raise ValueError(f"Unsupported toolname {wildcards.tool}")


def load_pf6_validation(tsv_fname: str) -> List[str]:
    with open(tsv_fname) as tsvfile:
        reader = DictReader(tsvfile, delimiter="\t")
        return list(map(lambda row: row["Sample"], reader))
