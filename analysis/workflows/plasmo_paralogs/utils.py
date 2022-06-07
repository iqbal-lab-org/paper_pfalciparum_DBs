from typing import List
from itertools import chain


def pp_get_all_paralog_regions(str_list_to_format):
    result = list()
    for str_to_format in str_list_to_format:
        for paralog_abbrev in config["paralog_abbrevs"]:
            result.append(str_to_format.format(paralog_abbrev =paralog_abbrev))
    return result


def pp_filter_to_paralogs(genes: List[str]) -> List[str]:
    genes_of_interest = set(chain.from_iterable(map(lambda gene: config["paralog_names"][gene], config["paralog_abbrevs"])))
    return list(filter(lambda gene: gene in genes_of_interest, genes))


def pp_get_all_vcfs(wildcards, mode):
    result = list()
    sample_names = _pp_load_sample_names(wildcards)
    for sample_name in sample_names:
        if mode == "sample_name":
            result.append(sample_name)
        elif mode == "vcf":
            result.append(_pp_get_one_vcf(wildcards, sample_name))
        else:
            raise ValueError(f"Unsupported mode: {mode}")
    return result

def _pp_load_sample_names(wildcards):
    if wildcards.sample_set_name == "pf6_analysis_set_fws95":
        result = ["ref"] + cu_get_sample_names_fws_matching(wildcards.sample_set_name)
    elif wildcards.sample_set_name == "pacb_ilmn_pf@pf6_analysis_set_fws95":
        result = cu_record_to_sample_names(cu_load_pacb_ilmn_pf(config["pacb_ilmn_pf_tsv"]))
    else:
        raise ValueError(f"Unsupported sample set name: {wildcards.sample_set_name}")
    return result


def _pp_get_one_vcf(wildcards, sample_name):
    GENE_LIST_NAME = "pf6_26_genes"
    if sample_name == "ref":
        return f'{config["input_data"]}/template.vcf.gz'
    if wildcards.tool == "malariaGEN":
        tool_name = "pf7" if "pacb" in wildcards.sample_set_name else "pf6"
        return f'{config["dl_output_dir"]}/vcfs/{tool_name}/combined_{GENE_LIST_NAME}_filterPASS.vcf.gz'
    elif wildcards.tool.startswith("gram_jointgeno"):
        elems = wildcards.sample_set_name.split("pf6_")
        return f'{config["jointgeno_dir"]}/{wildcards.tool}/{elems[0]}pf6/{elems[1]}/gapfiller/{GENE_LIST_NAME}_7_13/{sample_name}/final.vcf.gz'
    raise ValueError(f"Unsupported set of wildcards: {wildcards}")
