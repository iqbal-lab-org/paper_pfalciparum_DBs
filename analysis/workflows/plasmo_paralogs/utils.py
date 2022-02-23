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


def pp_get_one_vcf(wildcards):
    if wildcards.sample_name == "ref":
        return f'{config["input_data"]}/template.vcf.gz'
    if wildcards.sample_set_name == "pf6_analysis_set_fws95":
        return f'{config["jointgeno_dir"]}/gram_jointgeno_{GMTOOLS_COMMIT}/pf6/analysis_set_fws95/{wildcards.tool}/pf6_26_genes_7_13/{wildcards.sample_name}/final.vcf.gz'
    if wildcards.sample_set_name == "pacb_ilmn_pf@pf6_analysis_set_fws95":
        return f'{config["jointgeno_dir"]}/gram_jointgeno_{GMTOOLS_COMMIT}/pacb_ilmn_pf@pf6/analysis_set_fws95/{wildcards.tool}/pf6_26_genes_7_13/{wildcards.sample_name}/final.vcf.gz'
    raise ValueError(f"Unsupported set of wildcards: {wildcards}")
