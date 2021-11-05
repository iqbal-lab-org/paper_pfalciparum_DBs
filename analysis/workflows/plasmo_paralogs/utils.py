from typing import List
from itertools import chain

def pp_get_all_paralog_regions(str_list_to_format):
    result = list()
    for str_to_format in str_list_to_format:
        for paralog in config["paralog_abbrevs"]:
            for paralog_region in config["focal_regions"][paralog]:
                result.append(str_to_format.format(paralog=paralog,paralog_region=paralog_region))
    return result

def pp_filter_to_paralogs(genes: List[str]) -> List[str]:
    genes_of_interest = set(chain.from_iterable(map(lambda gene: config["paralog_names"][gene], config["paralog_abbrevs"])))
    return list(filter(lambda gene: gene in genes_of_interest, genes))

def pp_get_one_vcf(wildcards):
    pass
