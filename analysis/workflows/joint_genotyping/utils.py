def jg_get_gram_build(wildcards):
    """
    Converts a dataset name into a gramtools build output path, supporting pacb_ilmn_pf matching to pf6-built prgs
    """
    for recognised_name in ["pf6", "pvgv", "pacb_ilmn_pf"]:
        if wildcards.dataset_name.startswith(recognised_name):
            elems = wildcards.dataset_name.split("@")
            dataset_name = elems[-1]
            return f"{output_gram_build}/{dataset_name}/{wildcards.gene_list_name}_mml{wildcards.min_match_len}_k{wildcards.kmer_size}/kmers_stats"
    raise ValueError(f"Support for {wildcards.dataset_name} not implemented")

