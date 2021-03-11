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


def filter_to_validation(sample_names, dataset_name):
    validation_samples = set()
    if dataset_name.startswith("pf6"):
        validation_samples = set(
            map(
                lambda record: record.sample_name,
                load_pf6(config["pf6_validation_tsv"]),
            )
        )
    elif dataset_name.startswith("pvgv"):
        validation_samples = set(
            map(
                lambda record: record.sample_name,
                load_pvgv(config["pvgv_validation_tsv"]),
            )
        )
    elif dataset_name.startswith("pacb_ilmn_pf"):
        return sample_names
    else:
        raise ValueError(f"Support for {dataset_name} not implemented")
    return [sname for sname in sample_names if sname in validation_samples]
