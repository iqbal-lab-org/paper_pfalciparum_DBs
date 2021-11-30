GMTOOLS_COMMIT = cu_get_gmtools_commit(config["container_gramtools"])
gram_adju = f"gram_adju_{GMTOOLS_COMMIT}"
gram_jointgeno = f"gram_jointgeno_{GMTOOLS_COMMIT}"

def qc_get_tool_vcf(wildcards):
    varcall_tool = map(
        lambda target: wildcards.tool.startswith(target),
        ["cortex", "octopus", "gram_adju"],
    )
    if any(varcall_tool):
        dataset_name = "pf6" if wildcards.dataset_name.startswith("pf6") else wildcards.dataset_name
        tool_path = f'{config["varcall_dir"]}/{wildcards.tool}/{dataset_name}/{wildcards.sample_name}/'
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

def qc_get_expected_vcf_stats(wildcards):
    """
    Cf `qc_get_one_vcf` for how tool names in this function map to vcf output file names.
    """
    if wildcards.dataset_name.startswith("pf6"):
        tools = [gram_adju, f"{gram_jointgeno}__{wildcards.dataset_name}__7__13"]
        samples = cu_record_to_sample_names(cu_load_pf6(config["pf6_validation_tsv"]))
        #samples = cu_get_sample_names_fws_matching(wildcards.dataset_name)
    elif wildcards.dataset_name.startswith("pacb_ilmn_pf"):
        tools = [gram_adju, f"{gram_jointgeno}__{wildcards.dataset_name}@pf6_analysis_set_fws95__7__13"]
        samples = cu_record_to_sample_names(cu_load_pacb_ilmn_pf(config["pacb_ilmn_pf_tsv"]))
    else:
        raise ValueError(f"Unsupported dataset name: {wildcards.dataset_name}")
    return expand(
        f"{output_vcf_stats}/{wildcards.dataset_name}/{wildcards.gene_list}/{{tool}}/{{sample_name}}_vcf_stats.tsv",
        tool=tools,
        sample_name=samples,
    )


def qc_get_mapped_reads(wildcards):
    dataset_name = "pf6" if wildcards.dataset_name.startswith("pf6") else wildcards.dataset_name
    return f'{config["varcall_dir"]}/mapped_reads/{dataset_name}/{wildcards.sample_name}/mapped.bam'

def qc_get_expected_read_stats(wildcards):
    dataset_name = wildcards.dataset_name
    if wildcards.dataset_name.startswith("pf6"):
        samples = cu_record_to_sample_names(cu_load_pf6(config["pf6_validation_tsv"]))
        dataset_name = "pf6"
        #samples = cu_get_sample_names_fws_matching(wildcards.dataset_name)
    elif wildcards.dataset_name.startswith("pacb_ilmn_pf"):
        samples = cu_record_to_sample_names(cu_load_pacb_ilmn_pf(config["pacb_ilmn_pf_tsv"]))
    else:
        raise ValueError(f"Unsupported dataset name: {wildcards.dataset_name}")
    return expand(
        f"{output_read_stats}/per_sample/{dataset_name}/{{sample_name}}_read_stats.tsv",
        sample_name=samples,
    )
