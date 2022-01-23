def mp_get_one_vcf(wildcards):
    if wildcards.sample_name == "ref":
        return f'{config["input_data"]}/template.vcf.gz'
    if dataset_name.startswith("pf6"):
        ds_name = "pf6"
    elif dataset_name == "pvgv":
        ds_name = "pvgv"
    else:
        raise ValueError(f"Support for {dataset_name} not implemented")
    vcf_fname = wildcards.tool
    if wildcards.tool.startswith("gram_adju"):
        vcf_fname = "final"
    return f'{config["varcall_dir"]}/{wildcards.tool}/{ds_name}/{wildcards.sample_name}/{vcf_fname}.vcf.gz'

def mp_get_initial_vcf(wildcards):
    """
    Used to extract an initial vcf: e.g. gapfiller run on top of gram_adju's output
    """
    if wildcards.tool == "gapfiller":
        wildcards.tool = f"gram_adju_{GMTOOLS_COMMIT}"
    return mp_get_one_vcf(wildcards)

def mp_get_ref_genome(wildcards):
    if wildcards.dataset_name.startswith("pf6") and wildcards.tool == "gapfiller" and wildcards.sample_name != "ref":
        return f'{config["varcall_dir"]}/{wildcards.tool}/pf6/{wildcards.sample_name}/induced_ref.fa.gz'
    else:
        return cu_get_ref_genome(wildcards)
