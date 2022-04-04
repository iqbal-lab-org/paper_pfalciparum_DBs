from functools import lru_cache

@lru_cache
def mp_get_sample_inducers(dataset_name, tool, mode):
    result = list()
    # "ref" must be prepended to sample names so that we place a ref genome slice at the beginning of the msas/prgs.
    for sample_name in ["ref"] + sample_names:
        if mode == "vcf":
            result.append(_mp_get_one_vcf(sample_name, dataset_name, tool))
        elif mode == "ref_genome":
            result.append(_mp_get_ref_genome(sample_name, dataset_name, tool))
        elif mode == "initial_vcf":
            result.append(_mp_get_initial_vcf(sample_name, dataset_name, tool))
        elif mode == "initial_ref_genome":
            result.append(cu_get_ref_genome_no_wildcards(dataset_name))

        elif mode == "sample_name":
            if sample_name == "ref":
                sample_name = "sample" # to induce from template.vcf.gz
            result.append(sample_name)
        else:
            raise ValueError(f"Unsupported mode: {mode}")
    return result

def _mp_get_one_vcf(sample_name, dataset_name, tool):
    if sample_name == "ref":
        return f'{config["input_data"]}/template.vcf.gz'
    if dataset_name.startswith("pf6"):
        ds_name = "pf6"
    elif dataset_name == "pvgv":
        ds_name = "pvgv"
    else:
        raise ValueError(f"Support for {dataset_name} not implemented")
    vcf_fname = tool
    if vcf_fname.startswith("gram_adju"):
        vcf_fname = "final"
    return f'{config["varcall_dir"]}/{tool}/{ds_name}/{sample_name}/{vcf_fname}.vcf.gz'

def _mp_get_initial_vcf(sample_name, dataset_name, tool):
    """
    Used to extract an initial vcf: e.g. gapfiller run on top of gram_adju's output
    """
    if tool == "gapfiller":
        tool = f"gram_adju_{GMTOOLS_COMMIT}"
    return _mp_get_one_vcf(sample_name, dataset_name, tool)

def _mp_get_ref_genome(sample_name, dataset_name, tool):
    if dataset_name.startswith("pf6") and tool == "gapfiller" and sample_name != "ref":
        return f'{config["varcall_dir"]}/{tool}/pf6/{sample_name}/induced_ref.fa.gz'
    else:
        return cu_get_ref_genome_no_wildcards(dataset_name)
