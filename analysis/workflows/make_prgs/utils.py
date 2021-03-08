def get_one_vcf(wildcards):
    if wildcards.sample_name == "ref":
        return f'{config["input_data"]}/template.vcf.gz'
    if wildcards.dataset_name == "pf6_analysis_set":
        return (
            f'{config["varcall_dir"]}/cortex/pf6/{wildcards.sample_name}/cortex.vcf.gz'
        )
    elif wildcards.dataset_name == "pvgv":
        return (
            f'{config["varcall_dir"]}/cortex/pvgv/{wildcards.sample_name}/cortex.vcf.gz'
        )
    else:
        raise ValueError(f"Support for {wildcards.dataset_name} not implemented")
