def get_one_vcf(wildcards):
    if wildcards.sample_name == "ref":
        return f'{config["input_data"]}/template.vcf.gz'
    if dataset_name.startswith("pf6"):
        return (
            f'{config["varcall_dir"]}/cortex/pf6/{wildcards.sample_name}/cortex.vcf.gz'
        )
    elif dataset_name == "pvgv":
        return f'{config["varcall_dir"]}/cortex/pvgv/{sample_name}/cortex.vcf.gz'
    else:
        raise ValueError(f"Support for {dataset_name} not implemented")
