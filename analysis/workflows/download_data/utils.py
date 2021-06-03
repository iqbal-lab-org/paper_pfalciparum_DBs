def get_vcf_url(wildcards):
    if wildcards.dataset_name == "pf6":
        url_template = config["vcf"]["pf6_url"] + config["vcf"]["pf6_fname_template"]
        return url_template.format(chrom=wildcards.chrom)
    raise ValueError(f"Unsupported dataset name: {wildcards.dataset_name}")
