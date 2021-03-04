def get_one_vcf(wildcards):
    if wildcards.dataset_name == "pf6_analysis_set":
        if wildcards.sample_name == "ref":
            return f'{config["input_data"]}/template.vcf.gz'
        else:
            return f'{config["varcall_dir"]}/cortex/pf6/{wildcards.sample_name}/cortex.vcf.gz'
    else:
        raise ValueError(f"Support for {wildcards.dataset_name} not implemented")


def get_sample_names(dataset_name):
    if dataset_name == "pf6_analysis_set":
        pf6_analysis_samples = load_pf6(config["pf6_tsv"], use_analysis_set=True)
        return ["ref"] + [rec.sample_name for rec in pf6_analysis_samples]
    else:
        raise ValueError(f"Support for {dataset_name} not implemented")


def load_bed(gene_list_name):
    bed_fname = f'{config["gene_bed_dir"]}/{gene_list_name}.bed'
    result = list()
    with open(bed_fname) as fin:
        for line in fin:
            name_field = line.split("\t")[3]
            result.append(name_field.split(";")[0])
    return result
