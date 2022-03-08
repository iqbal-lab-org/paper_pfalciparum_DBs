def ct_get_reads(wildcards):
    wildcards.dataset_name = "clone_trees"
    return cu_get_reads(wildcards)

def ct_get_ref_genome(wildcards):
    if wildcards.tool.startswith("gram_jointgeno"):
        wildcards.dataset_name = "clone_trees"
        return cu_get_ref_genome(wildcards)
    elif wildcards.tool.startswith("gram_adju"):
        return f"{output_gram_joint_geno}/{wildcards.sample_name}/induced_ref.fa.gz", 
    else:
        raise ValueError(f"Unsupported tool: {wildcards.tool}")
