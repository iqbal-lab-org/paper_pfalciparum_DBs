import re

def ct_get_reads(wildcards):
    wildcards.dataset_name = "clone_trees"
    return cu_get_reads(wildcards)

_var_caller_regexp = re.compile("^gram_adju|^octopus|^cortex")
def ct_get_ref_genome(wildcards):
    if wildcards.tool.startswith("gram_jointgeno"):
        wildcards.dataset_name = "clone_trees"
        return cu_get_ref_genome(wildcards)
    elif wildcards.tool.startswith("gapfiller"):
        return f"{output_gram_adju}/{wildcards.sample_name}/induced_ref.fa.gz", 
    elif _var_caller_regexp.match(wildcards.tool):
        return f"{output_gram_joint_geno}/{wildcards.sample_name}/induced_ref.fa.gz", 
    else:
        raise ValueError(f"Unsupported tool: {wildcards.tool}")


def ct_get_translation_bed(wildcards):
    if wildcards.tool.startswith("gram_jointgeno"):
        return GENE_BED
    elif wildcards.tool.startswith("gapfiller"):
        return f"{output_gram_adju}/{wildcards.sample_name}/{GENE_LIST_NAME}_induced.bed",
    elif _var_caller_regexp.match(wildcards.tool):
        return f"{output_gram_joint_geno}/{wildcards.sample_name}/{GENE_LIST_NAME}_induced.bed",
    else:
        raise ValueError(f"Unsupported tool: {wildcards.tool}")

def ct_get_translation_bed_for_eval(wildcards):
    if wildcards.tool.startswith("gram_jointgeno"):
        return GENE_BED_FOR_EVAL
    elif wildcards.tool.startswith("gapfiller"):
        return f"{output_gram_adju}/{wildcards.sample_name}/{GENE_LIST_NAME}_for_eval_induced.bed",
    elif _var_caller_regexp.match(wildcards.tool):
        return f"{output_gram_joint_geno}/{wildcards.sample_name}/{GENE_LIST_NAME}_for_eval_induced.bed",
    else:
        raise ValueError(f"Unsupported tool: {wildcards.tool}")

def ct_get_expected_ir_stats(wildcards):
    return expand(
        f"{output_ir_stats_per_sample}/{{tool}}_{{sample_name}}.tsv",
        tool = [gram_joint_geno_path,"cortex","octopus",f"gram_adju_{GMTOOLS_COMMIT}","gapfiller"],
        sample_name = clone_tree_samples
        )
