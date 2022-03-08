import re

_var_caller_regexp = re.compile("^gram_adju|^octopus|^cortex")
def gs_get_ref_genome(wildcards):
    if wildcards.tool.startswith("gram_jointgeno"):
        return cu_get_ref_genome(wildcards)
    elif wildcards.tool.startswith("gapfiller"):
        return f"{output_gram_adju}/{wildcards.dataset_name}/{wildcards.sample_name}/induced_ref.fa.gz", 
    elif _var_caller_regexp.match(wildcards.tool):
        return f"{output_gram_joint_geno}/{wildcards.dataset_name}/{wildcards.sample_name}/induced_ref.fa.gz", 
    else:
        raise ValueError(f"Unsupported tool: {wildcards.tool}")


def gs_get_translation_bed(wildcards):
    if wildcards.tool.startswith("gram_jointgeno"):
        return GENE_BED
    elif wildcards.tool.startswith("gapfiller"):
        return f"{output_gram_adju}/{wildcards.dataset_name}/{wildcards.sample_name}/{GENE_LIST_NAME}_induced.bed",
    elif _var_caller_regexp.match(wildcards.tool):
        return f"{output_gram_joint_geno}/{wildcards.dataset_name}/{wildcards.sample_name}/{GENE_LIST_NAME}_induced.bed",
    else:
        raise ValueError(f"Unsupported tool: {wildcards.tool}")

def gs_get_translation_bed_for_eval(wildcards):
    if wildcards.tool.startswith("gram_jointgeno"):
        return GENE_BED_FOR_EVAL
    elif wildcards.tool.startswith("gapfiller"):
        return f"{output_gram_adju}/{wildcards.dataset_name}/{wildcards.sample_name}/{GENE_LIST_NAME}_for_eval_induced.bed",
    elif _var_caller_regexp.match(wildcards.tool):
        return f"{output_gram_joint_geno}/{wildcards.dataset_name}/{wildcards.sample_name}/{GENE_LIST_NAME}_for_eval_induced.bed",
    else:
        raise ValueError(f"Unsupported tool: {wildcards.tool}")

def gs_get_expected_ir_stats(wildcards):
    if wildcards.dataset_name == "clone_trees":
        return expand(
            f"{output_ir_stats_per_sample}/{wildcards.dataset_name}/{{tool}}_{{sample_name}}.tsv",
            tool = ALL_TOOLS,
            sample_name = clone_tree_samples
            )
    elif wildcards.dataset_name == "crosses":
        return expand(
            f"{output_ir_stats_per_sample}/{wildcards.dataset_name}/{{tool}}_{{sample_name}}.tsv",
            tool = ALL_TOOLS,
            sample_name = crosses_samples
            )
    else:
        raise ValueError(f"Unsupported dataset")
