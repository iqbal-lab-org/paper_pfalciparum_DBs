from pathlib import Path
from typing import List, Optional
from csv import DictReader
from subprocess import run as sp_run, PIPE
from glob import glob
from itertools import product
import re


def cu_get_gmtools_commit(container: str):
    """Get gramtools commit version through venv/singularity container it is installed in"""
    gmtools_commit_script = Path(config["scripts"]) / "gmtools_commit.py"
    gmtools_commit_script = str(gmtools_commit_script.resolve())
    try:  # Try in virtual environment first
        GMTOOLS_COMMIT = sp_run(
            ["python3", gmtools_commit_script],
            capture_output=True,
            text=True,
            check=True,
        )
    except:
        GMTOOLS_COMMIT = sp_run(
            [
                "singularity",
                "exec",
                str(Path(container).resolve()),
                "python3",
                gmtools_commit_script,
            ],
            capture_output=True,
            text=True,
            check=True,
        )
    GMTOOLS_COMMIT = GMTOOLS_COMMIT.stdout.strip()
    return GMTOOLS_COMMIT


def cu_mk_output_dirs(variables):
    """For each variable starting with 'output', makes the directory name it holds"""
    for variable in filter(lambda name: name.startswith("output"), variables):
        Path(eval(variable)).mkdir(exist_ok=True, parents=True)


def cu_gene_list_to_genome(wildcards):
    if "pf6" in wildcards.gene_list_name:
        return "Pfalciparum"
    elif "pvivax" in wildcards.gene_list_name:
        return "PvivaxP01"
    else:
        raise ValueError(f"wildcards.{gene_list_name} not in {{pf6, pvgv}}")


def cu_get_assembly(wildcards):
    assembly = glob(f'{config["assemblies_dir"]}/{wildcards.sample_name}*.fasta.gz')
    if len(assembly) != 1:
        raise ValueError(
            f'Error: expected one assembly for sample {wildcards.sample} in {config["assemblies_dir"]}, found {assembly}'
        )
    return assembly


def cu_get_reads(wildcards):
    templated_reads = [
        f'{config["dl_output_dir"]}/{{dataset_name}}/{wildcards.sample_name}/reads_{i}.final.fastq.gz'
        for i in [1, 2]
    ]
    for recognised_name in ["pf6", "pvgv", "pacb_ilmn_pf","clone_trees","crosses"]:
        if wildcards.dataset_name.startswith(recognised_name):
            return [
                read_file.format(dataset_name=recognised_name)
                for read_file in templated_reads
            ]
    raise ValueError(f"Support for {wildcards.dataset_name} not implemented")


def cu_get_ref_genome_no_wildcards(dataset_name):
    ds_to_ref = {
        "pf6": "Pfalciparum",
        "pf6_analysis_set": "Pfalciparum",
        "pvgv": "PvivaxP01",
        "pacb_ilmn_pf": "Pfalciparum",
        "clone_trees": "Pfalciparum",
        "crosses": "Pfalciparum",
    }
    for key, val in ds_to_ref.items():
        if dataset_name.startswith(key):
            return f'{config["dl_output_dir"]}/ref_genomes/{val}.genome.fasta.gz'
    raise ValueError(f"Support for {dataset_name} not implemented")


def cu_get_ref_genome(wildcards):
    return cu_get_ref_genome_no_wildcards(wildcards.dataset_name)

_PF_GENOME = "23.33Mb"
_PV_GENOME = "29.05Mb"

##################
# Load sample tsvs#
##################
class ENARecord:
    DATASETS = {"pf6": _PF_GENOME, "pvgv": _PV_GENOME, "pacb_ilmn_pf": _PF_GENOME, "clone_trees": _PF_GENOME, "crosses": _PF_GENOME}
    def __init__(self, dataset_name: str, sample_name: str, ena_IDs: str):
        if dataset_name not in self.DATASETS:
            raise ValueError(f"{dataset_name} not in {self.DATASETS}")
        self.dataset_name = dataset_name
        self.sample_name = sample_name
        self.ena_IDs = ena_IDs


def cu_get_genome_size(dataset_name: str):
    return ENARecord.DATASETS[dataset_name]


def _get_samples_below_fws(tsv_fname: str, fws_threshold: int):
    if fws_threshold <= 0 or fws_threshold >= 100:
        raise ValueError("User-provided Fws threshold must be 0<x<100")
    used_fws_threshold = fws_threshold / 100
    fws_fname = Path(tsv_fname).parent / "Pf_6_fws.tsv"
    result = set()
    with open(fws_fname) as tsvfile:
        reader = DictReader(tsvfile, delimiter="\t")
        for row in reader:
            if float(row["Fws"]) < used_fws_threshold:
                result.add(row["Sample"])
    return result


def cu_load_pf6(
    tsv_fname: str,
    use_analysis_set: bool = False,
    fws_threshold: Optional[float] = None,
) -> List[ENARecord]:
    ignored_samples = set()
    ignored_file = Path(tsv_fname).parent / "ignored_samples.tsv"
    with open(ignored_file) as tsvfile:
        reader = DictReader(tsvfile, delimiter="\t")
        for row in reader:
            ignored_samples.add(row["Sample"])
    if fws_threshold is not None:
        ignored_samples.update(_get_samples_below_fws(tsv_fname, fws_threshold))
    result = list()
    with open(tsv_fname) as tsvfile:
        reader = DictReader(tsvfile, delimiter="\t")
        for row in reader:
            # Skips the Indonesian samples, which were blocked from public release so have no ENA IDs,
            # and any other ignored samples (see the tsv for reason)
            if row["ENA"] == "" or row["Sample"] in ignored_samples:
                continue
            # Skip non analysis_set pf6 samples if asked for
            if use_analysis_set and row["Exclusion reason"] != "Analysis_set":
                continue
            result.append(ENARecord("pf6", row["Sample"], row["ENA"]))
    return result


def cu_load_pvgv(tsv_fname: str) -> List[ENARecord]:
    result = list()
    with open(tsv_fname) as tsvfile:
        reader = DictReader(tsvfile, delimiter="\t")
        for row in reader:
            ena_IDs = row["Sequencing Reads"].split(";")
            # This skips the Indonesian samples, which were blocked from public release so have no ENA IDs
            if all(map(lambda ID: ID == "", ena_IDs)):
                continue
            ena_IDs = ",".join(map(lambda s: s.strip(), ena_IDs))
            result.append(ENARecord("pvgv", row["Sample ID"], ena_IDs))
    return result


def cu_load_pacb_ilmn_pf(tsv_fname: str) -> List[ENARecord]:
    result = list()
    with open(tsv_fname) as tsvfile:
        reader = DictReader(tsvfile, delimiter="\t")
        for row in reader:
            # Filter out PfML01, as it requires too much memory for cortex to process (with mem_height=22)
            if row["sample_name"].startswith("#") or row["sample_name"] in {
                "Pf3D7",
                "PfML01",
            }:
                continue
            ena_IDs = row["ENA_sample_accession_ILMN"]
            result.append(ENARecord("pacb_ilmn_pf", row["sample_name"], ena_IDs))
    return result


def cu_load_clone_trees(tsv_fname: str) -> List[ENARecord]:
    result = list()
    with open(tsv_fname) as tsvfile:
        reader = DictReader(tsvfile, delimiter="\t")
        for row in reader:
            result.append(
                ENARecord("clone_trees", row["BAM_ID"], row["ENA_sample_accession"])
            )
    return result

def cu_load_crosses(tsv_fname: str) -> List[ENARecord]:
    result = list()
    with open(tsv_fname) as tsvfile:
        reader = DictReader(tsvfile, delimiter="\t")
        for row in reader:
            # This sample seg faults for octopus
            if row["Sample"] == "PG0017-C":
                continue
            result.append(
                ENARecord("crosses", row["Sample"], row["Accession"])
            )
    return result


def cu_record_to_sample_names(records: List[ENARecord]):
    return [rec.sample_name for rec in records]


def cu_get_sample_names_fws_matching(dataset_name):
    fws_pattern = "pf6.*fws([0-9]{1,2})"
    fws_threshold = re.match(fws_pattern, dataset_name)
    if fws_threshold is not None:
        fws_threshold = int(fws_threshold.groups()[0])
    sample_names = cu_get_sample_names(dataset_name, fws_threshold=fws_threshold)
    return sample_names


def cu_get_sample_names(dataset_name, fws_threshold=None):
    if dataset_name.startswith("pf6"):
        if fws_threshold is not None:
            supported_pattern = f"pf6.*fws{fws_threshold}"
            assert (
                re.match(supported_pattern, dataset_name) is not None
            ), f"Fws filtering only supported for dataset name matching pattern {supported_pattern}"
        use_analysis_set = "analysis_set" in dataset_name
        sample_tsv = config["pf6_tsv"]
        loaded_samples = cu_load_pf6(
            sample_tsv, use_analysis_set=use_analysis_set, fws_threshold=fws_threshold
        )
    elif dataset_name.startswith("pacb_ilmn_pf"):
        loaded_samples = cu_load_pacb_ilmn_pf(config["pacb_ilmn_pf_tsv"])
    elif dataset_name == "pvgv":
        loaded_samples = cu_load_pvgv(config["pvgv_tsv"])
    else:
        loaded_samples = None

    if loaded_samples is None:
        raise ValueError(f"Support for {dataset_name} not implemented")
    else:
        return cu_record_to_sample_names(loaded_samples)


def cu_load_bed(gene_list_name):
    bed_fname = f'{config["gene_bed_dir"]}/{gene_list_name}.bed'
    result = list()
    with open(bed_fname) as fin:
        for line in fin:
            name = line.split("\t")[3]
            result.append(name)
    return result
