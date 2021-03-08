from pathlib import Path
from typing import List
from csv import DictReader
from subprocess import run as sp_run, PIPE
from glob import glob
from itertools import product


def get_gmtools_commit(container: str):
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


def mk_output_dirs(variables):
    """For each variable starting with 'output', makes the directory name it holds"""
    for variable in filter(lambda name: name.startswith("output"), variables):
        Path(eval(variable)).mkdir(exist_ok=True, parents=True)


def gene_list_to_genome(wildcards):
    if "pf6" in wildcards.gene_list_name:
        return "Pfalciparum"
    elif "pvivax" in wildcards.gene_list_name:
        return "PvivaxP01"
    else:
        raise ValueError(f"wildcards.{gene_list_name} not in {{pf6, pvgv}}")


def get_assembly(wildcards):
    assembly = glob(f'{config["assemblies_dir"]}/{wildcards.sample_name}*.fasta.gz')
    if len(assembly) != 1:
        raise ValueError(
            f'Error: expected one assembly for sample {wildcards.sample} in {config["assemblies_dir"]}, found {assembly}'
        )
    return assembly


def get_reads(wildcards):
    templated_reads = [
        f'{config["dl_output_dir"]}/{{dataset_name}}/{wildcards.sample_name}/reads_{i}.final.fastq.gz'
        for i in [1, 2]
    ]
    for recognised_name in ["pf6", "pvgv", "pacb_ilmn_pf"]:
        if wildcards.dataset_name.startswith(recognised_name):
            # This allows joint genotyping to use pacb_ilmn_pf samples on a specifiable pf6 graph
            replacement = wildcards.dataset_name.split("@")[0]
            return [
                read_file.format(dataset_name=replacement)
                for read_file in templated_reads
            ]
    raise ValueError(f"Support for {wildcards.dataset_name} not implemented")


def get_ref_genome(wildcards):
    ds_to_ref = {
        "pf6": "Pfalciparum",
        "pf6_analysis_set": "Pfalciparum",
        "pvgv": "PvivaxP01",
        "pacb_ilmn_pf": "Pfalciparum",
    }
    for key, val in ds_to_ref.items():
        if wildcards.dataset_name.startswith(key):
            return f'{config["dl_output_dir"]}/ref_genomes/{val}.genome.fasta.gz'
    raise ValueError(f"Support for {wildcards.dataset_name} not implemented")


def get_sample_names(dataset_name):
    if dataset_name == "pf6_analysis_set":
        loaded_samples = load_pf6(config["pf6_tsv"], use_analysis_set=True)
    elif dataset_name.startswith("pacb_ilmn_pf"):
        loaded_samples = load_pacb_ilmn_pf(config["pacb_ilmn_pf_tsv"])
    elif dataset_name == "pvgv":
        loaded_samples = load_pvgv(config["pvgv_tsv"])
    else:
        loaded_samples = None

    if loaded_samples is None:
        raise ValueError(f"Support for {dataset_name} not implemented")
    else:
        return [rec.sample_name for rec in loaded_samples]


def load_bed(gene_list_name):
    bed_fname = f'{config["gene_bed_dir"]}/{gene_list_name}.bed'
    result = list()
    with open(bed_fname) as fin:
        for line in fin:
            name = line.split("\t")[3]
            result.append(name)
    return result


##################
# Load sample tsvs#
##################
class ENARecord:
    DATASETS = {"pf6": "23.33Mb", "pvgv": "29.05Mb", "pacb_ilmn_pf": "23.33Mb"}

    def __init__(self, dataset_name: str, sample_name: str, ena_IDs: str):
        if dataset_name not in self.DATASETS:
            raise ValueError(f"{dataset_name} not in {self.DATASETS}")
        self.dataset_name = dataset_name
        self.sample_name = sample_name
        self.ena_IDs = ena_IDs


def get_genome_size(dataset_name: str):
    return ENARecord.DATASETS[dataset_name]


def load_pf6(
    tsv_fname: str, ignored_pattern: str = "", use_analysis_set: bool = False
) -> List[ENARecord]:
    ignored_samples = set()
    ignored_file = Path(tsv_fname).parent / "ignored_samples.tsv"
    with open(ignored_file) as tsvfile:
        reader = DictReader(tsvfile, delimiter="\t")
        for row in reader:
            if ignored_pattern == "" or row["Reason"].startswith(ignored_pattern):
                ignored_samples.add(row["Sample"])
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


def load_pvgv(tsv_fname: str) -> List[ENARecord]:
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


def load_pacb_ilmn_pf(tsv_fname: str) -> List[ENARecord]:
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


def find_reads(wildcards) -> List[str]:
    reads_dir = f'{config["ilmn_reads_dir"]}/{wildcards.sample}'
    reads_files = glob(f"{reads_dir}/**/**/*.fastq.gz")
    reads_files += glob(f"{reads_dir}/**/*.fastq.gz")
    reads_files += glob(f"{reads_dir}/*.fastq.gz")
    if len(reads_files) == 0:
        raise FileNotFoundError(f"No reads files found in {reads_dir}")
    for read_file in reads_files:
        if " " in read_file:
            raise ValueError(
                f"file {read_file} has whitespace in it, this breaks the pipeline. rename the file or change the separator in the pipeline at {sys.argv[0]}"
            )
    return reads_files
