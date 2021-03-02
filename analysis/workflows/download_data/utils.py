from typing import List
from csv import DictReader
from pathlib import Path


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


def load_pf6(tsv_fname: str, ignored_pattern: str = "") -> List[ENARecord]:
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
            # This skips the Indonesian samples, which were blocked from public release so have no ENA IDs
            if row["ENA"] == "" or row["Sample"] in ignored_samples:
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
