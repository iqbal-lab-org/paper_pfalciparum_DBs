from typing import List
from csv import DictReader


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


def load_pf6(tsv_fname: str) -> List[ENARecord]:
    result = list()
    with open(tsv_fname) as tsvfile:
        reader = DictReader(tsvfile, delimiter="\t")
        for row in reader:
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
            if row["sample_name"].startswith("#") or '3D7' in row["sample_name"]:
                continue
            ena_IDs = row["ENA_sample_accession_ILMN"]
            result.append(ENARecord("pacb_ilmn_pf", row["sample_name"], ena_IDs))
    return result
