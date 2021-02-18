from typing import List
from csv import DictReader


class ENARecord:
    DATASET_NAMES = ["pf6", "pvgv", "pf_pacb_ilmn"]

    def __init__(self, dataset_name: str, sample_name: str, ena_IDs: str):
        if dataset_name not in self.DATASET_NAMES:
            raise ValueError(f"{dataset_name} not in {self.DATASET_NAMES}")
        self.dataset_name = dataset_name
        self.sample_name = sample_name
        self.ena_IDs = ena_IDs


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
