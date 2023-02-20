import re
from functools import lru_cache
from dataclasses import dataclass
from typing import List
from collections import defaultdict

from common_utils.metrics import MetricsRecorder
from common_utils.sqlite import drop_if_exists
from analysed_sequences.msa_coordinate_adjust import coordinate_adjust

@dataclass
class Donor:
    seq_ID: str
    start: int
    end: int

Target = str
DonorList = List[Donor]

class Donors:
    class RecombMetrics(MetricsRecorder):
        _headers = ["target", "donor", "start", "end", "end_is_recomb_breakpoint"]
        _required_headers = _headers

    def __init__(self):
        self._donor_map: Dict[Target, DonorList] = defaultdict(list)

    def __len__(self):
        return len(self._donors)

    def __iter__(self):
        return iter(self._donor_map.items())

    def __getitem__(self, key):
        return self._donor_map[key]

    def __contains__(self, key):
        return key in self._donor_map

    def keys(self):
        return self._donor_map.keys()

    def values(self):
        return self._donor_map.values()

    def items(self):
        return self._donor_map.items()

    def add_donors(self, target_ID: str, donors: List[Donor]):
        if target_ID in self._donor_map:
            raise ValueError(f"{target_ID} already populated")
        for donor in donors:
            self._donor_map[target_ID].append(donor)

    def position_match(self, target_ID, donor_ID, query_pos) -> bool:
        result = False
        matching_donors = [
            donor for donor in self._donor_map[target_ID] if donor.seq_ID == donor_ID
        ]
        for donor in matching_donors:
            if query_pos >= donor.start and query_pos <= donor.end:
                result = True
                break
        return result

    @lru_cache
    def id_match(self, target_ID, donor_ID) -> bool:
        result = False
        donors = self._donor_map.get(target_ID, None)
        if donors is not None:
            for donor in donors:
                if donor.seq_ID == donor_ID:
                    result = True
                    break
        return result

    def _to_tsv_metrics(self):
        result = list()
        for target_ID, donor_list in self:
            last_donor_idx = len(donor_list) - 1
            for i, donor in enumerate(donor_list):
                end_is_recomb_breakpoint = i != last_donor_idx
                result.append(
                    self.RecombMetrics(
                        target=target_ID,
                        donor=donor.seq_ID,
                        start=donor.start,
                        end=donor.end,
                        end_is_recomb_breakpoint=end_is_recomb_breakpoint,
                    )
                )
        return result

    def write_tsv(self, con, sqlite_db_fpath, table_name):
        db_cursor = con.cursor()
        table_name = f"mosaic_aligner_{table_name}"
        drop_if_exists(sqlite_db_fpath, db_cursor, table_name)
        db_cursor.execute(
            f"create table {table_name} (target,donor,start,end,end_is_recomb_breakpoint)",
        )
        for elem in self._to_tsv_metrics():
            db_cursor.execute(
                f"insert into {table_name} values (?,?,?,?,?)",
                (
                    elem["target"],
                    elem["donor"],
                    elem["start"],
                    elem["end"],
                    elem["end_is_recomb_breakpoint"],
                ),
            )
        con.commit()

def parse_mosaic_alignment_file(
    record_dict, mosaic_alignment_fname, get_names_only=False
) -> Donors:
    """
    `record_dict` must contain the {id: seq} of the MSA of all sequences present in
    `mosaic_alignment_fname`, so that the aligned sequences by mosaic aligner can be
    translated to MSA coordinate space.
    `get_names_only`: only the names of the targets and their donors are retrieved,
    not the position of breakpoints
    """
    target_regexp = re.compile(r"Target:\s*(?P<seq_ID>[\w-]+)", re.ASCII)
    donor_regexp = re.compile(
        r"\s*(?P<seq_ID>[\w-]+)\s*\[\s*(?P<start>[0-9]+)-\s*(?P<end>[0-9]+)\]", re.ASCII
    )
    result = Donors()
    donor_list = list()
    with open(mosaic_alignment_fname) as fin:
        alignment_text = fin.readlines()
    for line in alignment_text:
        target_match = target_regexp.match(line)
        if target_match is not None:
            if len(donor_list) > 0:
                result.add_donors(target_ID, (donor_list))
                donor_list = list()
            target_ID = target_match.group("seq_ID")
        donor_match = donor_regexp.match(line)
        if donor_match is not None:
            donor_ID = donor_match.group("seq_ID")
            if donor_ID != target_ID:
                if get_names_only:
                    new_donor = Donor(donor_ID, -1, -1)
                else:
                    donor_seq = record_dict[donor_ID]
                    donor_start, donor_end = coordinate_adjust(
                        donor_seq,
                        int(donor_match.group("start")),
                        int(donor_match.group("end")),
                    )
                    new_donor = Donor(donor_ID, donor_start, donor_end)
                donor_list.append(new_donor)
    if len(donor_list) > 0:
        result.add_donors(target_ID, (donor_list))
    return result
