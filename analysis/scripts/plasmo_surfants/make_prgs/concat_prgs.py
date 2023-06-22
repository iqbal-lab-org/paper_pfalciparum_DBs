import sys
import csv
import logging
from pathlib import Path
from typing import Dict, Tuple, List, TextIO
import itertools


def usage():
    print(f"usage: {sys.argv[0]} input_dir regions.bed output_file \n")
    print(
        f"dir is the directory containing the prgs, regions.bed lists the file prefixes in column 4"
    )
    exit(1)


PRG_Names = List[str]
PRG_Files = Dict[str, Path]
PRG_Ints = List[int]

ENDIANNESS = "little"
BYTES_PER_INT = 4


class PRGAggregationError(Exception):
    pass


class PRGDecodeError(Exception):
    pass


class Record:
    def __init__(self, translation: int, count: int):
        self.translation = translation
        self.count = count


class PRGAggregator:
    """
    Has counting for how many times a site (odd) marker is seen.
    make_prg produces sites as '5A6C6', so cannot be > 1.
    """

    def __init__(self):
        self.translations: Dict[str, Dict[int, Record]] = dict()
        self.next_allocated = 5

    def translate(self, ID: str, marker: int) -> int:
        if ID not in self.translations:
            self.translations[ID] = dict()

        if marker <= 4:
            raise PRGAggregationError(f"Marker {marker} is not >4")

        local_table = self.translations[ID]
        if marker % 2 == 0:
            site_ID = marker - 1
            if site_ID not in local_table:
                raise PRGAggregationError(
                    f"Error: {marker}'s site number {marker - 1} has never been seen"
                )
            return local_table[site_ID].translation + 1

        else:
            if marker in local_table:
                record = local_table[marker]
                if record.count >= 1:
                    raise PRGAggregationError(
                        f"Error: {marker} site number present >1 times in local PRG {ID}"
                    )
                else:
                    record.count += 1
                    return local_table[marker].translation + 1
            local_table[marker] = Record(self.next_allocated, 1)
            self.next_allocated += 2

            return local_table[marker].translation


def get_aggregated_prgs(agg: PRGAggregator, prg_files: PRG_Files) -> PRG_Ints:
    rescaled_prg_ints: PRG_Ints = list()
    for prg_name, prg_path in prg_files.items():
        logging.info(f"Processing: {prg_name}")
        with prg_path.open("rb") as f:
            all_bytes = f.read()
        for pos in range(0, len(all_bytes), BYTES_PER_INT):
            int_bytes = all_bytes[pos : pos + BYTES_PER_INT]

            decoded_int = int.from_bytes(int_bytes, ENDIANNESS)
            if decoded_int <= 0:
                raise PRGDecodeError(f"PRG marker {decoded_int} should be > 0")
            elif decoded_int <= 4:
                rescaled_prg_ints.append(decoded_int)
            else:
                rescaled_int = agg.translate(prg_name, decoded_int)
                rescaled_prg_ints.append(rescaled_int)
        logging.info(f"Cumulative len prg: {len(rescaled_prg_ints)}")
        logging.info(f"Cumulative num sites: {(agg.next_allocated - 3) // 2 - 1}\n")

    return rescaled_prg_ints


def load_prg_names(file_stream: TextIO) -> PRG_Names:
    prg_names: PRG_Names = list()
    reader = csv.reader(file_stream, delimiter="\t")
    for line in reader:
        prg_names.append(line[3])
    return prg_names


def get_file_names(base_dir: Path, prg_names: PRG_Names) -> PRG_Files:
    """Looks for, and checks we have, all the names in prg_names as files"""
    prg_files: PRG_Files = dict()
    non_var_dir = base_dir.parent / "nonvars"
    to_search = base_dir.iterdir()
    if non_var_dir.exists():
        to_search = itertools.chain(to_search, non_var_dir.iterdir())
    for child in to_search:
        if child.is_dir():
            continue
        if child.stem in prg_names:
            prg_files[child.stem] = child

    for prg_name in prg_names:
        if prg_name not in prg_files:
            raise FileNotFoundError(
                f"Error: {prg_name} required but file not found in {base_dir}\n"
            )

    reordered_dict = dict()
    for prg_name in prg_names:
        reordered_dict[prg_name] = prg_files[prg_name]

    return reordered_dict


def to_bytes(prg_ints: PRG_Ints) -> List[bytes]:
    output_bytes = map(
        lambda integer: integer.to_bytes(BYTES_PER_INT, ENDIANNESS), prg_ints
    )
    return b"".join(output_bytes)


def main():
    logging.basicConfig(level=logging.INFO)

    if len(sys.argv) != 4:
        usage()

    base_dir = Path(sys.argv[1]).resolve()
    if not base_dir.exists():
        print(f"{base_dir} not found")
        usage()

    regions = Path(sys.argv[2]).resolve()
    if not regions.exists():
        print(f"{regions} not found")
        usage()

    output_file = Path(sys.argv[3]).resolve()

    f_in = regions.open()
    prg_names = load_prg_names(f_in)
    f_in.close()

    prg_files: PRG_Files = get_file_names(base_dir, prg_names)

    agg = PRGAggregator()
    rescaled_prg_ints: PRG_Ints = get_aggregated_prgs(agg, prg_files)

    output_stream = output_file.open("wb")
    output_bytes = to_bytes(rescaled_prg_ints)
    output_stream.write(output_bytes)
    output_stream.close()


if __name__ == "__main__":
    main()
