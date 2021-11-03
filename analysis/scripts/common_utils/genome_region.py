from enum import Enum, auto
import re
from typing import List


class RegionMode(Enum):
    BED = "(\w+)\t([0-9]+)\t([0-9]+)\t(\w+).*"
    DASH_R = "(\w+):([0-9]+)-([0-9]+)"
    DASH_R_SLACKY = "([\w;]+):([\w;]+)-([\w;]+)"


class GenomeRegion:
    def __init__(self, region_str: str, mode: RegionMode):
        assert type(mode) is RegionMode, f"mode parameter must be of class {RegionMode}"
        match = re.match(mode.value, region_str)
        if match is None:
            raise ValueError(f"{region_str} does not match {mode.value}")
        groups = match.groups()
        self.gene = ""
        self.chrom = groups[0]
        self.start = groups[1]
        self.end = groups[2]
        if mode is RegionMode.BED:
            self.gene = groups[3]
            self.start = str(int(groups[1]) + 1)

    def to_dash_r(self):
        return f"{self.chrom}:{self.start}-{self.end}"

    def __len__(self):
        return int(self.end) - int(self.start) + 1


GenomeRegions = List[GenomeRegion]


def genome_regions_from_bed(bed_fname) -> GenomeRegions:
    result = list()
    with open(bed_fname) as fhandle_in:
        for line in fhandle_in:
            result.append(GenomeRegion(line, RegionMode.BED))
    return result
