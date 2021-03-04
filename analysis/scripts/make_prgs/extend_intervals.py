"""
Description
--------------
Border-constrained extension of bed intervals.
Provided an input bed and an integer flank_size, each bed interval is extended by flank_size either side,
under the condition that it cannot overlap other intervals.

Assumption
------------
The bed is sorted by feature (col1) then start pos (col2)
"""
import sys
from pathlib import Path
from typing import Dict, List, TextIO


class DisjointInterval:
    """
    Bed intervals are [0-based, 0-based). (Equivalently, [0-based, 1-based])
    So [0,2) and [2,3) do not overlap.
    """

    def __init__(self, start: int, end: int, name: str):
        self.start, self.end, self.name = start, end, name

    def left_intersects(self, flank: int, other: "DisjointInterval"):
        return (self.start - flank) < other.end

    def extend_left(self, flank: int, other: "DisjointInterval" = None):
        if other is not None and self.left_intersects(flank, other):
            self.start = other.end
        else:
            self.start = max(0, self.start - flank)

    def right_intersects(self, flank: int, other: "DisjointInterval" = None):
        return (self.end + flank) > other.start

    def extend_right(self, flank: int, other: "DisjointInterval"):
        if other is not None and self.right_intersects(flank, other):
            self.end = other.start
        else:
            self.end += flank

    def __repr__(self):
        return f"{self.start}\t{self.end}\t{self.name}"

    def __eq__(self, other: "DisjointInterval"):
        return (
            self.start == other.start
            and self.end == other.end
            and self.name == other.name
        )


Features = Dict[str, List[DisjointInterval]]


def load_existing_features(input_stream: TextIO) -> Features:
    features: Features = dict()
    for line in input_stream:
        elems = line.strip("\n").split("\t")
        interval = DisjointInterval(int(elems[1]), int(elems[2]), elems[3])
        feature = elems[0]
        if feature not in features:
            features[feature] = [interval]
        else:
            last_end = features[feature][-1].end
            if interval.start < last_end:
                raise ValueError(
                    f"New interval start {interval.start} overlaps with previous interval ending at {last_end}"
                )
            features[feature].append(interval)
    return features


def extend_features(features: Features, flank_size: int) -> None:
    """Extends intervals in-place"""
    for feature, intervals in features.items():
        last = len(intervals) - 1
        for i, interval in enumerate(intervals):
            prev_interval, next_interval = None, None
            if i > 0:
                prev_interval = intervals[i - 1]
            if i < last:
                next_interval = intervals[i + 1]

            interval.extend_left(flank_size, prev_interval)
            interval.extend_right(flank_size, next_interval)


def main():
    ## Parse args
    if not 4 <= len(sys.argv) <= 5:
        print(
            f"Usage: {sys.argv[0]} in.bed flank_size out.bed [--force]\n"
            "Intervals need to be sorted by start position and non-verlapping"
        )
        exit(1)

    input_bed = Path(sys.argv[1])
    if not input_bed.exists():
        raise FileNotFoundError(input_bed)

    try:
        flank_size = int(sys.argv[2])
    except ValueError:
        raise ValueError("flank_size should be an int") from None

    force = True if sys.argv[-1] == "--force" else False

    output_bed = Path(sys.argv[3])
    if output_bed.exists() and not force:
        raise ValueError(f"{output_bed} already exists. Run with --force to overwrite")

    ## Load features
    input_stream = input_bed.open()
    features = load_existing_features(input_stream)
    input_stream.close()

    extend_features(features, flank_size)

    ## Write extended features
    output_stream = output_bed.open("w")
    for feature, intervals in features.items():
        for interval in intervals:
            output_stream.write(f"{feature}\t{str(interval)}\n")
    output_stream.close()


if __name__ == "__main__":
    main()
