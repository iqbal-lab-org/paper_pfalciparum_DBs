"""
See main() docstring for explanation
"""
from collections import namedtuple
from typing import Dict, Tuple, List
from pathlib import Path

from pysam import AlignmentFile
import click


class MultipleAlignmentsError(Exception):
    pass


Scores = namedtuple("Scores", ["NM", "AS", "MAPQ"])
GeneScores = Dict[str, Scores]
BaselineScores: Dict[str, GeneScores] = dict()

MaskDict = Dict[str, float]


def load_mask_bed(mask_bed: click.Path) -> MaskDict:
    """
    Loads the amount of overlap between an aligned region and a mask
    of low-quality regions.
    The key to the returned dict is a composite of alignment query, target and tool.
    """
    result = dict()
    if mask_bed is not None:
        with Path(mask_bed).open() as fin:
            for line in fin:
                rows = line.strip().split("\t")
                new_key = f"{rows[4]}_{rows[0]}_{rows[3]}"
                assert new_key not in result
                result[new_key] = float(rows[8])
    return result


def get_sample_and_tool_name(sam_fname: Path) -> Tuple[str, str]:
    sample = sam_fname.stem
    tool = sam_fname.parent.name
    return sample, tool


def load_gene_lengths(bed_fname: Path) -> Dict[str, int]:
    result = dict()
    with bed_fname.open() as fin:
        for line in fin:
            cols = line.split("\t")
            gene_name = get_gene_name(cols[3].strip())
            result[gene_name] = int(cols[2]) - int(cols[1]) + 1
    return result


def get_gene_name(gene_string: str):
    """Converts 'Name=<x>;ID=<y>' into <x>"""
    return gene_string.split(";")[0].split("=")[1]


def get_scores(sam_fname: Path, gene_lengths: Dict[str, int]) -> GeneScores:
    result: GeneScores = dict()
    samfile = AlignmentFile(sam_fname, "r")
    for read in samfile.fetch(until_eof=True):
        gene_name = get_gene_name(read.query_name)
        if gene_name in result:
            raise MultipleAlignmentsError(
                f"Several alignments to gene {gene_name} in file {sam_fname}"
            )
        try:
            NM = read.get_tag("NM")
            NM = NM / gene_lengths[gene_name]
            # Convert to positive so that high is bad, low is good, like for NM.
            AS = read.get_tag("AS") * -1
            MAPQ = read.mapping_quality
            result[gene_name] = Scores(NM, AS, MAPQ)
        except KeyError:  # Case: read not aligned
            result[gene_name] = Scores("NA", "NA", "NA")

    return result


def compute_delta(scores_1: Scores, scores_2: Scores) -> Scores:
    newattrs = dict()
    for attr in Scores._fields:
        score_1 = getattr(scores_1, attr)
        score_2 = getattr(scores_2, attr)
        if score_1 == "NA" or score_2 == "NA":
            newattrs[attr] = "NA"
        else:
            newattrs[attr] = getattr(scores_2, attr) - getattr(scores_1, attr)
    return Scores(**newattrs)


def get_delta_scores(baseline_genes: GeneScores, query_genes: GeneScores) -> GeneScores:
    """Computes the scores of genes in `query_genes` minus the scores of genes in `baseline_genes`"""
    if baseline_genes.keys() != query_genes.keys():
        raise ValueError(
            f"Cannot get score deltas for dicts without same keys: {baseline_genes} vs {query_genes}"
        )
    result: GeneScores = dict()
    for gene_name, scores_1 in baseline_genes.items():
        scores_2 = query_genes[gene_name]
        result[gene_name] = compute_delta(scores_1, scores_2)
    return result


def write_stats(
    sam_file_list: List[Path], output_stats: Path, gene_lengths, mask_overlaps: MaskDict
):

    with output_stats.open("w") as stats_file:
        fieldnames = [
            "sample",
            "gene",
            "tool",
            "NM",
            "AS",
            "MAPQ",
            "mask_overlap",
        ]
        stats_file.write("\t".join(fieldnames) + "\n")

        for sam_fname in sam_file_list:
            sample, tool = get_sample_and_tool_name(sam_fname)
            scores = get_scores(sam_fname, gene_lengths)
            for gene in scores:
                mask_key = f"{gene}_{sample}_{tool}"
                mask_overlap = mask_overlaps.get(mask_key, 0)
                row = [
                    sample,
                    gene,
                    tool,
                    scores[gene].NM,
                    scores[gene].AS,
                    scores[gene].MAPQ,
                    mask_overlap,
                ]
                stats_file.write("\t".join(map(str, row)) + "\n")


@click.command()
@click.argument(
    "input_dir",
    type=click.Path(exists=True),
)
@click.argument("input_bed", type=click.Path(exists=True))
@click.argument("output_file", type=str)
@click.option("--mask_bed", type=click.Path(exists=True), default=None)
def main(input_dir, input_bed, output_file, mask_bed):
    """
    Computes stats from sam alignments and writes the stats to a tsv file. Computed stats include: NM (edit distance between query and target, scaled by gene length), AS
    (alignment score), MAPQ (probability that alignment is not unique).


    INPUT_DIR should be a directory containing other directories, each of which contain .sam files. Each sam-containing directory is listed as a tool in the output file.",
    """
    output_file = Path(output_file).resolve()
    output_file.parent.mkdir(exist_ok=True)

    if not output_file.exists():
        sam_file_list = list()
        for nextpath in Path(input_dir).glob("*"):
            if not nextpath.is_dir():
                if nextpath.suffix == ".sam":
                    sam_file_list.append(nextpath)
            else:
                sam_file_list.extend(Path(nextpath).glob(f"*.sam"))
        if len(sam_file_list) == 0:
            print(f"Error: no .sam files under {input_dir}")
            exit(1)

        gene_lengths = load_gene_lengths(Path(input_bed))
        mask_overlaps = load_mask_bed(mask_bed)

        write_stats(sam_file_list, output_file, gene_lengths, mask_overlaps)
    else:
        print(f"Found existing {output_file}, nothing to do. Delete it to regenerate")


if __name__ == "__main__":
    main()
