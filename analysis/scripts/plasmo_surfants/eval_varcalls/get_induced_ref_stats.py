import sys
from pathlib import Path
from typing import Tuple
import re

import click
from pysam import (
    AlignmentFile,
    stats as sam_stats,
)
from numpy import absolute, mean, std

from plasmo_surfants.common_utils.shift_to_induced_genome_coords import translate_bed
from plasmo_surfants.common_utils.pileup_utils import (
        get_mpileups, 
        pileup_matches,
        pileup_mismatches
)


def Round(value):
    return round(value, 5)

def print_header(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    print("\t".join(Stats.attributes), file=sys.stdout)
    ctx.exit()


def usage():
    print(
        f"usage: {sys.argv[0]} input.bed input.fa[.gz] input.vcf[.gz] input.bam output.tsv metadata\n"
        "metadata is of form <sample_name>:<tool_name>\n"
        f"To produce output tsv header only, run: {sys.argv[0]} --header"
    )
    exit(1)


class Stats:
    """
    Records induced_reference mapping and variant statistics
    """

    attributes = [
        "sample",
        "tool",
        "gene",
        "bases_mapped_cigar",
        "base_error_rate",
        "fraction_disagreeing_pileup_min1x",
        "fraction_disagreeing_pileup_min5x",
        "fraction_disagreeing_pileup_min10x",
        "fraction_positions_0x_or_less",
        "fraction_positions_4x_or_less",
        "fraction_positions_9x_or_less",
        "fraction_reads_paired",
        "fraction_reads_mapped_and_paired",
        "fraction_reads_paired_one_unmapped",
        "fraction_reads_mapped_and_paired_mapq_1_or_more",
        "fraction_reads_mapped_and_paired_mapq_10_or_more",
        "fraction_reads_mapped_and_paired_mapq_30_or_more",
        "fraction_reads_properly_paired_aligner",
        "fraction_reads_above_mean_ins_size_plus_two_std",
        "fraction_reads_below_mean_ins_size_minus_two_std",
        "max_insert_size",
        "mean_read_coverage",
        "std_read_coverage",
    ]

    def __init__(self, sample_name: str, tool_name: str, gene_name: str):
        self.gene = {attribute: "NA" for attribute in self.attributes}
        self.gene["sample"] = sample_name
        self.gene["tool"] = tool_name
        self.gene["gene"] = gene_name

    def __getitem__(self, key):
        return self.gene.get(key)

    def __setitem__(self, key, value):
        if key not in self.attributes:
            raise KeyError(f"{key} not in {self.attributes}")
        self.gene[key] = value

    @classmethod
    def print_header(cls):
        print("\t".join(cls.attributes))

    def __repr__(self):
        return "\t".join(map(str, [self.gene[key] for key in self.attributes]))


def extract_read_inserts_SN(input_bam, bed_fname) -> Tuple[float, float]:
    stats = sam_stats(str(input_bam), "-t",bed_fname, split_lines=True)
    mean, std = None, None
    for stat_line in stats:
        if not stat_line.startswith("SN"):
            continue
        if "insert size average" in stat_line:
            mean = float(stat_line.split()[4])
        if "insert size standard deviation" in stat_line:
            std = float(stat_line.split()[5])
        if mean is not None and std is not None:
            break
    return mean, std


def populate_stats_SN(record, input_bam, region) -> None:
    # Obtain mapping statistics in bed region via bam
    stats = sam_stats(str(input_bam), region, split_lines=True)
    for stat_line in stats:
        if not stat_line.startswith("SN"):
            continue
        if "bases mapped (cigar):" in stat_line:
            record["bases_mapped_cigar"] = int(stat_line.split()[4])
        elif "error rate" in stat_line:
            record["base_error_rate"] = float(stat_line.split()[3])
        else:
            continue



def majority_pileup_is_non_ref(pileup) -> bool:
    matches = sum(map(pileup_matches, pileup))
    mismatches = sum(map(pileup_mismatches, pileup))
    frac_matches = matches / (matches + mismatches)
    return frac_matches < 0.5


def get_fraction_leq_limit(values, limit):
    """
    leq: less than or equal to
    """
    num_below = len([value for value in values if value <= limit])
    return Round(num_below / len(values))


def populate_stats_pileup(record, input_bam, input_ref_genome, region) -> None:
    """
    Computes the fraction of positions that do not match, and fraction of positions under a given depth
    """
    pileups, depths, _ = get_mpileups(input_bam, input_ref_genome, region)
    record["mean_read_coverage"] = Round(mean(depths))
    record["std_read_coverage"] = Round(std(depths))
    if len(depths) > 0:
        for limit in [0, 4, 9]:
            record[f"fraction_positions_{limit}x_or_less"] = get_fraction_leq_limit(
                depths, limit
            )
    num_disagreeing = {1:0,5:0,10:0}
    num_total = 0
    for pileup, depth in zip(pileups, depths):
        if depth == 0:
            continue
        num_total += 1
        if majority_pileup_is_non_ref(pileup):
            for limit in num_disagreeing:
                if depth >= limit:
                    num_disagreeing[limit] += 1
    # Condition controls for no mapped reads in region
    if num_total == 0:
        return
    for limit in num_disagreeing:
        record[f"fraction_disagreeing_pileup_min{limit}x"] = Round(num_disagreeing[limit] / num_total)


def populate_stats_reads(
    record, input_bam, chrom, start, end, insert_mean: float, insert_std: float
) -> None:
    alignment_file = AlignmentFile(input_bam)
    read_stats = {metric: 0 for metric in Stats.attributes[11:19]}
    mapqs = list()
    tlens = list()
    max_tlen = 0
    num_reads = 0
    for read in alignment_file.fetch(chrom, start, end):
        num_reads += 1
        read_paired = read.is_paired
        read_mapped = not read.is_unmapped
        mate_mapped = not read.mate_is_unmapped
        if read_paired:
            read_stats["fraction_reads_paired"] += 1
        if read_paired and read_mapped and mate_mapped:
            read_stats["fraction_reads_mapped_and_paired"] += 1
            mapqs.append(read.mapping_quality)
            abs_tlen = absolute(read.tlen)
            if abs_tlen > max_tlen:
                max_tlen = abs_tlen
            tlens.append(abs_tlen)
        if read_paired and read.is_proper_pair and read_mapped:
            read_stats["fraction_reads_properly_paired_aligner"] += 1
        if read_paired and read_mapped and not mate_mapped:
            read_stats["fraction_reads_paired_one_unmapped"] += 1
    if num_reads == 0:
        return
    for metric in read_stats:
        read_stats[metric] = Round(read_stats[metric] / (num_reads + 1))
    for limit in [0, 9, 29]:
        metric_name = f"fraction_reads_mapped_and_paired_mapq_{limit + 1}_or_more"
        if len(mapqs) > 0:
            read_stats[metric_name] = Round(1 - get_fraction_leq_limit(mapqs, limit))
    for metric in read_stats:
        record[metric] = read_stats[metric]
    if insert_mean is not None:
        if max_tlen > 0:
            record["max_insert_size"] = max_tlen
        if len(tlens) > 0:
            record[
                "fraction_reads_above_mean_ins_size_plus_two_std"
            ] = Round(1 - get_fraction_leq_limit(tlens, insert_mean + 2 * insert_std))
            record[
                "fraction_reads_below_mean_ins_size_minus_two_std"
            ] = get_fraction_leq_limit(tlens, insert_mean - 2 * insert_std)


@click.command()
@click.argument("bed_fname", type=click.Path(exists=True))
@click.argument("ref_fname", type=click.Path(exists=True))
@click.argument("induced_ref_fname", type=click.Path(exists=True))
@click.argument("vcf_fname", type=click.Path(exists=True))
@click.argument("bam_fname", type=click.Path(exists=True))
@click.option("--bed_for_insert_size", "-bi", type=click.Path(exists=True))
@click.option("--sample_name", required=True)
@click.option("--tool_name", required=True)
@click.option(
    "--out_fname",
    default=sys.stdout,
    help="File to write to. Default: stdout",
)
@click.option(
    "--header_only",
    is_flag=True,
    help="Print the tsv header only",
    is_eager=True,
    callback=print_header,
)
def main(
    bed_fname,
    ref_fname,
    induced_ref_fname,
    vcf_fname,
    bam_fname,
    bed_for_insert_size,
    out_fname,
    header_only,
    sample_name,
    tool_name,
):
    if not Path(f"{vcf_fname}.csi").exists():
        raise Exception(f"bcftools index for {vcf_fname} not found")

    translated_bed = translate_bed(Path(bed_fname), ref_fname, vcf_fname)

    log_file = Path(str(out_fname)).with_suffix(".translated_bed.txt").open("w")
    out_fhandle = open(out_fname, "w")

    insert_mean, insert_std = None, None
    if bed_for_insert_size is not None:
        insert_mean, insert_std = extract_read_inserts_SN(
            bam_fname, bed_for_insert_size
        )

    for bed_line in translated_bed:
        log_file.write("\t".join(map(str,bed_line)))
        gene_name = bed_line[3].strip()
        reg_start = int(bed_line[1]) + 1
        reg_end = int(bed_line[2])
        record = Stats(sample_name, tool_name, gene_name)
        region = f"{bed_line[0]}:{reg_start}-{reg_end}"

        populate_stats_SN(record, bam_fname, region)
        populate_stats_pileup(record, bam_fname, induced_ref_fname, region)
        populate_stats_reads(
            record, bam_fname, bed_line[0], reg_start, reg_end, insert_mean, insert_std
        )
        print(record, file=out_fhandle)

    out_fhandle.close()


if __name__ == "__main__":
    main()
