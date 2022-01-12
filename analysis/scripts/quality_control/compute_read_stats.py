import click
import sys
from math import sqrt
from typing import List, Tuple, Iterable, Optional
import warnings

from pysam import (
    AlignmentFile,
    mpileup as sam_mpileup,
)
from numpy import absolute

from common_utils.qc import QCMeasure
from common_utils.genome_region import RegionMode, GenomeRegion, genome_regions_from_bed


def print_header(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    print(QCMeasure.get_header(), file=sys.stdout)
    ctx.exit()


@click.command()
@click.argument("bam_fname", type=click.Path(exists=True))
@click.option("--bed_fname", "-b", type=click.Path(exists=True))
@click.option("--sample_name", "-s", required=True)
@click.option(
    "--out_fname",
    type=click.File("w"),
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
@click.option(
    "--collapse_bed",
    is_flag=True,
    help="Compute the metric across all the bed features, instead of once for"
    " each feature",
)
@click.option(
    "--region",
    "-r",
    "dash_r_region",
    help="Specify a region (chrom:start-stop) instead of a bed file. "
    "Overrides any specified bed file input.",
)
def main(
    bam_fname,
    bed_fname,
    sample_name,
    out_fname,
    header_only,
    collapse_bed,
    dash_r_region,
):
    if bed_fname is None and dash_r_region is None:
        print("Must provide either bed filename or region", file=sys.stderr)
        exit(1)

    result = list()
    if dash_r_region is not None:
        region_list = [GenomeRegion(dash_r_region, RegionMode.DASH_R)]
    else:
        region_list = genome_regions_from_bed(bed_fname)
    if collapse_bed:
        chrom = ";".join(map(lambda reg: reg.chrom, region_list))
        start = ";".join(map(lambda reg: reg.start, region_list))
        end = ";".join(map(lambda reg: reg.end, region_list))
        region_for_naming = GenomeRegion(
            f"{chrom}:{start}-{end}", RegionMode.DASH_R_SLACKY
        )
        add_read_stats(result, region_list, sample_name, bam_fname, region_for_naming)
    else:
        for region in region_list:
            add_read_stats(result, [region], sample_name, bam_fname)
    for res in result:
        print(str(res), file=out_fname)


def mean_and_std(numbers: Iterable[int]) -> Tuple[float, float]:
    total = 0
    squared_total = 0
    region_length = 0
    for number in numbers:
        if number is None:
            continue
        region_length += 1
        total += number
        squared_total += number ** 2
    if region_length == 0:
        warnings.warn("Region is of size 0")
        return 0, 0
    result_mean = total / region_length
    result_stdv = sqrt((squared_total - (total ** 2) / region_length) / region_length)
    return result_mean, result_stdv


def get_pileup_depths(depth_strings: List[str]) -> Iterable[int]:
    for depth_string in depth_strings:
        depth = int(depth_string.split("\t")[3].strip())
        yield depth


def get_pileup_base_qual(depth_strings: List[str]) -> Iterable[int]:
    for depth_string in depth_strings:
        base_quals = depth_string.split("\t")[5].strip()
        for base_qual in base_quals:
            yield ord(base_qual) - 33


def add_read_stats(
    result_list: List[QCMeasure],
    regions: List[GenomeRegion],
    sample_name,
    bam_fname,
    region_for_naming: Optional[GenomeRegion] = None,
):
    """
    Extracts (mean, std) of read depth and base qualities, across all `regions`
    """
    pileup_out = list()
    read_lengths, tlens = list(), list()
    alignment_file = AlignmentFile(bam_fname)
    for region in regions:
        mpileup_command = ["-a", "-r", region.to_dash_r(), bam_fname]
        pileup_out += sam_mpileup(*mpileup_command, split_lines=True)
        for read in alignment_file.fetch(region.chrom, int(region.start), int(region.end)):
            read_lengths.append(read.infer_read_length())
            if read.tlen != 0:
                tlens.append(absolute(read.tlen))

    metrics = {
        "mapping_depth": get_pileup_depths,
        "base_quality": get_pileup_base_qual,
        "read_length": lambda _: read_lengths,
        "insert_size": lambda _: tlens,
        "insert_size_below_10000": lambda _: filter(lambda tlen: tlen < 10000, tlens)
    }
    if region_for_naming is None:
        region_for_naming = region
    for metric_category, metric_function in metrics.items():
        qc_fields = {
            "sample": sample_name,
            "start": region_for_naming.start,
            "end": region_for_naming.end,
            "chrom": region_for_naming.chrom,
            "metric_category": metric_category,
        }
        mean, stdv = mean_and_std(metric_function(pileup_out))
        for metric, value in dict(mean=mean, stdv=stdv).items():
            qc_fields["metric"] = metric
            qc_fields["value"] = str(round(value, 2))
            result_list.append(QCMeasure(**qc_fields))


if __name__ == "__main__":
    main()
