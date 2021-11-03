import click
import sys
from math import sqrt
from typing import List, Tuple, Iterable

from pysam import depth as samtools_depth

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
        region = GenomeRegion(dash_r_region, RegionMode.DASH_R)
        samtools_command = ["-r", dash_r_region, bam_fname]
        add_depth_results(result, region, len(region), sample_name, samtools_command)
    else:
        region_list = genome_regions_from_bed(bed_fname)
        if collapse_bed:
            chrom = ";".join(map(lambda reg: reg.chrom, region_list))
            start = ";".join(map(lambda reg: reg.start, region_list))
            end = ";".join(map(lambda reg: reg.end, region_list))
            region_length = sum(map(lambda reg: len(reg), region_list))
            region = GenomeRegion(f"{chrom}:{start}-{end}",RegionMode.DASH_R_SLACKY)
            samtools_command = ["-b", bed_fname, bam_fname]
            add_depth_results(result, region, region_length, sample_name, samtools_command)
        else:
            for region in region_list:
                samtools_command = ["-r", region.to_dash_r(), bam_fname]
                add_depth_results(result, region, len(region), sample_name, samtools_command)
    for res in result:
        print(str(res),file=out_fname)


def mean_and_std(numbers: Iterable[int], region_length) -> Tuple[float, float]:
    total = 0
    squared_total = 0
    for number in numbers:
        total += number
        squared_total += number ** 2
    result_mean = total / region_length
    result_stdv = sqrt((squared_total - (total ** 2) / region_length) / region_length)
    return result_mean, result_stdv


def parse_samtools_depth(depth_string: str) -> Iterable[int]:
    for line in depth_string.strip().split("\n"):
        depth = int(line.split("\t")[2].strip())
        yield depth

def add_depth_results(result_list: List[QCMeasure], region: GenomeRegion, region_length, sample_name, samtools_command: List[str]):
        samtools_out = samtools_depth(*samtools_command)
        mean, stdv = mean_and_std(parse_samtools_depth(samtools_out), region_length)
        qc_fields = {"sample":sample_name,"start":region.start,"end":region.end,
                "chrom":region.chrom,"metric_category": "mapping_depth"}
        for metric, value in dict(mean=mean,stdv=stdv).items():
            qc_fields["metric"] = metric
            qc_fields["value"] = str(round(value,2))
            result_list.append(QCMeasure(**qc_fields))

if __name__ == "__main__":
    main()
