import click
import sys
from typing import Set, List
from copy import deepcopy
import warnings

from pysam import VariantFile, VariantRecord

from common_utils.qc import QCMeasure
from common_utils.genome_region import RegionMode, GenomeRegion, genome_regions_from_bed

class VCFError(Exception):
    pass

@click.command()
@click.argument("vcf_fname", type=click.Path(exists=True))
@click.argument("bed_fname", type=click.Path(exists=True))
@click.option("--sample_name", "-s", required=True)
@click.option(
    "--out_fname",
    type=click.File("w"),
    default=sys.stdout,
    help="File to write to. Default: stdout",
)
@click.option(
    "--preset",
    type=click.Choice(["gramtools"]),
    help="Will extract known FORMAT fields for the supported genotyper."
)
@click.option(
    "--format_fields",
    "-f",
    help="The format fields to extract, in the form FIELD1:FIELD2:[...]"
)
def main(
    vcf_fname,
    bed_fname,
    sample_name,
    out_fname,
    preset,
    format_fields
):
    if preset is None and format_fields is None:
        print("Must provide either preset or format_fields", file=sys.stderr)
        exit(1)

    if preset is not None:
        format_fields = get_preset(preset)
    used_formats = set(format_fields.split(":"))
    if "COV" in used_formats:
        warnings.warn("COV field extraction will assume haploid genotyping (Reports first called allele)")
    check_headers(vcf_fname, used_formats, sample_name)

    result = list()
    region_list = genome_regions_from_bed(bed_fname)
    open_vcf = VariantFile(vcf_fname)

    if preset is None:
        preset = ""
    for region in region_list:
        for record in open_vcf.fetch(region.chrom, int(region.start), int(region.end)):
            QCTemplate = QCMeasure(sample=sample_name, chrom=region.chrom,
                    start=record.start, end=record.stop, gene=region.gene, tool=preset, 
                    metric_category="genotyping")
            result.extend(extract_genotyping_metrics(record, used_formats, QCTemplate,
                preset))

    for res in result:
        print(str(res),file=out_fname)

def get_preset(preset: str) -> str:
    supported_presets = {"gramtools": "DP:COV:GT_CONF_PERCENTILE"}
    if preset not in supported_presets:
        raise ValueError(f"{preset} preset is not supported")
    return supported_presets[preset]

def check_headers(vcf_fname: str, format_fields: Set[str], sample_name: str) -> None:
    open_vcf = VariantFile(vcf_fname)
    formats = set(list(open_vcf.header.formats))
    if not formats.issuperset(format_fields):
        raise VCFError(f"Requested format fields: {format_fields} are not all present " 
                f"in vcf file {vcf_fname}, which has fields: {formats}"
                )
    if sample_name not in open_vcf.header.samples:
        raise VCFError(f"{sample_name} not found in vcf")

def extract_genotyping_metrics(record: VariantRecord, used_formats: Set[str],
        QCTemplate: QCMeasure, preset = None) -> List[QCMeasure]:
    result = list()
    sample_call = record.samples[QCTemplate.sample]
    gt_idx = sample_call["GT"][0]
    for used_format in used_formats:
        new_measure = deepcopy(QCTemplate)
        try:
            value = sample_call[used_format]
        except:
            continue
        if used_format == "COV":
            value = sample_call["COV"][gt_idx]
            if preset == "gramtools":
                frs_value = sample_call["COV"][gt_idx] / sample_call["DP"]
                frs_measure = deepcopy(QCTemplate)
                frs_measure.update(metric="FRS", value=frs_value)
                result.append(frs_measure)
        new_measure.update(metric=used_format, value=value)
        result.append(new_measure)
    return result

if __name__ == "__main__":
    main()
