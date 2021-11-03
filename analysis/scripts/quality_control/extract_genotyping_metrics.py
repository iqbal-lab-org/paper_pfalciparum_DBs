import click
import sys
from typing import Set, List
from copy import deepcopy
import warnings

from pysam import VariantFile, VariantRecord
from edlib import align as edlib_align

from common_utils.qc import QCMeasure
from common_utils.genome_region import RegionMode, GenomeRegion, genome_regions_from_bed


class VCFError(Exception):
    pass


def print_header(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    print(QCMeasure.get_header(), file=sys.stdout)
    ctx.exit()


@click.command()
@click.option(
    "--header_only",
    is_flag=True,
    help="Print the tsv header only",
    is_eager=True,
    callback=print_header,
)
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
    "--preset", help="Will extract known FORMAT fields for the supported genotyper."
)
@click.option(
    "--format_fields",
    "-f",
    help="The format fields to extract, in the form FIELD1:FIELD2:[...]",
)
def main(
    header_only, vcf_fname, bed_fname, sample_name, out_fname, preset, format_fields
):
    if preset is None and format_fields is None:
        print("Must provide either preset or format_fields", file=sys.stderr)
        exit(1)

    if preset is not None:
        format_fields = get_preset(preset)
    used_formats = set(format_fields.split(":"))
    if "COV" in used_formats:
        warnings.warn(
            "COV field extraction will assume haploid genotyping (Reports first called allele)"
        )
    check_headers(vcf_fname, used_formats, sample_name)

    result = list()
    region_list = genome_regions_from_bed(bed_fname)
    open_vcf = VariantFile(vcf_fname)

    if preset is None:
        preset = ""
    for region in region_list:
        result.extend(extract_format_fields_in_region(open_vcf, region, sample_name, preset, used_formats))
        result.extend(extract_called_variation_in_region(open_vcf, region, sample_name, preset))

    for res in result:
        print(str(res), file=out_fname)


def get_preset(preset: str) -> str:
    supported_presets = {"gram": "DP:COV:GT_CONF_PERCENTILE"}
    used_preset = None
    for s_p in supported_presets:
        if preset.startswith(s_p):
            used_preset = s_p
    if used_preset is None:
        raise ValueError(f"{preset} preset is not supported")
    return supported_presets[used_preset]


def check_headers(vcf_fname: str, format_fields: Set[str], sample_name: str) -> None:
    open_vcf = VariantFile(vcf_fname)
    formats = set(list(open_vcf.header.formats))
    if not formats.issuperset(format_fields):
        raise VCFError(
            f"Requested format fields: {format_fields} are not all present "
            f"in vcf file {vcf_fname}, which has fields: {formats}"
        )
    if sample_name not in open_vcf.header.samples:
        raise VCFError(f"{sample_name} not found in vcf")

def extract_called_variation_in_region(open_vcf: VariantFile, region: GenomeRegion, sample_name: str, preset: str) -> List[QCMeasure]:
    result = list()
    QCTemplate = QCMeasure(
        sample=sample_name,
        chrom=region.chrom,
        start=region.start,
        end=region.end,
        gene=region.gene,
        tool=preset,
        metric_category="genotyping",
    )
    num_alt_calls, ref_bases_in_alt_calls = 0, 0
    call_edit_dist = 0
    for record in open_vcf.fetch(region.chrom, int(region.start), int(region.end)):
        sample_call = record.samples[QCTemplate.sample]
        gt_idx = sample_call["GT"][0]
        if gt_idx == 0 or gt_idx is None:
            continue
        ref_allele = record.alleles[0]
        alt_allele = record.alleles[gt_idx]
        if ref_allele == alt_allele:
            # Can happen in gramtools due to ambiguous graphs
            continue
        num_alt_calls += 1
        ref_bases_in_alt_calls += len(ref_allele)
        call_edit_dist += edlib_align(ref_allele, alt_allele)["editDistance"]
    frac_ref_bases_in_calls = round(ref_bases_in_alt_calls / len(region), 3)
    scaled_called_eddist = round(call_edit_dist / len(region), 3)
    for metric, value in dict(num_alt_calls=num_alt_calls, frac_ref_bases_in_calls=frac_ref_bases_in_calls, scaled_called_eddist=scaled_called_eddist).items():
        newQCMeasure = deepcopy(QCTemplate)
        newQCMeasure.update(metric=metric, value=value)
        result.append(newQCMeasure)
    return result


def extract_format_fields_in_region(open_vcf: VariantFile, region: GenomeRegion, sample_name: str, preset: str, used_formats) -> List[QCMeasure]:
    result = list()
    for record in open_vcf.fetch(region.chrom, int(region.start), int(region.end)):
        QCTemplate = QCMeasure(
            sample=sample_name,
            chrom=region.chrom,
            start=record.start,
            end=record.stop,
            gene=region.gene,
            tool=preset,
            metric_category="genotyping",
        )
        result.extend(
            extract_format_fields_in_record(record, used_formats, QCTemplate, preset)
        )
    return result

def extract_format_fields_in_record(
    record: VariantRecord, used_formats: Set[str], QCTemplate: QCMeasure, preset=None
) -> List[QCMeasure]:
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
            if preset == "gram":
                frs_value = sample_call["COV"][gt_idx] / sample_call["DP"]
                frs_measure = deepcopy(QCTemplate)
                frs_measure.update(metric="FRS", value=frs_value)
                result.append(frs_measure)
        new_measure.update(metric=used_format, value=value)
        result.append(new_measure)
    return result


if __name__ == "__main__":
    main()
