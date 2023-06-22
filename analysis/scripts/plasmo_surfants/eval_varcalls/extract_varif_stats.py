import click
import sys

from pysam import VariantFile, VariantRecord
from edlib import align as edlib_align

from plasmo_surfants.common_utils.metrics import MetricsRecorder
from plasmo_surfants.common_utils.genome_region import RegionMode, GenomeRegion, genome_regions_from_bed

class VarifStats(MetricsRecorder):
    _headers = [
            "chrom",
            "start",
            "stop",
            "gene",
            "sample",
            "tool",
            "eval_categ",
            "ed_dist_RA",
            "ed_dist_TR",
            "ed_dist_TA",
            "eval_categ_ed_dist_numerator",
            "eval_categ_ed_dist_denominator",
            ]

def print_header(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    print(VarifStats.get_header(), file=sys.stdout)
    ctx.exit()

def extract_sample_record(record, sample_name):
    try:
        sample_record = record.samples[sample_name]
    except KeyError:
        sample_record = record.samples["sample"]
    return sample_record

def extract_RA(record, sample_record):
    GT = sample_record["GT"][0]
    ref = record.alleles[0]
    alt = record.alleles[GT]
    return edlib_align(ref,alt)["editDistance"]

def extract_recall_metrics(record, sample_record):
    if "VFR_ED_TA" not in sample_record:
        return None
    result = dict(
            ed_dist_RA = extract_RA(record,sample_record),
            ed_dist_TA = int(sample_record["VFR_ED_TA"])
            )
    result["eval_categ_ed_dist_numerator"] = result["ed_dist_RA"] - result["ed_dist_TA"]
    result["eval_categ_ed_dist_denominator"] = result["ed_dist_RA"]
    return result

def extract_precision_metrics(record, sample_record):
    if "VFR_ED_TA" not in sample_record:
        return None
    result = dict(
            ed_dist_RA = extract_RA(record,sample_record),
            ed_dist_TA = int(sample_record["VFR_ED_TA"]),
            )
    try:
        ed_dist_TR = int(sample_record["VFR_ED_TR"])
    except KeyError:
        if result["ed_dist_TA"] == 0:
            ed_dist_TR = result["ed_dist_RA"]
        else:
            return None
    result["ed_dist_TR"] = ed_dist_TR
    result["eval_categ_ed_dist_numerator"] = result["ed_dist_TR"] - result["ed_dist_TA"]
    result["eval_categ_ed_dist_denominator"] = result["ed_dist_TR"]
    return result


@click.command()
@click.option(
    "--header_only",
    is_flag=True,
    help="Print the tsv header only",
    is_eager=True,
    callback=print_header,
)
@click.argument("recall_vcf_fname", type=click.Path(exists=True))
@click.argument("precision_vcf_fname", type=click.Path(exists=True))
@click.argument("bed_fname", type=click.Path(exists=True))
@click.option("--sample_name", "-s", required=True)
@click.option("--tool_name", "-t", required=True)
@click.option(
    "--out_fname",
    type=click.File("w"),
    default=sys.stdout,
    help="File to write to. Default: stdout",
)
def main(
    header_only, recall_vcf_fname, precision_vcf_fname, bed_fname, sample_name,
    tool_name, out_fname
):
    result = list()
    region_list = genome_regions_from_bed(bed_fname)
    for eval_categ, extracter_function, vcf_fname in zip(
            ["recall","precision"],[extract_recall_metrics, extract_precision_metrics],
            [recall_vcf_fname, precision_vcf_fname]
            ):
        open_vcf = VariantFile(vcf_fname)
        for region in region_list:
            for record in open_vcf.fetch(region.chrom, int(region.start), int(region.end)):
                stats = VarifStats(
                        chrom=region.chrom,
                        start=record.start,
                        stop=record.stop,
                        gene=region.gene,
                        sample=sample_name,
                        tool=tool_name,
                        eval_categ=eval_categ,
                        )
                sample_record = extract_sample_record(record, sample_name)
                if sample_record["VFR_RESULT"] == "FP_PROBE_UNMAPPED":
                    continue
                categ_specific_stats = extracter_function(record, sample_record)
                if categ_specific_stats is None:
                    continue
                stats.update(**categ_specific_stats)
                result.append(stats)
    for res in result:
        print(str(res), file=out_fname)

if __name__ == "__main__":
    main()
