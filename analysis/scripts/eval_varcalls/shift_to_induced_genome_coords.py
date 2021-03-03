"""
Expects an input bed and an input VCF and the reference genome fasta it refers to.

Produces an output bed whose features get mapped to the coordinates of the genome fasta
obtained by applying the calls in the vcf to the reference genome fasta.
"""
import sys
from pathlib import Path

from pysam import VariantFile

from gramtools.commands.genotype.seq_region_map import (
    SeqRegionMapper,
    SearchableSeqRegionsMap,
    BisectTarget,
    Chrom,
)
from gramtools.commands.common import load_fasta


def translate_pos(
    chrom: Chrom, pos: int, searchable_map: SearchableSeqRegionsMap
) -> int:
    region_idx = searchable_map.bisect(chrom, pos, BisectTarget.BASE_REF)
    region = searchable_map.get_region(chrom, region_idx)
    base_ref_offset = pos - region.base_ref_start
    translation = region.pers_ref_start + base_ref_offset
    return translation


def translate_bed(input_bed: Path, input_ref_genome: Path, input_vcf: Path):
    """
    Returns a list of list, where each inner list is a bed line with start and end positions translated to be expressed on the induced reference genome.
    """
    chrom_sizes = load_fasta(input_ref_genome, sizes_only=True)
    vcf_records = VariantFile(input_vcf).fetch()
    try:
        searchable_map = SearchableSeqRegionsMap(
            SeqRegionMapper(vcf_records, chrom_sizes).get_map()
        )
    # Catches case where vcf has no records
    except ValueError:
        searchable_map = None
    result = list()
    with input_bed.open("r") as bed_in:
        for line in bed_in:
            rows = line.split("\t")
            chrom = rows[0]
            start_pos = (
                int(rows[1]) + 1
            )  # Bed start is 0-based, SeqRegion coords are 1-based
            end_pos = int(rows[2])

            if searchable_map is not None:
                translated_start_pos = str(
                    translate_pos(chrom, start_pos, searchable_map) - 1
                )
                translated_end_pos = translate_pos(
                    chrom, end_pos, searchable_map
                ).__str__()
            else:
                translated_start_pos = start_pos
                translated_end_pos = end_pos
            result.append([chrom, translated_start_pos, translated_end_pos] + rows[3:])
    return result


def usage():
    print(f"usage: {sys.argv[0]} input.bed input.fa[.gz] input.vcf[.gz] output.bed")
    exit(1)


if __name__ == "__main__":
    if len(sys.argv) != 5:
        usage()

    input_bed = Path(sys.argv[1]).resolve()
    input_ref_genome = Path(sys.argv[2]).resolve()
    input_vcf = Path(sys.argv[3]).resolve()
    output_bed = Path(sys.argv[4]).resolve()

    for fname in (input_bed, input_vcf, input_ref_genome):
        if not fname.exists():
            print(f"Error: {fname} required but not found")
            usage()

    output_bed.parent.mkdir(exist_ok=True)

    translated_bed = translate_bed(input_bed, input_ref_genome, input_vcf)
    with output_bed.open("w") as bed_out:
        for line in translated_bed:
            bed_out.write("\t".join(line))
