import sys
from pathlib import Path
import re

from pysam import VariantFile, stats as sam_stats, mpileup as sam_mpileup

from eval_varcalls.shift_to_induced_genome_coords import translate_bed


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
        "reads_mapped_and_paired",
        "reads_properly_paired",
        "bases_mapped_cigar",
        "base_error_rate",
        "fraction_disagreeing_pileup",
        "fraction_positions_0x_or_less",
        "fraction_positions_9x_or_less",
    ]

    def __init__(self, sample_name: str, tool_name: str, gene_name: str):
        self.gene = {attribute: "None" for attribute in self.attributes}
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


def populate_stats(record, input_bam, region) -> None:
    # Obtain mapping statistics in bed region via bam
    stats = sam_stats(str(input_bam), region, split_lines=True)
    for stat_line in stats:
        if not stat_line.startswith("SN"):
            continue
        if "reads mapped and paired" in stat_line:
            record["reads_mapped_and_paired"] = int(stat_line.split()[5])
        elif "reads properly paired" in stat_line:
            record["reads_properly_paired"] = int(stat_line.split()[4])
        elif "bases mapped (cigar):" in stat_line:
            record["bases_mapped_cigar"] = int(stat_line.split()[4])
        elif "error rate" in stat_line:
            record["base_error_rate"] = float(stat_line.split()[3])
        else:
            continue


PILEUP_INDEL_RE = "[+-]([0-9]+)[ACGTNacgtn*#]+"
PILEUP_MATCH = ".,"
PILEUP_MISMATCH = "ATCGNatcgn*"

pileup_matches = lambda char: char in PILEUP_MATCH
pileup_mismatches = lambda char: char in PILEUP_MISMATCH


def remove_indels(pileup_string: str) -> str:
    """
    See man samtools-mpileup for what the pileup column can contain. Here I remove indels, making sure that no more than the indel is clipped out.
    E.g.
    Input:  '.,2aca,.' # the second a is a mismatch, not part of the insertion description
    Output: '.,a,.'
    """
    new_string = ""
    cur_pos = 0
    for match in re.finditer(PILEUP_INDEL_RE, pileup_string):
        indel_length = int(match.group(1))
        new_string += pileup_string[cur_pos : match.start(0)]
        cur_pos = match.end(1) + indel_length
    new_string += pileup_string[cur_pos:]
    return new_string


def remove_read_ends(pileup_string: str) -> str:
    new_string = ""
    cur_pos = 0
    for match in re.finditer(r"(\^.)|(\$)", pileup_string):
        new_string += pileup_string[cur_pos : match.start(0)]
        cur_pos = match.end(0)
    new_string += pileup_string[cur_pos:]
    return new_string


def majority_pileup_is_non_ref(pileup) -> bool:
    matches = sum(map(pileup_matches, pileup))
    mismatches = sum(map(pileup_mismatches, pileup))
    frac_matches = matches / (matches + mismatches)
    return frac_matches < 0.5


def get_fraction_below_depth(depths, limit):
    num_below = len([depth for depth in depths if depth <= limit])
    return round(num_below / len(depths), 3)


def populate_pileup(record, input_bam, input_ref_genome, region) -> None:
    """
    Computes the fraction of positions that do not match, and fraction of positions under a given depth
    """
    num_disagreeing, num_total = 0, 0
    mpileup_output = sam_mpileup(
        str(input_bam),
        "-r",
        region,
        "-f",
        str(input_ref_genome),
        "-a",
        split_lines=True,
    )
    pileups = map(
        lambda line: remove_read_ends(remove_indels(line.split("\t")[4])),
        mpileup_output,
    )
    depths = list(map(lambda line: int(line.split("\t")[3]), mpileup_output))
    for pileup, depth in zip(pileups, depths):
        if depth == 0:
            continue
        num_total += 1
        if majority_pileup_is_non_ref(pileup):
            num_disagreeing += 1
    record["fraction_disagreeing_pileup"] = round(num_disagreeing / num_total, 3)
    for limit in [0, 9]:
        record[f"fraction_positions_{limit}x_or_less"] = get_fraction_below_depth(
            depths, limit
        )


if __name__ == "__main__":
    try:
        if sys.argv[1] == "--header":
            Stats.print_header()
            exit(0)
    except IndexError:
        pass

    if len(sys.argv) != 7:
        usage()

    input_bed = Path(sys.argv[1]).resolve()
    input_ref_genome = Path(sys.argv[2]).resolve()
    input_vcf = Path(sys.argv[3]).resolve()
    input_bam = Path(sys.argv[4]).resolve()
    output_tsv = Path(sys.argv[5]).resolve()
    metadata = sys.argv[6]

    for fname in (input_bed, input_vcf, input_ref_genome, input_bam):
        if not fname.exists():
            print(f"Error: {fname} required but not found")
            usage()

    try:
        sample_name, tool_name = metadata.split(":")
    except ValueError:
        print("Metadata formatting error")
        usage()

    vcf_file = VariantFile(str(input_vcf))
    if not Path(f"{input_vcf}.csi").exists():
        raise Exception(f"bcftools index for {input_vcf} not found")

    translated_bed = translate_bed(input_bed, input_ref_genome, input_vcf)

    output_file = output_tsv.open("w")

    for bed_line in translated_bed:
        gene_name = bed_line[3]
        reg_start = int(bed_line[1]) + 1
        reg_end = int(bed_line[2])
        record = Stats(sample_name, tool_name, gene_name)
        region = f"{bed_line[0]}:{reg_start}-{reg_end}"

        populate_stats(record, input_bam, region)
        populate_pileup(record, input_bam, input_ref_genome, region)
        print(record, file=output_file)
    output_file.close()
