import sys
from pathlib import Path

from pysam import VariantFile, stats as sam_stats
from edlib import align as edalign

from shift_to_induced_genome_coords import translate_bed


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
        "num_alt_calls",
        "frac_ref_bases_in_alt_calls",
        "induced_ref_scaled_eddist",
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

    translated_bed = translate_bed(input_bed, input_ref_genome, input_vcf)

    vcf_file = VariantFile(str(input_vcf))
    output_file = output_tsv.open("w")

    for bed_line in translated_bed:
        gene_name = bed_line[3].split(";")[0].split("=")[1]
        reg_start = int(bed_line[1]) + 1
        reg_end = int(bed_line[2])
        record = Stats(sample_name, tool_name, gene_name)
        region = f"{bed_line[0]}:{reg_start}-{reg_end}"

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

        # Traverse vcf records in bed region
        vcf_records = vcf_file.fetch(bed_line[0], reg_start, reg_end)
        num_alt_calls, ref_bases_in_alt_calls = 0, 0
        edit_distance_induced_ref = 0
        for vcf_record in vcf_records:
            sample = vcf_record.samples[0]
            sample_gt = sample["GT"]
            if any(map(lambda el: el != 0 and el is not None, sample_gt)):
                num_alt_calls += 1
                ref_allele = vcf_record.alleles[0]
                ref_bases_in_alt_calls += len(ref_allele)
                alt_allele = [alt for alt in sample.alleles if alt != ref_allele][0]
                edit_distance_induced_ref += edalign(ref_allele, alt_allele)[
                    "editDistance"
                ]
        record["num_alt_calls"] = num_alt_calls
        region_length = reg_end - reg_start + 1
        record["frac_ref_bases_in_alt_calls"] = ref_bases_in_alt_calls / region_length
        record["induced_ref_scaled_eddist"] = edit_distance_induced_ref / region_length

        print(record, file=output_file)
    output_file.close()
