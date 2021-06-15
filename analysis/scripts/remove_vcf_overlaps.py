import sys
from pysam import VariantFile

def usage():
    print("purpose: Filters out variant records overlapping with previous record")
    print(f"usage: {sys.argv[0]} input_vcf output_vcf")
    exit(1)

def reset():
    return "", 0, ""

if __name__ == "__main__":
    if len(sys.argv) != 3:
        usage()
    input_vcf_fname = sys.argv[1]
    output_vcf_fname = sys.argv[2]

    input_vcf = VariantFile(input_vcf_fname)
    output_vcf = VariantFile(output_vcf_fname, "w", header = input_vcf.header)
    prev_chrom, prev_pos, prev_ref = reset()
    for record in input_vcf.fetch():
        if record.chrom != prev_chrom:
            prev_chrom, prev_pos, prev_ref = reset()
        if prev_pos + len(prev_ref) <= record.pos:
            output_vcf.write(record)
        prev_chrom = record.chrom
        prev_pos = record.pos
        prev_ref = record.alleles[0]
    output_vcf.close()

