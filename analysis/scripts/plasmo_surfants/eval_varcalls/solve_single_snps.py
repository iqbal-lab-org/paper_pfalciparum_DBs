from collections import Counter, defaultdict
from pathlib import Path

import click
import pandas as pd
from pybedtools import BedTool

from plasmo_surfants.common_utils.pileup_utils import get_mpileups, pileup_matches, pileup_mismatches


def get_induced_ref_file_names(induced_ref_dir, sname):
    bed = Path(induced_ref_dir) / f"{sname}_ir_stats.translated_bed.txt"
    bam = Path(induced_ref_dir) / sname / "mapped.bam"
    fasta = Path(induced_ref_dir) / sname / "induced_ref.fa.gz"
    result = {"bed": bed, "bam": bam, "fasta": fasta}
    for path in result.values():
        assert path.exists()
    return result


def get_gene_region_from_bed(bed_fname, gene_name) -> str:
    bed = BedTool(bed_fname)
    interval = [i for i in bed if i.name == gene_name]
    assert len(interval) == 1
    interval = interval[0]
    return f"{interval.chrom}:{interval.start+ 1}-{interval.end}"


def get_most_freq_base_and_its_frequency(mpileup_string):
    counter = Counter(mpileup_string.upper().replace(",","."))
    most_common = counter.most_common()[0]
    most_freq_base = most_common[0]
    freq = most_common[1] / sum(counter.values())
    return most_freq_base, freq


@click.command()
@click.argument("snp_solvable_tsv")
@click.argument("induced_ref_dir")
@click.argument("output_dir")
def main(snp_solvable_tsv, induced_ref_dir, output_dir):
    """
    Given a tsv containing the sample and gene names of sequences with a single SNP
    in induced ref pileup, solves the SNP, and outputs the corresponding gene sequence.

    This complements producing gene seqs with zero such differences
    """
    odir = Path(output_dir)
    odir.mkdir(exist_ok=True, parents=True)
    new_gene_seqs = defaultdict(list)

    df = pd.read_csv(snp_solvable_tsv, sep="\t", index_col=None)
    for row in df.iterrows():
        sname = row[1]["sample"]
        gene_name = row[1]["gene"]
        ir_files = get_induced_ref_file_names(induced_ref_dir, sname)
        gene_region = get_gene_region_from_bed(ir_files["bed"], gene_name)
        mpileups, depths, ref_seq = get_mpileups(
            ir_files["bam"], ir_files["fasta"], gene_region, output_ref_seq=True
        )
        assert all(map(lambda d: d >= 5, depths))
        total_non_ref_positions = 0
        result_seq = ""
        solvable = True
        for ref_base, mpileup in zip(ref_seq, mpileups):
            most_freq_base, freq = get_most_freq_base_and_its_frequency(mpileup)
            if pileup_mismatches(most_freq_base):
                total_non_ref_positions += 1
                if freq > 0.5:
                    result_seq += most_freq_base
                else:
                    solvable = False
                    break
            else:
                assert pileup_matches(most_freq_base)
                assert freq >= 0.5
                result_seq += ref_base
        if total_non_ref_positions > 1:
            print(
            f"Sample {sname}, gene {gene_name} has {total_non_ref_positions}"
            " non-reference positions; sequence not outputted."
            )
        elif not solvable:
            print(
            f"The SNP in sample {sname}, gene {gene_name} was not at freq > 0.5"
            " ;sequence not outputted."
            )
        else:
            new_gene_seqs[gene_name].append(
                    f">{gene_name}_{sname}\n"
                    f"{result_seq}\n"
            )
    for gene_name, seqs in new_gene_seqs.items():
        ofname = Path(odir / f"{gene_name}_single_SNP_solved.fasta").open("w")
        print(*seqs, sep="",file=ofname)

if __name__ == "__main__":
    main()
