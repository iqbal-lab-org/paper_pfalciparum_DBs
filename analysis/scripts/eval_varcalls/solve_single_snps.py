from collections import Counter, defaultdict
from pathlib import Path

import click
import pandas as pd
from pybedtools import BedTool

from common_utils.pileup_utils import get_mpileups, pileup_matches, pileup_mismatches


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


def get_majority_base(mpileup_string):
    result = Counter(mpileup_string.upper()).most_common()[0][0]
    return result


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
    df.sort_values("sample", inplace=True)
    prev_sname = None
    for row in df.iterrows():
        sname = row[1]["sample"]
        gene_name = row[1]["gene"]
        if sname != prev_sname:
            ir_files = get_induced_ref_file_names(induced_ref_dir, sname)
        gene_region = get_gene_region_from_bed(ir_files["bed"], gene_name)
        mpileups, depths, ref_seq = get_mpileups(
            ir_files["bam"], ir_files["fasta"], gene_region, output_ref_seq=True
        )
        assert all(map(lambda d: d >= 5, depths))
        total_snp_positions = 0
        result_seq = ""
        for base, mpileup in zip(ref_seq, mpileups):
            majority_base = get_majority_base(mpileup)
            if pileup_mismatches(majority_base):
                total_snp_positions += 1
                result_seq += majority_base
            else:
                assert pileup_matches(majority_base)
                result_seq += base
        assert total_snp_positions == 1
        new_gene_seqs[gene_name].append(
                f">{gene_name}_{sname}\n"
                f"{result_seq}\n"
        )
    for gene_name, seqs in new_gene_seqs.items():
        ofname = Path(odir / f"{gene_name}_single_SNP_solved.fasta").open("w")
        print(*seqs, file=ofname)

if __name__ == "__main__":
    main()
