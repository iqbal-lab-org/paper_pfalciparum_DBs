"""
Input: a dir containing Tandem Repeat Finder (trf) outputs
Output: a tsv with tandem repeat lengths for each input fasta
"""
from pathlib import Path
import re
import subprocess

import click
from pysam import AlignmentFile
from pyfaidx import Fasta

from plasmo_surfants.common_utils.metrics import MetricsRecorder


class RepeatStats(MetricsRecorder):
    _headers = [
        "gene",
        "sample",
        "tool",
        "repeat_consensus",
        "length",
        "copy_number",
        "total_size",
    ]


trf_seqname_re = re.compile("Sequence: (\w+)@(\w+)@(\w+)")
trf_copynumber_re = re.compile(r"Copynumber: (\d+\.\d+)")
trf_repeat_string_re = re.compile(r"Consensus pattern \(.*\):\s+([A-Z]+)")


def extract_seqs(ofname, sam_dir, ref_dir):
    sam_fnames = Path(sam_dir).glob("*.sam")
    output_multifasta = ofname.open("w")
    for sam_fname in sam_fnames:
        assemb_name = sam_fname.stem
        assemb_fname = list(Path(ref_dir).glob(f"{assemb_name}*.fasta.gz"))
        assert len(assemb_fname) == 1
        assemb_index = Fasta(str(assemb_fname[0]))
        samfile = AlignmentFile(sam_fname, "r")
        for read in samfile.fetch():
            chromo = read.reference_name
            gene_name = read.query_name
            if chromo is not None:
                truth_seq = assemb_index[chromo][
                    read.reference_start : read.reference_end
                ]
                truth_seq = str(truth_seq).upper()
                output_multifasta.write(f">{gene_name}@{assemb_name}@truth\n{truth_seq}\n")
            output_multifasta.write(
                f">{gene_name}@{assemb_name}@pipeline_inference\n{read.query_sequence}\n"
            )


@click.command()
@click.argument(
    "input_sam_dir",
    type=click.Path(exists=True),
)
@click.argument(
    "input_assembly_dir",
    type=click.Path(exists=True),
)
@click.argument("output_dir", type=str)
@click.option(
    "--skip_seqs",
    is_flag=True,
    help="Set this flag to avoid sequence extraction from aligned sams + truth assembs",
)
@click.option("--skip_trf", is_flag=True, help="Set this flag to avoid running TRF")
def main(input_sam_dir, input_assembly_dir, output_dir, skip_seqs, skip_trf):
    """
    @param input_sam_dir: directory containing tool sequence inference aligned to truth
    assemblies
    @param input_assembly_dir: directory containing fastas of truth assemblies, for gene
    sequence extraction.
    """
    output_multifasta = Path(output_dir) / "sequences.fasta"
    output_multifasta.parent.mkdir(exist_ok=True, parents=True)
    if not skip_seqs:
        extract_seqs(output_multifasta, input_sam_dir, input_assembly_dir)
    if not skip_trf:
        subprocess.run(
            f"trf409.linux64 {output_multifasta.name} 2 5 7 80 10 50 2000",
            cwd=output_dir,
            shell=True,
        )
    trf_text_files = Path(output_dir).glob("*txt.html")
    result = list()
    for ifname in trf_text_files:
        text = ifname.open().read()
        m = trf_seqname_re.search(text)
        gene_name = m.groups()[0]
        sample_name = m.groups()[1]
        tool_name = m.groups()[2]
        copy_numbers = trf_copynumber_re.findall(text)
        repeat_strings = trf_repeat_string_re.findall(text)
        for copy_number, repeat_string in zip(copy_numbers, repeat_strings):
            result.append(
                RepeatStats(
                    gene=gene_name,
                    sample=sample_name,
                    tool=tool_name,
                    repeat_consensus=repeat_string,
                    copy_number=copy_number,
                    length=len(repeat_string),
                    total_size=float(copy_number) * len(repeat_string),
                )
            )
    output_tsv = Path(output_dir) / "trf_table.tsv"
    with output_tsv.open("w") as ofile:
        ofile.write(RepeatStats.get_header() + "\n")
        for res in result:
            ofile.write(str(res) + "\n")


if __name__ == "__main__":
    main()
