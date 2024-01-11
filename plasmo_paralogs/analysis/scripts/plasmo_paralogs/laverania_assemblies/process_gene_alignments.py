import edlib
import click
from Bio import SeqIO
import pandas as pd

from plasmo_paralogs.common_utils.metrics import MetricsRecorder

FIELD_DELIM = "_"

class AlignmentMetrics(MetricsRecorder):
    _headers = [
        "gene_ID",
        "assembly_ID",
        "assembly_species",
        "hit_location",
        "hit_length",
        "hit_length_as_frac_of_3d7_length",
        "hit_percent_identity",
    ]
    _required_headers = _headers


def load_sequences(seq_ifname):
    result = {}
    with open(seq_ifname) as fhandle_in:
        for record in SeqIO.parse(fhandle_in, "fasta"):
            result[record.id.strip()] = str(record.seq)
    return result

def process_one_hit_sequence(hit_ID, hit_sequence, ref_sequences, species_names):
    hit_elems = hit_ID.split("::")
    hit_location = hit_elems[1]
    ID_fields = hit_elems[0].split(FIELD_DELIM)
    gene_ID = ID_fields[0]
    assembly_ID = ID_fields[2]
    species_name = species_names[assembly_ID]
    ref_sequence = ref_sequences[f'{gene_ID}_ref']
    # Mode 'HW' is used to discard 5' and 3' truncations of `ref_sequence` in `hit_sequence`,
    # that I don't want to include in the edit distance.
    alignment = edlib.align(hit_sequence, ref_sequence, mode="HW", task="path")
    hit_length = len(hit_sequence)
    hit_length_as_frac = round(hit_length / len(ref_sequence), 4)
    hit_perc_ID = round(alignment["editDistance"] / hit_length, 4)
    return AlignmentMetrics(
            gene_ID = gene_ID,
            assembly_ID = assembly_ID,
            assembly_species = species_name,
            hit_location = hit_location,
            hit_length = hit_length,
            hit_length_as_frac_of_3d7_length = hit_length_as_frac,
            hit_percent_identity = hit_perc_ID
            )


@click.command()
@click.argument("ref_seqs_ifname", type=click.Path(exists=True))
@click.argument("hit_seqs_ifname", type=click.Path(exists=True))
@click.argument("assembly_id_file", type=click.Path(exists=True))
@click.argument("tsv_ofname")
def main(ref_seqs_ifname, hit_seqs_ifname, assembly_id_file, tsv_ofname):
    ref_sequences = load_sequences(ref_seqs_ifname)
    hit_sequences = load_sequences(hit_seqs_ifname)
    laverania_names = pd.read_table(assembly_id_file, sep="\t", index_col=0)
    laverania_names = dict(zip(laverania_names["Abbreviation"], laverania_names["Species"]))
    with open(tsv_ofname, "w") as tsv_out:
        tsv_out.write(AlignmentMetrics.get_header())
        for hit_ID, hit_seq in hit_sequences.items():
            stats = process_one_hit_sequence(hit_ID, hit_seq, ref_sequences, laverania_names)
            tsv_out.write(str(stats))

if __name__ == "__main__":
    main()
