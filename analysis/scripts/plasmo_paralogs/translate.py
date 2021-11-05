from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import click

import sys

@click.command()
@click.option('--fstops', help=
        """file location for separately outputting
           translations with premature stop codons""", default=None)
@click.argument('fin')
@click.argument('fout')
def main(fstops,fin, fout):
    with open(fin) as fhandle_in:
        records = SeqIO.parse(fin,"fasta")
        unique = set()
        used_translations = list()
        separate_stops = list()
        total = 0
        num_stops = 0

        for record in records:
            total += 1
            translation = record.seq.translate()
            str_translation = str(translation)
            if str_translation.count("*") > 1:
                num_stops += 1
                if fstops is not None:
                    separate_stops.append(SeqRecord(translation,id=record.id,description=""))
            elif str_translation not in unique:
                unique.add(str_translation)
                used_translations.append(SeqRecord(translation,id=record.id,description=""))

    print(f"Made {len(unique)} translations from {total} input records")
    print(f"Found {num_stops} translations with >1 stop codons")

    with open(fout,"w") as fhandle_out:
        SeqIO.write(used_translations, fhandle_out, "fasta")
    if fstops is not None:
        with open(fstops,"w") as fhandle_out:
            SeqIO.write(separate_stops, fhandle_out, "fasta")

if __name__ == "__main__":
    main()
