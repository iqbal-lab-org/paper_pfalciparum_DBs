from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import click

import sys
import re

@click.command()
@click.option('--fstops', help=
        """file location for separately outputting
           translations with premature stop codons""", default=None)
@click.option('--aa_regexp',help=
        """AA pattern regexp to restrict the output translation to. this is useful for multi-exon genes, to look at only exon(s).""", default=None)
@click.argument('fin')
@click.argument('fout')
def main(fstops,aa_regexp,fin, fout):
    with open(fin) as fhandle_in:
        records = SeqIO.parse(fin,"fasta")
        used_translations = list()
        separate_stops = list()
        total = 0
        num_stops = 0

        for record in records:
            total += 1
            translation = record.seq.translate()
            str_translation = str(translation)
            if aa_regexp is not None:
                match = re.match(aa_regexp, str_translation)
                if match is not None:
                    str_translation = str_translation[match.start():match.end()]
            if str_translation.count("*") > 1:
                num_stops += 1
                if fstops is not None:
                    separate_stops.append(SeqRecord(Seq(str_translation),id=record.id,description=""))
            else:
                used_translations.append(SeqRecord(Seq(str_translation),id=record.id,description=""))

    print(f"Found {num_stops} translations with >1 stop codons")

    with open(fout,"w") as fhandle_out:
        SeqIO.write(used_translations, fhandle_out, "fasta")
    if fstops is not None:
        with open(fstops,"w") as fhandle_out:
            SeqIO.write(separate_stops, fhandle_out, "fasta")

if __name__ == "__main__":
    main()
