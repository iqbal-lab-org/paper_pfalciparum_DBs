import re

from pysam import mpileup as sam_mpileup

PILEUP_INDEL_RE = "[+-]([0-9]+)[ACGTNacgtn*#]+"
PILEUP_MATCH = ".,"
PILEUP_MISMATCH = "ATCGNatcgn*"

pileup_matches = lambda char: char in PILEUP_MATCH
pileup_mismatches = lambda char: char in PILEUP_MISMATCH


def remove_read_ends(pileup_string: str) -> str:
    new_string = ""
    cur_pos = 0
    for match in re.finditer(r"(\^.)|(\$)", pileup_string):
        new_string += pileup_string[cur_pos : match.start(0)]
        cur_pos = match.end(0)
    new_string += pileup_string[cur_pos:]
    return new_string


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


def get_mpileups(input_bam, input_ref_genome, region, output_ref_seq=False):
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
    if output_ref_seq:
        ref_seq = "".join(map(lambda line: line.split("\t")[2], mpileup_output))
        return pileups, depths, ref_seq
    else:
        return pileups, depths, None
