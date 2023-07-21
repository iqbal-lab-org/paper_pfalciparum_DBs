from typing import Tuple, Optional

from plasmo_paralogs.common_utils import REGION_DELIM
from plasmo_paralogs.common_utils.msa import MSA, ensure_valid_alignment_slice
from plasmo_paralogs.seq_stats.sharing import UniqueKmer, UniqueKmers

def extract_kmer(seq: str, start: int, stop: int) -> str:
    return seq[start:stop + 1]

def extract_kmers(alignment: MSA, start: int, stop: int, gene_id: str) -> UniqueKmers:
    """
    Generator producing all unique strings in an alignment slice
    """
    ensure_valid_alignment_slice(alignment, start, stop)
    for record in alignment:
        kmer = extract_kmer(str(record.seq), start, stop)
        sample_id = record.id
        new_kmer = UniqueKmer(kmer, feature_classes={gene_id}, instances={sample_id})
        yield new_kmer

def _compute_gap_sizes(kmer_size: int) -> Tuple[int, int]:
    pre_gap_size = int(kmer_size / 2)
    post_gap_size = pre_gap_size
    if kmer_size % 2 == 0:
        post_gap_size -= 1
    return pre_gap_size, post_gap_size

def _bounded_kmer_positions(start: int, stop: int, kmer_size: int) -> Tuple[int, int, int]:
    pre_gap_size, post_gap_size = _compute_gap_sizes(kmer_size)
    for center_position in range(start, stop):
        yield (center_position - pre_gap_size, center_position, center_position + post_gap_size)

def _all_kmer_positions(kmer_size: int, seqlen: int) -> Tuple[int,int,int]:
    """
    Returns all (start, center_pos, stop) tuples of kmer positions of size `kmer_size`
    """
    pre_gap_size, post_gap_size = _compute_gap_sizes(kmer_size)
    for center_position in range(pre_gap_size, seqlen - post_gap_size):
        yield (center_position - pre_gap_size, center_position, center_position + post_gap_size)

def get_pos_tuple(kmer_size, alignment, region:Optional[str]) -> Tuple[int, int, int]:
    if region is not None:
        region_spec = region.split(REGION_DELIM)
        pos_tuple = _bounded_kmer_positions(int(region_spec[0]), int(region_spec[1]), kmer_size)
    else:
        pos_tuple = _all_kmer_positions(kmer_size, len(alignment[0]))
    return pos_tuple
