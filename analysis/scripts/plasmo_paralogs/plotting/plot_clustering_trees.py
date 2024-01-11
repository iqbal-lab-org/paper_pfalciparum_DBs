"""
Produces a (Newick-formatted) hierarchical clustering tree
And annotation metadata in ITOL-readable format
For visualisation on: https://itol.embl.de/

A note on laverania:
    - All of the other `seq_stats` analyses use MSAs computed on P. falciparum only, not other laverania species. This is normal: sharing, sequence logos, conversion clusters are all defined with respect to the falciparum gene pool
    - We also want to add laverania sequences to the analysis, here, and thus expect an MSA with laverania sequences.
    - For plotting the metadata though (especially the fraction of shared kmers in each sequence), we load the original falciparum-only MSA and compute the stats on that.
"""
import sqlite3
from itertools import tee
from typing import List, Tuple, Dict
from pathlib import Path

from Bio.Align import MultipleSeqAlignment as MSA
from Bio import AlignIO
from scipy.cluster.hierarchy import linkage, to_tree
from scipy.spatial.distance import pdist
import click

from plasmo_paralogs.common_utils import VALID_SEQTYPES, ID_DELIM
from plasmo_paralogs.common_utils import (
    REGION_CLICK_HELP,
    REGION_DELIM,
    VALID_SQLITE_DATA,
)
from plasmo_paralogs.common_utils.msa import (
    split_alignment_by_gene,
    get_gene_name,
    get_sample_id,
)
from plasmo_paralogs.common_utils.sqlite import (
    get_sqlite_table_name,
    load_percent_identities,
)
from plasmo_paralogs.common_utils.metrics import MetricsRecorder
from plasmo_paralogs.common_utils.mosaic import Donors, parse_mosaic_alignment_file
from plasmo_paralogs.seq_stats.sharing import load_kmer_sharing_mapping
from plasmo_paralogs.common_utils.sequences import get_pos_tuple, extract_kmer


LAV_SPECIES = ["P. praefalciparum", "P. reichenowi", "P. billcollinsi", "P. falciparum", "P. adleri", "P. gaboni", "P. blacklocki"]


def get_species_name(record):
    species_map = {
        "PPRFG01": LAV_SPECIES[0],
        "PPRFG02": LAV_SPECIES[0],
        "PPRFG03": LAV_SPECIES[0],
        "PPRFG04": LAV_SPECIES[0],
        "PRG01": LAV_SPECIES[1],
        "PREICH001": LAV_SPECIES[1],
        "PBILCG01": LAV_SPECIES[2],
        "PADLG01": LAV_SPECIES[4],
        "PBLACG01": LAV_SPECIES[6],
        "PGABG01": LAV_SPECIES[5],
        "PGABG02": LAV_SPECIES[5],
    }
    sample_name = record.id.split(ID_DELIM)[-1].strip()
    return species_map.get(sample_name, LAV_SPECIES[3])


def load_MSA(msa_fname: str, start: int, end: int):
    with open(msa_fname) as fhandle_in:
        alignment = AlignIO.read(fhandle_in, "fasta")
        alignment = alignment[:, start - 1 : end]
        for record in alignment:
            record.id = record.id.split("::")[0]
    return alignment


def get_split_alignments(alignment, seq_names_to_include, SHARED_ID, dedup=True):
    if dedup:
        alignment = deduplicate(alignment, seq_names_to_include=seq_names_to_include)
    split_alignments = split_alignment_by_gene(alignment)
    split_alignments[SHARED_ID] = alignment
    return split_alignments


def pairwise(iterable):
    # pairwise('ABCDEFG') --> AB BC CD DE EF FG
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


SeqNames = List[str]
Newick = str
PROTEIN_CHARS = {
    "A",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "K",
    "L",
    "M",
    "N",
    "P",
    "Q",
    "R",
    "S",
    "T",
    "V",
    "W",
    "Y",
    "-",
}
PROT_TO_INT = {char: i for i, char in enumerate(PROTEIN_CHARS)}

convert_prot_to_int = lambda prot_seq: [PROT_TO_INT[char] for char in prot_seq]


class FracShared(MetricsRecorder):
    _headers = ["sample_ID", "frac_shared"]
    _required_headers = _headers


class NodeLinks(MetricsRecorder):
    _headers = ["sample_ID_1", "sample_ID_2", "width", "color", "style", "label"]
    _required_headers = _headers


def deduplicate(msa: MSA, seq_names_to_include: SeqNames = []) -> MSA:
    unique_seqs = set()
    seq_names_to_include = set(seq_names_to_include)
    result = MSA([])
    for record in msa:
        seq = str(record.seq)
        if record.id in seq_names_to_include and seq not in unique_seqs:
            unique_seqs.add(seq)
            result.append(record)
    for record in msa:
        seq = str(record.seq)
        if seq not in unique_seqs:
            unique_seqs.add(seq)
            result.append(record)
    return result


def tree2newick(node, newick, parentdist, leaf_names):
    """
    Taken from: https://gist.github.com/fbkarsdorp/a5d95f7ccc8c308a70c93aff38ac4860
    [Note] For the corresponding feature request in scipy, see https://github.com/scipy/scipy/issues/8274
    """
    if node.is_leaf():
        return f"{leaf_names[node.id]}:{parentdist - node.dist:.2f}{newick}"
    else:
        if len(newick) > 0:
            newick = f"):{parentdist - node.dist:.2f}{newick}"
        else:
            newick = ");"
        newick = tree2newick(node.get_left(), newick, node.dist, leaf_names)
        newick = tree2newick(node.get_right(), f",{newick}", node.dist, leaf_names)
        newick = f"({newick}"
        return newick


def newick_from_msa(msa: MSA) -> Newick:
    leaf_names = []
    int_seqs = []
    for record in msa:
        int_seq = convert_prot_to_int(str(record.seq))
        int_seqs.append(int_seq)
        leaf_names.append(record.id)
    distances = pdist(int_seqs, metric="hamming")
    linkages = linkage(distances, method="average")
    tree = to_tree(linkages)
    result = tree2newick(tree, "", tree.dist, leaf_names)
    return result


def compute_and_write_frac_shared(
    msa: MSA, kmer_sharing_mapping, kmer_size, region_start, gene_dir
) -> None:
    result = list()
    for record in msa:
        sharing_array = []
        seq = str(record.seq)
        for start, center_pos, stop in get_pos_tuple(kmer_size, [seq], region=None):
            kmer_seq = extract_kmer(seq, start, stop)
            kmer_mapping_key = f"{center_pos+region_start-1}{ID_DELIM}{kmer_seq}"
            sharing_array.append(int(kmer_sharing_mapping[kmer_mapping_key]))
        result.append(
            FracShared(
                sample_ID=record.id, frac_shared=sum(sharing_array) / len(sharing_array)
            )
        )
    with (gene_dir / f"frac_shared.txt").open("w") as fout:
        fout.write("DATASET_GRADIENT\n")
        fout.write("SEPARATOR TAB\n")
        fout.write(f"DATASET_LABEL\tFraction of shared {kmer_size}-mer peptides\t1\n")
        fout.write("COLOR\t#ff0000\n")
        fout.write("MARGIN\t10\n")
        fout.write("BORDER_WIDTH\t4\n")
        fout.write("BORDER_COLOR\t#000000\n")
        fout.write("COLOR_MIN\t#e02582\n")
        # fout.write("COLOR_MIN\t#cb0b0b\n")
        fout.write("COLOR_MAX\t#14ee14\n")
        fout.write("DATA\n")
        for entry in result:
            fout.write(str(entry))


def compute_and_write_recomb_links(donors: Donors, outdir, gene_name=None):
    result = list()
    for target_ID, donor_list in donors.items():
        for e1, e2 in pairwise(donor_list):
            result.append(
                NodeLinks(
                    sample_ID_1=e1.seq_ID,
                    sample_ID_2=e2.seq_ID,
                    width=4,
                    color="#8C99A6",
                    style="dashed",
                    label=target_ID,
                )
            )
    with (outdir / f"recomb_links.txt").open("w") as fout:
        fout.write("DATASET_CONNECTION\n")
        fout.write("SEPARATOR TAB\n")
        fout.write(f"DATASET_LABEL\tSample pairs that recombined\n")
        fout.write("COLOR\t#ff0000\n")
        fout.write("CENTER_CURVES\t1\n")
        fout.write("ALIGN_TO_LABELS\t1\n")
        fout.write("DRAW_ARROWS\t0\n")
        fout.write("MAXIMUM_LINE_WIDTH\t4\n")
        fout.write("DATA\n")
        for entry in result:
            fout.write(str(entry))


def write_gene_colours(msa: MSA, gene_dir):
    colour_map = {"DBLMSP": "#efbd24", "DBLMSP2": "#3d85c6"}
    all_colours = "\t".join(colour_map.values())
    all_labels = "\t".join(colour_map.keys())
    with (gene_dir / f"gene_colours.txt").open("w") as fout:
        fout.write("DATASET_COLORSTRIP\n")
        fout.write("SEPARATOR TAB\n")
        fout.write(f"DATASET_LABEL\tGene\n")
        fout.write("COLOR\t#ff0000\n")
        fout.write("BORDER_WIDTH\t4\n")
        fout.write("BORDER_COLOR\t#000000\n")
        fout.write(f"LEGEND_TITLE\tGene\n")
        fout.write(f"LEGEND_COLORS\t{all_colours}\n")
        fout.write(f"LEGEND_LABELS\t{all_labels}\n")
        fout.write(f"LEGEND_SHAPES\t1\t1\n")
        fout.write("DATA\n")
        for record in msa:
            gene_name = get_gene_name(record)
            colour = colour_map.get(gene_name, "#efbd24")
            fout.write(f"{record.id}\t{colour}\t{gene_name}\n")


def write_species_colours(msa: MSA, gene_dir):
    colour_map = {
        LAV_SPECIES[0]: "#0008e8",
        LAV_SPECIES[1]: "#0096ff",
        LAV_SPECIES[2]: "#00e8c4",
        LAV_SPECIES[3]: "#8a03ff",
        LAV_SPECIES[4]: "#d52d2d",
        LAV_SPECIES[5]: "#ffa246",
        LAV_SPECIES[6]: "#46ff74",
    }
    all_colours = "\t".join(colour_map.values())
    all_labels = "\t".join(colour_map.keys())
    with (gene_dir / f"species_colours.txt").open("w") as fout:
        fout.write("DATASET_COLORSTRIP\n")
        fout.write("SEPARATOR TAB\n")
        fout.write(f"DATASET_LABEL\tSpecies\n")
        fout.write("COLOR\t#ff0000\n")
        fout.write("MARGIN\t10\n")
        fout.write("BORDER_WIDTH\t4\n")
        fout.write("BORDER_COLOR\t#000000\n")
        fout.write(f"LEGEND_TITLE\tSpecies\n")
        fout.write(f"LEGEND_COLORS\t{all_colours}\n")
        fout.write(f"LEGEND_LABELS\t{all_labels}\n")
        fout.write(f"LEGEND_SHAPES"+"\t1"*len(colour_map.keys())+"\n")
        fout.write("DATA\n")
        for record in msa:
            species_name = get_species_name(record)
            colour = colour_map[species_name]
            fout.write(f"{record.id}\t{colour}\t{species_name}\n")


def write_premature_stops(msa: MSA, gene_dir):
    colour_map = {True: "#2B0000", False: "#FFCCAA"}
    all_colours = "\t".join(colour_map.values())
    all_labels = "\t".join(map(str,colour_map.keys()))
    with (gene_dir / f"premature_stops.txt").open("w") as fout:
        fout.write("DATASET_COLORSTRIP\n")
        fout.write("SEPARATOR TAB\n")
        fout.write(f"DATASET_LABEL\tPremature_stops\n")
        fout.write("COLOR\t#ff0000\n")
        fout.write("MARGIN\t10\n")
        fout.write("BORDER_WIDTH\t4\n")
        fout.write("BORDER_COLOR\t#000000\n")
        fout.write(f"LEGEND_TITLE\tProtein has premature stop codon(s)\n")
        fout.write(f"LEGEND_COLORS\t{all_colours}\n")
        fout.write(f"LEGEND_LABELS\t{all_labels}\n")
        fout.write(f"LEGEND_SHAPES\t1\t1\n")
        fout.write("DATA\n")
        for record in msa:
            premature_stop = "STOP" in record.id
            colour = colour_map[premature_stop]
            fout.write(f"{record.id}\t{colour}\t{premature_stop}\n")


def write_percent_identities(msa: MSA, percent_identities, gene_dir) -> None:
    with (gene_dir / f"percent_identities.txt").open("w") as fout:
        fout.write("DATASET_GRADIENT\n")
        fout.write("SEPARATOR TAB\n")
        fout.write(
            f"DATASET_LABEL\tFraction of identical codons between paralogs in sample\t1\n"
        )
        fout.write("COLOR\t#ff0000\n")
        fout.write("MARGIN\t10\n")
        fout.write("COLOR_MIN\t#0000ff\n")
        fout.write("COLOR_MAX\t#ff0000\n")
        fout.write("DATA\n")
        for record in msa:
            sample_ID = get_sample_id(record)
            if sample_ID in percent_identities:
                fout.write(f"{record.id}\t{percent_identities[sample_ID]}\n")


@click.command()
@click.argument("msa_fname")
@click.argument("dir_ofname")
@click.option("-s", "--sqlite_db_fpath", required=True)
@click.option("--seqtype", required=True, type=click.Choice(VALID_SEQTYPES))
@click.option("--seq_region", type=str, required=True, help=REGION_CLICK_HELP)
@click.option(
    "--kmer_sizes",
    type=str,
    help="kmer sizes to compute; comma-delimited",
    default="10",
    show_default=True,
)
@click.option("-m", "--mosaic_alignment_fname")
@click.option("--just_tree", is_flag=True,)
def main(
    msa_fname,
    dir_ofname,
    sqlite_db_fpath,
    seqtype,
    seq_region,
    kmer_sizes,
    mosaic_alignment_fname,
    just_tree
):
    SHARED_ID = "DBs"
    start, end = map(int, seq_region.split(REGION_DELIM))
    outdir = Path(dir_ofname)
    outdir.mkdir(exist_ok=True, parents=True)

    # Load MSA, with laverania
    alignment = load_MSA(msa_fname, start, end)

    if not just_tree:
        # Load MSA, without laverania
        msa_fname_no_laverania = msa_fname.replace("_and_laverania", "").replace("_and_premature_stops","")
        alignment_no_laverania = load_MSA(msa_fname_no_laverania, start, end)

        # Load and compute mosaic-recomb links
        if mosaic_alignment_fname is not None:
            recomb_donors = parse_mosaic_alignment_file({}, mosaic_alignment_fname, True)
            compute_and_write_recomb_links(recomb_donors, outdir)

    # Load percent identities
    if not just_tree:
        dna_msa_fname = msa_fname_no_laverania.replace(VALID_SEQTYPES[0], VALID_SEQTYPES[1])
        table_name = get_sqlite_table_name(
            dna_msa_fname,
            SHARED_ID,
            VALID_SEQTYPES[1],
            metric_info=VALID_SQLITE_DATA[5],
            seq_region=seq_region,
        )
        percent_identities = load_percent_identities(sqlite_db_fpath, table_name, SHARED_ID)
        percent_identities = {
            key: val for key, val in percent_identities.items() if val > 0.5
        }
        high_shared_samples = []
        for sample_ID in percent_identities:
            for gene_name in ["DBLMSP", "DBLMSP2"]:
                high_shared_samples.append(f"{gene_name}{ID_DELIM}{sample_ID}")

        seq_names_to_include = list(recomb_donors.keys()) + high_shared_samples
    else:
        seq_names_to_include = list()
    split_alignments = get_split_alignments(alignment, seq_names_to_include, SHARED_ID)
    # Write trees
    region_for_output = f"{start}{ID_DELIM}{end}"
    for gene_name, gene_alignment in split_alignments.items():
        gene_dir = Path(outdir) / f"{gene_name}{ID_DELIM}{region_for_output}"
        gene_dir.mkdir(exist_ok=True)
        tree = newick_from_msa(gene_alignment)
        with (gene_dir / "clustering_tree.nwk").open("w") as fout:
            fout.write(tree)

    assert len(kmer_sizes.split(",")) == 1

    ## Write additional data for trees
    for gene_name, gene_alignment in split_alignments.items():
        gene_dir = Path(outdir) / f"{gene_name}{ID_DELIM}{region_for_output}"
        if not just_tree:
            write_percent_identities(gene_alignment, percent_identities, gene_dir)
        if gene_name == SHARED_ID:
            write_species_colours(gene_alignment, gene_dir)
            write_gene_colours(gene_alignment, gene_dir)
            write_premature_stops(gene_alignment, gene_dir)

    if not just_tree:
        split_alignments_no_laverania = get_split_alignments(
            alignment_no_laverania, seq_names_to_include, SHARED_ID, dedup=False
        )
        ## Load kmer sharing stats
        con = sqlite3.connect(sqlite_db_fpath)
        input_table_name = get_sqlite_table_name(
            msa_fname_no_laverania, SHARED_ID, seqtype, metric_info=VALID_SQLITE_DATA[2]
        )
        kmer_sharing_mapping = load_kmer_sharing_mapping(input_table_name, con)
        con.close()
        for gene_name, gene_alignment in split_alignments_no_laverania.items():
            compute_and_write_frac_shared(
                gene_alignment, kmer_sharing_mapping, int(kmer_sizes), start, gene_dir
            )


if __name__ == "__main__":
    main()
