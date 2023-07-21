import sqlite3
from typing import List, Dict
from itertools import count
from pathlib import Path
from collections import defaultdict
from subprocess import run

import click
from Bio import SeqIO
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from plasmo_paralogs.common_utils import VALID_SEQTYPES, ID_DELIM
from plasmo_paralogs.common_utils.msa import get_gene_name, get_sample_id
from plasmo_paralogs.common_utils.metrics import MetricsRecorder
from plasmo_paralogs.common_utils.sqlite import get_sqlite_table_name
from plasmo_paralogs.common_utils.sequences import get_pos_tuple, extract_kmer
from plasmo_paralogs.common_utils.mosaic import parse_mosaic_alignment_file, DonorList
from plasmo_paralogs.common_utils import REGION_CLICK_HELP, REGION_DELIM, VALID_SQLITE_DATA
from plasmo_paralogs.plotting.svg import (
    SVG_HEADER,
    Coord,
    Grid,
    LineSpec,
    DEFAULT_RECT_HEIGHT,
    DEFAULT_RECT_WIDTH,
    DEFAULT_OPACITY,
    DEFAULT_TEXT_FILL,
)
from plasmo_paralogs.seq_stats.sharing import load_kmer_sharing_mapping


# Opacities
FOCAL_OP = 1
TEXT_MUTATION_OP = 0.8
NON_FOCAL_DONOR_OP = 0.5
BACKGROUND_OP = 0.2


def sort_records_by_sname(records):
    """
    Sorts records formatted as {gene_name}_{sample_name} first by {gene_name}, 
    then by {sample_name} - plus puts the samples called "ref" first
    """
    initial_len = len(records)
    snames = [key.split("_") for key in records if "ref" not in key]
    assert len(snames) == (initial_len - 2)
    sorted_snames = sorted(snames, key=lambda sname: [sname[0],sname[1]])
    prev_gene_name = None
    all_snames = list()
    for elem in sorted_snames:
        if elem[0] != prev_gene_name:
            all_snames.append("_".join([elem[0],"ref"]))
            prev_gene_name = elem[0]
        all_snames.append("_".join(elem))
    assert initial_len == len(all_snames)
    result = {}
    for sname in all_snames:
        result[sname] = records[sname]
    return result


def write_svg_grid(out_fname, svg_grid, height_margin=0):
    with open(out_fname, "w") as fout:
        width, height = svg_grid.size
        height += height_margin
        fout.write(SVG_HEADER)
        fout.write(f'<svg height="{height}" width="{width}">\n')
        fout.write(str(svg_grid))
        fout.write("</svg>")


def write_combined_grids(files_to_combine, dir_ofname, gene_name):
    snames = []
    for _fname in files_to_combine:
        snames.append(get_sample_id(Path(_fname).stem, already_ID=True))
    snames = "_".join(snames)
    combined_fname = f"{dir_ofname}/{gene_name}_{snames}.png"
    input_fnames = " ".join(files_to_combine)
    print(f"convert -append {input_fnames} {combined_fname}")
    run(f"convert -append {input_fnames} {combined_fname}", shell=True, check=True)


def write_combined_breakpoints(dir_ofname, all_vlines, ncol):
    # makes mini-grids showing all breakpoints
    for gene_name, vlines in all_vlines.items():
        for vline in vlines:
            vline.width = DEFAULT_RECT_WIDTH
        cumulative_svg_grid = Grid(
            Coord(0, 0), 1, ncol, vlines=vlines, rect_height=DEFAULT_RECT_HEIGHT * 3
        )
        for col_idx in range(ncol):
            cumulative_svg_grid[0][col_idx].add_properties(
                fill="#1f77b4", fill_opacity=FOCAL_OP
            )
        write_svg_grid(
            f"{dir_ofname}/all_breakpoints_{gene_name}.svg", cumulative_svg_grid
        )


def get_hlines(record_dict):
    result = list()
    previous = None
    for i, rec_id in enumerate(record_dict.keys()):
        gene_name = get_gene_name(rec_id, already_ID=True)
        if gene_name != previous:
            if previous != None:
                result.append(LineSpec(i - 1, "white", DEFAULT_RECT_HEIGHT * 1.25))
            previous = gene_name
    return result


def get_vlines(donor_list: DonorList):
    result = list()
    last_donor_idx = len(donor_list) - 1
    for i, donor in enumerate(donor_list):
        if i != last_donor_idx:
            result.append(LineSpec(donor.end, "red", DEFAULT_RECT_WIDTH / 2))
    return result


def load_sequence_colour_maps(con, seqtype):
    df_colours = pd.read_sql_query("SELECT * from seq_colour_schemes", con)
    df_colours = df_colours[df_colours["name"] == "seaview"]
    df_colours = df_colours[df_colours["seqtype"] == seqtype]
    colour_mapping = dict()
    for group, hex_colour in zip(df_colours["chars"], df_colours["rgb_hex"]):
        colour_mapping.update({char: hex_colour for char in group})
    return colour_mapping


def write_svg_sequence_grids(
    dir_ofname, colour_mapping, donors, record_dict, nrow, ncol, hlines
):
    Path(dir_ofname).mkdir(exist_ok=True, parents=True)

    made_svgs = {}
    all_vlines = defaultdict(list)
    best_mismatches = list()

    for target_ID, donor_list in donors:
        mismatch_dict = defaultdict(int)
        mismatch_dict[BEST_MATCH] = 0
        target_seq = record_dict[target_ID]
        vlines = get_vlines(donor_list)
        gene_name = get_gene_name(target_ID, already_ID=True)
        all_vlines[gene_name].extend(vlines)
        svg_grid = Grid(Coord(0, 0), nrow, ncol, vlines=vlines, hlines=hlines)
        for row_idx, rec_id, seq in zip(
            count(), record_dict.keys(), record_dict.values()
        ):
            for col_idx, char in enumerate(seq):
                text_fill = DEFAULT_TEXT_FILL
                text_specific_opacity = None
                if rec_id == target_ID:
                    opacity = FOCAL_OP
                elif donors.id_match(target_ID, rec_id):
                    opacity = NON_FOCAL_DONOR_OP
                    mismatch_found = char != target_seq[col_idx]
                    if donors.position_match(target_ID, rec_id, col_idx):
                        opacity = FOCAL_OP
                        if mismatch_found:
                            mismatch_dict[BEST_MATCH] += 1
                    if mismatch_found:
                        text_fill = "red"
                        text_specific_opacity = TEXT_MUTATION_OP
                        mismatch_dict[rec_id] += 1
                else:
                    opacity = BACKGROUND_OP
                svg_grid[row_idx][col_idx].add_properties(
                    fill=colour_mapping[char], fill_opacity=opacity
                )
                svg_grid[row_idx][col_idx].add_text(
                    char, opacity=opacity, fill=text_fill
                )
                if text_specific_opacity is not None:
                    svg_grid[row_idx][col_idx].add_text_properties(
                        opacity=text_specific_opacity
                    )

        no_recomb_vals = {
            mismatch_dict[key] for key in mismatch_dict if key != BEST_MATCH
        }
        best_mismatches.append(
            BestMismatches(
                target=target_ID,
                best_recomb=mismatch_dict[BEST_MATCH],
                best_no_recomb=min(no_recomb_vals),
                all_data=str(mismatch_dict),
            )
        )

        svg_fname = f"{dir_ofname}/{target_ID}.svg"
        write_svg_grid(svg_fname, svg_grid, height_margin=DEFAULT_RECT_HEIGHT * 2)
        made_svgs[target_ID] = svg_fname

    write_combined_breakpoints(dir_ofname, all_vlines, ncol)
    return made_svgs, best_mismatches


def write_svg_sharing_grids(
    dir_ofname,
    colour_mapping,
    kmer_sharing_mapping,
    donors,
    record_dict,
    nrow,
    ncol,
    hlines,
    region_start,
    kmer_sizes,
):
    Path(dir_ofname).mkdir(exist_ok=True, parents=True)
    made_svgs = {}
    for kmer_size in map(int, kmer_sizes.split(",")):
        for target_ID, donor_list in donors:
            target_seq = record_dict[target_ID]
            vlines = get_vlines(donor_list)
            rolling_shift = kmer_size // 2
            for vline in vlines:
                # Adjust recomb breakpoints because the grid is being depicted using
                # kmers centered (not starting at) each position
                vline.index -= rolling_shift
            svg_grid = Grid(Coord(0, 0), nrow, ncol - kmer_size + 1, vlines=vlines, hlines=hlines)
            for row_idx, rec_id, seq in zip(
                count(), record_dict.keys(), record_dict.values()
            ):
                for start, center_pos, stop in get_pos_tuple(
                    kmer_size, [seq], region=None
                ):
                    kmer_seq = extract_kmer(seq, start, stop)
                    kmer_mapping_key = (
                        f"{center_pos+region_start-1}{ID_DELIM}{kmer_seq}"
                    )
                    colour_key = get_gene_name(rec_id, already_ID=True)
                    if kmer_sharing_mapping[kmer_mapping_key]:
                        colour_key = "shared"
                    if rec_id == target_ID:
                        opacity = FOCAL_OP
                    elif donors.id_match(target_ID, rec_id):
                        opacity = NON_FOCAL_DONOR_OP
                        if donors.position_match(target_ID, rec_id, center_pos):
                            opacity = FOCAL_OP
                    else:
                        opacity = BACKGROUND_OP
                    svg_grid[row_idx][start].add_properties(
                        fill=colour_mapping[colour_key], fill_opacity=opacity
                    )

            svg_fname = f"{dir_ofname}/{target_ID}.svg"
            write_svg_grid(svg_fname, svg_grid, height_margin=DEFAULT_RECT_HEIGHT * 2)
            made_svgs[target_ID] = svg_fname
    return made_svgs

def write_best_mismatches(best_mismatches,dir_ofname):
    best_mismatches_ofname = f"{dir_ofname}/best_mismatches.tsv"
    with open(best_mismatches_ofname, "w") as fout:
        fout.write(BestMismatches.get_header())
        for r in best_mismatches:
            fout.write(str(r))

    df = pd.read_csv(best_mismatches_ofname,sep="\t")
    df["mosaic_edit_gain"] = df["best_no_recomb"] - df["best_recomb"]

    sns.set_theme()
    # Histogram
    fig = plt.figure()
    ax = fig.add_subplot()
    sns.histplot(df, x="mosaic_edit_gain",ax = ax,bins=15)
    ax.set_xlabel("Number of edits gained by \nallowing recombination")
    ax.set_ylabel("Number of aligned targets")
    fig.savefig(f"{dir_ofname}/best_mismatches_histo.pdf")

    # Dot plot
    fig = plt.figure()
    ax = fig.add_subplot()
    sns.scatterplot(df, y="best_recomb",x="best_no_recomb",ax = ax)
    ax.set_ylabel("Edit distance for mosaic alignment")
    ax.set_xlabel("Edit distance for best single-donor alignment")
    sns.lineplot(df, y="best_recomb",x="best_recomb",ax = ax,color="grey",linestyle="dotted")
    fig.savefig(f"{dir_ofname}/best_mismatches_scatter.pdf")

    df.to_csv(best_mismatches_ofname,sep="\t")

def to_pdf(made_svgs):
    dir_ofname = Path(list(made_svgs.values())[0]).parent
    dir_ofname = dir_ofname / "pdf"
    dir_ofname.mkdir(exist_ok=True)
    for target_ID, svg_fname in made_svgs.items():
        stem_name = Path(svg_fname).stem
        gene_name = get_gene_name(target_ID, already_ID=True)
        ofname = str(dir_ofname / f"{stem_name}.pdf")
        run(
            f'inkscape --without-gui --export-filename="{ofname}" {svg_fname}',
            shell=True,
            check=True,
        )
    with open(dir_ofname / f"file_order.txt","w") as fout:
        print(*made_svgs.values(), sep="\n", file=fout)

def rasterise_and_combine(made_svgs, seq_names_in_order):
    combined = []
    prev_gene_name = ""
    dir_ofname = Path(list(made_svgs.values())[0]).parent
    dir_ofname = dir_ofname / "rasterised"
    dir_ofname.mkdir(exist_ok=True)
    for target_ID in seq_names_in_order:
        svg_fname = made_svgs[target_ID]
        stem_name = Path(svg_fname).stem
        gene_name = get_gene_name(target_ID, already_ID=True)
        ofname = str(dir_ofname / f"{stem_name}.png")
        run(
            f'inkscape --without-gui --export-filename="{ofname}" --export-dpi 50 {svg_fname}',
            shell=True,
            check=True,
        )
        if len(combined) == 6 or (
            (gene_name != prev_gene_name) & (prev_gene_name != "")
        ):
            write_combined_grids(combined, dir_ofname, prev_gene_name)
            combined = []
        combined.append(ofname)
        prev_gene_name = gene_name
    write_combined_grids(combined, dir_ofname, prev_gene_name)


BEST_MATCH = "best_match"


class BestMismatches(MetricsRecorder):
    """
    Records number of mismatches for recombination-tolerant and -intolerant alignments
    """

    _headers = ["target", "best_recomb", "best_no_recomb", "all_data"]
    _required_headers = _headers


@click.command()
@click.argument("msa_fname")
@click.argument("mosaic_alignment_fname")
@click.argument("dir_ofname")
@click.option("-s", "--sqlite_db_fpath", required=True)
@click.option("--seqtype", required=True, type=click.Choice(VALID_SEQTYPES))
@click.option("--seq_region", type=str, required=True, help=REGION_CLICK_HELP)
@click.option("--rasterise", is_flag=True)
@click.option(
    "--kmer_sizes",
    type=str,
    help="kmer sizes to compute; comma-delimited",
    default="10",
    show_default=True,
)
def main(
    msa_fname,
    mosaic_alignment_fname,
    dir_ofname,
    sqlite_db_fpath,
    seqtype,
    seq_region,
    rasterise,
    kmer_sizes,
):
    con = sqlite3.connect(sqlite_db_fpath)

    # Load msa sequences
    start, end = map(int, seq_region.split(REGION_DELIM))
    donors = parse_mosaic_alignment_file(
        {}, mosaic_alignment_fname, get_names_only=True
    )
    all_targets = set(donors.keys())
    record_dict = dict()  # must be an ordered dict.
    for record in SeqIO.parse(msa_fname, "fasta"):
        record_ID = str(record.id)
        if record_ID in all_targets:
            record_dict[record_ID] = str(record.seq)[start - 1 : end]
    pre_names = set(record_dict)
    record_dict = sort_records_by_sname(record_dict)
    assert set(record_dict) == pre_names
    nrow = len(record_dict)
    ncol = len(next(iter(record_dict.values())))
    assert all(map(lambda seq: len(seq) == ncol, record_dict.values()))

    # Parse mosaic aligner file, populate sqlite table
    donors = parse_mosaic_alignment_file(record_dict, mosaic_alignment_fname)
    table_name = Path(mosaic_alignment_fname).stem
    donors.write_tsv(con, sqlite_db_fpath, table_name)

    hlines = get_hlines(record_dict)

    colour_mapping = load_sequence_colour_maps(con, seqtype)

    made_sequence_svgs, best_mismatches = write_svg_sequence_grids(
        Path(dir_ofname) / "sequence_grids",
        colour_mapping,
        donors,
        record_dict,
        nrow,
        ncol,
        hlines,
    )
    write_best_mismatches(best_mismatches, Path(dir_ofname) / "sequence_grids")

    input_table_name = get_sqlite_table_name(
        msa_fname, "DBs", seqtype, metric_info=VALID_SQLITE_DATA[2]
    )
    kmer_sharing_mapping = load_kmer_sharing_mapping(input_table_name, con)
    colour_mapping = {"DBLMSP": "#1B9E77", "DBLMSP2": "#D95F02", "shared": "#FFD92F"}
    made_sharing_svgs = write_svg_sharing_grids(
        Path(dir_ofname) / "sharing_grids",
        colour_mapping,
        kmer_sharing_mapping,
        donors,
        record_dict,
        nrow,
        ncol,
        hlines,
        start,
        kmer_sizes,
    )
    con.close()
    if rasterise:
        rasterise_and_combine(made_sequence_svgs, record_dict.keys())
        rasterise_and_combine(made_sharing_svgs, record_dict.keys())
    to_pdf(made_sequence_svgs)


if __name__ == "__main__":
    main()
