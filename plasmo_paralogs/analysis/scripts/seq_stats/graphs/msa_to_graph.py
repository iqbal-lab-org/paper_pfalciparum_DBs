"""
Generates cytoscape-viewable graph, in json format, from a MSA + recombination
breakpoints

The breakpoints give the size of the peptides that will be represented, and have been
parsed from mosaic aligner and loaded into sqlite db.
"""
import json
import sys
from dataclasses import dataclass, field
from collections import Counter, defaultdict
from typing import Set, Tuple, Mapping
from itertools import count, chain, tee
import sqlite3

import click
import pandas as pd

from common_utils.msa import split_alignment_by_gene, AlignIO, MSA, get_gene_name

LabelIDs = Mapping[str, int]


def load_recomb_breakpoints(db_fname, db_table_name):
    con = sqlite3.connect(db_fname)
    df_breakpoints = pd.read_sql_query(f"SELECT * from {db_table_name}", con)
    con.close()
    bpoints = dict()
    for row in df_breakpoints.iterrows():
        d = row[1]
        target = d.target
        donor = d.donor
        if not d.end_is_recomb_breakpoint:
            continue
        gene_name = get_gene_name(donor, already_ID=True)
        if target not in bpoints:
            bpoints[target] = (set([gene_name]), [d.end])
        else:
            bpoints[target][0].add(gene_name)
            bpoints[target][1].append(d.end)
    result = defaultdict(set)
    for key, val in bpoints.items():
        if len(val[0]) == 1:
            result[get_gene_name(key, already_ID=True)].update(val[1])
    return result


def pairwise(iterable):
    a, b = tee(iterable, 2)
    next(b, None)
    return zip(a, b)


def triplewise(iterable):
    a, b, c = tee(iterable, 3)
    next(b, None)
    next(c, None)
    next(c, None)
    return zip(a, b, c)


def get_freqs(alignment: MSA, start_pos: int, end_pos: int):
    sub_alignment = alignment[:, start_pos:end_pos]
    return Counter(map(lambda elem: str(elem.seq), sub_alignment))


@dataclass
class Node:
    ID: str
    x_pos: int
    label: str = ""
    y_pos: int = None
    gene_name: str = ""
    freq: int = 0

    def __hash__(self):
        return hash(self.ID)

    def __eq__(self, other):
        return self.ID == other.ID

    def __lt__(self, other):
        return self.ID < other.ID

    def to_dict(self):
        position_dict = dict(x=self.x_pos)
        if self.y_pos is not None:
            position_dict["y"] = self.y_pos
        else:
            print("warning, serialising a node with no y position", file=sys.stderr)
        return dict(
            data=dict(
                id=self.ID,
                label=self.label,
                gene_name=self.gene_name,
                frequency=self.freq,
            ),
            position=position_dict,
        )


@dataclass
class Edge:
    ID: str
    source: int
    target: int
    freq: int

    def to_dict(self):
        return dict(
            data=dict(
                id=self.ID, source=self.source, target=self.target, frequency=self.freq
            )
        )

    def __hash__(self):
        return hash(self.source + self.target)

    def __eq__(self, other):
        return self.source == other.source and self.target == other.target

    def __lt__(self, other):
        return self.ID < other.ID


Nodes = Set[Node]
Edges = Set[Edge]


@dataclass
class Graph:
    nodes: Nodes = field(default_factory=set)
    edges: Edges = field(default_factory=set)
    _node_ID: int = 0
    _edge_ID: int = 0
    _positional_ID_map = dict()

    def flush(self, pos: int):
        self._positional_ID_map = dict()

    def next_node_ID(self) -> str:
        result = f"n{self._node_ID}"
        self._node_ID += 1
        return result

    def next_edge_ID(self) -> str:
        result = f"e{self._edge_ID}"
        self._edge_ID += 1
        return result

    def add_node(self, node: Node):
        if node in self.nodes:
            raise ValueError("Node already in graph")
        self._positional_ID_map[f"{node.x_pos}:{node.label}{node.gene_name}"] = node.ID
        self.nodes.add(node)

    def retrieve_node_ID(self, x_pos: int, label: str) -> int:
        query = f"{x_pos}:{label}"
        return self._positional_ID_map[query]

    def add_edge(self, edge: Edge):
        if edge in self.edges:
            raise ValueError("Edge already exists")
        if (
            Node(edge.source, 0) not in self.nodes
            or Node(edge.target, 0) not in self.nodes
        ):
            raise ValueError("Edge refers to at least one non-existent node")
        self.edges.add(edge)

    def cytoscape_serialise(self):
        elems = list(
            map(
                lambda elem: elem.to_dict(),
                chain.from_iterable([self.nodes, self.edges]),
            )
        )
        for elem in elems:
            if "position" in elem:
                position_dict = elem["position"]
                position_dict["x"] *= 100
                position_dict["y"] *= 100
        result = dict(elements=list(elems))
        return json.dumps(result)


def make_graph(alignment: MSA, breakpoints_by_gene) -> Graph:
    sub_alignments = split_alignment_by_gene(alignment)
    graph = Graph()
    starting_y_pos = 0
    for gene_name, gene_alignment in sub_alignments.items():
        max_y_pos = starting_y_pos
        breakpoints = breakpoints_by_gene[gene_name]
        breakpoints.update(({0, len(gene_alignment[0]) - 1}))
        breakpoints = sorted(breakpoints)
        # Build the nodes
        for start, stop in pairwise(breakpoints):
            for xpos in range(start,stop):
                rolling_y_pos = starting_y_pos
                node_labels = get_freqs(gene_alignment, xpos, xpos+1)
                for label, freq in node_labels.items():
                    graph.add_node(
                        Node(
                            ID=graph.next_node_ID(),
                            label=label,
                            x_pos=xpos,
                            y_pos=rolling_y_pos,
                            gene_name=gene_name,
                            freq=freq,
                        )
                    )
                    rolling_y_pos += 1
                    max_y_pos = max(max_y_pos, rolling_y_pos)
                if xpos == start:
                    continue
                edge_labels = get_freqs(gene_alignment, xpos-1, xpos+1)
                for label, freq in edge_labels.items():
                    edge_source = graph.retrieve_node_ID(
                            xpos-1, label[0] + gene_name
                    )
                    edge_target = graph.retrieve_node_ID(
                        xpos, label[1]+ gene_name
                    )
                    graph.add_edge(
                        Edge(
                            ID=graph.next_edge_ID(),
                            source=edge_source,
                            target=edge_target,
                            freq=freq,
                        )
                    )
        starting_y_pos = max_y_pos + 1
    return graph


@click.command()
@click.option(
    "-r",
    "--region",
    help="0-based inclusive region to analyse, as" "'start:stop'",
    default=None,
)
@click.option(
    "--fout",
    type=click.File("w"),
    default=sys.stdout,
    help="File to write to. Default: stdout",
)
@click.argument("msa_fname")
@click.option("-s", "--sqlite_db_fpath", required=True)
@click.option("-n", "--sqlite_db_mosaic_table_name", required=True)
def main(region, fout, msa_fname, sqlite_db_fpath, sqlite_db_mosaic_table_name):
    with open(msa_fname) as fhandle_in:
        alignment = AlignIO.read(fhandle_in, "fasta")
    breakpoints_by_gene = load_recomb_breakpoints(
        sqlite_db_fpath, sqlite_db_mosaic_table_name
    )

    graph = make_graph(alignment, breakpoints_by_gene)
    print(graph.cytoscape_serialise(), file=fout)


if __name__ == "__main__":
    main()
