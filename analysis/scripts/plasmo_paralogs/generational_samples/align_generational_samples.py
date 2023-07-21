import csv
import re
from enum import Enum, auto
from typing import Dict
from copy import deepcopy
from pathlib import Path

from Bio import SeqIO
import edlib
import click
import pandas as pd

from plasmo_paralogs.common_utils.metrics import MetricsRecorder
from analysis_set.filter_genes import Filter, HARDCODED_FILTERS

# used_filters = re.compile(".*positions_0x.*|.*disagreeing_pileup.*")
# USED_FILTERS = [elem for elem in HARDCODED_FILTERS if used_filters.match(elem[0])]
USED_FILTERS = HARDCODED_FILTERS
FIELD_DELIM = "_"


class CloneTreeMetric(MetricsRecorder):
    _headers = [
        "Parent_Sample_ID",
        "Child_Sample_ID",
        "Parent_Tree_ID",
        "Child_Tree_ID",
        "Gene",
        "edit_distance",
        "cigar",
        "passes_all_filters",
    ]
    _required_headers = _headers


class CrossMetric(MetricsRecorder):
    _headers = [
        "Child_Sample_ID",
        "Cross",
        "Gene",
        "Parent_1_Sample_ID",
        "Parent_1_edit_distance",
        "Parent_1_cigar",
        "Parent_2_Sample_ID",
        "Parent_2_edit_distance",
        "Parent_2_cigar",
        "passes_all_filters",
    ]
    _required_headers = _headers


class FilteredSamples(MetricsRecorder):
    _headers = ["Sample_ID", "Gene", "Filter_reason"]
    _required_headers = _headers


class TreeNode:
    def __init__(self, Sample_ID, Tree_ID):
        self.Sample_ID = Sample_ID
        self.Tree_ID = Tree_ID

    def _align(self, other: "TreeNode"):
        return edlib.align(self.seq, other.seq, task="path")

    def assign_seq(self, seq):
        self.seq = seq


TreeDict = Dict[str, TreeNode]


class ParentNode(TreeNode):
    def add_child(self, other: TreeNode):
        if getattr(self, "children", None) is None:
            self.children = list()
        self.children.append(other)

    def align_to_children(self):
        for ch in self.children:
            if hasattr(ch,"seq"):
                alignment = self._align(ch)
                result_dict = {
                    "Parent_Sample_ID": self.Sample_ID,
                    "Child_Sample_ID": ch.Sample_ID,
                    "Parent_Tree_ID": self.Tree_ID,
                    "Child_Tree_ID": ch.Tree_ID,
                    "edit_distance": alignment["editDistance"],
                    "cigar": alignment["cigar"],
                }
                yield result_dict



class ChildNode(TreeNode):
    def add_parent(self, other: TreeNode):
        if getattr(self, "parents", None) is None:
            self.parents = list()
        self.parents.append(other)

    def align_to_parents(self):
        assert len(self.parents) == 2
        result_dict = {"Child_Sample_ID": self.Sample_ID}
        for i in range(2):
            alignment = self._align(self.parents[i])
            one_based = i + 1
            result_dict.update(
                {
                    f"Parent_{one_based}_Sample_ID": self.parents[i].Sample_ID,
                    f"Parent_{one_based}_edit_distance": alignment["editDistance"],
                    f"Parent_{one_based}_cigar": alignment["cigar"],
                }
            )
        result_dict["Cross"] = f"{self.parents[0].Tree_ID}_{self.parents[1].Tree_ID}"
        return result_dict



def get_filter_passing_samples(df_stats_ir, gene_name, tool_name):
    df_ir = Filter(df_stats_ir, "tool", tool_name, "==").apply(df_stats_ir)
    df_ir = Filter(df_ir, "gene", gene_name, "==").apply(df_ir)
    filters = [Filter(*(df_ir, *elem)) for elem in USED_FILTERS]
    for f in filters:
        df_ir = f.apply(df_ir)
    kept_samples = set(df_ir.index)
    return kept_samples


def load_sequences(seq_ifname):
    result = {}
    with open(seq_ifname) as fhandle_in:
        for record in SeqIO.parse(fhandle_in, "fasta"):
            result[record.id.strip()] = str(record.seq)
    return result


class GenMode(Enum):
    CLONE_TREE = auto()
    CROSSES = auto()

    @staticmethod
    def from_tsv(generational_tsv_fname):
        if "clone_tree" in generational_tsv_fname:
            return GenMode.CLONE_TREE
        elif "crosses" in generational_tsv_fname:
            return GenMode.CROSSES
        else:
            raise ValueError("Unsupported")


def load_generational_tsv(generational_tsv_fname) -> TreeDict:
    mode = GenMode.from_tsv(generational_tsv_fname)
    node_dict = dict()
    # First pass: make each sample a node in a tree
    with open(generational_tsv_fname) as fhandle_in:
        row_reader = csv.DictReader(fhandle_in, delimiter="\t")
        for row in row_reader:
            if mode is GenMode.CLONE_TREE:
                Sample_ID = row["BAM_ID"]
                Tree_ID = row["Tree_ID"]
            elif mode is GenMode.CROSSES:
                Sample_ID = row["Sample"]
                Tree_ID = row["Clone"]
            node_dict[Tree_ID] = TreeNode(Sample_ID, Tree_ID)
    # Second pass: create edges in the tree between parents and children
    if mode is GenMode.CROSSES:
        assign_parents(generational_tsv_fname, node_dict)
    elif mode is GenMode.CLONE_TREE:
        assign_children(generational_tsv_fname, node_dict)
    else:
        raise ValueError("Unknown mode")
    return node_dict


def assign_children(generational_tsv_fname, node_dict) -> None:
    with open(generational_tsv_fname) as fhandle_in:
        row_reader = csv.DictReader(fhandle_in, delimiter="\t")
        for row in row_reader:
            Tree_ID = row["Tree_ID"]
            Parent_Tree_ID = row["Parent_Tree_ID"]
            is_child = Parent_Tree_ID != "NA"
            if Tree_ID not in node_dict or not is_child:
                continue
            child_node = node_dict[Tree_ID]
            parent_node = node_dict[Parent_Tree_ID]
            if not isinstance(parent_node, ParentNode):
                parent_node = ParentNode(node_dict[Parent_Tree_ID].Sample_ID, Parent_Tree_ID)
                node_dict[Parent_Tree_ID] = parent_node
            parent_node.add_child(child_node)


def assign_parents(generational_tsv_fname, node_dict) -> None:
    with open(generational_tsv_fname) as fhandle_in:
        row_reader = csv.DictReader(fhandle_in, delimiter="\t")
        for row in row_reader:
            Tree_ID = row["Clone"]
            cross_name = row["Cross"]
            is_child = Tree_ID not in cross_name.upper()
            if not is_child:
                continue
            child_node = ChildNode(node_dict[Tree_ID].Sample_ID, Tree_ID)
            parent_IDs = cross_name.upper().split("_")
            assert len(parent_IDs) == 2
            if any(map(lambda par_ID: par_ID not in node_dict, parent_IDs)):
                continue
            for parent_ID in parent_IDs:
                child_node.add_parent(node_dict[parent_ID])
            node_dict[Tree_ID] = child_node


def get_alignments(stats_class, node_dict: TreeDict):
    result = []
    if stats_class is CrossMetric:
        for node in node_dict.values():
            if isinstance(node, ChildNode):
                alignment_dict = node.align_to_parents()
                if (
                    alignment_dict["Parent_1_edit_distance"] == 0
                    or alignment_dict["Parent_2_edit_distance"] == 0
                ):
                    continue
                result.append(stats_class(**alignment_dict))
    elif stats_class is CloneTreeMetric:
        for node in node_dict.values():
            if isinstance(node, ParentNode):
                for alignment_dict in node.align_to_children():
                    if alignment_dict["edit_distance"] == 0:
                        continue
                    result.append(stats_class(**alignment_dict))
    return result


def process_one_gene(
    df_stats_ir, sample_sequences, node_dict, mode, stats_class, gene_name, tool_name
):
    filtered_samples = list()
    samples_passing_filters = get_filter_passing_samples(
        df_stats_ir, gene_name, tool_name
    )
    # Filter
    filtered_samples = []
    new_node_dict = dict()
    for tree_ID, node in node_dict.items():
        sample_ID = node.Sample_ID
        sample_seq_ID = f"{gene_name}{FIELD_DELIM}{sample_ID}"
        if sample_seq_ID in sample_sequences:
            node.assign_seq(sample_sequences[sample_seq_ID])
            new_node_dict[tree_ID] = node
        if sample_seq_ID not in sample_sequences:
            filtered_samples.append(
                FilteredSamples(
                    Sample_ID=sample_ID,
                    Gene=gene_name,
                    Filter_reason="No sequence",
                )
            )
        if sample_ID not in samples_passing_filters:
            filtered_samples.append(
                FilteredSamples(
                    Sample_ID=sample_ID,
                    Gene=gene_name,
                    Filter_reason="Failing filters",
                )
            )

    stats_with_alignments = get_alignments(stats_class, new_node_dict)
    for elem in stats_with_alignments:
        if mode is GenMode.CLONE_TREE:
            queried_IDs = [elem["Parent_Sample_ID"], elem["Child_Sample_ID"]]
        elif mode is GenMode.CROSSES:
            queried_IDs = [elem["Parent_1_Sample_ID"], elem["Parent_2_Sample_ID"], elem["Child_Sample_ID"]]
        passing_filters = all(map(lambda query_ID: query_ID in samples_passing_filters,
            queried_IDs))
        elem.update(Gene=gene_name, passes_all_filters=passing_filters)
    return stats_with_alignments, filtered_samples


@click.command()
@click.argument("seq_ifname", type=click.Path(exists=True))
@click.argument("stats_tsv", type=click.Path(exists=True))
@click.argument("generational_tsv", type=click.Path(exists=True))
@click.argument("tsv_ofname")
@click.option("--gene_to_process",help="Process only this gene")
@click.option("--tool_name", required=True)
def main(seq_ifname, stats_tsv, generational_tsv, tsv_ofname, gene_to_process, tool_name):
    """ """
    df_stats_ir = pd.read_table(stats_tsv, sep="\t", index_col=0)
    sample_sequences = load_sequences(seq_ifname)
    node_dict = load_generational_tsv(generational_tsv)
    mode = GenMode.from_tsv(generational_tsv)
    filtered_ofname = tsv_ofname.replace(".tsv","_filtered_samples.tsv")
    if mode is GenMode.CROSSES:
        stats_class = CrossMetric
    elif mode is GenMode.CLONE_TREE:
        stats_class = CloneTreeMetric
    if gene_to_process is None:
        genes_to_process = set(df_stats_ir["gene"])
    else:
        genes_to_process = {gene_to_process}

    with open(tsv_ofname, "w") as tsv_out, open(filtered_ofname, "w") as filtered_out:
        tsv_out.write(stats_class.get_header())
        filtered_out.write(FilteredSamples.get_header())
        for gene_name in genes_to_process:
            stats, filtered_samples = process_one_gene(
                df_stats_ir,
                sample_sequences,
                node_dict,
                mode,
                stats_class,
                gene_name,
                tool_name
            )
            for elem in stats:
                tsv_out.write(str(elem))
            for elem in filtered_samples:
                filtered_out.write(str(elem))


if __name__ == "__main__":
    main()
