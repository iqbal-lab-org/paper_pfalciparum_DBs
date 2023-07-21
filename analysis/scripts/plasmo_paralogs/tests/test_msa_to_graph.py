from typing import List
from itertools import starmap

from test_utils import make_alignment, TEST_GENE_NAME
from plasmo_paralogs.seq_stats.msa_to_graph import Node, Edge, make_graph


def objects_are_equal(obj1, obj2):
    return obj1.__dict__ == obj2.__dict__

def lists_are_equal(list1: List, list2: List):
    if len(list1) != len(list2):
        return False
    return all(starmap(objects_are_equal, zip(sorted(list1),sorted(list2))))


class TestSimpleGraphs_SameGenes:
    def test_two_nodes_one_edge(self):
        alignment = make_alignment(["AA","AA"])
        expected_nodes = [
                Node(ID="n0",x_pos=0,y_pos=0,label="A",freq=2,gene_name=TEST_GENE_NAME),
                Node(ID="n1",x_pos=1,y_pos=0,label="A",freq=2,gene_name=TEST_GENE_NAME),
                ]

        expected_edges = [
                Edge(source="n0",target="n1",ID="e0",freq=2)
                ]

        result = make_graph(alignment, 0, len(alignment))
        assert lists_are_equal(result.nodes, expected_nodes)
        assert lists_are_equal(result.edges, expected_edges)

    def test_three_nodes_two_edges(self):
        alignment = make_alignment(["AT","AA"])
        result = make_graph(alignment, 0, len(alignment))

        expected_nodes = [
                Node(ID="n0",x_pos=0,y_pos=0,label="A",freq=2,gene_name=TEST_GENE_NAME),
                Node(ID="n1",x_pos=1,y_pos=0,label="T",freq=1,gene_name=TEST_GENE_NAME),
                Node(ID="n2",x_pos=1,y_pos=1,label="A",freq=1,gene_name=TEST_GENE_NAME),
                ]

        expected_edges = [
                Edge(source="n0",target="n1",ID="e0",freq=1),
                Edge(source="n0",target="n2",ID="e1",freq=1),
                ]

        result = make_graph(alignment, 0, len(alignment))
        assert lists_are_equal(result.nodes, expected_nodes)
        assert lists_are_equal(result.edges, expected_edges)

class TestSimpleGraphs_DifferentGenes:
    def test_four_nodes_two_edges(self):
        alignment = make_alignment(["TT","AA"],ids=["gene1_s1","gene2_s1"])
        result = make_graph(alignment, 0, len(alignment))

        expected_nodes = [
                Node(ID="n0",x_pos=0,y_pos=0,label="T",freq=1,gene_name="gene1"),
                Node(ID="n1",x_pos=0,y_pos=1,label="A",freq=1,gene_name="gene2"),
                Node(ID="n2",x_pos=1,y_pos=0,label="T",freq=1,gene_name="gene1"),
                Node(ID="n3",x_pos=1,y_pos=1,label="A",freq=1,gene_name="gene2"),
                ]
        expected_edges = [
                Edge(source="n0",target="n2",ID="e0",freq=1),
                Edge(source="n1",target="n3",ID="e1",freq=1),
                ]

        result = make_graph(alignment, 0, len(alignment))
        assert lists_are_equal(result.nodes, expected_nodes)
        assert lists_are_equal(result.edges, expected_edges)

