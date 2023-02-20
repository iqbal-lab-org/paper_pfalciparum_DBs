import operator
import sys

import click
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd

#GROUPS = {"DBs": ["DBLMSP_DBL", "DBLMSP2_DBL"]}
GROUPS = {"DBs": ["DBLMSP", "DBLMSP2"]}
for gene in ["AMA1","EBA175","MSP3","MSP6","MSP2","MSP1"]:
    GROUPS[gene] = [gene]
HARDCODED_FILTERS = (
        ("fraction_positions_4x_or_less", 0.0, "<="),
        ("fraction_disagreeing_pileup_min5x", 0.0, "<="),
        ("fraction_reads_above_mean_ins_size_plus_two_std", 0.15, "<"),
        )

def get_gene_name(gene_name: str) -> str:
    if "_DBL" in gene_name:
        return gene_name.replace("_DBL","")
    else:
        return gene_name

class Filter:
    OPERATORS = {"<": operator.lt, "<=": operator.le, "==": operator.eq}

    def __init__(
        self, df: pd.DataFrame, metric: str, value: float, operator: "OPERATORS"
    ):
        assert (
            metric in df.columns
        ), f"Error, {metric} not a valid metric, choose from {df.columns}"

        self.metric = metric
        self.value = value
        assert (
            operator in self.OPERATORS.keys()
        ), f"Error, valid operators are {self.OPERATORS}"
        self.operator = self.OPERATORS[operator]

    def apply(self, df: pd.DataFrame) -> pd.DataFrame:
        return df[self.operator(df[self.metric], self.value)]


@click.command(
        help=f"""
        Filters sequences based on stats computed by eval_varcalls workflow.\n
        Current filters are hardcoded:\n
        {HARDCODED_FILTERS}
        """
        )
@click.argument("seq_ifname", type=click.Path(exists=True))
@click.argument("stats_ifname", type=click.Path(exists=True))
@click.argument("seq_ofname")
@click.option(
    "--gene_group_name", required=True, type=click.Choice(list(GROUPS.keys()))
)
@click.option("--tool_name", required=True, help="e.g.: gapfiller")
def main(seq_ifname, stats_ifname, seq_ofname, gene_group_name, tool_name):
    ## Load sequences
    iseq_pool = dict()
    with open(seq_ifname) as ifhandle:
        for record in SeqIO.parse(ifhandle, "fasta"):
            iseq_pool[record.id] = record.seq
    isize = len(iseq_pool)

    oseq_pool = list()

    ## Perform filtering
    df = pd.read_table(stats_ifname, sep="\t", index_col=0)
    df = Filter(df, "tool", tool_name, "==").apply(df)
    filters = [Filter(*(df,*elem)) for elem in HARDCODED_FILTERS]
    for gene_name in GROUPS[gene_group_name]:
        try:
            # Always add ref seq
            ref_probe_name = f"{get_gene_name(gene_name)}_ref"
            oseq_pool.append(SeqRecord(iseq_pool[ref_probe_name], id=ref_probe_name, description=""))
        except KeyError:
            pass
        filtered_df = Filter(df, "gene", gene_name, "==").apply(df)
        for f in filters:
            before = len(filtered_df)
            filtered_df = f.apply(filtered_df)
            after = len(filtered_df)
            print(f'[Gene: {gene_name}; Metric: {f.metric}] removed {before - after} additional samples')
        for name_zip in zip(list(filtered_df["gene"]), list(filtered_df.index)):
            probe_name = get_gene_name(name_zip[0]) + "_" + name_zip[1]
            if probe_name in iseq_pool:
                oseq_pool.append(
                    SeqRecord(iseq_pool[probe_name], id=probe_name, description=""),
                )

    print(
        f"Filtered out {isize - len(oseq_pool)} of {isize} input records",
        file=sys.stdout,
    )
    with open(seq_ofname, "w") as ofhandle:
        SeqIO.write(oseq_pool, ofhandle, "fasta")


if __name__ == "__main__":
    main()
