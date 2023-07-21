"""
Given an alignment file (of both DBs), computes:
- for each sample the % identity of both genes (divergence)
- distribution of % identity for each gene between sample pairs (diversity)

Plotting:
- High inter-gene % ID are candidate gene conversion events. This is plotted.
- Two main conversion clusters appear- these are plotted by how much they appear 
across sample countries
"""
import itertools as it
import sqlite3
from pathlib import Path
from collections import defaultdict
import subprocess
import random
import operator as op

import click
from pyfaidx import Fasta
import numpy as np
import pandas as pd
import seaborn as sns
import plotly.express as px
from scipy.cluster.hierarchy import fcluster
import matplotlib.pyplot as plt
from Bio import AlignIO, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from plasmo_paralogs.common_utils import ID_DELIM, VALID_SQLITE_DATA, VALID_SEQTYPES
from plasmo_paralogs.common_utils import REGION_DELIM, REGION_CLICK_HELP
from plasmo_paralogs.common_utils.msa import get_gene_name, get_sample_id
from plasmo_paralogs.common_utils.sqlite import (
    get_sample_geospatial,
    get_DBL_coords,
    get_sqlite_table_name,
    drop_if_exists,
)
from plasmo_paralogs.common_utils.metrics import MetricsRecorder
from plasmo_paralogs.seq_stats.sharing.hmm_private_shared_logos import process_all_gap_columns
from plasmo_paralogs.analysed_sequences.msa_coordinate_adjust import coordinate_adjust

CODON_SIZE = 3
random.seed(42)


class CustomMetric(MetricsRecorder):
    _headers = ["sample_ID", "percent_identity", "gene_ID"]
    _required_headers = _headers


def protein_region_to_dna_region(representative_seq, seq_region, sqlite_db_fpath):
    """
    Converts protein (MSA) coordinates into dna (MSA) coordinates
    Uses as yardstick the DBL coordinates in MSA space, for both protein and dna, stored in `sqlite_db_fpath`
    """
    start, end = map(int, seq_region.split(REGION_DELIM))
    gapped_dna_start, gapped_dna_end = get_DBL_coords(sqlite_db_fpath, "dna")
    prot_start, prot_end = get_DBL_coords(sqlite_db_fpath, "protein")
    margin_start = prot_start - start
    margin_end = end - prot_end
    ungapped_dna_start, ungapped_dna_end = coordinate_adjust(
        representative_seq, gapped_dna_start, gapped_dna_end, reverse=True
    )
    adjusted_start = ungapped_dna_start - margin_start * 3
    adjusted_end = ungapped_dna_end + margin_end * 3
    result_start, result_end = coordinate_adjust(
        representative_seq, adjusted_start, adjusted_end
    )
    return (result_start, result_end)


def perc_ids_from_sample_pairs(seqs, used_sample_ids, gene_ids, seq_region):
    """
    gene_ID:
        - if a single gene ID, compares the same gene in randomly picked sample pairs
          (in this case the samples have to be different - else sequences are identical)
        - if two gene IDs, compares the two different genes in randomly picked sample pairs
          (in this case the samples can be the same)
    """
    assert 0 < len(gene_ids) <= 2
    multi_gene = len(gene_ids) == 2
    start, end = map(int, seq_region.split(REGION_DELIM))
    region_len = (end - start) // 3
    percent_identities = list()
    sample_IDs = list()
    while len(percent_identities) < len(used_sample_ids):
        sample_pair = random.choices(used_sample_ids, k=2)
        if not multi_gene and sample_pair[0] == sample_pair[1]:
            continue
        seq_pair = list()
        for i, sample_ID in enumerate(sample_pair):
            if multi_gene:
                fasta_ID = f"{gene_ids[i]}{ID_DELIM}{sample_ID}"
            else:
                fasta_ID = f"{gene_ids[0]}{ID_DELIM}{sample_ID}"
            seq_pair.append(np.array(codonify(seqs[fasta_ID][start:end].seq)))
        sample_IDs.append(ID_DELIM.join(sample_pair))
        percent_identities.append(sum(seq_pair[0] == seq_pair[1]) / region_len)
    return pd.DataFrame(
        {
            "sample_ID": sample_IDs,
            "percent_identity": percent_identities,
            "gene_ID": "DBs (randomly-assorted pairs)" if multi_gene else gene_ids[0],
        }
    )


def write_perc_id_table(
    sqlite_db_fpath, frac_shared_df, msa_fname, gene_ID, seq_region
):
    con = sqlite3.connect(sqlite_db_fpath)
    table_name = get_sqlite_table_name(
        msa_fname,
        gene_ID,
        VALID_SEQTYPES[1],
        metric_info=VALID_SQLITE_DATA[5],
        seq_region=seq_region,
    )
    frac_shared_df.to_sql(table_name, con, if_exists="replace", index=False)
    con.commit()
    con.close()


def codonify(seq: str):
    result = list()
    for i in range(0, len(seq), CODON_SIZE):
        result.append(seq[i : i + CODON_SIZE])
    return result


def get_conversion_clusters_by_country(sqlite_db_fpath, clustermap):
    sample_geospatial = get_sample_geospatial(sqlite_db_fpath)
    fl = fcluster(clustermap.dendrogram_row.linkage, 2, criterion="maxclust")
    clusters = {i + 1: defaultdict(int) for i in range(max(set(fl)))}
    clusters_samples = {i: [] for i in clusters.keys()}
    for sname, cluster in zip(clustermap.dendrogram_row.data.index, fl):
        country = sample_geospatial[sname].country
        clusters[cluster][country] += 1
        clusters_samples[cluster].append(sname)

    result = list()
    for cluster, country_dict in clusters.items():
        for country, count in country_dict.items():
            result.append([cluster, country, count])
    result = pd.DataFrame(
        result, columns=["Gene conversion cluster", "Country", "Sample count"]
    )
    return result, clusters_samples


@click.command()
@click.argument("msa_fname")
@click.argument("sqlite_db_fpath")
@click.argument("out_dirname")
@click.option(
    "-p",
    "protein_fasta",
    help="protein fasta for extracting conversion clusters",
    required=True,
)
@click.option("--seq_region", type=str, required=True, help=REGION_CLICK_HELP)
@click.option(
    "--plot",
    is_flag=True,
    help="Whether to make plots of (putative) gene conversion clusters",
)
@click.option(
    "--fold_changes_tsv","-f",
    required=True,
    help="Contains fold-changes computed from ir stats tsv",
)
def main(msa_fname, sqlite_db_fpath, out_dirname, protein_fasta, seq_region, plot, fold_changes_tsv):
    """
    `msa_fname`: dna-level MSA, as gene conversion is best observed at DNA level
    `seq_region`: in equivalent protein-level MSA space. conversion to DNA-space is done
    automatically inside this script
    """
    seqs = Fasta(msa_fname)
    representative_seq = seqs["DBLMSP_ref"][:].seq
    start, end = protein_region_to_dna_region(
        representative_seq, seq_region, sqlite_db_fpath
    )
    region_len = (end - start) // 3
    gene_ids = set()
    sample_ids = set()
    used_sample_ids = list()
    matrix = list()
    for fasta_id in seqs.keys():
        gene_ids.add(get_gene_name(fasta_id, True))
        sample_ids.add(get_sample_id(fasta_id, True))
    for sample_id in sample_ids:
        sample_seqs = list()
        for pre_fasta_id in it.product(gene_ids, [ID_DELIM], [sample_id]):
            fasta_id = "".join(pre_fasta_id)
            if fasta_id not in seqs:
                continue
            sample_seqs.append(np.array(codonify(seqs[fasta_id][start:end].seq)))
        if len(sample_seqs) != 2:
            continue
        used_sample_ids.append(sample_id)
        identities = sample_seqs[0] == sample_seqs[1]
        matrix.append(identities.astype(int))
    matrix = np.array(matrix)
    df = pd.DataFrame(matrix, index=used_sample_ids)
    print(f"Number of confidently-resolved pairs: {len(used_sample_ids)}")

    frac_shared = df.sum(axis=1) / region_len
    frac_shared_df = pd.DataFrame(
        {
            "sample_ID": frac_shared.index,
            "percent_identity": frac_shared,
            "gene_ID": "DBs",
        }
    )
    dna_seq_region = f"{start}{REGION_DELIM}{end}"
    # Identity level between same genes on randomly picked sample pairs
    for gene_ID in gene_ids:
        frac_shared_df = pd.concat(
            [
                frac_shared_df,
                perc_ids_from_sample_pairs(
                    seqs, used_sample_ids, [gene_ID], dna_seq_region
                ),
            ]
        )
    # Identity level between different genes on randomly picked sample pairs
    frac_shared_df = pd.concat(
        [
            frac_shared_df,
            perc_ids_from_sample_pairs(seqs, used_sample_ids, list(gene_ids), dna_seq_region),
        ]
    )
    write_perc_id_table(sqlite_db_fpath, frac_shared_df, msa_fname, "DBs", seq_region)

    if plot:
        input_start, input_end = map(int, seq_region.split(REGION_DELIM))
        out_dirname = Path(out_dirname) / seq_region.replace(REGION_DELIM, ID_DELIM)
        out_dirname.mkdir(exist_ok=True, parents=True)
        frac_shared = frac_shared.sort_values(ascending=False)
        fig = plt.figure()
        p = sns.histplot(frac_shared)
        p.set_ylabel("Number of samples")
        p.set_xlabel("Fraction of DBL shared (dna)")
        p.get_figure().savefig(f"{out_dirname}/fraction_shared_histogram.pdf")

        fold_changes = pd.read_csv(fold_changes_tsv, sep="\t")
        high_share_samples = frac_shared[frac_shared > 0.5]
        print(f"Number of confidently-resolved pairs with high frac_shared: {len(high_share_samples)}")
        h_s_snames = set(high_share_samples.index)
        DBs={"DBLMSP","DBLMSP2"}
        fold_changes = fold_changes[fold_changes["gene"].isin(DBs)]
        for comparator, fold_change in zip((op.lt,op.gt),(0.5,1.5)):
            putative_cnvs = fold_changes[comparator(fold_changes["mean_ratio"],fold_change)]
            filtered_out = set(h_s_snames.intersection(putative_cnvs["sample"]))
            print(f"Filtered out following putative gene conversion samples with putative CNV ({fold_change}-fold): {filtered_out}")
            high_share_samples = high_share_samples[~ high_share_samples.index.isin(filtered_out)]
        print(f"Number of confidently-resolved pairs with high frac_shared (post filtering): {len(high_share_samples)}")


        clustermap = sns.clustermap(
            df.loc[high_share_samples.index, :],
            col_cluster=False,
            yticklabels=False,
            cbar_pos=None,
            dendrogram_ratio=(0.2,0),
            figsize=(12,10)
        )
        clustermap.savefig(f"{out_dirname}/high_share_samples_clustered.pdf")

        conversion_df, clusters_samples = get_conversion_clusters_by_country(
            sqlite_db_fpath, clustermap
        )
        p = px.scatter_geo(
            conversion_df,
            facet_row="Gene conversion cluster",
            locations="Country",
            size="Sample count",
            locationmode="country names",
            color="Sample count",
        )
        p.write_image(f"{out_dirname}/gene_conv_clusters_map.pdf",width=800,height=600)

        prot_seqs = Fasta(protein_fasta)
        for cluster_num, snames in clusters_samples.items():
            with open(f"{out_dirname}/conv_cluster_{cluster_num}.fasta", "w") as fout:
                for sname in snames:
                    try:
                        for gene_ID in gene_ids:
                            seq_id = f"{gene_ID}_{sname}"
                            seq = prot_seqs[seq_id][input_start:input_end].seq
                            fout.write(f">{seq_id}\n{seq}\n")
                    except KeyError:
                        print(f"Not found: {sname}")
                        continue

        # Process all-gap columns
        for cluster_num in clusters_samples.keys():
            msa_fname = f"{out_dirname}/conv_cluster_{cluster_num}.fasta"
            recs = []
            for rec in AlignIO.read(open(msa_fname), "fasta"):
                recs.append(SeqRecord(Seq(rec.seq), id=rec.id, description=""))
            process_all_gap_columns({"_": recs})
            SeqIO.write(recs, msa_fname, "fasta")
            hmm_fname = msa_fname.replace("fasta", "hmm")
            subprocess.run(
                f"hmmbuild --fragthresh 0 --enone --symfrac 0 --pnone {hmm_fname} {msa_fname}",
                shell=True,
                check=True,
            )


if __name__ == "__main__":
    main()
