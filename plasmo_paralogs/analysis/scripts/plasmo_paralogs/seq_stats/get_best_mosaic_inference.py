import glob
import click
import sqlite3
import shutil
from io import StringIO
from pathlib import Path

import pandas as pd

from plasmo_paralogs.common_utils.metrics import MetricsRecorder
from plasmo_paralogs.common_utils import ID_DELIM
from plasmo_paralogs.common_utils.sqlite import drop_if_exists


class CustomMetric(MetricsRecorder):
    _headers = ["recomb_value", "seq_num", "log_likelihood"]
    _required_headers = _headers


@click.command()
@click.argument("input_dirname")
@click.argument("best_inference_fname")
@click.option("-s", "--sqlite_db_fpath", required=True)
def main(input_dirname, best_inference_fname, sqlite_db_fpath):
    result = list()
    file_dict = dict()
    for fname in glob.glob(f"{input_dirname}/*_log.txt"):
        recomb_value = float(Path(fname).name.split("_")[0])
        file_dict[recomb_value] = fname
        seq_num = 0
        for line in open(fname, "r").readlines():
            if "Maximum Log likelihood  =" in line:
                result.append(
                    CustomMetric(
                        recomb_value=recomb_value,
                        seq_num=seq_num,
                        log_likelihood=float(line.split("=")[1].strip()),
                    )
                )
                seq_num += 1
    result_string = CustomMetric.get_header()
    con = sqlite3.connect(sqlite_db_fpath)
    db_cursor = con.cursor()
    table_name = f"mosaic{ID_DELIM}recomb{ID_DELIM}{str(Path(best_inference_fname).stem)}{ID_DELIM}"
    drop_if_exists(sqlite_db_fpath, db_cursor, table_name)
    db_cursor.execute(
        f"create table {table_name} (recomb_value, seq_num, log_likelihood)",
    )
    for elem in result:
        elem_string = str(elem)
        db_cursor.execute(
            f"insert into {table_name} values (?,?,?)",
            elem_string.split("\t"),
        )
        result_string += elem_string
    con.commit()
    con.close()

    ## Find best recomb value, and copy the file it is linked to
    df = pd.read_csv(StringIO(result_string), sep="\t")
    df_sums = df.groupby(["recomb_value"]).sum()
    assert len(set(df_sums["seq_num"])) == 1
    ml_recomb_value = df_sums.idxmax()["log_likelihood"]
    best_file = file_dict[ml_recomb_value]
    best_file = best_file.replace("_log", "_align")
    shutil.copy(best_file, best_inference_fname)


if __name__ == "__main__":
    main()
