import sqlite3
from typing import Set, Dict, Optional, Tuple
from hashlib import md5
import re

from loguru import logger
import pandas as pd

from common_utils import VALID_SEQTYPES, VALID_TOOLNAMES, VALID_SQLITE_DATA, REGION_DELIM, ID_DELIM


def md5sum_of_file(fname: str):
    contents = open(fname,"rb").read()
    digest = md5(contents).hexdigest()
    return digest

def get_sqlite_table_name(fname: str, gene_name: str, seqtype: str, metric_info: str, tool_name: str = None, seq_region = None):
    info_regexp = re.compile("|".join(VALID_SQLITE_DATA))
    assert info_regexp.match(metric_info) is not None
    assert seqtype in VALID_SEQTYPES
    if tool_name is not None:
        assert tool_name in VALID_TOOLNAMES
        tool_name = f"{tool_name}_"
    else:
        tool_name = ""
    if seq_region is not None:
        seq_region = seq_region.replace(REGION_DELIM, ID_DELIM)
        seq_region = f"{seq_region}_"
    else:
        seq_region = ""
    digest = md5sum_of_file(fname)
    return f"{tool_name}{gene_name}_{seqtype}_{seq_region}{metric_info}_{digest}"

def drop_if_exists(sqlite_db_fpath, db_cursor,table_name):
    # TODO: if table already exists, log that it is being destroyed
    db_cursor.execute(f"drop table if exists {table_name}")
    logger.info(f"Filling in table {table_name} in {sqlite_db_fpath}")

class SampleGeospatial:
    TABLE_PF6 = "pf6_metadata"
    def __init__(self, country: str, year: int):
        self.country = country
        self.year = year

def get_sample_geospatial(sqlite_db_fpath) -> Dict[str, SampleGeospatial]:
    result = dict()
    con = sqlite3.connect(sqlite_db_fpath)
    con.row_factory = sqlite3.Row
    db_cursor = con.cursor()
    for row in db_cursor.execute(f"select * from {SampleGeospatial.TABLE_PF6}"):
        result[row["Sample"]] = SampleGeospatial(row["Country"],row["Year"])
    con.close()
    lab_strains = {"Pf7G8","PfCD01","PfDd2","PfGA01","PfGB4","PfGN01","PfHB3","PfIT","PfKE01","PfKH01","PfKH02","PfML01","PfSD01","PfSN01","PfTG01","ref"}
    for strain in lab_strains:
        result[strain] = SampleGeospatial("lab_strain","NA")
    return result

def get_sample_ids_matching_geo_region(sqlite_db_fpath, geo_region: Optional[str]) -> Set[str]:
    result = set()
    if geo_region is None:
        return result
    con = sqlite3.connect(sqlite_db_fpath)
    con.row_factory = sqlite3.Row
    db_cursor = con.cursor()
    geo_regexp = re.compile(geo_region)

    for row in db_cursor.execute(f"select * from {SampleGeospatial.TABLE_PF6}"):
        if (
            geo_regexp.match(row["Site"]) is not None
            or geo_regexp.match(row["Country"]) is not None
        ):
            result.add(row["Sample"])
    con.close()
    return result

def get_DBL_coords(sqlite_db_fpath, seqtype, margin: int = 0) -> Tuple[int, int]:
    table_name = "domain_coordinates"
    con = sqlite3.connect(sqlite_db_fpath)
    df = pd.read_sql_query(f"select * from {table_name}",con,index_col="Feature")
    con.close()
    column_name = f"{seqtype}_MSA_coords"
    result = df[column_name]["DBL_operational"]
    result = list(map(int,result.split(":")))
    return result[0] - margin, result[1] + margin


def load_percent_identities(sqlite_db_fpath, table_name, gene_ID = None):
    """
    Table computed in paralog_diversity_and_divergence.py
    """
    result = dict()
    con = sqlite3.connect(sqlite_db_fpath)
    con.row_factory = sqlite3.Row
    db_cursor = con.cursor()
    for row in db_cursor.execute(f"select * from {table_name}"):
        if gene_ID is not None:
            if row["gene_ID"] != gene_ID:
                continue
        result[row["sample_ID"]] = float(row["percent_identity"])
    return result
