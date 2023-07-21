import pandas as pd

def lav_load_assembly_names():
    df = pd.read_csv(config["input_tsvs"]["laverania_assemblies"], sep="\t", index_col=0)
    result = dict(zip(df["Abbreviation"],df["Accession"]))
    return result
