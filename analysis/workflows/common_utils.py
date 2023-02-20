import pandas as pd

TOOL_NAMES = ["gramtools","malariaGEN"]
SEQTYPES = ["dna","protein"]

def cu_mk_output_dirs(variables):
    """For each variable starting with 'output', makes the directory name it holds"""
    for variable in filter(lambda name: name.startswith("output"), variables):
        Path(eval(variable)).mkdir(exist_ok=True, parents=True)

def _cu_shift_by_window_size(region_str, seqtype):
    start, end = map(int, region_str.split(":"))
    window_size = config["domain_extra_amino_acid_window_size"]
    if seqtype == SEQTYPES[0]:
        window_size *= 3
    start -= window_size
    end += window_size
    return start, end

def cu_load_DBL_coords(seqtype):
    df = pd.read_csv(config["sqlite_fnames"]["domain_coordinates"], sep="\t", index_col=0)
    DBL_coords = {
        SEQTYPES[0]: df.loc["DBL_operational"]["dna_MSA_coords"],
        SEQTYPES[1]: df.loc["DBL_operational"]["protein_MSA_coords"],
    }
    region_string = DBL_coords[seqtype]
    return _cu_shift_by_window_size(region_string, seqtype)
