ID_DELIM = "_"
REGION_DELIM = ":"
SHARED_FEATURE_DELIM = "@"
REGION_CLICK_HELP = f"Region to focus on, format: 'start{REGION_DELIM}stop'."
VALID_SEQTYPES = ["protein", "dna"]
VALID_TOOLNAMES = ["gramtools", "malariaGEN"]

VALID_SQLITE_DATA = [
    "homozygosity",
    "sharing_def_by_geo",
    "sharing_assignment_by_sample",
    "num_used_samples_by_country",
    "positional_matching",
    "sample_identities",
]
