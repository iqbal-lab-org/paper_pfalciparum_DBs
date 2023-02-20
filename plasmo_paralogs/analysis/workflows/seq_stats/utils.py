def sqs_get_mosaic_seqtype(wildcards):
    if wildcards.seqtype == SEQTYPES[0]:
        return "nt"
    elif wildcards.seqtype == SEQTYPES[1]:
        return "aa"
    else:
        raise ValueError()
