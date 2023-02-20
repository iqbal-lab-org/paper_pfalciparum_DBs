def as_get_tool_name(wildcards):
    tn = wildcards.tool_name
    if tn == TOOL_NAMES[0]:
        return "gram_jointgeno_ebf7bcd5__pf6__analysis_set_fws95__gapfiller__7__13"
    elif tn == TOOL_NAMES[1]:
        return "pf6"
    else:
        raise ValueError(f"Unsupported tool name: {tn}")
