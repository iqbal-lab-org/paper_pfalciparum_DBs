from itertools import repeat, starmap

_qc_headers = ["sample", "chrom", "start", "end", "gene","tool","metric_category","metric","value"]
_required_qc_headers = set(_qc_headers[0:4] + _qc_headers[6:9])

class QCMeasure:
    _qc_headers_set = set(_qc_headers)
    def __init__(self, **fields):
        provided_keys = set(fields.keys())
        assert self._qc_headers_set.issuperset(provided_keys), (
                "One or more provided keys not in set of supported keys: "
                f"{_qc_headers}"
                )
        assert provided_keys.issuperset(_required_qc_headers), (
            f"Must initialise class with all of the following keys: {_required_qc_headers}"
            )
        for key in _qc_headers:
            if key in fields:
                setattr(self, key, fields[key])
            else:
                setattr(self, key, "")

    @classmethod
    def get_header(cls):
        return "\t".join(_qc_headers)

    def __repr__(self):
        values = starmap(getattr, zip(repeat(self),_qc_headers))
        return "\t".join(map(str,values))
