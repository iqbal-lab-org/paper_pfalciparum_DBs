from itertools import repeat, starmap

_qc_headers = ["sample", "chrom", "start", "end", "gene","tool","metric_category","metric","value"]
_required_qc_headers = set(_qc_headers[0:4] + _qc_headers[6:9])

class QCMeasure:
    def __init__(self, **fields):
        self.check_fields_are_supported(set(fields.keys()))
        for key in _qc_headers:
            if key in fields:
                setattr(self, key, fields[key])
            else:
                setattr(self, key, "")

    def check_fields_are_supported(self, field_names):
        assert set(_qc_headers).issuperset(field_names), (
                "One or more provided field names not in set of supported fields: "
                f"{_qc_headers}"
                )

    @classmethod
    def get_header(cls):
        return "\t".join(_qc_headers)

    def get_missing_fields(self):
        result = set()
        for required_field in _required_qc_headers:
            if getattr(self, required_field) == "":
                result.add(required_field)
        return result

    def update(self, **fields):
        self.check_fields_are_supported(set(fields.keys()))
        for key, val in fields.items():
            if getattr(self, key) != "":
                raise ValueError(f"{key} is already set")
            setattr(self, key, val)

    def __repr__(self):
        missing_fields = self.get_missing_fields()
        assert len(missing_fields) == 0, (
                f"Missing unset fields: {missing_fields}. Class can only be serialised with all of the following fields set: {_required_qc_headers}"
            )
        values = starmap(getattr, zip(repeat(self),_qc_headers))
        return "\t".join(map(str,values))
