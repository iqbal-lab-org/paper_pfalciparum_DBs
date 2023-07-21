from itertools import starmap, repeat

class PositionalMetric:
    _headers = ["pos", "feature_ID", "metric_ID", "metric_value"]

    def __init__(self, **fields):
        assert set(self._headers) == set(fields.keys()), (
            f"Must initialise class with exactly the following keys: {self._headers}"
            )
        for key, val in fields.items():
            setattr(self, key, val)

    @classmethod
    def get_header(cls):
        return "\t".join(cls._headers) + "\n"

    def __repr__(self):
        values = starmap(getattr, zip(repeat(self),self._headers))
        return "\t".join(map(str,values)) + "\n"


class MetricsRecorder:
    """
    This derived from metrics.py in plasmo_surfants repo [TODO: deduplicate somehow]
    """
    _headers = list()
    _required_headers = list()
    _UNSET = "NA"

    def __init__(self, **fields):
        self._check_fields_are_supported(set(fields.keys()))
        for key in self._headers:
            if key in fields:
                setattr(self, key, fields[key])
            else:
                setattr(self, key, self._UNSET)

    def _check_fields_are_supported(self, field_names):
        assert set(self._headers).issuperset(field_names), (
                "One or more provided field names not in set of supported fields: "
                f"{self._headers}"
                )

    @classmethod
    def get_header(cls):
        return "\t".join(cls._headers) + "\n"

    def _get_missing_fields(self):
        result = set()
        for required_field in self._required_headers:
            if getattr(self, required_field) == self._UNSET:
                result.add(required_field)
        return result

    def update(self, **fields):
        self._check_fields_are_supported(set(fields.keys()))
        for key, val in fields.items():
            if getattr(self, key) != self._UNSET:
                raise ValueError(f"{key} is already set")
            setattr(self, key, val)

    def __getitem__(self, key):
        return getattr(self, key)

    def __repr__(self):
        missing_fields = self._get_missing_fields()
        assert len(missing_fields) == 0, (
                f"Missing unset fields: {missing_fields}. Class can only be serialised with all of the following fields set: {self._required_headers}"
            )
        values = starmap(getattr, zip(repeat(self),self._headers))
        return "\t".join(map(str,values)) + "\n"
