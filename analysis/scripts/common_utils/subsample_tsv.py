import sys
from pathlib import Path
from csv import DictReader
from collections import defaultdict


def usage():
    print(
        f"usage: {sys.argv[0]} input_tsv num_samples output_tsv [inclusion_tsv]\n"
        "inclusion_tsv lists samples to include in output_tsv\n"
        "Maximises sampling diversity by sampling from as many countries as possible."
    )
    exit(0)


def main():
    try:
        input_tsv = Path(sys.argv[1])
        num_samples = int(sys.argv[2])
        output_tsv = Path(sys.argv[3])
        inclusion_tsv = None
        if len(sys.argv) > 4:
            inclusion_tsv = Path(sys.argv[4])
    except IndexError:
        usage()

    included_samples = set()
    if inclusion_tsv is not None:
        with inclusion_tsv.open() as fin:
            reader = DictReader(fin, delimiter="\t")
            for row in reader:
                included_samples.add(row["Sample"])

    if len(included_samples) > num_samples:
        raise ValueError(
            f"{num_samples} requested for sampling, but already {len(included_samples)} in {input_tsv}"
        )

    selection_samples = defaultdict(list)
    output_samples = list()
    with input_tsv.open() as fin:
        reader = DictReader(fin, delimiter="\t")
        header = "\t".join(reader.fieldnames)
        for row in reader:
            if row["Exclusion reason"] == "Analysis_set":
                line = "\t".join(row.values())
                if row["Sample"] in included_samples:
                    output_samples.append(line)
                else:
                    selection_samples[row["Country"]].append(line)

    if len(output_samples) != len(included_samples):
        raise ValueError(
            f"Error: {len(included_samples)} listed for inclusion in {inclusion_tsv} but found {len(output_samples)} in {input_tsv}"
        )

    num_sampled = len(output_samples)
    num_remaining = sum(map(len, selection_samples.values()))
    if num_sampled + num_remaining < num_samples:
        raise ValueError(
            f"{num_sampled + num_remaining} samples can be sampled, but {num_samples} samples requested"
        )

    while num_sampled < num_samples:
        for key in selection_samples:
            elems = selection_samples[key]
            if len(elems) > 0:
                output_samples.append(elems[0])
                selection_samples[key] = elems[1:]
                num_sampled += 1
            if num_sampled >= num_samples:
                break

    with output_tsv.open("w") as fout:
        fout.write(header + "\n")
        for line in output_samples:
            fout.write(line + "\n")


if __name__ == "__main__":
    main()
