import sys
from pathlib import Path

import cortex.calls as cortex


if not (4 <= len(sys.argv) <= 5):
    print(
        f"Usage: {sys.argv[0]} ref.fa[.gz] reads1.fq[.gz] [reads_n.fq[.gz]] out.vcf [sample_id]"
    )
    exit(0)

fasta_ref = Path(sys.argv[1]).resolve()
reads_files = sys.argv[2].split(" ")
out_vcf = Path(sys.argv[3])
try:
    sample_name = sys.argv[4]
except IndexError:
    sample_name = "sample"

cortex.run(
    fasta_ref,
    reads_files,
    out_vcf,
    sample_name=sample_name,
    mem_height=22,
    tmp_directory=out_vcf.parent,
    cleanup=False,
)
