import sys
from pathlib import Path

import cortex.calls as cortex


if len(sys.argv) != 4:
    print(f"Usage: {sys.argv[0]} ref.fa[.gz] reads1.fq[.gz] [reads_n.fq[.gz]] out.vcf")
    exit(1)

fasta_ref = Path(sys.argv[1]).resolve()
reads_files = sys.argv[2].split(" ")
out_vcf = Path(sys.argv[3])

cortex.run(fasta_ref, reads_files, out_vcf, mem_height=26)
