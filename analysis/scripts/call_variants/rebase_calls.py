"""
From a discovery VCF, the ref genome it refers to, and a genotyped VCF, express calls in discovery VCF against the reference of genotyped VCF.
(=rebasing in gramtools terminology)
"""

import sys
from pathlib import Path

from gramtools.commands.discover.discover import _rebase_vcf, _dump_rebased_vcf

def usage():
    print(f"usage: {sys.argv[0]} discov_vcf induced_ref.fa[.gz] genotyped.vcf[.gz] output_fname")
    exit(1)

class DiscoPaths:
    """
    Adapter to the class used in gramtools
    """
    def __init__(self, discov_vcf, genotype_vcf, induced_ref, output_vcf):
        self.discov_vcf = discov_vcf
        self.geno_vcf = genotype_vcf
        self.pers_ref = induced_ref
        self.final_vcf = output_vcf

if __name__ == "__main__":
    if len(sys.argv) != 5:
        usage()
    discov_vcf = Path(sys.argv[1]).resolve()
    induced_ref = Path(sys.argv[2]).resolve()
    genotype_vcf = Path(sys.argv[3]).resolve()
    output_vcf = Path(sys.argv[4]).resolve()

    disco_paths = DiscoPaths(discov_vcf, genotype_vcf, induced_ref, output_vcf)
    rebased_vcf_records = _rebase_vcf(disco_paths, check_records=False)
    _dump_rebased_vcf(rebased_vcf_records, disco_paths)

