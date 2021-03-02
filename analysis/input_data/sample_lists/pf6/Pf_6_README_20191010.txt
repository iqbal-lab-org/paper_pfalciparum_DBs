======================================================================
MalariaGEN Plasmodium falciparum Community Project - Pf 6 data release
======================================================================

Date: 2019-10-10

======================================================================
Description
======================================================================

Through an analysis of 7,113 parasite samples collected at 73 different locations in Africa, Asia, America and Oceania, we identified more than 6 million variant poisitons (SNPs and indels). This download includes genotyping data for samples contributed to the MalariaGEN Plasmodium falciparum Community Project.

Potential data users are asked to respect the legitimate interest of the Community Project and its partners by abiding any restrictions on the use of a data as described in the Terms of Use: http://www.malariagen.net/projects/parasite/pf/use-p-falciparum-community-project-data

For more information on the P. falciparum Community Project that generated these data, please visit: http://www.malariagen.net/projects/parasite/pf


The methods used to generate the data are described in detail in MalariaGEN Plasmodium falciparum Community Project, biorxiv (2019), DOI: ???.


======================================================================
Citation information
======================================================================

Publications using these data should acknowledge and cite the source of the data using the following format: "This publication uses data from the MalariaGEN Plasmodium falciparum Community Project as described in An open dataset of Plasmodium falciparum genome variation in 7,000 worldwide samples, biorxiv, 2019 (DOI: ???)."



======================================================================
Files in the release
======================================================================

Files in the release include:

    - Pf_6_samples.txt : sample metadata file in tab-delimited format
    - Pf_6_drug_resistance_marker_genotypes.txt : drug resistance marker genotypes file in tab-delimited format
    - Pf_6_inferred_resistance_status_classification.txt : inferred resistance status classification file in tab-delimited format
    - Pf_6_fws.txt : Fws values in tab-delimited format
    - Pf_6_mapping_genetic_markers_to_inferred_resistance_status_classification.docx : rules used to create Pf_6_inferred_resistance_status_classification.txt in Word format
    - Pf_6_studies.xlsx : details on partner studies in Excel format
    - Pf_6_genes_data.txt : genes file in tab-delimited format
    - Pf_6_vcf : directory containing vcf files, one per chromosome
    - Pf_6.zarr.zip : sample genotypes in zipped zarr format



File descriptions:

======================================================================
- Pf_6_samples.txt
======================================================================

This file includes sample metadata for all 7,113 samples collected from partners and details of sequence read data available at the European Nucleotide Archive (ENA, https://www.ebi.ac.uk/ena). It contains the following columns:

    - Sample                         Unique ID of each sample (which can be used to link to other sample information and genotypes)
    - Study                          Code of the partner study which collected this sample
    - Site                           Location the sample was collected (in the case of returning travellers this is the country visited)
    - Country                        Country in which the sample was collected (in the case of returning travellers this is the country visited)
    - Lat                            GPS coordinate of the latitude of the Site
    - Long                           GPS coordinate of the longitude of the Site
    - Year                           Year in which the sample was collected
    - ENA                            ENA run accession(s) for the sequencing read data. In some cases multiple runs of sequencing data were merged
    - All samples same individual    Identifies all the sample in the data set which were collected from the same individual, either at the same of different time points
    - Is returning traveller         Flag indicating whether the sample was from a returning traveller (True) or not (False)
    - Population                     Population to which we assigned the sample (SAM=South America, WAF=West Africa, CAF=Central Africa, EAF=East African, SAS=South Asia, WSEA=Western SE Asia, ESEA=Eastern SE Asia, OCE=Oceania)
    - % callable                     Percentage of the genome which has coverage of at least 5 reads and less than 10% of reads with mapping quality 0
    - QC pass                        Flag indicated whether the sample passed QC (True=passed QC, False=failed QC)
    - Exclusion reason               Reason samples failed QC


======================================================================
- Pf_6_drug_resistance_marker_genotypes.txt
======================================================================

This file contains genotypes at drug resistance markers for all 7,113 samples derived from analysis of sequence data. It contains the following columns:

    - Sample                         Unique ID of each sample (which can be used to link to other sample information and genotypes)
    - crt_76[K]                      Amino acid at crt position 76. For explanation see below.
    - crt_72-76[CVMNK]               Amino acids at crt positions 72 to 76. For explanation see below.
    - dhfr_51[N]                     Amino acid at dhfr position 51. For explanation see below.
    - dhfr_59[C]                     Amino acid at dhfr position 59. For explanation see below.
    - dhfr_108[S]                    Amino acid at dhfr position 108. For explanation see below.
    - dhfr_164[I]                    Amino acid at dhfr position 164. For explanation see below.
    - dhps_437[G]                    Amino acid at dhps position 437. For explanation see below.
    - dhps_540[K]                    Amino acid at dhps position 540. For explanation see below.
    - dhps_581[A]                    Amino acid at dhps position 581. For explanation see below.
    - dhps_613[A]                    Amino acid at dhps position 613. For explanation see below.
    - k13_class                      Genotype class of kelch13 gene (WT=wild type, MU=homozygous mutant, HE=heterozygous between wild type and mutant, MH=heterozygous between two different mutants, MI=missing call)
    - k13_allele                     All the mutations found in the "KBPD" region of kelch13 (amino acids 349-726). Homozygous mutations are shown in upper case and heterozygous in lower case. If multiple non-synonympus mutations are found, these will be separated by a comma ','. If all SNPs are hom ref, genotype will be empty string. If any SNP is missing, genotype will be '-'.
    - cn_mdr1                        Copy number of mdr1. Heterozygous duplications will have a copy number of 1.5. If copy number is undetermined will be shown as -1.0
    - cn_pm2                         Copy number of plasmepsin 2-3. Heterozygous duplications will have a copy number of 1.5. If copy number is undetermined will be shown as -1.0
    - cn_gch1                        Copy number of gch1. Heterozygous duplications will have a copy number of 1.5. If copy number is undetermined will be shown as -1.0
    - dup_mdr1                       1.0=mdr1 duplicated, 0.0=mdr1 not duplicated, -1.0=duplication status of mdr1 undetermined
    - dup_pm2                        1.0=plasmepsin 2-3 duplicated, 0.0=plasmepsin 2-3 not duplicated, -1.0=duplication status of plasmepsin 2-3 undetermined
    - dup_gch1                       1.0=gch1 duplicated, 0.0=gch1 not duplicated, -1.0=duplication status of gch1 undetermined
    - breakpoint_mdr1                Tandem duplication breakpoints around mdr1. Refer to Supplementary Table 9 for details.
    - breakpoint_pm2                 Tandem duplication breakpoints around plasmepsin 2-3. Refer to Supplementary Table 10 for details.
    - breakpoint_gch1                Tandem duplication breakpoints around gch1. Refer to Supplementary Table 8 for details.

Explanation of amino acid columns in crt, dhfr and dhps:

Each value can have a single haplotype if homozygous or two haplotypes separated by a comma if heterozygous
It is possible to have heterozygous calls where both amino acid haplotypes are the same. The heterozygosity here is at the nucleotide level. These could perhaps be considered homozygous alt.
- represents missing (missing genotype in at least one of the positions)
X represents an unphased het followed by another het. Because hets are unphased it is not possible to resolve the two haplotypes. These are perhaps best considered missing.
? represents two different values in the GATK PID field. Again it is not possible to resolve the two haplotypes in such cases. These are perhaps best considered missing.
! represents a frame-shift in the haplotype. These are perhaps best considered missing.


======================================================================
- Pf_6_inferred_resistance_status_classification.txt
======================================================================

This file includes sample phenotype data for 5,970 samples that passed QC derived from the data in Pf_6_drug_resistance_marker_genotypes.txt, using the rules outlined in "Pf6 mapping genetic markers to inferred resistance status classification.docx", together with deletion genotypes for HRP2 and HRP3 that can be used to determine resitance to parid diagnostic tests (RDTs). It contains the following columns:

    - Chloroquine                    Chloroquine resistance status. Resistant/Sensitive/Undetermined
    - Pyrimethamine                  Pyrimethamine resistance status. Resistant/Sensitive/Undetermined
    - Sulfadoxine                    Sulfadoxine resistance status. Resistant/Sensitive/Undetermined
    - Mefloquine                     Mefloquine resistance status. Resistant/Sensitive/Undetermined
    - Artemisinin                    Artemisinin resistance status. Resistant/Sensitive/Undetermined
    - Piperaquine                    Piperaquine resistance status. Resistant/Sensitive/Undetermined
    - SP (uncomplicated)             Sulfadoxine-Pyrimethamine treatment resistance status. Samples carrying the dhfr triple mutant, which is strongly associated with SP failure. Resistant/Sensitive/Undetermined
    - SP (IPTp)                      Sulfadoxine-Pyrimethamine intermittent preventive treatment in pregnancy resistance status. Samples carrying the dhfr/dhps sextuple mutant, which confers a higher level of SP resistance. Resistant/Sensitive/Undetermined
    - AS-MQ                          Artesunate-mefloquine resistance status. Resistant/Sensitive/Undetermined
    - DHA-PPQ                        Dihydroartemisinin-piperaquine resistance status. Resistant/Sensitive/Undetermined
    - HRP2                           Deletions at HRP2 associated with failure of rapid diagnostic tests. del=HRP2 deleted, nodel=HRP2 not deleted
    - HRP3                           Deletions at HRP3 associated with failure of rapid diagnostic tests. del=HRP3 deleted, nodel=HRP3 not deleted
    - HRP2 and HRP3                  Deletions at HRP2 and HRP3 associated with failure of rapid diagnostic tests. del=both HRP2 and HRP3 deleted, nodel=either HRP2, HRP3 or both not deleted


======================================================================
- Pf_6_fws.txt
======================================================================

This file includes Fws values for 5,970 samples that passed QC. It contains the following columns:

    - Sample                         Unique ID of each sample (which can be used to link to other sample information and genotypes)
    - Fws                            Fws value


======================================================================
- Pf_6_mapping_genetic_markers_to_inferred_resistance_status_classification.docx
======================================================================

This file describes rules used to create Pf_6_inferred_resistance_status_classification.txt


======================================================================
- Pf_6_genes_data.txt
======================================================================

This file contains gene information including global and local differentiation scores derived from analysis of SNP differentiation. It contains the following columns:

    - gene_id                        GeneDB ID (see https://www.genedb.org/)
    - gene_name                      GeneDB gene name
    - chrom                          Chromsome on which the gene lies in the reference genome
    - start                          Start position (bp) of the gene in the reference genome
    - end                            End position (bp) of the gene in the reference genome
    - global_differentiation_score   Differentiation score derived from analysis of Fst between populations
    - local_differentiation_score    Differentiation score derived from analysis of Fst within populations
    - distance_to_higher_local_diff_score    Distance (in bp) to a gene with a higher local_differentiation_score


======================================================================
- Pf_6_vcf
======================================================================

This directory contains vcf files, one per chromosome. Each file is in bgzip format (.vcf.gz) and has an associated tabix index file (.vcf.gz.tbi) and MD5 checksum (.vcf.gz.md5). There are sixteen files in total, fourteen for each of the autosomes (Pf3D7_01_v3 - Pf3D7_14_v3), one for the mitochondrial sequence (Pf_M76611) and one for the apicoplast sequence (Pf3D7_API_v3).

The files, once unzipped, are tab-separated text files, but may be too large to open in Excel.

The VCF format is described in https://github.com/samtools/hts-specs

Tools to assist in handling VCF files are freely available from
http://samtools.github.io/bcftools/

The VCF files contains details of 6,051,696 discovered variant genome positions.
These variants were discovered amongst all samples from the release.
3,168,721 of these variant positions are SNPs, with the remainder being either
short insertion/deletions (indels), or a combination of SNPs and indels. It is
important to note that many of these variants are considered low quality. Only
the variants for which the FILTER column is set to PASS should be considered of
reasonable quality. There are 3,110,990 such PASS variants of which 1,835,719
are SNPs and 1,275,271 indels.

The FILTER column is based on two types of information. Firstly certain regions
of the genome are considered "non-core". This includes sub-telomeric regions,
centromeres and internal VAR gene regions on chromosomes 4, 6, 7, 8 and 12. All
variants within non-core regions are considered to be low quality, and hence
will not have the FILTER column set to PASS. The regions which are core and
non-core can be found in the file
ftp://ngs.sanger.ac.uk/production/malaria/pf-crosses/1.0/regions-20130225.onebased.txt.

Secondly, variants are filtered out based on a quality score called VQSLOD. All
variants with a VQSLOD score below 0 are filtered out, i.e. will have a value of
Low_VQSLOD in the FILTER column, rather than PASS. The VQSLOD score for each
variant can be found in the INFO field of the VCF file. It is possible to use
the VQSLOD score to define a more or less stringent set of variants (see next
section for further details).

It is also important to note that many variants have more than two alleles. For
example, amongst the 3,110,990 PASS variants, 1,782,849 are biallelic. The
remaining 1,328,141 PASS variants have 3 or more alleles. The maximum number of
alternative alleles represented is 6. Note that some positions can in truth have
more than 6 alternative alleles, particularly those at the start of short tandem
repeats. In such cases, some true alternative alleles will be missing.

In addition to alleles representing SNPs and indels, some variants have an
alternative allele denoted by the * symbol. This is used to denote a "spanning
deletion". For samples that have this allele, the base at this position has been
deleted. Note that this is not the same as a missing call - the * denotes that
there are reads spanning across this position, but that the reads have this
position deleted yet map on either side of the deletion. For further details see
https://software.broadinstitute.org/gatk/guide/article?id=6926

In addition to the VQSLOD score mentioned above, The INFO field contains many
other variant-level metrics. The metrics QD, FS, SOR, DP are all measures
related to the quality of the variant. The VQSLOD score is derived from these
five metrics.

AC contains the number of non-reference alleles amongst the samples in the file.
Because the file contains diploid genotype calls, homozygous non-reference calls
will be counted as two non-reference alleles, whereas heterozygous calls will be
counted as one non-reference allele. Where a variant position has more than one
one non-reference allele, counts of each different non-reference allele are
given. AN contains the total number of called alleles, including reference
alleles. A simple non-reference allele frequency can be calculated as AC/AN.
AC and AN values are all specfic to the samples in the study the VCF was created
for.

Various functional annotations are held in the the SNPEFF variables of the INFO
field. Where appropriate, the amino acid change caused by the variant can be
found in SNPEFF_AMINO_ACID_CHANGE. Note that for multi-allelic variants, only
one annotation is given, and therefore this should not be relied on for non-
biallelic variants. SNPEFF_AMINO_ACID_CHANGE also does not take account of
nearby variants, so if two SNPs are present in the same codon, the
amino acid change given is likely to be wrong. Similarly, if two coding indels
are found in the same exon, the SNPEFF annotations are likely to be wrong. This
situation occurs at the CRT locus (see next section for further details).

Coding variants are identified using the CDS flag in the INFO field.

Columns 10 and onwards of the VCF contain the information for each sample.
The first component of this (GT) is always the diploid genotype call as
determined by GATK. A value of 0/0 indicates a homozygous reference call. A
value of 1/1 indicates a homozygous alternative allele call. 0/1 indicates a
heterozygous call. A value of 2 indicates the sample has the second alternative
allele, i.e. the second value in the ALT column. For example 2/2 would mean the
sample is homozygous the the second alternative allele, 0/2 would mean the
sample is heterozygous for the reference and second alternative alleles, and 1/2
would mean the sample is heterozygous for the first and second alternative
alleles. A value of ./. indicates a missing genotype call, usually because there
are no reads mapping at this position in that sample.


Recommendations regarding sets of variants to use in analyses
-------------------------------------------------------------

Variants are filtered using the VQSLOD metric. VQSLOD is log(p/q) where p is the
probability of being true and q is the probability of being false. Theoretically,
when VQSLOD > 0, p is greater than q, and therefore the variant is more likely
true than false. Conversely, when VQSLOD < 0, the variant is theoretically more
like false than true. This is why we have chosen 0 as the threshold to use to
declare that variants have passed the filters: all PASS variants are
theoretically more likely true than false. Of course, for variants where VQSLOD
is only slightly above 0, there is only a slightly greater probability of being
true than of being false. Therefore, for example, many of the variants with
values between 0 and 1 are likely to be false.

Empirically we have found that SNPs tend to be more accurate than indels, coding
variants tend to be more accurate than non-coding variants, and bi-allelic
variants tend to be more accurate than multi-allelic variants. If you require a
very reliable set of variants for genome-wide analysis, and don't mind if you
miss some real variants, we recommend using only bi-allelic coding SNPs in the
core genome with a VQSLOD score > 6. There are 83,168 such stringent SNPs in the
call set. We include a command below to create such a set of variants.

If instead you would like to know of all likely variation within a certain
region, even if this means including a few false variants, we recommend using
all PASS variants. Finally, if you want to ensure you miss as little as possible
of the true variation, at the risk of including large numbers of false positives,
you could ignore the FILTER column and use all variants in the VCF.

In general, we recommend caution in analysing indels. For any given sample, the
majority of differences from the reference genome are likely to be due to indels
in low-complexity non-coding regions, e.g. in length polymorphisms of short
tandem repeats (STRs), such as homopolymer runs or AT repeats. In general, it is
difficult to map short reads reliably in such regions, and this is compounded by
the fact that these regions tend to have high AT content, and in general we
typically have much lower coverage in high AT regions. Indels also tend to be
multi-allelic, making analysis much more challenging than for (typically
bi-allelic) SNPs.

Despite what is written above, it may often be important to analyse indels in
order to determine the true nature of variation at a locus. An example of this
is analysis of haplotypes of codons 72-76 of the chloroquine resistance
transporter gene PfCRT. These codons are translated from bases 403,612-403,626
on chromosome 7 (Pf3D7_07_v3:403612-403626). The vcf contains two indels here,
an insertion of a T after position 403,618, and a deletion of T after position
403,622. However, it turns out that every sample that has the insertion after
403,618 also has the deletion after 403,622. If we write out the full sequence
of the 3D7 reference (Ref), and a sample which has both the insertion and
deletion (Alt), we see the following:
Ref: TGTGTAAT-GAATAAA = TGTGTAATGAATAAA (amino acid sequence CVMNK)
Alt: TGTGTAATTGAA-AAA = TGTGTAATTGAAAAA (amino acid sequence CVIEK)
As can be seen from the above, the single base insertion and deletion are
equivalent to three SNPs at positions 403,620 (G/T), 403,621 (A/G) and
403,623 (T/A). If we had only chosen to analyse SNPs at this locus, we would
not have seen the CVIEK haplotype.


Extracting data from the VCF file
-----------------------------

We recommend the use of bcftools. To install bcftools, follow the instructions
at: https://github.com/samtools/bcftools/wiki/HOWTOs

The following are some commands which you might find useful for extracting data
from the vcf.gz files. We've used an example the vcf for chromosome 5
(Pf_60_public_Pf3D7_05_v3.final.vcf.gz), but similar commands should work on all vcf files.

To extract sample IDs and put into a file, one per line:
bcftools query --list-samples Pf_60_public_Pf3D7_05_v3.final.vcf.gz > samples.txt

To extract chromosome, position, reference allele, all alternate alleles,
filter value and VQSLOD for all variants into a tab-delimited file:
bcftools query -f \
'%CHROM\t%POS\t%REF\t%ALT{0}\t%ALT{1}\t%ALT{2}\t%ALT{3}\t%ALT{4}\t%ALT{5}\t%FILTER\t%VQSLOD\n' \
Pf_60_public_Pf3D7_05_v3.final.vcf.gz > all_variants.txt

To extract chromosome, position, reference allele, all alternate alleles and
VQSLOD for PASS SNPs only into a tab-delimited file:
bcftools query -f \
'%CHROM\t%POS\t%REF\t%ALT{0}\t%ALT{1}\t%ALT{2}\t%ALT{3}\t%ALT{4}\t%ALT{5}\t%VQSLOD\n' \
--include 'FILTER="PASS" && TYPE="snp"' \
Pf_60_public_Pf3D7_05_v3.final.vcf.gz > pass_snps.txt

To extract chromosome, position, reference allele, alternate allele and VQSLOD
for biallelic PASS SNPs only into a tab-delimited file:
bcftools query -f \
'%CHROM\t%POS\t%REF\t%ALT{0}\t%VQSLOD\n' \
--include 'FILTER="PASS" && TYPE="snp" && N_ALT=1' \
Pf_60_public_Pf3D7_05_v3.final.vcf.gz > biallelic_pass_snps.txt

To extract chromosome, position, reference allele, alternate allele and VQSLOD
for biallelic PASS SNPs that are segregating within the study into a
tab-delimited file:
bcftools query -f \
'%CHROM\t%POS\t%REF\t%ALT{0}\t%VQSLOD\n' \
--include 'FILTER="PASS" && TYPE="snp" && N_ALT=1 && AC>0' \
Pf_60_public_Pf3D7_05_v3.final.vcf.gz > biallelic_segregating_pass_snps.txt

To create a vcf file which contains only PASS bi-allelic coding SNPs with
VQSLOD > 6:
bcftools view \
--include 'FILTER="PASS" && N_ALT=1 && CDS==1 && TYPE="snp" && VQSLOD>6.0' \
--output-type z \
--output-file output_filename.vcf.gz \
Pf_60_public_Pf3D7_05_v3.final.vcf.gz
bcftools index --tbi output_filename.vcf.gz

To extract diploid genotype calls for biallelic PASS SNPs in gene MDR1 into a
tab-delimited text file, including the chromosome, position, ref and alt
alleles, VQSLOD score and amino acid substituion, and a header containing
sample names:
bcftools query \
-f '%CHROM\t%POS\t%REF\t%ALT{0}\t%VQSLOD\t%SNPEFF_AMINO_ACID_CHANGE[\t%GT]\n' \
--regions Pf3D7_05_v3:957890-962149 \
--include 'FILTER="PASS" && TYPE="snp" && N_ALT=1' \
--print-header \
Pf_60_public_Pf3D7_05_v3.final.vcf.gz > mdr1_genotypes.txt

To extract ref allele depths for biallelic PASS SNPs in gene MDR1 into a
tab-delimited text file, including the chromosome, position, ref and alt
alleles, VQSLOD score and amino acid substituion, and a header containing
sample names:
bcftools query \
-f '%CHROM\t%POS\t%REF\t%ALT{0}\t%VQSLOD\t%SNPEFF_AMINO_ACID_CHANGE[\t%AD{0}]\n' \
--regions Pf3D7_05_v3:957890-962149 \
--include 'FILTER="PASS" && TYPE="snp" && N_ALT=1' \
--print-header \
Pf_60_public_Pf3D7_05_v3.final.vcf.gz > mdr1_ref_allele_depth.txt

To extract alt allele depths for biallelic PASS SNPs in gene MDR1 into a
tab-delimited text file, including the chromosome, position, ref and alt
alleles, VQSLOD score and amino acid substituion, and a header containing
sample names:
bcftools query \
-f '%CHROM\t%POS\t%REF\t%ALT{0}\t%VQSLOD\t%SNPEFF_AMINO_ACID_CHANGE[\t%AD{1}]\n' \
--regions Pf3D7_05_v3:957890-962149 \
--include 'FILTER="PASS" && TYPE="snp" && N_ALT=1' \
--print-header \
Pf_60_public_Pf3D7_05_v3.final.vcf.gz > mdr1_alt_allele_depth.txt


======================================================================
    - Pf_6.zarr.zip
======================================================================

This file contains the information that is encoded in the VCF files, but in zipped zarr format.

We recommend analysing data using the scikit-allel package with the zarr file. For more details
on using scikit-allel, please see https://scikit-allel.readthedocs.io/en/stable/


======================================================================
Release notes:
======================================================================

Data excluded from release:
Sequence read data on samples collected in Indonesia cannot be made publically available because of national export restrictions.

