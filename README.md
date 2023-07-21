# In-depth analysis of sequence variation in *P. falciparum* genes DBLMSP and DBLMSP2 

author: Brice Letcher

## Summary

This repo contains the analysis of highly-variable, paralogous *P. falciparum* genes
DBLMSP (abbrev: DB1) and DBLMSP2 (abbrev: DB2) (collective abbrev: DBs)

It is a 'sister project' of project `plasmo_surfants`, used for genotyping
highly-variable *P. falciparum* genes. It is structured in the same way: snakemake
workflows + scripts under `analysis`.

Like its sister project, it supports the following publication:

> doi: https://doi.org/10.1101/2023.02.27.530215

Inputs: Input data is copied over from the outputs of `plasmo_surfants`, see below.

**Data availability**: The analysed sequences are all available on zenodo (see paper).

Documentation: I made a literature review of DBLMSP and DBLMSP2, [here](docs/dblmsp_lit_review/report.pdf).

## Workflows

### Requirements for running

1. A symlink needs to exist at the root of the project pointing to the repo `plasmo_surfants`
2. The appropriate input files, generated by `plasmo_surfants` repo, must be copied into
`analysis/input_data`.

For point 2, these are:

`$ find analysis/input_data`
```sh
analysis/input_data/
analysis/input_data/analysed_sequences
analysis/input_data/analysed_sequences/DBs_single_snp_solved_seqs.fa
analysis/input_data/analysed_sequences/assembly_sequences
analysis/input_data/analysed_sequences/assembly_sequences/dna
analysis/input_data/analysed_sequences/assembly_sequences/dna/AMA1_pf_assemblies.fa
analysis/input_data/analysed_sequences/assembly_sequences/dna/DBs_gram_joint_geno_ebf7bcd5__pacb_ilmn_pf@pf6__analysis_set_fws95__gapfiller__7__13.fa
analysis/input_data/analysed_sequences/assembly_sequences/dna/DBs_pf_assemblies.fa
analysis/input_data/analysed_sequences/assembly_sequences/dna/EBA175_pf_assemblies.fa
analysis/input_data/analysed_sequences/assembly_sequences/dna/MSP1_pf_assemblies.fa
analysis/input_data/analysed_sequences/assembly_sequences/dna/MSP2_pf_assemblies.fa
analysis/input_data/analysed_sequences/assembly_sequences/dna/MSP3_pf_assemblies.fa
analysis/input_data/analysed_sequences/assembly_sequences/dna/MSP6_pf_assemblies.fa
analysis/input_data/analysed_sequences/assembly_sequences/protein
analysis/input_data/analysed_sequences/assembly_sequences/protein/DBs_laverania_assemblies.fa
analysis/input_data/analysed_sequences/assembly_sequences/protein/DBs_laverania_assemblies_with_DBLMSP2_reichenowi.fa
analysis/input_data/analysed_sequences/ir_stats_all.tsv
analysis/input_data/analysed_sequences/ir_stats_fold_coverages.tsv
analysis/input_data/analysed_sequences/pf6
analysis/input_data/analysed_sequences/pf6/AMA1_gramtools_joint_geno_ebf7bcd5__pf6__analysis_set_fws95__gapfiller__7__13.fa
analysis/input_data/analysed_sequences/pf6/AMA1_malariaGEN_pf6.fa
analysis/input_data/analysed_sequences/pf6/DBLMSP2_gramtools_joint_geno_ebf7bcd5__pf6__analysis_set_fws95__gapfiller__7__13.fa
analysis/input_data/analysed_sequences/pf6/DBLMSP_gramtools_joint_geno_ebf7bcd5__pf6__analysis_set_fws95__gapfiller__7__13.fa
analysis/input_data/analysed_sequences/pf6/DBs_gramtools_joint_geno_ebf7bcd5__pf6__analysis_set_fws95__gapfiller__7__13.fa
analysis/input_data/analysed_sequences/pf6/DBs_malariaGEN_pf6.fa
analysis/input_data/analysed_sequences/pf6/EBA175_gramtools_joint_geno_ebf7bcd5__pf6__analysis_set_fws95__gapfiller__7__13.fa
analysis/input_data/analysed_sequences/pf6/EBA175_malariaGEN_pf6.fa
analysis/input_data/analysed_sequences/pf6/MSP1_gramtools_joint_geno_ebf7bcd5__pf6__analysis_set_fws95__gapfiller__7__13.fa
analysis/input_data/analysed_sequences/pf6/MSP1_malariaGEN_pf6.fa
analysis/input_data/analysed_sequences/pf6/MSP2_gramtools_joint_geno_ebf7bcd5__pf6__analysis_set_fws95__gapfiller__7__13.fa
analysis/input_data/analysed_sequences/pf6/MSP2_malariaGEN_pf6.fa
analysis/input_data/analysed_sequences/pf6/MSP3_gramtools_joint_geno_ebf7bcd5__pf6__analysis_set_fws95__gapfiller__7__13.fa
analysis/input_data/analysed_sequences/pf6/MSP3_malariaGEN_pf6.fa
analysis/input_data/analysed_sequences/pf6/MSP6_gramtools_joint_geno_ebf7bcd5__pf6__analysis_set_fws95__gapfiller__7__13.fa
analysis/input_data/analysed_sequences/pf6/MSP6_malariaGEN_pf6.fa
analysis/input_data/domain_coordinates.tsv
analysis/input_data/laverania_assemblies.tsv
analysis/input_data/sequence_colour_schemes.tsv
```

These files are copied over from the outputs of workflows in the `plasmo_surfants` repository.
They are available on zenodo.

### List and description of workflows

#### `analysed_sequences`

Filters the sequences described above based on how well-resolved they are (see paper for
definition). Translates those sequences to protein, makes MSAs, and clusters sequences.

#### `seq_stats`

Performs all the analyses shown in the paper: computing heterozygosity, defining
DBL-spanning domain, defining shared sequences, mosaic alignments, detecting gene
conversion.

This workflow relies on a common SQLite database that it creates and fills in, and is
later used for plotting figures.

## Miscellaneous

### DBL Coordinates

These are stored in `analysis/input_data/analysed_sequences/domain_coordinates.tsv`, 1-based, in both protein and dna space, and on both the 3D7 sequence and the MSA of DBs. I computed two sets of coordinates, as follows.

#### From known domain model
DBL domain from PFAM: [PF05424](http://pfam.xfam.org/family/PF05424) was downloaded and queried on DB1 and DB2 protein sequences from 3D7 using `hmmscan` from `hmmer` suite.

Best hits are:

```
DB1
>> Duffy_binding  Duffy binding domain
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   75.7   0.9   2.5e-25   2.5e-25       1     156 [.     185     310 ..     185     330 .. 0.82

DB2

>> Duffy_binding  Duffy binding domain
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   2 !   68.9   2.5   3.2e-23   3.2e-23       1     169 [.     201     341 ..     201     378 .. 0.77
```

Clear, single domain hit, at (1-based) coordinates: 185:301 (DB1) and 201:341 (DB2).

#### Operational definition

I also defined coordinates for DBL in DBs by taking the window inside which shared sequences exists (defined as the largest region in DBs alignment with no fully diverged nucleotides/AAs).
