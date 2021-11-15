# Rationale

This project genotypes malariaGEN samples using gramtools.

# Workflows

## Constraints

In each workflow I load a global configfile ("common.yaml") and a set of common python utilities ("common_utils.py"). "common.yaml" must be loaded before "common_utils.py", and "common_utils.py" must be loaded before other "utils.py".

## Conventions

Inside workflows, functions prefixed with `cu_` come from `analysis/workflows/common_utils.py`

Similarly, each workflow's own `analysis/workflows/<workflow_name>/utils.py` defines functions prefixed with a two-letter code: e.g. for `download_data`, functions will start with `dd_`.


## Running a workflow

To run a workflow on the cluster, call
```
sh analysis/cluster_submit.sh
```

## List and description of workflows

### download_data

#### Requirements

None

#### Function
3 main datasets:

  * pf6: downloads all p. falciparum read sets from malariaGEN
  * pvgv: downloads all p. vivax read sets from malariaGEN
  * pacb_ilmn_pf: downloads paired illumina reads and pacbio assemblies for 15 samples


### call_variants

#### Requirements

download_data

#### Function

2 main operations:

   * Runs cortex on all pf6 and pvgv samples
   * Runs samtools, cortex and gramtools (adjudicating samtools and cortex calls) on pacb_ilmn_pf samples. This allows `eval_varcalls` workflow to evaluate the calls for all three tools.

### eval_varcalls

#### Requirements

call_variants

#### Function


### make_prgs

#### Requirements 

call_variants

#### Function

* Makes a prg based on cortex calls in pf6 samples. Configurable parameters are:
       - Which pf6 samples to use for graph construction
       - Which genes to build the prg on
       - What `min_match_len` to use for `make_prg`
* Makes a prg based on cortex calls in pvgv samples.

### joint_genotyping

#### Requirements

make_prg

#### Function

Takes as input a genome graph made by make_prg, and runs gramtools genotyping on all specified samples.

### plasmo_paralogs

#### Requirements

joint_genotyping

#### Function

Produces sequences of the paralogs to study, makes MSAs, translates to protein, makes clusters using cd-hit.


## Development

### Testing without built container
To avoid rebuilding a container everytime a different version of the gramtools codebase is tested, it is installed locally:

```
. venv/bin/activate
pip install -e <path/to/gramtools>
```

This means the file pyrequirements.txt should not be produced by running `pip freeze > pyrequirements.txt`.


# Input data

We use Plasmodium data from the MalariaGEN projects: https://www.malariagen.net/ 
The following Sanger ftp ftp://ngs.sanger.ac.uk/production/malaria/ lists:
  - pf6 metadata/vcfs (under pfcommunityproject/Pf6)
  - pvivax vcfs (under pvgv)

## Reference genomes

For *P. falciparum*, I have used the 2018-11 version for variant calling and graph building (ftp://ftp.sanger.ac.uk/pub/project/pathogens/gff3/2018-11/). pf6 used the 2016-07 version in their release VCFs (see pf6 paper). However, I checked the following:

* The two genomes have exactly the same chromosome sizes
* The two annotation files (GFF) have exactly the same coordinates for pf6_26_genes

On that basis it is fine to use 2018-11 version for evaluation purposes (see eval_varcalls workflow)

## pf6

Some samples are duplicates with low coverage, estimated to be non-clonal (via Fws), 'continent mismatches'
(labeled from a continent but their genotypic data maps them to another); these are marked as not 'Analysis_set' in the
'Exclusion reason' column. 

We build graphs from and genotype only 'Analysis_set' samples.


# Genes

## P. falciparum

See the [documentation](docs/pf_genes.pdf) for a list and rationale for the genes we chose to analyse.

## Barry lab

Use GATK best practices pipeline with slight modifications [preprint link][myo_1]

### Myo's project

i) Characterising the global diversity of P. falciparum surface antigens. [preprint link][myo_1]
ii) Longitudinal study of children with asymptomatic/symptomatic malaria,
and spotting allelic switches favouring symptomatic. [preprint link](https://www.medrxiv.org/content/10.1101/2020.09.16.20196253v1)

In context of i), list of falciparum genes where GATK pipeline gave lots of null/low quality calls:

LSA1 (PF3D7_1036400)
LSA3 (PF3D7_0220000)
DBLMSP1, DBLMSP2
CSP (PF3D7_0304600) - NANP repeat region
MSP2 (PF3D7_0206800) 
MSP1 (PF3D7_0930300)


## Paolo's project

Designing P. vivax major allelic classes at surface antigens, express them and identify antibodies binding to them in serum from exposed humans. Allows defining antibodies for serotyping, to monitor prevalence.

List of vivax genes where GATK pipeline also gave null/low quality calls:

PVP01_1031000- pvmsp3a
PVP01_1031100- pvmsp3 hypothetical
PVP01_1031500- pvmsp3b 
PVP01_1031700- pvmsp3y
PVP01_0000130- fam-a hypothetical
PVP01_0102300- pv dbp2
PVP01_1220400- msp7a

Ref genome they used: PvO1.


[myo_1]: https://app.box.com/s/105k1utbcqskclnghhobu1c2cgpgcdmk
