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

Requirements: None

#### Function
3 main datasets:

  * pf6: downloads all p. falciparum read sets from malariaGEN
  * pvgv: downloads all p. vivax read sets from malariaGEN
  * pacb_ilmn_pf: downloads paired illumina reads and pacbio assemblies for 15 samples,
    from [Otto et al. (2018)][otto_2018]


### call_variants

Requirements: download_data

#### Function

* (Variant calling) Runs cortex and octopus on all pf6, pvgv and pacb_ilmn_pf samples
* (Genotyping) Adjudicates between cortex and octopus using gramtools
* (Variant calling) Runs gapfiller on top of gramtools adjudicated output

### eval_varcalls

Requirements: call_variants

#### Function


### make_prgs

Requirements: call_variants

#### Function

* Makes a PRG (Population Reference Graph; aka genome graph) based on end output of call_variants workflow. Configurable parameters are:
       - Which pf6 samples to use for graph construction
       - Which genes to build the prg on (for e.g., set of 26 hypervariable genes I have
         selected for study; see [Genes](#Genes) section below)
       - What `min_match_len` to use for `make_prg`

### joint_genotyping

Requirements: make_prg

#### Function

Takes as input a genome graph made by make_prg, and runs gramtools genotyping on all specified samples.

### plasmo_paralogs

Requirements: joint_genotyping

#### Function

Produces sequences of all the genes in the gramtools-built and genotyped PRG. So-named
because also concatenates paralog sequences together (e.g., DBLMSP and DBLMSP2)

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

Some samples are duplicates with low coverage, estimated to be mixed species, 'continent mismatches'
(labeled from a continent but their genotypic data maps them to another); these are marked as not 'Analysis_set' in the
'Exclusion reason' column. 

Further, some samples are mixed within-species infections. An estimate of clonality is
given by Fws- malariaGEN typically use >0.95 as threshold for clonality.

We build graphs from and genotype only 'Analysis_set' and fws>0.95 samples.
To get these, inside `analysis/input_data/sample_lists/pf6`:
```sh
grep -f <(awk '{if ($2 > 0.95){print "^"$1}}' Pf_6_fws.tsv) pf6_samples.tsv | grep "Analysis_set" | cat <(head -n1 pf6_samples.tsv) - > analysis_set_fws95.tsv
```
Note the code in `analysis/workflows/common_utils` also derives this (and other) sample subsets.

## Generational datasets

This includes clone trees and crosses datasets, literature refs and data sources are given below.

### Clone trees

* Claessens et al. (2014)[doi](): clone trees for parental lines 3D7, W2, Dd2, HB3. All
  parental lines are lab strains. 197 WGS samples in total.
* Hamilton et al. (2016)[doi](https://doi.org/10.1093/nar/gkw1259): clone trees for
  parental lines KH1-01 and KH2-01. These were culture-adapted from clinical samples
  from Cambodia. 87 WGS samples in total.

To build the input tsvs, I first took the initial excel files containing sample IDs, tree IDs
and ENA accessions from the supplementary materials of each paper, and combined them together. 
I then added a parent tree ID for each sample, using the clone tree diagrams in each
paper, also in the supplementary material. This is stored in analysis/input_data/sample_lists/generational_samples/clone_trees.tsv
For six samples, I found convincing evidence of sample mislabeling: across several
genes, I find large distances between the inferred sequences for the parent and the
child samples, but a distance of zero to another parent in a different clone tree. I
corrected these, and stored the new sample tsv in analysis/input_data/sample_lists/generational_samples/clone_trees_corrected.tsv

### Crosses

The crosses input tsv is produced from the input txt files of each of two papers as follows:
* Miles et al. (2016) table: ftp://ngs.sanger.ac.uk/production/malaria/pf-crosses/1.0/samples.txt, which I call `miles_etal_crosses.txt`
* Garimella et al. (2020) table: https://github.com/mcveanlab/Corticall/raw/master/manuscript/manifest.illumina.txt which I call `garimella_etal_crosses.txt`

```sh
echo -e "Sample\tCross\tClone\tAccession" > crosses.tsv
awk -F'\t' '{print $3"\t"$1"\t"$2"\t"$4}' miles_etal_crosses.txt | tail -n+2 >> crosses.tsv
awk -F'\t' '{print $2"\t"$1"\t"$3"\t"$4}' garimella_etal_crosses.txt | grep "803xGB4" | sed s/803xGB4/803_gb4/ >> crosses.tsv
```

Column `Clone` allows identifying the parent samples for each cross. All accessions are 
run IDs, but all samples were sequenced in single run (as per Miles et al. (2016) supplementary)

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
