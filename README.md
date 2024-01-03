# Genotyping highly-variable *P.falciparum* genes

author: Brice Letcher

## Summary

This project genotypes malariaGEN samples at a set of specific *P. falciparum* genes 
using a newly-developed genome-graph-based pipeline. See the [documentation](docs/pf_genes.pdf) for a list and rationale of
the genes picked for analysis.

Two of these, DBLMSP and DBLMSP2, were analysed in detail, supporting the following publication:

> doi: https://doi.org/10.1101/2023.02.27.530215

Some of the others are analysed a little bit more in my thesis:

> Genome-graph based genotyping with applications to highly variable genes in *P. falciparum*. https://doi.org/10.17863/CAM.93368

## Publication reproducibility

* To regenerate results, run the Snakemake workflows as detailed below. These rely on version-frozen software 
  through a singularity container, available in the [zenodo][zenodo] release associated with the publication above.
* To regenerate figures in the publication above, also refer to this [zenodo][zenodo] release.

## Project structure

The project is divided into Snakemake workflows located at `analysis/workflows`.
Code used by the workflows is in `analysis/scripts`.

Reproducibility is mediated by a singularity container whose definition is at
`reproducibility/container/singu_def.def`.

A built version of this container is expected by the workflows in
`reproducibility/container/built/singu.sif`. A built version is available on [zenodo][zenodo].

An additional repository is included in this project, as a git subtree, called
`plasmo_paralogs`. That sub-project runs all the sequence analyses of DBLMSP and
DBLMSP2 and is structured in the same way as this one (README.md, Snakemake workflows, scripts).

## Workflows

### Constraints & Conventions

In each workflow I load a global configfile ("common.yaml") and a set of common python
utilities ("common_utils.py"). "common.yaml" must be loaded before "common_utils.py",
and "common_utils.py" must be loaded before other "utils.py".

Inside workflows, functions prefixed with `cu_` come from `analysis/workflows/common_utils.py`

Similarly, each workflow's own `analysis/workflows/<workflow_name>/utils.py` defines
functions prefixed with a two-letter code: e.g. for `download_data`, functions will
start with `dd_`.


### Running a workflow

To run a workflow on the cluster, call
```
sh analysis/cluster_submit.sh
```

### List and description of workflows

#### `download_data`

Requirements: None

Downloads 3 main datasets, using input TSVs located at `analysis/input_data/sample_lists`
(and also released on [zenodo][zenodo]):

  * `pf6`: downloads all p. falciparum read sets from malariaGEN [pf6 release][pf6_release]
  * `pacb_ilmn_pf`: downloads paired illumina reads and pacbio assemblies for 15 samples from [Otto et al. (2018)][otto_2018]
  * `generational_samples`: downloads clone tree data from [Hamilton et al. (2017)][hamilton_2017] and crosses data (see paper for references)
  * `laverania_illumina`: downloads Illumina reads from non-pf laverania species from [Otto et al. (2018)][otto_2018b] (note this is a different paper reference than for `pacb_ilmn_pf` ;-))


#### `call_variants`

Requirements: download_data

Performs genotyping using a range of tools.

* (Variant calling) Runs cortex and octopus on all pf6 and pacb_ilmn_pf samples
* (Genotyping) Adjudicates between cortex and octopus using gramtools
* (Variant calling) Runs gapfiller on top of gramtools adjudicated output

#### `eval_varcalls`

Requirements: call_variants

Evaluates variant calls, using two orthogonal strategies (see paper).

#### `make_prgs`

Requirements: call_variants

* Makes a PRG (Population Reference Graph; aka genome graph) based on end output of call_variants workflow. Configurable parameters are:
       - Which pf6 samples to use for graph construction
       - Which genes to build the prg on (for e.g., set of 26 hypervariable genes I have
         selected for study; see [Genes](#Genes) section below)
       - What `min_match_len` to use for `make_prg`

#### `joint_genotyping`

Requirements: make_prg

Takes as input a genome graph made by make_prg, and runs gramtools genotyping on all specified samples.

#### `plasmo_paralogs`

Requirements: joint_genotyping

Produces sequences for use by `plasmo_paralogs` repository.

#### `non_pf6_samples`

Genotypes Illumina-sequenced samples not part of malariaGEN's [pf6][pf6_release]:
    - Generational samples: pf clones and crosses
    - Laverania samples: non-pf laverania samples

For the laverania samples, four datasets are available:

    - PPRFG01: the one with Pacb data as well, from which assembly was built
    - PPRFG02: DBLMSP, DBLMSP2, AMA1 were all fully resolved at the end of
      `non_pf6_samples` worfklow.
    - PPRFG03: almost no read coverage in DBLMSP, DBLMSP2, and AMA1 after joint genotyping on the pf6-graph
      This may be due to the high rate of host contamination (~86%) ([Otto et al. (2018)[otto_2018b] supp. table 1).
    - PPRFG04: has more diversity than the three above - this could be due to it being highly mixed with P. adleri. DBLMSP was almost fully resolved at the end of the `non_pf6_samples` workflow, not DBLMSP2.

Sanity check: I ran `bcftools view -R analysis/input_data/otto_2018_core_genome_def.bed -f "PASS" -v snps | wc -l` on each vcf made by octopus for:
    - The 15 `pacb_ilmn_pf` Pf samples
    - The four laverania samples above
And get ~18k SNPs for the Pf samples vs > 130k SNPs for the laverania samples; including for PPRFG02, which I initially found suspiciously not-diverged from 3D7, when looking at read alignments.

## Development

### Testing without built container
To avoid rebuilding a container everytime a different version of the gramtools codebase is tested, it is installed locally:

```
. venv/bin/activate
pip install -e <path/to/gramtools>
```

This means the file pyrequirements.txt should not be produced by running `pip freeze > pyrequirements.txt`.


## Input data

We use Plasmodium data from the MalariaGEN projects: https://www.malariagen.net/ 
The following Sanger ftp ftp://ngs.sanger.ac.uk/production/malaria/ lists:
  - pf6 metadata/vcfs (under pfcommunityproject/Pf6)
  - pvivax vcfs (under pvgv)

### Reference genomes

For *P. falciparum*, I have used the 2018-11 version for variant calling and graph
building (ftp://ftp.sanger.ac.uk/pub/project/pathogens/gff3/2018-11/). pf6 used the
2016-07 version in their release VCFs (see pf6 paper). However, I checked the following:

* The two genomes have exactly the same chromosome sizes
* The two annotation files (GFF) have exactly the same coordinates for pf6_26_genes

On that basis it is fine to use 2018-11 version for evaluation purposes (see eval_varcalls workflow)

### pf6

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

### Generational datasets

This includes clone trees and crosses datasets, literature refs and data sources are given below.

#### Clone trees

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

#### Crosses

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

## Genes

### P. falciparum

See the [documentation](docs/pf_genes.pdf) for a list and rationale for the genes we chose to analyse.

[pf6_release]: https://doi.org/10.12688/wellcomeopenres.16168.1
[hamilton_2017]: https://doi.org/10.1093/nar/gkw1259
[otto_2018]: https://doi.org/10.12688/wellcomeopenres.14571.1
[otto_2018b]: https://doi.org/10.1038/s41564-018-0162-2
[zenodo]: https://zenodo.org/doi/10.5281/zenodo.7677547
