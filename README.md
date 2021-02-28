# Input data

We use Plasmodium data from the MalariaGEN projects: https://www.malariagen.net/ 
The following Sanger ftp ftp://ngs.sanger.ac.uk/production/malaria/ lists:
  - pf6 metadata/vcfs (under pfcommunityproject/Pf6)
  - pvivax vcfs (under pvgv)

# Workflows

## Running a workflow

To run a workflow, ... [TODO]

## List of workflows

### download_data

### call_variants

### eval_varcalls


## Development

To avoid rebuilding a container everytime a different version of the gramtools codebase is tested, it is installed locally:

```
. venv/bin/activate
pip install -e <path/to/gramtools>
```

This means the file pyrequirements.txt should not be produced by running `pip freeze > pyrequirements.txt`.


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
