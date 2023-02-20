----------------------------------
# TODOs

## HMM logos

* Produce hmm from each split msa
  Command-line: `hmmbuild --fragthresh 0 --enone --symfrac 0 --pnone`

## Mosaic

[x] Automate production of pngs from svgs
  Command-line: `inkscape --without-gui --export-png="${basename}.png" --export-dpi 96 ${fname}`

## CNVs

[x] Checked for copy-number changes in all my studied isolates
    None of the high-sharing samples (putative gene convs, in paralog_diversity_and_divergence.py) have
    fold-cov < 0.5, one has fold-cov > 2 [compared to AMA1]: PV0287-C. 
    To be safe, I filtered out all three with fold-cov > 1.5: QG0337-C, PF0655-C, PV0287-C


## Repeatiness

Possible tools for quantification:

* http://repeatmasker.org/RepeatMasker/
* 'trf' (VNTRs; see HPRC Liao et al paper)

## QC

[x] For samples with evidence of conversion, measure i)Read coverage ii)Read paired-end distances along DBs and AMA1, and
look for any outlier samples by comparison of DB values (test) and AMA1 values
(control).

---------------------------

# Ideas

## Plots

* Matrix of DB1 (x) vs DB2 (y) at DBL domain, and each square coloured by LD/KL div
* Same as above but coloured by enrichment for same AA

## QC

* Still extract the VCF level FORMAT fields and can plot for e.g. the DP distribution
  (e.g. what proportion of variants have DP > 0) in DBL domain

*  XA tag in SAM alignment gives alternate mappings for each read


## Where do the diverged bits come from

* Larremore et al (2015) dataset, supp. fig. 4. Get a hold of that dataset?

## Coming up with a model

* Is the conversion biased
* Are the private forms possibly derived from other DBL domains, e.g. var genes

---------------------------


# Useful info

## Duplications and premature stop codons

* Tetteh et al. (2009) (and D. Conway, who is last author on that study, from personal
  communication) find 3/14 lab isolates. isolates with two DBLMSP haplotypes, one of which is
  identical across all three lab isolates, and different from all others. Called them
  'copy B', highlighted in blue in supp. fig. S5 (10.1371/journal.pone.0005568)

  Lab isolates: FCC2, D6 and Palo Alto.
  --> Could be DBLMSP2 seq? (PCR chimeric product)
      (or DBLMSP2 --> DBLMSP gene conversion in a subset of samples!)

* Van Tyne et al. (2011) find a SNP associated with halofantrine resistance in DBLMSP2
  (https://doi.org/10.1371/journal.pgen.1001383). And find that overexpression/>1
  copy-number of gene increases res to three antimalarials. And find that several
  culture-adapted isolates have >1 DBLMSP2.

## Direct evolution

* Of all cross and clone tree lab isolates that are also in 15 PacBio sequenced by Otto
  et al. (2018), PfKH01 looks the most gene-conversion-like, but not actually a lot. 
  The one that really looks like it has gene conv, but apparetnly hasn't been clone treed/crossed,
  is PfGN01.

  --> This established by running the `paralog_diversity_and_divergence.py` script 
  (the part making the codon identity matrix) on the 15 isolates.

## Gene conversion

* Origin of lab clones: 3D7 cloned from an airport malaria case in the
  netherlands; Dd2 from Southeast Asia; HB3 from Honduras (Tetteh et al 2009)
* Using longitudinal samples, or taking a malariaGEN subset from a given time and place,
  and showing the samples are close from rest of genome, but conversion seems to happen
  in DBL
* The following samples have at least one premature stop codon in at least one of the DBs, 
  and have evidence of gene conversion:
    PT0079-C
    PA0383-C
    QG0299-C
    QQ0119-C
    PF0946-C
    PF0769-C
    PF0185-C
    PF0427-C
    FP0093-C
    QC0275-C
  Have not been analysed so far
* The following samples have evidence of a conversion + flanking deletion:
  - PF0304-C (Ghana)
  - PC0302-C (Kenya)


## Discussion (from thesis)

New pipeline for antigens, performs v. well on DBs
    - Recap
    - Some measure of how well my pipeline does across all genes; appendix figure
    - Which genes are not resolved, and possible ideas why

Widespread recomb and gene conv. in DBs
    - Recap
    - Literature on prevalence of recomb and gene conv in Pf. Similarities/differences
    - Methods to study recomb/conversion. Mosaic; bockhorst. Lack of regular phylo-structure
    - Recomb rate estimation? Which clade did DBs appear in?
Functional bio of DBL co-evolution
    - Where did three forms come from? What is relationship to other DBL domains? (Larremore et
      al)
    - Biology: structural analyses: mapping polymorphism on structures, docking to human receptors
    - Biology: different life stages 
    - Co-evolution with human variation: Band et al type analysis
