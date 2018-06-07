## Synopsis

**Insight2** is a high performance implimentation of the [INSIGHT](http://compgen.cshl.edu/INSIGHT) statistical model for estimating number of genomic positions under selective pressure in a collection of intersperced loci. The INSIGHT model can be thought of as a generalization of the [McDonald-Kreitman](https://en.wikipedia.org/wiki/McDonald%E2%80%93Kreitman_test) framework to include non-coding DNA. INSIGHT is sensitive to strong purifying, weak purifying and adaptive selection, using polymorphism data from a collection of individuals in a species and divergence from neighboring outgroup species. Insight2 is approximately 10,000x faster than INSIGHT, but is limited to hg19 genomic positions in humans, with a set of fixed evolutionary parameters derived from the human-chimpanzee common ancestor, a set of polymorphism data from 54 relatively unrelated individuals (subset of [Complete Genomics data](http://www.completegenomics.com/public-data/69-genomes/)), and a fixed set of putative neutral positions consisting of approximately 1/2 of the human genome.

## Motivation
Selective pressure can be used as an indicator of the [potential for genomic function](https://www.biorxiv.org/content/early/2014/09/11/006825). When this indicator is used as part of a machine learning algorithm to analyze collections of genomic loci (as per [FitCons2](https://www.biorxiv.org/content/early/2018/05/09/317719)) rapid evaluation and large data set size is critical for realistic algorithm development. The original INSIGHT implementation was limited to approximately 10 mega-base regions each requiring tens of minutes to calculate. **Insight2** can process a monolithic data set of all 2.8 billion human genomic positions in approximately 90 seconds. In addition, Insight2 accepts fractional inclusion of individual positions, can apply Bayesian priors, and produce posterior probabilities for parameter estimates.

## Publication Checklist
### General	

 - Stand Alone Binary: Is available bin/Insight2.bz2
 - Database Needed to run the system: ADD LINK
 - Demonstration test: in test, see XXXX for instructions

### System Requirements 
 - Software Dependencies
 -- Primary terget is Linux (RedHat Enterprise 6.6 or above x64) for demos and tests. gcc 4.9.2+ std=c++11
 -- Windows 7 & 10 (x64) have been used as well, and Visual C++ 2015 project files are included.
 -- Code is C++ and uses C++2011 features, but no C++2014 features.
 -- No external dependencies or libraries are required.
  - Tested on 
 -- Linux (RedHat Enterprise 6.6 & 6.9, x64) as well as Windows 7 & 10 x64
 - Hardware Requirements:
-- 64 bit architecture and operating system.
 -- 8GB+ of RAM is recommended, along with 10GB of free hard drive space
 -- Database requires approximately 3GB to load. A 32 bit build is possible, but not supported.
### Installation Guide:
 - Instructions:
 -- Obtain and compile the [bedops](http://bedops.readthedocs.io/en/latest/index.html) package to obtain the [unstarch](http://bedops.readthedocs.io/en/latest/content/reference/file-management/compression/starch.html) decompression program 
 -- Download [Insight2 database](http://nextgen.cshl.edu/~bgulko/research/Insight2/db/FitCons2/Insight2DB.tar) (about 0.5GB) and expand files using unstarch and gunzip (expands to about 3.2GB) 
`wget http://nextgen.cshl.edu/~bgulko/research/Insight2/db/FitCons2/Insight2DB.tar`
`tar -xvf Insight2DB.tar`
`gunzip Insight2DB/monoDB.db.gz`
`unstarch Insight2DB/block.bedg.starch > Insight2DB/block.bedg`
`unstarch Insight2DB/poly.bedg.starch > Insight2DB/poly.bedg`
`unstarch Insight2DB/polyn.bedg.starch > Insight2DB/polyn.bedg`
the Insight2 database directory (`didb`) is now **`FitCons2/Insight2/Insight2DB`**
 -- Clone the [FitCons2 GitHub repository](https://github.com/CshlSiepelLab/FitCons2), this contains the Insight2 code, about 10Mb.
 `git clone https://github.com/CshlSiepelLab/FitCons2.git`
-- Generate the executable.. 
  `  cd FitCons2/Insight2`
  `  bunzip2 bin/Insight2.bz2`
  `  chmod u+x bin/Insight2`,
  or
 `make`
 -- No separate installation is needed, to view options type
 `bin/Insight2`
 - Total install time: after download, <5 minutes.
### Demo
 - Instructions
 -- Execute line
 - Expected Output:
 -- This will produce following output while running:
 -- And the following output files:XXX
 -- The files should match those in: XXX 
 - Expected Run-time:
 -- First invocation should take XXX, this includes building high-performance db cache  files 
 -- Second invocation should take YYY
### 	Instruction - how to run software
- Command Line
- Arguments 
 -- -database, -fin, -maxsampw, -qmap, -priors, -posteriors, -priorwt, -rho, -eta,  - Server Mode
-gamma
 - HardCoded hg19 priors
 - ML / MAP
 ### Additional:
 - License: See XXX
 - Link to repository: YYYY
 - Pseudocode: ZZZZ
	

