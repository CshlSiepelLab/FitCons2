## Synopsis

**Insight2** is a high performance implimentation of the [INSIGHT](http://compgen.cshl.edu/INSIGHT) statistical model for estimating number of genomic positions under selective pressure in a collection of intersperced loci. The INSIGHT model can be thought of as a generalization of the [McDonald-Kreitman](https://en.wikipedia.org/wiki/McDonald%E2%80%93Kreitman_test) framework to include non-coding DNA. INSIGHT is sensitive to strong purifying, weak purifying and adaptive selection, using polymorphism data from a collection of individuals in a species and divergence from neighboring outgroup species. Insight2 is approximately 10,000x faster than INSIGHT, but is limited to hg19 genomic positions in humans, with a set of fixed evolutionary parameters derived from the human-chimpanzee common ancestor, a set of polymorphism data from 54 relatively unrelated individuals (subset of [Complete Genomics data](http://www.completegenomics.com/public-data/69-genomes/)), and a fixed set of putative neutral positions consisting of approximately 1/2 of the human genome.

## Motivation
Selective pressure can be used as an indicator of the [potential for genomic function](https://www.biorxiv.org/content/early/2014/09/11/006825). When this indicator is used as part of a machine learning algorithm to analyze collections of genomic loci (as per [FitCons2](https://www.biorxiv.org/content/early/2018/05/09/317719)) rapid evaluation and large data set size is critical for realistic algorithm development. The original INSIGHT implementation was limited to approximately 10 mega-base regions each requiring tens of minutes to calculate. **Insight2** can process a monolithic data set of all 2.8 billion human genomic positions in approximately 90 seconds. In addition, Insight2 accepts fractional inclusion of individual positions, can apply Bayesian priors, and produce posterior probabilities for parameter estimates.

## Publication Checklist
### General	

 - Stand Alone Binary: For linux is available in repository as FitCons2/Insight2/bin/Insight2.bz2
 - Database needed to run Insight2: Download [here](http://nextgen.cshl.edu/~bgulko/research/Insight2/db/FitCons2/Insight2DB.tar), see installation instructions for use.
 - Demonstration test: in test, see **Demos** section, below, for instructions

### System Requirements 
 - Software Dependencies
   - Primary terget is Linux (RedHat Enterprise 6.6 or above x64) for demos and tests. gcc 4.9.2+ std=c++11
   - Windows 7 & 10 (x64) have been used as well, and Visual C++ 2015 project files are included.
   - Code is C++ and uses C++2011 features, but no C++2014 features.
   - No external dependencies or libraries are required.
  - Tested on 
    - Linux (RedHat Enterprise 6.6 & 6.9, x64) as well as Windows 7 & 10 x64
 - Hardware Requirements:
   - 64 bit architecture and operating system.
   - 8GB+ of RAM is recommended, along with 10GB of free hard drive space
   - Database requires approximately 3GB to load. A 32 bit build is possible, but not supported.
### Installation Guide:
 - Instructions:
   - Obtain and compile the [bedops](http://bedops.readthedocs.io/en/latest/index.html) package to obtain the [unstarch](http://bedops.readthedocs.io/en/latest/content/reference/file-management/compression/starch.html) decompression program 
   - Download [Insight2 database](http://nextgen.cshl.edu/~bgulko/research/Insight2/db/FitCons2/Insight2DB.tar) (about 0.5GB) and expand files using unstarch and gunzip (expands to about 3.2GB). Ihe Insight2 database directory (`didb`) will be **`FitCons2/Insight2/Insight2DB`** <br>
`wget http://nextgen.cshl.edu/~bgulko/research/Insight2/db/FitCons2/Insight2DB.tar`<br>
`tar -xvf Insight2DB.tar`<br>
`gunzip Insight2DB/monoDB.db.gz`<br>
`unstarch Insight2DB/block.bedg.starch > Insight2DB/block.bedg`<br>
`unstarch Insight2DB/poly.bedg.starch > Insight2DB/poly.bedg`<br>
`unstarch Insight2DB/polyn.bedg.starch > Insight2DB/polyn.bedg`<br> 
   - Clone the [FitCons2 GitHub repository](https://github.com/CshlSiepelLab/FitCons2), this contains the Insight2 code, about 10Mb.<br>
 `git clone https://github.com/CshlSiepelLab/FitCons2.git`
     - Decompress the pre-built executable.. <br>
  `  cd FitCons2/Insight2`<br>
  `  bunzip2 bin/Insight2.bz2`<br>
  `  chmod u+x bin/Insight2`<br>
     - or, build the software using<br>
 `make`
   - No separate installation is needed, to view options type<br>
 `bin/Insight2`
 - Total install time: after download, <5 minutes.
### Demos
#### Demo 0 - Help Menu
 - Execute Line:<br>`bin/Insight2 -h`
 - Expected Output
 ```
 INSIGHT2: Error Processing Input args:



Insight2 Version: 0.16e

Insight2 DirDB [-fin Fname] [args, become defaults in server mode]

May be invoked in 2 modes - single use and server.
        Single Use - requires specification of -fin argument. Reads db, processes command line and exits.
        Server     - without -fin, reads db, then awaits series of lines from stdin.
                        the word "done" on a line by itself terminates server mode and exits.
                        Line format:
                                InFileName # [list of command line args]

DirDB dierctory must contain 6 files
        block.bedg- string_chrom int_start int_end float_theta float_lambda
        poly.bedg - string_chrom int_start int_end char_freqLH float_priMaj float_priMin
        polyN.bedg- Poly file, but filtered to only contain positions in neutral loci.
        Monomorphism database, consistign of .tags .chroms .db
                monoDB.tags      int_ind str_tag  - a list of up to 255 tags, ind is in range 1-255, tag is any string (usually a string formatted float)
                monoDB.chroms    sorted lsit of chromosome extents, bed format. string_chromid 0 int_chromLen
                monoDB.db        binary file, one byte per position in .chroms, each byte is the index of an entry in .tags. 0 for missing data.

QuickCommands:
        -h            - Print this help menu.
        -v [level]    - Verbosity flag. 0 is default without flag. 1 is default with flag increasingly verbose up to ~5.
        -ddb dname    - override DirDB, but in the argument list. Has no effect when provided in server mode as DB is read before server loop is entered.
        -fin fname    - source file for positions of interest. Bed format. Do not include .bed suffix.
        -qmap [N]     - ML/MaP estimate for beta,rho,eta,gamma parameters. N must be an integer>=0. N=0 performs ML. N>0 (default:347) provides N pseudocoutns of hg19 distrib as prior.
        -qexp [N]     - Expectation for rho,eta,gamma parameters. N, if specified Must be >0, (default:347) hg19 pseudocount prior. Betas estimated via ML.
        -mkins fname  - Generate Insight1 input files from database and estimate betas.

DetailedCommands:
        -a int        - number of alleles, echoed in first line of output file. Default is 108.
        -mineta D     - Minimum value for eta in ML asessment. Generally 0. May be set lower to investigate ML stability. Less than -5 not reccomended.
        -datasum      - Primt summary of nucleotide positions and intersection with database.
        -postsum      - Summary of Insight1 Posterior information, over all sites innput file from refined model.
        -postdet      - Summary of Insight1 Posterior for each site with database information.
        -postall      - Insight2 Posterior allele distributions for each site with database information.
        -postcntfac I - Insight2 Counterfactual analysis of impact of novel allele at each site w/ database information.
                        I is a 2 digit number XY, with Y=[1,2] representing type of analysis and Y=[0,1,2,3] verbosity of output.
        -expval       - Generate expected values for parameters rho, eta, gamma.
        -expdist      - Generate prior, likelihood and posterior probabilities at point sampled grid and marginals for each parameter.
        -gres         - Grid resolution for posterior parameter asessment. Runtime grows cubically. 5 is fast, 30 is impossibly slow. 10-20 typ. default 10.
        -nthread      - Number of threads used in posterior parameter asessment. Inefficent, so if >1, use >=4. 4-12 typ.
        -nll          - Print total data likelihood under model for positions in database, in nats.

        -maxsampw D   - Maximum input position weight (number of samples). Defaults to 1 (unweighted case).

        -priorwt D    - Strength of the prior measured as D pasudocounts of data from whole autosome hg19 distribution.

        -ins1comp [N] - Insight1 compatibilty mode def N=1: 0=Insight2, 1=use only blocks w/ informative poly sites in fin (Insight1).

Parameter initialization:
        For each parameter listed below, provide default values, as well as a flag indicating if that parameter should be refined.
        Generally Initial values are the hg19, whole autosome expectations. Priors default to Init values when not provided.
        Format is [Init[,Refine[,Prior]]
                Init   - real value for initial search value in ML gradient descent.
                Refine - 0 - supress refinement of this variable, 1 - allow.
                Prior  - Defailts to Init. Overrides init value as prior value for this parameter.

                Example: -rho .65,1,.077
        Available parameters are:
                -rho, -eta, -gamma, -lambda, -lambdaN, -theta, -thetaN
                lambda and theta are only used in prior model, and are generally be inferred from actual data. N referes to Neutral Sites.

        -betas b1,b2,b3 - override default values for beta when not inferred from data. Values renormalized to 1. Defaults are based on whole autosome hg19.




Exiting.
```
 #### Demo 1 - HG19
 - Execute line
 ` bin/Insight2 Insight2DB -fin ./tests/1-hg19/hg19.bed -qmap`
 - Expected output while running (date stamps may differ):
```
Thu Jun  7 19:16:41 2018        Reading binary Insight databse from Insight2DB/
Thu Jun  7 19:17:30 2018        Intersecting database with positive loci from ./tests/1-hg19/hg19.bed
Thu Jun  7 19:17:39 2018        Estimating Beta values
Thu Jun  7 19:17:40 2018        Estimating MAP/ML parameters
Thu Jun  7 19:18:38 2018        Recalculating data likelihood.
Thu Jun  7 19:18:49 2018        Sourcefile: ./tests/1-hg19/hg19  Runtime: 127.8000      secs
Thu Jun  7 19:18:49 2018        Done.
 ```
 
 - Expected output output files (compare with files provided in `./tests/1-hg19/ref`)
 ```
-rw-r--r-- 1 bgulko siepellab  379 Jun  7 19:18 tests/1-hg19/hg19.insres
-rw-r--r-- 1 bgulko siepellab 2747 Jun  7 19:18 tests/1-hg19/hg19.model
```
 - Expected run-time: 90-200 secs.
 #### Demo 2 - Weighted positions (FitCons2 class ID 17,  genome-wide priors)
 - Execute lines
`unstarch ./tests/2-weighted/FC2-017.bed.starch > ./tests/2-weighted/FC2-017.bed`
`bin/Insight2 Insight2DB -fin ./tests/2-weighted/FC2-017.bed -qmap -maxsampw 115`
 - Expected output while running (date stamps may differ):
 ```
Thu Jun  7 19:53:16 2018        Reading binary Insight databse from Insight2DB/
Thu Jun  7 19:53:32 2018        Intersecting database with positive loci from ./tests/2-weighted/FC2-017.bed
Thu Jun  7 19:53:33 2018        Estimating Beta values
Thu Jun  7 19:53:33 2018        Estimating MAP/ML parameters
Thu Jun  7 19:53:34 2018        Recalculating data likelihood.
Thu Jun  7 19:53:34 2018        Sourcefile: ./tests/2-weighted/FC2-017   Runtime: 17.9300       secs
Thu Jun  7 19:53:34 2018        Done.
```
 - Expected output output files (compare with files provided in `./tests/2-weighted/ref`)
 ```
-rw-r--r-- 1 bgulko siepellab     370 Jun  7 19:53 tests/2-weighted/FC2-017.insres
-rw-r--r-- 1 bgulko siepellab    2741 Jun  7 19:53 tests/2-weighted/FC2-017.model
```
 #### Demo 3 - Posterior and expectation - (FitCons2 class ID 17, parental priors)
  - Execute lines (note this uses 16 processes on the host computer)
`unstarch ./tests/3-posterior/FC2-017.bed.starch > ./tests/3-posterior/FC2-017.bed`
 `bin/Insight2 Insight2DB -fin ./tests/3-posterior/FC2-017.bed -qmap -maxsampw 115 -rho 0.359788,1 -eta 0.788613,1 -gamma 1.277413,1 -qexp -nthread 16 -gres 10 -expval -expdist`
 - Expected output while running (date stamps may differ):
 ```
Thu Jun  7 20:22:00 2018        Reading binary Insight databse from Insight2DB/
Thu Jun  7 20:22:24 2018        Intersecting database with positive loci from ./tests/3-posterior/FC2-017.bed
Thu Jun  7 20:22:25 2018        Estimating Beta values
Thu Jun  7 20:22:25 2018        Estimating MAP/ML parameters
Thu Jun  7 20:22:25 2018        Estimating posterior parameters distributions.
Thu Jun  7 20:23:49 2018        Recalculating data likelihood.
Thu Jun  7 20:23:49 2018        Sourcefile: ./tests/3-posterior/FC2-017  Runtime: 1207.4200     secs
Thu Jun  7 20:23:49 2018        Done.
```
- Expected output output files (compare with files provided in `./tests/3-posterior/ref`)
```
-rw-r--r-- 1 bgulko siepellab     272 Jun  7 20:23 tests/3-posterior/FC2-017.insres
-rw-r--r-- 1 bgulko siepellab    2735 Jun  7 20:23 tests/3-posterior/FC2-017.model
-rw-r--r-- 1 bgulko siepellab     551 Jun  7 20:22 tests/3-posterior/FC2-017.model.map
-rw-r--r-- 1 bgulko siepellab 6772665 Jun  7 20:23 tests/3-posterior/FC2-017.ppost.full
-rw-r--r-- 1 bgulko siepellab    6183 Jun  7 20:23 tests/3-posterior/FC2-017.ppost.marg
```
### 	Instruction - how to run software
- Command Line: See **Demo 0**, above for invocation.
- Arguments: See **Demo 0**, above for list of parameters.
 ### Additional:
 - License: See LICENSE file in repository.
 - Link to repository: 
    - Repository Base	https://github.com/CshlSiepelLab/FitCons2
    - Repository Detail https://github.com/CshlSiepelLab/FitCons2/tree/master/Insight2
    - Repository download link and command: `git clone https://github.com/CshlSiepelLab/FitCons2.git`
 - Pseudocode: *In Preparation*
	
