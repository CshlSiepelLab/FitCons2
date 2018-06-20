# covtree

For an overview on how to easily and automatically download and use the Insight2 & covtree programs, please see the README.md file in the root directory of this repository.

## Overview

The covtree software is a high-performance decision tree inference tool (following a [CART](https://en.wikipedia.org/wiki/Decision_tree_learning) methodology). The covtree program utilizes tree-conditioned likelihood improvement to guide selection of single-covariate decision rules.

The covtree software is optimized to utilize [Insight2 software](https://github.com/CshlSiepelLab/FitCons2/tree/master/Insight2) to generate likelihoods. Insight2 is utilized by executing multiple instances of Insight2 in server mode as a sub-processes, generating input files representing exhaustive covariate bi-partitionings, providing the generated input files to Insight2, then recovering Insight likelihoods from Insight2 model files.


## Publication Checklist

### General	

 - Stand Alone Binary: For Linux is available in repository as [FitCons2/covtree/bin/covtree.bz2](https://github.com/CshlSiepelLab/FitCons2/blob/master/covtree/bin/covtree.bz2)
 - Database needed to run covtree: 
    - Full database available [here](https://github.com/CshlSiepelLab/FitCons2/blob/master/covtree/bin/covtree.bz2)
    - Demonstration subsets available via demo scripts in FitCons2/covtree/tests . See top level repository documentation  for instructions.
 - Demonstration test: in test, see **Demos** section, below, for instructions

### System Requirements 
 - Software Dependencies
   - Primary target is Linux (RedHat Enterprise 6.6 or above x64) for demos and tests. gcc 4.9.2+ std=c++14
   - Windows 7 & 10 (x64) have been used as well, and Visual C++ 2015 project files are included.
   - Code is C++ and uses C++2014 features.
   - No external dependencies or libraries other than Insight2 are required.
  - Tested on 
    - Linux (RedHat Enterprise 6.6 & 6.9, x64) as well as Windows 7 & 10 x64
 - Hardware Requirements:
   - 64-bit architecture and operating system. Requires only 1 core, but demos are optimized for 3 cores (recommended).
   - 20GB+ of RAM is recommended, along with 40GB of free hard drive space
   - Full database requires approximately 6GB to load, +4 GB for each Insight2 instance. Demos require less.
### Installation Guide:
 - Instructions:
   - covtree is best installed by downloading the entire FitCons2 repository (contains covtree) and installing using the ```./configure.sh``` command at the top level. This will configure covtree and Insight2,m as well as downloading the full compressed database. Running associated demos will download and decompress database subsets needed to operate software.
### Demos
Follow demonstration instructions listed in README.md in the [top level of the FitCons2 repository](https://github.com/CshlSiepelLab/FitCons2).

covtree arguments are as listed below, 
```
$ ../../bin/covtree -h
Tue Jun 19 20:59:51 2018 Main: Entering initialization.


covtree Version: 00.00.12 (Jun 15 2018 23:12:23)

covtree fCovDef fInsight2 dInsight2DB dOutBase [args]

Required Arguments:
        fnCovDef      - Path to covariate definition file, see documentation for details.
        insightExe    - Path to Insight2 executable.
        dInsight2DB   - Directory containing Insight2 database.
        dOutBase      - Directory at top of hierarchy containing results.

Options:
        -h                - Print this help menu and exit.
        -v [level]        - Verbosity flag. 0 is quiet. 1 is default. Higher values are increasingly verbose up to ~5.
        -iargs    [s]     - s is a quoted string containing any special Insight2 arguments. See Docs for Insight2. Defult is -qmap.
        -idb      [d]     - Local Insight2Database Directory. Overrides Insight2DB argument, but is not persistable.
        -iprocs   [u]     - Number of instances of Insight to run in parallel. Each uses 3-6 GB of memory. Defult is 1.
        -refgen   [f]     - f is a psth to a sorted bedfile defining refernce genome limits. Hg19 is default. Covaraites must fully parition this reference space.
        -bedmask  [f]     - f is pathname to a sorted bedfile defining subset of genomic positions. Only operate on these positions (not yet implimented).
        -covstart [s]     - Begin with covaraite subset. Serialized covsetDefSubset object. See Docs. Generally like ...Anno:{CvnA1;1,2};{CvnA2;3,4,5}:CTS:{CvnS1;1,2};{CvnS2;4,5,8}.
        -scts     [s]     - Comma seeperated list of cell-type names.
        -ipriors  [s]     - Comma seeperated list of floats rho,eta,gamma,N,DatNLL,PosteriorNLL any can be mising and will default to hg19. N is pseudocounts.
        -implprho [d]     - Alone or followed by 1 sets implicit-parent-Rho mode. ParRho = weighted sum of child rho. If 0 (default) use -iprior rho as parent rho, if it is provided, >0, and <1.
        -rhosplit         - if present, split on rho entropy, rather than default of Insight2NLL.
        -fout     [s]     - Nale of output file for external recursion. Written to dOutBase/[s]. [s].done ic completion semaphor file.
        -erecurse [u]     - write [u]=0, 1 or 2 children to fout for external recursions. Others are calculated serially in this process.
        -minbits  [d]     - Terminate recursion when best split yields fewer than this many bits. Def is 100,000. Must be positive. 100 is small, 1e6 is big.
        -hcpos    [u]     - Minimum number of positions in each child class needed to consider a split. Defaults to 0. Generally 300-10,000 is good. 1000 is naive est.
        -hcerr1   [d]     - Maximum Bonferroni Corrected probability of ordering error in childs. Defaults to 1.0 (disabled) 0.05 is good est.
        -hcerr2   [d]     - Maximum Bonferroni Corrected probability of distance error in childs (rho1-rho2)<(rho1_hat-rho2_hat)/2. Defaults to 1.0 (disabled) 0.05 is good est.
        -hcmndinf [d]     - If >0 use expected directed Information to determine ordering: E[DeltaInf(rho1,rho2)&(rho1>rho2)] Instead of point estimate DeltaInf(rho1_hat,rho2_hat). Value is min inf.
        -maxdepth [u]     - Terminate recursion no later than this depth. If 0, program does nothing and exits.
        -restart  [s]     - Restart recursion from previous covtree.out file in OutBase. Creats Hi/Lo sub dires. String Vals = Hi, LO or PR (both) children.
        -dwrk     [s]     - Write all physical files in this directory, however recursion logically proceeds on dirbase. Shell muse move files, so only useful when erecurse = 2.
        -debug    [d,[s]] - Debugging options d=1-skip singelton,2-skip partiton,3-skip both,0-skip none; s=process only single chromosme with tag s.
        -sarg     s       - If s=- print serialized args and exit, else desterilize next arg as arglist.

CovTree: Unable to initialize runtime elements.
        Init() returns: runtimeInit: Unable to parse user arguments

```

### Full covtree database
The full epi/genomic property database used to generate the reported FitCons2 results may be found in a compressed form in ```FitCons2D/covtree/db```
This database requires around __2.5GB__ in its compressed form. It can be expanded using the command 
```
$ cd FitCons2/covtree/db
$ chmod a+x unroll.simple.sh; ./unroll.simple.sh
```
This can require two hours to decompress the ```.starch``` files into ```.bedg``` files. The results occupies an additional 90GB of storage space.



### 	Instruction - how to run software
- Command Line: See **FitCons2/covtree/tests/runtest1.sh**,  for invocation.
- Arguments: See, above for list of parameters.
 ### Additional:
 - License: See LICENSE file in repository.
 - Link to repository: 
    - Repository Base	https://github.com/CshlSiepelLab/FitCons2
    - Repository Detail https://github.com/CshlSiepelLab/FitCons2/tree/master/covtree
    - Repository download link and command: `git clone https://github.com/CshlSiepelLab/FitCons2.git`
 - Pseudocode: see [FitCons2/covtree/docs/Covtree_FunctionalOverview.docx](https://github.com/CshlSiepelLab/FitCons2/blob/master/covtree/docs/Covtree_FunctionalOverview.pdf)
	
