# FitCons2

## Synopsis

**FitCons2** employs hominin selective pressure to identify patterns of epi/genomic features that characterize human genomic loci with the potential for biological function.

For details please see our [bioRxiv posting](https://www.biorxiv.org/content/early/2018/05/23/317719), in brief,
 - Under an assumption of random genomic mutation, variation that inhibits procreation is found at a reduced frequency. This can be thought of as selective constraint, and can be quantified as the reduction in observed variation from what would be expected at random.
 - Mutation is relatively rare, so identification of recent selective pressure requires identifying a depletion in the frequency of already-rare events, and thus it tends to be underpowered. We improve power by measuring depletion in variation over a collection of genomic positions that share epi/genomic features.
 - To identify selective pressure we use the [INSIGHT](http://compgen.cshl.edu/INSIGHT) model to identify selective constraint. A highly optimized implementation of this model is provided below as [Insight2](https://github.com/CshlSiepelLab/FitCons2/tree/master/Insight2).
 - To identify the most informative patterns of epi/genomic properties a decision tree model is implemented across a set of quantized epi/genomic features, each spanning the entire human genome, and available for all 127 [Roadmap](http://www.roadmapepigenomics.org/) cell-types. This tree deries from a [CART](https://en.wikipedia.org/wiki/Decision_tree_learning) methodology, with each decision tree node representing a binary split on a single covariate.
   - Each split is chosen via exhaustive search, so as to maximize the joint likelihood of all genomic variation across all cell-types under the INSIGHT model.  An optimized decision tree model, that invokes Insight2 as needed is provided as [CovTree](https://github.com/CshlSiepelLab/FitCons2/tree/master/covtree), below.

## Using the Code and Data
### Overview

***CAUTION***. This is not commercial-quality software, and was designed for a relatively limited investigational purpose.

The code in this repository includes source code, script to build/configure the code, and scripts to download and analyze precompiled databases. The databases include:
 - Insight2 reference database derived from primate primary sequences and human variation data.
 - A complete version of the epi/genomioc database used in the references bioRxiv pre-publication. Complete analysis of this data set is computationally demanding and is not described here.
 - Subsets of the epi/genomioc database suitable for realtime analysis as code samples.

### Build instructions
To examine the code in this repository, first clone the Git repository via ```git```.

```$ git clone https://github.com/CshlSiepelLab/FitCons2.git```

or manually download an extract the current repository via

```https://github.com/CshlSiepelLab/FitCons2/archive/master.zip```

The repository requires approximately __20MB__ of space, and is configured for Linux (thought Visual Studio 2015 project files and Windows7/10 x64 executables are included as well). There are no external code dependencies. At the top level directory configure (compile or extract prebuilt binaries if compile fails)

```$ cd FitCons2; chmod a+x ./configure.sh; ./configure.sh```

Alternatively, you can visit each of the FitCons2 sub directories
```bin```, ```Insight2```, and ```covtree``` and type 

```$ chmod a+x ./configure.sh; ./configure.sh```

This configuration takes approximately 4 min to complete. The process downloads compressed databases for each program, so the time to completion can depend on the speed of your internet connection. When configuration is complete the file space requirement is around __9 GB__. 

### Demo Instructions

#### Demo: Insight2
The Insight2 demo prints the help page
To run the Insight 2 demo, 
```
$ cd FitCons2/Insight2/tests; chmod a+x ./runtests.sh;
$ ( ./runtests.sh 2>&1 >./runtests.sh.log)
```
The demos require up to 8 cores and a total of 10GB or RAM. The tests complete in approximately 5 Mins and expands the storage requirements to __10 GB__.

This command places demo results to the ```runtests.sh.log```. The demos are as follows
0. Help Menu, this is the same as executing the command ```Insight2 -h```
1. Generate maximum likelihood values for all of hg19
2. Generate maximum a-posterior values for the FitCons2 class 17, using genome-wide priors. This demonstrates the use of positional weights. The corresponding input ```.bed``` file has a number 1-115 at each locus, representing the number of the 115 training cell types that locus in class 17.
3. Generate full distribution and posterior expectation over the three model parameters rho, eta and gamma via a gridded sampling of likelihood and prior for FitCons2 class 17.

#### Demo: covtree
To run the covtree demo,  
```
$ cd FitCons2D/covtree/tests
$ chmod a+x ./runtest1.sh; (./runtest1.sh 2>&1 ) > runtest1.sh.log
$ chmod a+x ./runtest2.sh; (./runtest2.sh 2>&1 ) > runtest2.sh.log
```

These demos unpack FitCons2 genomic feature databases. Each demo requires up to 4 cores and a total of 16GB of RAM, however by modifying the ```runtest``` script to reduce the number of concurrent instances of Insight2 from 3, to 1, the demos can be reconfigured to run more slowly using 1 core and 8GB of ram.

##### covtree Demo 1: 3 cell types, 3 genomic features, all monotonic
The first demo generates a FitCons covariate tree based on 3 cell types (E003-H1hESC, E116-GM12878, and E122-HUVEC) and 3 covariates (CDS, DNase-seq, and RNA-seq), to a maximum tree depth of 5 and a minimum information gain of 100 bits per split.
- The root of the result tree is found in FitCons2/covtree/tests/1-basic/res.
  - The child nodes of this root are in the ./hi and ./lo subdirectories.
  - The exhaustive covariate splits are found in ./part, 
  - The nonmonotonic covariate re-ordering data is in ./sing (empty in this demo, all covariates are monotonic)
-  Comparative results from the reference environment be found based in FitCons2/covtree/tests/1-basic/ref
-   The first demo requires approximately 32 min to run, and after the run the storage requirement is __34GB__

##### covtree Demo 2: 3 cell types, 4 genomic features, one non-monotonic
 
The second demo generates a FitCons covariate tree based on 3 cell types (E003-H1hESC, E116-GM12878, and E122-HUVEC) but uses __4__ covaraites (CDS, MeltMap, DNase-seq, and RNA-seq), to a maximum tree depth of 5 and a minimum information gain of 100 bits per split. This demonstrates the use of a nonmonotonic covariate (MeltMap).
 - The root of the result tree may be found in FitCons2/covtree/tests/2-nonmono/res.
   - The child nodes of this root are in the ./hi and ./lo subdirectories. 
   - The exhaustive covariate splits are found in ./part,
   - The nonmonotonic covariate re-ordering data is in ./sing (entries for the nonmonotonic MeltMap covariate). 
 - Comparitive results from the reference environment be found based in FitCons2/covtree/tests/2-nonmono/ref
 - The second demo requires approximately 50 min to run, and after the run the storage requirement is __65__ GB.

### Full covtree database
The full epi/genomic property database used to generate the reported FitCons2 results may be found in a compressed form in ```FitCons2D/covtree/db```
This database requires around __2.5GB__ in its compressed form. It can be expanded using the command 
```
$ cd FitCons2/covtree/db
$ chmod a+x unroll.simple.sh; ./unroll.simple.sh
```
This can require two hours to decompress the ```.starch``` files into ```.bedg``` files. The results occupies an additional 90GB of storage space.

