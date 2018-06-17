#!/bin/bash


dbase="../.."
source "${dbase}/bin/lib.configure.sh"
# insight exe in $fiex
# insight db  in $didb

fcex="${dbase}/covtree/bin/covtree"
# Confirm Covtree
if [[ -x "$fcex" ]]; then
  echo -e "`date` Found covtree executable,\n\t\t${fcex} \n"
else
  echo -e "`date` Did not find covtree executable, configuring covtree."
  pushd ..
  ./configure.sh
  popd
  if [[ ! -x "$fcex" ]]; then
    echo -e "`date` Unable to configure covtree. Error. Exiting."
    exit 1
  fi
fi

# confirm Insight
if [[ -x "$fiex" ]]; then
  echo -e "`date` Found Insight2 executable,\n\t\t $fiex \n"
else
  echo -e "`date` Did not find Insight2 executable, configuring Insight2."
  pushd ../../Insight2
  ./configure.sh
  popd
  if [[ ! -x "$fiex" ]]; then
    echo -e "`date` Unable to configure Insight2. Error. Exiting."
    exit 1
  fi
fi



# covtree test 2
# Simple decomposition
dt2="./2-nonmono"; mkdir -p "$dt2"; dt2=`readlink -f "$dt2"`
echo -e "\n\n==============================================================================================================================================\n"
echo -e "`date` Configuring covtree Demo 2: Decomposition using CDS, Melt, RNA-seq and DNase-seq covaraites, for E003,E116,E122 (h1hesc, gm12878, HUVEC).\n"
echo -e "\n\t Working directory:\n\t\t${dt2} \n"


# Get database
ddb="${dt2}/db"; mkdir -p "${ddb}"
fdb="${ddb}/covtree-test2-DB.tar"

#furl='http://nextgen.cshl.edu/~bgulko/research/Insight2/db/FitCons2/Insight2DB.tar'
#surl='http://nextgen.cshl.edu/~bgulko/research/covtree/db/FitCons2/1-test1-simple/1-test1-simple-db.tar'
surl='http://compgen.cshl.edu/~bgulko/research/covtree/db/FitCons2/2-test2-nonmonotonic/covtree-test2-DB.tar'

if [[ -s "$fdb" ]]; then
  echo -e "`date` Found database archive $fdb . "
else
  echo -e "`date` Fetching database archive $fdb from\n\t ${surl} . "
  wget "$surl" -O "$fdb"
fi
if [[ ! -s "$fdb" ]]; then
  echo -e "`date` Unable to fetch database archive. Error. Exiting."
  exit 1
fi



# extract the database 
if [[ -s "${ddb}/X-Annotations/E999/covs.bedg.starch" ]]; then
  echo -e "`date` Found extracted database files."
else
  echo -e "`date` Extracting database files. "
  pushd "${ddb}"
  tar -xvf "$fdb"
  popd
fi
if [[ ! -s "${ddb}/X-Annotations/E999/covs.bedg.starch" ]]; then
  echo -e "`date` unable to extrct database files. Error. Exiting."
  exit 1
fi



# decompressing
if [[ -s "${ddb}/X-Annotations/E999/covs.bedg" ]]; then
  echo -e "`date` Found decompressed files."
else
  echo -e "`date` Decompressing database files. "
  pushd "${ddb}"
  ./unroll.simple.sh "Y"
  popd
fi
if [[ ! -s "${ddb}/X-Annotations/E999/covs.bedg" ]]; then
  echo -e "`date` unable to decompress  database files. Error. Exiting."
  exit 1
fi



# Get reference results
dref="${dt2}/ref"; mkdir -p "${dref}"
fref="${dref}/covtree-test2-ref.tar"
#surl='http://nextgen.cshl.edu/~bgulko/research/covtree/db/FitCons2/1-test1-simple/1-test1-simple-db.tar'
#surl='http://nextgen.cshl.edu/~bgulko/research/covtree/db/FitCons2/1-test1-simple/1-test1-simple-ref.tar'
surl='http://compgen.cshl.edu/~bgulko/research/covtree/db/FitCons2/2-test2-nonmonotonic/covtree-test2-ref.tar'
if [[ -s "$fref" ]]; then
  echo -e "`date` Found reference archive $fref . "
else
  echo -e "`date` Fetching database archive $fref from\n\t ${surl} . "
  wget "$surl" -O "$fref"
fi
if [[ ! -s "$fref" ]]; then
  echo -e "`date` Unable to fetch refrence archive. Warning. Continuing."
  #exit 1
fi

# extract the reference
if [[ -s "${dref}/X-Annotations/E999/covs.bedg.starch" ]]; then
  echo -e "`date` Found extracted reference files."
else
  echo -e "`date` Extracting reference files. "
  pushd "${dref}"
  tar -xvf "$fref"
  popd
fi
if [[ ! -s "${dref}/X-Annotations/E999/covs.bedg.starch" ]]; then
  echo -e "`date` unable to extrct refrence files. Warning. Continuing."
  #exit 1
fi


# make outdir
dres="${dt2}/res"; mkdir -p "$dres"

# run Covtree
echo -e "\n\n\n"
echo -e "`date` Executing covtree, with three Insight2 treads. This requires 3 cores and approximately 20GB of RAM"
echo -e "\t\t To reduce memory footprint reduce -iproces from 3 to 1. This should reduce RAM needs to around 10GB, with 1 core.\n\n\n"
#ctargs="-v 1 -iprocs 1 -scts E003,E116,E122 -erecurse 0 -minbits 100 -hcpos 10000 -maxdepth 5 "  # one thread, low memory
#ctargs="-v 4 -iprocs 3 -scts E003,E116,E122 -erecurse 0 -minbits 100 -hcpos 10000 -maxdepth 5 "  # show details of nonomontonic ordering
ctargs="-v 1 -iprocs 3 -scts E003,E116,E122 -erecurse 0 -minbits 100 -hcpos 10000 -maxdepth 5 "  # default run

fdb="${ddb}/covdef.test2-nonmono.txt"

# execute covtree!
pushd "$dt2"
echo "${fcex} $fdb $fiex $didb $dres ${ctargs}"
"${fcex}" "$fdb" "$fiex" "$didb" "$dres" ${ctargs}
popd

echo -e "\n\n`date` covtree execution complete - to verify output, compare output directory (test-1/res) with reference output (test-1/ref)."
echo -e "\n\n\n"
exit



if [[ 1 == 0 ]]; then
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
        -minbits  [d]     - Terminate recursion when best split yeilds fewer than this many bits. Def is 100,000. Must be positive. 100 is small, 1e6 is big.
        -hcpos    [u]     - Minimum number of positions in each child class needed to consider a split. Defaults to 0. Generally 300-10,000 is good. 1000 is niave est.
        -hcerr1   [d]     - Maximum Bonferroni Corrected probability of ordering error in childs. Defaults to 1.0 (disabled) 0.05 is good est.
        -hcerr2   [d]     - Maximum Bonferroni Corrected probability of distance error in childs (rho1-rho2)<(rho1_hat-rho2_hat)/2. Defaults to 1.0 (disabled) 0.05 is good est.
        -hcmndinf [d]     - If >0 use expected directedInformation to determine ordering: E[DeltaInf(rho1,rho2)&(rho1>rho2)] Instead of point estimate DeltaInf(rho1_hat,rho2_hat). Value is min inf.
        -maxdepth [u]     - Terminate recursion no later than this depth. If 0, program does nothing and exits.
        -restart  [s]     - Restart recursion from previous covtree.out file in OutBase. Creats Hi/Lo sub dires. String Vals = Hi, LO or PR (both) children.
        -dwrk     [s]     - Write all physical files in this directory, however recursion logically proceeds on dirbase. Shell muse move files, so only useful when erecurse = 2.
        -debug    [d,[s]] - Debugging options d=1-skip singelton,2-skip partiton,3-skip both,0-skip none; s=process only single chromosme with tag s.
        -sarg     s       - If s=- print serialized args and exit, else deseralize next arg as arglist.
fi



