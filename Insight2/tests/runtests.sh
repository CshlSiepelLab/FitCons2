#!/bin/bash

sarg="$1"
dbase="../.."
source "${dbase}/bin/lib.configure.sh"
# Insight2db from lib
#insight2 executable from lib...
fbin="${fiex}"

echo -e "`date` Running Insight2 demos. \n\t directory `pwd` \n\t (directory `readlink -f .` )\n"

if [[ ! -x "$fbin" ]]; then
  echo -e "`date`\tUnable to find Insight2 executable ${fbin} . Please configure Insight2 first. Failure. Exiting.\n"
  exit 1
fi

echo -e "\n`date` Starting Demo 0: Command line argument list"
${fbin} -h
echo -e "\n`date` Demo 0 complete.\n\n"


echo -e "\n\n\n`date` Starting Demo 1: HG19"

dt="1-hg19"
pushd "$dt"

fin="hg19.bed"
echo -e "\t${fbin}  ${fidb} -fin $fin -qmap "
${fbin}  "${didb}" -fin "$fin" -qmap 

popd
echo -e "\n`date` Demo 1 complete.\n\n"



echo -e "\n\n\n`date` Starting Demo 2: Weighted positions (FitCons2 class ID 17,  genome-wide priors)."

dt="2-weighted"
pushd "$dt"

# expand input file if needed 
fin="FC2-017.bed"
echo -e "\t`date` Verifying input file $fin ."
expandFile "${fin}.starch"

echo -e "\t`date`\n\t${fbin}  ${fidb} -fin $fin -qmap -maxsampw 115"
${fbin}  "${didb}" -fin "$fin" -qmap -maxsampw 115

popd
echo -e "\n`date` Demo 2 complete.\n\n"



echo -e "\n\n\n`date` Starting Demo 3: Posterior and expectation - (FitCons2 class ID 17, parental priors)"
echo -e "\tNOTE: this uses 8 cores on the host computer. "

dt="3-posterior"
pushd "$dt"

# expand input file if needed 
fin="FC2-017.bed"
echo -e "\t`date` Verifying input file $fin ."
expandFile "${fin}.starch"

echo -e "\t`date`\n\t${fbin}  ${fidb} -fin $fin -qmap -maxsampw 115  -rho 0.359788,1 -eta 0.788613,1 -gamma 1.277413,1 -qexp -nthread 8 -gres 10 -expval -expdist"
${fbin}  "${didb}" -fin "$fin" -qmap -maxsampw 115 -rho 0.359788,1 -eta 0.788613,1 -gamma 1.277413,1 -qexp -nthread 8 -gres 10 -expval -expdist

popd
echo -e "\n`date` Demo 3 complete.\n\n"


echo -e "\n\nTo compare demo resuts to expected output, compare results in each directory (e.g. 1-hg19) with fiels in ref subdirectories (e.g. 1-ht19/ref).\n"

echo -e "\n`date` All Demos compelte. Exiting."

