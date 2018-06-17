#!/bin/bash

echo -e "\n\n`date` FitCons2 Master configuration file, configures top level bin dire, then Insight2 then covtree."

db="."; db=`readlink -f "${db}"`
echo -e "\n`date`\tConfiguration base directory:\n\t\t${db}\n"
echo -e "\n`date`\t Configuring ./bin "
pushd bin
./configure.sh
popd

echo -e "\n\n`date`\t Configuring ./Insight2 "
pushd ./Insight2
./configure.sh
popd

echo -e "\n\n`date`\t Configuring ./covtree "
pushd ./covtree
./configure.sh
popd


echo -e "\n\n`date`\t Configuration Complete."
echo -e "\t To run Insight2 demos\n\t\t cd ./Insight2/tests\n\t\t./runtests.sh\n\n"
echo -e "\t To run covtree  demos\n\t\t cd ./covtree/tests\n\t\t./runtest1.sh\n\t\truntest2.sh\n\n"

echo -e "Covtree demos download and operate on a subset of the FitCons2 database and may require 20GB of RAM, 100GB of storage and take an hour on three cores."
echo -e "Full FitCons2 database (2.5GB) is downloadded when covtree is configured but requries approximately 150GB of storage to decompress. This is not required for demos."

echo -e "\n\n`date` Done, exiting."
