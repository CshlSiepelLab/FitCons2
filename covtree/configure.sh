#!/bin/bash

sarg="$1"
if [[ "$sarg" != "clean" ]]; then sarg="make"; fi


echo -e "`date` Configuring covtree code\n\tmode: $sarg \n\t directory `pwd` \n\t (directory `readlink -f .` )\n"

dbase=".."
source "${dbase}/bin/lib.configure.sh"


if [[ ! -x "$fiex" ]]; then
  echo -e "\ncovtree requires Insight2, checking Insight2 configuration first...."
  pushd "../Insight2"
  chmod a+x ./configure.sh
  ./configure.sh
  popd
  if [[ ! -x "$fiex" ]]; then
    echo -e "\n Insight2 Configuration failed. Unable to find:\n\t ${fiex} \n\tError. Exiting. "
    exit 1
  fi
fi


# set up binary

fbin="./bin/covtree"

if [[ "$sarg" == "clean" ]]; then
  echo -e "Removing binary $fbin ."
  rm -f "$fbin"
else 
  if [[ -x "$fbin" ]]; then
    echo -e "\t Found binary $fbin ."
  else
    echo -e "\t Attempting to build binary $fbin ."
    make
    if [[ ! -x "$fbin" ]]; then
      echo -e "\t build failed - Attempting to expand precompiled binary $fbin."
      expandFile "./bin/covtree.bz2"
    fi
  fi
  
  if [[ -x "$fbin" ]]; then
    echo -e "\t Found binary $fbin ."
  else
    echo -e "\t Unable to generate bianry. Error. Exiting."
    exit 1
  fi
fi

echo " "

  

# fetching DB
echo "`date` Fetching covtree master database. Compressed archive requires approximately 2.5GB."

ddb="./db"; mkdir -p "${ddb}"; ddb=`readlink -f "${ddb}"`
fdb="covtreeDB.full.tar"
furl='http://compgen.cshl.edu/Insight2/db/FitCons2/covtreeDB.full.tar'

if [[ "$sarg" == "clean" ]]; then
  echo -e "`date`\tClearing  database in ./db ."
  rm -vf ./db/*
  echo " "
  echo -e "`date`\tRemoving downloaded database archive $fdb ."
  rm -vf "$fdb"
  echo " "
else
  if [[ -s "$fdb" ]] ; then
    echo -e "`date`\t Found database file $fdb . "
  else
    echo -e "`date`\t Retrieving database file $fdb from $furl ."
    wget "${furl}"
  fi

  if [[ ! -s "$fdb" ]] ; then
    echo -e "`date`\t Did not find database in $fdb - build failed. Exiting."
    exit 1
  fi
  
  echo -e "`date`\t Extracting database elements to ./db ."
  tar -xvf "$fdb" -C "./db"
  echo " "

  echo -e "`date`\t To Decompressing database elements in ./db, cd to ./db and execute unroll.simple.sh ."
  echo -e "`date`\t WARNING: uncompressed databased is approximately 150GB, cache fiels created durign execution extend this to about 200GB."
  # pushd ./db
  # ./unroll.simple.sh
  # popd
  echo -e "\n`date`\t Database Extracted, but must be decompressed before use."
fi

echo " "

echo -e "\n`date`\t Action Completed sucessfully. Exiting."

