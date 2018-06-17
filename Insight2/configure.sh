#!/bin/bash

sarg="$1"
if [[ "$sarg" != "clean" ]]; then sarg="make"; fi

echo -e "`date` Configuring Insight2 code\n\tmode: $sarg \n\t directory `pwd` \n\t (directory `readlink -f .` )\n"

dbase=".."
source "${dbase}/bin/lib.configure.sh"

if [[ ! -x "${dbase}/bin/unstarch" ]]; then
  echo -e "\n`date` Insight2 requires root configuration, configuring root."
  pushd "../bin"
  ./configure.sh
  popd
  if [[ ! -x "${dbase}/bin/unstarch" ]]; then
    echo -e "\n`date` root configuration failed. Error. Exiting."
    exit 1
  fi
fi

# set up binary

fbin="./bin/Insight2"

if [[ "$sarg" == "clean" ]]; then
  echo -e "Removing binary $fbin ."
  rm -f "$fbin"
else 
  if [[ -x "$fbin" ]]; then
    echo -e "\t Found binary $fbin ."
  else
    echo -e "\t Attempting to buind binary $fbin ."
    make
    if [[ ! -x "$fbin" ]]; then
      echo -e "\t build failed - Attempting to expand precompiled binary $fbin."
      expandFile "./bin/Insight2.bz2"
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

  

# fetchign DB
echo "`date` Configuring Insight2 database."

ddb="./db"; mkdir -p "${ddb}"; ddb=`readlink -f "${ddb}"`
fdb="Insight2DB.tar"
furl='http://compgen.cshl.edu/~bgulko/research/Insight2/db/FitCons2/Insight2DB.tar'

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

  echo -e "`date`\t Decompressing database elements in ./db ."
  flist="block.bedg.starch monoDB.db.gz poly.bedg.starch polyn.bedg.starch"
  for fin in ${flist}; do
    echo -e "`date`\t\t Expanding $fin "
    expandFile "./db/${fin}"
    fout="./db/${fin%.*}" 
    if [[ ! -s "$fout" ]]; then
      echo -e "`date`\t\t unable to find $fout . Error. Exiting. "
      exit 1
    fi
  done

  echo " "
  ls -l ./db/*
  echo -e "\n`date`\t Database present. "
fi

echo " "

echo -e "\n`date`\t Action Completed sucessfully. Exiting."

