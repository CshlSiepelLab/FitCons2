#!/bin/bash

source "lib.configure.sh"

scln="$1"

if [[ `exeFound "$egzu"` == "0" ]]; then
  echo -p "\n`date` Unable to find ${egzu}. This program is required for configuration. Exiting.\n"
  exit 1
fi

lparts="bzip2 starch unstarch bunzip2"
for fb in ${lparts} ; do
  if [[ "$scln" == "clean" ]]; then
    echo -ne "`date` cleaning  ${fb}.\n"
    rm -f "./${fb}"
  else 
    echo -ne "`date` Configuring ${fb}..."
    if [[ -x "./${fb}" ]]; then echo "found."; continue; fi
    ok=`exeSet "${fb}"`
    if [[ "$ok" == "1" ]]; then 
      echo "ok."
    else
      echo "Unable to configure. Exiting."
      exit 1
    fi
  fi
done

echo -e "`date` Sucessfully Configured.\n"



