#!/bin/bash

if [[ "$dbase" == "" ]]; then dbase=".."; fi

dbase=`readlink -f "$dbase"`
dbin="${dbase}/bin"

egzu="gunzip"
ebzu="${dbin}/bunzip2"
estu="${dbin}/unstarch"

#insight2 database directory
didb="${dbase}/Insight2/db"
fiex="${dbase}/Insight2/bin/Insight2"

# expanda compressed file, while keeping the original compressed version
function expandFile() {
  local fname="$1"
  local suf="${fname##*.}"
  local fout="${fname%.${suf}}"
  if [[ ! -s "$fout" ]]; then
    mv "${fname}"   "${fname}.1"
    cp "${fname}.1" "${fname}"
    case "$suf" in
      z|Z|gz|bgz)                     "${egzu}" "$fname"  ;;
      bz2)                            "${ebzu}" "$fname"  ;;
      starch)                         unstarch "$fname" > "$fout" ;;
      # bw|bigwig|bigWig|BigWig)        bigWigToBedGraph "$fname" /dev/stdout;;
      # bb)                             bigBedToBed "$fname" /dev/stdout;;
      *)         cat "$fname" > "$fout" ;;
    esac
    mv -f "${fname}.1" "${fname}"
  fi
}

# is there an exe file in the current path or corrent directory?
function exeFound() {
  local fb="$1"
  local fbp=`which "$fb" 2>&1 | grep -v "no ${fb}"`
  if [[ "$fbp" == "" ]]; then echo "0"; else echo "1"; fi
}


# If therre is no executable i nthe currentdirectory, link the the system onem or expand a local archive
function exeSet() {
  local fb="$1"
  local bok=`exeFound "${fb}"`

  if [[ -x "./${fb}" ]]; then echo "1"; return; fi

  if [[ `exeFound "${fb}"` == "1" ]]; then
    ln -s  `which "$fb" 2>&1` "./${fb}"
  else
    expandFile "./${fb}.gz"
  fi

  if [[ -x "./${fb}" ]]; then echo "1"; return; fi
  echo "0"
}


