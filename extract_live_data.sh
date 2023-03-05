#!/bin/bash

#
# Purpose: Extract data from 'macroscopic_vars.dat' which is written by JOREK during the code run.
#
# Date: 2011-03-30
# Author: Matthias Hoelzl, IPP Garching
#

function usage() {
  echo ""
  echo "Usage: `basename $0` <PREFIX> [<target_file>] [options]"
  echo ""
  echo "    -l          List all available prefixes"
  echo "    -f <file>   Use <file> instead of ./macroscopic_vars.dat"
  echo ""
  echo "Remark:"
  echo "  If <target_file> is omitted, the output goes to STDOUT"
  echo ""
  echo "Examples:"
  echo "  extract_live_data energies energies.dat"
  echo "  extract_live_data n_tor"
  echo ""
}

# --- Evaluate command line parameters
if [ $# -lt 1 ] || [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
  usage
  exit
fi

PREFIX=""
if [ "${1:0:1}" != "-" ]; then
  PREFIX="$1"
  shift
fi

TARGET=""
if [ "${1:0:1}" != "-" ]; then
  TARGET="$1"
  shift
fi

list=0
file="./macroscopic_vars.dat"
while [ $# -gt 0 ]; do
  if [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
    usage
    exit
  elif [ "$1" == "-l" ]; then
    list=1
    shift
  elif [ "$1" == "-f" ]; then
    file="$2"
    shift 2
  else
    echo ""
    echo "ERROR: Unkown option '$1'."
    usage
    exit 1
  fi
done

# --- List all available prefixes if -l option is specified
if [ $list -eq 1 ]; then
  echo ""
  tmp="/tmp/tmp_eld_$$"
  cat $file | sed -e "s/^@\([^ ]*\):.*/\1/" | sort | uniq -c | sed -e 's/^ *\([0-9]*\) *\([^ ]*\)/\2: \1 entries/'
  echo
  exit
fi

# --- Extract all data for the specified prefix
if [ "$PREFIX" == "" ]; then
  echo ""
  echo "ERROR: No prefix specified."
  usage
  exit 1
fi

function extract() { grep "@${PREFIX}:" $file | sed -e "s/^@${PREFIX}://" -e "s/^ *//" -e "s/ *$//"; }

if [ "$TARGET" == "-" ] || [ "$TARGET" == "" ]; then
  extract
else
  extract > ${TARGET}
fi
