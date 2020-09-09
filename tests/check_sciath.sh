#!/usr/bin/env sh

if ! python -c"import sciath"
then
  printf "SciATH not found. Either update your PYTHONPATH or\n"
  printf "\n  git clone https://github.com/sciath/sciath tests/sciath\n\n"
  exit 1
fi
