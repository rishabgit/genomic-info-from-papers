#!/usr/bin/env bash

if [[ $# -ne 1 ]]; then
    echo "Illegal number of parameters" >&2
    exit 2
fi
FILENAME=$1
VAREXTPATH=/hps/software/users/wormbase/variant-extraction/genomic-info-from-papers
VAREXTLOGS=/hps/scratch/flicek/wormbase/var_extraction
cp -f "${VAREXTPATH}/codon_launch_template.sh" "${VAREXTLOGS}/codon_launch.sh"
sed -i "s/INPUTFILE/${FILENAME}/g" "${VAREXTLOGS}/codon_launch.sh"

