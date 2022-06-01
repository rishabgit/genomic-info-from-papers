#!/usr/bin/env bash

#BSUB -J variant_extraction
#BSUB -q short
#BSUB -n 1
#BSUB -M 20000
#BSUB -R "rusage[mem=20GB]"
#BSUB -o /hps/scratch/flicek/wormbase/var_extraction/Output_%J.out
#BSUB -e /hps/scratch/flicek/wormbase/var_extraction/Error_%J.err


VAREXTPATH=/hps/software/users/wormbase/variant-extraction/genomic-info-from-papers
VAREXTLOGS=/hps/scratch/flicek/wormbase/var_extraction
export NLTK_DATA="${VAREXTPATH}/data/nltk"
cd "${VAREXTPATH}"
source venv/bin/activate
python process_papers.py -i "${VAREXTLOGS}/INPUTFILE" -o "${VAREXTLOGS}/INPUTFILE.csv" -m textpresso

