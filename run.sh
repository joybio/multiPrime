#!/bin/bsh
#snakemake --configfile multimulti-DegePrime.yaml -s multi-DegePrime.py --cores 10 --resources mem_mb=80000

#snakemake --configfile multiPrime-orinal.yaml -s multiPrime-orignal.py --cores 10 --resources mem_mb=80000

snakemake --configfile multiPrime.yaml -s multiPrime.py --cores 10 --resources mem_mb=80000

