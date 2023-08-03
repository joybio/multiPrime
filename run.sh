#!/bin/bsh
#snakemake --configfile multi-DegePrime.yaml -s multi-DegePrime.py --cores 10 --resources mem_mb=80000

#snakemake --configfile multiPrime-original.yaml -s multiPrime-original.py --cores 10 --resources mem_mb=80000

snakemake --configfile multiPrime.yaml -s multiPrime.py --cores 10 --resources mem_mb=80000

