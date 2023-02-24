#!/bin/bsh
#snakemake --configfile multiPrime.yaml -s multiPrime.py --cores 10 --resources mem_mb=80000

snakemake --configfile multiPrime2.yaml -s multiPrime2.py --cores 10 --resources mem_mb=80000


