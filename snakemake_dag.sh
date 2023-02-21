#!/bin/bash
# Note: Checkpoints result in inaccurate flow diagrams.

snakemake --configfile multiPrime.yaml -s multiPrime.py --dag | dot -Tpdf > dag.pdf

snakemake --configfile multiPrime2.yaml -s multiPrime2.py --dag | dot -Tpdf > dag.pdf

#snakemake --dag report.html | dot -Tsvg > dag.svg
#snakemake -s *** --cores 10
#snakemake --cores 10(-j 10)

