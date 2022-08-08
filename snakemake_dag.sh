#!/bin/bash
snakemake --configfile multi_PCR_primer_design_pipeline.yaml -s multi_PCR_primer_design_pipeline.py --dag | dot -Tpdf > dag.pdf

#snakemake --dag report.html | dot -Tsvg > dag.svg
#snakemake -s *** --cores 10
#snakemake --cores 10(-j 10)

