multiPrime: version 0.0.2
# multi PCR primers design processing pipeline

Scripts and pipelines provided in this repository aid to design multiplex PCR primer and return a maximal primerset for multi-PCR. It contains all scripts to allow a self-assembled processing and additionally provides pipeline scripts that run the entire processing automatically.

# Requirements

To run this pipeline, your computer requires **40 GB of available memory (RAM)** to process larger genomes (e.g. human or mouse). Moreover, snakemake was used to facilitate the automated execution of all analysis steps. The easiest way to make use of the pipeline is to set up a python3 virtual environment and run the pipeline is this environment. 
Download/Provide all necessary files:

DEGEPRIME-1.1.0: DOI: 10.1128/AEM.01403-14; "DegePrime, a program for degenerate primer design for broad-taxonomic-range PCR in microbial ecology studies"
		Links: https://github.com/EnvGen/DegePrime; please move this directory into scripts

biopython

mfeprimer-3.2.6: DOI: 10.1093/nar/gkz351; "MFEprimer-3.0: quality control for PCR primers."

blast+: Links: https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastNews

# snakemake
Snakemake is a workflow management system that helps to create and execute data processing pipelines. It requires python3 and can be most easily installed via the bioconda package of the python anaconda distribution.

conda create -n multiPCR-snakemake -c bioconda -c conda-forge --file requirement.txt

# Activate the environment
  ```bash
  source activate multiPCR-snakemake
  ```
  To exit the environment (after finishing the usage of the pipeline), just execute
  ```bash
  source deactivate
  ```
# Run the pipeline

# Configure input parameters

The working directory contains a file named `multiPrime.yaml`. It`s the central file in which all user settings, paramter values and path specifications are stored. During a run, all steps of the pipeline will retrieve their paramter values from this file. It follows the yaml syntax (find more information about yaml and it's syntax [here](http://www.yaml.org/)) what makes it easy to read and edit. The main principles are:
  - everything that comes after a `#` symbol is considered as comment and will not be interpreted
  - paramters are given as key-value pair, with `key` being the name and `value` the value of any paramter

Before starting the pipeline, open the `multiPrime.yaml` configuration file and set all options according as required. This should at least include:
  - **name of the input directory** - where are your input fastq files stored
	-input_dir: ["abs_path_to"]/multiPrime/test_data
  - **name of the output directory** - where should the pipeline store the output files (the direcotry is created if not existing)
	-results_dir: ["abs_path_to"]/multiPrime/test_data/results
  - **name of the log directory** - where should the pipeline store the log files
	-log_dir: ["abs_path_to"]/multiPrime/test_data/logs
  - **name of the scripts directory** - where should the pipeline store the scripts files
	-scripts_dir: ["abs_path_to"]/multiPrime/scripts
  - **name(s) of your input samples** - please note: If your sample is named `sample1.fa` then `sample1` will be kept as naming scheme throughout the entire run to indicate output files that belong to this input file, e.g. the pipeline will create a file called `sample1.3pSites.noIP.bed.gz`. If you have multiple input files, just follow the given pattern with one sample name per line (and a dash that indicates another list item).
	-virus:
		- test

# Start a run

Once you set up your configuration file, running the pipeline locally on your computer is as easy as invoking:
`sh run.sh`

# Output
logs: log file of the multiPrime.py 

results: results directory

	-Total_fa: genome file and cluster of genome file.

	-Clusters_fa: genome file split by each cluster.
		--*.fa: fasta of each cluster
		--*.tfa: top N {default: 500} fasta of each cluster
		--*.txt: accession id of each cluster
		--*.db: directory of database (for blastn).
		--*.blastout: output file of blastn of all paired primers.
		--*.number: number of fasta in each cluster

	-Clusters_msa: alginment by muscle
		--*.tmsa: muscle output of the top N {default: 500}

	-Clusters_trim_msa: trimmed alignment by degePrimer
		--*.trim.tmsa: trimmed muscle by degePrimer

	-Clusters_primer: get_degePrimer from degePrimer out
		--*.top.primer.out: paired primers designed by the top N {default: 500} fasta

	-Clusters_cprimer:
		--*.fa: candidate primers in fa format
		--*.txt: candidate primers in txt format (1 line)
		--*.nt.Check: tmp file; primers filter by blastn

	-Primers_set
		--candidate_primers_sets.txt: all candidate primers in each cluster
		--sort.candidate_primers_sets.txt: sort by the number of candidate primers in each line (cluster)
		--final_maxprimers_set.fa: fasta format of primers set for multiPCR
		--final_maxprimers_set.xls: primers information of primers set
		--Coverage_stast.xls: coverage of all primers in primer set
		--final_maxprimers_set.fa.dimer: dimer check by mfePrimer 
		--final_maxprimers_set.fa.hairpin: hairpin check by mfePrimer
		--PCR_product: PCR product for each primer


