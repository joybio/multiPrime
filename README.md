multiPrime: version 2.0.0
# multi PCR primers design processing pipeline
multiPrime was developed to design primers for large-scale sequences, identifying the core primers and optimized primer sets.

Scripts and pipelines provided in this repository aid to design multiplex PCR primer and return a minimal primerset for multi-PCR. It contains all scripts to allow a self-assembled processing and additionally provides pipeline scripts that run the entire processing automatically.

# Requirements

To run this pipeline, your computer requires **40 GB of available memory (RAM)** to process larger number of sequence (e.g. 100000). Moreover, snakemake was used to facilitate the automated execution of all analysis steps. The easiest way to make use of the pipeline is to set up a python 3.9 virtual environment and run the pipeline is this environment. 
Download/Provide all necessary files:

DEGEPRIME-1.1.0: DOI: 10.1128/AEM.01403-14; "DegePrime, a program for degenerate primer design for broad-taxonomic-range PCR in microbial ecology studies."
		Links: https://github.com/EnvGen/DegePrime; please move this directory into scripts

biopython: Not required in v1.0.1

mfeprimer-3.2.6: DOI: 10.1093/nar/gkz351; "MFEprimer-3.0: quality control for PCR primers."

blast+: Links: https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastNews

bowtie2: DOI:10.1038/nmeth.1923; "Fast gapped-read alignment with Bowtie 2."
		Links: https://www.nature.com/articles/nmeth.1923
# snakemake
Snakemake is a workflow management system that helps to create and execute data processing pipelines. It requires python3 and can be most easily installed via the bioconda package of the python anaconda distribution.

conda create -n multiPrime -c bioconda -c conda-forge --file requirement.txt

# Activate the environment
  ```bash
  source activate multiPrime
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
  - **name of the input directory** - where are your input fasta files stored
	-input_dir: ["abs_path_to_input_dir"]
  - **name of the output directory** - where should the pipeline store the output files (the direcotry is created if not existing)
	-results_dir: ["abs_path_to_results_dir"]
  - **name of the log directory** - where should the pipeline store the log files
	-log_dir: ["abs_path_to_log_dir"]
  - **name of the scripts directory** - where should the pipeline store the scripts files
	-scripts_dir: ["abs_path_to"]/multiPrime/scripts
  - **name(s) of your input samples** - please note: If your sample is named `sample1.fa` then `sample1` will be kept as naming scheme throughout the entire run to indicate output files that belong to this input file, e.g. the pipeline will create a file called `sample1.fa`. If you have multiple input files, just follow the given pattern with one sample name per line (and a dash that indicates another list item).
	-virus:
		- test

# Start a run

Once you set up your configuration file, running the pipeline locally on your computer is as easy as invoking:
  ```bash
  sh run.sh
  ```

# Start a run independently
If you want to run python script locally on your computer independently. It is as easy as invoking:
  ```bash
  python {path to script}/{target}.py --help
  ```
or 
  ```bash
  python {path to script}/{target}.py
  ```
For example:
  
  DPD (degenerate primer design):
  ```bash
  python scripts/multiPrime-core.py
  ```
  ```bash
  Usage: multiPrime-core.py -i input -o output -p 20
                 Options: { -l [18] -n [4] -d [10] -v [1] -g [0.2,0.7] -f [0.8] -c [4] -p [10] -a [4] }

  Options:
  -h, --help            show this help message and exit
  -i INPUT, --input=INPUT
                        Input file: multi-alignment output (muscle or others).
  -l PLEN, --plen=PLEN  Length of primer. Default: 18.
  -n DNUM, --dnum=DNUM  Number of degenerate. Default: 4.
  -d DEGENERACY, --degeneracy=DEGENERACY
                        degeneracy of primer. Default: 10.
  -v VARIATION, --variation=VARIATION
                        Max mismatch number of primer. Default: 1.
  -g GC, --gc=GC        Filter primers by GC content. Default [0.2,0.7].
  -s SIZE, --size=SIZE  Filter primers by mini PRODUCT size. Default 100.
  -f FRACTION, --fraction=FRACTION
                        Filter primers by match fraction. Default: 0.8.
  -c COORDINATE, --coordinate=COORDINATE
                        Mismatch index is not allowed to locate in start or
                        stop. otherwise, it won't be regard as the mis-
                        coverage. With this param, you can control the index
                        of Y-distance (number=variation and position of mismatch) when calculate
                        coverage with error.Default: 4.
  -p PROC, --proc=PROC  Number of process to launch. Default: 20.
  -a AWAY, --away=AWAY  Filter hairpin structure, which means distance of the
                        minimal paired bases. Default: 4. Example:(number of
                        X) AGCT[XXXX]AGCT. Primers should not have
                        complementary sequences (no consecutive 4 bp
                        complementarities),otherwise the primers themselves
                        will fold into hairpin structure.
  -o OUT, --out=OUT     Output file: candidate primers. e.g.
                        [*].candidate.primers.txt.
  ```
  Extract primers from ONT reads. Input file format can be "fasta", "fastq", "fa.gz", and "fq.gz":
  ```bash
  python scripts/FindONTprimerV3.py
  ```
  ```bash
  Usage: FindONTprimerV3.py -i [input] -s [primer set] -p [20] -l [primer length] -m [0.6] -f [fq] -o [output].

  Options:
  -h, --help            show this help message and exit
  -i INPUT, --input=INPUT
                        Input file: fastq or fasta or fq.gz or fa.gz.
  -s SET, --set=SET     primer set file.
  -p NPROC, --nproc=NPROC
                        Primer set file. option. Default: 10
  -l LEN, --len=LEN     Primer length. Default: 18
  -m MIN_IDENT, --min_ident=MIN_IDENT
                        min identity. Default: 0.6
  -f FORMAT, --format=FORMAT
                        Input format can be fasta, fastq, fa.gz and fq.gz. Default: fastq
  -o OUT, --out=OUT     Output file: candidate primers. e.g.
                        [*].candidate.primers.txt.
  ```
  Extract PCR products with perfect match:
  ```bash
  python scripts/extract_PCR_product.py
  ```
  ```
  Usage: extract_PCR_product.py -i [input] -r [reference] -p [10] -f [format] -o [output]

  Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -r REF, --ref=REF     reference file: template fasta or reference fasta.
  -i INPUT, --input=INPUT
                        Primer file. One of the followed three types:
                        final_maxprimers_set.xls   primer.fa
                        primer_F,primer_R.
  -f FORMAT, --format=FORMAT
                        Format of primer file: xls or fa or seq; default: xls.
                        xls: final_primer_set.xls, output of multiPrime.  fa:
                        fasta format.  seq: sequence format, comma seperate.
                        e.g. primer_F,Primer_R.
  -o OUT, --out=OUT     Output_dir. default: PCR_product.
  -p PROCESS, --process=PROCESS
                        Number of process to launch.  default: 10.
  -s STAST, --stast=STAST
                        Stast information: number of coverage and total.
                        default: Coverage.xls
  ```
  Get primer information of PCR products with mismatch
  ```bash
  python scripts/primer_coverage_confirmation_by_BWT_V3.py
  ```
  ```
  Usage: primer_coverage_confirmation_by_BWT_V3.py -i [input] -r [bowtie index] -l [150,2000] -p [10]-o [output]

  Options:
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input=INPUT_FILE
                        input file: primer.fa.
  -r REF, --ref=REF     reference file: bowtie index.
  -l LEN, --len=LEN     Length of primer, which is used for mapping. Default:
                        18
  -t TERM, --term=TERM  Position of mismatch is not allowed in the 3 term of
                        primer. Default: 4
  -s SIZE, --s=SIZE     Length of PCR product, default: 150,2000.
  -p PROC, --proc=PROC  Number of process. Default: 20
  -b BOWTIE, --bowtie=BOWTIE
                        bowtie or bowtie2 was employed for mapping. Default:
                        bowtie2
  -m SEEDMMS, --seedmms=SEEDMMS
                        Bowtie: Mismatches in seed (can be 0 - 3, default: -n
                        1).Bowtie2: Gap or mismatches in seed (can be 0 - 1,
                        default: -n 1).
  -o OUT, --out=OUT     Prodcut of PCR product with primers.

  ```
  ```
  ...
  ```

# Output
logs: log file of the multiPrime.py 

results: results directory

	-cluster.identities: identity of each sequence.

	-cluster.txt: cluster information. for example: Cluster_0_222.fa, 0 ==> cluster rank; 222 ==> sequence number.

	-Total_fa: genome file and cluster of genome file.

	-Clusters_fa: genome file split by each cluster.
		--*.fa: fasta of each cluster
		--*.tfa: top N {default: 500 randomly selected. Always contain the representative seq} fasta of each cluster
		--*.txt: accession id of each cluster
		--*.db: directory of database (for bowtie2).
		--*.number: number of fasta in each cluster

	-Clusters_msa: alginment by muscle
		--*.tmsa: muscle output of the top N {default: 500 randomly selected. Always contain the representative seq}

	-Clusters_trim_msa: trimmed alignment by degePrimer
		--*.trim.tmsa: trimmed muscle by degePrimer

	-Clusters_primer: get_degePrimer from degePrimer out
		--*.top.primer.out: paired primers designed by the top N {default: 500} fasta

	-Clusters_cprimer: candidate primers for each cluster.
		--*.bed: candidate PCR product (1 mismatch and mismatch position must 9bp away from 3'end at least.)
		--*.fa: candidate primers in fa format
		--*.txt: candidate primers in txt format (1 line)
		--*.Check: tmp file; primers filter by bowtie2 (1 mismatch and mismatch position must 9bp away from 3'end at least.)

	-Primers_set:
		--candidate_primers_sets.txt: all candidate primers in each cluster
		--candidate_primers_sets: directory contain all candidate primers in fasta
		--sort.candidate_primers_sets.txt: sorted by the number of candidate primers in each line (cluster)
		--final_maxprimers_set.fa: fasta format of primers set for multiPCR
		--final_maxprimers_set.xls: primers information of primers set
		--final_maxprimers_set.next.xls: primer set 2
		--Coverage_stast.xls: coverage of all primers in primer set
		--final_maxprimers_set.fa.dimer: dimer check by mfePrimer 
		--final_maxprimers_set.fa.hairpin: hairpin check by mfePrimer
		--PCR_product: perfect PCR product of each primer

	-Core_primers_set:
		--core_candidate_primers_sets.txt: core candidate primers in each cluster
		--core_candidate_primers_sets:  directory contain all core candidate primers in fasta
		--sort.core_candidate_primers_sets.txt: sorted by the number of core candidate primers in each line (cluster)
		--core_final_maxprimers_set.fa: fasta format of core primers set for multiPCR
		--core_final_maxprimers_set.xls: primers information of core primers set
		--core_final_maxprimers_set.next.xls: primer set 2
		--core_Coverage_stast.xls: coverage of all primers in core primer set
		--core_candidate_primers_sets.fa: core primer set fasta
		--core_candidate_primers_sets.number: candidate primer number of each core cluster
		--core_final_maxprimers_set.fa.dimer: dimer check by mfePrimer
		--core_final_maxprimers_set.fa.hairpin: hairpin check by mfePrimer
		--core_PCR_product: core PCR product of each primer



