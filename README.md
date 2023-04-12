multiPrime: version 2.0.3

# multi PCR primer pairs design processing pipeline
MultiPrime is a pipeline designed for broad-spectrum detection of target sequences using tNGS. It is implemented in Python and Snakemake and takes a FASTA format file as input. The pipeline has three main steps: classification by identity, primer design, and primer set combination. In the classification step, redundant sequences are removed and clusters are formed by identity. Rare sequence clusters are compared to others by average nucleotide identity, and if they are deemed similar enough, they are merged. In the primer design step, multi-alignment is performed using MUSCLE or MAFFT, and candidate primers are designed using the nearest-neighbor model. Primer pairs are selected based on PCR product length, melting temperature, dimer examination, coverage with errors, and other factors. Finally, a greedy algorithm is used to combine primer pairs into a minimal primer set according to dimer examination.

multiPrime1: Degenerate primer design by DEGEPRIME (MC-DPD).

mulitPrime2: Degenerate primer design by multiPrime-core (MC-EDPD or MC-DPD).

Scripts and pipelines provided in this repository aid to design multiplex PCR primer and return a minimal primerset for multi-PCR. It contains all scripts to allow a self-assembled processing and additionally provides pipeline scripts that run the entire processing automatically.

# Requirements

To run this pipeline, your computer requires **30 GB of available memory (RAM)** to process larger number of sequence (e.g. 1,000,000). We **don't suggest** that Input sequences contains those sequences whose **average length is greater than 100K**,  if necessary, you'd better set the Maxseq in yaml file as small as possible, but **do not smaller than 200**. Snakemake was used to facilitate the automated execution of all analysis steps. The easiest way to make use of the pipeline is to set up a python 3.9 virtual environment and run the pipeline is this environment. 

Download/Provide all necessary files:

DEGEPRIME-1.1.0: DOI: 10.1128/AEM.01403-14; "DegePrime, a program for degenerate primer design for broad-taxonomic-range PCR in microbial ecology studies."
		Links: https://github.com/EnvGen/DegePrime; please move this directory into scripts.

biopython: Not required in v1.0.1 and the subsequent version.

mfeprimer-3.2.6: DOI: 10.1093/nar/gkz351; Please cite: "MFEprimer-3.0: quality control for PCR primers." please move this it into scripts. 
Please add "execute" to mfeprimer-3.2.6

Programs we employed:

MUSCLE: It is already in the requirement.txt. version=v3.8.1551. http://www.drive5.com/muscle This software is donated to the public domain. Please cite: Edgar, R.C. Nucleic Acids Res 32(5), 1792-97.

MAFFT: It is already in the requirement.txt. version=v7.508 (2022/Sep/07). Please cite: "MAFFT multiple sequence alignment software version 7: improvements in performance and usability".

fastANI: It is already in the requirement.txt. version=version 1.33. Please cite: "FastANI, Mash and Dashing equally differentiate between Klebsiella species." 

blast+: It is already in the requirement.txt. version=BLAST 2.13.0+. Links: https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastNews.

bowtie2: It is already in the requirement.txt. version=version 2.2.5. DOI:10.1038/nmeth.1923; Please cite: "Fast gapped-read alignment with Bowtie 2."
		Links: https://www.nature.com/articles/nmeth.1923
# snakemake
Snakemake is a workflow management system that helps to create and execute data processing pipelines. It requires python3 and dependent environment (multiPrime == multiPrime2) can be most easily installed via the bioconda package of the python anaconda distribution. 

  ```bash
  conda create -n multiPrime -c bioconda -c conda-forge --file requirement.txt
  ```
  if conflict:
  ```bash
  conda create -n multiPrime -c bioconda -c conda-forge --file requirement2.txt
  ```

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

The working directory contains files named `multiPrime.yaml` and `multiPrime2.yaml`. These are the central file in which all user settings, paramter values and path specifications are stored. `multiPrime.yaml` employs DEGEPRIME-1.1.0 for maximum coverage degenerate primer design (MC-DPD), `multiPrime2.yaml` use multiPrime-core.py for MC-DPD or MC-DPD with error. During a run, all steps of the pipeline will retrieve their paramter values from these file. It follows the yaml syntax (find more information about yaml and it's syntax [here](http://www.yaml.org/)) what makes it easy to read and edit. The main principles are:
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
  - **identity** - threshold for classification. please note: If you set 1, multiPrime will design candidate primer pairs for each fasta in input files. Suggestion: 0.7-0.8. 

# Start a run

Once you set up your configuration file, running the pipeline locally on your computer is as easy as invoking:
  ```bash
  sh run.sh
  ```
  minimal degeneracy degenerate primer design (MC-DPD):
  ```bash
  snakemake --configfile multiPrime.yaml -s multiPrime.py --cores 10 --resources disk_mb=80000
  ```
  minimal degeneracy degenerate primer design with errors (MC-EDPD) or MC-DPD. It depends on the parameter in multiPrime2.yaml. MC-DPD when you set variation=0.
  ```bash
  snakemake --configfile multiPrime2.yaml -s multiPrime2.py --cores 10 --resources disk_mb=80000
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
or
  ```bash
  pip install multiPrime
  ```
  ```bash
  multiPrime --help
  ```
  multiPrime package contains primer design; primer pair selection and primer pair coverage statistics, more functions will be improved in the future. All manual instruction for multiPrime can be found in [here](https://pypi.org/project/multiPrime/).

For example:
  
  MC-DPD (--variation 0) or MC-EDPD (--variation 1 or 2, we do not suggest you set --variation greater than 2, because amplification efficiency was severely inhibited when there are 3 mismathes). 

  ```bash
  python scripts/multiPrime-core.py
  ```
  ```bash
  Usage: multiPrime-core.py -i input -o output -p 20
                 Options: { -l [18] -n [4] -d [10] -v [1] -e [3.6] -g [0.2,0.7] -f [0.8] -c [4] -p [10] -a [4] }

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
  -e ENTROPY, --entropy=ENTROPY
                        Entropy is actually a measure of disorder. This
                        parameter is used to judge the window is conservation or not. 
			Entropy of primer-length window. Default: 3.6.
  -g GC, --gc=GC        Filter primers by GC content. Default [0.2,0.7].
  -s SIZE, --size=SIZE  Filter primers by mini PRODUCT size. Default 100.
  -f FRACTION, --fraction=FRACTION
                        Filter primers by match fraction. Default: 0.8.
  -c COORDINATE, --coordinate=COORDINATE
                        Mismatch index is not allowed to locate in start or
                        stop region. otherwise, it won't be regard as the mis-
                        coverage. With this param, you can control the index
                        of Y-distance (number and position of mismatch) when calculate
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
  Get candidate degenerate primer with high (error) coverage:
  ```bash
  python scripts/get_multiPrime.py
  ```
  ```bash
  Usage: get_multiPrime.py -i input -r sequence.fa -o output
                 Options: {-f [0.6] -m [500] -n [200] -e [4] -p [9] -s [250,500] -g [0.2,0.7] -d [4] -a ","}.

  Options:
  -h, --help            show this help message and exit
  -i INPUT, --input=INPUT
                        Input file: degeprimer out.
  -r REF, --ref=REF     Reference sequence file: all the sequence in 1 fasta,
                        for example: (Cluster_96_171.fa).
  -g GC, --gc=GC        Filter primers by GC content. Default [0.2,0.7].
  -f FRACTION, --fraction=FRACTION
                        Filter primers by match fraction. Default: 0.6.
                        Sometimes you need a small fraction to get output.
  -e END, --end=END     Filter primers by degenerate base position. e.g. [-t
                        4] means I dont want degenerate base appear at the end
                        four bases when primer pre-filter. Default: 4.
  -p PROC, --proc=PROC  Number of process to launch.  default: 10.
  -s SIZE, --size=SIZE  Filter primers by PRODUCT size. Default [250,500].
  -d DIST, --dist=DIST  Filter param of hairpin, which means distance of the
                        minimal paired bases. Default: 4. Example:(number of
                        X) AGCT[XXXX]AGCT.
  -a ADAPTOR, --adaptor=ADAPTOR
                        Adaptor sequence, which is used for NGS next. Hairpin
                        or dimer detection for [adaptor--primer]. For example: 
                        TCTTTCCCTACACGACGCTCTTCCGATCT,TCTTTCCCTACACGACGCTCTTCCGATCT 
                        (Default). If you dont want adaptor,
                        use [","]
  -m MAXSEQ, --maxseq=MAXSEQ
                        Limit of sequence number. Default: 500. If 0, then all
                        sequence will take into account. This param should
                        consistent with [max_seq] in multi-alignment [muscle].
  -o OUT, --out=OUT     Output file: candidate primers. e.g.
                        [*].candidate.primers.txt.
  ```
  Extract primers from ONT reads. FindONTprimerV2.py = FindONTprimerV3.py:
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
  python scripts/primer_coverage_validation_by_BWT.py
  ```
  ```
  Usage: primer_coverage_validation_by_BWT.py -i [input] -r [bowtie index] -l [150,2000] -p [10]-o [output]

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
  Others ...
  ```

# Output
logs: log file of the multiPrime.py 

results: results directory

	-cluster.identities: identity of each sequence.

	-cluster.txt: cluster information. for example: Cluster_0_222.fa, 0 ==> cluster rank; 222 ==> sequence number.

	-history.txt:  history of clusters with rare sequence numbers are compared with others by average nucleotide identity.
	
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
		--*.out: paired primers designed by the top N {default: 500} fasta
		--*.gap_seq_id_json: Positions and non-contained sequences caused by gap.
		--*.non_coverage_seq_id_json: Positions and non-contained sequences caused by others.

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
		--Coverage_stast.xls: coverage of all primers in primer set (perfect match)
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
		--core_Coverage_stast.xls: coverage of all primers in core primer set (perfect match)
		--core_candidate_primers_sets.fa: core primer set fasta
		--core_candidate_primers_sets.number: candidate primer number of each core cluster
		--core_final_maxprimers_set.fa.dimer: dimer check by mfePrimer
		--core_final_maxprimers_set.fa.hairpin: hairpin check by mfePrimer
		--core_PCR_product: core PCR product of each primer


# Contact
Please send comments, suggestions, bug reports and bug fixes to 1806389316@pku.edu.cn

# Todo
This repository is updated frequently. Expect breaking changes. More functions will be improved in the future and a new version is comming.


