multiPrime: version 2.1.1

# Multi PCR primer pairs design processing pipeline

MultiPrime [https://multiPrime.cn](https://multiPrime.cn) is a pipeline designed for broad-spectrum detection of target sequences using tNGS. It is implemented in Python and Snakemake and takes a FASTA format file as input. The pipeline has three main steps: classification by identity, primer design, and primer set combination. In the classification step, redundant sequences are removed and clusters are formed by identity. Rare sequence clusters are compared to others by average nucleotide identity, and if they are deemed similar enough, they are merged. In the primer design step, multi-alignment is performed using MUSCLE or MAFFT, and candidate primers are designed using the nearest-neighbor model. Primer pairs are selected based on PCR product length, melting temperature, dimer examination, coverage with errors, and other factors. Finally, a greedy algorithm is used to combine primer pairs into a minimal primer set according to dimer examination.

If you only require primer design without the need for primer set combination, you may use the primer design module of MultiPrime, which is accessible through scripts/multiPrime-core.py or pip install multiPrime (version >=2.4.8) and utilize the DPrime function.

multi-DegePrime: Degenerate primer design by DEGEPRIME (MC-DPD).

multiPrime-original: Degenerate primer design by multiPrime-core (MC-EDPD or MC-DPD). It allows for avoidance of mismatches at 3'end region.

multiPrime: It is an update of original. New multiPrime allows for easy avoidance of mismatches at any position, making it flexible for experienced users.

Scripts and pipelines provided in this repository aid to design multiplex PCR primer and return a minimal primerset for multi-PCR. It contains all scripts to allow a self-assembled processing and additionally provides pipeline scripts that run the entire processing automatically.

We have provided a video tutorial at [here](https://figshare.com/articles/media/Installation_video_of_multiPrime/23904159) to assist you with the installation and usage of multiPrime.

# Why multiPrime

1) MultiPrime is a user-friendly and one-step tool for designing tNGS primer sets (cluster-specific primers or ultra multiplex PCR).
2) It integrates degenerate primer design theory with mismatch handling, resulting in improved accuracy and specificity in detecting broad spectrum sequences.
3) It outperformed conventional programs in terms of run time, primer number, and primer coverage.
4) The versatility and potential of multiPrime is highlighted by its potential application in detecting single or multiple genes, exons, antisense strands, RNA, or other specific DNA segments.

# Citation
If you find multiPrime helpful for your research or project, we kindly request that you cite the following publication:

Xia, Han et al. 2023. ["MultiPrime: A Reliable and Efficient Tool for Targeted Next-Generation Sequencing.” iMeta e143. https://doi.org/10.1002/imt2.143"](https://doi.org/10.1002/imt2.143).

Citing the publication will acknowledge the contribution of multiPrime to your work and help us in further development and improvement of the tool. Thank you for your support!

# Requirements

To run this pipeline, your computer requires **30 GB of available memory (RAM)** to process larger number of sequence (e.g. 1,000,000). **Note:** We **don't suggest** that Input sequences contains those sequences whose **length is greater than 100K**,  if necessary, you'd better set the Maxseq in yaml file as small as possible, but **do not smaller than 200**. Alternatively, you may consider using conserved genes/regions instead of whole genomes. 
Snakemake was used to facilitate the automated execution of all analysis steps. The easiest way to make use of the pipeline is to set up a python 3.9 virtual environment and run the pipeline in this environment. 

Download/Provide all necessary files:

Comparison:

DEGEPRIME-1.1.0 (multi-DegePrime): DOI: 10.1128/AEM.01403-14; Please cite: "DegePrime, a program for degenerate primer design for broad-taxonomic-range PCR in microbial ecology studies."
		Links: https://github.com/EnvGen/DegePrime; please move this directory into scripts.

mfeprimer-3.2.6: DOI: 10.1093/nar/gkz351; Please cite: "MFEprimer-3.0: quality control for PCR primers." please move this it into scripts. 
Please add "execute" to mfeprimer-3.2.6

Programs we employed:

biopython: Not required in v1.0.1 and the subsequent version.

The method for calculating Tm values in this study is a slightly modified version of primer3-py [here](https://github.com/libnano/primer3-py). Reference paper: Owczarzy et al., 2004; Owczarzy et al., 2008.

The method for calculating deltaG in this study is a slightly modified version of the approach proposed by Martin et al., 2020. "Base-Pairing and Base-Stacking Contributions to Double-Stranded DNA Formation."

The method for dimer examination in this study is a slightly modified version of the approach proposed by Xie et al., 2022. "Designing highly multiplex PCR primer sets with Simulated Annealing Design using Dimer Likelihood Estimation (SADDLE)"

MUSCLE: It is already in the requirement.txt. version=v3.8.1551. http://www.drive5.com/muscle This software is donated to the public domain. Please cite: Edgar, R.C. Nucleic Acids Res 32(5), 1792-97.

MAFFT: It is already in the requirement.txt. version=v7.508 (2022/Sep/07). Please cite: "MAFFT multiple sequence alignment software version 7: improvements in performance and usability".

fastANI: It is already in the requirement.txt. version=version 1.33. Please cite: "FastANI, Mash and Dashing equally differentiate between Klebsiella species." 

blast+: It is already in the requirement.txt. version=BLAST 2.13.0+. Links: https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastNews.

bowtie2: It is already in the requirement.txt. version=version 2.2.5. DOI:10.1038/nmeth.1923; Please cite: "Fast gapped-read alignment with Bowtie 2."
		Links: https://www.nature.com/articles/nmeth.1923

# Installation and Snakemake
Snakemake is a workflow management system that helps to create and execute data processing pipelines. It requires python3 and dependent environment (multi-DegePrime == multiPrime-original == multiPrime) can be most easily installed via the bioconda package of the python anaconda distribution. 

  ```bash
  conda create -n multiPrime -c bioconda -c conda-forge --file requirement.txt
  ```
  if conflicts:
  ```bash
  conda create -n multiPrime -c bioconda -c conda-forge --file requirement2.txt
  ```
  or
 Copy multiPrime.tar.gz files from the ENV directory to your Conda environment directory and then unpack them.
  ```
  cp ENV/multiPrime.tar.gz  ${/path/to/your/conda/env}
  tar -xzvf multiPrime.tar.gz 
  conda activate multiPrime
  ```

# Activate and exit the environment
  To activate the environment 
  ```bash
  source activate multiPrime
  ```
  To exit the environment (after finishing the usage of the pipeline), just execute
  ```bash
  source deactivate
  ```
# Run the pipeline

# Configure input parameters

The working directory contains files named `multi-DegePrime.yaml`, `multiPrime-original.yaml` and `multiPrime.yaml`. These are the central file in which all user settings, paramter values and path specifications are stored. `multi-DegePrime.yaml` employs DEGEPRIME-1.1.0 for maximum coverage degenerate primer design (MC-DPD), `multiPrime-orignal.yaml` and `multiPrime.yaml` use multiPrime-core.py for MC-DPD or MC-DPD with error. During a run, all steps of the pipeline will retrieve their paramter values from these file. It follows the yaml syntax (find more information about yaml and it's syntax [here](http://www.yaml.org/)) what makes it easy to read and edit. The main principles are:
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
  - **others** - for more information on the parameters, please refer to the YAML file.

# Start a run

Once you set up your configuration file, running the pipeline locally on your computer is as easy as invoking:
  ```bash
  sh run.sh
  ```
  maximal coverage degenerate primer design (MC-DPD). The approach employed DegePrime to design degenerate primers for the target sequence.
  ```bash
  snakemake --configfile multi-DegePrime.yaml -s multi-DegePrime.py --cores 10 --resources disk_mb=80000
  ```
  maximal coverage degenerate primer design with errors tolerant (MC-EDPD) or MC-DPD. MultiPrime-orignal is capable of avoiding mismatches that occur at the 3'end position. The approach used in multiPrime-orignal.yaml depends on the value of the "variation" parameter.

  If "variation" is set to 0, then multiPrime uses the MC-DPD approach to design degenerate primers for the target sequence. In this approach, the primer sequences are designed with prefect match (0-mismatch).

  If "variation" is set to a value greater than 0, then multiPrime uses the MC-EDPD approach to design degenerate primers with errors (mismatches) tolerance (1-mismatch or 2-mismatches). In this approach, the primer sequences are allowed to contain a limited number of errors (mismatches), which increases the probability of finding suitable primer sequences for the target sequence.
  ```bash
  snakemake --configfile multiPrime-orignal.yaml -s multiPrime-orignal.py --cores 10 --resources disk_mb=80000
  ```
  multiPrime is similiar to multiPrime-orignal, but it enables easy avoidance of mismatches at any position.
  ```bash
  snakemake --configfile multiPrime.yaml -s multiPrime.py --cores 10 --resources disk_mb=80000
  ```


# Start a run independently
Setting default parameters may not always be suitable for all conditions. If you want to design primers with more flexible parameters, you can install the multiPrime package through PyPI (Python Package Index) using pip.
  ```bash
  pip install multiPrime (make sure the version is latest)
  ```
  ```bash
  multiPrime --help
  ```
  The multiPrime package includes various functions for designing primers, selecting primer pairs, and calculating the coverage of these primer pairs. These features have been developed to assist researchers in performing PCR experiments efficiently and accurately. The package is continuously being improved and expanded upon, with additional functions expected to be added in the future. All manual instruction for multiPrime can be found in [here](https://pypi.org/project/multiPrime/).

If you have already cloned this repository and do not wish to install the multiPrime package from PyPI, you can use the scripts included in this repository to perform primer design. These scripts have the same functionality as the multiPrime package and can be run locally on your computer without the need for a separate installation.
  ```bash
  python {path to script}/{target}.py --help
  ```
or 
  ```bash
  python {path to script}/{target}.py
  ```

For example:
  
  Primer design with MC-DPD (--variation 0) or MC-EDPD (--variation 1 or 2. We do not recommend setting the --variation parameter greater than 2, as amplification efficiency can be severely inhibited when there are more than 2 mismatches between the primers and their targets). 

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
                        Entropy is a measurement of the degree of randomness or disorder in a system. 
                        This parameter is utilized to determine whether a window is conserved or not. 
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
  To get candidate degenerate primer pairs with high coverage.
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
			The default value for the limit of sequence number is set at 500. 
			However, if the value is set to 0, then all sequences will be taken into consideration. 
			It is important that this parameter remains consistent with the [max_seq] 
			parameter used in multi-alignment [muscle]. 
  -o OUT, --out=OUT     Output file: candidate primers. e.g.
                        [*].candidate.primers.txt.
  ```
  To extract primers of your amplicons from Oxford Nanopore Technology (ONT) reads. FindONTprimerV2.py = FindONTprimerV3.py. If your primers are degenerate, you can use the "FindONTexpandprimer.py" script included in the multiPrime package. This script is designed specifically to identify degenerate primer sequences from ONT reads and expand them into their full, non-degenerate forms.:
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
  To extract PCR products with perfect matches from your input FASTA file:
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
  To extract PCR products with mismatches from your input FASTA file
  ```bash
  python scripts/primer_coverage_validation_by_BWT.py
  ```
  ```
  Usage: primer_coverage_validation_by_BWT.py -i [input] -r [bowtie index] -l [150,2000] -p [10]-o [output]

  Options:
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input=INPUT_FILE
                        input file: primer.fa.
  -r REF, --ref=REF     reference file: template fasta or reference fasta.
  -l LEN, --len=LEN     Length of primer, which is used for mapping. If the length of the primer used for 
                        mapping is set to 0,the entire length of the primer will be utilized. Default: 0
  -t TERM, --term=TERM  Position of mismatch is not allowed in the 3 term of
                        primer. Default: 4
  -s SIZE, --s=SIZE     Length of PCR product, default: 150,2000.
  -p PROC, --proc=PROC  Number of process. Default: 20
  -b BOWTIE, --bowtie=BOWTIE
                        bowtie/ABS_path(bowtie) or bowtie2/ABS_path(bowtie2) was employed for mapping. 
                        Default: bowtie2
  -m SEEDMMS, --seedmms=SEEDMMS
                        Bowtie: Mismatches in seed (can be 0 - 3, default: -n
                        1).Bowtie2: Gap or mismatches in seed (can be 0 - 1,
                        default: -n 1).
  -d DICT, --dict=DICT
                        Dictionary of targets sequences, binary format. Default: None.
                        It can be obtained from prepare_fa_pickle.py (https://github.com/joybio/multiPrime).
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
		--Bowtie_DB: the Bowtie_DB is a pre-built index of a reference genome that is created from an input FASTA file. 
		--*.format.fa: reformat input fasta file. This process involves modifying the structure of the file to ensure compatibility with downstream analysis tools.
		--*.dict: binary dict of input fasta file. 
		--*.filtered.fa: filtered fasta. These fasta will not used in the downstream steps.
		--*.clstr: cluster information

	-Clusters_fa: genome file split by each cluster.
		--*.fa: fasta of each cluster, only acc_ID. For complete names, please refer to the files located in the "Clusters_target" directory.
		--*.tfa: top N {default: 500 randomly selected. Always contain the representative seq} fasta of each cluster
		--*.txt: accession id of each cluster

	-Clusters_msa: alginment by muscle
		--*.tmsa: muscle output of the top N {default: 500 randomly selected. Always contain the representative seq}

	-Clusters_trim_msa: trimmed alignment by degePrimer
		--*.trim.tmsa: trimmed muscle by degePrimer

	-Clusters_primer: get_degePrimer from degePrimer out
		--*.out: paired primers designed by the top N {default: 500} fasta
		--*.gap_seq_id_json: Positions and non-contained sequences caused by gap.
		--*.non_coverage_seq_id_json: Positions and non-contained sequences caused by others.

	-Clusters_cprimer: candidate primers for each cluster.
		--*.fa: candidate primers in fa format
		--*.txt: candidate primers in txt format (1 line)
		--*.xls: candidate primers in xls format

	-Clusters_target: information of each cluster.
		--*.txt: full name of each sequences in each cluster. You can get the full information of sequences from these files.

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
		--BWT_coverage: coverage of all primers in core primer set (up to 2-mismatch)
			--*.out: the start and stop positions of each accession and its corresponding primer pair.
			--*.pair.num: the target number of each primer.
			--*.total.acc.num: the total target number of final primer sets.
			--*.unmatched.fa: the sequences that were not captured by the core primer set with setting mismatches.
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
The project was conceptualized and scripted by Junbo Yang.Please send comments, suggestions, bug reports and bug fixes to 1806389316@pku.edu.cn / yang_junbo_hi@126.com.

# Todo
This repository is updated frequently. Expect breaking changes. More functions will be improved in the future and a new version is comming.


