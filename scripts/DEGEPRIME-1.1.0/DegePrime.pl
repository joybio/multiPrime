#!/usr/bin/perl -w


=head1 NAME

DegePrime.pl - Finds degenerate oligomers of as high coverage as possible at every position of a multiple sequence alignment. 


=head1 USAGE

perl degeprime.pl -i <ALIGNMENT_FILE> -l <OLIGOMER_LENGTH> -d <MAX_DEGENERACY> -o <OUTPUT_FILE> [-depth <MIN_DEPTH>] [-skip <SKIP_LENGTH>] [-iter <NUMBER_ITERATIONS>] [-taxfile <TAXONOMY_FILE>] [-taxlevel <TAXONOMY_LEVEL>] [-h]


=head1 POSITIONAL ARGUMENTS

-i <ALIGNMENT_FILE>			Specify alignment file. This has some special requirements, see README

-l <OLIGOMER_LENGTH>		Specify oligomer length (integer)

-d <MAX_DEGENERACY>			Specify maximum degeneracy of oligomer (integer)

-o <OUTPUT_FILE>			Specify output file name (overwrites existing file with same name)


=head1 OPIONAL ARGUMENTS

-depth <MIN_DEPTH>			Specify minimum number of spanning sequences without in/dels for a window position to be included in analysis, default 1

-skip <SKIP_LENGTH>			Specify number of bases in start and end of sequences that will not be considered in analysis, default 20

-iter <NUMBER_ITERATIONS>	Specify number of iterations during Weighted Randomised Merging, default 100

-taxfile <TAXONOMY_FILE>	Specify file with taxonomic annotations for the sequences in the alignment file. This is required if <TAXONOMY_LEVEL> is given

-taxlevel <TAXONOMY_LEVEL>	Specify what taxonomic level to use (a positive integer). This is required if <TAXONOMY_FILE> is given

-h							Prints this help message

[Press q to close this help message]

=cut

#use strict;
use Storable qw<dclone>;
use Getopt::Long;
use Time::HiRes qw(time);

### default parameter setting ###

$min_depth = 1; # minimum number of spanning sequences without in/dels to include in a window.
$skip_length = 20; # number of bases in both ends of sequence that will not be considered in analysis.
$number_iteratons = 100; # number of iterations in Weighted Randomised Merging
$aligned_short_file = undef;
$window_length = undef;
$max_deg = undef;
$taxonomy_file = undef;
$taxonomy_level = undef;

######################

&GetOptions('i=s' => \$aligned_short_file, 'l=i' => \$window_length, 'd=i' => \$max_deg, 'o=s' => \$outfile, 'depth=i' => \$min_depth, 'skip=i' => \$skip_length, 'iter=i' => \$number_iteratons, 'taxfile=s' => \$taxonomy_file, 'taxlevel=i' => \$taxonomy_level, 'h!' => \$help);
if (!$aligned_short_file or !$window_length or !$max_deg or !$outfile or $help) {
	system ('perldoc', $0);
	exit;
}
if ($taxonomy_file) {
	if (!$taxonomy_level) {
		print "\n  Error: No TAXONOMY_LEVEL given. This (positive integer) is needed when a TAXONOMY_FILE is given. \n";
		print "\n  perl degeprime.pl -h for more help \n\n";
		exit;
	}
}
if ($taxonomy_level) {
	if (!$taxonomy_file) {
		print "\n  Error: No TAXONOMY_FILE given. This is needed when a TAXONOMY_LEVEL is given. \n";
		print "\n  perl degeprime.pl -h for more help \n\n";
		exit;
	}
	if ($taxonomy_level < 1) { 
		print "\n  Error: TAXONOMY_LEVEL should be a positive integer. \n"; 
		print "\n  perl degeprime.pl -h for more help \n\n";
		exit;
	}
}

#$taxonomy_file = "rdp_10_7_taxonomy.txt";
#$aligned_short_file = "arch_align_rdp_10_7_short.fasta";
#$aligned_short_file = "rdp_download_138807seqs_v9_1200nt_good_aligned_short_r2000.fna";
#$aligned_short_file = "polyG_mouse_all_slim.fa";
#$aligned_short_file = "rdp_download_138807seqs_v9_1200nt_good_aligned_short.fna";
#$aligned_short_file = "eukaryotic/arb-silva.de_2011-06-01_id2791_0.90.fasta";

#############

$totaltime = time;

print"\nRunning DegePrime\n";
$changed_max_deg = undef;
while (&check_max_deg($max_deg) == 0) {
	$max_deg = $max_deg - 1;
	$changed_max_deg = 1;
}
if ($changed_max_deg) {
	print"Max degeneracy was not a valid degeneracy and has been changed to $max_deg\n";
}
&make_iupac;
if ($taxonomy_level) {
	&get_taxonomy;
}
&read_short_seqs;
&check_lengths;
&calc_coverage;
print"Finnished DegePrime run succesfully\n\n";

$totaltime = time - $totaltime;
print"Run time:\t$totaltime\n";

#############

sub check_max_deg {
	local($ok) = 0;
	local($a) = $_[0];
	local($b);
	local($c);
	local($n) = 0;
	while ( int($a/2**($n+1)) == $a/2**($n+1) ) {
		$n++;
	}
	$b = $a/2**$n;
	$n = 0;
	while (int($b/3**($n+1)) == $b/3**($n+1)) {
		$n++;
	}
	$c = $b/3**$n;
	$ok = 1 if ($c == 1);
	return($ok);
}

sub calc_coverage {
	print"Finding primers\n";
	open (OUT, ">$outfile");
	print OUT "Pos\tNumberSpanning\tUniqueMers\tEntropy\tPrimerDeg\tPrimerSeq\tNumberMatching\tFractionMatching";
	if ($taxonomy_level) {
		foreach $taxon (@taxa) {
			print OUT "\tSpanning $taxon\tMatching $taxon\tFraction $taxon";
		}
	}
    print OUT "\n";
	for ($pos = 0; $pos <= ($seq_length - $window_length); $pos++) {
	#for ($pos = 385; $pos <= 385; $pos++) {
	#for ($pos = 470; $pos <= 519; $pos++) {
		print"	finding primer in position $pos\n";
		$time0 = time;
		($total_spanning, $zero_insertions, $entropy) = &get_mer_ranking($pos);
		if ($zero_insertions >= $min_depth) {
			$unique = @sorted_mer;
			print OUT "$pos\t$total_spanning\t$unique\t$entropy";
			&weighted_randomised_combination($pos);
			#&hyden($pos);
			$time1 = time;
			$time = ($time1 - $time0);
			#print OUT "\t$time";
			print OUT "\n";	
		}
	}
	#print"\n\n  Elapsed Time: $time sec\n\n";
}

sub read_short_seqs {
	print"Reading sequences\n";
	open (INFILE, $aligned_short_file) || die ("Can't open");
	$seq = "";
	while (<INFILE>) {
		chomp;
		$row = $_;
		if (substr($row, 0, 1) eq ">") {
			if ($seq ne "") {
				$id_seq{$id} = $seq;
				$seq = "";
			}
			@fields = split(/\s+/, $row); # split id at space
			$id = $fields[0];
			substr($id, 0, 1) = "";
		} else {
			$seq = $seq.$row;
		}	
	}
	$id_seq{$id} = $seq;
	close (INFILE);
}

sub check_lengths {
	print"Checking sequence lengths\n";
	foreach $id (keys %id_seq) {
		$seq = $id_seq{$id};
		$seq_length = length($seq);
		$length_id{$seq_length} = $id;
		die if ($seq !~ m/^-*/); # ??
		$seq =~ m/(^\W*)/;
		#$seq =~ m/(^-*)/;
		$start{$id} = length($1); # position for first base
		die if ($seq !~ m/-*$/); # ??
		$seq =~ m/(\W*$)/;
		#$seq =~ m/(-*$)/;
		$end{$id} = length($seq) - length($1) - 1; # position for last base
	}
	if ((keys %length_id) > 1) {
		die ("Error: Not all aligned sequences have same length. \n");
	}
}

sub get_mer_ranking {
	local($startpos) = $_[0];
	local($mer);
	local($total_spanning) = 0;
	local($zero_gaps) = 0;
	local($entropy) = 0;
	local(%mer_count_entropy) = ();
	%id_spanning = ();
	%mer_n_id = (); 
	%mer_count = ();
	@sorted_mer = ();
	foreach $id (keys %id_seq) {
		if (($start{$id} + $skip_length) <= $startpos) {
			if (($end{$id} - $skip_length) >= $startpos + $window_length - 1) {
				$seq = $id_seq{$id};
				$mer = substr($seq, $startpos, $window_length);
				substr($mer, -1, 1) =~ tr/[a-z]/[A-Z]/;
				if ($mer !~ m/[^(A|T|C|G)]/) {
					$mer_count{$mer}++;
					$zero_gaps++;
					$mer_n_id{$mer}{$id} = 1;
				}
				$total_spanning++;
				$id_spanning{$id} = 1;
				# for entropy calc:
				$mer_count_entropy{$mer}++;	
			}
		}
	}
	#@sorted_mer = (sort {$mer_count{$b} <=> $mer_count{$a}} keys %mer_count); # original
	#@sorted_mer = (sort {$mer_count{$a} <=> $mer_count{$b}} keys %mer_count);
	@sorted_mer = keys %mer_count;
	# Calculate entropy in window
	if ($total_spanning > 0) {
		foreach $mer (keys %mer_count_entropy) {
			$entropy = $entropy - ($mer_count_entropy{$mer}/$total_spanning)*log($mer_count_entropy{$mer}/$total_spanning)/log(2);
		}
	}
	return($total_spanning, $zero_gaps, $entropy);
}

sub weighted_randomised_combination {
	local($pos);
	local($deg);
	local($matching);
	local($min_index);
	local($max_index);
	local(%pos_n_base);
	local(%temp_pos_n_base);
	local(%pos_n_primers);
	local(%already_chosen);
	local(@end_indices) = ();
	local(@temp_end_indices);
	local(@sorted_indices);
	local($bestmatch) = 0;
	local($bestdeg) = "";
	local($adjust_counts);
	local($adjust_index);
    local($fraction_match);
	# Make list from which to make random draws:
	$end_indices[0] = $mer_count{ $sorted_mer[0] };
	for ($i = 1; $i < @sorted_mer; $i++) {
		$end_indices[$i] = $end_indices[$i - 1] + $mer_count{ $sorted_mer[$i] };
	}
	# Iterate randomised summation procedure:
	$timeA = $timeB = $timeC = 0;
	for ($iter = 0; $iter < $number_iteratons; $iter++) {
		@temp_end_indices = @end_indices;
		@temp_sorted_mer = @sorted_mer;
		$deg = 0;
		$trial = 0;
		$matching = 0;
		%pos_n_base = ();
		%already_chosen = ();
		while ($deg < $max_deg) {
			$trial++;
			last if ($trial > 100);
			# Get random mer:
			$random = int(rand($temp_end_indices[-1]));
			$min_index = 0;
			$max_index = @temp_end_indices;
			$timeA = time;
			while ($max_index > $min_index) {
				$i = int( ($min_index + $max_index)/2 );
				if ($temp_end_indices[$i] < $random) {
					$min_index = $i;
				} elsif ($temp_end_indices[$i - 1] >= $random) {
					$max_index = $i;
				} else {
					last;					
				}
			}
			$timeB = $timeB + time - $timeA;
			$selected_index = $i;	
			$mer = $temp_sorted_mer[$selected_index];
			if (defined $already_chosen{$selected_index}) {
				$trial--;
				$fails++;
			} else {
				$already_chosen{$selected_index} = $selected_index;
				$fails = 0;
				# Calculate newdeg:
				$newdeg = 1;
				%temp_pos_n_base = %{dclone(\%pos_n_base)};
				for ($pos = 0; $pos < $window_length; $pos++) {
					$nt = substr($mer, $pos, 1);
					$temp_pos_n_base{$pos}{$nt} = 1;
					$newdeg = $newdeg * (keys %{$temp_pos_n_base{$pos}});
				}
				# Add mer if newdeg <= max_deg:
				if ($newdeg <= $max_deg) {
					$deg = $newdeg;
					$matching = $matching + $mer_count{$mer};
					%pos_n_base = %{dclone(\%temp_pos_n_base)};
				}
			}
			# Update @temp_end_indices & @temp_sorted_mer
			if ($fails > 4) { # 4		
				$timeA = time;
				@sorted_indices = (sort {$already_chosen{$a} <=> $already_chosen{$b}} keys %already_chosen);
				$adjust_counts = 0;
				$adjust_index = 0;				
				for ($i = $sorted_indices[0]; $i < (@temp_end_indices - @sorted_indices - 1); $i++) {	
					while (defined $already_chosen{$i + $adjust_index}) {
						$adjust_counts = $adjust_counts + $mer_count{ $temp_sorted_mer[$i + $adjust_index] };
						$adjust_index++;
					}		
					$temp_end_indices[$i] = $temp_end_indices[$i + $adjust_index] - $adjust_counts;
				}
				for ($i = 0; $i < @sorted_indices; $i++) {
					splice(@temp_end_indices, -1, 1);			
				}
				for ($i = (@sorted_indices - 1); $i >= 0; $i--) {
					splice(@temp_sorted_mer, $sorted_indices[$i], 1);
				}
				%already_chosen = ();
				$timeC = $timeC + time - $timeA;
 			}
			# Stop if all mers have been tried:
			last if (@temp_end_indices == 0);
		}
		
		# Count matches:
		%pos_n_primers = ();
		$pos_n_primers{-1}{""} = 1;
		for ($pos = 0; $pos < $window_length; $pos++) {
			foreach $primer (keys %{$pos_n_primers{$pos - 1}}) {
				foreach $nt (keys %{$pos_n_base{$pos}}) {		
					$new_primer = $primer.$nt;
					$pos_n_primers{$pos}{$new_primer} = 1;
				}
			}
		}
		$match = 0;
		foreach $primer (keys %{$pos_n_primers{$window_length - 1}}) {
			if (defined $mer_count{$primer}) { 
				$match = $match + $mer_count{$primer};
			}				
		}
		if ($match > $bestmatch) {
			$bestmatch = $match;
			$bestdeg = $deg;
			$deg_primer = "";
			for ($pos = 0; $pos < $window_length; $pos++) {
				$string = join("", keys %{$pos_n_base{$pos}});
				$deg_primer = $deg_primer.$deg{$string};
			}
			$bestprimer = $deg_primer;
			@best_primers = (keys %{$pos_n_primers{$window_length - 1}});
		}
	}
    if ($total_spanning > 0) {
        $fraction_match = $bestmatch/$total_spanning;
    } else {
        $fraction_match = "NA";
    }
	if ($taxonomy_level) {
		%taxon_counts_spanning = ();
		%taxon_counts_matching = ();
		foreach $id (keys %id_spanning) {
			$taxon_counts_spanning{$id_taxon{$id}}++;
		}
		foreach $primer (@best_primers) {
			foreach $id (keys %{$mer_n_id{$primer}}) {
				$taxon_counts_matching{$id_taxon{$id}}++;
			}
		}
        print OUT "\t$bestdeg\t$bestprimer\t$bestmatch\t$fraction_match";
		foreach $taxon (@taxa) {
			if (defined $taxon_counts_spanning{$taxon}) {
				print OUT "\t$taxon_counts_spanning{$taxon}";
			} else {
				print OUT "\t0";
			}
			if (defined $taxon_counts_matching{$taxon}) {
                $fraction_match = $taxon_counts_matching{$taxon}/$taxon_counts_spanning{$taxon};
				print OUT "\t$taxon_counts_matching{$taxon}\t$fraction_match";
			} else {
				print OUT "\t0\tNA";
			}
		}
	} else {
        print OUT "\t$bestdeg\t$bestprimer\t$bestmatch\t$fraction_match";
	}
}

sub make_iupac {
	$deg{"A"} = "A";
	$deg{"C"} = "C";
	$deg{"G"} = "G";
	$deg{"T"} = "T";

	$deg{"AT"} = "W";
	$deg{"TA"} = "W";

	$deg{"CG"} = "S";
	$deg{"GC"} = "S";

	$deg{"AC"} = "M";
	$deg{"CA"} = "M";

	$deg{"GT"} = "K";
	$deg{"TG"} = "K";

	$deg{"AG"} = "R";
	$deg{"GA"} = "R";

	$deg{"CT"} = "Y";
	$deg{"TC"} = "Y";

	$deg{"CGT"} = "B";
	$deg{"CTG"} = "B";
	$deg{"TCG"} = "B";
	$deg{"TGC"} = "B";
	$deg{"GTC"} = "B";
	$deg{"GCT"} = "B";

	$deg{"AGT"} = "D";
	$deg{"ATG"} = "D";
	$deg{"TAG"} = "D";
	$deg{"TGA"} = "D";
	$deg{"GAT"} = "D";
	$deg{"GTA"} = "D";

	$deg{"ACT"} = "H";
	$deg{"ATC"} = "H";
	$deg{"TAC"} = "H";
	$deg{"TCA"} = "H";
	$deg{"CTA"} = "H";
	$deg{"CAT"} = "H";

	$deg{"ACG"} = "V";
	$deg{"AGC"} = "V";
	$deg{"GAC"} = "V";
	$deg{"GCA"} = "V";
	$deg{"CAG"} = "V";
	$deg{"CGA"} = "V";

	$deg{"ACGT"} = "N";
	$deg{"ACTG"} = "N";
	$deg{"ATCG"} = "N";
	$deg{"ATGC"} = "N";
	$deg{"AGTC"} = "N";
	$deg{"AGCT"} = "N";

	$deg{"CAGT"} = "N";
	$deg{"CATG"} = "N";
	$deg{"CTAG"} = "N";
	$deg{"CTGA"} = "N";
	$deg{"CGAT"} = "N";
	$deg{"CGTA"} = "N";

	$deg{"GACT"} = "N";
	$deg{"GATC"} = "N";
	$deg{"GTAC"} = "N";
	$deg{"GTCA"} = "N";
	$deg{"GCTA"} = "N";
	$deg{"GCAT"} = "N";

	$deg{"TACG"} = "N";
	$deg{"TAGC"} = "N";
	$deg{"TGAC"} = "N";
	$deg{"TGCA"} = "N";
	$deg{"TCAG"} = "N";
	$deg{"TCGA"} = "N";

	$deg_n_nt{"A"}{"A"} = 1;
	$deg_n_nt{"C"}{"C"} = 1;
	$deg_n_nt{"G"}{"G"} = 1;
	$deg_n_nt{"T"}{"T"} = 1;

	$deg_n_nt{"W"}{"A"} = 1;
	$deg_n_nt{"W"}{"T"} = 1;

	$deg_n_nt{"S"}{"C"} = 1;
	$deg_n_nt{"S"}{"G"} = 1;

	$deg_n_nt{"M"}{"A"} = 1;
	$deg_n_nt{"M"}{"C"} = 1;

	$deg_n_nt{"K"}{"G"} = 1;
	$deg_n_nt{"K"}{"T"} = 1;

	$deg_n_nt{"R"}{"A"} = 1;
	$deg_n_nt{"R"}{"G"} = 1;

	$deg_n_nt{"Y"}{"C"} = 1;
	$deg_n_nt{"Y"}{"T"} = 1;

	$deg_n_nt{"B"}{"C"} = 1;
	$deg_n_nt{"B"}{"G"} = 1;
	$deg_n_nt{"B"}{"T"} = 1;

	$deg_n_nt{"D"}{"A"} = 1;
	$deg_n_nt{"D"}{"G"} = 1;
	$deg_n_nt{"D"}{"T"} = 1;

	$deg_n_nt{"H"}{"A"} = 1;
	$deg_n_nt{"H"}{"C"} = 1;
	$deg_n_nt{"H"}{"T"} = 1;

	$deg_n_nt{"V"}{"A"} = 1;
	$deg_n_nt{"V"}{"C"} = 1;
	$deg_n_nt{"V"}{"G"} = 1;

	$deg_n_nt{"N"}{"A"} = 1;
	$deg_n_nt{"N"}{"C"} = 1;
	$deg_n_nt{"N"}{"G"} = 1;
	$deg_n_nt{"N"}{"T"} = 1;
}

sub get_taxonomy {
	print"Reading taxonomy\n";
	open (INFILE, $taxonomy_file) || die ("Can't open $taxonomy_file");
	while (<INFILE>) {
		chomp;
		@fields = split(/\t/);
		@subfields = split(/;/, $fields[1]);
		if (@subfields < $taxonomy_level) {
			$taxon = $fields[1];
		} else {
			$taxon = join(';', @subfields[0..($taxonomy_level-1)]);
		}
		$id_taxon{$fields[0]} = $taxon;
		$taxon{$taxon} = 1;
	}
	close (INFILE);
	@taxa = (keys %taxon);
	@taxa = sort @taxa;
}

##################################################
##################################################
# The remainder is not used by Degeprime

sub hyden {
	$hyden_outfile = "hyden_out.txt";
	$hyden_infile = "hyden_in.fasta";
	local($startpos) = $_[0];
	local($deg) = undef;
	local($match) = 0;
	local($deg_primer) = "";
	local($new_primer) = "";
	local(%pos_n_primers) = ();
	local($pos) = 0;
	local($total_spanning) = 0;
	local($zero_gaps) = 0;
	open (OUT2, ">$hyden_infile") || die "Can't make $hyden_infile\n\n";
	foreach $id (keys %id_seq) {
		if (($start{$id} + $skip_length) <= $startpos) {
			if (($end{$id} - $skip_length) >= $startpos + $window_length - 1) {
				$seq = $id_seq{$id};
				$mer = substr($seq, $startpos, $window_length);
				substr($mer, -1, 1) =~ tr/[a-z]/[A-Z]/;
				if ($mer !~ m/[^(A|T|C|G)]/) {
					$mer = $mer."A";
					print OUT2 ">$id\n$mer\n";
					$zero_gaps++;
				}
				$total_spanning++;
			}
		}
	}
	$command = "hyden.exe -dna $hyden_infile -len5 $window_length -len3 1 -deg5 $max_deg -deg3 1 -from5 0 -to5 0 -from3 0 -to3 0 -mis5 0 -mis3 0 -nentropy $zero_gaps -mprimers 1 > $hyden_outfile";
	system($command);
	($deg, $match, $deg_primer) = &get_results;
	print OUT "\t$deg\t$match\t$deg_primer";
}

sub get_results {
	local($stop) = 0;
	local($primer) = undef;
	local($match) = undef;
	local($total) = undef;
	local($deg) = undef;
	open (INFILE, $hyden_outfile) || die ("Can't open $hyden_outfile");
	while (<INFILE>) {
		chomp;
		if ($stop == 1) {
			@fields = split(/\s+/);
			$primer = $fields[1];
			last;
		}
		if (substr($_, 0, 6) eq "Best 5") {
			@fields = split(/\s+/);
			$match = $fields[4];
			$total = $fields[8];
		}
		if (substr($_, 0, 7) eq "Primer:") {
			@fields = split(/\s+/);
			$deg = $fields[4];
			$stop = 1;
		}		
	}
	return($deg, $match, $primer);
	#print OUT "$i\t$deg\t$match\t$primer\n";
}

#####################################################

sub get_base_ranking {
	local($startpos) = $_[0];
	local($pos);
	local($total_spanning) = 0;
	local($zero_gaps) = 0;
	@nt_matr = ();
	@freq_matr = ();
	%pos_n_base = ();
	foreach $id (keys %id_seq) {
		if (($start{$id} + $skip_length) <= $startpos) {
			if (($end{$id} - $skip_length) >= $startpos + $window_length - 1) {
				$seq = $id_seq{$id};
				$mer = substr($seq, $startpos, $window_length);
				substr($mer, -1, 1) =~ tr/[a-z]/[A-Z]/; 
				if ($mer !~ m/[^(A|T|C|G)]/) {
					for ($pos = 0; $pos < $window_length; $pos++) {
						$nt = substr($mer, $pos, 1);
						$pos_n_base{$pos}{$nt}++;
					}
					$zero_gaps++;
				}
				$total_spanning++;
			}
		}
	}
	if ($zero_gaps > 0) {
		for ($pos = 0; $pos < $window_length; $pos++) {
			@sorted_nt = (sort {$pos_n_base{$pos}{$b} <=> $pos_n_base{$pos}{$a}} keys %{$pos_n_base{$pos}});
			for ($i = 0; $i < 4; $i++) {
				if ($pos_n_base{$pos}{$sorted_nt[$i]} > 0) {
					$nt_matr[$pos][$i] = $sorted_nt[$i];
					$freq_matr[$pos][$i] = $pos_n_base{$pos}{$sorted_nt[$i]} / $total;
				} else {
					$freq_matr[$pos][$i] = 0;
					$nt_matr[$pos][$i] = "X";
				}		
			}
		}
	}
	return($total_spanning, $zero_gaps);
}

sub expansion {
	local($pos);
	local($deg) = 1;
	local($matching);
	@rank = ();
	# set rank to most frequent base:
	for ($pos = 0; $pos < $window_length; $pos++) {
		$rank[$pos] = 0;
	}
	while ($deg < $max_deg) {
		$highest_freq = 0;
		for ($pos = 0; $pos < $window_length; $pos++) {
			if ($freq_matr[$pos][$rank[$pos] + 1] > $highest_freq) {
				$newdeg = $deg * ($rank[$pos] + 2) / ($rank[$pos] + 1);
				if ($newdeg <= $max_deg) {
					$highest_freq = $freq_matr[$pos][$rank[$pos] + 1];
					$highest_pos = $pos;
				}
			}
		}
		last if ($highest_freq == 0);
		$newdeg = $deg * ($rank[$highest_pos] + 2) / ($rank[$highest_pos] + 1);
		$deg = $newdeg;
		$rank[$highest_pos]++;
	}
	$matching = &check_matching;
	print"\t$deg\t$matching";
	#print OUT "\t$matching";
}

sub constriction {
	local($pos);
	local($deg) = 1;
	local($matching);
	@rank = ();
	# set rank to base with lowest freq (above zero)
	for ($pos = 0; $pos < $window_length; $pos++) {
		for ($rank = 0; $rank < 4; $rank++) {
			if ($freq_matr[$pos][$rank] > 0) {
				$rank[$pos] = $rank;
			}
		}
	}
	for ($pos = 0; $pos < $window_length; $pos++) {
		$deg = $deg * ($rank[$pos] + 1);
	}
	while ($deg > $max_deg) {
		$lowest_freq = 10;
		for ($pos = 0; $pos < $window_length; $pos++) {
			if ($rank[$pos] > 0) {
				if ($freq_matr[$pos][$rank[$pos]] < $lowest_freq) {
					$lowest_freq = $freq_matr[$pos][$rank[$pos]];
					$lowest_pos = $pos;
				}
			}
		}
		last if ($lowest_freq == 10);
		$newdeg = $deg * ($rank[$lowest_pos]) / ($rank[$lowest_pos] + 1);
		$deg = $newdeg;
		$rank[$lowest_pos]--;
		last if ($newdeg < $max_deg);
	}
	$matching = &check_matching;
	print"\t$deg\t$matching";
	#print OUT "\t$matching";
}

sub summation {
	local($pos);
	local($deg);
	local($matching) = 0;
	local(%pos_n_base) = (); #there is a global one as well!
	local(%temp_pos_n_base) = (); #there is a global one as well!
	local(%pos_n_primers) = ();
	#sum best primers while deg < max_deg:
	foreach $mer (@sorted_mer) {
		$deg = 1;
		%temp_pos_n_base = %{dclone(\%pos_n_base)};
		for ($pos = 0; $pos < $window_length; $pos++) {
			$nt = substr($mer, $pos, 1);
			$temp_pos_n_base{$pos}{$nt} = 1;
			$deg = $deg * (keys %{$temp_pos_n_base{$pos}});
		}
		if ($deg <= $max_deg) {
			%pos_n_base = %{dclone(\%temp_pos_n_base)};
		} else {
			last;
		}
	}
	#count matches:
	$pos_n_primers{-1}{""} = 1;
	for ($pos = 0; $pos < $window_length; $pos++) {
		foreach $primer (keys %{$pos_n_primers{$pos - 1}}) {
			foreach $nt (keys %{$pos_n_base{$pos}}) {		
				$new_primer = $primer.$nt;
				$pos_n_primers{$pos}{$new_primer} = 1;
			}
			
		}
	}
	$match = 0;
	foreach $primer (keys %{$pos_n_primers{$window_length - 1}}) {
		if (defined $mer_count{$primer}) { 
			$match = $match + $mer_count{$primer};
		}				
	}
	#make deg primer:
	$deg_primer = "";
	$deg = 1;
	for ($pos = 0; $pos < $window_length; $pos++) {
		$string = join("", keys %{$pos_n_base{$pos}});
		$deg_primer = $deg_primer.$deg{$string};
		$deg = $deg * (keys %{$pos_n_base{$pos}});
	}
	print"\t$deg\t$match\t$deg_primer";
	#print OUT "\t$deg\t$match";
}

sub print_primer {
	for ($pos = 0; $pos < $window_length; $pos++) {
		print"$pos\t";
		for ($rank = 0; $rank <= $rank[$pos]; $rank++) {
			print" $nt_matr[$pos][$rank]";
		}
		print"\n";
	}
}

sub check_matching {
	local($match) = 0;
	local($primer) = "";
	local($new_primer) = "";
	local(%pos_n_primers) = ();
	local($pos) = 0;	
	$pos_n_primers{-1}{""} = 1;
	for ($pos = 0; $pos < $window_length; $pos++) {
		foreach $primer (keys %{$pos_n_primers{$pos - 1}}) {
			for ($rank = 0; $rank <= $rank[$pos]; $rank++) {
				$nt = $nt_matr[$pos][$rank];
				$new_primer = $primer.$nt;
				$pos_n_primers{$pos}{$new_primer} = 1;
			}
		}
	}
	foreach $primer (keys %{$pos_n_primers{$window_length - 1}}) {
		if (defined $mer_count{$primer}) {
			$match = $match + $mer_count{$primer};
		}				
	}
	return ($match);
}

