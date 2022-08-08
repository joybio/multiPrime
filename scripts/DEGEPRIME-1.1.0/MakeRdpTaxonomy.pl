#!/usr/bin/perl -w

=head1 NAME

MakeRdpTaxonomy.pl - Extracts taxonomic information from an RDP (http://rdp.cme.msu.edu/) file in genbank format and outputs a tab-separated file with sequence id in first column and taxonomy in second column (with the different taxonoimc levels separated by semicolon).

=head1 USAGE

perl MakeRdpTaxonomy.pl -i RDP_GENBANK_FILE -o OUTPUT_FILE [-h]


=head1 POSITIONAL ARGUMENTS

-i RDP_GENBANK_FILE			Specify RDP file in genbank format

-o OUTPUT_FILE				Specify output file name (overwrites existing file with same name)


=head1 OPIONAL ARGUMENTS

-h							Print this help message

=cut

use Getopt::Long;

$infile = undef;
$outfile = undef;

&GetOptions('i=s' => \$infile, 'o=s' => \$outfile, 'h!' => \$help);
if (!$infile or !$outfile or $help) {
	system ('perldoc', $0);
	exit;
}

#####
print"\nRunning MakeRdpTaxonomy\n";
&extract_taxonomy_from_genbank;
print"Finnished MakeRdpTaxonomy succesfully\n\n";


#####

sub extract_taxonomy_from_genbank {
	$in = 0;
	open (INFILE, "$infile") || die ("Can't open $infile\n");
	open (OUT, ">$outfile");
	while (<INFILE>) {
		chomp;
		$row = $_;
		@fields = split(/\s+/, $row);
		if ($fields[0] eq "LOCUS") {
			$rdp_id = $fields[1];
		}
		if ($fields[0] eq "REFERENCE") {
			$in = 0;
		}
		if ($fields[0] eq "COMMENT") {
			$in = 0;
		}
		if ($in == 1) {	
			$row =~ tr/ //d;
			$tax = "$tax$row";
		}
		if (@fields > 1) {
			if ($fields[1] eq "ORGANISM") {
				$in = 1;
				$tax = "";
			}
		}
		if ($fields[0] eq "//") {
			$tax =~ m/(Root.*[^\.]\.)/;
			$tax = $1;
			print OUT "$rdp_id\t$tax\n";		
		}
	}
	close (INFILE);
	close(OUT);
}