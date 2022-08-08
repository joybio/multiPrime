#!/usr/bin/perl -w

=head1 NAME

MakeSilvaTaxonomy.pl - Extracts taxonomic information from a Silva (http://www.arb-silva.de/) file in fasta format and outputs a tab-separated file with sequence id in first column and taxonomy in second column (with the different taxonoimc levels separated by semicolon).

=head1 USAGE

perl MakeSilvaTaxonomy.pl -i SILVA_TAX_FASTA_FILE -o OUTPUT_FILE [-h]


=head1 POSITIONAL ARGUMENTS

-i RDP_GENBANK_FILE			Specify SILVA file in fasta format with taxonomies in sequence headers

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
print"\nRunning MakeSilvaTaxonomy\n";
&extract_taxonomy_from_silvafasta;
print"Finnished MakeSilvaTaxonomy succesfully\n\n";

#####

sub extract_taxonomy_from_silvafasta {
	open (INFILE, "$infile") || die ("Can't open $infile\n");
	open (OUT, ">$outfile");
	while (<INFILE>) {
		chomp;
		$row = $_;
		if (substr($row, 0, 1) eq ">") {
			@fields = split(/\s+/, $row, 2); # split at first space
			$id = $fields[0];
			substr($id, 0, 1) = "";
			$tax = $fields[1];
			print OUT "$id\t$tax\n"; 
		}	
	}
	close (INFILE);
	close(OUT);
}



