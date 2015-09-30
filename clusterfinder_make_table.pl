#!/usr/bin/perl

# Program that makes the input table expected by ClusterFinder.
# Takes the table output of an hmmscan run, using the Pfam-A database, and
# either a Prokka annotation table or a prodigal output table 

# For prokka - supply the .tbl file
# For prodigal, make sure to generate a protein fasta, and supply the fasta file

use strict;
use warnings;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);

my $gene_table;
my $hmmscan_table;
my $input_type = 'prokka'; # Other possible value 'prodigal'
my $sequencing_status = 'Draft'; # Other possible value 'Finished'
my $organism_name;
my $organism_id;
my $output_path;

GetOptions(	"gene_positions|g=s"	=>	\$gene_table,
			"hmmscan_table|h=s"		=>	\$hmmscan_table,
			"type|t=s"				=>	\$input_type,
			"status|s=s"			=>	\$sequencing_status,
			"organism_name|n=s"		=>	\$organism_name,
			"organism_id|i=s"		=>	\$organism_id,
			"output|o=s"			=>	\$output_path )
or die ("Error in command line arguments\n");

# First parse gene table
my $geneHashRef = getGeneInfo($gene_table);

# Now parse hmmscan table
open (HMMTABLE, "<$hmmscan_table") or die "Couldn't open $hmmscan_table\n";
open (OUTPUT, ">$output_path") or die "Couldn't open $output_path\n";
while (<HMMTABLE>)
{
	my $line = $_;
	next if substr($line, 0, 1) eq '#';

	chomp $line;
	my @lineArray = split (/\s+/, $line);

	my $locus_tag = $lineArray[3];
	my $pfam_accession = $lineArray[1];
	my @pfamArray = split (/\./, $pfam_accession);
	my $pfam_id = $pfamArray[0];

	my $pfam_template_start = $lineArray[15];
	my $pfam_template_end = $lineArray[16];
	my $pfam_start = $lineArray[19];
	my $pfam_end = $lineArray[20];
	my $pfam_escore = $lineArray[6];
	my $enzyme_id = 'n/a';

	my $output_line = 	$geneHashRef->{ $locus_tag }->{ 'id' } . "\t" .
						$sequencing_status . "\t" .
						$organism_name . "\t" .
						$geneHashRef->{ $locus_tag }->{ 'contig' } . "\t" .
						$organism_id . "\t" .
						$locus_tag . "\t" .
						$geneHashRef->{ $locus_tag }->{ 'start' } . "\t" .
						$geneHashRef->{ $locus_tag }->{ 'end' } . "\t" .
						$geneHashRef->{ $locus_tag }->{ 'strand' } . "\t" .
						$pfam_template_start . "\t" .
						$pfam_template_end . "\t" . 
						$pfam_start . "\t" .
						$pfam_end . "\t" . 
						$pfam_id . "\t" .
						$pfam_escore . "\t" . 
						$enzyme_id . "\n";

	print OUTPUT $output_line;
}

sub getGeneInfo
{
	my $input_table = $_[0];
	my $tempHashRef;
	if ($input_type eq 'prokka')
	{
		$tempHashRef = parseProkkaTable( $input_table );
	} elsif ( $input_type eq 'prodigal' ) {
		$tempHashRef = parseProdigalTable( $input_table );
	} else {
		print "Error, invalid input type $input_type\n";
		exit 1;
	}

	return $tempHashRef;
}

sub parseProkkaTable
{
	my $input_table = $_[0];
	my %tempHash; # Hash will be keyed by locus tag
	my $current_contig;
	my $gene_id = 0;
	my $gene_start;
	my $gene_end;
	my $gene_strand;

	open (TABLE, "<$input_table") or die "\nCouldn't open $input_table\n";
	while (<TABLE>)
	{
		my $line = $_;
		chomp $line;
		if ( substr($line, 0, 1) eq '>' )
		{
			my @linearray = split (' ', $line);
			$current_contig = $linearray[1];
		} else {
			my @linearray = split ("\t", $line);
			if (looks_like_number($linearray[0]) && looks_like_number($linearray[1]) && ($linearray[2] eq 'CDS'))
			{
				if ( $linearray[1] > $linearray[0] )
				{
					$gene_id++;
					$gene_start = $linearray[0];
					$gene_end = $linearray[1];
					$gene_strand = '+';
				} else {
					$gene_id++;
					$gene_start = $linearray[1];
					$gene_end = $linearray[0];
					$gene_strand = '-';
				}
			} elsif ( (exists $linearray[3]) && ($linearray[3] eq 'locus_tag') ) {
				# First check everything is defined
				unless (defined( $current_contig ) && defined( $gene_start ))
				{
					print "Some variables not initialized\n";
					exit 1;
				}

				my $locus_tag = $linearray[4];
				$tempHash{ $locus_tag } = {
					'contig'	=>	$current_contig,
					'start'		=>	$gene_start,
					'end'		=>	$gene_end,
					'strand'	=>	$gene_strand,
					'id'		=>	$gene_id,
				};
			}
		}
	}
	close TABLE;

	return ( \%tempHash );
}

sub parseProdigalTable
{
	my $input_table = [0];

	my %tempHash; # Hash will be keyed by locus tag
	my $current_contig;
	my $gene_id;
	my $gene_start;
	my $gene_end;
	my $gene_strand;

	open (TABLE, "<$input_table") or die "\nCouldn't open $input_table\n";
	while (<TABLE>)
	{
		my $line = $_;
		chomp $line;															

		next unless ( substr( $line, 0, 1 ) eq '>' );

		my @linearray = split (' # ', $line);
		$gene_start = $linearray[1];
		$gene_end = $linearray[2];

		if ( $linearray[3] == 1 )
		{
			$gene_strand = '+';
		} elsif ( $linearray[3] == -1 ) {
			$gene_strand = '-';
		}

		# Get locus tag and contig
		my $locus_tag = substr( $linearray[0], 1 );
		my @locusArray = split ( '_', $locus_tag );
		$current_contig = $locusArray[0];

		# Get gene id
		my @extraInfoArray = split ( ';', $linearray[4] );
		my @idArray = split ( '=', $extraInfoArray[0] );
		$gene_id = $idArray[1];

		$tempHash{ $locus_tag } = {
					'contig'	=>	$current_contig,
					'start'		=>	$gene_start,
					'end'		=>	$gene_end,
					'strand'	=>	$gene_strand,
					'id'		=>	$gene_id,
				};
	}
	close TABLE;

	return ( \%tempHash );
}