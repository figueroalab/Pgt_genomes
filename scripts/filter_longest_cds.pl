#! /usr/bin/perl -w

# filter_longest_cds.pl

# This script will take a FASTA sequence file and retain only the longest sequence for each transcript

# Cory Hirsch
# March 25, 2019

use strict;
use Getopt::Long();
use Bio::SeqIO;

my $usage = "\nUsage: $0 --fasta_in <gene seq file> --rep_out <longest representative output file> --help <help on usage of script>\n\n";

my ($fasta_in, $rep_out, $help);

Getopt::Long::GetOptions('fasta_in=s' => \$fasta_in,
                         'rep_out=s' => \$rep_out,
                         'h|help' => \$help
                         );

if (defined($help)) {
    die $usage;
}

# Use Bioperl module to call in sequences
my $inseq = Bio::SeqIO->new(
                            -file   => "<$fasta_in",
                            -format => 'Fasta',
                           );

my %longest;
my $short_id;
# Call in each sequence in the file 
while (my $seq = $inseq->next_seq) {
    # Get the sequence identifier, the sequence, the length of the sequence, and the description
    my $id = $seq->display_id;
    my $sequence = $seq->seq;
    my $len = $seq->length();
    my $desc = $seq->desc();
    print "$id\t$desc\n";

    # get rid of transcript id from sequence identifier
    #if ($id =~ /^PGT/) {
    #$id =~ /(PGT\w+)\-\d+/;
    #$short_id = $1;
    #print "$short_id\n";
    #}
    #if ($id !~ /^PGT/) {
    #$short_id = $id;
    #print "$short_id\n";
    #}
    
    # Retain information for longest transcript in the hash
    # gene -> chr, start, stop, length, sequence
    if (exists $longest{$desc}) {
        my $hash_length = $longest{$desc}{'LENGTH'};
        if ($hash_length < $len) {
            $longest{$desc}{'LENGTH'} = $len;
            $longest{$desc}{'SEQ'} = $sequence;
            $longest{$desc}{'GENE'} = $id;
        }
    }
    if (! (exists $longest{$desc}) ) {
        $longest{$desc}{'LENGTH'} = $len;
        $longest{$desc}{'SEQ'} = $sequence;
        $longest{$desc}{'GENE'} = $id;
    }
}

# Open Bioperl for longetst transcript output file
my $seqout = Bio::SeqIO->new(
                             -file   => ">$rep_out",
                             -format => 'Fasta',
                            );

foreach my $key (keys %longest) {
    my $seq = Bio::Seq->new(
                            -seq => $longest{$key}{'SEQ'},
                            -id  => $longest{$key}{'GENE'},
                           );
    
    $seqout->write_seq($seq);
}

exit;