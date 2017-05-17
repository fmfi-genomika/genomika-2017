#!/usr/bin/perl -w
use strict;

# (1) quit unless we have the correct number of command-line args
my $num_args = $#ARGV + 1;
if ($num_args != 2) {
    print "\nUsage: name.pl file.psl\n";
    exit;
}

# (2) we got one command line args, so assume it is the
# infile name and outfile
my $infile = $ARGV[0];
my $outfile = $ARGV[1];

open my $inf, $infile or die "Could not open $infile: $!";
open(my $outf, '>', $outfile) or die "Could not open file '$outfile' $!";
while(my $line = <$inf>) {
    if ($line =~ /^#.*/) {
        print $outf $line;
    } else {
        my @columns = split "\t", $line;
        if ($columns[1] ne "0" || $columns[9] ne $columns[13] || $columns[11] ne $columns[15] || $columns[12] ne $columns[16]) {
            print $outf join("\t", @columns);
        }
    }
}
close $inf;
close $outf;
