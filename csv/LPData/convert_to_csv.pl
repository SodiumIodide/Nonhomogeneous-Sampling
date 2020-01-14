#!/usr/bin/env perl

use warnings;
use strict;

sub main() {
    opendir my $dh, '.' or die "Cannot open local directory\n";
    my @files = grep { /.+\.ultra$/ and -f "./$_"} readdir $dh;
    foreach my $file (@files) {
        $file =~ /(.+)\.ultra/;
        my $filename = $1;
        open my $fh, $file or die "Cannot open $file\n";
        my $output_fh;
        print "Writing $filename files...";
        my $iterator = -1;
        while (my $line = <$fh>) {
            if ($line =~ /#/) {
                $iterator++;
                open $output_fh, '>', "$filename" . "_" . "$iterator.csv" or die "Cannot create $filename.csv\n";
                $line =~ s/#//;
                $line =~ s/psi/phi/;
            }
            $line =~ s/\s*(\S+)\s+(\S+)$/$1,$2/;
            print $output_fh "$line";
        }
        close $fh;
        close $output_fh;
        print "Done\n";
    }
}

main() unless caller;

__END__