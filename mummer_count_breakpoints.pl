#!/usr/bin/env perl

=head1 NAME

mummer_count_breakpoints.pl - Count syntenic breakpoints in Mummer-aligned genome pair

=head1 SYNOPSIS

perl mummer_count_breakpoints.pl --mummer mummer_output.coords [--plot] [--verbose]

perl mummer_count_breakpoints.pl --help

=head1 DESCRIPTION

This script counts the number of synteny breakpoints between two genomes that
have been aligned with Mummer. Two aligned genomes will usually share conserved,
aligned regions. A group of conserved regions with the same arrangement and
orientation between both genomes constitute a syntenic block. This script counts
the number of breakpoints between syntenic blocks, which are evidence for
genomic rearrangement. It accounts for fragmented genome assemblies (i.e. where
the genome assembles into several contigs) by not counting breakpoints that
touch contig boundaries, because these are ambiguous with respect to genomic
rearrangement. 

=head1 PREPARING DATA WITH MUMMER

The input to this script is a coordinates file from a pairwise alignment of two
genomes produced by Mummer. You can use either the Nucmer or Promer aligners to
perform the alignment. Nucmer aligns nucleotide sequences, whereas Promer uses
translated sequences, which helps to align genomes that are not as closely
related.

This script has been tested with Mummer 3. Other versions of Mummer may not
function as expected.

The following describes alignment of two nucleotide Fasta files, genome1.fasta
and genome2.fasta, with Promer, and how to convert the output to a coordinates
file for the breakpoint calculation. PREFIX is the filename prefix used for the
output files by Mummer.

=over 4

=item promer --prefix PREFIX genome1.fasta genome2.fasta

Perform the alignment with Promer aligner.

=item delta-filter -q -r PREFIX.delta > PREFIX.filter

Filter to retain only bidirectional best hits with the -q and -r options

=item show-coords -B -k PREFIX.filter > PREFIX.coords

Convert to genome coordinates in the btab format (option -B) and remove
alignments that overlap other better-scoring alignments (option -k).

=back

=head1 ARGUMENTS

=over 8

=item --mummer I<FILE>

Coordinates of aligned regions between two genomes, as produced by show-coords
from Mummer, in tab-separated btab format (-B option of show-coords).

=item --plot

Optional: Produce dot-plot of the two aligned genomes. Requires the R script
plot_segments.R.

=item --verbose

Optional: Report the number of breakpoints counted in a full English sentence.

=back

=head1 OUTPUT

By default reports input filename, count of synteny breakpoints, and count of
aligned segments with breakpoints as tab-separated fields to STDOUT.

When --verbose flag is used, the counts are reported in a sentence to STDOUT.

When --plot flag is used, a PDF file of the genome alignment dotplot is
produced with the same filename prefix as the input file.

=head1 COPYRIGHT AND LICENSE

Copyright 2017, Brandon Seah (kbseah@mpi-bremen.de)

LICENSE
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

=cut 

use strict;
use warnings;
use diagnostics;

use Getopt::Long;
use File::Basename;
use Pod::Usage;
use FindBin;

my $mummer_file;
my $file_prefix;
my $path_to_r_script = $FindBin::RealBin."/plot_segments.R";
my $do_plot;
my $verbose; # verbose output

pod2usage(-message=>"Invalid input arguments", -exitstatus=>2) if (! @ARGV); 

GetOptions("mummer=s"=>\$mummer_file,
           "plot"=>\$do_plot,
           "verbose"=>\$verbose,
           "help|h" => sub { pod2usage(-exitstatus=>2, -verbose=>2); },
           "man|m" => sub { pod2usage(-exitstatus=>0, -verbose=>2); }
           ) or pod2usage(-exitstatus=>2, -verbose=>2);

# Parse input file

$file_prefix = fileparse($mummer_file, (".coords"));

###############################################################################
# IDENTIFY AND COUNT BREAKPOINTS

my %segment_pos_hash;
my %ref_len_hash;
my %query_len_hash;
my %neighbor_hash; # Hash of neighboring contig IDs

hash_mummer($mummer_file, \%segment_pos_hash, \%ref_len_hash, \%query_len_hash);

# Hash references
foreach my $contig (keys %ref_len_hash) {
    getneighbors_ref(\%segment_pos_hash, \%neighbor_hash, $contig, "ref");
}
# Hash queries
foreach my $contig (keys %query_len_hash) {
    getneighbors_ref(\%segment_pos_hash, \%neighbor_hash, $contig, "query");
}

# Identify segments with breakpoints
my @breakpoints = get_breakpoints(\%neighbor_hash);

# Count total breakpoints
my $total_breakpoints = 0;
foreach my $id (keys %neighbor_hash) {
    $total_breakpoints += $neighbor_hash{$id}{"flag"};
}

# Print output

if ($verbose) {
    print "BREAKPOINTS FOUND:   ";
    print $total_breakpoints;
    print "\n";
    print "SEGMENTS WITH BREAKPOINTS:   ";
    print scalar @breakpoints;
    print "\n";
} else {
    my $breakpoint_segments = scalar @breakpoints;
    print join "\t", ($mummer_file,
                      $total_breakpoints,
                      $breakpoint_segments);
    print "\n";
}

###############################################################################
# DRAW PLOTS

if ($do_plot) {
    # Variables for making plots
    my @vlines; # Vertical lines marking contig boundaries in fasta1
    my @hlines; # Horizontal lines marking contig boundaries in fasta2
    my %reflen_offsets;
    my %querylen_offsets;
    
    # get offsets for contig boundaries
    get_offsets (\%ref_len_hash, \%reflen_offsets);
    get_offsets (\%query_len_hash, \%querylen_offsets);
    
    # contig boundary lines
    contig_boundaries (\%reflen_offsets, \@vlines);
    contig_boundaries (\%querylen_offsets, \@hlines);
    
    my %write_hash;
    parse_mummer_coords(\%neighbor_hash, \%reflen_offsets, \%querylen_offsets, \%write_hash);
    
    # Segments corresponding to non-breakpoint alignments
    write_segment_files ("x0_unflagged.list", "0", "x0", \%write_hash);
    write_segment_files ("x1_unflagged.list", "0", "x1", \%write_hash);
    write_segment_files ("y0_unflagged.list", "0", "y0", \%write_hash);
    write_segment_files ("y1_unflagged.list", "0", "y1", \%write_hash);
    
    # Segments corresponding to breakpoint alignments
    write_segment_files ("x0_flagged.list", "1", "x0", \%write_hash);
    write_segment_files ("x1_flagged.list", "1", "x1", \%write_hash);
    write_segment_files ("y0_flagged.list", "1", "y0", \%write_hash);
    write_segment_files ("y1_flagged.list", "1", "y1", \%write_hash);
    
    # Horizontal and vertical contig boundary lines
    open(OUT1, ">", "hlines.list") or die ("$!");
    print OUT1 join "\n", @hlines;
    close(OUT1);
    open(OUT1, ">", "vlines.list") or die ("$!");
    print OUT1 join "\n", @vlines;
    close(OUT1);
    
    # Make the plot with R script
    system("Rscript $path_to_r_script $file_prefix.diagonals.pdf");
}

## SUBROUTINES ################################################################

sub write_segment_files {
    my ($file, $flag, $val, $href) = @_;
    open(OUT, ">", $file) or die ("$!");
    print OUT join "\n", @{${href}->{$flag}{$val}};
    close(OUT);
}

sub parse_mummer_coords {
    my ($href, $href_offset_ref, $href_offset_query, $href_out) = @_;
    foreach my $id (keys %$href) {
        my $contig_ref = ${href}->{$id}{"contig"}{"ref"};
        my $contig_query = ${href}->{$id}{"contig"}{"query"};
        my $refstart = ${href}->{$id}{"start"}{"ref"} + ${href_offset_ref}->{$contig_ref}{"start"};
        my $refend = ${href}->{$id}{"end"}{"ref"} + ${href_offset_ref}->{$contig_ref}{"start"};
        my $querystart = ${href}->{$id}{"start"}{"query"} + ${href_offset_query}->{$contig_query}{"start"};
        my $queryend = ${href}->{$id}{"end"}{"query"} + ${href_offset_query}->{$contig_query}{"start"};
        my $flag = ${href}->{$id}{"flag"};
        push @{${href_out}->{$flag}{"x0"}}, $refstart;
        push @{${href_out}->{$flag}{"x1"}}, $refend;
        push @{${href_out}->{$flag}{"y0"}}, $querystart;
        push @{${href_out}->{$flag}{"y1"}}, $queryend;
    }
}

sub contig_boundaries {
    my ($href_offsets, $aref_lines) = @_;
    foreach my $contig (sort {$a cmp $b} keys %$href_offsets) {
        push @$aref_lines, ${href_offsets}->{$contig}{"end"};
    }
}

sub get_offsets {
    my ($href, $href_out) = @_;
    my $counter = 0;
    foreach my $contig (sort {$a cmp $b} keys %$href) {
        my $start = $counter + 1;
        my $end = $counter + ${href}->{$contig};
        $counter += ${href}->{$contig};
        ${href_out}->{$contig}{"start"} = $start;
        ${href_out}->{$contig}{"end"} = $end;
    }
}

sub hash_mummer {
    my ($file, $href, $rlen_href, $qlen_href) = @_;
    open(IN, "<", $file) or die ("$!");
    my $counter=0;
    while (<IN>) {
        chomp;
        my @splitline = split "\t";
        my ($ref, $query, $rlen, $qlen) = (@splitline)[5,0,18,2];
        my ($startref, $endref, $startquery, $endquery);
        if ($splitline[9] < $splitline[8]) { # Swap so that startref is always before endref
            ($startref, $endref, $startquery, $endquery) = (@splitline)[9,8,7,6]
        } else {
            ($startref, $endref, $startquery, $endquery) = (@splitline)[8,9,6,7];
        }
        ${href}->{$ref}{$counter}{"start"} = $startref;
        ${href}->{$ref}{$counter}{"end"} = $endref;
        ${href}->{$query}{$counter}{"start"} = $startquery;
        ${href}->{$query}{$counter}{"end"} = $endquery;
        # Hash lengths of reference sequences by name
        ${rlen_href}->{$ref} = $rlen unless defined ${rlen_href}->{$ref};
        ${qlen_href}->{$query} = $qlen unless defined ${qlen_href}->{$query};
        $counter++;
    }
    close(IN);
}

sub getneighbors_ref {
    my ($href, $href_out, $contig, $type) = @_;
    # $type specifies either "query" or "ref"
    my @sorted_ids = sort { ${${href}->{$contig}}{$a}{"start"} <=> ${${href}->{$contig}}{$b}{"start"} } keys %{${href}->{$contig}};
    my $last_index = $#sorted_ids;
    my $index = 0;
    foreach my $id (@sorted_ids) {
        # Check whether segment is reversed orientation
        my $rev = 0;
        if (${href}->{$contig}{$id}{"start"} > ${href}->{$contig}{$id}{"end"}) { # Deal with reverse direction
            $rev = 1;
        }
        # Check for first and last segment
        # TO DO: Implement a cutoff for max distance from contig start/end to count as a boundary
        if ($index ==0 && $last_index == 0) { # Single alignment on contig
            ${href_out}->{$id}{"leftneighbor"}{$type} = "BOUNDARY";
            ${href_out}->{$id}{"rightneighbor"}{$type} = "BOUNDARY";
        } elsif ($index==0 && $rev==0) {
            ${href_out}->{$id}{"leftneighbor"}{$type} = "BOUNDARY";
            ${href_out}->{$id}{"rightneighbor"}{$type} = $sorted_ids[$index+1];
        } elsif ($index==0 && $rev==1) {
            ${href_out}->{$id}{"rightneighbor"}{$type} = "BOUNDARY";
            ${href_out}->{$id}{"leftneighbor"}{$type} = $sorted_ids[$index+1];
        } elsif ($index==$last_index && $rev==0) {
            ${href_out}->{$id}{"rightneighbor"}{$type} = "BOUNDARY";
            ${href_out}->{$id}{"leftneighbor"}{$type} = $sorted_ids[$index-1];
        } elsif ($index==$last_index && $rev==1) {
            ${href_out}->{$id}{"leftneighbor"}{$type} = "BOUNDARY";
            ${href_out}->{$id}{"rightneighbor"}{$type} = $sorted_ids[$index-1];
        } else {
            # Non-boundary cases
            my ($leftindex, $rightindex);
            if ($rev==0) {
                $leftindex = $index - 1;
                $rightindex = $index + 1;
            } elsif ($rev==1) {
                $leftindex = $index + 1;
                $rightindex = $index - 1;
            }
            ${href_out}->{$id}{"leftneighbor"}{$type} = $sorted_ids[$leftindex];
            ${href_out}->{$id}{"rightneighbor"}{$type} = $sorted_ids[$rightindex];
        }
        ${href_out}->{$id}{"contig"}{$type} = $contig;
        ${href_out}->{$id}{"start"}{$type} = ${href}->{$contig}{$id}{"start"};
        ${href_out}->{$id}{"end"}{$type} = ${href}->{$contig}{$id}{"end"};
        # Test output
        #print join "\t", ($contig,
        #                  $id,
        #                  $rev,
        #                  $type,
        #                  ${href}->{$contig}{$id}{"start"},
        #                  ${href_out}->{$id}{"leftneighbor"}{$type},
        #                  ${href_out}->{$id}{"rightneighbor"}{$type}
        #                  );
        #print "\n";
        # Increment index
        $index++;
    }
}

sub get_breakpoints {
    # Find breakpoints by comparing left and right neighbors
    my ($href) = @_;
    my @breakpoint_ids;
    foreach my $id (sort {$a <=> $b} keys %$href) {
        my $flag = 0;
        if (${href}->{$id}{"leftneighbor"}{"ref"} ne ${href}->{$id}{"leftneighbor"}{"query"}) {
            # Left neighbor does not match between reference and query
            unless (${href}->{$id}{"leftneighbor"}{"ref"} eq "BOUNDARY" || ${href}->{$id}{"leftneighbor"}{"query"} eq "BOUNDARY") {
                $flag++;
            }
        }
        if (${href}->{$id}{"rightneighbor"}{"ref"} ne ${href}->{$id}{"rightneighbor"}{"query"}) {
            # Right neighbor does not match between reference and query
            unless (${href}->{$id}{"rightneighbor"}{"ref"} eq "BOUNDARY" || ${href}->{$id}{"rightneighbor"}{"query"} eq "BOUNDARY") {
                $flag++;
            }
        }
        #print join "\t", ($id,
        #                  ${href}->{$id}{"contig"}{"ref"},
        #                  ${href}->{$id}{"contig"}{"query"},
        #                  ${href}->{$id}{"start"}{"ref"},
        #                  ${href}->{$id}{"end"}{"ref"},
        #                  ${href}->{$id}{"start"}{"query"},
        #                  ${href}->{$id}{"end"}{"query"},
        #                  $flag
        #                  );
        #print "\n";
        if ($flag > 0) {
            push @breakpoint_ids, $id;
        }
        ${href}->{$id}{"flag"} = $flag;
    }
    return (@breakpoint_ids);
}
