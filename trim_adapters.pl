#!/usr/bin/perl -w
my $help_string = <<_HELP_STRING_;
#
#   trim_adapters.pl - Trim adapters out of Trinity Transcriptome assembly
#
# SYNOPSIS
#  Required command line arguments:
#   -infile       transcriptome assembly file made by Trinity and possibly
#                 reduced/merged from multiple libraries.  Fasta files must have
#                 one line for header, one line for sequence.
#   -scanfile     Output of blastn scan of UniVec with infile, see command below
#   -outfile      Name of output file
#   -adapters     Text, space delimited set of adapters to operate on, like "NGB00360:len NGB00362:len"
#
#  Optional command line arguments:
#   -pepfile      OPTIONAL - TransDecoder peptide file corresponding to infile.
#                   Only used to determine if trimming out the adapter fell within
#                   the CDS of the protein.  The Transdecoder cds file may also be used.
#   -minlen       Drop a transcript if it is trimmed below this length.
#   -drop_nopep   Drop a transcript if it does not have a protein, requires pepfile
#   -help         Print help
#
# Example trim run:
#   ./trim_adapters.pl \
#    -infile testing.fasta \
#    -outfile testing_out.fasta \
#    -scanfile univec_vs_testing.fmt6.out \
#    -pepfile testing.pep \
#    -adapters "NGB00360:58 NGB00362:61" \
#    -minlen 200 
#
# Example command to calculate scanfile (derived from NCBI TSA submission page):
#     nice blastn -task blastn -reward 1 -penalty -5 \
#        -gapopen 3 -gapextend 3 -dust yes -soft_masking true \
#        -evalue 700 -searchsp 1750000000000 -outfmt 6 \
#        -db \$PATH_TO_UNIVEC/UniVec -query testing.fasta \
#        -num_threads 40 -out testing.fmt6.out \
#        >testing.fmt6.log 2>&1 &
#
# Warning - the blastn command only gives one hit for an adapter even
#  if multiple copies are present.  So if a sequence looks like:
#    -------------------------->  transcript
#    -A> -B>                      two copies of the adapter
#  and blastn returns the A alignment the B adapter will still
#  be present.  If this happens run another cycle of this script.
#
# Warning - contaminants other than Illumina adapters may be present
#  and this tool will probably not be appropriate for removing them.
#
#
_HELP_STRING_
;
# Author
#   David Mathog <mathog\@caltech.edu>
#
# License
#   GPL 3
# Copyright (c) 2018, David Mathog & Caltech
#  All rights reserved.
#
#
# 2018-10-04 initial coding
## WARNING!  VERY LITTLE ERROR CHECKING!!!!!
#
#

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use List::Util  qw(min max);
use POSIX qw(strftime);
use FileHandle;

#command line parameters
my $infile;
my $scanfile;
my $pepfile;
my $outfile;
my $adapters;
my $help;
my $minlen=0;
my $drop_nopep;
my %match_set;

#globals
my %tr_lens;
my %tr_cds_ranges;
my %trim_set_lo;
my %trim_set_hi;
my $use_pep=0;
my $close_end=15;   #one end of alignment with the adapter must be within 15 bp of a transcript end.


GetOptions ("infile=s"     => \$infile,
            "outfile=s"    => \$outfile,
            "scanfile=s"   => \$scanfile,
            "pepfile=s"    => \$pepfile,
            "adapters=s"   => \$adapters,
            "minlen=i"     => \$minlen,
            "drop_nopep"   => \$drop_nopep,
            "help"         => \$help,
            );

if($help){
   print $help_string;
   exit;
}

select(STDOUT);             ##  
$| = 1; #disable buffering  ##  
select(STDERR);             ##  
$| = 1; #disable buffering  ##  

print_timestamp_message("Starting");

if(!defined($infile)){    die "-infile not specified";  }
if(!defined($outfile)){   die "-outfile not specified";  }
if(!defined($scanfile)){  die "-scanfile not specified";  }
if(!defined($adapters)){  die "-adapters not specified";  }

#############
print "Command line arguments were";
print " infile: $infile";
print " outfile: $outfile";
print " scanfile: $scanfile";
print " adapters: $adapters";
if(defined($pepfile)){
   print " pepfile: $pepfile";
   $use_pep=1;
}
print " minlen: $minlen";
print " drop_nopep: " . (defined($drop_nopep) ? "yes": "no");
print "\n";

#############
# store all adapter names to match in match_set
my @the_set = split(/\s+/, $adapters);
foreach(@the_set){
   my ($aname,$alen) = split(/\:/, $_);
   $match_set{$aname}=$alen;
}

#############
#get the lengths of all the transcripts
print_timestamp_message("Obtaining sequence lengths from: $infile");
open(IFILE,"<","$infile") or die "could not open $infile\n";
while(<IFILE>){
   my $line = $_;
   chomp($line);
   if(substr($line,0,1) eq '>'){
      my ($front_part, $len_part, $rest) = split(/\s+/, $line);
      my $tname = substr($front_part,1);
      my ($drop, $seqlen) = split(/\=/, $len_part);
      $tr_lens{$tname}=$seqlen;
   }
}
close(IFILE);

#############
#get the CDS positions from the pepfile
my $ranges=0;
my $multiples=0;
if(!$use_pep){
   print_timestamp_message("Not processing pepfile: (none specified)");
}
else {
   print_timestamp_message("Processing pepfile: $pepfile");
   open(PFILE,"<","$pepfile") or die "could not open $pepfile\n";
   while(<PFILE>){
      my $line = $_;
      chomp($line);
      if(substr($line,0,1) eq '>'){
         $ranges++;
         my ($front_piece, $remainder) = split(/\s+/, $line);
         my ($tname,$pnum) = split(/\./,substr($line,1));
         my @all_pieces = split(/\:/, $line);
         my $end_piece = $all_pieces[-1];
         my ($cds_range, $ignore) = split(/\(/, $end_piece);
         if(defined($tr_cds_ranges{$tname})){
            $tr_cds_ranges{$tname} .= " $cds_range";
            $multiples++;
         }
         else {
            $tr_cds_ranges{$tname}=$cds_range;
         }
      }
   } #while
   print_timestamp_message("CDS ranges:$ranges, 2nd or higher CDS range for transcript:$multiples");
   close(PFILE);
} #if $use_pep

#############
#process scan file

print_timestamp_message("Processing scanfile: $scanfile");
open(SFILE,"<","$scanfile") or die "could not open $scanfile\n";
while(<SFILE>) {
   my $line = $_;
   chomp($line);
   my ($tname, $vname, $ident, $len, $mismatch, $indel, $ts, $te, $vs, $ve, $ignore) = split(/\s+/, $line);
   my ($front_part, $back_part) = split(/\./, $vname);
   my $key = substr($front_part,3);
   #
   #  For Illumina adapters the end SHOULD be like this:
   #  ---------------> transcript
   #  --->       <---  adapters
   #  So if direction is flipped it is the high limit
   #  Adapters are accepted if they are within 10bp of one end of the transcript
   #     or the 3' end of the adapter is within 2bp of full length.
   #  Otherwise the adapter match is ignored (probably a weak similarity to a sequence).
   #
   my $klen = $match_set{$key};
   if(defined($klen)){
      if($ts<$close_end){ #end by position on transcript
         assign_lo($tname, max($ts,$te)+1); #coord is first base to KEEP, in 1->N
      }
      elsif($tr_lens{$tname} - $te <$close_end){ #end by position on transcript
         assign_hi($tname, min($ts,$te)-1); #coord is last base to KEEP, in 1->N
      }
      elsif($klen - 2 <= max($ve,$vs)){
         #somewhere internal, decide which side to trim by the orientation of the adapter
         if(($te-$ts)*($ve-$vs) < 0){
            assign_hi($tname, min($ts,$te)-1); #coord is last base to KEEP, in 1->N
         }
         else {
            assign_lo($tname, max($ts,$te)+1); #coord is first base to KEEP, in 1->N
         }
      }
   }
}
close(SFILE);

#############
#process input file
print_timestamp_message("Processing $infile");

my $keep_hi;
my $keep_lo;
my $newlen;
my $trimmed=0;
my $lt200=0;
my $gt100delta=0;
my $Nseqs=0;
my $cds_cut=0;
my $cds_drop=0;
my $len_drop=0;
my $noprot_drop=0;
my $drop_this=0;  #set if an entire CDS is removed or the transcript is too short
open(IFILE,"<","$infile") or die "could not open $infile\n";
open(OFILE,">","$outfile") or die "could not create output sequence file $outfile\n";
OUTER: while(<IFILE>) {
   my $line = $_;
   chomp($line);
   if(substr($line,0,1) eq '>'){
      $Nseqs++;
      #trinity transcript headers look like:
      #>SRR531950_TRINITY_DN15052_c0_g1_i5 len=297 path=[6:0-25 8:26-59 9:...
      #------fpart------------------------ -lpart- ---rest--->
      my ($front_part, $len_part, $rest) = split(/\s+/, $line);
      my $tname = substr($front_part,1);
      my ($drop, $seqlen) = split(/\=/, $len_part);
      $keep_lo = $trim_set_lo{$tname};
      $keep_hi = $trim_set_hi{$tname};
      
      if(!defined($keep_lo)){  $keep_lo=1;       }
      if(!defined($keep_hi)){  $keep_hi=$seqlen; }
      $newlen = $keep_hi - $keep_lo + 1;
      my $delta = $seqlen-$newlen;
      if($newlen < $minlen){
         $len_drop++;
         $drop_this=1;
         next OUTER;
      }
      if($use_pep){
         my $all_ranges = $tr_cds_ranges{$tname};
         if(defined($all_ranges)){ #a protein might not have been called for this transcript
            foreach my $this_range (split(/\s+/, $all_ranges)) {
               my ($cds_lo, $cds_hi) = split(/\-/,$this_range);
               if($cds_lo > $keep_hi || $cds_hi < $keep_lo){
                  $cds_drop++;
                  $drop_this=1;
                  next OUTER;
               }
               elsif($cds_lo < $keep_lo || $cds_hi > $keep_hi){
                  $cds_cut++;
               }
            }
         }
         elsif(defined($drop_nopep)){
            $noprot_drop++;
            $drop_this=1;
            next OUTER;
         }
      }
      # only counted for transcripts which were not dropped
      if($delta > 100){        $gt100delta++;    }
      if($seqlen != $newlen){  $trimmed++;       }
      if($newlen<200){         $lt200++;         }
      print OFILE "$front_part len=$newlen $rest\n";
   }
   else {
      if($drop_this){
         $drop_this=0;
         next OUTER;
      }
      print OFILE substr($line,$keep_lo-1, $newlen);
      print OFILE "\n";
   }
}
close(IFILE);
my $wrote=$Nseqs - $cds_drop - $len_drop - $noprot_drop;
print_timestamp_message("Done  Read:$Nseqs Wrote:$wrote Trimmed:$trimmed newlen<200:$lt200 deltaLen>100:$gt100delta cds_cut:$cds_cut cds_drop:$cds_drop len_drop:$len_drop noprot_drop:$noprot_drop");
exit;

sub assign_lo{
   my ($tname,$val) = @_;
   if(defined($trim_set_lo{$tname})){  #this is very rare
      if($trim_set_lo{$tname} > $val){ 
         $val = $trim_set_lo{$tname};
      }
   }
   $trim_set_lo{$tname}=$val;
}

sub assign_hi{
   my ($tname,$val) = @_;
   if(defined($trim_set_hi{$tname})){  #this is very rare
      if($trim_set_hi{$tname} < $val){ 
         $val = $trim_set_hi{$tname};
      }
   }
   $trim_set_hi{$tname}=$val;
}

sub print_timestamp_message{
   my ($message) = @_;
   my $now = strftime('%Y-%m-%d %H:%M:%S',localtime);
   print "$now $message\n";
}
