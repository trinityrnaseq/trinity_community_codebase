#!/usr/bin/perl -w
my $help_string = <<_HELP_STRING_;
#
#   trim_adapters.pl - Trim adapters out of Trinity Transcriptome assembly
#
# SYNOPSIS
#  Required command line arguments:
#   -infile       transcriptome assembly file made by Trinity and possibly
#                 reduced/merged from multiple libraries.  Fasta files must have
#                 one line for header, one line for sequence.  Unless -inranges
#                 is used the header must be verbatim as Trinity created it.  When
#                 -inranges is used it must have a name, space, and then string
#                 is allowed after that.
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
#   -inranges     If the infile is the result of previous processing there may be a ".ranges"
#                 file describing any previous edits.  When this is specified trims are applied
#                 to those preexisting values and a modifed ranges sent to -outranges. 
#                 Do NOT combine with -pepfile.  Output fasta file headers will be identical to the 
#                 input (other than those which are dropped) but the sequence will be trimmed
#                 if needed.  Trims will be reflected in -outranges start/stop fields.
#   -outranges    Required if -inranges is used, may be used without -inranges.
#   -ranges_fmt   Printf string for header fields: name start stop original_length remainder.  
#                 Default: '%-40s %6d %6d %6d %s'.  (No line feed character should be used!)
#   -relaxed      Normally the alignment found must extend to within the a small distance of 
#                 the end of the mRNA.  When this is set if the aligment is within 10bp of the
#                 3' end of the adapter an implied 5' end of the alignment, made by offsetting 
#                 by the adapter length is also tested.  This will find adapters in inaccurate
#                 mRNA sequence, where the sequence has been somewhat mangled.
#   -help         Print help
#
# Example trim runs:
#   ./trim_adapters.pl \
#    -infile testing.fasta \
#    -outfile testing_out.fasta \
#    -scanfile univec_vs_testing.fmt6.out \
#    -pepfile testing.pep \
#    -adapters "NGB00360:58 NGB00362:61" \
#    -minlen 200 
#
#   ./trim_adapters.pl \
#    -infile testing.fasta \
#    -outfile testing_out.fasta \
#    -inranges testing.ranges \
#    -ouranges testing.trimmed.ranges \
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
#  IMPORTANT!!!!! 
#  The UniVec database should be formatted with this or the equivalent:
#     formatdb -p F -i UniVec -o T
#     makeblastdb -in UniVec -dbtype nucl -parse_seqids 
#  So that the name for the entry name shown in the blast result has this format:
#    "uv|NGB00360.1:1-58"
#  Other formatting options MAY work as long as the first "uv" precedes the adapter
#  name by a single character.  If the "uv" is not present in the blastn output this
#  script will NOT work.
#
#  ranges data records are space delimited and look like:
#
#     name  start  stop original_length optional_data
#
#  Where start/stop define where in the original Trinity sequence the current sequence was located.
#  The length of the original sequence is the 4th value.  Positions are 1 -> original_length.
#  Optional data can be any string.  Any number of comments may be present at the beginning
#  of the file, these are indicated by a pound sign first character.
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
# 2019-01-24 updated,  added -relaxed, more conservative trim (by a few bases, in some cases)
# 2019-01-22 updated,  More robust handling of different blast database formats, support for "ranges" to keep
#     track of changes.
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
my $inranges;
my $outranges;
my $help;
my $minlen=0;
my $drop_nopep;
my %match_set;
my $ranges_fmt;
my $relaxed=0;

#globals
my %tr_original_slen;  # length of the transcript as made by Trinity
my %tr_current_start;  # after some previous processing, where this transcript starts (1->N)
my %tr_current_stop;   # after some previous processing, where this transcript starts (1->N)
my %tr_cds_ranges;     # used with pepfile
my %tr_remainder;      # used with inranges
my %trim_set_lo;
my %trim_set_hi;
my $process_type;
my $USE_NORMAL=0;
my $USE_PEP=1;
my $USE_RANGES=2;
my $close_mRNA_end=15;   #one end of alignment with the adapter must be within 15 bp of a transcript end.
my $close_adapter_end=10; #the alignment extends to 10 or less bp of 3' end of adapter, used with -relaxed
my $retained_ranges_comments="";  

GetOptions ("infile=s"     => \$infile,
            "outfile=s"    => \$outfile,
            "scanfile=s"   => \$scanfile,
            "pepfile=s"    => \$pepfile,
            "inranges=s"   => \$inranges,
            "outranges=s"  => \$outranges,
            "adapters=s"   => \$adapters,
            "ranges_fmt=s" => \$ranges_fmt,
            "minlen=i"     => \$minlen,
            "relaxed"      => \$relaxed,
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
if(defined($pepfile) && (defined($inranges) || defined($outranges))){
   die "-pepfile cannot be combined with -[in|out]ranges";
}
if(defined($inranges) && !defined($outranges)){
   die "-inranges requres -outranges (but not vice versa)";
}
if(!defined($ranges_fmt)){
   $ranges_fmt='%-40s %6d %6d %6d %s';
}

emit_parameters_summary();

get_adapters();

if($process_type == $USE_RANGES){
   get_range_data();
}
else {
   transcript_lengths_from_input();
}

if($process_type == $USE_PEP){
   get_CDS_from_pepfile();
}

process_scan_file();

process_infile();

exit;

sub emit_parameters_summary{
   print "Command line arguments were";
   print " infile: $infile";
   print " outfile: $outfile";
   print " scanfile: $scanfile";
   print " adapters: $adapters";
   if(defined($pepfile)){
      print " pepfile: $pepfile";
      $process_type = $USE_PEP;
   }
   elsif(defined($inranges)){
      print " inranges: $inranges";
      print " outranges: $outranges";
      $process_type=$USE_RANGES;
   }
   else {
      $process_type=$USE_NORMAL;
   }
   print " minlen: $minlen";
   print " relaxed: $relaxed";
   print " drop_nopep: " . (defined($drop_nopep) ? "yes": "no");
   print " ranges_fmt: $ranges_fmt";
   print "\n";
}



#############
# store all adapter names to match in match_set
sub get_adapters{
   my @the_set = split(/\s+/, $adapters);
   foreach(@the_set){
      my ($aname,$alen) = split(/\:/, $_);
      $match_set{$aname}=$alen;
   }
}

#############
# get range data from -inranges.
sub get_range_data{
   print_timestamp_message("Obtaining range data from: $inranges");
   open(RFILE,"<","$inranges") or die "could not open $inranges\n";
   my $idx=0;
   while(<RFILE>){
     my $line = $_;
     if(substr($line,0,1) eq '#'){  #keep all comments, 
        $retained_ranges_comments .= $line;
     }
     chomp($line);
     my ($tname, $rest) = split(/\s+/, $line, 2);
     ($tr_current_start{$tname}, $tr_current_stop{$tname},
      $tr_original_slen{$tname}, $tr_remainder{$tname}) = split(/\s+/, $rest, 4);
   }
   close(RFILE);
   
}

#############
#get the lengths of all the transcripts using the values stored in the
#header.  Header syntax must be straight from Trinity and claimed length
#must match actual length.
sub transcript_lengths_from_input{
   print_timestamp_message("Obtaining sequence lengths from: $infile");
   open(IFILE,"<","$infile") or die "could not open $infile\n";
   my $tname="";
   my $seqlen=0;
   while(<IFILE>){
     my $line = $_;
     chomp($line);
     if(substr($line,0,1) eq '>'){
        my ($front_part, $rest) = split(/\s+/, $line);
        $tname = substr($front_part,1);
        my $len_pos = index($line," len=");
        if($len_pos == -1){ die "fatal error: Trinity fasta header lacks len= component.  Line:\n$line"; }
        ($seqlen, my $drop) = split(/\s+/, substr($line,$len_pos+5));
        $tr_original_slen{$tname}=$seqlen;
        $tr_current_start{$tname}=1;
        $tr_current_stop{$tname}=$seqlen;
     }
     else {
        my $real_seqlen = length($line);
        if($real_seqlen != $seqlen){
           die "Fatal input error: $tname has Len=$seqlen but actual length is $real_seqlen\n";
        }
     }
   }
   close(IFILE);
}

#############
#get the CDS positions from the pepfile
sub get_CDS_from_pepfile{
  my $ranges=0;
  my $multiples=0;
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
}

#############
#process scan file
sub process_scan_file{
   print_timestamp_message("Processing scanfile: $scanfile");
   open(SFILE,"<","$scanfile") or die "could not open $scanfile\n";
   my $cut_lo_end=0;
   my $cut_lo_mid=0;
   my $cut_hi_end=0;
   my $cut_hi_mid=0;
   while(<SFILE>) {
      my $line = $_;
      chomp($line);
      my ($tname, $vname, $ident, $len, $mismatch, $indel, $ts, $te, $vs, $ve, $ignore) = split(/\s+/, $line);
      my ($front_part, $back_part) = split(/\./, $vname);
      #depending on how the database was formatted, it might see either of these forms for $front_part:
      #  gnl|uv|NGB00360.1:1-58  uv:NGB00360.1:1-58 uv|NGB00360.1:1-58 or maybe something else.  uv always seems to be last.
      my $first_uv = index($front_part,'uv');
      if($first_uv == -1 ){ die "fatal error: adapter name in blastn output does not look like \"<something>uv|name.<something>\" Is:\n$line"; }
      my $key = substr($front_part,$first_uv + 3);
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
      my $current_seqlen = $tr_current_stop{$tname} - $tr_current_start{$tname} + 1;
      my $dir = ((($te-$ts)*($ve-$vs) < 0) ? 1 : 0);  # alignment direction, 0 is forward, 1 is reversed
      if($relaxed && ($klen - max($ve,$vs) + 1 <= $close_adapter_end)){
         #extrapolate to new $ts, $te, as if the aligmment was perfect rather than just the smaller part from the scan
         if($dir){
            $ts = $ts - ($klen - $vs);
            $te = min($current_seqlen, $ts + $klen - 1);
         }
         else {
            $te = $te + ($klen - $ve);
            $ts = max(1,$te - $klen + 1);
         }
      }
      if(defined($klen)){
         if($ts<$close_mRNA_end){ #end by position on transcript
            assign_lo($tname, max($ts,$te)+1); #coord is first base to KEEP, in 1->N
            $cut_lo_end++;
         }
         elsif($current_seqlen - $te < $close_mRNA_end){ #end by position on transcript (or remaining part of transcript)
            assign_hi($tname, min($ts,$te)-1); #coord is last base to KEEP, in 1->N
            $cut_hi_end++;
         }
         elsif($klen - 2 <= max($ve,$vs)){
#print "DEBUG internal klen $klen vs $vs ve $ve ts $ts te $te seqlen $current_seqlen\n";
            #somewhere internal, decide which side to trim by the orientation of the adapter
            if($dir){
               assign_hi($tname, min($ts,$te)-1); #coord is last base to KEEP, in 1->N
               $cut_hi_mid++;
            }
            else {
               assign_lo($tname, max($ts,$te)+1); #coord is first base to KEEP, in 1->N
               $cut_lo_mid++;
            }
         }
      }
   }
   print_timestamp_message("adapter cuts at: end (lo $cut_lo_end hi $cut_hi_end ) mid (lo $cut_lo_mid hi $cut_hi_mid )"); 
   close(SFILE);
}

#############
#process input file

sub process_infile{
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
   my $trimmed_lo=0;
   my $trimmed_hi=0;
   my $trimmed_both=0;
   open(IFILE,"<","$infile") or die "could not open $infile\n";
   open(OFILE,">","$outfile") or die "could not create output sequence file $outfile\n";
   if($process_type == $USE_RANGES){
      open(RFILE,">","$outranges") or die "could not create output ranges file $outranges\n";
      print RFILE "$retained_ranges_comments";
   }
   OUTER: while(<IFILE>) {
      my $line = $_;
      chomp($line);
      if(substr($line,0,1) eq '>'){
         $Nseqs++;
         #trinity transcript headers look like:
         #>SRR531950_TRINITY_DN15052_c0_g1_i5 len=297 path=[6:0-25 8:26-59 9:...
         #------fpart------------------------ -lpart- ---rest--->
         my ($front_part, $ignore_lpart, $rest) = split(/\s+/, $line, 3);
         my $tname = substr($front_part,1);
         my $current_seqlen = $tr_current_stop{$tname} - $tr_current_start{$tname} + 1;
         $keep_lo = $trim_set_lo{$tname};
         $keep_hi = $trim_set_hi{$tname};

         if(!defined($keep_lo)){  $keep_lo=1;       }
         if(!defined($keep_hi)){  $keep_hi=$current_seqlen; }
         $newlen = $keep_hi - $keep_lo + 1;
         if($keep_lo > 1 ){              $trimmed_lo++;   }
         if($keep_hi < $current_seqlen){ $trimmed_hi++;   }
         if($keep_lo > 1  && $keep_hi < $current_seqlen){ 
                                         $trimmed_both++; }
         
         my $delta = $current_seqlen-$newlen;
         if($newlen < $minlen){
            $len_drop++;
            $drop_this=1;
            next OUTER;
         }
         if($process_type == $USE_PEP){
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
            print OFILE "$front_part len=$newlen $rest\n";
         }
         elsif($process_type == $USE_RANGES){ 
           # fasta header is not modified in any way
           # changes, if any are stored in new ranges file
            print OFILE "$line\n";
            my $new_start = $tr_current_start{$tname} + $keep_lo - 1;  # coordinates in 1->N, changes if trimmed
            my $new_stop  = $tr_current_start{$tname} + $keep_hi - 1;  # coordinates in 1->N, changes if trimmed
            printf RFILE $ranges_fmt,$tname,$new_start,$new_stop,$tr_original_slen{$tname},"$tr_remainder{$tname}\n";
         }
         else {
            print OFILE "$front_part len=$newlen $rest\n";
         }
         # only counted for transcripts which were not dropped
         if($delta > 100){                $gt100delta++;    }
         if($current_seqlen != $newlen){  $trimmed++;       }
         if($newlen<200){                 $lt200++;         }
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
   if($process_type == $USE_RANGES){
      close(RFILE);
   }
   my $wrote=$Nseqs - $cds_drop - $len_drop - $noprot_drop;
   print_timestamp_message("Done  Read:$Nseqs Wrote:$wrote Trimmed:$trimmed newlen<200:$lt200 deltaLen>100:$gt100delta cds_cut:$cds_cut cds_drop:$cds_drop len_drop:$len_drop noprot_drop:$noprot_drop trim_lo:$trimmed_lo trim_hi:$trimmed_hi trim_both:$trimmed_both");
}

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
