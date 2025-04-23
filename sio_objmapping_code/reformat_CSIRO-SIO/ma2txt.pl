#!/usr/bin/perl
use strict;

# my $version = "Version 1.1 Nov 23, 2009";
# my $version = "Version 1.2 May 20, 2011";
# Version 1.2: Modifications to qcdoc section to clarify
my $version = "Version 1.3 May 25, 2011";
# Version 1.3: Handles QC hexadecimals up to 8 hex digits
# 27nov2017 LL chop up to write only data I want...
#
print "gtspp2txt.pl - Converts MEDS-ASCII to Text or CSV\n";
print "$version\n";
######
sub qcdoc {
print <<EOFDOC;
----------
QC Surface Code Groups using standard tests & bit assignments:
The bit is set to "1" if the test has been performed,
or if the station has failed the test.

  MEDS (Using QC2_LIB.FOR):
    QCP\$: QC tests "Performed"
    QCF\$: QC tests "Fail"
  APDRC (Using qced):
    QP1\$: QC tests "Performed"
    QT1#: QC tests "Fail" for TEMP
    QP1#: QC tests "Fail" for PSAL
  NODC (Using qced):
    QNP\$: QC tests "Performed"
    QNF\$: QC tests "Fail"
    QTE#: QC tests "Fail" for TEMP
    QPS#: QC tests "Fail" for PSAL
  Other than NODC or Science QC Center, but using NODC's qced (IDL):
    QP9\$: QC tests "Performed"
    QT9#: QC tests "Fail" for TEMP
    QP9#: QC tests "Fail" for PSAL

    HEXADECIMAL:  0003FFFFFF (10 CHARACTERS = 40 BITS; ONLY LOW 26 BITS USED)
      BIT VALUE:  0000 0000 0000 0011 1111 1111 1111 1111 1111 1111
 HEX CHAR VALUE:  0    0    0    3    F    F    F    F    F    F
          BYTE#:  10   9    8    7    6    5    4    3    2    1
           BIT#:  4321 4321 4321 4321 4321 4321 4321 4321 4321 4321
SEQUENTIAL BIT#:  {... .... .*.. .}{. .... .... .... ...9 8765 432}
                  {14 BITS UNUSED }{26 BITS (IN SEQ R to L) USED  }
* Bit 8,3 is used by ISDM (MEDS); Meaning unknown

Test Stages and Bit Assignments:
(DECIMAL is the base-10 value for the bit, = 2^(SEQ-1))
  STAGE  TEST                           CHAR BIT  SEQ  (DECIMAL)
  -----  -----------------------------  ---- ---  ---  ----------
         Location and Identification
   1.1     Platform Identification      1    1    1    (1)
   1.2     Impossible Date/Time         1    2    2    (2)
   1.3     Impossible Location          1    3    3    (4)
   1.4     Position on Land             1    4    4    (8)
   1.5     Impossible Speed             2    1    5    (16)
   1.6     Impossible Sounding          2    2    6    (32)
         Profile Tests
   2.1     Global Impossible Values     2    3    7    (64)
   2.2     Regional Impossible Values   2    4    8    (128)
   2.3     Increasing Depth             3    1    9    (256)
   2.4     Profile Envelope             3    2   10    (512)
   2.5     Constant Profile             3    3   11    (1024)
   2.6     Freezing Point               3    4   12    (2048)
   2.7     Spike                        4    1   13    (4096)
   2.8     Top and Bottom Spike         4    2   14    (8192)
   2.9     Gradient                     4    3   15    (16384)
   2.10    Density Inversion            4    4   16    (32768)
   2.11    Bottom                       6    4   24    (8388608)
   2.12    Temperature Inversion        7    1   25    (16777216)
         Climatology Tests
   3.1     Levitus Seasonal Statistics  5    1   17    (65536)
   3.2     Emery and Dewar Climatology  5    2   18    (131072)
   3.3     Asheville Climatology        5    3   19    (262144)
   3.4     Levitus Monthly Climatology  5    4   20    (524288)
   3.5     Levitus Annual Climatology   7    2   26    (33554432)
         Profile Consistency Tests
   4.1     Waterfall                    6    1   21    (1048576)
         Visual Inspection
   5.1     Cruise Track                 6    2   22    (2097152)
   5.2     Profiles                     6    3   23    (4194304)
         Undocumented
   UNK     UNDEFINED                    7    3   27    (67108864)
   UNK     UNDEFINED                    7    4   28    (134217728)
   UNK     UNDEFINED                    8    1   29    (268435456)
   UNK     UNDEFINED                    8    2   30    (536870912)
   UNK     UNDEFINED But used by ISDM   8    3   31    (1073741824)
   UNK     UNDEFINED                    8    4   32    (2147483648)
----------
EOFDOC
exit;
}

######
sub usage {
print <<EOF1;
----------
Usage: gtspp2txt.pl [options]
 Options flags may be in any order, flags may be upper or lower case
 Space between flag and parameter is optional
 Input filename may be provided without the -i flag
 If gzip is installed, will read compressed input filename ending in .gz
  -i filename = Input MEDS-ASCII file (may be gzipped)
  -o filename = Output file
  -unw        = Text output with QC Surface Codes 'unwound' (Default)
  -txt        = Text output with no interpretation of surface codes
  -csv        = Comma Separated Values output
  -doc        = prints description of QC group expansion
  -h          = Help (Prints this message)
Examples:
  gtspp2txt.pl -i input.meds -o output.txt -unw
    (writes output to file as text with QC interpretation)
  gtspp2txt.pl -i input.meds -o output.txt -txt
    (writes output as text with no QC interpretation)
  gtspp2txt.pl -i input.meds -o output.txt -csv
    (writes output as comma separated values)
  gtspp2txt.pl input.meds.gz
    (reads gzipped file, writes text output to screen with QC interpretation)
EOF1
exit;
}
#
######
# Interpret Command Line
#
my $infile = '';
my $outfile = '';
my $format = 'unw';
CMD: for (my $n = 0; $n <= $#ARGV; $n++) {
  if ($ARGV[$n] =~ /^-i$/i ) {
    $infile = $ARGV[$n+1]; $ARGV[$n]=""; $ARGV[$n+1]=""; next CMD;}
  if ($ARGV[$n] =~ /^-o$/i ) {
    $outfile = $ARGV[$n+1]; $ARGV[$n]=""; $ARGV[$n+1]=""; next CMD;}
  #
  if ($ARGV[$n] =~ /^-i(.+)$/i) {$infile = $1; $ARGV[$n]=""; next CMD;}
  if ($ARGV[$n] =~ /^-o(.+)$/i) {$outfile = $1; $ARGV[$n]=""; next CMD;}
  #
  if ($ARGV[$n] =~ /^-unw$/i) {$format = 'unw'; $ARGV[$n]=""; next CMD;}
  if ($ARGV[$n] =~ /^-txt$/i) {$format = 'txt'; $ARGV[$n]=""; next CMD;}
  if ($ARGV[$n] =~ /^-csv$/i) {$format = 'csv'; $ARGV[$n]=""; next CMD;}
  if ($ARGV[$n] =~ /^-doc$/i) {&qcdoc;}
  if ($ARGV[$n] =~ /^-h$/i)   {&usage;}
}

for (@ARGV) {
  if (/^-/) {print "UNKNOWN COMMAND LINE SWITCH ($_)\n"; &usage;}
  if (length($_)) {
    if (length($infile)) {
      print "UNKNOWN COMMAND LINE ARGUMENT ($_)\n";
      &usage;
    }
    else {$infile = $_;}
  }
}

if (length($infile)==0) {
  print "INPUT FILE NAME REQUIRED\n";
  &usage;
}
#
######
## DATA STRUCTURE -- GLOBAL VARIABLES
######
my $CSDPDEPH; #Flag for CSIRO depth correction
my $cor_dep_found; #Flag for actually finding depth corrected data
my ($sta, $MKey);
my ($station, $seq, $totalsegs);
my ($One_Deg_sq, $Cruise_ID, $Obs_Year, $Obs_Month, 
    $Obs_Day, $Obs_Time, $Data_Type, $Iumsgno, $Stream_Source, 
    $Uflag, $MEDS_Sta, $Latitude, $Longitude, $Q_Pos, $Q_Date_Time,
    $Q_Record, $Up_Date, $Bul_Time, $Bul_Header, $Source_ID,
    $Stream_Ident, $QC_Version, $Data_Avail, $No_Prof, $Nparms,
    $Nsurfc, $Num_Hists, $loadepoch, $active, $dmode, $oceancode);
my (@No_Seg, @Prof_Type, @Dup_flag, @Digit_Code, @Standard, @Deep_Depth);
my (@Pcode, @Parm, @Q_Parm);
my (@SRFC_Code, @SRFC_Parm, @SRFC_Q_Parm);
my (@Ident_Code, @PRC_Code, @Version, @PRC_Date, @Act_Code, @Act_Parm, 
    @Aux_ID, @Previous_Val);
my (@profile, @Profile_Type, @Profile_Seg, @No_Depths, @D_P_Code);
my (@Depth_Press, @Depres_Q, @Prof_Parm, @Prof_Q_Parm);
######


## ADDITIONAL VARIABLES FOR MEDS-ASCII READ/WRITE
######
my ($rec, $seg, $line, $reccount);
my @Prof_Inf;
my @SPGp;
my @SCGp;
my @HGp;
my @dataline;
my (@MKey_r, @One_Deg_sq, @Cruise_ID, @Obs_Year, @Obs_Month, 
  @Obs_Day, @Obs_Time, @Data_Type, @Iumsgno);
######
my $qcn;
my (@qclab, @bin);
######
#
print "Converting MEDS-ASCII ";
if ($format =~ /unw/) {print "to text (With interpretation of QC)\n";}
if ($format =~ /txt/) {print "to text (No interpretation of QC)\n";}
if ($format =~ /csv/) {print "to Comma Separated Values (csv) format\n";}
#
# READ STATION FROM FILE AND PRINT TO STDOUT OR OUTPUT FILE
#  Open compressed or non-compressed files
if ((substr($infile, length($infile)-3,3) eq '.gz') || (-B $infile) ){
  # GZIP file
  open(INFILE, "gzip -dc $infile | ") || die "Cannot open file $infile: $!\n";
  print "  Opened GZIP compressed file: $infile\n";
}
else {
  #regular ASCII file
  open(INFILE, $infile) || die "Cannot open file $infile: $!\n";
  print "  Opened non-compressed file: $infile\n";
}
#  Open output file and select output; or remain as STDOUT
my $old_fh;
if (length($outfile)) {
  open(OUTFILE, ">".$outfile) || die "Can't open output file '$outfile': $!\n";
  print "  Opened file for output: $outfile\n";
  $old_fh = select(OUTFILE);
}
else {print "  Output to Standard Output\n";}
#  Loop through the stations in the file
my $countstations = 0;
while (!eof(INFILE)) {
  &readsta;
# LL
  printf("%3d \n", $countstations+1 ); 
  if ($format =~ /unw|txt/) {&printtxt;}
  if ($format =~ /csv/) {&printcsv;}
  $countstations++;
}
close INFILE;
if (length($outfile)) {select($old_fh);}
close OUTFILE;
print "  Converted $countstations stations\n";
print "Done\n";
# END OF MAIN

######
sub readsta {
  $_ = <INFILE>;
  s/[\012\015]+$//; # Remove EOL characters
  $rec = 0;
  $seg = 0;
  $line++;
  $reccount++;

  $MKey          = substr($_,   0,  8);
  $One_Deg_sq    = substr($_,   8,  8);
  $Cruise_ID     = substr($_,  16, 10);
  $Obs_Year      = substr($_,  26,  4);
  $Obs_Month     = substr($_,  30,  2);
  $Obs_Day       = substr($_,  32,  2);
  $Obs_Time      = substr($_,  34,  4);
  $Data_Type     = substr($_,  38,  2);
  $Iumsgno       = substr($_,  40, 12);
  $Stream_Source = substr($_,  52,  1);
  $Uflag         = substr($_,  53,  1);
  $MEDS_Sta      = substr($_,  54,  8);
  $Latitude      = substr($_,  62,  8);
  $Longitude     = substr($_,  70,  9);
  $Q_Pos         = substr($_,  79,  1);
  $Q_Date_Time   = substr($_,  80,  1);
  $Q_Record      = substr($_,  81,  1);
  $Up_Date       = substr($_,  82,  8);
  $Bul_Time      = substr($_,  90, 12);
  $Bul_Header    = substr($_, 102,  6);
  $Source_ID     = substr($_, 108,  4);
  $Stream_Ident  = substr($_, 112,  4);
  $QC_Version    = substr($_, 116,  4);
  $Data_Avail    = substr($_, 120,  1);
  $No_Prof       = substr($_, 121,  2);
  $Nparms        = substr($_, 123,  2);
  $Nsurfc        = substr($_, 125,  2);
  $Num_Hists     = substr($_, 127,  3);
#
#
  my $B = 130;
  $totalsegs = 0;
  for (my $n = 0; $n < $No_Prof; $n++) {
    $Prof_Inf[$n] = substr($_, $B, 14);
    $B += 14;
    $No_Seg[$n] = substr($Prof_Inf[$n], 0, 2);
    $Prof_Type[$n] = substr($Prof_Inf[$n], 2, 4);
    $Dup_flag[$n] = substr($Prof_Inf[$n], 6, 1);
    $Digit_Code[$n] = substr($Prof_Inf[$n], 7, 1);
    $Standard[$n] = substr($Prof_Inf[$n], 8, 1);
    $Deep_Depth[$n] = substr($Prof_Inf[$n], 9, 5);
    $totalsegs += $No_Seg[$n];
  }
#
  for (my $n = 0; $n < $Nparms; $n++) {
    $SPGp[$n] = substr($_, $B, 15);
    $B += 15;
    $Pcode[$n] = substr($SPGp[$n], 0, 4);
    $Parm[$n]  = substr($SPGp[$n], 4, 10);
    $Q_Parm[$n] = substr($SPGp[$n], 14, 1);
  }
#
  for (my $n = 0; $n < $Nsurfc; $n++) {
    $SCGp[$n] = substr($_, $B, 15);
    $B += 15;
    $SRFC_Code[$n] = substr($SCGp[$n], 0, 4);
    $SRFC_Parm[$n] = substr($SCGp[$n], 4, 10);
    $SRFC_Q_Parm[$n] = substr($SCGp[$n], 14, 1);
  }
#
  for (my $n = 0; $n < $Num_Hists; $n++) {
    $HGp[$n] = substr($_, $B, 42);
    $B += 42;
    $Ident_Code[$n] = substr($HGp[$n], 0, 2);
    $PRC_Code[$n] = substr($HGp[$n], 2, 4);
    $Version[$n] = substr($HGp[$n], 6, 4);
    $PRC_Date[$n] = substr($HGp[$n], 10, 8);
    $Act_Code[$n] = substr($HGp[$n], 18, 2);
    $Act_Parm[$n] = substr($HGp[$n], 20, 4);
    $Aux_ID[$n] = substr($HGp[$n], 24, 8);
    $Previous_Val[$n] = substr($HGp[$n], 32, 10);
  }
# 
  $totalsegs = 0;
  for (my $rec = 0; $rec < $No_Prof; $rec++) {
    for (my $seg = 0; $seg < $No_Seg[$rec]; $seg++) {
      $_ = <INFILE>;
      s/[\012\015]+$//; # Remove EOL characters
      $dataline[$totalsegs] = $_;
      $line++;
      $MKey_r[$totalsegs]       = substr($_,  0,  8);
      $One_Deg_sq[$totalsegs]   = substr($_,  8,  8);
      $Cruise_ID[$totalsegs]    = substr($_, 16, 10);
      $Obs_Year[$totalsegs]     = substr($_, 26,  4);
      $Obs_Month[$totalsegs]    = substr($_, 30,  2);
      $Obs_Day[$totalsegs]      = substr($_, 32,  2);
      $Obs_Time[$totalsegs]     = substr($_, 34,  4);
      $Data_Type[$totalsegs]    = substr($_, 38,  2);
      $Iumsgno[$totalsegs]      = substr($_, 40, 12);
      $Profile_Type[$totalsegs] = substr($_, 52,  4);
      $Profile_Seg[$totalsegs]  = substr($_, 56,  2);
      $No_Depths[$totalsegs]    = substr($_, 58,  4);
      $D_P_Code[$totalsegs]     = substr($_, 62,  1);
#
      $B = 63;
      for (my $d = 0; $d < $No_Depths[$totalsegs]; $d++) {
        $Depth_Press[$totalsegs][$d] = substr($_, $B, 6); $B += 6;
        $Depres_Q[$totalsegs][$d]    = substr($_, $B, 1); $B += 1;
        $Prof_Parm[$totalsegs][$d]   = substr($_, $B, 9); $B += 9;
        $Prof_Q_Parm[$totalsegs][$d] = substr($_, $B, 1); $B += 1;
      }
      $totalsegs++;
    }
  }
}

######
sub printtxt {
#LL
  print
  "$Obs_Year\n",
  "$Obs_Month\n",
  "$Obs_Day\n",
  "$Obs_Time\n",
  "$Latitude\n",
  "$Longitude\n",
  "$Nsurfc\n";
#
  $qcn = 0;
  for (my $n = 0; $n < $Nsurfc; $n++) {
    printf ("%s    %s\n", $SRFC_Code[$n], $SRFC_Parm[$n]);
  }

  #
  my $tsegs = 0;
  for (my $rec = 0; $rec < $No_Prof; $rec++) {
# i want seg 1    for (my $seg = 0; $seg < $No_Seg[$rec]; $seg++) {
    for (my $seg = 0; $seg < 1; $seg++) {
      print
      "$No_Depths[$tsegs]\n";
      #
      #
      for (my $d = 0; $d < $No_Depths[$tsegs]; $d++) {
        printf ("%s  %s %s\n", $Depth_Press[$tsegs][$d],
        $Prof_Parm[$tsegs][$d], $Prof_Q_Parm[$tsegs][$d]);
      }
      $tsegs++;
    }
  }
#
  print
  "$Num_Hists\n";
#
  for (my $n = 0; $n < $Num_Hists; $n++) {
    printf ("%s   %s %s\n", $Act_Code[$n], $Aux_ID[$n], $Previous_Val[$n]);
  }
  print "END\n";
}

#####
# USAGE OF interpretqc:
# Requires global $qcn = number of interpretable QC groups
# Requires global @qclab = array of labels (srfccode)
# Requires global @bin = array of QC codes as binary strings
# Prints chart with Stage number, Test description, HexChar/Bit, Values=0/1
# One column of values for each interpretable QC group
#
sub interpretqc {
  print "\n  INTERPRETATION OF QC SURFACE CODE GROUPS:\n";
  print "  STAGE  TEST                         BIT";
  for (my $q = 0; $q < $qcn; $q++) {print "  $qclab[$q]";}
  print "\n";
  print "  -----  ---------------------------  ---";
  for (my $q = 0; $q < $qcn; $q++) {print "  ----";}
  print "\n";
  #
  my @line = ("zero",
  "   1.1   Platform Identification      1,1  1",
  "   1.2   Impossible Date/Time         1,2  2",
  "   1.3   Impossible Location          1,3  3",
  "   1.4   Position on Land             1,4  4",
  "   1.5   Impossible Speed             2,1  5",
  "   1.6   Impossible Sounding          2,2  6",
  "   2.1   Global Impossible Values     2,3  7",
  "   2.2   Regional Impossible Values   2,4  8",
  "   2.3   Increasing Depth             3,1  9",
  "   2.4   Profile Envelope             3,2 10",
  "   2.5   Constant Profile             3,3 11",
  "   2.6   Freezing Point               3,4 12",
  "   2.7   Spike                        4,1 13",
  "   2.8   Top and Bottom Spike         4,2 14",
  "   2.9   Gradient                     4,3 15",
  "   2.10  Density Inversion            4,4 16",
  "   2.11  Bottom                       6,4 24",
  "   2.12  Temperature Inversion        7,1 25",
  "   3.1   Levitus Seasonal Statistics  5,1 17",
  "   3.2   Emery and Dewar Climatology  5,2 18",
  "   3.3   Asheville Climatology        5,3 19",
  "   3.4   Levitus Monthly Climatology  5,4 20",
  "   3.5   Levitus Annual Climatology   7,2 26",
  "   4.1   Waterfall                    6,1 21",
  "   5.1   Cruise Track                 6,2 22",
  "   5.2   Profiles                     6,3 23",
  "   UNK   UNDEFINED                    7,3 27",
  "   UNK   UNDEFINED                    7,4 28",
  "   UNK   UNDEFINED                    8,1 29",
  "   UNK   UNDEFINED                    8,2 30",
  "   UNK   UNDEFINED                    8,3 31",
  "   UNK   UNDEFINED                    8,4 32");
  #0123456789012345678901234567890123456789012345

  # Loop through descriptive lines and test bits
  for (my $ln = 1; $ln <=32; $ln++) {
    print substr($line[$ln], 0, 41);
    my $bit = substr($line[$ln], 42, 2);
    # Loop through interpretable groups (as binary strings)
    for (my $q = 0; $q < $qcn; $q++) {
      print "     ", substr($bin[$q], 32-$bit, 1);
    }
    print "\n";
  }
}

#####
sub printcsv {
  print "MKey,One_Deg_sq,Cruise_ID,Obs_Year,Obs_Month,Obs_Day,Obs_Time,",
    "Data_Type,Iumsgno,Stream_Source,Uflag,MEDS_Sta,Latitude,Longitude,",
    "Q_Pos,Q_Date_Time,Q_Record,Up_Date,Bul_Time,Bul_Header,Source_ID,",
    "Stream_Ident,QC_Version,Data_Avail,No_Prof,Nparms,Nsurfc,Num_Hists\n";
  my @vals = &trim($MKey,$One_Deg_sq,$Cruise_ID,$Obs_Year,$Obs_Month,
    $Obs_Day,$Obs_Time,$Data_Type,$Iumsgno,$Stream_Source,$Uflag,$MEDS_Sta,
    $Latitude,$Longitude,$Q_Pos,$Q_Date_Time,$Q_Record,$Up_Date,$Bul_Time,
    $Bul_Header,$Source_ID,$Stream_Ident,$QC_Version,$Data_Avail,$No_Prof,
    $Nparms,$Nsurfc,$Num_Hists);
  print join(",", @vals), "\n";
  #
  print "Prof_Info_Gp#,No_Seg,Prof_Type,Dup_flag,Digit_Code,Standard,Deep_Depth\n";
  for (my $n = 0; $n < $No_Prof; $n++) {
  @vals = &trim($n+1, $No_Seg[$n], $Prof_Type[$n], $Dup_flag[$n], $Digit_Code[$n], 
    $Standard[$n], $Deep_Depth[$n]);
    print join(",", @vals), "\n";
  }
  #
  print "Sfc_Parm_Gp#,Pcode,Parm,Q_Parm\n";
  for (my $n = 0; $n < $Nparms; $n++) {
    @vals = &trim($n+1, $Pcode[$n], $Parm[$n], $Q_Parm[$n]);
    print join(",", @vals), "\n";
  }
  #
  print "Sfc_Code_Gp#,SRFC_Code,SRFC_Parm,SRFC_Q_Parm\n";
  for (my $n = 0; $n < $Nsurfc; $n++) {
    @vals = &trim($n+1, $SRFC_Code[$n], $SRFC_Parm[$n], $SRFC_Q_Parm[$n]);
    print join(",", @vals), "\n";
  }
  #
  print "Hist_Gp#,Ident_Code,PRC_Code,Version,PRC_Date,Act_Code,Act_Parm,Aux_ID,Prev_Val\n";
  for (my $n = 0; $n < $Num_Hists; $n++) {
    @vals = &trim($n+1, $Ident_Code[$n], $PRC_Code[$n], $Version[$n], $PRC_Date[$n],
      $Act_Code[$n], $Act_Parm[$n], $Aux_ID[$n], $Previous_Val[$n]);
    print join(",", @vals), "\n";
  }
  #
  my $tsegs = 0;
  for (my $rec = 0; $rec < $No_Prof; $rec++) {
    for (my $seg = 0; $seg < $No_Seg[$rec]; $seg++) {
      print "Prof_Rec,Seg,MKey_r,One_Deg_sq_r,Cruise_ID_r,Obs_Year_r,Obs_Month_r,",
        "Obs_Day_r,Obs_Time_r,Data_Type_r,Iumsgno_r,Profile_Type_r,Profile_Seg,",
        "No_Depths,D_P_Code_r\n";
      @vals = &trim($rec+1, $seg+1,$MKey_r[$tsegs],$One_Deg_sq[$tsegs],
        $Cruise_ID[$tsegs],$Obs_Year[$tsegs],$Obs_Month[$tsegs],$Obs_Day[$tsegs],
        $Obs_Time[$tsegs],$Data_Type[$tsegs],$Iumsgno[$tsegs],$Profile_Type[$tsegs],
        $Profile_Seg[$tsegs],$No_Depths[$tsegs],$D_P_Code[$tsegs]);
      print join(",", @vals), "\n";
      #
      print "Level,D/Press,DPq,ProfParm,PPq\n";
      for (my $d = 0; $d < $No_Depths[$tsegs]; $d++) {
        @vals = &trim($d+1, $Depth_Press[$tsegs][$d], $Depres_Q[$tsegs][$d],
          $Prof_Parm[$tsegs][$d], $Prof_Q_Parm[$tsegs][$d]);
        print join(",", @vals), "\n";
      }
      $tsegs++;
    }
  }
}

#####
sub hextobin {
  my $hex = $_[0];
  $hex =~ s/^ +//;
  $hex =~ s/ +$//;
  my $bin32 = sprintf("%032b", oct("0x".$hex));
  return $bin32;
}

######
sub trim {
  my @out = @_;
  for (@out) {
    s/^\s+//; # trim left
    s/\s+$//; # trim right
  }
  return @out == 1
    ? $out[0] # Only one to return
    : @out;   # or many
}


