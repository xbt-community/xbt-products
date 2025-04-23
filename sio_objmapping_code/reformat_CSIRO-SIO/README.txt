3 routines convert the CSIRO cruise data file to the SIO file formats for mapping: ma2txt.pl perl script converts medsascii to text, matxt2qe.x creates one e and one q file for each profile (SIO's txt versions of the qc'ed 2m averaged data and full resolution data), mkstns.x creates our stations.dat (or pXXMMYY.dat) file that summarizes the cruise drops, one line for each profile. All these routines rely on the data being named pXXYYMM pXX=transect name eg p05, YYMM are the year and month the cruise started, a,b, c are appended if there is more than one cruise on the line that month; the routines are also setup to look for data in directory structure like pXX/YYMMProcess CSIRO data:

Notes from Lisa:

Download from Bec email <name>.MA files (meds-ascii), this example using
filename:
D5LS3_SB30022023.MA

Bec will tell you what line, or you can see here:
D5LS3 = callsign
SB = ship initials, SB=Seatrade Blue
30 = px30 (we call it p31)
02 = uhhh, 2nd cruise of 2023, or is it 2nd cruise of Seatrade Blue...?
2023 = year
.MA meds-ascii format

create a dir under p31, guess it's 2301 here.  (cd /data/xbt ; mkdir 2301)
 cp D5LS3_SB30022023.MA /data/xbt/p31/2301/.

Note if it's p28 you'd do 2301a - 1st SB transect of month 01
Note if it's p28 you'd do 2301b - 1st NB transect of month 01
Note if it's p28 you'd do 2301c - 2nd SB transect of month 01
Note if it's p28 you'd do 2301d - 2nd NB transect of month 01
Note CSIRO PX30 = SIO p31
Note CSIRO PX32 = SIO p34
Note CSIRO IX28 = SIO p28

*) run ma2txt.pl   
Original perl script written by Norm Hall, modified/renamed by Lisa

ma = medsascii
2 "to"
txt = text
.pl perl

-i input-filename
-o output-filename
-unw unwind

ma2txt.pl -i D5LS3_SB30022023.MA -o D5LS3_SB30022023.txt -unw

Outputs a human readable <file>.txt file, tells you number of profiles

'more' the file to see date of beginning of cruise. If you need to
move to a different named folder, do it now. Let's say you guessed 
it was 2302 and you put it there, then find out it should be 2301 then:
   cd /data/xbt/p31
   mv 2302 2301
You can also see the shipname/callsign,
Make sure it's in file /data/xbt/callsign.txt
		-> watch the format here, mapxbt3.x reads it, NO TABS!

*) Run matxt2qe.x
clever, eh? 
matxt (output of ma2txt)
2 to
qe - outputting q and e files
.x - Lisa names exe files .x

Input filename: D5LS3_SB30022023.txt
enter your cruise name: p312301

outputs: in this case, 24 e and 24 q files.

*) mkstns.x
create a stations.dat file to match above e/q files

*) wpsxbt	/data/xbt/prog/
to view qc
If any are less than, say 50m, consider marking as NG - to do that
edit stations.dat and put a '-1' in the edit column 9.
ALSO - check the NG's and near surface failures are -1 in edit column
Make a note in your cruise README 

If any are PL premature launch - I've been marking those bad too since I've not written
a program to deal with CSIRO's PL's... Generally they have a redrop.

*) cp stations.dat p312301.dat    (use your cruise name here, delete last line "ENDDATA")

*) run tenm3_qc.x
Note! this is DIFFERENT than tenm3_chgcoef.x!
CSIRO already has new coefficients, so do not change them again.
*) day.x to get julian day for xbtinfo file
*) add to xbtinfo.p?? 

*) continue normal cruise processing found in process-general - 
steps AFTER tenm3_chgcoef.x...
*) NOTE: for p28 - don't process past temperature grid. Not doing
salinity/etc there.
