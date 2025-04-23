! TODO: when you change the surface CS back to 1 from 5, you're not
! checking if it should actually be changed to a 2... TODO

! 13may2021 LL - erl script incorrectly converting longs to negative -  fix
! 04apr2019 LL - set drop number (filename) to loop counter.
! 19feb2019 LL - Deal with NaN from input file. How about if you find a
!                NaN, set temp=-999.0 and qc=4 ?
! 01feb2019 LL - Change your skip of "NG" to look if NG set at dep=0.67
!                (index=1), if so, then mark all NG
! 22jun2018 LL - add read of "MFD#" - looks like probe DoM
!              - add create control.dat if does not exist (just cruise
!              name)
! nov2017:
! read output of ma2txt.pl (modified gtspp2txt.pl (read meds-ascii *.MA))
! to output q and e format files.
       program matxt2qe
       parameter(ndepmax=3000)
       parameter(nchangemax=200)
!
       character*22 infile 
       character*7 cruisename 
       character*12 qfile , efile
       character*4 asurf 
       character*10 aasurf
       character*10 acall,acsid,apeq,arct,aser,adom
       character*2 achange(nchangemax)
       character*3 aend 
       character*2 aqc(ndepmax)                        ! qc code, ie HB,LE,
       character*5 adep
       character*85 hdr
!
       real d(ndepmax), t(ndepmax) , torig(ndepmax)
       integer iqc(ndepmax)
       real dc(nchangemax), tc(nchangemax)
       integer stat
!
       data infile/'JM3403hed2018.txt     '/
       data qfile /'p341803q.000'/ 
       data efile /'p341803e.000'/ 
!                   12345678901234
!
       infile = '                    '
       adom =   '          '
       write(*,*)'Enter input filename: (ie JM3403hed2018.txt)'
       read(5,'(a)') infile
       write(*,*)' Enter cruisename: (ie, p311802) (7 chars!)'
       read(5,'(a7)') qfile(1:7)
       efile(1:7) = qfile(1:7)
       write(*,*)'reading: ', infile
       write(*,*)'writing: ', qfile(1:7), efile(1:7)
!
! 22jun2018 if no control.dat exists, create basic cruise name one
       open(40,file='control.dat',status='old',err=5)
       close(40)
       go to 6      ! if we are here, control.dat exists
5      open(40,file='control.dat',status='unknown',form='formatted')
       write(40,'(a7)') qfile(1:7)
       close(40)
6      continue
!
       open(10,file=infile,status='old',form='formatted')
       open(33,file='error.out',status='unknown',form='formatted')
!
       do 200 k = 1, 500
! 04apr2019 set idrop (file number) to k, no inputdrop...
       idrop = k
! read basics:
       read(10,'(i3)',end=201,err=201) inputdrop
       write(33,*) 'inputdrop=', inputdrop
       read(10,'(i4)') iyear
       read(10,'(i2)') imo
       read(10,'(i2)') idy
       read(10,'(i2,i2)') ihr,imn
       read(10,'(f8.2)') xlat
       read(10,'(f9.2)') xlon
! the perl script is incorrectly reading longs as negative, fix:
       if(qfile(1:3).eq.'p09'.and.xlon.lt.0.0) xlon=-xlon  ! Szechuen 5/2021
       if(qfile(1:3).eq.'p31'.and.xlon.lt.0.0) xlon=-xlon
       if(qfile(1:3).eq.'p34'.and.xlon.lt.0.0) xlon=-xlon
       if(qfile(2:3).eq.'28'.and.xlon.lt.0.0) xlon=-xlon
       if(xlon.lt.0.0) xlon = 360.0 + xlon
       write(33,*) 'inputdrop=',inputdrop,year,imo,idy,ihr,imn,xlat,xlon
!
! read surface codes:
       read(10,'(i3)') isurf
       do 10, i = 1, isurf
          read(10,'(a4,4x,a10)') asurf, aasurf
          if(asurf(1:4).eq.'GCLL') acall(1:10)=aasurf(1:10)
          if(asurf(1:4).eq.'CSID') acsid(1:10)=aasurf(1:10)
          if(asurf(1:4).eq.'PEQ$') apeq(1:10) =aasurf(1:10)
          if(asurf(1:4).eq.'RCT$') arct(1:10) =aasurf(1:10)
          if(asurf(1:4).eq.'SER#') aser(1:10) =aasurf(1:10)
          if(asurf(1:4).eq.'MFD#') adom(1:10) =aasurf(1:10)
10     continue
       write(33,*) acall,acsid,apeq,arct,aser
!
! read station data d,t
       read(10,'(a)') adep
!       print *, 'adep=',adep
       read(adep,'(i5)') ndep
!       print *, 'ndep=',ndep
       write(33,*) 'ndep=', ndep
       if(ndep.gt.ndepmax) stop ' ndep>ndepmax '
       ifoundnan = 0             ! init to 0 - no NaN's found
       do 20, i = 1, ndep
          read(10,'(f6.2,f11.3,i2)') d(i),t(i),iqc(i)
          aqc(i) = '  '
! 19feb2019 test if t(i) is NaN:
          if(isnan(t(i))) then
            if(ifoundnan.eq.1) then
               aqc(i) = '  '        ! only set first NaN
            else
               aqc(i) = 'SP'        ! make SP up since I dunno why NaN
            endif
            t(i) = 0.0
            iqc(i) = 4
            ifoundnan = 1        ! set to 1 once found a NaN
          endif
!          write(33,'(f6.2,f11.3,i2)') d(i),t(i),iqc(i)
          torig(i) = t(i)    ! ??
          if(ifoundnan.eq.1) iqc(i) = 4  ! be sure all below NaN are qc=4
20     continue
!
! read station data d,t changes:
       read(10,'(i3)') nchange
       if(nchange.gt.nchangemax) stop ' nchange>nchangemax '
       if(nchange.eq.0) go to 31
       do 30, i = 1, nchange
          read(10,'(a2,f11.2,f11.3)') achange(i),dc(i),tc(i)
30     continue
31     continue
!
       read(10,'(a3)') aend
       if(aend.ne.'END') stop ' does not = END '
! now put changed data and qc codes back into data:
       do 40 i = 1, nchange
!
          if(achange(i).eq.'QC') go to 40   !skip
          if(achange(i).eq.'RE') go to 40   !skip
!
! 01feb2019          if(achange(i).eq.'NG') go to 40   !skip
!  Only set whole prog=NG if index=1, otherwise skip
          if(achange(1).eq.'NG') then
            iqc(1) = 4
            aqc(1) = 'NG'
            go to 40
          elseif(achange(i).eq.'NG') then
            go to 40 !skip
          endif
          if(achange(i).eq.'CS'.and.tc(i).ne.99.990) then       !put surface data back...
             do 35 j = 1, ndep
                if(dc(i).eq.d(j)) then      ! found matching depth
                   t(j) = tc(i)             ! replace t
! 8feb2022: BUT you should check if iqc should be 2 here!!!! TO DO.....
                   iqc(j) = 1               ! reset qc to 1 (from 5)
                   go to 40
                endif
35           continue
          endif
!
          if(achange(i).eq.'HB'.or.achange(i).eq.'WB') then
             do 36 j = 1, ndep
                if(dc(i).eq.d(j)) then      ! found matching depth
                   aqc(j) = achange(i)      ! set qc char code
                   go to 40
                endif
36           continue
          endif
!
! if we are here, it's another code:
          do 39 j = 1, ndep
             if(dc(i).eq.d(j)) then         ! found matching depth
! they have changed tems, yet still keep qc code=2! (ipa, interp areas)
! what is iqc of this d(j) ?:
                if(iqc(j).eq.5) then    ! 5 is changed data
                   aqc(j) = achange(i)      ! set qc char code
                   torig(j) = tc(i)
! this is for they qc=2, yet the tem has been changed:
                elseif((iqc(j).eq.2).and.(t(j).ne.tc(i))) then
                   aqc(j) = achange(i)
                   torig(j) = tc(i)
                   iqc(j) = 5
                else
                   aqc(j) = achange(i)      ! set qc char code
                endif
                go to 40
             endif
39        continue
! if we are here, I didn't deal with the changed code:
          if(achange(i).eq.'WB') then
             write(33,*)'WB deeper than data:',i,achange(i),dc(i),tc(i)
             write(33,*)'-->>>>> ignore it'
          else
           write(33,*) 'what to do with ', i, achange(i), dc(i),tc(i)
           stop
          endif

40     continue       
!
! 19apr2018 LL if aser(1:4)=Test, set iqc(1)=4 and aqc(1)='TP'
       if(aser(1:4).eq.'Test') then
          write(33,*)'change aqc(1) to TP'
          iqc(1) = 4
          aqc(1) = 'TP' 
! and since this is q file, do iqc(2:ndep) too:
          do 44 ii = 2, ndep
             iqc(ii) = 4
44        continue

       endif
! 19apr2018 LL if iqc=3, then change it to iqc(1)=4 & aqc(1)='NG' &
!      carry original aqc(1) code to aqc(2), this is because csiro
!      just labels bad drops as iqc(1)=3, and I want them NG'd
! since q file, fiddle with which depths to change so carries in e file:
       if(iqc(1).eq.3) then
          aqc(4) = aqc(1)
          iqc(4) = iqc(1)
          iqc(1) = 4
          iqc(2) = 4
          iqc(3) = 4
          aqc(1) = 'NG'
          write(33,*)'change aqc(1) to NG'
       endif
!
       write(33,*)'finish nchange'
!
! and here let's try writing a 'q' file:
       qfile(10:12) = '000'                    !data qfile /'p340811q.000'/ 
       efile(10:12) = '000'                    !data efile /'p340811e.000'/ 
       if(idrop.le.9) then
          write(qfile(12:12),'(i1)') idrop
       elseif(idrop.ge.10.and.idrop.le.99) then
          write(qfile(11:12),'(i2)') idrop
       else
          write(qfile(10:12),'(i3)') idrop
       endif
       print *, 'qfile=',qfile
       open(20,file=qfile,status='unknown',form='formatted')
!
! write q header:
       hdr(1:85) = ' '
       write(hdr(1:36),500)qfile(10:12),apeq(1:3),arct(1:2),
     $              acall(1:7),xlat,xlon
500    format(2x,a3,1x,a3,1x,a2,1x,a7,2f8.3)
       write(hdr(37:85),501) idy,imo,iyear,ihr,imn,acsid(1:8),aser(1:7),
     $                       adom(5:6),adom(7:8),adom(1:4)
501    format(1x,i2,'-',i2,'-',i4,2x,i2,':',i2,':00 ',a8,1x,a7,
     $        1x,a2,'-',a2,'-',a4)
        write(20,'(a36)',err=201) hdr(1:36)
        write(20,'(a49)',err=201) hdr(37:85)
!
! write q data, double check if have any extra iqc=5 without the torig:
       do 100 j = 1, ndep
          if(aqc(j).eq.'CS') aqc(j) = '  '    ! take care of extra CS
          if(iqc(j).ne.5) then
             write(20,502) nint(1000.0*d(j)),nint(1000.0*t(j)),
     $                     iqc(j),aqc(j)
502          format(i7,1x,i5,1x,i1,1x,a2)
          else
             write(20,503) nint(1000.0*d(j)),nint(1000.0*t(j)),iqc(j),
     $                        aqc(j),nint(1000.0*torig(j))
503          format(i7,1x,i5,1x,i1,1x,a2,1x,i5)
          endif
100    continue
       close(20)
       call Convertq2e(qfile,d,t,ndep,hdr,iqc,aqc,torig,12)
       close(21)
200    continue
201    continue
       stop
       end
!
!-----------------------------------------------
! pulled Convertq2e from srpedit.for:  SOME MODS FOR THIS PROGRAM!!!
!*****************************************************
        SUBROUTINE Convertq2e(f1,d,t,num,hdr,iclass,code,torig,
     $                         len_f1)
! Try to convert q file dep,tem,iclass,code,torig full resolution to
!                e file binned resolution
! Have to keep any surface codes at the surface, and then as work
! the way down in depth, use the deepest class& code of the bin ?
        parameter(ndepmax=3000)
        character*(*) f1, hdr*85
        real*4 d(ndepmax), t(ndepmax),torig(ndepmax)
        integer*4 iclass(ndepmax), iclassbin
        character*2 code(ndepmax), codebin
        character*12 f2

! create matching e filename:
!
!LL MOD
! p340811e.000
! 123456789012
        f2(1:12) = f1(1:12)
        f2(8:8) = 'e'
!        f2(1:len_f1+5)=f1(1:len_f1+5)
!        f2(len_f1+6:len_f1+6)=char(0)
!        f2(len_f1+1:len_f1+1)='e'

        write(33,*)'in Convertq2e', f2

        open(30, file=f2, form='formatted',status='unknown', err=200)
!!!!        write(30,'(a35)',err=201) hdr(1:35)
!!!!        write(30,'(a39)',err=201) hdr(36:74)
        write(30,'(a36)',err=201) hdr(1:36)
        write(30,'(a49)',err=201) hdr(37:85)

        ibins = 0
        nbin = 0
        avt = 0.0
        avtorig = 0.0
        iclassbin = 0
        codebin = '  '

        do 310 j=1,num
        temp=t(j)
        temporig=torig(j)
!!!!        depth=-d(j)
        depth=d(j)
        ibin=ifix(depth/2.)+1
!        write(33,*) 'temp,depth,ibin',temp, depth, ibin
        if(ibin.eq.ibins.or.ibins.eq.0)then
                nbin=nbin+1
                avt=avt+temp
                avtorig=avtorig+temporig
                if(iclassbin.ne.iclass(j)) iclassbin = iclass(j)
                if(codebin.eq.'  ') then   ! only overwrite upper code if it's blank
                   codebin = code(j)
                endif
        else
                avt=avt/float(nbin)
                avtorig=avtorig/float(nbin)
                iavt=nint(1000.*avt)
                iavtorig=nint(1000.*avtorig)
                idep=2*ibin-3
! 08jun2015 for idep>999 change idep=idep-1000:
                if(idep.gt.999) idep = idep-1000
                if(ibin.eq.1) then               ! if surface bin, write iclass&code of t(1)
                   write(30,500)idep,iavt,iclass(1),code(1)
                elseif(ibin.gt.1.and.iclassbin.ne.5) then
                   write(30,500)idep,iavt,iclassbin,codebin
                elseif(ibin.gt.1.and.iclassbin.eq.5) then
                   write(30,501)idep,iavt,iclassbin,codebin,iavtorig
                endif
500             format(i3,i6,1x,i1,1x,a2)
501             format(i3,i6,1x,i1,1x,a2,1x,i5)
                avt=temp
                avtorig=temporig
                nbin=1
                iclassbin = iclass(j)
                codebin = code(j)
        endif
                ibins=ibin
310     continue
        close(30)
200     continue
201     continue
        return
        END SUBROUTINE Convertq2e
!*****************************************************



