       program mkstns
       parameter(ndrp=350)
!
! compile with 
!     gfortran mkstns.f -o mkstns.x
!
!      (28mar2018 used to use s files, change to e files)
! create missing stations.dat from e files.
! this version just assumes all are good, iedt=1. You could search
! for NG's , another time...
       character statfilen*12, ayr*4, efilen*12
       character drpn*3, lath*1, lonh*1, depint*2, time*4
!                   123456789012
       data efilen/'p099709e.001'/
       data statfilen/'stations.dat'/
!
       write(*,*)'Enter cruisename, ie p311802 (7 chars)'
       read(5,'(a)') efilen(1:7)
       open(11,file=statfilen,status='unknown',form='formatted')
       do 100 i = 1, ndrp
          if(i.le.9) write(efilen(12:12),'(i1)') i
          if(i.ge.10.and.i.le.99) write(efilen(11:12),'(i2)') i
          if(i.ge.100) write(efilen(10:12),'(i3)') i
          iedt = 1
c read e file
          open(10,file=efilen,status='old',form='formatted',
     $         err=101)
          read(10,500,end=101) drpn, xlat, xlon
500       format(2x,a3,15x,2f8.3)
          read(10,501,end=101)  idy, imo, ayr, ihr, imin, isec
501       format(1x,i2,1x,i2,1x,a4,2x,i2,1x,i2,1x,i2)
c read in temps:
          do 50 j = 1, 1000
             read(10,502,end=51) idep, itemp, iclass
! if j=1 (1st depth) is iclass=4, it's no good, put "-1" in stations.dat
! edit col, do you want iclass=3 here too??
             if(j.eq.1) then
               if(iclass.eq.4) then
                iedt = -1
               else 
                iedt = 1
               endif
             endif
             if(idep.ge.697) go to 51
!502          format(i3,1x,i5)
502          format(i3,1x,i5,1x,i1)
50        continue
51        continue
          close(10)
          write(11,505) drpn, 0, real(itemp*.001), idy, imo, ayr(3:4),
     $            ihr,imin,isec, xlat, xlon, 0, iedt, 0
505          format(1x,a3,i6,f7.3,1x,i2,1h/,i2,1h/,a2,1x,i2,1h:,i2,1h:,
     $i2,f9.3,f9.3,i4,i5,i6)
100       continue
101       continue
       write(11,'(a7)') 'ENDDATA'
       close(11)
       write(*,*)
     $'Check for missed NG and shallow e/q files to modify iedt (col 9)'
       write(*,*)
       write(*,*)'If this is CSIRO data, use wpsxbt to view'
       
       stop
       end
