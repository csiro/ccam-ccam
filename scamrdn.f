      subroutine scamrdn
c reads the parameters for scam from scamfil into common block scamc

      include 'newmpar.h'
      include 'arrays.h'
      include 'filnames.h'  ! list of files, read in once only
      include 'nsibd.h'     ! ivegt,isoilm
      include 'parm.h'      ! id,jd,idjd
      include 'pbl.h'       ! tss
      include 'scamdim.h'   ! dimension of patches
      include 'scampar.h'   ! dimension of input,store,output
      include 'sigs.h'
      include 'soil.h'      ! land, zolnd passed to indata
      include 'soilsnow.h'  ! initial soil temperatures and moistures
      include 'soilv.h'

* --- physical
      parameter (r=287., g=9.806)
c     dimension ssat(9),swilt(9)                               ! sflux,rdnsib
c     data ssat/.398,.479,.482,.443,.426,.482,.420,.45,.479/   ! sflux,rdnsib
c     data swilt/.072,.216,.286,.135,.219,.283,.175,.395,.216/ ! sflux,rdnsib

      real z040(0:44)
      data z040/ 0.0001,
     &   1.8, 1.6, 1.41, 1.415,1.599,0.515,0.573,0.666,0.649,0.293,
     & 0.270,0.078,0.126,0.061,0.012,0.052,0.012,0.061,0.079,0.049,
     & 0.009,0.049,0.049,0.032,0.033,0.039,0.005,0.230,0.0001,0.0001,
     & 0.0001,
     & 1.6,1.1,0.83,0.67,0.57,0.32,0.22,0.14,0.01,0.18,0.001,0.06,0.01/
c deans veg types 0-31 (0 is water, etc)
c global veg types 32-43
c stand alone (for scam) veg type 44
c-----------------------------------------------------------------------
c input parameter file for scam surface scheme
      print*,'in scamrdn'
      open(unit=8,file=scamfile,status='old') 

      call comskp(8)
      do i=1,npara0                    ! gridpoint-indept input params
        read(8,110) param0(i),pname0(i)
        write(*,*) i,param0(i),pname0(i)
      end do
c      write(*,*)(param0(i),i=1,npara0)
      call comskp(8)
      do i=1,npara1                    ! gridpoint-dept input params
        read(8,110) param1(i),pname1(i)
        write(*,*) i,param1(i),pname1(i)
      end do
c      write(*,*)(param1(i),i=1,npara1)
      call comskp(8)
      do i=1,npara2                   ! gridpoint-dept input params
        read(8,110) param2(1,i),pname2(i)
        write(*,*) i,param2(1,i),pname2(i)
      end do
c      write(*,*)(param2(1,i),i=1,npara2)
      call comskp(8)
      do i=1,nstore                    ! initial storages
        read(8,110) store(1,i),sname(i)
        write(*,*) i,store(1,i),sname(i)
      end do
c      write(*,*)(store(1,i),i=1,nstore)
110   format(f10.3,a8)

      param0(2) = dt
      param0(1) = -r*t(il/2,jl/2,1)*log(sig(1))/g

      write(*,111) (pname0(i),param0(i),i=1,npara0)
111   format(/,4(1x,a8,1x,f8.2))

      print*,'initalisation of tgg,wb values'
      do i=1,ifull
         tgg(i,1) = tss(i,1) ! inital temperature at second layer(6.5.97 KF)
cc       wb(i,1)  = w(i,1)   ! top soil moisture
         wbice(i,1) = 0.0    ! top soil ice
         smass(i,1) = 0.0
         smass(i,2) = 0.0
         smass(i,3) = 0.0
         ssdn(i,1)  = 120.
         ssdn(i,2)  = 120.
         ssdn(i,3)  = 120.
         do k=2,ms
c           tgg(i,k) = t(i,1,2)   ! inital temperature at second layer(6.5.97 KF)
c           wb(i,k)  = w2(i,1)    ! layer initialisation of moisture
            wbice(i,k)  = 0.0
         enddo
c        tgg(i,k) = ts(i,1,1  !inital temperature from GCM runs with 3 layers 
c        tgg(i,k)   = ts(i,1,2) !inital temperature from GCM runs with 3 layers 
c       if( ntsur2.gt.2) then
c         tgg(i,2) = ts(i,1,2) !inital temperature from GCM runs with 3 layers
c         tgg(i,3) = ts(i,1,2) !inital temperature from GCM runs with 3 layers
c         do kk=4,ms
c          tgg(i,kk) = ts(i,1,1)
c         enddo
c        endif
         store(i,1) = 0.0      ! relative canopy water store initalisation (on leaves)
      enddo
        print*,'in scamrdn smass,ssdn',idjd,
     &  smass(idjd,1),smass(idjd,2),smass(idjd,3),
     &  ssdn(idjd,1),ssdn(idjd,2),ssdn(idjd,3)
 
      print *,'correct w and w2 here from Evas to Klaras'
c     *******    probably needs kstart or nrungcm switch
      do iq=1,ij
       istype = isoilm(iq)      ! Deans and global set
       if (istype.gt.0)then
cc       win=w(iq)
cc       w2in=w2(iq)
cc       w (iq)   = min( max(w (iq),swilt(istype)) ,ssat(istype))
cc       w2(iq)   = min( max(w2(iq),swilt(istype)) ,ssat(istype))
cc       if(iq.eq.idjd)print *,'iq,istype,swilt,ssat,win,w2in,
cc   .    w,w2 ',
cc   .    iq,istype,swilt(istype),ssat(istype),win,w2in,
cc   .    w(iq),w2(iq)
       endif! (istype.gt.0)then
      enddo

c     first test with grid averages for soil type and veg type
c     param2 dimensioned with il*jl
      idjd = id+il*(jd-1) 
      do iq=1,ifull
c        if(land(iq)) then   ! non-Dean's stuff
c          ivegt(iq)  = ivegt(iq) + 31  ! offset for global Veg type
c        endif
         param2(iq) = ivegt(iq)
         param2(iq,2) = isoilm(iq)
         param2(iq,3) = 1.              ! zsm
         param2(iq,4) = 1.              ! fraction of gridpoint
c        if (ivegt(i,j).gt.44) then
c           print*, 'ERROR in scamrdn: ivegt>44 at ',i,j
c           stop
c        endif
         zolnd(iq) = z040(ivegt(iq))
         if (iq.eq.idjd)then
           print *,'iq,land,ivegt,zolnd ',
     .              iq,land(iq),ivegt(iq),zolnd(iq)
           print *,'initial at id,jd', (store(iq,k),k=1,1)
         endif
c        initialistion of stability parameter (neutral stability)
c        zetai(iq) = 0.
      enddo
c check that global parameters are within possible ranges
c     print*,'check param0,1 input values'
c     isinf  = param0(5)   ! soil infiltration: 1,2 = MP84, FC96
c     isevap = param0(6)   ! soil evap: 1,2,3 = alfa,beta,threshold
c     niter  = param0(11)  ! number of iterations for za/L
c     zdlin  = param0(16)  ! height frac of d below which TL linear
c     fdrain = param0(20)  ! exfilt param: qexf=fdrain*cond2
c     print*,'ismois,isinf,isevap,zdlin,fdrain'
c     print*,ismois,isinf,isevap,zdlin,fdrain
      if (isinf .gt.2.or.isinf.lt.1)      stop 'illegal isinf'
      if (isevap.gt.2.or.isevap.lt.1)     stop 'illegal isevap'
      if (niter.gt.10)                    stop 'illegal niter'
      if (zdlin.gt.1.0.or.zdlin.lt.0.0)   stop 'illegal zdlin'
      if (fdrain.gt.1.0.or.fdrain.lt.0.0) stop 'illegal fdrain'
c this check has to be done with the patches !!!! so that the
c number of patches in grid cell i,j does not exceed npmax
ckf      if (npij.gt.npmax)                  stop 'illegal npij'

c Set names for parameters in array hold
      hname(1)  = '"rlai"  '
      hname(2)  = '"hc  "  '
      hname(3)  = '"disp" '
      hname(4)  = '"usuh"   '
      hname(5)  = '"coexp"  '
      hname(6)  = '"zruffs"  '
      hname(7)  = '"trans"  '
      hname(8)  = '"rt0us" '
      hname(9)  = '"rt1usa"'
      hname(10) = '"rt1usb"'

c       call readreal(albfile,alb,il*jl)
c       call readreal(rsmfile,rsmin,il*jl)
c       call readreal(zofile,zolnd,il*jl)
c       call readint(vegfile,ivegt,il*jl)
c
c index,npij and param2 arrays need to be assigned
c

      close(8)

c Open file
      print*,'open scam files'
      open(unit=2, file='scama.prn',status='unknown')   ! detailed time series file 1
      open(unit=3, file='scamb.prn',status='unknown')   ! detailed time series file 2
      open(unit=4, file='scamz.prn',status='unknown')   ! detailed time series file 3
      open(unit=7, file='scamo.prn',status='unknown')   ! input from darlam saved for stand alone

c write input data to detailed output file
      write(2,121) text,metfil
121   format(' "',a40,'"',/,' "metfile:"',2x,'"',a40,'"')
      write(2,122) (param0(j),pname0(j),j=1,npara0)
      write(2,123) (param1(j),pname1(j),j=1,npara1)
122   format(' "GLOBAL input parameters (array param0):"',/,
     &          (5(e11.4,1x,a8,1x)))
123   format(' "LOCAL input parameters (array param1):"',/,
     &          (5(e11.4,1x,a8,1x)))

      write(2,124) (store(idjd,j),sname(j),j=1,nstore)
124   format(' "initial STORAGE contents at i,j:"',/,
     &      (5(f10.5,1x,a8,1x)))

      return
99    print*,'error in opening scam input file',scamfile

      stop
      end
c#######################################################################

      SUBROUTINE COMSKP(IUNIT)
C MRR, 5-AUG-83
C SKIPS COMMENT LINES IN CONTROL DATA FILE, STOPS IF EOF FOUND
      CHARACTER*1 COM
1     READ(IUNIT,100,END=2) COM
100   FORMAT(A1)
      IF(COM.EQ.'C'.or.com.eq.'c') GOTO 1
      BACKSPACE IUNIT
      RETURN
2     STOP 'CONTROL FILE EOF'
      END

c#######################################################################

