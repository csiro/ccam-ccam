module cable_ccam
  ! rml 02/10/07 removed sib4 and other cbm/cable related subroutines from sflux
  ! rml 23/10/07 added cbm_pack and cbm_unpack to this module
  !
  ! This file contains the following subroutines:
  !   sib4,
  !   setco2for,
  !   cbmrdn,
  !   comskp,
  !   rlaiday,
  !   cbm_pack,
  !   cbm_unpack, and
  !   eva_output
  !
  include 'newmpar.h'
  integer, parameter :: CABLE = 4
  ! will need to have vertical dimension for interaction with radiation?
  real, dimension(ifull) :: atmco2

  character(len=80), save :: vegtypefile=''
  character(len=80), save :: vegparmfile=''
  logical, save :: vegparmnew=.true. ! old or new format for veg parm file
  integer, save :: CO2forcingtype=1  ! 1 constant, 2 prescribed 1900-2004,
                                     ! 3 interactive
  integer, save :: initcarbpools=999  ! 1 initialise from veg_parm file values
  namelist/cableinput/vegtypefile,vegparmfile,vegparmnew, &
                    & CO2forcingtype,initcarbpools

  ! arrays for vegetation names and types
  character(len=25), dimension(:), allocatable, save :: vegname
  character(len=10), dimension(:), allocatable, save :: vegtype
  logical, dimension(:), allocatable, save :: forest
  logical, dimension(:), allocatable, save :: deciduous
  logical, dimension(:), allocatable, save :: shrub
  logical, dimension(:), allocatable, save :: grass
  logical, dimension(:), allocatable, save :: crop
  logical, dimension(:), allocatable, save :: noveg  ! ice, possibly also barren
       

  contains
  ! ****************************************************************************

  subroutine sib4(nvegt)     ! new version of sib1 with soilsnowv
! BP added nvegt to be passed down to cbm (BP jan 2008)
!  subroutine sib4     ! new version of sib1 with soilsnowv

      use cc_mpi
      use zenith_m
      USE define_types, cbm_ms => ms
      USE air_module
      USE canopy_module
      USE carbon_module
      USE cbm_module
      USE soil_snow_module
      INTEGER, PARAMETER          :: mp_max = 99 ! max. no. geographic points
      TYPE (air_type)             :: air
      TYPE (bgc_pool_type)        :: bgc
      TYPE (canopy_type)          :: canopy
      TYPE (met_type)          :: met
      TYPE (balances_type)     :: bal
      TYPE (radiation_type)       :: rad
      TYPE (roughness_type)       :: rough
      type (soil_parameter_type)      :: soil       ! soil parameters
      TYPE (soil_snow_type)       :: ssoil
      TYPE (sum_flux_type)        :: sum_flux
      type (veg_parameter_type)       :: veg        ! vegetation parameters
      ! Save these so only have to do the allocation once.
      save air, bgc, canopy, met, bal, rad, rough, soil, ssoil, &
     &     sum_flux, veg


      include 'newmpar.h'
      include 'aalat.h'    ! slwa
      include 'arrays.h'
      include 'carbpools.h'
      include 'dates.h' ! ktime,kdate,timer,timeg,xg,yg
      include 'extraout.h'
      include 'filnames.h'
      include 'latlong.h'  ! rlatt,rlongg
      include 'map.h'      ! id,jd,idjd
      include 'morepbl.h'
      include 'nsibd.h'    ! rsmin,sigmf,tgf,ssdn,res,rmc,tsigmf
      include 'parm.h'
      include 'permsurf.h'
      include 'pbl.h'
      include 'prec.h'
      include 'screen.h'   ! tscrn etc
      include 'sigs.h'
      include 'soil.h'     ! ... zmin zolnd zolog sice alb
      include 'soilsnow.h' ! 
      include 'soilsnin.h'
      include 'soilv.h'    ! ssat, clay,..
      include 'vegpar.h' ! 
      include 'trcom2.h'   ! nstn,slat,slon,istn,jstn
!                     met forcing for CBM
      common/work2/dirad(ifull),dfgdt(ifull),degdt(ifull)  &
     & ,wetfac(ifull),degdw(ifull),cie(ifull)              &
     & ,factch(ifull),qsttg(ifull),rho(ifull),zo(ifull)    &
     & ,aft(ifull),fh(ifull),ri(ifull),theta(ifull)        &
     & ,gamm(ifull),rg(ifull),vmod(ifull),dummwk2(ifull)

      ! for calculation of zenith angle
      real fjd, r1, dlt, slag, dhr, coszro2(ifull),taudar2(ifull)

      integer imonth(12),iwrk(ifull)
      data imonth /31,28,31,30,31,30,31,31,30,31,30,31/
      dimension ndoy(12)   ! days from beginning of year (1st Jan is 0)
      data ndoy/ 0,31,59,90,120,151,181,212,243,273,304,334/
      real ssumcbmfl(14,20)
      logical, save :: cbm_allocated = .false.
      save ktauplus,iswt
      data iswt/0/
      integer nvegt

       if ( .not. cbm_allocated ) then
          ! These variables are all saved so only need to be allocated once
          ! in the run.
          cbm_allocated = .true.
          mp = ipland
          call alloc_cbm_var(air, mp)
          call alloc_cbm_var(bgc, mp)
          call alloc_cbm_var(canopy, mp)
          call alloc_cbm_var(met, mp)
          call alloc_cbm_var(bal, mp)
          call alloc_cbm_var(rad, mp)
          call alloc_cbm_var(rough, mp)
          call alloc_cbm_var(soil, mp)
          call alloc_cbm_var(ssoil, mp)
          call alloc_cbm_var(sum_flux, mp)
          call alloc_cbm_var(veg, mp)
       end if
!
!      set meteorological forcing
!

       jyear=kdate/10000
       jmonth=(kdate-jyear*10000)/100
       jday=kdate-jyear*10000-jmonth*100
       jhour=ktime/100
       jmin=ktime-jhour*100
       ! mins from start of year
       mstart=1440*(ndoy(jmonth)+jday-1) + 60*jhour + jmin
       nperday=nint(24.*3600./dt)
       ktauplus=0
       do k=1,jmonth-1
        ktauplus = ktauplus + imonth(k)*nperday
       enddo

       ! Initialise sum flux variables
       IF (ktau==1) THEN
          ! rml added intialisation of sumpn, sumrp, sumrs and sumrd
          sumpn = 0.
          sumrp = 0.
          sumrs = 0.
          sumrd = 0.
          sum_flux%sumrpw = 0.
          sum_flux%sumrpr = 0.
          sum_flux%dsumpn = 0.
          sum_flux%dsumrp = 0.
          sum_flux%dsumrd = 0.
       END IF

       print *,'jyear,jmonth',jyear,jmonth,ktauplus
       ! mtimer contains number of minutes since the start of the run.
       mins = mtimer + mstart
       bpyear = 0.
       fjd = float(mod(mins,525600))/1440.  ! 525600 = 1440*365
       call solargh(fjd,bpyear,r1,dlt,alp,slag)
       dhr = 1.e-6
       call zenith(fjd,r1,dlt,slag,rlatt, &
     &               rlongg,dhr,ifull,coszro2,taudar2)


       call setco2for(jyear)

       !  rml: these not moved to cbm_pack because mins, theta, vmod, 
       !      atmco2, coszro2 not accessible there
       met%doy = float(mod(mins,24*60*365))/(24.*60.)
       met%tk = pack(theta,land)
       met%tc = met%tk - 273.16
       met%ua = pack(vmod,land)
       where (met%ua < 1.) met%ua = 1.
       met%ca = 1.e-6 * pack(atmco2,land)
       met%coszen = pack(coszro2,land)   ! use instantaneous value
       where (met%coszen < 1.e-8) met%coszen=1.e-8

       kstart = 1
       call cbm_pack(air, bgc, canopy, met, bal, rad, &
     &          rough, soil, ssoil, sum_flux, veg )
 

!      rml 21/09/07 remove ktauplus+ktau due to change in cable_offline
       CALL cbm(ktau, kstart, ntau, dt, air, bgc, canopy, met, &
     &      bal, rad, rough, soil, ssoil, sum_flux, veg, nvegt, mxst)

       call cbm_unpack(ktauplus+ktau, air, bgc, canopy, met, bal, rad, &
     &          rough, soil, ssoil, sum_flux, veg)

        
      return
      end subroutine sib4

! *************************************************************************************
      subroutine setco2for(jyear)
!     set co2 forcing for cable
!     constant: atmospheric co2 = 360 ppm 
!     prescribed: atmospheric co2 follows prescribed trend from 1900-2004
!                 based on ice core and Mauna Loa/South Pole data
!     interactive: atmospheric co2 taken from tracer (usually cable+fos+ocean)
      integer, parameter :: constantCO2 = 1
      integer, parameter :: prescribedCO2 = 2
      integer, parameter :: interactiveCO2 = 3
      include 'tracers.h'

!     rml added values for 2000-2004 using scripps records from cdiac
!     0.75*mauna loa + 0.25*south pole
!     2005 from gv07, noaa flask record
      real co2for(0:105)
      data co2for/296.0049,296.3785,296.7731,297.1795,297.5887,297.9919 &
     &           ,298.3842,298.7654,299.1349,299.4925,299.838 ,300.1709 &
     &           ,300.491 ,300.801 ,301.106 ,301.4113,301.7205,302.0357 &
     &           ,302.3587,302.6915,303.036,303.3941,303.7677,304.1587 &
     &           ,304.569,304.9971,305.4388,305.8894,306.3444,306.7992 &
     &           ,307.2494,307.6902,308.117,308.521,308.8895,309.2135 &
     &           ,309.4877,309.7068,309.8658,309.9668,310.019,310.0358 &
     &           ,310.035,310.0345,310.0522,310.0989,310.1794,310.2977 &
     &           ,310.4581,310.6661,310.928,311.2503,311.6395,312.1015 &
     &           ,312.6341,313.227,313.8694,314.5506,315.2599,315.9866 &
     &           ,316.7167,317.4268,318.106,318.7638,319.4402,320.1901 &
     &           ,321.0398,321.9713,322.9779,324.045,325.1356,326.2445 &
     &           ,327.3954,328.5711,329.723,330.8865,332.1331,333.5012 &
     &           ,334.9617,336.4614,337.9588,339.4225,340.8724,342.3488 &
     &           ,343.8625,345.421,347.0481,348.7452,350.4397,352.0169 &
     &           ,353.4269,354.6917,355.8849,357.1348,358.5514,360.1441 &
     &           ,361.8578,363.6479,365.4682,367.2137 &
     &           ,368.87,370.35,372.49,374.93,376.69,379.09/

      select case (CO2forcingtype)
        case (constantCO2); atmco2 = 360.
        case (prescribedCO2) 
          ipco2 = jyear-1900
          if (ipco2.lt.0.or.ipco2.gt.105) stop 'year out of range for&
     & co2 forcing for cable'
          atmco2 = co2for(ipco2)
        case (interactiveCO2)
          write(6,*) 'need to replace with tracer sum'
          atmco2 = tr(1:ifull,1,1)
      end select

      return
      end subroutine setco2for


!*******************************************************************************
      subroutine cbmrdn(nveg)
!      subroutine cbmrdn(imonth)
!     BP removed the redundant parameter imonth and added nveg to be passed
!     upward to globpe.f (BP jan 2008)
!
!     reads the parameters required by land surface scheme 
!     called from indata

      use cc_mpi
      use define_dimensions, only : ncs, ncp
      include 'newmpar.h'
      include 'arrays.h'
      include 'carbpools.h'
      include 'dates.h'     ! dt,kdate,iyd
      include 'nsibd.h'     ! ivegt,isoilm
      include 'parm.h'      ! id,jd,idjd
      include 'permsurf.h'
      include 'map.h'      ! id,jd,idjd
      include 'pbl.h'       ! tss
      include 'sigs.h'
      include 'soil.h'      ! land, zolnd passed to indata
      include 'soilv.h'
      include 'soilsnow.h'  ! initial soil temperatures and moistures
      include 'vegpar.h'

      integer nveg
      integer jyear,jmonth
      character*80 comments
      character*2 chflag
      integer, dimension(ifull_g) :: ivegt_g
!     rml 28/09/07 igbp veg types and frac-c4 defined all gridpoints
!     c4frac defined in vegpar.h, c4frac_g defined here for mpi stuff
      logical igbp
      real, dimension(ifull_g) :: c4frac_g
!     rml 01/10/07 add alternative parameter file format
      real notused
      character*10 vegtypetmp
      character*25 vegnametmp
      character*41 vegtypecheck
      data vegtypecheck/'forestdeciduousshrubgrasscropnovegnotused'/

      if ( myid == 0) then
         open(unit=8,file=vegtypefile,status='old')
!     rml 28/09/07 add alternative read for igbp file
         read(8,*) comments
         print *,'read vegetation types'
         if (comments(1:4).eq.'IGBP') then
           igbp=.true.
           print *, comments(1:4)
         else
           igbp=.false.
           print *, 'CASA'
           c4frac = 0.
         endif
         do ii=1,ifull_g
            if (igbp) then
              read(8,*) iq,i,j,rl1,rl2,ivnew,c4
              c4frac_g(iq) =c4
            else
              read(8,*) iq,i,j,rl1,rl2,ivold,ivnew
            endif
            ivegt_g(iq)=ivnew
!           if(land(iq)) print *,iq,i,j,rl1,rl2,ivold,ivegt(iq)
         enddo                  ! ii
         close(8)
         call ccmpi_distribute(ivegt, ivegt_g)
         if (igbp) call ccmpi_distribute(c4frac, c4frac_g)
      else
         call ccmpi_distribute(ivegt)
         if (igbp) call ccmpi_distribute(c4frac)
      end if


      ! This is a small file, so simpler to let every processor read if
!     rml 01/10/07 old or new format given in cableinput namelist
      if (vegparmnew) then
!       read vegetation parameters from new format file
        open(unit=8,file=vegparmfile,status='old')
        read(8,*) comments
        if ( myid == 0 ) write(*,802) comments
        read(8,*) nveg
        if ( myid == 0 )  print *,'nveg=',nveg
        if (maxval(ivegt).gt.nveg) stop 'veg parameter file has less &
     & vegetation types than vegtype.dat'
        if (nveg.gt.mxvt) stop 'increase mxvt in newmpar.h'
        allocate(vegname(nveg),vegtype(nveg))
        do n=1,nveg
          read(8,*) jveg,vegtypetmp,vegnametmp
          if (jveg.gt.nveg) stop 'jveg out of range in paramgeter file'
          vegname(jveg) = vegnametmp
          vegtype(jveg) = vegtypetmp
          if (index(vegtypecheck,vegtypetmp).eq.0) write(6,*) &
     & 'WARNING: ',vegtypetmp,' is not a recognised vegetation type'
          read(8,*) hc(jveg),xfang(jveg),notused,dleaf(jveg)
          read(8,*)  ! rholeaf not used
          read(8,*)  ! tauleaf not used
          read(8,*) notused,notused,notused,xalbnir(jveg) ! rhosoil not used
          read(8,*) notused,wai(jveg),canst1(jveg),shelrb(jveg), &
     &              vegcf(jveg),extkn(jveg)
          read(8,*) vcmax(jveg),rp20(jveg),rpcoef(jveg),rs20(jveg)
          read(8,*) tminvj(jveg),tmaxvj(jveg),vbeta(jveg),rootbeta(jveg)
          read(8,*) tcplant(jveg,1:3),tcsoil(jveg,1:2)
!  rates not currently set to vary with veg type
          read(8,*) ratecp(1:3),ratecs(1:2) 
        enddo
        ejmax = 2.*vcmax
        if ( myid == 0 ) then
          print *, 'vegname',vegname(1:nveg)
          print *, 'vegtype',vegtype(1:nveg)
          print *, 'canst1',canst1(1:nveg)
          print *, 'dleaf',dleaf(1:nveg)
          print *, 'xalbnir',xalbnir(1:nveg)
          print *, 'vcmax',vcmax(1:nveg)
          print *, 'ejmax',ejmax(1:nveg)
          print *, 'hc',hc(1:nveg)
          print *, 'xfang',xfang(1:nveg)
          print *, 'rp20',rp20(1:nveg)
          print *, 'rpcoef',rpcoef(1:nveg)
          print *, 'rprs20',rs20(1:nveg)
          print *, 'wai',wai(1:nveg)
          print *, 'shelrb',shelrb(1:nveg)
          print *, 'vegcf',vegcf(1:nveg)
          print *, 'extkn',extkn(1:nveg)
          print *, 'tminvj',tminvj(1:nveg)
          print *, 'tmaxvj',tmaxvj(1:nveg)
          print *, 'rootbeta',rootbeta(1:nveg)
          print *,'cplant 1',tcplant(1:nveg,1)
          print *,'cplant 2',tcplant(1:nveg,2)
          print *,'cplant 3',tcplant(1:nveg,3)
          print *,'csoil 1',tcsoil(1:nveg,1)
          print *,'csoil 2',tcsoil(1:nveg,2)
          print *,'ratecp',(ratecp(1:3))
          print *,'ratecs',(ratecs(1:2))
        endif
        close(8)
      else
!       read vegetation parameters from old format file
        open(unit=8,file=vegparmfile,status='old')
        read(8,*) comments
        if ( myid == 0 ) write(*,802) comments
802     format(1x,a80)
        call comskp(8)
        read(8,*) nveg
!       rml 28/09/07 compare nveg and max veg type number
        if (maxval(ivegt).gt.nveg) stop 'veg_parm.txt has less vegetation &
     & types than vegtype.dat'
        if ( myid == 0 )  print *,'nveg=',nveg
        call comskp(8)
        read(8,*) (canst1(jveg),jveg=1,nveg)
        if ( myid == 0 )  print *, 'canst1',(canst1(jveg),jveg=1,nveg)
!      read(8,*) (cansto(jveg),jveg=1,nveg)
!      print *, 'cansto',(cansto(jveg),jveg=1,nveg)
        read(8,*) (dleaf(jveg),jveg=1,nveg)
        if ( myid == 0 ) print *, 'dleaf',(dleaf(jveg),jveg=1,nveg)
        read(8,*) (vcmax(jveg),jveg=1,nveg)
        if ( myid == 0 ) print *, 'vcmax',(vcmax(jveg),jveg=1,nveg)
        do jveg=1,nveg
          ejmax(jveg)=2.*vcmax(jveg)
        enddo
        if ( myid == 0 ) print *, 'ejmax',(ejmax(jveg),jveg=1,nveg)
        read(8,*) (hc(jveg),jveg=1,nveg)
        if ( myid == 0 ) print *, 'hc',(hc(jveg),jveg=1,nveg)
        read(8,*) (xfang(jveg),jveg=1,nveg)
        if ( myid == 0 ) print *, 'xfang',(xfang(jveg),jveg=1,nveg)
        read(8,*) (rp20(jveg),jveg=1,nveg)
        if ( myid == 0 ) print *, 'rp20',(rp20(jveg),jveg=1,nveg)
        read(8,*) (rpcoef(jveg),jveg=1,nveg)
        if ( myid == 0 ) print *, 'rpcoef',(rpcoef(jveg),jveg=1,nveg)
        read(8,*) (rs20(jveg),jveg=1,nveg)
        if ( myid == 0 ) print *, 'rprs20',(rs20(jveg),jveg=1,nveg)
        read(8,*) (shelrb(jveg),jveg=1,nveg)
        if ( myid == 0 ) print *, 'shelrb',(shelrb(jveg),jveg=1,nveg)
        read(8,*) (frac4(jveg),jveg=1,nveg)
        read(8,*) (tminvj(jveg),jveg=1,nveg)
        if ( myid == 0 ) print *, 'tminvj',(tminvj(jveg),jveg=1,nveg)
        read(8,*) (tmaxvj(jveg),jveg=1,nveg)
        if ( myid == 0 ) print *, 'tmaxvj',(tmaxvj(jveg),jveg=1,nveg)
        read(8,*) (vbeta(jveg),jveg=1,nveg)
        if ( myid == 0 ) print *, 'vbeta',(vbeta(jveg),jveg=1,nveg)
        read(8,*) (tcplant(jveg,1),jveg=1,nveg)
        read(8,*) (tcplant(jveg,2),jveg=1,nveg)
        read(8,*) (tcplant(jveg,3),jveg=1,nveg)
        if ( myid == 0) print *,'cplant 1',(tcplant(jveg,1),jveg=1,nveg)
        if ( myid == 0) print *,'cplant 2',(tcplant(jveg,2),jveg=1,nveg)
        if ( myid == 0) print *,'cplant 3',(tcplant(jveg,3),jveg=1,nveg)
        read(8,*) (tcsoil(jveg,1),jveg=1,nveg)
        read(8,*) (tcsoil(jveg,2),jveg=1,nveg)
        if ( myid == 0 ) print *,'csoil 1',(tcsoil(jveg,1),jveg=1,nveg)
        if ( myid == 0 ) print *,'csoil 2',(tcsoil(jveg,2),jveg=1,nveg)
        call comskp(8)
        read(8,*) (ratecp(j),j=1,ncp)
        if ( myid == 0 ) print *, 'ratecp',(ratecp(j),j=1,ncp)
        call comskp(8)
        read(8,*) (ratecs(j),j=1,ncs)
        if ( myid == 0 ) print *, 'ratecs',(ratecs(j),j=1,ncs)
        close(8)
!       set vegcf and xalbnir here as they are not in old format file;
!       later on, also true values for wai (not currently used)
        allocate(vegtype(nveg))
        vegtype(:) = 'others'   ! BP added initialization (1Feb2008)
        SELECT CASE (nveg)
          CASE (13)     ! CASA vegetation types
            vegcf = (/ 1.95, 1.5, 1.55, 0.91, 0.73, 2.8, 2.75, 0.0,  &
                     & 2.05, 0.6, 0.4, 2.8, 0.0, 0.0, 0.0, 0.0, 0.0 /)
            xalbnir = (/ 0.96, 1.0, 1.08, 0.79, 0.81, 1.02, 1.23, 1.0, &
                       & 1.20, 1.14, 1.15, 0.98, 1.00, 1.00, 1.00, 1.00, 1.00 /)
            vegtype(2)='deciduous'
            vegtype(5)='deciduous'
          CASE (16,17)  ! IGBP vegetation types without/with water bodies
!            vegcf = (/ 11.82, 13.06, 6.71, 11.34, 8.59, 0.6, 2.46, 10.0,  &
!                     & 15.93, 3.77, 0.0, 11.76, 0.0, 2.8, 10.0, 10.0, 0.0 /)
            vegcf = (/ 0.91, 1.95, 0.73, 1.50, 1.55, 0.6, 2.05, 2.80,  &
                     & 2.75, 2.75, 0.0, 2.80, 0.0, 2.80, 0.0, 0.4, 0.4 /)
            xalbnir = (/ 0.79, 0.96, 0.81, 1.0, 1.08, 1.14, 1.20, 1.02, &
                       & 1.23, 1.16, 0.89, 0.98, 1.10, 1.13, 1.00, 1.15, 1.00 /)
            vegtype(3)='deciduous'
            vegtype(4)='deciduous'
          CASE DEFAULT
            PRINT *,'Error! Dimension not compatible with CASA or IGBP types!'
            PRINT *,'Dimension = ', nveg
            PRINT *,'At the vegcf section.'
            STOP
        END SELECT
        wai = 0.
        extkn = 0.4
        rootbeta = 0.

      endif

      jyear=kdate/10000
      jmonth=(kdate-jyear*10000)/100

!     if( jmonth .eq. 1 .and. jyear .eq. inyear_carb) then
      if (initcarbpools == 1) then
      do iq=1,ifull
      if(land(iq))then
       cplant(iq,1)=tcplant(ivegt(iq),1) 
       cplant(iq,2)=tcplant(ivegt(iq),2) 
       cplant(iq,3)=tcplant(ivegt(iq),3) 
      
       csoil(iq,1)=tcsoil(ivegt(iq),1) 
       csoil(iq,2)=tcsoil(ivegt(iq),2) 
       cansto(iq) = 0.
      endif !land
      enddo !iq
      endif ! jmonth
       if ( myid==0 ) print *,'cbmrdn',cplant(13419,1),cplant(13419,2), &
     &       cplant(13419,3),csoil(13419,1),csoil(13419,2)
!      call comskp(8)
!      read(8,*) (froot(j),j=1,ms)
!     rml 01/10/07 old format reads rate constants here, new format 
!     includes with veg types
!      if (.not.vegparmnew) then   ! BP moved this block up 56 lines (1Feb2008)
!        call comskp(8)
!        read(8,*) (ratecp(j),j=1,ncp)
!        if ( myid == 0 ) print *, 'ratecp',(ratecp(j),j=1,ncp)
!        call comskp(8)
!        read(8,*) (ratecs(j),j=1,ncs)
!        close(8)
!        if ( myid == 0 ) print *, 'ratecs',(ratecs(j),j=1,ncs)
!      endif
 
      end subroutine cbmrdn

! ************************************************************************

      SUBROUTINE COMSKP(IUNIT)
! MRR, 5-AUG-83
! SKIPS COMMENT LINES IN CONTROL DATA FILE, STOPS IF EOF FOUND
      CHARACTER*1 COM
1     READ(IUNIT,100,END=2) COM
100   FORMAT(A1)
      IF(COM.EQ.'C'.or.com.eq.'c') GOTO 1
      BACKSPACE IUNIT
      RETURN
2     STOP 'CONTROL FILE EOF'
      END subroutine comskp

! ************************************************************************

      subroutine rlaiday
!     called from globpe.f, interpolates monthly lai to daily

      include 'newmpar.h'
      include 'nsibd.h'
      include 'dates.h'     !  kdate,ktime,timer,mtimer
      include 'parm.h'      ! id,jd
      include 'permsurf.h'  ! iperm etc
      include 'vegpar.h'  
      include 'mpif.h'
      integer, parameter, dimension(0:13) :: mdays = &
     &     (/ 31, 31,28,31,30,31,30,31,31,30,31,30,31, 31 /)
      integer iyr, imo, iday, ijd
      save iyr,imo,iday

      ijd = 9061
      print *,'rlaiday  called'
      iyr=kdate/10000
      imo=(kdate-10000*iyr)/100
      iday=kdate-10000*iyr-100*imo  +mtimer/(60*24)
      do while (iday>mdays(imo))
       iday=iday-mdays(imo)
       imo=imo+1
       if(imo>12)then
         imo=1               ! no leap years for the moment
         iyr=iyr+1
       endif
       enddo
!       print *,'iday,imo,iyr',iday,imo,iyr

       if(iday<mdays(imo)/2)then  ! 1st half of month
        rat1=(mdays(imo)-2.*iday)/(mdays(imo)+mdays(imo-1))
        rat2=(2.*iday+mdays(imo-1))/(mdays(imo)+mdays(imo-1))
        rlai=rat1*rlai123(:,1)+rat2*rlai123(:,2) 
       else                             ! 2nd half of month
        rat1=(mdays(imo+1)+2.*mdays(imo)-2.*iday)/ &
     &                               (mdays(imo+1)+mdays(imo))
        rat2=(2.*iday-mdays(imo))/(mdays(imo+1)+mdays(imo))
        rlai=rat1*rlai123(:,2)+rat2*rlai123(:,3) 
       endif
       if(iday==mdays(imo)/2) rlai=rlai123(:,2)
       do iq=1,ifull
         rlai(iq)=min(7.,max(0.0011,rlai(iq)))
         if (isoilm(iq).eq.9) rlai(iq)=0.0011
         vlai(iq)=rlai(iq)
       enddo

       print *,'month_imo,iday',imo,iday,rlai(ijd),rlai123(ijd,1), &
     &       rlai123(ijd,2),rlai123(ijd,3)
      return
      end subroutine rlaiday

  ! ****************************************************************************
  subroutine cbm_pack(air, bgc, canopy, met, bal, rad,  &
                    & rough, soil, ssoil, sum_flux, veg)

      use define_types, cbm_ms => ms
      use roughness_module
      use air_module
      use radiation_module
      use canopy_module
      use soil_snow_module             

      implicit none

      type (air_type)            :: air 
      type (canopy_type)         :: canopy     ! updated each timestep
      type (balances_type)       :: bal        ! energy/water bal variables
      type (met_type)            :: met ! met forcing from gcm each timestep
      type (soil_parameter_type) :: soil       ! soil parameters
      type (soil_snow_type)      :: ssoil      ! various data for land point
      type (sum_flux_type)       :: sum_flux   ! accumulated carbon fluxes
      type (radiation_type)      :: rad        ! updated each timestep
      type (roughness_type)      :: rough
      type (veg_parameter_type)  :: veg        ! vegetation parameters
      type (bgc_pool_type)       :: bgc        ! bgc parameters

      INTEGER(i_d) :: ip,iq,j,k  

      include 'newmpar.h' 
      include 'aalat.h'
      include 'arrays.h'
      include 'const_phys.h' ! grav
      include 'dates.h'
      include 'extraout.h'   ! sgsave, rgsave
      include 'morepbl.h' 
      include 'nsibd.h' 
      include 'carbpools.h'
      include 'permsurf.h'  
      include 'parm.h'  
      include 'sigs.h'       ! sig
      include 'soilbal.h'
      include 'soilsnow.h'
      include 'soilv.h' !zse,zshh,bch,css,rhos,cnsd,hyds,sucs,hsbh,sfc,swilt,ibp2,i2bp3
      include 'vegpar.h'
      include 'soil.h'

      REAL(r_1) :: totdepth,totroot(mxvt)   ! local variable for froot calc

!     these constant throughout simulation so only need to be set once
      if (ktau == 1) then
!       vegetation parameters
        veg%iveg   = pack(ivegt, land)
        veg%canst1 = canst1(veg%iveg)
        veg%ejmax  =  ejmax(veg%iveg)
        veg%hc     =     hc(veg%iveg)
        veg%frac4  =  frac4(veg%iveg)
!       rml 28/09/07 c4 fraction defined for each grid point, for
!       compatibility just overwrite values from parameter file if
!       find data in c4frac array
        if (sum(c4frac).gt.0) veg%frac4 = pack(c4frac,land)
        veg%tminvj = tminvj(veg%iveg)
        veg%tmaxvj = tmaxvj(veg%iveg)
        veg%vbeta  =  vbeta(veg%iveg)
        veg%rp20   =   rp20(veg%iveg)
        veg%rpcoef = rpcoef(veg%iveg)
        veg%shelrb = shelrb(veg%iveg)
        veg%vcmax  =  vcmax(veg%iveg)
        veg%xfang  =  xfang(veg%iveg)
        veg%dleaf  =  dleaf(veg%iveg)
        veg%xalbnir = xalbnir(veg%iveg)
        veg%wai    =    wai(veg%iveg)
        veg%vegcf  =  vegcf(veg%iveg)
        veg%extkn  =  extkn(veg%iveg)
        veg%deciduous = .false.
        where (vegtype(veg%iveg).eq.'deciduous') veg%deciduous=.true.

!       soil parameters
        soil%zse  = zse
        soil%zshh = zshh
        soil%zse = (/.022, .058, .154, .409, 1.085, 2.872/) ! soil layer thickness
        soil%zshh(1) = 0.5 * soil%zse(1)
        soil%zshh(ms+1) = 0.5 * soil%zse(ms)
        soil%zshh(2:ms) = 0.5 * (soil%zse(1:ms-1) + soil%zse(2:ms))

        soil%isoilm  =  pack(isoilm, land)
        soil%albsoil =  pack(albsoil, land) 
        soil%bch     =   bch(soil%isoilm)
        soil%css     =   css(soil%isoilm)
        soil%rhosoil =  rhos(soil%isoilm)
        soil%cnsd    =  cnsd(soil%isoilm)
        soil%hyds    =  hyds(soil%isoilm)
        soil%sucs    =  sucs(soil%isoilm)
        soil%hsbh    =  hsbh(soil%isoilm)
        soil%sfc     =   sfc(soil%isoilm)
        soil%ssat    =  ssat(soil%isoilm)
        soil%swilt   = swilt(soil%isoilm)
        soil%ibp2    =  ibp2(soil%isoilm)
        soil%i2bp3   = i2bp3(soil%isoilm)
        soil%rs20    =  rs20(veg%iveg)

        if (sum(rootbeta).eq.0) then
!        assume rootbeta not in parameter file and need to prescribe froot here
         IF (MAXVAL(ivegt) == 13) THEN   ! CASA vegetation types  
          froot(:,1) = (/.02,.04,.04,.04,.04,.05,.05,.05,.05,.05,.05,.05,.05,.01,.01,.01,.01/)
          froot(:,2) = (/.06,.11,.11,.11,.11,.15,.15,.10,.10,.10,.10,.15,.15,.01,.01,.01,.01/)
          froot(:,3) = (/.14,.20,.20,.20,.20,.35,.35,.35,.35,.35,.35,.34,.35,.01,.01,.01,.01/)
          froot(:,4) = (/.28,.26,.26,.26,.26,.39,.39,.35,.35,.35,.35,.38,.40,.01,.01,.01,.01/)
          froot(:,5) = (/.35,.24,.24,.24,.24,.05,.05,.10,.10,.10,.10,.06,.04,.01,.01,.01,.01/)
          froot(:,6) = (/.15,.15,.15,.15,.15,.01,.01,.05,.05,.05,.05,.02,.01,.01,.01,.01,.01/)
         ELSEIF (MAXVAL(ivegt) == 16 .OR. MAXVAL(ivegt) == 17) THEN
          ! IGBP vegetation types without/with water bodies
          froot(:,1) = (/.04,.02,.04,.04,.04,.05,.05,.05,.05,.05,.05,.05,.05,.05,.05,.05,.01/)
          froot(:,2) = (/.11,.06,.11,.11,.11,.10,.10,.15,.15,.15,.10,.15,.10,.15,.15,.10,.01/)
          froot(:,3) = (/.20,.14,.20,.20,.20,.35,.35,.35,.35,.35,.35,.34,.35,.34,.35,.35,.01/)
          froot(:,4) = (/.26,.28,.26,.26,.26,.35,.35,.39,.39,.39,.35,.38,.35,.38,.40,.35,.01/)
          froot(:,5) = (/.24,.35,.24,.24,.24,.10,.10,.05,.05,.05,.10,.06,.10,.06,.04,.10,.01/)
          froot(:,6) = (/.15,.15,.15,.15,.15,.05,.05,.01,.01,.01,.05,.02,.05,.02,.01,.05,.01/)
         ELSE
          PRINT *, 'Error! Dimension not compatible with CASA or IGBP types!'
          PRINT *, 'Dimension = ', MAXVAL(veg%iveg)
          PRINT *, 'Within cbm_pack.'
          STOP
         END IF
        else
!        preferred option
!        froot is now calculated from soil depth and the new parameter rootbeta 
!        according to Jackson et al. 1996, Oceologica, 108:389-411
         totroot(:) = 1.0-rootbeta(:)**(sum(soil%zse)*100.0)
         totdepth = 0.0
         do k=1,ms
           totdepth = totdepth + soil%zse(k)*100.0
           froot(:,k) = min(1.0,1.0-rootbeta(:)**totdepth)
         enddo
         do k = ms, 2, -1
           froot(:,k) = froot(:,k) - froot(:,k-1)
         enddo
        endif
        veg%froot = froot(veg%iveg,:)

        rad%latitude = pack(alat,land)

      endif  ! ktau=1

! rml 18/10/07 moved this code from sflux.f (and cleaned up)
      met%fld = -1. * pack(rgsave,land)        ! long wave down  
      met%qv = pack(qg(1:ifull,1),land)        ! specific humidity in kg/kg
      met%pmb = .01*pack(ps(1:ifull),land)     ! pressure in mb at ref height
      met%precip = pack(condx,land)
      ! name changed to precip_s (EK nov2007)
      met%precip_s = 0.0                        ! in mm not mm/sec
      where ( met%tc < 0.0 ) met%precip_s = met%precip
      do ip=1,ipland
        iq = iperm(ip)
          met%hod(ip)=(met%doy(ip)-int(met%doy(ip)))*24.0 + along(iq)/15.
          if (met%hod(ip).gt.24.0) met%hod(ip)=met%hod(ip)-24.0
          rough%za(ip) = -287.*t(iq,1)*log(sig(1))/grav   ! reference height
          met%fsd(ip) = sgsave(iq)/(1.-alb(iq))! short wave down (positive) W/m^2
      enddo
!     met%hod = (met%doy - int(met%doy))*24.0 + along(iperm)/15.
!     where (met%hod > 24.0) met%hod = met%hod - 24.0
!     met%fsd = sgsave(iperm)/(1.-alb(iperm))  ! short wave down (positive) W/m^2

!     rough%za = -287.*t(iperm,1)*log(sig(1))/grav   ! reference height
      write(6,*) along(iperm(1)),met%hod(1),met%fsd(1),rough%za(1)
      write(6,*) along(iperm(mp)),met%hod(mp),met%fsd(mp),rough%za(mp)

      ssoil%albsoilsn(:,1) = pack(albsoilsn(:,1), land)
      ssoil%albsoilsn(:,2) = pack(albsoilsn(:,2), land)
      ssoil%albsoilsn(:,3)   = 0.05
      do k = 1,ms
        ssoil%tgg(:,k) = pack(tgg(:,k), land)
        ssoil%wb(:,k) = pack(wb(:,k), land)
        ssoil%wbice(:,k) = pack(wbice(:,k), land)
      enddo
! rml check ssoil%wbtot calculation may be redone in cable_soilsnow anyway
! is bal%wbtot0 ever used?
      if (ktau == 1) then
        ssoil%wbtot = 0.
        DO j=1,ms
          ssoil%wbtot  = ssoil%wbtot + ssoil%wb(:,j) * soil%zse(j) * 1000.0
        END DO
        bal%wbtot0 = ssoil%wbtot
      endif

!     bal%wbtot0 = pack(wbtot0, land) !think this is redundant
!  are any of these balance variables used?
      bal%evap_tot = pack(tevap, land)
      bal%precip_tot = pack(tprecip, land)
      bal%ebal_tot = pack(totenbal, land)
      bal%rnoff_tot = pack(trnoff, land)

      ssoil%wbtot = pack(wbtot, land)
      ssoil%isflag = pack(isflag, land)
      do k = 1,3
        ssoil%tggsn(:,k) = pack(tggsn(:,k), land)
        ssoil%smass(:,k) = pack(smass(:,k), land)
        ssoil%ssdn(:,k) = pack(ssdn(:,k), land)
      enddo

      ssoil%ssdnn = pack(ssdnn, land)
      ssoil%snowd = pack(snowd, land)
      ssoil%osnowd = pack(osnowd, land)
      bal%osnowd0 = pack(osnowd0, land)
      ssoil%snage = pack(snage, land)
      ssoil%isflag = pack(isflag, land)
      ssoil%rtsoil = pack(rtsoil, land)

      canopy%ghflux = pack(gflux, land)
      canopy%sghflux = pack(sgflux, land)
      canopy%cansto = pack(cansto, land)
      veg%vlai = pack(vlai, land)

      sum_flux%sumpn = pack(sumpn, land)
      sum_flux%sumrp = pack(sumrp, land)
      sum_flux%sumrs = pack(sumrs, land)
      sum_flux%sumrd = pack(sumrd, land)
      sum_flux%sumrpw = pack(sumrpw, land)
      sum_flux%sumrpr = pack(sumrpr, land)
      sum_flux%dsumpn = pack(dsumpn, land)
      sum_flux%dsumrp = pack(dsumrp, land)
      sum_flux%dsumrs = pack(dsumrs, land)
      sum_flux%dsumrd = pack(dsumrd, land)
      do k=1,ncp
        bgc%cplant(:,k) = pack(cplant(:,k), land)
        bgc%ratecp(k) = ratecp(k)
      enddo
      do k=1,ncs
        bgc%csoil(:,k) = pack(csoil(:,k), land)
        bgc%ratecs(k) = ratecs(k)
      enddo

      end subroutine cbm_pack

! **************************************************************************
      subroutine cbm_unpack(ktauyear, air, bgc, canopy, met, bal, rad,  &
               rough, soil, ssoil, sum_flux, veg)

      use define_types, cbm_ms => ms
      use roughness_module
      use air_module
      use radiation_module
      use canopy_module
      use soil_snow_module

!  implicit none

      INTEGER(i_d), INTENT(IN)    :: ktauyear
      type (air_type)             :: air
      type (canopy_type)          :: canopy     ! updated each timestep
      type (balances_type)        :: bal        ! energy/water bal variables
      type (met_type)             :: met ! met forcing from gcm (timestep at a time)
      type (soil_parameter_type)  :: soil       ! soil parameters
      type (soil_snow_type)       :: ssoil        ! various data for land point
      type (sum_flux_type)        :: sum_flux
      type (radiation_type)       :: rad        ! updated each timestep
      type (roughness_type)       :: rough
      type (veg_parameter_type)   :: veg        ! vegetation parameters
      type (bgc_pool_type)        :: bgc        ! bgc parameters

      INTEGER(i_d) :: ip,iq,j,k

      include 'newmpar.h'
      include 'aalat.h'
      include 'arrays.h'
      include 'dates.h'
      include 'morepbl.h'
      include 'nsibd.h'
      include 'carbpools.h'
      include 'pbl.h'    ! new cduv for vetmix
      include 'permsurf.h'
      include 'soilbal.h'
      include 'soilsnow.h'
      include 'soilv.h'
      include 'screen.h'
      include 'vegpar.h'
      include 'soil.h'

      do k=1,2
        albsoilsn(:,k) = unpack(ssoil%albsoilsn(:,k), land, albsoilsn(:,k))
        albvisnir(:,k) = unpack(rad%albedo(:,k), land, albvisnir(:,k))
      enddo
      runoff= unpack(ssoil%runoff, land, runoff)
      rnof1= unpack(ssoil%rnof1, land, rnof1)
      rnof2= unpack(ssoil%rnof2, land, rnof2)
      wbtot= unpack(ssoil%wbtot, land, wbtot)
      tevap = unpack(bal%evap_tot, land, tevap)
      tprecip = unpack(bal%precip_tot, land, tprecip)
      totenbal = unpack(bal%ebal_tot, land, totenbal)
      trnoff = unpack(bal%rnoff_tot, land, trnoff)
      do k = 1,ms
        tgg(:,k)= unpack(ssoil%tgg(:,k), land, tgg(:,k))
        wb(:,k)= unpack(real(ssoil%wb(:,k),r_1), land, wb(:,k))
        wbice(:,k)= unpack(real(ssoil%wbice(:,k),r_1), land, wbice(:,k))
      enddo
      do k = 1,3
        tggsn(:,k)= unpack(ssoil%tggsn(:,k), land, tggsn(:,k))
        smass(:,k)= unpack(ssoil%smass(:,k), land, smass(:,k))
        ssdn(:,k)= unpack(ssoil%ssdn(:,k), land, ssdn(:,k))
      enddo

      ssdnn= unpack(ssoil%ssdnn, land, ssdnn)
      snowd= unpack(ssoil%snowd, land, snowd)
      osnowd= unpack(ssoil%osnowd, land, osnowd)
      osnowd0= unpack(bal%osnowd0, land, osnowd0)
      snage= unpack(ssoil%snage, land, snage)
      isflag= unpack(ssoil%isflag, land, isflag)
      rtsoil= unpack(ssoil%rtsoil, land, rtsoil)
      rnet = unpack(canopy%rnet,  land, rnet)
      fg = unpack(canopy%fh,  land, fg)
      eg = unpack(canopy%fe,  land, eg)
      epot = unpack(ssoil%potev,  land, epot)
      tss = unpack(rad%trad,  land, tss)
      tscrn = unpack(canopy%tscrn,  land, tscrn)
      qgscrn = unpack(canopy%qscrn,  land, qgscrn)
      uscrn = unpack(canopy%uscrn,  land, uscrn)
      cduv= unpack(canopy%cduv, land, cduv)
      cansto= unpack(canopy%cansto, land, cansto)
      vlai= unpack(veg%vlai, land, vlai)
      gflux = unpack(canopy%ghflux, land, gflux)
      sgflux = unpack(canopy%sghflux, land, sgflux)
      rtsoil = unpack(ssoil%rtsoil, land, rtsoil)
      fnee= unpack(canopy%fnee, land, fnee)
      fpn= unpack(canopy%fpn, land, fpn)
      frd= unpack(canopy%frday, land, frd)
      frp= unpack(canopy%frp, land, frp)
      frpw= unpack(canopy%frpw, land, frpw)
      frpr= unpack(canopy%frpr, land, frpr)
      frs= unpack(canopy%frs, land, frs)
      sumpn= unpack(sum_flux%sumpn, land, sumpn)
      sumrp= unpack(sum_flux%sumrp, land, sumrp)
      sumrpw= unpack(sum_flux%sumrpw, land, sumrpw)
      sumrpr= unpack(sum_flux%sumrpr, land, sumrpr)
      sumrs= unpack(sum_flux%sumrs, land, sumrs)
      sumrd= unpack(sum_flux%sumrd, land, sumrd)
! rml - if not bothering to unpack these here (or write them out), perhaps 
! should delete everywhere
!    dsumpn= unpack(sum_flux%dsumpn, land, dsumpn)
!    dsumrp= unpack(sum_flux%dsumrp, land, dsumrp)
!    dsumrs= unpack(sum_flux%dsumrs, land, dsumrs)
!    dsumrd= unpack(sum_flux%dsumrd, land, dsumrd)
      do k=1,3
        cplant(:,k)= unpack(bgc%cplant(:,k), land, cplant(:,k))
! need? rates don't change?
        ratecp(k)= bgc%ratecp(k)
      enddo
      do k=1,2
        csoil(:,k)= unpack(bgc%csoil(:,k), land, csoil(:,k))
! need? rates don't change?
        ratecs(k)= bgc%ratecs(k)
      enddo

      if ( nproc == 1 ) then
         call eva_output( ktauyear, dels, air, bgc, &
                        canopy, met, bal, &
                        rad, rough, soil, ssoil, sum_flux, veg)
      else
         print*, "Skipping eva_output - not yet implemented in parallel version"
      end if

  END Subroutine cbm_unpack

! ******************************************************************************
  SUBROUTINE eva_output(ktauyear, dels, air, bgc, &
          canopy, met, bal, &
          rad, rough, soil, ssoil, sum_flux, veg)
  USE define_types, cbm_ms => ms ! Rename to avoic conflict with newmpar.h value
  USE air_module
  USE canopy_module
  USE carbon_module
  USE soil_snow_module
  USE physical_constants
  IMPLICIT NONE
!  INTEGER(i_d), INTENT(IN)              :: ktau ! integration step number
!  INTEGER(i_d), INTENT(IN)              :: kend
  INTEGER(i_d), INTENT(IN)              :: ktauyear
  REAL(r_1), INTENT(IN)                 :: dels ! integration time setp (s)
  TYPE (air_type), INTENT(INOUT)                :: air
  TYPE (bgc_pool_type), INTENT(INOUT)   :: bgc
  TYPE (canopy_type), INTENT(INOUT)     :: canopy
  TYPE (met_type), INTENT(INOUT)     :: met
!  TYPE (gcmmet_type), INTENT(INOUT)     :: met
!  TYPE (misc_output_type), INTENT(INOUT)        :: mout
  type (balances_type), INTENT(INOUT)   :: bal        ! energy/water bal variables

  TYPE (radiation_type), INTENT(INOUT)  :: rad
  TYPE (roughness_type), INTENT(INOUT)  :: rough
  TYPE (soil_parameter_type), INTENT(INOUT)               :: soil
  TYPE (soil_snow_type), INTENT(INOUT)  :: ssoil
  TYPE (sum_flux_type), INTENT(INOUT)   :: sum_flux
!  TYPE (veg_type), INTENT(INOUT)                :: veg
  type (veg_parameter_type), INTENT(INOUT)       :: veg        ! vegetation parameters
  REAL(r_1), dimension(mp)                :: xx1,xx2,xx3
  INTEGER(i_d)                          :: u         ! output unit
  INTEGER(i_d)                          :: ijd        
  LOGICAL                               :: is_open   ! Is file open?
  REAL(r_1), DIMENSION(mp)              :: owbtot
  REAL(r_1), DIMENSION(mp)              :: swnet
  REAL(r_1), DIMENSION(mp)              :: lwnet
  INTEGER :: k ! do loop counter
  INTEGER(i_d),DIMENSION(10) :: listp ! 

  include 'newmpar.h'
  include 'morepbl.h'
  include 'parm.h'

  listp(1)=2726
  listp(2)=3002
  listp(3)=3009
  listp(4)=2468
  listp(5)=2560
  listp(6)=1547
  listp(7)=3688



  ijd = listp(1)

!  u = 102
!  INQUIRE (u, opened=is_open)
!  IF (.not. is_open) THEN
!     OPEN (u, file='f102.txt', status='replace')
!     WRITE(u, '(a6,a7,a6,5a7,4a7,2a7,5a9,a12,3a7,2a7)') &
!          'ktau', 'fsd', 'fld', 'fnv', 'fev', 'fevc', 'fevw', 'fhv', &
!          'fns', 'fes', 'fhs', 'ga', 'fe', 'fh', 'drybal', 'wetbal', &
!          'ebalt', 'ebal', 'ebal', 'totenbal', 'trad', 'tv', 'tss', 'vlaiw', 'transd'
!  END IF
!  WRITE (u, '(i6,f7.1,f6.1,5f7.1,4f7.1,2f7.1,5f9.1,f12.1,3f7.1,2f7.2)') &
!       ktau,met%fsd(ijd),met%fld(ijd), &
!       canopy%fnv(ijd),canopy%fev(ijd),canopy%fevc(ijd),canopy%fevw(ijd),canopy%fhv(ijd), &   !43-47
!!       canopy%fns(ijd),canopy%fes(ijd),canopy%fhs(ijd),canopy%ga(ijd), &
!       canopy%fe(ijd),canopy%fh(ijd), &
!       mout%drybal(ijd),mout%wetbal(ijd), &
!       sum(rad%qcan(ijd,:,1))+sum(rad%qcan(ijd,:,2))+rad%qssabs(ijd) &
!       +met%fld(ijd) - sboltz*emsoil*rad%trad(ijd)**4    &
!       -canopy%fev(ijd) - canopy%fes(ijd)*ssoil%cls(ijd) -canopy%fh(ijd) -canopy%ga(ijd), &
!       met%fsd(ijd)*(1.-(rad%albedo(ijd,1)+rad%albedo(ijd,2))/2) &
!       +met%fld(ijd)-sboltz*emleaf*canopy%tv(ijd)**4*(1.-rad%transd(ijd))-mout%flws(ijd)*rad%transd(ijd) &
!       -canopy%fev(ijd) - canopy%fes(ijd)*ssoil%cls(ijd) -canopy%fh(ijd) -canopy%ga(ijd), &
!       mout%ebal(ijd),ssoil%totenbal(ijd), &
!       rad%trad(ijd),canopy%tv(ijd),ssoil%tss(ijd),veg%vlaiw(ijd),rad%transd(ijd)
!  ! New file-- - -- -
!  u = 101
!  INQUIRE (u, opened=is_open)
!  IF (.not. is_open) THEN
!     OPEN (u, file='f101.txt', status='replace')
!     WRITE(u, '(a6,23(1x,a12))') &
!          'ktau', 'precip', 'cansto', 'delwc', 'rnof1', 'rnof2', 'fe*dels/rlam', &
!          'fev*dels/rlam', 'fevc*dels/rlam', 'fevw*dels/rlam', 'fes/cls*dels/rlam', &
!          'owbtot', 'wbtot-owbtot', 'delwat', 'tprecip', 'tevap', 'trnoff', &
!          'wbtot', 'tdelwat', 'totenbal'
!  END IF
!  WRITE(u, '(i6,12f13.6,11f13.5)') &
!       ktau,met%precip(ijd),veg%cansto(ijd),canopy%delwc(ijd), &
!       ssoil%rnof1(ijd),ssoil%rnof2(ijd),canopy%fe(ijd)*dels/air%rlam(ijd), &
!       canopy%fev(ijd)*dels/air%rlam(ijd),canopy%fevc(ijd)*dels/air%rlam(ijd), &
!       canopy%fevw(ijd)*dels/air%rlam(ijd),  &
!       canopy%fes(ijd)/ssoil%cls(ijd)*dels/air%rlam(ijd),owbtot(ijd), &
!       ssoil%wbtot(ijd)-owbtot(ijd),ssoil%delwat(ijd), &
!       ssoil%tprecip(ijd), ssoil%tevap(ijd),ssoil%trnoff(ijd),  &
!       ssoil%wbtot(ijd),ssoil%tdelwat(ijd),ssoil%totenbal(ijd)
!  IF(ktau==ntau) PRINT *,'wbt res=', &
!       ssoil%tprecip(ijd)-ssoil%tevap(ijd)-ssoil%trnoff(ijd)+(ssoil%wbtot0(ijd)-ssoil%wbtot(ijd)) &
!       +ssoil%osnowd0(ijd)-ssoil%snowd(ijd), &
!       ssoil%tprecip(ijd), ssoil%tevap(ijd),ssoil%trnoff(ijd),ssoil%totenbal(ijd)
! New file-- --- --
  u = 98
  INQUIRE (u, opened=is_open)
  IF (.not. is_open) THEN
     OPEN (u, file='f98.txt', status='replace')
     WRITE(u, '(a6,60a14)') &
          'ktau', 'hod', 'fsd', 'fld', 'precip', & !1-5
          'tk', 'pmb', 'ua', 'qv', 'ca', & !6-11
          'swnet', 'lwnet', &
          'fe', 'fh', & !12-14
          'fnv+fns', 'trad', 'tv', & !15-17
          'acond', 'wb(1,1)', 'wb(1,2)', & !18-20
          'wb(1,3)', 'wb(1,4)', 'wb(1,5)', 'wb(1,6)', & !21-24
          'xx1', 'xx2', 'xx3', & !25-27
          'tgg(1,1)', 'tgg(1,2)', 'tgg(1,3)', 'tgg(1,4)', & !28-31
          'tgg(1,5)', 'tgg(1,6)', & !32-33
          'fns', 'fes', 'fhs', 'ga', & !34-37
          'frp', 'frpw', 'frpr', 'fpn', 'frs', & !38-42
          'pfrs', &
          'fnv', 'fev', 'fevc', 'fevw', 'fhv', & !43-47
          'precis', 'vlai', 'ebal', 'ebalcan', 'ebalcan2', 'ebalsoil', & !48-53
          'runoff', 'rnof1', 'rnof2', 'delwc', 'snowd', & !54-58
          'delwat', 'albedo(:,1)', 'albedo(:,2)' !59-61
  END IF
  xx1=ssoil%tggsn(:,1)
  xx2=ssoil%tggsn(:,2)
  xx3=ssoil%tggsn(:,3)
  where (ssoil%isflag == 0 )
     xx1=ssoil%tgg(:,1)
     xx2=ssoil%tgg(:,2)
     xx3=ssoil%tgg(:,3)
  end where
  do k=1,mp
    swnet(k)=SUM(rad%qcan(k,:,1))+SUM(rad%qcan(k,:,2))+rad%qssabs(k)
    lwnet(k)=met%fld(k)-sboltz*emleaf*canopy%tv(k)**4 *(1-rad%transd(k))-rad%flws(k)*rad%transd(k)
  enddo
!
  do k=1,7
   ijd =listp(k)
     WRITE(u, '(i5,"xr",i7,65e14.5)') &
          ijd,ktauyear,met%hod(ijd),met%fsd(ijd),met%fld(ijd),met%precip(ijd), & ! 1-5
          met%tk(ijd), met%pmb(ijd),met%ua(ijd),met%qv(ijd),met%ca(ijd), & ! 6-11
          swnet(ijd),lwnet(ijd), &
          canopy%fe(ijd),canopy%fh(ijd),&                                 !12-14
          canopy%fnv(ijd)+canopy%fns(ijd), rad%trad(ijd),canopy%tv(ijd), & !15-17
          1./max(0.0001,rough%rt1(ijd)), ssoil%wb(ijd,1),ssoil%wb(ijd,2),& !18-20
          ssoil%wb(ijd,3),ssoil%wb(ijd,4), ssoil%wb(ijd,5),ssoil%wb(ijd,6),& !21-24
          xx1(ijd),xx2(ijd),xx3(ijd),      &   !25-27
                    !               0.0, 0.0, 0.0,    & !25-27
          ssoil%tgg(ijd,1),ssoil%tgg(ijd,2), ssoil%tgg(ijd,3),ssoil%tgg(ijd,4),& !28-31
          ssoil%tgg(ijd,5),ssoil%tgg(ijd,6),    & !32-33
          canopy%fns(ijd),canopy%fes(ijd),canopy%fhs(ijd),canopy%ga(ijd), & !34-37
          canopy%frp(ijd),canopy%frpw(ijd),canopy%frpr(ijd),canopy%fpn(ijd),canopy%frs(ijd),& !38-42
          canopy%fnv(ijd),canopy%fev(ijd),canopy%fevc(ijd),canopy%fevw(ijd),canopy%fhv(ijd), & !43-47
          canopy%precis(ijd),veg%vlai(ijd),bal%ebal(ijd), &
          canopy%fnv(ijd)-canopy%fhv(ijd)-canopy%fev(ijd), &
          bal%drybal(ijd)+bal%wetbal(ijd), &
          canopy%fns(ijd)-canopy%fhs(ijd)-canopy%fes(ijd)*ssoil%cls(ijd)-canopy%ga(ijd), &
          ssoil%runoff(ijd),ssoil%rnof1(ijd),ssoil%rnof2(ijd),canopy%delwc(ijd),ssoil%snowd(ijd), &
          bal%wbal(ijd),rad%albedo(ijd,1),rad%albedo(ijd,2), &
         1./max(0.001,ssoil%rtsoil(ijd)), canopy%us(ijd),canopy%cduv(ijd), pblh(ijd)
         enddo

!!! these came from carbon.f90; the end of carbon_pl, before the        bgc%cplant(:,1) = max statements

!  u = 100
!  INQUIRE (u, opened=is_open)
!  IF (.not. is_open) THEN
!     OPEN (u, file='f100.txt', status='replace')
!     WRITE(u, *) 'clitt ',  'frs ', 'cplant(:,1)', 'cplant(:,2)','cplant(:,3)', &
!          'csoil(:,1) ', 'csoil(:,2)' , 'coef_cd'
!  END IF
!  !     WRITE(u, '(11e13.5)') clitt, cfrts, cfwd, cfsf, cfrts, canopy%frs, &
!!  !             bgc%csoil(:,1), bgc%csoil(:,2) ,coef_cd, coef_cdnew
!  WRITE(u, '(11e13.5)') clitt(ijd),canopy%frs(ijd),bgc%cplant(ijd,1),bgc%cplant(ijd,2),bgc%cplant(ijd,3), &
!       bgc%csoil(ijd,1), bgc%csoil(ijd,2) ,coef_cd(ijd), coef_cdnew(ijd)
!!!  !! and these from the very end of SUBROUTINE soilcarb:
!
!  u = 103
!  INQUIRE (u, opened=is_open)
!  IF (.not. is_open) THEN
!     OPEN (u, file='f103.txt', status='replace')
!     WRITE(u, *) 'frs '
!  END IF
!  WRITE(u, '(9e13.5)') canopy%frs(ijd)



  END SUBROUTINE eva_output
!*******************************************************************************
end module cable_ccam

