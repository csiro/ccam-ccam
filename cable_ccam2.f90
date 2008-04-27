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
!  character(len=80), save :: vegtypefile=''
!  character(len=80), save :: vegparmfile=''
  logical, save :: vegparmnew=.false. ! old or new format for veg parm file
  integer, save :: CO2forcingtype=1   ! 1 constant, 2 prescribed 1900-2004,
                                      ! 3 interactive
  integer, save :: initcarbpools=999  ! 1 initialise from veg_parm file values
!  namelist/cableinput/vegtypefile,vegparmfile,vegparmnew, &
!                    & CO2forcingtype,initcarbpools
  ! arrays for vegetation names and types
!  character(len=25), dimension(:), allocatable, save :: vegname
  character(len=10), dimension(:), allocatable, save :: vegtype
!  logical, dimension(:), allocatable, save :: forest
!  logical, dimension(:), allocatable, save :: deciduous
!  logical, dimension(:), allocatable, save :: shrub
!  logical, dimension(:), allocatable, save :: grass
!  logical, dimension(:), allocatable, save :: crop
!  logical, dimension(:), allocatable, save :: noveg  ! ice, possibly also barren

  contains
  ! ****************************************************************************

  subroutine sib4(nvegt)     ! new version of sib1 with soilsnowv
! BP added nvegt to be passed down to cbm (BP jan 2008)

      use cc_mpi
      use zenith_m
      USE define_types, cbm_ms => ms
      USE air_module
      USE canopy_module
      USE carbon_module
      USE cbm_module
      USE soil_snow_module
      TYPE (air_type)             :: air
      TYPE (bgc_pool_type)        :: bgc
      TYPE (canopy_type)          :: canopy
      TYPE (met_type)             :: met
      TYPE (balances_type)        :: bal
      TYPE (radiation_type)       :: rad
      TYPE (roughness_type)       :: rough
      type (soil_parameter_type)  :: soil       ! soil parameters
      TYPE (soil_snow_type)       :: ssoil
      TYPE (sum_flux_type)        :: sum_flux
      type (veg_parameter_type)   :: veg        ! vegetation parameters
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
      subroutine cbmrdn2 ! sib hardwired dummy version.  Use cbmrdn from cable_ccam.f90 when avaliable.
!
!     reads the parameters required by land surface scheme 
!     called from indata

      use cc_mpi, only : myid
      use define_dimensions, only : ncs, ncp

      implicit none

      include 'newmpar.h'
      include 'carbpools.h'
      include 'nsibd.h'
      include 'latlong.h'
      include 'soil.h'
      include 'soilsnow.h'
      include 'soilv.h'
      include 'vegpar.h'
      
      if (myid == 0) print *,"Setting CABLE defaults (sib)"
      
      vegparmnew=.false.

      allocate(vegtype(mxvt))

      where ((ivegt.eq.9).and.(rlatt.ge.-30.).and.(rlatt.le.30.))
        c4frac=0.95
      else where ((ivegt.eq.8).and.(rlatt.ge.-30.).and.(rlatt.le.0.))
        c4frac=0.5
      else where ((ivegt.eq.8).and.(rlatt.ge.0.).and.(rlatt.le.20.))
        c4frac=0.8
      else where ((ivegt.eq.8).and.(rlatt.ge.20.).and.(rlatt.le.30.))
        c4frac=0.5
      else where (ivegt.eq.6)
        c4frac=0.75
      else where ((ivegt.eq.7).and.(rlatt.ge.-30.).and.(rlatt.le.-20.))
        c4frac=0.5
      else where ((ivegt.eq.7).and.(rlatt.ge.-20.).and.(rlatt.le.20.))
        c4frac=0.95
      else where ((ivegt.eq.7).and.(rlatt.ge.20.).and.(rlatt.le.35.))
        c4frac=0.5
      else where ((ivegt.eq.12).and.(rlatt.ge.0.).and.(rlatt.le.40.))
        c4frac=0.3
      else where
        c4frac=0.
      end where
      frac4=0.
      hc=(/ 35.,20.,20.,17.,17.,1.,0.5,0.6,0.5,4.,0.05,1.,0.01,0.01,0.01,0.01,0.01 /)
      xfang=(/ 0.1,0.25,0.125,0.01,0.01,-0.3,-0.3,0.2,0.01,0.2,0.01,-0.3,0.01,0.01,0.01,0.01,0.01 /)
      dleaf=(/ 0.075,0.12,0.07,0.028,0.02,0.16,0.16,0.155,0.0165,0.155,0.005,0.155,0.005,0.005,0.005,0.005,0.005 /)
      wai=(/ 1.6,1.2,1.,1.,0.6,0.5,0.,0.5,0.5,0.5,0.,0.,0.,0.,0.,0.,0. /)
      canst1=0.1
      shelrb=2.
      vegcf=(/ 1.95,1.5,1.55,0.91,0.73,2.8,2.75,0.,2.05,0.6,0.4,2.8,0.,0.,0.,0.,0. /)
      vegtype(:)='others'
      vegtype(2)='deciduous'
      vegtype(5)='deciduous'
      extkn=0.4
      vcmax=(/ 85.E-6,85.E-6,80.E-6,65.2E-6,70.E-6,10.1E-6,9.E-6,35.053E-6,39.E-6,38.053E-6,17.02E-6,120.E-6,1.E-6,1.E-6, &
               1.E-6,1.E-6,1.E-6 /)
      ejmax=2.*vcmax
      rp20=(/ 1.1342,1.4733,2.3704,3.3039,3.2879,1.0538,0.8037,12.032,1.2994,7.4175,2.8879,3.,0.1,0.1,0.1,0.1,0.1 /)
      rpcoef=0.0832
      rs20=(/ 1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,1.,0. /) ! ?????? wrong size
      tminvj=(/ 5.,5.,5.,2.,5.,10.,10.,0.,5.,-5.,5.,10.,-5.,-5.,-5.,-5.,-5. /)
      tmaxvj=(/ 15.,15.,10.,5.,10.,15.,15.,5.,10.,0.,10.,15.,0.,0.,0.,0.,0. /)
      vbeta=(/ 1.,1.,1.,1.,1.,4.,4.,4.,4.,1.,4.,4.,1.,1.,1.,1.,1. /)
      rootbeta=0.
      tcplant(:,1)=(/ 300.,300.,200.,200.,200.,250.,250.,250.,150.,150.,1.,150.,1.,1.,1.,1.,1. /)
      tcplant(:,2)=(/ 16833.,12000.,10217.,10217.,5967.,5247.,0.,5247.,5000.,5000.,0.,0.,0.,0.,0.,0.,0. /)
      tcplant(:,3)=(/ 1443.,1029.,876.,876.,511.,1124.,500.,1124.,500.,500.,1.,607.,1.,1.,1.,1.,1. /)
      tcsoil(:,1)=(/ 303.,216.,184.,184.,107.,275.,275.,275.,100.,100.,1.,149.,1.,1.,1.,1.,1. /)
      tcsoil(:,2)=(/ 606.,432.,367.,367.,214.,314.,314.,314.,250.,250.,1.,300.,1.,1.,1.,1.,1. /)
      ratecp(1)=1.
      ratecp(2)=0.03
      ratecp(3)=0.14
      ratecs(1)=2.
      ratecs(2)=5.

      if (all(cplant.eq.0.)) then
        if (myid == 0) print *,"Using default carbpools"
        where (land)
          cplant(:,1)=tcplant(ivegt(:),1)
          cplant(:,2)=tcplant(ivegt(:),2)
          cplant(:,3)=tcplant(ivegt(:),3)
          csoil(:,1)=tcsoil(ivegt(:),1)
          csoil(:,2)=tcsoil(ivegt(:),2)
          cansto(:)=0.
        end where
      else
        if (myid == 0) print *,"Loading carbpools from ifile"
      end if

      albsoil(:)=0.5*sum(albvisnir(:,:),2)
      where ((isoilm(:).eq.9).and.land)
        albsoil(:)=min(albsoil(:),0.61)
      else where ((isoilm(:).eq.1).and.land)
        albsoil(:)=min(albsoil(:),0.41)
      else where (land)
        albsoil(:)=min(albsoil(:),0.25)
      end where  
      if (all(albsoilsn.eq.0.)) then
        if (myid == 0) print *,"Using default albsoilsn"
        where (albsoil(:).le.0.14)
          albsoilsn(:,2)=2.*albsoil(:)/1.5
          albsoilsn(:,1)=0.5*albsoilsn(:,2)
        else where (albsoil(:).le.0.21)
          albsoilsn(:,2)=2.*albsoil(:)/1.62
          albsoilsn(:,1)=0.62*albsoilsn(:,2)
        else where
          albsoilsn(:,2)=2.*albsoil(:)/1.68
          albsoilsn(:,1)=0.68*albsoilsn(:,2)
        end where
      else
        if (myid == 0) print *,"Loading albsoilsn from ifile"
      end if
 
      end subroutine cbmrdn2

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
      real, dimension(mxvt,5) :: froot2

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
        veg%wai    =    wai(veg%iveg)
        veg%vegcf  =  vegcf(veg%iveg)
        veg%extkn  =  extkn(veg%iveg)
        veg%deciduous =(vegtype(veg%iveg).eq.'deciduous') ! MJT suggestion

!       soil parameters
        !soil%zse  = zse
        !soil%zshh = zshh
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
         IF (.not.vegparmnew) THEN   ! CASA vegetation types  
          froot2(:,1) = (/.02,.04,.04,.04,.04,.05,.05,.05,.05,.05,.05,.05,.05,.01,.01,.01,.01/)
          froot2(:,2) = (/.06,.11,.11,.11,.11,.15,.15,.10,.10,.10,.10,.15,.15,.01,.01,.01,.01/)
          froot2(:,3) = (/.14,.20,.20,.20,.20,.35,.35,.35,.35,.35,.35,.34,.35,.01,.01,.01,.01/)
          froot2(:,4) = (/.28,.26,.26,.26,.26,.39,.39,.35,.35,.35,.35,.38,.40,.01,.01,.01,.01/)
          froot2(:,5) = (/.35,.24,.24,.24,.24,.05,.05,.10,.10,.10,.10,.06,.04,.01,.01,.01,.01/)
          froot2(:,6) = (/.15,.15,.15,.15,.15,.01,.01,.05,.05,.05,.05,.02,.01,.01,.01,.01,.01/)
         ELSE
          ! IGBP vegetation types without/with water bodies
          froot2(:,1) = (/.04,.02,.04,.04,.04,.05,.05,.05,.05,.05,.05,.05,.05,.05,.05,.05,.01/)
          froot2(:,2) = (/.11,.06,.11,.11,.11,.10,.10,.15,.15,.15,.10,.15,.10,.15,.15,.10,.01/)
          froot2(:,3) = (/.20,.14,.20,.20,.20,.35,.35,.35,.35,.35,.35,.34,.35,.34,.35,.35,.01/)
          froot2(:,4) = (/.26,.28,.26,.26,.26,.35,.35,.39,.39,.39,.35,.38,.35,.38,.40,.35,.01/)
          froot2(:,5) = (/.24,.35,.24,.24,.24,.10,.10,.05,.05,.05,.10,.06,.10,.06,.04,.10,.01/)
          froot2(:,6) = (/.15,.15,.15,.15,.15,.05,.05,.01,.01,.01,.05,.02,.05,.02,.01,.05,.01/)
         END IF
        else
!        preferred option
!        froot is now calculated from soil depth and the new parameter rootbeta 
!        according to Jackson et al. 1996, Oceologica, 108:389-411
         totroot(:) = 1.0-rootbeta(:)**(sum(soil%zse)*100.0)
         totdepth = 0.0
         do k=1,ms
           totdepth = totdepth + soil%zse(k)*100.0
           froot2(:,k) = min(1.0,1.0-rootbeta(:)**totdepth)
         enddo
         do k = ms, 2, -1
           froot2(:,k) = froot2(:,k) - froot2(:,k-1)
         enddo
        endif
        veg%froot = froot2(veg%iveg,:)

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
          met%fsd(ip) = sgsave(iq)/(1.-albvisnir(iq,1))! short wave down (positive) W/m^2
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
      tscrn = unpack(canopy%tscrn,  land, tscrn)    ! clobbered by scrnout?
      qgscrn = unpack(canopy%qscrn,  land, qgscrn)  ! clobbered by scrnout?
      uscrn = unpack(canopy%uscrn,  land, uscrn)    ! clobbered by scrnout?
      cduv= unpack(canopy%cduv, land, cduv)         ! clobbered by sflux?
      cansto= unpack(canopy%cansto, land, cansto)
      vlai= unpack(veg%vlai, land, vlai)
      gflux = unpack(canopy%ghflux, land, gflux)
      sgflux = unpack(canopy%sghflux, land, sgflux)
      !rtsoil = unpack(ssoil%rtsoil, land, rtsoil)
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

  END Subroutine cbm_unpack

end module cable_ccam

