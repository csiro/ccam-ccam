 module pack_unpack_m
  public :: cbm_pack       
  contains

  subroutine cbm_pack(air, bgc, canopy, met, bal, rad,  &
               rough, soil, ssoil, sum_flux, veg)

  use define_types, cbm_ms => ms
  use roughness_module
  use air_module
  use radiation_module
  use canopy_module
!  use output_module 
  use soil_snow_module             
  implicit none
  type (air_type)                       :: air 
  type (canopy_type)                    :: canopy     ! updated each timestep
  type (balances_type)            :: bal        ! energy/water bal variables
!  type(misc_output_type)           :: mout       ! non-type output variables
!  type (gcmmet_type)             :: met ! met forcing from gcm (timestep at a time)
  type (met_type)             :: met ! met forcing from gcm (timestep at a time)
!  type (soil_type)               :: soil        ! various data for land point
  type (soil_parameter_type)      :: soil       ! soil parameters
  type (soil_snow_type)           :: ssoil        ! various data for land point
  type (sum_flux_type)            :: sum_flux
  type (radiation_type)                 :: rad        ! updated each timestep
  type (roughness_type)                 :: rough
!  type (veg_type)              :: veg        ! updated each timestep
  type (veg_parameter_type)       :: veg        ! vegetation parameters

!  type (ssoil_type), save               :: ssoil      ! updated each timestep
  type (bgc_pool_type)            :: bgc        ! bgc parameters
  INTEGER(i_d) :: ip,iq,j,k  
  INTEGER(i_d) ::  jyear, ipco2,jmonth,jday,jhour,jmin,mstart,ndoy,mins 
  dimension ndoy(12)   ! days from beginning of year (1st Jan is 0)
  data ndoy/ 0,31,59,90,120,151,181,212,243,273,304,334/

!      include 'cbmdim.h'
       include 'newmpar.h' 
       include 'aalat.h'
       include 'arrays.h'
       include 'dates.h'
!       include 'histave.h'
       include 'morepbl.h' 
       include 'nsibd.h' 
       include 'carbpools.h'
       include 'permsurf.h'  
       include 'parm.h'  
       include 'soilbal.h'
       include 'soilsnow.h'
!       include 'soilroot.h'
       include 'soilv.h'
       include 'vegpar.h'
       include 'soil.h'
       REAL(r_1), DIMENSION(mp) :: hour   ! 

!      define soil parameters
       soil%zse = (/.022, .058, .154, .409, 1.085, 2.872/) ! soil layer thickness
       soil%zshh(1) = 0.5 * soil%zse(1)
       soil%zshh(ms+1) = 0.5 * soil%zse(ms)
       soil%zshh(2:ms) = 0.5 * (soil%zse(1:ms-1) + soil%zse(2:ms))
       froot(:,1) = (/.02,.04,.04,.04,.04,.05,.05,.05,.05,.05,.05,.05,.05/)
       froot(:,2) = (/.06,.11,.11,.11,.11,.15,.15,.10,.10,.10,.10,.15,.15/)
       froot(:,3) = (/.14,.20,.20,.20,.20,.35,.35,.35,.35,.35,.35,.34,.35/)
       froot(:,4) = (/.28,.26,.26,.26,.26,.39,.39,.35,.35,.35,.35,.38,.40/)
       froot(:,5) = (/.35,.24,.24,.24,.24,.05,.05,.10,.10,.10,.10,.06,.04/)
       froot(:,6) = (/.15,.15,.15,.15,.15,.01,.01,.05,.05,.05,.05,.02,.01/)


       do ip=1,ipland  ! all land points
       iq=iperm(ip)
!         rml 21/09/07 changed soil%froot to veg%froot to match changes 
!         in cable_offline
          veg%froot(ip,1)= froot(ivegt(iq),1)
          veg%froot(ip,2)= froot(ivegt(iq),2)
          veg%froot(ip,3)= froot(ivegt(iq),3)
          veg%froot(ip,4)= froot(ivegt(iq),4)
          veg%froot(ip,5)= froot(ivegt(iq),5)
          veg%froot(ip,6)= froot(ivegt(iq),6)
          soil%bch(ip) = bch(isoilm(iq))
!         rml 21/09/07 c3 no longer needed in updated soilsnow
!         soil%c3(ip)= c3(isoilm(iq))
          soil%cnsd(ip) = cnsd(isoilm(iq))
          soil%css(ip) = css(isoilm(iq))
          soil%rhosoil(ip) = rhos(isoilm(iq))
          soil%hyds(ip) = hyds(isoilm(iq))
          soil%sucs(ip) = sucs(isoilm(iq))
          soil%hsbh(ip)  = hyds(isoilm(iq))*abs(sucs(isoilm(iq)))*bch(isoilm(iq))
          soil%ibp2(ip) = ibp2(isoilm(iq))
          soil%i2bp3(ip) =i2bp3(isoilm(iq))
          soil%sfc(ip) = sfc(isoilm(iq))
          soil%ssat(ip) = ssat(isoilm(iq))
          soil%swilt(ip) = swilt(isoilm(iq))
          soil%rs20(ip) = rs20(ivegt(iq))
       enddo !ip
       soil%albsoil = pack(albsoil, land) 
!    print *,'pack albsoil',soil%albsoil(2468),albsoil(7560)
       soil%isoilm = pack(isoilm, land)
!    print *,'pack',soil%sfc(soil%isoilm(2)),sfc(2)
!    print *,'pack',soil%rhosoil(soil%isoilm(2)),rhos(2)
       rad%latitude = pack(alat, land)
!!    soil%silt = pack(silt, land) 
!!    soil%clay = pack(clay, land)
!!    soil%sand = pack(sand, land)
       ssoil%albsoilsn(:,1) = pack(albsoilsn(:,1), land)
       ssoil%albsoilsn(:,2) = pack(albsoilsn(:,2), land)
       ssoil%albsoilsn(:,3) = 0.05
!     print *, 'pack7560',soil%albsoilsn(2468,1),soil%albsoilsn(2468,2),albsoilsn(7560,1),albsoilsn(7560,2)
!    soil%soildepth = pack(soildepth, land)
!    print *, 'pack',soil%albsoilsn,albvisnir
!      define soil variables 
!       ssoil%runoff = pack(runoff, land)
!       ssoil%rnof1 = pack(rnof1, land)
!       ssoil%rnof2 = pack(rnof2, land)
!       ssoil%wbtot0 = pack(wbtot0, land)
!       bal%tevap = pack(tevap, land)
!       bal%tprecip = pack(tprecip, land)
!       bal%totenbal = pack(totenbal, land)
!       bal%trnoff = pack(trnoff, land)
       bal%wbtot0 = pack(wbtot0, land)
       bal%evap_tot = pack(tevap, land)
       bal%precip_tot = pack(tprecip, land)
       bal%ebal_tot = pack(totenbal, land)
       bal%rnoff_tot = pack(trnoff, land)


       do k = 1,6
        ssoil%tgg(:,k) = pack(tgg(:,k), land)
        ssoil%wb(:,k) = pack(wb(:,k), land)
        ssoil%wbice(:,k) = pack(wbice(:,k), land)
        soil%zse(k) = zse(k)
       enddo
       if(ktau == 1) then
        bal%wbtot0 = 0.0
        ssoil%wbtot = 0.
       DO j=1,ms
        bal%wbtot0 = bal%wbtot0 + ssoil%wb(:,j) * soil%zse(j) * 1000.0
        ssoil%wbtot  = ssoil%wbtot + ssoil%wb(:,j) * soil%zse(j) * 1000.0
       END DO
       endif
       ssoil%wbtot = pack(wbtot, land)
!       ssoil%froot = froot
       ssoil%isflag = pack(isflag, land)
       veg%iveg = pack(ivegt, land)
!     print *, 'pack7560',veg%iveg(2468),ivegt(7560),veg%iveg(3688),ivegt(13419)
!    soil%soildepth = pack(soildepth, land)
!       do k = 1,6
!       ssoil%froot(:,k) = froot(veg%iveg(:),k)
!       enddo
!    print *,'cbm froot0',ivegt(13419),ssoil%froot(veg%iveg(2468),:)
!    print *,'cbm zse',zse
!    soil%zshh(6) = zshh(6)
       do k = 1,3
        ssoil%tggsn(:,k) = pack(tggsn(:,k), land)
        ssoil%smass(:,k) = pack(smass(:,k), land)
        ssoil%ssdn(:,k) = pack(ssdn(:,k), land)
       enddo
!    print *,'cbm tgg',ivegt(13419),ssoil%tgg(3688,:)
!    print *,'cbm tgg',ivegt(7560),ssoil%tgg(2468,:)
!    print *,'cbm tggsn',ivegt(7560),ssoil%tggsn(2468,:)

       ssoil%ssdnn = pack(ssdnn, land)
       ssoil%snowd = pack(snowd, land)
       ssoil%osnowd = pack(osnowd, land)
       bal%osnowd0 = pack(osnowd0, land)
       ssoil%snage = pack(snage, land)
       ssoil%isflag = pack(isflag, land)
!    soil%otgsoil = pack(otgsoil, land)
       ssoil%rtsoil = pack(rtsoil, land)

!    print *,'cbm snow',ivegt(7560),ssoil%snowd(2468),ssoil%isflag(2468)
!        canopy%ga = pack(oga, land)
        canopy%ghflux = pack(gflux, land)
        canopy%sghflux = pack(sgflux, land)
!    canopy%fhs = pack(fhs, land)
!    canopy%dgdtg = pack(dgdtg, land)
!    canopy%fevc = pack(fev, land)
!    canopy%fes = pack(fes, land)
!    canopy%precis = pack(condxpr, land)

!         define soil parameters
       do ip=1,ipland  ! all land points
       iq=iperm(ip)
         veg%canst1(ip) = canst1(ivegt(iq))
         veg%ejmax(ip) = ejmax(ivegt(iq))
         veg%frac4(ip) = frac4(ivegt(iq))
         veg%tminvj(ip) = tminvj(ivegt(iq))
         veg%tmaxvj(ip) = tmaxvj(ivegt(iq))
         veg%vbeta(ip) = vbeta(ivegt(iq))
         veg%hc(ip) = hc(ivegt(iq))
         veg%rp20(ip) = rp20(ivegt(iq))
         veg%rpcoef(ip) = rpcoef(ivegt(iq))
         veg%shelrb(ip) = shelrb(ivegt(iq))
         veg%vcmax(ip) = vcmax(ivegt(iq))
         veg%xfang(ip) = xfang(ivegt(iq))
         veg%dleaf(ip) = dleaf(ivegt(iq))
       enddo ! ip
       canopy%cansto = pack(cansto, land)
       veg%vlai = pack(vlai, land)
       veg%vlaimax = pack(rlaimax, land)

!      rml 28/09/07 c4 fraction defined for each grid point, for
!      compatibility just overwrite values from parameter file if
!      find data in c4frac array
       if (sum(c4frac).gt.0) veg%frac4 = pack(c4frac,land)
 
!       print *,' pack_cbm vlai',veg%vlai(3688),vlai(13419),veg%vlai(2468), &
!        veg%vlaimax(2468)
!    print *,'cbm veg',ivegt(13419),veg%vlai(3688),veg%vcmax(ivegt(13419)), &
!    veg%vlaimax(3688),veg%dleaf(ivegt(13419))
   
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
!       print *,'packbgc',bgc%cplant(3688,:),bgc%csoil(3688,:)
!       print *,'packbgc',bgc%cplant(2468,:),bgc%csoil(2468,:)
!      define met data
!       soil%latitude = pack(alat, land)
!       print *,'cbm_pack 2',soil%latitude(3688),soil%latitude(2468)
!       soil%longitude = pack(along, land)
! need a better way to calculate local day and hour
!       met%doy = 

        
      end subroutine cbm_pack

  subroutine cbm_unpack(ktauyear, air, bgc, canopy, met, bal, rad,  &
               rough, soil, ssoil, sum_flux, veg)

!    use define_types
!    use define_dimensions
!    use roughness_module
!    use air_module
!    use radiation_module
!    use canopy_module
!!    use output_module
!    use physical_constants
!    use soil_snow_module
!    use cbm_data
!!    implicit none
  use define_types, cbm_ms => ms
  use roughness_module
  use air_module
  use radiation_module
  use canopy_module
!  use output_module 
  use soil_snow_module
!  implicit none
  INTEGER(i_d), INTENT(IN)              :: ktauyear
  type (air_type)                       :: air
  type (canopy_type)                    :: canopy     ! updated each timestep
  type (balances_type)            :: bal        ! energy/water bal variables
!  type(misc_output_type)           :: mout       ! non-type output variables
!  type (gcmmet_type)             :: met ! met forcing from gcm (timestep at a time)
  type (met_type)             :: met ! met forcing from gcm (timestep at a time)
!  type (soil_type)               :: soil        ! various data for land point
  type (soil_parameter_type)      :: soil       ! soil parameters
  type (soil_snow_type)           :: ssoil        ! various data for land point
  type (sum_flux_type)            :: sum_flux
  type (radiation_type)                 :: rad        ! updated each timestep
  type (roughness_type)                 :: rough
!  type (veg_type)              :: veg        ! updated each timestep
  type (veg_parameter_type)       :: veg        ! vegetation parameters
!  type (ssoil_type), save               :: ssoil      ! updated each timestep
  type (bgc_pool_type)            :: bgc        ! bgc parameters
  INTEGER(i_d) :: ip,iq,j,k
  INTEGER(i_d) ::  jyear, ipco2,jmonth,jday,jhour,jmin,mstart,ndoy,mins
  dimension ndoy(12)   ! days from beginning of year (1st Jan is 0)
  data ndoy/ 0,31,59,90,120,151,181,212,243,273,304,334/


  include 'newmpar.h'
!    include 'cbmdim.h'
       include 'aalat.h'
       include 'arrays.h'
       include 'dates.h'
!       include 'histave.h'
       include 'morepbl.h'
       include 'nsibd.h'
       include 'carbpools.h'
       include 'pbl.h'    ! new cduv for vetmix
       include 'permsurf.h'
       include 'soilbal.h'
       include 'soilsnow.h'
!       include 'soilroot.h'
       include 'soilv.h'
       include 'screen.h'
       include 'vegpar.h'
       include 'soil.h'

    albsoilsn(:,1) = unpack(ssoil%albsoilsn(:,1), land, albsoilsn(:,1))
    albsoilsn(:,2) = unpack(ssoil%albsoilsn(:,2), land, albsoilsn(:,2))
    albvisnir(:,1) = unpack(rad%albedo(:,1), land, albvisnir(:,1))
    albvisnir(:,2) = unpack(rad%albedo(:,2), land, albvisnir(:,2))
!    print *,'unpack',albvisnir(:,1),albvisnir(:,2)
!    print *,'unpack',soil%refs(:,1),soil%refs(:,2)
!    print *,'unpack',rad%albedo,ssoil%snowd,veg%vlai,veg%vlaiw
    runoff= unpack(ssoil%runoff, land, runoff)
    rnof1= unpack(ssoil%rnof1, land, rnof1)
    rnof2= unpack(ssoil%rnof2, land, rnof2)
    wbtot= unpack(ssoil%wbtot, land, wbtot)
!   owbtot= unpack(ssoil%owbtot, land, owbtot)
!   wbtot0 = unpack(ssoil%wbtot0, land, wbtot0)
!    tevap = unpack(ssoil%tevap, land, tevap)
    tevap = unpack(bal%evap_tot, land, tevap)
!    tprecip = unpack(ssoil%tprecip, land, tprecip)
    tprecip = unpack(bal%precip_tot, land, tprecip)
!    totenbal = unpack(ssoil%totenbal, land, totenbal)
    totenbal = unpack(bal%ebal_tot, land, totenbal)
!    trnoff = unpack(ssoil%trnoff, land, trnoff)
    trnoff = unpack(bal%rnoff_tot, land, trnoff)
    do k = 1,6
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
!    print *,'cbm_unpack',rnet(13419),fg(13419),eg(13419),tss(13419),   & 
!       tscrn(13419),cduv(13419),vlai(13419)
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
!    dsumpn= unpack(sum_flux%dsumpn, land, dsumpn)
!    dsumrp= unpack(sum_flux%dsumrp, land, dsumrp)
!    dsumrs= unpack(sum_flux%dsumrs, land, dsumrs)
!    dsumrd= unpack(sum_flux%dsumrd, land, dsumrd)
    do k=1,3
    cplant(:,k)= unpack(bgc%cplant(:,k), land, cplant(:,k))
    ratecp(k)= bgc%ratecp(k)
    enddo
    do k=1,2
    csoil(:,k)= unpack(bgc%csoil(:,k), land, csoil(:,k))
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
!!$SUBROUTINE eva_ebal
!!$  REAL(r_1), DIMENSION(mp)           :: swnet
!!$  REAL(r_1), DIMENSION(mp)           :: lwnet
!!$  ! calculate energy balances:
!!$  ! sw absorbed + lw absorbed - (lh+sh+ghflux) should = 0
!!$  do k=1,mp
!!$     mout%ebal(k) = sum(rad%qcan(k,:,1))+sum(rad%qcan(k,:,2))+rad%qssabs(k) &
!!$          +met%fld(k)-sboltz*emleaf*canopy%tv(k)**4*(1-rad%transd(k))-mout%flws(k)*rad%transd(k) &
!!$          -canopy%fev(k) - canopy%fes(k)*ssoil%cls(k) -canopy%fh(k) -canopy%ga(k)
!!$     canopy%fe = canopy%fev + canopy%fes
!!$     swnet(k)=SUM(rad%qcan(k,:,1))+SUM(rad%qcan(k,:,2))+rad%qssabs(k)
!!$     lwnet(k)=met%fld(k)-sboltz*emleaf*canopy%tv(k)**4 *(1-rad%transd(k))-mout%flws(k) *rad%transd(k)
!!$  enddo
!!$  IF (ktau==1) THEN
!!$     ssoil%totenbal = 0.
!!$     mout%drybal = 0.
!!$  END IF
!!$  ssoil%totenbal=ssoil%totenbal+mout%ebal
!!$  IF (ktau==0) THEN
!!$     owbtot=ssoil%wbtot0
!!$     ssoil%tdelwat=0.
!!$     ssoil%tprecip=0.
!!$     ssoil%tevap=0.
!!$     ssoil%trnoff=0.
!!$     ssoil%osnowd0=ssoil%osnowd
!!$  END IF
!!$  ssoil%delwat=met%precip - canopy%delwc - ssoil%snowd+ssoil%osnowd  &
!!$       -ssoil%rnof1-ssoil%rnof2-(canopy%fev+canopy%fes/ssoil%cls)*dels/air%rlam &
!!$       - ssoil%wbtot + owbtot         ! soil water change
!!$  ssoil%tdelwat=ssoil%tdelwat+ssoil%delwat
!!$  ssoil%tprecip=ssoil%tprecip+met%precip
!!$  ssoil%tevap=ssoil%tevap+(canopy%fev+canopy%fes/ssoil%cls)*dels/air%rlam
!!$  ssoil%trnoff=ssoil%trnoff+ssoil%rnof1+ssoil%rnof2
!!$
!!$
!!$END SUBROUTINE eva_ebal


END MODULE pack_unpack_m
