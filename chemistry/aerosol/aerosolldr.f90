! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2025 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
! This file is part of the Conformal Cubic Atmospheric Model (CCAM)
!
! CCAM is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! CCAM is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with CCAM.  If not, see <http://www.gnu.org/licenses/>.

!------------------------------------------------------------------------------

! This is a prognostic aerosol model for CCAM based on the LDR scheme used in Mk3.6 (Rotstayn and Lohmann 2002)
    
! Plan to replace diagnosed sea-salt with prognostic version based on CTM.  Eventually employ GLOMAP for modes.

module aerosolldr

use aerosol_arrays, only : ndust, nsalt, naero, ndcls,     &
  itracdms, itracso2, itracso4, itracbc, itracoc, itracdu, &
  itracsa
  

implicit none

private
public aldrcalc
public enhanceu10, Ch_dust, zvolcemi
public dustreff
public so4mtn, carbmtn, saltsmallmtn, saltlargemtn

! options
integer, save :: enhanceu10 = 0                 ! Modify 10m wind speed for emissions (0=none, 1=quadrature, 2=linear)
real, parameter :: zmin     = 1.e-10            ! Minimum concentration tolerance

! physical constants
real, parameter :: grav      = 9.80616          ! Gravitation constant
real, parameter :: rdry      = 287.04           ! Specific gas const for dry air
real, parameter :: cp        = 1004.64          ! Heat capacity of air
real, parameter :: hl        = 2.5104e6         ! Latent heat of vaporisation
real, parameter :: vkar      = 0.4              ! von Karman constant
real, parameter :: rhos      = 100.             ! Assumed density of snow in kg/m^3

! emission and deposition constants
real, save :: zvolcemi       = 8.               ! Total emission from volcanoes (TgS/yr)
real, save :: Ch_dust        = 1.e-9            ! Transfer coeff for type natural source (kg*s2/m5)

! Indirect effect coefficients
! converts from mass (kg/m3) to number concentration (1/m3) for dist'n
real, save :: so4mtn = 1.24e17                  ! Penner et al (1998)
real, save :: carbmtn = 1.25e17                 ! Penner et al (1998)
real, save :: saltsmallmtn = 1.89e16            ! Nillson et al. number mode radius = 0.1 um, sd=2
real, save :: saltlargemtn = 1.1e14             ! Nillson et al. number mode radius = 0.5 um, sd=2
!real, save :: so4mtn = 1.69e17                 ! IPCC (2001) Table 5.1
!real, save :: carbmtn = 1.21e17                ! IPCC (2001) Table 5.1
!real, save :: carbmtn = 2.30e17                ! counts Aitken mode as well as accumulation mode carb aerosols
!real, save :: saltsmallmtn = 3.79e17           ! Herzog number mode radius = 0.035 um, sd=1.92, rho=2.165 g/cm3
!real, save :: saltlargemtn = 7.25e14           ! Herzog number mode radius = 0.35 um, sd=1.7, rho=2.165

! emission indices
integer, parameter :: iso2a1 = 1      ! SO2/SO4 Anthropogenic surface
integer, parameter :: iso2a2 = 2      ! SO2/SO4 Anthropogenic upper level
integer, parameter :: ibca1  = 3      ! BC Anthropogenic surface
integer, parameter :: ibca2  = 4      ! BC Anthropogenic upper level
integer, parameter :: ioca1  = 5      ! OC Anthropogenic surface
integer, parameter :: ioca2  = 6      ! OC Anthropogenic upper level
integer, parameter :: iso2b1 = 7      ! SO2/SO4 BiomassBurning surface
integer, parameter :: iso2b2 = 8      ! SO2/SO4 BiomassBurning upper level
integer, parameter :: ibcb1  = 9      ! BC BiomassBurning surface
integer, parameter :: ibcb2  = 10     ! BC BiomassBurning upper level
integer, parameter :: iocb1  = 11     ! OC BiomassBurning surface
integer, parameter :: iocb2  = 12     ! OC BiomassBurning upper level
integer, parameter :: idmso  = 13     ! DMS ocean
integer, parameter :: idmst  = 14     ! DMS land
integer, parameter :: iocna  = 15     ! Natural organic

! Dust coefficients
real, dimension(ndust), parameter :: dustden = (/ 2500., 2650., 2650., 2650. /)    ! Density of dust (kg/m3)
                                                                                   ! (Clay, small silt, small slit, small silt)
real, dimension(ndust), parameter :: dustreff = (/ 0.73e-6,1.4e-6,2.4e-6,4.5e-6 /) ! Main effective radius (m)
                                                                                   ! (Clay, small silt, small slit, small silt)

! Salt coefficients
real, dimension(nsalt), parameter :: saltden = (/ 2165., 2165. /)     ! density of salt
real, dimension(nsalt), parameter :: saltreff = (/ 0.1e-6, 0.5e-6 /)  ! radius of salt (um)


contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Main routine

subroutine aldrcalc(dt,sig,dz,wg,pblh,prf,ts,ttg,condc,snowd,taudar,fg,eg,v10m,                    &
                    ustar,zo,land,fracice,tsigmf,qvg,qlg,qfg,stratcloud,clcon,cldcon,pccw,rhoa,vt, &
                    pfprec,pfmelt,pfsnow,pfsubl,plambs,pmrate,pmaccr,                              &
                    pqfsedice,prscav,prfreeze,pfevap,zdayfac,kbsav,locean)

use aerosol_arrays, only : xtg, xtosav, duste, dustdd, dustwd, dust_burden,          &
  bce, bcdd, bcwd, bc_burden, oce, ocdd, ocwd, oc_burden, vso2,                      &
  dmse, so2e, so4e, so2dd, so4dd, so2wd, so4wd, dms_burden, so2_burden, so4_burden,  &
  salte, saltdd, saltwd, salt_burden, dmsso2o, so2so4o, erod, zoxidant_g, emissfield
use newmpar_m

integer, dimension(:), intent(in) :: kbsav  ! Bottom of convective cloud
real, intent(in) :: dt                      ! Time step
real, dimension(:), intent(in) :: sig       ! Sigma levels
real, dimension(:), intent(in) :: wg        ! Soil moisture fraction of field capacity
real, dimension(:), intent(in) :: prf       ! Surface pressure
real, dimension(:), intent(in) :: ts        ! Surface temperture
real, dimension(:), intent(in) :: pblh      ! Boundary layer height
real, dimension(:), intent(in) :: v10m      ! 10m wind speed
real, dimension(:), intent(in) :: condc     ! Convective rainfall
real, dimension(:), intent(in) :: snowd     ! Snow depth
real, dimension(:), intent(in) :: taudar    ! Fraction of time sunlit
real, dimension(:), intent(in) :: fg        ! Sensible heat flux
real, dimension(:), intent(in) :: eg        ! Latent heat flux
real, dimension(:), intent(in) :: ustar     ! Friction velocity
real, dimension(:), intent(in) :: zo        ! Roughness length
real, dimension(:), intent(in) :: fracice   ! Sea-ice fraction
real, dimension(:), intent(in) :: tsigmf    ! Vegetation fraction
real, dimension(:), intent(in) :: vt        ! transfer velocity
real, dimension(:), intent(in) :: zdayfac   ! scale factor for day length
real, dimension(:), intent(in) :: cldcon    ! Convective rainfall area fraction
real, dimension(:,:), intent(in) :: dz     ! thickness of vertical levels (m)
real, dimension(:,:), intent(in) :: ttg    ! Air temperature
real, dimension(:,:), intent(in) :: qvg    ! liquid water mixing ratio
real, dimension(:,:), intent(in) :: qlg    ! liquid water mixing ratio
real, dimension(:,:), intent(in) :: qfg    ! frozen water mixing ratio
real, dimension(:,:), intent(in) :: stratcloud  ! stratiform cloud fraction
real, dimension(:,:), intent(in) :: clcon  ! convective cloud fraction
real, dimension(:,:), intent(in) :: pccw
real, dimension(:,:), intent(in) :: rhoa   ! density of air (kg/m3)
real, dimension(:,:), intent(in) :: pfprec, pfmelt, pfsnow         ! from LDR prog cloud
real, dimension(:,:), intent(in) :: pfsubl, plambs, pmrate         ! from LDR prog cloud
real, dimension(:,:), intent(in) :: pmaccr, pqfsedice, prscav      ! from LDR prog cloud
real, dimension(:,:), intent(in) :: prfreeze                       ! from LDR prog cloud
real, dimension(:,:), intent(in) :: pfevap                         ! from LDR prog cloud
real, dimension(ifull) :: cgssnowd,bbem
real, dimension(ifull) :: veff,vefn,v10n
real, dimension(ifull) :: Vgust_free,Vgust_deep
real, dimension(ifull,kl,naero) :: xte, xtm1, xtu
real, dimension(ifull,kl) :: prhop1,ptp1,pclcon
real, dimension(ifull,kl) :: pclcover,pcfcover,pmlwc,pmiwc,pfconv
real, dimension(ifull,kl) :: aphp2
real, dimension(ifull) :: fracc
real, dimension(ifull) :: so2oh,so2h2,so2o3,dmsoh,dmsn3
real, parameter :: beta = 0.65
real thetav, Wstar3, rrate, qtot
real, dimension(imax,kl) :: aphp1
real, dimension(imax,kl) :: lrhoa, ldz, lttg
real, dimension(imax,kl,naero) :: lxte, lxtg
real, dimension(imax,kl,4) :: lzoxidant
real, dimension(imax,naero) :: lxtem
real, dimension(imax,15) :: lemissfield
real, dimension(imax,ndcls) :: lerod
real, dimension(imax,ndust) :: dcola,dcolb
real, dimension(imax,ndust) :: oldduste,ldustwd,lduste
real, dimension(imax) :: oldsalte
real, dimension(imax,naero) :: burden
logical, dimension(:), intent(in) :: land   ! land/water mask (t=land).  Water includes lakes and ocean
logical, dimension(:), intent(in) :: locean ! sea mask without lakes (t=ocean)
integer nt,k,iq,tile,js,je

if ( maxval(xtg(1:ifull,:,:))>2.e-3 ) then
  write(6,*) "xtg out-of-range at start of aldrcalc"
  write(6,*) "xtg maxval,maxloc ",maxval(xtg(1:ifull,:,:)),maxloc(xtg(1:ifull,:,:))
end if


do tile = 1,ntiles
  js = (tile-1)*imax + 1
  je = tile*imax

  cgssnowd(js:je) = 1.E-3*snowd(js:je)

  ! Calculate sub-grid Vgust
  do iq = js,je
    v10n(iq) = ustar(iq)*log(10./zo(iq))/vkar ! neutral wind speed
    ! Mesoscale enhancement follows Redelsperger et al. (2000), J. Climate 13, 402-421.
    ! Equation numbers follow Fairall et al. 1996, JGR 101, 3747-3764.
    ! Calculate convective scaling velocity (Eq.17) and gustiness velocity (Eq.16)
    thetav = ttg(iq,1)*(1.+0.61*qvg(iq,1))
    Wstar3 = max(0.,(grav*pblh(iq)/thetav)*(fg(iq)/cp+0.61*ttg(iq,1)*eg(iq)/hl)/rhoa(iq,1))
    Vgust_free(iq) = beta*Wstar3**(1./3.)
    ! Calculate the Redelsperger-based Vgust_deep if deep convection is present.
    ! Note that Redelspreger gives two other parameterizations, based on
    ! the updraft or downdraught mass fluxes respectively.
    rrate = 8640.*condc(iq)/dt   !Rainfall rate in cm/day
    Vgust_deep(iq) = (19.8*rrate**2/(1.5+rrate+rrate**2))**0.4
  end do
end do

! Calculate effective 10m wind (Eq. 15)
! These can plausibly be added in quadrature, or linearly, the latter giving a much larger
! effect (Lunt & Valdes, JGR, 2002).
select case(enhanceu10)
  case(0)
    do tile = 1,ntiles
      js = (tile-1)*imax + 1
      je = tile*imax
      veff(js:je) = v10m(js:je)
      vefn(js:je) = v10n(js:je)
    end do
  case(1)
    do tile = 1,ntiles
      js = (tile-1)*imax + 1
      je = tile*imax
      veff(js:je) = sqrt( v10m(js:je)**2 + Vgust_free(js:je)**2 + Vgust_deep(js:je)**2 )
      vefn(js:je) = sqrt( v10n(js:je)**2 + Vgust_free(js:je)**2 + Vgust_deep(js:je)**2 )
    end do
  case(2)
    do tile = 1,ntiles
      js = (tile-1)*imax + 1
      je = tile*imax
      veff(js:je) = v10m(js:je) + Vgust_free(js:je) + Vgust_deep(js:je)
      vefn(js:je) = v10n(js:je) + Vgust_free(js:je) + Vgust_deep(js:je)
    end do
end select

! Emission and dry deposition (sulfur cycle and carbonaceous aerosols)
do tile = 1,ntiles
  js = (tile-1)*imax + 1
  je = tile*imax
  lrhoa = rhoa(js:je,:)
  ldz = dz(js:je,:)
  lemissfield = emissfield(js:je,:)
  lxtg = xtg(js:je,:,:)
  call xtemiss(dt, sig, lrhoa, ts(js:je), fracice(js:je), vefn(js:je),            & !Inputs
               land(js:je), tsigmf(js:je), cgssnowd(js:je), wg(js:je), ldz,       & !Inputs
               lxte, lxtem, bbem(js:je),                                          & !Outputs
               lemissfield,vso2(js:je),dmse(js:je),so2e(js:je),so4e(js:je),       & !Outputs
               bce(js:je),oce(js:je),lxtg,so2dd(js:je),so4dd(js:je),bcdd(js:je),  & !Output
               ocdd(js:je))                                                         !Inputs
  !xtem(js:je,:) = lxtem
  do nt = 1,naero
    do k = 1,kl
      xtg(js:je,k,nt) = max( xtg(js:je,k,nt)+lxte(:,k,nt)*dt, 0. )
    end do
  end do
end do


if ( maxval(xtg(1:ifull,:,:))>2.e-3 ) then
  write(6,*) "xtg out-of-range after xtemiss"
  write(6,*) "xtg maxval,maxloc ",maxval(xtg(1:ifull,:,:)),maxloc(xtg(1:ifull,:,:))
end if


do tile = 1,ntiles
  js = (tile-1)*imax + 1
  je = tile*imax
  ! Emission and dry deposition of dust
  do k = 1,kl
    ! calculate air pressure
    aphp1(:,k) = prf(js:je)*sig(k)*0.01 ! hPa
  end do
  lrhoa = rhoa(js:je,:)
  lxtg = xtg(js:je,:,:)
  ldz = dz(js:je,:)
  lttg = ttg(js:je,:)
  
  lerod = erod(js:je,:)
  lduste = duste(js:je,1:ndust)
  ! Calculate integrated column dust loading before settling and deposition
  oldduste(:,1:ndust) = lduste(:,1:ndust) ! duste is cumulative dust emissions
  dcola(:,1:ndust) = 0.
  do nt = 1,ndust
    do k = 1,kl
      dcola(:,nt) = dcola(:,nt) + lrhoa(:,k)*lxtg(:,k,itracdu+nt-1)*ldz(:,k)
    end do
  end do  
  ! Calculate the settling of large dust particles
  call dsettling(dt,lrhoa,lttg,ldz,aphp1,lxtg)
  call dustem(dt,lrhoa(:,1),wg(js:je),veff(js:je),ldz(:,1),vt(js:je),snowd(js:je), &
              lerod,lduste,lxtg)
  dcolb(:,1:ndust) = 0.
  do nt = 1,ndust
    do k = 1,kl
      ! Calculate integrated column dust after settling
      dcolb(:,nt) = dcolb(:,nt) + lrhoa(:,k)*lxtg(:,k,itracdu+nt-1)*ldz(:,k)
    end do
    ! Calculate deposition flux to surface
    dustdd(js:je,nt) = dustdd(js:je,nt) + (dcola(:,nt)-dcolb(:,nt))/dt + lduste(:,nt) - oldduste(:,nt)
  end do
  duste(js:je,:) = lduste

  call xtsink(dt,lxtg)

  oldsalte(:) = salte(js:je) ! salte is cumulative salt emissions
  dcola(:,1) = 0.
  do nt = 1,nsalt
    do k = 1,kl
      dcola(:,1) = dcola(:,1) + lrhoa(:,k)*lxtg(:,k,itracsa+nt-1)*ldz(:,k)
    end do
  end do
  ! Calculate the settling of large salt particles
  call ssettling(dt,lrhoa,lttg,ldz,aphp1,lxtg)
  ! Calculate salt emission and turbulent dry deposition at the surface
  call seasaltem(dt,veff(js:je),vt(js:je),lrhoa(:,1),ldz(:,1),salte(js:je), &
                 lxtg,locean(js:je))
  dcolb(:,1) = 0.
  do nt = 1,nsalt
    do k = 1,kl
      ! Calculate integrated column dust after settling
      dcolb(:,1) = dcolb(:,1) + lrhoa(:,k)*lxtg(:,k,itracsa+nt-1)*ldz(:,k)
    end do
  end do  
  ! Calculate deposition flux to surface
  saltdd(js:je) = saltdd(js:je) + (dcola(:,1)-dcolb(:,1))/dt + salte(js:je) - oldsalte(:)  
  xtg(js:je,:,:) = lxtg
end do


if ( maxval(xtg(1:ifull,:,:))>2.e-3 ) then
  write(6,*) "xtg out-of-range after settling, xtsink and em"
  write(6,*) "xtg maxval,maxloc ",maxval(xtg(1:ifull,:,:)),maxloc(xtg(1:ifull,:,:))
end if


do tile = 1,ntiles
  js = (tile-1)*imax + 1
  je = tile*imax
  
  ! Aerosol chemistry and wet deposition
  ! Need to invert vertical levels for ECHAM code... Don't you hate that?
  do nt = 1,naero
    do k = 1,kl
      xtm1(js:je,kl+1-k,nt) = xtg(js:je,k,nt)
      ! Convert from aerosol concentration outside convective cloud (used by CCAM)
      ! to aerosol concentration inside convective cloud
      xtu(js:je,kl+1-k,nt) = max(xtg(js:je,k,nt)-(1.-clcon(js:je,k))*xtosav(js:je,k,nt),0.)/max(clcon(js:je,k),1.E-8)
    end do
  end do
  do k = 1,kl
    do iq = js,je
      aphp2(iq,kl+1-k)  = rhoa(iq,k)*dz(iq,k)                       ! density * thickness
      prhop1(iq,kl+1-k) = rhoa(iq,k)                                ! air density
      ptp1(iq,kl+1-k)   = ttg(iq,k)                                 ! air temperature
      pclcon(iq,kl+1-k) = min(max(clcon(iq,k),0.),1.)               ! convective cloud fraction
      qtot = qlg(iq,k) + qfg(iq,k)                                  ! total liquid and ice mixing ratio
      pclcover(iq,kl+1-k) = stratcloud(iq,k)*qlg(iq,k)/max(qtot,1.E-8) ! Liquid-cloud fraction
      pcfcover(iq,kl+1-k) = stratcloud(iq,k)*qfg(iq,k)/max(qtot,1.E-8) ! Ice-cloud fraction
      pmlwc(iq,kl+1-k) = qlg(iq,k)
      pmiwc(iq,kl+1-k) = qfg(iq,k)
      if ( k<=kbsav(iq) ) then
        pfconv(iq,kl+1-k) = condc(iq)/dt
      else
        pfconv(iq,kl+1-k) = 0.
      end if
      fracc(iq) = cldcon(iq)
    end do
  end do
end do

call xtchemie(imax, 2, dt, zdayfac, aphp2, pmrate, pfprec,           & !Inputs
              pclcover, pmlwc, prhop1, ptp1, taudar, xtm1,           & !Inputs
              pfsnow,pfsubl,pcfcover,pmiwc,pmaccr,pfmelt,            & !Inputs
              pqfsedice,plambs,prscav,prfreeze,pfevap,pclcon,fracc,  & !Inputs
              pccw,pfconv,xtu,                                       & !Inputs
              xte, so2oh, so2h2, so2o3, dmsoh, dmsn3,                & !Output
              zoxidant_g,so2wd,so4wd,bcwd,ocwd,dustwd,saltwd)

do tile = 1,ntiles
  js = (tile-1)*imax + 1
  je = tile*imax  
  do nt = 1,naero
    do k = 1,kl
      xtg(js:je,k,nt) = max( xtg(js:je,k,nt)+xte(js:je,kl+1-k,nt)*dt, 0. )
    end do
  enddo
  dmsso2o(js:je) = dmsso2o(js:je) + dmsoh(js:je) + dmsn3(js:je)          ! oxidation of DMS to SO2
  so2so4o(js:je) = so2so4o(js:je) + so2oh(js:je) + so2h2(js:je) + so2o3(js:je)  ! oxidation of SO2 to SO4
end do


if ( maxval(xtg(1:ifull,:,:))>2.e-3 ) then
  write(6,*) "xtg out-of-range after xtchemie"
  write(6,*) "xtg maxval,maxloc ",maxval(xtg(1:ifull,:,:)),maxloc(xtg(1:ifull,:,:))
end if


do tile = 1,ntiles
  js = (tile-1)*imax + 1
  je = tile*imax

  burden(:,:) = 0.
  do k = 1,kl
    do nt = 1,naero
      do iq = 1,imax
        burden(iq,nt) = burden(iq,nt) + xtg(iq+js-1,k,nt)*rhoa(iq+js-1,k)*dz(iq+js-1,k)
      end do
    end do
  end do

  do nt = 1,ndust
    dust_burden(js:je,nt) = dust_burden(js:je,nt) + burden(:,nt+itracdu-1)
  end do

  do nt = itracbc,itracbc+1
    bc_burden(js:je) = bc_burden(js:je) + burden(:,nt)
  end do

  do nt = itracoc,itracoc+1
    oc_burden(js:je) = oc_burden(js:je) + burden(:,nt)
  end do

  dms_burden(js:je) = dms_burden(js:je) + burden(:,itracdms)
  so2_burden(js:je) = so2_burden(js:je) + burden(:,itracso2)
  so4_burden(js:je) = so4_burden(js:je) + burden(:,itracso4)

  do nt = itracsa,itracsa+nsalt-1
    salt_burden(js:je) = salt_burden(js:je) + burden(:,nt)
  end do

end do

return
end subroutine aldrcalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! xt emiss

SUBROUTINE XTEMISS(ztmst, sig, rhoa, TSM1M, SEAICEM, ZZSPEED,                    & !Inputs
                   LOLAND, PFOREST, PSNOW, WSM1M, dz,                            & !Inputs
                   XTE, PXTEMS, bbem,                                            & !Outputs
                   emissfield,vso2,dmse,so2e,so4e,bce,oce,xtg,so2dd,so4dd,bcdd,ocdd)

!
!    THIS ROUTINE CALCULATES THE LOWER BOUNDARY CONDITIONS
!    FOR VDIFF DEPENDING ON THE SURFACE EMISSION AND THE
!    DRY DEPOSITION FLUX.
!
!    JOHANN FEICHTER          UNI/HAMBURG         08/91
!    MODIFIED  U. SCHLESE    DKRZ-HAMBURG        JAN-95
!    Adapted for CSIRO GCM by Leon Rotstayn, 12/99
!    Adapted for CCAM by Marcus Thatcher 2012
!
!    PURPOSE
!   ---------
!    THE LOWER BOUNDARY CONDITION FOR CALCULATING THE
!    TURBULENT EXCHANGE IN THE BOUNDARY LAYER IS
!    DETERMINED BY THE EMISSION AND THE DRY DEPOSITION FLUX.
!

! Argument list
real, intent(in) :: ztmst                           !Timestep [s]
real, dimension(:), intent(in) :: sig
real, dimension(:,:), intent(in) :: rhoa            !Density of air
real, dimension(:), intent(in) :: TSM1M             !Surface temp
real, dimension(:), intent(in) :: SEAICEM           !Sea-ice fraction
real, dimension(:), intent(in) :: ZZSPEED           !10m wind (corrected to neutral for Nightingale scheme)
real, dimension(:,:), intent(in) :: dz              ! layer thickness [m]
real, dimension(:), intent(in) :: PFOREST           !Fractional vegetation cover
real, dimension(:), intent(in) :: PSNOW             !Snow depth [m]
! Land-surface details needed to specify dry deposition velocity
real, dimension(:), intent(in) :: WSM1M             !surface wetness [vol fraction for CSIRO GCM, not m]
real, dimension(:,:,:), intent(out) :: XTE          !Tracer tendencies (kg/kg/s)
real, dimension(:,:), intent(out) :: PXTEMS         !Sfc. flux of tracer passed to vertical mixing [kg/m2/s]
logical, dimension(:), intent(in) :: LOLAND         !Land flag
! Some diagnostics
real, dimension(:), intent(out) :: bbem

real, dimension(:,:), intent(in) :: emissfield
real, dimension(:), intent(in) :: vso2
real, dimension(:), intent(inout) :: dmse
real, dimension(:), intent(inout) :: so2e
real, dimension(:), intent(inout) :: so4e
real, dimension(:), intent(inout) :: bce
real, dimension(:), intent(inout) :: oce
real, dimension(:,:,:), intent(in) :: xtg
real, dimension(:), intent(inout) :: so2dd
real, dimension(:), intent(inout) :: so4dd
real, dimension(:), intent(inout) :: bcdd
real, dimension(:), intent(inout) :: ocdd
real zdmscon, zsst, scdms, vpco2, zvdms, vpco2liss
real wtliss, gdp, zdmsemiss
real zhilbco, zhilbcy, zhiloco, zhilocy, zhilso2, zhilso4
real zvd2nof, zvd4nof
real zvdrd1, zvdrd2
integer jk2, jk3, jk4, jk5, jk6, jk8, jk9, k
integer jk, iq, imax, kl

!     M WATER EQUIVALENT  CRITICAL SNOW HEIGHT (FROM *SURF*)
real, parameter :: ZSNCRI = 0.025
!     COEFFICIENTS FOR ZVDRD = FUNCTION OF SOIL MOISTURE
real, parameter :: ZVWC2 = (0.8E-2 - 0.2E-2)/(1. - 0.9)
real, parameter :: ZVW02 = ZVWC2-0.8E-2
real, parameter :: ZVWC4 = (0.2E-2 - 0.025E-2)/(1. - 0.9)
real, parameter :: ZVW04 = ZVWC4-0.2E-2
real, parameter :: tmelt = 273.05
!     Dry deposition
real, parameter :: ZVDPHOBIC = 0.025E-2
!     DMS emissions
real, parameter :: ScCO2     = 600.
!real, parameter :: a_vpco2   = 0.222 ! nightingale (2000)
!real, parameter :: b_vpco2   = 0.333 ! nightingale (2000)
real, parameter :: a_vpco2  = 0.166 ! approx Liss and Merlivat (see nightingale 2000)
real, parameter :: b_vpco2  = 0.133 ! approx Liss and Merlivat (see nightingale 2000)

! Start code : ----------------------------------------------------------

imax = size(rhoa,1)
kl   = size(rhoa,2)

! MJT - define injection levels

! scale to CSIRO9 18 levels
jk2=2 ! jk2 is top of level=1, bottom of level=2
jk3=3
jk4=4
jk5=5
jk6=6
jk8=7
jk9=8
do k = 3,kl
  if ( sig(k)>=0.967 ) then ! 300m
    jk3 = k ! jk3 is top of level=2, bottom of level=3
  end if
  if ( sig(k)>=0.930 ) then ! 600m
    jk4 = k ! jk4 is top of level=3, bottom of level=4
  end if
  if ( sig(k)>=0.882 ) then ! 1,000m
    jk5 = k ! jk5 is top of level=4, bottom of level=5
  end if
  if ( sig(k)>=0.800 ) then ! 1,800m
    jk6 = k ! jk6 is top of level=5, bottom of level=6
  end if
  if ( sig(k)>=0.530 ) then ! 5,000m
    jk8 = k ! jk8 is top of level=7, bottom of level=8
  end if
  if ( sig(k)>=0.350 ) then ! 8,000m
    jk9 = k ! jk9 is top of level=8, bottom of level=9
  end if
end do

pxtems(:,:) = 0.
xte(:,:,:) = 0.

! --------------------------------------------------------------
!
!*     1.   SURFACE EMISSION.
!           ------- --------
!
!   CALCULATE DMS EMISSIONS FOLLOWING LISS+MERLIVAT
!   DMS SEAWATER CONC. FROM KETTLE ET AL.

do iq = 1,imax
  ZDMSCON = EMISSFIELD(iq,idmso)*(1.-SEAICEM(iq))**2
  ZSST = min( TSM1M(iq)-273.15, 45. )   ! Even Saltzman Sc formula has trouble over 45 deg C
  ! The formula for ScDMS from Saltzman et al (1993) is given by Kettle & Andreae (ref below)
  ScDMS = 2674. - 147.12*ZSST + 3.726*ZSST**2 - 0.038*ZSST**3 !Sc for DMS (Saltzman et al.)
  ! Nightingale (2000) scheme (J. Biogeochem. Cycles, 14, 373-387)
  ! For this scheme, zzspeed is the 10m wind adjusted to neutral stability.
  VpCO2 = a_vpco2*zzspeed(iq)**2 + b_vpco2*zzspeed(iq) !Nightingale et al
  !  ZZSPEED:  10-M WINDS
  if ( ZZSPEED(iq)<3.6 ) then
    zVdms = VpCO2*(ScCO2/ScDMS)**(2./3.)
  else if ( zzspeed(iq)<20. ) then
    ! Phase in Liss & Merlivat from 13 to 18 m/s, since Nightingale is doubtful for high windspeeds,
    ! due to limited data.
    VpCO2liss = 5.9*ZZSPEED(iq) - 49.3
    wtliss = min( max( (zzspeed(iq)-13.)/5., 0. ), 1. )
    VpCO2 = wtliss*VpCO2liss + (1.-wtliss)*VpCO2        
    zVdms = VpCO2*sqrt(ScCO2/ScDMS)
  else
    ! limit wind speed to 20 m/s for emissions - MJT suggestion  
    VpCO2liss = 5.9*20. - 49.3
    wtliss = 1.
    VpCO2 = VpCO2liss
    zVdms = VpCO2*sqrt(ScCO2/ScDMS)
  end if
  if ( loland(iq) ) then
    zdmsemiss = (1./1.938)*emissfield(iq,idmst) !kgS/m2/s
  else
    zdmsemiss = ZDMSCON*ZVDMS*32.06e-11/3600.
    ! NANOMOL/LTR*CM/HOUR --> KG/M**2/SEC
  end if
  gdp = 1./(rhoa(iq,1)*dz(iq,1))
  xte(iq,1,itracdms) = xte(iq,1,itracdms) + zdmsemiss*gdp
end do

if ( any(xtg(1:imax,1,itracdms)+xte(1:imax,1,itracdms)*ztmst>2.e-3) ) then
  write(6,*) "xtg out-of-range after DMS emission (xtemiss)"
  write(6,*) "xtg maxval,maxloc ",maxval(xtg(1:imax,1,itracdms)+xte(1:imax,1,itracdms)*ztmst), &
                                  maxloc(xtg(1:imax,1,itracdms)+xte(1:imax,1,itracdms)*ztmst)
end if

! Other biomass emissions of SO2 are done below (with the non-surface S emissions)
PXTEMS(:,ITRACSO2)  =(EMISSFIELD(:,iso2a1)+EMISSFIELD(:,iso2b1))*0.97
PXTEMS(:,ITRACSO4)  =(EMISSFIELD(:,iso2a1)+EMISSFIELD(:,iso2b1))*0.03
! Apply these here as a tendency (XTE), rather than as a surface flux (PXTEMS) via vertmix.
do iq = 1,imax
  gdp=1./(rhoa(iq,1)*dz(iq,1))
  xte(iq,1,itracso2)  =xte(iq,1,itracso2)  +pxtems(iq,itracso2)*gdp
  xte(iq,1,itracso4)  =xte(iq,1,itracso4)  +pxtems(iq,itracso4)*gdp
end do

!  EMISSION OF ANTHROPOGENIC SO2 IN THE NEXT HIGHER LEVEL PLUS BIOMASS BURNING
do jk=jk2,jk3-1
  do iq = 1,imax
    gdp=1./(rhoa(iq,jk)*dz(iq,jk))/real(jk3-jk2)
    XTE(iq,JK,ITRACSO2)  =XTE(iq,JK,ITRACSO2)  +0.97*EMISSFIELD(iq,iso2a2)*gdp !100% of the "above 100m" SO2 emission
    XTE(iq,JK,ITRACSO4)  =XTE(iq,JK,ITRACSO4)  +0.03*EMISSFIELD(iq,iso2a2)*gdp !100% of the "above 100m" SO4 emission
  end do
end do
do jk=jk3,jk4-1
  do iq = 1,imax
    gdp=1./(rhoa(iq,jk)*dz(iq,jk))/real(jk4-jk3)
    xte(iq,jk,ITRACSO2)=xte(iq,jk,ITRACSO2)+0.3*emissfield(iq,iso2b2)*gdp
  end do
end do
do jk=jk4,jk5-1
  do iq = 1,imax
    gdp=1./(rhoa(iq,jk)*dz(iq,jk))/real(jk5-jk4)
    xte(iq,jk,ITRACSO2)=xte(iq,jk,ITRACSO2)+0.4*emissfield(iq,iso2b2)*gdp
  end do
end do
do jk=jk5,jk6-1
  do iq = 1,imax
    gdp=1./(rhoa(iq,jk)*dz(iq,jk))/real(jk6-jk5)
    xte(iq,jk,ITRACSO2)=xte(iq,jk,ITRACSO2)+0.3*emissfield(iq,iso2b2)*gdp
  end do
end do
  
!    VOLCANIC BACKGROUND EMISSIONS 
!
!   3 EMISSION LEVELS: 
!    1. PRE-INTRA ERUPTION IN LEVEL IVOLC-HEIGHT (=TOP OF VOLCANO)
!    2. POST-EXTRA ERUPTION IN LEVEL 15 -16 (CA 550-1736M)
!    3. EXPLOSIVE ERUPTION IN LEVEL 10 - 11 (CA 5000-7900M)
do iq = 1,imax
  gdp=1./(rhoa(iq,1)*dz(iq,1))
  XTE(iq,1,ITRACSO2)=XTE(iq,1,ITRACSO2)+ZVOLCEMI*0.36*vso2(iq)*gdp
end do
do jk=jk4,jk6-1
  do iq = 1,imax
    gdp=1./(rhoa(iq,jk)*dz(iq,jk))/real(jk6-jk4)
    XTE(iq,jk,ITRACSO2)=XTE(iq,jk,ITRACSO2)+ZVOLCEMI*0.36*vso2(iq)*gdp
  end do
end do
do jk=jk8,jk9-1
  do iq = 1,imax
    gdp=1./(rhoa(iq,jk)*dz(iq,jk))/real(jk9-jk8)
    XTE(iq,jk,ITRACSO2)=XTE(iq,jk,ITRACSO2)+ZVOLCEMI*0.28*vso2(iq)*gdp
  end do
end do


!Do carbonaceous aerosols
! Inject the low-level fossil-fuel and natural SOA emissions into layer 1
! Assume BC 80% hydrophobic, OC 50%.
PXTEMS(:,ITRACBC)  =0.8*EMISSFIELD(:,ibca1)
PXTEMS(:,ITRACBC+1)=0.2*EMISSFIELD(:,ibca1)
PXTEMS(:,ITRACOC)  =0.5*(EMISSFIELD(:,ioca1)+EMISSFIELD(:,iocna))
PXTEMS(:,ITRACOC+1)=0.5*(EMISSFIELD(:,ioca1)+EMISSFIELD(:,iocna))
! Apply these here as a tendency (XTE), rather than as a surface flux (PXTEMS) via vertmix.
do iq = 1,imax
  gdp=1./(rhoa(iq,1)*dz(iq,1))
  xte(iq,1,itracbc)  =xte(iq,1,itracbc)  +pxtems(iq,itracbc)*gdp
  xte(iq,1,itracbc+1)=xte(iq,1,itracbc+1)+pxtems(iq,itracbc+1)*gdp
  xte(iq,1,itracoc)  =xte(iq,1,itracoc)  +pxtems(iq,itracoc)*gdp
  xte(iq,1,itracoc+1)=xte(iq,1,itracoc+1)+pxtems(iq,itracoc+1)*gdp
end do
! Inject the upper-level fossil-fuel emissions into layer 2
! Assume BC 80% hydrophobic, OC 50%.
PXTEMS(:,ITRACBC)  =0.8*EMISSFIELD(:,ibca2)
PXTEMS(:,ITRACBC+1)=0.2*EMISSFIELD(:,ibca2)
PXTEMS(:,ITRACOC)  =0.5*EMISSFIELD(:,ioca2)
PXTEMS(:,ITRACOC+1)=0.5*EMISSFIELD(:,ioca2)
! Apply these here as a tendency (XTE), rather than as a surface flux (PXTEMS) via vertmix.
do jk=jk2,jk3-1
  do iq = 1,imax
    gdp=1./(rhoa(iq,jk)*dz(iq,jk))/real(jk3-jk2)
    xte(iq,jk,itracbc)  =xte(iq,jk,itracbc)  +pxtems(iq,itracbc)*gdp
    xte(iq,jk,itracbc+1)=xte(iq,jk,itracbc+1)+pxtems(iq,itracbc+1)*gdp
    xte(iq,jk,itracoc)  =xte(iq,jk,itracoc)  +pxtems(iq,itracoc)*gdp
    xte(iq,jk,itracoc+1)=xte(iq,jk,itracoc+1)+pxtems(iq,itracoc+1)*gdp
  end do
end do
! Inject the lower-level biomass emissions into layer 2 (NB: Doesn't include biofuel any more)
! Assume BC and OC both 50% hydrophobic.
PXTEMS(:,ITRACBC)  =0.5*EMISSFIELD(:,ibcb1)
PXTEMS(:,ITRACBC+1)=0.5*EMISSFIELD(:,ibcb1)
PXTEMS(:,ITRACOC)  =0.5*EMISSFIELD(:,iocb1)
PXTEMS(:,ITRACOC+1)=0.5*EMISSFIELD(:,iocb1)
! Apply these here as a tendency (XTE)
do jk=jk2,jk3-1
  do iq = 1,imax
    gdp=1./(rhoa(iq,jk)*dz(iq,jk))/real(jk3-jk2)
    xte(iq,jk,itracbc)  =xte(iq,jk,itracbc)  +pxtems(iq,itracbc)*gdp
    xte(iq,jk,itracbc+1)=xte(iq,jk,itracbc+1)+pxtems(iq,itracbc+1)*gdp
    xte(iq,jk,itracoc)  =xte(iq,jk,itracoc)  +pxtems(iq,itracoc)*gdp
    xte(iq,jk,itracoc+1)=xte(iq,jk,itracoc+1)+pxtems(iq,itracoc+1)*gdp
  end do
end do
! Inject the upper-level biomass emissions into layers 3-5 (30%, 40%, 30%)
! Assume BC and OC both 50% hydrophobic.
PXTEMS(:,ITRACBC)  =0.5*EMISSFIELD(:,ibcb2)
PXTEMS(:,ITRACBC+1)=0.5*EMISSFIELD(:,ibcb2)
PXTEMS(:,ITRACOC)  =0.5*EMISSFIELD(:,iocb2)
PXTEMS(:,ITRACOC+1)=0.5*EMISSFIELD(:,iocb2)
! Apply these here as a tendency (XTE)
do jk=jk3,jk4-1
  do iq = 1,imax
    gdp=1./(rhoa(iq,jk)*dz(iq,jk))/real(jk4-jk3)
    xte(iq,jk,itracbc)  =xte(iq,jk,itracbc)  +0.3*pxtems(iq,itracbc)*gdp
    xte(iq,jk,itracbc+1)=xte(iq,jk,itracbc+1)+0.3*pxtems(iq,itracbc+1)*gdp
    xte(iq,jk,itracoc)  =xte(iq,jk,itracoc)  +0.3*pxtems(iq,itracoc)*gdp
    xte(iq,jk,itracoc+1)=xte(iq,jk,itracoc+1)+0.3*pxtems(iq,itracoc+1)*gdp
  end do
end do
do jk=jk4,jk5-1
  do iq = 1,imax
    gdp=1./(rhoa(iq,jk)*dz(iq,jk))/real(jk5-jk4)
    xte(iq,jk,itracbc)  =xte(iq,jk,itracbc)  +0.4*pxtems(iq,itracbc)*gdp
    xte(iq,jk,itracbc+1)=xte(iq,jk,itracbc+1)+0.4*pxtems(iq,itracbc+1)*gdp
    xte(iq,jk,itracoc)  =xte(iq,jk,itracoc)  +0.4*pxtems(iq,itracoc)*gdp
    xte(iq,jk,itracoc+1)=xte(iq,jk,itracoc+1)+0.4*pxtems(iq,itracoc+1)*gdp
  end do
end do
do jk=jk5,jk6-1
  do iq = 1,imax
    gdp=1./(rhoa(iq,jk)*dz(iq,jk))/real(jk6-jk5)
    xte(iq,jk,itracbc)  =xte(iq,jk,itracbc)  +0.3*pxtems(iq,itracbc)*gdp
    xte(iq,jk,itracbc+1)=xte(iq,jk,itracbc+1)+0.3*pxtems(iq,itracbc+1)*gdp
    xte(iq,jk,itracoc)  =xte(iq,jk,itracoc)  +0.3*pxtems(iq,itracoc)*gdp
    xte(iq,jk,itracoc+1)=xte(iq,jk,itracoc+1)+0.3*pxtems(iq,itracoc+1)*gdp
  end do
end do

!   --------------------------------------------------------------
!
!*      2.    DRY DEPOSITION.
!             --- ----------

! Sulfur emission diagnostic (hard-coded for 3 sulfur variables)
do jk=1,kl
  dmse=dmse+xte(:,jk,ITRACDMS)*rhoa(:,jk)*dz(:,jk)   !Above surface
  so2e=so2e+xte(:,jk,ITRACSO2)*rhoa(:,jk)*dz(:,jk)   !Above surface
  so4e=so4e+xte(:,jk,ITRACSO4)*rhoa(:,jk)*dz(:,jk)   !Above surface

  ! Assume that BC and OC emissions are passed in through xte()
  bce=bce+(xte(:,jk,ITRACBC)+xte(:,jk,ITRACBC+1))*rhoa(:,jk)*dz(:,jk)
  oce=oce+(xte(:,jk,ITRACOC)+xte(:,jk,ITRACOC+1))*rhoa(:,jk)*dz(:,jk)
enddo

! Total biomass burning primary emissions (note 1.3 for organic carbon)
bbem=emissfield(:,ibcb1)+emissfield(:,ibcb2)+1.3*(emissfield(:,iocb1)+emissfield(:,iocb2))

do iq = 1,imax

!      DRY DEPOSITION OF SO2, SO4

  !         -  SNOW/NO SNOW -
  if ( PSNOW(iq)>ZSNCRI .and. tsm1m(iq)>=tmelt ) then
  !            - MELTING/NOT MELTING SNOW -
    ZVD2NOF=0.8E-2
    ZVD4NOF=0.2E-2
  else if ( PSNOW(iq)>ZSNCRI ) then
    ZVD2NOF=0.1E-2
    ZVD4NOF=0.025E-2
  else if ( tsm1m(iq)<=tmelt ) then
  !           -  FROZEN SOIL -
    ZVD2NOF=0.2E-2
    ZVD4NOF=0.025E-2
  else
  !           - PARTLY WET -
    ZVD2NOF=max(min(ZVWC2*WSM1M(iq)-ZVW02,0.8E-2),0.2E-2)
    ZVD4NOF=max(min(ZVWC4*WSM1M(iq)-ZVW04,0.2E-2),0.025E-2)
  end if

  ! ZVDRD   DRY DEPOSITION VELOCITY IN M/S
  ! ZVDRD1  FOR SO2 GAS
  ! ZVDRD2  FOR AEROSOLS

  !     -  SEA -
  if ( .NOT.LOLAND(iq) .and. tsm1m(iq)>=(tmelt-0.1) ) then
  !         - MELTING SEA ICE -
    ZVDRD1=(1.-SEAICEM(iq))*0.8E-2+SEAICEM(iq)*0.8E-2 !So leads agree with ocean
    ZVDRD2=(1.-SEAICEM(iq))*0.2E-2+SEAICEM(iq)*0.2E-2
  else if ( .not.loland(iq) ) then
  !         - NOT MELTING SEA ICE -
    ZVDRD1=(1.-SEAICEM(iq))*0.8E-2+SEAICEM(iq)*0.1E-2 !So leads agree with ocean
    ZVDRD2=(1.-SEAICEM(iq))*0.2E-2+SEAICEM(iq)*0.025E-2
  else
  !      - LAND -
    ZVDRD1=PFOREST(iq)*0.8E-2+(1.-PFOREST(iq))*ZVD2NOF
    ZVDRD2=PFOREST(iq)*0.2E-2+(1.-PFOREST(iq))*ZVD4NOF
  end if


  gdp=1./(rhoa(iq,1)*dz(iq,1))

  zhilso2=(xtg(iq,1,itracso2)+xte(iq,1,itracso2)*ztmst)   &
         *(1.-exp(-ztmst*zvdrd1/dz(iq,1)))/(ztmst*gdp)
  xte(iq,1,ITRACSO2)  =xte(iq,1,ITRACSO2)  -zhilso2*gdp
  
  zhilso4=(xtg(iq,1,itracso4)+xte(iq,1,itracso4)*ztmst)   &
         *(1.-exp(-ztmst*zvdrd2/dz(iq,1)))/(ztmst*gdp)
  xte(iq,1,ITRACSO4)  =xte(iq,1,ITRACSO4)  -zhilso4*gdp

  ZHILBCO=(xtg(iq,1,ITRACBC)+xte(iq,1,itracbc)*ztmst)     &
         *(1.-exp(-ztmst*ZVDPHOBIC/dz(iq,1)))/(ztmst*gdp)
  xte(iq,1,itracbc)  =xte(iq,1,itracbc)    -zhilbco*gdp

  ZHILBCY=(xtg(iq,1,ITRACBC+1)+xte(iq,1,itracbc+1)*ztmst) &
         *(1.-exp(-ztmst*ZVDRD2/dz(iq,1)))/(ztmst*gdp)
  xte(iq,1,itracbc+1)=xte(iq,1,itracbc+1)  -zhilbcy*gdp

  ZHILOCO=(xtg(iq,1,ITRACOC)+xte(iq,1,itracoc)*ztmst)     &
         *(1.-exp(-ztmst*ZVDPHOBIC/dz(iq,1)))/(ztmst*gdp)
  xte(iq,1,itracoc)  =xte(iq,1,itracoc)    -zhiloco*gdp

  ZHILOCY=(xtg(iq,1,ITRACOC+1)+xte(iq,1,itracoc+1)*ztmst) &
         *(1.-exp(-ztmst*ZVDRD2/dz(iq,1)))/(ztmst*gdp)
  xte(iq,1,itracoc+1)=xte(iq,1,itracoc+1)  -zhilocy*gdp

  so2dd(iq)=so2dd(iq)+zhilso2
  so4dd(iq)=so4dd(iq)+zhilso4
  bcdd(iq)=bcdd(iq)+zhilbco+zhilbcy
  ocdd(iq)=ocdd(iq)+zhiloco+zhilocy

end do


return
end subroutine xtemiss

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! xt sink

SUBROUTINE XTSINK(PTMST,xtg)

!
!   *XTSINK*  CALCULATES THE DECREASE OF TRACER CONCENTRATION
!             FOR  A GIVEN HALF-LIFE-TIME
!
!   JOHANN FEICHTER               UNI-HAMBURG    08-91
!
!   PURPOSE
!  ---------------
!   THE MASS MIXING-RATIO OF TRACERS IS MULTIPLIED WITH
!   EXP(LOG(0.5)*TIME-STEP/HALF-LIFE-TIME).
!

REAL, intent(in) :: PTMST
real, dimension(:,:,:), intent(inout) :: xtg
real zdxtdt
real pqtmst,zfac,zdecay
integer jk, jl, imax, kl

! Start code : ----------------------------------------------------------

imax = size(xtg,1)
kl   = size(xtg,2)

PQTMST=1./PTMST
ZFAC=LOG(0.5)*PTMST

ZDECAY=EXP(ZFAC/86400.) ! 1 day
do jk = 1,kl
  do jl = 1,imax
    ZDXTDT=xtg(JL,JK,ITRACBC)*(ZDECAY-1.)
    xtg(JL,JK,ITRACBC)   = xtg(JL,JK,ITRACBC)   + ZDXTDT
    xtg(JL,JK,ITRACBC+1) = xtg(JL,JK,ITRACBC+1) - ZDXTDT
  end do
end do

ZDECAY=EXP(ZFAC/86400.) ! 1 day
do jk = 1,kl
  do jl = 1,imax
    ZDXTDT=xtg(JL,JK,ITRACOC)*(ZDECAY-1.)
    xtg(JL,JK,ITRACOC)   = xtg(JL,JK,ITRACOC)   + ZDXTDT
    xtg(JL,JK,ITRACOC+1) = xtg(JL,JK,ITRACOC+1) - ZDXTDT
  end do
end do

RETURN
END subroutine xtsink

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! xt chemie

SUBROUTINE XTCHEMIE(imax, KTOP, PTMST, zdayfac, rhodz, PMRATEP, PFPREC,              & !Inputs
                    PCLCOVER, PMLWC, PRHOP1, PTP1, taudar, xtm1,                     & !Inputs
                    pfsnow,pfsubl,pcfcover,pmiwc,pmaccr,pfmelt,                      & !Inputs
                    pqfsedice,plambs,prscav,prfreeze,pfevap,pclcon,fracc,            & !Inputs
                    pccw,pfconv,xtu,                                                 & !Inputs
                    xte,so2oh,so2h2,so2o3,dmsoh,dmsn3,                               & !Outputs
                    zoxidant,so2wd,so4wd,bcwd,ocwd,dustwd,saltwd)                                                           !Inputs

! Inputs
! ktop: top level for aerosol processes (set to 1, counting downwards from top)
! ptmst: timestep (seconds; tdt in main program)
! rhodz: density * thickness (si units)
! pmratep: precip formation rate (kg/kg/s)
! pfprec: rainfall flux (entering from above) (kg/m2/s)
! pclcover: liquid-water cloud fraction (input; don't pass in cfrac though)
! pmlwc: liquid-water mixing ratio (kg/kg)
! prhop1: density of air (kg/m3)
! ptp1: temperature (k)
! taudar: fraction of time sunlit (used to determine if daytime)
! xtm1: tracer mixing ratio (kg/kg)
! pfsnow: snowfall flux (entering from above) (kg/m2/s)
! pfsubl: snowfall flux evaporating in layer k (kg/m2/s)
! pcfcover: ice cloud fraction
! pmiwc: ice mixing ratio (kg/kg)
! pmaccr: accretion rate (kg/kg/s)
! pfmelt: snowfall flux melting in layer k (kg/m2/s)
! pfstayice: snowfall flux staying in layer k (kg/m2/s)
! pqfsedice: fractional ice sedimentation in timestep
! pfevap: rainfall flux evaporating in layer k (kg/m2/s)
! plambs: slope (lambda) for snow crystal size distribution (m**-1)
! prscav: fractional rain scavenging rate in time step (needs to be mult. by coll. eff.)
! pclcon: convective cloud fraction
! fracc: Convective rain fraction
! pccw: convective cloud water mixing ratio (kg/kg)
! pfconv: convective rainfall flux (kg/m2/s)
! xtu: tracer mixing ratio in convective updraught (kg/kg)

! Outputs
! xte: tracer tendency (kg/kg/s)
! so2oh: oxidation of SO2 by OH (diagnostic)
! so2h2: oxidation of SO2 by H2O2 (diagnostic)
! so2o3: oxidation of SO2 by O3 (diagnostic)
! dmsoh: oxidation of DMS by OH (diagnostic)
! dmsn3: oxidation of DMS by NO3 (diagnostic)

!**** *XTCHEMIE*  CALCULATES DRY AND WET CHEMISTRY
!
!      J. FEICHTER             UNI HAMBURG    30/06/92
!
!      PURPOSE
!      ---------
!      THIS ROUTINE COMPUTES THE OXIDATION AND THE WET SCAVENGING
!      OF CHEMICAL SPECIES.
!
!**    INTERFACE
!      -----------
!      *XTCHEMIE*   IS CALLED FROM  progcld in CSIRO GCM
!
!      EXTERNALS
!      ------------
!          *XTWETDEP*  CALCULATES THE WET DEPOSITION

! Argument list
integer, intent(in) :: KTOP, imax
real, intent(in) :: PTMST
real, dimension(:,:,:), intent(inout) :: XTM1
real, dimension(:,:,:), intent(inout) :: xtu
real, dimension(:,:), intent(in) :: rhodz
real, dimension(:,:), intent(in) :: PMRATEP
real, dimension(:,:), intent(in) :: PFPREC
real, dimension(:,:), intent(in) :: PCLCOVER
real, dimension(:,:), intent(in) :: PMLWC
real, dimension(:,:), intent(in) :: PRHOP1
real, dimension(:,:), intent(in) :: PTP1
real, dimension(:,:), intent(in) :: pfsnow
real, dimension(:,:), intent(in) :: pfconv
real, dimension(:,:), intent(in) :: pfsubl
real, dimension(:,:), intent(in) :: pcfcover
real, dimension(:,:), intent(in) :: pmiwc
real, dimension(:,:), intent(in) :: pmaccr
real, dimension(:,:), intent(in) :: pfmelt
real, dimension(:,:), intent(in) :: pqfsedice
real, dimension(:,:), intent(in) :: plambs
real, dimension(:,:), intent(in) :: prscav
real, dimension(:,:), intent(in) :: prfreeze
real, dimension(:,:), intent(in) :: pfevap
real, dimension(:,:), intent(in) :: pclcon
real, dimension(:,:), intent(in) :: pccw
real, dimension(:), intent(in) :: taudar
real, dimension(:), intent(in) :: fracc
real, dimension(:), intent(in) :: zdayfac
real, dimension(:,:,:), intent(out) :: xte
real, dimension(:), intent(out) :: dmsoh, dmsn3, so2oh, so2h2, so2o3 !Diagnostic output

real, dimension(:,:,:), intent(in) :: zoxidant
real, dimension(:), intent(inout) :: so2wd
real, dimension(:), intent(inout) :: so4wd
real, dimension(:), intent(inout) :: bcwd
real, dimension(:), intent(inout) :: ocwd
real, dimension(:,:), intent(inout) :: dustwd
real, dimension(:), intent(inout) :: saltwd

! Local work arrays and variables
integer jt, jk, jn
integer jl, ifull, kl, iq, ntiles, tile, js, je
integer i, ktrac
real, dimension(imax,size(rhodz,2),naero) :: xto
real, dimension(imax,size(rhodz,2)) :: ZHENRY, ZSO4, ZSO4i, ZSO4C, ZHENRYC 
real, dimension(imax,size(rhodz,2),2:naero) :: ZXTP10, ZXTP1C, ZXTP1CON, zsolub
real, dimension(imax,size(rhodz,2)) :: ZZOH, ZZH2O2, ZZO3, ZZNO2
real, dimension(imax,size(rhodz,2)) :: zlwcic, ziwcic
real, dimension(imax,2:naero) :: wd
real x,pqtmst
real zlwcl, zlwcv, zhp, zqtp1, zrk, zrke
real zh_so2, zpfac, zp_so2, zf_so2, zh_h2o2, zp_h2o2, zf_h2o2
real ZRKH2O2
real ze1,ze2,ze3,zfac1,zrkfac
real zza,za21,za22,zph_o3,zf_o3,zdt
real zh2o2m,zso2m,zso4m,zsumh2o2,zsumo3
real zq,zso2mh,zdso2h,zso2l,zso4l
real zzb,zzp,zzq,zzp2,zqhp,za2,zheneff
real zrko3,zso2mo,zdso2o,zdso2tot,zfac
real zxtp1dms,zso2,ztk23b
real zhil,zexp,zm,zdms,t,ztk1,zqt,zqt3
real zrhoair,zkno2o3,zkn2o5aq,zrx1,zrx12
real zkno2no3,ztk3,ztk2,zkn2o5
real zno3,zxtp1so2


!    REACTION RATE SO2-OH
real, parameter :: ZK2I=2.0E-12
real, parameter :: ZK2=4.0E-31
real, parameter :: ZK2F=0.45

!   REACTION RATE DMS-NO3
real, parameter :: ZK3=1.9E-13

!   MOLECULAR WEIGHTS IN G
real, parameter :: ZMOLGS=32.064
real, parameter :: ZMOLGAIR=28.84

real, parameter :: ZHPBASE=2.5E-06
real, parameter :: ZE1K=1.1E-02
real, parameter :: ZE1H=2300.
real, parameter :: ZE2K=1.23
real, parameter :: ZE2H=3020.
real, parameter :: ZE3K=1.2E-02
real, parameter :: ZE3H=2010.
real, parameter :: ZQ298=1./298.
real, parameter :: ZRGAS=8.2E-02

real, parameter :: ZAVO=6.022E+23
real, parameter :: ZNAMAIR=1.E-03*ZAVO/ZMOLGAIR

real, dimension(imax,2:naero) :: ZDEPS, ZDEPR
real, dimension(imax) :: PDEP, ZMTOF, ZCLR0
real, dimension(imax) :: zilcscav, xdep, plambda
real, dimension(imax) :: zbcscav, xbcscav, zstay_t
real, dimension(imax) :: xstay, zmelt, xmelt
real, dimension(imax) :: zicscav, xicscav, ziicscav
real, dimension(imax) :: xfreeze, zfreeze, zcollefc
real, dimension(imax) :: frc, ZDXTE, zxtp1

logical, parameter :: assume_convliq = .true. ! assume convective rainfall is liquid

! Below-cloud collection eff. for rain
real, dimension(naero), parameter :: zcollefr = (/0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.10,0.20,0.40,0.05,0.10/)
! Below-cloud collection eff. for snow
real, dimension(naero), parameter :: zcollefs = (/0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.02,0.04,0.08,0.01,0.02/)
! Retention coeff. on riming
real, dimension(naero), parameter :: Rcoeff = (/1.00,0.62,1.00,0.00,1.00,0.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00/)
! Allow in-cloud scavenging in ice clouds for hydrophobic BC and OC, and dust
real, dimension(naero), parameter :: Ecols = (/0.00,0.00,0.00,0.05,0.00,0.05,0.00,0.05,0.05,0.05,0.05,0.05,0.05/)
! wet deposition coefficients
! Relative re-evaporation rate
real, dimension(naero), parameter :: Evfac = (/0.25,1.00,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25/)


! Start code : ----------------------------------------------------------

ifull = size(rhodz,1)
kl    = size(rhodz,2)
ntiles = ifull/imax

!    CONSTANTS
PQTMST=1./PTMST

if ( maxval(xtm1(1:ifull,:,:))>6.5e-5 ) then
  write(6,*) "xtg is out-of-range at start of xtchemie"
  write(6,*) "xtg maxval,maxloc ",maxval(xtm1(1:ifull,:,:)), &
                                  maxloc(xtm1(1:ifull,:,:))
end if

!$acc data create(prhop1,ptp1,rhodz,xtm1,xte,zoxidant)
!$acc update device(prhop1,ptp1,rhodz,xtm1,zoxidant)

!$acc parallel loop copy(so2wd,so4wd,bcwd,ocwd,dustwd,saltwd)                             &
!$acc   copyin(pclcover,pcfcover,pmlwc,pmiwc,pccw,xtu)                                    &
!$acc   copyin(pclcon,pqfsedice,pmaccr,plambs,pfsnow,pfsubl,pfmelt,pmratep,prscav)        &
!$acc   copyin(pfevap,pfprec,prfreeze,fracc,pfconv)                                       &
!$acc   copyout(so2h2,so2o3)                                                              &
!$acc   private(xto,zlwcic,ziwcic,zhenry,zhenryc,zso4,zso4c,zso4i,zxtp10,zxtp1c,zxtp1con) &
!$acc   private(zzh2o2,zzo3,zsolub,wd,zdepr,zdeps,pdep,zclr0,zmtof,zilcscav)              &
!$acc   private(xdep,plambda,zbcscav,xbcscav,xstay,zmelt,xmelt,zicscav,xicscav,ziicscav)  &
!$acc   private(xfreeze,zfreeze,zcollefc,frc,zdxte,zxtp1)                                 &
!$acc   present(prhop1,ptp1,rhodz,xtm1,xte,zoxidant)
do tile = 1,ntiles
  js = (tile-1)*imax + 1
  je = tile*imax
  
  !$acc loop vector
  do jl = js,je
    so2h2(jl)=0.
    so2o3(jl)=0.
  end do

  !$acc loop collapse(2) seq
  do jt = 1,naero
    do jk = 1,kl
      !$acc loop vector  
      do jl = 1,imax
        iq = jl + js - 1
        xte(iq,jk,jt)=0.
        ! Calculate xto, tracer mixing ratio outside convective updraughts
        ! Assumes pclcon < 1, but this shouldn't be a problem.
        xto(jl,jk,jt)=(xtm1(iq,jk,jt)-pclcon(iq,jk)*xtu(iq,jk,jt))/(1.-pclcon(iq,jk))
        xto(jl,jk,jt)=max(0.,xto(jl,jk,jt))
      end do  
    end do
  end do

  !$acc loop seq
  do jk = 1,kl
    !$acc loop vector  
    do jl = 1,imax
      iq = jl + js - 1  

      ! Calculate in-cloud ql
      zlwcic(jl,jk)=0.
      if ( pclcover(iq,jk)>1.e-8 ) then
        zlwcic(jl,jk)=pmlwc(iq,jk)/pclcover(iq,jk)
      end if
      ziwcic(jl,jk)=0.
      if ( pcfcover(iq,jk)>1.e-8 ) then
        ziwcic(jl,jk)=pmiwc(iq,jk)/pcfcover(iq,jk)
      end if

      !  OXIDANT CONCENTRATIONS IN MOLECULE/CM**3
      ! -- levels are already inverted --
      ZZH2O2(jl,jk) = ZOXIDANT(iq,jk,2)*PRHOP1(iq,jk)*1.e-3
      ZZO3(jl,jk)   = ZOXIDANT(iq,jk,3)*PRHOP1(iq,jk)*1.e-3
  
      zhenry(jl,jk)=0.
      zhenryc(jl,jk)=0.

      !   PROCESSES WHICH ARE DIFERENT INSIDE AND OUTSIDE OF CLOUDS
      ZSO4(jl,jk)=max(XTO(jl,jk,ITRACSO4),0.)
      ZSO4i(jl,jk)=max(XTO(jl,jk,ITRACSO4),0.)
      ZXTP1CON(jl,jk,ITRACSO2)=max(XTU(iq,jk,ITRACSO2),0.)
      ZSO4C(jl,jk)            =max(XTU(iq,jk,ITRACSO4),0.)
    end do  
  end do

  !$acc loop seq
  do jk = ktop,kl
    !$acc loop vector  
    do jl = 1,imax
      iq = jl + js - 1
      !   CALCULATE THE REACTION-RATES FOR SO2-H2O2
      if ( zlwcic(jl,jk)>zmin ) then
        ZLWCL=ZLWCIC(JL,JK)*PRHOP1(iq,JK)*1.E-06
        ZLWCV=ZLWCIC(JL,JK)*PRHOP1(iq,JK)*1.E-03
        ZHP=ZHPBASE+ZSO4(JL,JK)*1000./(ZLWCIC(JL,JK)*ZMOLGS)
        ZQTP1=1./PTP1(iq,JK)-ZQ298
        ZRK=8.E+04*EXP(-3650.*ZQTP1)/(0.1+ZHP)
        ZRKE=ZRK/(ZLWCL*ZAVO)

        ZH_SO2=ZE2K*EXP(ZE2H*ZQTP1)
        ZPFAC=ZRGAS*ZLWCV*PTP1(iq,JK)
        ZP_SO2=ZH_SO2*ZPFAC
        ZF_SO2=ZP_SO2/(1.+ZP_SO2)

        ZH_H2O2=9.7E+04*EXP(6600.*ZQTP1)
        ZP_H2O2=ZH_H2O2*ZPFAC
        ZF_H2O2=ZP_H2O2/(1.+ZP_H2O2)

        ZRKH2O2=ZRKE*ZF_SO2*ZF_H2O2
      else
        ZRKH2O2=0.
      end if

      !   HETEROGENEOUS CHEMISTRY
      ZXTP1(jl)              = XTO(JL,JK,ITRACSO2)
      ZXTP10(JL,JK,ITRACSO2) = XTO(JL,JK,ITRACSO2)
      ZXTP1C(JL,JK,ITRACSO2) = XTO(JL,JK,ITRACSO2)
      IF ( ZXTP1(jl)>ZMIN .AND. ZLWCIC(JL,JK)>ZMIN ) THEN

        ZQTP1=1./PTP1(iq,JK)-ZQ298
        ZE1=ZE1K*EXP(ZE1H*ZQTP1)
        ZE2=ZE2K*EXP(ZE2H*ZQTP1)
        ZE3=ZE3K*EXP(ZE3H*ZQTP1)

        ZLWCL=ZLWCIC(JL,JK)*PRHOP1(iq,JK)*1.E-06
        !    ZLWCL = LWC IN L/CM**3
        ZLWCV=ZLWCIC(JL,JK)*PRHOP1(iq,JK)*1.E-03
        !   ZLWCV = LWC IN VOL/VOL
        ZFAC1=1./(ZLWCL*ZAVO)
        !   ZFAC1 CALCULATES MOLECULES PER CM**3 TO MOLE PER LTR H2O
        ZRKFAC=ZRGAS*PTP1(iq,JK)*ZLWCV
        !   ZRKFAC CALCULATES DIMENSIONLESS HENRY-COEFF.
        ZZA=ZE2*ZRKFAC
        ZA21=4.39E+11*EXP(-4131./PTP1(iq,JK))
        ZA22=2.56E+03*EXP(-966./PTP1(iq,JK)) !926 corrected to 966 here
        ZPH_O3=ZE1*ZRKFAC
        ZF_O3=ZPH_O3/(1.+ZPH_O3)
        ZDT=PTMST/5.

        ZH2O2M=ZZH2O2(JL,JK)
        ZSO2M=ZXTP1(jl)*PRHOP1(iq,JK)*6.022E+20/ZMOLGS
        ZSO4M=ZSO4(JL,JK)*PRHOP1(iq,JK)*6.022E+20/ZMOLGS

        ZSUMH2O2=0.
        ZSUMO3=0.

        !$acc loop seq
        DO JN=1,5
          ZQ=ZRKH2O2*ZH2O2M
          ZSO2MH=ZSO2M*EXP(-ZQ*ZDT)

          ZDSO2H=ZSO2M-ZSO2MH
          ZH2O2M=ZH2O2M-ZDSO2H
          ZH2O2M=MAX(0.,ZH2O2M)
          ZSUMH2O2=ZSUMH2O2+ZDSO2H

          ZSO4M=ZSO4M+ZDSO2H
          !   CALCULATE THE PH VALUE
          ZSO2L=ZSO2MH*ZFAC1
          ZSO4L=ZSO4M*ZFAC1
          ZZB=ZHPBASE+ZSO4L
          ZZP=(ZZA*ZE3-ZZB-ZZA*ZZB)/(1.+ZZA)
          ZZQ=-ZZA*ZE3*(ZZB+ZSO2L)/(1.+ZZA)
          ZZP=0.5*ZZP
          ZZP2=ZZP*ZZP
          ZHP=-ZZP+SQRT(ZZP2-ZZQ)
          ZQHP=1./ZHP

          !   CALCULATE THE REACTION RATE FOR SO2-O3
          ZA2=(ZA21+ZA22*ZQHP)*ZFAC1
          ZHENEFF=1.+ZE3*ZQHP
          ZP_SO2=ZZA*ZHENEFF
          ZF_SO2=ZP_SO2/(1.+ZP_SO2)
          ZRKO3=ZA2*ZF_O3*ZF_SO2

          ZQ=ZZO3(JL,JK)*ZRKO3
          ZSO2MO=ZSO2MH*EXP(-ZQ*ZDT)
          ZDSO2O=ZSO2MH-ZSO2MO
          ZSO4M=ZSO4M+ZDSO2O
          ZSO2M=ZSO2MO
          ZSUMO3=ZSUMO3+ZDSO2O
        end do  !End of iteration loop

        ZDSO2TOT=ZXTP1(jl)-ZSO2M*ZMOLGS/(6.022E+20*PRHOP1(iq,JK))
        ZDSO2TOT=MIN(ZDSO2TOT,ZXTP1(jl))
        ZXTP1C(JL,JK,ITRACSO2)=ZXTP1(jl)-ZDSO2TOT
        ZSO4(JL,JK)=ZSO4(JL,JK)+ZDSO2TOT

        ZHENRY(JL,JK)=ZF_SO2
        ! Diagnostic only...
        ZFAC=PQTMST*PCLCOVER(iq,JK)*ZMOLGS/(6.022E+20*PRHOP1(iq,JK))
        ZFAC1=ZFAC*rhodz(iq,JK)
        so2h2(iq)=so2h2(iq)+ZSUMH2O2*ZFAC1
        so2o3(iq)=so2o3(iq)+ZSUMO3*ZFAC1
      END IF
    end do
  end do


  ! Repeat the aqueous oxidation calculation for ice clouds.

!******************************************************************************
!   CALCULATE THE REACTION-RATES FOR SO2-H2O2
!DO JK=KTOP,KL
!  DO JL=1,ifull
!    IF(ziwcic(JL,JK).GT.ZMIN) THEN
!      ZLWCL(jl)=ziwcic(JL,JK)*PRHOP1(JL,JK)*1.E-06
!      ZLWCV(jl)=ziwcic(JL,JK)*PRHOP1(JL,JK)*1.E-03
!      ZHP(jl)=ZHPBASE+ZSO4i(JL,JK)*1000./(ziwcic(JL,JK)*ZMOLGS)
!      ZQTP1(jl)=1./PTP1(JL,JK)-ZQ298
!      ZRK(jl)=8.E+04*EXP(-3650.*ZQTP1(jl))/(0.1+ZHP(jl))
!      ZRKE(jl)=ZRK(jl)/(ZLWCL(jl)*ZAVO)
!
!      ZH_SO2(jl)=ZE2*EXP(ZE2H*ZQTP1(jl))
!      ZPFAC(jl)=ZRGAS*ZLWCV(jl)*PTP1(JL,JK)
!      ZP_SO2(jl)=ZH_SO2(jl)*ZPFAC(jl)
!      ZF_SO2(jl)=ZP_SO2(jl)/(1.+ZP_SO2(jl))
!
!      ZH_H2O2(jl)=9.7E+04*EXP(6600.*ZQTP1(jl))
!      ZP_H2O2(jl)=ZH_H2O2(jl)*ZPFAC(jl)
!      ZF_H2O2(jl)=ZP_H2O2(jl)/(1.+ZP_H2O2(jl))
!
!      ZRKH2O2(JL,JK)=ZRKE(jl)*ZF_SO2(jl)*ZF_H2O2(jl)
!    ELSE
!      ZRKH2O2(JL,JK)=0.
!    ENDIF
!  ENDDO
!
!!   HETEROGENEOUS CHEMISTRY
!  DO JL=1,ifull
!    ZXTP1(jl)=XTO(JL,JK,ITRACSO2)
!    IF(ZXTP1(jl)>ZMIN.AND.ziwcic(JL,JK)>ZMIN) THEN
!      X=PRHOP1(JL,JK)
!
!      ZQTP1(jl)=1./PTP1(JL,JK)-ZQ298
!      ZE1=ZE1K*EXP(ZE1H*ZQTP1(jl))
!      ZE2=ZE2K*EXP(ZE2H*ZQTP1(jl))
!      ZE3=ZE3K*EXP(ZE3H*ZQTP1(jl))
!
!      ZLWCL(jl)=ziwcic(JL,JK)*PRHOP1(JL,JK)*1.E-06
!!    ZLWCL = LWC IN L/CM**3
!      ZLWCV(jl)=ziwcic(JL,JK)*PRHOP1(JL,JK)*1.E-03
!!   ZLWCV = LWC IN VOL/VOL
!      ZFAC1=1./(ZLWCL(jl)*ZAVO)
!!   ZFAC1 CALCULATES MOLECULES PER CM**3 TO MOLE PER LTR H2O
!      ZRKFAC=ZRGAS*PTP1(JL,JK)*ZLWCV(jl)
!!   ZRKFAC CALCULATES DIMENSIONLESS HENRY-COEFF.
!      ZZA=ZE2*ZRKFAC
!      ZA21=4.39E+11*EXP(-4131./PTP1(JL,JK))
!      ZA22=2.56E+03*EXP(-966./PTP1(JL,JK)) !926 corrected to 966 here
!      ZPH_O3=ZE1*ZRKFAC
!      ZF_O3=ZPH_O3/(1.+ZPH_O3)
!      ZDT=PTMST/5.
!
!      ZH2O2M=ZZH2O2(JL,JK)
!      ZSO2M=ZXTP1(jl)*X*6.022E+20/ZMOLGS
!      ZSO4M=ZSO4i(JL,JK)*X*6.022E+20/ZMOLGS
!
!      ZSUMH2O2=0.
!      ZSUMO3=0.
!
!      DO JN=1,5
!        ZQ=ZRKH2O2(JL,JK)*ZH2O2M
!        ZSO2MH=ZSO2M*EXP(-ZQ*ZDT)
!
!        ZDSO2H=ZSO2M-ZSO2MH
!        ZH2O2M=ZH2O2M-ZDSO2H
!        ZH2O2M=MAX(0.,ZH2O2M)
!        ZSUMH2O2=ZSUMH2O2+ZDSO2H
!
!        ZSO4M=ZSO4M+ZDSO2H
!!   CALCULATE THE PH VALUE
!        ZSO2L=ZSO2MH*ZFAC1
!        ZSO4L=ZSO4M*ZFAC1
!        ZZB=ZHPBASE+ZSO4L
!       ZZP=(ZZA*ZE3-ZZB-ZZA*ZZB)/(1.+ZZA)
!        ZZQ=-ZZA*ZE3*(ZZB+ZSO2L)/(1.+ZZA)
!        ZZP=0.5*ZZP
!        ZZP2=ZZP*ZZP
!        ZHP(jl)=-ZZP+SQRT(ZZP2-ZZQ)
!        ZQHP=1./ZHP(jl)
!
!!   CALCULATE THE REACTION RATE FOR SO2-O3
!        ZA2=(ZA21+ZA22*ZQHP)*ZFAC1
!        ZHENEFF=1.+ZE3*ZQHP
!        ZP_SO2(jl)=ZZA*ZHENEFF
!        ZF_SO2(jl)=ZP_SO2(jl)/(1.+ZP_SO2(jl))
!        ZRKO3=ZA2*ZF_O3*ZF_SO2(jl)
!
!        ZQ=ZZO3(JL,JK)*ZRKO3
!        ZSO2MO=ZSO2MH*EXP(-ZQ*ZDT)
!        ZDSO2O=ZSO2MH-ZSO2MO
!        ZSO4M=ZSO4M+ZDSO2O
!        ZSO2M=ZSO2MO
!        ZSUMO3=ZSUMO3+ZDSO2O
!      ENDDO  !End of iteration loop
!
!      ZDSO2TOT=ZXTP1(jl)-ZSO2M*ZMOLGS/(6.022E+20*X)
!      ZDSO2TOT=MIN(ZDSO2TOT,ZXTP1(jl))
!
!      ZXTP10(JL,JK,ITRACSO2)=ZXTP1(jl)-ZDSO2TOT*pcfcover(jl,jk)/(1.-pclcover(jl,jk))
!      ZSO4i(JL,JK)=ZSO4i(JL,JK)+ZDSO2TOT*pcfcover(jl,jk)/(1.-pclcover(jl,jk))
!      ZHENRY(JL,JK)=ZF_SO2(jl)
!! Diagnostic only...
!      ZFAC=PQTMST*pcfcover(jl,jk)*ZMOLGS/(6.022E+20*X)
!      ZFAC1=ZFAC*rhodz(JL,JK)
!      so2h2(JL)=so2h2(JL)+ZSUMH2O2*ZFAC1
!      so2o3(JL)=so2o3(JL)+ZSUMO3*ZFAC1
!    ENDIF
!  ENDDO
!ENDDO
!******************************************************************************


  ! Repeat the aqueous oxidation calculation for convective clouds.
  
  !$acc loop seq
  do jk = ktop,kl
    !$acc loop vector
    do jl = 1,imax
      iq = jl + js - 1
      !   CALCULATE THE REACTION-RATES FOR SO2-H2O2
      if ( PCCW(iq,JK)>ZMIN ) then
        ZLWCL=PCCW(iq,JK)*PRHOP1(iq,JK)*1.E-06
        ZLWCV=PCCW(iq,JK)*PRHOP1(iq,JK)*1.E-03
        ZHP=ZHPBASE+ZSO4C(JL,JK)*1000./(PCCW(iq,JK)*ZMOLGS)
        ZQTP1=1./PTP1(iq,JK)-ZQ298
        ZRK=8.E+04*EXP(-3650.*ZQTP1)/(0.1+ZHP)
        ZRKE=ZRK/(ZLWCL*ZAVO)

        ZH_SO2=ZE2K*EXP(ZE2H*ZQTP1)
        ZPFAC=ZRGAS*ZLWCV*PTP1(iq,JK)
        ZP_SO2=ZH_SO2*ZPFAC
        ZF_SO2=ZP_SO2/(1.+ZP_SO2)

        ZH_H2O2=9.7E+04*EXP(6600.*ZQTP1)
        ZP_H2O2=ZH_H2O2*ZPFAC
        ZF_H2O2=ZP_H2O2/(1.+ZP_H2O2)

        ZRKH2O2=ZRKE*ZF_SO2*ZF_H2O2
      else
        ZRKH2O2=0.
      end if

      !   HETEROGENEOUS CHEMISTRY
      ZXTP1(jl)=XTU(iq,JK,ITRACSO2)
      IF(ZXTP1(jl)>ZMIN.AND.PCCW(iq,JK)>ZMIN) THEN

        ZQTP1=1./PTP1(iq,JK)-ZQ298
        ZE1=ZE1K*EXP(ZE1H*ZQTP1)
        ZE2=ZE2K*EXP(ZE2H*ZQTP1)
        ZE3=ZE3K*EXP(ZE3H*ZQTP1)

        ZLWCL=PCCW(iq,JK)*PRHOP1(iq,JK)*1.E-06
        !    ZLWCL = LWC IN L/CM**3
        ZLWCV=PCCW(iq,JK)*PRHOP1(iq,JK)*1.E-03
        !   ZLWCV = LWC IN VOL/VOL
        ZFAC1=1./(ZLWCL*ZAVO)
        !   ZFAC1 CALCULATES MOLECULES PER CM**3 TO MOLE PER LTR H2O
        ZRKFAC=ZRGAS*PTP1(iq,JK)*ZLWCV
        !   ZRKFAC CALCULATES DIMENSIONLESS HENRY-COEFF.
        ZZA=ZE2*ZRKFAC
        ZA21=4.39E+11*EXP(-4131./PTP1(iq,JK))
        ZA22=2.56E+03*EXP(-966./PTP1(iq,JK)) !926 corrected to 966 here
        ZPH_O3=ZE1*ZRKFAC
        ZF_O3=ZPH_O3/(1.+ZPH_O3)
        ZDT=PTMST/5.

        ZH2O2M=ZZH2O2(JL,JK)
        ZSO2M=ZXTP1(jl)*PRHOP1(iq,JK)*6.022E+20/ZMOLGS
        ZSO4M=ZSO4C(JL,JK)*PRHOP1(iq,JK)*6.022E+20/ZMOLGS

        ZSUMH2O2=0.
        ZSUMO3=0.

        !$acc loop seq
        DO JN=1,5
          ZQ=ZRKH2O2*ZH2O2M
          ZSO2MH=ZSO2M*EXP(-ZQ*ZDT)

          ZDSO2H=ZSO2M-ZSO2MH
          ZH2O2M=ZH2O2M-ZDSO2H
          ZH2O2M=MAX(0.,ZH2O2M)
          ZSUMH2O2=ZSUMH2O2+ZDSO2H

          ZSO4M=ZSO4M+ZDSO2H
!   CALCULATE THE PH VALUE
          ZSO2L=ZSO2MH*ZFAC1
          ZSO4L=ZSO4M*ZFAC1
          ZZB=ZHPBASE+ZSO4L
          ZZP=(ZZA*ZE3-ZZB-ZZA*ZZB)/(1.+ZZA)
          ZZQ=-ZZA*ZE3*(ZZB+ZSO2L)/(1.+ZZA)
          ZZP=0.5*ZZP
          ZZP2=ZZP*ZZP
          ZHP=-ZZP+SQRT(ZZP2-ZZQ)
          ZQHP=1./ZHP

!   CALCULATE THE REACTION RATE FOR SO2-O3
          ZA2=(ZA21+ZA22*ZQHP)*ZFAC1
          ZHENEFF=1.+ZE3*ZQHP
          ZP_SO2=ZZA*ZHENEFF
          ZF_SO2=ZP_SO2/(1.+ZP_SO2)
          ZRKO3=ZA2*ZF_O3*ZF_SO2
!
          ZQ=ZZO3(JL,JK)*ZRKO3
          ZSO2MO=ZSO2MH*EXP(-ZQ*ZDT)
          ZDSO2O=ZSO2MH-ZSO2MO
          ZSO4M=ZSO4M+ZDSO2O
          ZSO2M=ZSO2MO
          ZSUMO3=ZSUMO3+ZDSO2O
        END DO  !End of iteration loop

        ZDSO2TOT=ZXTP1(jl)-ZSO2M*ZMOLGS/(6.022E+20*PRHOP1(iq,JK))
        ZDSO2TOT=MIN(ZDSO2TOT,ZXTP1(jl))
        ZXTP1CON(JL,JK,ITRACSO2)=ZXTP1CON(JL,JK,ITRACSO2)-ZDSO2TOT
        ZSO4C(JL,JK)=ZSO4C(JL,JK)+ZDSO2TOT
        ZHENRYC(JL,JK)=ZF_SO2
        ! Diagnostic only...
        ZFAC=PQTMST*pclcon(iq,jk)*ZMOLGS/(6.022E+20*PRHOP1(iq,JK))
        ZFAC1=ZFAC*rhodz(iq,JK)
        so2h2(iq)=so2h2(iq)+ZSUMH2O2*ZFAC1
        so2o3(iq)=so2o3(iq)+ZSUMO3*ZFAC1
      END IF
    END DO
  END DO

  !*******************************************************************************
  !
  !    CALCULATE THE WET DEPOSITION
  !    (True for all except DMS)
  !

  !$acc loop seq
  do jk = 1,kl
    !$acc loop vector  
    do jl = 1,imax
      iq = jl + js - 1  
      zsolub(jl,jk,ITRACSO2)=zhenry(jl,jk)
      zxtp10(jl,jk,ITRACSO4)=zso4i(jl,jk)
      zxtp1c(jl,jk,ITRACSO4)=zso4(jl,jk)    
      zxtp1con(jl,jk,ITRACSO4)=zso4c(jl,jk)
      zsolub(jl,jk,ITRACSO4)=0.6
      zxtp10(jl,jk,ITRACbc)=xto(jl,jk,ITRACbc)
      zxtp1c(jl,jk,ITRACbc)=xto(jl,jk,ITRACbc)
      zxtp1con(jl,jk,ITRACbc)=xtu(iq,jk,ITRACbc)
      zsolub(jl,jk,ITRACbc)=0.
      zxtp10(jl,jk,ITRACBC+1)=xto(jl,jk,ITRACBC+1)
      zxtp1c(jl,jk,ITRACBC+1)=xto(jl,jk,ITRACBC+1)
      zxtp1con(jl,jk,ITRACBC+1)=xtu(iq,jk,ITRACBC+1)
      zsolub(jl,jk,ITRACBC+1)=0.2
      zxtp10(jl,jk,ITRACOC)=xto(jl,jk,ITRACOC)
      zxtp1c(jl,jk,ITRACOC)=xto(jl,jk,ITRACOC)
      zxtp1con(jl,jk,ITRACOC)=xtu(iq,jk,ITRACOC)
      zsolub(jl,jk,ITRACOC)=0.
      zxtp10(jl,jk,ITRACOC+1)=xto(jl,jk,ITRACOC+1)
      zxtp1c(jl,jk,ITRACOC+1)=xto(jl,jk,ITRACOC+1)
      zxtp1con(jl,jk,ITRACOC+1)=xtu(iq,jk,ITRACOC+1)
      zsolub(jl,jk,ITRACOC+1)=0.2
    end do
  end do
  !$acc loop vector
  do jl = 1,imax
    wd(jl,ITRACSO2) = 0.
    wd(jl,ITRACSO4) = 0.
    wd(jl,ITRACbc) = 0.
    wd(jl,ITRACbc+1) = 0.
    wd(jl,ITRACOC) = 0.
    wd(jl,ITRACOC+1) = 0.
  end do  
  !$acc loop seq
  DO JT=ITRACDU,ITRACDU+NDUST-1
    !$acc loop seq
    do jk = 1,kl
      !$acc loop vector  
      do jl = 1,imax
        iq = jl + js - 1  
        zxtp10(jl,jk,jt)=xto(jl,jk,jt)
        zxtp1c(jl,jk,jt)=xto(jl,jk,jt)
        zxtp1con(jl,jk,jt)=xtu(iq,jk,jt)
        zsolub(jl,jk,jt)=0.05
      end do
    end do
    !$acc loop vector
    do jl = 1,imax
      wd(jl,jt) = 0.
    end do
  end do
  !$acc loop seq
  DO JT=ITRACSA,ITRACSA+NSALT-1
    !$acc loop seq  
    do jk = 1,kl
      !$acc loop vector
      do jl = 1,imax
        iq = jl + js - 1  
        zxtp10(jl,jk,jt)=xto(jl,jk,jt)
        zxtp1c(jl,jk,jt)=xto(jl,jk,jt)
        zxtp1con(jl,jk,jt)=xtu(iq,jk,jt)
        zsolub(jl,jk,jt)=0.05
      end do
    end do
    !$acc loop vector
    do jl = 1,imax
      wd(jl,jt) = 0.
    end do
  end do

  !$acc loop seq
  do ktrac = 2,naero
    !$acc loop vector  
    do i = 1,imax
      zdepr(i,ktrac) = 0.
      zdeps(i,ktrac) = 0.
    end do
  end do

  ! Search for convective cloud base
  !kbase(:) = kl+1
  !do jk = ktop,kl
  !  where ( pclcon(:,jk)>zmin )
  !    kbase(:) = k
  !  end where
  !enddo

  !     BEGIN OF VERTICAL LOOP
  !$acc loop seq
  do JK = KTOP,kl
    do ktrac = 2,naero
        
      pdep(1:imax) = 0.  
      
      ! zdepr(i) = zdepr(i) + zdepr_save(i,jk)
      ! zdeps(i) = zdeps(i) + zdeps_save(i,jk)
      ! zdepr_save(i,jk) = 0.
      ! zdeps_save(i,jk) = 0.

      !ZCLEAR = 1. - PCLCOVER(i,JK) - pcfcover(i,jk) - pclcon(i,jk)
      ZCLR0(:) = max( 1. - PCLCOVER(js:je,jk) - pclcon(js:je,jk), 0. ) !Clear air or ice cloud (applies to zxtp10)
      ZMTOF(:) = rhodz(js:je,jk)*pqtmst
      ZXTP1C(:,JK,ktrac) = MAX( 0., ZXTP1C(:,JK,ktrac) )
      ZXTP10(:,JK,ktrac) = MAX( 0., ZXTP10(:,JK,ktrac) )

      ! In-cloud ice scavenging (vertical redistribution when snow falls into a layer).
      where ( zclr0(:)>zmin )
        ziicscav(:) = Ecols(ktrac)*pqfsedice(js:je,jk) !qfsedice is the fractional sedimentation in dt
        ziicscav(:) = max( min( ziicscav(:), 1. ), 0. )
        xdep(:) = max( zxtp10(1:imax,jk,ktrac)*ziicscav(:), 0.)
        pdep(:) = pdep(:) + xdep(:)*pcfcover(js:je,jk)
        zxtp10(:,jk,ktrac) = zxtp10(:,jk,ktrac) - xdep(:)*pcfcover(js:je,jk)/zclr0(:) ! MJT suggestion
        zdeps(:,ktrac) = zdeps(:,ktrac) + xdep(:)*pcfcover(js:je,jk)*zmtof(:)
      end where  

      ! This loop does riming (accretion of liquid water by falling snow)
      where ( pmlwc(js:je,jk)>zmin )
        zilcscav(:) = Rcoeff(ktrac)*zsolub(:,jk,ktrac)*pmaccr(js:je,jk)*ptmst/pmlwc(js:je,jk)
        zilcscav(:) = max( min( zilcscav(:), 1. ), 0. )
        xdep(:) = max( zxtp1c(:,jk,ktrac)*zilcscav, 0. )
        pdep(:) = pdep(:) + xdep(:)*pclcover(js:Je,jk)
        zxtp1c(:,jk,ktrac) = zxtp1c(:,jk,ktrac) - xdep(:)
        zdeps(:,ktrac) = zdeps(:,ktrac) + xdep(:)*pclcover(js:Je,jk)*zmtof(:)
      end where

      ! Below-cloud scavenging by snow
      plambda(:) = min( plambs(js:Je,jk), 8.e3 ) !Cut it off at about -30 deg. C
      zbcscav(:) = zcollefs(ktrac)*plambda(:)*pfsnow(js:je,jk)*ptmst/(2.*rhos)
      zbcscav(:) = max( min( 1., zbcscav(:)/(1.+0.5*zbcscav(:)) ), 0. ) !Time-centred
      xbcscav(:) = max( zbcscav(:)*zxtp10(:,jk,ktrac), 0. )
      pdep(:) = pdep(:) + xbcscav(:)*zclr0(:)
      zxtp10(:,jk,ktrac) = zxtp10(:,jk,ktrac) - xbcscav(:)
      zdeps(:,ktrac) = zdeps(:,ktrac) + xbcscav(:)*zclr0(:)*zmtof(:)

      ! Redistribution by snow that evaporates
      where ( pfsubl(js:je,jk)>zmin .and. pfsnow(js:je,jk)>zmin .and. zclr0(:)>zmin )
        !zstay_t(:) = (pfsubl(js:je,jk)+pfstayice(js:je,jk))/pfsnow(js:je,jk)
        zstay_t(:) = pfsubl(js:je,jk)/pfsnow(js:je,jk) ! MJT suggestion
        zstay_t(:) = max( min( 1., zstay_t(:) ), 0. )
        xstay(:) = max( zdeps(:,ktrac)*zstay_t(:)/zmtof(:), 0. )
        !limit sublimation to prevent crash - MJT suggestion
        xstay(:) = max( min( xstay(:), 6.e-6/(1.-pclcover(js:Je,jk)-pclcon(js:je,jk)) - zxtp10(:,jk,ktrac) ), 0. )
        pdep(:) = pdep(:) - xstay(:)*zclr0(:)
        zxtp10(:,jk,ktrac) = zxtp10(:,jk,ktrac) + xstay(:)
        zdeps(:,ktrac) = zdeps(:,ktrac) - xstay(:)*zclr0(:)*zmtof(:)
        zdeps(:,ktrac) = max( 0., zdeps(:,ktrac) )
      end where
    
      ! Redistribution by snow that stays in layer
      !if ( pfstayice(i,jk)>zmin .and. zclr0>zmin ) then
      !  zstay_t = pfstayice(i,jk)/(pfsnow(i,jk)+pfstayice(i,jk))
      !  zstay_t = max( min( 1., zstay_t ), 0. )
      !  xstay = max( zdeps(i,ktrac)*zstay_t/zmtof, 0. )
      !  pdep = pdep - xstay*zclr0
      !  zdeps_save(i,jk,ktrac) = zdeps_save(i,jk,ktrac) + xstay*zclr0*zmtof
      !  zdeps(i,ktrac) = zdeps(i,ktrac) - xstay*zclr0*zmtof
      !  zdeps(i,ktrac) = max( 0., zdeps(i,ktrac) )
      !end if    

      ! Melting of snow... 
      zmelt(:) = pfmelt(js:je,jk)/max(pfsnow(js:je,jk)+pfmelt(js:je,jk),zmin) 
      zmelt(:) = max( min( 1., zmelt(:) ), 0. )
      xmelt(:) = zmelt(:)*zdeps(:,ktrac)
      zdepr(:,ktrac) = zdepr(:,ktrac) + xmelt(:)
      zdeps(:,ktrac) = zdeps(:,ktrac) - xmelt(:)
      zdeps(:,ktrac) = max( 0., zdeps(:,ktrac) )
  
      !  In-cloud scavenging by warm-rain processes (autoconversion and collection)
      where ( pmlwc(js:je,jk)>zmin ) ! MJT suggestion
        zicscav(:) = zsolub(:,jk,ktrac)*pmratep(js:Je,jk)*ptmst/pmlwc(js:je,jk)
        zicscav(:) = max( min( zicscav(:), 1. ), 0. )
        xicscav(:) = max( zxtp1c(:,jk,ktrac)*zicscav(:), 0. )
        pdep(:) = pdep(:) + xicscav(:)*pclcover(js:je,jk)
        zxtp1c(:,jk,ktrac) = zxtp1c(:,jk,ktrac) - xicscav(:)
        zdepr(:,ktrac) = zdepr(:,ktrac) + xicscav(:)*pclcover(js:je,jk)*zmtof(:)
      end where
 
      ! Below-cloud scavenging by stratiform rain (conv done below)
      zbcscav(:) = zcollefr(ktrac)*prscav(js:je,jk)
      zbcscav(:) = max( min( 1., zbcscav(:)/(1.+0.5*zbcscav(:)) ), 0. ) !Time-centred
      xbcscav(:) = max( zbcscav(:)*zxtp10(:,jk,ktrac), 0. )
      pdep(:) = pdep(:) + xbcscav(:)*zclr0(:)
      zxtp10(:,jk,ktrac) = zxtp10(:,jk,ktrac) - xbcscav(:) 
      zdepr(:,ktrac) = zdepr(:,ktrac) + xbcscav(:)*zclr0(:)*zmtof(:)

      ! MJT - suggestion (only include evaporation)
      ! Redistribution by rain that evaporates or stays in layer
      where ( pfevap(js:Je,jk)>zmin .and. pfprec(js:je,jk)>zmin .and. zclr0(:)>zmin )
        !zstay_t(:) = (pfevap(js:Je,jk)+pfstayliq(js:je,jk))/pfprec(js:je,jk)  
        zstay_t(:) = pfevap(js:je,jk)/pfprec(js:je,jk) ! MJT suggestion
        zstay_t(:) = max( min( 1., zstay_t(:) ), 0. )
        xstay(:) = max( zdepr(:,ktrac)*zstay_t(:)*evfac(ktrac)/zmtof(:), 0. )
        !limit sublimation to prevent crash - MJT suggestion
        xstay(:) = max( min( xstay(:), 6.e-6/(1.-pclcover(js:je,jk)-pclcon(js:Je,jk)) - zxtp10(:,jk,ktrac) ), 0. )
        pdep(:) = pdep(:) - xstay(:)*zclr0(:)
        zxtp10(:,jk,ktrac) = zxtp10(:,jk,ktrac) + xstay(:)
        zdepr(:,ktrac) = zdepr(:,ktrac) - xstay(:)*zclr0(:)*zmtof(:)
        zdepr(:,ktrac) = max( 0., zdepr(:,ktrac) )
      end where
    
      ! Redistribution by rain that evaporates or stays in layer
      !if ( pfstayliq(i,jk)>zmin .and. zclr0>zmin ) then
      !  zstay_t = pfstayliq(i,jk)/(pfprec(i,jk)+pfstayliq(i,jk))
      !  zstay_t = max( min( 1., zstay_t ), 0. )
      !  xstay = max( zdepr(i,ktrac)*zstay_t/zmtof, 0. )
      !  pdep = pdep - xstay*zclr0
      !  zdepr_save(i,jk,ktrac) = zdepr_save(i,jk,ktrac) + xstay*zclr0*zmtof
      !  zdepr(i,ktrac) = zdepr(i,ktrac) - xstay*zclr0*zmtof
      !  zdepr(i,ktrac) = max( 0., zdepr(i,ktrac) )
      !end if

      ! Freezing of rain... 
      zfreeze(:) = prfreeze(js:Je,jk)/max(pfprec(js:Je,jk)+prfreeze(js:je,jk),zmin) 
      zfreeze(:) = max( min( 1., zfreeze(:) ), 0. )
      xfreeze(:) = zfreeze(:)*zdepr(:,ktrac)
      zdeps(:,ktrac) = zdeps(:,ktrac) + xfreeze(:)
      zdepr(:,ktrac) = zdepr(:,ktrac) - xfreeze(:)
      zdepr(:,ktrac) = max( 0., zdepr(:,ktrac) )

      ! Now do the convective below-cloud bit...
      ! In-cloud convective bit was done in convjlm.

      ! Use collection efficiencies for rain below melting level, snow above

      ! MJT notes - Assume rain for JLM convection
      where ( ptp1(js:je,jk)>273.15 .or. assume_convliq )
        zcollefc(:) = zcollefr(ktrac)
      else where
        zcollefc(:) = zcollefs(ktrac)
      end where

      ! Below-cloud scavenging by convective precipitation
      where ( fracc(js:je)>zmin )
        Frc(:) = max( 0., pfconv(js:je,jk-1)/fracc(js:je) )
        zbcscav(:) = zcollefc(:)*fracc(js:je)*0.24*ptmst*sqrt(Frc(:)*sqrt(Frc(:)))
        !zbcscav(:) = min( 1., zbcscav(:)/(1.+0.5*zbcscav(:)) ) !Time-centred
        zbcscav(:) = max( min( 1., zbcscav(:) ), 0. ) ! MJT suggestion
        xbcscav(:) = max( zbcscav(:)*zxtp10(:,jk,ktrac), 0. )
        pdep(:) = pdep(:) + xbcscav(:)*zclr0(:)
        zxtp10(:,jk,ktrac) = zxtp10(:,jk,ktrac) - xbcscav(:)
        !conwd(js:je,ktrac) = conwd(js:je,ktrac) + xbcscav(:)*zclr0(:)*zmtof(:)
      end where

      ! Below-cloud reevaporation of convective rain
      ! This never triggers for JLM convection because pcevap=0.
      ! lmask(:) = jk>kbase(:) .and. pfconv(:,jk-1)>zmin .and. zclr0(:)>zmin
      ! where ( lmask(:) )
      !   pcevap = pfconv(:,jk-1) - pfconv(:,jk)
      !   zevap = pcevap/pfconv(:,jk-1)
      ! elsewhere
      !   zevap(:)=0.
      ! end where
      ! where ( lmask(:) .and. zevap<1. )
      !   zevap = Evfac(ktrac)*zevap
      ! end where
      ! where ( lmask(:) )
      !   zevap = max( 0., min( 1., zevap ) )
      !   xevap = max( conwd(:,ktrac)*zevap/zmtof(:), 0. ) !xevap is the grid-box-mean m.r. change
      !   pdep = pdep - xevap*zclr0
      !   zxtp10(:,jk,ktrac) = zxtp10(:,jk,ktrac) + xevap
      !   conwd(:,ktrac) = conwd(:,ktrac) - xevap*zclr0*zmtof(:)
      !   conwd(:,ktrac) = max( 0., conwd(:,ktrac) )
      ! end where

      ZXTP1(:)=(1.-pclcover(js:je,jk)-pclcon(js:je,jk))*ZXTP10(:,JK,ktrac)+ &
                PCLCOVER(js:je,JK)*ZXTP1C(:,JK,ktrac)+                      &
                pclcon(js:je,jk)*zxtp1con(:,jk,ktrac)
      zxtp1(:)=max(zxtp1(:),0.)
      ZDXTE(:)=(ZXTP1(:)-XTM1(js:je,JK,ktrac))*PQTMST  !Total tendency (Dep + chem)
      !    CHANGE THE TOTAL TENDENCIES
      xte(js:je,jk,ktrac) = xte(js:je,jk,ktrac) + zdxte(:)
      wd(:,ktrac) = wd(:,ktrac) + pqtmst*pdep(:)*rhodz(js:je,jk)

    end do ! ktrac = 2,naero
  end do   ! jk = top,kl

  !$acc loop vector
  do jl = 1,imax
    iq = jl + js - 1
    so2wd(iq) = so2wd(iq) + wd(jl,ITRACSO2) 
    so4wd(iq) = so4wd(iq) + wd(jl,ITRACSO4)
    bcwd(iq) = bcwd(iq) + wd(jl,ITRACBC) + wd(jl,ITRACBC+1)
    ocwd(iq) = ocwd(iq) + wd(jl,ITRACOC) + wd(jl,ITRACOC+1)
  end do
  !$acc loop seq
  DO JT=ITRACDU,ITRACDU+NDUST-1
    !$acc loop vector  
    do jl = 1,imax
      iq = jl + js - 1
      dustwd(iq,jt-itracdu+1) = dustwd(iq,jt-itracdu+1) + wd(jl,jt)
    end do  
  end do
  !$acc loop seq
  do jt = ITRACSA,ITRACSA+NSALT-1
    !$acc loop vector  
    do jl = 1,imax
      iq = jl + js - 1
      saltwd(iq) = saltwd(iq) + wd(jl,jt)
    end do  
  end do
  
end do ! tile = 1,ntiles
!$acc end parallel loop

!!$acc update self(xte)
!if ( maxval(xtm1(1:ifull,:,:)+xte(1:ifull,:,:)*PTMST)>6.5e-5 ) then
!  write(6,*) "xtg out-of-range after xtwetdep"
!  write(6,*) "xtg maxval,maxloc ",maxval(xtm1(1:ifull,:,:)+xte(1:ifull,:,:)*PTMST), &
!                                  maxloc(xtm1(1:ifull,:,:)+xte(1:ifull,:,:)*PTMST)
!end if

!$acc parallel loop copyin(taudar,zdayfac) copyout(dmsoh,dmsn3,so2oh)          &
!$acc   private(zzoh,zzo3,zzno2) present(prhop1,ptp1,rhodz,xtm1,xte,zoxidant)
do tile = 1,ntiles
  js = (tile-1)*imax + 1
  je = tile*imax

  !$acc loop vector
  do jl = js,je
    dmsoh(jl)=0.
    dmsn3(jl)=0.
    so2oh(jl)=0.    
  end do  
  
  !$acc loop seq
  do jk = 1,kl
    !$acc loop vector  
    do jl = 1,imax
      iq = jl + js - 1  

      !  OXIDANT CONCENTRATIONS IN MOLECULE/CM**3
      ! -- levels are already inverted --
      ZZOH(jl,jk) = ZOXIDANT(iq,jk,1)
      ZZO3(jl,jk) = ZOXIDANT(iq,jk,3)*PRHOP1(iq,jk)*1.e-3
      ZZNO2(jl,jk) = ZOXIDANT(iq,jk,4)*PRHOP1(iq,jk)*1.e-3
  
    end do  
  end do
  
  
  !$acc loop seq
  do jk = 1,kl
    !$acc loop vector  
    do jl = 1,imax
      iq = jl + js - 1
      X=PRHOP1(iq,JK)      
      IF ( taudar(iq)>0.5 ) THEN
        !   DAY-TIME CHEMISTRY        
        ZXTP1SO2=XTM1(iq,JK,ITRACSO2)+XTE(iq,JK,ITRACSO2)*PTMST
        ZTK2=ZK2*(PTP1(iq,JK)/300.)**(-3.3)
        ZM=X*ZNAMAIR
        ZHIL=ZTK2*ZM/ZK2I
        ZEXP=LOG10(ZHIL)
        ZEXP=1./(1.+ZEXP*ZEXP)
        ZTK23B=ZTK2*ZM/(1.+ZHIL)*ZK2F**ZEXP
        ZSO2=ZXTP1SO2*ZZOH(JL,JK)*ZTK23B*ZDAYFAC(iq)
        ZSO2=MIN(ZSO2,ZXTP1SO2*PQTMST)
        ZSO2=MAX(ZSO2,0.)
        XTE(iq,JK,ITRACSO2)=XTE(iq,JK,ITRACSO2)-ZSO2
        XTE(iq,JK,ITRACSO4)=XTE(iq,JK,ITRACSO4)+ZSO2
        so2oh(iq)=so2oh(iq)+zso2*rhodz(iq,jk)

        ZXTP1DMS=XTM1(iq,JK,ITRACDMS)+XTE(iq,JK,ITRACDMS)*PTMST
        T=PTP1(iq,JK)
!       ZTK1=(T*EXP(-234./T)+8.46E-10*EXP(7230./T)+ &
!            2.68E-10*EXP(7810./T))/(1.04E+11*T+88.1*EXP(7460./T)) !Original
        ztk1=1.646e-10-1.850e-12*t+8.151e-15*t**2-1.253e-17*t**3 !Cubic fit good enough
        ztk1=max(ztk1,5.e-12) !Because cubic falls away for T > 300 K
        ! MJT notes - LDR employs 1.5 here, but CCAM's DMS burden is too high
        ! so we revert to the original factor of 2 used for rotstayn and lohmann 2002
        ztk1=2.*ztk1          !This is the fudge factor to account for other oxidants
        !ztk1=1.5*ztk1        !This is the fudge factor to account for other oxidants
        ZDMS=ZXTP1DMS*ZZOH(JL,JK)*ZTK1*ZDAYFAC(iq)
        ZDMS=MIN(ZDMS,ZXTP1DMS*PQTMST)
        XTE(iq,JK,ITRACDMS)=XTE(iq,JK,ITRACDMS)-ZDMS
        XTE(iq,JK,ITRACSO2)=XTE(iq,JK,ITRACSO2)+ZDMS
        dmsoh(iq)=dmsoh(iq)+zdms*rhodz(iq,jk)
      ELSE

        !   NIGHT-TIME CHEMISTRY
        ZXTP1DMS=XTM1(iq,JK,ITRACDMS)+XTE(iq,JK,ITRACDMS)*PTMST
        ZTK3=ZK3*EXP(520./PTP1(iq,JK))
        !    CALCULATE THE STEADY STATE CONCENTRATION OF NO3
        ZQT=1./PTP1(iq,JK)
        ZQT3=300.*ZQT
        ZRHOAIR=PRHOP1(iq,JK)*ZNAMAIR
        ZKNO2O3=1.2E-13*EXP(-2450.*ZQT)
        ZKN2O5AQ=0.1E-04
        ZRX1=2.2E-30*ZQT3**3.9*ZRHOAIR
        !ZRX2=1.5E-12*ZQT3**0.7
        ZRX12=1.467e-18*ZQT3**3.2*ZRHOAIR !=ZRX1/ZRX2
        ZKNO2NO3=ZRX1/(1.+ZRX12)*0.6**(1./(1.+(LOG10(ZRX12))**2))
        !ZEQN2O5=4.E-27*EXP(10930.*ZQT)
        !ZKN2O5=ZKNO2NO3/ZEQN2O5
        ZKN2O5=5.5E-4*ZQT3**3.9*ZRHOAIR*EXP(-10930.*ZQT)/(1.+ZRX12)*0.6**(1./(1.+(LOG10(ZRX12))**2))

        ZNO3=ZKNO2O3*(ZKN2O5+ZKN2O5AQ)*ZZNO2(JL,JK)*ZZO3(JL,JK)
        ZZQ=ZKNO2NO3*ZKN2O5AQ*ZZNO2(JL,JK)+(ZKN2O5+ZKN2O5AQ)*ZTK3*ZXTP1DMS*X*6.022E+20/ZMOLGS
        ZNO3=0.
        IF(ZZQ>0.) THEN
          ZNO3=ZNO3/ZZQ
        ENDIF
        ZDMS=ZXTP1DMS*ZNO3*ZTK3
        ZDMS=MIN(ZDMS,ZXTP1DMS*PQTMST)
        XTE(iq,JK,ITRACDMS)=XTE(iq,JK,ITRACDMS)-ZDMS
        XTE(iq,JK,ITRACSO2)=XTE(iq,JK,ITRACSO2)+ZDMS
        dmsn3(iq)=dmsn3(iq)+zdms*rhodz(iq,jk)
      ENDIF
    end do
  end do
  
end do ! tile = 1,ntiles
!$acc end parallel loop

!$acc update self(xte)
!$acc end data

if ( maxval(xtm1(1:ifull,:,:)+xte(1:ifull,:,:)*PTMST)>6.5e-5 ) then
  write(6,*) "xtg is out-of-range at end of xtchemie"
  write(6,*) "xtg maxval,maxloc ",maxval(xtm1(1:ifull,:,:)+xte(1:ifull,:,:)*PTMST), &
                                  maxloc(xtm1(1:ifull,:,:)+xte(1:ifull,:,:)*PTMST)
end if

RETURN
END subroutine xtchemie


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Dust settling

subroutine dsettling(tdt,rhoa,tmp,delz,prf,xtg)

!     Inputs:
real, intent(in) :: tdt                  !timestep (s)
real, dimension(:,:), intent(in) :: rhoa !air density (kg/m3)
real, dimension(:,:), intent(in) :: tmp  !temperature (K)
real, dimension(:,:), intent(in) :: delz !Layer thickness (m)
real, dimension(:,:), intent(in) :: prf  !Pressure (hPa)
real, dimension(:,:,:), intent(inout) :: xtg

! Local work arrays and variables
real, dimension(size(rhoa,1)) :: dfall
real c_stokes, corr, c_cun, vd_cor, b, newxtg
integer nt, k, iq, imax, kl

! Start code : ----------------------------------------------------------

imax = size(rhoa,1)
kl   = size(rhoa,2)

! Settling velocity (m/s) for each soil classes (Stokes Law)
! DUSTDEN     soil class density             (kg/m3)
! DUSTREFF    effective radius of soil class (m)
! grav        gravity                        (m/s2)
! 0.5         upper limit with temp correction (already incorporated with dzmin_gbl - MJT)

do nt = 1, NDUST

  ! Solve at the model top    
  do iq = 1,imax
    ! Dynamic viscosity
    C_Stokes = 1.458E-6*TMP(iq,kl)**1.5/(TMP(iq,kl)+110.4) 
    ! Cuningham correction
    Corr = 6.6E-8*prf(iq,kl)/1013.*TMP(iq,kl)/293.15
    C_Cun = 1. + 1.249*corr/dustreff(nt)
    ! Settling velocity
    Vd_cor = 2./9.*grav*dustden(nt)*dustreff(nt)**2/C_Stokes*C_Cun
  
    ! Update mixing ratio
    b = tdt*VD_cor/DELZ(iq,kl)
    newxtg = xtg(iq,kl,nt+itracdu-1)*exp(-b)
    newxtg = max( newxtg, 0. )
    dfall(iq) = max( xtg(iq,kl,nt+itracdu-1) - newxtg, 0. )
    xtg(iq,kl,nt+itracdu-1) = newxtg
  end do

  ! Solve each vertical layer successively (layer k)
  do k = kl-1,1,-1
    do iq = 1,imax
      ! Dynamic viscosity
      C_Stokes = 1.458E-6*TMP(iq,k)**1.5/(TMP(iq,k)+110.4) 
      ! Cuningham correction
      Corr = 6.6E-8*prf(iq,k)/1013.*TMP(iq,k)/293.15
      C_Cun = 1. + 1.249*corr/dustreff(nt)
      ! Settling velocity
      Vd_cor = 2./9.*grav*dustden(nt)*dustreff(nt)**2/C_Stokes*C_Cun
      
      ! Update mixing ratio
      b = tdt*Vd_cor/DELZ(iq,k)
      dfall(iq) = dfall(iq)*delz(iq,k+1)*rhoa(iq,k+1)/(delz(iq,k)*rhoa(iq,k))
      ! Fout = 1.-exp(-b)
      ! Fthru = 1.-Fout/b
      newxtg = xtg(iq,k,nt+itracdu-1)*exp(-b) + dfall(iq)*(1.-exp(-b))/b
      newxtg = max( newxtg, 0. )
      dfall(iq) = max( xtg(iq,k,nt+itracdu-1) + dfall(iq) - newxtg, 0. )
      xtg(iq,k,nt+itracdu-1) = newxtg
    end do
  end do
  
end do

return
end subroutine dsettling

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Dust emissions

subroutine dustem(tdt,rhoa,wg,w10m,dz1,vt,snowd,erod,duste,xtg)

real, intent(in) :: tdt                      !Leapfrog timestep (s) (substep and long step)
real, dimension(:), intent(in) :: rhoa       !air density (kg/m3)
real, dimension(:), intent(in) :: wg         !ground wetness (fraction of field capacity)
real, dimension(:), intent(in) :: w10m       !10m windspeed (m/s)
real, dimension(:), intent(in) :: dz1        !Lowest layer thickness (m)
real, dimension(:), intent(in) :: vt         !Transfer velocity at surface for dry deposition (m/s)
real, dimension(:), intent(in) :: snowd      !Snow depth (mm equivalent water)
real, dimension(:,:), intent(in) :: erod
real, dimension(:,:), intent(inout) :: duste
real, dimension(:,:,:), intent(inout) :: xtg
real, dimension(size(xtg,1)) :: snowa         !Estimated snow areal coverage
real, dimension(size(xtg,1)) :: airmas, airden
real g, den, diam
real u_ts0, u_ts, srce, dsrc, a, b, veff
integer n, m, iq, imax, kl

real, parameter :: w_dust = 15. ! maxmium wind speed for dust emissions

! This frac_s array gives fraction of source in each size bin.
! Does not quite add to 1, because larger sizes (omitted) account for some too.
! All source is in first four bins, even when using eight bins, since next four are hydrophilic.
real, dimension(ndust), parameter :: frac_s = (/ 0.1, 0.25, 0.25, 0.25 /)
integer, dimension(ndust), parameter :: ipoint = (/ 3, 2, 2, 2 /)                  ! Pointer used for dust classes
                                                                                   ! (sand=1, silt=2, clay=3)

! Start code : ----------------------------------------------------------

imax = size(xtg,1)
kl   = size(xtg,2)

g = grav*1.e2
airden = rhoa*1.e-3

! Convert snow depth to estimated areal coverage (see Zender et al 2003, JGR 108, D14, 4416)
! 0.05 m is the geometrical snow thickness for 100% areal coverage.
snowa = min( 1., snowd/5. )
airmas = dz1 * rhoa  ! kg/m2

do n = 1, ndust
  do iq = 1,imax
  
    ! Threshold velocity as a function of the dust density and the diameter
    den = dustden(n)*1.e-3
    diam = 2.*dustreff(n)*1.e2
    ! Pointer to the 3 classes considered in the source data files
    m = ipoint(n)
      
    ! Following is from Ginoux et al (2004) Env. Modelling & Software.
    u_ts0 = 0.13*1.e-2*sqrt(den*g*diam/airden(iq))*sqrt(1.+0.006/den/g/diam**2.5)/ &
            sqrt(1.928*(1331.*diam**1.56+0.38)**0.092-1.)
  
    ! Case of surface dry enough to erode
    if ( wg(iq)<0.1 ) then
      !Tuning suggested for Asian source by P. Ginoux
      u_ts = max( 0., u_ts0*(1.2+0.2*log10(max( 1.e-3, wg(iq) ))) )
    else
      ! Case of wet surface, no erosion
      u_ts = 100.
    end if
  
    ! MJT notes - erod should be zero for ocean points
    
    !srce = frac_s(n)*erod(iq,m)*dxy(iq) ! (m2)
    srce = frac_s(n)*erod(iq,m) ! (fraction) - MJT suggestion
    if ( w10m(iq) < w_dust ) then
      dsrc = (1.-snowa(iq))*Ch_dust*srce*W10m(iq)**2*(W10m(iq)-u_ts) ! (kg/s/m2)
    else
      ! limit maximum wind speed to w_dust m/s for emissions - MJT sugestion  
      dsrc = (1.-snowa(iq))*Ch_dust*srce*w_dust**2*(w_dust-u_ts) ! (kg/s/m2)  
    end if
    dsrc = max( 0., dsrc )

    ! Calculate dust mixing ratio tendency at first model level.
    a = dsrc / max(airmas(iq),0.1)
    duste(iq,n) = duste(iq,n) + dsrc ! Diagnostic
      
    ! Calculate turbulent dry deposition at surface
    ! Use full layer thickness for CSIRO model (should be correct if Vt is relative to mid-layer)
    veff = max( Vt(iq)*(wg(iq)+(1.-wg(iq))*exp(-max( 0., w10m(iq)-u_ts0 ))), 0. )
    b = Veff / dz1(iq)

    ! Update mixing ratio
    ! Write in form dx/dt = a - b*x (a = source term, b = drydep term)
    ! solution is x = a/b + (X0-a/b)*exp(-b*tdt).  However, in split form
    ! x = X0 + a*tdt, and x = X0*exp(-b*tdt), or combined
    ! x = (X0 + a*tdt)*exp(-b*tdt)
    xtg(iq,1,n+itracdu-1) = (xtg(iq,1,n+itracdu-1)+a*tdt)*exp(-min(b*tdt,50.))
    xtg(iq,1,n+itracdu-1) = max( 0., xtg(iq,1,n+itracdu-1) )
  end do
end do

return
end subroutine dustem

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! dustage

!subroutine dustage(tdt,rhg) !Inputs
!
!!     Inputs:
!real tdt              !Leapfrog timestep
!real rhg(imax,kl)    !RH (fraction 0 to 1)
!
!! Local work arrays and variables
!real rrate(imax,kl)
!real vvso2,rk,xd,dx
!integer k,mg,nt
!
!! Start code : ----------------------------------------------------------
!
!! Reference for this scheme is
!! Fan et al. (2004): Impact of air pollution on wet deposition of mineral dust aerosols, GRL 31,
!! L02104, doi:10.1029/2003GL018501.
!
!! Loop over 4 dust size bins
!! Reaction rate is a linear function of SO2 concentration
!
!do k=1,kl
!  do mg=1,imax
!    vvso2 = xtg(mg,k,itracso2)*29./64. !Vol. mixing ratio of SO2
!    rk = 0.01 * max (0.1, dim(rhg(mg,k),0.5))
!    xd = 0.
!    do nt = itracdu, itracdu+ndust-1
!      xd = xd + xtg(mg,k,nt)
!    enddo
!    rrate(mg,k) = rk * vvso2 / max (1.e-12, xd)
!  enddo
!enddo
!
!! Reduce hydrophobic dust and increase hydrophilic dust
!do nt = itracdu, itracdu+ndust-1
!  do k=1,kl
!    do mg=1,imax
!      dx = min (xtg(mg,k,nt), rrate(mg,k)*tdt*xtg(mg,k,nt))
!      xtg(mg,k,nt) = xtg(mg,k,nt) - dx
!      xtg(mg,k,nt+ndust) = xtg(mg,k,nt+ndust) + dx
!    enddo
!  enddo
!enddo
!
!return
!end subroutine dustage

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! A simple diagnostic treatment of seasalt aerosol (LDR 3/02)
!
!subroutine seasalt(land,fracice,zmid,pblh,v10m,ssn,imax,kl) !Inputs
!
!integer, intent(in) :: imax, kl
!logical, dimension(imax), intent(in) :: land   !True for land points
!real, dimension(imax), intent(in) :: fracice   !Sea-ice fraction
!real, dimension(imax,kl), intent(in) :: zmid   !Height of full level (m)
!real, dimension(imax), intent(in) :: pblh      !PBL height (m)
!real, dimension(imax), intent(in) :: v10m      !10m windpseed, including effect of sub-grid gustiness (m/s)
!real, dimension(imax,2) :: ssn_base
!real, dimension(imax,kl,2), intent(inout) :: ssn
!integer k
!
!! Calculate number and mass concentration of seasalt within the marine BL.
!! Set seasalt conc. to zero elsewhere.
!! Height of BL taken from ncarpbl scheme, so needs this turned on (although we 
!! set it to 2000m in hvertmx if ncarpbl=F)
!! The first mode is the "film-drop" mode, and the second is the "jet-drop" mode.
!! Smaller number mode radii are given by Nilsson et al. (2001) cf. O'Dowd's.
!!
!! References: O'Dowd et al. (1997) Atmos. Environ. 31, 73-80
!!             Jones et al. (2001) JGR 106, 20293-20310.
!!             Nilsson et al. (2001) JGR 106, 32139-32154.
!
!! Jones et al. give different windspeed relations for v10m < 2, 2 <= v10m <= 17.5,
!! and v10m > 17.5.
!where ( v10m<2. )
!  ssn_base(:,1)=3.856e6*(1.-exp(-0.736*v10m))
!  ssn_base(:,2)=0.671e6*(1.-exp(-1.351*v10m))
!elsewhere ( v10m<17.5 )
!  ssn_base(:,1)=10.**(0.0950*v10m+6.2830) 
!  ssn_base(:,2)=10.**(0.0422*v10m+5.7122) 
!elsewhere
!  ssn_base(:,1)=1.5e8*(1.-97.874*exp(-0.313*v10m))
!  ssn_base(:,2)=3.6e6*(1.-103.926*exp(-0.353*v10m))
!end where
!  
!! Number-to-mass conversion factors are derived from the parameters of the two
!! lognormal modes given by Nilsson, together with rhosalt=2.0e3 kg/m3.
!do k=1,kl
!  where ( .not.land .and. zmid(:,k)<pblh )
!    ssn(:,k,1) = ssn_base(:,1)
!    ssn(:,k,2) = ssn_base(:,2)
!  elsewhere ( .not.land )
!    ssn(:,k,1) = 1.e7*exp(-zmid(:,k)/3000.)
!    ssn(:,k,2) = 1.e6*exp(-zmid(:,k)/1500.)
!  elsewhere
!    ssn(:,k,1) = 0.
!    ssn(:,k,2) = 0.
!  end where
!  ! Reduce over sea ice...
!  ssn(:,k,1) = (1.-fracice)*ssn(:,k,1)
!  ssn(:,k,2) = (1.-fracice)*ssn(:,k,2)
!end do
!
!! These relations give ssm in kg/m3 based on ssn in m^{-3}...
!!ssm(:,:,1)=ssn(:,:,1)/saltsmallmtn
!!ssm(:,:,2)=ssn(:,:,2)/saltlargemtn
!
!return
!end subroutine seasalt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! sea salt settling
subroutine ssettling(tdt,rhoa,tmp,delz,prf,xtg)

real, intent(in) :: tdt                  !timestep (s)
real, dimension(:,:), intent(in) :: rhoa !air density (kg/m3)
real, dimension(:,:), intent(in) :: tmp  !temperature (K)
real, dimension(:,:), intent(in) :: delz !Layer thickness (m)
real, dimension(:,:), intent(in) :: prf  !Pressure (hPa)
real, dimension(:,:,:), intent(inout) :: xtg
real, dimension(size(rhoa,1)) :: dfall
real vd_cor, newxtg, b, c_stokes, corr, c_cun
integer nt, k, iq, imax, kl

imax = size(rhoa,1)
kl   = size(rhoa,2)

do nt = 1,nsalt
    
  ! Solve at the model top  
  do iq = 1,imax
    ! Dynamic viscosity
    C_Stokes = 1.458E-6*TMP(iq,kl)**1.5/(TMP(iq,kl)+110.4) 
    ! Cuningham correction
    Corr = 6.6E-8*prf(iq,kl)/1013.*TMP(iq,kl)/293.15
    C_Cun = 1. + 1.249*corr/saltreff(nt)
    ! Settling velocity
    Vd_cor = 2./9.*grav*saltden(nt)*saltreff(nt)**2/C_Stokes*C_Cun
  
    ! Update mixing ratio
    b = tdt*VD_cor/DELZ(iq,kl)
    newxtg = xtg(iq,kl,nt+itracsa-1)*exp(-b)
    newxtg = max( newxtg, 0. )
    dfall(iq) = max( xtg(iq,kl,nt+itracsa-1) - newxtg, 0. )
    xtg(iq,kl,nt+itracsa-1) = newxtg
  end do
  
  ! Solve each vertical layer successively (layer k)
  do k = kl-1,1,-1
    do iq = 1,imax
      ! Dynamic viscosity
      C_Stokes = 1.458E-6*TMP(iq,k)**1.5/(TMP(iq,k)+110.4) 
      ! Cuningham correction
      Corr = 6.6E-8*prf(iq,k)/1013.*TMP(iq,k)/293.15
      C_Cun = 1. + 1.249*corr/saltreff(nt)
      ! Settling velocity
      Vd_cor = 2./9.*grav*saltden(nt)*saltreff(nt)**2/C_Stokes*C_Cun
      
      ! Update mixing ratio
      b = tdt*Vd_cor/DELZ(iq,k)
      dfall(iq) = dfall(iq)*delz(iq,k+1)*rhoa(iq,k+1)/(delz(iq,k)*rhoa(iq,k))
      ! Fout = 1.-exp(-b)
      ! Fthru = 1.-Fout/b
      newxtg = xtg(iq,k,nt+itracsa-1)*exp(-b) + dfall(iq)*(1.-exp(-b))/b
      newxtg = max( newxtg, 0. )
      dfall(iq) = max( xtg(iq,k,nt+itracsa-1) + dfall(iq) - newxtg, 0. )
      xtg(iq,k,nt+itracsa-1) = newxtg
    end do
  end do
  
end do

return
end subroutine ssettling

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! sea salt emissions
subroutine seasaltem(tdt,v10m,vt,rhoa,dz1,salte,xtg,locean)

integer n, iq, imax, kl
real, intent(in) :: tdt
real df0dd, dfdd, ftw, fsw, a, b, veff, diam
real, dimension(:), intent(in) :: v10m    ! 10m wind speed
real, dimension(:), intent(in) :: vt      ! transfer velocity
real, dimension(:), intent(in) :: dz1     ! layer thickness
real, dimension(:), intent(in) :: rhoa    ! air density
real, dimension(:), intent(inout) :: salte
real, dimension(:,:,:), intent(inout) :: xtg
real, dimension(size(xtg,1)) :: wu10
real, dimension(nsalt) :: mtnfactor, saltrange
logical, dimension(:), intent(in) :: locean

imax = size(xtg,1)
kl   = size(xtg,2)

! Follows Sofiev et al 2011 for emissions
wu10(:) = 3.84e-6*v10m(:)**3.41

mtnfactor(1) = saltsmallmtn
mtnfactor(2) = saltlargemtn
saltrange = (/ 0.4e-6, 3.5e-6 /) ! 0.1-0.5um and 0.5-4um

do n = 1,nsalt
  do iq = 1,imax
    ! Calculate emissions  
    diam = 2.*saltreff(n)*1.e6
    
    !df0dr = 3.6e5*(1.+0.057*diam**1.05)/(diam**3)*10**(1.19*exp(-((0.38-log(diam))/0.65)**2))
  
    ! Follows Sofiev et al 2011 for emissions (accounts for <0.1 microns)
    df0dd = 1.e6*exp(-0.09/(diam+3.e-3))/(2.+exp(-5./diam))*(1.+0.05*diam**1.05)/(diam**3) &
           *10**(1.05*exp(-((0.27-log(diam))/1.1)**2))
    
    ftw = 0.48*diam**(-0.36)   ! 15C
    !ftw = 0.15*diam**(-0.88)  ! 5C
    !ftw = 0.092*diam**(-0.96) ! -2C
  
    fsw = 0.12*diam**(-0.71)
    !fsw = 5.85*e-5*diam**-1.7 ! sal = 0.

    dfdd = wu10(iq)*df0dd*ftw*fsw  ! number/micro/m^2/s
  
    if ( locean(iq) ) then
      a = dfdd*2.*saltrange(n)*1.e6/(rhoa(iq)*dz1(iq)) ! number/kg/s
    else
      a = 0.
    end if 
  
    a = a/mtnfactor(n) ! kg/kg/s  
  
    salte(iq) = salte(iq) + a*rhoa(iq)*dz1(iq) ! Diagnostic
  
    ! Calculate turbulent dry deposition at surface
    veff = Vt(iq)
    b = Veff / dz1(iq)

    ! Update mixing ratio
    ! Write in form dx/dt = a - b*x (a = source term, b = drydep term) 
    ! solution is x = a/b + (X0-a/b)*exp(-b*tdt).  However, in split form
    ! x = X0 + a*tdt, and x = X0*exp(-b*tdt), or combined
    ! x = (X0 + a*tdt)*exp(-b*tdt)
    xtg(iq,1,n+itracsa-1) = (xtg(iq,1,n+itracsa-1)+a*tdt)*exp(-min(b*tdt,50.))
    xtg(iq,1,n+itracsa-1) = max( 0., xtg(iq,1,n+itracsa-1) )
  end do
end do    

return
end subroutine seasaltem

end module aerosolldr
