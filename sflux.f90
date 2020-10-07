! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2020 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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

! CCAM interface for surface flux routines. Includes standard land-surface scheme,
! prescribed SSTs and sea-ice, CABLE interface, urban interface and MLO interface.
      
! nsib=3              Standard land-surface scheme with SIB and Gratz data
! nsib=5              Standard land-surface scheme with MODIS data
! nsib=6              CABLE land-surface scheme with CABLE diagnostics
! nsib=7              CABLE land-surface scheme with CCAM diagnostics
! nmlo=0              Prescriped SSTs and sea-ice with JLM skin enhancement
! nmlo>0 and mlo<=9   KPP ocean mixing
! nmlo>9              Use external PCOM ocean model
! nurban>0            Use urban scheme
! nriver=1            Use river routing (automatically enabled for abs(nmlo)>1)
    
module sflux_m

implicit none

private
public sflux, sflux_init

real, save :: ri_max, zologbgin, ztv, z1onzt, chnsea
real, save :: srcp
real, dimension(:), allocatable, save :: taftfh, taftfhg

integer, parameter :: ntest=0   ! ntest= 0 for diags off; ntest= 1 for diags on
real, parameter :: bprm=5.,cms=5.,chs=2.6,vkar=.4
real, parameter :: fmroot=.57735     ! was .4 till 7 Feb 1996

real, dimension(:), allocatable :: vmag,azmin,uav,vav,rho
real, dimension(:), allocatable :: oldrunoff,oldsnowmelt,oldevspsbl,oldsbl
real, dimension(:), allocatable :: af,aft,ri
real, dimension(:), allocatable :: fg_ocn, fg_ice, eg_ocn, eg_ice
real, dimension(:), allocatable :: taux_ocn, taux_ice, tauy_ocn, tauy_ice
real, dimension(:), allocatable :: river_inflow
real, dimension(:), allocatable :: tss_save

#ifdef csircoupled
real, dimension(:), allocatable :: dumrg, dumx
real, dimension(:), allocatable :: dumtaux_ocn, dumtauy_ocn, dumtaux_ice, dumtauy_ice
real, dimension(:), allocatable :: zonx, zony, zonz, costh, sinth
#endif

contains

subroutine sflux_init

use const_phys                     ! Physical constants
use newmpar_m, only : ifull        ! Grid parameters
use parm_m                         ! Model configuration
use permsurf_m                     ! Fixed surface arrays
use sigs_m                         ! Atmosphere sigma levels
use soil_m                         ! Soil and surface data

implicit none

if ( nsib==3 .or. nsib==5 ) then
  allocate( taftfh(ifull), taftfhg(ifull) )
  taftfh(:)=.05        ! just a diag default for sea points
  taftfhg(:)=7.e-4     ! just a diag default for sea points
end if

ri_max=(1./fmroot -1.)/bprm    ! i.e. .14641
zologbgin=log(zmin/zobgin)     ! pre-calculated for all except snow points
ztv=exp(vkar/sqrt(chn10))/10.  ! proper inverse of ztsea
z1onzt=300.*rdry*(1.-sig(1))*ztv/grav
chnsea=(vkar/log(z1onzt))**2   ! should give .00085 for csiro9
srcp = sig(1)**(rdry/cp)

allocate( vmag(ifull), azmin(ifull), uav(ifull), vav(ifull), rho(ifull) )
allocate( oldrunoff(ifull), oldsnowmelt(ifull), oldevspsbl(ifull), oldsbl(ifull) )
allocate( af(ifull), aft(ifull), ri(ifull) )
allocate( fg_ocn(ifull), fg_ice(ifull), eg_ocn(ifull), eg_ice(ifull) )
allocate( taux_ocn(ifull), taux_ice(ifull), tauy_ocn(ifull), tauy_ice(ifull) )
allocate( river_inflow(ifull) )
allocate( tss_save(ifull) )
#ifdef csircoupled
allocate( dumrg(ifull), dumx(ifull) )
allocate( dumtaux_ocn(ifull), dumtauy_ocn(ifull) )
allocate( dumtaux_ice(ifull), dumtauy_ice(ifull) )
allocate( zonx(ifull), zony(ifull), zonz(ifull) )
allocate( costh(ifull), sinth(ifull) )
#endif

return
end subroutine sflux_init


subroutine sflux
      
use arrays_m                       ! Atmosphere dyamics prognostic arrays
use cable_ccam, only : sib4        ! CABLE interface
use cc_mpi                         ! CC MPI routines
use cc_omp, only : imax, ntiles    ! CC OpenMP routines
use const_phys                     ! Physical constants
use diag_m                         ! Diagnostic routines
use estab                          ! Liquid saturation function
use extraout_m                     ! Additional diagnostics
use morepbl_m                      ! Additional boundary layer diagnostics
use newmpar_m                      ! Grid parameters
use nharrs_m                       ! Non-hydrostatic atmosphere arrays
use nsibd_m                        ! Land-surface arrays
use parm_m                         ! Model configuration
use pbl_m                          ! Boundary layer arrays
use prec_m                         ! Precipitation
use raddiag_m                      ! Radiation diagnostic
use riverarrays_m                  ! River data
use savuvt_m                       ! Saved dynamic arrays
use screen_m                       ! Screen level diagnostics
use sigs_m                         ! Atmosphere sigma levels
use soil_m                         ! Soil and surface data
use soilsnow_m                     ! Soil, snow and surface data
use work2_m                        ! Diagnostic arrays
use work3_m                        ! Mk3 land-surface diagnostic arrays
      
#ifdef csircoupled
use parmgeom_m
use vcom_ccam
use vecsuv_m
use xyzinfo_m
#endif

implicit none
    
integer is, ie, tile, iq, k
real, dimension(imax) :: rg_error, qsat_tmp

!     stability dependent drag coefficients using Louis (1979,blm) f'
!     n.b. cduv, cdtq are returned as drag coeffs mult by vmod
!          (cduv=cduv*vmod; cdtq=cdtq*vmod)

!     t, u, v, qg are current values
!     tss is surface temperature
!     dw is soil wetness availability (1. for ocean) - not needed here
!     fg is sensible heat flux (was h0)
!     eg is latent heat flux (was wv)
!     dfgdt is dfgdt (was csen in surfupa/b)
!     degdt is degdt (was ceva in surfupa/b)
  
! allocated arrays
!$omp do schedule(static) private(is,ie)
do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax  
  oldrunoff(is:ie) = runoff(is:ie)
  oldsnowmelt(is:ie) = snowmelt(is:ie)
  oldevspsbl(is:ie) = evspsbl(is:ie)
  oldsbl(is:ie) = sbl(is:ie)
  fg_ocn(is:ie)=0.      ! dummy value
  fg_ice(is:ie)=0.      ! dummy value
  eg_ocn(is:ie)=0.      ! dummy value
  eg_ice(is:ie)=0.      ! dummy value
  taux_ocn(is:ie)=0.    ! dummy value
  taux_ice(is:ie)=0.    ! dummy value
  tauy_ocn(is:ie)=0.    ! dummy value
  tauy_ice(is:ie)=0.    ! dummy value
!     using av_vmod (1. for no time averaging)
  azmin(is:ie) = bet(1)*t(is:ie,1)/grav
  rho(is:ie) = ps(is:ie)/(rdry*tss(is:ie))
  uav(is:ie) = av_vmod*u(is:ie,1) + (1.-av_vmod)*savu(is:ie,1)   
  vav(is:ie) = av_vmod*v(is:ie,1) + (1.-av_vmod)*savv(is:ie,1)  
  vmag(is:ie) = max( sqrt(uav(is:ie)**2+vav(is:ie)**2), vmodmin )    ! vmag used to calculate ri
  tss_save(is:ie) = tss(is:ie)

  taux(is:ie) = 0.      ! dummy value
  tauy(is:ie) = 0.      ! dummy value  
  zo(is:ie) = 0.        ! dummy value
  zoh(is:ie) = 0.       ! dummy value
  zoq(is:ie) = 0.       ! dummy value
  cduv(is:ie) = 0.      ! dummy value
  cdtq(is:ie) = 0.      ! dummy value
  ga(is:ie) = 0.        ! for ocean points in ga_ave diagnostic
  theta(is:ie) = t(is:ie,1)/srcp
  if ( ntsur==7 ) then
    vmod(is:ie) = sqrt(uav(is:ie)**2+vav(is:ie)**2)  
  else
    vmod(is:ie) = vmag(is:ie)  
  end if
end do
!$omp end do nowait

if (diag.or.ntest==1) then
!$omp barrier
!$omp single    
  if (mydiag) then
    if (land(idjd)) then
      write(6,*) 'entering sflux ktau,nsib,ivegt,isoilm,land ',ktau,nsib,ivegt(idjd),isoilm(idjd),land(idjd)
      write(6,*) 'idjd,id,jd,slwa,sgsave ',idjd,id,jd,slwa(idjd),sgsave(idjd)
      write(6,*) 'snowd,sicedep,condx ',snowd(idjd),sicedep(idjd),condx(idjd)
      write(6,*) 't1,tss ',t(idjd,1),tss(idjd)
      write(6,*) 'wb ',(wb(idjd,k),k=1,ms)
      write(6,*) 'tgg ',(tgg(idjd,k),k=1,ms)
    endif
  end if
  call maxmin(t,' t',ktau,1.,kl)
!$omp end single
endif

!--------------------------------------------------------------
!$omp master
call START_LOG(sfluxwater_begin)
!$omp end master
if ( nmlo==0 ) then ! prescribed SSTs
!$omp barrier
!$omp single
  call sflux_sea                                                                                 ! sea
!$omp end single
elseif (abs(nmlo)>=1.and.abs(nmlo)<=9) then ! prognostic SSTs                                    ! sea
  call sflux_mlo                                                                                 ! sea
end if                                                                                           ! sea
!$omp do schedule(static) private(is,ie)
do tile = 1,ntiles                                                                               ! sea
  is = (tile-1)*imax + 1                                                                         ! sea
  ie = tile*imax                                                                                 ! sea
  call nantest("after sflux_water",is,ie)                                                        ! sea
end do                                                                                           ! sea
!$omp end do nowait
!$omp master
call END_LOG(sfluxwater_end)
!$omp end master


!--------------------------------------------------------------      
!$omp master
call START_LOG(sfluxland_begin)
!$omp end master
select case(nsib)                                                                                ! land
  case(3,5)                                                                                      ! land
!$omp barrier
!$omp single
    call sflux_land                                                                              ! land
!$omp end single
  case(7)                                                                                        ! land
    call sib4 ! call cable                                                                       ! land
end select                                                                                       ! land
! update remaining diagnostic arrays                                                             ! land
!$omp do schedule(static) private(is,ie,iq)
do tile = 1,ntiles                                                                               ! land
  is = (tile-1)*imax + 1                                                                         ! land
  ie = tile*imax                                                                                 ! land
  qsat_tmp(1:imax) = qsat(ps(is:ie),tss(is:ie),imax)                                             ! land
  do iq = is,ie                                                                                  ! land
    if ( land(iq) ) then                                                                         ! land
      qsttg(iq) = qsat_tmp(iq-is+1)                                                              ! land
      taux(iq) = rho(iq)*cduv(iq)*u(iq,1)                                                        ! land
      tauy(iq) = rho(iq)*cduv(iq)*v(iq,1)                                                        ! land
      sno(iq) = sno(iq) + conds(iq)                                                              ! land
      grpl(iq) = grpl(iq) + condg(iq)                                                            ! land
    end if                                                                                       ! land
  end do                                                                                         ! land
  call nantest("after sflux_land",is,ie)                                                         ! land
end do                                                                                           ! land
!$omp end do nowait
!$omp master
call END_LOG(sfluxland_end)
!$omp end master
!----------------------------------------------------------


!$omp master
call START_LOG(sfluxurban_begin)
!$omp end master
if (nurban/=0) then                                                                              ! urban
  call sflux_urban                                                                               ! urban
end if                                                                                           ! urban
!$omp do schedule(static) private(is,ie)
do tile = 1,ntiles                                                                               ! urban
  is = (tile-1)*imax + 1                                                                         ! urban
  ie = tile*imax                                                                                 ! urban
  call nantest("after sflux_urban",is,ie)                                                        ! urban
end do                                                                                           ! urban
!$omp end do nowait
!$omp master
call END_LOG(sfluxurban_end)
!$omp end master
! ----------------------------------------------------------------------


#ifdef csircoupled
!$omp barrier
!$omp single
if ( nmlo==0 ) then                                                                              ! VCOM
  call START_LOG(sfluxwater_begin)                                                               ! VCOM
                                                                                                 ! VCOM
  dumrg(:)=-rgsave(:)                                                                            ! VCOM
  dumx(:)=condx(:)/dt ! total precip                                                             ! VCOM
                                                                                                 ! VCOM
  zonx=real(                             -sin(rlat0*pi/180.)*y(1:ifull))                         ! VCOM  
  zony=real(sin(rlat0*pi/180.)*x(1:ifull)+cos(rlat0*pi/180.)*z(1:ifull))                         ! VCOM
  zonz=real(-cos(rlat0*pi/180.)*y(1:ifull)                       )                               ! VCOM
  costh= (zonx*ax(1:ifull)+zony*ay(1:ifull)+zonz*az(1:ifull)) &                                  ! VCOM
        /sqrt( max(zonx**2+zony**2+zonz**2,1.e-7) )                                              ! VCOM
  sinth=-(zonx*bx(1:ifull)+zony*by(1:ifull)+zonz*bz(1:ifull)) &                                  ! VCOM
        /sqrt( max(zonx**2+zony**2+zonz**2,1.e-7) )                                              ! VCOM 
                                                                                                 ! VCOM             
  dumtaux_ocn(:) = costh*taux_ocn - sinth*tauy_ocn                                               ! VCOM
  dumtauy_ocn(:) = sinth*taux_ocn + costh*tauy_ocn                                               ! VCOM
  dumtaux_ice(:) = costh*taux_ice - sinth*tauy_ice                                               ! VCOM
  dumtauy_ice(:) = sinth*taux_ice + costh*tauy_ice                                               ! VCOM
                                                                                                 ! VCOM
  call vcom_ccam_physics(sgdn,dumrg,condx,river_inflow,     &                                    ! VCOM
                 dumtaux_ocn,dumtauy_ocn,fg_ocn,eg_ocn,     &                                    ! VCOM
                 dumtaux_ice,dumtauy_ice,fg_ice,eg_ice,     &                                    ! VCOM
                 tgg(:,1),tggsn(:,1),fracice,sicedep)                                            ! VCOM
  tss(:) = (1.-fracice(:))*tgg(:,1) + fracice(:)*tggsn(:,1)                                      ! VCOM
                                                                                                 ! VCOM  
  !call nantest("after VCOM",1,ifull)                                                            ! VCOM
  call END_LOG(sfluxwater_end)                                                                   ! VCOM
end if                                                                                           ! VCOM
!$omp end single
#endif
! ----------------------------------------------------------------------


!$omp do schedule(static) private(is,ie)
do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax  

  ! scrnout is the standard CCAM screen level diagnostics.
  ! autoscrn contains the newer diagnostic calculation
  if (nmlo==0.and.(nsib==3.or.nsib==5).and.rescrn==0) then
    call scrnout(zo(is:ie),ustar(is:ie),zoh(is:ie),wetfac(is:ie),qsttg(is:ie),      &
                 qgscrn(is:ie),tscrn(is:ie),uscrn(is:ie),u10(is:ie),rhscrn(is:ie),  &
                 af(is:ie),aft(is:ie),ri(is:ie),vmod(is:ie),bprm,cms,chs,chnsea,    &
                 is,ie)
  else
    call autoscrn(is,ie)
  end if

  ! ----------------------------------------------------------------------
  ! now use evspsbl to account for snow and ice
  !evap(is:ie) = evap(is:ie) + dt*eg(is:ie)/hl        !time integ value in mm (wrong for snow)
  
  ! Update runoff for river routing
  if ( abs(nriver)==1 ) then
    watbdy(is:ie) = watbdy(is:ie) + runoff(is:ie) - oldrunoff(is:ie) ! runoff in mm
  end if
  
  ! correct longwave radiation due to change in tss
  if ( odcalc ) then
    rg_error(1:imax) = stefbo*( tss_save(is:ie)**4 - tss(is:ie)**4 )
    rt(is:ie) = rt(is:ie) + rg_error(1:imax)
    rtclr(is:ie) = rtclr(is:ie) + rg_error(1:imax)
    rgn(is:ie) = rgn(is:ie) - rg_error(1:imax)
    rgclr(is:ie) = rgclr(is:ie) - rg_error(1:imax)
  end if  
  
  
  !***  end of surface updating loop
  ! ----------------------------------------------------------------------
end do
!$omp end do nowait
  
if(diag.or.ntest==1)then
!$omp barrier
!$omp single
  if ( mydiag ) then
    write(6,*) 'at end of sflux'
    write(6,*) 'eg,fg ',eg(idjd),fg(idjd)
    write(6,*) 'tscrn,cduv,zolnd ',tscrn(idjd),cduv(idjd),zolnd(idjd)
    write(6,*) 'snowd,sicedep ',snowd(idjd),sicedep(idjd)
    write(6,*) 'u1,v1,qg1 ',u(idjd,1),v(idjd,1),qg(idjd,1)
    write(6,*) 'w,w2,condx ',wb(idjd,1),wb(idjd,ms),condx(idjd)
    write(6,*) 't1,tss,tgg_2,tgg_ms ',t(idjd,1),tss(idjd),tgg(idjd,2),tgg(idjd,ms)
  end if
  call maxmin(tscrn,'tc',ktau,1.,1)
!$omp end single
endif
if(ntest==4.and.ktau==10)then
!$omp barrier
!$omp single
  do iq=1,ifull
    if(.not.land(iq))then
      write(45,'(2g13.4)') sqrt(u(iq,1)**2+v(iq,1)**2),fg(iq)
      write(46,'(2g13.4)') sqrt(u(iq,1)**2+v(iq,1)**2),eg(iq)
    endif
    write(47,'(2g13.4)') sqrt(u(iq,1)**2+v(iq,1)**2),eg(iq)
  enddo
!$omp end single
endif

return
end subroutine sflux


subroutine sflux_sea

use arrays_m                        ! Atmosphere dyamics prognostic arrays
use cc_mpi                          ! CC MPI routines
use const_phys                      ! Physical constants
use estab                           ! Liquid saturation function
use extraout_m                      ! Additional diagnostics
use morepbl_m                       ! Additional boundary layer diagnostics
use newmpar_m, only : ifull         ! Grid parameters
use parm_m                          ! Model configuration
use pbl_m                           ! Boundary layer arrays
use prec_m, only : evspsbl, sbl     ! Precipitation
use riverarrays_m                   ! River data
use screen_m                        ! Screen level diagnostics
use soil_m                          ! Soil and surface data
use soilsnow_m                      ! Soil, snow and surface data
use work2_m                         ! Diagnostic arrays
use work3_m                         ! Mk3 land-surface diagnostic arrays

implicit none

integer iq, it
real, dimension(ifull) :: charnck,fgf,rgg,fev,dirad,dfgdt,degdt
real, dimension(ifull) :: cie,gamm,fh,factch,evap_sea,evap_ice
real afrootpan,es,constz,xx,afroot,fm,consea,con,daf,den,dden,dfm
real dtsol,con1,conw,zminlog,drst,ri_ice,factchice,zoice,zologice
real conh,epotice,qtgnet,qtgair,eg2,gbot,deltat,b1,deg,eg1,denha
real denma,root,dcs,fq,afq,facqch,denqa

integer, parameter :: ntss_sh=0 ! 0 for original, 3 for **3, 4 for **4
! Beljaars parameters
! ..... high wind speed - rough sea
real, parameter :: zcom1 = 0.018     ! Charnock's constant
real, parameter :: zcoh1 = 0.0       ! Beljaars 1994 values
real, parameter :: zcoq1 = 0.0
! ..... low wind speed - smooth sea
real, parameter :: gnu   = 1.5e-5
real, parameter :: zcom2 = 0.11
real, parameter :: zcoh2 = 0.40
real, parameter :: zcoq2 = 0.62

gamm=3.471e+05 ! dummy value

if(ntest==2.and.mydiag)write(6,*) 'before sea loop'                                              ! sea
! from June '03 use basic sea temp from tgg1 (so leads is sensible)                              ! sea
! all sea points in this loop; also open water of leads                                          ! sea
if( charnock>0. )then                                                                            ! sea
  charnck(:)=charnock                                                                            ! sea
elseif ( charnock<-10. ) then ! Beljaars                                                         ! sea
  charnck(:)=zcom1                                                                               ! sea
elseif(charnock<-1.)then  ! zo like Moon (2004)                                                  ! sea
  charnck(:)=max(.0000386*u10(:),.000085*u10(:)-.00058)                                          ! sea
else                      ! like Makin (2002)                                                    ! sea
  charnck(:)=.008+3.e-4*(u10(:)-9.)**2/(1.+(.006+.00008*u10(:))*u10(:)**2)                       ! sea
endif                                                                                            ! sea
do iq=1,ifull                                                                                    ! sea
  if(.not.land(iq))then                                                                          ! sea
    wetfac(iq)=1.                                                                                ! sea
    ! tpan holds effective sea for this loop                                                     ! sea
    if(ntss_sh==0)then                                                                           ! sea
      dtsol=tss_sh*.01*sgsave(iq)/(1.+.25*vmod(iq)**2)  ! solar heating                          ! sea
      tpan(iq)=tgg(iq,1)+min(dtsol,8.)                  ! of ssts                                ! sea
    elseif(ntss_sh==3)then                                                                       ! sea
      dtsol=tss_sh*.01*sgsave(iq)/(1.+.035*vmod(iq)**3) ! solar heating                          ! sea
      tpan(iq)=tgg(iq,1)+min(dtsol,8.)                  ! of ssts                                ! sea
    elseif(ntss_sh==4)then                                                                       ! sea
      dtsol=tss_sh*.01*sgsave(iq)/(1.+vmod(iq)**4/81.) ! solar heating                           ! sea
      tpan(iq)=tgg(iq,1)+min(dtsol,8.)                 ! of ssts                                 ! sea
    endif   ! (ntss_sh==0) .. else ..                                                            ! sea
    if(ntsea==1.and.condx(iq)>.1)tpan(iq)=t(iq,2)                                                ! sea
    if(ntsea==2.and.condx(iq)>.1)tpan(iq)=t(iq,1)                                                ! sea
    if(ntsea==3.and.condx(iq)>.1)tpan(iq)=.5*(t(iq,2)+tgg(iq,1))                                 ! sea
    if(ntsea==4.and.condx(iq)>.1)tpan(iq)=.5*(t(iq,1)+tgg(iq,1))                                 ! sea
  endif  ! (.not.land(iq))                                                                       ! sea
enddo   ! iq loop                                                                                ! sea
                                                                                                 ! sea
! here calculate fluxes for sea point, and nominal pan points                                    ! sea
afrootpan=vkar/log(zmin/panzo)                                                                   ! sea
do iq=1,ifull                                                                                    ! sea
  ! drag coefficients  for momentum cduv                                                         ! sea
  ! for heat and moisture  cdtq                                                                  ! sea
  es = establ(tpan(iq))                                                                          ! sea
  constz=max(ps(iq)-es,0.1)                                                                      ! sea
  qsttg(iq)= .98*.622*es/constz  ! with Zeng 1998 for sea water                                  ! sea
  xx=grav*zmin*(1.-tpan(iq)*srcp/t(iq,1))                                                        ! sea
  ri(iq)=min(xx/vmag(iq)**2 , ri_max)                                                            ! sea
  ! this is in-line ocenzo using latest coefficient, i.e. .018                                   ! sea
  if(land(iq))then                                                                               ! sea
    zo(iq)=panzo                                                                                 ! sea
    af(iq)=afrootpan**2                                                                          ! sea
    aft(iq)=chnsea                                                                               ! sea
    afq=chnsea                                                                                   ! sea
    factch(iq)=1.                                                                                ! sea
    facqch=1.                                                                                    ! sea
  else                                                                                           ! sea
    if ( charnock<-10. ) then ! Beljarrs                                                         ! sea
      zo(iq)=.001    ! .0005 better first guess                                                  ! sea
      if(ri(iq)>0.)then             ! stable sea points                                          ! sea
        do it=1,4                                                                                ! sea
          afroot=vkar/log(zmin/zo(iq))                                                           ! sea
          af(iq)=afroot**2                                                                       ! sea
          daf=2.*af(iq)*afroot/(vkar*zo(iq))                                                     ! sea
          fm=1./(1.+bprm*ri(iq))**2                                                              ! sea
          consea=zcom1*vmod(iq)**2*af(iq)*fm/grav+zcom2*gnu/(vmod(iq)*sqrt(fm*af(iq)))           ! sea
          dcs=(zcom1*vmod(iq)**2/grav                                       &                    ! sea
              -0.5*zcom2*gnu/(vmod(iq)*sqrt(fm*af(iq))*fm*af(iq)))*(fm*daf)                      ! sea
          zo(iq)=zo(iq)-(zo(iq)-consea)/(1.-dcs)                                                 ! sea
          zo(iq)=max(1.5e-5,zo(iq))                                                              ! sea
          zo(iq)=min(zo(iq),6.)                                                                  ! sea
        enddo    ! it=1,4                                                                        ! sea
        afroot=vkar/log(zmin/zo(iq))                                                             ! sea
        af(iq)=afroot**2                                                                         ! sea
      else                        ! unstable sea points                                          ! sea
        do it=1,4                                                                                ! sea
          afroot=vkar/log(zmin/zo(iq))                                                           ! sea
          af(iq)=afroot**2                                                                       ! sea
          daf=2.*af(iq)*afroot/(vkar*zo(iq))                                                     ! sea
          con=cms*2.*bprm*sqrt(-ri(iq)*zmin/zo(iq))                                              ! sea
          den=1.+af(iq)*con                                                                      ! sea
          fm=1.-2.*bprm*ri(iq)/den                                                               ! sea
          dfm=2.*bprm*ri(iq)*(con*daf+af(iq)*cms*bprm*sqrt(-ri(iq))*zmin &                       ! sea
              /(sqrt(zmin/zo(iq))*zo(iq)**2))/(den**2)                                           ! sea
          consea=zcom1*vmod(iq)**2*af(iq)*fm/grav+zcom2*gnu/(vmod(iq)*sqrt(fm*af(iq)))           ! sea
          dcs=(zcom1*vmod(iq)**2/grav                                              &             ! sea
              -0.5*zcom2*gnu/(vmod(iq)*sqrt(fm*af(iq))*fm*af(iq)))*(fm*daf-dfm*af(iq))           ! sea
          zo(iq)=zo(iq)-(zo(iq)-consea)/(1.-dcs)                                                 ! sea
          zo(iq)=max(1.5e-5,zo(iq))                                                              ! sea
          zo(iq)=min(zo(iq),6.)                                                                  ! sea
        enddo  ! it=1,4                                                                          ! sea
        afroot=vkar/log(zmin/zo(iq))                                                             ! sea
        af(iq)=afroot**2                                                                         ! sea
      endif    ! (xx>0.) .. else..                                                               ! sea        
      zoh(iq)=max(zcoh1+zcoh2*gnu/(vmod(iq)*sqrt(fm*af(iq))),1.5e-7)                             ! sea
      zoq(iq)=max(zcoq1+zcoq2*gnu/(vmod(iq)*sqrt(fm*af(iq))),1.5e-7)                             ! sea
      aft(iq)=vkar**2/(log(zmin/zo(iq))*log(zmin/zoh(iq)))                                       ! sea
      afq=vkar**2/(log(zmin/zo(iq))*log(zmin/zoq(iq)))                                           ! sea
      factch(iq)=sqrt(zo(iq)/zoh(iq))                                                            ! sea
      facqch=sqrt(zo(iq)/zoq(iq))                                                                ! sea
    else if(charnock<0.)then  ! Moon (2004) over sea                                             ! sea
      zo(iq)=charnck(iq)                                                                         ! sea
      afroot=vkar/log(zmin/zo(iq))                                                               ! sea
      af(iq)=afroot**2                                                                           ! sea
      aft(iq)=chnsea                                                                             ! sea
      afq=chnsea                                                                                 ! sea
      if ( newztsea==0 ) then ! 0 for original, 1 for different zt over sea                      ! sea
        ! enhanced formula used in Feb '92 Air-Sea conference follows:                           ! sea
        ! factch=sqrt(zo*exp(vkar*vkar/(chnsea*log(zmin/zo)))/zmin)                              ! sea
        factch(iq)=1. ! factch is sqrt(zo/zt) only for use in unstable fh                        ! sea
      else                                                                                       ! sea
        factch(iq)=sqrt(zo(iq)*ztv) ! for use in unstable fh                                     ! sea
      end if                                                                                     ! sea
      facqch=factch(iq)                                                                          ! sea
    else            ! usual charnock method over sea                                             ! sea
      consea=vmod(iq)*charnck(iq)/grav  ! usually charnock=.018                                  ! sea  
      zo(iq)=.001    ! .0005 better first guess                                                  ! sea
      if(ri(iq)>0.)then             ! stable sea points                                          ! sea
        fm=vmod(iq) /(1.+bprm*ri(iq))**2 ! N.B. this is vmod*fm                                  ! sea
        con=consea*fm                                                                            ! sea
        do it=1,3                                                                                ! sea
          afroot=vkar/log(zmin/zo(iq))                                                           ! sea
          af(iq)=afroot**2                                                                       ! sea
          daf=2.*af(iq)*afroot/(vkar*zo(iq))                                                     ! sea
          zo(iq)=max(1.5e-5,zo(iq)-(zo(iq)-con*af(iq))/(1.-con*daf))                             ! sea
          zo(iq)=min(zo(iq),9.) ! JLM fix                                                        ! sea
        enddo    ! it=1,3                                                                        ! sea
        afroot=vkar/log(zmin/zo(iq))                                                             ! sea
        af(iq)=afroot**2                                                                         ! sea
      else                        ! unstable sea points                                          ! sea
        do it=1,3                                                                                ! sea
          afroot=vkar/log(zmin/zo(iq))                                                           ! sea
          af(iq)=afroot**2                                                                       ! sea
          daf=2.*af(iq)*afroot/(vkar*zo(iq))                                                     ! sea
          con1=cms*2.*bprm*sqrt(-ri(iq)*zmin/zo(iq))                                             ! sea
          den=1.+af(iq)*con1                                                                     ! sea
          dden=con1*(daf-.5*af(iq)/zo(iq))                                                       ! sea
          fm=vmod(iq)-vmod(iq)*2.*bprm *ri(iq)/den                                               ! sea
          dfm=2.*bprm*ri(iq)*dden/den**2                                                         ! sea
          zo(iq)=max(1.5e-5,zo(iq)-(zo(iq)-consea*af(iq)*fm)/(1.-consea*(daf*fm+af(iq)*dfm)))    ! sea
          zo(iq)=min(zo(iq),6.) ! JLM fix                                                        ! sea
        enddo  ! it=1,3                                                                          ! sea
      endif    ! (xx>0.) .. else..                                                               ! sea
      aft(iq)=chnsea                                                                             ! sea
      afq=chnsea                                                                                 ! sea
      if ( newztsea==0 ) then ! 0 for original, 1 for different zt over sea                      ! sea
        ! enhanced formula used in Feb '92 Air-Sea conference follows:                           ! sea
        ! factch=sqrt(zo*exp(vkar*vkar/(chnsea*log(zmin/zo)))/zmin)                              ! sea
        factch(iq)=1. ! factch is sqrt(zo/zt) only for use in unstable fh                        ! sea
      else                                                                                       ! sea
        factch(iq)=sqrt(zo(iq)*ztv) ! for use in unstable fh                                     ! sea
      end if                                                                                     ! sea
      facqch=factch(iq)                                                                          ! sea
    endif     ! (charnock<-1.) .. else ..                                                        ! sea
  endif      ! (land(iq)) .. else ..                                                             ! sea

  ! Having settled on zo & af now do actual fh and fm calcs                                      ! sea
  if(ri(iq)>0.)then                                                                              ! sea
    fm=vmod(iq)/(1.+bprm*ri(iq))**2  ! no zo contrib for stable                                  ! sea
    fh(iq)=fm                                                                                    ! sea
    fq=fm                                                                                        ! sea
  else        ! ri is -ve                                                                        ! sea
    root=sqrt(-ri(iq)*zmin/zo(iq))                                                               ! sea
    ! First do momentum                                                                          ! sea
    denma=1.+cms*2.*bprm*af(iq)*root                                                             ! sea
    fm=vmod(iq)-vmod(iq)*2.*bprm *ri(iq)/denma                                                   ! sea
    ! n.b. fm denotes ustar**2/(vmod(iq)*af)                                                     ! sea
    ! Now heat ; allow for smaller zo via aft and factch                                         ! sea
    ! N.B. for newztsea=1, zo contrib cancels in factch*root,                                    ! sea
    ! so eg (& epan) and fg  (also aft) then indept of zo                                        ! sea
    denha=1.+chs*2.*bprm*factch(iq)*aft(iq)*root                                                 ! sea
    fh(iq)=vmod(iq)-vmod(iq)*2.*bprm *ri(iq)/denha                                               ! sea
    denqa=1.+chs*2.*bprm*facqch*afq*root                                                         ! sea
    fq=vmod(iq)-vmod(iq)*2.*bprm *ri(iq)/denqa                                                   ! sea
  endif                                                                                          ! sea
                                                                                                 ! sea
  conh=rho(iq)*aft(iq)*cp                                                                        ! sea
  conw=rho(iq)*afq*hl                                                                            ! sea
  fg(iq)=conh*fh(iq)*(tpan(iq)-theta(iq))                                                        ! sea
  eg(iq)=conw*fq*(qsttg(iq)-qg(iq,1))                                                            ! sea
  evap_sea(iq)=eg(iq)/hl                                                                         ! sea
  rnet(iq)=sgsave(iq)-rgsave(iq)-stefbo*tpan(iq)**4                                              ! sea
  zoh(iq)=zo(iq)/(factch(iq)**2)                                                                 ! sea
  zoq(iq)=zo(iq)/(facqch**2)                                                                     ! sea
  ! cduv is now drag coeff *vmod                                                                 ! sea
  cduv(iq) =af(iq)*fm                                                                            ! sea
  cdtq(iq) =aft(iq)*fh(iq)                                                                       ! sea
  ustar(iq) = sqrt(vmod(iq)*cduv(iq))                                                            ! sea
  ! Surface stresses taux, tauy: diagnostic only - unstaggered now                               ! sea
  taux(iq)=rho(iq)*cduv(iq)*u(iq,1)                                                              ! sea
  tauy(iq)=rho(iq)*cduv(iq)*v(iq,1)                                                              ! sea
  fg_ocn(iq)=fg(iq)                                                                              ! sea
  eg_ocn(iq)=eg(iq)                                                                              ! sea
  taux_ocn(iq)=taux(iq)                                                                          ! sea
  tauy_ocn(iq)=tauy(iq)                                                                          ! sea
  ! note that iq==idjd  can only be true on the correct processor                                ! sea
  if(ntest==1.and.iq==idjd.and.mydiag)then                                                       ! sea
    write(6,*) 'in sea-type loop for iq,idjd: ',iq,idjd                                          ! sea
    write(6,*) 'zmin,zo,factch ',zmin,zo(iq),factch(iq)                                          ! sea
    write(6,*) 'ri,ustar,es ',ri(iq),ustar(iq),es                                                ! sea
    write(6,*) 'af,aft ',af(iq),aft(iq)                                                          ! sea
    write(6,*) 'tpan,tss,theta ',tpan(iq),tss(iq),theta(iq)                                      ! sea
    write(6,*) 'chnsea,rho,t1 ',chnsea,rho(iq),t(iq,1)                                           ! sea
    write(6,*) 'fm,fh,conh ',fm,fh(iq),conh                                                      ! sea
    write(6,*) 'vmod,cduv,fg ',vmod(iq),cduv(iq),fg(iq)                                          ! sea
  endif                                                                                          ! sea
enddo     ! iq loop                                                                              ! sea
epot(:) = eg(:)                                                                                  ! sea
epan(:) = eg(:)                                                                                  ! sea
! section to update pan temperatures                                                             ! sea
do iq=1,ifull                                                                                    ! sea
  if(land(iq))then                                                                               ! sea
    rgg(iq)=5.67e-8*tpan(iq)**4                                                                  ! sea
    ! assume gflux = 0                                                                           ! sea
    ! note pan depth=.254 m, spec heat water=4186 joule/kg K                                     ! sea
    ! and change in heat supplied=spec_heatxmassxdelta_T                                         ! sea
    ga(iq)=-slwa(iq)-rgg(iq)-panfg*fg(iq)                                                        ! sea
    tpan(iq)=tpan(iq)+ga(iq)*dt/(4186.*.254*1000.)                                               ! sea
  else                                                                                           ! sea
    sno(iq)=sno(iq)+conds(iq)                                                                    ! sea
    grpl(iq)=grpl(iq)+condg(iq)                                                                  ! sea
  endif  ! (land(iq))                                                                            ! sea
enddo   ! iq loop                                                                                ! sea
                                                                                                 ! sea
if(nmaxpr==1.and.mydiag)then                                                                     ! sea
  iq=idjd                                                                                        ! sea
  write (6,"('after sea loop fg,tpan,epan,ri,fh,vmod',9f9.4)")           &                       ! sea
      fg(idjd),tpan(idjd),epan(idjd),ri(idjd),fh(idjd),vmod(idjd)                                ! sea
  write (6,"('u10,ustar,charnck,zo,cd',3f9.4,2f9.6)")                    &                       ! sea
      u10(idjd),ustar(idjd),charnck(idjd),zo(idjd),cduv(idjd)/vmod(idjd)                         ! sea
  if(ri(iq)>0.)then                                                                              ! sea
    fm=vmod(iq)/(1.+bprm*ri(iq))**2                                                              ! sea
  else                                                                                           ! sea
    root=sqrt(-ri(iq)*zmin/zo(iq))                                                               ! sea
    denma=1.+cms*2.*bprm*af(iq)*root                                                             ! sea
    fm=vmod(iq)-vmod(iq)*2.*bprm *ri(iq)/denma                                                   ! sea
    denha=1.+chs*2.*bprm*factch(iq)*aft(iq)*root                                                 ! sea
  endif                                                                                          ! sea
  write (6,"('after sea loop af,aft,factch,root,denha,denma,fm',9f9.4)") &                       ! sea
      af(iq),aft(iq),factch(iq),root,denha,denma,fm                                              ! sea
endif                                                                                            ! sea
zminlog=log(zmin)                                                                                ! sea

fgf=0.                                                                                           ! sice
fev=0.                                                                                           ! sice
evap_ice=0.                                                                                      ! sice
do iq=1,ifull                                                                                    ! sice
  if(sicedep(iq)>0.)then                                                                         ! sice
    ! non-leads for sea ice points                                                               ! sice
    ! N.B. tggsn( ,1) holds tice                                                                 ! sice
    es = establ(tggsn(iq,1))                                                                     ! sice
    constz=max(ps(iq)-es,0.1)                                                                    ! sice
    qsttg(iq)= .622*es/constz                                                                    ! sice
    drst=qsttg(iq)*ps(iq)*hlars/(tggsn(iq,1)**2*constz)                                          ! sice
    xx=grav*zmin*(1.-tggsn(iq,1)*srcp/t(iq,1))                                                   ! sice
    ri_ice=min(xx/vmag(iq)**2 , ri_max)                                                          ! sice
    !factchice=1. ! factch is sqrt(zo/zt) for use in unstable fh                                 ! sice
    factchice=sqrt(7.4) ! same as land from 27/4/99                                              ! sice
    zoice=.001                                                                                   ! sice
    zologice=zminlog-log(zoice)   !   i.e. log(zmin/zo(iq))                                      ! sice
    af(iq)=(vkar/zologice)**2                                                                    ! sice
    ! MJT notes - following line does not agree with factchice                                   ! sice
    aft(iq)=vkar**2/(zologice*zologice )                                                         ! sice
    wetfac(iq)=1+.008*(tggsn(iq,1)-273.16)                                                       ! sice
                                                                                                 ! sice
    ! now do fh and fm calcs for sice                                                            ! sice
    if(ri_ice>0.)then                                                                            ! sice
      fm=vmod(iq)/(1.+bprm*ri_ice)**2                                                            ! sice
      fh(iq)=fm                                                                                  ! sice
    else                                                                                         ! sice
      root=sqrt(-ri_ice*zmin/zoice)                                                              ! sice
      ! First do momentum                                                                        ! sice
      denma=1.+cms*2.*bprm*af(iq)*root                                                           ! sice
      fm=vmod(iq)-vmod(iq)*2.*bprm *ri_ice/denma                                                 ! sice
      ! n.b. fm denotes ustar**2/(vmod(iq)*af)                                                   ! sice
      ! Now heat ; allow for smaller zo via aft and factch                                       ! sice
      denha=1.+chs*2.*bprm*factchice*aft(iq)*root                                                ! sice
      fh(iq)=vmod(iq)-(2.*bprm *ri_ice)/denha                                                    ! sice
    endif                                                                                        ! sice
                                                                                                 ! sice
    conh=rho(iq)*aft(iq)*cp                                                                      ! sice
    conw=rho(iq)*aft(iq)*hl                                                                      ! sice
    ! fgice & egice renamed as fgf and fev from Aug 05 to aid diags                              ! sice	
    fgf(iq)=conh*fh(iq)*(tggsn(iq,1)-theta(iq))                                                  ! sice
    dfgdt(iq)=conh*fh(iq)                                                                        ! sice
    if(ntest==1.and.iq==idjd.and.mydiag)then                                                     ! sice
      write(6,*) 'in sice loop'                                                                  ! sice
      write(6,*) 'zmin,zo,wetfac ',zmin,zoice,wetfac(iq)                                         ! sice
      write(6,*) 'ri_ice,es ',ri_ice,es                                                          ! sice
      write(6,*) 'af,aft,ustar ',af(iq),aft(iq),ustar(iq)                                        ! sice
      write(6,*) 'chnsea,rho ',chnsea,rho(iq)                                                    ! sice
      write(6,*) 'fm,fh,conh ',fm,fh(iq),conh                                                    ! sice
    endif                                                                                        ! sice
                                                                                                 ! sice
    if(nalpha==1)then    ! beta scheme         sice here                                         ! sice
      epotice=conw*fh(iq)*(qsttg(iq)-qg(iq,1))                                                   ! sice
      fev(iq)=wetfac(iq)*epotice                                                                 ! sice
      degdt(iq)=wetfac(iq)*conw*fh(iq)*drst                                                      ! sice
    else                   ! alpha scheme                                                        ! sice
      ! following trick reduces -ve evap (dew) to 1/10th value                                   ! sice
      qtgnet=qsttg(iq)*wetfac(iq) -qg(iq,1)                                                      ! sice
      qtgair=qsttg(iq)*wetfac(iq)-max(qtgnet,.1*qtgnet)                                          ! sice
      eg2=-conw*fh(iq)*qtgair                                                                    ! sice
      eg1=conw*fh(iq)*qsttg(iq)                                                                  ! sice
      fev(iq) = eg1*wetfac(iq) + eg2                                                             ! sice
      epotice = conw*fh(iq)*(qsttg(iq)-qg(iq,1))                                                 ! sice
      deg=wetfac(iq)*conw*fh(iq)*drst                                                            ! sice
      ! following reduces degdt by factor of 10 for dew                                          ! sice
      degdt(iq)=.55*deg+sign(.45*deg,qtgnet)                                                     ! sice
    endif                                                                                        ! sice
    evap_ice(iq)=fev(iq)/hls                                                                     ! sice
                                                                                                 ! sice
    ! section to update sea ice surface temperature;                                             ! sice
    ! specified sea-ice thickness                                                                ! sice
    ! over sea ice, set a minimum depth for this experiment of .1                                ! sice
    ! sicedep(iq) = max(sicedep(iq) , 0.1)  fixed in indata/nestin from Jan 06                   ! sice
    ! no snow on the ice assumed for now                                                         ! sice
    gamm(iq) = 3.471e+05                                                                         ! sice
    cie(iq) = 2.04/sicedep(iq)                                                                   ! sice
    rgg(iq) = 5.67e-8*tggsn(iq,1)**4                                                             ! sice
    ! gflux here is flux from ice to water, +ve downwards                                        ! sice
    gflux(iq)=cie(iq)*(tggsn(iq,1)-271.2)                                                        ! sice
    ga(iq)=-slwa(iq)-rgg(iq)-fev(iq)-fgf(iq)-gflux(iq)                                           ! sice
    dirad(iq)=4.*5.67e-8*tggsn(iq,1)**3                                                          ! sice
    b1=dirad(iq)+degdt(iq)+dfgdt(iq)+cie(iq)                                                     ! sice
    gbot=(gamm(iq)/dt)+b1                                                                        ! sice
    deltat=ga(iq)/gbot                                                                           ! sice
    tggsn(iq,1)=tggsn(iq,1)+deltat                                                               ! sice
    tggsn(iq,1)=min(tggsn(iq,1),271.2)   ! jlm fix Tue  05-30-2000                               ! sice
    fgf(iq) =fgf(iq) +deltat*dfgdt(iq)                                                           ! sice
    fev(iq) =fev(iq) +deltat*degdt(iq)                                                           ! sice
    es = establ(tggsn(iq,1))                                                                     ! sice
    constz=max(ps(iq)-es,0.1)                                                                    ! sice
    qsttg(iq)=.622*es/constz                                                                     ! sice
                                                                                                 ! sice
    ! combine ice and leads contributions here                                                   ! sice
    eg(iq) =fracice(iq)*fev(iq) + (1.-fracice(iq))*eg(iq)                                        ! sice
    fg(iq) =fracice(iq)*fgf(iq) + (1.-fracice(iq))*fg(iq)                                        ! sice
    ri(iq) =fracice(iq)*ri_ice  + (1.-fracice(iq))*ri(iq)                                        ! sice
    zo(iq) =fracice(iq)*zoice   + (1.-fracice(iq))*zo(iq)                                        ! sice
    factch(iq) = fracice(iq)*factchice + (1.-fracice(iq))*factch(iq)                             ! sice
    zoh(iq)=zo(iq)/(factch(iq)*factch(iq))                                                       ! sice
    zoq(iq)=zoh(iq)                                                                              ! sice
    cduv(iq) =fracice(iq)*af(iq)*fm + (1.-fracice(iq))*cduv(iq)                                  ! sice
    cdtq(iq) =fracice(iq)*aft(iq)*fh(iq)+(1.-fracice(iq))*cdtq(iq)                               ! sice
    ustar(iq) = sqrt(vmod(iq)*cduv(iq))                                                          ! sice
    ! N.B. potential evaporation is now eg+eg2                                                   ! sice
    epot(iq) =fracice(iq)*epotice + (1.-fracice(iq))*epot(iq)                                    ! sice
    tss(iq) = fracice(iq)*tggsn(iq,1)+(1.-fracice(iq))*tpan(iq)                                  ! sice
    rnet(iq)=sgsave(iq)-rgsave(iq)-stefbo*tss(iq)**4                                             ! sice
    ! Surface stresses taux, tauy: diagnostic only - unstag now                                  ! sice
    taux(iq)=rho(iq)*cduv(iq)*u(iq,1)                                                            ! sice
    tauy(iq)=rho(iq)*cduv(iq)*v(iq,1)                                                            ! sice
    fg_ice(iq)=fgf(iq)                                                                           ! sice
    eg_ice(iq)=fev(iq)                                                                           ! sice
    taux_ice(iq)=rho(iq)*af(iq)*fm*u(iq,1)                                                       ! sice
    tauy_ice(iq)=rho(iq)*af(iq)*fm*v(iq,1)                                                       ! sice
  endif  ! (sicedep(iq)>0.)                                                                      ! sice
enddo       ! iq loop                                                                            ! sice
where (.not.land)                                                                                ! sice
  snowd=0.                                                                                       ! sice
  evspsbl = evspsbl + dt*(fracice*evap_ice + (1.-fracice)*evap_sea)                              ! sice
  sbl = sbl + dt*fracice*evap_ice                                                                ! sice
end where                                                                                        ! sice
if(mydiag.and.nmaxpr==1)then                                                                     ! sice
  write(6,*) 'after sice loop'                                                                   ! sice
  iq=idjd                                                                                        ! sice
  if(sicedep(iq)>0.)then                                                                         ! sice
    b1=dirad(iq)+degdt(iq)+dfgdt(iq)+cie(iq)                                                     ! sice
    gbot=(gamm(iq)/dt)+b1                                                                        ! sice
    deltat=ga(iq)/gbot                                                                           ! sice
    write(6,*) 'ri,vmag,vmod,cduv ',ri(iq),vmag(iq),vmod(iq),cduv(iq)                            ! sice
    write(6,*) 'fh,tss,tpan,tggsn1 ',fh(iq),tss(iq),tpan(iq),tggsn(iq,1)                         ! sice
    write(6,*) 'theta,t1,deltat ',theta(iq),t(iq,1),deltat                                       ! sice
    write(6,*) 'b1,ga,gbot,af,aft ',b1,ga(iq),gbot,af(iq),aft(iq)                                ! sice
    write(6,*) 'fg,fgice,factch ',fg(iq),fgf(iq),factch(iq)                                      ! sice
    write(6,*) 'cie ',cie(iq)                                                                    ! sice
    write(6,*) 'eg,egice(fev),ustar ',eg(iq),fev(iq),ustar(iq)                                   ! sice
  end if ! sicedep(iq)>0.                                                                        ! sice
endif    ! (mydiag.and.nmaxpr==1)                                                                ! sice

if ( abs(nriver)==1 ) then
  where ( .not.land(1:ifull) )
    river_inflow(1:ifull) = watbdy(1:ifull)  
    watbdy(1:ifull) = 0. ! water enters ocean and is removed from rivers
  end where
end if  

return
end subroutine sflux_sea

subroutine sflux_mlo

use arrays_m                       ! Atmosphere dyamics prognostic arrays
use cc_mpi                         ! CC MPI routines
use cc_omp                         ! CC OpenMP routines
use const_phys                     ! Physical constants
use estab                          ! Liquid saturation function
use extraout_m                     ! Additional diagnostics
use map_m                          ! Grid map arrays
use mlo                            ! Ocean physics and prognostic arrays
use mlodynamicsarrays_m            ! Ocean dynamics data
use morepbl_m                      ! Additional boundary layer diagnostics
use newmpar_m                      ! Grid parameters
use parm_m                         ! Model configuration
use pbl_m                          ! Boundary layer arrays
use permsurf_m                     ! Fixed surface arrays
use prec_m, only : evspsbl, sbl    ! Precipitation
use raddiag_m                      ! Radiation diagnostic
use riverarrays_m                  ! River data
use soil_m                         ! Soil and surface data
use soilsnow_m                     ! Soil, snow and surface data
use work2_m                        ! Diagnostic arrays
use work3_m                        ! Mk3 land-surface diagnostic arrays

implicit none

integer tile, is, ie, k
real, dimension(imax) :: lwatbdy
logical, dimension(imax) :: loutflowmask

!$omp do schedule(static) private(is,ie,k),                                         &
!$omp private(lwatbdy,loutflowmask)
do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax

  if ( nriver/=0 ) then
    lwatbdy      = watbdy(is:ie)
    loutflowmask = outflowmask(is:ie)
  end if
  
  call sflux_mlo_work(ri(is:ie),srcp,vmag(is:ie),ri_max,bprm,chs,ztv,chnsea,rho(is:ie),azmin(is:ie),     &
                      uav(is:ie),vav(is:ie),                                                             &
                      ps(is:ie),t(is:ie,1),qg(is:ie,1),sgdn(is:ie),sgsave(is:ie),rgsave(is:ie),          &
                      swrsave(is:ie),fbeamvis(is:ie),fbeamnir(is:ie),taux(is:ie),tauy(is:ie),            &
                      ustar(is:ie),f(is:ie),water_g(tile),wpack_g(:,tile),wfull_g(tile),depth_g(tile),   &
                      dgice_g(tile),dgscrn_g(tile),dgwater_g(tile),ice_g(tile),                          &
                      tpan(is:ie),epan(is:ie),rnet(is:ie),condx(is:ie),                                  &
                      conds(is:ie),condg(is:ie),fg(is:ie),eg(is:ie),evspsbl(is:ie),sbl(is:ie),           &
                      epot(is:ie),tss(is:ie),cduv(is:ie),cdtq(is:ie),lwatbdy,loutflowmask,land(is:ie),   &
                      fracice(is:ie),sicedep(is:ie),snowd(is:ie),sno(is:ie),grpl(is:ie),qsttg(is:ie),    &
                      vmod(is:ie),zo(is:ie),wetfac(is:ie),                                               &
                      zoh(is:ie),zoq(is:ie),theta(is:ie),ga(is:ie),turb_g(tile))
  
  if ( nriver/=0 ) then
    watbdy(is:ie) = lwatbdy
  end if
  
  do k = 1,ms 
    call mloexport(0,tgg(is:ie,k),k,0,water_g(tile),wpack_g(:,tile),wfull_g(tile)) 
    where ( tgg(is:ie,k)<100. )     
      tgg(is:ie,k) = tgg(is:ie,k) + wrtemp 
    end where                 
  end do                      
  do k = 1,3                   
    call mloexpice(tggsn(is:ie,k),k,0,ice_g(tile),wpack_g(:,tile),wfull_g(tile)) 
  end do 

end do
!$omp end do nowait

return
end subroutine sflux_mlo

subroutine sflux_mlo_work(ri,srcp,vmag,ri_max,bprm,chs,ztv,chnsea,rho,azmin,uav,vav,      &
                          ps,t,qg,sgdn,sgsave,rgsave,swrsave,fbeamvis,fbeamnir,taux,tauy, &
                          ustar,f,water,wpack,wfull,depth,dgice,dgscrn,dgwater,ice,       &
                          tpan,epan,rnet,condx,conds,condg,fg,eg,evspsbl,sbl,epot,        &
                          tss,cduv,cdtq,watbdy,outflowmask,land,                          &
                          fracice,sicedep,snowd,sno,grpl,qsttg,vmod,zo,wetfac,            &
                          zoh,zoq,theta,ga,turb)
use cc_mpi                           ! CC MPI routines
use cc_omp                           ! CC OpenMP routines
use const_phys                       ! Physical constants
use estab                            ! Liquid saturation function
use mlo, only : waterdata,icedata, & ! Ocean physics and prognostic arrays
  dgwaterdata,dgicedata,dgscrndata,& 
  depthdata,mloeval,               & 
  mloexport,mloimport,mloextra,    &
  mloexpice,turbdata
use newmpar_m                        ! Grid parameters
use parm_m                           ! Model configuration
use soil_m, only : zmin              ! Soil and surface data

implicit none

integer iq
integer, intent(in) :: wfull
real root, denha, esatf
real, intent(in) :: srcp, ri_max, bprm, chs, ztv, chnsea
real, dimension(imax), intent(in) :: t, qg
real, dimension(imax), intent(inout) :: ri, taux, tauy, ustar, tpan, epan, rnet
real, dimension(imax), intent(inout) :: fg, eg, evspsbl, sbl, epot, tss, cduv, cdtq, watbdy, fracice, sicedep
real, dimension(imax), intent(inout) :: snowd, sno, grpl, qsttg, zo, wetfac, zoh, zoq, ga
real, dimension(imax), intent(in) :: vmag, rho, azmin, uav, vav, ps, sgdn, sgsave, rgsave, swrsave
real, dimension(imax), intent(in) :: fbeamvis, fbeamnir, f, condx, conds, condg, vmod, theta
real, dimension(imax) :: neta, oflow, dumw, dumrg, dumx, dums, fhd
real, dimension(imax) :: umod, umag, levspsbl, lsbl
logical, dimension(imax), intent(in) :: wpack, outflowmask, land
type(waterdata), intent(inout) :: water
type(dgicedata), intent(inout) :: dgice
type(dgscrndata), intent(inout) :: dgscrn
type(dgwaterdata), intent(inout) :: dgwater
type(icedata), intent(inout) :: ice
type(depthdata), intent(in) :: depth
type(turbdata), intent(inout) :: turb

umod = vmod ! this is vmod, but accounts for moving surface                                    ! MLO
                                                                                               ! MLO
if ( abs(nmlo)==1 ) then                                                                       ! MLO
  ! Single column                                                                              ! MLO
  ! set free surface to zero when water is not conserved                                       ! MLO
  neta=0.                                                                                      ! MLO
  call mloimport(4,neta,0,0,water,depth,wpack,wfull)                                           ! MLO
end if                                                                                         ! MLO
                                                                                               ! MLO
! inflow and outflow model for rivers                                                          ! MLO
if ( abs(nmlo)>=2 ) then                                                                       ! MLO
  dumw(1:imax) = 0.                                                                            ! MLO
  where ( .not.land(1:imax) )                                                                  ! MLO
    dumw(1:imax) = watbdy(1:imax)/dt                                                           ! MLO
    watbdy(1:imax) = 0.                                                                        ! MLO
  end where                                                                                    ! MLO
  neta(1:imax) = 0.                                                                            ! MLO
  call mloexport(4,neta,0,0,water,wpack,wfull)                                                 ! MLO
  where ( outflowmask(1:imax) )                                                                ! MLO
    oflow(:) = max( neta(1:imax), 0. )                                                         ! MLO
    watbdy(1:imax) = watbdy(1:imax) + 1000.*oflow(:)                                           ! MLO
    neta(1:imax) = neta(1:imax) - oflow(:)                                                     ! MLO
  end where                                                                                    ! MLO
  call mloimport(4,neta,0,0,water,depth,wpack,wfull)                                           ! MLO
else                                                                                           ! MLO
  dumw(1:imax) = 0.                                                                            ! MLO
end if                                                                                         ! MLO
                                                                                               ! MLO
! Ocean mixing                                                                                 ! MLO
where (.not.land(1:imax))                                                                      ! MLO
  rnet(:)=sgsave(:)-rgsave(:)-stefbo*tss(:)**4 ! use tss as should be tss(t=tau) for MLO       ! MLO
end where                                                                                      ! MLO
dumrg(:)=-rgsave(:)                                                                            ! MLO
dumx(:)=condx(:)/dt ! total precip                                                             ! MLO
dums(:)=(conds(:)+condg(:))/dt  ! ice, snow and graupel precip                                 ! MLO
call mloeval(tss,zo,cduv,cdtq,umod,fg,eg,levspsbl,lsbl,wetfac,epot,epan,fracice,sicedep,     & ! MLO
             snowd,dt,azmin,azmin,sgdn,dumrg,dumx,dums,uav,vav,t(1:imax),qg(1:imax),         & ! MLO
             ps(1:imax),f(1:imax),swrsave,fbeamvis,fbeamnir,dumw,0,.true.,                   & ! MLO
             depth,dgice,dgscrn,dgwater,ice,water,wfull,wpack,turb)                            ! MLO
call mloextra(0,zoh,azmin,0,dgwater,dgice,ice,wpack,wfull)                                     ! MLO
call mloextra(3,zoq,azmin,0,dgwater,dgice,ice,wpack,wfull)                                     ! MLO
call mloextra(1,taux,azmin,0,dgwater,dgice,ice,wpack,wfull)                                    ! MLO
call mloextra(2,tauy,azmin,0,dgwater,dgice,ice,wpack,wfull)                                    ! MLO
                                                                                               ! MLO
where ( .not.land(1:imax) )                                                                    ! MLO
  tpan = tss(:) ! assume tss_sh=0.                                                             ! MLO
  evspsbl = evspsbl + dt*levspsbl ! levspsbl should be zero over land                          ! MLO
  sbl = sbl + dt*lsbl ! lsbl should be zero over land                                          ! MLO
end where                                                                                      ! MLO
                                                                                               ! MLO
! pan evaporation diagnostic                                                                   ! MLO
umag = max( umod, vmodmin )                                                                    ! MLO
qsttg = qsat(ps,tpan,imax)                                                                     ! MLO
ri = min(grav*zmin*(1.-tpan*srcp/t)/umag**2,ri_max)                                            ! MLO
                                                                                               ! MLO
! stuff to keep tpan over land working                                                         ! MLO
where (ri>0.)                                                                                  ! MLO
  fhd = umod/(1.+bprm*ri)**2                                                                   ! MLO
elsewhere                                                                                      ! MLO
  fhd = umod-umod*2.*bprm*ri/(1.+chs*2.*bprm*chnsea*sqrt(-ri*zmin/panzo))                      ! MLO
end where                                                                                      ! MLO
qsttg = qsat(ps,tpan,imax)                                                                     ! MLO
                                                                                               ! MLO
where ( .not.land(1:imax) )                                                                    ! MLO
  snowd = snowd*1000.                                                                          ! MLO
  ga = 0.                                                                                      ! MLO
  sno = sno + conds                                                                            ! MLO
  grpl = grpl + condg                                                                          ! MLO
  ! This cduv accounts for a moving surface                                                    ! MLO
  cduv = cduv*umod                                                                             ! MLO
  cdtq = cdtq*umod                                                                             ! MLO
  ustar = sqrt(cduv*umod)                                                                      ! MLO
elsewhere                                                                                      ! MLO
  ! assume gflux = 0                                                                           ! MLO
  ! note pan depth=.254 m, spec heat water=4186 joule/kg K                                     ! MLO
  ! and change in heat supplied=spec_heatxmassxdelta_T                                         ! MLO
  fg=rho*chnsea*cp*fhd*(tpan-theta)                                                            ! MLO
  ga=sgsave-rgsave-5.67e-8*tpan**4-panfg*fg                                                    ! MLO
  tpan=tpan+ga*dt/(4186.*0.254*1000.)                                                          ! MLO
  epan=rho*chnsea*hl*fhd*(qsttg-qg)                                                            ! MLO
end where                                                                                      ! MLO
qsttg = qsat(ps,tss,imax)                                                                      ! MLO

return
end subroutine sflux_mlo_work
    
subroutine sflux_urban

use arrays_m                       ! Atmosphere dyamics prognostic arrays
use ateb                           ! Urban
use cc_mpi                         ! CC MPI routines
use cc_omp                         ! CC OpenMP routines
use const_phys                     ! Physical constants
use estab                          ! Liquid saturation function
use extraout_m                     ! Additional diagnostics
use morepbl_m                      ! Additional boundary layer diagnostics
use newmpar_m                      ! Grid parameters
use parm_m                         ! Model configuration
use parmgeom_m                     ! Coordinate data
use pbl_m                          ! Boundary layer arrays
use prec_m, only : evspsbl, sbl    ! Precipitation
use raddiag_m                      ! Radiation diagnostic
use soilsnow_m                     ! Soil, snow and surface data
use vecsuv_m                       ! Map to cartesian coordinates
use work2_m                        ! Diagnostic arrays
use xyzinfo_m                      ! Grid coordinate arrays

implicit none

integer :: tile, is, ie
real, dimension(imax) :: luzon, lvmer, costh, sinth, zonx, zony, zonz

!$omp do schedule(static) private(is,ie),                                 &
!$omp private(luzon, lvmer, costh, sinth, zonx, zony, zonz) 
do tile=1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax

#ifdef scm
    luzon = uav(is:ie) ! zonal wind
    lvmer = vav(is:ie) ! meridonal wind 
#else
    zonx=real(                       -sin(rlat0*pi/180.)*y(is:ie))      
    zony=real(sin(rlat0*pi/180.)*x(is:ie)+cos(rlat0*pi/180.)*z(is:ie))  
    zonz=real(-cos(rlat0*pi/180.)*y(is:ie)                       )      
    costh= (zonx*ax(is:ie)+zony*ay(is:ie)+zonz*az(is:ie)) &
          /sqrt( max(zonx**2+zony**2+zonz**2,1.e-7) )      
    sinth=-(zonx*bx(is:ie)+zony*by(is:ie)+zonz*bz(is:ie)) &
          /sqrt( max(zonx**2+zony**2+zonz**2,1.e-7) )   
    luzon= costh*uav(is:ie)-sinth*vav(is:ie)  ! zonal wind
    lvmer= sinth*uav(is:ie)+costh*vav(is:ie)  ! meridonal wind
#endif

  call sflux_urban_work(azmin(is:ie),oldrunoff(is:ie),rho(is:ie),vmag(is:ie),                                        &
                        oldsnowmelt(is:ie),f_g(tile),                                                                &
                        f_intm(tile),f_road(tile),f_roof(tile),f_slab(tile),f_wall(tile),intm_g(:,tile),p_g(:,tile), &
                        rdhyd_g(:,tile),rfhyd_g(:,tile),rfveg_g(tile),road_g(:,tile),roof_g(:,tile),room_g(:,tile),  &
                        slab_g(:,tile),walle_g(:,tile),wallw_g(:,tile),cnveg_g(tile),intl_g(tile),                   &
                        luzon,lvmer,cdtq(is:ie),cduv(is:ie),conds(is:ie),                                            &
                        condg(is:ie),condx(is:ie),eg(is:ie),fg(is:ie),ps(is:ie),qg(is:ie,1),                         &
                        qsttg(is:ie),rgsave(is:ie),rnet(is:ie),runoff(is:ie),sgdn(is:ie),sgsave(is:ie),              &
                        snowmelt(is:ie),t(is:ie,1),taux(is:ie),tauy(is:ie),tss(is:ie),u(is:ie,1),                    &
                        ustar(is:ie),v(is:ie,1),vmod(is:ie),wetfac(is:ie),zo(is:ie),zoh(is:ie),zoq(is:ie),           &
                        evspsbl(is:ie),oldevspsbl(is:ie),sbl(is:ie),oldsbl(is:ie),anthropogenic_flux(is:ie),         &
                        urban_ts(is:ie),urban_wetfac(is:ie),urban_zom(is:ie),                                        &
                        urban_zoh(is:ie),urban_zoq(is:ie),urban_emiss(is:ie),urban_storage_flux(is:ie),              &
                        urban_elecgas_flux(is:ie),urban_heating_flux(is:ie),urban_cooling_flux(is:ie),               &
                        upack_g(:,tile),ufull_g(tile))

end do
!$omp end do nowait

return
end subroutine sflux_urban

subroutine sflux_urban_work(azmin,oldrunoff,rho,vmag,oldsnowmelt,fp,fp_intm,fp_road,fp_roof,                      &
                            fp_slab,fp_wall,intm,pd,rdhyd,rfhyd,rfveg,road,roof,room,slab,walle,wallw,cnveg,intl, &
                            uzon,vmer,cdtq,cduv,                                                                  &
                            conds,condg,condx,eg,fg,ps,qg,qsttg,rgsave,rnet,runoff,sgdn,sgsave,snowmelt,          &
                            t,taux,tauy,tss,u,ustar,v,vmod,wetfac,zo,zoh,zoq,evspsbl,oldevspsbl,sbl,oldsbl,       &
                            anthropogenic_flux,urban_ts,urban_wetfac,urban_zom,urban_zoh,urban_zoq,urban_emiss,   &
                            urban_storage_flux,urban_elecgas_flux,urban_heating_flux,urban_cooling_flux,          &
                            upack,ufull)

use ateb                           ! Urban
use cc_mpi                         ! CC MPI routines
use cc_omp                         ! CC OpenMP routines
use const_phys                     ! Physical constants
use estab                          ! Liquid saturation function
use newmpar_m                      ! Grid parameters
use parm_m                         ! Model configuration
use parmgeom_m                     ! Coordinate data

implicit none

integer, intent(in) :: ufull
real, dimension(imax), intent(in) :: uzon,vmer
real, dimension(imax) :: dumrg,dumx,dums
real, dimension(imax) :: newsnowmelt
real, dimension(imax) :: u_fg, u_eg, u_rn, u_evspsbl, u_sbl
real, dimension(imax) :: zo_work, zoh_work, zoq_work, u_sigma
real, dimension(imax) :: cduv_work, cdtq_work
real, dimension(imax), intent(in) :: azmin, oldrunoff, rho, vmag, oldsnowmelt, oldevspsbl, oldsbl
real, dimension(imax), intent(inout) :: anthropogenic_flux, urban_ts, urban_wetfac
real, dimension(imax), intent(inout) :: urban_zom, urban_zoh, urban_zoq, urban_emiss
real, dimension(imax), intent(inout) :: urban_storage_flux,urban_elecgas_flux
real, dimension(imax), intent(inout) :: urban_heating_flux,urban_cooling_flux
real, dimension(imax), intent(in) :: conds, condg, condx
real, dimension(imax), intent(in) :: ps, rgsave, sgdn, sgsave, vmod
real, dimension(imax), intent(inout) :: cdtq, cduv
real, dimension(imax), intent(inout) :: eg, fg, qsttg
real, dimension(imax), intent(inout) :: rnet, runoff, snowmelt
real, dimension(imax), intent(inout) :: taux, tauy, tss, ustar, wetfac
real, dimension(imax), intent(inout) :: zo, zoh, zoq, evspsbl, sbl
real, dimension(imax), intent(in) :: qg, t, u, v
logical, dimension(imax), intent(in) :: upack
type(facetparams), intent(in) :: fp_intm, fp_road, fp_roof, fp_slab, fp_wall
type(hydrodata), dimension(nfrac), intent(inout) :: rdhyd, rfhyd
type(vegdata), intent(inout) :: rfveg, cnveg
type(facetdata), dimension(nfrac), intent(inout) :: intm, road, roof, room, slab, walle, wallw
type(intldata), intent(in) :: intl
type(fparmdata), intent(in) :: fp
type(pdiagdata), dimension(nfrac), intent(inout) :: pd

! set-up radiative fluxes and precipitation                                                      ! urban
dumrg = -rgsave                                                                                  ! urban
dumx = condx/dt                                                                                  ! urban
dums = (conds+condg)/dt                                                                          ! urban
! default fluxes                                                                                 ! urban
u_fg = 0.                                                                                        ! urban
u_eg = 0.                                                                                        ! urban
urban_ts = 0.                                                                                    ! urban
urban_wetfac = 0.                                                                                ! urban
u_rn = 0.                                                                                        ! urban
u_evspsbl = 0.                                                                                   ! urban
u_sbl = 0.                                                                                       ! urban
! call aTEB                                                                                      ! urban
call atebcalc(u_fg,u_eg,urban_ts,urban_wetfac,u_rn,u_evspsbl,u_sbl,dt,azmin,sgdn,dumrg,     &    ! urban
              dumx,dums,rho,t(1:imax),qg(1:imax),ps(1:imax),uzon,vmer,vmodmin,              &    ! urban
              fp,fp_intm,fp_road,fp_roof,fp_slab,fp_wall,intm,pd,rdhyd,rfhyd,rfveg,road,    &    ! urban
              roof,room,slab,walle,wallw,cnveg,intl,upack,ufull,0,raw=.true.)                    ! urban
u_sigma = unpack(fp%sigmau,upack,0.)                                                             ! urban
where ( u_sigma>0. )                                                                             ! urban
  fg = (1.-u_sigma)*fg + u_sigma*u_fg                                                            ! urban
  eg = (1.-u_sigma)*eg + u_sigma*u_eg                                                            ! urban
  tss = (1.-u_sigma)*tss + u_sigma*urban_ts                                                      ! urban
  wetfac = (1.-u_sigma)*wetfac + u_sigma*urban_wetfac                                            ! urban
  ! since ateb will blend non-urban and urban runoff, it is                                      ! urban
  ! easier to remove the new runoff and add it again after the                                   ! urban
  ! urban scheme has been updated                                                                ! urban
  runoff = oldrunoff + (1.-u_sigma)*(runoff-oldrunoff) + u_sigma*u_rn                            ! urban
  evspsbl = oldevspsbl + (1.-u_sigma)*(evspsbl-oldevspsbl) + u_sigma*dt*u_evspsbl                ! urban
  sbl = oldsbl + (1.-u_sigma)*(sbl-oldsbl) + u_sigma*dt*u_sbl                                    ! urban
end where                                                                                        ! urban
! default urban roughness lengths                                                                ! urban
urban_zom = 1.e-10                                                                               ! urban
urban_zoh = 1.e-10                                                                               ! urban
urban_zoq = 1.e-10                                                                               ! urban
call atebzo(urban_zom,urban_zoh,urban_zoq,0,pd,fp,upack,ufull,raw=.true.)                        ! urban
! blend zo with the urban part                                                                   ! urban
zo_work  = sqrt((1.-u_sigma)/log(azmin/zo)**2+u_sigma/log(azmin/urban_zom)**2)                   ! urban
zoh_work = (1.-u_sigma)/(log(azmin/zo)*log(azmin/zoh))                                  &        ! urban
           +u_sigma/(log(azmin/urban_zom)*log(azmin/urban_zoh))                                  ! urban
zoh_work = zoh_work/zo_work                                                                      ! urban
zoq_work = (1.-u_sigma)/(log(azmin/zo)*log(azmin/zoq))                                  &        ! urban
           +u_sigma/(log(azmin/urban_zom)*log(azmin/urban_zoq))                                  ! urban
zoq_work = zoq_work/zo_work                                                                      ! urban
where ( u_sigma>0. )                                                                             ! urban
  zo  = azmin*exp(-1./zo_work)                                                                   ! urban
  zoh = azmin*exp(-1./zoh_work)                                                                  ! urban
  zoq = azmin*exp(-1./zoq_work)                                                                  ! urban
end where                                                                                        ! urban
! calculate ustar                                                                                ! urban
cduv_work = cduv/vmag                                                                            ! urban
cdtq_work = cdtq/vmag                                                                            ! urban
call atebcd(cduv_work,cdtq_work,0,pd,fp,upack,ufull)                                             ! urban
where ( u_sigma>0. )                                                                             ! urban
  cduv = cduv_work*vmag                                                                          ! urban
  cdtq = cdtq_work*vmag                                                                          ! urban
  ustar = sqrt(vmod*cduv)                                                                        ! urban
end where                                                                                        ! urban
! calculate snowmelt                                                                             ! urban
newsnowmelt = snowmelt - oldsnowmelt                                                             ! urban
newsnowmelt = newsnowmelt/dt                                                                     ! urban
call atebhydro(newsnowmelt,"snowmelt",0,pd,fp,upack,ufull)                                       ! urban
newsnowmelt = newsnowmelt*dt                                                                     ! urban
where ( u_sigma>0. )                                                                             ! urban
  snowmelt = oldsnowmelt + newsnowmelt                                                           ! urban
end where                                                                                        ! urban
! call urban fluxes                                                                              ! urban
anthropogenic_flux = 0.                                                                          ! urban
urban_elecgas_flux = 0.                                                                          ! urban
urban_heating_flux = 0.                                                                          ! urban
urban_cooling_flux = 0.                                                                          ! urban
urban_storage_flux = 0.                                                                          ! urban
call atebenergy(anthropogenic_flux,"anthropogenic",0,fp,pd,rdhyd,rfhyd,upack,ufull)              ! urban
call atebenergy(urban_elecgas_flux,"elecgas",0,fp,pd,rdhyd,rfhyd,upack,ufull)                    ! urban
call atebenergy(urban_heating_flux,"heating",0,fp,pd,rdhyd,rfhyd,upack,ufull)                    ! urban
call atebenergy(urban_cooling_flux,"cooling",0,fp,pd,rdhyd,rfhyd,upack,ufull)                    ! urban
call atebenergy(urban_storage_flux,"storage",0,fp,pd,rdhyd,rfhyd,upack,ufull)                    ! urban
where ( u_sigma>0. )                                                                             ! urban
  qsttg(1:imax) = qsat(ps,tss,imax)                                                              ! urban
  rnet(1:imax) = sgsave(1:imax) - rgsave(1:imax) - stefbo*tss(1:imax)**4                         ! urban
  taux(1:imax) = rho(1:imax)*cduv(1:imax)*u(1:imax)                                              ! urban
  tauy(1:imax) = rho(1:imax)*cduv(1:imax)*v(1:imax)                                              ! urban
end where                                                                                        ! urban
! calculate emissivity                                                                           ! urban
urban_emiss = 0.                                                                                 ! urban
call atebmisc(urban_emiss,"emissivity",0,fp,pd,upack,ufull)                                      ! urban

end subroutine sflux_urban_work

    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine sflux_land

use arrays_m                       ! Atmosphere dyamics prognostic arrays
use cc_mpi                         ! CC MPI routines
use const_phys                     ! Physical constants
use diag_m                         ! Diagnostic routines
use estab                          ! Liquid saturation function
use extraout_m                     ! Additional diagnostics
use gdrag_m                        ! Gravity wave drag
use morepbl_m                      ! Additional boundary layer diagnostics
use newmpar_m                      ! Grid parameters
use nsibd_m                        ! Land-surface arrays
use parm_m                         ! Model configuration
use pbl_m                          ! Boundary layer arrays
use soil_m                         ! Soil and surface data
use soilsnow_m                     ! Soil, snow and surface data
use work2_m                        ! Diagnostic arrays

implicit none

integer iq
real, dimension(ifull) :: fh,taftfhg_temp,factch
real zologx,xx,fhbg,es,afroot,fm,root,denma,denha
real zobg,zologbg,afland,aftlandg,rootbg,denhabg,thnew,thgnew,thnewa,aftland
real thgnewa,ri_tmp,fh_tmp

integer, parameter :: nblend=0  ! 0 for original non-blended, 1 for blended af
                                                                                               ! land
do iq = 1,ifull                                                                                ! land
  if ( land(iq) ) then                                                                         ! land 
    ! fh itself was only used outside this loop in sib0 (jlm)                                  ! land
    zobg=zobgin                                                                                ! land
    es = establ(tss(iq))                                                                       ! land
    qsttg(iq)= .622*es/max(ps(iq)-es,0.1)  ! prim for scrnout, bur recalc end sib3             ! land
    ! factch is sqrt(zo/zt) for land use in unstable fh                                        ! land
    factch(iq)=sqrt(7.4)                                                                       ! land
    if(snowd(iq)>0.)then                                                                       ! land
      ! reduce zo over snow;                                                                   ! land
      zobg=max(zobgin -snowd(iq)*0.00976/12., 0.00024)                                         ! land
      zologbg=log(zmin/zobg)                                                                   ! land
      ! following line is bit simpler than csiro9                                              ! land
      zo(iq)=max(zolnd(iq) -.001*snowd(iq), .01)                                               ! land
      zologx=log(zmin/zo(iq))                                                                  ! land
    else  ! land but not snow                                                                  ! land
      zo(iq)=zolnd(iq)                                                                         ! land
      zologbg=zologbgin                                                                        ! land
      zologx=zolog(iq)                                                                         ! land
    endif     ! (snowd(iq)>0.)                                                                 ! land
    if(nblend==1)then  ! blended zo for momentum                                               ! land
      ! note that Dorman & Sellers zo is already an average,                                   ! land
      ! accounting for sigmf, so may not wish to further blend zo                              ! land
      afland=(vkar/((1.-sigmf(iq))*zologbg+sigmf(iq)*zologx))**2                               ! land
    else    ! non-blended zo for momentum                                                      ! land
      afland=(vkar/zologx)**2                                                                  ! land
    endif   ! (nblend==1)                                                                      ! land
    aftland=vkar**2/( zologx * (2.+zologx) )                                                   ! land
    aft(iq)=aftland                                                                            ! land
    aftlandg=vkar**2/( zologbg * (2.+zologbg) )                                                ! land
    ! lgwd>0 enhances cduv (momentum) over orog under (stable & unst) condns                   ! land
    af(iq)=afland+helo(iq)       ! jlm special gwd4b                                           ! land
                                                                                               ! land
    ! Having settled on zo (and thus af) now do actual fh and fm calcs                         ! land
    xx=grav*zmin*(1.-tss(iq)*srcp/t(iq,1))                                                     ! land
    ri(iq)=min(xx/vmag(iq)**2 , ri_max)                                                        ! land
    if(ri(iq)>0.)then                                                                          ! land
      fm=vmod(iq)/(1.+bprm*ri(iq))**2                                                          ! land
      fh(iq)=fm                                                                                ! land
      fhbg=fh(iq)                                                                              ! land
    else                                                                                       ! land
      root=sqrt(-ri(iq)*zmin/zo(iq))  ! ignoring blending here                                 ! land
      ! First do momentum                                                                      ! land
      denma=1.+cms*2.*bprm*af(iq)*root                                                         ! land 
      fm=vmod(iq)-vmod(iq)*2.*bprm *ri(iq)/denma                                               ! land
      ! n.b. fm denotes ustar**2/(vmod(iq)*af)                                                 ! land
      ! Now heat ; allow for smaller zo via aft and factch                                     ! land
      denha=1.+chs*2.*bprm*factch(iq)*aft(iq)*root                                             ! land
      fh(iq)=vmod(iq)-vmod(iq)*2.*bprm *ri(iq)/denha                                           ! land
      rootbg=sqrt(-ri(iq)*zmin/zobg)                                                           ! land
      denhabg=1.+chs*2.*bprm*factch(iq)*aftlandg*rootbg                                        ! land
      fhbg=vmod(iq)-vmod(iq)*2.*bprm *ri(iq)/denhabg                                           ! land
    endif                                                                                      ! land
    taftfhg_temp(iq)=aftlandg*fhbg ! uses fmroot above, for sib3                               ! land
    taftfhg(iq)=aftlandg*fhbg ! value used for ntaft=3 (may need improving)                    ! land
                                                                                               ! land
    zoh(iq)=zo(iq)/(factch(iq)*factch(iq))                                                     ! land
    zoq(iq)=zoh(iq)                                                                            ! land
    ! cduv is now drag coeff *vmod                                                             ! land
    cduv(iq) =af(iq)*fm                                                                        ! land
    cdtq(iq) =aft(iq)*fh(iq)                                                                   ! land
    ustar(iq) = sqrt(vmod(iq)*cduv(iq))                                                        ! land
    ! Surface stresses taux, tauy: diagnostic only - unstaggered now                           ! land
    taux(iq)=rho(iq)*cduv(iq)*u(iq,1)                                                          ! land
    tauy(iq)=rho(iq)*cduv(iq)*v(iq,1)                                                          ! land
    if(ntest==1.and.iq==idjd.and.mydiag)then                                                   ! land
      write(6,*) 'in main land loop'                                                           ! land
      write(6,*) 'zmin,zobg,zobgin,snowd ',zmin,zobg,zobgin,snowd(iq)                          ! land
      write(6,*) 'afland,aftland,zologbg ',afland,aftland,zologbg                              ! land
      write(6,*) 'af,vmag,vmod,es ',af(iq),vmag(iq),vmod(iq),es                                ! land
      write(6,*) 'tss,theta,t1 ',tss(iq),theta(iq),t(iq,1)                                     ! land
      write(6,*) 'aft,fm,fh,rho ',aft(iq),fm,fh(iq),rho(iq)                                    ! land
      write(6,*) 'ri,vmod,cduv,fg ',ri(iq),vmod(iq),cduv(iq),fg(iq)                            ! land
    endif  ! (ntest==1.and.iq==idjd)                                                           ! land
  end if ! if (land(iq))                                                                       ! land  
enddo    ! iq=1,ifull                                                                          ! land
                                                                                               ! land
if(ntaft==0.or.ktau==1)then                                                                    ! land
  do iq=1,ifull  ! will only use land values                                                   ! land
    if(land(iq))then                                                                           ! land
      taftfh(iq)=aft(iq)*fh(iq) ! uses fmroot above                                            ! land
      taftfhg(iq)=taftfhg_temp(iq)                                                             ! land
    endif                                                                                      ! land
  enddo                                                                                        ! land
elseif(ntaft==1)then                                                                           ! land
  do iq=1,ifull         ! will only use land values                                            ! land
    if(land(iq))then                                                                           ! land
      thnew=aft(iq)*fh(iq) ! uses fmroot above                                                 ! land
      thgnew=taftfhg_temp(iq)                                                                  ! land
      if(thnew>2.*taftfh(iq).or.thnew<.5*taftfh(iq))then                                       ! land
        taftfh(iq)=.5*(thnew+taftfh(iq))                                                       ! land
      else                                                                                     ! land
        taftfh(iq)=thnew                                                                       ! land
      endif                                                                                    ! land
      if(thgnew>2.*taftfhg(iq).or.thgnew<.5*taftfhg(iq))then                                   ! land
        taftfhg(iq)=.5*(thgnew+taftfhg(iq))                                                    ! land
      else                                                                                     ! land
        taftfhg(iq)=thgnew                                                                     ! land
      endif                                                                                    ! land
    endif                                                                                      ! land
  enddo                                                                                        ! land
elseif(ntaft==2)then    ! preferred faster option                                              ! land
  do iq=1,ifull         ! will only use land values                                            ! land
    if(land(iq))then                                                                           ! land
      thnew=aft(iq)*fh(iq) ! uses fmroot above                                                 ! land
      thgnew=taftfhg_temp(iq)                                                                  ! land
      thnewa=min(thnew,max(2.*taftfh(iq),.5*(thnew+taftfh(iq))))                               ! land
      taftfh(iq)=max(thnewa,min(.5*taftfh(iq),.5*(thnew+taftfh(iq))))                          ! land
      thgnewa=min(thgnew,max(2.*taftfhg(iq),.5*(thgnew+taftfhg(iq))))                          ! land
      taftfhg(iq)=max(thgnewa,min(.5*taftfhg(iq),.5*(thgnew+taftfhg(iq))))                     ! land
    endif                                                                                      ! land
  enddo                                                                                        ! land
elseif(ntaft==3)then                                                                           ! land
  ! do vegetation calulation for taftfh                                                        ! land
  do iq=1,ifull                                                                                ! land
    if(land(iq))then                                                                           ! land
      xx=grav*zmin*(1.-tgf(iq)*srcp/t(iq,1)) ! actually otgf                                   ! land
      ri_tmp=min(xx/vmag(iq)**2 , ri_max)                                                      ! land
      if(ri_tmp>0.)then                                                                        ! land
        fh_tmp=vmod(iq)/(1.+bprm*ri_tmp)**2                                                    ! land
      else                                                                                     ! land
        root=sqrt(-ri_tmp*zmin/zo(iq))  ! ignoring blending                                    ! land
        ! Now heat ; allow for smaller zo via aft and factch                                   ! land
        denha=1.+chs*2.*bprm*factch(iq)*aft(iq)*root                                           ! land
        fh_tmp=vmod(iq)-vmod(iq)*2.*bprm *ri_tmp/denha                                         ! land
      endif                                                                                    ! land
      taftfh(iq)=aft(iq)*fh_tmp ! uses fmroot above, for sib3                                  ! land
    endif                                                                                      ! land
  enddo                                                                                        ! land
endif  ! (ntaft==0.or.ktau==1)  .. else ..                                                     ! land
if(ntest>0.and.mydiag)then                                                                     ! land
  write(6,*) 'before sib3 zo,zolnd,af ',zo(idjd),zolnd(idjd),af(idjd)                          ! land
  write(6,*) 'av_vmod,u,v',av_vmod,u(idjd,1),v(idjd,1)                                         ! land
endif                                                                                          ! land
                                                                                               ! land
call sib3                 ! for nsib=3, 5                                                      ! land
                                                                                               ! land
if(diag.or.ntest>0)then                                                                        ! land
  if (mydiag) write(6,*) 'before call scrnout'                                                 ! land
  call maxmin(t,' t',ktau,1.,kl)                                                               ! land
endif                                                                                          ! land
! MJT notes - This clobbers the af, ri, cduv, ustar, taux and                                  ! land
! tauy from the sice calculation above.                                                        ! land
if(ntsur/=5)then    ! ntsur=6 is default from Mar '05                                          ! land
  ! preferred option to recalc cduv, ustar (gives better uscrn, u10)                           ! land
  do iq=1,ifull                                                                                ! land
    afroot=vkar/log(zmin/zo(iq))! land formula is bit different above                          ! land
    af(iq)=afroot**2+helo(iq)                                                                  ! land
    xx=grav*zmin*(1.-tss(iq)*srcp/t(iq,1))                                                     ! land
    ri(iq)=min(xx/vmag(iq)**2 , ri_max)                                                        ! land
    if(ri(iq)>0.)then                                                                          ! land
      fm=vmod(iq)/(1.+bprm*ri(iq))**2         ! Fm * vmod                                      ! land
    else                                                                                       ! land
      root=sqrt(-ri(iq)*zmin/zo(iq))                                                           ! land
      denma=1.+cms*2.*bprm*af(iq)*root                                                         ! land
      fm=vmod(iq)-vmod(iq)*2.*bprm *ri(iq)/denma   ! Fm * vmod                                 ! land
      ! n.b. fm denotes ustar**2/(vmod(iq)*af)                                                 ! land
    endif                                                                                      ! land
    ! cduv is now drag coeff *vmod                                                             ! land
    cduv(iq) =af(iq)*fm                       ! Cd * vmod                                      ! land
    ustar(iq) = sqrt(vmod(iq)*cduv(iq))                                                        ! land
    ! Surface stresses taux, tauy: diagnostic only                                             ! land   
    taux(iq)=rho(iq)*cduv(iq)*u(iq,1)                                                          ! land
    tauy(iq)=rho(iq)*cduv(iq)*v(iq,1)                                                          ! land
  enddo                                                                                        ! land
endif  ! (ntsur==6)                                                                            ! land
                                                                                               ! land
return
end subroutine sflux_land
    
subroutine sib3

! This is the standard land-surface scheme
      
use arrays_m                     ! Atmosphere dyamics prognostic arrays
use cc_mpi                       ! CC MPI routines
use const_phys                   ! Physical constants
use dates_m                      ! Date data
use estab                        ! Liquid saturation function
use extraout_m                   ! Additional diagnostics
use latlong_m                    ! Lat/lon coordinates
use liqwpar_m                    ! Cloud water mixing ratios
use morepbl_m                    ! Additional boundary layer diagnostics
use newmpar_m                    ! Grid parameters
use nsibd_m                      ! Land-surface arrays
use parm_m                       ! Model configuration
use pbl_m                        ! Boundary layer arrays
use permsurf_m                   ! Fixed surface arrays
use prec_m, only : evspsbl, sbl  ! Precipitation
use screen_m                     ! Screen level diagnostics
use sigs_m                       ! Atmosphere sigma levels
use soil_m                       ! Soil and surface data
use soilsnow_m                   ! Soil, snow and surface data
use soilv_m                      ! Soil parameters
use vegpar_m                     ! Vegetation arrays
use work2_m                      ! Diagnostic arrays
use work3_m                      ! Mk3 land-surface diagnostic arrays

implicit none
      
integer iq,k,iveg,layer,isoil,ip,icount
real xxx,tgss,esattg,tgss2,fle,frac,conw_fh,qtgnet
real qtgair,eg1,eg2,deg,egg_alph1,sstar,ff,rsi,den
real wbav,f1,f2,f3,f4,esatf,qsatgf,beta,etr,betetrdt
real prz,dirad1,devf,ewwwa,delta_t0
real delta_t,deltat,es,tsoil
real, dimension(ifull) :: airr,cc,ccs,condxg,condsg,delta_tx,evapfb1,evapfb2,evapfb3,evapfb4
real, dimension(ifull) :: evapfb5,evapfb1a,evapfb2a,evapfb3a,evapfb4a,evapfb5a,otgf,rmcmax
real, dimension(ifull) :: tgfnew,evapfb,dqsttg,tstom,omc
real, dimension(ifull) :: ftsoil,rlai,srlai,res,tsigmf,fgf,egg
real, dimension(ifull) :: evapxf,ewww,fgg,rdg,rgg,residf,fev
real, dimension(ifull) :: extin,dirad,dfgdt,degdt,cls

integer, parameter :: ntest=0 ! ntest= 0 for diags off; ntest= 1 for diags on
!                                      2 for ewww diags      
!                     N.B. may need vsafe for correct diags
integer, parameter :: itnmeth=5
integer, parameter :: newfgf=0    ! 0 for original; 1 with tscrn; 2 with aft too

!     fle(isoil,w)=(w-swilt(isoil))/(sfc(isoil)-swilt(isoil))           !0 Eva's
!     fle(isoil,w)= w/ssat(isoil)                                       !1 simplest bare
!     fle(isoil,w)= (w-frac*max( w,swilt(isoil) ))/                     !2 jlm sugg. bare
!    .               (ssat(isoil)*(1.-frac))                            !2 jlm sugg. bare
!     fle(isoil,w)=10.*((w-ssoil(isoil))/ssat(isoil)+.1))               ! an old one
!     fle(isoil,w)= w/sfc(isoil)                                        ! jlm for PIRCS
!     fle(isoil,w)=(w-frac*swilt(isoil))/(sfc(isoil)-frac*swilt(isoil)) ! jlm special
     
do iq=1,ifull
  if(land(iq))then
    iveg=ivegt(iq)
    ! evaluate seasonal impact and the snow depth
    ! impact on the fractional vegetation cover
    tstom(iq)=298.
    if(iveg==6+31)tstom(iq)=302.
    if(iveg>=10.and.iveg<=21.and.abs(rlatt(iq)*180./pi)<25.)tstom(iq)=302.
    tsoil=min(tstom(iq), .5*(.3333*tgg(iq,2)+.6667*tgg(iq,3)+.95*tgg(iq,4) + .05*tgg(iq,5)))
    ftsoil(iq)=max(0.,1.-.0016*(tstom(iq)-tsoil)**2)
    ! which is same as:  ftsoil=max(0.,1.-.0016*(tstom-tsoil)**2)
    ! if( tsoil >= tstom ) ftsoil=1.
  endif ! (land)
enddo

if(nsib==3)then
  do iq=1,ifull
    if(land(iq))then
      iveg=ivegt(iq)
      rlai(iq)=max(.1,rlaim44(iveg)-slveg44(iveg)*(1.-ftsoil(iq)))
      srlai(iq)=rlai(iq)+rlais44(iveg)    ! nsib=3  leaf area index
      rsmin(iq) = rsunc44(iveg)/rlai(iq)  ! nsib=3  
      tsigmf(iq)=max(.001, sigmf(iq)-scveg44(iveg)*(1.-ftsoil(iq)))
      vlai(iq)=rlai(iq)
    else
      vlai(iq)=0.
    endif ! (land)
  enddo
else     ! i.e. nsib=5
  where (land)
    rlai=max(.1,vlai)
    srlai=rlai                  ! nsib=5 leaf area index
    tsigmf=max(.001,sigmf)
  end where
endif  !(nsib==3) .. else ..

if(ktau==1)then
  !if(mydiag)write(6,*) 'ipland,ipsice,ipsea in sflux: ',ipland,ipsice,ipsea
  do iq=1,ifull     ! gives default over sea too
    tgf(iq)=t(iq,1)  ! was tss(iq)
    cansto(iq)=0.
  enddo 
  do iq=1,ifull
    if(land(iq))then  ! following gives agreement on restarts
      tscrn(iq)=theta(iq)  ! first guess, needed for newfgf=1
      if(nrungcm==3)then
        do layer=2,ms
          wb(iq,layer)=wb(iq,1)   ! w, w2 and wb all same initially 
        enddo
      endif  ! (nrungcm==3)
    endif  ! (land)
  enddo
  if ( mydiag ) then
    if ( land(idjd) ) then ! MJT bugfix
      iveg=ivegt(idjd)
      isoil = isoilm(idjd)
      tsoil=0.5*(0.3333*tgg(idjd,2)+0.6667*tgg(idjd,3)+0.95*tgg(idjd,4) +  0.05*tgg(idjd,5))
      write(6,*) 'nsib,iveg,isoil,nalpha,newfgf,nsigmf,tsigmf ',nsib,iveg,isoil,nalpha,newfgf,nsigmf,tsigmf(idjd)
      write(6,*) 'ftsoil,scveg44,sigmf ',ftsoil(idjd),scveg44(iveg),sigmf(idjd)
      write(6,*) 'swilt,sfc,wb1-6 ',swilt(isoil),sfc(isoil),(wb(idjd,k),k=1,ms)
      write(6,*) 'srlai,rsmin ',srlai(idjd),rsmin(idjd)
    endif
  endif
endif           ! (ktau==1)

if(ntest==1.and.mydiag) then
  iq=idjd
  iveg=ivegt(iq)
  write(6,*) 'in sib3a iq,iveg ',iq,iveg
  write(6,*) 'snowd,zo,zolnd,tstom ',snowd(iq),zo(iq),zolnd(iq),tstom(iq)
  write(6,*) 'in sib3b iq,idjd,iveg ',iq,idjd,iveg
  write(6,*) 'iveg,sigmf(iq),tsigmfa ',iveg,sigmf(iq),tsigmf(iq)
  tsoil=0.5*(0.3333*tgg(idjd,2)+0.6667*tgg(idjd,3)+0.95*tgg(idjd,4) +  0.05*tgg(idjd,5))
  write(6,*) 'rlaim44,tsoil,ftsoil ',rlaim44(iveg),tsoil,ftsoil(iq)
  write(6,*) 'scveg44,snowd,zo,zolnd,tstom ',scveg44(iveg),snowd(iq),zo(iq),zolnd(iq),tstom(iq)
  write(6,*) 'w2,rlai ',wb(iq,ms),rlai(iq)
endif ! ntest
 
do iq=1,ifull
  if ( land(iq) ) then
    tsigmf(iq)=(1.-snowd(iq)/(snowd(iq)+5.*100.*zo(iq)))*tsigmf(iq)
    ! extin(iq)=exp(-0.6*max(1.,rlai(iq)))  ! good approx uses next 2 (jlm)
    xxx=.6*max(1.,rlai(iq))
    extin(iq)=1.-xxx/(1. +.5*xxx +xxx*xxx/12.) 
    if(ntest==1.and.iq==idjd.and.mydiag) then
      write(6,*) 'in sib3c ip,iq,idjd,iveg ',ip,iq,idjd,ivegt(iq)
      write(6,*) 'iveg,sigmf(iq),tsigmf ',ivegt(iq),sigmf(iq),tsigmf(iq)
      write(6,*) 'scveg44,snowd,zo,zolnd,tstom ',scveg44(iveg),snowd(iq),zo(iq),zolnd(iq),tstom(iq)
      write(6,*) 'alb,sgsave ',albvisnir(iq,1),sgsave(iq)
      write(6,*) 'w2,rlai,extin ',wb(iq,ms),rlai(iq),extin(iq)
    endif ! ntest
    ! bare ground calculation
    tgss=isflag(iq)*tggsn(iq,1) + (1-isflag(iq))*tgg(iq,1)
    esattg=establ(tgss)
    qsttg(iq)=.622*esattg/max(ps(iq)-esattg,0.1)
    tgss2=tgss*tgss
    dqsttg(iq)=qsttg(iq)*ps(iq)*hlars/((ps(iq)-esattg)*tgss2)
    rgg(iq) =  stefbo*tgss2**2   ! i.e. stefbo*tgss**4
    dirad(iq)=4.*rgg(iq)/tgss
    ! sensible heat flux
    dfgdt(iq)=taftfhg(iq)*rho(iq)*cp
    fgg(iq)=dfgdt(iq)*(tgss-theta(iq))
  end if      ! if land(iq)
enddo         ! iq=1,ifull

select case(nbarewet)
  case(0) ! original Eva's, same as NCAR
    do iq=1,ifull  ! all land points in this nsib=3 loop
      if ( land(iq) ) then
        isoil = isoilm(iq)
        fle=(wb(iq,1)-swilt(isoil))/(sfc(isoil)-swilt(isoil))          
        wetfac(iq)=max( 0.,min(1.,fle) )
      end if ! if land(iq)  
    enddo   ! iq=1,ifull

  case(1) ! simplest bare
    do iq=1,ifull  ! all land points in this nsib=3 loop
      if ( land(iq) ) then
        isoil = isoilm(iq)
        fle= wb(iq,1)/ssat(isoil)                                   
        wetfac(iq)=max( 0.,min(1.,fle) )
      end if ! if land(iq)  
    enddo   ! iq=1,ifull

  case(2) ! jlm suggestion
    do iq=1,ifull  ! all land points in this nsib=3 loop
      if ( land(iq) ) then
        isoil = isoilm(iq)
        frac=max(.01,tsigmf(iq))     ! jlm special for fle
        fle= (wb(iq,1)-frac*max( wb(iq,1),swilt(isoil) ))/(ssat(isoil)*(1.-frac))                         
        wetfac(iq)=max( 0.,min(1.,fle) )
      end if ! if land(iq)  
    enddo   ! iq=1,ifull

  case(3) ! jlm suggestion
    do iq=1,ifull  ! all land points in this nsib=3 loop
      if ( land(iq) ) then
        isoil = isoilm(iq)
        frac=max(.01,tsigmf(iq))     ! jlm special for fle
        fle= (wb(iq,1)-frac*swilt(isoil) )/(sfc(isoil)-frac*swilt(isoil))                    
        wetfac(iq)=max( 0.,min(1.,fle) )
      end if ! if land(iq)  
    enddo   ! iq=1,ifull

  case(4)  ! jlm, similar to Noilhan & Planton
    do iq=1,ifull  ! all land points in this nsib=3 loop
      if ( land(iq) ) then
        isoil = isoilm(iq)
        fle=min( 1.,wb(iq,1)/sfc(isoil) )         
        wetfac(iq)=fle*fle*(3.-2.*fle)
      end if ! if land(iq)  
    enddo   ! iq=1,ifull

  case(5)  ! jlm, similar to Noilhan & Planton with swilt
    do iq=1,ifull  ! all land points in this nsib=3 loop
      if ( land(iq) ) then
        isoil = isoilm(iq)
        fle=max( 0.,min( 1.,(wb(iq,1)-swilt(isoil))/(sfc(isoil)-swilt(isoil)) ) )
        wetfac(iq)=fle*fle*(3.-2.*fle)
      end if ! if land(iq)  
    enddo   ! iq=1,ifull

  case(6) ! newer jlm
    do iq=1,ifull  ! all land points in this nsib=3 loop
      if ( land(iq) ) then
        isoil = isoilm(iq)
        fle= max( 0.,min(1.,wb(iq,1)/ssat(isoil)) )
        wetfac(iq)=fle*fle*(2.2-1.2*fle)  ! .4 for fle=.5
      end if ! if land(iq)  
    enddo   ! iq=1,ifull

  case(7) ! newest piecewise jlm (use with nalpha=1, beta)
    do iq=1,ifull  ! all land points in this nsib=3 loop
      if ( land(iq) ) then
        isoil = isoilm(iq)
        wetfac(iq)=.2*(wb(iq,1)/swilt(isoil))**2
        if(wb(iq,1)>swilt(isoil))then
          wetfac(iq)=(.2*(sfc(isoil)-wb(iq,1))+.8*(wb(iq,1)-swilt(isoil)))/(sfc(isoil)-swilt(isoil))
        endif
        if(wb(iq,1)>sfc(isoil))then
          wetfac(iq)=(.8*(ssat(isoil)-wb(iq,1))+(wb(iq,1)-sfc(isoil)))/(ssat(isoil)-sfc(isoil))
        endif
      end if ! if land(iq)  
    enddo   ! iq=1,ifull

  case(8) ! like NCAR but uses ssat
    do iq=1,ifull  ! all land points in this nsib=3 loop
      if ( land(iq) ) then
        isoil = isoilm(iq)
        fle=(wb(iq,1)-swilt(isoil))/(ssat(isoil)-swilt(isoil))        
        wetfac(iq)=max( 0.,min(1.,fle) )
      end if ! if land(iq)  
    enddo   ! iq=1,ifull

end select
      
if ( nalpha==1 ) then    ! beta scheme
  do iq=1,ifull  ! all land points in this nsib=3 loop
    if ( land(iq) ) then
      isoil = isoilm(iq)
      conw_fh = rho(iq)*taftfhg(iq)*hl    ! taftfhg(iq)=aftlandg*fhbg
      epot(iq) = conw_fh*(qsttg(iq)-qg(iq,1))
      egg(iq) = wetfac(iq)*epot(iq)
      degdt(iq) = wetfac(iq)*conw_fh*dqsttg(iq)
    end if ! if land(iq)
  end do   ! iq=1,ifull
else
  ! following is alpha scheme
  do iq=1,ifull  ! all land points in this nsib=3 loop
    if ( land(iq) ) then
      isoil = isoilm(iq)
      conw_fh = rho(iq)*taftfhg(iq)*hl    ! taftfhg(iq)=aftlandg*fhbg
      qtgnet = qsttg(iq)*wetfac(iq) - qg(iq,1)
      qtgair = qsttg(iq)*wetfac(iq) - max(qtgnet, .1*qtgnet)
      eg2 = -conw_fh*qtgair
      eg1 =  conw_fh*qsttg(iq)
      ! evaporation from the bare ground
      egg(iq) = eg1*wetfac(iq) +eg2
      epot(iq) = conw_fh*(qsttg(iq)-qg(iq,1))
      deg = wetfac(iq)*conw_fh*dqsttg(iq)
      ! following reduces degdt by factor of 10 for dew
      degdt(iq) = .55*deg + sign(.45*deg, qtgnet)
    end if ! if land(iq)
  end do   ! iq=1,ifull
end if    ! (nalpha==1) .. else ..
if((ntest==1.or.diag).and.mydiag ) then
  if (land(idjd))then ! MJT bugfix
    iq=idjd
    write(6,*) 'epot,egg,tgg1,snowd ',epot(iq),egg(iq),tgg(iq,1),snowd(iq)
    isoil = isoilm(iq)
    conw_fh=rho(iq)*taftfhg(iq)*hl    ! taftfhg(iq)=aftlandg*fhbg
    qtgnet=  qsttg(iq)*wetfac(iq) -qg(iq,1)
    qtgair=  qsttg(iq)*wetfac(iq)-max(qtgnet,.1*qtgnet)
    eg2=   -conw_fh*qtgair
    eg1=    conw_fh*qsttg(iq)
    ! evaporation from the bare ground
    egg_alph1=wetfac(iq)*epot(iq)
    write(6,*) 'then iq,isoil,conw_fh,qsttg,qtgair ',iq,isoil,conw_fh,qsttg(iq),qtgair
    write(6,*) 'eg1,eg2,wetfac ',eg1,eg2,wetfac(iq)
    write(6,*) 'epot,egg,egg_alph1 ',epot(iq),egg(iq),egg_alph1
  endif
endif  ! (ntest==1)

do iq=1,ifull  ! all land points in this nsib=3 loop
  if ( land(iq) ) then
    if(snowd(iq)>1.)then
      egg(iq)=epot(iq)
      wetfac(iq)=1.   ! added jlm 18/3/04 to fix qgscrn inconsistency
      cls(iq)=1.+hlf/hl
    else
      egg(iq)=min(egg(iq),wb(iq,1)*zse(1)*1000.*hl/dt)
      cls(iq)=1.
    endif  ! (snowd(iq)>1.)
  end if ! if land(iq)
enddo    ! iq=1,ifull

select case(nsigmf)
    
  case(0)   ! original
    do iq=1,ifull  ! all land points in this nsib=3 loop
      if ( land(iq) ) then
        ga(iq)=-slwa(iq)-rgg(iq)-fgg(iq)-cls(iq)*egg(iq)       
        dgdtg(iq)=-dirad(iq)-dfgdt(iq)-cls(iq)*degdt(iq)
      end if ! if land(iq)
    enddo    ! iq=1,ifull
   
  case(1)   ! jlm preferred
    ! spreads bare-soil flux across whole grid square      
    do iq=1,ifull  ! all land points in this nsib=3 loop
      if ( land(iq) ) then
        ga(iq)=(-slwa(iq)-rgg(iq)-fgg(iq)-cls(iq)*egg(iq))*(1.-tsigmf(iq))
        ! dgtdg is used in soilsnow
        dgdtg(iq)=-(dirad(iq)+dfgdt(iq)+cls(iq)*degdt(iq))*(1.-tsigmf(iq))
      end if ! if land(iq)
    enddo    ! iq=1,ifull
    
end select

if(ntest==1.and.mydiag)then
  iq=idjd
  write(6,*)'dgdtg,dirad,dfgdt,cls,degdt,tsigmf',dgdtg(iq),dirad(iq),dfgdt(iq),cls(iq),degdt(iq),tsigmf(iq) 
endif

! ----------------------------------------------
do iq=1,ifull  ! all land points in this nsib=3 loop
  if ( land(iq) ) then
    isoil = isoilm(iq)
    iveg=ivegt(iq)
    ! components of the stomatal resistance
    sstar=90.+sign(60.,.5-zo(iq))  ! equiv to above 2 lines
    ff= 1.1*sgsave(iq)/(rlai(iq)*sstar)
    rsi = rsmin(iq) * rlai(iq)
    f1= (1.+ff)/(ff+rsi/5000.)
    den=sfc(isoil)-swilt(isoil)                          ! sib3
    wbav=(max(0.,froot(1)*(wb(iq,1)-swilt(isoil)))+  &
          max(0.,froot(2)*(wb(iq,2)-swilt(isoil)))+  &
          max(0.,froot(3)*(wb(iq,3)-swilt(isoil)))+  &
          max(0.,froot(4)*(wb(iq,4)-swilt(isoil)))+  &
          max(0.,froot(5)*(wb(iq,5)-swilt(isoil)))   )/den
    f2=max(1. , .5/ max( wbav,1.e-7)) ! N.B. this is equiv to next 2 (jlm)
    f4=max(1.-.0016*(tstom(iq)-t(iq,1))**2 , .05) ! zero for delta_t=25
    airr(iq) = 1./taftfh(iq)
    cc(iq) =min(condx(iq),           4./(1440. *60./dt))  ! jlm speedup for 4 mm/day
    ccs(iq)=min(conds(iq)+condg(iq), 4./(1440. *60./dt))
    ! depth of the reservoir of water on the canopy
    rmcmax(iq) = max(0.5,srlai(iq)) * .1
    omc(iq) = cansto(iq)  ! value from previous timestep as starting guess
    f3=max(1.-.00025*(establ(t(iq,1))-qg(iq,1)*ps(iq)/.622),.05)
    res(iq)=max(30.,rsmin(iq)*f1*f2/(f3*f4))
    if(ntest==1.and.iq==idjd.and.mydiag)then
      write(6,*) 'rlai,srlai,wbav,den ',rlai(iq),srlai(iq),wbav,den
      write(6,*) 'f1,f2,f3,f4 ',f1,f2,f3,f4
      write(6,*) 'ff,f124,rsi,res ',ff,f1*f2/f4,rsi,res(iq)
      write(6,*) 'qg,qfg,qlg ',qg(iq,1),qfg(iq,1),qlg(iq,1)
    endif
    otgf(iq)=tgf(iq)
    tgfnew(iq)=tgf(iq)
    delta_tx(iq)=10.   ! just to supply max change of 5 deg first time
    evapfb1a(iq)=max(0.,wb(iq,1)-swilt(isoil)) *zse(1)*1000.
    evapfb2a(iq)=max(0.,wb(iq,2)-swilt(isoil)) *zse(2)*1000.
    evapfb3a(iq)=max(0.,wb(iq,3)-swilt(isoil)) *zse(3)*1000.
    evapfb4a(iq)=max(0.,wb(iq,4)-swilt(isoil)) *zse(4)*1000.
    evapfb5a(iq)=max(0.,wb(iq,5)-swilt(isoil)) *zse(5)*1000.
  end if ! if land(iq)
enddo    ! iq=1,ifull

do icount=1,itnmeth     ! jlm new iteration
  ! transpiration
  do iq=1,ifull  ! all land points in this nsib=3 loop
    if ( land(iq) ) then
      esatf = establ(tgfnew(iq))
      qsatgf=.622*esatf/(ps(iq)-esatf)
      ! wet evaporation
      ewwwa = rho(iq) *(qsatgf-qg(iq,1))/airr(iq) ! in W/m**2 /hl
      ! max available dewfall is 
      ! -(qsatgf-qg1)*dsig*1000*ps/grav in mm (mult by hl/dt for W/m**2)
      ewwwa=max(ewwwa,-abs((qsatgf-qg(iq,1))*dsig(1)*ps(iq))/(grav*dt))
      ewww(iq)=min(dt*ewwwa,omc(iq),dt*ewwwa*omc(iq)/rmcmax(iq))
      ! dew(-ve), no_dew ,        no_dew
      ! cansto is reservoir on leaf
      cansto(iq)=(omc(iq)-ewww(iq)) +cc(iq)
      ewww(iq)=ewww(iq)/dt  ! these changes on 19/1/06 jlm

      ! precipitation reaching the ground under the canopy
      ! water interception by the canopy
      condxg(iq)=max(condx(iq)          -cc(iq) +max(0.,cansto(iq)-rmcmax(iq)),0.) ! keep
      condsg(iq)=max(conds(iq)+condg(iq)-ccs(iq)+max(0.,cansto(iq)-rmcmax(iq)),0.)
      cansto(iq) = min( max(0.,cansto(iq)), rmcmax(iq))
      beta =      cansto(iq)/rmcmax(iq)
      Etr=rho(iq)*max(0.,qsatgf-qg(iq,1))/(airr(iq) +res(iq))  ! jlm
      betetrdt =(1.-beta)*Etr*dt*tsigmf(iq)   ! fixed 23/5/01
      evapfb1(iq)=min(betetrdt*froot(1),evapfb1a(iq))
      evapfb2(iq)=min(betetrdt*froot(2),evapfb2a(iq))
      evapfb3(iq)=min(betetrdt*froot(3),evapfb3a(iq))
      evapfb4(iq)=min(betetrdt*froot(4),evapfb4a(iq))
      evapfb5(iq)=min(betetrdt*froot(5),evapfb5a(iq))
      evapfb(iq)=(evapfb1(iq)+evapfb2(iq)+evapfb3(iq)+evapfb4(iq)+evapfb5(iq))/tsigmf(iq)
      evapxf(iq) = (evapfb(iq)/dt + ewww(iq))*hl  ! converting to W/m**2
      prz = rho(iq)*cp*taftfh(iq)
      if(newfgf==0)fgf(iq) = prz*(tgfnew(iq)-theta(iq))  ! original/usual
      if(newfgf==1)fgf(iq) = prz*(tgfnew(iq)-tscrn(iq))
      if(newfgf==2)fgf(iq)=rho(iq)*aft(iq)*cp*(tgfnew(iq)-tscrn(iq))
      ! limit extreme fgf to avoid undue tgf oscillations  June '04
      fgf(iq)=max(-1000.,min(fgf(iq),1000.))
      rdg(iq) =  stefbo*tgfnew(iq)**4
      residf(iq) = -slwa(iq) - rdg(iq) - fgf(iq) - evapxf(iq)
      dirad1 = 4.*rdg(iq)/300.
      ! next 2 expressions can be approximated without effects
      ! dqg=qsatgf*hlars/300.**2
      ! devf= hl*rho(iq)*dqg*( (1.-beta)/(airr(iq) + res(iq))+beta/airr(iq) ) ! re-factored by jlm
      ! according to jlm prints, devf has only small effect
      devf= (hl*hlars/300.**2)*qsatgf*(1.-beta)/res(iq)
      delta_t0=residf(iq)/(dirad1 + devf + prz)
      delta_t=sign(min(abs(delta_t0),.5*abs(delta_tx(iq))),delta_t0)
      tgfnew(iq)=tgfnew(iq)+delta_t
      delta_tx(iq)=tgfnew(iq)-otgf(iq)
    end if ! if land(iq)
  enddo    ! iq=1,ifull
  if((ntest==1.or.diag).and.mydiag)then 
    if(land(idjd))then
      iq=idjd
      write(6,*) 'ktau,icount,iq,omc,cc ',ktau,icount,iq,omc(iq),cc(iq)
      write(6,*) 'rmc,rmcmax,ewww ',cansto(iq),rmcmax(iq),ewww(iq)
      esatf = establ(tgfnew(iq))  ! value for next itn
      qsatgf=.622*esatf/(ps(iq)-esatf)
      ewwwa = rho(iq) *(qsatgf-qg(iq,1))/airr(iq)
      write(6,*) 'esatf,qsatgf,ewwwa ',esatf,qsatgf,ewwwa
      prz = rho(iq)*cp*taftfh(iq)
      beta =      cansto(iq)/rmcmax(iq)
      devf= (hl*hlars/300.**2)*qsatgf*(1.-beta)/res(iq)
      dirad1 = 4.*rdg(iq)/300.
      write(6,*) 'beta,airr,res ',beta,airr(iq),res(iq)
      write(6,*) 'dirad1,devf,prz ',dirad1,devf,prz
      write(6,*) 'theta,tscrn,slwa ',theta(iq),tscrn(iq),slwa(iq)
      write(6,*) 'taftfh,condxg ',taftfh(iq),condxg(iq)
      write(6,*) 'rdg,fgf,evapxf,evapfb ',rdg(iq),fgf(iq),evapxf(iq),evapfb(iq)
      write(6,*) 'delta_tx ',delta_tx(iq)
      write(6,*) 'otgf,tgfnew,residf ',otgf(iq),tgfnew(iq),residf(iq)
    endif  ! (land(idjd))
  endif   ! ((ntest==2.or.diag).and.mydiag)
enddo     !  icount=1,5

do iq=1,ifull  ! all land points in this nsib=3 loop
  if ( land(iq) ) then
    if(cansto(iq)<1.e-10)cansto(iq)=0.  ! to avoid underflow 24/1/06
    if(tsigmf(iq) <= .0101) then
      condxpr(iq)=condx(iq)
      condspr(iq)=conds(iq)+condg(iq)
      evapfb(iq) = 0.
      evapxf(iq) = egg(iq)
      fgf(iq)  = fgg(iq)
      rdg(iq)=rgg(iq)
      tgf(iq) = tss(iq)
    else
      tgf(iq)=tgfnew(iq)
      wb(iq,1)=wb(iq,1)-evapfb1(iq)/(zse(1)*1000.)
      wb(iq,2)=wb(iq,2)-evapfb2(iq)/(zse(2)*1000.)
      wb(iq,3)=wb(iq,3)-evapfb3(iq)/(zse(3)*1000.)
      wb(iq,4)=wb(iq,4)-evapfb4(iq)/(zse(4)*1000.)
      wb(iq,5)=wb(iq,5)-evapfb5(iq)/(zse(5)*1000.)
      condxpr(iq)=(1.-tsigmf(iq))*condx(iq)+tsigmf(iq)*condxg(iq)
      condspr(iq)=(1.-tsigmf(iq))*(conds(iq)+condg(iq))+tsigmf(iq)*condsg(iq)
      if(ntest==1.and.abs(residf(iq))>10.) then
        write(6,*) 'iq,otgf(iq),tgf,delta_tx,residf ',iq,otgf(iq),tgf(iq),delta_tx(iq),residf(iq)
      end if
    endif          ! tsigmf <= .01   ... else ...
    fev(iq)=evapfb(iq)/dt*hl*tsigmf(iq) ! passed to soilsnow to update wb
    fes(iq)=(1.-tsigmf(iq))*egg(iq)*cls(iq)  ! also passed to soilsnow
    otgsoil(iq)=isflag(iq)*tggsn(iq,1) + (1-isflag(iq))*tgg(iq,1)
  end if ! if land(iq)
enddo    ! iq=1,ifull
if((ntest==1.or.diag).and.mydiag) then
  if(land(idjd))then ! MJT bugfix
    iq=idjd
    isoil = isoilm(iq)
    iveg=ivegt(iq)
    write(6,*) 'in sib3 before soilsnowv'
    write(6,*) 'evapxf,epot,egg,fev,wetfac ',evapxf(iq),epot(iq),egg(iq),fev(iq),wetfac(iq)
    write(6,*) 'fgf,fgg,fes ',fgf(iq),fgg(iq),fes(iq)
    write(6,*) 'isoil,ssat,tsigmf,rlai ',isoil,ssat(isoil),tsigmf(iq),rlai(iq)
    write(6,*) 'tgg1,t1,theta,tscrn ',tgg(iq,1),t(iq,1),theta(iq),tscrn(iq)
    write(6,*) 'qg1,qsttg ',qg(iq,1),qsttg(iq)
    write(6,*) 'dfgdt,taftfhg,rho ',dfgdt(iq),taftfhg(iq),rho(iq)
    write(6,*) 'rmc,rmcmax(iq) ',cansto(iq),rmcmax(iq)
    if (abs(tgf(iq)-otgf(iq))>4.9) then
      write(6,"('ktau,iq,otgf,tgf,dtgf,t1,t2',i4,i6,5f8.2)") ktau,iq,otgf(iq),tgf(iq),tgf(iq)-otgf(iq),t(iq,1),t(iq,2)
    end if
  endif
endif
!-------------------------------------

call soilsnowv

do iq=1,ifull  ! all land points in this nsib=3 loop
  if ( land(iq) ) then
    if(isflag(iq)==0) then
      deltat=tgg(iq,1)-otgsoil(iq)
      fgg(iq)=fgg(iq)+deltat*dfgdt(iq)
      egg(iq)=egg(iq)+deltat*degdt(iq)
      egg(iq)=min(egg(iq),wb(iq,1)*zse(1)*1000.*hl/dt)
      rgg(iq)=rgg(iq)+deltat*dirad(iq)
    else
      deltat=tggsn(iq,1)-otgsoil(iq)
      fgg(iq)=fgg(iq)+deltat*dfgdt(iq)
      egg(iq)=egg(iq)+deltat*degdt(iq)
      rgg(iq)=rgg(iq)+deltat*dirad(iq)
    endif
    ! combined fluxes
    if(snowd(iq)>1.)then
      eg(iq)=tsigmf(iq)*evapxf(iq) + egg(iq)
      evspsbl(iq) = evspsbl(iq) + dt*(tsigmf(iq)*evapxf(iq)/hl + egg(iq)/(hl*cls(iq)))
    else
      eg(iq) = tsigmf(iq)*evapxf(iq) + (1. - tsigmf(iq))*egg(iq)
      evspsbl(iq) = evspsbl(iq) + dt*(tsigmf(iq)*evapxf(iq)/hl + (1. - tsigmf(iq))*egg(iq)/(hl*cls(iq)))
    endif
    if ( cls(iq)>1.01 ) then
      sbl(iq) = sbl(iq) + dt*egg(iq)/(hl*cls(iq))  
    end if    
    if(nsigmf==2)then
      fg(iq)=tsigmf(iq)*fgf(iq)+fgg(iq)
    else
      fg(iq)=tsigmf(iq)*fgf(iq)+(1.-tsigmf(iq))*fgg(iq)
    endif
    rnet(iq)=-slwa(iq)-(1.-tsigmf(iq))*rgg(iq)-tsigmf(iq)*rdg(iq)

    tgss=isflag(iq)*tggsn(iq,1) + (1-isflag(iq))*tgg(iq,1)  ! jlm
    if(tsigmf(iq)<= .01) then
      tss(iq) = tgss
      tgf(iq) = tgss
    else
      tss(iq)=tsigmf(iq)*tgf(iq)+(1.-tsigmf(iq))*tgss
    endif       ! tsigmf<= .01
    es = establ(tss(iq))     !  from 27/12/05
    qsttg(iq)= .622*es/max(ps(iq)-es,0.1)  ! recal for scrnout, esp. snow    

  end if ! if land(iq)
enddo    ! iq=1,ifull
      
! Calculate fraction of canopy which is wet
fwet=0.
where (land)
  fwet=cansto/rmcmax
end where

if((ntest==1.or.diag).and.mydiag) then
  if (land(idjd))then ! MJT bugfix
    iq=idjd
    write(6,*) 'even further down sib3 after soilsnowv'
    write(6,*) 'tgg ',(tgg(iq,k),k=1,ms)
    write(6,*) 'wb ',(wb(iq,k),k=1,ms)
    write(6,*) 'isflag,snowd ',isflag(iq),snowd(iq)
    write(6,*) 'evapfb,fev,ewww ',evapfb(iq),fev(iq),ewww(iq)
    write(6,*) 'tsigmf,evapxf,egg ',tsigmf(iq),evapxf(iq),egg(iq)
    write(6,*) 'deltat,degdt,wb,zse ',tgg(iq,1)-otgsoil(iq),degdt(iq),wb(iq,1),zse(1)
    write(6,*) 'eg,fg ',eg(iq),fg(iq)
  endif
endif

return
end subroutine sib3

end module sflux_m
