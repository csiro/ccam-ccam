! Conformal Cubic Atmospheric Model
    
! Copyright 2015 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
! This module is the Rotstayn 1997 cloud microphysics parameterisation

! The scheme has been modifed by MJT for max/rnd cloud overlap and to include prognostic rainfall.  There is
! also an optional prognostic cloud fraction option (see cloudmod.f90).

! ldr    = 0    Diagnosed cloud scheme (depreciated)
! ldr   /= 0    Prognostic cloud condensate (different ice fall speed options)
    
! ncloud = 0    Standard LDR cloud microphysics with water vapour, liquid cloud and ice cloud
! ncloud = 1    Use newer LDR autoconvection from Mk3.6
! ncloud = 2    Same as ncloud=1, but with prognostic rain
! ncloud = 3    Same as ncloud=2, but with prognostic graupel and snow
! ncloud = 4    Use prognostic cloud fraction based on Tiedtke from GFDL-CM3, but autoconversion from ncloud=0
! ncloud = 5    Same as ncloud=4, but convective sources are included in prognostic cloud fraction

! Currently we are developing the graupel and snow components, based on Lin et al 1983 and GFDL-AM3.
   
!                            Water vapour (qg)
!
!   Cloud water (qlg,cfrac)                      Cloud ice (qfg,cfrac)
!
!   Rain (qrg,rfrac)                             Snow (qsg,sfrac)         graupel (qgrg,gfrac)

! qg, qlg, qfg, qrg, qsg and qgrg are mixing ratios (g/g) and cfrac, rfrac, sfrac, gfrac are area cover
! fractions
    
module leoncld_mod
    
private
public leoncld

real, parameter :: maxlintime = 150. ! time-step for Lin et al 83 microphysics

contains
    
! This subroutine is the interface for the LDR cloud microphysics
subroutine leoncld
      
use aerointerface                 ! Aerosol interface
use arrays_m                      ! Atmosphere dyamics prognostic arrays
use cc_mpi, only : mydiag         ! CC MPI routines
use cfrac_m                       ! Cloud fraction
use cloudmod                      ! Prognostic cloud fraction
use diag_m                        ! Diagnostic routines
use estab                         ! Liquid saturation function
use kuocomb_m                     ! JLM convection
use latlong_m                     ! Lat/lon coordinates
use liqwpar_m                     ! Cloud water mixing ratios
use morepbl_m                     ! Additional boundary layer diagnostics
use nharrs_m                      ! Non-hydrostatic atmosphere arrays 
use prec_m                        ! Precipitation
use sigs_m                        ! Atmosphere sigma levels
use soil_m                        ! Soil and surface data
use work3f_m                      ! Grid work arrays
      
implicit none
      
include 'newmpar.h'               ! Grid parameters
include 'const_phys.h'            ! Physical constants
include 'cparams.h'               ! Cloud scheme parameters
include 'kuocom.h'                ! Convection parameters
include 'parm.h'                  ! Model configuration
      
! Local variables
integer k
integer, dimension(ifull) :: kbase,ktop                   !Bottom and top of convective cloud 

real, dimension(ifull,kl) :: prf                          !Pressure on full levels (hPa)
real, dimension(ifull,kl) :: dprf                         !Pressure thickness (hPa)
real, dimension(ifull,kl) :: rhoa                         !Air density (kg/m3)
real, dimension(ifull,kl) :: dz                           !Layer thickness (m)
real, dimension(ifull,kl) :: cdso4                        !Cloud droplet conc (#/m3)
real, dimension(ifull,kl) :: ccov                         !Cloud cover (may differ from cloud frac if vertically subgrid)
real, dimension(ifull,kl) :: cfa                          !Cloud fraction in which autoconv occurs (option in newrain.f)
real, dimension(ifull,kl) :: qca                          !Cloud water mixing ratio in cfa(:,:)    (  "    "     "     )
real, dimension(ifull,kl) :: clcon                        !Convective cloud fraction in layer 
real, dimension(ifull,kl) :: qsatg                        !Saturation mixing ratio
real, dimension(ifull,kl) :: qcl                          !Vapour mixing ratio inside convective cloud
real, dimension(ifull,kl) :: qenv                         !Vapour mixing ratio outside convective cloud
real, dimension(ifull,kl) :: tenv                         !Temperature outside convective cloud
real, dimension(ifull,kl) :: tnhs                         !Non-hydrostatic temperature adjusement
real, dimension(ifull,kl) :: tv                           !Virtual air temperature
real, dimension(ifull) :: precs                           !Amount of stratiform precipitation in timestep (mm)
real, dimension(ifull) :: preci                           !Amount of stratiform snowfall in timestep (mm)
real, dimension(ifull) :: precg                           !Amount of stratiform graupel in timestep (mm)
real, dimension(ifull) :: wcon                            !Convective cloud water content (in-cloud, prescribed)

real, dimension(ifull,kl) :: qevap, qsubl, qauto, qcoll, qprog, qaccr, qaccf
real, dimension(ifull,kl) :: fluxr, fluxi, fluxs, fluxg, fluxmelt, pqfsed
real, dimension(ifull,kl) :: pfstayice, pfstayliq, slopes, prscav
real, dimension(ifull) :: prf_temp, fl, qtot


! Non-hydrostatic terms
tnhs(1:ifull,1) = phi_nh(:,1)/bet(1)
do k = 2,kl
  ! representing non-hydrostatic term as a correction to air temperature
  tnhs(1:ifull,k) = (phi_nh(:,k)-phi_nh(:,k-1)-betm(k)*tnhs(:,k-1))/bet(k)
end do

! meterological fields
do k = 1,kl
  prf(1:ifull,k)    = 0.01*ps(1:ifull)*sig(k)    !ps is SI units
  prf_temp(1:ifull) = 100.*prf(:,k)
  dprf(1:ifull,k)   = -0.01*ps(1:ifull)*dsig(k)  !dsig is -ve
  qtot(1:ifull)     = qg(1:ifull,k)+qlg(1:ifull,k)+qrg(1:ifull,k)+qfg(1:ifull,k)+qsng(1:ifull,k)+qgrg(1:ifull,k)
  tv(1:ifull,k)     = t(1:ifull,k)*(1.+1.61*qg(1:ifull,k)-qtot(:))                                   ! virtual temperature
  rhoa(1:ifull,k)   = prf_temp(:)/(rdry*tv(1:ifull,k))                                               ! air density
  qsatg(1:ifull,k)  = qsat(prf_temp(:),t(1:ifull,k))                                                 ! saturated mixing ratio
  dz(1:ifull,k)     = 100.*dprf(1:ifull,k)/(rhoa(1:ifull,k)*grav)*(1.+tnhs(1:ifull,k)/tv(1:ifull,k)) ! level thickness in metres
enddo
 
! Calculate droplet concentration from aerosols (for non-convective faction of grid-box)
call aerodrop(1,ifull,cdso4,rhoa,outconv=.true.)

! default values
kbase(1:ifull) = 0  ! default
ktop(1:ifull)  = 0  ! default
precs(1:ifull) = 0. ! rain
preci(1:ifull) = 0. ! snow
precg(1:ifull) = 0. ! hail

!     Set up convective cloud column
call convectivecloudfrac(clcon)
where ( ktsav(1:ifull)<kl-1 )
  ktop(1:ifull)   = ktsav(:)
  kbase(1:ifull)  = kbsav(:)+1
  wcon(1:ifull)   = wlc
elsewhere
  wcon(1:ifull)   = 0.
end where

if ( nmaxpr==1 .and. mydiag ) then
  if ( ktau == 1 ) then
    write(6,*)'in leoncloud acon,bcon,Rcm ',acon,bcon,Rcm
  end if
  write(6,*) 'entering leoncld'
  write(6,"('qv  ',9f8.3/4x,9f8.3)") qg(idjd,:)
  write(6,"('qf  ',9f8.3/4x,9f8.3)") qfg(idjd,:)
  write(6,"('ql  ',9f8.3/4x,9f8.3)") qlg(idjd,:)
  write(6,"('qr  ',9f8.3/4x,9f8.3)") qrg(idjd,:)
  write(6,"('qs  ',9f8.3/4x,9f8.3)") qsng(idjd,:)
  write(6,"('qg  ',9f8.3/4x,9f8.3)") qgrg(idjd,:)
endif

! Calculate convective cloud fraction and adjust moisture variables before calling newcloud
if ( ncloud <= 4 ) then

  ! diagnose cloud fraction (ncloud<=3) or prognostic strat. cloud but diagnostic conv. cloud (ncloud==4)
  if ( nmr > 0 ) then
    ! Max/Rnd cloud overlap
    do k = 1,kl
      where ( clcon(1:ifull,k)>0. )
        !ccw=wcon(:)/rhoa(:,k)  !In-cloud l.w. mixing ratio
        qccon(1:ifull,k) = clcon(:,k)*wcon(:)/rhoa(:,k)
        qcl(1:ifull,k)   = max(qsatg(:,k),qg(1:ifull,k))  ! jlm
        qenv(1:ifull,k)  = max(1.e-8,qg(1:ifull,k)-clcon(:,k)*qcl(:,k))/(1.-clcon(:,k))
        qcl(1:ifull,k)   = (qg(1:ifull,k)-(1.-clcon(:,k))*qenv(1:ifull,k))/clcon(:,k)
        qlg(1:ifull,k)   = qlg(1:ifull,k)/(1.-clcon(:,k))
        qfg(1:ifull,k)   = qfg(1:ifull,k)/(1.-clcon(:,k))
        qrg(1:ifull,k)   = qrg(1:ifull,k)/(1.-clcon(:,k))
        qsng(1:ifull,k)  = qsng(1:ifull,k)/(1.-clcon(:,k))
        qgrg(1:ifull,k)  = qgrg(1:ifull,k)/(1.-clcon(:,k))
        rfrac(1:ifull,k) = rfrac(1:ifull,k)/(1.-clcon(:,k))
        sfrac(1:ifull,k) = sfrac(1:ifull,k)/(1.-clcon(:,k))
        gfrac(1:ifull,k) = gfrac(1:ifull,k)/(1.-clcon(:,k))
      elsewhere
        clcon(1:ifull,k) = 0.
        qccon(1:ifull,k) = 0.
        qcl(1:ifull,k)   = qg(1:ifull,k)
        qenv(1:ifull,k)  = qg(1:ifull,k)
      end where
    end do
  else
    ! usual random cloud overlap
    do k = 1,kl
      where ( clcon(1:ifull,k)>0. )
        !ccw=wcon(iq)/rhoa(iq,k)  !In-cloud l.w. mixing ratio
        qccon(1:ifull,k) = clcon(:,k)*wcon(1:ifull)/rhoa(1:ifull,k)
        qcl(1:ifull,k)   = max(qsatg(1:ifull,k),qg(1:ifull,k))  ! jlm
        qenv(1:ifull,k)  = max(1.e-8,qg(1:ifull,k)-clcon(:,k)*qcl(1:ifull,k))/(1.-clcon(:,k))
        qcl(1:ifull,k)   = (qg(1:ifull,k)-(1.-clcon(:,k))*qenv(1:ifull,k))/clcon(:,k)
        qlg(1:ifull,k)   = qlg(1:ifull,k)/(1.-clcon(:,k))
        qfg(1:ifull,k)   = qfg(1:ifull,k)/(1.-clcon(:,k))
        qrg(1:ifull,k)   = qrg(1:ifull,k)/(1.-clcon(:,k))
        qsng(1:ifull,k)  = qsng(1:ifull,k)/(1.-clcon(:,k))
        qgrg(1:ifull,k)  = qgrg(1:ifull,k)/(1.-clcon(:,k))
        rfrac(1:ifull,k) = rfrac(1:ifull,k)/(1.-clcon(:,k))
        sfrac(1:ifull,k) = sfrac(1:ifull,k)/(1.-clcon(:,k))
        gfrac(1:ifull,k) = gfrac(1:ifull,k)/(1.-clcon(:,k))        
      elsewhere
        clcon(1:ifull,k) = 0.
        qccon(1:ifull,k) = 0.
        qcl(1:ifull,k)   = qg(1:ifull,k)
        qenv(1:ifull,k)  = qg(1:ifull,k)
      end where
    enddo
  end if

else
  ! prognostic strat. and conv. cloud fraction (ncloud>=5)
  ! MJT notes - no rescaling is performed because the prognostic cloud fraction scheme
  ! also accounts for convection when ncloud=5
  clcon(1:ifull,1:kl) = 0.
  qccon(1:ifull,1:kl) = 0.
  qcl(1:ifull,1:kl)   = qg(1:ifull,1:kl)
  qenv(1:ifull,1:kl)  = qg(1:ifull,1:kl)
end if
      
tenv(1:ifull,:) = t(1:ifull,:) !Assume T is the same in and out of convective cloud

if ( nmaxpr==1 .and. mydiag ) then
  write(6,*) 'before newcloud',ktau
  write(6,"('t   ',9f8.2/4x,9f8.2)") t(idjd,:)
  write(6,"('qv  ',9f8.3/4x,9f8.3)") qg(idjd,:)
  write(6,"('qf  ',9f8.3/4x,9f8.3)") qfg(idjd,:)
  write(6,"('ql  ',9f8.3/4x,9f8.3)") qlg(idjd,:)
  write(6,"('qr  ',9f8.3/4x,9f8.3)") qrg(idjd,:)
  write(6,"('qs  ',9f8.3/4x,9f8.3)") qsng(idjd,:)
  write(6,"('qg  ',9f8.3/4x,9f8.3)") qgrg(idjd,:)
  write(6,"('qnv ',9f8.3/4x,9f8.3)") qenv(idjd,:)
  write(6,"('qsat',9f8.3/4x,9f8.3)") qsatg(idjd,:)
  write(6,"('qcl ',9f8.3/4x,9f8.3)") qcl(idjd,:)
  write(6,"('clc ',9f8.3/4x,9f8.3)") clcon(idjd,:)
  write(6,*) 'kbase,ktop ',kbase(idjd),ktop(idjd)
endif

!     Calculate cloud fraction and cloud water mixing ratios
call newcloud(dt,land,prf,rhoa,cdso4,tenv,qenv,qlg,qfg,cfrac,ccov,cfa,qca)

if ( nmaxpr==1 .and. mydiag ) then
  write(6,*) 'after newcloud',ktau
  write (6,"('tnv ',9f8.2/4x,9f8.2)") tenv(idjd,:)
  write (6,"('qv0 ',9f8.3/4x,9f8.3)") qg(idjd,:)
  write (6,"('qf  ',9f8.3/4x,9f8.3)") qfg(idjd,:)
  write (6,"('ql  ',9f8.3/4x,9f8.3)") qlg(idjd,:)
  write (6,"('qr  ',9f8.3/4x,9f8.3)") qrg(idjd,:)
  write (6,"('qs  ',9f8.3/4x,9f8.3)") qsng(idjd,:)
  write (6,"('qg  ',9f8.3/4x,9f8.3)") qgrg(idjd,:)
  write (6,"('qnv ',9f8.3/4x,9f8.3)") qenv(idjd,:) ! really new qg
endif

!     Weight output variables according to non-convective fraction of grid-box            
do k = 1,kl
  t(1:ifull,k)  = clcon(:,k)*t(1:ifull,k)+(1.-clcon(:,k))*tenv(:,k)
  qg(1:ifull,k) = clcon(:,k)*qcl(:,k)+(1.-clcon(:,k))*qenv(:,k)
  where ( k>=kbase(:) .and. k<=ktop(:) )
    cfrac(1:ifull,k) = cfrac(:,k)*(1.-clcon(:,k))
    rfrac(1:ifull,k) = rfrac(1:ifull,k)*(1.-clcon(:,k))
    sfrac(1:ifull,k) = sfrac(1:ifull,k)*(1.-clcon(:,k))
    gfrac(1:ifull,k) = gfrac(1:ifull,k)*(1.-clcon(:,k))
    ccov(1:ifull,k)  = ccov(:,k)*(1.-clcon(:,k))              
    qlg(1:ifull,k)   = qlg(1:ifull,k)*(1.-clcon(:,k))
    qfg(1:ifull,k)   = qfg(1:ifull,k)*(1.-clcon(:,k))
    qrg(1:ifull,k)   = qrg(1:ifull,k)*(1.-clcon(:,k))
    qsng(1:ifull,k)  = qsng(1:ifull,k)*(1.-clcon(:,k))
    qgrg(1:ifull,k)  = qgrg(1:ifull,k)*(1.-clcon(:,k))
    cfa(1:ifull,k)   = cfa(:,k)*(1.-clcon(:,k))
    qca(1:ifull,k)   = qca(:,k)*(1.-clcon(:,k))              
  end where
enddo

if ( nmaxpr==1 .and. mydiag ) then
  write(6,*) 'before newsnowrain',ktau
  write (6,"('t   ',9f8.2/4x,9f8.2)") t(idjd,:)
  write (6,"('qv  ',9f8.3/4x,9f8.3)") qg(idjd,:)
  write (6,"('qf  ',9f8.3/4x,9f8.3)") qfg(idjd,:)
  write (6,"('ql  ',9f8.3/4x,9f8.3)") qlg(idjd,:)
  !write (6,"('qr  ',9f8.3/4x,9f8.3)") qrg(idjd,:)
  !write (6,"('qs  ',9f8.3/4x,9f8.3)") qsng(idjd,:)
  !write (6,"('qg  ',9f8.3/4x,9f8.3)") qgrg(idjd,:)
endif
if ( diag ) then
  call maxmin(t,' t',ktau,1.,kl)
  call maxmin(qg,'qv',ktau,1.e3,kl)
  call maxmin(qfg,'qf',ktau,1.e3,kl)
  call maxmin(qlg,'ql',ktau,1.e3,kl)
  call maxmin(qrg,'qr',ktau,1.e3,kl)
  call maxmin(qsng,'qs',ktau,1.e3,kl)
  call maxmin(qgrg,'qg',ktau,1.e3,kl)
endif

! Add convective cloud water into fields for radiation
! cfrad replaced by updating cfrac Oct 2005
! Moved up 16/1/06 (and ccov,cfrac NOT UPDATED in newrain)
! done because sometimes newrain drops out all qlg, ending up with 
! zero cloud (although it will be rediagnosed as 1 next timestep)
do k = 1,kl
  fl(1:ifull)      = max(0., min(1., (t(1:ifull,k)-ticon)/(273.15-ticon) ) )
  qlrad(1:ifull,k) = qlg(1:ifull,k)+qrg(1:ifull,k)+fl(:)*qccon(:,k)
  qfrad(1:ifull,k) = qfg(1:ifull,k)+qsng(1:ifull,k)+qgrg(1:ifull,k)+(1.-fl(:))*qccon(:,k)
enddo

!     Calculate precipitation and related processes
call newsnowrain(dt,rhoa,dz,prf,cdso4,cfa,qca,t,qlg,qfg,qrg,qsng,qgrg,            &
                 precs,qg,cfrac,rfrac,sfrac,gfrac,ccov,preci,precg,qevap,qsubl,   &
                 qauto,qcoll,qprog,qaccr,qaccf,fluxr,fluxi,fluxs,fluxg,fluxmelt,  &
                 pfstayice,pfstayliq,pqfsed,slopes,prscav)

if ( nmaxpr==1 .and. mydiag ) then
  write(6,*) 'after newsnowrain',ktau
  write (6,"('t   ',9f8.2/4x,9f8.2)") t(idjd,:)
  write (6,"('qv  ',9f8.3/4x,9f8.3)") qg(idjd,:)
  write (6,"('qf  ',9f8.3/4x,9f8.3)") qfg(idjd,:)
  write (6,"('ql  ',9f8.3/4x,9f8.3)") qlg(idjd,:)
  write (6,"('qr  ',9f8.3/4x,9f8.3)") qrg(idjd,:)
  write (6,"('qs  ',9f8.3/4x,9f8.3)") qsng(idjd,:)
  write (6,"('qg  ',9f8.3/4x,9f8.3)") qgrg(idjd,:)
end if
if ( diag ) then
  call maxmin(t,' t',ktau,1.,kl)
  call maxmin(qg,'qv',ktau,1.e3,kl)
  call maxmin(qfg,'qf',ktau,1.e3,kl)
  call maxmin(qlg,'ql',ktau,1.e3,kl)
  call maxmin(qrg,'qr',ktau,1.e3,kl)
  call maxmin(qsng,'qs',ktau,1.e3,kl)
  call maxmin(qgrg,'qg',ktau,1.e3,kl)
endif

!--------------------------------------------------------------
! Store data needed by prognostic aerosol scheme
! MJT notes - invert levels for aerosol code
if ( abs(iaero)>=2 ) then
  ppfprec(:,1) = 0. !At TOA
  ppfmelt(:,1) = 0. !At TOA
  ppfsnow(:,1) = 0. !At TOA
  do k = 1,kl-1
    ppfprec(:,kl+1-k) = (fluxr(:,k+1)+fluxmelt(:,k))/dt                           !flux *entering* layer k
    ppfmelt(:,kl+1-k) = fluxmelt(:,k)/dt                                          !flux melting in layer k
    ppfsnow(:,kl+1-k) = (fluxi(:,k+1)+fluxs(:,k+1)+fluxg(:,k+1)-fluxmelt(:,k))/dt !flux *entering* layer k
  end do
  do k = 1,kl
    ppfevap(:,kl+1-k)    = qevap(:,k)*rhoa(:,k)*dz(:,k)/dt
    ppfsubl(:,kl+1-k)    = qsubl(:,k)*rhoa(:,k)*dz(:,k)/dt !flux sublimating or staying in k
    pplambs(:,kl+1-k)    = slopes(:,k)
    ppmrate(:,kl+1-k)    = (qauto(:,k)+qcoll(:,k)+qprog(:,k))/dt
    ppmaccr(:,kl+1-k)    = qaccr(:,k)/dt
    ppfstayice(:,kl+1-k) = pfstayice(:,k)
    ppfstayliq(:,kl+1-k) = pfstayliq(:,k)
    ppqfsed(:,kl+1-k)    = pqfsed(:,k)
    pprscav(:,kl+1-k)    = prscav(:,k)
  end do
end if
!--------------------------------------------------------------

! Add convective cloud water into fields for radiation
! cfrad replaced by updating cfrac Oct 2005
! Moved up 16/1/06 (and ccov,cfrac NOT UPDATED in newrain)
! done because sometimes newrain drops out all qlg, ending up with 
! zero cloud (although it will be rediagnosed as 1 next timestep)
cfrac(:,1:kl) = min(1.,ccov(:,1:kl)+clcon(:,1:kl))

!========================= Jack's diag stuff =========================
!if ( ncfrp==1 ) then  ! from here to near end; Jack's diag stuff
!  do iq = 1,icfrp
!    tautot(iq)  = 0.
!    cldmax(iq)  = 0.
!    ctoptmp(iq) = 0.
!    ctoppre(iq) = 0.
!    do k = 1,kl
!      fice(iq,k) = 0.
!    enddo
!    kcldfmax(iq) = 0.
!  enddo
!!      cfrp data
!  do k = 1,kl-1
!    do iq = 1,icfrp
!      taul(iq,k) = 0.
!      taui(iq,k) = 0.
!      Reffl = 0.
!      if ( cfrac(iq,k)>0. ) then
!        tau_sfac = 1.
!        fice(iq,k) = qfrad(iq,k)/(qfrad(iq,k)+qlrad(iq,k)) ! 16/1/06
!!            Liquid water clouds
!        if ( qlg(iq,k)>1.0e-8 ) then
!          Wliq = rhoa(iq,k)*qlg(iq,k)/(cfrac(iq,k)*(1-fice(iq,k))) !kg/m^3
!          if ( .not.land(iq) ) then !sea
!            rk = 0.8
!          else            !land
!            rk = 0.67
!          endif
!! Reffl is the effective radius at the top of the cloud (calculated following
!! Martin etal 1994, JAS 51, 1823-1842) due to the extra factor of 2 in the
!! formula for reffl. Use mid cloud value of Reff for emissivity.
!          Reffl = (3*2*Wliq/(4*rhow*pi*rk*cdso4(iq,k)))**(1./3)
!          qlpath = Wliq*dz(iq,k)
!          taul(iq,k) = tau_sfac*1.5*qlpath/(rhow*Reffl)
!        endif ! qlg
!! Ice clouds
!        if ( qfg(iq,k)>1.0e-8 ) then
!          Wice = rhoa(iq,k)*qfg(iq,k)/(cfrac(iq,k)*fice(iq,k)) !kg/m**3
!          sigmai = aice*Wice**bice !visible ext. coeff. for ice
!          taui(iq,k) = sigmai*dz(iq,k) !visible opt. depth for ice
!          taui(iq,k) = tau_sfac*taui(iq,k)
!        endif ! qfg
!      endif !cfrac
!    enddo ! iq
!  enddo ! k
!! Code to get vertically integrated value...
!! top down to get highest level with cfrac=cldmax (kcldfmax)
!  do k = kl-1,1,-1
!    do iq = 1,icfrp
!      tautot(iq) = tautot(iq)+cfrac(iq,k)*(fice(iq,k)*taui(iq,k)+(1.-fice(iq,k))*taul(iq,k))
!      if ( cfrac(iq,k)>cldmax(iq) ) kcldfmax(iq) = k
!      cldmax(iq) = max(cldmax(iq),cfrac(iq,k))
!    enddo ! iq
!  enddo ! k
!
!  do iq = 1,icfrp
!    if ( cldmax(iq)>1.e-10 ) then
!      tautot(iq) = tautot(iq)/cldmax(iq)
!
!      cfd = 0.
!      do k = kl,kcldfmax(iq),-1
!        fcf = max(0.,cfrac(iq,k)-cfd) ! cld frac. from above
!        ctoptmp(iq) = ctoptmp(iq)+fcf*t(iq,k)/cldmax(iq)
!        ctoppre(iq) = ctoppre(iq)+fcf*prf(iq,k)/cldmax(iq)
!        cfd = max(cfrac(iq,k),cfd)
!      enddo ! k=kl,kcldfmax(iq),-1
!
!    endif ! (cldmax(iq).gt.1.e-10) then
!  enddo   ! iq
!endif    ! ncfrp.eq.1
!========================= end of Jack's diag stuff ======================

condx(1:ifull)  = condx(1:ifull)+precs(1:ifull)
conds(1:ifull)  = conds(1:ifull)+preci(1:ifull)
condg(1:ifull)  = condg(1:ifull)+precg(1:ifull)
precip(1:ifull) = precip(1:ifull)+precs(1:ifull)

return
end subroutine leoncld


! from arguments
!      ttg - temperature (K)
!      qtg - water vapour mixing ratio (kg/kg) - called qenv in leoncld
!      qlg - cloud liquid water mixing ratio (kg/kg)
!      qfg - cloud ice mixing ratio (kg/kg)
!
! Output:
!
! from arguments
!      cfrac - cloudy fraction of grid box
!      ccov - cloud cover looking from above (currently = cloud fraction)
! 
!******************************************************************************

 subroutine newcloud(tdt,land,prf,rhoa,cdrop,ttg,qtg,qlg,qfg,cfrac,ccov,cfa,qca)

! This routine is part of the prognostic cloud water scheme

use cc_mpi, only : mydiag
use cloudmod
use diag_m      
use estab, only : esdiffx, qsati
use map_m
use sigs_m

implicit none

! Global parameters
include 'newmpar.h'
include 'const_phys.h' ! Input physical constants
include 'cparams.h'    ! Input cloud scheme parameters
include 'kuocom.h'     ! Input cloud scheme parameters rcrit_l & rcrit_s
include 'parm.h'

! Argument list
real, intent(in) :: tdt
real, dimension(ifull,kl), intent(in) :: prf
real, dimension(ifull,kl), intent(in) :: rhoa
real, dimension(ifull,kl), intent(in) :: cdrop
real, dimension(ifull,kl), intent(inout) :: ttg
real, dimension(ifull,kl), intent(inout) :: qtg
real, dimension(ifull+iextra,kl), intent(inout) :: qlg
real, dimension(ifull+iextra,kl), intent(inout) :: qfg
real, dimension(ifull,kl), intent(inout) :: cfrac
real, dimension(ifull,kl), intent(inout) :: ccov
real, dimension(ifull,kl), intent(inout) :: cfa
real, dimension(ifull,kl), intent(inout) :: qca
logical, dimension(ifull), intent(in) :: land

! Local work arrays and variables
real, dimension(ifull,kl) :: qsl, qsw
real, dimension(ifull,kl) :: qcg, qtot, tliq
real, dimension(ifull,kl) :: fice, qcold, rcrit
real, dimension(ifull,kl) :: qsi, qfnew
real, dimension(ifull) :: pk_v, deles_v
real, dimension(ifull) :: tk, fl

integer k, mg

real al, alf, aprpr, bprpr, cice, cm0, crate, decayfac
real deles, delq, dqsdt, es, fd, hlrvap, pk, qc
real qfdep, qi0, qs, rhoic
real qcic, qcrit, qc2, qto, wliq, r3c, r6c, eps, beta6

! Start code : ----------------------------------------------------------

if ( diag.and.mydiag ) then
  write(6,*) 'entering newcloud'
  write(6,'(a,30f10.3)') 'prf ',(prf(idjd,k),k=1,kl)
  write(6,'(a,30f10.3)') 'ttg ',(ttg(idjd,k),k=1,kl)
  write(6,*) 'qtg ',(qtg(idjd,k),k=1,kl)
  write(6,*) 'qlg ',(qlg(idjd,k),k=1,kl)
  write(6,*) 'qfg ',(qfg(idjd,k),k=1,kl)
end if

! First melt cloud ice or freeze cloud water to give correct ice fraction fice.
! Then calculate the cloud conserved variables qtot and tliq.
! Note that qcg is the total cloud water (liquid+frozen)

where ( ttg(1:ifull,:)>=tfrz )
  fice(1:ifull,:) = 0.
elsewhere ( ttg(1:ifull,:)>=tice.and.qfg(1:ifull,:)>1.e-12 )
  fice(1:ifull,:) = min(qfg(1:ifull,:)/(qfg(1:ifull,:)+qlg(1:ifull,:)),1.)
elsewhere( ttg(1:ifull,:)>=tice )
  fice(1:ifull,:) = 0.
elsewhere
  fice(1:ifull,:) = 1.
end where
qcg(1:ifull,:)   = qlg(1:ifull,:)+qfg(1:ifull,:)
qcold(1:ifull,:) = qcg(:,:)
qfnew(1:ifull,:) = fice(:,:)*qcg(:,:)
ttg(1:ifull,:)   = ttg(1:ifull,:)+hlfcp*(qfnew-qfg(1:ifull,:)) !Release L.H. of fusion
qfg(1:ifull,:)   = fice(:,:)*qcg(:,:)
qlg(1:ifull,:)   = max(0.,qcg(:,:)-qfg(1:ifull,:))

if ( diag.and.mydiag ) then
  write(6,*) 'within newcloud'
  write(6,*) 'ttg ',ttg(idjd,:)
  write(6,*) 'qcold ',qcold(idjd,:)
  write(6,*) 'qcg ',qcg(idjd,:)
  write(6,*) 'qlg ',qlg(idjd,:)
  write(6,*) 'qfg ',qfg(idjd,:)
  write(6,*) 'fice ',fice(idjd,:)
end if

! Precompute the array of critical relative humidities 
if ( nclddia == -3 ) then
  do k = 1,kl
    where ( land(1:ifull) )
      rcrit(:,k) = max( rcrit_l, (1.-16.*(1.-sig(k))**3) )
    elsewhere
      rcrit(:,k) = max( rcrit_s, (1.-16.*(1.-sig(k))**3) )
    end where
  enddo
else if ( nclddia < 0 ) then
  do k = 1,kl
    where ( land(1:ifull) )
      rcrit(:,k) = max( rcrit_l, (1.-4.*(1.-sig(k))**2) )
    elsewhere
      rcrit(:,k) = max( rcrit_s, (1.-4.*(1.-sig(k))**2) )
    end where
  enddo
else if ( nclddia == 1 ) then
  do k = 1,kl
    where ( land(1:ifull) )
      rcrit(:,k) = max( rcrit_l, sig(k)**3 )
    elsewhere
      rcrit(:,k) = max( rcrit_s, sig(k)**3 )
    end where
  enddo
else if ( nclddia == 2 ) then
  do k = 1,kl
    where ( land(1:ifull) )
      rcrit(:,k) = rcrit_l
    elsewhere
      rcrit(:,k) = rcrit_s
    end where
  enddo
else if ( nclddia == 3 ) then
  do k = 1,kl
    where ( land(1:ifull) )
      rcrit(:,k) = max( rcrit_l, sig(k)**2 )          ! .75 for R21 Mk2
    elsewhere
      rcrit(:,k) = max( rcrit_s, sig(k)**2 )          ! .85 for R21 Mk2
    end where
  enddo
else if ( nclddia == 4 ) then
  do k = 1,kl
    where ( land(1:ifull) )
      rcrit(:,k) = max( rcrit_l, sig(k) )             ! .75 for test Mk2/3
    elsewhere
      rcrit(:,k) = max( rcrit_s, sig(k) )             ! .9  for test Mk2/3
    end where
  enddo
else if ( nclddia == 5 ) then  ! default till May 08
  do k = 1,kl
    where ( land(1:ifull) )
      rcrit(:,k) = max( rcrit_l, min( .99, sig(k) ) )    ! .75 for same as T63
    elsewhere
      rcrit(:,k) = max( rcrit_s, min( .99, sig(k) ) )    ! .85 for same as T63
    end where
  enddo
else if ( nclddia == 6 ) then
  do k = 1,kl
    where ( land(1:ifull) )
      rcrit(:,k) = max( rcrit_l*(1.-.15*sig(k)), sig(k)**4 )
    elsewhere
      rcrit(:,k) = max( rcrit_s*(1.-.15*sig(k)), sig(k)**4 )
    end where
  enddo
else if ( nclddia == 7 ) then
  do k = 1,kl
    where ( land(1:ifull) )
      rcrit(:,k) = max( rcrit_l*(1.-.2*sig(k)), sig(k)**4 )
    elsewhere
      rcrit(:,k) = max( rcrit_s*(1.-.2*sig(k)), sig(k)**4 )
    end where
  enddo
else if ( nclddia > 7 ) then  ! e.g. 12    JLM
  ! MJT notes - Lopez (2002) "Implementation and validation of a new pronostic large-scale cloud
  ! and precipitation scheme for climate and data-assimilation purposes" Q J R Met Soc 128, 229-257,
  ! has a useful discussion of the dependence of RHcrit on grid spacing
  do k = 1,kl  ! typically set rcrit_l=.75,  rcrit_s=.85
    tk(1:ifull) = ds/(em(1:ifull)*208498.) ! MJT suggestion
    fl(1:ifull) = (1.+real(nclddia))*tk(1:ifull)/(1.+real(nclddia)*tk(1:ifull))
    ! for rcit_l=.75 & nclddia=12 get rcrit=(0.751, 0.769, .799, .901, .940, .972, .985) for (200, 100, 50, 10, 5, 2, 1) km
    where ( land(1:ifull) )
      rcrit(1:ifull,k) = max( 1.-fl(1:ifull)*(1.-rcrit_l), sig(k)**3 )        
    elsewhere
      rcrit(1:ifull,k) = max( 1.-fl(1:ifull)*(1.-rcrit_s), sig(k)**3 )         
    end where
  end do
end if  ! (nclddia<0)  .. else ..


if ( ncloud <= 3 ) then
  ! usual diagnostic cloud fraction
      
  ! Calculate cloudy fraction of grid box (cfrac) and gridbox-mean cloud water
  ! using the triangular PDF of Smith (1990)

  do k = 1,kl
    do mg = 1,ifull
      hlrvap=(hl+fice(mg,k)*hlf)/rvap
      qtot(mg,k)=qtg(mg,k)+qcg(mg,k)
      tliq(mg,k)=ttg(mg,k)-hlcp*qcg(mg,k)-hlfcp*qfg(mg,k)

      ! Calculate qs and gam=(L/cp)*dqsdt,  at temperature tliq
      pk=100.0*prf(mg,k)
      qsi(mg,k)=qsati(pk,tliq(mg,k))           !Ice value
      deles=esdiffx(tliq(mg,k))                ! MJT suggestion
      qsl(mg,k)=qsi(mg,k)+epsil*deles/pk       !qs over liquid
      qsw(mg,k)=fice(mg,k)*qsi(mg,k)+(1.-fice(mg,k))*qsl(mg,k) !Weighted qs at temperature Tliq
      qs=qsw(mg,k)
      dqsdt=qs*hlrvap/tliq(mg,k)**2
      !qvc(mg,k)=qs !Vapour mixing ratio in cloud

      al=1./(1.+(hlcp+fice(mg,k)*hlfcp)*dqsdt)    !Smith's notation
      qc=qtot(mg,k)-qs

      delq=(1.-rcrit(mg,k))*qs      !UKMO style (equivalent to above)
      cfrac(mg,k)=1.
      qcg(mg,k)=al*qc
      if ( qc<delq ) then
        cfrac(mg,k)=max(1.e-6 , 1.-.5*((qc-delq)/delq)**2)     ! for roundoff
        qcg(mg,k)=max(1.e-8,al*(qc-(qc-delq)**3/(6.*delq**2))) ! for roundoff
      end if
      if ( qc<=0. ) then
        cfrac(mg,k)=max(1.e-6 , .5*((qc+delq)/delq)**2)    ! for roundoff
        qcg(mg,k)=max(1.e-8, al*(qc+delq)**3/(6.*delq**2)) ! for roundoff
      end if
      if ( qc<=-delq ) then
        cfrac(mg,k)=0.
        qcg(mg,k)=0.
      end if

      ! Calculate the cloud fraction (cfa) in which ql exceeds qcrit, and
      ! the corresponding gridbox-mean cloud water mixing ratio qca. 
      ! This (qca) is the cloud-water mixing ratio inside cfa times cfa.
      ! The new variable qc2 is like qc above, but is used for integration limits
      ! only, not the integrand

      if ( cfrac(mg,k)>0. ) then
        qcic=qcg(mg,k)/cfrac(mg,k) !Mean in cloud value

        ! Following few lines are for Yangang Liu's new scheme (2004: JAS, GRL)
        ! Need to do first-order estimate of qcrit using mean in-cloud qc (qcic)

        Wliq = max( 1.e-10, 1000. * qcic * rhoa(mg,k)) !g/m3
        R6c = 4.09e-4 * ( 1.15e23*1.e-6*cdrop(mg,k) / Wliq**2 ) ** (1./6.)
        eps = 1. - 0.7 * exp(-0.003e-6*cdrop(mg,k)) !mid range
        beta6 = ((1.+3.*eps**2)*(1.+4.*eps**2)*(1.+5.*eps**2) / ((1.+eps**2)*(1.+2.*eps**2)) )**(1./6.)
        R3c = 1.e-6*R6c/beta6 !in metres
        qcrit=(4.*pi/3.)*rhow*R3c**3*cdrop(mg,k)/rhoa(mg,k) !New qcrit

        qc2=qtot(mg,k)-qs-qcrit/al
        cfa(mg,k)=1.
        qca(mg,k)=al*qc
        if ( qc2<delq ) then
          cfa(mg,k)=1.-0.5*((qc2-delq)/delq)**2
          qto=(qtot(mg,k)-delq+2.*(qs+qcrit/al))/3.
          qca(mg,k)=al*(qtot(mg,k) - qto + cfa(mg,k)*(qto-qs))
        end if
        if ( qc2<=0. ) then
          cfa(mg,k)=0.5*((qc2+delq)/delq)**2
          qca(mg,k)=cfa(mg,k)*(al/3.)*(2.*qcrit/al + qc+delq)
        end if
        if ( qc2<=-delq ) then
          cfa(mg,k)=0.
          qca(mg,k)=0.
        end if
      else
        cfa(mg,k)=0.
        qca(mg,k)=0.
      end if

    end do
  end do

  if ( diag .and. mydiag ) then
    write(6,*) 'rcrit ',rcrit(idjd,:)
    write(6,*) 'qtot ',qtot(idjd,:)
    write(6,*) 'qsi',qsi(idjd,:)
    write(6,*) 'tliq',tliq(idjd,:)
    write(6,*) 'qsl ',qsl(idjd,:)
    write(6,*) 'qsw ',qsw(idjd,:)
    write(6,*) 'cfrac ',cfrac(idjd,:)
    write(6,*) 'qc  ',  qtot(idjd,:)-qsw(idjd,:)
    write(6,*) 'qcg ',qcg(idjd,:)
    write(6,*) 'delq ', (1.-rcrit(idjd,:))*qsw(idjd,:)
  endif

  ! Assume condensation or evaporation retains ice fraction fice.
  ! Introduce a time-decay factor for cirrus (as suggested by results of Khvorostyanov & Sassen,
  ! JAS, 55, 1822-1845, 1998). Their suggested range for the time constant is 0.5 to 2 hours.
  ! The grid-box-mean values of qtg and ttg are adjusted later on (below).
  decayfac = exp ( -tdt/7200. )              ! Try 2 hrs
  !decayfac = 0.                             ! Instant adjustment (old scheme)
  where( ttg(1:ifull,:) >= Tice )
    qfg(1:ifull,:) = fice*qcg
    qlg(1:ifull,:) = qcg - qfg(1:ifull,:)
  elsewhere                                    ! Cirrus T range
    qfg(1:ifull,:) = qcold*decayfac + qcg*(1.-decayfac)
    qlg(1:ifull,:) = 0.
    qcg(1:ifull,:) = qfg(1:ifull,:)
  end where
  
else
  
  ! Tiedtke prognostic cloud fraction model
  ! MJT notes - we use ttg instead of tliq
  qtot(1:ifull,1:kl)=qtg(1:ifull,1:kl)+qcg(1:ifull,1:kl)
  tliq(1:ifull,1:kl)=ttg(1:ifull,1:kl)-hlcp*qcg(1:ifull,1:kl)-hlfcp*qfg(1:ifull,1:kl)
  do k=1,kl
    pk_v=100.*prf(1:ifull,k)
    qsi(1:ifull,k)=qsati(pk_v(1:ifull),ttg(1:ifull,k))      ! Ice value
    deles_v=esdiffx(ttg(1:ifull,k))
    qsl(1:ifull,k)=qsi(1:ifull,k)+epsil*deles_v/pk_v ! Liquid value
  end do
  qsw(:,:)=fice*qsi+(1.-fice)*qsl        ! Weighted qs at temperature Tliq
  call progcloud(cfrac,qcg,qtot,prf,rhoa,fice,qsw,ttg,rcrit)
        
  ! Use 'old' autoconversion with prognostic cloud
  cfa(:,:)=0.
  qca(:,:)=0.

  where(ttg(1:ifull,:)>=Tice)
    qfg(1:ifull,:) = fice(:,:)*qcg(:,:)
    qlg(1:ifull,:) = qcg(:,:) - qfg(1:ifull,:)
  elsewhere
    qfg(1:ifull,:) = qcg(:,:)
    qlg(1:ifull,:) = 0.
    qcg(1:ifull,:) = qfg(1:ifull,:)
  end where
  
end if ! ncloud<=3 ..else..


! Do the vapour deposition calculation in mixed-phase clouds:
! Calculate deposition on cloud ice, assuming es(T) is the weighted value of the 
! liquid and ice values.
do k=1,kl
  do mg=1,ifull
    if ( cfrac(mg,k)>0. ) then
      Tk(mg)=tliq(mg,k)+hlcp*(qlg(mg,k)+qfg(mg,k))/cfrac(mg,k) !T in liq cloud
      !fl(mg)=qlg(mg,k)/max(qfg(mg,k)+qlg(mg,k),1.e-30)
      if ( Tk(mg)<tfrz.and.qlg(mg,k)>1.e-8 ) then
        pk=100*prf(mg,k)
        qs=qsati(pk,Tk(mg))
        es=qs*pk/0.622 !ice value
        Aprpr=hl/(rKa*Tk(mg))*(hls/(rvap*Tk(mg))-1.)
        Bprpr=rvap*Tk(mg)/((Dva/pk)*es)
        deles=(1.-fice(mg,k))*esdiffx(Tk(mg))
        Cice=1.e3*exp(12.96*deles/es - 0.639) !Meyers et al 1992
        cm0=1.e-12 !Initial crystal mass
        qi0=cm0*Cice/rhoa(mg,k) !Initial ice mixing ratio
        ! Next 2 lines are for assumption of fully mixed ql and qf (also a line further down).
        qi0=max(qi0, qfg(mg,k)/cfrac(mg,k)) !Assume all qf and ql are mixed
        fd=1.     !Fraction of cloud in which deposition occurs
        ! fd=fl   !Or, use option of adjacent ql,qf
        alf=1./3.
        rhoic=700.
        Crate=7.8*((Cice/rhoa(mg,k))**2/rhoic)**(1./3.)*deles/((Aprpr+Bprpr)*es)
        qfdep=fd*cfrac(mg,k)*sqrt(((2./3.)*Crate*tdt+qi0**(2./3.))**3)
        ! Also need this line for fully-mixed option...
        qfdep = qfdep - qfg(mg,k)
        qfdep=min(qfdep,qlg(mg,k))
        qlg(mg,k)=qlg(mg,k)-qfdep
        qfg(mg,k)=qfg(mg,k)+qfdep
      end if
      fice(mg,k)=qfg(mg,k)/max(qfg(mg,k)+qlg(mg,k),1.e-30)
    end if
  end do
end do  

! Calculate new values of vapour mixing ratio and temperature
qtg(1:ifull,:) = qtot(1:ifull,:) - qcg(1:ifull,:)
ttg(1:ifull,:) = tliq(1:ifull,:) + hlcp*qcg(1:ifull,:) + hlfcp*qfg(1:ifull,:)
ccov(1:ifull,:) = cfrac(1:ifull,:) !Do this for now

! Vertically sub-grid cloud
where ( cfrac(1:ifull,2:kl-1)>1.e-2 .and. cfrac(1:ifull,3:kl)==0. .and. cfrac(1:ifull,1:kl-2)==0. )
  ccov(1:ifull,2:kl-1) = sqrt(cfrac(1:ifull,2:kl-1))
end where
     
if ( diag .and. mydiag ) then
   write(6,*) 'at end of newcloud'
   write(6,*) 'ttg ',ttg(idjd,:)
   write(6,*) 'qcg ',qcg(idjd,:)
   write(6,*) 'qlg ',qlg(idjd,:)
   write(6,*) 'qfg ',qfg(idjd,:)
   write(6,*) 'qtg ',qtg(idjd,:)
end if

return
 end subroutine newcloud

! This routine is part of the prognostic cloud scheme. It calculates rainfall
! and the evaporation of rain, and also does the frozen precipitation. It is
! called by progcld.
!
! INPUT/OUTPUT
!
! Input:
!
! see also include files physparams.h (model physical constants)
!                        cparams.h    (cloud scheme parameters)
!
! from arguments
!      tdt - leapfrog timestep (seconds)
!      rhoa - air density (kg/m**3)
!      dz - layer thicknes (m)
!      prf - pressure at full levels (in hPa. NB: not SI units)
!
! In/Out:
!
! from arguments
!      ttg - temperature (K)
!      qlg - cloud liquid water mixing ratio (kg/kg)
!      qfg - cloud ice mixing ratio (kg/kg)
!      qrg - falling rain (kg/kg)
!      qsng - falling snow (kg/kg)
!      qgrg - falling graupel (kg/kg)
!      precs - amount of stratiform precipitation in timestep (mm)
!      qtg - water vapour mixing ratio (kg/kg) - called qg in C-CAM
!      cfrac - stratiform cloud fraction
!      cfrainfall - falling rain fraction
!      cfsnowfall - falling snow fraction
!      cfgraupelfall - falling graupel fraction
!      ccov - stratiform cloud *cover* looking from above (currently = cfrac)
!
! Output:
!
! from arguments
!      preci - amount of stratiform snowfall in timestep (mm)
!      precg - amount of stratiform graupel in timestep (mm)
!      qevap - evaporation of rainfall (kg/kg)
!      qsubl - sublimation of snowfall (kg/kg)
!      qauto - autoconversion of cloud liquid water (kg/kg)
!      qcoll - collection by rain of cloud liquid water (kg/kg)
 !     qprog - existing prognostic rain
!      qaccr - accretion by snow of cloud liquid water (kg/kg)
!
!**************************************************************************

subroutine newsnowrain(tdt_in,rhoa,dz,prf,cdrop,cfa,qca,ttg,qlg,qfg,qrg,qsng,qgrg,precs,qtg,cfrac,cfrainfall,     &
                       cfsnowfall,cfgraupelfall,ccov,preci,precg,qevap,qsubl,qauto,qcoll,qprog,qaccr,qaccf,fluxr, &
                       fluxi,fluxs,fluxg,fluxmelt,pfstayice,pfstayliq,pqfsed,pslopes,prscav)

use cc_mpi, only : mydiag
use estab, only : esdiffx, qsati, pow75
use kuocomb_m
use morepbl_m

implicit none

! Global parameters
include 'newmpar.h'
include 'const_phys.h' !Input physical constants
include 'cparams.h'    !Input cloud scheme parameters
include 'kuocom.h'     !acon,bcon,Rcm,ktsav,nevapls
include 'parm.h'

! Argument list
real, intent(in) :: tdt_in
real, dimension(ifull,kl), intent(in) :: rhoa
real, dimension(ifull,kl), intent(in) :: dz
real, dimension(ifull,kl), intent(in) :: prf
real, dimension(ifull,kl), intent(in) :: cdrop
real, dimension(ifull+iextra,kl), intent(inout) :: ttg
real, dimension(ifull+iextra,kl), intent(inout) :: qlg
real, dimension(ifull+iextra,kl), intent(inout) :: qfg
real, dimension(ifull+iextra,kl), intent(inout) :: qrg
real, dimension(ifull+iextra,kl), intent(inout) :: qsng
real, dimension(ifull+iextra,kl), intent(inout) :: qgrg
real, dimension(ifull+iextra,kl), intent(inout) :: qtg
real, dimension(ifull), intent(inout) :: precs
real, dimension(ifull), intent(inout) :: preci
real, dimension(ifull), intent(inout) :: precg
real, dimension(ifull,kl), intent(in) :: cfrac
real, dimension(ifull+iextra,kl), intent(inout) :: cfrainfall
real, dimension(ifull+iextra,kl), intent(inout) :: cfsnowfall
real, dimension(ifull+iextra,kl), intent(inout) :: cfgraupelfall
real, dimension(ifull,kl), intent(in) :: ccov
real, dimension(ifull,kl), intent(out) :: qevap
real, dimension(ifull,kl), intent(out) :: qsubl
real, dimension(ifull,kl), intent(out) :: qauto
real, dimension(ifull,kl), intent(out) :: qcoll
real, dimension(ifull,kl), intent(out) :: qprog
real, dimension(ifull,kl), intent(out) :: qaccr
real, dimension(ifull,kl), intent(out) :: qaccf
real, dimension(ifull,kl), intent(out) :: pqfsed
real, dimension(ifull,kl), intent(out) :: pfstayice
real, dimension(ifull,kl), intent(out) :: pfstayliq
real, dimension(ifull,kl), intent(out) :: pslopes
real, dimension(ifull,kl), intent(out) :: prscav
real, dimension(ifull,kl), intent(in) :: cfa
real, dimension(ifull,kl), intent(in) :: qca
real, dimension(ifull,kl), intent(out) :: fluxr
real, dimension(ifull,kl), intent(out) :: fluxi
real, dimension(ifull,kl), intent(out) :: fluxs
real, dimension(ifull,kl), intent(out) :: fluxg
real, dimension(ifull,kl), intent(out) :: fluxmelt

! Local work arrays and variables
real, dimension(ifull,kl-1) :: fthruliq,foutliq,fthruice,foutice
real, dimension(ifull,kl-1) :: fthrusnow,foutsnow,fthrugraupel,foutgraupel
real, dimension(ifull,kl-1) :: rhor,gam
real, dimension(ifull,kl) :: cfrain,fluxauto,vl2
real, dimension(ifull,kl) :: vi2,rhoi
real, dimension(ifull,kl) :: qprecipitation,cfsnow,fluxprecipitation,vs2,rhos
real, dimension(ifull,kl) :: qautograupel,cfgraupel,fluxautograupel,vg2,rhog
real, dimension(ifull,kl) :: clfr,cifr,qsatg,cfmelt
real, dimension(ifull) :: slopes
real, dimension(ifull) :: fluxice,fluxsnow,fluxgraupel,fluxrain
real, dimension(ifull) :: rhoiin,rhoiout,rhorin,rhorout
real, dimension(ifull) :: rhosin,rhosout,rhogin,rhogout
real, dimension(ifull) :: cffluxin,cffluxout
real, dimension(ifull) :: clfra,cifra,csfra,cgfra
real, dimension(ifull) :: mxclfrliq,rdclfrliq,mxclfrice,rdclfrice
real, dimension(ifull) :: mxclfrsnow,rdclfrsnow,mxclfrgraupel,rdclfrgraupel
real, dimension(ifull) :: fsclr_g,fsclr_s,fsclr_i,frclr
real, dimension(ifull) :: caccr_g,caccr_s,caccr_i
real, dimension(ifull) :: caccf_g,caccf_s,caccf_i
real, dimension(ifull) :: sublflux,dqf,dql,qif,dttg,csb,bf
real, dimension(ifull) :: dqs,ql,qf,qsn,qgr,qrn,cdt
real, dimension(ifull) :: rhodz,evap,qpf,clrevap,fr
real, dimension(ifull) :: mxovr,rdovr,fcol,coll,alph
real, dimension(ifull) :: alphaf,tk,pk,es,aprpr,bprpr
real, dimension(ifull) :: curly,Csbsav
real, dimension(ifull) :: n0s,lambdadum,rica
real, dimension(ifull) :: cftmp, xwgt
real, dimension(ifull) :: rhodum_g, rhodum_s, rhodum_i, rhodum_r
real, dimension(1) :: cgfr
real, dimension(3) :: cac

real, parameter :: rnzr = 8.e6
real, parameter :: rnzs = 3.e6
real, parameter :: rnzg = 4.e6
real, parameter :: rho_r = 1.0e3 ! rain density
real, parameter :: rho_s = 0.1e3 ! snow density
real, parameter :: rho_g = 0.4e3 ! grauple density
real, parameter :: qr0_crt = 2.e-4 ! rain -> snow or graupel density threshold
real, parameter :: qi0_crt = 8.e-5 ! ice -> snow density threshold
real, parameter :: qs0_crt = 6.e-3 ! snow -> graupel density threshold
real, parameter :: c_piacr = 0.1   ! accretion rate of rain -> ice
real, parameter :: c_psaut = 1.e-3 ! autoconversion rate of ice -> snow
real, parameter :: tau_i = 5.   ! (sec) cloud ice melt
real, parameter :: tau_s = 90.  ! (sec) cloud snow melt
real, parameter :: tau_g = 180. ! (sec) cloud graupel melt

integer n, k, mg, ncount

real tdt
real apr,bpr,bl,cev,crate,dqsdt,frb,qcic,qcrit,ql1,ql2,R6c,R3c,beta6,eps
real satevap,selfcoll,Wliq,cfla,dqla,qla,qsl
real craci,cracs,csacr,cgacw,cgacr,cgacs,cgaci,csacw

craci = pi*rnzr*842.*4.694155/(4.*(pi*rnzr*rho_r)**0.95)
cracs = pi*pi*rnzr*rnzs*rho_s
csacr = pi*pi*rnzr*rnzs*rho_r
cgacw = pi*rnzg*3.323363*40.74/(4.*(pi*rnzg*rho_g)**0.875)
cgacr = pi*pi*rnzr*rnzg*rho_r
cgacs = pi*pi*rnzg*rnzs*rho_s
cgaci = 0.1*cgacw
csacw = pi*rnzs*4.8*2.54925/(4.*(pi*rnzs*rho_s)**0.8125)

do k = 1,kl
  fluxr(1:ifull,k)             = 0.
  fluxi(1:ifull,k)             = 0.
  fluxs(1:ifull,k)             = 0.
  fluxg(1:ifull,k)             = 0.
  fluxmelt(1:ifull,k)          = 0.  
  qevap(1:ifull,k)             = 0.
  qauto(1:ifull,k)             = 0.
  qprecipitation(1:ifull,k)    = 0.
  qautograupel(1:ifull,k)      = 0.
  qcoll(1:ifull,k)             = 0.
  qprog(1:ifull,k)             = qrg(1:ifull,k)
  qsubl(1:ifull,k)             = 0.
  qaccr(1:ifull,k)             = 0.
  qaccf(1:ifull,k)             = 0.
  pqfsed(1:ifull,k)            = 0.
  pfstayice(1:ifull,k)         = 0.
  prscav(1:ifull,k)            = 0.
  pfstayliq(1:ifull,k)         = 0.
  pslopes(1:ifull,k)           = 0.
end do

! Use sub timestep if required
if ( ncloud >= 3 ) then
  ncount = int(tdt_in/(maxlintime+0.01)) + 1
  tdt    = tdt_in/real(ncount)
else
  ncount = 1
  tdt = tdt_in
end if

do n = 1,ncount

  do k = 1,kl
    fluxauto(1:ifull,k)          = 0.
    fluxprecipitation(1:ifull,k) = 0.
    fluxautograupel(1:ifull,k)   = 0.
    pk(1:ifull)                  = 100.*prf(1:ifull,k)
    qsatg(1:ifull,k)             = qsati(pk(1:ifull),ttg(1:ifull,k))
    cifr(1:ifull,k)              = cfrac(1:ifull,k)*qfg(1:ifull,k)/max(qlg(1:ifull,k)+qfg(1:ifull,k),1.E-30)
    clfr(1:ifull,k)              = cfrac(1:ifull,k)*qlg(1:ifull,k)/max(qlg(1:ifull,k)+qfg(1:ifull,k),1.E-30)
    cfrain(1:ifull,k)            = 0.
    cfsnow(1:ifull,k)            = 0.
    cfgraupel(1:ifull,k)         = 0.
  end do

  !**************** Cut here to insert new auto scheme ********************            
  if ( ncloud>0 .and. ncloud<=3 ) then

    ! Using new (subgrid) autoconv scheme... 
    do k = kl-1,1,-1
      do mg = 1,ifull
        cfrain(mg,k) = 0.
        rhodz(mg) = rhoa(mg,k)*dz(mg,k)
        if ( clfr(mg,k) > 0. ) then
          ql = qlg(mg,k)
          cfla = 0.
          dqla = 0.
          if ( cfa(mg,k) > 0. ) then
            cfla = cfa(mg,k)*clfr(mg,k)/(clfr(mg,k)+cifr(mg,k))
            qla = qca(mg,k)/cfa(mg,k)
            ! Following few lines are for Yangang Liu's new scheme (2004: JAS, GRL)
            Wliq = max(1.e-10, 1000. * qla * rhoa(mg,k)) !g/m3
            R6c = 4.09e-4 * ( 1.15e23*1.e-6*cdrop(mg,k) / Wliq**2 )**(1./6.)
            eps = 1. - 0.7 * exp(-0.003e-6*cdrop(mg,k)) !mid range
            beta6 = ((1.+3.*eps**2)*(1.+4.*eps**2)*(1.+5.*eps**2) / ((1.+eps**2)*(1.+2.*eps**2)) )**(1./6.)
            R3c = 1.e-6*R6c/beta6 !in metres
            qcrit = (4.*pi/3.)*rhow*R3c**3*Cdrop(mg,k)/rhoa(mg,k) !New qcrit
            if ( qla <= qcrit ) then
              ql2 = qla
            else
              ! Following is Liu & Daum (JAS, 2004)
              Crate = 1.9e17*(0.75*rhoa(mg,k)/(pi*rhow))**2*beta6**6/cdrop(mg,k)
              ql1 = qla/sqrt(1.+2.*crate*qla**2*tdt)
              ql1 = max(ql1, qcrit) !Intermediate qlg after auto
              Frb = dz(mg,k)*rhoa(mg,k)*(qla-ql1)/tdt
              cdt(mg) = tdt*0.5*Ecol*0.24*pow75(Frb)
              selfcoll = min(ql1,ql1*cdt(mg))
              ql2 = ql1 - selfcoll
            end if
            dqla = cfla*(qla-ql2)
            ql(mg) = max( 1.e-20, qlg(mg,k)-dqla )
          end if
          dql(mg) = max( qlg(mg,k)-ql(mg), 0. )
          cfrain(mg,k) = max( cfla*dql(mg)/qlg(mg,k), 0.)
          qauto(mg,k) = qauto(mg,k) + dql(mg)
          qlg(mg,k) = qlg(mg,k) - dql(mg)
          fluxauto(mg,k) = dql(mg)*rhodz(mg)
        end if
      end do
    end do

  ! Or, using old autoconv scheme... also used by prognostic cloud scheme
  else

    do k = kl-1,1,-1
      do mg = 1,ifull
        cfrain(mg,k) = 0.0
        rhodz(mg) = rhoa(mg,k)*dz(mg,k)
        if ( clfr(mg,k) > 0. ) then
          qcrit = (4.*pi/3.)*rhow*Rcm**3*cdrop(mg,k)/rhoa(mg,k)
          qcic = qlg(mg,k)/clfr(mg,k) !In cloud value
          if ( qcic < qcrit ) then
            ql = qlg(mg,k)
          else
            Crate = Aurate*rhoa(mg,k)*(rhoa(mg,k)/(cdrop(mg,k)*rhow))**(1./3.)
            ql1 = 1./pow75(qcic**(-4./3.)+(4./3.)*Crate*tdt)
            ql1 = max( ql1, qcrit ) !Intermediate qlg after auto
            Frb = dz(mg,k)*rhoa(mg,k)*(qcic-ql1)/tdt
            cdt(mg) = tdt*0.5*Ecol*0.24*pow75(Frb) !old
            selfcoll = min( ql1, ql1*cdt(mg) )
            ql2 = ql1 - selfcoll
            ql(mg) = clfr(mg,k)*ql2
          end if
          dql(mg) = max( qlg(mg,k)-ql(mg), 0. )
          cfrain(mg,k) = max( clfr(mg,k)*dql(mg)/qlg(mg,k), 0. )
          qauto(mg,k) = qauto(mg,k) + dql(mg)
          qlg(mg,k) = qlg(mg,k) - dql(mg)
          fluxauto(mg,k) = dql(mg)*rhodz(mg)
        end if
      end do
    end do

  end if ! ( ncloud>0 .and. ncloud<=3 ) ..else..


  ! calculate rate of precipitation of frozen cloud water to snow
  if ( ncloud <= 2 ) then

    ! LDR97 includes snow (qsng) in cloud ice (qfg).
    qprecipitation(1:ifull,:)    = 0.
    qautograupel(1:ifull,:)      = 0.

  else

    do k = 1,kl

      rhodz(1:ifull) = rhoa(1:ifull,k)*dz(1:ifull,k)
      
      ! autoconversion of snow to graupel (from Lin et al 1983)
      where ( ttg(1:ifull,k)<tfrz .and. qsng(1:ifull,k)>1.e-10 )
        qf(:)  = max( qsng(1:ifull,k)-qs0_crt/rhoa(:,k), 0. )
        cdt(:) = tdt*1.e-3*exp(0.09*(ttg(1:ifull,k)-tfrz))
        dqf(:) = max( min( qf(1:ifull), qf(:)*cdt(:)/(1.+cdt(:)) ), 0. )
        qautograupel(1:ifull,k)    = qautograupel(1:ifull,k) + dqf(:)
        cfgraupel(1:ifull,k)       = max( cfsnowfall(1:ifull,k)*dqf(:)/qsng(1:ifull,k), 0. )
        cfsnowfall(1:ifull,k)      = max( cfsnowfall(1:ifull,k)*(1.-dqf(:)/qsng(1:ifull,k)), 0. )
        qsng(1:ifull,k)            = qsng(1:ifull,k) - dqf(:)
        fluxautograupel(1:ifull,k) = dqf(:)*rhodz(:)
      end where

      ! autoconversion of ice to snow (from Lin et al 1983)
      ! Threshold from WSM6 scheme, Hong et al 2004, Eq(13) : qi0_crt ~8.e-5
      where ( ttg(1:ifull,k)<tfrz .and. qfg(1:ifull,k)>1.e-10 )
        qf(:)  = max( qfg(1:ifull,k)-qi0_crt/rhoa(:,k), 0. )
        cdt(:) = tdt*c_psaut*exp(0.025*(ttg(1:ifull,k)-tfrz))
        dqf(:) = max( min( qf(1:ifull), qf(:)*cdt(:)/(1.+cdt(:)) ), 0.)
        qprecipitation(1:ifull,k)    = qprecipitation(1:ifull,k) + dqf(:)
        cfsnow(1:ifull,k)            = max( cifr(1:ifull,k)*dqf(:)/qfg(1:ifull,k), 0. )
        !cifr is updated below
        qfg(1:ifull,k)               = qfg(1:ifull,k) - dqf(:)
        fluxprecipitation(1:ifull,k) = dqf(:)*rhodz(:)
      end where
      
    end do

  end if ! ( ncloud<=2 ) ..else..


  ! update cloud liquid and frozen water fractions
  cifr(1:ifull,1:kl)    = cfrac(1:ifull,1:kl)*qfg(1:ifull,1:kl)/max( qlg(1:ifull,1:kl)+qfg(1:ifull,1:kl), 1.e-30 )
  clfr(1:ifull,1:kl)    = max( cfrac(1:ifull,1:kl)-cifr(1:ifull,1:kl), 0.)
  rhoi(1:ifull,1:kl)    = qfg(1:ifull,1:kl)*rhoa(1:ifull,1:kl)
  cfmelt(1:ifull,1:kl)  = 0.
  vi2(1:ifull,kl)       = 0.1 ! Assume no cloud at top level
  fluxice(1:ifull)      = 0.
  cifra(1:ifull)        = 0.
  mxclfrice(1:ifull)    = 0. ! max overlap ice fraction
  rdclfrice(1:ifull)    = 0. ! rnd overlap ice fraction
  rica(1:ifull)         = 0. ! backward compatibility for ncloud<=2

  ! update prognostic rain
  rhor(1:ifull,1:kl-1)   = qrg(1:ifull,:)*rhoa(1:ifull,:)
  ! max overlap autoconversion and rainfall from previous time step
  cfrain(1:ifull,1:kl-1) = max( cfrain(1:ifull,1:kl-1), cfrainfall(1:ifull,1:kl-1) ) 
  vl2(1:ifull,kl)        = 0.
  clfra(1:ifull)         = 1.e-6
  fluxrain(1:ifull)      = 0.
  mxclfrliq(1:ifull)     = 0. ! max overlap rain fraction
  rdclfrliq(1:ifull)     = 0. ! rnd overlap rain fraction

  ! Set up snow fields
  rhos(1:ifull,:)     = qsng(1:ifull,:)*rhoa(1:ifull,:)
  cfsnow(1:ifull,:)   = max( cfsnow(1:ifull,:), cfsnowfall(1:ifull,:) ) 
  vs2(1:ifull,kl)     = 0.1
  fluxsnow(1:ifull)   = 0.
  csfra(1:ifull)      = 0.
  mxclfrsnow(1:ifull) = 0. ! max overlap snow fraction
  rdclfrsnow(1:ifull) = 0. ! rnd overlap snow fraction

  ! Set up graupel fields
  rhog(1:ifull,:)        = qgrg(1:ifull,:)*rhoa(1:ifull,:)
  cfgraupel(1:ifull,:)   = max( cfgraupel(1:ifull,:), cfgraupelfall(1:ifull,:) ) 
  vg2(1:ifull,kl)        = 0.1
  fluxgraupel(1:ifull)   = 0.
  cgfra(1:ifull)         = 0.
  mxclfrgraupel(1:ifull) = 0. ! max overlap graupel fraction
  rdclfrgraupel(1:ifull) = 0. ! rnd overlap graupel fraction


  if ( diag .and. mydiag ) then
    write(6,*) 'cfrac     ',cfrac(idjd,:)
    write(6,*) 'cifr      ',cifr(idjd,:)
    write(6,*) 'clfr      ',clfr(idjd,:)
    write(6,*) 'cfrain    ',cfrain(idjd,:)
    write(6,*) 'cfsnow    ',cfsnow(idjd,:)
    write(6,*) 'cfgraupel ',cfgraupel(idjd,:)
    write(6,*) 'qlg ',qlg(idjd,:)
    write(6,*) 'qfg ',qfg(idjd,:)
    write(6,*) 'qrg ',qrg(idjd,:)
    write(6,*) 'qsng',qsng(idjd,:)
    write(6,*) 'qgrg',qgrg(idjd,:)
  end if  ! (diag.and.mydiag)

  ! Now work down through the levels...
  do k = kl-1,1,-1
  
    ! misc fields
    pk(:)    = 100.*prf(:,k)
    rhodz(:) = rhoa(:,k)*dz(:,k)  
    slopes(1:ifull) = 1.6e3*10**(-0.023*(ttg(1:ifull,k)-tfrz))
    pslopes(1:ifull,k) = pslopes(1:ifull,k) + slopes(1:ifull)
  
    ! default fall velocities
    vg2(1:ifull,k) = vg2(1:ifull,k+1)
    vs2(1:ifull,k) = vs2(1:ifull,k+1)
    vi2(1:ifull,k) = vi2(1:ifull,k+1)
    vl2(1:ifull,k) = vl2(1:ifull,k+1)

  
    if ( ncloud >= 3 ) then

      ! Misc ------------------------------------------------------------------------------
      
      ! Accretion of cloud ice by rain to produce snow or grauple (from Lin et al 1983 - praci)
      ! (Neglected in UM and ACCESS 1.3)
      qf(1:ifull)       = fluxice(1:ifull)/rhodz(1:ifull)  
      rhodum_r(1:ifull) = (fluxrain(1:ifull)+evap(1:ifull))/dz(1:ifull,k)
      where ( fluxrain(1:ifull)+evap(1:ifull)>0. .and. qf(1:ifull)>1.e-10 .and. ttg(1:ifull,k)<tfrz .and. ncloud>=3 )
        cdt(1:ifull)           = tdt*craci*rhodum_r(1:ifull)**0.95/sqrt(rhoa(1:ifull,k))
        xwgt(1:ifull)          = (rhodum_r(1:ifull)-0.995*qr0_crt)/(0.01*qr0_crt) ! MJT suggestion to switch from snow to graupel
        xwgt(1:ifull)          = max( min( xwgt(1:ifull), 1. ), 0. )
        dqf(1:ifull)           = max( min( clfra(1:ifull)*qf(1:ifull), qf(1:ifull)*cdt(1:ifull)/(1.+cdt(1:ifull)) ), 0. )
        fluxice(1:ifull)       = fluxice(:)     - rhodz(:)*dqf(:)
        fluxgraupel(1:ifull)   = fluxgraupel(:) + rhodz(:)*dqf(:)*xwgt(:)
        fluxsnow(1:ifull)      = fluxsnow(:)    + rhodz(:)*dqf(:)*(1.-xwgt(:))
        rdclfrice(1:ifull)     = max( rdclfrice(:)*(1.-dqf(:)/qf(:)), 0. )
        mxclfrice(1:ifull)     = max( mxclfrice(:)*(1.-dqf(:)/qf(:)), 0. )
        cftmp(1:ifull)         = max( mxclfrice(:) + rdclfrice(:) - mxclfrice(:)*rdclfrice(:), 0. )  
        cfgraupel(1:ifull,k)   = max( cfgraupel(1:ifull,k), max( cifra(:)-cftmp(:)*xwgt(:), 0. ) )
        cfsnow(1:ifull,k)      = max( cfsnow(1:ifull,k), max( cifra(:)-cftmp(:)*(1.-xwgt(:)), 0. ) )
        cifra(1:ifull)         = cftmp(:)
      end where  
  
  
      ! Graupel ---------------------------------------------------------------------------
      alphaf(1:ifull)   = hls*qsatg(1:ifull,k)/(rvap*ttg(1:ifull,k)**2)
      gam(1:ifull,k)    = hlscp*alphaf(:) !(L/cp)*dqsdt (HBG notation)
      sublflux(1:ifull) = 0.
      caccr_g(1:ifull)  = 0.
      caccf_g(1:ifull)  = 0.

      ! Set up the rate constant for graupel sublimation
      ! MJT notes - curly and Csbsav depend on vg2(:,k+1), so vg2(:,k) can be updated below
      Tk(1:ifull)    = ttg(1:ifull,k)
      es(1:ifull)    = qsatg(1:ifull,k)*pk(:)/epsil
      Aprpr(1:ifull) = (hls/(rKa*Tk(1:ifull)))*(hls/(rvap*Tk(1:ifull))-1.)
      Bprpr(1:ifull) = rvap*Tk(1:ifull)/((Dva/pk(1:ifull))*es(1:ifull))
      where ( nevapls==-1 .or. (nevapls==-2.and.condx(:)>0..and.k<=ktsav(:)) )
        curly(1:ifull) = 0.
      elsewhere
        ! MJT notes - follow Lin83 and UM and treat graupel as spheres with a much higher density and should
        ! be treated more like like raindrops (0.31*0.493/0.44 = 0.347)
        !curly(1:ifull) = 0.65*slopes(:)**2 + 0.493*slopes(:)*sqrt(slopes(:)*vg2(:,k+1)*rhoa(:,k)/um)
        curly(1:ifull) = 0.78*slopes(:)**2 + 0.347*slopes(:)*sqrt(slopes(:)*vg2(:,k+1)*rhoa(:,k)/um)
      end where
      ! Define the rate constant for sublimation of graupel, omitting factor rhog
      Csbsav(1:ifull) = 4.*curly(:)/(rhoa(:,k)*qsatg(1:ifull,k)*(Aprpr(:)+Bprpr(:))*pi*vg2(:,k+1)*rho_g)
    
      ! The following flag detects max/random overlap clouds
      ! that are separated by a clear layer
      where ( cfgraupel(1:ifull,k)<1.e-10 .or. nmr==0 )
        ! combine max overlap from last cloud with net random overlap
        rdclfrgraupel(1:ifull) = rdclfrgraupel(:) + mxclfrgraupel(:) - rdclfrgraupel(:)*mxclfrgraupel(:)
        mxclfrgraupel(1:ifull) = 0.
      end where
    
      fluxgraupel(:) = max( fluxgraupel(:) + fluxautograupel(:,k), 0. )
    
      ! graupel fall speed (from Lin et al 1983 - see GFDL AM3)
      where ( cfgraupel(1:ifull,k) >= 1.e-10 )
        rhodum_g(:) = max( fluxgraupel(:)/dz(:,k), 0. )
        vg2(1:ifull,k) = max( 0.1, 87.2382675*(rhodum_g(:)/cfgraupel(:,k)/5026548245.74367)**0.125 )
      elsewhere
        vg2(1:ifull,k) = vg2(1:ifull,k+1)
      end where
    
      ! Set up the parameters for the flux-divergence calculation
      alph(:)           = tdt*vg2(:,k)/dz(:,k)
      foutgraupel(:,k)  = 1. - exp(-alph(:))             !analytical
      fthrugraupel(:,k) = 1. - foutgraupel(:,k)/alph(:)  !analytical
    
      ! Melt graupel if > 2 deg C  (based on Lin et al 83)
      qgr(1:ifull)      = fluxgraupel(1:ifull)/rhodz(1:ifull)
      rhodum_g(1:ifull) = qgr(1:ifull)*rhoa(:,k)
      where ( ttg(1:ifull,k)>tfrz+2. .and. qgr(1:ifull)>1.e-10 )
        cdt(1:ifull)           = tdt*((2.*pi*2.36e-2*rnzg/hlf)*(ttg(1:ifull,k)-tfrz)/rhoa(:,k)           &
                                      -2.*pi*2.11e-5*rnzg*hl/hlf*(qsatg(:,k)-qtg(1:ifull,k)))            &
                                     *((0.78/sqrt(pi*rnzg*rho_g))*sqrt(rhodum_g(:))                      &
                                      +(0.31*(1.259e-5/2.11e-5)**(1./3.)*1.608355*sqrt(40.74/1.259e-5)   &
                                      /(pi*rnzg*rho_g)**0.6875)*rhodum_g(:)**0.6875/rhoa(:,k)**0.25)
        qif(1:ifull)           = max( min( qgr(:), qgr(:)*cdt(:)/(1.+cdt(:)) ), 0. ) !Mixing ratio of graupel
        dttg(1:ifull)          = -hlfcp*qif(:)
        ttg(1:ifull,k)         = ttg(1:ifull,k) + dttg(:)
        qsatg(1:ifull,k)       = qsatg(1:ifull,k) + gam(:,k)*dttg(:)/hlscp
        fluxmelt(1:ifull,k)    = fluxmelt(:,k)  + qif(:)*rhodz(:)
        fluxgraupel(1:ifull)   = fluxgraupel(:) - qif(:)*rhodz(:)
        rdclfrgraupel(1:ifull) = max( rdclfrgraupel(:)*(1.-qif(:)/qgr(:)), 0. )
        mxclfrgraupel(1:ifull) = max( mxclfrgraupel(:)*(1.-qif(:)/qgr(:)), 0. )
        cftmp(1:ifull)         = max( mxclfrgraupel(:) + rdclfrgraupel(:) - mxclfrgraupel(:)*rdclfrgraupel(:), 0. )
        cfmelt(1:ifull,k)      = max( cfmelt(:,k), max( cgfra(:)-cftmp(:), 0. ) )
        cgfra(1:ifull)         = cftmp(:)
      end where

      ! Melt falling graupel if > 2 deg C  (based on Lin et al 83)      
      qgr(1:ifull) = fluxgraupel(1:ifull)/rhodz(1:ifull) !Mixing ratio of graupel
      where ( ttg(1:ifull,k)>tfrz+2. .and. qgr(1:ifull)>1.e-10 )
        qif(1:ifull)           = min( qgr(:), min(tdt/tau_g,1.)*(ttg(1:ifull,k)-tfrz-2.)/hlfcp )  
        dttg(1:ifull)          = -hlfcp*qif(:)
        ttg(1:ifull,k)         = ttg(1:ifull,k) + dttg(:)
        qsatg(1:ifull,k)       = qsatg(1:ifull,k) + gam(:,k)*dttg(:)/hlscp
        fluxmelt(1:ifull,k)    = fluxmelt(:,k)  + qif(:)*rhodz(:)
        fluxgraupel(1:ifull)   = fluxgraupel(:) - qif(:)*rhodz(:)
        rdclfrgraupel(1:ifull) = max( rdclfrgraupel(:)*(1.-qif(:)/qgr(:)), 0. )
        mxclfrgraupel(1:ifull) = max( mxclfrgraupel(:)*(1.-qif(:)/qgr(:)), 0. )
        cftmp(1:ifull)         = max( mxclfrgraupel(:) + rdclfrgraupel(:) - mxclfrgraupel(:)*rdclfrgraupel(:), 0. )
        cfmelt(1:ifull,k)      = max( cfmelt(:,k), max( cgfra(:)-cftmp(:), 0. ) )
        cgfra(1:ifull)         = cftmp(:)
      end where
      
      ! Sublimation of graupel is neglected in the UM and ACCESS 1.3.
      ! (Currently treated the same as LDR97 ice sublimation)
      fsclr_g(:) = max( (1.-cifr(:,k)-clfr(:,k))*fluxgraupel(:), 0. )
      where ( fluxgraupel(:)>0. .and. qtg(1:ifull,k)<qsatg(1:ifull,k) ) ! sublime graupel
        Csb(1:ifull)         = Csbsav(:)*fluxgraupel(:)/tdt
        bf(1:ifull)          = 1. + 0.5*Csb(:)*tdt*(1.+gam(:,k))
        dqs(1:ifull)         = max( 0., tdt*(Csb(:)/bf(:))*(qsatg(1:ifull,k)-qtg(1:ifull,k)) )
        dqs(1:ifull)         = min( dqs(:), (qsatg(1:ifull,k)-qtg(1:ifull,k))/(1.+gam(:,k)) ) !Don't supersat.
        sublflux(1:ifull)    = min( dqs(:)*rhodz(:), fsclr_g(:) )
        fluxgraupel(1:ifull) = fluxgraupel(:) - sublflux(:)
        fsclr_g(1:ifull)     = fsclr_g(:) - sublflux(:)
        dqs(1:ifull)         = sublflux(:)/rhodz(:)
        qsubl(1:ifull,k)     = qsubl(:,k) + dqs(:)
        qtg(1:ifull,k)       = qtg(1:ifull,k) + dqs(:)
        dttg(1:ifull)        = -hlscp*dqs(:)
        ttg(1:ifull,k)       = ttg(1:ifull,k) + dttg(:)
        qsatg(1:ifull,k)     = qsatg(1:ifull,k) + gam(:,k)*dttg(:)/hlscp
      end where
    
      ! Save flux for the wet deposition scheme.
      pfstayice(:,k) = max( pfstayice(:,k) + fluxgraupel(:)*(1.-fthrugraupel(:,k))/tdt, 0. )

      ! Accretion of cloud liquid by falling graupel (from Lin et al 1983 - pgacw)
      ! This calculation uses the incoming fluxgraupel without subtracting sublimation
      ! (since subl occurs only outside cloud), so add sublflux back to fluxgraupel.
      ql(1:ifull)       = qlg(1:ifull,k)
      rhodum_g(1:ifull) = (fluxgraupel(:)+sublflux(:))/dz(:,k)
      where ( fluxgraupel(:)+sublflux(:)>0. .and. ql(1:ifull)>1.e-10 )
        cdt(1:ifull)         = tdt*cgacw*rhodum_g(:)**0.875/sqrt(rhoa(:,k))
        dql(1:ifull)         = max( min( cgfra(:)*ql(:), ql(:)*cdt(:)/(1.+cdt(:)) ), 0. )
        qlg(1:ifull,k)       = qlg(1:ifull,k)   - dql(:)
        qaccr(1:ifull,k)     = qaccr(1:ifull,k) + dql(:)
        fluxgraupel(1:ifull) = fluxgraupel(:) + rhodz(:)*dql(:)
        dttg(1:ifull)        = hlfcp*dql(:)
        ttg(1:ifull,k)       = ttg(1:ifull,k) + dttg(:)
        qsatg(1:ifull,k)     = qsatg(1:ifull,k) + gam(:,k)*dttg(:)/hlscp
        cftmp(1:ifull)       = max( clfr(:,k)*dql(:)/ql(1:ifull), 0. )
        clfr(1:ifull,k)      = max( clfr(:,k)*(1.-dql(:)/ql(1:ifull)), 0. )
        caccr_g(1:ifull)     = max( caccr_g(1:ifull), cftmp(1:ifull) )
      end where
    
      ! Accretion of rain by falling graupel (from Lin et al 1983 - pgacr)
      ! (Neglected in UM and ACCESS 1.3)
      cac(1) = 5./((pi*rnzr*rho_r)**1.5*(pi*rnzg*rho_g)**0.25) ! modified from acco(1,3)
      cac(2) = 2./((pi*rnzr*rho_r)**1.25*(pi*rnzg*rho_g)**0.5)
      cac(3) = 0.5/((pi*rnzr*rho_r)*(pi*rnzg*rho_g)**0.75)
      qrn(1:ifull)      = fluxrain(1:ifull)/rhodz(:)
      rhodum_r(1:ifull) = qrn(1:ifull)*rhoa(1:ifull,k)
      rhodum_g(1:ifull) = (fluxgraupel(:)+sublflux(:))/dz(:,k)
      where ( fluxgraupel(:)+sublflux(:)>0. .and. qrn(1:ifull)>1.e-10 .and. ttg(1:ifull,k)<tfrz )
        cdt(1:ifull)         = tdt*cgacr*abs(vg2(:,k)-vl2(:,k))*qrn(:)*rhodum_g(:)**0.25*           &
                               (cac(1)*sqrt(rhodum_r(:))+cac(2)*rhodum_r(:)**0.25*rhodum_g(:)**0.25 &
                               +cac(3)*sqrt(rhodum_g(:)))    
        dql(1:ifull)         = max( min( cgfra(:)*qrn(:), qrn(:)*cdt(:)/(1.+cdt(:)) ), 0. )
        qaccr(1:ifull,k)     = qaccr(:,k) + dql(:)  
        fluxrain(1:ifull)    = fluxrain(:)    - rhodz(:)*dql(:)
        fluxgraupel(1:ifull) = fluxgraupel(:) + rhodz(:)*dql(:)
        dttg(1:ifull)        = hlfcp*dql(:)
        ttg(1:ifull,k)       = ttg(1:ifull,k) + dttg(:)
        qsatg(1:ifull,k)     = qsatg(1:ifull,k) + gam(:,k)*dttg(:)/hlscp 
        rdclfrliq(1:ifull)   = max( rdclfrliq(:)*(1.-dql(:)/qrn(:)), 0. )
        mxclfrliq(1:ifull)   = max( mxclfrliq(:)*(1.-dql(:)/qrn(:)), 0. )
        cftmp(1:ifull)       = max( mxclfrliq(:) + rdclfrliq(:) - mxclfrliq(:)*rdclfrliq(:), 0. )      
        caccr_g(1:ifull)     = max( caccr_g(:), max( clfra(:)-cftmp(:), 0. ) )
        clfra(1:ifull)       = cftmp(:)
      end where   
    
      ! Accretion of cloud ice by falling graupel (from Lin et al 1983 - pgaci)
      ! (Neglected in UM and ACCESS 1.3)
      qf(1:ifull)       = fluxice(1:ifull)/rhodz(:)
      rhodum_g(1:ifull) = (fluxgraupel(:)+sublflux(:))/dz(:,k)
      where ( fluxgraupel(:)+sublflux(:)>0. .and. qf(1:ifull)>1.e-10 .and. ttg(1:ifull,k)<tfrz )
        cdt(1:ifull)         = tdt*cgaci*rhodum_g(:)**0.875/sqrt(rhoa(:,k))
        dqf(1:ifull)         = max( min( cgfra(:)*qf(:), qf(:)*cdt(:)/(1.+cdt(:)) ), 0. )
        qaccf(1:ifull,k)     = qaccf(:,k) + dqf(:)      
        fluxice(1:ifull)     = fluxice(:)     - rhodz(:)*dqf(:)
        fluxgraupel(1:ifull) = fluxgraupel(:) + rhodz(:)*dqf(:)
        rdclfrice(1:ifull)   = max( rdclfrice(:)*(1.-dqf(:)/qf(:)), 0. )
        mxclfrice(1:ifull)   = max( mxclfrice(:)*(1.-dqf(:)/qf(:)), 0. )
        cftmp(1:ifull)       = max( mxclfrice(:) + rdclfrice(:) - mxclfrice(:)*rdclfrice(:), 0. )  
        caccf_g(1:ifull)     = max( caccf_g(:), max( cifra(:)-cftmp(:), 0. ) )
        cifra(1:ifull)       = cftmp(:)
      end where

      ! Accretion of snow by falling graupel (from Lin et al 1983 - pgacs )
      cac(1) = 5./((pi*rnzs*rho_s)**1.5*(pi*rnzg*rho_g)**0.25)
      cac(2) = 2./((pi*rnzs*rho_s)**1.25*(pi*rnzg*rho_g)**0.5)
      cac(3) = 0.5/((pi*rnzs*rho_s)*(pi*rnzg*rho_g)**0.75)
      qsn(1:ifull)      = fluxsnow(1:ifull)/rhodz(1:ifull)
      rhodum_s(1:ifull) = qsn(1:ifull)*rhoa(1:ifull,k)
      rhodum_g(1:ifull) = (fluxgraupel(:)+sublflux(:))/dz(:,k)
      where ( fluxgraupel(:)+sublflux(:)>0. .and. qsn(1:ifull)>1.e-10 .and. ttg(1:ifull,k)<tfrz )
        cdt(1:ifull)         = tdt*cgacs*abs(vg2(:,k)-vs2(:,k))*qsn(:)*rhodum_g(:)**0.25*           &
                               (cac(1)*sqrt(rhodum_s(:))+cac(2)*rhodum_s(:)**0.25*rhodum_g(:)**0.25 &
                               +cac(3)*sqrt(rhodum_g(:)))
        dqf(1:ifull)         = max( min( cgfra(:)*qsn(:), qsn(:)*cdt(:)/(1.+cdt(:)) ), 0. )
        qaccf(1:ifull,k)     = qaccf(:,k) + dqf(:)
        fluxsnow(1:ifull)    = fluxsnow(:)    - rhodz(:)*dqf(:)
        fluxgraupel(1:ifull) = fluxgraupel(:) + rhodz(:)*dqf(:)
        rdclfrsnow(1:ifull)  = max( rdclfrsnow(:)*(1.-dqf(:)/qsn(:)), 0. )
        mxclfrsnow(1:ifull)  = max( mxclfrsnow(:)*(1.-dqf(:)/qsn(:)), 0. )
        cftmp(1:ifull)       = max( mxclfrsnow(:) + rdclfrsnow(:) - mxclfrsnow(:)*rdclfrsnow(:), 0. )  
        caccf_g(1:ifull)     = max( caccf_g(1:ifull), max( csfra(:)-cftmp(:), 0. ) )
        csfra(1:ifull)       = cftmp(:)
      end where

  
      ! Snow ------------------------------------------------------------------------------
      alphaf(1:ifull)   = hls*qsatg(1:ifull,k)/(rvap*ttg(1:ifull,k)**2)
      gam(1:ifull,k)    = hlscp*alphaf(:) !(L/cp)*dqsdt (HBG notation)
      sublflux(1:ifull) = 0.
      caccr_s(1:ifull)  = 0.
      caccf_s(1:ifull)  = 0.

      ! Set up the rate constant for snow sublimation
      ! MJT notes - curly and Csbsav depend on vs2(:,k+1), so vs2(:,k) can be updated below
      Tk(1:ifull)    = ttg(1:ifull,k)
      es(1:ifull)    = qsatg(1:ifull,k)*pk(:)/epsil
      Aprpr(1:ifull) = (hls/(rKa*Tk(1:ifull)))*(hls/(rvap*Tk(1:ifull))-1.)
      Bprpr(1:ifull) = rvap*Tk(1:ifull)/((Dva/pk(1:ifull))*es(1:ifull))
      where ( nevapls==-1 .or. (nevapls==-2.and.condx(:)>0..and.k<=ktsav(:)) )
        curly(1:ifull) = 0.
      elsewhere
        curly(1:ifull) = 0.65*slopes(:)**2 + 0.493*slopes(:)*sqrt(slopes(:)*vs2(:,k+1)*rhoa(:,k)/um)
      end where
      ! Define the rate constant for sublimation of snow, omitting factor rhos
      Csbsav(1:ifull) = 4.*curly(:)/(rhoa(:,k)*qsatg(1:ifull,k)*(Aprpr(:)+Bprpr(:))*pi*vs2(:,k+1)*rho_s)
        
      ! The following flag detects max/random overlap clouds
      ! that are separated by a clear layer
      where ( cfsnow(1:ifull,k)<1.e-10 .or. nmr==0 )
        ! combine max overlap from last cloud with net random overlap
        rdclfrsnow(1:ifull) = rdclfrsnow(:) + mxclfrsnow(:) - rdclfrsnow(:)*mxclfrsnow(:)
        mxclfrsnow(1:ifull) = 0.
      end where
      
      fluxsnow(:) = max( fluxsnow(:) + fluxprecipitation(:,k), 0. )
  
      ! Snow fall speed (from Lin et al 1983 - see GFDL AM3)
      where ( cfsnow(1:ifull,k) >= 1.e-10 )
        rhodum_s(:) = max( fluxsnow(:)/dz(:,k), 0. )
        vs2(1:ifull,k) = max( 0.1, 6.63*(rhodum_s(:)/cfsnow(:,k)/942477796.)**0.0625 )
      elsewhere
        vs2(1:ifull,k) = vs2(1:ifull,k+1)
      end where

      ! Set up the parameters for the flux-divergence calculation
      alph(:)        = tdt*vs2(:,k)/dz(:,k)
      foutsnow(:,k)  = 1. - exp(-alph(:))          !analytical
      fthrusnow(:,k) = 1. - foutsnow(:,k)/alph(:)  !analytical

      ! Melt falling snow if > 1 deg C due to rain accretion (based on Lin et al 83, but using 0.65 and 0.493 coeffs
      ! following the UM approach)
      qsn(1:ifull)      = fluxsnow(1:ifull)/rhodz(1:ifull)
      rhodum_s(1:ifull) = qsn(1:ifull)*rhoa(1:ifull,k)
      where ( ttg(1:ifull,k)>tfrz+1. .and. qsn(1:ifull)>1.e-10 )
        cdt(1:ifull)        = tdt*((2.*pi*2.36e-2*rnzs/hlf)*(ttg(1:ifull,k)-tfrz)/rhoa(:,k)         &
                                      -2.*pi*2.11e-5*rnzs*hl/hlf*(qsatg(:,k)-qtg(1:ifull,k)))       &
                                  *((0.65/sqrt(pi*3.e6*0.1e3))*sqrt(rhodum_s(:))                    &
                                   +(0.493*(1.259e-5/2.11e-5)**(1./3.)*1.456943*sqrt(4.8/1.259e-5)  &
                                   /(pi*3.e6*1.e2)**0.65625)*rhodum_s(:)**0.65625)
        qif(1:ifull)        = max( min( qsn(:), qsn(:)*cdt(:)/(1.+cdt(:)) ), 0. ) !Mixing ratio of snow
        dttg(1:ifull)       = -hlfcp*qif(:)
        ttg(1:ifull,k)      = ttg(1:ifull,k) + dttg(:)
        qsatg(1:ifull,k)    = qsatg(1:ifull,k) + gam(:,k)*dttg(:)/hlscp
        fluxmelt(1:ifull,k) = fluxmelt(:,k) + qif(:)*rhodz(:)
        fluxsnow(1:ifull)   = fluxsnow(:)   - qif(:)*rhodz(:)
        rdclfrsnow(1:ifull) = max( rdclfrsnow(:)*(1.-qif(:)/qsn(:)), 0. )
        mxclfrsnow(1:ifull) = max( mxclfrsnow(:)*(1.-qif(:)/qsn(:)), 0. )
        cftmp(1:ifull)      = max( mxclfrsnow(:) + rdclfrsnow(:) - mxclfrsnow(:)*rdclfrsnow(:), 0. )
        cfmelt(1:ifull,k)   = max( cfmelt(:,k), max( csfra(:)-cftmp(:), 0. ) )
        csfra(1:ifull)      = cftmp(:)      
      end where

      ! Melt falling snow if > 2 deg C (based on Lin et al 83)
      qsn(1:ifull) = fluxsnow(:)/rhodz(:) !Mixing ratio of snow
      where ( ttg(1:ifull,k)>tfrz+2. .and. qsn(1:ifull)>1.e-10 )
        qif(1:ifull)        = min( qsn(:), min(tdt/tau_s,1.)*(ttg(1:ifull,k)-tfrz-2.)/hlfcp )
        dttg(1:ifull)       = -hlfcp*qif(:)
        ttg(1:ifull,k)      = ttg(1:ifull,k) + dttg(:)
        qsatg(1:ifull,k)    = qsatg(1:ifull,k) + gam(:,k)*dttg(:)/hlscp
        fluxmelt(1:ifull,k) = fluxmelt(:,k) + qif(:)*rhodz(:)
        fluxsnow(1:ifull)   = fluxsnow(:)   - qif(:)*rhodz(:)
        rdclfrsnow(1:ifull) = max( rdclfrsnow(:)*(1.-qif(:)/qsn(:)), 0. )
        mxclfrsnow(1:ifull) = max( mxclfrsnow(:)*(1.-qif(:)/qsn(:)), 0. )
        cftmp(1:ifull)      = max( mxclfrsnow(:) + rdclfrsnow(:) - mxclfrsnow(:)*rdclfrsnow(:), 0. )
        cfmelt(1:ifull,k)   = max( cfmelt(:,k), max( csfra(:)-cftmp(:), 0. ) )
        csfra(1:ifull)      = cftmp(:)      
      end where
      
      ! Compute the sublimation of snow falling from level k+1 into level k
      ! (Currently treated the same as LDR97 ice sublimation - see UM and ACCESS 1.3)
      fsclr_s(:) = max( (1.-cifr(:,k)-clfr(:,k))*fluxsnow(:), 0. )
      where ( fluxsnow(:)>0. .and. qtg(1:ifull,k)<qsatg(1:ifull,k) ) ! sublime snow
        Csb(1:ifull)      = Csbsav(:)*fluxsnow(:)/tdt
        bf(1:ifull)       = 1. + 0.5*Csb(:)*tdt*(1.+gam(:,k))
        dqs(1:ifull)      = max( 0., tdt*(Csb(:)/bf(:))*(qsatg(1:ifull,k)-qtg(1:ifull,k)) )
        dqs(1:ifull)      = min( dqs(:), (qsatg(1:ifull,k)-qtg(1:ifull,k))/(1.+gam(:,k)) ) !Don't supersat.
        sublflux(1:ifull) = min( dqs(:)*rhodz(:), fsclr_s(:) )
        fluxsnow(1:ifull) = fluxsnow(:) - sublflux(:)
        fsclr_s(1:ifull)  = fsclr_s(:) - sublflux(:)
        dqs(1:ifull)      = sublflux(:)/rhodz(:)
        qsubl(1:ifull,k)  = qsubl(:,k) + dqs(:)
        qtg(1:ifull,k)    = qtg(1:ifull,k) + dqs(:)
        dttg(1:ifull)     = -hlscp*dqs(:)
        ttg(1:ifull,k)    = ttg(1:ifull,k) + dttg(:)
        qsatg(1:ifull,k)  = qsatg(1:ifull,k) + gam(:,k)*dttg(:)/hlscp
      end where

      ! Save flux for the wet deposition scheme.
      pfstayice(:,k) = max( pfstayice(:,k) + fluxsnow(:)*(1.-fthrusnow(:,k))/tdt, 0. )

      ! Accretion of cloud liquid by falling snow (from Lin et al 1983 - psacw)
      ql(1:ifull)       = qlg(1:ifull,k)
      rhodum_s(1:ifull) = (fluxsnow(:)+sublflux(:))/dz(:,k)
      where ( fluxsnow(:)+sublflux(:)>0. .and. qlg(1:ifull,k)>1.e-10 .and. ttg(1:ifull,k)<tfrz )
        cdt(1:ifull)      = tdt*csacw*rhodum_s(:)**0.8125/sqrt(rhoa(:,k))
        dql(1:ifull)      = max( min( csfra(:)*ql(:), ql(:)*cdt(:)/(1.+cdt(:)) ), 0. )
        qlg(1:ifull,k)    = qlg(1:ifull,k) - dql(:)
        qaccr(1:ifull,k)  = qaccr(:,k)     + dql(:)
        fluxsnow(1:ifull) = fluxsnow(:) + rhodz(:)*dql(:)
        dttg(1:ifull)     = hlfcp*dql(:)
        ttg(1:ifull,k)    = ttg(1:ifull,k) + dttg(:)
        qsatg(1:ifull,k)  = qsatg(1:ifull,k) + gam(:,k)*dttg(:)/hlscp
        cftmp(1:ifull)    = max( clfr(:,k)*dql(:)/ql(1:ifull), 0. )
        clfr(1:ifull,k)   = max( clfr(:,k)*(1.-dql(:)/ql(1:ifull)), 0. )
        caccr_s(1:ifull)  = max( caccr_s(1:ifull), cftmp(1:ifull) )
      end where

      ! Accretion of rain by falling snow to form snow (from Lin et al 1983 - psacr)
      cac(1) = 5./((pi*rnzr*rho_r)**1.5*(pi*rnzs*rho_s)**0.25)
      cac(2) = 2./((pi*rnzr*rho_r)**1.25*(pi*rnzs*rho_s)**0.5)
      cac(3) = 0.5/((pi*rnzr*rho_r)*(pi*rnzs*rho_s)**0.75)
      qrn(1:ifull) = fluxrain(1:ifull)/rhodz(1:ifull)
      rhodum_r(1:ifull) = qrn(1:ifull)*rhoa(1:ifull,k)
      rhodum_s(1:ifull) = (fluxsnow(:)+sublflux(:))/dz(:,k)
      where ( fluxsnow(:)+sublflux(:)>0. .and. qrn(1:ifull)>1.e-10 .and. ttg(1:ifull,k)<tfrz )
        cdt(1:ifull)       = tdt*csacr*abs(vs2(:,k)-vl2(:,k))*qrn(:)*rhodum_s(:)**0.25*           &
                             (cac(1)*sqrt(rhodum_r(:))+cac(2)*rhodum_r(:)**0.25*rhodum_s(:)**0.25 &
                             +cac(3)*sqrt(rhodum_s(:)))
        dql(1:ifull)       = max( min( clfra(:)*qrn(:), qrn(:)*cdt(:)/(1.+cdt(:)) ), 0. )
        qaccr(1:ifull,k)   = qaccr(:,k) + dql(:)
        fluxrain(1:ifull)  = fluxrain(:) - rhodz(:)*dql(:)
        fluxsnow(1:ifull)  = fluxsnow(:) + rhodz(:)*dql(:)
        dttg(1:ifull)      = hlfcp*dql(:)
        ttg(1:ifull,k)     = ttg(1:ifull,k) + dttg(:)
        qsatg(1:ifull,k)   = qsatg(1:ifull,k) + gam(:,k)*dttg(:)/hlscp      
        rdclfrliq(1:ifull) = max( rdclfrliq(:)*(1.-dql(:)/qrn(:)), 0. )
        mxclfrliq(1:ifull) = max( mxclfrliq(:)*(1.-dql(:)/qrn(:)), 0. )
        cftmp(1:ifull)     = max( mxclfrliq(:) + rdclfrliq(:) - mxclfrliq(:)*rdclfrliq(:), 0. )      
        caccr_s(1:ifull)   = max( caccr_s(:), max( clfra(:)-cftmp(:), 0. ) )
        clfra(1:ifull)     = cftmp(:)      
      end where

      ! Accretion of rain by falling snow to form graupel (neglected in Lin83 but included in UM)   
    
      ! Accretion of cloud ice by falling snow (from Lin et al 1983 - psaci)
      qf(1:ifull) = fluxice(1:ifull)/rhodz(1:ifull)
      rhodum_s(1:ifull) = (fluxsnow(:)+sublflux(:))/dz(:,k)
      where ( fluxsnow(:)+sublflux(:)>0. .and. qf(1:ifull)>1.e-10 .and. ttg(1:ifull,k)<tfrz )
        n0s(1:ifull)       = 2.e6*exp(-0.12*max(ttg(1:ifull,k)-tfrz,-100.))
        lambdadum(1:ifull) = rhodum_s(:)/(pi*rho_s*n0s(:))
        cdt(1:ifull)       = tdt*27.737*n0s(:)*exp(0.05*max(ttg(1:ifull,k)-tfrz,-100.))*lambdadum(:)**0.8525/sqrt(rhoa(:,k))
        dqf(1:ifull)       = max( min( csfra(:)*qf(:), qf(:)*cdt(:)/(1.+cdt(:)) ), 0. )
        qaccf(1:ifull,k)   = qaccf(:,k) + dqf(:)
        fluxice(1:ifull)   = fluxice(:)  - rhodz(:)*dqf(:)
        fluxsnow(1:ifull)  = fluxsnow(:) + rhodz(:)*dqf(:)
        rdclfrice(1:ifull) = max( rdclfrice(:)*(1.-dqf(:)/qf(:)), 0. )
        mxclfrice(1:ifull) = max( mxclfrice(:)*(1.-dqf(:)/qf(:)), 0. )
        cftmp(1:ifull)     = max( mxclfrice(:) + rdclfrice(:) - mxclfrice(:)*rdclfrice(:), 0. )  
        caccf_s(1:ifull)   = max( caccf_s(:), max( cifra(:)-cftmp(:), 0. ) )
        cifra(1:ifull)     = cftmp(:)
      end where

    
    end if
  
  
    ! Ice ---------------------------------------------------------------------------------
    alphaf(1:ifull)   = hls*qsatg(1:ifull,k)/(rvap*ttg(1:ifull,k)**2)
    gam(1:ifull,k)    = hlscp*alphaf(:) !(L/cp)*dqsdt (HBG notation)
    sublflux(1:ifull) = 0.
    caccr_i(1:ifull)  = 0.
    caccf_i(1:ifull)  = 0.

    ! Set up the rate constant for ice sublimation
    ! MJT notes - curly and Csbsav depend on vi2(:,k+1), so vi2(:,k) can be updated below
    Tk(1:ifull)    = ttg(1:ifull,k)
    es(1:ifull)    = qsatg(1:ifull,k)*pk(:)/epsil
    Aprpr(1:ifull) = (hls/(rKa*Tk(1:ifull)))*(hls/(rvap*Tk(1:ifull))-1.)
    Bprpr(1:ifull) = rvap*Tk(1:ifull)/((Dva/pk(1:ifull))*es(1:ifull))
    where ( nevapls==-1 .or. (nevapls==-2.and.condx(:)>0..and.k<=ktsav(:)) )
      curly(1:ifull) = 0.
    elsewhere
      curly(1:ifull) = 0.65*slopes(:)**2 + 0.493*slopes(:)*sqrt(slopes(:)*vi2(:,k+1)*rhoa(:,k)/um) !Factor in curly brackets
    end where
    ! Define the rate constant for sublimation of snow, omitting factor rhoi
    Csbsav(1:ifull) = 4.*curly(:)/(rhoa(:,k)*qsatg(1:ifull,k)*(Aprpr(:)+Bprpr(:))*pi*vi2(:,k+1)*rho_s)
    
    ! The following flag detects max/random overlap clouds
    ! that are separated by a clear layer
    where ( cifr(1:ifull,k)<1.e-10 .or. nmr==0 )
      ! combine max overlap from last cloud with net random overlap
      rdclfrice(1:ifull) = rdclfrice(:) + mxclfrice(:) - rdclfrice(:)*mxclfrice(:)
      mxclfrice(1:ifull) = 0.
    end where
  
    ! Set up ice fall speed field
    if ( ncloud >= 3 ) then
      ! Ice fall speed from Lin et al 1983
      where ( cifr(1:ifull,k) >= 1.e-10 )
        rhodum_i(:) = max( fluxice(:)/dz(:,k), 0. )
        vi2(1:ifull,k) = max( 0.1, 3.29*(rhodum_i(:)/cifr(:,k))**0.16 )
      end where
    else
      select case(abs(ldr))
        case(1)
          ! Ice fall speed from LDR 1997
          where ( cifr(1:ifull,k) >= 1.e-10 )
            vi2(1:ifull,k) = 3.23*(rhoi(:,k)/cifr(:,k))**0.17
          end where
        case(2)
          where ( cifr(:,k) >= 1.e-10 )
            vi2(:,k) = 0.9*3.23*(rhoi(:,k)/cifr(:,k))**0.17
          end where
        case(3)
          where ( cifr(1:ifull,k) >= 1.e-10 )
            vi2(1:ifull,k) = max( 0.1, 2.05+0.35*log10(rhoi(:,k)/rhoa(:,k)/cifr(1:ifull,k)) )
          end where
        case(4)
          where ( cifr(1:ifull,k) >= 1.e-10 )
            vi2(1:ifull,k) = 1.4*3.23*(rhoi(:,k)/cifr(1:ifull,k))**0.17
          end where
        case(11)
          ! following are alternative slightly-different versions of above
          ! used for I runs from 29/4/05 till 30/8/05
          ! for given qfg, large cifr implies small ice crystals, 
          ! with a small fall speed. 
          ! Note that for very small qfg, cifr is small.
          ! But rhoi is like qfg, so ratio should also be small and OK.
          vi2(1:ifull,k) = max( vi2(1:ifull,k+1),3.23*(rhoi(:,k)/max(cifr(1:ifull,k),1.e-30))**0.17 )
        case(22)
          vi2(1:ifull,k) = max( vi2(1:ifull,k+1),0.9*3.23*(rhoi(:,k)/max(cifr(1:ifull,k),1.e-30))**0.17 )
        case(33)
          ! following max gives vi2=.1 for qfg=cifr=0
          vi2(1:ifull,k) = max( vi2(1:ifull,k+1),2.05+0.35*log10(max(rhoi(:,k)/rhoa(:,k),2.68e-36)/max(cifr(1:ifull,k),1.e-30)) )
      end select
    end if

    ! Set up the parameters for the flux-divergence calculation
    alph(:)       = tdt*vi2(:,k)/dz(:,k)
    foutice(:,k)  = 1. - exp(-alph(:))         !analytical
    fthruice(:,k) = 1. - foutice(:,k)/alph(:)  !analytical
  
    ! Melt falling ice if > 0 deg C
    qif(1:ifull) = fluxice(:)/rhodz(:)      !Mixing ratio of ice
    where ( ttg(1:ifull,k)>tfrz .and. fluxice(:)>0. .and. ncloud<=2 )
      dttg(1:ifull)       = -hlfcp*qif(:)
      ttg(1:ifull,k)      = ttg(1:ifull,k) + dttg(:)
      qsatg(1:ifull,k)    = qsatg(1:ifull,k) + gam(:,k)*dttg(:)/hlscp
      fluxmelt(1:ifull,k) = fluxmelt(:,k) + fluxice(:)
      fluxice(1:ifull)    = 0.
      rdclfrice(1:ifull)  = 0.
      mxclfrice(1:ifull)  = 0.
      cfmelt(1:ifull,k)   = max( cfmelt(:,k), cifra(:) )
      cifra(1:ifull)      = 0.
    end where

    ! Melt falling ice if > 2 deg C ( following Lin et al 83)
    qsn(1:ifull) = fluxice(:)/rhodz(:)      !Mixing ratio of ice
    where ( ttg(1:ifull,k)>tfrz+2. .and. qsn(1:ifull)>1.e-10 .and. ncloud>=3 )
      qif(1:ifull)        = min( qsn(:), min(tdt/tau_i,1.)*(ttg(1:ifull,k)-tfrz-2.)/hlfcp )  
      dttg(1:ifull)       = -hlfcp*qif(:)
      ttg(1:ifull,k)      = ttg(1:ifull,k) + dttg(:)
      qsatg(1:ifull,k)    = qsatg(1:ifull,k) + gam(:,k)*dttg(:)/hlscp
      fluxmelt(1:ifull,k) = fluxmelt(:,k) + qif(:)*rhodz(:)
      fluxice(1:ifull)    = fluxice(:)    - qif(:)*rhodz(:)
      rdclfrice(1:ifull)  = max( rdclfrice(:)*(1.-qif(:)/qsn(:)), 0. )
      mxclfrice(1:ifull)  = max( mxclfrice(:)*(1.-qif(:)/qsn(:)), 0. )
      cftmp(1:ifull)      = max( mxclfrice(:) + rdclfrice(:) - mxclfrice(:)*rdclfrice(:), 0. )
      cfmelt(1:ifull,k)   = max( cfmelt(:,k), max( cifra(:)-cftmp(:), 0. ) )
      cifra(1:ifull)      = cftmp(:)
    end where    
    
    ! Compute the sublimation of ice falling from level k+1 into level k
    fsclr_i(:) = max( (1.-cifr(:,k)-clfr(:,k))*fluxice(:), 0. )
    where ( fluxice(:)>0. .and. qtg(1:ifull,k)<qsatg(1:ifull,k) ) ! sublime ice
      Csb(1:ifull)      = Csbsav(:)*fluxice(:)/tdt
      bf(1:ifull)       = 1. + 0.5*Csb(:)*tdt*(1.+gam(:,k))
      dqs(1:ifull)      = max( 0., tdt*(Csb(:)/bf(:))*(qsatg(1:ifull,k)-qtg(1:ifull,k)) )
      dqs(1:ifull)      = min( dqs(:), (qsatg(1:ifull,k)-qtg(1:ifull,k))/(1.+gam(:,k)) ) !Don't supersat.
      sublflux(1:ifull) = min( dqs(:)*rhodz(:), fsclr_i(:) )
      fluxice(1:ifull)  = fluxice(:) - sublflux(:)
      fsclr_i(1:ifull)  = fsclr_i(:) - sublflux(:)
      dqs(1:ifull)      = sublflux(:)/rhodz(:)
      qsubl(1:ifull,k)  = qsubl(:,k) + dqs(:)
      qtg(1:ifull,k)    = qtg(1:ifull,k) + dqs(:)
      dttg(1:ifull)     = -hlscp*dqs(:)
      ttg(1:ifull,k)    = ttg(1:ifull,k) + dttg(:)
      qsatg(1:ifull,k)  = qsatg(1:ifull,k) + gam(:,k)*dttg(:)/hlscp
    end where

    ! Save flux for the wet deposition scheme.
    pfstayice(:,k) = max( pfstayice(:,k) + fluxice(:)*(1.-fthruice(:,k))/tdt, 0. )
  
    ! Accretion of cloud liquid by falling ice (neglected in Lin et al 1983, but
    ! included in UM and ACCESS 1.3 as piacw)
    ! This calculation uses the incoming fluxice without subtracting sublimation
    ! (since subl occurs only outside cloud), so add sublflux back to fluxice.
    ql(1:ifull) = qlg(1:ifull,k)
    where ( fluxice(:)+sublflux(:)>0. .and. ql(1:ifull)>1.e-10 .and. ncloud<=2 )
      cdt(1:ifull)     = Eac*slopes(:)*(fluxice(:)+sublflux(:))/(2.*rhosno)
      dql(1:ifull)     = max( min( cifra(:)*ql(:), ql(:)*cdt(:)/(1.+0.5*cdt(:)) ), 0. )
      qlg(1:ifull,k)   = qlg(1:ifull,k) - dql(:)
      qaccr(1:ifull,k) = qaccr(:,k)     + dql(:)
      fluxice(1:ifull) = fluxice(:) + rhodz(:)*dql(:)
      dttg(1:ifull)    = hlfcp*dql(:)
      ttg(1:ifull,k)   = ttg(1:ifull,k) + dttg(:)
      qsatg(1:ifull,k) = qsatg(1:ifull,k) + gam(:,k)*dttg(:)/hlscp
      cftmp(1:ifull)   = max( clfr(:,k)*dql(:)/ql(1:ifull), 0. )
      clfr(1:ifull,k)  = max( clfr(:,k)*(1.-dql(:)/ql(1:ifull)), 0.)
      caccr_i(1:ifull) = max( caccr_i(1:ifull), cftmp(1:ifull) )
    end where
  
    ! Accretion of rain by falling ice to produce ice (from Lin et al 1983 - piacr)
    ! (see UM and ACCESS 1.3 piacr-c for an alternate formulation)
    qrn(1:ifull) = fluxrain(1:ifull)/rhodz(1:ifull)
    qf(1:ifull)  = (fluxice(1:ifull)+sublflux(1:ifull))/rhodz(1:ifull)
    where ( fluxice(:)+sublflux(:)>0. .and. qrn(1:ifull)>1.e-10 .and. ncloud>=3 )
      cdt(1:ifull)      = tdt*c_piacr*qf(:)/sqrt(rhoa(:,k))
      dql(1:ifull)      = max( min( cifra(:)*qrn(:), qrn(:)*cdt(:)/(1.+cdt(:)) ), 0. )
      qaccr(1:ifull,k)  = qaccr(:,k) + dql(:)
      fluxrain(1:ifull) = fluxrain(:) - rhodz(:)*dql(:)
      fluxice(1:ifull)  = fluxice(:)  + rhodz(:)*dql(:)
      dttg(1:ifull)     = hlfcp*dql(:)
      ttg(1:ifull,k)    = ttg(1:ifull,k) + dttg(:)
      qsatg(1:ifull,k)  = qsatg(1:ifull,k) + gam(:,k)*dttg(:)/hlscp
      rdclfrliq(1:ifull) = max( rdclfrliq(:)*(1.-dql(:)/qrn(:)), 0. )
      mxclfrliq(1:ifull) = max( mxclfrliq(:)*(1.-dql(:)/qrn(:)), 0. )
      cftmp(1:ifull)     = max( mxclfrliq(:) + rdclfrliq(:) - mxclfrliq(:)*rdclfrliq(:), 0. )      
      caccr_i(1:ifull)   = max( caccr_i(:), max( clfra(:)-cftmp(:), 0. ) )
      clfra(1:ifull)     = cftmp(:)    
    end where

    ! Accretion of rain by falling ice to produce graupel (Neglected in Lin et al 1983)
    ! (see UM and ACCESS 1.3 piacr-g for an alternate formulation)

  
    ! Rain --------------------------------------------------------------------------------
    evap(:) = 0.

    ! The following flag detects maximum/random overlap clouds
    ! that are separated by a clear layer
    where ( cfrain(1:ifull,k)<1.e-10 .or. nmr==0 )
      ! combine max overlap from above cloud with net random overlap
      rdclfrliq(:) = rdclfrliq(:) + mxclfrliq(:) - rdclfrliq(:)*mxclfrliq(:)
      mxclfrliq(:) = 0.
    end where

    ! Add flux of melted snow to fluxrain
    fluxrain(:) = max( fluxrain(:) + fluxmelt(:,k) + fluxauto(:,k), 0. )
  
    ! Calculate rain fall speed (MJT)
    if ( ncloud > 1 ) then
      Fr(:)         = max( fluxrain(:)/tdt/max( clfra(:), 1.e-15 ), 0. )
      vl2(:,k)      = 11.3*Fr(:)**(1./9.)/sqrt(rhoa(:,k))  !Actual fall speed
      vl2(:,k)      = max( vl2(:,k), 0.1 )
      alph(:)       = tdt*vl2(:,k)/dz(:,k)
      foutliq(:,k)  = 1. - exp(-alph(:))         !analytical
      fthruliq(:,k) = 1. - foutliq(:,k)/alph(:)  !analytical
    else
      foutliq(:,k)  = 1.
      fthruliq(:,k) = 1.
    end if
  
    ! Evaporation of rain
    qpf(:)     = max( fluxrain(:)/rhodz(:), 0. ) !Mix ratio of rain which falls into layer
    clrevap(:) = max( (1.-clfr(:,k))*qpf(:), 0. )
    do mg = 1,ifull
      if ( fluxrain(mg) > 0. ) then
        qsatg(mg,k) = qsati(pk(mg),ttg(mg,k))
        if ( ttg(mg,k)<tfrz .and. ttg(mg,k)>=tice ) then
          qsl = qsatg(mg,k) + epsil*esdiffx(ttg(mg,k))/pk(mg)
        else
          qsl = qsatg(mg,k)
        end if             !qsl is qs value over liquid surface
        Tk(mg)   = ttg(mg,k)
        es(mg)   = qsl*pk(mg)/epsil 
        Apr      = (hl/(rKa*Tk(mg)))*(hl/(rvap*Tk(mg))-1.)
        Bpr      = rvap*Tk(mg)/((Dva/pk(mg))*es(mg))
        Fr(mg)   = max( fluxrain(mg)/tdt/max( clfra(mg), 1.e-15 ), 0. )
        Cev      = clfra(mg)*3.8e2*sqrt(max( Fr(mg)/rhoa(mg,k), 0. ))/(qsl*(Apr+Bpr))
        dqsdt    = hl*qsl/(rvap*ttg(mg,k)**2)
        bl       = 1. + 0.5*Cev*tdt*(1.+hlcp*dqsdt)
        evap(mg) = tdt*(Cev/bl)*(qsl-qtg(mg,k))
        satevap  = (qsl-qtg(mg,k))/(1.+hlcp*dqsdt) !Evap to saturate
        ! Vl2=11.3*Fr**(1./9.)/sqrt(rhoa(mg,k))    !Actual fall speed
        ! Vl2=5./sqrt(rhoa(mg,k))                  !Nominal fall speed
        evap(mg) = max( 0., min( evap(mg), satevap, qpf(mg), clrevap(mg) ) )
      end if
    end do
    if ( nevapls == -1 ) then
      evap(1:ifull) = 0.
    else if ( nevapls == -2 ) then
      where ( k<=ktsav(1:ifull) .and. condx(1:ifull)>0. )
        evap(1:ifull) = 0.
      end where
    else if ( nevapls==-3 ) then
      evap(1:ifull) = 0.5*evap(:)
    else if ( nevapls==-4 ) then
      where ( k<=ktsav(1:ifull) .and. condx(1:ifull)>0. )
        evap(1:ifull) = 0.5*evap(:) ! usual
      end where
    end if
    qevap(1:ifull,k) = qevap(:,k) + evap(:)
    qtg(1:ifull,k)   = qtg(:,k) + evap(:)
    ttg(1:ifull,k)   = ttg(1:ifull,k) - hlcp*evap(:)
    frclr(1:ifull)   = rhodz(:)*(clrevap(:)-evap(:)) !over tdt

    ! Freezing rain to produce graupel (pgfr)
    ! (Neglected in UM and ACCESS 1.3)
    cgfr(1) = 20.e2*pi*pi*rnzr*rho_r/(pi*rnzr*rho_r)**1.75
    qrn(1:ifull) = (fluxrain(1:ifull)+evap(1:ifull))/rhodz(1:ifull)
    rhodum_r(1:ifull) = qrn(1:ifull)*rhoa(1:ifull,k)
    where ( qrn(1:ifull)>1.e-10 .and. ttg(1:ifull,k)<tfrz .and. ncloud>=3 )
      ! MJT notes - limit temperature to -100 C to avoid overflow with single precision
      cdt(1:ifull)         = tdt*cgfr(1)*(exp(-0.66*max( ttg(1:ifull,k)-tfrz, -100. ))-1.)  &
                             *(rhodum_r(:))**1.75/rhoa(:,k)
      dql(1:ifull)         = max( min( qrn(1:ifull), qrn(1:ifull)*cdt(:)/(1.+cdt(:)) ), 0. )
      fluxrain(1:ifull)    = fluxrain(:)    - rhodz(:)*dql(:)
      fluxgraupel(1:ifull) = fluxgraupel(:) + rhodz(:)*dql(:)
      dttg(1:ifull)        = hlfcp*dql(:)
      ttg(1:ifull,k)       = ttg(1:ifull,k) + dttg(:)
      qsatg(1:ifull,k)     = qsatg(1:ifull,k) + gam(:,k)*dttg(:)/hlscp
      rdclfrliq(1:ifull)   = max( rdclfrliq(:)*(1.-dql(:)/qrn(:)), 0. )
      mxclfrliq(1:ifull)   = max( mxclfrliq(:)*(1.-dql(:)/qrn(:)), 0. )
      cftmp(1:ifull)       = max( mxclfrliq(:) + rdclfrliq(:) - mxclfrliq(:)*rdclfrliq(:), 0. )
      cfgraupel(1:ifull,k) = max( cfgraupel(:,k), max( clfra(:)-cftmp(:), 0. ) )
      clfra(1:ifull)       = cftmp(:)      
    end where
  
    ! store liquid flux for aerosols
    pfstayliq(:,k) = max( fluxrain(:)*(1.-fthruliq(:,k))/tdt, 0. )
  
    ! Now do the collection of liquid cloud by rain term (cf. pracc in Lin83).
    where ( fluxrain(:) > 0. )
      Fr(1:ifull)       = max( fluxrain(:)/tdt/max( clfra(:), 1.e-15 ), 0. )
      mxovr(1:ifull)    = min( mxclfrliq(:), clfr(:,k) )          ! max overlap
      mxovr(1:ifull)    = max( cfrain(1:ifull,k), mxovr(:) )
      rdovr(1:ifull)    = rdclfrliq(:)*clfr(:,k)                  ! rnd overlap
      cfrain(1:ifull,k) = mxovr(:) + rdovr(:) - mxovr(:)*rdovr(:) ! combine collection
    elsewhere
      Fr(1:ifull) = 0.
    end where
    ! The collection term comprises collection by stratiform rain falling from
    ! above (Fr), stratiform rain released in this grid box (Frb).
    ! Frb term now done above.
    fcol(1:ifull)     = min( 1., mxclfrliq(:)/(1.e-20+clfr(:,k)) )    !max overlap
    fcol(1:ifull)     = fcol(:) + rdclfrliq(:) - fcol(:)*rdclfrliq(:) !rnd overlap
    cdt(1:ifull)      = tdt*Ecol*0.24*fcol(:)*pow75(Fr(:))
    prscav(1:ifull,k) = prscav(1:ifull,k) + tdt*0.24*fcol(:)*pow75(Fr(:))  !Strat only
    coll(1:ifull)     = min( qlg(1:ifull,k), qlg(1:ifull,k)*cdt(:)/(1.+0.5*cdt(:)) )
    qcoll(1:ifull,k)  = qcoll(:,k) + coll(:)
    qlg(1:ifull,k)    = qlg(1:ifull,k) - coll(:)
    fluxrain(1:ifull) = fluxrain(:) + coll(:)*rhodz(:)
  
    ! Accretion of cloud snow by rain to produce rain (from Lin et al 1983 - pracs)
    cac(1) = 5./((pi*rnzs*rho_s)**1.5*(pi*rnzr*rho_r)**0.25)
    cac(2) = 2./((pi*rnzs*rho_s)**1.25*(pi*rnzr*rho_r)**0.5)
    cac(3) = 0.5/((pi*rnzs*rho_s)*(pi*rnzr*rho_r)**0.75)
    qsn(1:ifull)      = fluxsnow(1:ifull)/rhodz(1:ifull)
    rhodum_s(1:ifull) = fluxsnow(1:ifull)/dz(1:ifull,k)
    rhodum_r(1:ifull) = (fluxrain(1:ifull)+evap(1:ifull))/dz(1:ifull,k)
    where ( fluxrain(1:ifull)+evap(1:ifull)>0. .and. qsn(1:ifull)>1.e-10 .and. ttg(1:ifull,k)>tfrz+1. .and. ncloud>=3 )
      cdt(1:ifull)         = tdt*cracs*abs(vl2(:,k)-vs2(:,k))*qsn(:)*rhodum_r(:)**0.25*           &
                             (cac(1)*sqrt(rhodum_s(:))+cac(2)*rhodum_s(:)**0.25*rhodum_r(:)**0.25 &
                             +cac(3)*sqrt(rhodum_r(:)))/sqrt(rhoa(:,k))
      dqf(1:ifull)         = max( min( clfra(:)*qsn(:), qsn(:)*cdt(:)/(1.+cdt(:)) ), 0. )
      fluxsnow(1:ifull)    = fluxsnow(:) - rhodz(:)*dqf(:)
      fluxrain(1:ifull)    = fluxrain(:) + rhodz(:)*dqf(:)
      dttg(1:ifull)        = hlfcp*dqf(:)
      ttg(1:ifull,k)       = ttg(1:ifull,k) - dttg(:)
      qsatg(1:ifull,k)     = qsatg(1:ifull,k) - gam(:,k)*dttg(:)/hlscp      
      rdclfrsnow(1:ifull)  = max( rdclfrsnow(:)*(1.-dqf(:)/qsn(:)), 0. )
      mxclfrsnow(1:ifull)  = max( mxclfrsnow(:)*(1.-dqf(:)/qsn(:)), 0. )
      cftmp(1:ifull)       = max( mxclfrsnow(:) + rdclfrsnow(:) - mxclfrsnow(:)*rdclfrsnow(:), 0. )  
      csfra(1:ifull)       = cftmp(:)
    end where
  
    ! subtract evaporated rain
    fluxrain(:) = max( fluxrain(:)-rhodz(:)*evap(:), 0. )
  
  
    ! Liquid ------------------------------------------------------------------------------
    ! (Currently cloud droplet settling is negected, although included in UM and ACCESS 1.3)


  
    ! Update fluxes and area fractions for graupel, snow, ice and rain
  
    if ( ncloud >= 3 ) then

      
      ! Grauple
      ! calculate maximum and random overlap for falling graupel
      cfgraupel(1:ifull,k) = min( 1., cfgraupel(1:ifull,k)+caccr_g(1:ifull)-cfgraupel(1:ifull,k)*caccr_g(1:ifull) )  ! rnd overlap
      cfgraupel(1:ifull,k) = max( cfgraupel(1:ifull,k), caccf_g(1:ifull) )                                           ! max overlap
      where ( fluxgraupel(:) < 1.e-15 )
        rdclfrgraupel(1:ifull) = 0.
        mxclfrgraupel(1:ifull) = 0.
      end where
      !max overlap
      mxclfrgraupel(1:ifull) = max( mxclfrgraupel(:), cfgraupel(1:ifull,k) )                
      !rnd overlap the mx and rd ice fractions
      cgfra(1:ifull)         = max( 1.e-15, mxclfrgraupel(:)+rdclfrgraupel(:)-mxclfrgraupel(:)*rdclfrgraupel(:) ) 
      ! Compute fluxes into the box
      cffluxin(:) = cgfra(:) - cfgraupel(:,k)
      rhogin(:)   = fluxgraupel(:)/dz(:,k)
      ! Compute the fluxes of snow leaving the box
      cffluxout(:) = cfgraupel(:,k)*foutgraupel(:,k)
      rhogout(:)   = rhog(:,k)*foutgraupel(:,k)
      ! Update the rhos and cfsnow fields
      cfgraupelfall(1:ifull,k) = max( cfgraupel(1:ifull,k) - cffluxout(:) + cffluxin(:)*(1.-fthrugraupel(:,k)), 0. )
      rhog(1:ifull,k)          = max( rhog(:,k) - rhogout(:) + rhogin(:)*(1.-fthrugraupel(:,k)), 0. )
      fluxgraupel(1:ifull)     = max( rhogout(:)*dz(:,k) + fluxgraupel(:)*fthrugraupel(:,k), 0. )
      ! Now fluxgraupel is flux leaving layer k
      fluxg(1:ifull,k)         = max( fluxg(:,k) + fluxgraupel(:), 0. )
   
    
      ! Snow
      ! calculate maximum and random overlap for falling snow
      cfsnow(1:ifull,k) = min( 1., cfsnow(1:ifull,k)+caccr_s(1:ifull)-cfsnow(1:ifull,k)*caccr_s(1:ifull) ) ! rnd overlap
      cfsnow(1:ifull,k) = max( cfsnow(1:ifull,k), caccf_s(1:ifull) )                                       ! max overlap
      where ( fluxsnow(:) < 1.e-15 )
        rdclfrsnow(1:ifull) = 0.
        mxclfrsnow(1:ifull) = 0.
      end where
      !max overlap
      mxclfrsnow(1:ifull) = max( mxclfrsnow(:), cfsnow(1:ifull,k) )
      !rnd overlap the mx and rd snow fractions
      csfra(1:ifull)      = max( 1.e-15, mxclfrsnow(:)+rdclfrsnow(:)-mxclfrsnow(:)*rdclfrsnow(:) ) 
      ! Compute fluxes into the box
      cffluxin(:) = csfra(:) - cfsnow(:,k)
      rhosin(:)   = fluxsnow(:)/dz(:,k)
      ! Compute the fluxes of snow leaving the box
      cffluxout(:) = cfsnow(:,k)*foutsnow(:,k)
      rhosout(:)   = rhos(:,k)*foutsnow(:,k)
      ! Update the rhos and cfsnow fields
      cfsnowfall(1:ifull,k) = max( cfsnow(1:ifull,k) - cffluxout(:) + cffluxin(:)*(1.-fthrusnow(:,k)), 0. )
      rhos(1:ifull,k)       = max( rhos(:,k) - rhosout(:) + rhosin(:)*(1.-fthrusnow(:,k)), 0. )
      fluxsnow(1:ifull)     = max( rhosout(:)*dz(:,k) + fluxsnow(:)*fthrusnow(:,k), 0. )
      ! Now fluxsnow is flux leaving layer k
      fluxs(1:ifull,k)      = max( fluxs(:,k) + fluxsnow(:), 0. )

    
    end if ! ncloud>=3

  
    ! Ice
    ! calculate maximum and random overlap for falling ice
    cifr(1:ifull,k) = min( 1., cifr(1:ifull,k)+caccr_i(:)-cifr(1:ifull,k)*caccr_i(:) )  ! rnd overlap
    where ( fluxice(:) < 1.e-15 )
      rdclfrice(1:ifull) = 0.
      mxclfrice(1:ifull) = 0.
    end where
    mxclfrice(1:ifull) = max( mxclfrice(:), cifr(1:ifull,k) )                             !max overlap
    cifra(1:ifull)     = max( 0.01, mxclfrice(:)+rdclfrice(:)-mxclfrice(:)*rdclfrice(:) ) !rnd overlap the mx and rd ice fractions
    ! Compute fluxes into the box
    if ( ncloud <= 2 ) then
      where ( rica(:) > 0. )
        rhoiin(:)   = max( fluxice(:)/dz(:,k), 0. )
        cffluxin(:) = min( rhoiin(:)/rica(:), 1. )
      elsewhere
        rhoiin(:)   = 0.
        cffluxin(:) = 0.
      end where
      where ( cifr(:,k) >= 1.e-10 ) 
        rica(:) = max( rhoi(:,k)/cifr(:,k), 0. ) ! in cloud rhoi
      end where
    else
      rhoiin(:)   = max( fluxice(:)/dz(:,k), 0. )
      cffluxin(:) = cifra(:) - cifr(:,k)
    end if
    ! Compute the fluxes of ice leaving the box
    where ( cifr(:,k) >= 1.e-10 )
      cffluxout(:) = cifr(:,k)*foutice(:,k)
      rhoiout(:)   = rhoi(:,k)*foutice(:,k)
    elsewhere
      cffluxout(:) = 0.
      rhoiout(:)   = 0.
    end where
    ! Update the rhoi and cifr fields
    cifr(1:ifull,k)  = max( min( 1.-clfr(1:ifull,k), cifr(1:ifull,k) - cffluxout(:) + cffluxin(:)*(1.-fthruice(:,k)) ), 0. )
    rhoi(1:ifull,k)  = max( rhoi(:,k) - rhoiout(:) + rhoiin(:)*(1.-fthruice(:,k)), 0. )
    fluxice(1:ifull) = max( rhoiout(:)*dz(:,k) + fluxice(:)*fthruice(:,k), 0. )
    ! Now fluxice is flux leaving layer k
    fluxi(1:ifull,k) = max( fluxi(:,k) + fluxice(:), 0. )

  
    ! Rain
    ! Calculate the raining cloud cover down to this level, for stratiform (clfra).
    cfrain(1:ifull,k) = min( 1., cfrain(1:ifull,k)+cfmelt(:,k)-cfrain(1:ifull,k)*cfmelt(:,k) ) ! rnd overlap
    where ( fluxrain(:) < 1.e-15 )
      rdclfrliq(:) = 0.
      mxclfrliq(:) = 0.
    end where
    mxclfrliq(:) = max( mxclfrliq(:), cfrain(1:ifull,k) )                             !max overlap
    clfra(:)     = max( 1.e-15, rdclfrliq(:)+mxclfrliq(:)-rdclfrliq(:)*mxclfrliq(:) ) !rnd overlap the mx and rd rain fractions
    ! Compute fluxes into the box
    cffluxin(:) = clfra(:) - cfrain(:,k)
    rhorin(:)   = fluxrain(:)/dz(:,k)
    ! Compute the fluxes of rain leaving the box
    ! Use the flux-divergent form as in Rotstayn (QJRMS, 1997)
    cffluxout(:) = cfrain(:,k)*foutliq(:,k)
    rhorout(:)   = rhor(:,k)*foutliq(:,k)
    ! Update the rhor and cfrainfall fields
    cfrainfall(1:ifull,k) = max( cfrain(1:ifull,k) - cffluxout(:) + cffluxin(:)*(1.-fthruliq(:,k)), 0. )
    rhor(1:ifull,k)       = max( rhor(:,k) - rhorout(:) + rhorin(:)*(1.-fthruliq(:,k)), 0. )
    fluxrain(1:ifull)     = max( rhorout(:)*dz(:,k) + fluxrain(:)*fthruliq(:,k), 0. )
    ! Now fluxrain is flux leaving layer k
    fluxr(1:ifull,k)      = max( fluxr(:,k) + fluxrain(:), 0. )
 
  
    ! Save sedimentation rate for aerosol scheme
    pqfsed(:,k) = pqfsed(:,k) + (fluxgraupel(:)*foutgraupel(:,k)+fluxsnow(:)*foutsnow(:,k)+max(fluxice(:),1.e-30)*foutice(:,k)) &
                                /max( fluxgraupel(:)+fluxsnow(:)+fluxice(:), 1.e-30 )
  
  end do


  ! Re-create qrg, qfg, qsng and agrg fields
  qrg(1:ifull,kl)     = 0.
  qrg(1:ifull,1:kl-1) = rhor(1:ifull,1:kl-1)/rhoa(1:ifull,1:kl-1)
  qfg(1:ifull,1:kl)   = rhoi(1:ifull,1:kl)/rhoa(1:ifull,1:kl)
  qsng(1:ifull,1:kl)  = rhos(1:ifull,1:kl)/rhoa(1:ifull,1:kl)
  qgrg(1:ifull,1:kl)  = rhog(1:ifull,1:kl)/rhoa(1:ifull,1:kl)


  ! store precip, snow and hail
  precs(1:ifull) = precs(1:ifull) + fluxr(1:ifull,1) + fluxi(1:ifull,1) + fluxs(1:ifull,1) + fluxg(1:ifull,1)
  preci(1:ifull) = preci(1:ifull) + fluxi(1:ifull,1) + fluxs(1:ifull,1)
  precg(1:ifull) = precg(1:ifull) + fluxg(1:ifull,1)


  ! Remove small amounts of cloud
  where ( qlg(1:ifull,1:kl)<1.e-10 .or. clfr(1:ifull,1:kl)<1.e-5 )
    qtg(1:ifull,1:kl)  = qtg(1:ifull,1:kl) + qlg(1:ifull,1:kl)
    ttg(1:ifull,1:kl)  = ttg(1:ifull,1:kl) - hlcp*qlg(1:ifull,1:kl)
    qlg(1:ifull,1:kl)  = 0.
    clfr(1:ifull,1:kl) = 0.
  end where
  where ( qfg(1:ifull,1:kl)<1.e-10 .or. cifr(1:ifull,1:kl)<1.e-5 )
    qtg(1:ifull,1:kl)  = qtg(1:ifull,1:kl) + qfg(1:ifull,1:kl)
    ttg(1:ifull,1:kl)  = ttg(1:ifull,1:kl) - hlscp*qfg(1:ifull,1:kl)
    qfg(1:ifull,1:kl)  = 0.
    cifr(1:ifull,1:kl) = 0.
  end where
  where ( qrg(1:ifull,1:kl)<1.e-10 .or. cfrainfall(1:ifull,1:kl)<1.e-5 )
    qtg(1:ifull,1:kl)        = qtg(1:ifull,1:kl) + qrg(1:ifull,1:kl)
    ttg(1:ifull,1:kl)        = ttg(1:ifull,1:kl) - hlcp*qrg(1:ifull,1:kl)
    qrg(1:ifull,1:kl)        = 0.
    cfrainfall(1:ifull,1:kl) = 0.
  end where
  where ( qsng(1:ifull,1:kl)<1.e-10 .or. cfsnowfall(1:ifull,1:kl)<1.e-5 )
    qtg(1:ifull,1:kl)        = qtg(1:ifull,1:kl) + qsng(1:ifull,1:kl)
    ttg(1:ifull,1:kl)        = ttg(1:ifull,1:kl) - hlscp*qsng(1:ifull,1:kl)
    qsng(1:ifull,1:kl)       = 0.
    cfsnowfall(1:ifull,1:kl) = 0.
  end where
  where ( qgrg(1:ifull,1:kl)<1.e-10 .or. cfgraupelfall(1:ifull,1:kl)<1.e-5 )
    qtg(1:ifull,1:kl)           = qtg(1:ifull,1:kl) + qgrg(1:ifull,1:kl)
    ttg(1:ifull,1:kl)           = ttg(1:ifull,1:kl) - hlscp*qgrg(1:ifull,1:kl)
    qgrg(1:ifull,1:kl)          = 0.
    cfgraupelfall(1:ifull,1:kl) = 0.
  end where

  
end do


! estimate average sedimentation rate for aerosol scheme
pqfsed(:,:) = pqfsed(:,:)/real(ncount)

! estimate average slopes
pslopes(:,:) = pslopes(:,:)/real(ncount)


!      Adjust cloud fraction (and cloud cover) after precipitation
if ( nmaxpr==1 .and. mydiag ) then
  write(6,*) 'diags from newrain for idjd ',idjd
  write (6,"('cfrac         ',9f8.3/6x,9f8.3)") cfrac(idjd,:)
  write (6,"('cfrainfall    ',9f8.3/6x,9f8.3)") cfrainfall(idjd,:)
  write (6,"('cfsnowfall    ',9f8.3/6x,9f8.3)") cfsnowfall(idjd,:)
  write (6,"('cfgraupelfall ',9f8.3/6x,9f8.3)") cfgraupelfall(idjd,:)
  write (6,"('cftemp        ',9f8.3/6x,9f8.3)") cifr(idjd,:)+clfr(idjd,:)
  write (6,"('ccov_in       ',9f8.3/6x,9f8.3)") ccov(idjd,:)
end if

! Diagnostics for debugging
if (diag .and. mydiag) then  ! JLM
  write(6,*) 'vi2',vi2(idjd,:)
  write(6,*) 'cfraci ',cfrac(idjd,:)
  write(6,*) 'cifr',cifr(idjd,:)
  write(6,*) 'clfr',clfr(idjd,:)
  write(6,*) 'ttg',ttg(idjd,:)
  write(6,*) 'qsatg',qsatg(idjd,:)         
  write(6,*) 'qlg',qlg(idjd,:)
  write(6,*) 'qfg',qfg(idjd,:)
  write(6,*) 'qrg',qrg(idjd,:)
  write(6,*) 'qsng',qsng(idjd,:)
  write(6,*) 'qgrg',qgrg(idjd,:)
  write(6,*) 'qsubl',qsubl(idjd,:)
  write(6,*) 'rhoa',rhoa(idjd,:)
  write(6,*) 'rhos',rhos(idjd,:)
  write(6,*) 'fluxs ',fluxs(idjd,:)
  write(6,*) 'gam',gam(idjd,:)
  write(6,*) 'foutice',foutice(idjd,:)
  write(6,*) 'fthruice',fthruice(idjd,:)
  write(6,*) 'pqfsed',pqfsed(idjd,:)
  write(6,*) 'cfmelt',cfmelt(idjd,:)
  write(6,*) 'fluxmelt',fluxmelt(idjd,:)
  write(6,*) 'cifra,fluxsnow',cifra(idjd),fluxsnow(idjd)
end if  ! (diag.and.mydiag)

return
end subroutine newsnowrain

end module leoncld_mod
