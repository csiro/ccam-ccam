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

! The scheme has been modifed by MJT for max/rnd cloud overlap and to include prognostic rainfall, snow and
! graupel.  The snow and graupel components are based on Lin et al 1983, with modifications to the slope
! and intercept to be consistent with Rotstayn 97. There is also an optional prognostic cloud fraction option
! based on Tiedtke (see cloudmod.f90).

! ldr    = 0    Diagnosed cloud scheme (depreciated)
! ldr   /= 0    Prognostic cloud condensate (different ice fall speed options)
    
! ncloud = 0    Standard LDR cloud microphysics with water vapour, liquid cloud and ice cloud
! ncloud = 1    Use newer LDR autoconvection from Mk3.6
! ncloud = 2    Same as ncloud=1, but with prognostic rain
! ncloud = 3    Same as ncloud=2, but with prognostic graupel and snow
! ncloud = 4    Use prognostic cloud fraction based on Tiedtke from GFDL-CM3, but autoconversion from ncloud=0
! ncloud = 5    Same as ncloud=4, but convective sources are included in prognostic cloud fraction
   
!                            Water vapour (qg)
!
!   Cloud water (qlg,cfrac)                      Cloud ice (qfg,cfrac)
!
!   Rain (qrg,rfrac)                             Snow (qsg,sfrac)         Graupel (qgrg,gfrac)

! qg, qlg, qfg, qrg, qsg and qgrg are mixing ratios (g/g) and cfrac, rfrac, sfrac, gfrac are area cover
! fractions
    
module leoncld_mod
    
private
public leoncld

real, parameter :: maxlintime = 120. ! time-step for Lin et al 83 cloud microphysics

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

real, dimension(ifull,kl) :: qevap, qsubl, qauto, qcoll, qaccr, qaccf
real, dimension(ifull,kl) :: fluxr, fluxi, fluxs, fluxg, fluxm, pqfsed
real, dimension(ifull,kl) :: pfstayice, pfstayliq, pslopes, prscav
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
  qtot(1:ifull)     = qg(1:ifull,k) + qlg(1:ifull,k) + qrg(1:ifull,k) + qfg(1:ifull,k) &
                    + qsng(1:ifull,k) + qgrg(1:ifull,k)
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
precg(1:ifull) = 0. ! graupel

!     Set up convective cloud column
call convectivecloudfrac(clcon)
where ( ktsav(1:ifull)<kl-1 )
  ktop(1:ifull)   = ktsav(:)
  kbase(1:ifull)  = kbsav(:) + 1
  wcon(1:ifull)   = wlc
elsewhere
  wcon(1:ifull)   = 0.
end where

if ( nmaxpr==1 .and. mydiag ) then
  if ( ktau==1 ) then
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
if ( ncloud<=4 ) then

  ! diagnose cloud fraction (ncloud<=3) or prognostic strat. cloud but diagnostic conv. cloud (ncloud==4)
  if ( nmr>0 ) then
    ! Max/Rnd cloud overlap
    do k = 1,kl
      where ( clcon(1:ifull,k)>0. )
        !ccw = wcon(:)/rhoa(:,k)  !In-cloud l.w. mixing ratio
        qccon(1:ifull,k) = clcon(:,k)*wcon(:)/rhoa(:,k)
        qcl(1:ifull,k)   = max( qsatg(:,k), qg(1:ifull,k) )  ! qv in the convective fraction of the grid box
        qenv(1:ifull,k)  = max( 1.e-8, qg(1:ifull,k)-clcon(:,k)*qcl(:,k))/(1.-clcon(:,k) ) ! qv in the non-convective
                                                                                           ! fraction of the grid box
        qcl(1:ifull,k)   = (qg(1:ifull,k)-(1.-clcon(:,k))*qenv(1:ifull,k))/clcon(:,k)
        qlg(1:ifull,k)   = qlg(1:ifull,k)/(1.-clcon(:,k))    ! ql in the non-covective fraction of the grid box
        qfg(1:ifull,k)   = qfg(1:ifull,k)/(1.-clcon(:,k))    ! qf in the non-covective fraction of the grid box
        qrg(1:ifull,k)   = qrg(1:ifull,k)/(1.-clcon(:,k))    ! qr in the non-covective fraction of the grid box
        qsng(1:ifull,k)  = qsng(1:ifull,k)/(1.-clcon(:,k))   ! qs in the non-covective fraction of the grid box
        qgrg(1:ifull,k)  = qgrg(1:ifull,k)/(1.-clcon(:,k))   ! qg in the non-covective fraction of the grid box
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
        qcl(1:ifull,k)   = max( qsatg(1:ifull,k), qg(1:ifull,k) )  ! jlm
        qenv(1:ifull,k)  = max( 1.e-8, qg(1:ifull,k)-clcon(:,k)*qcl(1:ifull,k))/(1.-clcon(:,k) )
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
      
tenv(1:ifull,:) = t(1:ifull,:) ! Assume T is the same in and out of convective cloud

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
  t(1:ifull,k)  = clcon(:,k)*t(1:ifull,k) + (1.-clcon(:,k))*tenv(:,k)
  qg(1:ifull,k) = clcon(:,k)*qcl(:,k) + (1.-clcon(:,k))*qenv(:,k)
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
  fl(1:ifull)      = max(0., min(1., (t(1:ifull,k)-ticon)/(273.15-ticon)))
  qlrad(1:ifull,k) = qlg(1:ifull,k) + fl(:)*qccon(:,k)
  qfrad(1:ifull,k) = qfg(1:ifull,k) + (1.-fl(:))*qccon(:,k)
enddo

!     Calculate precipitation and related processes
call newsnowrain(dt,rhoa,dz,prf,cdso4,cfa,qca,t,qlg,qfg,qrg,qsng,qgrg,            &
                 precs,qg,cfrac,rfrac,sfrac,gfrac,ccov,preci,precg,qevap,qsubl,   &
                 qauto,qcoll,qaccr,qaccf,fluxr,fluxi,fluxs,fluxg,fluxm,           &
                 pfstayice,pfstayliq,pqfsed,pslopes,prscav)

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
    ppfprec(:,kl+1-k) = (fluxr(:,k+1)+fluxm(:,k))/dt                           !flux *entering* layer k
    ppfmelt(:,kl+1-k) = fluxm(:,k)/dt                                          !flux melting in layer k
    ppfsnow(:,kl+1-k) = (fluxi(:,k+1)+fluxs(:,k+1)+fluxg(:,k+1)-fluxm(:,k))/dt !flux *entering* layer k
  end do
  do k = 1,kl
    ppfevap(:,kl+1-k)    = qevap(:,k)*rhoa(:,k)*dz(:,k)/dt
    ppfsubl(:,kl+1-k)    = qsubl(:,k)*rhoa(:,k)*dz(:,k)/dt !flux sublimating or staying in k
    pplambs(:,kl+1-k)    = pslopes(:,k)
    ppmrate(:,kl+1-k)    = (qauto(:,k)+qcoll(:,k))/dt
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
cfrac(:,1:kl) = min(1., ccov(:,1:kl)+clcon(:,1:kl))

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

condx(1:ifull)  = condx(1:ifull) + precs(1:ifull)
conds(1:ifull)  = conds(1:ifull) + preci(1:ifull)
condg(1:ifull)  = condg(1:ifull) + precg(1:ifull)
precip(1:ifull) = precip(1:ifull) + precs(1:ifull)

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
real, dimension(ifull) :: pk_v, deles_v, qs_v, es_v
real, dimension(ifull) :: tk, fl, aprpr, bprpr, cice
real, dimension(ifull) :: qi0, fd, crate, qfdep

integer k, mg

real al, alf, cm0, decayfac
real deles, delq, dqsdt, hlrvap, pk, qc
real qs, rhoic
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
elsewhere ( ttg(1:ifull,:)>=tice .and. qfg(1:ifull,:)>1.e-12 )
  fice(1:ifull,:) = min(qfg(1:ifull,:)/(qfg(1:ifull,:)+qlg(1:ifull,:)), 1.)
elsewhere( ttg(1:ifull,:)>=tice )
  fice(1:ifull,:) = 0.
elsewhere
  fice(1:ifull,:) = 1.
end where
qcg(1:ifull,:)   = qlg(1:ifull,:) + qfg(1:ifull,:)
qcold(1:ifull,:) = qcg(:,:)
qfnew(1:ifull,:) = fice(:,:)*qcg(:,:)
ttg(1:ifull,:)   = ttg(1:ifull,:) + hlfcp*(qfnew(:,:)-qfg(1:ifull,:)) !Release L.H. of fusion
qfg(1:ifull,:)   = qfnew(:,:)
qlg(1:ifull,:)   = max(0., qcg(:,:)-qfg(1:ifull,:))

qtot(1:ifull,:) = qtg(1:ifull,:) + qcg(1:ifull,:)
tliq(1:ifull,:) = ttg(1:ifull,:) - hlcp*qcg(1:ifull,:) - hlfcp*qfg(1:ifull,:)  

if ( diag .and. mydiag ) then
  write(6,*) 'within newcloud'
  write(6,*) 'ttg ',ttg(idjd,:)
  write(6,*) 'qcold ',qcold(idjd,:)
  write(6,*) 'qcg ',qcg(idjd,:)
  write(6,*) 'qlg ',qlg(idjd,:)
  write(6,*) 'qfg ',qfg(idjd,:)
  write(6,*) 'fice ',fice(idjd,:)
end if

! Precompute the array of critical relative humidities 
if ( nclddia==-3 ) then
  do k = 1,kl
    where ( land(1:ifull) )
      rcrit(:,k)=max( rcrit_l, (1.-16.*(1.-sig(k))**3) )
    elsewhere
      rcrit(:,k)=max( rcrit_s, (1.-16.*(1.-sig(k))**3) )
    end where
  enddo
else if ( nclddia<0 ) then
  do k = 1,kl
    where ( land(1:ifull) )
      rcrit(:,k)=max( rcrit_l, (1.-4.*(1.-sig(k))**2) )
    elsewhere
      rcrit(:,k)=max( rcrit_s, (1.-4.*(1.-sig(k))**2) )
    end where
  enddo
else if ( nclddia==1 ) then
  do k = 1,kl
    where ( land(1:ifull) )
      rcrit(:,k)=max( rcrit_l, sig(k)**3 )
    elsewhere
      rcrit(:,k)=max( rcrit_s, sig(k)**3 )
    end where
  enddo
else if ( nclddia==2 ) then
  do k = 1,kl
    where ( land(1:ifull) )
      rcrit(:,k)=rcrit_l
    elsewhere
      rcrit(:,k)=rcrit_s
    end where
  enddo
else if ( nclddia==3 ) then
  do k = 1,kl
    where ( land(1:ifull) )
      rcrit(:,k)=max( rcrit_l, sig(k)**2 )          ! .75 for R21 Mk2
    elsewhere
      rcrit(:,k)=max( rcrit_s, sig(k)**2 )          ! .85 for R21 Mk2
    end where
  enddo
else if ( nclddia==4 ) then
  do k = 1,kl
    where ( land(1:ifull) )
      rcrit(:,k)=max( rcrit_l, sig(k) )             ! .75 for test Mk2/3
    elsewhere
      rcrit(:,k)=max( rcrit_s, sig(k) )             ! .9  for test Mk2/3
    end where
  enddo
else if ( nclddia==5 ) then  ! default till May 08
  do k = 1,kl
    where ( land(1:ifull) )
      rcrit(:,k)=max( rcrit_l, min(.99,sig(k)) )    ! .75 for same as T63
    elsewhere
      rcrit(:,k)=max( rcrit_s, min(.99,sig(k)) )    ! .85 for same as T63
    end where
  enddo
else if ( nclddia==6 ) then
  do k = 1,kl
    where ( land(1:ifull) )
      rcrit(:,k)=max(rcrit_l*(1.-.15*sig(k)),sig(k)**4)
    elsewhere
      rcrit(:,k)=max(rcrit_s*(1.-.15*sig(k)),sig(k)**4)
    end where
  enddo
else if ( nclddia==7 ) then
  do k = 1,kl
    where ( land(1:ifull) )
      rcrit(:,k)=max(rcrit_l*(1.-.2*sig(k)),sig(k)**4)
    elsewhere
      rcrit(:,k)=max(rcrit_s*(1.-.2*sig(k)),sig(k)**4)
    end where
  enddo
else if ( nclddia>7 ) then  ! e.g. 12    JLM
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


if ( ncloud<=3 ) then
  ! usual diagnostic cloud fraction
      
  ! Calculate cloudy fraction of grid box (cfrac) and gridbox-mean cloud water
  ! using the triangular PDF of Smith (1990)

  do k = 1,kl
    do mg = 1,ifull
      hlrvap = (hl+fice(mg,k)*hlf)/rvap

      ! Calculate qs and gam=(L/cp)*dqsdt,  at temperature tliq
      pk = 100.0*prf(mg,k)
      qsi(mg,k) = qsati(pk,tliq(mg,k))           !Ice value
      deles = esdiffx(tliq(mg,k))                ! MJT suggestion
      qsl(mg,k) = qsi(mg,k) + epsil*deles/pk     !qs over liquid
      qsw(mg,k) = fice(mg,k)*qsi(mg,k) + (1.-fice(mg,k))*qsl(mg,k) !Weighted qs at temperature Tliq
      qs = qsw(mg,k)
      dqsdt = qs*hlrvap/tliq(mg,k)**2
      !qvc(mg,k) = qs !Vapour mixing ratio in cloud

      al = 1./(1.+(hlcp+fice(mg,k)*hlfcp)*dqsdt)    !Smith's notation
      qc = qtot(mg,k) - qs

      delq = (1.-rcrit(mg,k))*qs      !UKMO style (equivalent to above)
      cfrac(mg,k) = 1.
      qcg(mg,k) = al*qc
      if ( qc<delq ) then
        cfrac(mg,k) = max( 1.e-6, 1.-.5*((qc-delq)/delq)**2 )       ! for roundoff
        qcg(mg,k) = max( 1.e-8, al*(qc-(qc-delq)**3/(6.*delq**2)) ) ! for roundoff
      end if
      if ( qc<=0. ) then
        cfrac(mg,k) = max( 1.e-6, .5*((qc+delq)/delq)**2 )     ! for roundoff
        qcg(mg,k) = max( 1.e-8, al*(qc+delq)**3/(6.*delq**2) ) ! for roundoff
      end if
      if ( qc<=-delq ) then
        cfrac(mg,k) = 0.
        qcg(mg,k) = 0.
      end if

      ! Calculate the cloud fraction (cfa) in which ql exceeds qcrit, and
      ! the corresponding gridbox-mean cloud water mixing ratio qca. 
      ! This (qca) is the cloud-water mixing ratio inside cfa times cfa.
      ! The new variable qc2 is like qc above, but is used for integration limits
      ! only, not the integrand

      if ( cfrac(mg,k)>0. ) then
        qcic = qcg(mg,k)/cfrac(mg,k) !Mean in cloud value

        ! Following few lines are for Yangang Liu's new scheme (2004: JAS, GRL)
        ! Need to do first-order estimate of qcrit using mean in-cloud qc (qcic)

        Wliq = max( 1.e-10, 1000. * qcic * rhoa(mg,k)) !g/m3
        R6c = 4.09e-4 * ( 1.15e23*1.e-6*cdrop(mg,k) / Wliq**2 ) ** (1./6.)
        eps = 1. - 0.7 * exp(-0.003e-6*cdrop(mg,k)) !mid range
        beta6 = ((1.+3.*eps**2)*(1.+4.*eps**2)*(1.+5.*eps**2) / ((1.+eps**2)*(1.+2.*eps**2)) )**(1./6.)
        R3c = 1.e-6*R6c/beta6 !in metres
        qcrit = (4.*pi/3.)*rhow*R3c**3*cdrop(mg,k)/rhoa(mg,k) !New qcrit

        qc2 = qtot(mg,k) - qs - qcrit/al
        cfa(mg,k) = 1.
        qca(mg,k) = al*qc
        if ( qc2<delq ) then
          cfa(mg,k) = 1. - 0.5*((qc2-delq)/delq)**2
          qto = (qtot(mg,k)-delq+2.*(qs+qcrit/al))/3.
          qca(mg,k) = al*(qtot(mg,k) - qto + cfa(mg,k)*(qto-qs))
        end if
        if ( qc2<=0. ) then
          cfa(mg,k) = 0.5*((qc2+delq)/delq)**2
          qca(mg,k) = cfa(mg,k)*(al/3.)*(2.*qcrit/al + qc+delq)
        end if
        if ( qc2<=-delq ) then
          cfa(mg,k) = 0.
          qca(mg,k) = 0.
        end if
      else
        cfa(mg,k) = 0.
        qca(mg,k) = 0.
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
  decayfac = exp ( -tdt/7200. )             ! Try 2 hrs
  !decayfac = 0.                            ! Instant adjustment (old scheme)
  where( ttg(1:ifull,:)>=Tice )
    qfg(1:ifull,:) = fice*qcg(:,:)
    qlg(1:ifull,:) = qcg(:,:) - qfg(1:ifull,:)
  elsewhere                                    ! Cirrus T range
    qfg(1:ifull,:) = qcold(:,:)*decayfac + qcg(:,:)*(1.-decayfac)
    qlg(1:ifull,:) = 0.
    qcg(1:ifull,:) = qfg(1:ifull,:)
  end where
  
else
  
  ! Tiedtke prognostic cloud fraction model
  ! MJT notes - we use ttg instead of tliq
  do k = 1,kl
    pk_v = 100.*prf(1:ifull,k)
    qsi(1:ifull,k) = qsati(pk_v(1:ifull),ttg(1:ifull,k)) ! Ice value
    deles_v = esdiffx(ttg(1:ifull,k))
    qsl(1:ifull,k) = qsi(1:ifull,k) + epsil*deles_v/pk_v ! Liquid value
  end do
  qsw(:,:) = fice*qsi + (1.-fice)*qsl        ! Weighted qs at temperature Tliq
  call progcloud(cfrac,qcg,qtot,prf,rhoa,fice,qsw,ttg,rcrit)
        
  ! Use 'old' autoconversion with prognostic cloud
  cfa(:,:) = 0.
  qca(:,:) = 0.

  where( ttg(1:ifull,:)>=Tice )
    qfg(1:ifull,:) = fice(:,:)*qcg(:,:)
    qlg(1:ifull,:) = qcg(:,:) - qfg(1:ifull,:)
  elsewhere
    qfg(1:ifull,:) = qcg(:,:)
    qlg(1:ifull,:) = 0.
    qcg(1:ifull,:) = qfg(1:ifull,:)
  end where
  
end if ! ncloud<4 ..else..


! Do the vapour deposition calculation in mixed-phase clouds:
! Calculate deposition on cloud ice, assuming es(T) is the weighted value of the 
! liquid and ice values.
alf = 1./3.
rhoic = 700.
cm0 = 1.e-12 !Initial crystal mass
do k = 1,kl
  where ( cfrac(1:ifull,k)>0. )
    Tk(:) = tliq(:,k) + hlcp*(qlg(1:ifull,k)+qfg(1:ifull,k))/cfrac(1:ifull,k) !T in liq cloud
    !fl(:) = qlg(1:ifull,k)/max(qfg(1:ifull,k)+qlg(1:ifull,k),1.e-30)
  end where
  where ( cfrac(1:ifull,k)>0. .and. Tk(:)<tfrz .and. qlg(1:ifull,k)>1.e-8 )
    pk_v(:)    = 100.*prf(:,k)
    qs_v(:)    = qsati(pk_v(:),Tk(:))
    es_v(:)    = qs_v(:)*pk_v(:)/0.622 !ice value
    Aprpr(:)   = hl/(rKa*Tk(:))*(hls/(rvap*Tk(:))-1.)
    Bprpr(:)   = rvap*Tk(:)/((Dva/pk_v(:))*es_v(:))
    deles_v(:) = (1.-fice(:,k))*esdiffx(Tk(:))
    Cice(:)    = 1.e3*exp(12.96*deles_v(:)/es_v(:) - 0.639) !Meyers et al 1992
    qi0(:)     = cm0*Cice(:)/rhoa(:,k) !Initial ice mixing ratio
    ! Next 2 lines are for assumption of fully mixed ql and qf (also a line further down).
    qi0(:)     = max(qi0(:), qfg(1:ifull,k)/cfrac(1:ifull,k)) !Assume all qf and ql are mixed
    fd(:)      = 1.       !Fraction of cloud in which deposition occurs
    !fd(:)     = fl(:)   !Or, use option of adjacent ql,qf
    Crate(:)   = 7.8*((Cice(:)/rhoa(:,k))**2/rhoic)**(1./3.)*deles_v(:)/((Aprpr(:)+Bprpr(:))*es_v(:))
    qfdep(:)   = fd(:)*cfrac(1:ifull,k)*sqrt(((2./3.)*Crate(:)*tdt+qi0(:)**(2./3.))**3)
    ! Also need this line for fully-mixed option...
    qfdep(:)   = qfdep(:) - qfg(1:ifull,k)
    qfdep(:)   = min(qfdep(:), qlg(1:ifull,k))
    qlg(1:ifull,k) = qlg(1:ifull,k) - qfdep(:)
    qfg(1:ifull,k) = qfg(1:ifull,k) + qfdep(:)
  end where
  fice(1:ifull,k) = qfg(1:ifull,k)/max(qfg(1:ifull,k)+qlg(1:ifull,k),1.e-30)
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
!      qaccr - accretion by snow of cloud liquid water (kg/kg)
!
!**************************************************************************

subroutine newsnowrain(tdt_in,rhoa,dz,prf,cdrop,cfa,qca,ttg,qlg,qfg,qrg,qsng,qgrg,precs,qtg,cfrac,cfrainfall,    &
                       cfsnowfall,cfgraupelfall,ccov,preci,precg,qevap,qsubl,qauto,qcoll,qaccr,qaccf,fluxr,      &
                       fluxi,fluxs,fluxg,fluxm,pfstayice,pfstayliq,pqfsed,pslopes,prscav)

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
real, dimension(ifull,kl), intent(inout) :: cfrac
real, dimension(ifull+iextra,kl), intent(inout) :: cfrainfall
real, dimension(ifull+iextra,kl), intent(inout) :: cfsnowfall
real, dimension(ifull+iextra,kl), intent(inout) :: cfgraupelfall
real, dimension(ifull,kl), intent(in) :: ccov
real, dimension(ifull,kl), intent(out) :: qevap
real, dimension(ifull,kl), intent(out) :: qsubl
real, dimension(ifull,kl), intent(out) :: qauto
real, dimension(ifull,kl), intent(out) :: qcoll
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
real, dimension(ifull,kl), intent(out) :: fluxm

! Local work arrays and variables
real, dimension(ifull,kl-1) :: fthruliq,foutliq,fthruice,foutice
real, dimension(ifull,kl-1) :: fthrusnow,foutsnow,fthrugraupel,foutgraupel
real, dimension(ifull,kl-1) :: rhor,gam
real, dimension(ifull,kl) :: cfrain,fluxauto,vl2
real, dimension(ifull,kl) :: vi2,rhoi
real, dimension(ifull,kl) :: cfsnow,fluxprecipitation,vs2,rhos
real, dimension(ifull,kl) :: cfgraupel,fluxautograupel,vg2,rhog
real, dimension(ifull,kl) :: clfr,cifr,qsatg
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
real, dimension(ifull) :: n0s, rica
real, dimension(ifull) :: cftmp, xwgt, cfmelt, fluxmelt
real, dimension(ifull) :: rhodum_r
real, dimension(ifull) :: slopes_i, slopes_s, slopes_g, slopes_r
real, dimension(ifull) :: denfac, esi, qsl, apr, bpr, cev
real, dimension(ifull) :: dqsdt, bl, satevap

real, parameter :: n0r = 8.e6        ! intercept for rain
real, parameter :: n0g = 4.e6        ! intercept for graupel
real, parameter :: rho_r = 1.0e3     ! rain density
real, parameter :: rho_s = 0.1e3     ! snow density
real, parameter :: rho_g = 0.4e3     ! grauple density
real, parameter :: qr0_crt = 2.e-4   ! rain -> snow or graupel density threshold
real, parameter :: qi0_crt = 8.e-5   ! ice -> snow density threshold
real, parameter :: qs0_crt = 6.e-3   ! snow -> graupel density threshold
real, parameter :: c_piacr = 0.1     ! accretion rate of rain -> ice
real, parameter :: c_psaut = 1.e-3   ! autoconversion rate of ice -> snow
real, parameter :: c_pgacs = 1.e-3   ! snow -> graupel "accretion" eff
real, parameter :: sfcrho = 1.2      ! reference density rho_0
real, parameter :: vdifu = 2.11e-5
real, parameter :: tcond = 2.36e-2
real, parameter :: visk = 1.259e-5
real, parameter :: gam263 = 1.456943 ! gamma function for 2.63
real, parameter :: gam275 = 1.608355 ! gamma function for 2.75
real, parameter :: gam325 = 2.54925  ! gamma function for 3.25
real, parameter :: gam350 = 3.323363 ! gamma function for 3.5
real, parameter :: gam380 = 4.694155 ! gamma function for 3.8
real, parameter :: alin = 842.
real, parameter :: clin = 4.8
real, parameter :: gcon = 44.628 ! = 40.74*sqrt(sfcrho)
real, parameter :: tau_s = 90.   ! (sec) snow melt
real, parameter :: tau_g = 180.  ! (sec) graupel melt

integer k, mg, n, njumps

real crate,frb,qcic,qcrit,ql1,ql2,R6c,R3c,beta6,eps
real selfcoll,Wliq,cfla,dqla,qla
real scm3, tdt

scm3 = (visk/vdifu)**(1./3.)

do k = 1,kl
  fluxr(1:ifull,k)             = 0.
  fluxi(1:ifull,k)             = 0.
  fluxs(1:ifull,k)             = 0.
  fluxg(1:ifull,k)             = 0.
  fluxm(1:ifull,k)             = 0.  
  fluxauto(1:ifull,k)          = 0.
  fluxprecipitation(1:ifull,k) = 0.
  fluxautograupel(1:ifull,k)   = 0.
  qevap(1:ifull,k)             = 0.
  qauto(1:ifull,k)             = 0.
  qcoll(1:ifull,k)             = 0.
  qsubl(1:ifull,k)             = 0.
  qaccr(1:ifull,k)             = 0.
  qaccf(1:ifull,k)             = 0.
  pqfsed(1:ifull,k)            = 0.
  prscav(1:ifull,k)            = 0.  
  pfstayice(1:ifull,k)         = 0.  
  pfstayliq(1:ifull,k)         = 0. 
  pslopes(1:ifull,k)           = 0.
  pk(1:ifull)                  = 100.*prf(1:ifull,k)
  qsatg(1:ifull,k)             = qsati(pk(1:ifull),ttg(1:ifull,k))
  cifr(1:ifull,k)              = cfrac(1:ifull,k)*qfg(1:ifull,k)/max(qlg(1:ifull,k)+qfg(1:ifull,k),1.E-30)
  clfr(1:ifull,k)              = cfrac(1:ifull,k)*qlg(1:ifull,k)/max(qlg(1:ifull,k)+qfg(1:ifull,k),1.E-30)
  cfrain(1:ifull,k)            = 0.
  cfsnow(1:ifull,k)            = 0.
  cfgraupel(1:ifull,k)         = 0.
end do


! Use full timestep for autoconversion
!njumps = 1
tdt = tdt_in


!**************** Cut here to insert new auto scheme ********************            
if ( ncloud>0 .and. ncloud<=3 ) then

  ! Using new (subgrid) autoconv scheme... 
  do k = kl-1,1,-1
    rhodz(1:ifull) = rhoa(1:ifull,k)*dz(1:ifull,k)
    do mg = 1,ifull
      if ( clfr(mg,k)>0. ) then
        if ( cfa(mg,k)>0. ) then
          cfla = cfa(mg,k)*clfr(mg,k)/(clfr(mg,k)+cifr(mg,k))
          qla = qca(mg,k)/cfa(mg,k)
          ! Following few lines are for Yangang Liu's new scheme (2004: JAS, GRL)
          Wliq = max( 1.e-10, 1000. * qla * rhoa(mg,k) ) !g/m3
          R6c = 4.09e-4 * ( 1.15e23*1.e-6*cdrop(mg,k) / Wliq**2 )**(1./6.)
          eps = 1. - 0.7 * exp(-0.003e-6*cdrop(mg,k)) !mid range
          beta6 = ((1.+3.*eps**2)*(1.+4.*eps**2)*(1.+5.*eps**2) / ((1.+eps**2)*(1.+2.*eps**2)) )**(1./6.)
          R3c = 1.e-6*R6c/beta6 !in metres
          qcrit = (4.*pi/3.)*rhow*R3c**3*Cdrop(mg,k)/rhoa(mg,k) !New qcrit
          if ( qla<=qcrit ) then
            ql2 = qla
          else
            ! Following is Liu & Daum (JAS, 2004)
            Crate = 1.9e17*(0.75*rhoa(mg,k)/(pi*rhow))**2*beta6**6/cdrop(mg,k)
            ql1 = qla/sqrt(1.+2.*crate*qla**2*tdt)
            ql1 = max( ql1, qcrit ) !Intermediate qlg after auto
            Frb = dz(mg,k)*rhoa(mg,k)*(qla-ql1)/tdt
            cdt(mg) = tdt*0.5*Ecol*0.24*pow75(Frb)
            selfcoll = min( ql1, ql1*cdt(mg) )
            ql2 = ql1 - selfcoll
          end if
          dqla = cfla*(qla-ql2)
          ql(mg) = max( 1.e-10, qlg(mg,k)-dqla )
          dql(mg) = max( qlg(mg,k)-ql(mg), 0. )
        else
          cfla = 0.
          dqla = 0.
          ql(mg) = max( 1.e-10, qlg(mg,k) )
          dql(mg) = 0.
        end if
        cfrain(mg,k) = cfla*dql(mg)/qlg(mg,k)
        qauto(mg,k) = qauto(mg,k) + dql(mg)
        qlg(mg,k) = qlg(mg,k) - dql(mg)
        fluxauto(mg,k) = dql(mg)*rhodz(mg)
      end if
    end do
  end do

! Or, using old autoconv scheme... also used by prognostic cloud scheme
else

  do k = kl-1,1,-1
    rhodz(1:ifull) = rhoa(1:ifull,k)*dz(1:ifull,k)
    do mg = 1,ifull
      if ( clfr(mg,k)>0. ) then
        qcrit = (4.*pi/3.)*rhow*Rcm**3*cdrop(mg,k)/rhoa(mg,k)
        qcic = qlg(mg,k)/clfr(mg,k) !In cloud value
        if ( qcic<qcrit ) then
          ql(mg) = qlg(mg,k)
          dql(mg) = 0.
        else
          Crate = Aurate*rhoa(mg,k)*(rhoa(mg,k)/(cdrop(mg,k)*rhow))**(1./3.)
          ql1 = 1./pow75(qcic**(-4./3.)+(4./3.)*Crate*tdt)
          ql1 = max( ql1, qcrit ) !Intermediate qlg after auto
          Frb = dz(mg,k)*rhoa(mg,k)*(qcic-ql1)/tdt
          cdt(mg) = tdt*0.5*Ecol*0.24*pow75(Frb) !old
          selfcoll = min( ql1, ql1*cdt(mg) )
          ql2 = ql1 - selfcoll
          ql(mg) = clfr(mg,k)*ql2
          dql(mg) = max( qlg(mg,k)-ql(mg), 0. )
        end if
        cfrain(mg,k) = clfr(mg,k)*dql(mg)/qlg(mg,k)        
        qauto(mg,k) = qauto(mg,k) + dql(mg)
        qlg(mg,k) = qlg(mg,k) - dql(mg)
        fluxauto(mg,k) = dql(mg)*rhodz(mg)
      end if
    end do
  end do

end if ! ( ncloud>0 .and. ncloud<=3 ) ..else..


! calculate rate of precipitation of frozen cloud water to snow
if ( ncloud>=3 ) then

  do k = 1,kl

    rhodz(1:ifull) = rhoa(1:ifull,k)*dz(1:ifull,k)
      
    ! autoconversion of ice to snow (from Lin et al 1983)
    ! Threshold from WSM6 scheme, Hong et al 2004, Eq(13) : qi0_crt ~8.e-5
    where ( ttg(1:ifull,k)<tfrz .and. qfg(1:ifull,k)*rhoa(1:ifull,k)>qi0_crt )
      qf(:)  = max( qfg(1:ifull,k)-qi0_crt/rhoa(:,k), 0. )
      cdt(:) = tdt*c_psaut*exp(0.025*(ttg(1:ifull,k)-tfrz))
      dqf(:) = max( min( qfg(1:ifull,k), qf(:)*cdt(:) ), 0. )
      cfsnow(1:ifull,k)            = cifr(1:ifull,k)*dqf(:)/qfg(1:ifull,k)
      cfrac(1:ifull,k)             = cfrac(1:ifull,k) - cfsnow(1:ifull,k)
      qfg(1:ifull,k)               = qfg(1:ifull,k) - dqf(:)
      fluxprecipitation(1:ifull,k) = dqf(:)*rhodz(:)
    end where
    
    ! autoconversion of snow to graupel (from Lin et al 1983)
    where ( ttg(1:ifull,k)<tfrz .and. qsng(1:ifull,k)*rhoa(1:ifull,k)>qs0_crt )
      qf(:)  = max( qsng(1:ifull,k)-qs0_crt/rhoa(:,k), 0. )
      cdt(:) = tdt*1.e-3*exp(0.09*(ttg(1:ifull,k)-tfrz))
      dqf(:) = max( min( qsng(1:ifull,k), qf(:)*cdt(:) ), 0.) 
      cfgraupel(1:ifull,k)       = cfsnowfall(1:ifull,k)*dqf(:)/qsng(1:ifull,k)
      cfsnowfall(1:ifull,k)      = cfsnowfall(1:ifull,k) - cfgraupel(1:ifull,k)
      qsng(1:ifull,k)            = qsng(1:ifull,k) - dqf(:)
      fluxautograupel(1:ifull,k) = dqf(:)*rhodz(:)
    end where

  end do

end if ! ( ncloud>=3 )


! update cloud liquid and frozen water fractions
cifr(1:ifull,1:kl) = cfrac(1:ifull,1:kl)*qfg(1:ifull,1:kl)/max(qlg(1:ifull,1:kl)+qfg(1:ifull,1:kl),1.e-30 )
clfr(1:ifull,1:kl) = max( cfrac(1:ifull,1:kl)-cifr(1:ifull,1:kl), 0. )
rhoi(1:ifull,1:kl) = qfg(1:ifull,1:kl)*rhoa(1:ifull,1:kl)
vi2(1:ifull,kl)    = 0.1 ! Assume no cloud at top level

! combine autoconversion and prognostic rain
! max overlap autoconversion and rainfall from previous time step
cfrain(1:ifull,1:kl-1) = max( cfrain(1:ifull,1:kl-1), cfrainfall(1:ifull,1:kl-1) ) 
rhor(1:ifull,1:kl-1)   = qrg(1:ifull,1:kl-1)*rhoa(1:ifull,1:kl-1)
vl2(1:ifull,kl)        = 0.

! Set up snow fields
cfsnow(1:ifull,:) = max( cfsnow(1:ifull,:), cfsnowfall(1:ifull,:) ) 
rhos(1:ifull,:)   = qsng(1:ifull,:)*rhoa(1:ifull,:)
vs2(1:ifull,kl)   = 0.1

! Set up graupel fields
cfgraupel(1:ifull,:) = max( cfgraupel(1:ifull,:), cfgraupelfall(1:ifull,:) ) 
rhog(1:ifull,:)      = qgrg(1:ifull,:)*rhoa(1:ifull,:)
vg2(1:ifull,kl)      = 0.1


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
endif  ! (diag.and.mydiag)


! Use sub time-step if required
if ( ncloud >= 3 ) then
  njumps = int(tdt_in/(maxlintime+0.01)) + 1
  tdt    = tdt_in/real(njumps)
else
  njumps = 1
  tdt = tdt_in
end if

do n = 1,njumps


  cfmelt(1:ifull) = 0.

  fluxgraupel(1:ifull)   = 0.
  mxclfrgraupel(1:ifull) = 0. ! max overlap graupel fraction
  rdclfrgraupel(1:ifull) = 0. ! rnd overlap graupel fraction
  cgfra(1:ifull)         = 0.

  fluxsnow(1:ifull)   = 0.
  mxclfrsnow(1:ifull) = 0. ! max overlap snow fraction
  rdclfrsnow(1:ifull) = 0. ! rnd overlap snow fraction
  csfra(1:ifull)      = 0.

  fluxice(1:ifull)   = 0.
  mxclfrice(1:ifull) = 0. ! max overlap ice fraction
  rdclfrice(1:ifull) = 0. ! rnd overlap ice fraction
  cifra(1:ifull)     = 0.
  rica(1:ifull)      = 0. ! backward compatibility for ncloud<=2

  fluxrain(1:ifull)  = 0.
  mxclfrliq(1:ifull) = 0. ! max overlap rain fraction
  rdclfrliq(1:ifull) = 0. ! rnd overlap rain fraction
  clfra(1:ifull)     = 1.e-6


  ! Now work down through the levels...
  do k = kl-1,1,-1
  
    ! misc fields
    pk(:)              = 100.*prf(:,k)
    rhodz(:)           = rhoa(:,k)*dz(:,k)  
    denfac(1:ifull)    = sqrt(sfcrho/rhoa(:,k))
    n0s(1:ifull)       = 2.e6*exp(-0.12*max(ttg(1:ifull,k)-tfrz,-200.)) ! intercept
    fluxmelt(1:ifull)  = 0.
  
    ! default fall velocities
    vg2(1:ifull,k) = vg2(1:ifull,k+1)
    vs2(1:ifull,k) = vs2(1:ifull,k+1)
    vi2(1:ifull,k) = vi2(1:ifull,k+1)
    vl2(1:ifull,k) = vl2(1:ifull,k+1)

    ! default slopes
    slopes_g(1:ifull) = ( max( fluxgraupel(:), 0. )/dz(:,k)/(pi*n0g*rho_g))**0.25 ! from Lin et al 83
    slopes_s(1:ifull) = ( max( fluxsnow(:), 0. )/dz(:,k)/(pi*rho_s*n0s(:)))**0.25 ! from Lin et al 83
    slopes_i(1:ifull)  = 1.6e3*10**(-0.023*(ttg(1:ifull,k)-tfrz))     ! from HDC 04
    slopes_r(1:ifull) = (( max( fluxrain(:), 0. )/max( clfra(:),1.e-15 )/tdt)**0.22)/714. ! from LDR97
    
    pslopes(1:ifull,k) = pslopes(1:ifull,k) + slopes_i(1:ifull)*tdt/tdt_in  
    
    if ( ncloud>=3 ) then
      
  
      ! Graupel ---------------------------------------------------------------------------
      alphaf(1:ifull)   = hls*qsatg(1:ifull,k)/(rvap*ttg(1:ifull,k)**2)
      gam(1:ifull,k)    = hlscp*alphaf(:) !(L/cp)*dqsdt (HBG notation)
      sublflux(1:ifull) = 0.
      caccr_g(1:ifull)  = 0.
      caccf_g(1:ifull)  = 0.

      ! The following flag detects max/random overlap clouds
      ! that are separated by a clear layer
      where ( cfgraupel(1:ifull,k)<1.e-10 .or. nmr==0 )
        ! combine max overlap from last cloud with net random overlap
        rdclfrgraupel(1:ifull) = rdclfrgraupel(:) + mxclfrgraupel(:) - rdclfrgraupel(:)*mxclfrgraupel(:)
        mxclfrgraupel(1:ifull) = 0.
      end where
    
      fluxgraupel(:) = fluxgraupel(:) + fluxautograupel(:,k)*tdt/tdt_in
 
      ! graupel fall speed (from Lin et al 1983 - see GFDL AM3)
      where ( cfgraupel(1:ifull,k)>=1.e-10 )
        vg2(1:ifull,k) = max( 0.1, 87.2382675*(max( rhog(:,k), 0. )/cfgraupel(:,k)/5026548245.74367)**0.125 )
      elsewhere
        vg2(1:ifull,k) = vg2(1:ifull,k+1)
      end where
    
      ! Set up the parameters for the flux-divergence calculation
      alph(:)           = tdt*vg2(:,k)/dz(:,k)
      foutgraupel(:,k)  = 1. - exp(-alph(:))             !analytical
      fthrugraupel(:,k) = 1. - foutgraupel(:,k)/alph(:)  !analytical

      ! Melt falling graupel (based on Lin et al 83)
!      qgr(1:ifull)      = fluxgraupel(1:ifull)/rhodz(1:ifull)
!      slopes_g(1:ifull) = ( max( fluxgraupel(:), 0. )/dz(:,k)/(pi*n0g*rho_g))**0.25
!      where ( ttg(1:ifull,k)>tfrz .and. qgr(1:ifull)>1.e-10 )
!        cdt(1:ifull)           = tdt*2.*pi*n0g/hlf*(tcond*(ttg(1:ifull,k)-tfrz)/rhoa(:,k)-vdifu*hl*(qsatg(:,k)-qtg(1:ifull,k)))  &
!                                    *(0.78*slopes_g(:)**2+0.31*scm3*gam275*sqrt(gcon/visk)*slopes_g(:)**2.75*sqrt(denfac(:)))
!        qif(1:ifull)           = max( min( qgr(:), qgr(:)*cdt(:) ), 0. ) !Mixing ratio of graupel
!        dttg(1:ifull)          = -hlfcp*qif(:)
!        ttg(1:ifull,k)         = ttg(1:ifull,k) + dttg(:)
!        qsatg(1:ifull,k)       = qsatg(1:ifull,k) + gam(:,k)*dttg(:)/hlscp
!        fluxmelt(1:ifull)      = fluxmelt(:)    + qif(:)*rhodz(:)
!        fluxgraupel(1:ifull)   = fluxgraupel(:) - qif(:)*rhodz(:)
!        rdclfrgraupel(1:ifull) = rdclfrgraupel(:)*(1.-qif(:)/qgr(:))
!        mxclfrgraupel(1:ifull) = mxclfrgraupel(:)*(1.-qif(:)/qgr(:))
!        cftmp(1:ifull)         = mxclfrgraupel(:) + rdclfrgraupel(:) - mxclfrgraupel(:)*rdclfrgraupel(:)
!        cfmelt(1:ifull)        = max( cfmelt(:), max( cgfra(:)-cftmp(:), 0. ) )
!        cgfra(1:ifull)         = cftmp(:)
!      end where
    
      ! Melt falling graupel if > 2 deg C  (based on Lin et al 83)
      qgr(1:ifull) = fluxgraupel(1:ifull)/rhodz(1:ifull)
      where ( ttg(1:ifull,k)>tfrz+2. .and. qgr(1:ifull)>1.e-10 )
        cdt(1:ifull)           = min( 1., tdt/tau_g )*(ttg(1:ifull,k)-tfrz-2.)/hlfcp
        qif(1:ifull)           = min( qgr(:), cdt(:) ) !Mixing ratio of graupel
        dttg(1:ifull)          = -hlfcp*qif(:)
        ttg(1:ifull,k)         = ttg(1:ifull,k) + dttg(:)
        qsatg(1:ifull,k)       = qsatg(1:ifull,k) + gam(:,k)*dttg(:)/hlscp
        fluxmelt(1:ifull)      = fluxmelt(:)    + qif(:)*rhodz(:)
        fluxgraupel(1:ifull)   = fluxgraupel(:) - qif(:)*rhodz(:)
        rdclfrgraupel(1:ifull) = rdclfrgraupel(:)*(1.-qif(:)/qgr(:))
        mxclfrgraupel(1:ifull) = mxclfrgraupel(:)*(1.-qif(:)/qgr(:))
        cftmp(1:ifull)         = mxclfrgraupel(:) + rdclfrgraupel(:) - mxclfrgraupel(:)*rdclfrgraupel(:)
        cfmelt(1:ifull)        = max( cfmelt(:), max( cgfra(:)-cftmp(:), 0. ) )
        cgfra(1:ifull)         = cftmp(:)
      end where

      ! Sublimation of graupel is neglected in the UM and ACCESS 1.3.
      ! (Currently treated the same as LDR97 ice sublimation)
      fsclr_g(:)        = max( (1.-cifr(:,k)-clfr(:,k))*fluxgraupel(:), 0. )
      slopes_g(1:ifull) = ( max( fluxgraupel(:), 0. )/dz(:,k)/(pi*n0g*rho_g))**0.25
      where ( fluxgraupel(:)>0. .and. qtg(1:ifull,k)<qsatg(1:ifull,k) ) ! sublime graupel
        cdt(1:ifull)         = 2.*pi*vdifu*tcond*rvap*n0g*ttg(1:ifull,k)**2                                            &
                              *(0.78*slopes_g(:)**2+0.31*scm3*gam275*sqrt(gcon/visk)*slopes_g(:)**2.75*sqrt(denfac(:))) &
                              /(tcond*rvap*ttg(1:ifull,k)**2+hls**2*vdifu*qsatg(:,k)*rhoa(1:ifull,k))
        dqs(1:ifull)         = tdt*cdt(:)*(qsatg(:,k)-qtg(1:ifull,k))
        dqs(1:ifull)         = min( dqs(:), (qsatg(1:ifull,k)-qtg(1:ifull,k))/(1.+gam(:,k)) ) !Don't supersat.
        sublflux(1:ifull)    = min( dqs(:)*rhodz(:), fsclr_g(:) )
        fluxgraupel(1:ifull) = fluxgraupel(:) - sublflux(:)
        fsclr_g(1:ifull)     = fsclr_g(:)     - sublflux(:)
        dqs(1:ifull)         = sublflux(:)/rhodz(:)
        qsubl(1:ifull,k)     = qsubl(:,k)     + dqs(:)
        qtg(1:ifull,k)       = qtg(1:ifull,k) + dqs(:)
        dttg(1:ifull)        = -hlscp*dqs(:)
        ttg(1:ifull,k)       = ttg(1:ifull,k) + dttg(:)
        qsatg(1:ifull,k)     = qsatg(1:ifull,k) + gam(:,k)*dttg(:)/hlscp
      end where
    
      ! Save flux for the wet deposition scheme.
      pfstayice(:,k) = pfstayice(:,k) + fluxgraupel(:)*(1.-fthrugraupel(:,k))/tdt_in

      ! Accretion of cloud liquid by falling graupel (from Lin et al 1983 - pgacw)
      ! This calculation uses the incoming fluxgraupel without subtracting sublimation
      ! (since subl occurs only outside cloud), so add sublflux back to fluxgraupel.
      ql(1:ifull)       = qlg(1:ifull,k)
      slopes_g(1:ifull) = ( max( fluxgraupel(:), 0. )/dz(:,k)/(pi*n0g*rho_g))**0.25
      where ( fluxgraupel(:)+sublflux(:)>0. .and. ql(1:ifull)>1.e-10 .and. ttg(1:ifull,k)<tfrz )
        cdt(1:ifull)         = tdt*pi*n0g*gam350*gcon/4.0*slopes_g(:)**3.5/sqrt(rhoa(:,k))
        dql(1:ifull)         = max( min( cgfra(:)*ql(:), ql(:)*cdt(:)/(1.+0.5*cdt(:)) ), 0. )
        qlg(1:ifull,k)       = qlg(1:ifull,k)   - dql(:)
        qaccr(1:ifull,k)     = qaccr(1:ifull,k) + dql(:)
        fluxgraupel(1:ifull) = fluxgraupel(:) + rhodz(:)*dql(:)
        dttg(1:ifull)        = hlfcp*dql(:)
        ttg(1:ifull,k)       = ttg(1:ifull,k) + dttg(:)
        qsatg(1:ifull,k)     = qsatg(1:ifull,k) + gam(:,k)*dttg(:)/hlscp
        cftmp(1:ifull)       = clfr(:,k)*dql(:)/ql(1:ifull)
        clfr(1:ifull,k)      = clfr(:,k) - cftmp(:)
        caccr_g(1:ifull)     = max( caccr_g(1:ifull), cftmp(1:ifull) )
      end where
    
      ! Accretion of rain by falling graupel (from Lin et al 1983 - pgacr)
      ! (Neglected in UM and ACCESS 1.3)
      qrn(1:ifull)      = rhor(:,k)/rhoa(:,k)
      slopes_g(1:ifull) = ( max( fluxgraupel(:), 0. )/dz(:,k)/(pi*n0g*rho_g))**0.25
      slopes_r(1:ifull) = (( max( rhor(:,k), 0. )*rhodz(:)/max( cfrain(:,k),1.e-15 )/tdt)**0.22)/714.
      where ( fluxgraupel(:)+sublflux(:)>0. .and. qrn(1:ifull)>1.e-10 .and. ttg(1:ifull,k)<tfrz )
        cdt(1:ifull)         = tdt*pi*pi*n0g*n0r*abs(vg2(:,k)-vl2(:,k))*qsn(:)*(rho_r/rhoa(:,k))     &
                               *(5.*slopes_r(:)**6*slopes_g(:)+2.*slopes_r(:)**5*slopes_g(:)**2      &
                                +0.5*slopes_r(:)**4*slopes_g(:)**3)          
        dql(1:ifull)         = max( min( cgfra(:)*qrn(:), qrn(:)*cdt(:)/(1.+0.5*cdt(:)) ), 0. )
        rhor(1:ifull,k)      = rhor(:,k) - rhoa(:,k)*dql(:)
        fluxgraupel(1:ifull) = fluxgraupel(:) + rhodz(:)*dql(:)
        dttg(1:ifull)        = hlfcp*dql(:)
        ttg(1:ifull,k)       = ttg(1:ifull,k) + dttg(:)
        qsatg(1:ifull,k)     = qsatg(1:ifull,k) + gam(:,k)*dttg(:)/hlscp 
        cftmp(1:ifull)       = cfrain(:,k)*dql(:)/qrn(:)
        cfrain(1:ifull,k)    = cfrain(:,k) - cftmp(:)
        caccr_g(1:ifull)     = max( caccr_g(:), cftmp(:) )        
      end where     
    
      ! Accretion of cloud ice by falling graupel (from Lin et al 1983 - pgaci)
      ! (Neglected in UM and ACCESS 1.3)
      qf(1:ifull)       = rhoi(:,k)/rhoa(:,k)
      slopes_g(1:ifull) = ( max( fluxgraupel(:), 0. )/dz(:,k)/(pi*n0g*rho_g))**0.25
      where ( fluxgraupel(:)+sublflux(:)>0. .and. qf(1:ifull)>1.e-10 .and. ttg(1:ifull,k)<tfrz )
        cdt(1:ifull)         = tdt*0.1*pi*n0g*gam350*gcon/4.*slopes_g(:)**3.5/sqrt(rhoa(:,k))
        dqf(1:ifull)         = max( min( cgfra(:)*qf(:), qf(:)*cdt(:)/(1.+0.5*cdt(:)) ), 0. )
        qaccf(1:ifull,k)     = qaccf(:,k) + dqf(:)      
        rhoi(1:ifull,k)      = rhoi(:,k) - rhoa(:,k)*dqf(:)
        fluxgraupel(1:ifull) = fluxgraupel(:) + rhodz(:)*dqf(:)
        cftmp(1:ifull)       = cifr(:,k)*dqf(:)/qf(:)
        cifr(1:ifull,k)      = cifr(:,k) - cftmp(:)
        caccf_g(1:ifull)     = max( caccf_g(:), cftmp(:) )
      end where

      ! Accretion of snow by falling graupel (from Lin et al 1983 - pgacs )
      qsn(1:ifull) = rhos(:,k)/rhoa(:,k)
      slopes_s(1:ifull) = (max( rhos(:,k), 0. )*rhoa(:,k)/(pi*rho_s*n0s(:)))**0.25
      slopes_g(1:ifull) = (max( fluxgraupel(:), 0.)/dz(:,k)/(pi*n0g*rho_g))**0.25
      where ( fluxgraupel(:)+sublflux(:)>0. .and. qsn(1:ifull)>1.e-10 .and. ttg(1:ifull,k)<tfrz )
        cdt(1:ifull)         = tdt*pi*pi*n0g*n0s(:)*abs(vg2(:,k)-vs2(:,k))*qsn(:)*(rho_s/rhoa(:,k))  &
                               *(5.*slopes_s(:)**6*slopes_g(:)+2.*slopes_s(:)**5*slopes_g(:)**2      &
                                +0.5*slopes_s(:)**4*slopes_g(:)**3)        
        dqf(1:ifull)         = max( min( cgfra(:)*qsn(:), qsn(:)*cdt(:)/(1.+0.5*cdt(:)) ), 0. )
        qaccf(1:ifull,k)     = qaccf(:,k) + dqf(:)
        rhos(1:ifull,k)      = rhos(:,k) - rhoa(:,k)*dqf(:)
        fluxgraupel(1:ifull) = fluxgraupel(:) + rhodz(:)*dqf(:)
        cftmp(1:ifull)       = cfsnow(:,k)*dqf(:)/qsn(:)
        cfsnow(1:ifull,k)    = cfsnow(:,k) - cftmp(:)
        caccf_g(1:ifull)     = max( caccf_g(1:ifull), cftmp(:) )
      end where

  
      ! Snow ------------------------------------------------------------------------------
      alphaf(1:ifull)   = hls*qsatg(1:ifull,k)/(rvap*ttg(1:ifull,k)**2)
      gam(1:ifull,k)    = hlscp*alphaf(:) !(L/cp)*dqsdt (HBG notation)
      sublflux(1:ifull) = 0.
      caccr_s(1:ifull)  = 0.
      caccf_s(1:ifull)  = 0.

      ! The following flag detects max/random overlap clouds
      ! that are separated by a clear layer
      where ( cfsnow(1:ifull,k)<1.e-10 .or. nmr==0 )
        ! combine max overlap from last cloud with net random overlap
        rdclfrsnow(1:ifull) = rdclfrsnow(:) + mxclfrsnow(:) - rdclfrsnow(:)*mxclfrsnow(:)
        mxclfrsnow(1:ifull) = 0.
      end where
      
      fluxsnow(:) = fluxsnow(:) + fluxprecipitation(:,k)*tdt/tdt_in
  
      ! Snow fall speed (from Lin et al 1983 - see GFDL AM3)
      where ( cfsnow(1:ifull,k)>=1.e-10 )
        vs2(1:ifull,k) = max( 0.1, 6.63*(max( rhos(:,k), 0. )/cfsnow(:,k)/942477796.)**0.0625 )
      elsewhere
        vs2(1:ifull,k) = vs2(1:ifull,k+1)
      end where

      ! Set up the parameters for the flux-divergence calculation
      alph(:)        = tdt*vs2(:,k)/dz(:,k)
      foutsnow(:,k)  = 1. - exp(-alph(:))          !analytical
      fthrusnow(:,k) = 1. - foutsnow(:,k)/alph(:)  !analytical

      ! Melt falling snow if > 0 deg C due to rain accretion
      ! (based on Lin et al 83, but using 0.65 and 0.44 coeffs following the UM approach)
!      qsn(1:ifull)      = fluxsnow(1:ifull)/rhodz(1:ifull)
!      slopes_s(1:ifull) = ( max( fluxsnow(:), 0. )/dz(:,k)/(pi*rho_s*n0s(:)))**0.25
!      where ( ttg(1:ifull,k)>tfrz .and. qsn(1:ifull)>1.e-10 )
!        cdt(1:ifull)        = tdt*2.*pi*n0s(:)/hlf*(tcond*(ttg(1:ifull,k)-tfrz)/rhoa(:,k)-vdifu*hl*(qsatg(:,k)-qtg(1:ifull,k))) &
!                                 *(0.65*slopes_s(:)**2+0.44*scm3*gam263*sqrt(clin/visk)*slopes_s(:)**2.63*sqrt(denfac(:)))
!        qif(1:ifull)        = max( min( qsn(:), qsn(:)*cdt(:) ), 0. ) !Mixing ratio of snow
!        dttg(1:ifull)       = -hlfcp*qif(:)
!        ttg(1:ifull,k)      = ttg(1:ifull,k) + dttg(:)
!        qsatg(1:ifull,k)    = qsatg(1:ifull,k) + gam(:,k)*dttg(:)/hlscp
!        fluxmelt(1:ifull)   = fluxmelt(:) + qif(:)*rhodz(:)
!        fluxsnow(1:ifull)   = fluxsnow(:) - qif(:)*rhodz(:)
!        rdclfrsnow(1:ifull) = rdclfrsnow(:)*(1.-qif(:)/qsn(:))
!        mxclfrsnow(1:ifull) = mxclfrsnow(:)*(1.-qif(:)/qsn(:))
!        cftmp(1:ifull)      = mxclfrsnow(:) + rdclfrsnow(:) - mxclfrsnow(:)*rdclfrsnow(:)
!        cfmelt(1:ifull)     = max( cfmelt(:), max( csfra(:)-cftmp(:), 0. ) )
!        csfra(1:ifull)      = cftmp(:)      
!      end where
    
      ! Melt falling snow if > 2 deg C
      qsn(1:ifull) = fluxsnow(1:ifull)/rhodz(1:ifull)
      where ( ttg(1:ifull,k)>tfrz+2. .and. qsn(1:ifull)>1.e-10 )
        cdt(1:ifull)        = min( 1., tdt/tau_s )*(ttg(1:ifull,k)-tfrz-2.)/hlfcp
        qif(1:ifull)        = min( qsn(:), cdt(:) ) !Mixing ratio of graupel
        dttg(1:ifull)       = -hlfcp*qif(:)
        ttg(1:ifull,k)      = ttg(1:ifull,k) + dttg(:)
        qsatg(1:ifull,k)    = qsatg(1:ifull,k) + gam(:,k)*dttg(:)/hlscp
        fluxmelt(1:ifull)   = fluxmelt(:) + qif(:)*rhodz(:)
        fluxsnow(1:ifull)   = fluxsnow(:) - qif(:)*rhodz(:)
        rdclfrsnow(1:ifull) = rdclfrsnow(1:ifull)*(1.-qif(1:ifull)/qsn(1:ifull))
        mxclfrsnow(1:ifull) = mxclfrsnow(1:ifull)*(1.-qif(1:ifull)/qsn(1:ifull))
        cftmp(1:ifull)      = mxclfrsnow(:) + rdclfrsnow(:) - mxclfrsnow(:)*rdclfrsnow(:)
        cfmelt(1:ifull)     = max( cfmelt(:), max( csfra(:)-cftmp(:), 0. ) )
        csfra(1:ifull)      = cftmp(:)      
      end where
    
      ! Compute the sublimation of snow falling from level k+1 into level k
      ! (Currently treated the same as LDR97 ice sublimation - see UM and ACCESS 1.3)
      fsclr_s(1:ifull)  = max( (1.-cifr(:,k)-clfr(:,k))*fluxsnow(:), 0. )
      slopes_s(1:ifull) = ( max( fluxsnow(:), 0. )/dz(:,k)/(pi*rho_s*n0s(:)))**0.25
      where ( fluxsnow(:)>0. .and. qtg(1:ifull,k)<qsatg(1:ifull,k) ) ! sublime snow
        cdt(1:ifull)      = 2.*pi*vdifu*tcond*rvap*n0s(:)*ttg(1:ifull,k)**2                                          &
                           *(0.65*slopes_s(:)**2+0.44*scm3*gam263*sqrt(clin/visk)*slopes_s(:)**2.63*sqrt(denfac(:))) &
                           /(tcond*rvap*ttg(1:ifull,k)**2+hls**2*vdifu*qsatg(:,k)*rhoa(1:ifull,k))
        dqs(1:ifull)      = tdt*cdt(:)*(qsatg(:,k)-qtg(1:ifull,k))
        dqs(1:ifull)      = min( dqs(:), (qsatg(1:ifull,k)-qtg(1:ifull,k))/(1.+gam(:,k)) ) !Don't supersat.
        sublflux(1:ifull) = min( dqs(:)*rhodz(:), fsclr_s(:) )
        dqs(1:ifull)      = sublflux(:)/rhodz(:)        
        fluxsnow(1:ifull) = fluxsnow(:) - sublflux(:)
        fsclr_s(1:ifull)  = fsclr_s(:) - sublflux(:)
        qsubl(1:ifull,k)  = qsubl(:,k) + dqs(:)
        qtg(1:ifull,k)    = qtg(1:ifull,k) + dqs(:)
        dttg(1:ifull)     = -hlscp*dqs(:)
        ttg(1:ifull,k)    = ttg(1:ifull,k) + dttg(:)
        qsatg(1:ifull,k)  = qsatg(1:ifull,k) + gam(:,k)*dttg(:)/hlscp
      end where

      ! Save flux for the wet deposition scheme.
      pfstayice(:,k) = pfstayice(:,k) + fluxsnow(:)*(1.-fthrusnow(:,k))/tdt_in

      ! Accretion of cloud liquid by falling snow (from Lin et al 1983 - psacw)
      ql(1:ifull)       = qlg(1:ifull,k)
      slopes_s(1:ifull) = ( max(fluxsnow(:), 0. )/dz(:,k)/(pi*rho_s*n0s(:)))**0.25
      where ( fluxsnow(:)+sublflux(:)>0. .and. qlg(1:ifull,k)>1.e-10 .and. ttg(1:ifull,k)<tfrz )
        cdt(1:ifull)      = tdt*denfac(:)*pi*clin*gam325*n0s(:)/4.*slopes_s(:)**3.25
        dql(1:ifull)      = max( min( csfra(:)*ql(:), ql(:)*cdt(:)/(1.+0.5*cdt(:)) ), 0. )
        qlg(1:ifull,k)    = qlg(1:ifull,k) - dql(:)
        qaccr(1:ifull,k)  = qaccr(:,k)     + dql(:)
        fluxsnow(1:ifull) = fluxsnow(:) + rhodz(:)*dql(:)
        dttg(1:ifull)     = hlfcp*dql(:)
        ttg(1:ifull,k)    = ttg(1:ifull,k) + dttg(:)
        qsatg(1:ifull,k)  = qsatg(1:ifull,k) + gam(:,k)*dttg(:)/hlscp
        cftmp(1:ifull)    = clfr(:,k)*dql(:)/ql(1:ifull)
        clfr(1:ifull,k)   = clfr(:,k) - cftmp(:)
        caccr_s(1:ifull)  = max( caccr_s(1:ifull), cftmp(1:ifull) )
      end where

      ! Accretion of rain by falling snow to form snow (from Lin et al 1983 - psacr)
      qrn(1:ifull)      = rhor(:,k)/rhoa(:,k)
      slopes_r(1:ifull) = (( max( rhor(:,k), 0. )*rhodz(:)/max( cfrain(:,k),1.e-15 )/tdt)**0.22)/714.
      slopes_s(1:ifull) = ( max( fluxsnow(:), 0. )/dz(:,k)/(pi*rho_s*n0s(:)))**0.25
      where ( fluxsnow(:)+sublflux(:)>0. .and. qrn(1:ifull)>1.e-10 .and. ttg(1:ifull,k)<tfrz )
        cdt(1:ifull)       = tdt*pi*pi*n0r*n0s(:)*abs(vs2(:,k)-vl2(:,k))*qrn(:)*(rho_r/rhoa(:,k))  &
                             *(5.*slopes_r(:)**6*slopes_s(:)+2.*slopes_r(:)**5*slopes_s(:)**2      &
                              +0.5*slopes_r(:)**4*slopes_s(:)**3)
        dql(1:ifull)       = max( min( clfra(:)*qrn(:), qrn(:)*cdt(:)/(1.+0.5*cdt(:)) ), 0. )
        rhor(1:ifull,k)    = rhor(:,k) - rhoa(:,k)*dql(:)
        fluxsnow(1:ifull)  = fluxsnow(:) + rhodz(:)*dql(:)
        dttg(1:ifull)      = hlfcp*dql(:)
        ttg(1:ifull,k)     = ttg(1:ifull,k) + dttg(:)
        qsatg(1:ifull,k)   = qsatg(1:ifull,k) + gam(:,k)*dttg(:)/hlscp  
        cftmp(1:ifull)     = cfrain(:,k)*dql(:)/qrn(:)
        cfrain(1:ifull,k)  = cfrain(:,k) - cftmp(:)
        caccr_s(1:ifull)   = max( caccr_s(:), cftmp(:) )
      end where

      ! Accretion of rain by falling snow to form graupel (neglected in Lin83 but included in UM)   
    
      ! Accretion of cloud ice by falling snow (from HDC 2004 - psaci)
      qf(1:ifull) = rhoi(:,k)/rhoa(:,k)
      slopes_s(1:ifull) = ( max( fluxsnow(:), 0. )/dz(:,k)/(pi*rho_s*n0s(:)))**0.25
      where ( fluxsnow(:)+sublflux(:)>0. .and. qf(1:ifull)>1.e-10 .and. ttg(1:ifull,k)<tfrz )
        esi(1:ifull)       = exp(0.05*max(ttg(1:ifull,k)-tfrz,-100.))       ! efficiency
        cdt(1:ifull)       = tdt*denfac(:)*27.737*n0s(:)*esi(:)*slopes_s(:)**3.41
        dqf(1:ifull)       = max( min( csfra(:)*qf(:), qf(:)*cdt(:)/(1.+0.5*cdt(:)) ), 0. )
        qaccf(1:ifull,k)   = qaccf(:,k) + dqf(:)
        rhoi(1:ifull,k)    = rhoi(:,k) - rhoa(:,k)*dqf(:)
        fluxsnow(1:ifull)  = fluxsnow(:) + rhodz(:)*dqf(:)
        cftmp(1:ifull)     = cifr(:,k)*dqf(:)/qf(:)
        cifr(1:ifull,k)    = cifr(:,k) - cftmp(:)
        caccf_s(1:ifull)   = max( caccf_s(:), cftmp(:) )
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
      curly(1:ifull) = 0.65*slopes_i(:)**2+0.493*slopes_i(:)*sqrt(slopes_i(:)*vi2(:,k+1)*rhoa(:,k)/um) !Factor in curly brackets
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
  
    ! Set up snow fall speed field
    select case(abs(ldr))
      case(1)
        where ( cifr(1:ifull,k)>=1.e-10 )
          vi2(1:ifull,k) = max( 0.1, 3.23*(max( rhoi(:,k), 0. )/cifr(:,k))**0.17 )  ! Ice fall speed from LDR 1997
          !vi2(1:ifull,k) = max( 0.1, 3.29*(max( rhoi(:,k), 0. )/cifr(:,k))**0.16 ) ! from Lin et al 1983
        end where
      case(2)
        where ( cifr(:,k)>=1.e-10 )
          vi2(:,k) = 0.9*3.23*(rhoi(:,k)/cifr(:,k))**0.17
        end where
      case(3)
        where ( cifr(1:ifull,k)>=1.e-10 )
          vi2(1:ifull,k) = max( 0.1, 2.05+0.35*log10(rhoi(:,k)/rhoa(:,k)/cifr(1:ifull,k)) )
        end where
      case(4)
        where ( cifr(1:ifull,k)>=1.e-10 )
          vi2(1:ifull,k) = 1.4*3.23*(rhoi(:,k)/cifr(1:ifull,k))**0.17
        end where
      case(11)
        ! following are alternative slightly-different versions of above
        ! used for I runs from 29/4/05 till 30/8/05
        ! for given qfg, large cifr implies small ice crystals, 
        ! with a small fall speed. 
        ! Note that for very small qfg, cifr is small.
        ! But rhoi is like qfg, so ratio should also be small and OK.
        vi2(1:ifull,k) = max( vi2(1:ifull,k+1), 3.23*(rhoi(:,k)/max(cifr(1:ifull,k),1.e-30))**0.17 )
      case(22)
        vi2(1:ifull,k) = max( vi2(1:ifull,k+1), 0.9*3.23*(rhoi(:,k)/max(cifr(1:ifull,k),1.e-30))**0.17 )
      case(33)
        ! following max gives vi2=.1 for qfg=cifr=0
        vi2(1:ifull,k) = max( vi2(1:ifull,k+1), 2.05+0.35*log10(max(rhoi(:,k)/rhoa(:,k),2.68e-36)/max(cifr(1:ifull,k),1.e-30)) )
    end select

    ! Set up the parameters for the flux-divergence calculation
    alph(:)       = tdt*vi2(:,k)/dz(:,k)
    foutice(:,k)  = 1. - exp(-alph(:))         !analytical
    fthruice(:,k) = 1. - foutice(:,k)/alph(:)  !analytical
  
    ! Melt falling ice if > 0 deg C
    where ( ttg(1:ifull,k)>tfrz .and. fluxice(:)>0. )
      qif(1:ifull)        = fluxice(:)/rhodz(:)      !Mixing ratio of ice
      dttg(1:ifull)       = -hlfcp*qif(:)
      ttg(1:ifull,k)      = ttg(1:ifull,k) + dttg(:)
      qsatg(1:ifull,k)    = qsatg(1:ifull,k) + gam(:,k)*dttg(:)/hlscp
      fluxmelt(1:ifull)   = fluxmelt(:) + fluxice(:)
      cfmelt(1:ifull)     = max( cfmelt(:), cifra(:) )
      fluxice(1:ifull)    = 0.
      cifra(1:ifull)      = 0.
      rdclfrice(1:ifull)  = 0.
      mxclfrice(1:ifull)  = 0.
    end where

    ! Compute the sublimation of ice falling from level k+1 into level k
    fsclr_i(:) = (1.-cifr(:,k)-clfr(:,k))*fluxice(:)
    where ( fluxice(:)>0. .and. qtg(1:ifull,k)<qsatg(1:ifull,k) ) ! sublime ice
      Csb(1:ifull)      = Csbsav(:)*fluxice(:)/tdt
      bf(1:ifull)       = 1. + 0.5*Csb(:)*tdt*(1.+gam(:,k))
      dqs(1:ifull)      = max( 0., tdt*(Csb(:)/bf(:))*(qsatg(1:ifull,k)-qtg(1:ifull,k)) )
      dqs(1:ifull)      = min( dqs(:), (qsatg(1:ifull,k)-qtg(1:ifull,k))/(1.+gam(:,k)) ) !Don't supersat.
      sublflux(1:ifull) = min( dqs(:)*rhodz(:), fsclr_i(:) )
      dqs(1:ifull)      = sublflux(:)/rhodz(:)      
      fluxice(1:ifull)  = fluxice(:) - sublflux(:)
      fsclr_i(1:ifull)  = fsclr_i(:) - sublflux(:)
      qsubl(1:ifull,k)  = qsubl(:,k) + dqs(:)
      qtg(1:ifull,k)    = qtg(1:ifull,k) + dqs(:)
      dttg(1:ifull)     = -hlscp*dqs(:)
      ttg(1:ifull,k)    = ttg(1:ifull,k) + dttg(:)
      qsatg(1:ifull,k)  = qsatg(1:ifull,k) + gam(:,k)*dttg(:)/hlscp
    end where

    ! Save flux for the wet deposition scheme.
    pfstayice(:,k) = pfstayice(:,k) + fluxice(:)*(1.-fthruice(:,k))/tdt_in
  
    ! Accretion of cloud liquid by falling ice (neglected in Lin et al 1983, but
    ! included in UM and ACCESS 1.3 as piacw)
    ! This calculation uses the incoming fluxice without subtracting sublimation
    ! (since subl occurs only outside cloud), so add sublflux back to fluxice.
    ql(1:ifull)      = qlg(1:ifull,k)
    where ( fluxice(:)+sublflux(:)>0. .and. qlg(1:ifull,k)>1.e-10 )
      cdt(1:ifull)     = Eac*slopes_i(:)*(fluxice(:)+sublflux(:))/(2.*rhosno)
      dql(1:ifull)     = max( min( cifra(:)*ql(:), ql(:)*cdt(:)/(1.+0.5*cdt(:)) ), 0. )
      qlg(1:ifull,k)   = qlg(1:ifull,k) - dql(:)
      qaccr(1:ifull,k) = qaccr(:,k) + dql(:)
      fluxice(1:ifull) = fluxice(:) + rhodz(:)*dql(:)
      dttg(1:ifull)    = hlfcp*dql(:)
      ttg(1:ifull,k)   = ttg(1:ifull,k) + dttg(:)
      qsatg(1:ifull,k) = qsatg(1:ifull,k) + gam(:,k)*dttg(:)/hlscp
      cftmp(1:ifull)   = clfr(:,k)*dql(:)/ql(:)
      clfr(1:ifull,k)  = clfr(:,k) - cftmp(:)
      caccr_i(1:ifull) = max( caccr_i(:), cftmp(:) )
    end where
  
    ! Accretion of rain by falling ice to produce ice (from Lin et al 1983 - piacr)
    ! (see UM and ACCESS 1.3 piacr-c for an alternate formulation)
    qrn(1:ifull) = rhor(:,k)/rhoa(:,k)
    qf(1:ifull)  = (fluxice(1:ifull)+sublflux(1:ifull))/rhodz(1:ifull)
    where ( fluxice(:)+sublflux(:)>0. .and. qrn(1:ifull)>1.e-10 .and. ttg(1:ifull,k)<tfrz .and. ncloud>=3 )
      cdt(1:ifull)       = tdt*denfac(:)*c_piacr*qf(:)/sqrt(rhoa(:,k))
      dql(1:ifull)       = max( min( cifra(:)*qrn(:), qrn(:)*cdt(:)/(1.+0.5*cdt(:)) ), 0. )
      rhor(1:ifull,k)    = rhor(1:ifull,k) - rhoa(:,k)*dql(:)
      fluxice(1:ifull)   = fluxice(:)  + rhodz(:)*dql(:)
      dttg(1:ifull)      = hlfcp*dql(:)
      ttg(1:ifull,k)     = ttg(1:ifull,k) + dttg(:)
      qsatg(1:ifull,k)   = qsatg(1:ifull,k) + gam(:,k)*dttg(:)/hlscp
      cftmp(1:ifull)     = cfrain(:,k)*dql(:)/qrn(:)
      cfrain(1:ifull,k)  = cfrain(:,k) - cftmp(:)
      caccr_i(1:ifull)   = max( caccr_i(:), cftmp(:) )
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
    fluxrain(:) = fluxrain(:) + fluxmelt(:) + fluxauto(:,k)*tdt/tdt_in
    
    ! Calculate rain fall speed (MJT)
    if ( ncloud>1 ) then
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
    qpf(:)     = fluxrain(:)/rhodz(:) !Mix ratio of rain which falls into layer
    clrevap(:) = (1.-clfr(:,k)-cifr(:,k))*qpf(:)
    Tk(:)      = ttg(1:ifull,k)
    qsatg(:,k) = qsati(pk(:),Tk(:))
    where ( Tk(:)<tfrz .and. Tk(:)>=tice )
      qsl(:) = qsatg(:,k) + epsil*esdiffx(Tk(:))/pk(:)
    elsewhere
      qsl(:) = qsatg(:,k)
    end where
    where ( fluxrain(:)>0. .and. clfra(:)>0. )
      es(:)      = qsl(:)*pk(:)/epsil 
      Apr(:)     = (hl/(rKa*Tk(:)))*(hl/(rvap*Tk(:))-1.)
      Bpr(:)     = rvap*Tk(:)/((Dva/pk(:))*es(:))
      Fr(:)      = fluxrain(:)/tdt/clfra(:)
      Cev(:)     = clfra(:)*3.8e2*sqrt(Fr(:)/rhoa(:,k))/(qsl(:)*(Apr(:)+Bpr(:)))
      dqsdt(:)   = hl*qsl(:)/(rvap*Tk(:)**2)
      bl(:)      = 1. + 0.5*Cev(:)*tdt*(1.+hlcp*dqsdt(:))
      evap(:)    = tdt*(Cev(:)/bl(:))*(qsl(:)-qtg(:,k))
      satevap(:) = (qsl(:)-qtg(:,k))/(1.+hlcp*dqsdt(:)) !Evap to saturate
      ! Vl2=11.3*Fr**(1./9.)/sqrt(rhoa(mg,k))    !Actual fall speed
      ! Vl2=5./sqrt(rhoa(mg,k))                  !Nominal fall speed
      evap(:) = max( 0., min( evap(:), satevap(:), clrevap(:) ) )
    end where
    if ( nevapls==-1 ) then
      evap(1:ifull) = 0.
    else if ( nevapls==-2 ) then
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

    ! store liquid flux for aerosols
    pfstayliq(:,k) = pfstayliq(:,k) + fluxrain(:)*(1.-fthruliq(:,k))/tdt_in
  
    ! Now do the collection of liquid cloud by rain term (cf. pracc in Lin83).
    where ( fluxrain(:)>0. )
      Fr(1:ifull)       = fluxrain(:)/tdt/max( clfra(:), 1.e-15 )
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
    fcol(1:ifull)     = min( 1., mxclfrliq(:)/(1.e-20+clfr(:,k)) )     !max overlap
    fcol(1:ifull)     = fcol(:) + rdclfrliq(:) - fcol(:)*rdclfrliq(:)  !rnd overlap
    cdt(1:ifull)      = tdt*Ecol*0.24*fcol(:)*pow75(Fr(:))
    prscav(1:ifull,k) = tdt*0.24*fcol(:)*pow75(Fr(:))                  !Strat only
    coll(1:ifull)     = max( min( qlg(1:ifull,k), qlg(1:ifull,k)*cdt(:)/(1.+0.5*cdt(:)) ), 0. )
    qcoll(1:ifull,k)  = qcoll(:,k) + coll(:)
    qlg(1:ifull,k)    = qlg(1:ifull,k) - coll(:)
    fluxrain(1:ifull) = fluxrain(:) + coll(:)*rhodz(:)
  
    ! subtract evaporated rain
    fluxrain(:) = max( fluxrain(:) - rhodz(:)*evap(:), 0. ) !To avoid roundoff -ve's
    
    ! Freezing rain to produce graupel (pgfr)
    ! (Neglected in UM and ACCESS 1.3)
    qrn(1:ifull)      = fluxrain(1:ifull)/rhodz(1:ifull)
    slopes_r(1:ifull) = (( max( fluxrain(:), 0. )/max( clfra(:),1.e-15 )/tdt)**0.22)/714. ! from LDR97
    where ( qrn(1:ifull)>1.e-10 .and. ttg(1:ifull,k)<tfrz .and. ncloud>=3 )
      ! MJT notes - limit temperature to -100 C to avoid overflow with single precision
      cdt(1:ifull)         = tdt*20.e2*pi*pi*n0r*(rho_r/rhoa(:,k))*slopes_r(:)**7 &
                             *(exp(-0.66*max( ttg(1:ifull,k)-tfrz, -100. ))-1.)
      dql(1:ifull)         = max( min( qrn(1:ifull), qrn(1:ifull)*cdt(:)/(1.+0.5*cdt(:)) ), 0. )
      fluxrain(1:ifull)    = fluxrain(:)    - rhodz(:)*dql(:)
      fluxgraupel(1:ifull) = fluxgraupel(:) + rhodz(:)*dql(:)
      dttg(1:ifull)        = hlfcp*dql(:)
      ttg(1:ifull,k)       = ttg(1:ifull,k) + dttg(:)
      qsatg(1:ifull,k)     = qsatg(1:ifull,k) + gam(:,k)*dttg(:)/hlscp
      rdclfrliq(1:ifull)   = rdclfrliq(:)*(1.-dql(:)/qrn(:))
      mxclfrliq(1:ifull)   = mxclfrliq(:)*(1.-dql(:)/qrn(:))
      cftmp(1:ifull)       = mxclfrliq(:) + rdclfrliq(:) - mxclfrliq(:)*rdclfrliq(:)
      cftmp(1:ifull)       = clfra(:) - cftmp(:)
      clfra(1:ifull)       = clfra(:) - cftmp(:)
      cfgraupel(1:ifull,k) = cfgraupel(:,k) + cftmp(:) - cfgraupel(:,k)*cftmp(:)
    end where
    
    ! Accretion of cloud snow by rain (from Lin et al 1983 - pracs)
    qsn(1:ifull)      = rhos(:,k)/rhoa(:,k)
    slopes_s(1:ifull) = ( max( rhos(:,k), 0. )*rhoa(:,k)/(pi*rho_s*n0s(:)))**0.25    
    slopes_r(1:ifull) = (( max( fluxrain(:), 0. )/max( clfra(:),1.e-15 )/tdt)**0.22)/714. ! from LDR97
    where ( fluxrain(1:ifull)>0. .and. qsn(1:ifull)>1.e-10 .and. ttg(1:ifull,k)>tfrz+1. .and. ncloud>=3 )
      cdt(1:ifull)         = tdt*pi*pi*n0r*n0s(:)*abs(vl2(:,k)-vs2(:,k))*qsn(:)*(rho_s/rhoa(:,k))   &
                             *(5.*slopes_s(:)**6*slopes_r(:)+2.*slopes_s(:)**5*slopes_r(:)**2       &
                              +0.5*slopes_s(:)**4*slopes_r(:)**3)
      dqf(1:ifull)         = max( min( clfra(:)*qsn(:), qsn(:)*cdt(:)/(1.+0.5*cdt(:)) ), 0. )
      rhos(1:ifull,k)      = rhos(:,k) - rhoa(:,k)*dqf(:)
      fluxrain(1:ifull)    = fluxrain(:) + rhodz(:)*dqf(:)
      dttg(1:ifull)        = hlfcp*dqf(:)
      ttg(1:ifull,k)       = ttg(1:ifull,k) - dttg(:)
      qsatg(1:ifull,k)     = qsatg(1:ifull,k) - gam(:,k)*dttg(:)/hlscp      
      cftmp(1:ifull)       = cfsnow(:,k)*dqf(:)/qsn(:)
      cfsnow(1:ifull,k)    = cfsnow(:,k) - cftmp(:)
      cfrain(1:ifull,k)    = cfrain(:,k) + cftmp(:) - cfrain(:,k)*cftmp(:)
    end where

    
    ! Liquid ------------------------------------------------------------------------------
    ! (Currently cloud droplet settling is negected, although included in UM and ACCESS 1.3)


    ! Misc ------------------------------------------------------------------------------
      
    ! Accretion of cloud ice by rain to produce snow or grauple (from Lin et al 1983 - praci)
    ! (Neglected in UM and ACCESS 1.3)
    qf(1:ifull)       = rhoi(:,k)/rhoa(:,k)
    rhodum_r(1:ifull) = fluxrain(1:ifull)/dz(:,k)
    slopes_r(1:ifull) = (( max( fluxrain(:), 0. )/max( clfra(:),1.e-15 )/tdt)**0.22)/714. ! from LDR97
    where ( fluxrain(1:ifull)>0. .and. qf(1:ifull)>1.e-10 .and. ttg(1:ifull,k)<tfrz .and. ncloud>=3 )
      cdt(1:ifull)           = tdt*pi*n0r*alin*gam380/4.*slopes_r(:)**3.8*denfac(:)
      xwgt(1:ifull)          = (rhodum_r(1:ifull)-0.995*qr0_crt)/(0.01*qr0_crt) ! MJT suggestion to switch from snow to graupel
      xwgt(1:ifull)          = max( min( xwgt(1:ifull), 1. ), 0. )
      dqf(1:ifull)           = max( min( clfra(1:ifull)*qf(1:ifull), qf(1:ifull)*cdt(1:ifull)/(1.+0.5*cdt(1:ifull)) ), 0. )
      rhoi(1:ifull,k)        = rhoi(:,k) - rhoa(:,k)*dqf(:)
      fluxgraupel(1:ifull)   = fluxgraupel(:) + rhodz(:)*dqf(:)*xwgt(:)
      fluxsnow(1:ifull)      = fluxsnow(:)    + rhodz(:)*dqf(:)*(1.-xwgt(:))
      cftmp(1:ifull)         = cifr(:,k)*dqf(:)/qf(:)
      cifr(1:ifull,k)        = cifr(:,k) - cftmp(:)
      cfgraupel(1:ifull,k)   = max( cfgraupel(1:ifull,k), cftmp(:)*xwgt(:) )
      cfsnow(1:ifull,k)      = max( cfsnow(1:ifull,k), cftmp(:)*(1.-xwgt(:)) )
    end where 
  
    
    ! Update fluxes and area fractions for graupel, snow, ice and rain
  
    fluxm(:,k) = fluxm(:,k) + fluxmelt(:)
  
    if ( ncloud>=3 ) then

      
      ! Grauple
      ! calculate maximum and random overlap for falling graupel
      cfgraupel(1:ifull,k) = min( 1., cfgraupel(1:ifull,k)+caccr_g(1:ifull)-cfgraupel(1:ifull,k)*caccr_g(1:ifull) )  ! rnd overlap
      cfgraupel(1:ifull,k) = max( cfgraupel(1:ifull,k), caccf_g(1:ifull) )                                           ! max overlap
      where ( fluxgraupel(:)<=0. )
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
      cfgraupelfall(1:ifull,k) = cfgraupel(1:ifull,k) - cffluxout(:) + cffluxin(:)*(1.-fthrugraupel(:,k))
      rhog(1:ifull,k)          = rhog(:,k) - rhogout(:) + rhogin(:)*(1.-fthrugraupel(:,k))
      fluxgraupel(1:ifull)     = max( rhogout(:)*dz(:,k) + fluxgraupel(:)*fthrugraupel(:,k), 0. )
      ! Now fluxgraupel is flux leaving layer k
      fluxg(1:ifull,k)         = fluxg(:,k) + fluxgraupel(:)      
   
    
      ! Snow
      ! calculate maximum and random overlap for falling snow
      cfsnow(1:ifull,k) = min( 1., cfsnow(1:ifull,k)+caccr_s(1:ifull)-cfsnow(1:ifull,k)*caccr_s(1:ifull) ) ! rnd overlap
      cfsnow(1:ifull,k) = max( cfsnow(1:ifull,k), caccf_s(1:ifull) )                                       ! max overlap
      where ( fluxsnow(:)<=0. )
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
      cfsnowfall(1:ifull,k) = cfsnow(1:ifull,k) - cffluxout(:) + cffluxin(:)*(1.-fthrusnow(:,k))
      rhos(1:ifull,k)       = rhos(:,k) - rhosout(:) + rhosin(:)*(1.-fthrusnow(:,k))
      fluxsnow(1:ifull)     = max( rhosout(:)*dz(:,k) + fluxsnow(:)*fthrusnow(:,k), 0. )
      ! Now fluxsnow is flux leaving layer k
      fluxs(1:ifull,k)      = fluxs(:,k) + fluxsnow(:)

    
    end if ! ncloud>=3

  
    ! Ice
    ! calculate maximum and random overlap for falling ice
    cifr(1:ifull,k) = min( 1., cifr(1:ifull,k)+caccr_i(:)-cifr(1:ifull,k)*caccr_i(:) )  ! rnd overlap
    where ( fluxice(:)<=0. )
      rdclfrice(1:ifull) = 0.
      mxclfrice(1:ifull) = 0.
    end where
    mxclfrice(1:ifull) = max( mxclfrice(:), cifr(1:ifull,k) )                             !max overlap
    cifra(1:ifull)     = max( 0.01, mxclfrice(:)+rdclfrice(:)-mxclfrice(:)*rdclfrice(:) ) !rnd overlap the mx and rd ice fractions
    ! Save sedimentation rate for aerosol scheme
    pqfsed(:,k) = pqfsed(:,k) + foutice(:,k)*tdt/tdt_in
    ! Compute fluxes into the box
    if ( ncloud<=2 ) then
      where ( rica(:)>0. )
        rhoiin(:)   = max( fluxice(:)/dz(:,k), 0. )
        cffluxin(:) = min( rhoiin(:)/rica(:), 1. )
      elsewhere
        rhoiin(:)   = 0.
        cffluxin(:) = 0.
      end where
      where ( cifr(:,k)>=1.e-10 ) 
        rica(:) = max( rhoi(:,k)/cifr(:,k), 0. ) ! in cloud rhoi
      end where
    else
      rhoiin(:)   = max( fluxice(:)/dz(:,k), 0. )
      cffluxin(:) = cifra(:) - cifr(:,k)
    end if
    ! Compute the fluxes of ice leaving the box
    where ( cifr(:,k)>=1.e-10 )
      cffluxout(:) = cifr(:,k)*foutice(:,k)
      rhoiout(:)   = rhoi(:,k)*foutice(:,k)
    elsewhere
      cffluxout(:) = 0.
      rhoiout(:)   = 0.
    end where
    ! Update the rhoi and cifr fields
    cifr(1:ifull,k)  = min( 1.-clfr(1:ifull,k), cifr(1:ifull,k)-cffluxout(:)+cffluxin(:)*(1.-fthruice(:,k)) )
    rhoi(1:ifull,k)  = rhoi(:,k) - rhoiout(:) + rhoiin(:)*(1.-fthruice(:,k))
    fluxice(1:ifull) = max( rhoiout(:)*dz(:,k) + fluxice(:)*fthruice(:,k), 0. )
    ! Now fluxice is flux leaving layer k
    fluxi(1:ifull,k) = fluxi(:,k) + fluxice(:)

  
    ! Rain
    ! Calculate the raining cloud cover down to this level, for stratiform (clfra).
    cfrain(1:ifull,k) = min( 1., cfrain(1:ifull,k)+cfmelt(:)-cfrain(1:ifull,k)*cfmelt(:) ) ! rnd overlap
    where ( fluxrain(:)<=0. )
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
    cfrainfall(1:ifull,k) = cfrain(1:ifull,k) - cffluxout(:) + cffluxin(:)*(1.-fthruliq(:,k))
    rhor(1:ifull,k)       = rhor(:,k) - rhorout(:) + rhorin(:)*(1.-fthruliq(:,k))
    fluxrain(1:ifull)     = max( rhorout(:)*dz(:,k) + fluxrain(:)*fthruliq(:,k), 0. )
    ! Now fluxrain is flux leaving layer k
    fluxr(1:ifull,k)      = fluxr(:,k) + fluxrain(:)

  end do ! k
  
end do   ! n


! Re-create qrg, qfg, qsng and qgrg fields
qrg(1:ifull,kl)     = 0.
qrg(1:ifull,1:kl-1) = rhor(1:ifull,1:kl-1)/rhoa(1:ifull,1:kl-1)
qfg(1:ifull,1:kl)   = rhoi(1:ifull,1:kl)/rhoa(1:ifull,1:kl)
qsng(1:ifull,1:kl)  = rhos(1:ifull,1:kl)/rhoa(1:ifull,1:kl)
qgrg(1:ifull,1:kl)  = rhog(1:ifull,1:kl)/rhoa(1:ifull,1:kl)


! store precip, snow and graupel
precs(1:ifull) = precs(1:ifull) + fluxr(1:ifull,1) + fluxi(1:ifull,1) + fluxs(1:ifull,1) + fluxg(1:ifull,1)
preci(1:ifull) = preci(1:ifull) + fluxi(1:ifull,1) + fluxs(1:ifull,1)
precg(1:ifull) = precg(1:ifull) + fluxg(1:ifull,1)


! Remove small amounts of cloud and precip
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

if ( ncloud>=3 ) then
  cfrac(1:ifull,1:kl) = clfr(1:ifull,1:kl) + cifr(1:ifull,1:kl)
end if

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
if ( diag .and. mydiag ) then  ! JLM
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
  write(6,*) 'fluxm',fluxm(idjd,:)
  write(6,*) 'cifra,fluxsnow',cifra(idjd),fluxsnow(idjd)
end if  ! (diag.and.mydiag)

return
end subroutine newsnowrain

end module leoncld_mod
