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
    
! ncloud = 0    Standard LDR cloud microphysics
! ncloud = 1    Use newer LDR autoconvection from Mk3.6
! ncloud = 2    Same as ncloud=1, but with prognostic rain
! ncloud = 3    RESERVED for grauple, snow and ice
! ncloud = 4    Use prognostic cloud fraction based on Tiedtke from GFDL-CM3, but autoconversion from ncloud=0
! ncloud = 5    Same as ncloud=4, but convective sources are included in prognostic cloud fraction

! Currently we are developing the grauple, snow and ice components, based on Lin et al 1983 and GFDL-AM3.
   
!                            Water vapour (qg)
!
!   Cloud water (qlg,cfrac)                      Cloud ice (qfg,cfrac)
!
!   Rain (qrg,rfrac)                             Snow (qsg,sfrac)         Grauple (qgrg,gfrac)

! qg, qlg, qfg, qrg, qsg and qgrg are mixing ratios (g/g) and cfrac, rfrac, sfrac, gfrac are area cover
! fractions
    
module leoncld_mod
    
private
public leoncld

contains
    
! This subroutine is the interface for the LDR cloud microphysics
subroutine leoncld
      
use aerointerface                 ! Aerosol interface
use arrays_m                      ! Atmosphere dyamics prognostic arrays
use cc_mpi, only : mydiag, myid   ! CC MPI routines
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
integer iq,k,ncl
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
real, dimension(ifull) :: precg                           !Amount of stratiform grauple in timestep (mm)
real, dimension(ifull) :: wcon                            !Convective cloud water content (in-cloud, prescribed)

real, dimension(ifull,kl) :: qevap, qsubl, qauto, qcoll, qaccr
real, dimension(ifull,kl) :: fluxr, fluxi, fluxs, fluxg, fluxmelt, pqfsed
real, dimension(ifull,kl) :: pfstayice, pfstayliq, slopes, prscav
real, dimension(ifull) :: prf_temp, fl, qtot


! Non-hydrostatic terms
tnhs(:,1) = phi_nh(:,1)/bet(1)
do k = 2,kl
  ! representing non-hydrostatic term as a correction to air temperature
  tnhs(:,k) = (phi_nh(:,k)-phi_nh(:,k-1)-betm(k)*tnhs(:,k-1))/bet(k)
end do

! meterological fields
do k = 1,kl
  prf(:,k)    = 0.01*ps(1:ifull)*sig(k)    !ps is SI units
  prf_temp(:) = 100.*prf(:,k)
  dprf(:,k)   = -0.01*ps(1:ifull)*dsig(k)  !dsig is -ve
  qtot(:)     = qg(1:ifull,k)+qlg(1:ifull,k)+qrg(1:ifull,k)+qfg(1:ifull,k)+qsg(1:ifull,k)+qgrg(1:ifull,k)
  tv(:,k)     = t(1:ifull,k)*(1.+1.61*qg(1:ifull,k)-qtot(:))                                   ! virtual temperature
  rhoa(:,k)   = prf_temp(:)/(rdry*tv(1:ifull,k))                                               ! air density
  qsatg(:,k)  = qsat(prf_temp(:),t(1:ifull,k))                                                 ! saturated mixing ratio
  dz(:,k)     = 100.*dprf(1:ifull,k)/(rhoa(1:ifull,k)*grav)*(1.+tnhs(1:ifull,k)/tv(1:ifull,k)) ! level thickness in metres
enddo
 
! Calculate droplet concentration from aerosols (for non-convective faction of grid-box)
call aerodrop(1,ifull,cdso4,rhoa,outconv=.true.)

! default values
kbase(:)    = 0  ! default
ktop(:)     = 0  ! default
precs(:)    = 0. ! rain
preci(:)    = 0. ! snow
precg(:)    = 0. ! hail

!     Set up convective cloud column
call convectivecloudfrac(clcon)
where ( ktsav<kl-1 )
  ktop   = ktsav
  kbase  = kbsav+1
  wcon   = wlc
elsewhere
  wcon   = 0.
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
  write(6,"('qs  ',9f8.3/4x,9f8.3)") qsg(idjd,:)
  write(6,"('qg  ',9f8.3/4x,9f8.3)") qgrg(idjd,:)
endif

! Calculate convective cloud fraction and adjust moisture variables before calling newcloud
if ( ncloud<=4 ) then

  ! diagnose cloud fraction (ncloud<=3) or prognostic strat. cloud but diagnostic conv. cloud (ncloud==4)
  if ( nmr>0 ) then
    ! Max/Rnd cloud overlap
    do k = 1,kl
      where ( clcon(:,k)>0. )
        !ccw=wcon(:)/rhoa(:,k)  !In-cloud l.w. mixing ratio
        qccon(:,k)       = clcon(:,k)*wcon(:)/rhoa(:,k)
        qcl(:,k)         = max(qsatg(:,k),qg(1:ifull,k))  ! jlm
        qenv(1:ifull,k)  = max(1.e-8,qg(1:ifull,k)-clcon(:,k)*qcl(:,k))/(1.-clcon(:,k))
        qcl(:,k)         = (qg(1:ifull,k)-(1.-clcon(:,k))*qenv(1:ifull,k))/clcon(:,k)
        qlg(1:ifull,k)   = qlg(1:ifull,k)/(1.-clcon(:,k))
        qfg(1:ifull,k)   = qfg(1:ifull,k)/(1.-clcon(:,k))
        qrg(1:ifull,k)   = qrg(1:ifull,k)/(1.-clcon(:,k))
        qsg(1:ifull,k)   = qsg(1:ifull,k)/(1.-clcon(:,k))
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
        qccon(1:ifull,k) = clcon(1:ifull,k)*wcon(1:ifull)/rhoa(1:ifull,k)
        qcl(1:ifull,k)   = max(qsatg(1:ifull,k),qg(1:ifull,k))  ! jlm
        qenv(1:ifull,k)  = max(1.e-8,qg(1:ifull,k)-clcon(1:ifull,k)*qcl(1:ifull,k))/(1.-clcon(1:ifull,k))
        qcl(1:ifull,k)   = (qg(1:ifull,k)-(1.-clcon(1:ifull,k))*qenv(1:ifull,k))/clcon(1:ifull,k)
        qlg(1:ifull,k)   = qlg(1:ifull,k)/(1.-clcon(1:ifull,k))
        qfg(1:ifull,k)   = qfg(1:ifull,k)/(1.-clcon(1:ifull,k))
        qrg(1:ifull,k)   = qrg(1:ifull,k)/(1.-clcon(1:ifull,k))
        qsg(1:ifull,k)   = qsg(1:ifull,k)/(1.-clcon(1:ifull,k))
        qgrg(1:ifull,k)  = qgrg(1:ifull,k)/(1.-clcon(1:ifull,k))
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
  clcon(:,:)      = 0.
  qccon(:,:)      = 0.
  qcl(1:ifull,:)  = qg(1:ifull,:)
  qenv(1:ifull,:) = qg(1:ifull,:)
end if
      
tenv(:,:) = t(1:ifull,:) !Assume T is the same in and out of convective cloud

if ( nmaxpr==1 .and. mydiag ) then
  write(6,*) 'before newcloud',ktau
  write(6,"('t   ',9f8.2/4x,9f8.2)") t(idjd,:)
  write(6,"('qv  ',9f8.3/4x,9f8.3)") qg(idjd,:)
  write(6,"('qf  ',9f8.3/4x,9f8.3)") qfg(idjd,:)
  write(6,"('ql  ',9f8.3/4x,9f8.3)") qlg(idjd,:)
  write(6,"('qr  ',9f8.3/4x,9f8.3)") qrg(idjd,:)
  write(6,"('qs  ',9f8.3/4x,9f8.3)") qsg(idjd,:)
  write(6,"('qg  ',9f8.3/4x,9f8.3)") qgrg(idjd,:)
  write(6,"('qnv ',9f8.3/4x,9f8.3)") qenv(idjd,:)
  write(6,"('qsat',9f8.3/4x,9f8.3)") qsatg(idjd,:)
  write(6,"('qcl ',9f8.3/4x,9f8.3)") qcl(idjd,:)
  write(6,"('clc ',9f8.3/4x,9f8.3)") clcon(idjd,:)
  write(6,*) 'kbase,ktop ',kbase(idjd),ktop(idjd)
endif

!     Calculate cloud fraction and cloud water mixing ratios
call newcloud(dt,land,prf,kbase,ktop,rhoa,cdso4,tenv,qenv,qlg,qfg,cfrac,ccov,cfa,qca)

if ( nmaxpr==1 .and. mydiag ) then
  write(6,*) 'after newcloud',ktau
  write (6,"('tnv ',9f8.2/4x,9f8.2)") tenv(idjd,:)
  write (6,"('qv0 ',9f8.3/4x,9f8.3)") qg(idjd,:)
  write (6,"('qf  ',9f8.3/4x,9f8.3)") qfg(idjd,:)
  write (6,"('ql  ',9f8.3/4x,9f8.3)") qlg(idjd,:)
  write (6,"('qr  ',9f8.3/4x,9f8.3)") qrg(idjd,:)
  write (6,"('qs  ',9f8.3/4x,9f8.3)") qsg(idjd,:)
  write (6,"('qg  ',9f8.3/4x,9f8.3)") qgrg(idjd,:)
  write (6,"('qnv ',9f8.3/4x,9f8.3)") qenv(idjd,:) ! really new qg
endif

!     Weight output variables according to non-convective fraction of grid-box            
do k = 1,kl
  t(1:ifull,k)  = clcon(:,k)*t(1:ifull,k)+(1.-clcon(:,k))*tenv(:,k)
  qg(1:ifull,k) = clcon(:,k)*qcl(:,k)+(1.-clcon(:,k))*qenv(:,k)
  where ( k>=kbase .and. k<=ktop )
    cfrac(:,k)      = cfrac(:,k)*(1.-clcon(:,k))
    rfrac(1:ifull,k) = rfrac(1:ifull,k)*(1.-clcon(:,k))
    sfrac(1:ifull,k) = sfrac(1:ifull,k)*(1.-clcon(:,k))
    gfrac(1:ifull,k) = gfrac(1:ifull,k)*(1.-clcon(:,k))
    ccov(:,k)       = ccov(:,k)*(1.-clcon(:,k))              
    qlg(1:ifull,k)  = qlg(1:ifull,k)*(1.-clcon(:,k))
    qfg(1:ifull,k)  = qfg(1:ifull,k)*(1.-clcon(:,k))
    qrg(1:ifull,k)  = qrg(1:ifull,k)*(1.-clcon(:,k))
    qsg(1:ifull,k)  = qsg(1:ifull,k)*(1.-clcon(:,k))
    qgrg(1:ifull,k) = qgrg(1:ifull,k)*(1.-clcon(:,k))
    cfa(:,k)        = cfa(:,k)*(1.-clcon(:,k))
    qca(:,k)        = qca(:,k)*(1.-clcon(:,k))              
  end where
enddo

if ( nmaxpr==1 .and. mydiag ) then
  write(6,*) 'before newrain',ktau
  write (6,"('t   ',9f8.2/4x,9f8.2)") t(idjd,:)
  write (6,"('qv  ',9f8.3/4x,9f8.3)") qg(idjd,:)
  write (6,"('qf  ',9f8.3/4x,9f8.3)") qfg(idjd,:)
  write (6,"('ql  ',9f8.3/4x,9f8.3)") qlg(idjd,:)
  write (6,"('qr  ',9f8.3/4x,9f8.3)") qrg(idjd,:)
  write (6,"('qs  ',9f8.3/4x,9f8.3)") qsg(idjd,:)
  write (6,"('qg  ',9f8.3/4x,9f8.3)") qgrg(idjd,:)
endif
if ( diag ) then
  call maxmin(t,' t',ktau,1.,kl)
  call maxmin(qg,'qv',ktau,1.e3,kl)
  call maxmin(qfg,'qf',ktau,1.e3,kl)
  call maxmin(qlg,'ql',ktau,1.e3,kl)
  call maxmin(qrg,'qr',ktau,1.e3,kl)
  call maxmin(qsg,'qs',ktau,1.e3,kl)
  call maxmin(qgrg,'qg',ktau,1.e3,kl)
endif

! Add convective cloud water into fields for radiation
! cfrad replaced by updating cfrac Oct 2005
! Moved up 16/1/06 (and ccov,cfrac NOT UPDATED in newrain)
! done because sometimes newrain drops out all qlg, ending up with 
! zero cloud (although it will be rediagnosed as 1 next timestep)
do k = 1,kl
  fl         = max(0., min(1., (t(1:ifull,k)-ticon)/(273.15-ticon) ) )
  qlrad(:,k) = qlg(1:ifull,k)+qrg(1:ifull,k)+fl*qccon(:,k)
  qfrad(:,k) = qfg(1:ifull,k)+qsg(1:ifull,k)+qgrg(1:ifull,k)+(1.-fl)*qccon(:,k)
enddo

!     Calculate precipitation and related processes
call newicerain(land,dt,rhoa,dz,prf,cdso4,cfa,qca,t,qlg,qfg,qrg,                 &
                precs,qg,cfrac,rfrac,ccov,preci,precg,qevap,qsubl,qauto,qcoll,   &
                qaccr,fluxr,fluxi,fluxs,fluxg,fluxmelt,pfstayice,pfstayliq,      &
                pqfsed,slopes,prscav)

if ( nmaxpr==1 .and. mydiag ) then
  write(6,*) 'after newrain',ktau
  write (6,"('t   ',9f8.2/4x,9f8.2)") t(idjd,:)
  write (6,"('qv  ',9f8.3/4x,9f8.3)") qg(idjd,:)
  write (6,"('qf  ',9f8.3/4x,9f8.3)") qfg(idjd,:)
  write (6,"('ql  ',9f8.3/4x,9f8.3)") qlg(idjd,:)
  write (6,"('qr  ',9f8.3/4x,9f8.3)") qrg(idjd,:)
  write (6,"('qs  ',9f8.3/4x,9f8.3)") qsg(idjd,:)
  write (6,"('qg  ',9f8.3/4x,9f8.3)") qgrg(idjd,:)
end if
if ( diag ) then
  call maxmin(t,' t',ktau,1.,kl)
  call maxmin(qg,'qv',ktau,1.e3,kl)
  call maxmin(qfg,'qf',ktau,1.e3,kl)
  call maxmin(qlg,'ql',ktau,1.e3,kl)
  call maxmin(qrg,'qr',ktau,1.e3,kl)
  call maxmin(qsg,'qs',ktau,1.e3,kl)
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

 subroutine newcloud(tdt,land,prf,kbase,ktop,rhoa,cdrop,ttg,qtg,qlg,qfg,cfrac,ccov,cfa,qca)

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
integer, dimension(ifull), intent(in) :: kbase, ktop

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

where(ttg>=tfrz)
  fice=0.
elsewhere(ttg>=tice.and.qfg(1:ifull,:)>1.e-12)
  fice=min(qfg(1:ifull,:)/(qfg(1:ifull,:)+qlg(1:ifull,:)),1.)
elsewhere(ttg>=tice)
  fice=0.
elsewhere
  fice=1.
end where
qcg(:,:)=qlg(1:ifull,:)+qfg(1:ifull,:)
qcold(:,:)=qcg(:,:)
qfnew=fice(:,:)*qcg(:,:)
ttg(1:ifull,:)=ttg(1:ifull,:)+hlfcp*(qfnew-qfg(1:ifull,:)) !Release L.H. of fusion
qfg(1:ifull,:)=fice(:,:)*qcg(:,:)
qlg(1:ifull,:)=max(0.,qcg(:,:)-qfg(1:ifull,:))

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
if ( nclddia==-3 ) then
  do k=1,kl
    where ( land(1:ifull) )
      rcrit(:,k)=max( rcrit_l, (1.-16.*(1.-sig(k))**3) )
    elsewhere
      rcrit(:,k)=max( rcrit_s, (1.-16.*(1.-sig(k))**3) )
    end where
  enddo
else if ( nclddia<0 ) then
  do k=1,kl
    where ( land(1:ifull) )
      rcrit(:,k)=max( rcrit_l, (1.-4.*(1.-sig(k))**2) )
    elsewhere
      rcrit(:,k)=max( rcrit_s, (1.-4.*(1.-sig(k))**2) )
    end where
  enddo
else if ( nclddia==1 ) then
  do k=1,kl
    where ( land(1:ifull) )
      rcrit(:,k)=max( rcrit_l, sig(k)**3 )
    elsewhere
      rcrit(:,k)=max( rcrit_s, sig(k)**3 )
    end where
  enddo
else if ( nclddia==2 ) then
  do k=1,kl
    where ( land(1:ifull) )
      rcrit(:,k)=rcrit_l
    elsewhere
      rcrit(:,k)=rcrit_s
    end where
  enddo
else if ( nclddia==3 ) then
  do k=1,kl
    where ( land(1:ifull) )
      rcrit(:,k)=max( rcrit_l, sig(k)**2 )          ! .75 for R21 Mk2
    elsewhere
      rcrit(:,k)=max( rcrit_s, sig(k)**2 )          ! .85 for R21 Mk2
    end where
  enddo
else if ( nclddia==4 ) then
  do k=1,kl
    where ( land(1:ifull) )
      rcrit(:,k)=max( rcrit_l, sig(k) )             ! .75 for test Mk2/3
    elsewhere
      rcrit(:,k)=max( rcrit_s, sig(k) )             ! .9  for test Mk2/3
    end where
  enddo
else if ( nclddia==5 ) then  ! default till May 08
  do k=1,kl
    where ( land(1:ifull) )
      rcrit(:,k)=max( rcrit_l, min(.99,sig(k)) )    ! .75 for same as T63
    elsewhere
      rcrit(:,k)=max( rcrit_s, min(.99,sig(k)) )    ! .85 for same as T63
    end where
  enddo
else if ( nclddia==6 ) then
  do k=1,kl
    where ( land(1:ifull) )
      rcrit(:,k)=max(rcrit_l*(1.-.15*sig(k)),sig(k)**4)
    elsewhere
      rcrit(:,k)=max(rcrit_s*(1.-.15*sig(k)),sig(k)**4)
    end where
  enddo
else if ( nclddia==7 ) then
  do k=1,kl
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
  do k=1,kl  ! typically set rcrit_l=.75,  rcrit_s=.85
    tk(:)=ds/(em(1:ifull)*208498.) ! MJT suggestion
    fl(:)=(1.+real(nclddia))*tk(:)/(1.+real(nclddia)*tk(:))
    ! for rcit_l=.75 & nclddia=12 get rcrit=(0.751, 0.769, .799, .901, .940, .972, .985) for (200, 100, 50, 10, 5, 2, 1) km
    where ( land(1:ifull) )
      rcrit(:,k)=max(1.-fl*(1.-rcrit_l),sig(k)**3)        
    elsewhere
      rcrit(:,k)=max(1.-fl*(1.-rcrit_s),sig(k)**3)         
    end where
  end do
end if  ! (nclddia<0)  .. else ..


if ( ncloud<=3 ) then
  ! usual diagnostic cloud fraction
      
  ! Calculate cloudy fraction of grid box (cfrac) and gridbox-mean cloud water
  ! using the triangular PDF of Smith (1990)

  do k=1,kl
    do mg=1,ifull
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
    write(6,*) 'qsi',qsi(mg,k)
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
  where( ttg(1:ifull,:)>=Tice )
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
  qtot(:,:)=qtg(:,:)+qcg(:,:)
  tliq(:,:)=ttg(:,:)-hlcp*qcg(:,:)-hlfcp*qfg(1:ifull,:)
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

  where(ttg>=Tice)
    qfg(1:ifull,:) = fice*qcg
    qlg(1:ifull,:) = qcg - qfg(1:ifull,:)
  elsewhere
    qfg(1:ifull,:) = qcg
    qlg(1:ifull,:) = 0.
    qcg(1:ifull,:) = qfg(1:ifull,:)
  end where
  
end if ! ncloud<4 ..else..


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
qtg(1:ifull,:)=qtot(1:ifull,:)-qcg(1:ifull,:)
ttg(1:ifull,:)=tliq(1:ifull,:)+hlcp*qcg(1:ifull,:)+hlfcp*qfg(1:ifull,:)
ccov(1:ifull,:)=cfrac(1:ifull,:) !Do this for now

! Vertically sub-grid cloud
where ( cfrac(1:ifull,2:kl-1)>1.e-2.and.cfrac(1:ifull,3:kl)==0..and.cfrac(1:ifull,1:kl-2)==0. )
  ccov(1:ifull,2:kl-1)=sqrt(cfrac(1:ifull,2:kl-1))
end where
     
if(diag.and.mydiag)then
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
!      land - logical variable for surface type ( = T for land points)
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
!      precs - amount of stratiform precipitation in timestep (mm)
!      qtg - water vapour mixing ratio (kg/kg) - called qg in C-CAM
!      cfrac - stratiform cloud fraction
!      ccov - stratiform cloud *cover* looking from above (currently = cfrac)
!
! Output:
!
! from arguments
!      preci - amount of stratiform snowfall in timestep (mm)
!      precg - amount of stratiform grauple in timestep (mm)
!      qevap - evaporation of rainfall (kg/kg)
!      qsubl - sublimation of snowfall (kg/kg)
!      qauto - autoconversion of cloud liquid water (kg/kg)
!      qcoll - collection by rain of cloud liquid water (kg/kg)
!      qaccr - accretion by snow of cloud liquid water (kg/kg)
!
!**************************************************************************

subroutine newicerain(land,tdt,rhoa,dz,prf,cdrop,cfa,qca,ttg,qlg,qfg,qrg,precs,qtg,cfrac,cffall,ccov,              &
                      preci,precg,qevap,qsubl,qauto,qcoll,qaccr,fluxr,fluxi,fluxs,fluxg,fluxmelt,pfstayice,        &
                      pfstayliq,pqfsed,slopes,prscav)

use cc_mpi, only : mydiag
use estab, only : esdiffx, qsati, pow75
use kuocomb_m
use morepbl_m  !condx  

implicit none

! Global parameters
include 'newmpar.h'
include 'const_phys.h' !Input physical constants
include 'cparams.h'    !Input cloud scheme parameters
include 'kuocom.h'     !acon,bcon,Rcm,ktsav,nevapls
include 'parm.h'

! Argument list
logical, dimension(ifull), intent(in) :: land
real, intent(in) :: tdt
real, dimension(ifull,kl), intent(in) :: rhoa
real, dimension(ifull,kl), intent(in) :: dz
real, dimension(ifull,kl), intent(in) :: prf
real, dimension(ifull,kl), intent(in) :: cdrop
real, dimension(ifull+iextra,kl), intent(inout) :: ttg
real, dimension(ifull+iextra,kl), intent(inout) :: qlg
real, dimension(ifull+iextra,kl), intent(inout) :: qfg
real, dimension(ifull+iextra,kl), intent(inout) :: qrg
real, dimension(ifull+iextra,kl), intent(inout) :: qtg
real, dimension(ifull), intent(inout) :: precs
real, dimension(ifull), intent(inout) :: preci
real, dimension(ifull), intent(inout) :: precg
real, dimension(ifull,kl), intent(inout) :: cfrac
real, dimension(ifull+iextra,kl), intent(inout) :: cffall
real, dimension(ifull,kl), intent(in) :: ccov
real, dimension(ifull,kl), intent(out) :: qevap
real, dimension(ifull,kl), intent(out) :: qsubl
real, dimension(ifull,kl), intent(out) :: qauto
real, dimension(ifull,kl), intent(out) :: qcoll
real, dimension(ifull,kl), intent(out) :: qaccr
real, dimension(ifull,kl), intent(out) :: pqfsed
real, dimension(ifull,kl), intent(out) :: pfstayice
real, dimension(ifull,kl), intent(out) :: pfstayliq
real, dimension(ifull,kl), intent(out) :: slopes
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
real, dimension(ifull,kl-1) :: rhor,gam
real, dimension(ifull,kl) :: clfr,cifr,qsatg,cfrain,cfmelt,fluxauto
real, dimension(ifull,kl) :: rhoi,vi2,vl2
real, dimension(ifull) :: clfra,mxclfrliq,rdclfrliq,fluxrain
real, dimension(ifull) :: mxclfrice,rdclfrice,rica,fluxice,cifra
real, dimension(ifull) :: sublflux,fsclr,caccr,dqf,qif,dttg,csb,bf
real, dimension(ifull) :: dqs,ql,cdt,rhoiin,cffluxin,rhoiout
real, dimension(ifull) :: cffluxout,rhodz,evap,qpf,clrevap,fr
real, dimension(ifull) :: mxovr,rdovr,fcol,coll,alph,rhorin
real, dimension(ifull) :: rhorout,tc,alphaf,tk,pk,es,aprpr,bprpr
real, dimension(ifull) :: curly,frclr,Csbsav

integer k,mg

real apr,bpr,bl,cev,crate,dql,dqsdt,frb,qcic,qcrit,ql1,ql2,R6c,R3c,beta6,eps
real satevap,selfcoll,Wliq,cfla,dqla,qla,viin,qsl

!real, parameter :: normg = 5026548245.74367
!real, parameter :: norms = 942477796.076938
!real, parameter :: vcong = 87.2382675
!real, parameter :: vcons = 6.6280504

do k=1,kl
  fluxr(1:ifull,k)=0.
  fluxi(1:ifull,k)=0.
  fluxs(1:ifull,k)=0.
  fluxg(1:ifull,k)=0.
  fluxmelt(1:ifull,k)=0.  
  fluxauto(1:ifull,k)=0.
  qevap(1:ifull,k)=0.
  qauto(1:ifull,k)=0.
  qcoll(1:ifull,k)=0.
  cfrain(1:ifull,k)=0.
  pk(1:ifull)=100.0*prf(1:ifull,k)
  qsatg(1:ifull,k)=qsati(pk(1:ifull),ttg(1:ifull,k))
  cifr(1:ifull,k)=cfrac(1:ifull,k)*qfg(1:ifull,k)/max(qlg(1:ifull,k)+qfg(1:ifull,k),1.E-30)
  clfr(1:ifull,k)=cfrac(1:ifull,k)*qlg(1:ifull,k)/max(qlg(1:ifull,k)+qfg(1:ifull,k),1.E-30)
end do

!**************** Cut here to insert new auto scheme ********************            
if ( ncloud>0 .and. ncloud<4 ) then

  ! Using new (subgrid) autoconv scheme... 
  do k = kl-1,1,-1
    do mg = 1,ifull
      cfrain(mg,k)=0.0
      rhodz(mg)=rhoa(mg,k)*dz(mg,k)
      if ( clfr(mg,k)>0. ) then
        ql=qlg(mg,k)
        cfla=0.
        dqla=0.
        if ( cfa(mg,k)>0. ) then
          cfla=cfa(mg,k)*clfr(mg,k)/(clfr(mg,k)+cifr(mg,k))
          qla=qca(mg,k)/cfa(mg,k)
          ! Following few lines are for Yangang Liu's new scheme (2004: JAS, GRL)
          Wliq = max(1.e-10, 1000. * qla * rhoa(mg,k)) !g/m3
          R6c = 4.09e-4 * ( 1.15e23*1.e-6*cdrop(mg,k) / Wliq**2 )**(1./6.)
          eps = 1. - 0.7 * exp(-0.003e-6*cdrop(mg,k)) !mid range
          beta6 = ((1.+3.*eps**2)*(1.+4.*eps**2)*(1.+5.*eps**2) / ((1.+eps**2)*(1.+2.*eps**2)) )**(1./6.)
          R3c = 1.e-6*R6c/beta6 !in metres
          qcrit=(4.*pi/3.)*rhow*R3c**3*Cdrop(mg,k)/rhoa(mg,k) !New qcrit
          if ( qla<=qcrit ) then
            ql2=qla
          else
            ! Following is Liu & Daum (JAS, 2004)
            Crate=1.9e17*(0.75*rhoa(mg,k)/(pi*rhow))**2*beta6**6/cdrop(mg,k)
            ql1=qla/sqrt(1.+2.*crate*qla**2*tdt)
            ql1=max(ql1, qcrit) !Intermediate qlg after auto
            Frb=dz(mg,k)*rhoa(mg,k)*(qla-ql1)/tdt
            cdt(mg)=tdt*0.5*Ecol*0.24*pow75(Frb)
            selfcoll=min(ql1,ql1*cdt(mg))
            ql2=ql1-selfcoll
            cfrain(mg,k)=cfla
          end if
          dqla=cfla*(qla-ql2)
          ql(mg)=max(1.e-10,qlg(mg,k)-dqla)
        end if
        dql=qlg(mg,k)-ql(mg)
        qauto(mg,k)=qauto(mg,k)+dql
        qlg(mg,k)=qlg(mg,k)-dql
        fluxauto(mg,k)=dql*rhodz(mg)
      end if
    end do
  end do

! Or, using old autoconv scheme... also used by prognostic cloud scheme
else

  do k=kl-1,1,-1
    do mg=1,ifull
      cfrain(mg,k)=0.0
      rhodz(mg)=rhoa(mg,k)*dz(mg,k)
      if ( clfr(mg,k)>0. ) then
        qcrit=(4.*pi/3.)*rhow*Rcm**3*cdrop(mg,k)/rhoa(mg,k)
        qcic=qlg(mg,k)/clfr(mg,k) !In cloud value
        if ( qcic<qcrit ) then
          ql=qlg(mg,k)
        else
          Crate=Aurate*rhoa(mg,k)*(rhoa(mg,k)/(cdrop(mg,k)*rhow))**(1./3.)
          ql1=1./pow75(qcic**(-4./3.)+(4./3.)*Crate*tdt)
          ql1=max(ql1, qcrit) !Intermediate qlg after auto
          Frb=dz(mg,k)*rhoa(mg,k)*(qcic-ql1)/tdt
          cdt(mg)=tdt*0.5*Ecol*0.24*pow75(Frb) !old
          selfcoll=min(ql1,ql1*cdt(mg))
          ql2=ql1-selfcoll
          ql(mg)=clfr(mg,k)*ql2
          cfrain(mg,k)=clfr(mg,k)
        end if
        dql=qlg(mg,k)-ql(mg)
        qauto(mg,k)=qauto(mg,k)+dql
        qlg(mg,k)=qlg(mg,k)-dql
        fluxauto(mg,k)=dql*rhodz(mg)
      end if
    end do
  end do

end if ! ( ncloud>0 .and. ncloud<4 )

! Set up prognostic rain - MJT
! The following has been modified according to LDR's flux divergence calculation.
! LDR's original scheme can be recovered by setting foutliq=1 and fthruliq=1 or using
! ncloud<=1.

! combine autoconversion and prognostic rain
qauto(1:ifull,1:kl-1)=qauto(1:ifull,1:kl-1)+qrg(1:ifull,1:kl-1)
rhor(1:ifull,1:kl-1)=qrg(1:ifull,1:kl-1)*rhoa(1:ifull,1:kl-1)
! max overlap autoconversion and rain from previous time step
cfrain(1:ifull,1:kl-1)=max(cfrain(1:ifull,1:kl-1),cffall(1:ifull,1:kl-1)) 

! Set up frozen fields
! Convert from mixing ratio to density of ice, and work out ice cloud fraction
qfg(1:ifull,1:kl)=max( qfg(1:ifull,1:kl), 0. )  
rhoi(1:ifull,1:kl)=rhoa(1:ifull,1:kl)*qfg(1:ifull,1:kl) 
cifr(1:ifull,1:kl)=cfrac(1:ifull,1:kl)*qfg(1:ifull,1:kl)/max( qlg(1:ifull,1:kl)+qfg(1:ifull,1:kl), 1.e-30 )
clfr(1:ifull,1:kl)=max( cfrac(1:ifull,1:kl)-cifr(1:ifull,1:kl), 0. )
qsubl(1:ifull,1:kl)=0.
qaccr(1:ifull,1:kl)=0.
cfmelt(1:ifull,1:kl)=0.

if ( diag .and. mydiag ) then
  write(6,*) 'cfrac  ',cfrac(idjd,:)
  write(6,*) 'cifr   ',cifr(idjd,:)
  write(6,*) 'clfr   ',clfr(idjd,:)
  write(6,*) 'cfrain ',cfrain(idjd,:)
  write(6,*) 'qlg ',qlg(idjd,:)
  write(6,*) 'qfg ',qfg(idjd,:)
  write(6,*) 'qrg ',qrg(idjd,:)
endif  ! (diag.and.mydiag)

! The following has been modified to track the random overlap rain fraction (rdclfrliq)
! and the max overlap rain fraction (mxclfrliq) so than both random overlaped and
! maximum/random overlaped clouds are supported - MJT
vl2(1:ifull,kl)=0.
prscav(1:ifull,kl)=0.
pfstayliq(1:ifull,kl)=0.
clfra(1:ifull)=1.e-6
fluxrain(1:ifull)=0.
mxclfrliq(1:ifull)=0. ! max overlap rain fraction
rdclfrliq(1:ifull)=0. ! rnd overlap rain fraction

! The following has been modified to track the random overlap ice fraction (rdclfrice)
! and the max overlap ice fraction (mxclfrice) so than both random overlapped and
! max/random overlapped clouds are supported - MJT
vi2(1:ifull,kl)=0.1 ! Assume no cloud at top level
slopes(1:ifull,kl)=0.
pqfsed(1:ifull,kl)=0.
pfstayice(1:ifull,kl)=0.
fluxice(1:ifull)=0.
cifra(1:ifull)=0.
rica(1:ifull)=0.
mxclfrice(1:ifull)=0. ! max overlap ice fraction
rdclfrice(1:ifull)=0. ! rnd overlap ice fraction

! Now work down through the levels...
do k=kl-1,1,-1

  pk(:)=100.*prf(:,k)
    
  ! Frozen ----------------------------------------------------------------------------
  tc(:)=ttg(1:ifull,k)-tfrz
  slopes(:,k)=1.6e3*10**(-0.023*tc(:))
  alphaf(:)=hls*qsatg(1:ifull,k)/(rvap*ttg(1:ifull,k)**2)
  gam(:,k)=hlscp*alphaf(:) !(L/cp)*dqsdt (HBG notation)
  sublflux(:)=0.
  caccr(:)=0.
  dqf(:)=0.  

  ! Set up the Rate constant for snow sublimation
  ! MJT notes - curly and Csbsav depend on vi2(:,k+1), so vi2(:,k) can be updated below
  Tk(:)=ttg(:,k)
  es(:)=qsatg(1:ifull,k)*pk(:)/epsil
  Aprpr(:)=(hls/(rKa*Tk(:)))*(hls/(rvap*Tk(:))-1.)
  Bprpr(:)=rvap*Tk(:)/((Dva/pk(:))*es(:))
  where ( nevapls==-1 .or. (nevapls==-2.and.condx(:)>0..and.k<=ktsav(:)) )
    curly=0.
  elsewhere
    curly(:)=0.65*slopes(:,k)**2+0.493*slopes(:,k)*sqrt(slopes(:,k)*vi2(:,k+1)*rhoa(:,k)/um)
  end where

  ! Define the rate constant for sublimation of snow, omitting factor rhoi
  Csbsav(:)=4.*curly(:)/(rhoa(:,k)*qsatg(1:ifull,k)*(Aprpr(:)+Bprpr(:))*pi*vi2(:,k+1)*rhosno)
    
  ! The following flag detects max/random overlap clouds
  ! that are separated by a clear layer
  where ( cifr(:,k)<1.e-10.or.nmr==0 )
    ! combine max overlap from last cloud with net random overlap
    rdclfrice(:)=rdclfrice(:)+mxclfrice(:)-rdclfrice(:)*mxclfrice(:)
    mxclfrice(:)=0.
  end where

  ! Melt falling ice if > 0 deg C
  where ( ttg(1:ifull,k)>tfrz.and.fluxice(:)>0. )
    qif(:)=fluxice(:)/(dz(:,k)*rhoa(:,k))      !Mixing ratio of ice
    dttg(:)=-hlfcp*qif(:)
    ttg(1:ifull,k)=ttg(1:ifull,k)+dttg(:)
    qsatg(1:ifull,k)=qsatg(1:ifull,k)+gam(:,k)*dttg(:)/hlscp
    fluxmelt(:,k)=fluxmelt(:,k)+fluxice(:)
    cfmelt(:,k)=cifra(:)
    fluxice(:)=0.
    cifra(:)=0.
    rdclfrice(:)=0.
    mxclfrice(:)=0.
  end where

  ! Compute the sublimation of ice falling from level k+1 into level k
  fsclr(:)=(1.-cifr(:,k)-clfr(:,k))*fluxice(:)
  where ( fluxice(:)>0..and.qtg(1:ifull,k)<qsatg(1:ifull,k) ) ! sublime snow
    Csb(:)=Csbsav(:)*fluxice(:)/tdt !LDR
    bf(:)=1.+0.5*Csb(:)*tdt*(1.+gam(:,k))
    dqs(:)=max(0.,tdt*(Csb(:)/bf(:))*(qsatg(1:ifull,k)-qtg(1:ifull,k)))
    dqs(:)=min(dqs(:),(qsatg(1:ifull,k)-qtg(1:ifull,k))/(1.+gam(:,k))) !Don't supersat.
    sublflux(:)=min(dqs(:)*rhoa(:,k)*dz(:,k),fsclr(:))
    fluxice(:)=fluxice(:)-sublflux(:)
    fsclr(:)=fsclr(:)-sublflux(:)
    dqs(:)=sublflux(:)/(rhoa(:,k)*dz(:,k))
    qsubl(:,k)=qsubl(:,k)+dqs(:)
    qtg(1:ifull,k)=qtg(1:ifull,k)+dqs(:)
    dttg(:)=-hlscp*dqs(:)
    ttg(1:ifull,k)=ttg(1:ifull,k)+dttg(:)
    qsatg(1:ifull,k)=qsatg(1:ifull,k)+gam(:,k)*dttg(:)/hlscp
  end where

  ! Accretion of cloud water by falling ice
  ! This calculation uses the incoming fluxice without subtracting sublimation
  ! (since subl occurs only outside cloud), so add sublflux back to fluxice.
  where ( fluxice(:)+sublflux(:)>0.and.qlg(1:ifull,k)>1.e-10 )
    ql(:)=qlg(1:ifull,k)
    cdt(:)=Eac*slopes(:,k)*(fluxice(:)+sublflux(:))/(2.*rhosno)
    dqf(:)=min(ql(:),cifra(:)*ql(:),ql(:)*cdt(:)/(1.+0.5*cdt(:)))
    clfr(:,k)=clfr(:,k)*(1.-dqf(:)/qlg(1:ifull,k))
    caccr(:)=clfr(:,k)*dqf(:)/qlg(1:ifull,k)
    qlg(1:ifull,k)=qlg(1:ifull,k)-dqf(:)
    qaccr(:,k)=qaccr(:,k)+dqf(:)
    fluxice(:)=fluxice(:)+rhoa(:,k)*dz(:,k)*dqf(:)
    dttg(:)=hlfcp*dqf(:)
    ttg(1:ifull,k)=ttg(1:ifull,k)+dttg(:)
    qsatg(1:ifull,k)=qsatg(1:ifull,k)+gam(:,k)*dttg(:)/hlscp
  end where

  where ( fsclr(:)<1.e-15 )
    rdclfrice(:)=0.
    mxclfrice(:)=0.
  end where
  mxclfrice(:)=max(mxclfrice(:),cifr(:,k)+caccr(:))                      !max overlap
  cifra(:)=max(0.01,mxclfrice(:)+rdclfrice(:)-mxclfrice(:)*rdclfrice(:)) !rnd overlap the mx and rd ice fractions

  ! Set up ice fall speed field
  ! MJT notes - currently no feedback as rhoi is updated below
  select case(abs(ldr))
    case(1)  ! 1 for R21 runs, like prev lw=22
      where ( cifr(1:ifull,k)>=1.e-10 )
        vi2(1:ifull,k)=3.23*(rhoi(1:ifull,k)/cifr(1:ifull,k))**0.17
      elsewhere
        vi2(1:ifull,k)=vi2(1:ifull,k+1)
      end where
    case(2)
      vi2(:,k)=vi2(:,k+1)
      where ( cifr(:,k)>=1.e-10 )
        vi2(:,k)=0.9*3.23*(rhoi(:,k)/cifr(:,k))**0.17
      end where
    case(3)
      vi2(1:ifull,k)=vi2(1:ifull,k+1)
      where ( cifr(1:ifull,k)>=1.e-10 )
        vi2(1:ifull,k)=max(0.1,2.05+0.35*log10(qfg(1:ifull,k)/cifr(1:ifull,k)))
      end where
    case(4)
      vi2(1:ifull,k)=vi2(1:ifull,k+1)
      where ( cifr(1:ifull,k)>=1.e-10 )
        vi2(1:ifull,k)=1.4*3.23*(rhoi(1:ifull,k)/cifr(1:ifull,k))**0.17
      end where
    case(11) ! 1 for R21 runs, like prev lw=22
      ! following are alternative slightly-different versions of above
      ! used for I runs from 29/4/05 till 30/8/05
      ! for given qfg, large cifr implies small ice crystals, 
      ! with a small fall speed. 
      ! Note that for very small qfg, cifr is small.
      ! But rhoi is like qfg, so ratio should also be small and OK.         
      vi2(1:ifull,k)=max( vi2(1:ifull,k+1),3.23*(rhoi(1:ifull,k)/max(cifr(1:ifull,k),1.e-30))**0.17 )
    case(22)
      vi2(1:ifull,k)=max( vi2(1:ifull,k+1),.9*3.23*(rhoi(1:ifull,k)/max(cifr(1:ifull,k),1.e-30))**0.17 )
    case(33)
      ! following max gives vi2=.1 for qfg=cifr=0
      vi2(1:ifull,k)=max( vi2(1:ifull,k+1),2.05 +0.35*log10(max(qfg(1:ifull,k),2.68e-36)/max(cifr(1:ifull,k),1.e-30)) )
  end select

  ! Snow fall speed (from Lin et al 1983 - see GFDL AM3)
  !vs2(:,k)=max(0.1, vcons*rhof*(qs(k)*den(k)/norms)**0.0625)
    
  ! Grauple fall speed (from Lin et al 1983 - see GFDL AM3)
  !vg2(:,k)=max(0.1, vg_fac*vcong*rhof*(qg(k)*den(k)/normg)**0.125)  
    
  ! Set up the parameters for the flux-divergence calculation
  alph(:)=tdt*vi2(:,k)/dz(:,k)
  foutice(:,k)=1.-exp(-alph(:))          !analytical
  fthruice(:,k)=1.-foutice(:,k)/alph(:)  !analytical

  ! Save sedimentation rate for aerosol scheme
  pqfsed(:,k)=foutice(:,k)

  ! Save this for the wet deposition scheme.
  pfstayice(:,k)=fluxice(:)*(1.-fthruice(:,k))/tdt !Flux staying in layer k
  
  ! Compute fluxes into the box
  where ( fluxice(:)>0. )
    rhoiin(:)=fluxice(:)/dz(:,k)
    cffluxin(:)=min(1.,rhoiin(:)/rica(:))
  elsewhere
    rhoiin(:)=0.
    cffluxin(:)=0.
  end where

  ! Compute the fluxes of ice and cloud amount leaving the box
  where ( cifr(:,k)>=1.e-10 )
    ! Use the flux-divergent form as in Rotstayn (QJRMS, 1997)
    rhoiout(:)=rhoi(:,k)*foutice(:,k)
    cffluxout(:)=cifr(:,k)*foutice(:,k)
    rica(:)=rhoi(:,k)/cifr(:,k) !in cloud rhoi above
  elsewhere
    ! Keep value of rica from above
    rhoiout(:)=0.
    cffluxout(:)=0.
  end where
            
  ! Update the rhoi and cifr fields
  cifr(:,k)=min(1.-clfr(:,k),(cifr(:,k)-cffluxout(:))+cffluxin(:)*(1.-fthruice(:,k)))
  rhoi(:,k)=rhoi(:,k)-rhoiout(:)+rhoiin(:)*(1.-fthruice(:,k))
  fluxice(:)=rhoiout(:)*dz(:,k)+fluxice(:)*fthruice(:,k) 
  ! Now fluxice is flux leaving layer k
  fluxi(:,k)=fluxi(:,k)+fluxice(:)


  ! Liquid ----------------------------------------------------------------------------
  rhodz(:)=rhoa(:,k)*dz(:,k)
  evap(:)=0.

  ! The following flag detects maximum/random overlap clouds
  ! that are separated by a clear layer
  where ( (clfr(:,k)<1.e-10.and.cfrain(1:ifull,k)<1.e-10).or.nmr==0 )
    ! combine max overlap from above cloud with net random overlap
    rdclfrliq(:)=rdclfrliq(:)+mxclfrliq(:)-rdclfrliq(:)*mxclfrliq(:)
    mxclfrliq(:)=0.
  end where

  ! Add flux of melted snow to fluxrain
  fluxrain(:)=fluxrain(:)+fluxmelt(:,k)+fluxauto(:,k)

  ! Evaporation of rain
  qpf(:)=fluxrain(:)/rhodz(:) !Mix ratio of rain which falls into layer  ! MJT suggestion
  clrevap(:)=(1.-clfr(:,k))*qpf(:)                                       ! MJT suggestion
  do mg=1,ifull
    if ( fluxrain(mg)>0. ) then
      qsatg(mg,k)=qsati(pk(mg),ttg(mg,k))
      if ( ttg(mg,k)<tfrz .and. ttg(mg,k)>=tice ) then
        qsl=qsatg(mg,k)+epsil*esdiffx(ttg(mg,k))/(100.*prf(mg,k))    ! MJT suggestion
      else
        qsl=qsatg(mg,k)
      end if             !qsl is qs value over liquid surface
      Tk(mg)=ttg(mg,k)
      es(mg)=qsl*pk(mg)/epsil 
      Apr=(hl/(rKa*Tk(mg)))*(hl/(rvap*Tk(mg))-1.)
      Bpr=rvap*Tk(mg)/((Dva/pk(mg))*es(mg))
      Fr(mg)=fluxrain(mg)/tdt/clfra(mg)
      Cev=clfra(mg)*3.8e2*sqrt(Fr(mg)/rhoa(mg,k))/(qsl*(Apr+Bpr))
      dqsdt=hl*qsl/(rvap*ttg(mg,k)**2)
      bl=1.+0.5*Cev*tdt*(1.+hlcp*dqsdt)
      evap(mg)=tdt*(Cev/bl)*(qsl-qtg(mg,k))
      satevap=(qsl-qtg(mg,k))/(1.+hlcp*dqsdt) !Evap to saturate
!     Vl2=11.3*Fr**(1./9.)/sqrt(rhoa(mg,k))  !Actual fall speed
!     Vl2=5./sqrt(rhoa(mg,k))                !Nominal fall speed

      evap(mg)=max(0., min(evap(mg),satevap,qpf(mg),clrevap(mg)))
      if(nevapls==-1)evap(mg)=0.
      if(nevapls==-2.and.k<=ktsav(mg).and.condx(mg)>0.)evap(mg)=0.
      if(nevapls==-3)evap(mg)=.5*evap(mg)
      if(nevapls==-4.and.k<=ktsav(mg).and.condx(mg)>0.)evap(mg)=.5*evap(mg) ! usual
      qevap(mg,k)=qevap(mg,k)+evap(mg)
      qtg(mg,k)=qtg(mg,k)+evap(mg)
      ttg(mg,k)=ttg(mg,k)-hlcp*evap(mg)
    end if
  end do
  frclr(:)=rhodz(:)*(clrevap(:)-evap(:)) !over tdt ! MJT suggestion

  ! Now do the collection term.
  where ( fluxrain(:)>0. )
    Fr(:)=fluxrain(:)/clfra(:)/tdt
    mxovr(:)=min(mxclfrliq(:),clfr(:,k))                      ! max overlap
    mxovr(:)=max(cfrain(1:ifull,k),mxovr(:))
    rdovr(:)=rdclfrliq(:)*clfr(:,k)                           ! rnd overlap
    cfrain(1:ifull,k)=mxovr(:)+rdovr(:)-mxovr(:)*rdovr(:)     ! combine collection
  elsewhere
    Fr(:)=0.
  end where

  ! The collection term comprises collection by stratiform rain falling from
  ! above (Fr), stratiform rain released in this grid box (Frb).
  ! Frb term now done above.
  fcol(:)=min(1.,mxclfrliq(:)/(1.e-20+clfr(:,k)))          !max overlap
  fcol(:)=fcol(:)+rdclfrliq(:)-fcol(:)*rdclfrliq(:)        !rnd overlap
  cdt(:)=tdt*Ecol*0.24*fcol(:)*pow75(Fr(:))
  prscav(:,k)=tdt*0.24*fcol(:)*pow75(Fr(:))                !Strat only

  coll(:)=min(qlg(1:ifull,k),qlg(1:ifull,k)*cdt(:)/(1.+0.5*cdt(:)))
  qcoll(:,k)=qcoll(:,k)+coll(:)
  qlg(1:ifull,k)=qlg(1:ifull,k)-coll(:)
  fluxrain(:)=fluxrain(:)+coll(:)*rhodz(:)

  ! subtract evaporated rain
  fluxrain(:)=fluxrain(:)-rhodz(:)*evap(:)
  fluxrain(:)=max(fluxrain(:),0.) !To avoid roundoff -ve's

  ! Calculate the raining cloud cover down to this level, for stratiform (clfra).
  cfrain(1:ifull,k)=min(1.,cfrain(1:ifull,k)+cfmelt(:,k)-cfrain(1:ifull,k)*cfmelt(:,k)) ! rnd overlap
  where ( frclr(:)<1.e-15 )
    rdclfrliq(:)=0.
    mxclfrliq(:)=0.
  end where
  mxclfrliq(:)=max(mxclfrliq(:),cfrain(1:ifull,k)) !max overlap
  clfra(:)=max(1.e-15,rdclfrliq(:)+mxclfrliq(:)-rdclfrliq(:)*mxclfrliq(:)) !rnd overlap the mx and rd rain fractions

  ! Calculate rain fall speed (MJT)
  if ( ncloud>1 ) then
    Fr(:)=fluxrain(:)/tdt/clfra(:)
    vl2(:,k)=11.3*Fr(:)**(1./9.)/sqrt(rhoa(:,k))  !Actual fall speed
    vl2(:,k)=max(vl2(:,k),0.1)
    alph(:)=tdt*vl2(:,k)/dz(:,k)
    foutliq(:,k)=1.-exp(-alph(:))          !analytical
    fthruliq(:,k)=1.-foutliq(:,k)/alph(:)  !analytical
  else
    foutliq(:,k)=1.
    fthruliq(:,k)=1.
  end if
            
  pfstayliq(:,k)=fluxrain(:)*(1.-fthruliq(:,k))/tdt

  ! Compute fluxes into the box
  cffluxin(:)=clfra(:)-cfrain(:,k)
  rhorin(:)=fluxrain(:)/dz(:,k)

  ! Compute the fluxes of rain leaving the box
  ! Use the flux-divergent form as in Rotstayn (QJRMS, 1997)
  cffluxout(:)=cfrain(:,k)*foutliq(:,k)
  rhorout(:)=rhor(:,k)*foutliq(:,k)
            
  ! Update the rhor and cffall fields
  cffall(1:ifull,k)=cfrain(1:ifull,k)-cffluxout(:)+cffluxin(:)*(1.-fthruliq(:,k))
  rhor(:,k)=rhor(:,k)-rhorout(:)+rhorin(:)*(1.-fthruliq(:,k))
  fluxrain(:)=rhorout(:)*dz(:,k)+fluxrain(:)*fthruliq(:,k) 
  ! Now fluxrain is flux leaving layer k
  fluxr(:,k)=fluxr(:,k)+fluxrain(:)
    
end do

! Re-create qrg field
qrg(1:ifull,kl)=0.
qrg(1:ifull,1:kl-1)=rhor(1:ifull,1:kl-1)/rhoa(1:ifull,1:kl-1)

! Re-create qfg field
qfg(1:ifull,1:kl)=rhoi(1:ifull,1:kl)/rhoa(1:ifull,1:kl)

! Factor 0.5 here accounts for leapfrog scheme
precs(1:ifull)=precs(1:ifull)+fluxr(1:ifull,1)+fluxi(1:ifull,1)+fluxs(1:ifull,1)+fluxg(1:ifull,1)
preci(1:ifull)=preci(1:ifull)+fluxi(1:ifull,1)+fluxs(1:ifull,1)
precg(1:ifull)=precg(1:ifull)+fluxg(1:ifull,1)

! Remove small amounts of cloud
where ( qlg(1:ifull,1:kl)<1.e-10.or.clfr(1:ifull,1:kl)<1.e-5 )
  qtg(1:ifull,1:kl)=qtg(1:ifull,1:kl)+qlg(1:ifull,1:kl)
  ttg(1:ifull,1:kl)=ttg(1:ifull,1:kl)-hlcp*qlg(1:ifull,1:kl)
  qlg(1:ifull,1:kl)=0.
  clfr(1:ifull,1:kl)=0.
end where
where ( qfg(1:ifull,1:kl)<1.e-10.or.cifr(1:ifull,1:kl)<1.e-5 )
  qtg(1:ifull,1:kl)=qtg(1:ifull,1:kl)+qfg(1:ifull,1:kl)
  ttg(1:ifull,1:kl)=ttg(1:ifull,1:kl)-hlscp*qfg(1:ifull,1:kl)
  qfg(1:ifull,1:kl)=0.
  cifr(1:ifull,1:kl)=0.
end where
where ( qrg(1:ifull,1:kl)<1.e-10.or.cffall(1:ifull,1:kl)<1.e-5 )
  qtg(1:ifull,1:kl)=qtg(1:ifull,1:kl)+qrg(1:ifull,1:kl)
  ttg(1:ifull,1:kl)=ttg(1:ifull,1:kl)-hlcp*qrg(1:ifull,1:kl)
  qrg(1:ifull,1:kl)=0.
  cffall(1:ifull,1:kl)=0.
end where

!      Adjust cloud fraction (and cloud cover) after precipitation
if(nmaxpr==1.and.mydiag)then
  write(6,*) 'diags from newrain for idjd ',idjd
  write (6,"('cfrac ',9f8.3/6x,9f8.3)") cfrac(idjd,:)
  write (6,"('cftemp',9f8.3/6x,9f8.3)") cifr(idjd,:)+clfr(idjd,:)
  write (6,"('ccov_in',9f8.3/6x,9f8.3)") ccov(idjd,:)
endif

! Diagnostics for debugging
if (diag.and.mydiag) then  ! JLM
  write(6,*) 'vi2',vi2(idjd,:)
  write(6,*) 'cfraci ',cfrac(idjd,:)
  write(6,*) 'cifr',cifr(idjd,:)
  write(6,*) 'clfr',clfr(idjd,:)
  write(6,*) 'ttg',ttg(idjd,:)
  write(6,*) 'qsg',qsatg(idjd,:)         
  write(6,*) 'qlg',qlg(idjd,:)
  write(6,*) 'qfg',qfg(idjd,:)
  write(6,*) 'qsubl',qsubl(idjd,:)
  write(6,*) 'rhoa',rhoa(idjd,:)
  write(6,*) 'rhoi',rhoi(idjd,:)
  write(6,*) 'fluxi ',fluxi(idjd,:)
  write(6,*) 'gam',gam(idjd,:)
  write(6,*) 'foutice',foutice(idjd,:)
  write(6,*) 'fthruice',fthruice(idjd,:)
  write(6,*) 'pqfsed',pqfsed(idjd,:)
  write(6,*) 'cfmelt',cfmelt(idjd,:)
  write(6,*) 'fluxmelt',fluxmelt(idjd,:)
  write(6,*) 'cifra,fluxice,rica',cifra(idjd),fluxice(idjd),rica(idjd)
endif  ! (diag.and.mydiag)

return
end subroutine newicerain

end module leoncld_mod
