! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2022 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
module cloudmod

! This module is the Rotstayn 1997 cloud microphysics parameterisation for cloud fraction.
! prognostic cloud fraction scheme based on Tiedtke from GFDL-CM3.

! ncloud = 0    Standard LDR cloud microphysics with water vapour, liquid cloud and ice cloud
! ncloud = 2    Same as ncloud=0, but with prognostic rain and modified cfrac
! ncloud = 3    Same as ncloud=2, but with prognostic graupel and snow, as well as modified cfrac
! ncloud = 4    Use prognostic cloud fraction based on Tiedtke from GFDL-CM3
! ncloud = 10   Same as ncloud=0 with Tiedtke from GFDL-CM3
! ncloud = 12   Same as ncloud=2 with Tiedtke from GFDL-CM3
! ncloud = 13   Same as ncloud=3 with Tiedtke from GFDL-CM3 (i.e., same as ncloud=4)
! ncloud = 20   Same as ncloud=3 with MG cloud
! ncloud = 21   Same as ncloud=3 with 2nd moment condensate
! ncloud = 22   Same as ncloud=3 with MG cloud and 2nd momement condensate
! ncloud = 100  Use Lin et al 2nd moment microphysics
! ncloud = 110  Same as ncloud=100 with Tiedtke from GFDL-CM3
! ncloud = 120  Same as ncloud=100 with MG cloud fraction
    
implicit none
    
private
public convectivecloudfrac, convectivecloudarea
public update_cloud_fraction

! Physical constants
real, parameter :: Dva=2.21    !Diffusivity of qv in air (0 deg. and 1 Pa)
real, parameter :: rKa=2.4e-2  !Thermal conductivity of air (0 deg)

! Tunable parameters for qcloud scheme
real, parameter :: ti = -40.               ! Min T for liquid water clouds in Celsius
real, parameter :: tice=273.15+ti          !Convert ti to Kelvin

! Parameters related to diagnosed convective cloud
real, parameter :: wlc=0.2e-3   !LWC of deep conv cloud (kg/m**3)
real, parameter :: ticon=238.15 !Temp at which conv cloud becomes ice

contains

! This subroutine is the interface for the LDR cloud microphysics
subroutine update_cloud_fraction(cfrac,kbsav,ktsav,land,                            &
                    ps,qccon,qfg,qfrad,qg,qlg,qlrad,t,                              &
                    dpsldt,nettend,stratcloud,clcon,cdrop,em,pblh,idjd,mydiag,      &
                    ncloud,nclddia,ldr,rcrit_l,rcrit_s,rcm,cld_decay,               &
                    vdeposition_mode,tiedtke_form,rkmsave,rkhsave,imax,kl)

use const_phys                    ! Physical constants
use estab                         ! Liquid saturation function
use mgcloud_m , only : mg_progcld
                                  ! MG cloud microphysics
use parm_m, only : nmaxpr, dt
                                  ! Model configuration
use sigs_m                        ! Atmosphere sigma levels

implicit none

integer, intent(in) :: idjd, ncloud, nclddia, ldr, vdeposition_mode
integer, intent(in) :: tiedtke_form
integer, intent(in) :: imax, kl
integer, dimension(imax), intent(in) :: kbsav
integer, dimension(imax), intent(in) :: ktsav
real, dimension(imax,kl), intent(inout) :: cfrac
real, dimension(imax,kl), intent(inout) :: qg, qlg, qfg
real, dimension(imax,kl), intent(inout) :: qlrad, qfrad
real, dimension(imax,kl), intent(inout) :: t
real, dimension(imax,kl), intent(inout) :: nettend
real, dimension(imax,kl), intent(inout) :: stratcloud, clcon, cdrop
real, dimension(imax,kl), intent(out) :: qccon
real, dimension(imax,kl), intent(in) :: dpsldt, rkmsave, rkhsave
real, dimension(imax), intent(in) :: ps
real, dimension(imax), intent(in) :: em, pblh
real, intent(in) :: rcrit_l, rcrit_s, rcm, cld_decay
logical, intent(in) :: mydiag
logical, dimension(imax), intent(in) :: land

integer, dimension(imax) :: kbase,ktop  !Bottom and top of convective cloud
real, dimension(imax,kl) :: prf      !Pressure on full levels (hPa)
real, dimension(imax,kl) :: dprf     !Pressure thickness (hPa)
real, dimension(imax,kl) :: rhoa     !Air density (kg/m3)
real, dimension(imax,kl) :: dz       !Layer thickness (m)
real, dimension(imax,kl) :: ccov     !Cloud cover (may differ from cloud frac if vertically subgrid)
real, dimension(imax,kl) :: qsatg    !Saturation mixing ratio
real, dimension(imax,kl) :: qcl      !Vapour mixing ratio inside convective cloud
real, dimension(imax,kl) :: qenv     !Vapour mixing ratio outside convective cloud
real, dimension(imax,kl) :: tenv     !Temperature outside convective cloud
real, dimension(imax) :: precs       !Amount of stratiform precipitation in timestep (mm)
real, dimension(imax) :: preci       !Amount of stratiform snowfall in timestep (mm)
real, dimension(imax) :: precg       !Amount of stratiform graupel in timestep (mm)
real, dimension(imax) :: wcon        !Convective cloud water content (in-cloud, prescribed)

integer k, iq
real, dimension(imax) :: prf_temp, fl
real, dimension(imax) :: diag_temp
real invdt

! meterological fields
do k = 1,kl
  prf_temp(:) = ps*sig(k)
  prf(:,k)    = 0.01*prf_temp    !ps is SI units
  dprf(:,k)   = -0.01*ps*dsig(k) !dsig is -ve
  rhoa(:,k)   = prf_temp/(rdry*t(:,k))             ! air density
  qsatg(:,k)  = qsat(prf_temp,t(:,k),imax)         ! saturated mixing ratio
  dz(:,k)     = -rdry*dsig(k)*t(:,k)/(grav*sig(k)) ! level thickness in metres
  dz(:,k)     = min( max(dz(:,k), 1.), 2.e4 )
end do

! default values
kbase(:) = 0  ! default
ktop(:)  = 0  ! default

!     Set up convective cloud column
where ( ktsav(:)<kl-1 )
  ktop(:)  = ktsav(:)
  kbase(:) = kbsav(:) + 1
  wcon(:)  = wlc
elsewhere
  wcon(:)  = 0.
end where

#ifdef debug
if ( nmaxpr==1 .and. mydiag ) then
  write(6,*) 'entering leoncld'
  diag_temp(:) = qg(idjd,:)
  write(6,"('qv  ',9f8.3/4x,9f8.3)") diag_temp(:)
  diag_temp(:) = qfg(idjd,:)
  write(6,"('qf  ',9f8.3/4x,9f8.3)") diag_temp(:)
  diag_temp(:) = qlg(idjd,:)
  write(6,"('ql  ',9f8.3/4x,9f8.3)") diag_temp(:)
  diag_temp(:) = qrg(idjd,:)
  write(6,"('qr  ',9f8.3/4x,9f8.3)") diag_temp(:)
  diag_temp(:) = qsng(idjd,:)
  write(6,"('qs  ',9f8.3/4x,9f8.3)") diag_temp(:)
  diag_temp(:) = qgrg(idjd,:)
  write(6,"('qg  ',9f8.3/4x,9f8.3)") diag_temp(:)
endif
#endif


! Calculate convective cloud fraction and adjust moisture variables before calling newcloud
do k = 1,kl
  where ( clcon(:,k)>0. )
    !ccw=wcon(iq)/rhoa(iq,k)  !In-cloud l.w. mixing ratio
    qccon(:,k)  = clcon(:,k)*wcon(:)/rhoa(:,k)
    qenv(:,k)   = max( 1.e-8, (qg(:,k)-clcon(:,k)*max(qsatg(:,k),qg(:,k)))/(1.-clcon(:,k)) )
    qcl(:,k)    = (qg(:,k)-(1.-clcon(:,k))*qenv(:,k))/clcon(:,k)
    qlg(:,k)    = qlg(:,k)/(1.-clcon(:,k))
    qfg(:,k)    = qfg(:,k)/(1.-clcon(:,k))
    stratcloud(:,k) = stratcloud(:,k)/(1.-clcon(:,k))
  elsewhere
    qccon(:,k)  = 0.
    qcl(:,k)    = qg(:,k)
    qenv(:,k)   = qg(:,k)
  end where
  tenv(:,k)   = t(:,k) ! Assume T is the same in and out of convective cloud
end do

#ifdef debug
if ( nmaxpr==1 .and. mydiag ) then
  write(6,*) 'before newcloud'
  diag_temp(:) = t(idjd,:)
  write(6,"('t   ',9f8.2/4x,9f8.2)") diag_temp
  diag_temp(:) = qg(idjd,:)
  write(6,"('qv  ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qfg(idjd,:)
  write(6,"('qf  ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qlg(idjd,:)
  write(6,"('ql  ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qrg(idjd,:)
  write(6,"('qr  ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qsng(idjd,:)
  write(6,"('qs  ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qgrg(idjd,:)
  write(6,"('qg  ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qenv(idjd,:)
  write(6,"('qnv ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qsatg(idjd,:)
  write(6,"('qsat',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qcl(idjd,:)
  write(6,"('qcl ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = clcon(idjd,:)
  write(6,"('clc ',9f8.3/4x,9f8.3)") diag_temp
  write(6,*) 'kbase,ktop ',kbase(idjd),ktop(idjd)
endif
#endif


!     Calculate cloud fraction and cloud water mixing ratios
if ( ncloud==20 .or. ncloud==22 .or. ncloud==120 ) then
  call mg_progcld
else
  call newcloud(dt,land,prf,rhoa,tenv,qenv,qlg,qfg,       &
                dpsldt,nettend,stratcloud,em,pblh,idjd,   &
                mydiag,ncloud,nclddia,rcrit_l,rcrit_s,    &
                cld_decay,vdeposition_mode,tiedtke_form,  &
                rkmsave,rkhsave,imax,kl)
end if

! Vertically sub-grid cloud
ccov(1:imax,1:kl) = stratcloud(1:imax,1:kl)
do k = 2,kl-1
  where ( stratcloud(:,k-1)<1.e-10 .and. stratcloud(:,k)>1.e-2 .and. stratcloud(:,k+1)<1.e-10 )
    ccov(:,k) = sqrt(stratcloud(:,k))
  end where
end do


#ifdef debug
if ( nmaxpr==1 .and. mydiag ) then
  write(6,*) 'after newcloud'
  diag_temp(:) = tenv(idjd,:)
  write (6,"('tnv ',9f8.2/4x,9f8.2)") diag_temp
  diag_temp(:) = qg(idjd,:)
  write (6,"('qv0 ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qfg(idjd,:)
  write (6,"('qf  ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qlg(idjd,:)
  write (6,"('ql  ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qrg(idjd,:)
  write (6,"('qr  ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qsng(idjd,:)
  write (6,"('qs  ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qgrg(idjd,:)
  write (6,"('qg  ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qenv(idjd,:) ! really new qg
  write (6,"('qnv ',9f8.3/4x,9f8.3)") diag_temp
endif
#endif


!     Weight output variables according to non-convective fraction of grid-box
do k = 1,kl
  t(:,k)  = clcon(:,k)*t(:,k) + (1.-clcon(:,k))*tenv(:,k)
  qg(:,k) = clcon(:,k)*qcl(:,k) + (1.-clcon(:,k))*qenv(:,k)
  where ( k>=kbase(:) .and. k<=ktop(:) )
    stratcloud(:,k) = stratcloud(:,k)*(1.-clcon(:,k))
    ccov(:,k) = ccov(:,k)*(1.-clcon(:,k))
    qlg(:,k)  = qlg(:,k)*(1.-clcon(:,k))
    qfg(:,k)  = qfg(:,k)*(1.-clcon(:,k))
  end where
end do

if ( nmaxpr==1 .and. mydiag ) then
  write(6,*) 'before newsnowrain'
  diag_temp(:) = t(idjd,:)
  write (6,"('t   ',9f8.2/4x,9f8.2)") diag_temp
  diag_temp(:) = qg(idjd,:)
  write (6,"('qv  ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qfg(idjd,:)
  write (6,"('qf  ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qlg(idjd,:)
  write (6,"('ql  ',9f8.3/4x,9f8.3)") diag_temp
endif
!if ( diag .and. ntiles==1 ) then
!  call maxmin(t,' t',ktau,1.,kl)
!  call maxmin(qg,'qv',ktau,1.e3,kl)
!  call maxmin(qfg,'qf',ktau,1.e3,kl)
!  call maxmin(qlg,'ql',ktau,1.e3,kl)
!  call maxmin(qrg,'qr',ktau,1.e3,kl)
!  call maxmin(qsng,'qs',ktau,1.e3,kl)
!  call maxmin(qgrg,'qg',ktau,1.e3,kl)
!endif


! Add convective cloud water into fields for radiation
! done because sometimes newrain drops out all qlg, ending up with
! zero cloud (although it will be rediagnosed as 1 next timestep)
do k = 1,kl
  fl(:)      = max(0., min(1., (t(:,k)-ticon)/(273.15-ticon)))
  qlrad(:,k) = qlg(:,k) + fl(:)*qccon(:,k)
  qfrad(:,k) = qfg(:,k) + (1.-fl(:))*qccon(:,k)
  cfrac(:,k) = min( 1., ccov(:,k)+clcon(:,k) ) ! original
end do

return
end subroutine update_cloud_fraction


subroutine newcloud(tdt,land,prf,rhoa,ttg,qtg,qlg,qfg,        &
                    dpsldt,nettend,stratcloud,em,pblh,idjd,   &
                    mydiag,ncloud,nclddia,rcrit_l,rcrit_s,    &
                    cld_decay,vdeposition_mode,tiedtke_form,  &
                    rkmsave,rkhsave,imax,kl)

! This routine is part of the prognostic cloud water scheme

use const_phys                    ! Physical constants
use estab                         ! Liquid saturation function
use parm_m, only : diag, ds       ! Model configuration
use sigs_m                        ! Atmosphere sigma levels

implicit none

! Argument list
integer, intent(in) :: idjd, ncloud, nclddia, vdeposition_mode
integer, intent(in) :: tiedtke_form
integer, intent(in) :: imax, kl
real, dimension(imax,kl), intent(in) :: prf
real, dimension(imax,kl), intent(in) :: rhoa
real, dimension(imax,kl), intent(inout) :: ttg
real, dimension(imax,kl), intent(inout) :: qtg
real, dimension(imax,kl), intent(inout) :: qlg
real, dimension(imax,kl), intent(inout) :: qfg
real, dimension(imax,kl), intent(in) :: dpsldt, rkmsave, rkhsave
real, dimension(imax,kl), intent(inout) :: nettend
real, dimension(imax,kl), intent(inout) :: stratcloud
real, dimension(imax), intent(in) :: em, pblh
real, dimension(kl) :: rfull
real, intent(in) :: tdt
real, intent(in) :: rcrit_l, rcrit_s, cld_decay
logical, intent(in) :: mydiag
logical, dimension(imax), intent(in) :: land

! Local work arrays and variables
real, dimension(imax,kl) :: qsw
real, dimension(imax,kl) :: qcg, qtot, tliq
real, dimension(imax,kl) :: fice, qcold, rcrit
real, dimension(imax) :: tk
real, dimension(imax) :: pk, deles
real, dimension(imax) :: qsi, qsl
real, dimension(kl) :: diag_temp
real decayfac
real es, Aprpr, Bprpr, Cice
real qi0, fd, Crate, Qfdep
real fl, hlrvap, qs, dqsdt
real al, qc, delq, qfnew

integer k, iq

real, parameter :: rhoic = 700.
real, parameter :: cm0 = 1.e-12 !Initial crystal mass

! Start code : ----------------------------------------------------------


#ifdef debug
if ( diag.and.mydiag ) then
  write(6,*) 'entering newcloud'
  diag_temp(:) = prf(idjd,:)
  write(6,'(a,30f10.3)') 'prf ',diag_temp
  diag_temp(:) = ttg(idjd,:)
  write(6,'(a,30f10.3)') 'ttg ',diag_temp
  diag_temp(:) = qtg(idjd,:)
  write(6,*) 'qtg ',diag_temp
  diag_temp(:) = qlg(idjd,:)
  write(6,*) 'qlg ',diag_temp
  diag_temp(:) = qfg(idjd,:)
  write(6,*) 'qfg ',diag_temp
end if
#endif


! First melt cloud ice or freeze cloud water to give correct ice fraction fice.
! Then calculate the cloud conserved variables qtot and tliq.
! Note that qcg is the total cloud water (liquid+frozen)

do k = 1,kl
  do iq = 1,imax
    if ( ttg(iq,k)>=tfrz ) then
      fice(iq,k) = 0.
    else if ( ttg(iq,k)>=tice .and. qfg(iq,k)>1.e-12 ) then
      fice(iq,k) = min(qfg(iq,k)/(qfg(iq,k)+qlg(iq,k)), 1.)
    else if ( ttg(iq,k)>=tice ) then
      fice(iq,k) = 0.
    else
      fice(iq,k) = 1.
    end if
    qcg(iq,k)   = qlg(iq,k) + qfg(iq,k)
    qcold(iq,k) = qcg(iq,k)
    qfnew       = fice(iq,k)*qcg(iq,k)
    ttg(iq,k)   = ttg(iq,k) + hlfcp*(qfnew-qfg(iq,k)) !Release L.H. of fusion
    qfg(iq,k)   = qfnew
    qlg(iq,k)   = max(0., qcg(iq,k)-qfg(iq,k))
    qtot(iq,k)  = qtg(iq,k) + qcg(iq,k)
    tliq(iq,k)  = ttg(iq,k) - hlcp*qcg(iq,k) - hlfcp*qfg(iq,k)
  end do
end do

#ifdef debug
if ( diag .and. mydiag ) then
  write(6,*) 'within newcloud'
  diag_temp = ttg(idjd,:)
  write(6,*) 'ttg ',diag_temp
  diag_temp = qcold(idjd,:)
  write(6,*) 'qcold ',diag_temp
  diag_temp = qcg(idjd,:)
  write(6,*) 'qcg ',diag_temp
  diag_temp = qlg(idjd,:)
  write(6,*) 'qlg ',diag_temp
  diag_temp = qfg(idjd,:)
  write(6,*) 'qfg ',diag_temp
  diag_temp = fice(idjd,:)
  write(6,*) 'fice ',diag_temp
end if
#endif


! Precompute the array of critical relative humidities
if ( nclddia==-3 ) then
  do k = 1,kl
    where ( land(:) )
      rcrit(:,k)=max( rcrit_l, (1.-16.*(1.-sig(k))**3) )
    elsewhere
      rcrit(:,k)=max( rcrit_s, (1.-16.*(1.-sig(k))**3) )
    end where
  enddo
else if ( nclddia<0 ) then
  do k = 1,kl
    where ( land(:) )
      rcrit(:,k)=max( rcrit_l, (1.-4.*(1.-sig(k))**2) )
    elsewhere
      rcrit(:,k)=max( rcrit_s, (1.-4.*(1.-sig(k))**2) )
    end where
  enddo
else if ( nclddia==1 ) then
  do k = 1,kl
    where ( land(:) )
      rcrit(:,k)=max( rcrit_l, sig(k)**3 )
    elsewhere
      rcrit(:,k)=max( rcrit_s, sig(k)**3 )
    end where
  enddo
else if ( nclddia==2 ) then
  do k = 1,kl
    where ( land(:) )
      rcrit(:,k)=rcrit_l
    elsewhere
      rcrit(:,k)=rcrit_s
    end where
  enddo
else if ( nclddia==3 ) then  !   JLM  Feb '21
! rcmid=.85   ! rcrit value at sigma=sigmidd (expect >= rcrit_s)
! rfmid=1.15  ! rh for full cloud at (or below) sigmidd
! rctop=0.9   ! rcrit_top  Lohmann et al (2008) show clear sky for cirrus RH up to 140%
! rftop=1.2   ! rh for full cloud at sig=0
! sigmidd=.4
  do k=1,kl   ! typically set rcrit_l=.85,  rcrit_s=.85
    rfull(k)=1.15    ! full cloud RH below sig=.sigmidd=.4 say
    if (sig(k)<.4     ) then
      rfull(k)=1.2  +sig(k)*(1.15 -1.2  )/.4
    end if
    where ( land(1:imax) )
      rcrit(:,k)=rcrit_l+(.85  -rcrit_l)*(1.-sig(k))/(1.-.4     ) ! gives rcmid at sigmidd
    elsewhere
      rcrit(:,k)=rcrit_s+(.85  -rcrit_l)*(1.-sig(k))/(1.-.4     ) ! gives rcmid at sigmidd
    end where
    if (sig(k)<.4     ) then
      rcrit(:,k)=.9   +sig(k)*(.85  -.9   )/.4
    end if
    where ( sig(k)>(1.-pblh(:)*9.8/(287.*300.)) )
      rcrit(:,k)=.98
    end where
  end do
else if ( nclddia==4 ) then
  do k = 1,kl
    where ( land(:) )
      rcrit(:,k)=max( rcrit_l, sig(k) )             ! .75 for test Mk2/3
    elsewhere
      rcrit(:,k)=max( rcrit_s, sig(k) )             ! .9  for test Mk2/3
    end where
  enddo
else if ( nclddia==5 ) then  ! default till May 08
  do k = 1,kl
    where ( land(:) )
      rcrit(:,k)=max( rcrit_l, min(.99,sig(k)) )    ! .75 for same as T63
    elsewhere
      rcrit(:,k)=max( rcrit_s, min(.99,sig(k)) )    ! .85 for same as T63
    end where
  enddo
else if ( nclddia==6 ) then
  do k = 1,kl
    where ( land(:) )
      rcrit(:,k)=max(rcrit_l*(1.-.15*sig(k)),sig(k)**4)
    elsewhere
      rcrit(:,k)=max(rcrit_s*(1.-.15*sig(k)),sig(k)**4)
    end where
  enddo
else if ( nclddia==7 ) then
  do k = 1,kl
    where ( land(:) )
      rcrit(:,k)=max(rcrit_l*(1.-.2*sig(k)),sig(k)**4)
    elsewhere
      rcrit(:,k)=max(rcrit_s*(1.-.2*sig(k)),sig(k)**4)
    end where
  enddo
else if ( nclddia>7 ) then  ! e.g. 12    JLM
  ! MJT notes - Lopez (2002) "Implementation and validation of a new pronostic large-scale cloud
  ! and precipitation scheme for climate and data-assimilation purposes" Q J R Met Soc 128, 229-257,
  ! has a useful discussion of the dependence of RHcrit on grid spacing
  do k = 1,kl ! typically set rcrit_l=.75,  rcrit_s=.85
    do iq = 1,imax
      al = ds/(em(iq)*208498.)
      fl = (1.+real(nclddia))*al/(1.+real(nclddia)*al)
      ! for rcit_l=.75 & nclddia=12 get rcrit=(0.751, 0.769, .799, .901, .940, .972, .985) for (200, 100, 50, 10, 5, 2, 1) km
      ! for rcit_l=.75 % nclddia=120 get rcrit=(0.750, 0.752, .756, .785, .813, .865, .907) for (200, 100, 50, 10, 5, 2, 1) km
      if ( land(iq) ) then
        rcrit(iq,k) = max( 1.-fl*(1.-rcrit_l), sig(k)**3 )
      else
        rcrit(iq,k) = max( 1.-fl*(1.-rcrit_s), sig(k)**3 )
      end if
    end do
  end do
end if  ! (nclddia<0)  .. else ..

if ( (ncloud/=4 .and. ncloud<10) .or. ncloud==100 ) then
  ! usual diagnostic cloud fraction

  ! Calculate cloudy fraction of grid box (stratcloud) and gridbox-mean cloud water
  ! using the triangular PDF of Smith (1990)

  do k = 1,kl
    ! Calculate qs at temperature tliq
    pk(:) = 100.0*prf(:,k)
    qsi(:) = qsati(pk,tliq(:,k),imax)   !Ice value
    deles(:) = esdiffx(tliq(:,k),imax)  ! MJT suggestion
    qsl(:) = qsi(:) + epsil*deles/pk    !qs over liquid
    qsw(:,k) = fice(:,k)*qsi(:) +    &
               (1.-fice(:,k))*qsl(:)    !Weighted qs at temperature Tliq
    do iq = 1,imax
      ! Calculate gam=(L/cp)*dqsdt at temperature tliq
      hlrvap = (hl+fice(iq,k)*hlf)/rvap
      qs = qsw(iq,k)
      dqsdt = qs*hlrvap/tliq(iq,k)**2
      al = 1./(1.+(hlcp+fice(iq,k)*hlfcp)*dqsdt)  !Smith's notation
      qc = qtot(iq,k) - qs
      delq = (1.-rcrit(iq,k))*qs     !UKMO style (equivalent to above)
      if ( qc<=-delq ) then
        stratcloud(iq,k) = 0.
        qcg(iq,k) = 0.
      else if ( qc<=0. ) then
        stratcloud(iq,k) = max( 1.e-6, 0.5*((qc+delq)/delq)**2 )  ! for roundoff
        qcg(iq,k) = max( 1.e-8, al*(qc+delq)**3/(6.*delq**2) )    ! for roundoff
      else if ( qc<delq ) then
        stratcloud(iq,k) = max( 1.e-6, 1.-0.5*((qc-delq)/delq)**2 ) ! for roundoff
        qcg(iq,k) = max( 1.e-8, al*(qc-(qc-delq)**3/(6.*delq**2)) ) ! for roundoff
      else
        stratcloud(iq,k) = 1.
        qcg(iq,k) = al*qc
      end if
    end do ! iq loop
  end do   ! k loop

  if ( nclddia==3 ) then
    do k = 1,kl
      do iq = 1,imax
        stratcloud(iq,k)=max(0.,min(1.,(qtot(iq,k)/qsw(iq,k)-rcrit(iq,k))/(rfull(k)-rcrit(iq,k))))
      end do
    end do
  end if

#ifdef debug
  if ( diag .and. mydiag ) then
    diag_temp(:) = rcrit(idjd,:)
    write(6,*) 'rcrit ',diag_temp
    diag_temp(:) = qtot(idjd,:)
    write(6,*) 'qtot ',diag_temp
    !diag_temp(:) = qsi(idjd,:)
    !write(6,*) 'qsi',diag_temp
    diag_temp(:) = tliq(idjd,:)
    write(6,*) 'tliq',diag_temp
    !diag_temp(:) = qsl(idjd,:)
    !write(6,*) 'qsl ',diag_temp
    diag_temp(:) = qsw(idjd,:)
    write(6,*) 'qsw ',diag_temp
    diag_temp(:) = stratcloud(idjd,:)
    write(6,*) 'stratcloud',diag_temp
    diag_temp(:) = qtot(idjd,:)-qsw(idjd,:)
    write(6,*) 'qc  ',diag_temp
    diag_temp(:) = qcg(idjd,:)
    write(6,*) 'qcg ',diag_temp
    diag_temp(:) = (1.-rcrit(idjd,:))*qsw(idjd,:)
    write(6,*) 'delq ',diag_temp
  endif
#endif


  ! Assume condensation or evaporation retains ice fraction fice.
  ! Introduce a time-decay factor for cirrus (as suggested by results of Khvorostyanov & Sassen,
  ! JAS, 55, 1822-1845, 1998). Their suggested range for the time constant is 0.5 to 2 hours.
  ! The grid-box-mean values of qtg and ttg are adjusted later on (below).
  decayfac = exp ( -tdt/cld_decay )  ! Try cld_decay=2 hrs (decayfac=0. is instant adjustment for the old scheme)
  do k = 1,kl
    where( ttg(:,k)>=Tice )
      qfg(:,k) = fice(:,k)*qcg(:,k)
      qlg(:,k) = qcg(:,k) - qfg(:,k)
    elsewhere                                 ! Cirrus T range
      qfg(:,k) = qcold(:,k)*decayfac + qcg(:,k)*(1.-decayfac)
      qlg(:,k) = 0.
      qcg(:,k) = qfg(:,k)
    end where
  end do

else

  ! Tiedtke prognostic cloud fraction model
  ! MJT notes - we use ttg instead of tliq
  do k = 1,kl
    pk = 100.*prf(:,k)
    qsi(:) = qsati(pk,ttg(:,k),imax) ! Ice value
    deles = esdiffx(ttg(:,k),imax)
    qsl(:) = qsi(:) + epsil*deles/pk ! Liquid value
    qsw(:,k) = fice(:,k)*qsi(:) + (1.-fice(:,k))*qsl(:)        ! Weighted qs at temperature Tliq
    qcg(:,k) = qlg(:,k) + qfg(:,k)
  end do

  call progcloud(tdt,qcg,qtot,prf,rhoa,fice,qsw,ttg,rcrit,  &
                 dpsldt,nettend,stratcloud,tiedtke_form,    &
                 rkmsave,rkhsave,imax,kl)

  decayfac = exp ( -tdt/cld_decay )  ! Try 2 hrs
  !decayfac = 0.                     ! Instant adjustment (old scheme)
  do k = 1,kl
    where( ttg(:,k)>=Tice )
      qfg(:,k) = fice(:,k)*qcg(:,k)
      qlg(:,k) = qcg(:,k) - qfg(:,k)
    elsewhere                                 ! Cirrus T range
      qfg(:,k) = qcold(:,k)*decayfac + qcg(:,k)*(1.-decayfac)
      qlg(:,k) = 0.
      qcg(:,k) = qfg(:,k)
    end where
  end do

end if ! ncloud/=4 .and. ncloud<10 ..else..

! Do the vapour deposition calculation in mixed-phase clouds:
! Calculate deposition on cloud ice, assuming es(T) is the weighted value of the
! liquid and ice values.
if ( vdeposition_mode==0 ) then
  do k = 1,kl
    Tk(:) = tliq(:,k) + hlcp*(qlg(:,k)+qfg(:,k))/max(stratcloud(:,k),1.e-10) !T in liq cloud
    pk(:) = 100.*prf(:,k)
    qsi(:) = qsati(pk,Tk,imax)
    deles(:) = (1.-fice(:,k))*esdiffx(Tk,imax)
    do iq = 1,imax
      if ( stratcloud(iq,k)>0. .and. Tk(iq)<tfrz .and. qlg(iq,k)>1.e-8 ) then
        qs        = qsi(iq)
        es        = qs*pk(iq)/0.622 !ice value
        Aprpr     = hl/(rKa*Tk(iq))*(hls/(rvap*Tk(iq))-1.)
        Bprpr     = rvap*Tk(iq)/((Dva/pk(iq))*es)
        Cice      = 1.e3*exp(12.96*deles(iq)/es - 0.639) !Meyers et al 1992
        qi0       = cm0*Cice/rhoa(iq,k) !Initial ice mixing ratio
        ! Next 2 lines are for assumption of fully mixed ql and qf (also a line further down).
        qi0       = max(qi0, qfg(iq,k)/stratcloud(iq,k)) !Assume all qf and ql are mixed
        fd        = 1.       !Fraction of cloud in which deposition occurs
        Crate     = 7.8*((Cice/rhoa(iq,k))**2/rhoic)**(1./3.)*deles(iq)/((Aprpr+Bprpr)*es)
        qfdep     = fd*stratcloud(iq,k)*sqrt(((2./3.)*Crate*tdt+qi0**(2./3.))**3)
        ! Also need this line for fully-mixed option...
        qfdep     = qfdep - qfg(iq,k)
        qfdep     = min(qfdep, qlg(iq,k))
        qlg(iq,k) = qlg(iq,k) - qfdep
        qfg(iq,k) = qfg(iq,k) + qfdep
        fice(iq,k) = qfg(iq,k)/max(qfg(iq,k)+qlg(iq,k),1.e-30)
      end if
    end do
  end do
else
! use ql/(qf+ql) for deposition
  do k = 1,kl
    Tk(:) = tliq(:,k) + hlcp*(qlg(:,k)+qfg(:,k))/max(stratcloud(:,k),1.e-10) !T in liq cloud
    pk(:) = 100.*prf(:,k)
    qsi(:) = qsati(pk,Tk,imax)
    deles(:) = (1.-fice(:,k))*esdiffx(Tk,imax)
    do iq = 1,imax
      if ( stratcloud(iq,k)>0. .and. Tk(iq)<tfrz .and. qlg(iq,k)>1.e-8 ) then
        fl        = qlg(iq,k)/max(qfg(iq,k)+qlg(iq,k),1.e-30)
        qs        = qsi(iq)
        es        = qs*pk(iq)/0.622 !ice value
        Aprpr     = hl/(rKa*Tk(iq))*(hls/(rvap*Tk(iq))-1.)
        Bprpr     = rvap*Tk(iq)/((Dva/pk(iq))*es)
        Cice      = 1.e3*exp(12.96*deles(iq)/es - 0.639) !Meyers et al 1992
        qi0       = cm0*Cice/rhoa(iq,k) !Initial ice mixing ratio
        qi0       = max(qi0, qfg(iq,k)/stratcloud(iq,k)) !Assume all qf and ql are mixed
        fd        = fl      !Or, use option of adjacent ql,qf
        Crate     = 7.8*((Cice/rhoa(iq,k))**2/rhoic)**(1./3.)*deles(iq)/((Aprpr+Bprpr)*es)
        qfdep     = fd*stratcloud(iq,k)*sqrt(((2./3.)*Crate*tdt+qi0**(2./3.))**3)
        ! Also need this line for fully-mixed option...
        qfdep     = qfdep - qfg(iq,k)
        qfdep     = min(qfdep, qlg(iq,k))
        qlg(iq,k) = qlg(iq,k) - qfdep
        qfg(iq,k) = qfg(iq,k) + qfdep
        fice(iq,k) = qfg(iq,k)/max(qfg(iq,k)+qlg(iq,k),1.e-30)
      end if
    end do
  end do

end if

! Calculate new values of vapour mixing ratio and temperature
do k = 1,kl
  qtg(:,k) = qtot(:,k) - qcg(:,k)
  ttg(:,k) = tliq(:,k) + hlcp*qcg(:,k) + hlfcp*qfg(:,k)
end do

#ifdef debug
if ( diag .and. mydiag ) then
   write(6,*) 'at end of newcloud'
   diag_temp(:) = ttg(idjd,:)
   write(6,*) 'ttg ',diag_temp
   diag_temp(:) = qcg(idjd,:)
   write(6,*) 'qcg ',diag_temp
   diag_temp(:) = qlg(idjd,:)
   write(6,*) 'qlg ',diag_temp
   diag_temp(:) = qfg(idjd,:)
   write(6,*) 'qfg ',diag_temp
   diag_temp(:) = qtg(idjd,:)
   write(6,*) 'qtg ',diag_temp
end if
#endif

return
end subroutine newcloud

subroutine progcloud(dt,qc,qtot,press,rho,fice,qs,t,rcrit, &
                     dpsldt,nettend,stratcloud,tiedtke_form, &
                     rkmsave,rkhsave,imax,kl)

use const_phys                    ! Physical constants
use parm_m, only : qgmin          ! Model configuration

implicit none

integer, intent(in) :: tiedtke_form
integer, intent(in) :: imax, kl
integer k
real, dimension(imax,kl), intent(inout) :: qc ! condensate = qf + ql
real, dimension(imax,kl), intent(in) :: qtot, rho, fice, qs, t, rcrit, press
real, dimension(imax,kl), intent(in) :: dpsldt, rkmsave, rkhsave
real, dimension(imax,kl), intent(inout) :: nettend
real, dimension(imax,kl), intent(inout) :: stratcloud
real, dimension(imax,kl) :: xf
real, dimension(imax) :: aa, bb, cc, a_dt, b_dt, cf1, cfeq, cfbar
real, dimension(imax) :: qv, omega, hlrvap, dqsdT, gamma, dqs
real, dimension(imax) :: delq, da
real, intent(in) :: dt
real, dimension(imax,kl) :: erosion_scale
real, parameter :: erosion_scale_d = 1.e-6
real, parameter :: erosion_scale_t = 5.e-5
real, parameter :: diff_threshold_t = 0.1
real, parameter :: u00 = 0.8

! background erosion scale in 1/secs
erosion_scale(:,:) = erosion_scale_d

!Convection is treated independently for now
!if ( ncloud==15 ) then
!  ! convert convective mass flux from half levels to full levels
!  do k = 1,kl-1
!    cmflx(:,k) = rathb(k)*fluxtot(:,k)+ratha(k)*fluxtot(:,k+1)
!  end do
!  cmflx(:,kl) = rathb(kl)*fluxtot(:,kl)
!else ! ncloud==4 .or. ncloud==10 .or. ncloud==12 .or. ncloud==13
!  ! use convective area fraction in leoncld.f, instead of convective mass flux
!  cmflx = 0.
!end if

! Turbulence
do k = 1,kl
  where ( 0.5*(rkmsave(:,k)+rkhsave(:,k)) > diff_threshold_t )
    erosion_scale(:,k) = erosion_scale_t
  end where
end do

! Aerosols
!do k = 1,kl
!  edum = drop(:,k)*0.4/150. + 0.8
!  erosion_scale(:,k) = edum*erosion_scale(:,k)
!end do


! calculate dqs = ((omega + grav*Mc)/(cp*rho)+nettend)*dqsdT*dt
!                 -------------------------------------------------------
!                 1 + (stratcloud + 0.5*da)*gamma
! MJT notes - GFDL AM adds (stratcloud+0.5*at*da)*gamma term

! Change in saturated volume fraction
! da = -0.5*(1.-cf)^2*dqs/(qs-qv)
! MJT notes - Tiedtke 93 does not use 0.5

! gamma = L/cp*dqsdT

! Follow GFDL AM approach since da=da(dqs), hence need to solve the above
! quadratic equation for dqs if da/=0

! dqs*dqs*AA + dqs*BB + CC = 0
! AA = 0.25*gamma*(1-cf)^2/(qs-qv)
! BB = -(1+gamma*cf)
! CC = ((omega + grav*mflx)/(cp*rho)+netten)*dqsdT*dt

select case(tiedtke_form)
  case(1) ! Original GFDL-AM4
    do k = 1,kl
      hlrvap = (hl+fice(:,k)*hlf)/rvap
      dqsdT = qs(:,k)*hlrvap/(t(:,k)**2)
      gamma = (hlcp+fice(:,k)*hlfcp)*dqsdT
      omega = press(:,k)*dpsldt(:,k)
      ! 1st estimate with da=0.
      dqs = (omega/(cp*rho(:,k))+nettend(:,k)*dqsdt*dt) &
           /(1.+stratcloud(:,k)*gamma)
      qv = qtot(:,k) - qc(:,k)
      where ( qv>u00*qs(:,k) .and. dqs<0. .and. stratcloud(:,k)<1. )
        xf(:,k) = 1.
      elsewhere
        xf(:,k) = 0.
      end where
    end do
  case default
    write(6,*) "ERROR: Unsupported tiedtke_form = ",tiedtke_form
    stop
end select

do k = 1,kl
  stratcloud(:,k) = max( min( stratcloud(:,k), 1. ), 0. )

  qv = qtot(:,k) - qc(:,k)
  ! calculate vertical velocity, dqs/dT and gamma
  omega = press(:,k)*dpsldt(:,k)
  hlrvap = (hl+fice(:,k)*hlf)/rvap
  dqsdT = qs(:,k)*hlrvap/(t(:,k)**2)
  gamma = (hlcp+fice(:,k)*hlfcp)*dqsdT

  !cc = ((omega + grav*cmflx(:,k))/(cp*rho(:,k))+nettend(:,k))*dt*dqsdT
  cc = (omega/(cp*rho(:,k))+nettend(:,k))*dt*dqsdT ! neglect cmflx
  aa = 0.5*(1.-stratcloud(:,k))**2/max( qs(:,k)-qv, 1.e-10 )
  bb = 1.+gamma*stratcloud(:,k)
  !dqs = ( bb - sqrt( bb*bb - 2.*gamma*xf*aa*cc ) ) / ( gamma*xf*aa ) ! GFDL style
  !dqs = min( dqs, cc/(1. + 0.5*bb) )                                 ! GFDL style
  dqs = 2.*cc/( bb + sqrt( bb**2 - 2.*gamma*xf(:,k)*aa*cc ) ) ! alternative form of quadratic equation
                                                              ! note that aa and bb have been multipled by 2 and -1, respectively.
  da = -aa*dqs

  ! Large scale cloud formation via condensation (A)
  ! Large scale cloud destruction via erosion (B)
  where ( da>0. )
    !a_dt = xf(:,k)*da/max(1.-stratcloud(:,k),1.e-10)
    a_dt = -xf(:,k)*dqs*0.5*(1.-stratcloud(:,k))/max( qs(:,k)-qv, 1.e-10 )
    b_dt = 0.
  elsewhere ( qc(:,k)>1.e-10 )
    a_dt = 0.
    b_dt = stratcloud(:,k)*erosion_scale(:,k)*dt*max(qs(:,k)-qv,0.)/qc(:,k)
  elsewhere
    a_dt = 0.
    b_dt = 0.
  end where

  ! Integrate
  !   dcf/dt = (1-cf)*A - cf*B
  ! to give (use cf' = A-cf*(A+B))
  !   cf(t=1) = cfeq + (cf(t=0) - cfeq)*exp(-(A+B)*dt)
  !   cfeq = A/(A+B)
  ! Average cloud fraction over the interval t=tau to t=tau+1
  !   cfbar = cfeq - (cf(t=1) - cf(t=0))/((A+B)*dt)
  ! cfeq is the equilibrum cloud fraction that is approached with
  ! a time scale of 1/(A+B)
  where ( a_dt>1.e-7 .or. b_dt>1.e-7 )
    cfeq  = a_dt/(a_dt+b_dt)
    cf1   = cfeq + (stratcloud(:,k) - cfeq)*exp(-a_dt-b_dt)
    cfbar = cfeq + (stratcloud(:,k) - cf1 )/(a_dt+b_dt)
  elsewhere
    cfeq  = stratcloud(:,k)
    cf1   = stratcloud(:,k)
    cfbar = stratcloud(:,k)
  end where

  ! Change in condensate
  ! dqc = -dqs*(stratcloud+0.5*da)
  qc(:,k) = max(min( qc(:,k) - dqs*(cf1+0.5*da), qtot(:,k)-qgmin), 0.)
  !qc(:,k) = max(min( qc(:,k) - cfbar*dqs, qtot(:,k)-qgmin ), 0. )

  stratcloud(:,k) = min( cf1, 1. )

  ! Change in cloud fraction
  where ( qc(:,k)>1.e-10 )
    stratcloud(:,k) = max( stratcloud(:,k), 1.e-10 )
  elsewhere
    ! MJT notes - cloud fraction is maintained (da=0.) while condesate evaporates (dqc<0.) until
    ! the condesate dissipates
    stratcloud(:,k) = 0.
    qc(:,k) = 0.
  end where

  ! Reset tendency and mass flux for next time-step
  nettend(:,k) = 0.

end do

return
end subroutine progcloud
    
subroutine convectivecloudfrac(clcon,kbsav,ktsav,condc,cldcon)

use parm_m           ! Model configuration

implicit none

include 'kuocom.h'   ! Convection parameters

integer k, kl
real, dimension(:,:), intent(out) :: clcon
real, dimension(:), intent(out), optional :: cldcon
real, dimension(size(clcon,1)) :: cldcon_temp
real, dimension(size(clcon,1)) :: n, cldcon_local
integer, dimension(:), intent(in) :: kbsav
integer, dimension(:), intent(in) :: ktsav
real, dimension(:), intent(in) :: condc

! MJT notes - This is an old parameterisation from NCAR.  acon and
! bcon represent shallow and deep convection, respectively.  It can
! be argued that acon should be zero in CCAM to avoid discontinuous
! evolution of the model.  Furthermore, a much better fit can be
! obtained from the mass flux, rather than rainfall.  It also
! should be noted that acon and bcon are likely to depend on the
! spatial resolution.

kl = size(clcon,2)

cldcon_temp = 0. ! for cray compiler
call convectivecloudarea(cldcon_temp,ktsav,condc)
if ( present(cldcon) ) then
  cldcon = cldcon_temp
end if

! Impose cloud overlap assumption
if ( nmr>=1 ) then
  cldcon_local = cldcon_temp               ! maximum overlap
else
  n = 1./real(max(ktsav-kbsav,1))
  cldcon_local = 1. - (1.-cldcon_temp)**n  ! random overlap
end if

do k = 1,kl
  where( k<kbsav+1 .or. k>ktsav )
    clcon(:,k) = 0.
  elsewhere
    clcon(:,k) = cldcon_local
  end where
end do

return
end subroutine convectivecloudfrac

pure subroutine convectivecloudarea(cldcon,ktsav,condc)

use newmpar_m        ! Grid parameters
use parm_m           ! Model configuration

implicit none

include 'kuocom.h'   ! Convection parameters

integer, dimension(:), intent(in) :: ktsav
real, dimension(:), intent(in) :: condc
real, dimension(:), intent(out) :: cldcon

where ( ktsav<kl-1 )
  cldcon = min( acon+bcon*log(1.+condc*86400./dt), 0.8 ) !NCAR
elsewhere
  cldcon = 0.
end where

return
end subroutine convectivecloudarea

end module cloudmod
