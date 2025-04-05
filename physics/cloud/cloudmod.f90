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
    
module cloudmod

! This module is the Rotstayn 1997 cloud microphysics parameterisation for cloud fraction.
! prognostic cloud fraction scheme based on Tiedtke from GFDL-CM3.

implicit none
    
private
public convectivecloudfrac
public update_cloud_fraction
public cloud_ice_method

! Physical constants
real, parameter :: Dva = 2.21    !Diffusivity of qv in air (0 deg. and 1 Pa)
real, parameter :: rKa = 2.4e-2  !Thermal conductivity of air (0 deg)

! Tunable parameters for qcloud scheme
real, parameter :: ti = -40.               ! Min T for liquid water clouds in Celsius
real, parameter :: tice = 273.15+ti        ! Convert ti to Kelvin

! Parameters related to diagnosed convective cloud
real, parameter :: wlc = 0.2e-3   !LWC of deep conv cloud (kg/m**3)
real, parameter :: ticon = 238.15 !Temp at which conv cloud becomes ice

integer, save :: cloud_ice_method = 0   ! method for 0C to -40C (0=liq, 1=linear)

contains

! This subroutine is the interface for the LDR cloud microphysics
subroutine update_cloud_fraction(cfrac,land,                             &
                    ps,qccon,qfg,qfrad,qg,qlg,qlrad,t,                   &
                    dpsldt,rad_tend,trb_tend,trb_qend,stratcloud,clcon,  &
                    em,pblh,idjd,mydiag,nclddia,rcrit_l,                 &
                    rcrit_s,rcm,cld_decay,vdeposition_mode,              &
                    tiedtke_form,rkmsave,rkhsave,cmode)

use const_phys                    ! Physical constants
use estab                         ! Liquid saturation function
use parm_m, only : nmaxpr, dt     ! Model configuration
use sigs_m                        ! Atmosphere sigma levels

implicit none

real, dimension(:,:), intent(inout) :: t
integer, intent(in) :: idjd, nclddia, vdeposition_mode
integer, intent(in) :: tiedtke_form
real, dimension(size(t,1),size(t,2)), intent(inout) :: qg, qlg, qfg
real, dimension(size(t,1),size(t,2)), intent(out) :: qlrad, qfrad
real, dimension(size(t,1),size(t,2)), intent(inout) :: rad_tend, trb_tend, trb_qend
real, dimension(size(t,1),size(t,2)), intent(inout) :: stratcloud
real, dimension(size(t,1),size(t,2)), intent(out) :: qccon, cfrac
real, dimension(size(t,1),size(t,2)), intent(in) :: dpsldt, rkmsave, rkhsave
real, dimension(size(t,1),size(t,2)), intent(in) :: clcon
real, dimension(size(t,1)), intent(in) :: ps
real, dimension(size(t,1)), intent(in) :: em, pblh
real, intent(in) :: rcrit_l, rcrit_s, rcm, cld_decay
logical, intent(in) :: mydiag
logical, dimension(size(t,1)), intent(in) :: land
character(len=*), intent(in) :: cmode

!integer, dimension(size(t,1)) :: kbase,ktop  !Bottom and top of convective cloud
real, dimension(size(t,1),size(t,2)) :: prf      !Pressure on full levels (hPa)
real, dimension(size(t,1),size(t,2)) :: rhoa     !Air density (kg/m3)
real, dimension(size(t,1),size(t,2)) :: ccov     !Cloud cover (may differ from cloud frac if vertically subgrid)
real, dimension(size(t,1),size(t,2)) :: qcl      !Vapour mixing ratio inside convective cloud
real, dimension(size(t,1),size(t,2)) :: qenv     !Vapour mixing ratio outside convective cloud
real, dimension(size(t,1),size(t,2)) :: tenv     !Temperature outside convective cloud
real, dimension(size(t,1),size(t,2)) :: fice
real, dimension(size(t,1),size(t,2)) :: qtot
real, dimension(size(t,1),size(t,2)) :: qcold
real, dimension(size(t,1),size(t,2)) :: qcg
real, dimension(size(t,1),size(t,2)) :: tliq

integer imax, kl, k, iq
real fl, prf_temp
real qsatg !Saturation mixing ratio
real wcon  !Convective cloud water content (in-cloud, prescribed)
real, dimension(size(t,2)) :: diag_temp


imax = size(t,1)
kl = size(t,2)

! meterological fields
do k = 1,kl
  do iq = 1,imax
    prf_temp = ps(iq)*sig(k)
    !dprf(iq,k) = -0.01*ps(iq)*dsig(k)  !dsig is -ve
    prf(iq,k) = 0.01*prf_temp           !ps is SI units
    rhoa(iq,k) = prf_temp/(rdry*t(iq,k))     ! air density
  end do
end do

#ifdef debug
if ( nmaxpr==1 .and. mydiag ) then
  write(6,*) 'entering leoncld'
  diag_temp(:) = qg(idjd,:)
  write(6,"('qv  ',9f8.3/4x,9f8.3)") diag_temp(:)
  diag_temp(:) = qfg(idjd,:)
  write(6,"('qf  ',9f8.3/4x,9f8.3)") diag_temp(:)
  diag_temp(:) = qlg(idjd,:)
  write(6,"('ql  ',9f8.3/4x,9f8.3)") diag_temp(:)
endif
#endif


! Calculate convective cloud fraction and adjust moisture variables before calling newcloud
tenv(:,:) = t(:,:) ! Assume T is the same in and out of convective cloud
do k = 1,kl
  do iq = 1,imax
    if ( clcon(iq,k)>1.e-10 ) then
      wcon = wlc
      !ccw=wcon(iq)/rhoa(iq,k)  !In-cloud l.w. mixing ratio
      prf_temp = ps(iq)*sig(k)
      qsatg = qsat(prf_temp,t(iq,k)) ! saturated mixing ratio
      qccon(iq,k)  = clcon(iq,k)*wcon/rhoa(iq,k)
      qcl(iq,k)    = max(qsatg,qg(iq,k))
      qenv(iq,k)   = max( 0., (qg(iq,k)-clcon(iq,k)*qcl(iq,k))/(1.-clcon(iq,k)) )
      !qcl(iq,k)   = (qg(iq,k)-(1.-clcon(iq,k))*qenv(iq,k))/clcon(iq,k)
      qlg(iq,k)    = qlg(iq,k)/(1.-clcon(iq,k))
      qfg(iq,k)    = qfg(iq,k)/(1.-clcon(iq,k))
      stratcloud(iq,k) = stratcloud(iq,k)/(1.-clcon(iq,k))
    else
      qccon(iq,k)  = 0.
      qcl(iq,k)    = qg(iq,k)
      qenv(iq,k)   = qg(iq,k)
    end if
  end do
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
  diag_temp(:) = qenv(idjd,:)
  write(6,"('qnv ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = qcl(idjd,:)
  write(6,"('qcl ',9f8.3/4x,9f8.3)") diag_temp
  diag_temp(:) = clcon(idjd,:)
  write(6,"('clc ',9f8.3/4x,9f8.3)") diag_temp
endif
#endif


do k = 1,kl
  do iq = 1,imax
    qtot(iq,k)  = qenv(iq,k) + qlg(iq,k) + qfg(iq,k)
    qcold(iq,k) = qlg(iq,k) + qfg(iq,k)
  end do
end do

! Update ice fraction
fice = calc_fice(tenv,qlg,qfg)

do k = 1,kl
  do iq = 1,imax  
    qcg(iq,k)  = qlg(iq,k) + qfg(iq,k)
    tliq(iq,k) = tenv(iq,k) - (hlcp+hlfcp*fice(iq,k))*qcg(iq,k)
  end do  
end do

!     Calculate cloud fraction and cloud water mixing ratios
call newcloud(dt,land,ps,prf,rhoa,tliq,qtot,qcg,fice,   &
              dpsldt,rad_tend,trb_tend,trb_qend,        &
              stratcloud,em,pblh,idjd,mydiag,           &
              nclddia,rcrit_l,rcrit_s,tiedtke_form,     &
              rkmsave,rkhsave,imax,kl,cmode)

! Update condensate
call saturation_adjustment(dt,cld_decay,vdeposition_mode, &
                           tliq,qtot,qcg,qcold,fice,      &
                           stratcloud,prf,rhoa,           &
                           tenv,qenv,qlg,qfg)


#ifdef debug
if ( nmaxpr==1 .and. mydiag ) then
   write(6,*) 'after saturation_adjustment'
   diag_temp(:) = t(idjd,:)
   write(6,*) 'ttg ',diag_temp
   diag_temp(:) = qcg(idjd,:)
   write(6,*) 'qcg ',diag_temp
   diag_temp(:) = qlg(idjd,:)
   write(6,*) 'qlg ',diag_temp
   diag_temp(:) = qfg(idjd,:)
   write(6,*) 'qfg ',diag_temp
   diag_temp(:) = qg(idjd,:)
   write(6,*) 'qtg ',diag_temp
end if
#endif


! Vertically sub-grid cloud
ccov(1:imax,1:kl) = stratcloud(1:imax,1:kl)

do k = 2,kl-1
  do iq = 1,imax
    if ( stratcloud(iq,k-1)<1.e-10 .and. stratcloud(iq,k)>1.e-10 .and. stratcloud(iq,k+1)<1.e-10 ) then
      ccov(iq,k) = sqrt(stratcloud(iq,k))
    end if
  end do
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
  diag_temp(:) = qenv(idjd,:) ! really new qg
  write (6,"('qnv ',9f8.3/4x,9f8.3)") diag_temp
endif
#endif


!     Weight output variables according to non-convective fraction of grid-box
do k = 1,kl
  do iq = 1,imax
    t(iq,k)  = clcon(iq,k)*t(iq,k) + (1.-clcon(iq,k))*tenv(iq,k)
    qg(iq,k) = clcon(iq,k)*qcl(iq,k) + (1.-clcon(iq,k))*qenv(iq,k)
    stratcloud(iq,k) = stratcloud(iq,k)*(1.-clcon(iq,k))
    ccov(iq,k) = ccov(iq,k)*(1.-clcon(iq,k))
    qlg(iq,k)  = qlg(iq,k)*(1.-clcon(iq,k))
    qfg(iq,k)  = qfg(iq,k)*(1.-clcon(iq,k))
  end do
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
!endif


! Add convective cloud water into fields for radiation
! done because sometimes newrain drops out all qlg, ending up with
! zero cloud (although it will be rediagnosed as 1 next timestep)
do k = 1,kl
  do iq = 1,imax
    fl = max(0., min(1., (t(iq,k)-ticon)/(273.15-ticon)))
    qlrad(iq,k) = qlg(iq,k) + fl*qccon(iq,k)
    qfrad(iq,k) = qfg(iq,k) + (1.-fl)*qccon(iq,k)
    cfrac(iq,k) = max( 0., min( 1., ccov(iq,k)+clcon(iq,k) ) ) ! original
  end do
end do

return
end subroutine update_cloud_fraction

function calc_fice(ttg,qlg,qfg) result(fice)

use const_phys                    ! Physical constants

implicit none

integer k, kl
real, dimension(:,:), intent(in) :: ttg
real, dimension(:,:), intent(inout) :: qlg, qfg
real, dimension(size(ttg,1),size(ttg,2)) :: fice
real, dimension(size(ttg,1)) :: qcg

kl = size(ttg,2)

! First melt cloud ice or freeze cloud water to give correct ice fraction fice.
! Then calculate the cloud conserved variables qtot and tliq.
! Note that qcg is the total cloud water (liquid+frozen)

select case(cloud_ice_method)
  case(0) ! original
    do k = 1,kl
      where ( ttg(:,k)>=tfrz )
        fice(:,k) = 0.
      else where ( ttg(:,k)>=tice .and. qfg(:,k)>0. )
        fice(:,k) = qfg(:,k)/(qfg(:,k)+qlg(:,k))
      else where ( ttg(:,k)>=tice )
        fice(:,k) = 0.
      elsewhere
        fice(:,k) = 1.
      end where
    end do
  !case(1) ! to be depreciated
  !  do k = 1,kl
  !    where ( ttg(:,k)>=tfrz )
  !      fice(:,k) = 0.
  !    else where ( ttg(:,k)>=tice .and. qlg(:,k)+qfg(:,k)>0. )
  !      fice(:,k) = qfg(:,k)/(qfg(:,k)+qlg(:,k))
  !    else where ( ttg(:,k)>=tice )
  !      fice(:,k) = 1. - (ttg(:,k)-tice)/(tfrz-tice) ! MJT suggestion
  !    elsewhere
  !      fice(:,k) = 1.
  !    end where
  !  end do
  case default
    write(6,*) "ERROR: Invalid cloud_ice_method ",cloud_ice_method
    stop
end select
  
return
end function calc_fice
  
subroutine newcloud(tdt,land,ps,prf,rhoa,tliq,qtot,qcg,fice,  &
                    dpsldt,rad_tend,trb_tend,trb_qend,        &
                    stratcloud,em,pblh,idjd,mydiag,           &
                    nclddia,rcrit_l,rcrit_s,tiedtke_form,     &
                    rkmsave,rkhsave,imax,kl,cmode)

! This routine is part of the prognostic cloud water scheme

use const_phys                    ! Physical constants
use estab                         ! Liquid saturation function
use parm_m, only : diag, ds       ! Model configuration
use sigs_m                        ! Atmosphere sigma levels

implicit none

! Argument list
integer, intent(in) :: idjd, nclddia
integer, intent(in) :: tiedtke_form
integer, intent(in) :: imax, kl
real, dimension(imax,kl), intent(in) :: prf
real, dimension(imax,kl), intent(in) :: rhoa
real, dimension(imax,kl), intent(in) :: tliq
real, dimension(imax,kl), intent(in) :: qtot
real, dimension(imax,kl), intent(in) :: fice
real, dimension(imax,kl), intent(inout) :: qcg
real, dimension(imax,kl), intent(in) :: dpsldt, rkmsave, rkhsave
real, dimension(imax,kl), intent(inout) :: rad_tend, trb_tend, trb_qend
real, dimension(imax,kl), intent(inout) :: stratcloud
real, dimension(imax), intent(in) :: em, pblh, ps
real, dimension(kl) :: rfull
real, intent(in) :: tdt
real, intent(in) :: rcrit_l, rcrit_s
logical, intent(in) :: mydiag
logical, dimension(imax), intent(in) :: land
character(len=*), intent(in) :: cmode

! Local work arrays and variables
real, dimension(imax,kl) :: qsw
real, dimension(imax,kl) :: rcrit
real, dimension(imax,kl) :: pk, deles
real, dimension(imax,kl) :: qsi, qsl
real, dimension(kl) :: diag_temp
real steepness
real es, Aprpr, Bprpr, Cice
real qi0, fd, Crate, Qfdep
real fl, hlrvap, qs, dqsdt
real al, qc, delq, qfnew
real tk
real decayfac

integer k, iq

real, parameter :: rhoic = 700.
real, parameter :: cm0 = 1.e-12 !Initial crystal mass

! Start code : ----------------------------------------------------------

#ifdef debug
if ( diag.and.mydiag ) then
  write(6,*) 'entering newcloud'
  diag_temp(:) = prf(idjd,:)
  write(6,'(a,30f10.3)') 'prf ',diag_temp
  write(6,*) 'within newcloud'
  diag_temp = qcg(idjd,:)
  write(6,*) 'qcg ',diag_temp
  diag_temp = fice(idjd,:)
  write(6,*) 'fice ',diag_temp
end if
#endif


! -------------------------------------------------------------------
! Update cloud fraction


! Precompute the array of critical relative humidities
if ( nclddia==-4 ) then
  do k = 1,kl ! typically set rcrit_l=.805,  rcrit_s=.805
    do iq = 1,imax
      al = ds/(em(iq)*208498.)
      fl = (1.+0.526)*al/(1.+0.526*al) ! nclddia=0.526
      ! for rcit_l=.805 & nclddia=0.526 get rcrit=(0.805, 0.882, .934, .986, .993, .997, .999) for (200, 100, 50, 10, 5, 2, 1) km
      if ( land(iq) ) then
        rcrit(iq,k) = 1.-fl*(1.-rcrit_l)
      else
        rcrit(iq,k) = 1.-fl*(1.-rcrit_s)  
      end if    
      if ( sig(k) >= 0.8 ) then
        steepness = 65.
        rcrit(iq,k) = rcrit(iq,k) + (1. - rcrit(iq,k)) * 0.5 * (1.0 + tanh(steepness * (sig(k) - 0.9)))
      end if
    end do
  end do
else if ( nclddia==-3 ) then
  do k = 1,kl
    where ( land(:) )
      rcrit(:,k)=max( rcrit_l, (1.-16.*(1.-sig(k))**3) )
    elsewhere
      rcrit(:,k)=max( rcrit_s, (1.-16.*(1.-sig(k))**3) )
    end where
  end do
else if ( nclddia<0 ) then
  do k = 1,kl
    where ( land(:) )
      rcrit(:,k)=max( rcrit_l, (1.-4.*(1.-sig(k))**2) )
    elsewhere
      rcrit(:,k)=max( rcrit_s, (1.-4.*(1.-sig(k))**2) )
    end where
  end do
else if ( nclddia==1 ) then
  do k = 1,kl
    where ( land(:) )
      rcrit(:,k)=max( rcrit_l, sig(k)**3 )
    elsewhere
      rcrit(:,k)=max( rcrit_s, sig(k)**3 )
    end where
  end do
else if ( nclddia==2 ) then
  do k = 1,kl
    where ( land(:) )
      rcrit(:,k)=rcrit_l
    elsewhere
      rcrit(:,k)=rcrit_s
    end where
  end do
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
      ! for rcit_l=.825 & nclddia=8 get rcrit=(0.826, 0.844, .871, .945, .968, .986, .993) for (200, 100, 50, 10, 5, 2, 1) km
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

select case(cmode)
case("SMITH")
  ! usual diagnostic cloud fraction

  ! Calculate cloudy fraction of grid box (stratcloud) and gridbox-mean cloud water
  ! using the triangular PDF of Smith (1990)

  ! Calculate qs at temperature tliq
  pk(:,:) = 100.*prf(:,:)
  qsi(:,:) = qsati(pk(:,:),tliq(:,:))             !Ice value
  deles(:,:) = esdiffx(tliq(:,:))                 ! MJT suggestion
  qsl(:,:) = qsi(:,:) + epsil*deles(:,:)/pk(:,:)  !qs over liquid
  qsw(:,:) = fice(:,:)*qsi(:,:) +    &
             (1.-fice(:,:))*qsl(:,:)    !Weighted qs at temperature Tliq
  do k = 1,kl
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
        stratcloud(iq,k) = max( 1.e-10, 0.5*((qc+delq)/delq)**2 )     ! for roundoff
        qcg(iq,k) = al*(qc+delq)**3/(6.*delq**2)
      else if ( qc<delq ) then
        stratcloud(iq,k) = max( 1.e-10, 1.-0.5*((qc-delq)/delq)**2 )  ! for roundoff
        qcg(iq,k) = al*(qc-(qc-delq)**3/(6.*delq**2))
      else
        stratcloud(iq,k) = 1.
        qcg(iq,k) = al*qc
      end if
    end do ! iq loop
  end do   ! k loop

  if ( nclddia==3 ) then
    do k = 1,kl
      do iq = 1,imax
        tk = min(1.,ds/(em(iq)*40000.)) ! used to increase rcrit linearly for ds<40 km
        rcrit(iq,k) = rfull(k)-tk*(rfull(k)-rcrit(iq,k))
        stratcloud(iq,k) = max(0.,min(1.,(qtot(iq,k)/qsw(iq,k)-rcrit(iq,k))/(rfull(k)-rcrit(iq,k)+0.001)))
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

case("TIEDTKE")

  ! Tiedtke prognostic cloud fraction model
  pk(:,:) = 100.*prf(:,:)
  qsi(:,:) = qsati(pk(:,:),tliq(:,:))             ! Ice value
  deles(:,:) = esdiffx(tliq(:,:))
  qsl(:,:) = qsi(:,:) + epsil*deles(:,:)/pk(:,:) ! Liquid value
  qsw(:,:) = fice(:,:)*qsi(:,:) + (1.-fice(:,:))*qsl(:,:)  ! Weighted qs at temperature Tliq

  call progcloud(tdt,qcg,qtot,ps,prf,rhoa,fice,qsw,tliq,rcrit,  &
                 dpsldt,rad_tend,trb_tend,trb_qend,stratcloud,  &
                 tiedtke_form,rkmsave,rkhsave,imax,kl)

case default
  write(6,*) "ERROR: Invalid cmode in cloudmod ",trim(cmode)
  write(6,*) "Likely to be an invalid option for ncloud"
  stop -1
  
end select

return
end subroutine newcloud

subroutine saturation_adjustment(tdt,cld_decay,vdeposition_mode,    &
                                 tliq,qtot,qcg_in,qcold,fice,       &
                                 stratcloud,prf,rhoa,               &
                                 ttg,qtg,qlg,qfg)

use const_phys                    ! Physical constants
use estab                         ! Liquid saturation function
use parm_m, only : diag           ! Model configuration

implicit none

integer, intent(in) :: vdeposition_mode
integer k, iq, imax, kl
real, dimension(:,:), intent(in) :: tliq, qtot, qcg_in, qcold, fice, stratcloud
real, dimension(:,:), intent(in) :: prf, rhoa
real, dimension(:,:), intent(out) :: ttg, qtg, qlg, qfg
real, dimension(size(tliq,1),size(tliq,2)) :: pk, qsi, deles, qcg
real, intent(in) :: tdt, cld_decay
real decayfac, Tk, qs, es, aprpr, bprpr, cice, qi0, fd, crate, qfdep
real fl
real, parameter :: rhoic = 700.
real, parameter :: cm0 = 1.e-12 !Initial crystal mass

imax = size(tliq,1)
kl = size(tliq,2)

qcg(:,:) = qcg_in(:,:)

! Assume condensation or evaporation retains ice fraction fice.
! Introduce a time-decay factor for cirrus (as suggested by results of Khvorostyanov & Sassen,
! JAS, 55, 1822-1845, 1998). Their suggested range for the time constant is 0.5 to 2 hours.
! The grid-box-mean values of qtg and ttg are adjusted later on (below).
decayfac = exp ( -tdt/cld_decay )  ! Try cld_decay=2 hrs (decayfac=0. is instant adjustment for the old scheme)
do k = 1,kl
  do iq = 1,imax
    Tk = tliq(iq,k) + hlcp*qcg(iq,k) + hlfcp*fice(iq,k)*qcg(iq,k)  
    if ( Tk>=Tice ) then
      qfg(iq,k) = fice(iq,k)*qcg(iq,k)
      qlg(iq,k) = qcg(iq,k) - qfg(iq,k)
    else                                 ! Cirrus T range
      qfg(iq,k) = qcold(iq,k)*decayfac + qcg(iq,k)*(1.-decayfac)
      qlg(iq,k) = 0.
      qcg(iq,k) = qfg(iq,k)
    end if
  end do
end do


! Do the vapour deposition calculation in mixed-phase clouds:
! Calculate deposition on cloud ice, assuming es(T) is the weighted value of the
! liquid and ice values.
if ( vdeposition_mode==0 ) then
  do k = 1,kl
    do iq = 1,imax
      Tk = tliq(iq,k) + hlcp*(qlg(iq,k)+qfg(iq,k))/max(stratcloud(iq,k),1.e-10) !T in liq cloud
      if ( stratcloud(iq,k)>1.e-10 .and. Tk<tfrz .and. qlg(iq,k)>0. ) then
        pk(iq,k) = 100.*prf(iq,k)
        qsi(iq,k) = qsati(pk(iq,k),Tk)
        deles(iq,k) = (1.-fice(iq,k))*esdiffx(Tk)
        qs        = qsi(iq,k)
        es        = qs*pk(iq,k)/0.622 !ice value
        Aprpr     = hl/(rKa*Tk)*(hls/(rvap*Tk)-1.)
        Bprpr     = rvap*Tk/((Dva/pk(iq,k))*es)
        Cice      = 1.e3*exp(12.96*deles(iq,k)/es - 0.639) !Meyers et al 1992
        qi0       = cm0*Cice/rhoa(iq,k) !Initial ice mixing ratio
        ! Next 2 lines are for assumption of fully mixed ql and qf (also a line further down).
        qi0       = max(qi0, qfg(iq,k)/stratcloud(iq,k)) !Assume all qf and ql are mixed
        fd        = 1.       !Fraction of cloud in which deposition occurs
        Crate     = 7.8*((Cice/rhoa(iq,k))**2/rhoic)**(1./3.)*deles(iq,k)/((Aprpr+Bprpr)*es)
        qfdep     = fd*stratcloud(iq,k)*sqrt(((2./3.)*Crate*tdt+qi0**(2./3.))**3)
        ! Also need this line for fully-mixed option...
        qfdep     = qfdep - qfg(iq,k)
        qfdep     = min(qfdep, qlg(iq,k))
        qlg(iq,k) = qlg(iq,k) - qfdep
        qfg(iq,k) = qfg(iq,k) + qfdep
      end if
    end do
  end do
  
else if ( vdeposition_mode==1 ) then
  ! use ql/(qf+ql) for deposition (see fl and fd)
  do k = 1,kl
    do iq = 1,imax
      Tk = tliq(iq,k) + hlcp*(qlg(iq,k)+qfg(iq,k))/max(stratcloud(iq,k),1.e-10) !T in liq cloud
      if ( stratcloud(iq,k)>1.e-10 .and. Tk<tfrz .and. qlg(iq,k)>0. ) then
        pk(iq,k) = 100.*prf(iq,k)
        qsi(iq,k) = qsati(pk(iq,k),Tk)
        deles(iq,k) = (1.-fice(iq,k))*esdiffx(Tk)
        fl        = qlg(iq,k)/(qfg(iq,k)+qlg(iq,k))
        qs        = qsi(iq,k)
        es        = qs*pk(iq,k)/0.622 !ice value
        Aprpr     = hl/(rKa*Tk)*(hls/(rvap*Tk)-1.)
        Bprpr     = rvap*Tk/((Dva/pk(iq,k))*es)
        Cice      = 1.e3*exp(12.96*deles(iq,k)/es - 0.639) !Meyers et al 1992
        qi0       = cm0*Cice/rhoa(iq,k) !Initial ice mixing ratio
        qi0       = max(qi0, qfg(iq,k)/stratcloud(iq,k)) !Assume all qf and ql are mixed
        fd        = fl      !Or, use option of adjacent ql,qf
        Crate     = 7.8*((Cice/rhoa(iq,k))**2/rhoic)**(1./3.)*deles(iq,k)/((Aprpr+Bprpr)*es)
        qfdep     = fd*stratcloud(iq,k)*sqrt(((2./3.)*Crate*tdt+qi0**(2./3.))**3)
        ! Also need this line for fully-mixed option...
        qfdep     = qfdep - qfg(iq,k)
        qfdep     = min(qfdep, qlg(iq,k))
        qlg(iq,k) = qlg(iq,k) - qfdep
        qfg(iq,k) = qfg(iq,k) + qfdep
      end if
    end do
  end do

else
  write(6,*) "ERROR: Unknow vdeposition_mode ",vdeposition_mode
  stop
end if ! vdeposition_mode==0 ..else..


! Calculate new values of vapour mixing ratio and temperature
do k = 1,kl
  do iq = 1,imax
    qtg(iq,k) = max( qtot(iq,k) - qcg(iq,k), 0. )
    ttg(iq,k) = tliq(iq,k) + hlcp*qcg(iq,k) + hlfcp*qfg(iq,k)
  end do
end do

return
end subroutine saturation_adjustment


subroutine progcloud(dt,qc,qtot,ps,press,rho,fice,qs,tliq,rcrit,   &
                     dpsldt,rad_tend,trb_tend,trb_qend,stratcloud, &
                     tiedtke_form,rkmsave,rkhsave,imax,kl)

use const_phys                    ! Physical constants
use estab                         ! Liquid saturation function
use sigs_m                        ! Atmosphere sigma levels

implicit none

integer, intent(in) :: tiedtke_form
integer, intent(in) :: imax, kl
integer iq, k
real, dimension(imax,kl), intent(inout) :: qc ! condensate = qf + ql
real, dimension(imax,kl), intent(in) :: qtot, rho, fice, qs, tliq, rcrit, press
real, dimension(imax,kl), intent(in) :: dpsldt, rkmsave, rkhsave
real, dimension(imax,kl), intent(inout) :: rad_tend, trb_tend, trb_qend
real, dimension(imax,kl), intent(inout) :: stratcloud
real, dimension(imax,kl) :: new_stratcloud, new_qcg
real, dimension(imax), intent(in) :: ps
real, dimension(imax) :: dqsdt
real, dimension(imax) :: dqs_adiabatic, dqs_radiation, dqs_turbulence
real, intent(in) :: dt
real, dimension(imax,kl) :: erosion_scale
real, parameter :: erosion_scale_d = 1.e-6
real, parameter :: erosion_scale_t = 5.e-5
real, parameter :: diff_threshold_t = 0.1
real a_dt, b_dt, cf1, cfeq, cfbar, dqs, gam, qv
real al, delq, qc_local, omega, hlrvap, aa, bb, cc

! calculate diagnostic cloud fraction using Smith.  This is uses for cloud initialisation
! (0% cloud cover) or eroding from 100% cloud cover.


! Calculate cloudy fraction of grid box (stratcloud) and gridbox-mean cloud water
! using the triangular PDF of Smith (1990)
do k = 1,kl
  do iq = 1,imax
    hlrvap = (hl+fice(iq,k)*hlf)/rvap
    dqsdt(iq) = qs(iq,k)*hlrvap/tliq(iq,k)**2
    al = 1./(1.+(hlcp+fice(iq,k)*hlfcp)*dqsdt(iq))  !Smith's notation
    qc_local = qtot(iq,k) - qs(iq,k)     ! qc_local can be negative
    delq = (1.-rcrit(iq,k))*qs(iq,k)     !UKMO style (equivalent to above)
    if ( qc_local<=-delq ) then
      new_stratcloud(iq,k) = 0.
      new_qcg(iq,k) = 0.
    else if ( qc_local<=0. ) then
      new_stratcloud(iq,k) = max( 1.e-10, 0.5*((qc_local+delq)/delq)**2 )  ! for roundoff
      new_qcg(iq,k) = al*(qc_local+delq)**3/(6.*delq**2)
    else if ( qc_local<delq ) then
      new_stratcloud(iq,k) = max( 1.e-10, 1.-0.5*((qc_local-delq)/delq)**2 )  ! for roundoff
      new_qcg(iq,k) = al*(qc_local-(qc_local-delq)**3/(6.*delq**2))
    else
      new_stratcloud(iq,k) = 1.
      new_qcg(iq,k) = al*qc_local
    end if
  end do ! iq loop
end do   ! k loop

! use diagnostic cloud as seed for new cloud, following PC2
where ( new_stratcloud(:,:)>1.e-10 .and. stratcloud(:,:)<=1.e-10 )
  stratcloud(:,:) = new_stratcloud(:,:)
  qc(:,:) = new_qcg(:,:)
end where


select case(tiedtke_form)
  case(1)
    ! erosion scale in 1/secs
    do k = 1,kl
      do iq = 1,imax  
        qv = qtot(iq,k) - qc(iq,k)        
        if ( 0.5*(rkmsave(iq,k)+rkhsave(iq,k)) > diff_threshold_t ) then
          ! Turbulence          
          erosion_scale(iq,k) = erosion_scale_t*stratcloud(iq,k)*max(qs(iq,k)-qv,0.)/max(qc(iq,k),1.e-20)
        else
          ! Background
          erosion_scale(iq,k) = erosion_scale_d*stratcloud(iq,k)*max(qs(iq,k)-qv,0.)/max(qc(iq,k),1.e-20)
        end if
      end do  
    end do
    ! Aerosols
    !do k = 1,kl
    !  edum = drop(:,k)*0.4/150. + 0.8
    !  erosion_scale(:,k) = edum*erosion_scale(:,k)
    !end do

  case(2)  
    ! PC2 version  

    ! da_erosion = -0.5*(1-stratcloud)^2 * Qc * (1/bs)*dbsdt / aL*(qsat - qv )
    ! Qc = aL*(qt-qsat)
    ! (1/bs)*dbsdt = -2.25e-5*exp(-3.1*Qc/(aL*qsat)
    ! aL = 1/( 1 + dqsdT*L/Cp ) = 1/(1+gamma)  
    do k = 1,kl
      do iq = 1,imax  
        qv = qtot(iq,k) - qc(iq,k)   
        erosion_scale(iq,k) = 0.5*(1.-stratcloud(iq,k))**2/max(stratcloud(iq,k),1.e-10)  &
            *(qtot(iq,k)-qs(iq,k))*2.25e-5*exp(-3.1*(qtot(iq,k)-qs(iq,k))/qs(iq,k))      &
            /max(qs(iq,k)-qv,1.e-8)
        erosion_scale(iq,k) = max( erosion_scale(iq,k), 0. )
      end do  
    end do
    
  case default
    write(6,*) "ERROR: Unknown option for tiedtke_form = ",tiedtke_form
    stop
end select      

  
! Tiedtke scheme

! calculate dqs = ( (omega+grav*Mc)/(cp*rho) + nettend )*dqsdT*dt
!                 -------------------------------------------------------
!                 1 + (stratcloud + 0.5*da)*gamma
! MJT notes - GFDL AM adds (stratcloud+0.5*at*da)*gamma term

! Change in saturated volume fraction
! da = -0.5*(1.-stratcloud)^2*dqs/(qs-qv)
! MJT notes - Tiedtke 93 does not use 0.5

! gamma = L/cp*dqsdT

! Follow GFDL AM approach since da=da(dqs), hence need to solve the above
! quadratic equation for dqs if da/=0

! dqs**2*AA + dqs*BB + CC = 0
! AA = 0.25*gamma*(1-stratcloud)^2/(qs-qv)
! BB = -(1+gamma*stratcloud)
! CC = ((omega + grav*mflx)/(cp*rho)+netten)*dqsdT*dt

do k = 1,kl

  do iq = 1,imax
    stratcloud(iq,k) = max( min( stratcloud(iq,k), 1. ), 0. )
    ! calculate dqs/dT
    hlrvap = (hl+fice(iq,k)*hlf)/rvap
    dqsdT(iq) = qs(iq,k)*hlrvap/(tliq(iq,k)**2)
  end do  
  
  ! tiedtke_form==1
  do iq = 1,imax
    omega = ps(iq)*dpsldt(iq,k)   
    dqs_adiabatic(iq) = omega/(cp*rho(iq,k))*dqsdT(iq)
    dqs_radiation(iq) = rad_tend(iq,k)*dqsdT(iq)
    dqs_turbulence(iq) = trb_tend(iq,k)*dqsdT(iq)
  end do    
  if ( tiedtke_form==2 ) then
    do iq = 1,imax  
      dqs_turbulence(iq) = dqs_turbulence(iq) - trb_qend(iq,k)
    end do  
  end if  
  !dqs_convection(:) = grav*cmflx(:,k)/(cp*rho(:,k))*dqsdT
  
  do iq = 1,imax
    qv = qtot(iq,k) - qc(iq,k)  
    gam = (hlcp+fice(iq,k)*hlfcp)*dqsdT(iq)  
      
    !cc = ((omega + grav*cmflx(:,k))/(cp*rho(:,k))+nettend(:,k))*dqsdT
    aa = 0.25*gam*(1.-stratcloud(iq,k))**2/max( qs(iq,k)-qv, 1.e-8 )
    bb = -(1.+gam*stratcloud(iq,k))  
    cc = dqs_adiabatic(iq) + dqs_radiation(iq) + dqs_turbulence(iq) ! neglect dqs_convection
    dqs = 2.*cc/( -bb + sqrt( bb**2 - 4.*aa*cc ) ) ! alternative form of quadratic equation
    
    !da = -2.*aa*dqs/gam
    !da(:) = -0.5*dqs*(1.-stratcloud(:,k))**2/max( qs(:,k)-qv, 1.e-8 )

    ! Large scale cloud formation via condensation (A)
    ! Large scale cloud destruction via erosion (B)
    if ( qc(iq,k)>1.e-8 ) then
      !a_dt = da/max(1.-stratcloud(:,k),1.e-10)
      a_dt = max( -0.5*dqs*(1.-stratcloud(iq,k))/max( qs(iq,k)-qv, 1.e-8 ), 0. )
      !a_dt = da_erosion/max(stratcloud(:,k),1.e-10)
    else
      a_dt = 0.  
    end if
    b_dt = erosion_scale(iq,k)

    ! Integrate
    !   dcf/dt = (1-cf)*A - cf*B
    ! to give (use cf' = A-cf*(A+B))
    !   cf(t=1) = cfeq + (cf(t=0) - cfeq)*exp(-(A+B)*dt)
    !   cfeq = A/(A+B)
    ! Average cloud fraction over the interval t=tau to t=tau+1
    !   cfbar = cfeq - (cf(t=1) - cf(t=0))/((A+B)*dt)
    ! cfeq is the equilibrum cloud fraction that is approached with
    ! a time scale of 1/(A+B)
    if ( a_dt+b_dt>1.e-20 ) then
      cfeq  = a_dt/(a_dt+b_dt)
      cf1   = cfeq + (stratcloud(iq,k) - cfeq)*exp(-(a_dt+b_dt)*dt)
      cfbar = cfeq + (stratcloud(iq,k) - cf1)/((a_dt+b_dt)*dt)
    else
      cfeq  = stratcloud(iq,k)
      cf1   = stratcloud(iq,k)
      cfbar = stratcloud(iq,k)
    end if
  
    ! Change in condensate
    ! dqc = -dqs*(stratcloud+0.5*da)
    !qc(:,k) = min( max( qc(:,k) - dqs*(cf1+0.5*da)*dt, 0.), qtot(:,k) )
    qc(iq,k) = max( qc(iq,k) - cfbar*dqs, 0. )

    stratcloud(iq,k) = max( min( cf1, 1. ), 0. )
  
  end do ! i=1,imax

  ! Change in cloud fraction
  where ( qc(:,k)>0. )
    stratcloud(:,k) = max( stratcloud(:,k), 1.e-10 )
  elsewhere
    ! MJT notes - cloud fraction is maintained (da=0.) while condesate evaporates (dqc<0.) until
    ! the condesate dissipates
    stratcloud(:,k) = 0.
    qc(:,k) = 0.
  end where

end do ! k = 1,kl

return
end subroutine progcloud

subroutine convectivecloudfrac(clcon,kbsav,ktsav,condc,acon,bcon,cldcon)

use parm_m           ! Model configuration

implicit none

integer k, kl, iq
real, intent(in) :: acon, bcon
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
where ( ktsav<kl-1 )
  cldcon_temp = min( acon+bcon*log(1.+condc*86400./dt), 0.8 ) !NCAR
elsewhere
  cldcon_temp = 0.
end where
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
  do iq = 1,size(clcon,1)
    if ( k<kbsav(iq)+1 .or. k>ktsav(iq) ) then
      clcon(iq,k) = 0.
    else
      clcon(iq,k) = cldcon_local(iq)
    end if
  end do
end do

return
end subroutine convectivecloudfrac

end module cloudmod
