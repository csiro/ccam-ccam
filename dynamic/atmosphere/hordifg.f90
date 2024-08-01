! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2024 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
! This subroutine calculates horizontal diffusion for atmospheric fields

!     this horizontal diffusion routine allows compensation for
!     being on sigma surfaces through reductions of kh over orography
!     this is called for -nhor.ge.50; the value of nhor (in metres) gives the
!     scaling value for deltaz
!     in namelist khdif is fujio's a**2, e.g. 4.
    
!     nhorx used in hordif  ! previous code effectively has nhorx=0
!           = 1 u, v, T, qg  diffusion reduced near mountains
!           = 4              diffusion not reduced near mountains    
!           = 7 u, v, T, qg  diffusion reduced near mountains (bottom 2/3 only)

!     khor: if non-zero increases the diffusion in the
!     upper model levels   (never used these days)
!         khor > 0 progessively doubles for khor upper levels
!         khor < 0 linearly increases above sig = 0.2 (up to factor of 1-khor)
    
!     It has -ve nhorps option:
!        nhorps=0  does only T, qg, TKE, u & v            horiz diff.
!        nhorps=-1 does only T, qg & TKE                  horiz diff.
!        nhorps=-2 does only u &  v                       horiz diff.
!        nhorps=-3 does only qg                           horiz diff.
!        nhorps=-4 does only T, qg, cloud, TKE & aerosols horiz diff.
!        nhorps=-5 does only T                            horiz diff.
!        nhorps=-6 does only T & qg                       horiz diff.
  
module hordifg_m

private
public hordifgt

contains
    
subroutine hordifgt

use aerointerface
use arrays_m
use cc_mpi
use cfrac_m
use const_phys
use dpsdt_m
use indices_m
use kuocom_m
use liqwpar_m
use map_m
use newmpar_m
use nharrs_m
use parm_m
use parmdyn_m
use parmhdff_m
use savuvt_m
use sigs_m
use tkeeps, only : tke,eps,shear,mintke,mineps,cm0,minl,maxl, &
                   u_ema,v_ema,w_ema,update_ema
use vecsuv_m
use vvel_m

implicit none

real, dimension(ifull+iextra,kl) :: uc, vc, wc
real, dimension(ifull+iextra,kl) :: uav, vav
real, dimension(ifull+iextra,kl) :: xfact, yfact, t_kh
real, dimension(ifull,kl) :: ww, dwdx, dwdy
real, dimension(ifull) :: emi
real, dimension(ifull,kl) :: zg
real, dimension(ifull) :: zgh_a, zgh_b
real, dimension(ifull) :: ptemp, tx_fact, ty_fact
real dudx, dudy, dvdx, dvdy, dudz, dvdz
real r1, r2, cc, base
real ucc, vcc, wcc
real delphi, hdif
real tv
integer iq, k, nhora, nhorx
integer nstart, nend, nt, ntr, kmax
integer, parameter :: nf=2

!     set up topography reduction factors for each type of location
!     expect power nf to be about 1 or 2 (see data statement)

nhorx = 0
delphi = 1.e6  ! turns off reduction (can also use nhorx=4)
if ( abs(nhor)>=50 ) then
  nhora = 10*(abs(nhor)/10)    ! e.g. 150  for nhor=-157
  nhorx = abs(nhor) - nhora    ! e.g.   7  for nhor=-157
  delphi = nhora*grav
endif

if ( nhorx==1 ) then
  kmax = kl
else if ( nhorx>=7 ) then
  kmax = 1
  do while ( sig(kmax)>=0.25 .and. kmax<kl )
    kmax = kmax + 1
  end do
else
  kmax = 0
end if

do k = 1,kl
  t_kh(:,k) = 0.
  xfact(:,k) = 0.
  yfact(:,k) = 0.
  uav(:,k) = 0.
  vav(:,k) = 0.
  uc(:,k) = 0.
  vc(:,k) = 0.
  wc(:,k) = 0.  
end do
emi(1:ifull) = ps(1:ifull)/em(1:ifull)
!ptemp(1:ifull) = ps(1:ifull)**.286
ptemp(1:ifull) = 26.915348*exp(0.286*psl(1:ifull))
tx_fact(1:ifull) = 1./(1.+(abs(zs(ie)-zs(1:ifull))/delphi)**nf)
ty_fact(1:ifull) = 1./(1.+(abs(zs(in)-zs(1:ifull))/delphi)**nf)

if ( diag .and. mydiag ) then
  write(6,*) 'hordifgt u ',(u(idjd,k),k=1,kl)
  write(6,*) 'hordifgt v ',(v(idjd,k),k=1,kl)
  write(6,*) 'ax,ay,az ',ax(idjd),ay(idjd),az(idjd)
  write(6,*) 'bx,by,bz ',bx(idjd),by(idjd),bz(idjd)
  write(6,*) 'zs,zse,zsn ',zs(idjd),zs(ie(idjd)),zs(in(idjd))
  write(6,*) 'tx_,ty_ ',tx_fact(idjd),ty_fact(idjd)
endif

! weighted horizontal velocities
do k = 1,kl        
  uav(1:ifull,k) = av_vmod*u(1:ifull,k) + (1.-av_vmod)*savu(1:ifull,k)
  vav(1:ifull,k) = av_vmod*v(1:ifull,k) + (1.-av_vmod)*savv(1:ifull,k)
end do

! Calculate shear for tke
if ( nvmix==6 .or. nvmix==9 ) then

  ! calculate vertical velocity in m/s
  ! omega = ps*dpsldt
  ! ww = -R/g * (T+Tnhs) * dpsldt/sig
  do k = 1,kl
    do iq = 1,ifull  
      tv = t(iq,k)*(1.+0.61*qg(iq,k)-qlg(iq,k)-qfg(iq,k))
      !tnhs(iq,k) = (phi_nh(iq,k)-phi_nh(iq,k-1)-betm(k)*tnhs(iq,k-1))/bet(k)
      ww(iq,k) = (dpsldt(iq,k)/sig(k)-dpsdt(iq)/(864.*ps(iq))) &
                 *(-rdry/grav)*tv
    end do  
  end do

  ! calculate height on full levels (hydrostatic terms)
  zg(:,1) = (zs(1:ifull)+bet(1)*t(1:ifull,1))/grav
  do k = 2,kl
    zg(:,k) = zg(:,k-1) + (bet(k)*t(1:ifull,k)+betm(k)*t(1:ifull,k-1))/grav
  end do ! k  loop
  !zg = zg + phi_nh(1:ifull,:)/grav

  ! time averaging of source terms for tke-eps
  call update_ema(ww,w_ema,dt)
  call update_ema(uav,u_ema,dt)
  call update_ema(vav,v_ema,dt)
  
  call bounds(w_ema)
  do k = 1,kl
    do iq = 1,ifull  
      dwdx(iq,k) = 0.5*(w_ema(ie(iq),k)-w_ema(iw(iq),k))*em(iq)/ds
      dwdy(iq,k) = 0.5*(w_ema(in(iq),k)-w_ema(is(iq),k))*em(iq)/ds
    end do
  end do
  do k = 1,kmax
    do iq = 1,ifull  
      dwdx(iq,k) = dwdx(iq,k)*tx_fact(iq)
      dwdy(iq,k) = dwdy(iq,k)*ty_fact(iq)
    end do
  end do
  
  ! calculate vertical gradients
  do iq = 1,ifull
    zgh_b(iq) = ratha(1)*zg(iq,2) + rathb(1)*zg(iq,1) ! upper half level
    r1 = u_ema(iq,1)
    r2 = ratha(1)*u_ema(iq,2) + rathb(1)*u_ema(iq,1)          
    dudz = (r2-r1)/(zgh_b(iq)-zg(iq,1))
    r1 = v_ema(iq,1)
    r2 = ratha(1)*v_ema(iq,2) + rathb(1)*v_ema(iq,1)          
    dvdz = (r2-r1)/(zgh_b(iq)-zg(iq,1))
    shear(iq,1) = (dudz+dwdx(iq,1))**2 + (dvdz+dwdy(iq,1))**2
  end do
  do k = 2,kl-1
    do iq = 1,ifull
      zgh_a(iq) = zgh_b(iq) ! lower half level
      zgh_b(iq) = ratha(k)*zg(iq,k+1) + rathb(k)*zg(iq,k)     ! upper half level
      r1 = ratha(k-1)*u_ema(iq,k) + rathb(k-1)*u_ema(iq,k-1)
      r2 = ratha(k)*u_ema(iq,k+1) + rathb(k)*u_ema(iq,k)          
      dudz = (r2-r1)/(zgh_b(iq)-zgh_a(iq))
      r1 = ratha(k-1)*v_ema(iq,k) + rathb(k-1)*v_ema(iq,k-1)
      r2 = ratha(k)*v_ema(iq,k+1) + rathb(k)*v_ema(iq,k)          
      dvdz = (r2-r1)/(zgh_b(iq)-zgh_a(iq))
      shear(iq,k) = (dudz+dwdx(iq,k))**2 + (dvdz+dwdy(iq,k))**2
    end do
  end do
  do iq = 1,ifull
    zgh_a(iq) = zgh_b(iq) ! lower half level
    r1 = ratha(kl-1)*u_ema(iq,kl) + rathb(kl-1)*u_ema(iq,kl-1)
    r2 = u_ema(iq,kl)          
    dudz = (r2-r1)/(zg(iq,kl)-zgh_a(iq))
    r1 = ratha(kl-1)*v_ema(iq,kl) + rathb(kl-1)*v_ema(iq,kl-1)
    r2 = v_ema(iq,kl)          
    dvdz = (r2-r1)/(zg(iq,kl)-zgh_a(iq))
    shear(iq,kl) = (dudz+dwdx(iq,kl))**2 + (dvdz+dwdy(iq,kl))**2
  end do
  
end if ! nvmix=6 .or. nvmix==9
      
! usual deformation for nhorjlm=1 or nhorjlm=2
if ( nhorjlm==1 .or. nhorjlm==2 .or. nhorps==0 .or. nhorps==-2 ) then 
  do k = 1,kl
    ! in hordifgt, need to calculate Cartesian components 
    uc(1:ifull,k) = ax(1:ifull)*u(1:ifull,k) + bx(1:ifull)*v(1:ifull,k)
    vc(1:ifull,k) = ay(1:ifull)*u(1:ifull,k) + by(1:ifull)*v(1:ifull,k)
    wc(1:ifull,k) = az(1:ifull)*u(1:ifull,k) + bz(1:ifull)*v(1:ifull,k)
  end do
  call bounds(uc)
  call bounds(vc)
  call bounds(wc)
end if

! apply horz diffusion
select case(nhorjlm)
  case(0)
    ! This is based on 2D Smagorinsky closure
    ! uses (dv/dx+du/dy)**2 + .5*(du/dx)**2 + .5*(dv/dy)**2
    ! following Kikuchi et al. 1981      now Smag. Wed  04-30-1997
    ! N.B. original Smag. had m on top (in D formulae) and khdif=3.2
    ! More recently (21/9/00) I think original Smag has khdif=0.8
    ! Smag's actual diffusion also differentiated Dt and Ds
    ! t_kh is kh at t points
    call boundsuv(uav,vav,allvec=.true.)
    do k = 1,kl
      hdif = dt*hdiff(k) ! N.B.  hdiff(k)=khdif*.1  
      do iq = 1,ifull
        dudx = 0.5*(uav(ieu(iq),k)-uav(iwu(iq),k))*em(iq)/ds
        dudy = 0.5*(uav(inu(iq),k)-uav(isu(iq),k))*em(iq)/ds
        dvdx = 0.5*(vav(iev(iq),k)-vav(iwv(iq),k))*em(iq)/ds
        dvdy = 0.5*(vav(inv(iq),k)-vav(isv(iq),k))*em(iq)/ds
        r1 = (dvdx+dudy)**2 + 0.5*dudx**2 + 0.5*dvdy**2
        t_kh(iq,k) = sqrt(r1)*hdif*emi(iq)
      end do
    end do
             
  case(1)
    ! jlm deformation scheme using 3D uc, vc, wc
    do k = 1,kl
      hdif = dt*hdiff(k)/ds ! N.B.  hdiff(k)=khdif*.1  
      do iq = 1,ifull
        cc = (uc(ie(iq),k)-uc(iw(iq),k))**2 + (uc(in(iq),k)-uc(is(iq),k))**2 + &
             (vc(ie(iq),k)-vc(iw(iq),k))**2 + (vc(in(iq),k)-vc(is(iq),k))**2 + &
             (wc(ie(iq),k)-wc(iw(iq),k))**2 + (wc(in(iq),k)-wc(is(iq),k))**2
        ! N.B. using double grid length
        t_kh(iq,k)= 0.5*sqrt(cc)*hdif*ps(iq) ! this one without em in D terms
      end do
    end do

  case(2)
    ! jlm deformation scheme using 3D uc, vc, wc and omega (1st rough scheme)
    do k = 1,kl
      do iq = 1,ifull
        hdif = dt*hdiff(k)/ds ! N.B.  hdiff(k)=khdif*.1
        cc = (uc(ie(iq),k)-uc(iw(iq),k))**2 + (uc(in(iq),k)-uc(is(iq),k))**2 + &
             (vc(ie(iq),k)-vc(iw(iq),k))**2 + (vc(in(iq),k)-vc(is(iq),k))**2 + &
             (wc(ie(iq),k)-wc(iw(iq),k))**2 + (wc(in(iq),k)-wc(is(iq),k))**2 + &
             .01*(dpsldt(ie(iq),k)*ps(ie(iq))-dpsldt(iw(iq),k)*ps(iw(iq)))**2+ &
             .01*(dpsldt(in(iq),k)*ps(in(iq))-dpsldt(is(iq),k)*ps(is(iq)))**2 
        ! approx 1 Pa/s = .1 m/s     
        ! N.B. using double grid length
        t_kh(iq,k)= 0.5*sqrt(cc)*hdif*ps(iq) ! this one without em in D terms
      end do
    enddo

  case(3)
    ! K-eps model + Smag
    call boundsuv(uav,vav,allvec=.true.)
    do k = 1,kl
      do iq = 1,ifull
        hdif = dt*hdiff(k) ! N.B.  hdiff(k)=khdif*.1
        dudx = 0.5*(uav(ieu(iq),k)-uav(iwu(iq),k))*em(iq)/ds
        dudy = 0.5*(uav(inu(iq),k)-uav(isu(iq),k))*em(iq)/ds
        dvdx = 0.5*(vav(iev(iq),k)-vav(iwv(iq),k))*em(iq)/ds
        dvdy = 0.5*(vav(inv(iq),k)-vav(isv(iq),k))*em(iq)/ds
        r1 = (dudx-dvdy)**2+(dvdx+dudy)**2
        t_kh(iq,k) = sqrt(r1)*hdif*emi(iq)
        t_kh(iq,k)=max( t_kh(iq,k), tke(iq,k)**2/eps(iq,k)*dt*cm0*emi(iq) )
      end do
    end do

  case DEFAULT
    write(6,*) "ERROR: Unknown option nhorjlm=",nhorjlm
    call ccmpi_abort(-1)
end select


call bounds(t_kh,nehalf=.true.)
do k = 1,kl
  xfact(1:ifull,k) = (t_kh(ie,k)+t_kh(1:ifull,k))*.5
  yfact(1:ifull,k) = (t_kh(in,k)+t_kh(1:ifull,k))*.5
end do  
do k = 1,kmax
  xfact(1:ifull,k) = xfact(1:ifull,k)*tx_fact(1:ifull)
  yfact(1:ifull,k) = yfact(1:ifull,k)*ty_fact(1:ifull)
end do    
call boundsuv(xfact,yfact,stag=-9) ! MJT - can use stag=-9 option that will
                                   ! only update iwu and isv values


! perform diffusion ---------------------------------------------------

if ( nhorps==0 .or. nhorps==-1 .or. nhorps==-4 .or. nhorps==-5 .or. nhorps==-6 ) then
  do k = 1,kl
    t(1:ifull,k) = t(1:ifull,k)/ptemp(1:ifull)
  end do
  call bounds(t)
end if
if ( nhorps==0 .or. nhorps==-1 .or. nhorps==-3 .or. nhorps==-4 .or. nhorps==-6 ) then
  call bounds(qg)
end if
if ( nhorps==-4 .and. ldr/=0 ) then  
  call bounds(qlg)
  call bounds(qfg)
  call bounds(stratcloud)
  if ( ncloud>=100 .and. ncloud<200 ) then
    call bounds(ni)
  end if
end if
if ( (nhorps==0.or.nhorps==-1.or.nhorps==-4) .and. (nvmix==6.or.nvmix==9) ) then
  call bounds(eps)
  call bounds(tke)
end if
if ( nhorps==-4 .and. abs(iaero)>=2 ) then
  call bounds(xtg)  
end if


!$omp parallel
!$omp sections

! momentum U, V, W - bounds updated above
!$omp section
if ( nhorps==0 .or. nhorps==-2 ) then ! for nhorps=-1,-3,-4 don't diffuse u,v
  call hordifgt_work(uc,xfact,yfact,emi)
end if
!$omp section
if ( nhorps==0 .or. nhorps==-2 ) then ! for nhorps=-1,-3,-4 don't diffuse u,v
  call hordifgt_work(vc,xfact,yfact,emi)
end if
!$omp section
if ( nhorps==0 .or. nhorps==-2 ) then ! for nhorps=-1,-3,-4 don't diffuse u,v
  call hordifgt_work(wc,xfact,yfact,emi)
end if  

! potential temperture and water vapour
!$omp section
if ( nhorps==0 .or. nhorps==-1 .or. nhorps==-4 .or. nhorps==-5 .or. nhorps==-6 ) then
  call hordifgt_work(t,xfact,yfact,emi)
end if
!$omp section
if ( nhorps==0 .or. nhorps==-1 .or. nhorps==-3 .or. nhorps==-4 .or. nhorps==-6 ) then
  call hordifgt_work(qg,xfact,yfact,emi)
end if  

! cloud liquid & frozen water plus cloud fraction
!$omp section
if ( nhorps==-4 .and. ldr/=0 ) then  
  call hordifgt_work(qlg,xfact,yfact,emi)
end if
!$omp section
if ( nhorps==-4 .and. ldr/=0 ) then  
  call hordifgt_work(qfg,xfact,yfact,emi)
end if
!$omp section
if ( nhorps==-4 .and. ldr/=0 ) then  
  call hordifgt_work(stratcloud,xfact,yfact,emi)
end if
!$omp section
if ( nhorps==-4 .and. ldr/=0 ) then  
  if ( ncloud>=100 .and. ncloud<200 ) then
    call hordifgt_work(ni,xfact,yfact,emi)
  end if
end if

! tke and eps
!$omp section
if ( (nhorps==0.or.nhorps==-1.or.nhorps==-4) .and. (nvmix==6.or.nvmix==9) ) then
  call hordifgt_work(eps,xfact,yfact,emi)
end if
!$omp section
if ( (nhorps==0.or.nhorps==-1.or.nhorps==-4) .and. (nvmix==6.or.nvmix==9) ) then
  call hordifgt_work(tke,xfact,yfact,emi)
end if

!$omp end sections nowait

! prgnostic aerosols (disabled by default)
if ( nhorps==-4 .and. abs(iaero)>=2 ) then
  !$omp do schedule(static) private(ntr)
  do ntr = 1,naero
    call hordifgt_work(xtg(:,:,ntr),xfact,yfact,emi)
  end do
  !$omp end do nowait
end if  ! (nhorps==-4.and.abs(iaero)>=2)  

!$omp end parallel


if ( nhorps==0 .or. nhorps==-2 ) then ! for nhorps=-1,-3,-4 don't diffuse u,v
  do k = 1,kl
    u(1:ifull,k) = ax(1:ifull)*uc(1:ifull,k) &
                 + ay(1:ifull)*vc(1:ifull,k) &
                 + az(1:ifull)*wc(1:ifull,k)
    v(1:ifull,k) = bx(1:ifull)*uc(1:ifull,k) &
                 + by(1:ifull)*vc(1:ifull,k) &
                 + bz(1:ifull)*wc(1:ifull,k)
  end do
end if  
if ( nhorps==0 .or. nhorps==-1 .or. nhorps==-4 .or. nhorps==-5 .or. nhorps==-6 ) then
  do k = 1,kl
    t(1:ifull,k) = t(1:ifull,k)*ptemp(1:ifull)  
  end do
end if

if ( diag .and. mydiag ) then
  do k = 1,kl
    write(6,*) 'k,id,jd,idjd ',k,id,jd,idjd
    write(6,*) 'k, xfact, xfactw ',k,xfact(idjd,k),xfact(iwu(idjd),k)
    write(6,*) 'k, yfact, yfacts ',k,yfact(idjd,k),yfact(isv(idjd),k)
    write(6,*) 'k, uc,uce,ucw,ucn,ucs ',k,uc(idjd,k),uc(ie(idjd),k),uc(iw(idjd),k),uc(in(idjd),k),uc(is(idjd),k)
    write(6,*) 'k,u,v ',k,u(idjd,k),v(idjd,k)
  end do
endif

return
end subroutine hordifgt

subroutine hordifgt_work(work,xfact,yfact,emi)

use indices_m
use newmpar_m

implicit none

integer k, iq
real, dimension(ifull+iextra,kl), intent(in) :: xfact, yfact
real, dimension(ifull), intent(in) :: emi
real, dimension(ifull+iextra,kl), intent(inout) :: work
real, dimension(ifull) :: ans
real base, xfact_iwu, yfact_isv

do k = 1,kl
  do iq = 1,ifull  
    xfact_iwu = xfact(iwu(iq),k)
    yfact_isv = yfact(isv(iq),k)
    base = emi(iq)+xfact(iq,k)+xfact_iwu  &
                  +yfact(iq,k)+yfact_isv
    ans(iq) = ( emi(iq)*work(iq,k) +               &
                xfact(iq,k)*work(ie(iq),k) +       &
                xfact_iwu*work(iw(iq),k) +         &
                yfact(iq,k)*work(in(iq),k) +       &
                yfact_isv*work(is(iq),k) )         &
             / base 
  end do
  do iq = 1,ifull
    work(iq,k) = ans(iq)
  end do
end do

return
end subroutine hordifgt_work

end module hordifg_m
