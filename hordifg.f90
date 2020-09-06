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
    
subroutine hordifgt

use aerosolldr
use arrays_m
use cc_mpi
use cfrac_m
use const_phys
use dpsdt_m
use indices_m
use liqwpar_m
use map_m
use newmpar_m
use nharrs_m
use parm_m
use parmdyn_m
use parmhdff_m
use savuvt_m
use sigs_m
use tkeeps, only : tke,eps,shear,mintke,mineps,cm0,minl,maxl
use vecsuv_m
use vvel_m

implicit none

include 'kuocom.h'

real, dimension(ifull+iextra,kl,3) :: work
real, dimension(ifull+iextra,kl) :: uc, vc, wc
real, dimension(ifull+iextra,kl) :: uav, vav
real, dimension(ifull+iextra,kl) :: xfact, yfact, t_kh
real, dimension(ifull) :: emi
real, dimension(ifull,kl) :: zg
real, dimension(ifull) :: zgh_a, zgh_b
real, dimension(ifull) :: ptemp, tx_fact, ty_fact
real dudx, dudy, dvdx, dvdy, dudz, dvdz
real r1, r2, cc, base
real ucc, vcc, wcc
real delphi, hdif
integer iq, k, nhora, nhorx
integer nstart, nend, nt, ntr
integer, save :: kmax=-1
integer, parameter :: nf=2

nhorx = 0

if ( kmax<0 ) then
  kmax = 1
  do while ( sig(kmax)>=0.25 .and. kmax<kl )
    kmax = kmax + 1
  end do
end if

!     set up topography reduction factors for each type of location
!     expect power nf to be about 1 or 2 (see data statement)

delphi = 1.e6  ! turns off reduction (can also use nhorx=4)
if ( abs(nhor)>=50 ) then
  nhora = 10*(abs(nhor)/10)    ! e.g. 150  for nhor=-157
  nhorx = abs(nhor) - nhora    ! e.g.   7  for nhor=-157
  delphi = nhora*grav
endif

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
ptemp(1:ifull) = ps(1:ifull)**.286
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

do k = 1,kl        
  ! weighted horizontal velocities
  uav(1:ifull,k) = av_vmod*u(1:ifull,k) + (1.-av_vmod)*savu(1:ifull,k)
  vav(1:ifull,k) = av_vmod*v(1:ifull,k) + (1.-av_vmod)*savv(1:ifull,k)
end do

! Calculate shear for tke
if ( nvmix==6 ) then

  ! calculate height on full levels
  zg(:,1) = (zs(1:ifull)+bet(1)*t(1:ifull,1))/grav
  do k = 2,kl
    zg(:,k) = zg(:,k-1) + (bet(k)*t(1:ifull,k)+betm(k)*t(1:ifull,k-1))/grav
  end do ! k  loop
        
  ! calculate vertical gradients
  do concurrent (iq = 1:ifull)
    zgh_b(iq) = ratha(1)*zg(iq,2) + rathb(1)*zg(iq,1) ! upper half level
    r1 = uav(iq,1)
    r2 = ratha(1)*uav(iq,2) + rathb(1)*uav(iq,1)          
    dudz = (r2-r1)/(zgh_b(iq)-zg(iq,1))
    r1 = vav(iq,1)
    r2 = ratha(1)*vav(iq,2) + rathb(1)*vav(iq,1)          
    dvdz = (r2-r1)/(zgh_b(iq)-zg(iq,1))
    shear(iq,1) = dudz**2 + dvdz**2
  end do
  do k = 2,kl-1
    do concurrent (iq = 1:ifull)
      zgh_a(iq) = zgh_b(iq) ! lower half level
      zgh_b(iq) = ratha(k)*zg(iq,k+1) + rathb(k)*zg(iq,k) ! upper half level
      r1 = ratha(k-1)*uav(iq,k) + rathb(k-1)*uav(iq,k-1)
      r2 = ratha(k)*uav(iq,k+1) + rathb(k)*uav(iq,k)          
      dudz = (r2-r1)/(zgh_b(iq)-zgh_a(iq))
      r1 = ratha(k-1)*vav(iq,k) + rathb(k-1)*vav(iq,k-1)
      r2 = ratha(k)*vav(iq,k+1) + rathb(k)*vav(iq,k)          
      dvdz = (r2-r1)/(zgh_b(iq)-zgh_a(iq))
      shear(iq,k) = dudz**2 + dvdz**2
    end do
  end do
  do concurrent (iq = 1:ifull)
    zgh_a(iq) = zgh_b(iq) ! lower half level
    r1 = ratha(kl-1)*uav(iq,kl) + rathb(kl-1)*uav(iq,kl-1)
    r2 = uav(iq,kl)          
    dudz = (r2-r1)/(zg(iq,kl)-zgh_a(iq))
    r1 = ratha(kl-1)*vav(iq,kl) + rathb(kl-1)*vav(iq,kl-1)
    r2 = vav(iq,kl)          
    dvdz = (r2-r1)/(zg(iq,kl)-zgh_a(iq))
    shear(iq,k) = dudz**2 + dvdz**2
  end do
end if ! nvmix=6
      
! usual deformation for nhorjlm=1 or nhorjlm=2
if ( nhorjlm==1 .or. nhorjlm==2 .or. nhorps==0 .or. nhorps==-2 ) then 
  do k = 1,kl
    ! in hordifgt, need to calculate Cartesian components 
    work(1:ifull,k,1) = ax(1:ifull)*u(1:ifull,k) + bx(1:ifull)*v(1:ifull,k)
    work(1:ifull,k,2) = ay(1:ifull)*u(1:ifull,k) + by(1:ifull)*v(1:ifull,k)
    work(1:ifull,k,3) = az(1:ifull)*u(1:ifull,k) + bz(1:ifull)*v(1:ifull,k)
  end do
  call bounds(work(:,:,1:3))
  uc(:,:) = work(:,:,1)
  vc(:,:) = work(:,:,2)
  wc(:,:) = work(:,:,3)
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
    !$acc parallel loop collapse(2) copyin(hdiff,uav,vav,em,emi) copyout(t_kh(1:ifull,1:kl)) &
    !$acc   present(ieu,iwu,inu,isu,iev,iwv,inv,isv)
    do concurrent (k = 1:kl)
      do concurrent (iq = 1:ifull)
        hdif=dt*hdiff(k) ! N.B.  hdiff(k)=khdif*.1
        dudx = 0.5*(uav(ieu(iq),k)-uav(iwu(iq),k))*em(iq)/ds
        dudy = 0.5*(uav(inu(iq),k)-uav(isu(iq),k))*em(iq)/ds
        dvdx = 0.5*(vav(iev(iq),k)-vav(iwv(iq),k))*em(iq)/ds
        dvdy = 0.5*(vav(inv(iq),k)-vav(isv(iq),k))*em(iq)/ds
        r1 = (dvdx+dudy)**2 + 0.5*dudx**2 + 0.5*dvdy**2
        t_kh(iq,k)=sqrt(r1)*hdif*emi(iq)
      end do
    end do
    !$acc end parallel loop
             
  case(1)
    ! jlm deformation scheme using 3D uc, vc, wc
    !$acc parallel loop collapse(2) copyin(hdiff,uc,vc,wc,ps) copyout(t_kh(1:ifull,1:kl)) &
    !$acc   present(ie,iw,in,is)
    do concurrent (k = 1:kl)
      do concurrent (iq = 1:ifull)
        hdif = dt*hdiff(k)/ds ! N.B.  hdiff(k)=khdif*.1
        cc = (uc(ie(iq),k)-uc(iw(iq),k))**2 + (uc(in(iq),k)-uc(is(iq),k))**2 + &
             (vc(ie(iq),k)-vc(iw(iq),k))**2 + (vc(in(iq),k)-vc(is(iq),k))**2 + &
             (wc(ie(iq),k)-wc(iw(iq),k))**2 + (wc(in(iq),k)-wc(is(iq),k))**2
        ! N.B. using double grid length
        t_kh(iq,k)= .5*sqrt(cc)*hdif*ps(iq) ! this one without em in D terms
      end do
    end do
    !$acc end parallel loop

  case(2)
    ! jlm deformation scheme using 3D uc, vc, wc and omega (1st rough scheme)
    !$acc parallel loop collapse(2) copyin(hdiff,uc,vc,wc,ps,dpsldt) copyout(t_kh(1:ifull,1:kl)) &
    !$acc   present(ie,iw,in,is)
    do concurrent (k = 1:kl)
      do concurrent (iq = 1:ifull)
        hdif = dt*hdiff(k)/ds ! N.B.  hdiff(k)=khdif*.1
        cc = (uc(ie(iq),k)-uc(iw(iq),k))**2 + (uc(in(iq),k)-uc(is(iq),k))**2 + &
             (vc(ie(iq),k)-vc(iw(iq),k))**2 + (vc(in(iq),k)-vc(is(iq),k))**2 + &
             (wc(ie(iq),k)-wc(iw(iq),k))**2 + (wc(in(iq),k)-wc(is(iq),k))**2 + &
             .01*(dpsldt(ie(iq),k)*ps(ie(iq))-dpsldt(iw(iq),k)*ps(iw(iq)))**2+ &
             .01*(dpsldt(in(iq),k)*ps(in(iq))-dpsldt(is(iq),k)*ps(is(iq)))**2 
        ! approx 1 Pa/s = .1 m/s     
        ! N.B. using double grid length
        t_kh(iq,k)= .5*sqrt(cc)*hdif*ps(iq) ! this one without em in D terms
      end do
    enddo

  case(3)
    ! K-eps model + Smag
    call boundsuv(uav,vav,allvec=.true.)
    !$acc parallel loop collapse(2) copyin(hdiff,uav,vav,em,emi,tke,eps) copyout(t_kh(1:ifull,1:kl)) &
    !$acc   present(ieu,iwu,inu,isu,iev,iwv,inv,isv)
    do concurrent (k = 1:kl)
      do concurrent (iq = 1:ifull)
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
    !$acc end parallel loop      

  case DEFAULT
    write(6,*) "ERROR: Unknown option nhorjlm=",nhorjlm
    call ccmpi_abort(-1)
end select


call bounds(t_kh,nehalf=.true.)
do k = 1,kl
  xfact(1:ifull,k) = (t_kh(ie,k)+t_kh(1:ifull,k))*.5
  yfact(1:ifull,k) = (t_kh(in,k)+t_kh(1:ifull,k))*.5
end do    
if ( nhorx==1 ) then
  do k = 1,kl
    xfact(1:ifull,k) = xfact(1:ifull,k)*tx_fact(1:ifull)
    yfact(1:ifull,k) = yfact(1:ifull,k)*ty_fact(1:ifull)
  end do    
else if ( nhorx>=7 ) then
  do k = 1,kmax
    xfact(1:ifull,k) = xfact(1:ifull,k)*tx_fact(1:ifull)
    yfact(1:ifull,k) = yfact(1:ifull,k)*ty_fact(1:ifull)
  end do
end if
call boundsuv(xfact,yfact,stag=-9) ! MJT - can use stag=-9 option that will
                                   ! only update iwu and isv values


! momentum (disabled by default)
if ( nhorps==0 .or. nhorps==-2 ) then ! for nhorps=-1,-3,-4 don't diffuse u,v
  do k = 1,kl
    do concurrent (iq = 1:ifull)
      base = emi(iq)+xfact(iq,k)+xfact(iwu(iq),k)  &
                    +yfact(iq,k)+yfact(isv(iq),k)  
      ucc = ( emi(iq)*uc(iq,k) +              &
              xfact(iq,k)*uc(ie(iq),k) +      &
              xfact(iwu(iq),k)*uc(iw(iq),k) + &
              yfact(iq,k)*uc(in(iq),k) +      &
              yfact(isv(iq),k)*uc(is(iq),k) ) &
            / base
      vcc = ( emi(iq)*vc(iq,k) +              &
              xfact(iq,k)*vc(ie(iq),k) +      &
              xfact(iwu(iq),k)*vc(iw(iq),k) + &
              yfact(iq,k)*vc(in(iq),k) +      &
              yfact(isv(iq),k)*vc(iq,k) )     &
            / base
      wcc = ( emi(iq)*wc(iq,k) +              &
              xfact(iq,k)*wc(ie(iq),k) +      &
              xfact(iwu(iq),k)*wc(iw(iq),k) + &
              yfact(iq,k)*wc(in(iq),k) +      &
              yfact(isv(iq),k)*wc(is(iq),k) ) &
            / base
      u(iq,k) = ax(iq)*ucc + ay(iq)*vcc + az(iq)*wcc
      v(iq,k) = bx(iq)*ucc + by(iq)*vcc + bz(iq)*wcc
    end do
  end do
end if   ! nhorps==0 .or. nhorps==-2

if ( diag .and. mydiag ) then
  do k = 1,kl
    write(6,*) 'k,id,jd,idjd ',k,id,jd,idjd
    write(6,*) 'k, xfact, xfactw ',k,xfact(idjd,k),xfact(iwu(idjd),k)
    write(6,*) 'k, yfact, yfacts ',k,yfact(idjd,k),yfact(isv(idjd),k)
    write(6,*) 'k, uc,uce,ucw,ucn,ucs ',k,uc(idjd,k),uc(ie(idjd),k),uc(iw(idjd),k),uc(in(idjd),k),uc(is(idjd),k)
    write(6,*) 'k,u,v ',k,u(idjd,k),v(idjd,k)
  end do
endif

! do t diffusion based on potential temperature ff
! for nhorps=-3 don't diffuse T or cloud; only qg

! t + qg enabled by default
if ( nhorps==0 .or. nhorps==-1 .or. nhorps==-4 .or. nhorps==-6 ) then  
  do k = 1,kl
    work(1:ifull,k,1) = t(1:ifull,k)/ptemp(1:ifull) ! watch out for Chen!
    work(1:ifull,k,2) = qg(1:ifull,k)
  end do
  call bounds(work(:,:,1:2))
  !$acc parallel loop collapse(2) copyin(ptemp,emi,xfact,yfact,work(:,:,1:2)) copyout(t,qg)
  do k = 1,kl
    do concurrent (iq = 1:ifull)
      base = emi(iq)+xfact(iq,k)+xfact(iwu(iq),k)  &
                    +yfact(iq,k)+yfact(isv(iq),k)  
      t(iq,k) = ( emi(iq)*work(iq,k,1) +              &
                  xfact(iq,k)*work(ie(iq),k,1) +      &
                  xfact(iwu(iq),k)*work(iw(iq),k,1) + &
                  yfact(iq,k)*work(in(iq),k,1) +      &
                  yfact(isv(iq),k)*work(is(iq),k,1) ) &
                / base
      t(iq,k) = ptemp(iq)*t(iq,k)
      qg(iq,k) = ( emi(iq)*work(iq,k,2) +               &
                   xfact(iq,k)*work(ie(iq),k,2) +       &
                   xfact(iwu(iq),k)*work(iw(iq),k,2) +  &
                   yfact(iq,k)*work(in(iq),k,2) +       &
                   yfact(isv(iq),k)*work(is(iq),k,2) )  &
                 / base
    end do
  end do
  !$acc end parallel loop

else if ( nhorps==-5 ) then  
  do k = 1,kl
    work(1:ifull,k,1) = t(1:ifull,k)/ptemp(1:ifull) ! watch out for Chen!
  end do
  call bounds(work(:,:,1))
  do k = 1,kl
    do concurrent (iq = 1:ifull)
      base = emi(iq)+xfact(iq,k)+xfact(iwu(iq),k)  &
                    +yfact(iq,k)+yfact(isv(iq),k)  
      t(iq,k) = ( emi(iq)*work(iq,k,1) +              &
                  xfact(iq,k)*work(ie(iq),k,1) +      &
                  xfact(iwu(iq),k)*work(iw(iq),k,1) + &
                  yfact(iq,k)*work(in(iq),k,1) +      &
                  yfact(isv(iq),k)*work(is(iq),k,1) ) &
                / base
      t(iq,k) = ptemp(iq)*t(iq,k)
    end do
  end do

else if ( nhorps==-3 ) then  
  work(1:ifull,1:kl,1) = qg(1:ifull,1:kl)
  call bounds(work(:,:,1))
  do k = 1,kl      
    do concurrent (iq = 1:ifull)
      base = emi(iq)+xfact(iq,k)+xfact(iwu(iq),k)  &
                    +yfact(iq,k)+yfact(isv(iq),k)  
      qg(iq,k) = ( emi(iq)*work(iq,k,1) +               &
                   xfact(iq,k)*work(ie(iq),k,1) +       &
                   xfact(iwu(iq),k)*work(iw(iq),k,1) +  &
                   yfact(iq,k)*work(in(iq),k,1) +       &
                   yfact(isv(iq),k)*work(is(iq),k,1) )  &
                 / base
    end do
  end do
end if

! cloud microphysics (disabled by default) 
if ( nhorps==-4 .and. ldr/=0 ) then  
  work(1:ifull,1:kl,1) = qlg(1:ifull,1:kl)
  work(1:ifull,1:kl,2) = qfg(1:ifull,1:kl)
  work(1:ifull,1:kl,3) = stratcloud(1:ifull,1:kl)
  call bounds(work(:,:,1:3))
  do k = 1,kl
    do concurrent (iq = 1:ifull)
      base = emi(iq)+xfact(iq,k)+xfact(iwu(iq),k)  &
                      +yfact(iq,k)+yfact(isv(iq),k)  
      qlg(iq,k) = ( emi(iq)*work(iq,k,1) +               &
                    xfact(iq,k)*work(ie(iq),k,1) +       &
                    xfact(iwu(iq),k)*work(iw(iq),k,1) +  &
                    yfact(iq,k)*work(in(iq),k,1) +       &
                    yfact(isv(iq),k)*work(is(iq),k,1) )  &
                  / base
      qfg(iq,k) = ( emi(iq)*work(iq,k,2) +               &
                    xfact(iq,k)*work(ie(iq),k,2) +       &
                    xfact(iwu(iq),k)*work(iw(iq),k,2) +  &
                    yfact(iq,k)*work(in(iq),k,2) +       &
                    yfact(isv(iq),k)*work(is(iq),k,2) )  &
                  / base
      stratcloud(iq,k) = ( emi(iq)*work(iq,k,3) +               & 
                           xfact(iq,k)*work(ie(iq),k,3) +       &
                           xfact(iwu(iq),k)*work(iw(iq),k,3) +  &
                           yfact(iq,k)*work(in(iq),k,3) +       &
                           yfact(isv(iq),k)*work(is(iq),k,3) )  &
                         / base
    end do
  end do
end if       ! (ldr/=0.and.nhorps==-4)

! apply horizontal diffusion to TKE and EPS terms (disabled by default)
if ( (nhorps==0.or.nhorps==-1.or.nhorps==-4) .and. nvmix==6 ) then
  work(1:ifull,1:kl,1) = tke(1:ifull,1:kl)
  work(1:ifull,1:kl,2) = eps(1:ifull,1:kl)
  call bounds(work(:,:,1:2))
  do k = 1,kl
    do concurrent (iq = 1:ifull)
      base = emi(iq)+xfact(iq,k)+xfact(iwu(iq),k)  &
                    +yfact(iq,k)+yfact(isv(iq),k)  
      tke(iq,k) = ( emi(iq)*work(iq,k,1) +              &
                    xfact(iq,k)*work(ie(iq),k,1) +      &
                    xfact(iwu(iq),k)*work(iw(iq),k,1) + &
                    yfact(iq,k)*work(in(iq),k,1) +      &
                    yfact(isv(iq),k)*work(is(iq),k,1) ) &
                  / base
      eps(iq,k) = ( emi(iq)*work(iq,k,2) +               &
                    xfact(iq,k)*work(ie(iq),k,2) +       &
                    xfact(iwu(iq),k)*work(iw(iq),k,2) +  &
                    yfact(iq,k)*work(in(iq),k,2) +       &
                    yfact(isv(iq),k)*work(is(iq),k,2) )  &
                  / base
    end do
  end do
end if   ! (nvmix==6)

! prgnostic aerosols (disabled by default)
if ( nhorps==-4 .and. abs(iaero)>=2 ) then
  do nstart = 1,naero,3
    nend = min( nstart+2, naero)
    nt = nend - nstart + 1
    work(1:ifull,1:kl,1:nt) = xtg(1:ifull,1:kl,nstart:nend)
    call bounds(work(:,:,1:nt))
    do ntr = 1,nt
      do k = 1,kl
        do concurrent (iq = 1:ifull)
          base = emi(iq)+xfact(iq,k)+xfact(iwu(iq),k)  &
                        +yfact(iq,k)+yfact(isv(iq),k)  
          xtg(iq,k,nstart+ntr-1) = ( emi(iq)*work(iq,k,nt) +               &
                                     xfact(iq,k)*work(ie(iq),k,nt) +       &
                                     xfact(iwu(iq),k)*work(iw(iq),k,nt) +  &
                                     yfact(iq,k)*work(in(iq),k,nt) +       &
                                     yfact(isv(iq),k)*work(is(iq),k,nt) )  &
                                   / base
        end do
      end do
    end do
  end do
end if  ! (nhorps==-4.and.abs(iaero)>=2)  

return
end subroutine hordifgt
