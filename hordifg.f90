! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2019 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
use cloudmod
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

real, dimension(ifull+iextra,kl,nagg) :: work
real, dimension(ifull+iextra,kl) :: uc, vc, wc
real, dimension(ifull+iextra,kl) :: uav, vav, ww
real, dimension(ifull+iextra,kl) :: xfact, yfact, t_kh
real, dimension(ifull,kl) :: xfact_iwu, yfact_isv
real, dimension(ifull,kl) :: dudx, dudy, dudz
real, dimension(ifull,kl) :: dvdx, dvdy, dvdz
real, dimension(ifull,kl) :: dwdx, dwdy, dwdz
real, dimension(ifull,kl) :: emi, base
real, dimension(ifull,kl) :: zg
real, dimension(ifull) :: zgh_a, zgh_b
real, dimension(ifull) :: ptemp, tx_fact, ty_fact
real, dimension(ifull) :: sx_fact, sy_fact
real, dimension(ifull) :: r1, r2, cc
real, dimension(ifull) :: ucc, vcc, wcc
real, dimension(ifull) :: zs_n, zs_s, zs_e, zs_w
real, dimension(ifull) :: ww_n, ww_s, ww_e, ww_w
real, dimension(ifull) :: uc_n, uc_s, uc_e, uc_w
real, dimension(ifull) :: vc_n, vc_s, vc_e, vc_w
real, dimension(ifull) :: wc_n, wc_s, wc_e, wc_w
real, dimension(ifull) :: t_kh_n, t_kh_e
real, dimension(ifull) :: work_n, work_s, work_e, work_w
real delphi, hdif
integer k, nhora, nhorx, ntr
integer nstart, nend
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
  nhorx = abs(nhor)-nhora      ! e.g.   7  for nhor=-157
  delphi = nhora*grav
endif

t_kh(:,:) = 0.
xfact(:,:) = 0.
yfact(:,:) = 0.
uav(:,:) = 0.
vav(:,:) = 0.
ww(:,:) = 0.
uc(:,:) = 0.
vc(:,:) = 0.
wc(:,:) = 0.
do k = 1,kl
  ! note the following line is the same on all levels
  ! but helps with the matrix multiplication below
  emi(1:ifull,k) = ps(1:ifull)/em(1:ifull)
end do
ptemp(1:ifull) = ps(1:ifull)**.286
call unpack_nsew(zs,zs_n,zs_s,zs_e,zs_w)
tx_fact(1:ifull) = 1./(1.+(abs(zs_e-zs(1:ifull))/delphi)**nf)
ty_fact(1:ifull) = 1./(1.+(abs(zs_n-zs(1:ifull))/delphi)**nf)
! sx_fact and sy_fact are unstaggered versions of tx_fact and ty_fact
sx_fact(1:ifull) = 1./(1.+(0.5*abs(zs_e-zs_w)/delphi)**nf)
sy_fact(1:ifull) = 1./(1.+(0.5*abs(zs_n-zs_s)/delphi)**nf)

if ( diag .and. mydiag ) then
  write(6,*) 'hordifgt u ',(u(idjd,k),k=1,kl)
  write(6,*) 'hordifgt v ',(v(idjd,k),k=1,kl)
  write(6,*) 'ax,ay,az ',ax(idjd),ay(idjd),az(idjd)
  write(6,*) 'bx,by,bz ',bx(idjd),by(idjd),bz(idjd)
  write(6,*) 'zs,zse,zsn ',zs(idjd),zs(ie(idjd)),zs(in(idjd))
  write(6,*) 'tx_,ty_ ',tx_fact(idjd),ty_fact(idjd)
endif

! This option is for Smag diffusion or prognostic TKE
if ( nhorjlm==0 .or. nhorjlm==3 .or. nvmix==6 ) then
  ! Calculate du/dx,dv/dx,du/dy,dv/dy, etc 

  ! calculate height on full levels
  zg(:,1) = (zs(1:ifull)+bet(1)*t(1:ifull,1))/grav
  do k = 2,kl
    zg(:,k) = zg(:,k-1) + (bet(k)*t(1:ifull,k)+betm(k)*t(1:ifull,k-1))/grav
  end do ! k  loop

  do k = 1,kl        
    ! weighted horizontal velocities
    uav(1:ifull,k) = av_vmod*u(1:ifull,k) + (1.-av_vmod)*savu(1:ifull,k)
    vav(1:ifull,k) = av_vmod*v(1:ifull,k) + (1.-av_vmod)*savv(1:ifull,k)
        
    ! calculate vertical velocity in m/s
    ! omega=ps*dpsldt
    ! ww = -R/g * T / (ps*sig) * ( omega - sig*dpsdt )
    ww(1:ifull,k) = (-rdry/grav)*t(1:ifull,k)/(ps(1:ifull)*sig(k))* &
                    ( ps(:ifull)*dpsldt(:,k) - sig(k)*dpsdt(:)/864. ) ! dpsdt is in hPa/day and need Pa/s
  end do
  call boundsuv(uav,vav,allvec=.true.)
  call bounds(ww)
  do k = 1,kl
    call unpack_nveu(uav(:,k),vav(:,k),vc_n,uc_e)
    call unpack_svwu(uav(:,k),vav(:,k),vc_s,uc_w)
    call unpack_nsew(ww(:,k),ww_n,ww_s,ww_e,ww_w)  
    dudx(:,k) = 0.5*(uc_e-uc_w)*em(1:ifull)/ds
    dudy(:,k) = 0.5*(uav(inu,k)-uav(isu,k))*em(1:ifull)/ds
    dvdx(:,k) = 0.5*(vav(iev,k)-vav(iwv,k))*em(1:ifull)/ds
    dvdy(:,k) = 0.5*(vc_n-vc_s)*em(1:ifull)/ds
    dwdx(:,k) = 0.5*(ww_e-ww_w)*em(1:ifull)/ds
    dwdy(:,k) = 0.5*(ww_n-ww_s)*em(1:ifull)/ds
  end do
        
  ! calculate vertical gradients
  zgh_b = ratha(1)*zg(:,2) + rathb(1)*zg(:,1) ! upper half level
  r1 = uav(1:ifull,1)
  r2 = ratha(1)*uav(1:ifull,2) + rathb(1)*uav(1:ifull,1)          
  dudz(1:ifull,1)=(r2-r1)/(zgh_b-zg(1:ifull,1))
  r1 = vav(1:ifull,1)
  r2 = ratha(1)*vav(1:ifull,2) + rathb(1)*vav(1:ifull,1)          
  dvdz(1:ifull,1) = (r2-r1)/(zgh_b-zg(1:ifull,1))
  r1 = ww(1:ifull,1)
  r2 = ratha(1)*ww(1:ifull,2) + rathb(1)*ww(1:ifull,1)          
  dwdz(1:ifull,1) = (r2-r1)/(zgh_b-zg(1:ifull,1))
  do k = 2,kl-1
    zgh_a = zgh_b ! lower half level
    zgh_b = ratha(k)*zg(:,k+1) + rathb(k)*zg(:,k) ! upper half level
    r1 = ratha(k-1)*uav(1:ifull,k) + rathb(k-1)*uav(1:ifull,k-1)
    r2 = ratha(k)*uav(1:ifull,k+1) + rathb(k)*uav(1:ifull,k)          
    dudz(1:ifull,k) = (r2-r1)/(zgh_b-zgh_a)
    r1 = ratha(k-1)*vav(1:ifull,k) + rathb(k-1)*vav(1:ifull,k-1)
    r2 = ratha(k)*vav(1:ifull,k+1) + rathb(k)*vav(1:ifull,k)          
    dvdz(1:ifull,k) = (r2-r1)/(zgh_b-zgh_a)
    r1 = ratha(k-1)*ww(1:ifull,k) + rathb(k-1)*ww(1:ifull,k-1)
    r2 = ratha(k)*ww(1:ifull,k+1) + rathb(k)*ww(1:ifull,k)          
    dwdz(1:ifull,k) = (r2-r1)/(zgh_b-zgh_a)
  end do
  zgh_a = zgh_b ! lower half level
  r1 = ratha(kl-1)*uav(1:ifull,kl) + rathb(kl-1)*uav(1:ifull,kl-1)
  r2 = uav(1:ifull,kl)          
  dudz(1:ifull,kl) = (r2-r1)/(zg(1:ifull,kl)-zgh_a)
  r1 = ratha(kl-1)*vav(1:ifull,kl) + rathb(kl-1)*vav(1:ifull,kl-1)
  r2 = vav(1:ifull,kl)          
  dvdz(1:ifull,kl) = (r2-r1)/(zg(1:ifull,kl)-zgh_a)
  r1 = ratha(kl-1)*ww(1:ifull,kl) + rathb(kl-1)*ww(1:ifull,kl-1)
  r2 = ww(1:ifull,kl)          
  dwdz(1:ifull,kl) = (r2-r1)/(zg(1:ifull,kl)-zgh_a)
        
end if   ! nhorjlm==0.or.nvmix==6
      
! usual deformation for nhorjlm=1 or nhorjlm=2
if ( nhorjlm==1 .or. nhorjlm==2 .or. nhorps==0 .or. nhorps==-2 ) then 
  do k = 1,kl
    ! in hordifgt, need to calculate Cartesian components 
    work(1:ifull,k,1) = ax(1:ifull)*u(1:ifull,k) + bx(1:ifull)*v(1:ifull,k)
    work(1:ifull,k,2) = ay(1:ifull)*u(1:ifull,k) + by(1:ifull)*v(1:ifull,k)
    work(1:ifull,k,3) = az(1:ifull)*u(1:ifull,k) + bz(1:ifull)*v(1:ifull,k)
  end do
  call bounds(work(:,:,1:3))
  uc = work(:,:,1)
  vc = work(:,:,2)
  wc = work(:,:,3)
end if

! apply horz diffusion
select case(nhorjlm)
  case(0)
    ! This is based on 2D Smagorinsky closure
    do k = 1,kl
      hdif=dt*hdiff(k) ! N.B.  hdiff(k)=khdif*.1
      r1(:)=(dudx(:,k)-dvdy(:,k))**2+(dvdx(:,k)+dudy(:,k))**2
      t_kh(1:ifull,k)=sqrt(r1(:))*hdif*emi(:,k)
    end do
    call bounds(t_kh,nehalf=.true.)
    do k = 1,kl
      call unpack_ne(t_kh(:,k),t_kh_n,t_kh_e)  
      xfact(1:ifull,k) = (t_kh_e+t_kh(1:ifull,k))*.5
      yfact(1:ifull,k) = (t_kh_n+t_kh(1:ifull,k))*.5
    end do
             
  case(1)
    ! jlm deformation scheme using 3D uc, vc, wc
    do k = 1,kl
      hdif = dt*hdiff(k)/ds ! N.B.  hdiff(k)=khdif*.1
      call unpack_nsew(uc(:,k),uc_n,uc_s,uc_e,uc_w)
      call unpack_nsew(vc(:,k),vc_n,vc_s,vc_e,vc_w)
      call unpack_nsew(wc(:,k),wc_n,wc_s,wc_e,wc_w)
      cc = (uc_e-uc_w)**2 + (uc_n-uc_s)**2 + &
           (vc_e-vc_w)**2 + (vc_n-vc_s)**2 + &
           (wc_e-wc_w)**2 + (wc_n-wc_s)**2
      ! N.B. using double grid length
      t_kh(1:ifull,k)= .5*sqrt(cc)*hdif*ps(1:ifull) ! this one without em in D terms
    end do
    call bounds(t_kh,nehalf=.true.)
    do k = 1,kl
      call unpack_ne(t_kh(:,k),t_kh_n,t_kh_e)  
      xfact(1:ifull,k) = (t_kh_e+t_kh(1:ifull,k))*.5
      yfact(1:ifull,k) = (t_kh_n+t_kh(1:ifull,k))*.5
    end do

  case(2)
    ! jlm deformation scheme using 3D uc, vc, wc and omega (1st rough scheme)
    do k=1,kl
      hdif=dt*hdiff(k)/ds ! N.B.  hdiff(k)=khdif*.1
      call unpack_nsew(uc(:,k),uc_n,uc_s,uc_e,uc_w)
      call unpack_nsew(vc(:,k),vc_n,vc_s,vc_e,vc_w)
      call unpack_nsew(wc(:,k),wc_n,wc_s,wc_e,wc_w)
      cc = (uc_e-uc_w)**2 + (uc_n-uc_s)**2 + &
           (vc_e-vc_w)**2 + (vc_n-vc_s)**2 + &
           (wc_e-wc_w)**2 + (wc_n-wc_s)**2 + &
           .01*(dpsldt(ie,k)*ps(ie)-dpsldt(iw,k)*ps(iw))**2+ &
           .01*(dpsldt(in,k)*ps(in)-dpsldt(is,k)*ps(is))**2 
      ! approx 1 Pa/s = .1 m/s     
      ! N.B. using double grid length
      t_kh(1:ifull,k)= .5*sqrt(cc)*hdif*ps(1:ifull) ! this one without em in D terms
    enddo
    call bounds(t_kh,nehalf=.true.)
    do k=1,kl
      call unpack_ne(t_kh(:,k),t_kh_n,t_kh_e)  
      xfact(1:ifull,k) = (t_kh_e+t_kh(1:ifull,k))*.5
      yfact(1:ifull,k) = (t_kh_n+t_kh(1:ifull,k))*.5
    end do

  case(3)
    ! Combine 2D Smagorinsky closure with K-eps model
    do k=1,kl
      hdif=dt*hdiff(k) ! N.B.  hdiff(k)=khdif*.1
      r1(:)=(dudx(:,k)-dvdy(:,k))**2+(dvdx(:,k)+dudy(:,k))**2
      t_kh(1:ifull,k)=sqrt(r1(:))*hdif*emi(:,k)
    end do      
    if (nvmix==6) then
      hdif=dt*cm0
      do k=1,kl
        tke(1:ifull,k)=max( tke(1:ifull,k), mintke )
        r1(:)=(cm0**0.75)*tke(1:ifull,k)*sqrt(tke(1:ifull,k))
        eps(1:ifull,k)=min( eps(1:ifull,k), r1(:)/minl )
        eps(1:ifull,k)=max( eps(1:ifull,k), r1(:)/maxl, mineps )
        t_kh(1:ifull,k)=max( t_kh(1:ifull,k),                                                    &
                             tke(1:ifull,k)*tke(1:ifull,k)/eps(1:ifull,k)*hdif*emi(1:ifull,k) )
      end do
    end if
    call bounds(t_kh,nehalf=.true.)
    do k=1,kl
      call unpack_ne(t_kh(:,k),t_kh_n,t_kh_e)  
      xfact(1:ifull,k) = (t_kh_e+t_kh(1:ifull,k))*.5
      yfact(1:ifull,k) = (t_kh_n+t_kh(1:ifull,k))*.5
    end do    

  case DEFAULT
    write(6,*) "ERROR: Unknown option nhorjlm=",nhorjlm
    call ccmpi_abort(-1)
end select
       
! Calculate horizontal diffusion based on prognostic TKE
! This can be combined with the diffusion coefficents above
! so as to operate over a large range of grid leng th scales
if (nvmix==6) then
  if (nhorx==1) then
    do k=1,kl
      shear(:,k) = 2.*(dwdz(:,k)**2                               &
                 + (dudx(:,k)*sx_fact)**2+(dvdy(:,k)*sy_fact)**2) &
                 + (dudy(:,k)*sy_fact+dvdx(:,k)*sx_fact)**2       &
                 + (dudz(:,k)+dwdx(:,k)*sx_fact)**2               &
                 + (dvdz(:,k)+dwdy(:,k)*sy_fact)**2
    end do
  else if (nhorx>=7) then
    do k=1,kmax
      shear(:,k) = 2.*(dwdz(:,k)**2                               &
                 + (dudx(:,k)*sx_fact)**2+(dvdy(:,k)*sy_fact)**2) &
                 + (dudy(:,k)*sy_fact+dvdx(:,k)*sx_fact)**2       &
                 + (dudz(:,k)+dwdx(:,k)*sx_fact)**2               &
                 + (dvdz(:,k)+dwdy(:,k)*sy_fact)**2
    end do
    do k=kmax+1,kl
      shear(:,k) = 2.*(dwdz(:,k)**2+dudx(:,k)**2+dvdy(:,k)**2)    &
                 + (dudy(:,k)+dvdx(:,k))**2                       &
                 + (dudz(:,k)+dwdx(:,k))**2                       &
                 + (dvdz(:,k)+dwdy(:,k))**2
    end do
  else
    do k = 1,kl
      shear(:,k) = 2.*(dwdz(:,k)**2+dudx(:,k)**2+dvdy(:,k)**2)    &
                 + (dudy(:,k)+dvdx(:,k))**2                       &
                 + (dudz(:,k)+dwdx(:,k))**2                       &
                 + (dvdz(:,k)+dwdy(:,k))**2
    end do
  end if
end if

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

do k = 1,kl
  call unpack_svwu(xfact(:,k),yfact(:,k),yfact_isv(:,k),xfact_iwu(:,k))  
  base(1:ifull,k) = emi(1:ifull,k) + xfact(1:ifull,k) + xfact_iwu(1:ifull,k) &
                                   + yfact(1:ifull,k) + yfact_isv(1:ifull,k)
  !emi(1:ifull,k) = emi(1:ifull,k)/base(1:ifull)
  !xfact(1:ifull,k) = xfact(1:ifull,k)/base(1:ifull)
  !yfact(1:ifull,k) = yfact(1:ifull,k)/base(1:ifull)
  !xfact_iwu(1:ifull,k) = xfact_iwu(1:ifull,k)/base(1:ifull)
  !yfact_isv(1:ifull,k) = yfact_isv(1:ifull,k)/base(1:ifull)
end do

if ( nhorps==0 .or. nhorps==-2 ) then ! for nhorps=-1,-3,-4 don't diffuse u,v
  do k = 1,kl
    call unpack_nsew(uc(:,k),uc_n,uc_s,uc_e,uc_w)  
    ucc = ( emi(1:ifull,k)*uc(1:ifull,k) +     &
            xfact(1:ifull,k)*uc_e +            &
            xfact_iwu(1:ifull,k)*uc_w +        &
            yfact(1:ifull,k)*uc_n +            &
            yfact_isv(1:ifull,k)*uc_s ) / base(1:ifull,k)
    call unpack_nsew(vc(:,k),vc_n,vc_s,vc_e,vc_w)
    vcc = ( emi(1:ifull,k)*vc(1:ifull,k) +     &
            xfact(1:ifull,k)*vc_e +            &
            xfact_iwu(1:ifull,k)*vc_w +        &
            yfact(1:ifull,k)*vc_n +            &
            yfact_isv(1:ifull,k)*vc_s ) / base(1:ifull,k)
    call unpack_nsew(wc(:,k),wc_n,wc_s,wc_e,wc_w)
    wcc = ( emi(1:ifull,k)*wc(1:ifull,k) +     &
            xfact(1:ifull,k)*wc_e +            &
            xfact_iwu(1:ifull,k)*wc_w +        &
            yfact(1:ifull,k)*wc_n +            &
            yfact_isv(1:ifull,k)*wc_s ) / base(1:ifull,k)
    u(1:ifull,k) = ax(1:ifull)*ucc + ay(1:ifull)*vcc + az(1:ifull)*wcc
    v(1:ifull,k) = bx(1:ifull)*ucc + by(1:ifull)*vcc + bz(1:ifull)*wcc
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
if ( nhorps==0 .or. nhorps==-1 .or. nhorps==-4 .or. nhorps==-5 .or. nhorps==-6 ) then  
  do k = 1,kl
    t(1:ifull,k) = t(1:ifull,k)/ptemp(1:ifull) ! watch out for Chen!
  end do
  work(1:ifull,1:kl,1) = t(1:ifull,1:kl)
  call bounds(work(:,:,1))
  do k = 1,kl
    call unpack_nsew(work(:,k,1),work_n,work_s,work_e,work_w)  
    t(1:ifull,k) = ( emi(1:ifull,k)*work(1:ifull,k,1) +   &
                     xfact(1:ifull,k)*work_e +            &
                     xfact_iwu(1:ifull,k)*work_w +        &
                     yfact(1:ifull,k)*work_n +            &
                     yfact_isv(1:ifull,k)*work_s ) / base(1:ifull,k)
    t(1:ifull,k) = ptemp(1:ifull)*t(1:ifull,k)
  end do
end if  
  
if ( nhorps==0 .or. nhorps==-1 .or. nhorps==-3 .or. nhorps==-4 .or. nhorps==-6 ) then  
  work(1:ifull,1:kl,1) = qg(1:ifull,1:kl)
  call bounds(work(:,:,1))
  do k = 1,kl      
    call unpack_nsew(work(:,k,1),work_n,work_s,work_e,work_w)
    qg(1:ifull,k) = ( emi(1:ifull,k)*work(1:ifull,k,1) +  &
                      xfact(1:ifull,k)*work_e +           &
                      xfact_iwu(1:ifull,k)*work_w +       &
                      yfact(1:ifull,k)*work_n +           &
                      yfact_isv(1:ifull,k)*work_s ) / base(1:ifull,k)
  end do
end if
  
if ( nhorps==-4 .and. ldr/=0 ) then  
  ! cloud microphysics
  work(1:ifull,1:kl,1) = qlg(1:ifull,1:kl)
  work(1:ifull,1:kl,2) = qfg(1:ifull,1:kl)
  work(1:ifull,1:kl,3) = stratcloud(1:ifull,1:kl)
  call bounds(work(:,:,1:3))
  do k = 1,kl
    call unpack_nsew(work(:,k,1),work_n,work_s,work_e,work_w)  
    qlg(1:ifull,k) = ( emi(1:ifull,k)*work(1:ifull,k,2) +  &
                       xfact(1:ifull,k)*work_e +           &
                       xfact_iwu(1:ifull,k)*work_w +       &
                       yfact(1:ifull,k)*work_n +           &
                       yfact_isv(1:ifull,k)*work_s ) / base(1:ifull,k)
  end do
  do k = 1,kl
    call unpack_nsew(work(:,k,2),work_n,work_s,work_e,work_w)
    qfg(1:ifull,k) = ( emi(1:ifull,k)*work(1:ifull,k,2) + &
                       xfact(1:ifull,k)*work_e +          &
                       xfact_iwu(1:ifull,k)*work_w +      &
                       yfact(1:ifull,k)*work_n +          &
                       yfact_isv(1:ifull,k)*work_s ) / base(1:ifull,k)
  end do
  do k = 1,kl
    call unpack_nsew(work(:,k,1),work_n,work_s,work_e,work_w)  
    stratcloud(1:ifull,k) = ( emi(1:ifull,k)*work(1:ifull,k,2) + &
                        xfact(1:ifull,k)*work_e +                &
                        xfact_iwu(1:ifull,k)*work_w +            &
                        yfact(1:ifull,k)*work_n +                &
                        yfact_isv(1:ifull,k)*work_s ) / base(1:ifull,k)
  end do
end if       ! (ldr/=0.and.nhorps==-4)

! apply horizontal diffusion to TKE and EPS terms
if ( (nhorps==0.or.nhorps==-1.or.nhorps==-4) .and. nvmix==6 ) then
  work(1:ifull,1:kl,1) = tke(1:ifull,1:kl)
  work(1:ifull,1:kl,2) = eps(1:ifull,1:kl)
  call bounds(work(:,:,1:2))
  do k = 1,kl
    call unpack_nsew(work(:,k,1),work_n,work_s,work_e,work_w) 
    tke(1:ifull,k) = ( emi(1:ifull,k)*work(1:ifull,k,1) +    &
                              xfact(1:ifull,k)*work_e +      &
                              xfact_iwu(1:ifull,k)*work_w +  &
                              yfact(1:ifull,k)*work_n +      &
                              yfact_isv(1:ifull,k)*work_s ) / base(1:ifull,k)
  end do
  do k = 1,kl
    call unpack_nsew(work(:,k,2),work_n,work_s,work_e,work_w)
    eps(1:ifull,k) = ( emi(1:ifull,k)*work(1:ifull,k,2) +    &
                              xfact(1:ifull,k)*work_e +      &
                              xfact_iwu(1:ifull,k)*work_w +  &
                              yfact(1:ifull,k)*work_n +      &
                              yfact_isv(1:ifull,k)*work_s ) / base(1:ifull,k)
  end do
end if   ! (nvmix==6)

! prgnostic aerosols
if ( nhorps==-4 .and. abs(iaero)>=2 ) then
  do nstart = 1,naero,nagg
    nend = min(nstart+nagg-1, naero)
    work(1:ifull,1:kl,1:nagg) = xtg(1:ifull,1:kl,nstart:nend)
    call bounds(work(:,:,1:nagg))
    do ntr = 1,nagg
      do k = 1,kl  
        call unpack_nsew(work(:,k,ntr),work_n,work_s,work_e,work_w)  
        xtg(1:ifull,k,nstart+ntr-1) = ( emi(1:ifull,k)*work(1:ifull,k,ntr) +    &
                               xfact(1:ifull,k)*work_e +                        &
                               xfact_iwu(1:ifull,k)*work_w +                    &
                               yfact(1:ifull,k)*work_n +                        &
                               yfact_isv(1:ifull,k)*work_s ) / base(1:ifull,k)
      end do
    end do
  end do
end if  ! (nhorps==-4.and.abs(iaero)>=2)  

return
end
