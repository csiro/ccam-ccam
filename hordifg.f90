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
    
! This subroutine calculates horizontal diffusion for atmospheric fields
      
subroutine hordifgt   !  globpea version    N.B. k loop in here

use aerosolldr
use arrays_m
use cc_mpi
use cfrac_m
use cloudmod
use indices_m
use liqwpar_m
use map_m
use nharrs_m
use parmhdff_m
use savuvt_m
use sigs_m
use tkeeps, only : tke,eps,shear,mintke,mineps,cm0,minl,maxl
use vecsuv_m
use vvel_m
implicit none
!     called from globpe (now not tendencies),
!     called for -ve nhor
!     for +ve nhor see hordifg 
!     It has -ve nhorps option:
!        nhorps=0  does everything
!        nhorps=-1 does only T & qg        horiz diff.
!        nhorps=-2 does only u &  v        horiz diff.
!        nhorps=-3 does only qg            horiz diff.
!        nhorps=-4 does only T, qg & cloud horiz diff.
!     and u,v have same options as T (e.g.nhor=-157)
!     this one has got map factors
!     has jlm nhorx option as last digit of nhor, e.g. -157
include 'newmpar.h'
include 'const_phys.h'
include 'kuocom.h'
include 'parm.h'
include 'parmdyn.h'

real, dimension(ifull+iextra,kl,nagg) :: ff
real, dimension(ifull+iextra,kl) :: uc, vc, wc
real, dimension(ifull+iextra,kl) :: ww, uav, vav
real, dimension(ifull+iextra,kl) :: xfact, yfact, t_kh
real, dimension(ifull,kl) :: xfact_iwu, yfact_isv
real, dimension(ifull,kl) :: dudx, dudy, dudz
real, dimension(ifull,kl) :: dvdx, dvdy, dvdz
real, dimension(ifull,kl) :: dwdx, dwdy, dwdz
real, dimension(ifull,kl) :: base, emi
real, dimension(ifull,kl) :: zg, tnhs, tv
real, dimension(ifull,2) :: zgh      
real, dimension(ifull) :: ptemp, tx_fact, ty_fact
real, dimension(ifull) :: sx_fact, sy_fact
real, dimension(ifull) :: r1, r2, cc
real, dimension(ifull) :: ucc, vcc, wcc
real delphi, hdif
integer k, nhora, nhorx, ntr
integer, save :: kmax=-1
integer, parameter :: nf=2


!     nhorx used in hordif  ! previous code effectively has nhorx=0
!           = 1 u, v, T, qg  diffusion reduced near mountains
!           = 4              diffusion not reduced near mountains    
!           = 7 u, v, T, qg  diffusion reduced near mountains (bottom 2/3 only)

!     khor (set up in darlam.f): if non-zero increases the diffusion in the
!     upper model levels   (never used these days)
!         khor > 0 progessively doubles for khor upper levels
!         khor < 0 linearly increases above sig = 0.2 (up to factor of 1-khor)

!     this horizontal diffusion routine allows compensation for
!     being on sigma surfaces through reductions of kh over orography
!     this is called for -nhor.ge.50; the value of nhor (in metres) gives the
!     scaling value for deltaz
!     in namelist khdif is fujio's a**2, e.g. 4.
!     if(nhorps.gt.1)stop 'nhorps > 1 not permitted in hordifgt'

nhorx=0

if (kmax<0) then
  kmax=1
  do while (sig(kmax)>=0.25.and.kmax<kl)
    kmax=kmax+1
  end do
end if

!     set up topography reduction factors for each type of location
!     expect power nf to be about 1 or 2 (see data statement)

delphi=1.e6  ! turns off reduction (can also use nhorx=4)
if (abs(nhor)>=50) then
  nhora=10*(abs(nhor)/10)    ! e.g. 150  for nhor=-157
  nhorx=abs(nhor)-nhora      ! e.g.   7  for nhor=-157
  delphi=nhora*grav
endif

do k=1,kl
  ! note the following line is the same on all levels
  ! but helps with the matrix multiplication below
  emi(1:ifull,k)=ps(1:ifull)/em(1:ifull)
end do
ptemp(1:ifull)=ps(1:ifull)**.286
tx_fact(1:ifull)=1./(1.+(abs(zs(ie)-zs(1:ifull))/delphi)**nf)
ty_fact(1:ifull)=1./(1.+(abs(zs(in)-zs(1:ifull))/delphi)**nf)
! sx_fact and sy_fact are unstaggered versions of tx_fact and ty_fact
sx_fact(1:ifull)=1./(1.+(0.5*abs(zs(ie)-zs(iw))/delphi)**nf)
sy_fact(1:ifull)=1./(1.+(0.5*abs(zs(in)-zs(is))/delphi)**nf)

#ifdef debug
if(diag.and.mydiag)then
  write(6,*) 'hordifgt u ',(u(idjd,k),k=1,kl)
  write(6,*) 'hordifgt v ',(v(idjd,k),k=1,kl)
  write(6,*) 'ax,ay,az ',ax(idjd),ay(idjd),az(idjd)
  write(6,*) 'bx,by,bz ',bx(idjd),by(idjd),bz(idjd)
endif
#endif

! This option is for Smag diffusion or prognostic TKE
if ( nhorjlm==0 .or. nhorjlm==3 .or. nvmix==6 ) then
  ! Calculate du/dx,dv/dx,du/dy,dv/dy, etc 

  ! calculate height on full levels and non-hydrostatic temp correction
  tv(:,:) = t(1:ifull,:)*(1.+0.61*qg(1:ifull,:)-qlg(1:ifull,:)-qfg(1:ifull,:) &
                         -qrg(1:ifull,:)-qsng(1:ifull,:)-qgrg(1:ifull,:))
  tnhs(:,1)=phi_nh(:,1)/bet(1)
  zg(:,1) = (zs(1:ifull)+bet(1)*tv(:,1))/grav
  do k=2,kl
    tnhs(:,k)=(phi_nh(:,k)-phi_nh(:,k-1)-betm(k)*tnhs(:,k-1))/bet(k)
    zg(:,k) = zg(:,k-1) + (bet(k)*tv(:,k)+betm(k)*tv(:,k-1))/grav
  end do ! k  loop
  zg(:,:) = zg(:,:) + phi_nh(:,:)/grav

  do k=1,kl        
    ! weighted horizontal velocities
    uav(1:ifull,k)=av_vmod*u(1:ifull,k)+(1.-av_vmod)*savu(1:ifull,k)
    vav(1:ifull,k)=av_vmod*v(1:ifull,k)+(1.-av_vmod)*savv(1:ifull,k)
        
    ! calculate vertical velocity in m/s
    ! omega=ps*dpsldt
    ! ww = -R/g * (T+Tnhs) * dpsldt/sig
    ww(1:ifull,k)=(dpsldt(:,k)/sig(k))*(-rdry/grav)*(tv(:,k)+tnhs(:,k))
  end do
        
  call boundsuv_allvec(uav,vav)
  call bounds(ww)

  do k=1,kl
    dudx(:,k)=0.5*(uav(ieu,k)-uav(iwu,k))*em(1:ifull)/ds
    dudy(:,k)=0.5*(uav(inu,k)-uav(isu,k))*em(1:ifull)/ds
    dvdx(:,k)=0.5*(vav(iev,k)-vav(iwv,k))*em(1:ifull)/ds
    dvdy(:,k)=0.5*(vav(inv,k)-vav(isv,k))*em(1:ifull)/ds
    dwdx(:,k)=0.5*(ww(ie,k)-ww(iw,k))*em(1:ifull)/ds
    dwdy(:,k)=0.5*(ww(in,k)-ww(is,k))*em(1:ifull)/ds
  end do
        
  ! calculate vertical gradients
  zgh(:,2)=ratha(1)*zg(:,2)+rathb(1)*zg(:,1) ! upper half level
  r1=uav(1:ifull,1)
  r2=ratha(1)*uav(1:ifull,2)+rathb(1)*uav(1:ifull,1)          
  dudz(1:ifull,1)=(r2-r1)/(zgh(1:ifull,2)-zg(1:ifull,1))
  r1=vav(1:ifull,1)
  r2=ratha(1)*vav(1:ifull,2)+rathb(1)*vav(1:ifull,1)
  dvdz(1:ifull,1)=(r2-r1)/(zgh(1:ifull,2)-zg(1:ifull,1))
  r1=ww(1:ifull,1)
  r2=ratha(1)*ww(1:ifull,2)+rathb(1)*ww(1:ifull,1)          
  dwdz(1:ifull,1)=(r2-r1)/(zgh(1:ifull,2)-zg(1:ifull,1))
  do k=2,kl-1
    zgh(:,1)=zgh(:,2) ! lower half level
    zgh(:,2)=ratha(k)*zg(:,k+1)+rathb(k)*zg(:,k) ! upper half level
    r1=ratha(k-1)*uav(1:ifull,k)+rathb(k-1)*uav(1:ifull,k-1)
    r2=ratha(k)*uav(1:ifull,k+1)+rathb(k)*uav(1:ifull,k)          
    dudz(1:ifull,k)=(r2-r1)/(zgh(1:ifull,2)-zgh(1:ifull,1))
    r1=ratha(k-1)*vav(1:ifull,k)+rathb(k-1)*vav(1:ifull,k-1)
    r2=ratha(k)*vav(1:ifull,k+1)+rathb(k)*vav(1:ifull,k)          
    dvdz(1:ifull,k)=(r2-r1)/(zgh(1:ifull,2)-zgh(1:ifull,1))
    r1=ratha(k-1)*ww(1:ifull,k)+rathb(k-1)*ww(1:ifull,k-1)
    r2=ratha(k)*ww(1:ifull,k+1)+rathb(k)*ww(1:ifull,k)          
    dwdz(1:ifull,k)=(r2-r1)/(zgh(1:ifull,2)-zgh(1:ifull,1))
  end do
  zgh(:,1)=zgh(:,2) ! lower half level
  r1=ratha(kl-1)*uav(1:ifull,kl)+rathb(kl-1)*uav(1:ifull,kl-1)
  r2=uav(1:ifull,kl)          
  dudz(1:ifull,kl)=(r2-r1)/(zg(1:ifull,kl)-zgh(1:ifull,1))
  r1=ratha(kl-1)*vav(1:ifull,kl)+rathb(kl-1)*vav(1:ifull,kl-1)
  r2=vav(1:ifull,kl)          
  dvdz(1:ifull,kl)=(r2-r1)/(zg(1:ifull,kl)-zgh(1:ifull,1))
  r1=ratha(kl-1)*ww(1:ifull,kl)+rathb(kl-1)*ww(1:ifull,kl-1)
  r2=ww(1:ifull,kl)          
  dwdz(1:ifull,kl)=(r2-r1)/(zg(1:ifull,kl)-zgh(1:ifull,1))
        
end if   ! nhorjlm==0.or.nvmix==6
      
! usual deformation for nhorjlm=1 or nhorjlm=2
if ( nhorjlm==1 .or. nhorjlm==2 .or. nhorps==0 .or. nhorps==-2 ) then 
  do k=1,kl
    ! in hordifgt, need to calculate Cartesian components 
    ff(1:ifull,k,1) = ax(1:ifull)*u(1:ifull,k) + bx(1:ifull)*v(1:ifull,k)
    ff(1:ifull,k,2) = ay(1:ifull)*u(1:ifull,k) + by(1:ifull)*v(1:ifull,k)
    ff(1:ifull,k,3) = az(1:ifull)*u(1:ifull,k) + bz(1:ifull)*v(1:ifull,k)
  end do
  call bounds(ff(:,:,1:3))
  uc(:,:)=ff(:,:,1)
  vc(:,:)=ff(:,:,2)
  wc(:,:)=ff(:,:,3)
end if

! apply horz diffusion
select case(nhorjlm)
  case(0)
    ! This is based on 2D Smagorinsky closure
    do k=1,kl
      hdif=dt*hdiff(k) ! N.B.  hdiff(k)=khdif*.1
      r1(:)=(dudx(:,k)-dvdy(:,k))**2+(dvdx(:,k)+dudy(:,k))**2
      t_kh(1:ifull,k)=sqrt(r1(:))*hdif*emi(:,k)
    end do
    call bounds(t_kh,nehalf=.true.)
    do k=1,kl
      xfact(1:ifull,k) = (t_kh(ie,k)+t_kh(1:ifull,k))*.5
      yfact(1:ifull,k) = (t_kh(in,k)+t_kh(1:ifull,k))*.5
    end do
             
  case(1)
    ! jlm deformation scheme using 3D uc, vc, wc
    do k=1,kl
      hdif=dt*hdiff(k)/ds ! N.B.  hdiff(k)=khdif*.1
      cc = (uc(ie,k)-uc(iw,k))**2 + (uc(in,k)-uc(is,k))**2 + &
           (vc(ie,k)-vc(iw,k))**2 + (vc(in,k)-vc(is,k))**2 + &
           (wc(ie,k)-wc(iw,k))**2 + (wc(in,k)-wc(is,k))**2
      ! N.B. using double grid length
      t_kh(1:ifull,k)= .5*sqrt(cc)*hdif*ps(1:ifull) ! this one without em in D terms
    enddo
    call bounds(t_kh,nehalf=.true.)
    do k=1,kl
      xfact(1:ifull,k) = (t_kh(ie,k)+t_kh(1:ifull,k))*.5
      yfact(1:ifull,k) = (t_kh(in,k)+t_kh(1:ifull,k))*.5
    end do

  case(2)
    ! jlm deformation scheme using 3D uc, vc, wc and omega (1st rough scheme)
    do k=1,kl
      hdif=dt*hdiff(k)/ds ! N.B.  hdiff(k)=khdif*.1
      cc = (uc(ie,k)-uc(iw,k))**2 + (uc(in,k)-uc(is,k))**2 + &
           (vc(ie,k)-vc(iw,k))**2 + (vc(in,k)-vc(is,k))**2 + &
           (wc(ie,k)-wc(iw,k))**2 + (wc(in,k)-wc(is,k))**2 + &
           .01*(dpsldt(ie,k)*ps(ie)-dpsldt(iw,k)*ps(iw))**2+ &
           .01*(dpsldt(in,k)*ps(in)-dpsldt(is,k)*ps(is))**2 
      ! approx 1 Pa/s = .1 m/s     
      ! N.B. using double grid length
      t_kh(1:ifull,k)= .5*sqrt(cc)*hdif*ps(1:ifull) ! this one without em in D terms
    enddo
    call bounds(t_kh,nehalf=.true.)
    do k=1,kl
      xfact(1:ifull,k) = (t_kh(ie,k)+t_kh(1:ifull,k))*.5
      yfact(1:ifull,k) = (t_kh(in,k)+t_kh(1:ifull,k))*.5
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
      xfact(1:ifull,k) = (t_kh(ie,k)+t_kh(1:ifull,k))*.5
      yfact(1:ifull,k) = (t_kh(in,k)+t_kh(1:ifull,k))*.5
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
      shear(:,k)=2.*(dwdz(:,k)**2                                   &
                +(dudx(:,k)*sx_fact)**2+(dvdy(:,k)*sy_fact)**2)     &
                +(dudy(:,k)*sy_fact+dvdx(:,k)*sx_fact)**2           &
                +(dudz(:,k)+dwdx(:,k)*sx_fact)**2                   &
                +(dvdz(:,k)+dwdy(:,k)*sy_fact)**2
    end do
  else if (nhorx>=7) then
    do k=1,kmax
      shear(:,k)=2.*(dwdz(:,k)**2                                   &
                +(dudx(:,k)*sx_fact)**2+(dvdy(:,k)*sy_fact)**2)     &
                +(dudy(:,k)*sy_fact+dvdx(:,k)*sx_fact)**2           &
                +(dudz(:,k)+dwdx(:,k)*sx_fact)**2                   &
                +(dvdz(:,k)+dwdy(:,k)*sy_fact)**2
    end do
    do k=kmax+1,kl
      shear(:,k)=2.*(dudx(:,k)**2+dvdy(:,k)**2+dwdz(:,k)**2)        &
                +(dudy(:,k)+dvdx(:,k))**2+(dudz(:,k)+dwdx(:,k))**2  &
                +(dvdz(:,k)+dwdy(:,k))**2
    end do
  else
    shear(:,:)=2.*(dudx(:,:)**2+dvdy(:,:)**2+dwdz(:,:)**2)          &
              +(dudy(:,:)+dvdx(:,:))**2+(dudz(:,:)+dwdx(:,:))**2    &
              +(dvdz(:,:)+dwdy(:,:))**2
  end if
end if

if (nhorx==1) then
  do k=1,kl
    xfact(1:ifull,k) = xfact(1:ifull,k)*tx_fact(1:ifull)
    yfact(1:ifull,k) = yfact(1:ifull,k)*ty_fact(1:ifull)
  end do    
else if (nhorx>=7) then
  do k=1,kmax
    xfact(1:ifull,k) = xfact(1:ifull,k)*tx_fact(1:ifull)
    yfact(1:ifull,k) = yfact(1:ifull,k)*ty_fact(1:ifull)
  end do
end if
call boundsuv(xfact,yfact,stag=-9) ! MJT - can use stag=-9 option that will
                                   ! only update iwu and isv values

xfact_iwu(:,:) = xfact(iwu,:)
yfact_isv(:,:) = yfact(isv,:)
base(:,:)=emi(:,:)+xfact(1:ifull,:)+xfact_iwu(:,:)+yfact(1:ifull,:)+yfact_isv(:,:)

if ( nhorps==0 .or. nhorps==-2 ) then ! for nhorps=-1,-3,-4 don't diffuse u,v
  do k=1,kl
    ucc = ( uc(1:ifull,k)*emi(1:ifull,k) +     &
            xfact(1:ifull,k)*uc(ie,k) +        &
            xfact_iwu(1:ifull,k)*uc(iw,k) +    &
            yfact(1:ifull,k)*uc(in,k) +        &
            yfact_isv(1:ifull,k)*uc(is,k) ) /  &
          base(1:ifull,k)
    vcc = ( vc(1:ifull,k)*emi(1:ifull,k) +     &
            xfact(1:ifull,k)*vc(ie,k) +        &
            xfact_iwu(1:ifull,k)*vc(iw,k) +    &
            yfact(1:ifull,k)*vc(in,k) +        &
            yfact_isv(1:ifull,k)*vc(is,k) ) /  &
          base(1:ifull,k)
    wcc = ( wc(1:ifull,k)*emi(1:ifull,k) +     &
            xfact(1:ifull,k)*wc(ie,k) +        &
            xfact_iwu(1:ifull,k)*wc(iw,k) +    &
            yfact(1:ifull,k)*wc(in,k) +        &
            yfact_isv(1:ifull,k)*wc(is,k) ) /  &
          base(1:ifull,k)
    u(1:ifull,k) = ax(1:ifull)*ucc + ay(1:ifull)*vcc + az(1:ifull)*wcc
    v(1:ifull,k) = bx(1:ifull)*ucc + by(1:ifull)*vcc + bz(1:ifull)*wcc
  end do
endif   ! nhorps.ge.0

#ifdef debug
if(diag.and.mydiag)then
  do k=1,kl
    write(6,*) 'k,id,jd,idjd ',k,id,jd,idjd
    write(6,*) 'k, xfact, xfactw ',k,xfact(idjd,k),xfact(iwu(idjd),k)
    write(6,*) 'k, yfact, yfacts ',k,yfact(idjd,k),yfact(isv(idjd),k)
    write(6,*) 'k, uc,uce,ucw,ucn,ucs ',k,uc(idjd,k),uc(ie(idjd),k),uc(iw(idjd),k),uc(in(idjd),k),uc(is(idjd),k)
    write(6,*) 'k,u,v ',k,u(idjd,k),v(idjd,k)
  end do
endif
#endif

if ( nhorps/=-2 ) then   ! for nhorps=-2 don't diffuse T, qg or cloud
  ! do t diffusion based on potential temperature ff
  if(nhorps/=-3)then  ! for nhorps=-3 don't diffuse T or cloud; only qg
    do k=1,kl
      ff(1:ifull,k,1)=t(1:ifull,k)/ptemp(1:ifull) ! watch out for Chen!
    end do
    ff(1:ifull,:,2)=qg(1:ifull,:)
    call bounds(ff(:,:,1:2))
    t(1:ifull,:)= ( ff(1:ifull,:,1)*emi(1:ifull,:) +    &
                   xfact(1:ifull,:)*ff(ie,:,1) +        &
                   xfact_iwu(1:ifull,:)*ff(iw,:,1) +    &
                   yfact(1:ifull,:)*ff(in,:,1) +        &
                   yfact_isv(1:ifull,:)*ff(is,:,1) ) /  &
                 base(1:ifull,:)
    qg(1:ifull,:) = ( ff(1:ifull,:,2)*emi(1:ifull,:) +  &
                   xfact(1:ifull,:)*ff(ie,:,2) +        &
                   xfact_iwu(1:ifull,:)*ff(iw,:,2) +    &
                   yfact(1:ifull,:)*ff(in,:,2) +        &
                   yfact_isv(1:ifull,:)*ff(is,:,2) ) /  &
                 base(1:ifull,:)
    do k=1,kl
      t(1:ifull,k) = ptemp(1:ifull) * t(1:ifull,k) 
    end do
  else
    ff(1:ifull,:,1)=qg(1:ifull,:)
    call bounds(ff(:,:,1:1))
    qg(1:ifull,:) = ( ff(1:ifull,:,1)*emi(1:ifull,:) +  &
                   xfact(1:ifull,:)*ff(ie,:,1) +        &
                   xfact_iwu(1:ifull,:)*ff(iw,:,1) +    &
                   yfact(1:ifull,:)*ff(in,:,1) +        &
                   yfact_isv(1:ifull,:)*ff(is,:,1) ) /  &
                 base(1:ifull,:)
  endif  ! (nhorps/=-3) ..else..
  
  ! cloud microphysics
  if ( ldr/=0 .and. nhorps==-4 ) then
    ff(1:ifull,:,1)=qlg(1:ifull,:)
    ff(1:ifull,:,2)=qfg(1:ifull,:)
    if ( ncloud>=2 ) then
      ff(1:ifull,:,3)=qrg(1:ifull,:)
      ff(1:ifull,:,4)=rfrac(1:ifull,:)    
      if ( ncloud>=3 ) then
        ff(1:ifull,:,5)=qsng(1:ifull,:)
        ff(1:ifull,:,6)=qgrg(1:ifull,:)
        ff(1:ifull,:,7)=sfrac(1:ifull,:)
        ff(1:ifull,:,8)=gfrac(1:ifull,:)
        if ( ncloud>=4 ) then
          ff(1:ifull,:,9)=stratcloud(1:ifull,:)
          call bounds(ff(:,:,1:9))
          stratcloud(1:ifull,:) = ( ff(1:ifull,:,9)*emi(1:ifull,:) +  &
                   xfact(1:ifull,:)*ff(ie,:,9) +                      &
                   xfact_iwu(1:ifull,:)*ff(iw,:,9) +                  &
                   yfact(1:ifull,:)*ff(in,:,9) +                      &
                   yfact_isv(1:ifull,:)*ff(is,:,9) ) /                &
                  base(1:ifull,:)
        else
          ! MJT notes - cfrac is not advected as it will be diagnosed later in leoncld.f90
          call bounds(ff(:,:,1:8))
        end if
        qsng(1:ifull,:) = ( ff(1:ifull,:,5)*emi(1:ifull,:) +   &
               xfact(1:ifull,:)*ff(ie,:,5) +                   &
               xfact_iwu(1:ifull,:)*ff(iw,:,5) +               &
               yfact(1:ifull,:)*ff(in,:,5) +                   &
               yfact_isv(1:ifull,:)*ff(is,:,5) ) /             &
               base(1:ifull,:)
        qgrg(1:ifull,:) = ( ff(1:ifull,:,6)*emi(1:ifull,:) +   &
               xfact(1:ifull,:)*ff(ie,:,6) +                   &
               xfact_iwu(1:ifull,:)*ff(iw,:,6) +               &
               yfact(1:ifull,:)*ff(in,:,6) +                   &
               yfact_isv(1:ifull,:)*ff(is,:,6) ) /             &
               base(1:ifull,:)
        sfrac(1:ifull,:) = ( ff(1:ifull,:,7)*emi(1:ifull,:) +  &
               xfact(1:ifull,:)*ff(ie,:,7) +                   &
               xfact_iwu(1:ifull,:)*ff(iw,:,7) +               &
               yfact(1:ifull,:)*ff(in,:,7) +                   &
               yfact_isv(1:ifull,:)*ff(is,:,7) ) /             &
               base(1:ifull,:)
        gfrac(1:ifull,:) = ( ff(1:ifull,:,8)*emi(1:ifull,:) +  &
               xfact(1:ifull,:)*ff(ie,:,8) +                   &
               xfact_iwu(1:ifull,:)*ff(iw,:,8) +               &
               yfact(1:ifull,:)*ff(in,:,8) +                   &
               yfact_isv(1:ifull,:)*ff(is,:,8) ) /             &
               base(1:ifull,:)          
      else
        call bounds(ff(:,:,1:4))          
      end if
      qrg(1:ifull,:) = ( ff(1:ifull,:,3)*emi(1:ifull,:) +    &
             xfact(1:ifull,:)*ff(ie,:,3) +                   &
             xfact_iwu(1:ifull,:)*ff(iw,:,3) +               &
             yfact(1:ifull,:)*ff(in,:,3) +                   &
             yfact_isv(1:ifull,:)*ff(is,:,3) ) /             &
             base(1:ifull,:)
      rfrac(1:ifull,:) = ( ff(1:ifull,:,4)*emi(1:ifull,:) +  &
             xfact(1:ifull,:)*ff(ie,:,4) +                   &
             xfact_iwu(1:ifull,:)*ff(iw,:,4) +               &
             yfact(1:ifull,:)*ff(in,:,4) +                   &
             yfact_isv(1:ifull,:)*ff(is,:,4) ) /             &
             base(1:ifull,:)
    else
      call bounds(ff(:,:,1:2))    
    end if
    qlg(1:ifull,:) = ( ff(1:ifull,:,1)*emi(1:ifull,:) +    &
           xfact(1:ifull,:)*ff(ie,:,1) +                   &
           xfact_iwu(1:ifull,:)*ff(iw,:,1) +               &
           yfact(1:ifull,:)*ff(in,:,1) +                   &
           yfact_isv(1:ifull,:)*ff(is,:,1) ) /             &
           base(1:ifull,:)
    qfg(1:ifull,:) = ( ff(1:ifull,:,2)*emi(1:ifull,:) +    &
           xfact(1:ifull,:)*ff(ie,:,2) +                   &
           xfact_iwu(1:ifull,:)*ff(iw,:,2) +               &
           yfact(1:ifull,:)*ff(in,:,2) +                   &
           yfact_isv(1:ifull,:)*ff(is,:,2) ) /             &
           base(1:ifull,:)
  end if       ! (ldr/=0.and.nhorps==-4)
end if         ! (nhorps/=-2)

if ( nhorps==-4 ) then
  ! apply horizontal diffusion to TKE and EPS terms
  if ( nvmix==6 ) then
    ff(1:ifull,:,1)=tke(1:ifull,:)
    ff(1:ifull,:,2)=eps(1:ifull,:)
    call bounds(ff(:,:,1:2))
    tke(1:ifull,:) = ( ff(1:ifull,:,1)*emi(1:ifull,:) +   &
                    xfact(1:ifull,:)*ff(ie,:,1) +         &
                    xfact_iwu(1:ifull,:)*ff(iw,:,1) +     &
                    yfact(1:ifull,:)*ff(in,:,1) +         &
                    yfact_isv(1:ifull,:)*ff(is,:,1) ) /   &
                  base(1:ifull,:)
    eps(1:ifull,:) = ( ff(1:ifull,:,2)*emi(1:ifull,:) +   &
                    xfact(1:ifull,:)*ff(ie,:,2) +         &
                    xfact_iwu(1:ifull,:)*ff(iw,:,2) +     &
                    yfact(1:ifull,:)*ff(in,:,2) +         &
                    yfact_isv(1:ifull,:)*ff(is,:,2) ) /   &
                  base(1:ifull,:)
  end if   ! (nvmix==6)
       
  ! prgnostic aerosols
  if ( abs(iaero)==2 ) then
    ff(1:ifull,:,1:naero)=xtg(1:ifull,:,:)
    call bounds(ff(:,:,1:naero))
    do ntr=1,naero
      xtg(1:ifull,:,ntr) =                       &
        ( ff(1:ifull,:,ntr)*emi(1:ifull,:) +     &
        xfact(1:ifull,:)*ff(ie,:,ntr) +          &
        xfact_iwu(1:ifull,:)*ff(iw,:,ntr) +      &
        yfact(1:ifull,:)*ff(in,:,ntr) +          &
        yfact_isv(1:ifull,:)*ff(is,:,ntr) ) /    &
        base(1:ifull,:)
    end do
  end if  ! (abs(iaero)==2)  
end if    ! (nhorps==-4)

return
end
