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
    
subroutine upglobal      ! globpea version   use ritchie 103

use aerosolldr
use arrays_m
use cc_mpi
use cfrac_m
use cloudmod
use diag_m
use epst_m
use indices_m
use liqwpar_m  ! ifullw
use map_m
use nharrs_m
use nlin_m
use sbar_m
use sigs_m
use staguvmod
use tkeeps, only : tke,eps
use tracers_m
use unn_m
use vadv
use vecsuv_m
use vvel_m     ! sdot
use work3f_m
use xarrs_m
use xyzinfo_m

implicit none

include 'newmpar.h'
include 'const_phys.h'
include 'kuocom.h'   ! ldr
include 'parm.h'
include 'parmdyn.h'  
include 'parmhor.h'  ! mhint, m_bs, nt_adv

integer, parameter :: ntest=0       ! ~8+ for diagnostic stability tests
integer ii,intsch, iq, jj,k, kk, ntr, ierr
integer l, idjdd, nstart, nend, ntot
integer, save :: num_hight = 0, numunstab = 0
integer, dimension(ifull) :: nits, nvadh_pass
real, save, allocatable, dimension(:,:):: tnsav,unsav,vnsav ! for npex=-1
real, dimension(ifull+iextra,kl,10) :: duma
real, dimension(ifull+iextra,kl) :: uc, vc, wc, dd
real, dimension(ifull+iextra) :: aa
real, dimension(ifull,kl) :: theta
real, dimension(ifull) :: sdmx, tempry
real, dimension(ifull) :: denb, vecdot, vdot1, vdot2
real, dimension(ifull) :: vec1x, vec1y, vec1z
real, dimension(ifull) :: vec2x, vec2y, vec2z
real, dimension(ifull) :: vec3x, vec3y, vec3z
real, dimension(kl) :: factr
real(kind=8), dimension(ifull,kl) :: x3d,y3d,z3d

call START_LOG(upglobal_begin)

intsch=mod(ktau,2)

do k=1,kl
  ! finish off RHS terms; this coriolis term was once in nonlin
  ux(1:ifull,k)=ux(1:ifull,k)+.5*dt*(1.-epsf)*f(1:ifull)*v(1:ifull,k) ! end of Eq. 129
  vx(1:ifull,k)=vx(1:ifull,k)-.5*dt*(1.-epsf)*f(1:ifull)*u(1:ifull,k) ! end of Eq. 130
enddo      ! k loop

call depts1(x3d,y3d,z3d)
      
!     calculate factr for choice of nt_adv, as usually used
select case(nt_adv)
  case(0)
    factr(:)=0.
  case(3)  ! 1. up to sig=.3
    do k=1,kl
      factr(k)=stdlapse*(rdry*nritch_t/grav)
      if(sig(k)<0.3)factr(k)=0.
    end do
  case(4) ! (1, .5, 0) for sig=(1, .75, .5)
    do k=1,kl
      factr(k)=max(2.*sig(k)-1., 0.)*stdlapse*(rdry*300./grav)
    end do
  case(5) ! 1 to 0 for sig=1 to 0
    do k=1,kl
      factr(k)=sig(k)*stdlapse*(rdry*nritch_t/grav)
    end do
  case(6) ! (1, .5625, 0) for sig=(1, .5, .2)
    do k=1,kl
      factr(k)=max(0.,1.25*(sig(k)-.2)*(2.-sig(k)))*stdlapse*(rdry*nritch_t/grav)
    end do
  case(7) ! 1 up to .4, then lin decr. to .2, then 0
    do k=1,kl
      factr(k)=max(0.,min(1.,(sig(k)-.2)/(.4-.2)))*stdlapse*(rdry*nritch_t/grav)
    end do
  case(8) ! .8 up to .4, then lin decr. to .2, then 0
    do k=1,kl
      factr(k)=.8*max(0.,min(1.,(sig(k)-.2)/(.4-.2)))*stdlapse*(rdry*nritch_t/grav)
    end do
  case(9) ! (1,1,.84375,.5,.15625,0) for sig=(1,.6,.5,.4,.3,.2)  
    do k=1,kl
      factr(k)=3.*((sig(k)-.2)/.4)**2 -2.*((sig(k)-.2)/.4)**3*stdlapse*(rdry*nritch_t/grav)
      if(sig(k)>.6)factr(k)=stdlapse*(rdry*nritch_t/grav)
      if(sig(k)<.2)factr(k)=0.
    end do
  case(10) ! (1,1,.741,.259,0) for sig=(1,.5,.4,.3,.2) 
    do k=1,kl
      factr(k)=3.*((sig(k)-.2)/.3)**2 -2.*((sig(k)-.2)/.3)**3*stdlapse*(rdry*nritch_t/grav)
      if(sig(k)>.5)factr(k)=stdlapse*(rdry*nritch_t/grav)
      if(sig(k)<.2)factr(k)=0.
    end do
end select

if ( mydiag ) then
  if(tx(idjd,kl)>264.)then  !cb
    write(6,*)
    write(6,*) 'upglobal a',ktau,id,jd,tx(idjd,kl)
  endif    ! (tx(idjd,kl)>264.)
end if

if(num_hight<100)then
  do iq=1,ifull
    if(tx(iq,kl)>264.)then  !cb
      write(6,*) 'upglobal ktau,myid,iq,large_tx  ',ktau,myid,iq,tx(iq,kl)
      write (6,"('sdot_iq',9f7.3/7x,9f7.3)") sdot(iq,1:kl)
      num_hight=num_hight+1
    endif
  enddo
endif 

aa(1:ifull)=zs(1:ifull)/(rdry*nritch_t)    ! save zs/(r*t) for nt_adv schemes 
do k=1,kl   
  dd(1:ifull,k)=aa(1:ifull)
end do     ! k loop

!-------------------------moved up here May 06---------------------------
! N.B. this moved one is doing vadv on just extra pslx terms      
sdmx(:) = maxval(abs(sdot),2)
nits(:)=1+sdmx(:)/2
nvadh_pass(:)=2*nits(:) ! use - for nvadu
call vadvtvd(tx,ux,vx,nvadh_pass,nits)
if( (diag.or.nmaxpr==1) .and. mydiag )then
  write(6,*) 'in upglobal after vadv1'
  write (6,"('qg  ',3p9f8.3/4x,9f8.3)")   qg(idjd,:)
endif
!------------------------------------------------------------------

do k=1,kl   
  ! N.B. [D + dsigdot/dsig] saved in adjust5 (or updps) as pslx
  pslx(1:ifull,k)=psl(1:ifull)-pslx(1:ifull,k)*dt*.5*(1.-epst(:))
  pslx(1:ifull,k)=pslx(1:ifull,k)+aa(1:ifull)
  tx(1:ifull,k)=tx(1:ifull,k)+aa(1:ifull)*factr(k)   !cy  
end do   ! k

if(nmaxpr==1.and.nproc==1)then
  write(6,*) 'pslx_3p before advection'
  write (6,"('pslx_b',3p9f8.4)") pslx(idjd,:)
  write (6,"(i6,8i8)") (ii,ii=id-4,id+4)
  write (6,"(3p9f8.4)") ((pslx(max(min(ii+jj*il,ifull),1),nlv),ii=idjd-4,idjd+4),jj=2,-2,-1)
endif

if ( mup/=0 ) then
  call ints_bl(dd,intsch,nface,xg,yg)  ! advection on all levels
  if ( nh/=0 ) then
    ! non-hydrostatic version
    duma(1:ifull,:,1)=pslx(1:ifull,:)
    duma(1:ifull,:,2)=h_nh(1:ifull,:)
    call ints(2,duma,intsch,nface,xg,yg,1)
    pslx(1:ifull,:)=duma(1:ifull,:,1)
    h_nh(1:ifull,:)=duma(1:ifull,:,2)
  else
    ! hydrostatic version
    call ints(1,pslx,intsch,nface,xg,yg,1)
  end if ! nh/=0
  call ints(1,tx,intsch,nface,xg,yg,3)
endif    ! mup/=0

do k=1,kl
  pslx(1:ifull,k)=pslx(1:ifull,k)-dd(1:ifull,k)      
  tx(1:ifull,k) = tx(1:ifull,k)  -dd(1:ifull,k)*factr(k)
end do
!------------------------------------------------------------------
if(nmaxpr==1.and.nproc==1)then
  write(6,*) 'pslx_3p & dd after advection'
  write (6,"('pslx_a',3p9f8.4)") pslx(idjd,:)
  write (6,"('aa#',3p9f8.4)") ((aa(ii+jj*il),ii=idjd-1,idjd+1),jj=-1,1)
  write (6,"('dd1#',3p9f8.4)") ((dd(ii+jj*il,1),ii=idjd-1,idjd+1),jj=-1,1)
  write (6,"('dd_a',3p9f8.4)") dd(idjd,:)
  write (6,"('nface',18i4)") nface(idjd,:)
  write (6,"('xg',9f8.4)") xg(idjd,:)
  write (6,"('yg',9f8.4)") yg(idjd,:)
  write (6,"(i6,8i8)") (ii,ii=id-4,id+4)
  idjdd=max(5+2*il,min(idjd,ifull-4-2*il))  ! for following prints
  write (6,"(3p9f8.4)") ((pslx(ii+jj*il,nlv),ii=idjdd-4,idjdd+4),jj=2,-2,-1)
  uc(1:ifull,1)=-pslx(1:ifull,1)*dsig(1) 
  do k=2,kl
    uc(1:ifull,1)=uc(1:ifull,1)-pslx(1:ifull,k)*dsig(k)
  enddo
  write(6,*) 'integ pslx after advection'
  write (6,"(i6,8i8)") (ii,ii=id-4,id+4)
  write (6,"(3p9f8.4)") ((uc(ii+jj*il,1),ii=idjdd-4,idjdd+4),jj=2,-2,-1)
  write(6,*) 'corresp integ ps after advection'
  write (6,"(i6,8i8)") (ii,ii=id-4,id+4)
  write (6,"(-2p9f8.2)") ((1.e5*exp(uc(ii+jj*il,1)),ii=idjdd-4,idjdd+4),jj=2,-2,-1)
endif
!     now comes ux & vx section
if(diag)then
  if ( mydiag ) then
    write(6,*) 'unstaggered now as uavx and vavx: globpea uses ux & vx'
    write(6,*) 'ux ',(ux(idjd,kk),kk=1,nlv)
  end if
  call printa('uavx',ux,ktau,nlv,ia,ib,ja,jb,0.,1.)
  if (mydiag) write(6,*) 'vx ',(vx(idjd,kk),kk=1,nlv)
  call printa('vavx',vx,ktau,nlv,ia,ib,ja,jb,0.,1.)
  if ( mydiag ) write(6,*)'unstaggered u and v as uav and vav: globpea uses u & v'
  call printa('uav ',u,ktau,nlv,ia,ib,ja,jb,0.,1.)
  call printa('vav ',v,ktau,nlv,ia,ib,ja,jb,0.,1.)
endif

!      convert uavx, vavx to cartesian velocity components
do k=1,kl
  uc(1:ifull,k)=ax(1:ifull)*ux(1:ifull,k) + bx(1:ifull)*vx(1:ifull,k)
  vc(1:ifull,k)=ay(1:ifull)*ux(1:ifull,k) + by(1:ifull)*vx(1:ifull,k)
  wc(1:ifull,k)=az(1:ifull)*ux(1:ifull,k) + bz(1:ifull)*vx(1:ifull,k)
enddo
if(diag)then
  if ( mydiag ) then
    write(6,*) 'uc,vc,wc before advection'
    write (6,'(a,18e20.10)') 'uc,vc,wc ',uc(idjd,nlv),vc(idjd,nlv),wc(idjd,nlv)
  end if
  call printa('uc  ',uc,ktau,nlv,ia,ib,ja,jb,0.,1.)
  call printa('vc  ',vc,ktau,nlv,ia,ib,ja,jb,0.,1.)
  call printa('wc  ',wc,ktau,nlv,ia,ib,ja,jb,0.,1.)
  call printa('xg  ',xg,ktau,nlv,ia,ib,ja,jb,0.,1.)
  call printa('yg  ',yg,ktau,nlv,ia,ib,ja,jb,0.,1.)
  if ( mydiag ) write(6,*) 'nface ',nface(idjd,:)
endif
if(mup/=0)then
  duma(1:ifull,:,1)=uc(1:ifull,:)
  duma(1:ifull,:,2)=vc(1:ifull,:)
  duma(1:ifull,:,3)=wc(1:ifull,:)
  call ints(3,duma,intsch,nface,xg,yg,2)
  uc(1:ifull,:)=duma(1:ifull,:,1)
  vc(1:ifull,:)=duma(1:ifull,:,2)
  wc(1:ifull,:)=duma(1:ifull,:,3)
endif
if(diag)then
  if ( mydiag ) write(6,*) 'uc,vc,wc after advection'
  call printa('uc  ',uc,ktau,nlv,ia,ib,ja,jb,0.,1.)
  call printa('vc  ',vc,ktau,nlv,ia,ib,ja,jb,0.,1.)
  call printa('wc  ',wc,ktau,nlv,ia,ib,ja,jb,0.,1.)
endif

! rotate wind vector to arrival point
do k=1,kl
  ! the following normalization may be done, but has ~zero effect
  ! dena=sqrt(x3d(iq)**2+y3d(iq)**2+z3d(iq)**2)
  ! x3d(iq)=x3d(iq)/dena
  ! y3d(iq)=y3d(iq)/dena
  ! z3d(iq)=z3d(iq)/dena
  ! cross product n1xn2 into vec1
  vec1x(1:ifull) = y3d(1:ifull,k)*z(1:ifull) - y(1:ifull)*z3d(1:ifull,k)
  vec1y(1:ifull) = z3d(1:ifull,k)*x(1:ifull) - z(1:ifull)*x3d(1:ifull,k)
  vec1z(1:ifull) = x3d(1:ifull,k)*y(1:ifull) - x(1:ifull)*y3d(1:ifull,k)
  denb(1:ifull) = vec1x(1:ifull)**2 + vec1y(1:ifull)**2 + vec1z(1:ifull)**2
  ! N.B. rotation formula is singular for small denb,
  ! but the rotation is unnecessary in this case
  where (denb>1.e-4)
    vecdot(1:ifull) = x3d(1:ifull,k)*x(1:ifull) + y3d(1:ifull,k)*y(1:ifull) + z3d(1:ifull,k)*z(1:ifull)
    vec2x(1:ifull) = x3d(1:ifull,k)*vecdot(1:ifull) - x(1:ifull)
    vec2y(1:ifull) = y3d(1:ifull,k)*vecdot(1:ifull) - y(1:ifull)
    vec2z(1:ifull) = z3d(1:ifull,k)*vecdot(1:ifull) - z(1:ifull)
    vec3x(1:ifull) = x3d(1:ifull,k) - vecdot(1:ifull)*x(1:ifull)
    vec3y(1:ifull) = y3d(1:ifull,k) - vecdot(1:ifull)*y(1:ifull)
    vec3z(1:ifull) = z3d(1:ifull,k) - vecdot(1:ifull)*z(1:ifull)
    vdot1(1:ifull) = (vec1x(1:ifull)*uc(1:ifull,k) + vec1y(1:ifull)*vc(1:ifull,k) + vec1z(1:ifull)*wc(1:ifull,k))/denb(1:ifull)
    vdot2(1:ifull) = (vec2x(1:ifull)*uc(1:ifull,k) + vec2y(1:ifull)*vc(1:ifull,k) + vec2z(1:ifull)*wc(1:ifull,k))/denb(1:ifull)
    uc(1:ifull,k) = vdot1(1:ifull)*vec1x(1:ifull) + vdot2(1:ifull)*vec3x(1:ifull)
    vc(1:ifull,k) = vdot1(1:ifull)*vec1y(1:ifull) + vdot2(1:ifull)*vec3y(1:ifull)
    wc(1:ifull,k) = vdot1(1:ifull)*vec1z(1:ifull) + vdot2(1:ifull)*vec3z(1:ifull)
  end where ! (denb>1.e-4)
end do ! k
if(diag)then
  if ( mydiag )then
    iq=idjd
    k=nlv
    vec1x(iq) = y3d(iq,k)*z(iq) - y(iq)*z3d(iq,k)
    vec1y(iq) = z3d(iq,k)*x(iq) - z(iq)*x3d(iq,k)
    vec1z(iq) = x3d(iq,k)*y(iq) - x(iq)*y3d(iq,k)
    denb(iq) = vec1x(iq)**2 + vec1y(iq)**2 + vec1z(iq)**2
    write(6,*) 'uc,vc,wc after nrot; denb = ',denb
  endif
  call printa('uc  ',uc,ktau,nlv,ia,ib,ja,jb,0.,1.)
  call printa('vc  ',vc,ktau,nlv,ia,ib,ja,jb,0.,1.)
  call printa('wc  ',wc,ktau,nlv,ia,ib,ja,jb,0.,1.)
endif

! convert back to conformal-cubic velocity components (unstaggered)
! globpea: this can be sped up later
do k=1,kl
  ux(1:ifull,k) = ax(1:ifull)*uc(1:ifull,k) + ay(1:ifull)*vc(1:ifull,k) + az(1:ifull)*wc(1:ifull,k)
  vx(1:ifull,k) = bx(1:ifull)*uc(1:ifull,k) + by(1:ifull)*vc(1:ifull,k) + bz(1:ifull)*wc(1:ifull,k)
end do   ! k loop
if(diag.and.k==nlv)then
  if ( mydiag ) write(6,*) 'after advection in upglobal; unstaggered ux and vx:'
  call printa('ux  ',ux,ktau,nlv,ia,ib,ja,jb,0.,1.)
  call printa('vx  ',vx,ktau,nlv,ia,ib,ja,jb,0.,1.)
endif

if ( mspec==1 .and. mup/=0 ) then   ! advect qg after preliminary step
  if ( ldr/=0 ) then
    duma(1:ifull,:,1)=qg(1:ifull,:)
    duma(1:ifull,:,2)=qlg(1:ifull,:)
    duma(1:ifull,:,3)=qfg(1:ifull,:)
    duma(1:ifull,:,4)=qrg(1:ifull,:)
    duma(1:ifull,:,5)=rfrac(1:ifull,:)
    duma(1:ifull,:,6)=qsg(1:ifull,:)
    duma(1:ifull,:,7)=qgrg(1:ifull,:)
    duma(1:ifull,:,8)=sfrac(1:ifull,:)
    duma(1:ifull,:,9)=gfrac(1:ifull,:)
    if ( ncloud>=4 ) then
      ! prognostic cloud fraction and condensate version
      duma(1:ifull,:,10)=stratcloud(1:ifull,:)
      call ints(10,duma,intsch,nface,xg,yg,4)
      stratcloud(1:ifull,:)=min(max(duma(1:ifull,:,10),0.),1.)
    else
      ! prognostic cloud condesate version
      if ( ncloud>=3 ) then
        ntot=9
      else if ( ncloud>=2 ) then
        ntot=5
      else
        ntot=3
      end if
      call ints(ntot,duma(:,:,1:ntot),intsch,nface,xg,yg,4)
    end if
    qg(1:ifull,:)   =duma(1:ifull,:,1)
    qlg(1:ifull,:)  =duma(1:ifull,:,2)
    qfg(1:ifull,:)  =duma(1:ifull,:,3)
    qrg(1:ifull,:)  =duma(1:ifull,:,4)
    rfrac(1:ifull,:)=min(max(duma(1:ifull,:,5),0.),1.)       
    qsg(1:ifull,:)  =duma(1:ifull,:,6)
    qgrg(1:ifull,:) =duma(1:ifull,:,7)
    sfrac(1:ifull,:)=min(max(duma(1:ifull,:,8),0.),1.)
    gfrac(1:ifull,:)=min(max(duma(1:ifull,:,9),0.),1.)        
  else
    call ints(1,qg,intsch,nface,xg,yg,3)
  end if    ! ldr/=0
  if ( ngas>0.or.nextout>=4 ) then
    if ( nmaxpr==1 .and. mydiag ) then
      write (6,"('xg#',9f8.2)") diagvals(xg(:,nlv))
      write (6,"('yg#',9f8.2)") diagvals(yg(:,nlv))
      write (6,"('nface#',9i8)") diagvals(nface(:,nlv))
      write (6,"('xlat#',9f8.2)") diagvals(tr(:,nlv,ngas+1))
      write (6,"('xlon#',9f8.2)") diagvals(tr(:,nlv,ngas+2))
      write (6,"('xpre#',9f8.2)") diagvals(tr(:,nlv,ngas+3))
    endif
    if ( ngas>0 ) then
      do nstart = 1, ngas, nagg
        nend = min(nstart+nagg-1,ngas)
        ntot = nend - nstart + 1
        call ints(ntot,tr(:,:,nstart:nend),intsch,nface,xg,yg,5)
      end do
    end if
    if ( nmaxpr==1 .and. mydiag ) then
      write (6,"('ylat#',9f8.2)") diagvals(tr(:,nlv,ngas+1))
      write (6,"('ylon#',9f8.2)") diagvals(tr(:,nlv,ngas+2))
      write (6,"('ypre#',9f8.2)") diagvals(tr(:,nlv,ngas+3))
    endif
  endif  ! (ngas>0.or.nextout>=4)
  if ( nvmix==6 ) then
    duma(1:ifull,:,1)=tke(1:ifull,:)
    duma(1:ifull,:,2)=eps(1:ifull,:)
    call ints(2,duma,intsch,nface,xg,yg,3)
    tke(1:ifull,:)=duma(1:ifull,:,1)
    eps(1:ifull,:)=duma(1:ifull,:,2)
  endif                 ! nvmix==6
  if ( abs(iaero)==2 ) then
    call ints(naero,xtg,intsch,nface,xg,yg,5)
  end if
end if     ! mspec==1


sdot(:,2:kl)=sbar(:,:)
if (mod(ktau,nmaxpr)==0.and.mydiag) then
  write(6,*) 'upglobal ktau,sdmx,nits,nvadh_pass ',ktau,sdmx(idjd),nits(idjd),nvadh_pass(idjd)
endif
if ( (diag.or.nmaxpr==1) .and. mydiag ) then
  write(6,*) 'in upglobal before vadv2'
  write (6,"('qg  ',3p9f8.3/4x,9f8.3)")   qg(idjd,:)
endif
call vadvtvd(tx,ux,vx,nvadh_pass,nits)
if ( (diag.or.nmaxpr==1) .and. mydiag ) then
  write(6,*) 'in upglobal after vadv2'
  write (6,"('qg  ',3p9f8.3/4x,9f8.3)")   qg(idjd,:)
endif

! adding later (after 2nd vadv) than for npex=1
ux(1:ifull,:) = ux(1:ifull,:)+0.5*dt*un(1:ifull,:) ! dyn contrib
vx(1:ifull,:) = vx(1:ifull,:)+0.5*dt*vn(1:ifull,:) ! dyn contrib
      
! second part of usual m=6 coriolis treatment (after 2nd vadv)
do k=1,kl
  ! incorporate coriolis terms (done here as for m=6 instead of in adjust5)
  tempry(1:ifull) = ux(1:ifull,k)+.5*dt*(1.+epsf)*f(1:ifull)*vx(1:ifull,k) ! Eq. 133
  vx(1:ifull,k)   = vx(1:ifull,k)-.5*dt*(1.+epsf)*f(1:ifull)*ux(1:ifull,k) ! Eq. 134
  ux(1:ifull,k)   = tempry(1:ifull)
enddo

tx(1:ifull,:)=tx(1:ifull,:)+.5*dt*tn(1:ifull,:) 

!     now interpolate ux,vx to the staggered grid
call staguv(ux,vx,ux,vx)

if( diag) then
  if(mydiag)then
    write(6,*) 'near end of upglobal staggered ux and vx:'
    write(6,*) 'un_u ',un(idjd,:)
    write(6,*) 'vn_u ',vn(idjd,:)
    write(6,*) 'tn_u ',tn(idjd,:)
    write (6,"('tx_u1',9f8.2/5x,9f8.2)") tx(idjd,:)
  endif
  call printa('ux  ',ux,ktau,nlv,ia,ib,ja,jb,0.,1.)
  call printa('vx  ',vx,ktau,nlv,ia,ib,ja,jb,0.,1.)
endif

if(ntest>4)then
  ! diagnostic check for unstable layers
  do iq=1,ifull
    do k=1,kl
      theta(iq,k)=tx(iq,k)*sig(k)**(-rdry/cp)
    enddo
    do k=ntest,kl   ! e.g. 8,kl
      if(theta(iq,k)<theta(iq,k-1))then  ! based on tx
        write(6,*)"unstable layer in upglobal for ktau,iq,k's,del ",ktau,iq,k-1,k,theta(iq,k-1)-theta(iq,k)
        write (6,"('theta',9f7.2/5x,9f7.2)") theta(iq,:)
        write (6,"('sdot',9f7.3/4x,9f7.3)")  sdot(iq,1:kl)
        numunstab=numunstab+1
      endif
    enddo
  enddo
endif

if( ( diag.or.nmaxpr==1) .and. mydiag ) then
  write(6,*) 'near end of upglobal for ktau= ',ktau
  write (6,"('tx_u2',9f8.2/5x,9f8.2)")  tx(idjd,:)
  write (6,"('qg_u',3p9f8.3/4x,9f8.3)") qg(idjd,:)
  write (6,"('ql_u',3p9f8.3/4x,9f8.3)") qlg(idjd,:)
  write (6,"('qf_u',3p9f8.3/4x,9f8.3)") qfg(idjd,:)
endif 

call END_LOG(upglobal_end)
      
return
end subroutine upglobal
