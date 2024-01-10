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
    
subroutine updps(iadj)    

use arrays_m
use cc_mpi, only : boundsuv,bounds,mydiag
use const_phys
use diag_m             ! for calls to maxmin
use indices_m
use newmpar_m
use map_m
use nlin_m             ! savs
use parm_m
use parmdyn_m
use parmhor_m
use savuvt_m
use savuv1_m
use sigs_m
use staguvmod
use vecsuv_m
use vvel_m
use xarrs_m
use xyzinfo_m

implicit none

integer, save :: num = 0
integer iq, k, ind, iadj
real, dimension(ifull+iextra,kl) :: cc, dd
real, dimension(ifull,kl) :: d, sdotin, pslxin, omgfin
real, dimension(ifull) :: derpsl
real, dimension(kl) :: pslx_k

! Always called for first time step
! Called every step for mup/=1
! Usual is mup=1, using simple centred (only first step)
!          mup=2, using simple centred (every time step)
!          mup=4 to use prior pslx (no longer here)
!          mup=5 (somewhat similar to 1) & 6 calc via staguv
!          mup<0 various others
      
if ( mup<=-4 ) then
  if ( ktau<3 ) then
    return
  else
    sdotin(:,1:kl) = sdot(:,1:kl)
    omgfin(:,1:kl) = dpsldt(:,1:kl)
    pslxin(:,1:kl) = pslx(1:ifull,1:kl)
  end if
end if  ! (mup<=-4)
      
if ( mup>=5 ) then
  ! MJT notes - make staguv a module so that assumed shape
  ! arrays can be used to avoid tmpu, tmpv, tmpcc and tmpdd
  call staguv(u,v,cc,dd)
  do k = 1,kl
    cc(1:ifull,k) = cc(1:ifull,k)/emu(1:ifull)
    dd(1:ifull,k) = dd(1:ifull,k)/emv(1:ifull)
  enddo
  call boundsuv(cc,dd)
  do k = 1,kl
    ! N.B. this div is calculated on the staggered grid
    d(1:ifull,k) = (cc(1:ifull,k)-cc(iwu,k) +dd(1:ifull,k)-dd(isv,k))*em(1:ifull)**2/ds
  enddo    ! k  loop
endif       ! (mup>=5) 
      
if ( mup<5 ) then       ! ************usual***********
  call boundsuv(u,v)
  do k = 1,kl
    ! N.B. this div is calculated on the non-staggered grid
    ! but in adjust5 it is on the staggered grid
    d(1:ifull,k) = (u(ieu,k)/em(ie)-u(iwu,k)/em(iw)+v(inv,k)/em(in)-v(isv,k)/em(is))*em(1:ifull)**2/(2.*ds)
  enddo    ! k  loop
  if ( mydiag .and. nmaxpr==1 ) then
    write(6,*) "d,em,ds ",d(idjd,nlv),em(idjd),ds
    write(6,*) "u ",u(ieu(idjd),nlv),u(iwu(idjd),nlv)
    write(6,*) "v ",v(inv(idjd),nlv),v(isv(idjd),nlv)
    write(6,*) "1em ",em(ie(idjd)),em(iw(idjd))
    write(6,*) "2em ",em(in(idjd)),em(is(idjd))
  end if
endif       ! (mup<5) 

if ( mup==-1 .or. mup<=-4 .or. (mup==-3.and.iadj==0) ) then
  tx(:,1) = zs(:)/(rdry*300.)
  tx(:,2) = psl(:) + tx(:,1)
  call bounds(tx)
  do k = 1,kl
    pslx(1:ifull,k) = -em(1:ifull)*(u(1:ifull,k)*(tx(ie,2)-tx(iw,2))+v(1:ifull,k)*(tx(in,2)-tx(is,2)))/(2.*ds)
    where ( u(1:ifull,k)>0. )
      pslx(1:ifull,k) = pslx(1:ifull,k) + em(1:ifull)*(tx(ie,1)-tx(1:ifull,1))*u(1:ifull,k)/ds
    elsewhere
      pslx(1:ifull,k) = pslx(1:ifull,k) + em(1:ifull)*(tx(1:ifull,1)-tx(iw,1))*u(1:ifull,k)/ds
    end where 
    where ( v(1:ifull,k)>0. )
      pslx(1:ifull,k) = pslx(1:ifull,k) + em(1:ifull)*(tx(in,1)-tx(1:ifull,1))*v(1:ifull,k)/ds
    elsewhere
      pslx(1:ifull,k) = pslx(1:ifull,k) + em(1:ifull)*(tx(1:ifull,1)-tx(is,1))*v(1:ifull,k)/ds
    end where
  enddo    ! k  loop
elseif ( mup==-2 .or. (mup==-3.and.iadj==1) ) then
  tx(:,1) = zs(:)/(rdry*300.)
  tx(:,2) = psl(:) + tx(:,1)
  call bounds(tx)
  do k = 1,kl
    pslx(1:ifull,k) = -em(1:ifull)*(u(1:ifull,k)*(tx(ie,2)-tx(iw,2))+v(1:ifull,k)*(tx(in,2)-tx(is,2)))/(2.*ds)
    where ( u(1:ifull,k)<0. )
      pslx(1:ifull,k) = pslx(1:ifull,k) + em(1:ifull)*(tx(ie,1)-tx(1:ifull,1))*u(1:ifull,k)/ds
    elsewhere
      pslx(1:ifull,k) = pslx(1:ifull,k) + em(1:ifull)*(tx(1:ifull,1)-tx(iw,1))*u(1:ifull,k)/ds
    end where
    where ( v(1:ifull,k)<0. )
      pslx(1:ifull,k) = pslx(1:ifull,k) + em(1:ifull)*(tx(in,1)-tx(1:ifull,1))*v(1:ifull,k)/ds
    elsewhere
      pslx(1:ifull,k) = pslx(1:ifull,k) + em(1:ifull)*(tx(1:ifull,1)-tx(is,1))*v(1:ifull,k)/ds
    end where 
  enddo    ! k  loop
elseif ( mup<4 .or. num==0 ) then  ! ************usual***********
  num = 1
  ! put -D(ln(psl))/Dt +d(ln(psl))/dt into pslx (i.e. Nps)
  do k = 1,kl
    pslx(1:ifull,k) = -em(1:ifull)*(u(1:ifull,k)*(psl(ie)-psl(iw))+v(1:ifull,k)*(psl(in)-psl(is)))/(2.*ds)
  enddo    ! k  loop
  if ( mydiag .and. nmaxpr==1 ) then
    write(6,*) "pslx,em,ds ",pslx(idjd,nlv),em(idjd),ds
    write(6,*) "u,v ",u(idjd,nlv),v(idjd,nlv)
    write(6,*) "psl ",psl(ie(idjd)),psl(iw(idjd)),psl(in(idjd)),psl(is(idjd))
  end if
  if ( diag ) then
    call bounds(ps)
    if ( mydiag ) then
      write(6,*) 'after bounds in updps'  
      iq = idjd
      write(6,*) 'in updps'
      write(6,*) 'dhatA ',(d(iq,k)-pslx(iq,k),k=1,kl)
      do k = 1,kl
        pslx_k(k) = -em(iq)*(u(iq,k)*(ps(ie(iq))-ps(iw(iq)))+v(iq,k)*(ps(in(iq))-ps(is(iq))))/(2.*ds*ps(iq))
      enddo    ! k  loop
      write(6,*) 'dhatAA ',(d(iq,k)-pslx_k(k),k=1,kl)
      write(6,*) 'part1AA ',(d(iq,k),k=1,kl)
      write(6,*) 'part2AA ',(-pslx_k(k),k=1,kl)
    endif
  endif
endif       ! (mup<4.or.num==0)

! integrate vertically {0 to 1} to get d(ln(psl))/dt
derpsl(:) = 0.
do k = 1,kl
  derpsl(:) = derpsl(:) - dsig(k)*(pslx(1:ifull,k)-d(:,k))
enddo      ! k  loop
! put -D(ln(psl))/Dt into pslx, i.e. to equal D+d(sdot)/d(sig)
do k = 1,kl
  pslx(1:ifull,k) = -derpsl(:) + pslx(1:ifull,k)
enddo      ! k  loop
if ( mydiag .and. nmaxpr==1 ) then
  write(6,*) "pslx,derpsl ",pslx(idjd,nlv),derpsl(idjd)
  write(6,*) "dsig,d ",dsig(nlv),d(idjd,nlv)
end if

! calculate sdot (at level k-.5) by vert. integ. {0 to sig(k-.5)} 
sdot(:,1) = 0.
sdot(:,kl+1) = 0.
do k = kl,2,-1
  sdot(:,k) = sdot(:,k+1) - dsig(k)*(pslx(1:ifull,k)-d(:,k))
enddo      ! k  loop

! full-level omega/ps into omgf (equivalenced to dpsldt)
do k = 1,kl
  dpsldt(:,k) = rata(k)*sdot(:,k+1) + ratb(k)*sdot(:,k) - sig(k)*pslx(1:ifull,k)
enddo      ! k  loop
      
!------------------------------------------------------------------      
! convert sdot (at level k-.5) into units of grid-steps/timestep
! half-level sigma-dot is used for vertical advection
do k = 2,kl
  sdot(:,k) = sdot(:,k)*dt/(sig(k)-sig(k-1))
enddo    ! k  loop

if ( mup<=-4 ) then
  do iq = 1,ifull
    ind = 0
    do k = 2,kl
     !if ( (sign(1.,sdot(iq,k)).ne.sign(1.,savs(iq,k))) .and. (sign(1.,savs1(iq,k)).ne.sign(1.,savs(iq,k))) ) ind = 1
     if ( sdot(iq,k)*savs(iq,k)<0. .and. savs1(iq,k)*savs(iq,k)<0. ) ind = 1
    enddo
    if ( ind==0 ) then
      do k = 1,kl
        sdot(iq,k) = sdotin(iq,k)
        dpsldt(iq,k) = omgfin(iq,k)
        pslx(iq,k) = pslxin(iq,k)
      enddo
    else
      do k = 1,kl
        t(iq,k) = t(iq,k) + dt*tn(iq,k)
        tn(iq,k) = 0.
        if ( mup==-5 ) then
          u(iq,k) = u(iq,k) + dt*un(iq,k)
          un(iq,k) = 0.
          v(iq,k) = v(iq,k) + dt*vn(iq,k)
          vn(iq,k) = 0.
        endif
      enddo
    endif  ! (ind==0) .. else ..
  enddo   ! iq loop
endif     ! (mup<=-4)

if ( diag .or. nmaxpr==1 ) then
  if ( mydiag ) then
    write(6,*) 'in updps'
    write (6,"('div5p',10f8.2)") d(idjd,:)
    iq = idjd
    k = nlv
    write(6,*) 'iq,ie,iw,in,is ',iq,ie(iq),iw(iq),in(iq),is(iq)
    write(6,*) 'em_iq,ie,iw,in,is ',em(iq),em(ie(iq)),em(iw(iq)),em(in(iq)),em(is(iq))
    write(6,*) 'iq,ieu,iwu,inv,isv ',iq,ieu(iq),iwu(iq),inv(iq),isv(iq)
    write(6,*) 'u_iq,ieu,iwu ',u(iq,k),u(ieu(iq),k),u(iwu(iq),k)
    write(6,*) 'v_iq,inv,isv ',v(iq,k),v(inv(iq),k),v(isv(iq),k)
  endif
  call printa('div5',d,ktau,nlv,ia,ib,ja,jb,0.,1.e5)
  call printa('u',u,ktau,nlv,ia,ib,ja,jb,0.,1.)
  call printa('v',u,ktau,nlv,ia,ib,ja,jb,0.,1.)
endif

return
end subroutine updps
