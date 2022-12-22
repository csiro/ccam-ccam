! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2021 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
module staguvmod

private
public staguv, unstaguv

integer, parameter :: itnmax=3

contains

    
subroutine staguv(u,v,uout,vout)

! stripped down version just with nstag=nstagu=3, mstagpt=-3      
! also now includes nstag=nstagu=4 & 5
! staguv    may be called from adjust5, upglobal
! unstaguv  may be called from adjust5,  nonlin
! nstag now in parm.h  ! for nstag=3-5 staguv: jlm reversible 
!                      ! -ve switches every abs(nstag) time steps
! nstagu now in parm.h ! same but for unstaguv
! N.B. staguv & unstaguv previously required 2D arrays as input
! - no longer, as transferred here to uin and vin
! unstaggered u & v as input; staggered as output

use cc_mpi
use indices_m
use map_m
use newmpar_m
use parm_m
use parmdyn_m
use vecsuv_m

implicit none

real, dimension(:,:), intent(inout)  :: u, v ! in case u=uout and v=vout
real, dimension(:,:), intent(out) :: uout, vout
real, dimension(ifull+iextra,size(u,2)) :: ua, va, ud, vd, uin, vin
real, dimension(ifull,size(u,2)) :: ug, vg
integer :: itn, kx, iq, k

call START_LOG(stag_begin)

kx = size(u,2)

#ifdef debug
if(nmaxpr==1.and.mydiag)then
  write(6,*) '  stag_ktau,nstag,nstagu',ktau,nstag,nstagu
endif
#endif

! Copying could be avoided if input arrays were dimensioned ifull+iextra
uin(1:ifull,1:kx) = u(1:ifull,1:kx)
vin(1:ifull,1:kx) = v(1:ifull,1:kx)
      
if (abs(nstag)<3) then
  call boundsuv(uin,vin,stag=2)
  do k = 1,kx
    do iq = 1,ifull
      uout(iq,k)=(9.*(uin(ieu(iq),k)+uin(iq,k))-uin(iwu(iq),k)-uin(ieeu(iq),k))/16.
      vout(iq,k)=(9.*(vin(inv(iq),k)+vin(iq,k))-vin(isv(iq),k)-vin(innv(iq),k))/16.
    end do
  end do  
  return
endif  ! (nstag==0)

if ( nstag==3 ) then
  call boundsuv(uin,vin,stag=1) ! inv, innv, ieu, ieeu
         
  !$omp parallel

  ! precalculate rhs terms with iwwu2 & issv2
  !$omp do schedule(static) private(k,iq)
  do k = 1,kx
    do iq = 1,ifull  
      ud(iq,k)=uin(iq,k)/2.+uin(ieu(iq),k)+uin(ieeu(iq),k)/10.
      vd(iq,k)=vin(iq,k)/2.+vin(inv(iq),k)+vin(innv(iq),k)/10.
    end do
  end do
  !$omp end do

  !$omp master
  call boundsuv(ud,vd,stag=-10) ! inv, ieu
  !$omp end master
  !$omp barrier
  
  !$omp do schedule(static) private(k,iq)
  do k = 1,kx
    do iq = 1,ifull  
      ua(iq,k)=ud(iq,k)-ud(ieu(iq),k)/2. ! 1st guess
      va(iq,k)=vd(iq,k)-vd(inv(iq),k)/2. ! 1st guess
      ug(iq,k)=ua(iq,k)
      vg(iq,k)=va(iq,k)
    end do
  end do
  !$omp end do

  do itn=1,itnmax-1        ! each loop is a double iteration
    !$omp master
    call boundsuv(ua,va,stag=2) ! isv, inv, innv, iwu, ieu, ieeu
    !$omp end master
    !$omp barrier
    !$omp do schedule(static) private(k,iq)
    do k = 1,kx
      do iq = 1,ifull  
        uin(iq,k)=(ug(iq,k)-ua(iwu(iq),k)/10. +ua(ieeu(iq),k)/4.)/.95
        vin(iq,k)=(vg(iq,k)-va(isv(iq),k)/10. +va(innv(iq),k)/4.)/.95
      end do
    end do
    !$omp end do
    !$omp master
    call boundsuv(uin,vin,stag=2) ! isv, inv, innv, iwu, ieu, ieeu
    !$omp end master
    !$omp barrier
    !$omp do schedule(static) private(k,iq)
    do k = 1,kx
      do iq = 1,ifull  
        ua(iq,k)=(ug(iq,k)-uin(iwu(iq),k)/10. +uin(ieeu(iq),k)/4.)/.95
        va(iq,k)=(vg(iq,k)-vin(isv(iq),k)/10. +vin(innv(iq),k)/4.)/.95
      end do
    end do
    !$omp end do
  end do                  ! itn=1,itnmax
  !$omp master
  call boundsuv(ua,va,stag=2) ! isv, inv, innv, iwu, ieu, ieeu
  !$omp end master
  !$omp barrier
  !$omp do schedule(static) private(k,iq)
  do k = 1,kx
    do iq = 1,ifull  
      uin(iq,k)=(ug(iq,k)-ua(iwu(iq),k)/10. +ua(ieeu(iq),k)/4.)/.95
      vin(iq,k)=(vg(iq,k)-va(isv(iq),k)/10. +va(innv(iq),k)/4.)/.95
    end do
  end do  
  !$omp end do
  !$omp master
  call boundsuv(uin,vin,stag=2) ! isv, inv, innv, iwu, ieu, ieeu
  !$omp end master
  !$omp barrier
  !$omp do schedule(static) private(k,iq)
  do k = 1,kx
    do iq = 1,ifull  
      uout(iq,k)=(ug(iq,k)-uin(iwu(iq),k)/10. +uin(ieeu(iq),k)/4.)/.95
      vout(iq,k)=(vg(iq,k)-vin(isv(iq),k)/10. +vin(innv(iq),k)/4.)/.95
    end do
  end do
  !$omp end do
  !$omp end parallel

else !if ( nstag==4 ) then
  call boundsuv(uin,vin,stag=3) ! issv, isv, inv, iwwu, iwu, ieu

  !$omp parallel
  !$omp do schedule(static) private(k,iq)
  do k = 1,kx
    do iq = 1,ifull  
      ua(iq,k)=-0.05*uin(iwwu(iq),k)-0.4*uin(iwu(iq),k)+0.75*uin(iq,k)+0.5*uin(ieu(iq),k) ! 1st guess
      va(iq,k)=-0.05*vin(issv(iq),k)-0.4*vin(isv(iq),k)+0.75*vin(iq,k)+0.5*vin(inv(iq),k) ! 1st guess
      ug(iq,k)=ua(iq,k)
      vg(iq,k)=va(iq,k)
    end do
  end do
  !$omp end do

  do itn=1,itnmax-1        ! each loop is a double iteration
    !$omp master
    call boundsuv(ua,va,stag=3) ! issv, isv, inv, iwwu, iwu, ieu
    !$omp end master
    !$omp barrier
    !$omp do schedule(static) private(k,iq)
    do k = 1,kx
      do iq = 1,ifull  
        uin(iq,k)=(ug(iq,k)-ua(ieu(iq),k)/10. +ua(iwwu(iq),k)/4.)/.95
        vin(iq,k)=(vg(iq,k)-va(inv(iq),k)/10. +va(issv(iq),k)/4.)/.95
      end do
    end do
    !$omp end do
    !$omp master
    call boundsuv(uin,vin,stag=3) ! issv, isv, inv, iwwu, iwu, ieu
    !$omp end master
    !$omp barrier
    !$omp do schedule(static) private(k,iq)
    do k = 1,kx
      do iq = 1,ifull  
        ua(iq,k)=(ug(iq,k)-uin(ieu(iq),k)/10. +uin(iwwu(iq),k)/4.)/.95
        va(iq,k)=(vg(iq,k)-vin(inv(iq),k)/10. +vin(issv(iq),k)/4.)/.95
      end do
    end do  
    !$omp end do
  end do                 ! itn=1,itnmax
  !$omp master
  call boundsuv(ua,va,stag=3) ! issv, isv, inv, iwwu, iwu, ieu
  !$omp end master
  !$omp barrier
  !$omp do schedule(static) private(k,iq)
  do k = 1,kx
    do iq = 1,ifull      
      uin(iq,k)=(ug(iq,k)-ua(ieu(iq),k)/10. +ua(iwwu(iq),k)/4.)/.95
      vin(iq,k)=(vg(iq,k)-va(inv(iq),k)/10. +va(issv(iq),k)/4.)/.95
    end do
  end do  
  !$omp end do
  !$omp master
  call boundsuv(uin,vin,stag=3) ! issv, isv, inv, iwwu, iwu, ieu
  !$omp end master
  !$omp barrier
  !$omp do schedule(static) private(k,iq)
  do k = 1,kx
    do iq = 1,ifull      
      uout(iq,k)=(ug(iq,k)-uin(ieu(iq),k)/10. +uin(iwwu(iq),k)/4.)/.95
      vout(iq,k)=(vg(iq,k)-vin(inv(iq),k)/10. +vin(issv(iq),k)/4.)/.95
    end do
  end do  
  !$omp end do
  !$omp end parallel
 
end if

call END_LOG(stag_end)

return
end subroutine staguv


subroutine unstaguv(u,v,uout,vout)

!     staggered u & v as input; unstaggered as output

use cc_mpi
use indices_m
use map_m
use newmpar_m
use parm_m
use parmdyn_m
use vecsuv_m

implicit none

real, dimension(:,:), intent(inout)  :: u, v ! in case u=uout and v=vout
real, dimension(:,:), intent(out) :: uout, vout
real, dimension(ifull+iextra,size(u,2)) :: ua, va, ud, vd, uin, vin
real, dimension(ifull,size(u,2)) :: ug, vg
integer :: itn, kx, iq, k

call START_LOG(stag_begin)

kx = size(u,2)

#ifdef debug
if(nmaxpr==1.and.mydiag)then
  write(6,*) 'unstag_ktau,nstag,nstagu',ktau,nstag,nstagu
endif
#endif

uin(1:ifull,1:kx) = u(1:ifull,1:kx)
vin(1:ifull,1:kx) = v(1:ifull,1:kx)

if (abs(nstagu)<3) then
  call boundsuv(uin,vin,stag=3)
  do k = 1,kx
    do iq = 1,ifull
      uout(iq,k)=(9.*(uin(iwu(iq),k)+uin(iq,k))-uin(iwwu(iq),k)-uin(ieu(iq),k))/16.
      vout(iq,k)=(9.*(vin(isv(iq),k)+vin(iq,k))-vin(issv(iq),k)-vin(inv(iq),k))/16.
    end do
  end do  
  return
endif  ! (nstagu==0)

if ( nstagu==3 ) then
  call boundsuv(uin,vin,stag=5) ! issv, isv, iwwu, iwu
  ! precalculate rhs terms with iwwu2 & issv2
  !$omp parallel
  !$omp do schedule(static) private(k,iq)
  do k = 1,kx
    do iq = 1,ifull  
      ud(iq,k)=uin(iq,k)/2.+uin(iwu(iq),k)+uin(iwwu(iq),k)/10.
      vd(iq,k)=vin(iq,k)/2.+vin(isv(iq),k)+vin(issv(iq),k)/10.
    end do
  end do  
  !$omp end do

  !$omp master
  call boundsuv(ud,vd,stag=-9) ! isv, iwu
  !$omp end master
  !$omp barrier
  !$omp do schedule(static) private(k,iq)
  do k = 1,kx
    do iq = 1,ifull  
      ua(iq,k)=ud(iq,k)-ud(iwu(iq),k)/2. ! 1st guess
      va(iq,k)=vd(iq,k)-vd(isv(iq),k)/2. ! 1st guess
      ug(iq,k)=ua(iq,k)
      vg(iq,k)=va(iq,k)
    end do
  end do  
  !$omp end do

  do itn=1,itnmax-1        ! each loop is a double iteration
    !$omp master  
    call boundsuv(ua,va,stag=3) ! issv, isv, inv, iwwu, iwu, ieu
    !$omp end master
    !$omp barrier
    !$omp do schedule(static) private(k,iq)    
    do k = 1,kx
      do iq = 1,ifull  
        uin(iq,k)=(ug(iq,k)-ua(ieu(iq),k)/10. +ua(iwwu(iq),k)/4.)/.95
        vin(iq,k)=(vg(iq,k)-va(inv(iq),k)/10. +va(issv(iq),k)/4.)/.95
      end do
    end do
    !$omp end do
    !$omp master
    call boundsuv(uin,vin,stag=3) ! issv, isv, inv, iwwu, iwu, ieu
    !$omp end master
    !$omp barrier
    !$omp do schedule(static) private(k,iq)
    do k = 1,kx
      do iq = 1,ifull  
        ua(iq,k)=(ug(iq,k)-uin(ieu(iq),k)/10. +uin(iwwu(iq),k)/4.)/.95
        va(iq,k)=(vg(iq,k)-vin(inv(iq),k)/10. +vin(issv(iq),k)/4.)/.95
      end do
    end do
    !$omp end do
  end do                 ! itn=1,itnmax
  !$omp master
  call boundsuv(ua,va,stag=3) ! issv, isv, inv, iwwu, iwu, ieu
  !$omp end master
  !$omp barrier
  !$omp do schedule(static) private(k,iq)
  do k = 1,kx
    do iq = 1,ifull  
      uin(iq,k)=(ug(iq,k)-ua(ieu(iq),k)/10. +ua(iwwu(iq),k)/4.)/.95
      vin(iq,k)=(vg(iq,k)-va(inv(iq),k)/10. +va(issv(iq),k)/4.)/.95
    end do
  end do  
  !$omp end do
  !$omp master
  call boundsuv(uin,vin,stag=3) ! issv, isv, inv, iwwu, iwu, ieu
  !$omp end master
  !$omp barrier
  !$omp do schedule(static) private(k,iq)
  do k = 1,kx
    do iq = 1,ifull  
      uout(iq,k)=(ug(iq,k)-uin(ieu(iq),k)/10. +uin(iwwu(iq),k)/4.)/.95
      vout(iq,k)=(vg(iq,k)-vin(inv(iq),k)/10. +vin(issv(iq),k)/4.)/.95
    end do
  end do
  !$omp end do
  !$omp end parallel

else !if ( nstagu==4 ) then
  call boundsuv(uin,vin,stag=2) ! isv, inv, innv, iwu, ieu, ieeu

  !$omp parallel
  !$omp do schedule(static) private(k,iq)
  do k = 1,kx
    do iq = 1,ifull  
      ua(iq,k)=-0.05*uin(ieeu(iq),k)-0.4*uin(ieu(iq),k)+0.75*uin(iq,k)+0.5*uin(iwu(iq),k) ! 1st guess
      va(iq,k)=-0.05*vin(innv(iq),k)-0.4*vin(inv(iq),k)+0.75*vin(iq,k)+0.5*vin(isv(iq),k) ! 1st guess
      ug(iq,k)=ua(iq,k)
      vg(iq,k)=va(iq,k)
    end do
  end do  
  !$omp end do

  do itn=1,itnmax-1        ! each loop is a double iteration
    !$omp master  
    call boundsuv(ua,va,stag=2) ! isv, inv, innv, iwu, ieu, ieeu
    !$omp end master
    !$omp barrier
    !$omp do schedule(static) private(k,iq)
    do k = 1,kx
      do iq = 1,ifull  
        uin(iq,k)=(ug(iq,k)-ua(iwu(iq),k)/10. +ua(ieeu(iq),k)/4.)/.95
        vin(iq,k)=(vg(iq,k)-va(isv(iq),k)/10. +va(innv(iq),k)/4.)/.95
      end do
    end do  
    !$omp end do
    !$omp master
    call boundsuv(uin,vin,stag=2) ! isv, inv, innv, iwu, ieu, ieeu
    !$omp end master
    !$omp barrier
    !$omp do schedule(static) private(k,iq)
    do k = 1,kx
      do iq = 1,ifull  
        ua(iq,k)=(ug(iq,k)-uin(iwu(iq),k)/10. +uin(ieeu(iq),k)/4.)/.95
        va(iq,k)=(vg(iq,k)-vin(isv(iq),k)/10. +vin(innv(iq),k)/4.)/.95
      end do
    end do  
    !$omp end do
  enddo                  ! itn=1,itnmax
  !$omp master
  call boundsuv(ua,va,stag=2) ! isv, inv, innv, iwu, ieu, ieeu
  !$omp end master
  !$omp barrier
  !$omp do schedule(static) private(k,iq)
  do k = 1,kx
    do iq = 1,ifull      
      uin(iq,k)=(ug(iq,k)-ua(iwu(iq),k)/10. +ua(ieeu(iq),k)/4.)/.95
      vin(iq,k)=(vg(iq,k)-va(isv(iq),k)/10. +va(innv(iq),k)/4.)/.95
    end do
  end do  
  !$omp end do
  !$omp master
  call boundsuv(uin,vin,stag=2) ! isv, inv, innv, iwu, ieu, ieeu
  !$omp end master
  !$omp barrier
  !$omp do schedule(static) private(k,iq)
  do k = 1,kx
    do iq = 1,ifull  
      uout(iq,k)=(ug(iq,k)-uin(iwu(iq),k)/10. +uin(ieeu(iq),k)/4.)/.95
      vout(iq,k)=(vg(iq,k)-vin(isv(iq),k)/10. +vin(innv(iq),k)/4.)/.95
    end do
  end do  
  !$omp end do
  !$omp end parallel
      
end if

call END_LOG(stag_end)

return
end subroutine unstaguv

end module staguvmod