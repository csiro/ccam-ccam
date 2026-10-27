! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2026 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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

interface staguv
  module procedure staguv3, staguv3x
end interface

interface unstaguv
  module procedure unstaguv3, unstaguv3x
end interface

contains

    
subroutine staguv3(u,v,uout,vout)

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

use atmosphere_arrays_m, only : nstag, nstagu
use cc_mpi
use indices_m
use map_m
use newmpar_m
use parm_m
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
do k = 1,kx
  uin(1:ifull,k) = u(1:ifull,k)
  vin(1:ifull,k) = v(1:ifull,k)
end do
      
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
         
  ! precalculate rhs terms with iwwu2 & issv2
  do k = 1,kx
    do iq = 1,ifull  
      ud(iq,k) = uin(iq,k)/2.+uin(ieu(iq),k)+uin(ieeu(iq),k)/10.
      vd(iq,k) = vin(iq,k)/2.+vin(inv(iq),k)+vin(innv(iq),k)/10.
    end do
  end do
  
  call boundsuv(ud,vd,stag=-10) ! inv, ieu
  
  do k = 1,kx
    do iq = 1,ifull  
      ua(iq,k) = ud(iq,k)-ud(ieu(iq),k)/2. ! 1st guess
      va(iq,k) = vd(iq,k)-vd(inv(iq),k)/2. ! 1st guess
      ug(iq,k) = ua(iq,k)
      vg(iq,k) = va(iq,k)
    end do
  end do

  do itn=1,itnmax        ! each loop is a double iteration
    call boundsuv(ua,va,stag=2) ! isv, inv, innv, iwu, ieu, ieeu
    do k = 1,kx
      do iq = 1,ifull  
        uin(iq,k) = (ug(iq,k)-ua(iwu(iq),k)/10.+ua(ieeu(iq),k)/4.)/.95
        vin(iq,k) = (vg(iq,k)-va(isv(iq),k)/10.+va(innv(iq),k)/4.)/.95
      end do
    end do
    call boundsuv(uin,vin,stag=2) ! isv, inv, innv, iwu, ieu, ieeu
    do k = 1,kx
      do iq = 1,ifull  
        ua(iq,k) = (ug(iq,k)-uin(iwu(iq),k)/10. +uin(ieeu(iq),k)/4.)/.95
        va(iq,k) = (vg(iq,k)-vin(isv(iq),k)/10. +vin(innv(iq),k)/4.)/.95
      end do
    end do
  end do                  ! itn=1,itnmax

else !if ( nstag==4 ) then
  call boundsuv(uin,vin,stag=3) ! issv, isv, inv, iwwu, iwu, ieu

  do k = 1,kx
    do iq = 1,ifull  
      ua(iq,k) = -0.05*uin(iwwu(iq),k)-0.4*uin(iwu(iq),k)+0.75*uin(iq,k)+0.5*uin(ieu(iq),k) ! 1st guess
      va(iq,k) = -0.05*vin(issv(iq),k)-0.4*vin(isv(iq),k)+0.75*vin(iq,k)+0.5*vin(inv(iq),k) ! 1st guess
      ug(iq,k) = ua(iq,k)
      vg(iq,k) = va(iq,k)
    end do
  end do

  do itn=1,itnmax        ! each loop is a double iteration
    call boundsuv(ua,va,stag=3) ! issv, isv, inv, iwwu, iwu, ieu
    do k = 1,kx
      do iq = 1,ifull  
        uin(iq,k) = ( ug(iq,k)-ua(ieu(iq),k)/10.+ua(iwwu(iq),k)/4. )/.95
        vin(iq,k) = ( vg(iq,k)-va(inv(iq),k)/10.+va(issv(iq),k)/4. )/.95
      end do
    end do
    call boundsuv(uin,vin,stag=3) ! issv, isv, inv, iwwu, iwu, ieu
    do k = 1,kx
      do iq = 1,ifull  
        ua(iq,k) = (ug(iq,k)-uin(ieu(iq),k)/10. +uin(iwwu(iq),k)/4.)/.95
        va(iq,k) = (vg(iq,k)-vin(inv(iq),k)/10. +vin(issv(iq),k)/4.)/.95
      end do
    end do
  end do                 ! itn=1,itnmax
 
end if

do k = 1,kx
  uout(1:ifull,k) = ua(1:ifull,k)
  vout(1:ifull,k) = va(1:ifull,k)
end do

call END_LOG(stag_end)

return
end subroutine staguv3


subroutine unstaguv3(u,v,uout,vout)

!     staggered u & v as input; unstaggered as output

use atmosphere_arrays_m, only : nstag, nstagu
use cc_mpi
use indices_m
use map_m
use newmpar_m
use parm_m
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

do k = 1,kx
  uin(1:ifull,k) = u(1:ifull,k)
  vin(1:ifull,k) = v(1:ifull,k)
end do

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
  do k = 1,kx
    do iq = 1,ifull  
      ud(iq,k) = uin(iq,k)/2.+uin(iwu(iq),k)+uin(iwwu(iq),k)/10.
      vd(iq,k) = vin(iq,k)/2.+vin(isv(iq),k)+vin(issv(iq),k)/10.
    end do
  end do  

  call boundsuv(ud,vd,stag=-9) ! isv, iwu
  
  do k = 1,kx
    do iq = 1,ifull  
      ua(iq,k) = ud(iq,k)-ud(iwu(iq),k)/2. ! 1st guess
      va(iq,k) = vd(iq,k)-vd(isv(iq),k)/2. ! 1st guess
      ug(iq,k) = ua(iq,k)
      vg(iq,k) = va(iq,k)
    end do
  end do

  do itn=1,itnmax        ! each loop is a double iteration
    call boundsuv(ua,va,stag=3) ! issv, isv, inv, iwwu, iwu, ieu
    do k = 1,kx
      do iq = 1,ifull  
        uin(iq,k) = (ug(iq,k)-ua(ieu(iq),k)/10. +ua(iwwu(iq),k)/4.)/.95
        vin(iq,k) = (vg(iq,k)-va(inv(iq),k)/10. +va(issv(iq),k)/4.)/.95
      end do
    end do
    call boundsuv(uin,vin,stag=3) ! issv, isv, inv, iwwu, iwu, ieu
    do k = 1,kx
      do iq = 1,ifull  
        ua(iq,k) = (ug(iq,k)-uin(ieu(iq),k)/10. +uin(iwwu(iq),k)/4.)/.95
        va(iq,k) = (vg(iq,k)-vin(inv(iq),k)/10. +vin(issv(iq),k)/4.)/.95
      end do
    end do
  end do                 ! itn=1,itnmax

else !if ( nstagu==4 ) then
  call boundsuv(uin,vin,stag=2) ! isv, inv, innv, iwu, ieu, ieeu

  do k = 1,kx
    do iq = 1,ifull  
      ua(iq,k) = -0.05*uin(ieeu(iq),k)-0.4*uin(ieu(iq),k)+0.75*uin(iq,k)+0.5*uin(iwu(iq),k) ! 1st guess
      va(iq,k) = -0.05*vin(innv(iq),k)-0.4*vin(inv(iq),k)+0.75*vin(iq,k)+0.5*vin(isv(iq),k) ! 1st guess
      ug(iq,k) = ua(iq,k)
      vg(iq,k) = va(iq,k)
    end do
  end do

  do itn=1,itnmax        ! each loop is a double iteration
    call boundsuv(ua,va,stag=2) ! isv, inv, innv, iwu, ieu, ieeu
    do k = 1,kx
      do iq = 1,ifull  
        uin(iq,k) = (ug(iq,k)-ua(iwu(iq),k)/10. +ua(ieeu(iq),k)/4.)/.95
        vin(iq,k) = (vg(iq,k)-va(isv(iq),k)/10. +va(innv(iq),k)/4.)/.95
      end do
    end do  
    call boundsuv(uin,vin,stag=2) ! isv, inv, innv, iwu, ieu, ieeu
    do k = 1,kx
      do iq = 1,ifull  
        ua(iq,k) = (ug(iq,k)-uin(iwu(iq),k)/10. +uin(ieeu(iq),k)/4.)/.95
        va(iq,k) = (vg(iq,k)-vin(isv(iq),k)/10. +vin(innv(iq),k)/4.)/.95
      end do
    end do
  enddo                  ! itn=1,itnmax
      
end if

do k = 1,kx
  uout(1:ifull,k) = ua(1:ifull,k)
  vout(1:ifull,k) = va(1:ifull,k)
end do

call END_LOG(stag_end)

return
end subroutine unstaguv3

subroutine staguv3x(u1,v1,u2,v2,uout1,vout1,uout2,vout2)

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

use atmosphere_arrays_m, only : nstag, nstagu
use cc_mpi
use indices_m
use map_m
use newmpar_m
use parm_m
use vecsuv_m

implicit none

real, dimension(:,:), intent(inout)  :: u1, v1 ! in case u=uout and v=vout
real, dimension(:,:), intent(inout)  :: u2, v2 ! in case u=uout and v=vout
real, dimension(:,:), intent(out) :: uout1, vout1
real, dimension(:,:), intent(out) :: uout2, vout2
real, dimension(ifull+iextra,size(u1,2)) :: ua1, va1, ud1, vd1, uin1, vin1
real, dimension(ifull,size(u1,2)) :: ug1, vg1
real, dimension(ifull+iextra,size(u2,2)) :: ua2, va2, ud2, vd2, uin2, vin2
real, dimension(ifull,size(u2,2)) :: ug2, vg2
integer :: itn, kx1, kx2, iq, k

call START_LOG(stag_begin)

kx1 = size(u1,2)
kx2 = size(u2,2)

#ifdef debug
if(nmaxpr==1.and.mydiag)then
  write(6,*) '  stag_ktau,nstag,nstagu',ktau,nstag,nstagu
endif
#endif

! Copying could be avoided if input arrays were dimensioned ifull+iextra
do k = 1,kx1
  uin1(1:ifull,k) = u1(1:ifull,k)
  vin1(1:ifull,k) = v1(1:ifull,k)
end do
do k = 1,kx2
  uin2(1:ifull,k) = u2(1:ifull,k)
  vin2(1:ifull,k) = v1(1:ifull,k)
end do
      
if (abs(nstag)<3) then
  call boundsuv(u1=uin1,v1=vin1,u2=uin2,v2=vin2,stag=2)
  do k = 1,kx1
    do iq = 1,ifull
      uout1(iq,k)=(9.*(uin1(ieu(iq),k)+uin1(iq,k))-uin1(iwu(iq),k)-uin1(ieeu(iq),k))/16.
      vout1(iq,k)=(9.*(vin1(inv(iq),k)+vin1(iq,k))-vin1(isv(iq),k)-vin1(innv(iq),k))/16.
    end do
  end do  
  do k = 1,kx2
    do iq = 1,ifull
      uout2(iq,k)=(9.*(uin2(ieu(iq),k)+uin2(iq,k))-uin2(iwu(iq),k)-uin2(ieeu(iq),k))/16.
      vout2(iq,k)=(9.*(vin2(inv(iq),k)+vin2(iq,k))-vin2(isv(iq),k)-vin2(innv(iq),k))/16.
    end do
  end do  
  return
endif  ! (nstag==0)

if ( nstag==3 ) then
  call boundsuv(u1=uin1,v1=vin1,u2=uin2,v2=vin2,stag=1) ! inv, innv, ieu, ieeu
         
  ! precalculate rhs terms with iwwu2 & issv2
  do k = 1,kx1
    do iq = 1,ifull  
      ud1(iq,k) = uin1(iq,k)/2.+uin1(ieu(iq),k)+uin1(ieeu(iq),k)/10.
      vd1(iq,k) = vin1(iq,k)/2.+vin1(inv(iq),k)+vin1(innv(iq),k)/10.
    end do
  end do
  do k = 1,kx2
    do iq = 1,ifull  
      ud2(iq,k) = uin2(iq,k)/2.+uin2(ieu(iq),k)+uin2(ieeu(iq),k)/10.
      vd2(iq,k) = vin2(iq,k)/2.+vin2(inv(iq),k)+vin2(innv(iq),k)/10.
    end do
  end do
  
  call boundsuv(u1=ud1,v1=vd1,u2=ud2,v2=vd2,stag=-10) ! inv, ieu
  
  do k = 1,kx1
    do iq = 1,ifull  
      ua1(iq,k) = ud1(iq,k)-ud1(ieu(iq),k)/2. ! 1st guess
      va1(iq,k) = vd1(iq,k)-vd1(inv(iq),k)/2. ! 1st guess
      ug1(iq,k) = ua1(iq,k)
      vg1(iq,k) = va1(iq,k)
    end do
  end do
  do k = 1,kx2
    do iq = 1,ifull  
      ua2(iq,k) = ud2(iq,k)-ud2(ieu(iq),k)/2. ! 1st guess
      va2(iq,k) = vd2(iq,k)-vd2(inv(iq),k)/2. ! 1st guess
      ug2(iq,k) = ua2(iq,k)
      vg2(iq,k) = va2(iq,k)
    end do
  end do

  do itn=1,itnmax        ! each loop is a double iteration
    call boundsuv(u1=ua1,v1=va1,u2=ua2,v2=va2,stag=2) ! isv, inv, innv, iwu, ieu, ieeu
    do k = 1,kx1
      do iq = 1,ifull  
        uin1(iq,k) = (ug1(iq,k)-ua1(iwu(iq),k)/10.+ua1(ieeu(iq),k)/4.)/.95
        vin1(iq,k) = (vg1(iq,k)-va1(isv(iq),k)/10.+va1(innv(iq),k)/4.)/.95
      end do
    end do
    do k = 1,kx2
      do iq = 1,ifull  
        uin2(iq,k) = (ug2(iq,k)-ua2(iwu(iq),k)/10.+ua2(ieeu(iq),k)/4.)/.95
        vin2(iq,k) = (vg2(iq,k)-va2(isv(iq),k)/10.+va2(innv(iq),k)/4.)/.95
      end do
    end do
    call boundsuv(u1=uin1,v1=vin1,u2=uin2,v2=vin2,stag=2) ! isv, inv, innv, iwu, ieu, ieeu
    do k = 1,kx1
      do iq = 1,ifull  
        ua1(iq,k) = (ug1(iq,k)-uin1(iwu(iq),k)/10. +uin1(ieeu(iq),k)/4.)/.95
        va1(iq,k) = (vg1(iq,k)-vin1(isv(iq),k)/10. +vin1(innv(iq),k)/4.)/.95
      end do
    end do
    do k = 1,kx2
      do iq = 1,ifull  
        ua2(iq,k) = (ug2(iq,k)-uin2(iwu(iq),k)/10. +uin2(ieeu(iq),k)/4.)/.95
        va2(iq,k) = (vg2(iq,k)-vin2(isv(iq),k)/10. +vin2(innv(iq),k)/4.)/.95
      end do
    end do
  end do                  ! itn=1,itnmax

else !if ( nstag==4 ) then
  call boundsuv(u1=uin1,v1=vin1,u2=uin2,v2=vin2,stag=3) ! issv, isv, inv, iwwu, iwu, ieu

  do k = 1,kx1
    do iq = 1,ifull  
      ua1(iq,k) = -0.05*uin1(iwwu(iq),k)-0.4*uin1(iwu(iq),k)+0.75*uin1(iq,k)+0.5*uin1(ieu(iq),k) ! 1st guess
      va1(iq,k) = -0.05*vin1(issv(iq),k)-0.4*vin1(isv(iq),k)+0.75*vin1(iq,k)+0.5*vin1(inv(iq),k) ! 1st guess
      ug1(iq,k) = ua1(iq,k)
      vg1(iq,k) = va1(iq,k)
    end do
  end do
  do k = 1,kx2
    do iq = 1,ifull  
      ua2(iq,k) = -0.05*uin2(iwwu(iq),k)-0.4*uin2(iwu(iq),k)+0.75*uin2(iq,k)+0.5*uin2(ieu(iq),k) ! 1st guess
      va2(iq,k) = -0.05*vin2(issv(iq),k)-0.4*vin2(isv(iq),k)+0.75*vin2(iq,k)+0.5*vin2(inv(iq),k) ! 1st guess
      ug2(iq,k) = ua2(iq,k)
      vg2(iq,k) = va2(iq,k)
    end do
  end do

  do itn=1,itnmax        ! each loop is a double iteration
    call boundsuv(u1=ua1,v1=va1,u2=ua2,v2=va2,stag=3) ! issv, isv, inv, iwwu, iwu, ieu
    do k = 1,kx1
      do iq = 1,ifull  
        uin1(iq,k) = ( ug1(iq,k)-ua1(ieu(iq),k)/10.+ua1(iwwu(iq),k)/4. )/.95
        vin1(iq,k) = ( vg1(iq,k)-va1(inv(iq),k)/10.+va1(issv(iq),k)/4. )/.95
      end do
    end do
    do k = 1,kx2
      do iq = 1,ifull  
        uin2(iq,k) = ( ug2(iq,k)-ua2(ieu(iq),k)/10.+ua2(iwwu(iq),k)/4. )/.95
        vin2(iq,k) = ( vg2(iq,k)-va2(inv(iq),k)/10.+va2(issv(iq),k)/4. )/.95
      end do
    end do
    call boundsuv(u1=uin1,v1=vin1,u2=uin2,v2=vin2,stag=3) ! issv, isv, inv, iwwu, iwu, ieu
    do k = 1,kx1
      do iq = 1,ifull  
        ua1(iq,k) = (ug1(iq,k)-uin1(ieu(iq),k)/10. +uin1(iwwu(iq),k)/4.)/.95
        va1(iq,k) = (vg1(iq,k)-vin1(inv(iq),k)/10. +vin1(issv(iq),k)/4.)/.95
      end do
    end do
    do k = 1,kx2
      do iq = 1,ifull  
        ua2(iq,k) = (ug2(iq,k)-uin2(ieu(iq),k)/10. +uin2(iwwu(iq),k)/4.)/.95
        va2(iq,k) = (vg2(iq,k)-vin2(inv(iq),k)/10. +vin2(issv(iq),k)/4.)/.95
      end do
    end do
  end do                 ! itn=1,itnmax
 
end if

do k = 1,kx1
  uout1(1:ifull,k) = ua1(1:ifull,k)
  vout1(1:ifull,k) = va1(1:ifull,k)
end do
do k = 1,kx1
  uout2(1:ifull,k) = ua2(1:ifull,k)
  vout2(1:ifull,k) = va2(1:ifull,k)
end do

call END_LOG(stag_end)

return
end subroutine staguv3x

subroutine unstaguv3x(u1,v1,u2,v2,uout1,vout1,uout2,vout2)

!     staggered u & v as input; unstaggered as output

use atmosphere_arrays_m, only : nstag, nstagu
use cc_mpi
use indices_m
use map_m
use newmpar_m
use parm_m
use vecsuv_m

implicit none

real, dimension(:,:), intent(inout)  :: u1, v1 ! in case u=uout and v=vout
real, dimension(:,:), intent(inout)  :: u2, v2 ! in case u=uout and v=vout
real, dimension(:,:), intent(out) :: uout1, vout1
real, dimension(:,:), intent(out) :: uout2, vout2
real, dimension(ifull+iextra,size(u1,2)) :: ua1, va1, ud1, vd1, uin1, vin1
real, dimension(ifull,size(u1,2)) :: ug1, vg1
real, dimension(ifull+iextra,size(u2,2)) :: ua2, va2, ud2, vd2, uin2, vin2
real, dimension(ifull,size(u2,2)) :: ug2, vg2
integer :: itn, kx1, kx2, iq, k

call START_LOG(stag_begin)

kx1 = size(u1,2)
kx2 = size(u2,2)

#ifdef debug
if(nmaxpr==1.and.mydiag)then
  write(6,*) 'unstag_ktau,nstag,nstagu',ktau,nstag,nstagu
endif
#endif

do k = 1,kx1
  uin1(1:ifull,k) = u1(1:ifull,k)
  vin1(1:ifull,k) = v1(1:ifull,k)
end do
do k = 1,kx2
  uin2(1:ifull,k) = u2(1:ifull,k)
  vin2(1:ifull,k) = v2(1:ifull,k)
end do

if (abs(nstagu)<3) then
  call boundsuv(u1=uin1,v1=vin1,u2=uin2,v2=vin2,stag=3)
  do k = 1,kx1
    do iq = 1,ifull
      uout1(iq,k)=(9.*(uin1(iwu(iq),k)+uin1(iq,k))-uin1(iwwu(iq),k)-uin1(ieu(iq),k))/16.
      vout1(iq,k)=(9.*(vin1(isv(iq),k)+vin1(iq,k))-vin1(issv(iq),k)-vin1(inv(iq),k))/16.
    end do
  end do  
  do k = 1,kx2
    do iq = 1,ifull
      uout2(iq,k)=(9.*(uin2(iwu(iq),k)+uin2(iq,k))-uin2(iwwu(iq),k)-uin2(ieu(iq),k))/16.
      vout2(iq,k)=(9.*(vin2(isv(iq),k)+vin2(iq,k))-vin2(issv(iq),k)-vin2(inv(iq),k))/16.
    end do
  end do  
  return
endif  ! (nstagu==0)

if ( nstagu==3 ) then
  call boundsuv(u1=uin1,v1=vin1,u2=uin2,v2=vin2,stag=5) ! issv, isv, iwwu, iwu

  ! precalculate rhs terms with iwwu2 & issv2
  do k = 1,kx1
    do iq = 1,ifull  
      ud1(iq,k) = uin1(iq,k)/2.+uin1(iwu(iq),k)+uin1(iwwu(iq),k)/10.
      vd1(iq,k) = vin1(iq,k)/2.+vin1(isv(iq),k)+vin1(issv(iq),k)/10.
    end do
  end do  
  do k = 1,kx2
    do iq = 1,ifull  
      ud2(iq,k) = uin2(iq,k)/2.+uin2(iwu(iq),k)+uin2(iwwu(iq),k)/10.
      vd2(iq,k) = vin2(iq,k)/2.+vin2(isv(iq),k)+vin2(issv(iq),k)/10.
    end do
  end do  

  call boundsuv(u1=ud1,v1=vd1,u2=ud2,v2=vd2,stag=-9) ! isv, iwu
  
  do k = 1,kx1
    do iq = 1,ifull  
      ua1(iq,k) = ud1(iq,k)-ud1(iwu(iq),k)/2. ! 1st guess
      va1(iq,k) = vd1(iq,k)-vd1(isv(iq),k)/2. ! 1st guess
      ug1(iq,k) = ua1(iq,k)
      vg1(iq,k) = va1(iq,k)
    end do
  end do
  do k = 1,kx2
    do iq = 1,ifull  
      ua2(iq,k) = ud2(iq,k)-ud2(iwu(iq),k)/2. ! 1st guess
      va2(iq,k) = vd2(iq,k)-vd2(isv(iq),k)/2. ! 1st guess
      ug2(iq,k) = ua2(iq,k)
      vg2(iq,k) = va2(iq,k)
    end do
  end do

  do itn=1,itnmax        ! each loop is a double iteration
    call boundsuv(u1=ua1,v1=va1,u2=ua2,v2=va2,stag=3) ! issv, isv, inv, iwwu, iwu, ieu
    do k = 1,kx1
      do iq = 1,ifull  
        uin1(iq,k) = (ug1(iq,k)-ua1(ieu(iq),k)/10. +ua1(iwwu(iq),k)/4.)/.95
        vin1(iq,k) = (vg1(iq,k)-va1(inv(iq),k)/10. +va1(issv(iq),k)/4.)/.95
      end do
    end do
    do k = 1,kx2
      do iq = 1,ifull  
        uin2(iq,k) = (ug2(iq,k)-ua2(ieu(iq),k)/10. +ua2(iwwu(iq),k)/4.)/.95
        vin2(iq,k) = (vg2(iq,k)-va2(inv(iq),k)/10. +va2(issv(iq),k)/4.)/.95
      end do
    end do
    call boundsuv(u1=uin1,v1=vin1,u2=uin2,v2=vin2,stag=3) ! issv, isv, inv, iwwu, iwu, ieu
    do k = 1,kx1
      do iq = 1,ifull  
        ua1(iq,k) = (ug1(iq,k)-uin1(ieu(iq),k)/10. +uin1(iwwu(iq),k)/4.)/.95
        va1(iq,k) = (vg1(iq,k)-vin1(inv(iq),k)/10. +vin1(issv(iq),k)/4.)/.95
      end do
    end do
    do k = 1,kx2
      do iq = 1,ifull  
        ua2(iq,k) = (ug2(iq,k)-uin2(ieu(iq),k)/10. +uin2(iwwu(iq),k)/4.)/.95
        va2(iq,k) = (vg2(iq,k)-vin2(inv(iq),k)/10. +vin2(issv(iq),k)/4.)/.95
      end do
    end do
  end do                 ! itn=1,itnmax

else !if ( nstagu==4 ) then
  call boundsuv(u1=uin1,v1=vin1,u2=uin2,v2=vin2,stag=2) ! isv, inv, innv, iwu, ieu, ieeu

  do k = 1,kx1
    do iq = 1,ifull  
      ua1(iq,k) = -0.05*uin1(ieeu(iq),k)-0.4*uin1(ieu(iq),k)+0.75*uin1(iq,k)+0.5*uin1(iwu(iq),k) ! 1st guess
      va1(iq,k) = -0.05*vin1(innv(iq),k)-0.4*vin1(inv(iq),k)+0.75*vin1(iq,k)+0.5*vin1(isv(iq),k) ! 1st guess
      ug1(iq,k) = ua1(iq,k)
      vg1(iq,k) = va1(iq,k)
    end do
  end do
  do k = 1,kx2
    do iq = 1,ifull  
      ua2(iq,k) = -0.05*uin2(ieeu(iq),k)-0.4*uin2(ieu(iq),k)+0.75*uin2(iq,k)+0.5*uin2(iwu(iq),k) ! 1st guess
      va2(iq,k) = -0.05*vin2(innv(iq),k)-0.4*vin2(inv(iq),k)+0.75*vin2(iq,k)+0.5*vin2(isv(iq),k) ! 1st guess
      ug2(iq,k) = ua2(iq,k)
      vg2(iq,k) = va2(iq,k)
    end do
  end do

  do itn=1,itnmax        ! each loop is a double iteration
    call boundsuv(u1=ua1,v1=va1,u2=ua2,v2=va2,stag=2) ! isv, inv, innv, iwu, ieu, ieeu
    do k = 1,kx1
      do iq = 1,ifull  
        uin1(iq,k) = (ug1(iq,k)-ua1(iwu(iq),k)/10. +ua1(ieeu(iq),k)/4.)/.95
        vin1(iq,k) = (vg1(iq,k)-va1(isv(iq),k)/10. +va1(innv(iq),k)/4.)/.95
      end do
    end do  
    do k = 1,kx2
      do iq = 1,ifull  
        uin2(iq,k) = (ug2(iq,k)-ua2(iwu(iq),k)/10. +ua2(ieeu(iq),k)/4.)/.95
        vin2(iq,k) = (vg2(iq,k)-va2(isv(iq),k)/10. +va2(innv(iq),k)/4.)/.95
      end do
    end do  
    call boundsuv(u1=uin1,v1=vin1,u2=uin2,v2=vin2,stag=2) ! isv, inv, innv, iwu, ieu, ieeu
    do k = 1,kx1
      do iq = 1,ifull  
        ua1(iq,k) = (ug1(iq,k)-uin1(iwu(iq),k)/10. +uin1(ieeu(iq),k)/4.)/.95
        va1(iq,k) = (vg1(iq,k)-vin1(isv(iq),k)/10. +vin1(innv(iq),k)/4.)/.95
      end do
    end do
    do k = 1,kx2
      do iq = 1,ifull  
        ua2(iq,k) = (ug2(iq,k)-uin2(iwu(iq),k)/10. +uin2(ieeu(iq),k)/4.)/.95
        va2(iq,k) = (vg2(iq,k)-vin2(isv(iq),k)/10. +vin2(innv(iq),k)/4.)/.95
      end do
    end do
  enddo                  ! itn=1,itnmax
      
end if

do k = 1,kx1
  uout1(1:ifull,k) = ua1(1:ifull,k)
  vout1(1:ifull,k) = va1(1:ifull,k)
end do
do k = 1,kx2
  uout2(1:ifull,k) = ua2(1:ifull,k)
  vout2(1:ifull,k) = va2(1:ifull,k)
end do

call END_LOG(stag_end)

return
end subroutine unstaguv3x

end module staguvmod