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
use vecsuv_m

implicit none

include 'newmpar.h'
include 'parm.h'
include 'parmdyn.h'

real, dimension(:,:), intent(inout)  :: u, v ! in case u=uout and v=vout
real, dimension(:,:), intent(out) :: uout, vout
real, dimension(ifull+iextra,kl) :: ua, va, ud, vd, uin, vin
real, dimension(ifull,kl) :: ug, vg
integer, parameter :: ntest=0    ! usually 0, 1 for test prints
integer :: itn

call START_LOG(stag_begin)

#ifdef debug
if(nmaxpr==1.and.mydiag)then
  write(6,*) '  stag_ktau,nstag,nstagu',ktau,nstag,nstagu
endif
#endif

! Copying could be avoided if input arrays were dimensioned ifull+iextra
uin(1:ifull,1:kl)=u(1:ifull,1:kl)
vin(1:ifull,1:kl)=v(1:ifull,1:kl)
      
if (abs(nstag)<3) then
  call boundsuv(uin,vin,stag=2)
  uout(1:ifull,1:kl)=(9.*(uin(ieu,1:kl)+uin(1:ifull,1:kl))-uin(iwu,1:kl)-uin(ieeu,1:kl))/16.
  vout(1:ifull,1:kl)=(9.*(vin(inv,1:kl)+vin(1:ifull,1:kl))-vin(isv,1:kl)-vin(innv,1:kl))/16.
  return
endif  ! (nstag==0)

if ( nstag==3 ) then
  call boundsuv(uin,vin,stag=1) ! inv, innv, ieu, ieeu
#ifdef debug  
  if(ntest==1)then
    write(6,*) 'staguv diags'
    write (6,"(2x,4i8,6x,4i8)") (i,i=1,4),(i,i=1,4)
    do j=93,96
      write (6,"(i4,4f8.3,6x,4f8.3)") j,(u(i+(j-1)*il,4),i=1,4),(v(i+(j-1)*il,4),i=1,4)
    enddo          
    write (6,"(2x,4i8,6x,4i8)") (i,i=1,4),(i,i=1,4)
    do j=189,192
      write (6,"(i4,4f8.3,6x,4f8.3)") j,(u(i+(j-1)*il,4),i=1,4),(v(i+(j-1)*il,4),i=1,4)
    enddo          
    write (6,"(2x,4i8,6x,4i8)") (i,i=1,4),(i,i=1,4)
    do j=285,288
      write (6,"(i4,4f8.3,6x,4f8.3)") j,(u(i+(j-1)*il,4),i=1,4),(v(i+(j-1)*il,4),i=1,4)
    enddo          
    do j=95,288,96
      do i=1,2
        iq=i+(j-1)*il
        write (6,"('i,j,uin(ieu),uin(ieeu) ',2i4,2f8.3)") i,j,uin(ieu(iq),4),uin(ieeu(iq),4) 
        write (6,"('i,j,uin(iwu),uin(iwwu) ',2i4,2f8.3)") i,j,uin(iwu(iq),4),uin(iwwu(iq),4) 
        write (6,"('i,j,vin(inv),vin(innv) ',2i4,2f8.3)") i,j,vin(inv(iq),4),vin(innv(iq),4) 
        write (6,"('i,j,vin(isv),vin(issv) ',2i4,2f8.3)") i,j,vin(isv(iq),4),vin(issv(iq),4) 
      enddo
    enddo
  endif  ! (ntest==1)
#endif
         
  ! precalculate rhs terms with iwwu2 & issv2
  ud(1:ifull,1:kl)=uin(1:ifull,1:kl)/2.+uin(ieu,1:kl)+uin(ieeu,1:kl)/10.
  vd(1:ifull,1:kl)=vin(1:ifull,1:kl)/2.+vin(inv,1:kl)+vin(innv,1:kl)/10.

  call boundsuv(ud,vd,stag=-10) ! inv, ieu
  ua(1:ifull,1:kl)=ud(1:ifull,1:kl)-ud(ieu,1:kl)/2. ! 1st guess
  va(1:ifull,1:kl)=vd(1:ifull,1:kl)-vd(inv,1:kl)/2. ! 1st guess
  ug(1:ifull,1:kl)=ua(1:ifull,1:kl)
  vg(1:ifull,1:kl)=va(1:ifull,1:kl)

  do itn=1,itnmax-1        ! each loop is a double iteration
    call boundsuv(ua,va,stag=2) ! isv, inv, innv, iwu, ieu, ieeu
    uin(1:ifull,1:kl)=(ug(1:ifull,1:kl)-ua(iwu,1:kl)/10. +ua(ieeu,1:kl)/4.)/.95
    vin(1:ifull,1:kl)=(vg(1:ifull,1:kl)-va(isv,1:kl)/10. +va(innv,1:kl)/4.)/.95

    call boundsuv(uin,vin,stag=2) ! isv, inv, innv, iwu, ieu, ieeu
    ua(1:ifull,1:kl)=(ug(1:ifull,1:kl)-uin(iwu,1:kl)/10. +uin(ieeu,1:kl)/4.)/.95
    va(1:ifull,1:kl)=(vg(1:ifull,1:kl)-vin(isv,1:kl)/10. +vin(innv,1:kl)/4.)/.95
  end do                  ! itn=1,itnmax
  call boundsuv(ua,va,stag=2) ! isv, inv, innv, iwu, ieu, ieeu
  uin(1:ifull,1:kl)=(ug(1:ifull,1:kl)-ua(iwu,1:kl)/10. +ua(ieeu,1:kl)/4.)/.95
  vin(1:ifull,1:kl)=(vg(1:ifull,1:kl)-va(isv,1:kl)/10. +va(innv,1:kl)/4.)/.95
  call boundsuv(uin,vin,stag=2) ! isv, inv, innv, iwu, ieu, ieeu
  uout(1:ifull,1:kl)=(ug(1:ifull,1:kl)-uin(iwu,1:kl)/10. +uin(ieeu,1:kl)/4.)/.95
  vout(1:ifull,1:kl)=(vg(1:ifull,1:kl)-vin(isv,1:kl)/10. +vin(innv,1:kl)/4.)/.95

else !if ( nstag==4 ) then
  call boundsuv(uin,vin,stag=3) ! issv, isv, inv, iwwu, iwu, ieu

  ua(1:ifull,1:kl)=-0.05*uin(iwwu,1:kl)-0.4*uin(iwu,1:kl)+0.75*uin(1:ifull,1:kl)+0.5*uin(ieu,1:kl) ! 1st guess
  va(1:ifull,1:kl)=-0.05*vin(issv,1:kl)-0.4*vin(isv,1:kl)+0.75*vin(1:ifull,1:kl)+0.5*vin(inv,1:kl) ! 1st guess
  ug(1:ifull,1:kl)=ua(1:ifull,1:kl)
  vg(1:ifull,1:kl)=va(1:ifull,1:kl)

  do itn=1,itnmax-1        ! each loop is a double iteration
    call boundsuv(ua,va,stag=3) ! issv, isv, inv, iwwu, iwu, ieu
    uin(1:ifull,1:kl)=(ug(1:ifull,1:kl)-ua(ieu,1:kl)/10. +ua(iwwu,1:kl)/4.)/.95
    vin(1:ifull,1:kl)=(vg(1:ifull,1:kl)-va(inv,1:kl)/10. +va(issv,1:kl)/4.)/.95

    call boundsuv(uin,vin,stag=3) ! issv, isv, inv, iwwu, iwu, ieu
    ua(1:ifull,1:kl)=(ug(1:ifull,1:kl)-uin(ieu,1:kl)/10. +uin(iwwu,1:kl)/4.)/.95
    va(1:ifull,1:kl)=(vg(1:ifull,1:kl)-vin(inv,1:kl)/10. +vin(issv,1:kl)/4.)/.95
  end do                 ! itn=1,itnmax
  call boundsuv(ua,va,stag=3) ! issv, isv, inv, iwwu, iwu, ieu
  uin(1:ifull,1:kl)=(ug(1:ifull,1:kl)-ua(ieu,1:kl)/10. +ua(iwwu,1:kl)/4.)/.95
  vin(1:ifull,1:kl)=(vg(1:ifull,1:kl)-va(inv,1:kl)/10. +va(issv,1:kl)/4.)/.95
  call boundsuv(uin,vin,stag=3) ! issv, isv, inv, iwwu, iwu, ieu
  uout(1:ifull,1:kl)=(ug(1:ifull,1:kl)-uin(ieu,1:kl)/10. +uin(iwwu,1:kl)/4.)/.95
  vout(1:ifull,1:kl)=(vg(1:ifull,1:kl)-vin(inv,1:kl)/10. +vin(issv,1:kl)/4.)/.95
 
end if

call END_LOG(stag_end)
return
end subroutine staguv


subroutine unstaguv(u,v,uout,vout)

!     staggered u & v as input; unstaggered as output

use cc_mpi
use indices_m
use map_m
use vecsuv_m

implicit none

include 'newmpar.h'
include 'parm.h'
include 'parmdyn.h'

real, dimension(:,:), intent(inout)  :: u, v ! in case u=uout and v=vout
real, dimension(:,:), intent(out) :: uout, vout
real, dimension(ifull+iextra,kl) :: ua, va, ud, vd, uin, vin
real, dimension(ifull,kl) :: ug, vg
integer :: itn

call START_LOG(stag_begin)

#ifdef debug
if(nmaxpr==1.and.mydiag)then
  write(6,*) 'unstag_ktau,nstag,nstagu',ktau,nstag,nstagu
endif
#endif

uin(1:ifull,1:kl) = u(1:ifull,1:kl)
vin(1:ifull,1:kl) = v(1:ifull,1:kl)
      
if (abs(nstagu)<3) then
  call boundsuv(uin,vin,stag=3)
  uout(1:ifull,1:kl)=(9.*(uin(iwu,1:kl)+uin(1:ifull,1:kl))-uin(iwwu,1:kl)-uin(ieu,1:kl))/16.
  vout(1:ifull,1:kl)=(9.*(vin(isv,1:kl)+vin(1:ifull,1:kl))-vin(issv,1:kl)-vin(inv,1:kl))/16.
! uout(1:ifull,1:kl)=.5*(uin(iwu,1:kl)+uin(1:ifull,1:kl))  ! for linear tests
! vout(1:ifull,1:kl)=.5*(vin(isv,1:kl)+vin(1:ifull,1:kl))  ! for linear tests
  return
endif  ! (nstagu==0)

if ( nstagu==3 ) then
  call boundsuv(uin,vin,stag=5) ! issv, isv, iwwu, iwu
  ! precalculate rhs terms with iwwu2 & issv2
  ud(1:ifull,1:kl)=uin(1:ifull,1:kl)/2.+uin(iwu,1:kl)+uin(iwwu,1:kl)/10.
  vd(1:ifull,1:kl)=vin(1:ifull,1:kl)/2.+vin(isv,1:kl)+vin(issv,1:kl)/10.

  call boundsuv(ud,vd,stag=-9) ! isv, iwu
  ua(1:ifull,1:kl)=ud(1:ifull,1:kl)-ud(iwu,1:kl)/2. ! 1st guess
  va(1:ifull,1:kl)=vd(1:ifull,1:kl)-vd(isv,1:kl)/2. ! 1st guess
  ug(1:ifull,1:kl)=ua(1:ifull,1:kl)
  vg(1:ifull,1:kl)=va(1:ifull,1:kl)

  do itn=1,itnmax-1        ! each loop is a double iteration
    call boundsuv(ua,va,stag=3) ! issv, isv, inv, iwwu, iwu, ieu
    uin(1:ifull,1:kl)=(ug(1:ifull,1:kl)-ua(ieu,1:kl)/10. +ua(iwwu,1:kl)/4.)/.95
    vin(1:ifull,1:kl)=(vg(1:ifull,1:kl)-va(inv,1:kl)/10. +va(issv,1:kl)/4.)/.95

    call boundsuv(uin,vin,stag=3) ! issv, isv, inv, iwwu, iwu, ieu
    ua(1:ifull,1:kl)=(ug(1:ifull,1:kl)-uin(ieu,1:kl)/10. +uin(iwwu,1:kl)/4.)/.95
    va(1:ifull,1:kl)=(vg(1:ifull,1:kl)-vin(inv,1:kl)/10. +vin(issv,1:kl)/4.)/.95
  end do                 ! itn=1,itnmax
  call boundsuv(ua,va,stag=3) ! issv, isv, inv, iwwu, iwu, ieu
  uin(1:ifull,1:kl)=(ug(1:ifull,1:kl)-ua(ieu,1:kl)/10. +ua(iwwu,1:kl)/4.)/.95
  vin(1:ifull,1:kl)=(vg(1:ifull,1:kl)-va(inv,1:kl)/10. +va(issv,1:kl)/4.)/.95
  call boundsuv(uin,vin,stag=3) ! issv, isv, inv, iwwu, iwu, ieu
  uout(1:ifull,1:kl)=(ug(1:ifull,1:kl)-uin(ieu,1:kl)/10. +uin(iwwu,1:kl)/4.)/.95
  vout(1:ifull,1:kl)=(vg(1:ifull,1:kl)-vin(inv,1:kl)/10. +vin(issv,1:kl)/4.)/.95

else !if ( nstagu==4 ) then
  call boundsuv(uin,vin,stag=2) ! isv, inv, innv, iwu, ieu, ieeu

  ua(1:ifull,1:kl)=-0.05*uin(ieeu,1:kl)-0.4*uin(ieu,1:kl)+0.75*uin(1:ifull,1:kl)+0.5*uin(iwu,1:kl) ! 1st guess
  va(1:ifull,1:kl)=-0.05*vin(innv,1:kl)-0.4*vin(inv,1:kl)+0.75*vin(1:ifull,1:kl)+0.5*vin(isv,1:kl) ! 1st guess
  ug(1:ifull,1:kl)=ua(1:ifull,1:kl)
  vg(1:ifull,1:kl)=va(1:ifull,1:kl)

  do itn=1,itnmax-1        ! each loop is a double iteration
    call boundsuv(ua,va,stag=2) ! isv, inv, innv, iwu, ieu, ieeu
    uin(1:ifull,1:kl)=(ug(1:ifull,1:kl)-ua(iwu,1:kl)/10. +ua(ieeu,1:kl)/4.)/.95
    vin(1:ifull,1:kl)=(vg(1:ifull,1:kl)-va(isv,1:kl)/10. +va(innv,1:kl)/4.)/.95

    call boundsuv(uin,vin,stag=2) ! isv, inv, innv, iwu, ieu, ieeu
    ua(1:ifull,1:kl)=(ug(1:ifull,1:kl)-uin(iwu,1:kl)/10. +uin(ieeu,1:kl)/4.)/.95
    va(1:ifull,1:kl)=(vg(1:ifull,1:kl)-vin(isv,1:kl)/10. +vin(innv,1:kl)/4.)/.95
  enddo                  ! itn=1,itnmax
  call boundsuv(ua,va,stag=2) ! isv, inv, innv, iwu, ieu, ieeu
  uin(1:ifull,1:kl)=(ug(1:ifull,1:kl)-ua(iwu,1:kl)/10. +ua(ieeu,1:kl)/4.)/.95
  vin(1:ifull,1:kl)=(vg(1:ifull,1:kl)-va(isv,1:kl)/10. +va(innv,1:kl)/4.)/.95
  call boundsuv(uin,vin,stag=2) ! isv, inv, innv, iwu, ieu, ieeu
  uout(1:ifull,1:kl)=(ug(1:ifull,1:kl)-uin(iwu,1:kl)/10. +uin(ieeu,1:kl)/4.)/.95
  vout(1:ifull,1:kl)=(vg(1:ifull,1:kl)-vin(isv,1:kl)/10. +vin(innv,1:kl)/4.)/.95
      
end if

call END_LOG(stag_end)

return
end subroutine unstaguv

end module staguvmod