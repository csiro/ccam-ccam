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
use newmpar_m
use parm_m
use parmdyn_m
use vecsuv_m

implicit none

real, dimension(:,:), intent(inout)  :: u, v ! in case u=uout and v=vout
real, dimension(:,:), intent(out) :: uout, vout
real, dimension(ifull+iextra,size(u,2)) :: ua, va, ud, vd, uin, vin
real, dimension(ifull,size(u,2)) :: ug, vg
integer :: itn, kx
#ifdef debug
integer, parameter :: ntest=0    ! usually 0, 1 for test prints
#endif

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
  uout(1:ifull,1:kx)=(9.*(uin(ieu,1:kx)+uin(1:ifull,1:kx))-uin(iwu,1:kx)-uin(ieeu,1:kx))/16.
  vout(1:ifull,1:kx)=(9.*(vin(inv,1:kx)+vin(1:ifull,1:kx))-vin(isv,1:kx)-vin(innv,1:kx))/16.
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
  ud(1:ifull,1:kx)=uin(1:ifull,1:kx)/2.+uin(ieu,1:kx)+uin(ieeu,1:kx)/10.
  vd(1:ifull,1:kx)=vin(1:ifull,1:kx)/2.+vin(inv,1:kx)+vin(innv,1:kx)/10.

  call boundsuv(ud,vd,stag=-10) ! inv, ieu
  ua(1:ifull,1:kx)=ud(1:ifull,1:kx)-ud(ieu,1:kx)/2. ! 1st guess
  va(1:ifull,1:kx)=vd(1:ifull,1:kx)-vd(inv,1:kx)/2. ! 1st guess
  ug(1:ifull,1:kx)=ua(1:ifull,1:kx)
  vg(1:ifull,1:kx)=va(1:ifull,1:kx)

  do itn=1,itnmax-1        ! each loop is a double iteration
    call boundsuv(ua,va,stag=2) ! isv, inv, innv, iwu, ieu, ieeu
    uin(1:ifull,1:kx)=(ug(1:ifull,1:kx)-ua(iwu,1:kx)/10. +ua(ieeu,1:kx)/4.)/.95
    vin(1:ifull,1:kx)=(vg(1:ifull,1:kx)-va(isv,1:kx)/10. +va(innv,1:kx)/4.)/.95

    call boundsuv(uin,vin,stag=2) ! isv, inv, innv, iwu, ieu, ieeu
    ua(1:ifull,1:kx)=(ug(1:ifull,1:kx)-uin(iwu,1:kx)/10. +uin(ieeu,1:kx)/4.)/.95
    va(1:ifull,1:kx)=(vg(1:ifull,1:kx)-vin(isv,1:kx)/10. +vin(innv,1:kx)/4.)/.95
  end do                  ! itn=1,itnmax
  call boundsuv(ua,va,stag=2) ! isv, inv, innv, iwu, ieu, ieeu
  uin(1:ifull,1:kx)=(ug(1:ifull,1:kx)-ua(iwu,1:kx)/10. +ua(ieeu,1:kx)/4.)/.95
  vin(1:ifull,1:kx)=(vg(1:ifull,1:kx)-va(isv,1:kx)/10. +va(innv,1:kx)/4.)/.95
  call boundsuv(uin,vin,stag=2) ! isv, inv, innv, iwu, ieu, ieeu
  uout(1:ifull,1:kx)=(ug(1:ifull,1:kx)-uin(iwu,1:kx)/10. +uin(ieeu,1:kx)/4.)/.95
  vout(1:ifull,1:kx)=(vg(1:ifull,1:kx)-vin(isv,1:kx)/10. +vin(innv,1:kx)/4.)/.95

else !if ( nstag==4 ) then
  call boundsuv(uin,vin,stag=3) ! issv, isv, inv, iwwu, iwu, ieu

  ua(1:ifull,1:kx)=-0.05*uin(iwwu,1:kx)-0.4*uin(iwu,1:kx)+0.75*uin(1:ifull,1:kx)+0.5*uin(ieu,1:kx) ! 1st guess
  va(1:ifull,1:kx)=-0.05*vin(issv,1:kx)-0.4*vin(isv,1:kx)+0.75*vin(1:ifull,1:kx)+0.5*vin(inv,1:kx) ! 1st guess
  ug(1:ifull,1:kx)=ua(1:ifull,1:kx)
  vg(1:ifull,1:kx)=va(1:ifull,1:kx)

  do itn=1,itnmax-1        ! each loop is a double iteration
    call boundsuv(ua,va,stag=3) ! issv, isv, inv, iwwu, iwu, ieu
    uin(1:ifull,1:kx)=(ug(1:ifull,1:kx)-ua(ieu,1:kx)/10. +ua(iwwu,1:kx)/4.)/.95
    vin(1:ifull,1:kx)=(vg(1:ifull,1:kx)-va(inv,1:kx)/10. +va(issv,1:kx)/4.)/.95

    call boundsuv(uin,vin,stag=3) ! issv, isv, inv, iwwu, iwu, ieu
    ua(1:ifull,1:kx)=(ug(1:ifull,1:kx)-uin(ieu,1:kx)/10. +uin(iwwu,1:kx)/4.)/.95
    va(1:ifull,1:kx)=(vg(1:ifull,1:kx)-vin(inv,1:kx)/10. +vin(issv,1:kx)/4.)/.95
  end do                 ! itn=1,itnmax
  call boundsuv(ua,va,stag=3) ! issv, isv, inv, iwwu, iwu, ieu
  uin(1:ifull,1:kx)=(ug(1:ifull,1:kx)-ua(ieu,1:kx)/10. +ua(iwwu,1:kx)/4.)/.95
  vin(1:ifull,1:kx)=(vg(1:ifull,1:kx)-va(inv,1:kx)/10. +va(issv,1:kx)/4.)/.95
  call boundsuv(uin,vin,stag=3) ! issv, isv, inv, iwwu, iwu, ieu
  uout(1:ifull,1:kx)=(ug(1:ifull,1:kx)-uin(ieu,1:kx)/10. +uin(iwwu,1:kx)/4.)/.95
  vout(1:ifull,1:kx)=(vg(1:ifull,1:kx)-vin(inv,1:kx)/10. +vin(issv,1:kx)/4.)/.95
 
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
integer :: itn, kx

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
  uout(1:ifull,1:kx)=(9.*(uin(iwu,1:kx)+uin(1:ifull,1:kx))-uin(iwwu,1:kx)-uin(ieu,1:kx))/16.
  vout(1:ifull,1:kx)=(9.*(vin(isv,1:kx)+vin(1:ifull,1:kx))-vin(issv,1:kx)-vin(inv,1:kx))/16.
  return
endif  ! (nstagu==0)

if ( nstagu==3 ) then
  call boundsuv(uin,vin,stag=5) ! issv, isv, iwwu, iwu
  ! precalculate rhs terms with iwwu2 & issv2
  ud(1:ifull,1:kx)=uin(1:ifull,1:kx)/2.+uin(iwu,1:kx)+uin(iwwu,1:kx)/10.
  vd(1:ifull,1:kx)=vin(1:ifull,1:kx)/2.+vin(isv,1:kx)+vin(issv,1:kx)/10.

  call boundsuv(ud,vd,stag=-9) ! isv, iwu
  ua(1:ifull,1:kx)=ud(1:ifull,1:kx)-ud(iwu,1:kx)/2. ! 1st guess
  va(1:ifull,1:kx)=vd(1:ifull,1:kx)-vd(isv,1:kx)/2. ! 1st guess
  ug(1:ifull,1:kx)=ua(1:ifull,1:kx)
  vg(1:ifull,1:kx)=va(1:ifull,1:kx)

  do itn=1,itnmax-1        ! each loop is a double iteration
    call boundsuv(ua,va,stag=3) ! issv, isv, inv, iwwu, iwu, ieu
    uin(1:ifull,1:kx)=(ug(1:ifull,1:kx)-ua(ieu,1:kx)/10. +ua(iwwu,1:kx)/4.)/.95
    vin(1:ifull,1:kx)=(vg(1:ifull,1:kx)-va(inv,1:kx)/10. +va(issv,1:kx)/4.)/.95

    call boundsuv(uin,vin,stag=3) ! issv, isv, inv, iwwu, iwu, ieu
    ua(1:ifull,1:kx)=(ug(1:ifull,1:kx)-uin(ieu,1:kx)/10. +uin(iwwu,1:kx)/4.)/.95
    va(1:ifull,1:kx)=(vg(1:ifull,1:kx)-vin(inv,1:kx)/10. +vin(issv,1:kx)/4.)/.95
  end do                 ! itn=1,itnmax
  call boundsuv(ua,va,stag=3) ! issv, isv, inv, iwwu, iwu, ieu
  uin(1:ifull,1:kx)=(ug(1:ifull,1:kx)-ua(ieu,1:kx)/10. +ua(iwwu,1:kx)/4.)/.95
  vin(1:ifull,1:kx)=(vg(1:ifull,1:kx)-va(inv,1:kx)/10. +va(issv,1:kx)/4.)/.95
  call boundsuv(uin,vin,stag=3) ! issv, isv, inv, iwwu, iwu, ieu
  uout(1:ifull,1:kx)=(ug(1:ifull,1:kx)-uin(ieu,1:kx)/10. +uin(iwwu,1:kx)/4.)/.95
  vout(1:ifull,1:kx)=(vg(1:ifull,1:kx)-vin(inv,1:kx)/10. +vin(issv,1:kx)/4.)/.95

else !if ( nstagu==4 ) then
  call boundsuv(uin,vin,stag=2) ! isv, inv, innv, iwu, ieu, ieeu

  ua(1:ifull,1:kx)=-0.05*uin(ieeu,1:kx)-0.4*uin(ieu,1:kx)+0.75*uin(1:ifull,1:kx)+0.5*uin(iwu,1:kx) ! 1st guess
  va(1:ifull,1:kx)=-0.05*vin(innv,1:kx)-0.4*vin(inv,1:kx)+0.75*vin(1:ifull,1:kx)+0.5*vin(isv,1:kx) ! 1st guess
  ug(1:ifull,1:kx)=ua(1:ifull,1:kx)
  vg(1:ifull,1:kx)=va(1:ifull,1:kx)

  do itn=1,itnmax-1        ! each loop is a double iteration
    call boundsuv(ua,va,stag=2) ! isv, inv, innv, iwu, ieu, ieeu
    uin(1:ifull,1:kx)=(ug(1:ifull,1:kx)-ua(iwu,1:kx)/10. +ua(ieeu,1:kx)/4.)/.95
    vin(1:ifull,1:kx)=(vg(1:ifull,1:kx)-va(isv,1:kx)/10. +va(innv,1:kx)/4.)/.95

    call boundsuv(uin,vin,stag=2) ! isv, inv, innv, iwu, ieu, ieeu
    ua(1:ifull,1:kx)=(ug(1:ifull,1:kx)-uin(iwu,1:kx)/10. +uin(ieeu,1:kx)/4.)/.95
    va(1:ifull,1:kx)=(vg(1:ifull,1:kx)-vin(isv,1:kx)/10. +vin(innv,1:kx)/4.)/.95
  enddo                  ! itn=1,itnmax
  call boundsuv(ua,va,stag=2) ! isv, inv, innv, iwu, ieu, ieeu
  uin(1:ifull,1:kx)=(ug(1:ifull,1:kx)-ua(iwu,1:kx)/10. +ua(ieeu,1:kx)/4.)/.95
  vin(1:ifull,1:kx)=(vg(1:ifull,1:kx)-va(isv,1:kx)/10. +va(innv,1:kx)/4.)/.95
  call boundsuv(uin,vin,stag=2) ! isv, inv, innv, iwu, ieu, ieeu
  uout(1:ifull,1:kx)=(ug(1:ifull,1:kx)-uin(iwu,1:kx)/10. +uin(ieeu,1:kx)/4.)/.95
  vout(1:ifull,1:kx)=(vg(1:ifull,1:kx)-vin(isv,1:kx)/10. +vin(innv,1:kx)/4.)/.95
      
end if

call END_LOG(stag_end)

return
end subroutine unstaguv

end module staguvmod