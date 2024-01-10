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
    
module betts1_m

implicit none

private
public la,itb,jtb,klh,ktm
public thl,deta,aeta,eta,hbm2,htm,res,dtq2,dum3a
public pl,rdq,rdth,rdp,rdthe
public qs0,sqs,the0,sthe,ptbl,ttbl
public sm,pd,prec,cuprec,cldefi,t,q
public betts1_init,betts1_end

integer, save :: ktm
integer, parameter :: la = 25
integer, parameter :: itb = 76
integer, parameter :: jtb = 171
real, save :: dtq2,pl,rdq,rdth,rdp,rdthe
real, save :: thl = 210.

integer, dimension(:), allocatable, save :: klh
real, dimension(:), allocatable, save :: deta,aeta,eta
real, dimension(:), allocatable, save :: hbm2,sm,res
real, dimension(:), allocatable, save :: qs0,sqs,the0,sthe
real, dimension(:), allocatable, save :: pd,prec,cuprec,cldefi
real, dimension(:,:), allocatable, save :: htm
real, dimension(:,:), allocatable, save :: ptbl,ttbl
real, dimension(:,:), allocatable, save :: t,q,dum3a

contains

subroutine betts1_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(deta(kl),aeta(kl),eta(kl+1))
allocate(hbm2(ifull),htm(ifull,kl),sm(ifull),klh(ifull),res(ifull))
allocate(qs0(jtb),sqs(jtb),the0(itb),sthe(itb),ptbl(itb,jtb),ttbl(jtb,itb))
allocate(pd(ifull),prec(ifull),cuprec(ifull),cldefi(ifull))
allocate(t(ifull,kl),q(ifull,kl),dum3a(ifull,kl))

return
end subroutine betts1_init

subroutine betts1_end

implicit none

deallocate(deta,aeta,eta)
deallocate(hbm2,htm,sm,klh,res)
deallocate(qs0,sqs,the0,sthe,ptbl,ttbl)
deallocate(pd,prec,cuprec,cldefi)
deallocate(t,q,dum3a)

return
end subroutine betts1_end

end module betts1_m