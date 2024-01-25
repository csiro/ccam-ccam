! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2024 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
module mlostag

implicit none

private
public mlostaguv, mlounstaguv
public mstagf, koff, nstagoffmlo

integer, save      :: nstagoffmlo = 0       ! staggering offset
integer, save      :: mstagf      = 30      ! alternating staggering (0=off left, -1=off right, >0 alternating)
integer, parameter :: koff        = 1       ! time split stagger relative to A-grid (koff=0) or C-grid (koff=1)
integer, parameter :: itnmax      = 3       ! number of interations for reversible staggering

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Stagger u and v
! Modified to include zero velocity next to land points
subroutine mlostaguv(u,v,uout,vout)

use cc_mpi
use indices_m
use mlo
use mlodynamicsarrays_m
use newmpar_m
use parm_m

implicit none

integer i,k,itn,kx,kn,iq
real, dimension(:,:), intent(in) :: u, v
real, dimension(:,:), intent(out) :: uout,vout
real, dimension(ifull+iextra,size(u,2)) :: uin,vin
real, dimension(ifull+iextra,size(u,2)) :: ua,va
real, dimension(ifull+iextra,size(u,2)) :: ud,vd
real, dimension(:,:,:), allocatable, save :: wtul,wtvl
real, dimension(:,:,:), allocatable, save :: wtur,wtvr
real, dimension(:,:,:), allocatable, save :: dtul,dtvl
real, dimension(:,:,:), allocatable, save :: dtur,dtvr
real, dimension(:,:), allocatable :: nwtu,nwtv
real, dimension(:,:), allocatable, save :: stul,stvl
real, dimension(:,:), allocatable, save :: stur,stvr
logical euetest,euwtest,evntest,evstest
logical euewtest,evnstest,eutest,evtest
logical ltest

call START_LOG(ocnstag_begin)

if (.not.allocated(wtul)) then

  allocate(wtul(ifull+iextra,wlev,0:3),wtvl(ifull+iextra,wlev,0:3))
  allocate(wtur(ifull+iextra,wlev,0:3),wtvr(ifull+iextra,wlev,0:3))
  allocate(dtul(ifull,wlev,3),dtvl(ifull,wlev,3))
  allocate(dtur(ifull,wlev,3),dtvr(ifull,wlev,3))
  allocate(stul(ifull,wlev),stvl(ifull,wlev))
  allocate(stur(ifull,wlev),stvr(ifull,wlev))

  allocate( nwtu(ifull,0:3), nwtv(ifull,0:3) )
  
  ! assign land arrays
  do k = 1,wlev
    do iq = 1,ifull
      eutest = ee(iq,k)*ee(ie(iq),k)>0.5
      evtest = ee(iq,k)*ee(in(iq),k)>0.5
      euetest = eutest .and. ee(ie(iq),k)*ee(iee(iq),k)>0.5
      euwtest = ee(iw(iq),k)*ee(iq,k)>0.5 .and. eutest
      evntest = evtest .and. ee(in(iq),k)*ee(inn(iq),k)>0.5
      evstest = ee(is(iq),k)*ee(iq,k)>0.5 .and. evtest
      euewtest = euetest .and. euwtest
      evnstest = evntest .and. evstest

      ! assign weights (left)

      ! |   *   | X E   |  EE  |     unstaggered
      ! W       * X     E            staggered
      if ( euewtest ) then
        wtul(iq,k,0)=1.
        wtul(iq,k,1)=-0.5
        wtul(iq,k,2)=-0.1
        wtul(iq,k,3)=0.
        dtul(iq,k,1)=0.1
        dtul(iq,k,2)=1.
        dtul(iq,k,3)=0.5

      ! #   *   | X E   |  EE  |     unstaggered
      ! 0       * X     E            staggered
      else if ( euetest ) then
        wtul(iq,k,0)=1.
        wtul(iq,k,1)=-0.5
        wtul(iq,k,2)=0.
        wtul(iq,k,3)=0.
        dtul(iq,k,1)=0.1
        dtul(iq,k,2)=1.
        dtul(iq,k,3)=0.5

      ! |   *   | X E   #  ##  #     unstaggered
      !         * X     0  ##  #     staggered
      else if ( euwtest ) then
        wtul(iq,k,0)=1.
        wtul(iq,k,1)=0.
        wtul(iq,k,2)=0.
        wtul(iq,k,3)=0.
        dtul(iq,k,1)=0.
        dtul(iq,k,2)=1.
        dtul(iq,k,3)=1./3.

      ! #   *   |   E   #  ##  #     unstaggered
      ! #       *       #  ##  #     staggered
      else if ( eutest ) then
        wtul(iq,k,0)=1.
        wtul(iq,k,1)=0.
        wtul(iq,k,2)=0.
        wtul(iq,k,3)=0.
        dtul(iq,k,1)=0.
        dtul(iq,k,2)=0.5
        dtul(iq,k,3)=0.5

      else
        wtul(iq,k,0)=0.
        wtul(iq,k,1)=0.
        wtul(iq,k,2)=0.
        wtul(iq,k,3)=0.
        dtul(iq,k,1)=0.
        dtul(iq,k,2)=0.
        dtul(iq,k,3)=0.            
      end if

      if ( evnstest ) then
        wtvl(iq,k,0)=1.
        wtvl(iq,k,1)=-0.5
        wtvl(iq,k,2)=-0.1
        wtvl(iq,k,3)=0.
        dtvl(iq,k,1)=0.1
        dtvl(iq,k,2)=1.
        dtvl(iq,k,3)=0.5
      else if ( evntest ) then
        wtvl(iq,k,0)=1.
        wtvl(iq,k,1)=-0.5
        wtvl(iq,k,2)=0.
        wtvl(iq,k,3)=0.
        dtvl(iq,k,1)=0.1
        dtvl(iq,k,2)=1.
        dtvl(iq,k,3)=0.5
      else if ( evstest ) then
        wtvl(iq,k,0)=1.
        wtvl(iq,k,1)=0.
        wtvl(iq,k,2)=0.
        wtvl(iq,k,3)=0.
        dtvl(iq,k,1)=0. 
        dtvl(iq,k,2)=1.
        dtvl(iq,k,3)=1./3.
      else if ( evtest ) then
        wtvl(iq,k,0)=1.
        wtvl(iq,k,1)=0.
        wtvl(iq,k,2)=0.
        wtvl(iq,k,3)=0.
        dtvl(iq,k,1)=0.
        dtvl(iq,k,2)=0.5
        dtvl(iq,k,3)=0.5
      else
        wtvl(iq,k,0)=0.
        wtvl(iq,k,1)=0.
        wtvl(iq,k,2)=0.
        wtvl(iq,k,3)=0.
        dtvl(iq,k,1)=0.
        dtvl(iq,k,2)=0.
        dtvl(iq,k,3)=0.
      end if
    end do
  end do
    
  ! Apply JLM's preconditioner
  do i = 0,3
    call boundsuv(wtul(:,:,i),wtvl(:,:,i),stag=-10)
  end do  
  do k = 1,wlev
    where (abs(wtul(ieu,k,0))>1.E-4.and.abs(wtul(1:ifull,k,1))>1.E-4)
      stul(1:ifull,k)=-wtul(1:ifull,k,1)/wtul(ieu,k,0)
      nwtu(1:ifull,0)=wtul(1:ifull,k,0)+stul(:,k)*wtul(ieu,k,2)
      nwtu(1:ifull,1)=wtul(1:ifull,k,1)+stul(:,k)*wtul(ieu,k,0)
      nwtu(1:ifull,2)=wtul(1:ifull,k,2)
      nwtu(1:ifull,3)=wtul(1:ifull,k,3)-stul(:,k)*wtul(ieu,k,1)
    elsewhere
      stul(1:ifull,k)=0.
      nwtu(1:ifull,0)=wtul(1:ifull,k,0)
      nwtu(1:ifull,1)=wtul(1:ifull,k,1)
      nwtu(1:ifull,2)=wtul(1:ifull,k,2)
      nwtu(1:ifull,3)=wtul(1:ifull,k,3)
    end where
    where (abs(wtvl(inv,k,0))>1.E-4.and.abs(wtvl(1:ifull,k,1))>1.E-4)
      stvl(1:ifull,k)=-wtvl(1:ifull,k,1)/wtvl(inv,k,0)
      nwtv(1:ifull,0)=wtvl(1:ifull,k,0)+stvl(:,k)*wtvl(inv,k,2)
      nwtv(1:ifull,1)=wtvl(1:ifull,k,1)+stvl(:,k)*wtvl(inv,k,0)
      nwtv(1:ifull,2)=wtvl(1:ifull,k,2)
      nwtv(1:ifull,3)=wtvl(1:ifull,k,3)-stvl(:,k)*wtvl(inv,k,1)
    elsewhere
      stvl(1:ifull,k)=0.
      nwtv(1:ifull,0)=wtvl(1:ifull,k,0)
      nwtv(1:ifull,1)=wtvl(1:ifull,k,1)
      nwtv(1:ifull,2)=wtvl(1:ifull,k,2)
      nwtv(1:ifull,3)=wtvl(1:ifull,k,3)
    end where
    where (abs(wtul(1:ifull,k,0))<1.E-4)
      wtul(1:ifull,k,0)=1.
      wtul(1:ifull,k,1)=0.
      wtul(1:ifull,k,2)=0.
      wtul(1:ifull,k,3)=0.
    elsewhere
      wtul(1:ifull,k,0)=nwtu(1:ifull,0)
      wtul(1:ifull,k,1)=nwtu(1:ifull,1)
      wtul(1:ifull,k,2)=nwtu(1:ifull,2)
      wtul(1:ifull,k,3)=nwtu(1:ifull,3)
    end where
    where (abs(wtvl(1:ifull,k,0))<1.E-4)
      wtvl(1:ifull,k,0)=1.
      wtvl(1:ifull,k,1)=0.
      wtvl(1:ifull,k,2)=0.
      wtvl(1:ifull,k,3)=0.
    elsewhere
      wtvl(1:ifull,k,0)=nwtv(1:ifull,0)
      wtvl(1:ifull,k,1)=nwtv(1:ifull,1)
      wtvl(1:ifull,k,2)=nwtv(1:ifull,2)
      wtvl(1:ifull,k,3)=nwtv(1:ifull,3)
    end where
  end do
    
  ! normalise
  do i = 0,3
    call boundsuv(wtul(:,:,i),wtvl(:,:,i),stag=-10)
  end do  
  do k = 1,wlev
    do i=1,3
      dtul(1:ifull,k,i)=dtul(1:ifull,k,i)/wtul(1:ifull,k,0)
      dtvl(1:ifull,k,i)=dtvl(1:ifull,k,i)/wtvl(1:ifull,k,0)
      wtul(1:ifull,k,i)=wtul(1:ifull,k,i)/wtul(1:ifull,k,0)
      wtvl(1:ifull,k,i)=wtvl(1:ifull,k,i)/wtvl(1:ifull,k,0)
    end do
    stul(1:ifull,k)=stul(1:ifull,k)*wtul(ieu,k,0)/wtul(1:ifull,k,0)
    stvl(1:ifull,k)=stvl(1:ifull,k)*wtvl(inv,k,0)/wtvl(1:ifull,k,0)
    wtul(1:ifull,k,0)=1.
    wtvl(1:ifull,k,0)=1.
  end do  

  ! assign weights (right)
  do k = 1,wlev
    do iq = 1,ifull
      eutest = ee(iq,k)*ee(ie(iq),k)>0.5
      evtest = ee(iq,k)*ee(in(iq),k)>0.5
      euetest = eutest     .and. ee(ie(iq),k)*ee(iee(iq),k)>0.5
      euwtest = ee(iw(iq),k)*ee(iq,k)>0.5 .and. eutest
      evntest = evtest     .and. ee(in(iq),k)*ee(inn(iq),k)>0.5
      evstest = ee(is(iq),k)*ee(iq,k)>0.5 .and. evtest
      euewtest = euetest .and. euwtest
      evnstest = evntest .and. evstest
 
      ! |   W   |   * X |  E   |     unstaggered
      !         W     X *      E     staggered
      if ( euewtest ) then
        wtur(iq,k,0)=1.
        wtur(iq,k,1)=-0.1
        wtur(iq,k,2)=-0.5
        wtur(iq,k,3)=0.
        dtur(iq,k,1)=0.1
        dtur(iq,k,2)=1.
        dtur(iq,k,3)=0.5

      ! |   W   |   * X |  E   #     unstaggered
      !         W     X *      0     staggered
      else if ( euwtest ) then
        wtur(iq,k,0)=1.
        wtur(iq,k,1)=0.
        wtur(iq,k,2)=-0.5
        wtur(iq,k,3)=0.
        dtur(iq,k,1)=0.1
        dtur(iq,k,2)=1.
        dtur(iq,k,3)=0.5

      ! #  ##   #   * X |   E   |     unstaggered
      ! #  ##   0     X *             staggered
      else if ( euetest ) then
        wtur(iq,k,0)=1.
        wtur(iq,k,1)=0.
        wtur(iq,k,2)=0.
        wtur(iq,k,3)=0.
        dtur(iq,k,1)=0.
        dtur(iq,k,2)=1.
        dtur(iq,k,3)=1./3.

      ! #  ##   #   *   |  E   #     unstaggered
      ! #  ##   #       *      #     staggered
      else if ( eutest ) then
        wtur(iq,k,0)=1.
        wtur(iq,k,1)=0.
        wtur(iq,k,2)=0.
        wtur(iq,k,3)=0.
        dtur(iq,k,1)=0.
        dtur(iq,k,2)=0.5
        dtur(iq,k,3)=0.5

      else
        wtur(iq,k,0)=0.
        wtur(iq,k,1)=0.
        wtur(iq,k,2)=0.
        wtur(iq,k,3)=0.
        dtur(iq,k,1)=0.
        dtur(iq,k,2)=0.
        dtur(iq,k,3)=0.
      end if

      if ( evnstest ) then
        wtvr(iq,k,0)=1.
        wtvr(iq,k,1)=-0.1
        wtvr(iq,k,2)=-0.5
        wtvr(iq,k,3)=0.
        dtvr(iq,k,1)=0.1
        dtvr(iq,k,2)=1.
        dtvr(iq,k,3)=0.5
      else if ( evstest ) then
        wtvr(iq,k,0)=1.
        wtvr(iq,k,1)=0.
        wtvr(iq,k,2)=-0.5
        wtvr(iq,k,3)=0.
        dtvr(iq,k,1)=0.1
        dtvr(iq,k,2)=1.
        dtvr(iq,k,3)=0.5
      else if ( evntest ) then
        wtvr(iq,k,0)=1.
        wtvr(iq,k,1)=0.
        wtvr(iq,k,2)=0.
        wtvr(iq,k,3)=0.
        dtvr(iq,k,1)=0.
        dtvr(iq,k,2)=1.
        dtvr(iq,k,3)=1./3.
      else if ( evtest ) then
        wtvr(iq,k,0)=1.
        wtvr(iq,k,1)=0.
        wtvr(iq,k,2)=0.
        wtvr(iq,k,3)=0.
        dtvr(iq,k,1)=0.
        dtvr(iq,k,2)=0.5
        dtvr(iq,k,3)=0.5
      else
        wtvr(iq,k,0)=0.
        wtvr(iq,k,1)=0.
        wtvr(iq,k,2)=0.
        wtvr(iq,k,3)=0.
        dtvr(iq,k,1)=0.
        dtvr(iq,k,2)=0.
        dtvr(iq,k,3)=0.
      end if
    end do
  end do  

  ! Apply JLM's preconditioner
  do i = 0,3
    call boundsuv(wtur(:,:,i),wtvr(:,:,i),stag=-9)
  end do  
  do k = 1,wlev
    where (abs(wtur(iwu,k,0))>1.E-4.and.abs(wtur(1:ifull,k,2))>1.E-4)
      stur(1:ifull,k)=-wtur(1:ifull,k,2)/wtur(iwu,k,0)
      nwtu(1:ifull,0)=wtur(1:ifull,k,0)+stur(:,k)*wtur(iwu,k,1)
      nwtu(1:ifull,1)=wtur(1:ifull,k,1)
      nwtu(1:ifull,2)=wtur(1:ifull,k,2)+stur(:,k)*wtur(iwu,k,0)
      nwtu(1:ifull,3)=wtur(1:ifull,k,3)-stur(:,k)*wtur(iwu,k,2)
    elsewhere
      stur(1:ifull,k)=0.
      nwtu(1:ifull,0)=wtur(1:ifull,k,0)
      nwtu(1:ifull,1)=wtur(1:ifull,k,1)
      nwtu(1:ifull,2)=wtur(1:ifull,k,2)
      nwtu(1:ifull,3)=wtur(1:ifull,k,3)
    end where
    where (abs(wtvr(isv,k,0))>1.E-4.and.abs(wtvr(1:ifull,k,2))>1.E-4)
      stvr(1:ifull,k)=-wtvr(1:ifull,k,2)/wtvr(isv,k,0)
      nwtv(1:ifull,0)=wtvr(1:ifull,k,0)+stvr(:,k)*wtvr(isv,k,1)
      nwtv(1:ifull,1)=wtvr(1:ifull,k,1)
      nwtv(1:ifull,2)=wtvr(1:ifull,k,2)+stvr(:,k)*wtvr(isv,k,0)
      nwtv(1:ifull,3)=wtvr(1:ifull,k,3)-stvr(:,k)*wtvr(isv,k,2)
    elsewhere
      stvr(1:ifull,k)=0.
      nwtv(1:ifull,0)=wtvr(1:ifull,k,0)
      nwtv(1:ifull,1)=wtvr(1:ifull,k,1)
      nwtv(1:ifull,2)=wtvr(1:ifull,k,2)
      nwtv(1:ifull,3)=wtvr(1:ifull,k,3)
    end where
    where (abs(wtur(1:ifull,k,0))<1.E-4)
      wtur(1:ifull,k,0)=1.
      wtur(1:ifull,k,1)=0.
      wtur(1:ifull,k,2)=0.
      wtur(1:ifull,k,3)=0.
    elsewhere
      wtur(1:ifull,k,0)=nwtu(1:ifull,0)
      wtur(1:ifull,k,1)=nwtu(1:ifull,1)
      wtur(1:ifull,k,2)=nwtu(1:ifull,2)
      wtur(1:ifull,k,3)=nwtu(1:ifull,3)
    end where
    where (abs(wtvr(1:ifull,k,0))<1.E-4)
      wtvr(1:ifull,k,0)=1.
      wtvr(1:ifull,k,1)=0.
      wtvr(1:ifull,k,2)=0.
      wtvr(1:ifull,k,3)=0.
    elsewhere
      wtvr(1:ifull,k,0)=nwtv(1:ifull,0)
      wtvr(1:ifull,k,1)=nwtv(1:ifull,1)
      wtvr(1:ifull,k,2)=nwtv(1:ifull,2)
      wtvr(1:ifull,k,3)=nwtv(1:ifull,3)
    end where
  end do  

  ! normalise
  do i = 0,3
    call boundsuv(wtur(:,:,i),wtvr(:,:,i),stag=-9)
  end do  
  do k = 1,wlev
    do i = 1,3
      dtur(1:ifull,k,i)=dtur(1:ifull,k,i)/wtur(1:ifull,k,0)
      dtvr(1:ifull,k,i)=dtvr(1:ifull,k,i)/wtvr(1:ifull,k,0)
      wtur(1:ifull,k,i)=wtur(1:ifull,k,i)/wtur(1:ifull,k,0)
      wtvr(1:ifull,k,i)=wtvr(1:ifull,k,i)/wtvr(1:ifull,k,0)
    end do
    stur(1:ifull,k)=stur(1:ifull,k)*wtur(iwu,k,0)/wtur(1:ifull,k,0)
    stvr(1:ifull,k)=stvr(1:ifull,k)*wtvr(isv,k,0)/wtvr(1:ifull,k,0)
    wtur(1:ifull,k,0)=1.
    wtvr(1:ifull,k,0)=1.
  end do  

  deallocate( nwtu, nwtv )

end if

kx = size(u,2)
kn = min(kx,wlev)

do k = 1,kn
  uin(1:ifull,k)=u(1:ifull,k)*ee(1:ifull,k)
  vin(1:ifull,k)=v(1:ifull,k)*ee(1:ifull,k)
end do
do k = kn+1,kx
  uin(1:ifull,k)=u(1:ifull,k)*ee(1:ifull,1)
  vin(1:ifull,k)=v(1:ifull,k)*ee(1:ifull,1)
end do

if ( mstagf==0 ) then
  ltest = .true.
else if ( mstagf<0 ) then
  ltest = .false.
else
  ! using ktau-1 ensures that the staggering is relative to the C grid
  ltest = mod(ktau-koff-nstagoffmlo,2*mstagf)<mstagf
end if

if ( ltest ) then

  call boundsuv(uin,vin,stag=1)
  !$omp parallel
  !$omp do schedule(static) private(k,iq)
  do k = 1,kn
    do iq = 1,ifull  
      ud(iq,k)=dtul(iq,k,1)*uin(ieeu(iq),k)+dtul(iq,k,2)*uin(ieu(iq),k)+dtul(iq,k,3)*uin(iq,k)
      vd(iq,k)=dtvl(iq,k,1)*vin(innv(iq),k)+dtvl(iq,k,2)*vin(inv(iq),k)+dtvl(iq,k,3)*vin(iq,k)
    end do  
  end do
  !$omp end do nowait
  !$omp do schedule(static) private(k,iq)
  do k = kn+1,kx
    do iq = 1,ifull  
      ud(iq,k)=dtul(iq,1,1)*uin(ieeu(iq),k)+dtul(iq,1,2)*uin(ieu(iq),k)+dtul(iq,1,3)*uin(iq,k)
      vd(iq,k)=dtvl(iq,1,1)*vin(innv(iq),k)+dtvl(iq,1,2)*vin(inv(iq),k)+dtvl(iq,1,3)*vin(iq,k)
    end do  
  end do
  !$omp end do

  call boundsuv(ud,vd,stag=-10)
  !$omp barrier
  !$omp do schedule(static) private(k,iq)  
  do k = 1,kn
    ! Apply JLM's preconditioner
    ua(1:ifull,k) = ud(1:ifull,k)-stul(1:ifull,k)*ud(ieu,k)
    va(1:ifull,k) = vd(1:ifull,k)-stvl(1:ifull,k)*vd(inv,k)
    ud(1:ifull,k) = ua(1:ifull,k)
    vd(1:ifull,k) = va(1:ifull,k)
  end do
  !$omp end do nowait
  !$omp do schedule(static) private(k,iq)  
  do k = kn+1,kx
    ! Apply JLM's preconditioner
    ua(1:ifull,k) = ud(1:ifull,k)-stul(1:ifull,1)*ud(ieu,k)
    va(1:ifull,k) = vd(1:ifull,k)-stvl(1:ifull,1)*vd(inv,k)
    ud(1:ifull,k) = ua(1:ifull,k)
    vd(1:ifull,k) = va(1:ifull,k)
  end do
  !$omp end do

  ! There are many ways to handle staggering near coastlines.
  ! This version supports the following properties:
  ! - Staggering is exactly reversible for the staggered (C grid) frame
  ! - Zero flux at land boundaries
  ! - The wave amplitude should be preserved

  do itn = 1,itnmax        ! each loop is a double iteration
    call boundsuv(ua,va,stag=2)
    !$omp barrier
    !$omp do schedule(static) private(k,iq)     
    do k = 1,kn
      do iq = 1,ifull  
        uin(iq,k)=ud(iq,k)+wtul(iq,k,1)*ua(ieu(iq),k)+wtul(iq,k,2)*ua(iwu(iq),k)+wtul(iq,k,3)*ua(ieeu(iq),k)
        vin(iq,k)=vd(iq,k)+wtvl(iq,k,1)*va(inv(iq),k)+wtvl(iq,k,2)*va(isv(iq),k)+wtvl(iq,k,3)*va(innv(iq),k)
      end do  
    end do
    !$omp end do nowait
    !$omp do schedule(static) private(k,iq)      
    do k = kn+1,kx
      do iq = 1,ifull  
        uin(iq,k)=ud(iq,k)+wtul(iq,1,1)*ua(ieu(iq),k)+wtul(iq,1,2)*ua(iwu(iq),k)+wtul(iq,1,3)*ua(ieeu(iq),k)
        vin(iq,k)=vd(iq,k)+wtvl(iq,1,1)*va(inv(iq),k)+wtvl(iq,1,2)*va(isv(iq),k)+wtvl(iq,1,3)*va(innv(iq),k)
      end do  
    end do
    !$omp end do
    call boundsuv(uin,vin,stag=2)
    !$omp barrier
    !$omp do schedule(static) private(k,iq)        
    do k = 1,kn
      do iq = 1,ifull  
        ua(iq,k)=ud(iq,k)+wtul(iq,k,1)*uin(ieu(iq),k)+wtul(iq,k,2)*uin(iwu(iq),k)+wtul(iq,k,3)*uin(ieeu(iq),k)
        va(iq,k)=vd(iq,k)+wtvl(iq,k,1)*vin(inv(iq),k)+wtvl(iq,k,2)*vin(isv(iq),k)+wtvl(iq,k,3)*vin(innv(iq),k)
      end do  
    end do
    !$omp end do nowait
    !$omp do schedule(static) private(k,iq)        
    do k = kn+1,kx
      do iq = 1,ifull  
        ua(iq,k)=ud(iq,k)+wtul(iq,1,1)*uin(ieu(iq),k)+wtul(iq,1,2)*uin(iwu(iq),k)+wtul(iq,1,3)*uin(ieeu(iq),k)
        va(iq,k)=vd(iq,k)+wtvl(iq,1,1)*vin(inv(iq),k)+wtvl(iq,1,2)*vin(isv(iq),k)+wtvl(iq,1,3)*vin(innv(iq),k)
      end do  
    end do
    !$omp end do
  end do                 ! itn=1,itnmax
  
  !$omp end parallel
 
else

  call boundsuv(uin,vin)
  !$omp parallel 
  !$omp do schedule(static) private(k,iq)      
  do k = 1,kn
    do iq = 1,ifull  
      ud(iq,k)=dtur(iq,k,1)*uin(iwu(iq),k)+dtur(iq,k,2)*uin(iq,k)+dtur(iq,k,3)*uin(ieu(iq),k)
      vd(iq,k)=dtvr(iq,k,1)*vin(isv(iq),k)+dtvr(iq,k,2)*vin(iq,k)+dtvr(iq,k,3)*vin(inv(iq),k)
    end do  
  end do
  !$omp end do nowait
  !$omp do schedule(static) private(k,iq)     
  do k = kn+1,kx
    do iq = 1,ifull  
      ud(iq,k)=dtur(iq,1,1)*uin(iwu(iq),k)+dtur(iq,1,2)*uin(iq,k)+dtur(iq,1,3)*uin(ieu(iq),k)
      vd(iq,k)=dtvr(iq,1,1)*vin(isv(iq),k)+dtvr(iq,1,2)*vin(iq,k)+dtvr(iq,1,3)*vin(inv(iq),k)
    end do  
  end do
  !$omp end do

  call boundsuv(ud,vd,stag=-9)
  !$omp barrier
  !$omp do schedule(static) private(k,iq)   
  do k = 1,kn
    ! Apply JLM's preconditioner
    ua(1:ifull,k) = ud(1:ifull,k)-stur(1:ifull,k)*ud(iwu,k)
    va(1:ifull,k) = vd(1:ifull,k)-stvr(1:ifull,k)*vd(isv,k)
    ud(1:ifull,k) = ua(1:ifull,k)
    vd(1:ifull,k) = va(1:ifull,k)    
  end do
  !$omp end do nowait
  !$omp do schedule(static) private(k,iq)      
  do k = kn+1,kx
    ! Apply JLM's preconditioner
    ua(1:ifull,k) = ud(1:ifull,k)-stur(1:ifull,1)*ud(iwu,k)
    va(1:ifull,k) = vd(1:ifull,k)-stvr(1:ifull,1)*vd(isv,k)
    ud(1:ifull,k) = ua(1:ifull,k)
    vd(1:ifull,k) = va(1:ifull,k)
  end do
  !$omp end do

  do itn = 1,itnmax        ! each loop is a double iteration
    call boundsuv(ua,va,stag=3)
    !$omp barrier
    !$omp do schedule(static) private(k,iq)       
    do k = 1,kn
      do iq = 1,ifull  
        uin(iq,k)=ud(iq,k)+wtur(iq,k,1)*ua(ieu(iq),k)+wtur(iq,k,2)*ua(iwu(iq),k)+wtur(iq,k,3)*ua(iwwu(iq),k)
        vin(iq,k)=vd(iq,k)+wtvr(iq,k,1)*va(inv(iq),k)+wtvr(iq,k,2)*va(isv(iq),k)+wtvr(iq,k,3)*va(issv(iq),k)
      end do  
    end do
    !$omp end do nowait
    !$omp do schedule(static) private(k,iq)      
    do k = kn+1,kx
      do iq = 1,ifull  
        uin(iq,k)=ud(iq,k)+wtur(iq,1,1)*ua(ieu(iq),k)+wtur(iq,1,2)*ua(iwu(iq),k)+wtur(iq,1,3)*ua(iwwu(iq),k)
        vin(iq,k)=vd(iq,k)+wtvr(iq,1,1)*va(inv(iq),k)+wtvr(iq,1,2)*va(isv(iq),k)+wtvr(iq,1,3)*va(issv(iq),k)
      end do  
    end do
    !$omp end do
    call boundsuv(uin,vin,stag=3)
    !$omp barrier
    !$omp do schedule(static) private(k,iq)      
    do k = 1,kn
      do iq = 1,ifull  
        ua(iq,k)=ud(iq,k)+wtur(iq,k,1)*uin(ieu(iq),k)+wtur(iq,k,2)*uin(iwu(iq),k)+wtur(iq,k,3)*uin(iwwu(iq),k)
        va(iq,k)=vd(iq,k)+wtvr(iq,k,1)*vin(inv(iq),k)+wtvr(iq,k,2)*vin(isv(iq),k)+wtvr(iq,k,3)*vin(issv(iq),k)
      end do  
    end do
    !$omp end do nowait
    !$omp do schedule(static) private(k,iq)     
    do k = kn+1,kx
      do iq = 1,ifull  
        ua(iq,k)=ud(iq,k)+wtur(iq,1,1)*uin(ieu(iq),k)+wtur(iq,1,2)*uin(iwu(iq),k)+wtur(iq,1,3)*uin(iwwu(iq),k)
        va(iq,k)=vd(iq,k)+wtvr(iq,1,1)*vin(inv(iq),k)+wtvr(iq,1,2)*vin(isv(iq),k)+wtvr(iq,1,3)*vin(issv(iq),k)
      end do  
    end do
    !$omp end do
  end do                 ! itn=1,itnmax
  
  !$omp end parallel

end if

do k = 1,kn
  uout(1:ifull,k)=ua(1:ifull,k)*eeu(1:ifull,k)
  vout(1:ifull,k)=va(1:ifull,k)*eev(1:ifull,k)
end do
do k = kn+1,kx
  uout(1:ifull,k)=ua(1:ifull,k)*eeu(1:ifull,1)
  vout(1:ifull,k)=va(1:ifull,k)*eev(1:ifull,1)
end do

call END_LOG(ocnstag_end)

return
end subroutine mlostaguv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Unstagger u and v
! Modified to include zero velocity over land points
subroutine mlounstaguv(u,v,uout,vout,toff)

use cc_mpi
use indices_m
use mlo
use mlodynamicsarrays_m
use newmpar_m
use parm_m

implicit none

integer, intent(in), optional :: toff
integer i,k,itn,kx,kn,zoff,iq
real, dimension(:,:), intent(in) :: u,v
real, dimension(:,:), intent(out) :: uout,vout
real, dimension(ifull+iextra,size(u,2)) :: uin,vin
real, dimension(ifull+iextra,size(u,2)) :: ua,va
real, dimension(ifull+iextra,size(u,2)) :: ud,vd
real, dimension(:,:), allocatable :: nwtu,nwtv
real, dimension(:,:,:), allocatable, save :: wtul,wtvl
real, dimension(:,:,:), allocatable, save :: wtur,wtvr
real, dimension(:,:,:), allocatable, save :: dtul,dtvl
real, dimension(:,:,:), allocatable, save :: dtur,dtvr
real, dimension(:,:), allocatable, save :: stul,stvl
real, dimension(:,:), allocatable, save :: stur,stvr
logical eutest,evtest
logical eetest,ewtest,entest,estest
logical euetest,euwtest,evntest,evstest
logical euewtest,evnstest
logical eeetest,ewwtest,enntest,esstest
logical ltest

call START_LOG(ocnstag_begin)

if (.not.allocated(wtul)) then
  allocate(wtul(ifull+iextra,wlev,0:3),wtvl(ifull+iextra,wlev,0:3))
  allocate(wtur(ifull+iextra,wlev,0:3),wtvr(ifull+iextra,wlev,0:3))
  allocate(dtul(ifull,wlev,3),dtvl(ifull,wlev,3))
  allocate(dtur(ifull,wlev,3),dtvr(ifull,wlev,3))
  allocate(stul(ifull,wlev),stvl(ifull,wlev))
  allocate(stur(ifull,wlev),stvr(ifull,wlev))

  allocate( nwtu(ifull,0:3), nwtv(ifull,0:3) )

  do k = 1,wlev
    do iq = 1,ifull
      ! assign land arrays
      eetest = ee(iq,k)*ee(ie(iq),k)>0.5
      entest = ee(iq,k)*ee(in(iq),k)>0.5

      eutest = ee(iw(iq),k)*ee(iq,k)>0.5
      evtest = ee(is(iq),k)*ee(iq,k)>0.5
      euetest = eutest .and. ee(iq,k)*ee(ie(iq),k)>0.5
      evntest = evtest .and. ee(iq,k)*ee(in(iq),k)>0.5
      euwtest = eutest .and. ee(iww(iq),k)*ee(iw(iq),k)>0.5
      evstest = evtest .and. ee(iss(iq),k)*ee(is(iq),k)>0.5
      euewtest = euetest .and. euwtest
      evnstest = evntest .and. evstest
      eeetest = eetest .and. ee(iee(iq),k)>0.5
      enntest = entest .and. ee(inn(iq),k)>0.5

      ! assign weights (left)

      !  |   W   | X *   |  E   |     unstaggered
      ! WW       W X     *            staggered
      if ( euewtest ) then
      wtul(iq,k,0)=1.
      wtul(iq,k,1)=-0.1
      wtul(iq,k,2)=-0.5
      wtul(iq,k,3)=0.
      dtul(iq,k,1)=0.1
      dtul(iq,k,2)=1.
      dtul(iq,k,3)=0.5

      ! ##   W   | X *   |  E   |     unstaggered
      ! #0       W X     *            staggered
      else if ( euetest ) then
        wtul(iq,k,0)=1.
        wtul(iq,k,1)=-0.1
        wtul(iq,k,2)=-0.5
        wtul(iq,k,3)=0.
        dtul(iq,k,1)=0.
        dtul(iq,k,2)=1.
        dtul(iq,k,3)=0.5
      
      !  |   W   | X *   #  ##  #     unstaggered
      !          W X     0  ##  #     staggered
      else if ( euwtest ) then
        wtul(iq,k,0)=1.
        wtul(iq,k,1)=0.
        wtul(iq,k,2)=-1./3.
        wtul(iq,k,3)=0.
        dtul(iq,k,1)=0.
        dtul(iq,k,2)=1.
        dtul(iq,k,3)=0.

      ! ##   ##  #   * X |   E   |     unstaggered
      ! ##   ##  0     X *             staggered
      else if ( eeetest ) then
        wtul(iq,k,0)=1.
        wtul(iq,k,1)=-1./3.
        wtul(iq,k,2)=0.
        wtul(iq,k,3)=0.
        dtul(iq,k,1)=0.
        dtul(iq,k,2)=0.
        dtul(iq,k,3)=1.

      ! ##   ##  #   *   |   E    #     unstaggered
      ! ##   ##  #       *        #     staggered
      else if ( eetest ) then
        wtul(iq,k,0)=1.
        wtul(iq,k,1)=0.
        wtul(iq,k,2)=0.
        wtul(iq,k,3)=0.
        dtul(iq,k,1)=0.
        dtul(iq,k,2)=0.
        dtul(iq,k,3)=1.

      ! ##   W   |   *   #  ##  #     unstaggered
      ! ##       W       #  ##  #     staggered
      else if ( eutest ) then
        wtul(iq,k,0)=1.
        wtul(iq,k,1)=0.
        wtul(iq,k,2)=0.
        wtul(iq,k,3)=0.
        dtul(iq,k,1)=0.
        dtul(iq,k,2)=1.
        dtul(iq,k,3)=0.

      else
        wtul(iq,k,0)=0.
        wtul(iq,k,1)=0.
        wtul(iq,k,2)=0.
        wtul(iq,k,3)=0.
        dtul(iq,k,1)=0.
        dtul(iq,k,2)=0.
        dtul(iq,k,3)=0.    
      end if

      if ( evnstest ) then
        wtvl(iq,k,0)=1.
        wtvl(iq,k,1)=-0.1
        wtvl(iq,k,2)=-0.5
        wtvl(iq,k,3)=0.
        dtvl(iq,k,1)=0.1
        dtvl(iq,k,2)=1.
        dtvl(iq,k,3)=0.5
      else if ( evntest ) then
        wtvl(iq,k,0)=1.
        wtvl(iq,k,1)=-0.1
        wtvl(iq,k,2)=-0.5
        wtvl(iq,k,3)=0.
        dtvl(iq,k,1)=0.
        dtvl(iq,k,2)=1.
        dtvl(iq,k,3)=0.5
      else if ( evstest ) then
        wtvl(iq,k,0)=1.
        wtvl(iq,k,1)=0.
        wtvl(iq,k,2)=-1./3.
        wtvl(iq,k,3)=0.
        dtvl(iq,k,1)=0.
        dtvl(iq,k,2)=1.
        dtvl(iq,k,3)=0.
      else if ( enntest ) then
        wtvl(iq,k,0)=1.
        wtvl(iq,k,1)=-1./3.
        wtvl(iq,k,2)=0.
        wtvl(iq,k,3)=0.
        dtvl(iq,k,1)=0.
        dtvl(iq,k,2)=0.
        dtvl(iq,k,3)=1.
      else if ( entest ) then
        wtvl(iq,k,0)=1.
        wtvl(iq,k,1)=0.
        wtvl(iq,k,2)=0.
        wtvl(iq,k,3)=0.
        dtvl(iq,k,1)=0.
        dtvl(iq,k,2)=0.
        dtvl(iq,k,3)=1.
      else if ( evtest ) then
        wtvl(iq,k,0)=1.
        wtvl(iq,k,1)=0.
        wtvl(iq,k,2)=0.
        wtvl(iq,k,3)=0.
        dtvl(iq,k,1)=0.
        dtvl(iq,k,2)=1.
        dtvl(iq,k,3)=0.
      else
        wtvl(iq,k,0)=0.
        wtvl(iq,k,1)=0.
        wtvl(iq,k,2)=0.
        wtvl(iq,k,3)=0.
        dtvl(iq,k,1)=0.
        dtvl(iq,k,2)=0.
        dtvl(iq,k,3)=0.
      end if
    end do
  end do  
    
  ! Apply JLM's preconditioner
  do i = 0,3
    call boundsuv(wtul(:,:,i),wtvl(:,:,i),stag=-9)
  end do   
  do k = 1,wlev
    where (abs(wtul(iwu,k,0))>1.E-4.and.abs(wtul(1:ifull,k,2))>1.E-4)
      stul(1:ifull,k)=-wtul(1:ifull,k,2)/wtul(iwu,k,0)
      nwtu(1:ifull,0)=wtul(1:ifull,k,0)+stul(:,k)*wtul(iwu,k,1)
      nwtu(1:ifull,1)=wtul(1:ifull,k,1)
      nwtu(1:ifull,2)=wtul(1:ifull,k,2)+stul(:,k)*wtul(iwu,k,0)
      nwtu(1:ifull,3)=wtul(1:ifull,k,3)-stul(:,k)*wtul(iwu,k,2)
    elsewhere
      stul(1:ifull,k)=0.
      nwtu(1:ifull,0)=wtul(1:ifull,k,0)
      nwtu(1:ifull,1)=wtul(1:ifull,k,1)
      nwtu(1:ifull,2)=wtul(1:ifull,k,2)
      nwtu(1:ifull,3)=wtul(1:ifull,k,3)
    end where
    where (abs(wtvl(isv,k,0))>1.E-4.and.abs(wtvl(1:ifull,k,2))>1.E-4)
      stvl(1:ifull,k)=-wtvl(1:ifull,k,2)/wtvl(isv,k,0)
      nwtv(1:ifull,0)=wtvl(1:ifull,k,0)+stvl(:,k)*wtvl(isv,k,1)
      nwtv(1:ifull,1)=wtvl(1:ifull,k,1)
      nwtv(1:ifull,2)=wtvl(1:ifull,k,2)+stvl(:,k)*wtvl(isv,k,0)
      nwtv(1:ifull,3)=wtvl(1:ifull,k,3)-stvl(:,k)*wtvl(isv,k,2)
    elsewhere
      stvl(1:ifull,k)=0.
      nwtv(1:ifull,0)=wtvl(1:ifull,k,0)
      nwtv(1:ifull,1)=wtvl(1:ifull,k,1)
      nwtv(1:ifull,2)=wtvl(1:ifull,k,2)
      nwtv(1:ifull,3)=wtvl(1:ifull,k,3)
    end where
    where (abs(wtul(1:ifull,k,0))<1.E-4)
      wtul(1:ifull,k,0)=1.
      wtul(1:ifull,k,1)=0.
      wtul(1:ifull,k,2)=0.
      wtul(1:ifull,k,3)=0.
    elsewhere
      wtul(1:ifull,k,0)=nwtu(1:ifull,0)
      wtul(1:ifull,k,1)=nwtu(1:ifull,1)
      wtul(1:ifull,k,2)=nwtu(1:ifull,2)
      wtul(1:ifull,k,3)=nwtu(1:ifull,3)
    end where
    where (abs(wtvl(1:ifull,k,0))<1.E-4)
      wtvl(1:ifull,k,0)=1.
      wtvl(1:ifull,k,1)=0.
      wtvl(1:ifull,k,2)=0.
      wtvl(1:ifull,k,3)=0.
    elsewhere
      wtvl(1:ifull,k,0)=nwtv(1:ifull,0)
      wtvl(1:ifull,k,1)=nwtv(1:ifull,1)
      wtvl(1:ifull,k,2)=nwtv(1:ifull,2)
      wtvl(1:ifull,k,3)=nwtv(1:ifull,3)
    end where
  end do  

  ! normalise
  do i = 0,3
    call boundsuv(wtul(:,:,i),wtvl(:,:,i),stag=-9)
  end do  
  do k = 1,wlev
    do i = 1,3
      wtul(1:ifull,k,i)=wtul(1:ifull,k,i)/wtul(1:ifull,k,0)
      wtvl(1:ifull,k,i)=wtvl(1:ifull,k,i)/wtvl(1:ifull,k,0)
      dtul(1:ifull,k,i)=dtul(1:ifull,k,i)/wtul(1:ifull,k,0)
      dtvl(1:ifull,k,i)=dtvl(1:ifull,k,i)/wtvl(1:ifull,k,0)
    end do
    stul(1:ifull,k)=stul(1:ifull,k)*wtul(iwu,k,0)/wtul(1:ifull,k,0)
    stvl(1:ifull,k)=stvl(1:ifull,k)*wtvl(isv,k,0)/wtvl(1:ifull,k,0)
    wtul(1:ifull,k,0)=1.
    wtvl(1:ifull,k,0)=1.
  end do  

  ! assign weights (right) 
  do k = 1,wlev
    do iq = 1,ifull
      
      ! assign land arrays  
      ewtest = ee(iw(iq),k)*ee(iq,k)>0.5
      estest = ee(is(iq),k)*ee(iq,k)>0.5
      
      eutest = ee(iq,k)*ee(ie(iq),k)>0.5
      evtest = ee(iq,k)*ee(in(iq),k)>0.5
      euetest = eutest .and. ee(ie(iq),k)*ee(iee(iq),k)>0.5
      evntest = evtest .and. ee(in(iq),k)*ee(inn(iq),k)>0.5
      euwtest = eutest .and. ee(iw(iq),k)*ee(iq,k)>0.5
      evstest = evtest .and. ee(is(iq),k)*ee(iq,k)>0.5
      euewtest = euetest .and. euwtest
      evnstest = evntest .and. evstest
      ewwtest = ewtest .and. ee(iww(iq),k)>0.5
      esstest = estest .and. ee(iss(iq),k)>0.5
  
  
      !  |   W   |   * X |  E   |     unstaggered
      !          W     X *      E     staggered
      if ( euewtest ) then
        wtur(iq,k,0)=1.
        wtur(iq,k,1)=-0.5
        wtur(iq,k,2)=-0.1
        wtur(iq,k,3)=0.
        dtur(iq,k,1)=0.1
        dtur(iq,k,2)=1.
        dtur(iq,k,3)=0.5
        
      !  |   W   |   * X |  E   #     unstaggered
      !          W     X *      0     staggered
      else if ( euwtest ) then
        wtur(iq,k,0)=1.
        wtur(iq,k,1)=-0.5
        wtur(iq,k,2)=-0.1
        wtur(iq,k,3)=0.
        dtur(iq,k,1)=0.
        dtur(iq,k,2)=1.
        dtur(iq,k,3)=0.5
      
      !  #   ##  #   * X |   E   |     unstaggered
      !  #   ##  0     X *             staggered
      else if ( euetest ) then
        wtur(iq,k,0)=1.
        wtur(iq,k,1)=-1./3.
        wtur(iq,k,2)=0.
        wtur(iq,k,3)=0.
        dtur(iq,k,1)=0.
        dtur(iq,k,2)=1.
        dtur(iq,k,3)=0.

      !  |   W   | X *   #  ##  #     unstaggered
      !          W X     0  ##  #     staggered
      else if ( ewwtest ) then
        wtur(iq,k,0)=1.
        wtur(iq,k,1)=0.
        wtur(iq,k,2)=-1./3.
        wtur(iq,k,3)=0.
        dtur(iq,k,1)=0.
        dtur(iq,k,2)=0.
        dtur(iq,k,3)=1.

      !  #   W   |   *   #  ##  #     unstaggered
      !  #       W       #  ##  #     staggered
      else if ( ewtest ) then
        wtur(iq,k,0)=1.
        wtur(iq,k,1)=0.
        wtur(iq,k,2)=0.
        wtur(iq,k,3)=0.
        dtur(iq,k,1)=0.
        dtur(iq,k,2)=0.
        dtur(iq,k,3)=1.

      !  #   ##  #   *   |  E   #     unstaggered
      !  #   ##  #       *      #     staggered
      else if ( eutest ) then
        wtur(iq,k,0)=1.
        wtur(iq,k,1)=0.
        wtur(iq,k,2)=0.
        wtur(iq,k,3)=0.
        dtur(iq,k,1)=0.
        dtur(iq,k,2)=1.
        dtur(iq,k,3)=0.
      else
        wtur(iq,k,0)=0.
        wtur(iq,k,1)=0.
        wtur(iq,k,2)=0.
        wtur(iq,k,3)=0.
        dtur(iq,k,1)=0.
        dtur(iq,k,2)=0.
        dtur(iq,k,3)=0.
      end if

      if ( evnstest ) then
        wtvr(iq,k,0)=1.
        wtvr(iq,k,1)=-0.5
        wtvr(iq,k,2)=-0.1
        wtvr(iq,k,3)=0.
        dtvr(iq,k,1)=0.1
        dtvr(iq,k,2)=1.
        dtvr(iq,k,3)=0.5
      else if ( evstest ) then
        wtvr(iq,k,0)=1.
        wtvr(iq,k,1)=-0.5
        wtvr(iq,k,2)=-0.1
        wtvr(iq,k,3)=0.
        dtvr(iq,k,1)=0.
        dtvr(iq,k,2)=1.
        dtvr(iq,k,3)=0.5
      else if ( evntest ) then
        wtvr(iq,k,0)=1.
        wtvr(iq,k,1)=-1./3.
        wtvr(iq,k,2)=0.
        wtvr(iq,k,3)=0.
        dtvr(iq,k,1)=0.
        dtvr(iq,k,2)=1.
        dtvr(iq,k,3)=0.
      else if ( esstest ) then
        wtvr(iq,k,0)=1.
        wtvr(iq,k,1)=0.
        wtvr(iq,k,2)=-1./3.
        wtvr(iq,k,3)=0.
        dtvr(iq,k,1)=0.
        dtvr(iq,k,2)=0.
        dtvr(iq,k,3)=1.
      else if ( estest ) then
        wtvr(iq,k,0)=1.
        wtvr(iq,k,1)=0.
        wtvr(iq,k,2)=0.
        wtvr(iq,k,3)=0.
        dtvr(iq,k,1)=0.
        dtvr(iq,k,2)=0.
        dtvr(iq,k,3)=1.
      else if ( evtest ) then
        wtvr(iq,k,0)=1.
        wtvr(iq,k,1)=0.
        wtvr(iq,k,2)=0.
        wtvr(iq,k,3)=0.
        dtvr(iq,k,1)=0.
        dtvr(iq,k,2)=1.
        dtvr(iq,k,3)=0.
      else
        wtvr(iq,k,0)=0.
        wtvr(iq,k,1)=0.
        wtvr(iq,k,2)=0.
        wtvr(iq,k,3)=0.
        dtvr(iq,k,1)=0.
        dtvr(iq,k,2)=0.
        dtvr(iq,k,3)=0.
      end if
    end do
  end do  

  ! Apply JLM's preconditioner
  do i = 0,3
    call boundsuv(wtur(:,:,i),wtvr(:,:,i),stag=-10)
  end do  
  do k = 1,wlev
    where (abs(wtur(ieu,k,0))>1.E-4.and.abs(wtur(1:ifull,k,1))>1.E-4)
      stur(1:ifull,k)=-wtur(1:ifull,k,1)/wtur(ieu,k,0)
      nwtu(1:ifull,0)=wtur(1:ifull,k,0)+stur(:,k)*wtur(ieu,k,2)
      nwtu(1:ifull,1)=wtur(1:ifull,k,1)+stur(:,k)*wtur(ieu,k,0)
      nwtu(1:ifull,2)=wtur(1:ifull,k,2)
      nwtu(1:ifull,3)=wtur(1:ifull,k,3)-stur(:,k)*wtur(ieu,k,1)
    elsewhere
      stur(1:ifull,k)=0.
      nwtu(1:ifull,0)=wtur(1:ifull,k,0)
      nwtu(1:ifull,1)=wtur(1:ifull,k,1)
      nwtu(1:ifull,2)=wtur(1:ifull,k,2)
      nwtu(1:ifull,3)=wtur(1:ifull,k,3)
    end where
    where (abs(wtvr(inv,k,0))>1.E-4.and.abs(wtvr(1:ifull,k,1))>1.E-4)
      stvr(1:ifull,k)=-wtvr(1:ifull,k,1)/wtvr(inv,k,0)
      nwtv(1:ifull,0)=wtvr(1:ifull,k,0)+stvr(:,k)*wtvr(inv,k,2)
      nwtv(1:ifull,1)=wtvr(1:ifull,k,1)+stvr(:,k)*wtvr(inv,k,0)
      nwtv(1:ifull,2)=wtvr(1:ifull,k,2)
      nwtv(1:ifull,3)=wtvr(1:ifull,k,3)-stvr(:,k)*wtvr(inv,k,1)
    elsewhere
      stvr(1:ifull,k)=0.
      nwtv(1:ifull,0)=wtvr(1:ifull,k,0)
      nwtv(1:ifull,1)=wtvr(1:ifull,k,1)
      nwtv(1:ifull,2)=wtvr(1:ifull,k,2)
      nwtv(1:ifull,3)=wtvr(1:ifull,k,3)
    end where
    where (abs(wtur(1:ifull,k,0))<1.E-4)
      wtur(1:ifull,k,0)=1.
      wtur(1:ifull,k,1)=0.
      wtur(1:ifull,k,2)=0.
      wtur(1:ifull,k,3)=0.
    elsewhere
      wtur(1:ifull,k,0)=nwtu(1:ifull,0)
      wtur(1:ifull,k,1)=nwtu(1:ifull,1)
      wtur(1:ifull,k,2)=nwtu(1:ifull,2)
      wtur(1:ifull,k,3)=nwtu(1:ifull,3)
    end where
    where (abs(wtvr(1:ifull,k,0))<1.E-4)
      wtvr(1:ifull,k,0)=1.
      wtvr(1:ifull,k,1)=0.
      wtvr(1:ifull,k,2)=0.
      wtvr(1:ifull,k,3)=0.
    elsewhere
      wtvr(1:ifull,k,0)=nwtv(1:ifull,0)
      wtvr(1:ifull,k,1)=nwtv(1:ifull,1)
      wtvr(1:ifull,k,2)=nwtv(1:ifull,2)
      wtvr(1:ifull,k,3)=nwtv(1:ifull,3)
    end where
  end do  

  ! normalise
  do i = 0,3
    call boundsuv(wtur(:,:,i),wtvr(:,:,i),stag=-10)
  end do  
  do k = 1,wlev
    do i = 1,3
      wtur(1:ifull,k,i)=wtur(1:ifull,k,i)/wtur(1:ifull,k,0)
      wtvr(1:ifull,k,i)=wtvr(1:ifull,k,i)/wtvr(1:ifull,k,0)
      dtur(1:ifull,k,i)=dtur(1:ifull,k,i)/wtur(1:ifull,k,0)
      dtvr(1:ifull,k,i)=dtvr(1:ifull,k,i)/wtvr(1:ifull,k,0)
    end do
    stur(1:ifull,k)=stur(1:ifull,k)*wtur(ieu,k,0)/wtur(1:ifull,k,0)
    stvr(1:ifull,k)=stvr(1:ifull,k)*wtvr(inv,k,0)/wtvr(1:ifull,k,0)
    wtur(1:ifull,k,0)=1.
    wtvr(1:ifull,k,0)=1.
  end do  

  deallocate( nwtu, nwtv )

end if

kx = size(u,2)
kn = min(kx,wlev)

zoff=0
if (present(toff)) then
  if (toff==1) zoff=koff
end if

do k = 1,kn
  uin(1:ifull,k)=u(1:ifull,k)*eeu(1:ifull,k)
  vin(1:ifull,k)=v(1:ifull,k)*eev(1:ifull,k)
end do
do k = kn+1,kx
  uin(1:ifull,k)=u(1:ifull,k)*eeu(1:ifull,1)
  vin(1:ifull,k)=v(1:ifull,k)*eev(1:ifull,1)
end do

if (mstagf==0) then
  ltest = .true.
else if (mstagf<0) then
  ltest = .false.
else
  ltest = mod(ktau-zoff-nstagoffmlo,2*mstagf)<mstagf
end if

if (ltest) then
  
  call boundsuv(uin,vin,stag=5)
  !$omp parallel 
  !$omp do schedule(static) private(k,iq)   
  do k = 1,kn
    do iq = 1,ifull  
      ud(iq,k)=dtul(iq,k,1)*uin(iwwu(iq),k)+dtul(iq,k,2)*uin(iwu(iq),k)+dtul(iq,k,3)*uin(iq,k)
      vd(iq,k)=dtvl(iq,k,1)*vin(issv(iq),k)+dtvl(iq,k,2)*vin(isv(iq),k)+dtvl(iq,k,3)*vin(iq,k)
    end do  
  end do
  !$omp end do nowait
  !$omp do schedule(static) private(k,iq)       
  do k = kn+1,kx
    do iq = 1,ifull  
      ud(iq,k)=dtul(iq,1,1)*uin(iwwu(iq),k)+dtul(iq,1,2)*uin(iwu(iq),k)+dtul(iq,1,3)*uin(iq,k)
      vd(iq,k)=dtvl(iq,1,1)*vin(issv(iq),k)+dtvl(iq,1,2)*vin(isv(iq),k)+dtvl(iq,1,3)*vin(iq,k)
    end do  
  end do
  !$omp end do

  call boundsuv(ud,vd,stag=-9)
  !$omp barrier
  !$omp do schedule(static) private(k,iq)    
  do k=1,kn
    ! Apply JLM's preconditioner
    ua(1:ifull,k) = ud(1:ifull,k)-stul(1:ifull,k)*ud(iwu,k)
    va(1:ifull,k) = vd(1:ifull,k)-stvl(1:ifull,k)*vd(isv,k)
    ud(1:ifull,k) = ua(1:ifull,k)
    vd(1:ifull,k) = va(1:ifull,k)    
  end do
  !$omp end do nowait
  !$omp do schedule(static) private(k,iq)   
  do k = kn+1,kx
    ! Apply JLM's preconditioner
    ua(1:ifull,k) = ud(1:ifull,k)-stul(1:ifull,1)*ud(iwu,k)
    va(1:ifull,k) = vd(1:ifull,k)-stvl(1:ifull,1)*vd(isv,k)
    ud(1:ifull,k) = ua(1:ifull,k)
    vd(1:ifull,k) = va(1:ifull,k)    
  end do
  !$omp end do

  do itn = 1,itnmax        ! each loop is a double iteration
    call boundsuv(ua,va,stag=3)
    !$omp barrier
    !$omp do schedule(static) private(k,iq)      
    do k = 1,kn
      do iq = 1,ifull  
        uin(iq,k)=ud(iq,k)+wtul(iq,k,1)*ua(ieu(iq),k)+wtul(iq,k,2)*ua(iwu(iq),k)+wtul(iq,k,3)*ua(iwwu(iq),k)
        vin(iq,k)=vd(iq,k)+wtvl(iq,k,1)*va(inv(iq),k)+wtvl(iq,k,2)*va(isv(iq),k)+wtvl(iq,k,3)*va(issv(iq),k)
      end do  
    end do
    !$omp end do nowait
    !$omp do schedule(static) private(k,iq)      
    do k = kn+1,kx
      do iq = 1,ifull  
        uin(iq,k)=ud(iq,k)+wtul(iq,1,1)*ua(ieu(iq),k)+wtul(iq,1,2)*ua(iwu(iq),k)+wtul(iq,1,3)*ua(iwwu(iq),k)
        vin(iq,k)=vd(iq,k)+wtvl(iq,1,1)*va(inv(iq),k)+wtvl(iq,1,2)*va(isv(iq),k)+wtvl(iq,1,3)*va(issv(iq),k)
      end do  
    end do
    !$omp end do
    call boundsuv(uin,vin,stag=3)
    !$omp barrier
    !$omp do schedule(static) private(k,iq)     
    do k = 1,kn
      do iq = 1,ifull  
        ua(iq,k)=ud(iq,k)+wtul(iq,k,1)*uin(ieu(iq),k)+wtul(iq,k,2)*uin(iwu(iq),k)+wtul(iq,k,3)*uin(iwwu(iq),k)
        va(iq,k)=vd(iq,k)+wtvl(iq,k,1)*vin(inv(iq),k)+wtvl(iq,k,2)*vin(isv(iq),k)+wtvl(iq,k,3)*vin(issv(iq),k)
      end do  
    end do
    !$omp end do nowait
    !$omp do schedule(static) private(k,iq)     
    do k = kn+1,kx
      do iq = 1,ifull  
        ua(iq,k)=ud(iq,k)+wtul(iq,1,1)*uin(ieu(iq),k)+wtul(iq,1,2)*uin(iwu(iq),k)+wtul(iq,1,3)*uin(iwwu(iq),k)
        va(iq,k)=vd(iq,k)+wtvl(iq,1,1)*vin(inv(iq),k)+wtvl(iq,1,2)*vin(isv(iq),k)+wtvl(iq,1,3)*vin(issv(iq),k)
      end do  
    end do
    !$omp end do
  end do                  ! itn=1,itnmax
  
  !$omp end parallel
  
else

  call boundsuv(uin,vin)
  !$omp parallel 
  !$omp do schedule(static) private(k,iq)   
  do k = 1,kn
    do iq = 1,ifull  
      ud(iq,k)=dtur(iq,k,1)*uin(ieu(iq),k)+dtur(iq,k,2)*uin(iq,k)+dtur(iq,k,3)*uin(iwu(iq),k)
      vd(iq,k)=dtvr(iq,k,1)*vin(inv(iq),k)+dtvr(iq,k,2)*vin(iq,k)+dtvr(iq,k,3)*vin(isv(iq),k)
    end do  
  end do
  !$omp end do nowait
  !$omp do schedule(static) private(k,iq)   
  do k = kn+1,kx
    do iq = 1,ifull  
      ud(iq,k)=dtur(iq,1,1)*uin(ieu(iq),k)+dtur(iq,1,2)*uin(iq,k)+dtur(iq,1,3)*uin(iwu(iq),k)
      vd(iq,k)=dtvr(iq,1,1)*vin(inv(iq),k)+dtvr(iq,1,2)*vin(iq,k)+dtvr(iq,1,3)*vin(isv(iq),k)
    end do  
  end do
  !$omp end do

  call boundsuv(ud,vd,stag=-10)
  !$omp barrier
  !$omp do schedule(static) private(k,iq)     
  do k = 1,kn
    ! Apply JLM's preconditioner
    ua(1:ifull,k) = ud(1:ifull,k)-stur(1:ifull,k)*ud(ieu,k)
    va(1:ifull,k) = vd(1:ifull,k)-stvr(1:ifull,k)*vd(inv,k)
    ud(1:ifull,k) = ua(1:ifull,k)
    vd(1:ifull,k) = va(1:ifull,k)     
  end do
  !$omp end do nowait
  !$omp do schedule(static) private(k,iq)    
  do k = kn+1,kx
    ! Apply JLM's preconditioner
    ua(1:ifull,k) = ud(1:ifull,k)-stur(1:ifull,1)*ud(ieu,k)
    va(1:ifull,k) = vd(1:ifull,k)-stvr(1:ifull,1)*vd(inv,k)
    ud(1:ifull,k) = ua(1:ifull,k)
    vd(1:ifull,k) = va(1:ifull,k)    
  end do
  !$omp end do

  do itn = 1,itnmax        ! each loop is a double iteration
    call boundsuv(ua,va,stag=2)
    !$omp barrier
    !$omp do schedule(static) private(k,iq)      
    do k = 1,kn
      do iq = 1,ifull  
        uin(iq,k)=ud(iq,k)+wtur(iq,k,1)*ua(ieu(iq),k)+wtur(iq,k,2)*ua(iwu(iq),k)+wtur(iq,k,3)*ua(ieeu(iq),k)
        vin(iq,k)=vd(iq,k)+wtvr(iq,k,1)*va(inv(iq),k)+wtvr(iq,k,2)*va(isv(iq),k)+wtvr(iq,k,3)*va(innv(iq),k)
      end do  
    end do
    !$omp end do nowait
    !$omp do schedule(static) private(k,iq)       
    do k = kn+1,kx
      do iq = 1,ifull  
        uin(iq,k)=ud(iq,k)+wtur(iq,1,1)*ua(ieu(iq),k)+wtur(iq,1,2)*ua(iwu(iq),k)+wtur(iq,1,3)*ua(ieeu(iq),k)
        vin(iq,k)=vd(iq,k)+wtvr(iq,1,1)*va(inv(iq),k)+wtvr(iq,1,2)*va(isv(iq),k)+wtvr(iq,1,3)*va(innv(iq),k)
      end do  
    end do
    !$omp end do
    call boundsuv(uin,vin,stag=2)
    !$omp barrier
    !$omp do schedule(static) private(k,iq)   
    do k = 1,kn
      do iq = 1,ifull  
        ua(iq,k)=ud(iq,k)+wtur(iq,k,1)*uin(ieu(iq),k)+wtur(iq,k,2)*uin(iwu(iq),k)+wtur(iq,k,3)*uin(ieeu(iq),k)
        va(iq,k)=vd(iq,k)+wtvr(iq,k,1)*vin(inv(iq),k)+wtvr(iq,k,2)*vin(isv(iq),k)+wtvr(iq,k,3)*vin(innv(iq),k)
      end do  
    end do
    !$omp end do nowait
    !$omp do schedule(static) private(k,iq)     
    do k = kn+1,kx
      do iq = 1,ifull  
        ua(iq,k)=ud(iq,k)+wtur(iq,1,1)*uin(ieu(iq),k)+wtur(iq,1,2)*uin(iwu(iq),k)+wtur(iq,1,3)*uin(ieeu(iq),k)
        va(iq,k)=vd(iq,k)+wtvr(iq,1,1)*vin(inv(iq),k)+wtvr(iq,1,2)*vin(isv(iq),k)+wtvr(iq,1,3)*vin(innv(iq),k)
      end do  
    end do
    !$omp end do
  end do                  ! itn=1,itnmax
  
  !$omp end parallel

end if

do k = 1,kn
  uout(1:ifull,k)=ua(1:ifull,k)*ee(1:ifull,k)
  vout(1:ifull,k)=va(1:ifull,k)*ee(1:ifull,k)
end do
do k = kn+1,kx
  uout(1:ifull,k)=ua(1:ifull,k)*ee(1:ifull,1)
  vout(1:ifull,k)=va(1:ifull,k)*ee(1:ifull,1)
end do

call END_LOG(ocnstag_end)

return
end subroutine mlounstaguv

end module mlostag
