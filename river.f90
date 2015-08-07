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
    
! This is the river routing which links to mlo.f90 and mlodynamics.f90
! ocean/lake model

! This version currently assumes that the orography is sufficently resolved
! to determine the river flow.  Plan to read in effective gradients between
! grid boxes from high resolution river flow datasets.    
    
module river

implicit none

private
public rvrinit, rvrrouter, watbdy

real, dimension(:), allocatable, save :: watbdy, ee
real, parameter :: leakrate    = 192.     ! E-folding time for leaking water into soil (hrs)
real, parameter :: maxwaterlvl = 1000.    ! Target maximum water level (mm)
integer, parameter :: basinmd = 3         ! basin mode (0=soil, 2=pile-up, 3=leak)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine initialises river arrays
!
subroutine rvrinit

use cc_mpi
use soil_m

implicit none

include 'newmpar.h'

! river water height
allocate(watbdy(ifull+iextra))
watbdy(1:ifull+iextra)=0.

! prep land-sea mask
allocate(ee(ifull+iextra))
ee(1:ifull+iextra)=0.
where(.not.land(1:ifull))
  ee(1:ifull)=1.
end where
call bounds(ee,nrows=2)

return
end subroutine rvrinit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculates the river routing.
!
subroutine rvrrouter

use arrays_m
use cable_ccam 
use cc_mpi
use indices_m
use map_m
use mlo
use nsibd_m
use soil_m
use soilsnow_m

implicit none

include 'newmpar.h'
include 'const_phys.h'
include 'parm.h'
include 'soilv.h'

integer i,k
integer nit
integer, dimension(ifull,4) :: xp
real, dimension(ifull+iextra) :: netflx,newwat,zsadj
real, dimension(ifull) :: vel,soilsink
real, dimension(ifull) :: tmpry,deltmpry,ll
real, dimension(ifull,4) :: idp,slope,flow
real, dimension(ifull,4) :: fta,ftb,ftx,fty
real rate

! To speed up the code, we use a (semi-)implicit solution rather than an iterative approach
! This avoids additional MPI calls.

! Note that unlike Miller et al (1994), this scheme currently has multiple outflows from a given
! grid box (i.e., following Gordon in Mk3.5).  So if the slope of the grid box is divergent, then
! on average water leaves the grid box, whereas if the slope of the grid box is convergent, then
! on average water enters the grid box.

! setup indices and grid spacing
xp(:,1)=in
xp(:,2)=ie
xp(:,3)=is
xp(:,4)=iw
idp(:,1)=emv(1:ifull)/ds
idp(:,2)=emu(1:ifull)/ds
idp(:,3)=emv(isv)/ds
idp(:,4)=emu(iwu)/ds
tmpry(:)=0. ! for cray compiler

!--------------------------------------------------------------------
! update processor boundaries
call bounds(watbdy)
newwat(1:ifull+iextra)=watbdy(1:ifull+iextra)

!--------------------------------------------------------------------
! predictor-corrector for water level and salinity
do nit=1,2
  
  if (nit==2) call bounds(newwat)

  ! calculate slopes
  do i=1,4
    zsadj(1:ifull+iextra)=0.1*max(0.5*(newwat(1:ifull+iextra)+watbdy(1:ifull+iextra))-maxwaterlvl,0.) ! increase at 100x
    where ( (ee(1:ifull)*ee(xp(:,i)))>0.5 )
      slope(:,i)=0. ! no orographic slope within ocean bounds
    elsewhere
      slope(:,i)=(zs(1:ifull)/grav+zsadj(1:ifull)-zs(xp(:,i))/grav-zsadj(xp(:,i)))*idp(:,i)
    end where
  end do

  newwat(1:ifull)=watbdy(1:ifull)

  ! Basic expression

  ! m = mass/area
  ! flow = m * vel / dx
  ! m(t+1)-m(t) = dt*sum(inflow)-dt*sum(outflow)

  ! outflow
  ! compute net outgoing flux for a grid box so that total water is conserved
  do i=1,4
    vel=min(0.35*sqrt(max(slope(:,i),0.)/0.00005),5.) ! from Miller et al (1994)
    fta(:,i)=-dt*vel*idp(:,i)     ! outgoing flux
  end do
  netflx(1:ifull)=sum(abs(fta),2) ! MJT notes - this will never trigger for sensible values of dt
  call bounds(netflx)
  
  ! water outflow
  do i=1,4
    where (netflx(1:ifull)>1.E-10)
      ftx(1:ifull,i)=-fta(1:ifull,i)/netflx(1:ifull) ! max fraction of total outgoing flux
      flow(1:ifull,i)=watbdy(1:ifull)*min(fta(1:ifull,i),ftx(1:ifull,i)) ! (kg/m^2)
    elsewhere
      flow(1:ifull,i)=0.
    end where
  end do
  newwat(1:ifull)=newwat(1:ifull)+sum(flow(1:ifull,1:4),2)

  ! inflow
  ! water inflow
  do i=1,4
    vel=min(0.35*sqrt(max(-slope(:,i),0.)/0.00005),5.) ! from Miller et al (1994)
    ftb(:,i)=dt*vel*idp(:,i)            ! incomming flux
    where (netflx(xp(1:ifull,i))>1.E-10)
      fty(1:ifull,i)=ftb(1:ifull,i)/netflx(xp(1:ifull,i)) ! max fraction of flux from outgoing cel
      flow(1:ifull,i)=watbdy(xp(1:ifull,i))*min(ftb(1:ifull,i),fty(1:ifull,i)) ! (kg/m^2)
      flow(1:ifull,i)=flow(1:ifull,i)*(em(1:ifull)/em(xp(1:ifull,i)))**2 ! change in gridbox area
    elsewhere
      flow(1:ifull,i)=0.
    end where
  end do
  newwat(1:ifull)=newwat(1:ifull)+sum(flow(1:ifull,1:4),2)

end do

watbdy(1:ifull)=max(newwat(1:ifull),0.)

!--------------------------------------------------------------------
! Water losses over land basins

  
! Method for land basins
select case(basinmd)
case(0)
    ! add water to soil moisture 
    ! estimate rate that water leaves river into soil
    rate=min(dt/(8.*3600.),1.) ! MJT suggestion
    if (nsib==6.or.nsib==7) then
      ! CABLE
      tmpry(1:ifull)=watbdy(1:ifull)
      call cableinflow(tmpry,rate)
      soilsink(1:ifull)=(tmpry(1:ifull)-watbdy(1:ifull))*(1.-sigmu(1:ifull))
      newwat(1:ifull)=newwat(1:ifull)+soilsink(1:ifull)
    else
      ! Standard land surface model
      deltmpry=0.
      do k=1,ms
        where (land(1:ifull))
          ll(:)=max(sfc(isoilm(:))-wb(:,k),0.)*1000.*zse(k)
          ll(:)=ll(:)*rate
          ll(:)=min(tmpry(:),ll(:))
          wb(:,k)=wb(:,k)+ll(:)/(1000.*zse(k))
          deltmpry(:)=deltmpry(:)-ll(:)
        end where
      end do
      soilsink=deltmpry*(1.-sigmu(:))
      newwat(1:ifull)=newwat(1:ifull)+soilsink
    end if
  case(2)
    ! pile-up water
  case(3)
    ! leak
    ! estimate rate that water leaves river into soil
    rate=dt/(leakrate*3600.) ! MJT suggestion
    if (nsib==6.or.nsib==7) then
      ! CABLE
      tmpry(1:ifull)=watbdy(1:ifull)
      call cableinflow(tmpry,rate)
      soilsink(1:ifull)=(tmpry(1:ifull)-watbdy(1:ifull))*(1.-sigmu(1:ifull))
      newwat(1:ifull)=newwat(1:ifull)+soilsink(1:ifull)
    else
      ! Standard land surface model
      deltmpry=0.
      do k=1,ms
        where (land(1:ifull))
          ll(1:ifull)=max(sfc(isoilm(1:ifull))-wb(1:ifull,k),0.)*1000.*zse(k)
          ll(1:ifull)=ll(1:ifull)*rate
          ll(1:ifull)=min(watbdy(1:ifull),ll(1:ifull))
          wb(1:ifull,k)=wb(1:ifull,k)+ll(1:ifull)/(1000.*zse(k))
          deltmpry(1:ifull)=deltmpry(1:ifull)-ll(1:ifull)
        end where
      end do
      soilsink(1:ifull)=deltmpry(1:ifull)*(1.-sigmu(1:ifull))
      newwat(1:ifull)=newwat(1:ifull)+soilsink(1:ifull)
    end if
  case default
    write(6,*) "ERROR: Unsupported basinmd ",basinmd
    call ccmpi_abort(-1)
end select

watbdy(1:ifull)=max(newwat(1:ifull),0.)

! MLO (or other ocean model) will remove watbdy from ocean points when it updates its
! river inflows.

return
end subroutine rvrrouter

end module river
