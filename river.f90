! These subroutines handle dynamics for the Mixed-Layer-Ocean model

! - Horizontal diffusion
! - River routing
! - Ocean dynamics
! - Ice dynamics

! This module also links to mlo.f90 which solves for 1D physics of
! ocean and ice processes.  Currently the code assumes the hydrostatic
! approximation which is reasonably valid to 1km resolution.

! Ocean and sea-ice dynamics are based on the R-grid used by CCAM

module river

implicit none

private
public rvrinit,rvrrouter,watbdy,salbdy

real, dimension(:), allocatable, save :: watbdy,salbdy,ee
real, parameter :: inflowbias = 10.       ! Bias height for inflow into ocean.  Levels above this will flow back onto land.
integer, parameter :: basinmd = 3         ! basin mode (0=soil, 1=redistribute, 2=pile-up, 3=leak)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine initialises mlo dynamical arrays
!
subroutine rvrinit

use cc_mpi
use soil_m

implicit none

include 'newmpar.h'

! river water height
allocate(watbdy(ifull+iextra),salbdy(ifull+iextra))
watbdy=0.
salbdy=0.

! prep land-sea mask
allocate(ee(ifull+iextra))
ee=0.
where(.not.land)
  ee(1:ifull)=1.
end where
call bounds(ee,nrows=2)

return
end subroutine rvrinit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculates the river routing.
! This version currently assumes that the orography is sufficently resolved
! to determine the river flow.  Plan to read in effective gradients between
! grid boxes from high resolution river flow datasets.
subroutine rvrrouter

use arrays_m
use cable_ccam 
use cc_mpi
use indices_m
use map_m
use mlo
use mlodynamics
use nsibd_m
use soil_m
use soilsnow_m

implicit none

include 'newmpar.h'
include 'const_phys.h'
include 'parm.h'
include 'soilv.h'

integer i,ii,iq,ierr,k
integer nit
integer, dimension(ifull,4) :: xp
real, dimension(ifull,wlev) :: sallvl
real, dimension(ifull+iextra) :: neta,netflx,cc
real, dimension(ifull+iextra) :: newwat
real, dimension(ifull+iextra,2) :: dum
real, dimension(ifull) :: newsal,cover,vel
real, dimension(ifull) :: deta,sal,salin,depdum
real, dimension(ifull) :: tmpry,ll
real, dimension(ifull,4) :: idp,slope,flow
real, dimension(ifull,4) :: fta,ftb,ftx,fty
real rate

! To speed up the code, we use a (semi-)implicit solution rather than an iterative approach
! This avoids additional MPI calls.

! Note that unlike Miller et al (1994), this scheme currently has multiple outflows from a given
! grid box (i.e., following Gordon in Mk3.5).  So if the slope of the grid box is divergent, then
! on average water leaves the grid box, whereas if the slope of the grid box is convergent, then
! on average water enters the grid box.

! This version supports salinity in rivers and overflow from lakes

! setup indices and grid spacing
xp(:,1)=in
xp(:,2)=ie
xp(:,3)=is
xp(:,4)=iw
idp(:,1)=emv(1:ifull)/ds
idp(:,2)=emu(1:ifull)/ds
idp(:,3)=emv(isv)/ds
idp(:,4)=emu(iwu)/ds

!--------------------------------------------------------------------
! extract data from ocean model if required
if (abs(nmlo)<=9) then
  neta=0.
  salin=0.
  sallvl=0.
  call mloexport(4,neta(1:ifull),0,0)
  neta(1:ifull)=neta(1:ifull)-inflowbias*ee(1:ifull)
  call mloexport3d(1,sallvl,0)
  do ii=1,wlev
    ! could modify the following to only operate above neta=0
    salin=salin+godsig(ii)*sallvl(:,ii)
  end do
  where (ee(1:ifull)>0.5)
    depdum=max(neta(1:ifull)+dd(1:ifull),minwater)
    ! collect any left over water over ocean
    salin=(salin*depdum+salbdy(1:ifull)*watbdy(1:ifull)*0.001)/(depdum+watbdy(1:ifull)*0.001)
    neta(1:ifull)=neta(1:ifull)+0.001*watbdy(1:ifull)
    ! move water from neta to watbdy for transport
    deta=1000.*max(neta(1:ifull),0.)
    salbdy(1:ifull)=salin
    watbdy(1:ifull)=deta
    neta(1:ifull)=min(neta(1:ifull),0.)
  elsewhere
    deta=0.
    depdum=0.
  end where
  salbdy(1:ifull)=salbdy(1:ifull)*watbdy(1:ifull) ! rescale salinity to PSU*mm for advection
end if

!--------------------------------------------------------------------
! update processor boundaries
dum(1:ifull,1)=watbdy(1:ifull)
dum(1:ifull,2)=salbdy(1:ifull)
call bounds(dum(:,1:2))
watbdy(ifull+1:ifull+iextra)=dum(ifull+1:ifull+iextra,1)
salbdy(ifull+1:ifull+iextra)=dum(ifull+1:ifull+iextra,2)

newwat(1:ifull+iextra)=watbdy(1:ifull+iextra)
newsal=salbdy(1:ifull)

!--------------------------------------------------------------------
! predictor-corrector for water level and salinity
do nit=1,2
  
  if (nit==2) call bounds(newwat)

  ! calculate slopes
  do i=1,4
    where ( (ee(1:ifull)*ee(xp(:,i)))>0.5 )
      slope(:,i)=0. ! no orographic slope within ocean bounds
    elsewhere
      !slope(:,i)=(zs(1:ifull)-zs(xp(:,i)))/(grav*dp(:,i))         ! basic
      slope(:,i)=(zs(1:ifull)/grav+0.001*newwat(1:ifull) &
                 -zs(xp(:,i))/grav-0.001*newwat(xp(:,i)))*idp(:,i) ! flood
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
  netflx(1:ifull)=sum(abs(fta),2)
  call bounds(netflx)
  
  ! water outflow
  do i=1,4
    where (netflx(1:ifull)>1.E-10)
      ftx(:,i)=-fta(:,i)/netflx(1:ifull) ! max fraction of total outgoing flux
      flow(:,i)=watbdy(1:ifull)*min(fta(:,i),ftx(:,i)) ! (kg/m^2)
    elsewhere
      flow(:,i)=0.
    end where
  end do
  newwat(1:ifull)=newwat(1:ifull)+sum(flow,2)

  ! inflow
  ! water inflow
  do i=1,4
    vel=min(0.35*sqrt(max(-slope(:,i),0.)/0.00005),5.) ! from Miller et al (1994)
    ftb(:,i)=dt*vel*idp(:,i)     ! incomming flux
    where (netflx(xp(:,i))>1.E-10)
      fty(:,i)=ftb(:,i)/netflx(xp(:,i)) ! max fraction of flux from outgoing cel
      flow(:,i)=watbdy(xp(:,i))*min(ftb(:,i),fty(:,i)) ! (kg/m^2)
      flow(:,i)=flow(:,i)*em(1:ifull)*em(1:ifull)/(em(xp(:,i))*em(xp(:,i))) ! change in gridbox area
    elsewhere
      flow(:,i)=0.
    end where
  end do
  newwat(1:ifull)=newwat(1:ifull)+sum(flow,2)

end do

! salinity outflow
do i=1,4
  where (netflx(1:ifull)>1.E-10)
    flow(:,i)=salbdy(1:ifull)*min(fta(:,i),ftx(:,i))
  elsewhere
    flow(:,i)=0.
  end where
end do
newsal=newsal+sum(flow,2)
! salinity inflow
do i=1,4
  where (netflx(xp(:,i))>1.E-10)
    flow(:,i)=salbdy(xp(:,i))*min(ftb(:,i),fty(:,i))
    flow(:,i)=flow(:,i)*em(1:ifull)*em(1:ifull)/(em(xp(:,i))*em(xp(:,i))) ! change in gridbox area
  elsewhere
    flow(:,i)=0.
  end where
end do
newsal=newsal+sum(flow,2)

watbdy(1:ifull)=max(newwat,0.)
salbdy(1:ifull)=max(newsal,0.)

!--------------------------------------------------------------------
! Water losses over land basins

! estimate grid box area covered by water
cover=min(0.001*watbdy(1:ifull)/minwater,1.)
! estimate rate that water leaves river into soil
rate=dt/(8.*3600.) ! MJT suggestion
  
! Method for land basins
select case(basinmd)
  case(0)
    ! add water to soil moisture 
    where (any(slope>=-1.E-10,2).or..not.land(1:ifull))
      cover=0.
    end where
    if (nsib==6.or.nsib==7) then
      ! CABLE
      tmpry=watbdy(1:ifull)
      call cableinflow(tmpry,cover,rate)
      newwat(1:ifull)=newwat(1:ifull)+(tmpry-watbdy(1:ifull))*(1.-sigmu(:))
    else
      ! Standard land surface model
      tmpry=watbdy(1:ifull)
      do k=1,ms
        where (land(1:ifull))
          ll(:)=max(sfc(isoilm(:))-wb(:,k),0.)*1000.*zse(k)
          ll(:)=ll(:)*rate*cover(:)
          ll(:)=min(tmpry(:),ll(:))
          wb(:,k)=wb(:,k)+ll(:)/(1000.*zse(k))
          tmpry(:)=tmpry(:)-ll(:)
        end where
      end do
      newwat(1:ifull)=newwat(1:ifull)+(tmpry-watbdy(1:ifull))*(1.-sigmu(:))
    end if
  case(2)
    ! pile-up water
  case(3)
    ! leak
    where (.not.land(1:ifull))
      cover=0.
    end where
    if (nsib==6.or.nsib==7) then
      ! CABLE
      tmpry=watbdy(1:ifull)
      call cableinflow(tmpry,cover,rate)
      newwat(1:ifull)=newwat(1:ifull)+(tmpry-watbdy(1:ifull))*(1.-sigmu(:))
    else
      ! Standard land surface model
      tmpry=watbdy(1:ifull)
      do k=1,ms
        where (land(1:ifull))
          ll(:)=max(sfc(isoilm(:))-wb(:,k),0.)*1000.*zse(k)
          ll(:)=ll(:)*rate*cover(:)
          ll(:)=min(tmpry(:),ll(:))
          wb(:,k)=wb(:,k)+ll(:)/(1000.*zse(k))
          tmpry(:)=tmpry(:)-ll(:)
        end where
      end do
      newwat(1:ifull)=newwat(1:ifull)+(tmpry-watbdy(1:ifull))*(1.-sigmu(:))
    end if
  case default
    write(6,*) "ERROR: Unsupported basinmd ",basinmd
    call ccmpi_abort(-1)
end select

watbdy(1:ifull)=max(newwat,0.)
salbdy(1:ifull)=max(newsal,0.)

! case when water is absorbed by land-surface, but there is
! non-zero salinity left behind.
if (any(salbdy(1:ifull)>0..and.watbdy(1:ifull)==0.)) then
  write(6,*) "WARN: Patch river salinity"
end if

!--------------------------------------------------------------------
! update ocean
if (abs(nmlo)<=9) then
  where (ee(1:ifull)>0.5)
    depdum=max(dd(1:ifull)+neta(1:ifull),minwater)
    ! salbdy is already multiplied by watbdy
    sal=(salin*depdum+0.001*salbdy(1:ifull))/(depdum+0.001*watbdy(1:ifull))
    neta(1:ifull)=neta(1:ifull)+0.001*watbdy(1:ifull)
    watbdy(1:ifull)=0.
    salbdy(1:ifull)=sal   ! ocean default
  elsewhere (watbdy(1:ifull)>1.E-10)
    salbdy(1:ifull)=salbdy(1:ifull)/watbdy(1:ifull) ! rescale salinity back to PSU
    sal=0.
  elsewhere
    salbdy(1:ifull)=0.    ! land default
    sal=0.
  end where

  ! import ocean data
  neta(1:ifull)=neta(1:ifull)+inflowbias*ee(1:ifull)
  call mloimport(4,neta(1:ifull),0,0)
  do ii=1,wlev
    where (ee(1:ifull)>0.5.and.salin>0.)
      ! rescale salinity profile
      sallvl(:,ii)=sallvl(:,ii)*sal/salin
    elsewhere (ee(1:ifull)>0.5)
      ! set new uniform salinity profile
      sallvl(:,ii)=sal
    end where
  end do
  call mloimport3d(1,sallvl,0)
end if

return
end subroutine rvrrouter

end module river
