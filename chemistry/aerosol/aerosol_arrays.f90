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

module aerosol_arrays

private
public naero
public aldrinit,aldrend,aldrloademiss,aldrloaderod
public xtg,xtgsav,xtosav
public itracdu,ndust
public dustdd,dustwd,duste,dust_burden
public itracbc,bce,bcdd,bcwd,bc_burden
public itracoc,oce,ocdd,ocwd,oc_burden
public itracdms,itracso2,itracso4
public dmse,dmsso2o,so2e,so2so4o,so2dd,so2wd,so4e,so4dd,so4wd
public dms_burden,so2_burden,so4_burden
public itracsa,nsalt,salte,saltdd,saltwd,salt_burden
public xtg_solub,zoxidant_g,erod,ndcls,emissfield,vso2

real, dimension(:,:,:), allocatable, save :: xtg            ! prognostic aerosols (see indexing below)
real, dimension(:,:,:), allocatable, save :: xtgsav         ! save for mass conservation in semi-Lagrangian models
real, dimension(:,:,:), allocatable, save :: xtosav         ! aerosol mixing ratio outside convective cloud
real, dimension(:,:,:), allocatable, save :: xtg_solub      ! aerosol mixing ratio that is dissolved in rain
real, dimension(:,:), allocatable, save :: erod             ! sand, clay and silt fraction that can erode
real, dimension(:,:), allocatable, save :: emissfield       ! non-volcanic emissions
real, dimension(:,:,:), allocatable, save :: zoxidant_g     ! oxidant fields
real, dimension(:), allocatable, save :: vso2               ! volcanic emissions
real, dimension(:,:), allocatable, save :: duste            ! Diagnostic - dust emissions
real, dimension(:,:), allocatable, save :: dustdd           ! Diagnostic - dust dry deposition
real, dimension(:,:), allocatable, save :: dustwd           ! Diagnostic - dust wet deposition
real, dimension(:,:), allocatable, save :: dust_burden      ! Diagnostic - dust burden
real, dimension(:), allocatable, save :: bce                ! Diagnostic - black carbon emissions
real, dimension(:), allocatable, save :: bcdd               ! Diagnostic - black carbon dry deposition
real, dimension(:), allocatable, save :: bcwd               ! Diagnostic - black carbon wet deposition
real, dimension(:), allocatable, save :: bc_burden          ! Diagnostic - black carbon burden
real, dimension(:), allocatable, save :: oce                ! Diagnostic - organic carbon emissions
real, dimension(:), allocatable, save :: ocdd               ! Diagnostic - organic carbon dry deposition
real, dimension(:), allocatable, save :: ocwd               ! Diagnostic - organic carbon wet deposition
real, dimension(:), allocatable, save :: oc_burden          ! Diagnostic - organic carbon burden
real, dimension(:), allocatable, save :: dmse               ! Diagnostic - DMS emissions
real, dimension(:), allocatable, save :: dmsso2o            ! Diagnostic - DMS->so2 oxidation
real, dimension(:), allocatable, save :: so2e               ! Diagnostic - so2 emissions
real, dimension(:), allocatable, save :: so2so4o            ! Diagnostic - so2->so4 oxidation
real, dimension(:), allocatable, save :: so2dd              ! Diagnostic - so2 dry deposition
real, dimension(:), allocatable, save :: so2wd              ! Diagnostic - so2 wet deposition
real, dimension(:), allocatable, save :: so4e               ! Diagnostic - so4 emissions
real, dimension(:), allocatable, save :: so4dd              ! Diagnostic - so4 dry deposition
real, dimension(:), allocatable, save :: so4wd              ! Diagnostic - so4 wet deposition
real, dimension(:), allocatable, save :: salte              ! Diagnostic - salt emission
real, dimension(:), allocatable, save :: saltdd             ! Diagnostic - salt dry deposition
real, dimension(:), allocatable, save :: saltwd             ! Diagnostic - salt wet deposition
real, dimension(:), allocatable, save :: dms_burden         ! Diagnostic - DMS burden
real, dimension(:), allocatable, save :: so2_burden         ! Diagnostic - so2 burden
real, dimension(:), allocatable, save :: so4_burden         ! Diagnostic - so4 burden
real, dimension(:), allocatable, save :: salt_burden        ! Diagnostic - salt burden

! tracers
integer, parameter :: nsulf = 3
integer, parameter :: ncarb = 4
integer, parameter :: ndust = 4
integer, parameter :: nsalt = 2
integer, parameter :: naero = nsulf+ncarb+ndust+nsalt ! Tracers: DMS, SO2, SO4, BCO, BCI, OCO, OCI, DUST(4), SALT(2)
integer, parameter :: itracdms = 1                  ! Index for DMS tracer (1)
integer, parameter :: itracso2 = 2                  ! Index for SO2 tracer (2)
integer, parameter :: itracso4 = 3                  ! Index for SO4 tracer (3)
integer, parameter :: itracbc = nsulf+1             ! Index for BC tracer (hydrophobic (4), hydrophillic (5))
integer, parameter :: itracoc = nsulf+3             ! Index for OC tracer (hydrophobic (6), hydrophillic (7))
integer, parameter :: itracdu = nsulf+ncarb+1       ! Index for dust tracer (8-11)
integer, parameter :: itracsa = nsulf+ncarb+ndust+1 ! Index for salt tracer (12-13)
integer, parameter :: ndcls = 3                     ! Number of dust emission classes (sand, silt, clay)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialisation

subroutine aldrinit(ifull,iextra,kl,sig)

implicit none

integer, intent(in) :: ifull,iextra,kl
integer k
real, dimension(kl), intent(in) :: sig

allocate(xtg(ifull+iextra,kl,naero),xtgsav(ifull,kl,naero))
allocate(xtosav(ifull,kl,naero),vso2(ifull))
allocate(emissfield(ifull,15))
allocate(zoxidant_g(ifull,kl,4),erod(ifull,ndcls))
allocate(duste(ifull,ndust),dustdd(ifull,ndust),dustwd(ifull,ndust),dust_burden(ifull,ndust))
allocate(bce(ifull),bcdd(ifull),bcwd(ifull))
allocate(bc_burden(ifull))
allocate(oce(ifull),ocdd(ifull),ocwd(ifull))
allocate(oc_burden(ifull))
allocate(dmse(ifull),dmsso2o(ifull))
allocate(so2e(ifull),so2so4o(ifull),so2dd(ifull),so2wd(ifull))
allocate(so4e(ifull),so4dd(ifull),so4wd(ifull))
allocate(dms_burden(ifull),so2_burden(ifull),so4_burden(ifull))
allocate(salte(ifull),saltdd(ifull),saltwd(ifull),salt_burden(ifull))

xtg=0.
xtgsav=0.
xtosav=0.
vso2=0.
emissfield=0.
zoxidant_g=0.
erod=0.
duste=0.
dustdd=0.
dustwd=0.
dust_burden=0.
bce=0.
bcdd=0.
bcwd=0.
bc_burden=0.
oce=0.
ocdd=0.
ocwd=0.
oc_burden=0.
dmse=0.
dmsso2o=0.
so2e=0.
so2so4o=0.
so2dd=0.
so2wd=0.
so4e=0.
so4dd=0.
so4wd=0.
salte=0.
saltdd=0.
saltwd=0.
dms_burden=0.
so2_burden=0.
so4_burden=0.
salt_burden=0.

return
end subroutine aldrinit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! End

subroutine aldrend

implicit none

deallocate(xtg,xtgsav,xtosav)
deallocate(vso2)
deallocate(emissfield)
deallocate(zoxidant_g,erod)
deallocate(duste,dustdd,dustwd,dust_burden)
deallocate(bce,bcdd,bcwd)
deallocate(bc_burden)
deallocate(oce,ocdd,ocwd)
deallocate(oc_burden)
deallocate(dmse,dmsso2o)
deallocate(so2e,so2so4o,so2dd,so2wd)
deallocate(so4e,so4dd,so4wd)
deallocate(dms_burden,so2_burden,so4_burden)
deallocate(salte,saltdd,saltwd,salt_burden)

return
end subroutine aldrend

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Load emission arrays

subroutine aldrloademiss(index,aa)

implicit none

integer, intent(in) :: index
real, dimension(:), intent(in) :: aa

if ( index<16 ) then
  emissfield(:,index)=aa(1:size(emissfield,1)) ! Then follow SO2, BC and OC from anthro (a) and biomass-burning (b) levels 1 and 2
elseif ( index==16 ) then
  vso2(:)=aa(1:size(vso2))             ! volcanic
else
  write(6,*) "ERROR: index out-of-range for aldrloademiss"
  stop
end if

#ifdef debug
if ( maxval(emissfield(:,14))>1. ) then
  write(6,*) "ERROR out-of-range dmst in aldrloademiss ",maxval(emissfield(:,14))
  stop -1
end if
#endif

return
end subroutine aldrloademiss

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Load soil data

subroutine aldrloaderod(inda,aa)

implicit none

integer, intent(in) :: inda
real, dimension(:), intent(in) :: aa

! EROD is the soil fraction of Sand (inda=1), Silt (inda=2) and Clay (inda=3) that can erode
erod(:,inda) = aa(1:size(erod,1))

return
end subroutine aldrloaderod

end module aerosol_arrays