module module_aux_rad

implicit none

private
public liqradmethod, iceradmethod
public cloud3

integer, save :: liqradmethod = 0     ! Method for calculating radius of liquid droplets
                                      ! (0=Martin Mid/NBer, 1=Martin Low/NBer, 2=Martin High/NBer
                                      !  3=Martin Mid/Ber,  4=Martin Low/Ber,  5=Martin High/Ber )
integer, save :: iceradmethod = 1     ! Method for calculating radius of ice droplets
                                      ! (0=Lohmann, 1=Donner smooth, 2=Fu, 3=Donner orig)

contains

! This subroutine is based on cloud2.f
subroutine cloud3(Rdrop,Rice,conl,coni,cfrac,qlrad,qfrad,prf,ttg,cdrop,imax,kl)

use cc_mpi              ! CC MPI routines
use const_phys          ! Physical constants
use parm_m              ! Model configuration

implicit none

integer, intent(in) :: imax, kl
integer iq, k, kr
integer, dimension(imax) :: tpos
real, dimension(imax,kl), intent(in) :: cfrac, qlrad, qfrad, prf, ttg
real, dimension(imax,kl), intent(in) :: cdrop
real(kind=8), dimension(imax,kl), intent(out) :: Rdrop, Rice, conl, coni
real, dimension(imax,kl) :: reffl, reffi, Wliq, rhoa
real, dimension(imax,kl) :: eps, rk, Wice
real, dimension(imax) :: basesize, tstore, tfrac
real :: liqparm
! parameters
logical :: do_brenguier   ! Adjust effective radius for vertically stratified cloud
real :: scale_factor      ! account for the plane-parallel homogenous cloud bias  (e.g. Cahalan effect)
! data
! MJT notes - add one extra index for overflow
integer, parameter :: n_donner = 8
real, dimension(n_donner+1), parameter :: temp_donner = &
    (/ 215.66, 220.66, 225.66, 230.66, 235.66, 240.66, 245.66, 250.66, 250.66 /)
real, dimension(n_donner+1), parameter :: rad_donner =  &
    (/   20.2,   21.6,   39.9,   42.5,   63.9,   93.5,   80.8,  100.6,  100.6 /)

real, parameter :: rhow     = 1000.            ! Density of water (kg/m^3) 

scale_factor = 1.
liqparm = -0.003e-6
do_brenguier = .false.

rhoa(:,:) = prf(:,:)/(rdry*ttg(:,:))

select case(liqradmethod)
  case(0)
    liqparm = -0.003e-6 ! mid range
    do_brenguier = .false.
    scale_factor = 1.
  case(1)
    liqparm = -0.001e-6 ! lower bound
    do_brenguier = .false.
    scale_factor = 1.
  case(2)
    liqparm = -0.008e-6 ! upper bound
    do_brenguier = .false.
    scale_factor = 1.
  case(3)
    liqparm = -0.003e-6 ! mid range
    do_brenguier = .true.
    scale_factor = 0.85
  case(4)
    liqparm = -0.001e-6 ! lower bound
    do_brenguier = .true.
    scale_factor = 0.85
  case(5)
    liqparm = -0.008e-6 ! upper bound
    do_brenguier = .true.
    scale_factor = 0.85
  case default
    write(6,*) "ERROR: Unknown option for liqradmethod ",liqradmethod
    call ccmpi_abort(-1)
end select


! Reffl is the effective radius calculated following
! Martin etal 1994, JAS 51, 1823-1842
where ( qlrad(:,:)>1.E-8 .and. cfrac(:,:)>1.E-4 )
  Wliq(:,:) = rhoa(:,:)*qlrad(:,:)/cfrac(:,:) !kg/m^3
  ! This is the Liu and Daum scheme for relative dispersion (Nature, 419, 580-581 and pers. comm.)
  eps(:,:) = 1. - 0.7*exp(liqparm*cdrop(:,:))
  rk(:,:)  = (1.+eps(:,:)**2)/(1.+2.*eps(:,:)**2)**2

  ! k_ratio = rk**(-1./3.)
  ! GFDL        k_ratio (land) 1.143 (water) 1.077
  ! mid range   k_ratio (land) 1.393 (water) 1.203
  ! lower bound k_ratio (land) 1.203 (water) 1.050

  ! Martin et al 1994
  reffl(:,:) = (2.*3.*Wliq(:,:)/(4.*pi*rhow*rk(:,:)*cdrop(:,:)))**(1./3.)
  !qlpath = Wliq*dz(iq,k)
  !taul(iq,k) = tau_sfac*1.5*qlpath/(rhow*Reffl)
elsewhere
  reffl(:,:) = 0.
  Wliq(:,:) = 0.
end where


! (GFDL NOTES)
!    for single layer liquid or mixed phase clouds it is assumed that
!    cloud liquid is vertically stratified within the cloud.  under
!    such situations for observed stratocumulus clouds it is found
!    that the cloud mean effective radius is between 80 and 100% of
!    the cloud top effective radius. (Brenguier et al., Journal of
!    Atmospheric Sciences, vol. 57, pp. 803-821 (2000))  for linearly
!    stratified cloud in liquid specific humidity, the cloud top
!    effective radius is greater than the effective radius of the
!    cloud mean specific humidity by a factor of 2**(1./3.).
!    this correction, 0.9*(2**(1./3.)) = 1.134, is applied only to
!    single layer liquid or mixed phase clouds.
if ( do_brenguier ) then
  if ( nmr>=1 ) then
    ! Max-Rnd overlap
    where ( cfrac(:,2)<1.e-4 )
      reffl(:,1) = reffl(:,1)*1.134
    end where
    do k = 2,kl-1
      where ( cfrac(:,k-1)<1.e-4 .and. cfrac(:,k+1)<1.e-4 )
        reffl(:,k) = reffl(:,k)*1.134
      end where
    end do
    where ( cfrac(:,kl-1)<1.e-4 )
      reffl(:,kl) = reffl(:,kl)*1.134
    end where
  else
    ! Rnd overlap
    reffl(:,:) = reffl(:,:)*1.134
  end if
end if

select case(iceradmethod)
  case(0)
    !Lohmann et al.(1999)
    where ( qfrad(:,:)>1.E-8 .and. cfrac(:,:)>1.E-4 )
      Wice(:,:) = rhoa(:,:)*qfrad(:,:)/cfrac(:,:) !kg/m**3
      reffi(:,:) = 0.5*min(150.e-6, 3.73e-4*Wice(:,:)**0.216)
    elsewhere
      Wice(:,:) = 0.
      reffi(:,:) = 0.
    end where

  case(1)
    !Donner et al (1997)
    ! linear interpolation by MJT
    do k = 1,kl
      tstore(:) = (ttg(:,k) - temp_donner(1))/(temp_donner(2) - temp_donner(1)) + 1.
      tstore(:) = min( max( tstore(:), 1. ), real(n_donner) )
      tfrac(:) = tstore(:) - aint(tstore(:))
      tpos(:) = int(tstore)
      basesize(:) = (1.-tfrac(:))*rad_donner(tpos(:)) + tfrac(:)*rad_donner(tpos(:)+1)
      where ( qfrad(:,k)>1.e-8 .and. cfrac(:,k)>1.e-4 )
        Wice(:,k) = rhoa(:,k)*qfrad(:,k)/cfrac(:,k) ! kg/m**3
        reffi(:,k) = 1.e-6*basesize(:) ! 5.e-7*basesize(:)
      elsewhere
        Wice(:,k) = 0.
        reffi(:,k) = 0.
      end where
    end do

  case(2)
    ! Fu 2007
    where ( qfrad(:,:)>1.E-8 .and. cfrac(:,:)>1.E-4 )
      Wice(:,:) = rhoa(:,:)*qfrad(:,:)/cfrac(:,:) !kg/m**3
      reffi(:,:) = 5.E-7*(47.05+0.6624*(ttg(:,:)-273.16)+0.001741*(ttg(:,:)-273.16)**2)
    elsewhere
      Wice(:,:) = 0.
      reffi(:,:) = 0.
    end where

  case(3)
    do k = 1,kl
      do iq = 1,imax
        if ( qfrad(iq,k)>1.E-8 .and. cfrac(iq,k)>1.E-4 ) then
          Wice(iq,k) = rhoa(iq,k)*qfrad(iq,k)/cfrac(iq,k) ! kg/m**3
          if ( ttg(iq,k)>248.16 ) then
            reffi(iq,k) = 5.E-7*100.6
          elseif ( ttg(iq,k)>243.16 ) then
            reffi(iq,k) = 5.E-7*80.8
          elseif ( ttg(iq,k)>238.16 ) then
            reffi(iq,k) = 5.E-7*93.5
          elseif ( ttg(iq,k)>233.16 ) then
            reffi(iq,k) = 5.E-7*63.9
          elseif ( ttg(iq,k)>228.16 ) then
            reffi(iq,k) = 5.E-7*42.5
          elseif ( ttg(iq,k)>223.16 ) then
            reffi(iq,k) = 5.E-7*39.9
          elseif ( ttg(iq,k)>218.16 ) then
            reffi(iq,k) = 5.E-7*21.6
          else
            reffi(iq,k) = 5.E-7*20.2
          end if
        else
          reffi(iq,k) = 0.
          Wice(iq,k) = 0.
        end if
      end do
    end do

end select

do k = 1,kl
  kr = kl + 1 - k
  Rdrop(:,kr) = real(1.E6*reffl(:,k), 8) ! convert to diameter and microns
  Rice(:,kr)  = real(1.E6*reffi(:,k), 8)
  conl(:,kr)  = real(1000.*scale_factor*Wliq(:,k), 8) !g/m^3
  coni(:,kr)  = real(1000.*scale_factor*Wice(:,k), 8)
end do


Rdrop(:,:) = min(max(Rdrop(:,:), 8.4_8), 33.2_8) ! constrain diameter to acceptable range (see microphys_rad.f90)
Rice(:,:) = min(max(Rice(:,:), 18.6_8), 130.2_8)

return
end subroutine cloud3
end module module_aux_rad 































































































