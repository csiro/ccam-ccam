module o3amip_m

implicit none

private
public gdat,glat,dp,gpri
public mo,jg,kg,lg
public o3amip_init,o3amip_end

integer, parameter :: mo=12    ! Months
integer, parameter :: jg=64    ! Latitudes
integer, parameter :: kg=59    ! Levels
integer, parameter :: lg=kg+1  ! Layer Interfaces
real, dimension(:), allocatable, save :: glat,dp,gpri
real, dimension(:,:,:), allocatable, save :: gdat

contains

subroutine o3amip_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(glat(jg),dp(kg),gpri(lg),gdat(jg,kg,mo))

return
end subroutine o3amip_init

subroutine o3amip_end

implicit none

deallocate(gdat,glat,dp,gpri)

return
end subroutine o3amip_end

end module o3amip_m