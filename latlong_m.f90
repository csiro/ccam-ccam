module latlong_m

implicit none

private
public rlatt_g,rlongg_g
public rlatt, rlongg
public latlong_init,latlong_end

real, dimension(:), allocatable, save :: rlatt_g,rlongg_g
real, dimension(:), allocatable, save :: rlatt,rlongg

contains

subroutine latlong_init(ifull_g,ifull,iextra,myid)

implicit none

integer, intent(in) :: ifull_g,ifull,iextra,myid

if (myid==0) then
  allocate(rlatt_g(ifull_g),rlongg_g(ifull_g))
end if
allocate(rlatt(ifull),rlongg(ifull))

return
end subroutine latlong_init

subroutine latlong_end

implicit none

if (allocated(rlatt_g)) deallocate(rlatt_g)
if (allocated(rlongg_g)) deallocate(rlongg_g)
deallocate(rlatt,rlongg)

return
end subroutine latlong_end

end module latlong_m