module soil_m

implicit none

private
public zoland,zmin,zolnd,zolog,albsav,so4t,albnirsav
public albvisdif,albnirdif,albvisdir,albnirdir
public land
public soil_init,soil_end

real, dimension(:), allocatable, save :: zolnd,zolog,albsav,so4t,albnirsav
real, dimension(:), allocatable, save :: albvisdif,albnirdif,albvisdir,albnirdir
real, save :: zoland,zmin
logical, dimension(:), allocatable, save :: land

contains

subroutine soil_init(ifull,iextra,kl,iaero,nsib)

implicit none

integer, intent(in) :: ifull,iextra,kl,iaero,nsib

allocate(zolnd(ifull),albsav(ifull),albnirsav(ifull))
allocate(albvisdif(ifull),albnirdif(ifull),albvisdir(ifull),albnirdir(ifull))
allocate(land(ifull))
if (iaero.ne.0) then
  allocate(so4t(ifull))
end if
if (nsib.eq.3.or.nsib.eq.5) then
  allocate(zolog(ifull))
end if
zoland=0.16

return
end subroutine soil_init

subroutine soil_end

implicit none

deallocate(zolnd,zolog,albsav,albnirsav)
deallocate(albvisdif,albnirdif,albvisdir,albnirdir)
deallocate(land)
if (allocated(so4t)) deallocate(so4t)
if (allocated(zolog)) deallocate(zolog)

return
end subroutine soil_end

end module soil_m
