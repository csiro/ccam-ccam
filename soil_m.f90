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

subroutine soil_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(zolnd(ifull),zolog(ifull),albsav(ifull),so4t(ifull),albnirsav(ifull))
allocate(albvisdif(ifull),albnirdif(ifull),albvisdir(ifull),albnirdir(ifull))
allocate(land(ifull))
zoland=0.16

return
end subroutine soil_init

subroutine soil_end

implicit none

deallocate(zolnd,zolog,albsav,so4t,albnirsav)
deallocate(albvisdif,albnirdif,albvisdir,albnirdir)
deallocate(land)

return
end subroutine soil_end

end module soil_m