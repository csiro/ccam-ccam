module vegpar_m

implicit none

private
public cansto,vlai,fwet
public vegpar_init,vegpar_end

real, dimension(:), allocatable, save :: cansto,vlai,fwet

contains

subroutine vegpar_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(cansto(ifull),vlai(ifull),fwet(ifull))
cansto=0.
vlai=0.
fwet=0.

return
end subroutine vegpar_init

subroutine vegpar_end

implicit none

deallocate(cansto,vlai,fwet)

return
end subroutine vegpar_end

end module vegpar_m