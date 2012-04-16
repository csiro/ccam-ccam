module vegpar_m

implicit none

private
public cansto,vlai
public vegpar_init,vegpar_end

real, dimension(:), allocatable, save :: cansto,vlai

contains

subroutine vegpar_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(cansto(ifull),vlai(ifull))

return
end subroutine vegpar_init

subroutine vegpar_end

implicit none

deallocate(cansto,vlai)

return
end subroutine vegpar_end

end module vegpar_m