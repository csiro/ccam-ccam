module kuocomb_m

implicit none

private
public kbsav,ktsav
public convpsav
public kuocomb_init,kuocomb_end

integer, dimension(:), allocatable, save :: kbsav,ktsav
real, dimension(:), allocatable, save :: convpsav

contains

subroutine kuocomb_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(kbsav(ifull),ktsav(ifull))
allocate(convpsav(ifull))

return
end subroutine kuocomb_init

subroutine kuocomb_end

implicit none

deallocate(kbsav,ktsav)
deallocate(convpsav)

return
end subroutine kuocomb_end

end module kuocomb_m