module arrays_m

implicit none

private
public t,u,v,qg
public psl,ps,zs
public arrays_init,arrays_end

real, dimension(:,:), allocatable, save :: t,u,v,qg
real, dimension(:), allocatable, save :: psl,ps,zs

contains

subroutine arrays_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(t(ifull+iextra,kl),u(ifull+iextra,kl),v(ifull+iextra,kl),qg(ifull+iextra,kl))
allocate(psl(ifull+iextra),ps(ifull+iextra),zs(ifull+iextra))

return
end subroutine arrays_init

subroutine arrays_end

implicit none

deallocate(t,u,v,qg)
deallocate(psl,ps,zs)

return
end subroutine arrays_end

end module arrays_m