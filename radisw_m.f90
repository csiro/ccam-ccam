module radisw_m

implicit none

private
public ktop,kbtm,nclds,ktopsw,kbtmsw,emcld
public temp,temp2,press,press2,rh2o,qo3
public camt,cuvrf,cirrf,cirab,coszro,taudar
public rrco2,ssolar,rrvco2,rrvch4,rrvn2o
public rrvf11,rrvf12,rrvf113,rrvf22
public radisw_init,radisw_end

integer, dimension(:), allocatable, save :: nclds
integer, dimension(:,:), allocatable, save :: ktop,kbtm,ktopsw,kbtmsw
real, dimension(:), allocatable, save :: coszro,taudar
real, dimension(:,:), allocatable, save :: emcld,camt,cuvrf,cirrf,cirab
real, dimension(:,:), allocatable, save :: temp,temp2,press,press2,rh2o,qo3
real, save :: rrco2,ssolar,rrvco2,rrvch4,rrvn2o
real, save :: rrvf11,rrvf12,rrvf113,rrvf22

contains

subroutine radisw_init(ifull,iextra,kl,imax)

implicit none

integer, intent(in) :: ifull,iextra,kl,imax

allocate(ktop(imax,kl+1),kbtm(imax,kl+1),nclds(imax),ktopsw(imax,kl+1),kbtmsw(imax,kl+1),emcld(imax,kl+1))
allocate(temp(imax,kl+1),temp2(imax,kl+1),press(imax,kl+1),press2(imax,kl+1),rh2o(imax,kl),qo3(imax,kl))
allocate(camt(imax,kl+1),cuvrf(imax,kl+1),cirrf(imax,kl+1),cirab(imax,kl+1),coszro(imax),taudar(imax))

return
end subroutine radisw_init

subroutine radisw_end

implicit none

deallocate(ktop,kbtm,nclds,ktopsw,kbtmsw,emcld)
deallocate(temp,temp2,press,press2,rh2o,qo3)
deallocate(camt,cuvrf,cirrf,cirab,coszro,taudar)

return
end subroutine radisw_end

end module radisw_m