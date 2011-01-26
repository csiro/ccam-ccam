module tabcom_m

implicit none

private
public em1,em1wde,table1,table2,table3,em3
public source,dsrce,ind,indx2,kmaxv,kmaxvm
public tabcom_init,tabcom_end

real, dimension(:,:), allocatable, save :: em1,em1wde,table1,table2,table3,em3
real, dimension(:,:), allocatable, save :: source,dsrce
integer, dimension(:), allocatable, save :: ind,indx2,kmaxv
integer, save :: kmaxvm

contains

subroutine tabcom_init(ifull,iextra,kl,imax,nbly)

implicit none

integer :: lp1v

integer, intent(in) :: ifull,iextra,kl,imax,nbly

lp1v=(kl+1)*(1+2*kl/2)
allocate(em1(28,180),em1wde(28,180),table1(28,180),table2(28,180),table3(28,180),em3(28,180))
allocate(source(28,nbly),dsrce(28,nbly),ind(imax),indx2(lp1v),kmaxv(kl+1))

return
end subroutine tabcom_init

subroutine tabcom_end

implicit none

deallocate(em1,em1wde,table1,table2,table3,em3)
deallocate(source,dsrce,ind,indx2,kmaxv)

return
end subroutine tabcom_end

end module tabcom_m