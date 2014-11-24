subroutine mslp(pmsl,psl,zs,t)

use cc_mpi, only : mydiag
use sigs_m

implicit none
! this one will ignore negative zs (i.e. over the ocean)

include 'newmpar.h'
include 'const_phys.h'
include 'parm.h'

integer, parameter :: meth=1 ! 0 for original, 1 for other jlm - always now
integer, save :: lev = -1
real c,conr,con
real, dimension(ifull), intent(out) :: pmsl
real, dimension(ifull), intent(in) :: psl,zs
real, dimension(ifull) :: phi1,tsurf,tav,dlnps
real, dimension(:,:), intent(in) :: t
      
c=grav/stdlapse
conr=c/rdry
if (lev<0) then
  lev=1
  do while (sig(lev+1)<=0.9)
    lev=lev+1
  end do
end if
con=sig(lev)**(rdry/c)/c
      
if(meth==1)then
  phi1(:)=t(1:ifull,lev)*rdry*(1.-sig(lev))/sig(lev) ! phi of sig(lev) above sfce
  tsurf(:)=t(1:ifull,lev)+phi1(:)*stdlapse/grav
  tav(:)=tsurf(:)+zs(1:ifull)*.5*stdlapse/grav
  dlnps(:)=zs(1:ifull)/(rdry*tav(:))
  pmsl(:)=1.e5*exp(psl(:)+dlnps(:))
endif  ! (meth==1)
      
if(nmaxpr==1.and.mydiag)then
  write(6,*) 'meth,lev,sig(lev) ',meth,lev,sig(lev)
  write(6,*) 'zs,t_lev,psl,pmsl ',zs(idjd),t(idjd,lev),psl(idjd),pmsl(idjd)
endif
      
return
end subroutine mslp
      
