 subroutine gettin(n)

!  gettina reads back preceding u,v,t,psl for n=2 (for nmi, assuming mex=1)
!  saves and re-reads initial arrays of t and psl
!  called only by darlam and vmodes
    
use arrays_m
use savuvt_m

implicit none

include 'newmpar.h'
      
integer, intent(in) :: n
      
if ( n==0 ) then
  savt(:,:)=t(1:ifull,:)
  savpsl(:)=psl(1:ifull)
else if ( n==2 ) then
  t(1:ifull,:)=savt(:,:)   ! for n=1, n=2
  psl(1:ifull)=savpsl(:)   ! for n=1, n=2
  u(1:ifull,:)=savu(:,:)  ! only for n=2 (VMI init)
  v(1:ifull,:)=savv(:,:)
else    ! for n=1
  t(1:ifull,:)=savt(:,:)   ! for n=1, n=2
  psl(1:ifull)=savpsl(:)   ! for n=1, n=2
endif

return
end subroutine gettin
