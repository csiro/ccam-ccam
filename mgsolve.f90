module mgsolve

! This module solves for the atmosphere helmholtz equation and the ocean free
! surface equation using a multi-grid approach.

! Currently we simply aggregate data as we upscale the grid, untill eventially
! all processors solve for the same coarse grid (thereby avoiding MPI_Allreduce)
! when solving the coarse grid.

implicit none

private
public mghelm,mgmlo,mgsor_init,mgzz_init

integer, save :: mg_maxlevel
integer, save :: mg_maxsize
integer, dimension(3), save :: mg_ifullc
integer, dimension(:,:), allocatable, save :: col_iq,col_iqn,col_iqe,col_iqs,col_iqw
logical, save :: sorfirst=.true.
logical, save :: zzfirst=.true.
logical, save :: mlofirst=.true.

type mgtype
  integer ifull,ifull_fine,ifull_coarse
  integer merge_len,merge_row,ipan
  integer comm,comm_mlo
  integer, dimension(:,:,:), allocatable :: fproc
  integer, dimension(:,:), allocatable :: merge_list
  integer, dimension(:), allocatable :: merge_pos
  integer, dimension(:), allocatable :: in,ie,is,iw,ine,inw,ise,isw
  integer, dimension(:), allocatable :: coarse_a,coarse_b,coarse_c,coarse_d
  integer, dimension(:), allocatable :: fine
  real, dimension(:), allocatable :: zzn,zze,zzs,zzw,zz
  logical globgath
end type mgtype

type mgbndtype
  integer(kind=4) len
  integer(kind=4) rlen,rlenx
  integer(kind=4) slen,slenx
  integer(kind=4), dimension(:), allocatable :: send_list
  integer(kind=4), dimension(:), allocatable :: unpack_list
  integer(kind=4), dimension(:), allocatable :: request_list
  real, dimension(:), allocatable :: rbuf
  real, dimension(:), allocatable :: sbuf
end type mgbndtype

type(mgtype), dimension(:), allocatable, save :: mg
type(mgbndtype), dimension(:,:), allocatable, save :: mg_bnds

integer, parameter :: itr_max=300 ! maximum number of iterations

contains

! This version is for the atmosphere
subroutine mghelm(izz,izzn,izze,izzw,izzs,ihelm,iv,irhs)

use cc_mpi
use indices_m

implicit none

include 'newmpar.h'
include 'parm.h'
include 'parmdyn.h'

integer, dimension(mg_maxsize) :: iqaa,iqbb,iqcc,iqdd
integer, dimension(kl) :: iters
integer itrc,itr,ng,ng4,g,k,jj
integer klimc,knew,klim,ir,ic
integer nc,ifc,n,iq_a,iq_b,iq_c,iq_d
!integer itsave2,itsave1,itstest,itc
real, dimension(ifull+iextra,kl), intent(inout) :: iv
real, dimension(ifull,kl), intent(in) :: ihelm,irhs
real, dimension(ifull), intent(in) :: izz,izzn,izze,izzw,izzs
real, dimension(ifull,maxcolour,kl) :: helmc,rhsc
real, dimension(ifull,maxcolour) :: zzc,zznc,zzec,zzwc,zzsc
real, dimension(ifull,3) :: zzcu,zzncu,zzecu,zzwcu,zzscu
real, dimension(mg_maxsize,kl,mg_maxlevel) :: v
real, dimension(mg_maxsize,kl,mg_maxlevel) :: rhs,helm
real, dimension(mg_maxsize,kl) :: w,dsol
real, dimension(mg_maxsize) :: vnew
real, dimension(kl) :: smax,smax_g,dsolmax,dsolmax_g
real, dimension(kl) :: smin,smin_g
real gd,ci,itserr2,itserr1

call start_log(helm_begin)

if (sorfirst.or.zzfirst) then
  write(6,*) "ERROR: mghelm requires mgsor_init and mgzz_init to be called first"
  call ccmpi_abort(-1)
end if


! Prepare input arrays
klim=kl
do k=1,kl
  smax(k)=maxval(iv(1:ifull,k))
  smin(k)=minval(iv(1:ifull,k))
end do
call ccmpi_allreduce(smax(1:kl),smax_g(1:kl),"max",comm_world)
call ccmpi_allreduce(smin(1:kl),smin_g(1:kl),"min",comm_world)

iters=0

!itsave2=0
!itserr2=9.E9
!itstest=1
!itc=0

rhs(1:ifull,:,1)=irhs(1:ifull,:)
helm(1:ifull,:,1)=ihelm(1:ifull,:)
! For when the inital grid cannot be upscaled
if (mg(1)%merge_len>1) then
  call mgcollect(1,rhs(:,:,1))
  call mgcollect(1,helm(:,:,1))
end if

do g=1,mg_maxlevel-1
  ng4=mg(g)%ifull_fine
  iqaa(1:ng4)=          mg(g)%fine    
  iqbb(1:ng4)= mg(g)%in(mg(g)%fine)
  iqcc(1:ng4)= mg(g)%ie(mg(g)%fine)
  iqdd(1:ng4)=mg(g)%ine(mg(g)%fine)
  helm(1:ng4,1:kl,g+1)=0.25*(helm(iqaa(1:ng4),1:kl,g)+helm(iqbb(1:ng4),1:kl,g) &
                            +helm(iqcc(1:ng4),1:kl,g)+helm(iqdd(1:ng4),1:kl,g))
  call mgcollect(g+1,helm(:,:,g+1))
end do

do nc=1,maxcolour
  ifc=ifullx(nc)
  zznc(1:ifc,nc) =izzn(iqx(1:ifc,nc))
  zzwc(1:ifc,nc) =izzw(iqx(1:ifc,nc))
  zzec(1:ifc,nc) =izze(iqx(1:ifc,nc))
  zzsc(1:ifc,nc) =izzs(iqx(1:ifc,nc))
  zzc(1:ifc,nc)  =izz(iqx(1:ifc,nc))
  do k=1,kl
    helmc(1:ifc,nc,k)=ihelm(iqx(1:ifc,nc),k)
    rhsc(1:ifc,nc,k) =irhs(iqx(1:ifc,nc),k)
  end do
end do

g=mg_maxlevel
do nc=1,3
  ifc=mg_ifullc(nc)
  zzncu(1:ifc,nc)=mg(g)%zzn(col_iq(1:ifc,nc))
  zzscu(1:ifc,nc)=mg(g)%zzs(col_iq(1:ifc,nc))
  zzecu(1:ifc,nc)=mg(g)%zze(col_iq(1:ifc,nc))
  zzwcu(1:ifc,nc)=mg(g)%zzw(col_iq(1:ifc,nc))
  zzcu(1:ifc,nc) =mg(g)%zz(col_iq(1:ifc,nc))
end do

call bounds(iv, klim=klim)


! Main loop
do itr=1,itr_max

  ! residual
  do k=1,klim
    w(1:ifull,k)=irhs(:,k)+(ihelm(:,k)-izz)*iv(1:ifull,k)     &
             -izze*iv(ie,k)-izzw*iv(iw,k)                     &
             -izzn*iv(in,k)-izzs*iv(is,k)
  end do

  ! fine grid
  g=1
  ng=mg(g)%ifull

  ! For when the inital grid cannot be upscaled
  call mgcollect(g,w(:,:),klim=klim)

  ! restriction
  ! (since this always operates within a panel, then ine = ien is always true)
  ng4=mg(g)%ifull_fine
  iqaa(1:ng4)=          mg(g)%fine    
  iqbb(1:ng4)= mg(g)%in(mg(g)%fine)
  iqcc(1:ng4)= mg(g)%ie(mg(g)%fine)
  iqdd(1:ng4)=mg(g)%ine(mg(g)%fine)
  rhs(1:ng4,1:klim,g+1)=0.25*(w(iqaa(1:ng4),1:klim)+w(iqbb(1:ng4),1:klim) &
                             +w(iqcc(1:ng4),1:klim)+w(iqdd(1:ng4),1:klim))

  ! upscale grid
  do g=2,mg_maxlevel-1
  
    ng=mg(g)%ifull

    ! merge grids if insufficent points on this processor
    call mgcollect(g,rhs(:,:,g),klim=klim)
    
    ! update
    ! (here we avoid using colours as it requires additional bounds calls)

    ! assume zero for first guess of residual (also avoids additional bounds call)
    !v(:,1:klim,g)=0.
    do k=1,klim
      v(1:ng,k,g)=-rhs(1:ng,k,g)/(helm(1:ng,k,g)-mg(g)%zz)
    end do
    
    ! residual
    call mgbounds(g,v(:,:,g),klim=klim)
    do k=1,klim
      w(1:ng,k)=rhs(1:ng,k,g)+(helm(1:ng,k,g)-mg(g)%zz)*v(1:ng,k,g) &
           -mg(g)%zze*v(mg(g)%ie,k,g)-mg(g)%zzw*v(mg(g)%iw,k,g)   &
           -mg(g)%zzn*v(mg(g)%in,k,g)-mg(g)%zzs*v(mg(g)%is,k,g)
    end do


    ! restriction
    ! (calculate finer grid before mgcollect as the messages sent/recv are shorter)
    ng4=mg(g)%ifull_fine
    iqaa(1:ng4)=          mg(g)%fine
    iqbb(1:ng4)= mg(g)%in(mg(g)%fine)
    iqcc(1:ng4)= mg(g)%ie(mg(g)%fine)
    iqdd(1:ng4)=mg(g)%ine(mg(g)%fine)
    rhs(1:ng4,1:klim,g+1)=0.25*(w(iqaa(1:ng4),1:klim)+w(iqbb(1:ng4),1:klim) &
                               +w(iqcc(1:ng4),1:klim)+w(iqdd(1:ng4),1:klim))

  end do

  ! solve coarse grid
  g=mg_maxlevel
  ng=mg(g)%ifull

  ! ensure all processors have a copy of the coarse grid
  v(1:ng,1:klim,g)=0.
  call mgcollect(g,rhs(:,:,g),klim=klim)

  klimc=klim
  do itrc=1,itr_max
    do nc=1,3 ! always three colours for coarse multigrid
      ifc=mg_ifullc(nc)
      do k=1,klimc
        ! update
        dsol(col_iq(1:ifc,nc),k)=(zzecu(1:ifc,nc)*v(col_iqe(1:ifc,nc),k,g)+zzwcu(1:ifc,nc)*v(col_iqw(1:ifc,nc),k,g)        &
                                 +zzncu(1:ifc,nc)*v(col_iqn(1:ifc,nc),k,g)+zzscu(1:ifc,nc)*v(col_iqs(1:ifc,nc),k,g)        &
                                 -rhs(col_iq(1:ifc,nc),k,g)-(helm(col_iq(1:ifc,nc),k,g)-zzcu(1:ifc,nc))*v(col_iq(1:ifc,nc),k,g)) &
                                /(helm(col_iq(1:ifc,nc),k,g)-zzcu(1:ifc,nc))
        v(col_iq(1:ifc,nc),k,g)=v(col_iq(1:ifc,nc),k,g)+dsol(col_iq(1:ifc,nc),k)
      end do
      ! no call to bounds since all points are on this processor
    end do
    do k=1,klimc  
      smax(k)=maxval(v(1:ng,k,g))
      smin(k)=minval(v(1:ng,k,g))
      dsolmax(k)=maxval(abs(dsol(1:ng,k)))
    end do    ! convergence test
    knew=klimc
    do k=klimc,1,-1
      if (dsolmax(k)>=restol*(smax(k)-smin(k))) exit
      knew=k-1      
    end do
    klimc=knew
    if (klimc<1) exit
  end do
  
  ! interpolation
  ! all points are present on this processor, so no need for MPI comms
  ! (there maybe issues interpolating near vertices, where we could
  !  possibly use a three point interpolation)
  ng4=mg(g)%ifull_coarse
  w(1:ng4,1:klim)=0.5625*v(mg(g)%coarse_a,1:klim,g) + 0.1875*v(mg(g)%coarse_b,1:klim,g) &
                + 0.1875*v(mg(g)%coarse_c,1:klim,g) + 0.0625*v(mg(g)%coarse_d,1:klim,g)


  ! downscale grid
  do g=mg_maxlevel-1,2,-1

    ng=mg(g)%ifull

    ! extension
    v(1:ng,1:klim,g)=v(1:ng,1:klim,g)+w(1:ng,1:klim)
    
    call mgbounds(g,v(:,:,g),klim=klim)

    ! post smoothing
    do k=1,klim
      vnew(1:ng)=(mg(g)%zze*v(mg(g)%ie,k,g)+mg(g)%zzw*v(mg(g)%iw,k,g) &
                 +mg(g)%zzn*v(mg(g)%in,k,g)+mg(g)%zzs*v(mg(g)%is,k,g) &
                 -rhs(1:ng,k,g))/(helm(1:ng,k,g)-mg(g)%zz)
      v(1:ng,k,g)=vnew(1:ng)
    end do
    

    call mgbounds(g,v(:,:,g),klim=klim,corner=.true.)
    
    
    ! interpolation
    ! (there maybe issues interpolating near vertices, where we could
    !  possibly use a three point interpolation)
    ! JLM suggests using a three point (plane) fit instead of the 4
    ! point bi-linear fit.  This treats the vertices the same as the
    ! other points.
    ng4=mg(g)%ifull_coarse
    w(1:ng4,1:klim)=0.5625*v(mg(g)%coarse_a,1:klim,g) + 0.1875*v(mg(g)%coarse_b,1:klim,g) &
                  + 0.1875*v(mg(g)%coarse_c,1:klim,g) + 0.0625*v(mg(g)%coarse_d,1:klim,g)
  end do


  ! fine grid
  g=1
  ng=mg(g)%ifull
  
  ! convert back to local grid if required
  if (mg(g)%merge_len>1) then
    do n=1,npan
      ir=mod(mg(g)%merge_pos(n)-1,mg(g)%merge_row)+1   ! index for proc row
      ic=(mg(g)%merge_pos(n)-1)/mg(g)%merge_row+1      ! index for proc col
      do k=1,klim
        do jj=1,jpan
          iq_a=1+(jj-1)*ipan+(n-1)*ipan*jpan
          iq_b=jj*ipan+(n-1)*ipan*jpan
          iq_c=1+(ir-1)*ipan+(jj-1+(ic-1)*jpan)*jpan+(n-1)*ipan*jpan*mg(g)%merge_len
          iq_d=ir*ipan+(jj-1+(ic-1)*jpan)*jpan+(n-1)*ipan*jpan*mg(g)%merge_len
          vnew(iq_a:iq_b)=w(iq_c:iq_d,k)
        end do
        w(1:ifull,k)=vnew(1:ifull)
      end do
    end do
  end if
  
  ! extension
  iv(1:ifull,1:klim)=iv(1:ifull,1:klim)+w(1:ifull,1:klim)
  
  call bounds(iv,klim=klim)

  ! post smoothing
  do nc=1,maxcolour
    ifc=ifullx(nc)
    do k=1,klim
      dsol(iqx(1:ifc,nc),k)=( zznc(1:ifc,nc)*iv(iqn(1:ifc,nc),k) + zzwc(1:ifc,nc)*iv(iqw(1:ifc,nc),k)  &
                            + zzec(1:ifc,nc)*iv(iqe(1:ifc,nc),k) + zzsc(1:ifc,nc)*iv(iqs(1:ifc,nc),k)  &
                            +(zzc(1:ifc,nc)-helmc(1:ifc,nc,k))*iv(iqx(1:ifc,nc),k) - rhsc(1:ifc,nc,k)) &
                           /(helmc(1:ifc,nc,k)-zzc(1:ifc,nc))
      iv(iqx(1:ifc,nc),k) = iv(iqx(1:ifc,nc),k) + dsol(iqx(1:ifc,nc),k)
    end do
    call bounds(iv, klim=klim, colour=nc)
  end do

  ! test for convergence
  !if (itr>=itstest) then
  
    do k=1,klim
      dsolmax(k)=maxval(abs(dsol(1:ifull,k)))
    end do    ! convergence test
    call ccmpi_allreduce(dsolmax(1:klim),dsolmax_g(1:klim),"max",comm_world)
    knew=klim
    do k=klim,1,-1
      iters(k)=itr
      if (dsolmax_g(k)>=restol*(smax_g(k)-smin_g(k))) exit
      knew=k-1
    end do
    klim=knew
    if (klim<1) exit
  
  !  ! MJT - reduce collective MPI calls by anticipating convergence
  !  itsave1=itsave2
  !  itsave2=itr
  !  itserr1=itserr2
  !  itserr2=log10(dsolmax_g(1))
  !  gd=(itserr2-itserr1)/real(itsave2-itsave1)
  !  ci=itserr2-gd*real(itsave2)
  !  if (gd/=0.) then
  !    itstest=nint((log10(restol*(smax_g(1)-smin_g(1)))-ci)/gd)
  !    itstest=max(itstest,itr+1)
  !  else
  !    itstest=itr+1
  !  end if
  !  if (myid==0.and.nmaxpr==1) then
  !    write(6,*) "itr,itstest ",itr,itstest
  !  end if
  !  itc=itc+1
  !end if

end do

if (myid==0) then
  if (ktau<6.or.iters(1)>itr_max-5) then
    do k=1,kl
      write(6,*) "mg ktau,k,iter ",ktau,k,iters(k)
    end do
    !write(6,*) "itc ",itc
  end if
end if

call end_log(helm_end)

return
end subroutine mghelm

! This version is for the ocean and ice
subroutine mgmlo

use cc_mpi

implicit none

if (sorfirst) then
  write(6,*) "ERROR: mgsormlo requires mgsor_init to be called first"
  call ccmpi_abort(-1)
end if

! Use
! DIV s - H s - J s^2 = R
! s=g+e

! DIV g - H g - J g^2 = R - E
! DIV e - (H + D) e - J e^2 = E
! D = 2 J g


return
end subroutine mgmlo

! this routine allows multi-grid bounds updates
! This is based on cc_mpi bounds routines, but
! accomodates the g-th multi-grid
subroutine mgbounds(g,vdat,klim,corner,gmode)

use cc_mpi

implicit none

include 'mpif.h'
include 'newmpar.h'

integer, intent(in) :: g
integer, intent(in), optional :: klim, gmode
integer :: kx, lmode, send_len, recv_len, iq
integer :: iproc, iq_b, iq_e
integer(kind=4) :: ierr, itag=0, rproc, sproc, ltype
integer(kind=4), dimension(MPI_STATUS_SIZE,2*nproc) :: status
real, dimension(:,:), intent(inout) :: vdat
logical, intent(in), optional :: corner
logical extra

!call start_log(mgbounds_begin)

extra=.false.
if (present(corner)) extra=corner
kx=size(vdat,2)
if (present(klim)) kx=klim
lmode=0
if (present(gmode)) lmode=gmode

#ifdef r8i8
ltype=MPI_DOUBLE_PRECISION
#else
ltype=MPI_REAL
#endif   

!     Set up the buffers to send
nreq = 0
do iproc = 1,nproc-1  !
  rproc = modulo(myid-iproc,nproc)  ! Recv from
  if ( extra ) then
    recv_len = mg_bnds(g,rproc)%rlenx
  else
    recv_len = mg_bnds(g,rproc)%rlen
  end if
  if ( lmode == 1 ) then
    recv_len=recv_len*bnds(rproc)%mlomsk
  end if
  if ( recv_len > 0 ) then
    nreq = nreq + 1
    recv_len=recv_len*kx
    call MPI_IRecv( mg_bnds(g,rproc)%rbuf(1), recv_len, &
                    ltype, rproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
  end if
end do
do iproc = 1,nproc-1
  sproc = modulo(myid+iproc,nproc)  ! Send to
  if ( extra ) then
    send_len = mg_bnds(g,sproc)%slenx
  else
    send_len = mg_bnds(g,sproc)%slen
  end if
  if (lmode==1) then
    send_len=send_len*bnds(sproc)%mlomsk
  end if
  if ( send_len > 0 ) then
!cdir nodep
    do iq=1,send_len
      iq_b = 1+(iq-1)*kx
      iq_e = iq*kx
      mg_bnds(g,sproc)%sbuf(iq_b:iq_e) = vdat(mg_bnds(g,sproc)%send_list(iq),:)
    end do
    nreq = nreq + 1
    send_len=send_len*kx
    call MPI_ISend( mg_bnds(g,sproc)%sbuf(1), send_len, &
                    ltype, sproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
  end if
end do

! Finally see if there are any points on my own processor that need
! to be fixed up. This will only be in the case when nproc < npanels.
if ( extra ) then
  recv_len = mg_bnds(g,myid)%rlenx
else
  recv_len = mg_bnds(g,myid)%rlen
end if
!cdir nodep
do iq=1,recv_len
  ! request_list is same as send_list in this case
  vdat(mg(g)%ifull+mg_bnds(g,myid)%unpack_list(iq),:) = vdat(mg_bnds(g,myid)%request_list(iq),:)
end do

if ( nreq > 0 ) then
  call start_log(mpiwait_begin)
  call MPI_Waitall(nreq,ireq,status,ierr)
  call end_log(mpiwait_end)
end if

do iproc = 1,nproc-1
  rproc = modulo(myid-iproc,nproc)  ! Recv from
  if ( extra ) then
    recv_len = mg_bnds(g,rproc)%rlenx
  else
    recv_len = mg_bnds(g,rproc)%rlen
  end if
  if (lmode==1) then
    recv_len=recv_len*bnds(rproc)%mlomsk
  end if
  if ( recv_len > 0 ) then
!cdir nodep
    do iq=1,recv_len
      iq_b = 1+(iq-1)*kx
      iq_e = iq*kx
      vdat(mg(g)%ifull+mg_bnds(g,rproc)%unpack_list(iq),:) = mg_bnds(g,rproc)%rbuf(iq_b:iq_e)
    end do
  end if
end do

!call end_log(mgbounds_end)

return
end subroutine mgbounds

! This subroutine merges datasets when upscaling with the multi-grid solver
subroutine mgcollect(g,vdat,klim,gmode)

implicit none

include 'newmpar.h'

integer, intent(in) :: g
integer, intent(in), optional :: klim,gmode
integer nmax,npanx,kx,lmode
integer msg_len,lmsg
real, dimension(:,:), intent(inout) :: vdat

!call start_log(mggather_begin)

! merge length
nmax=mg(g)%merge_len
if (nmax<=1) return

kx=size(vdat,2)
if (present(klim)) kx=klim
lmode=0
if (present(gmode)) lmode=gmode

! The following trick allows a multi-grid global gather
! without needing a separate subroutine
npanx=npan
if (mg(g)%globgath) npanx=1

msg_len=mg(g)%ifull/(nmax*npanx) ! message unit size
lmsg =msg_len*kx

if (npanx==1) then ! usually face_decomp

  call mgcollect_face(g,vdat,kx,lmode,nmax,msg_len,lmsg)

else ! usually npanx==6 for uniform_decomp

  call mgcollect_uniform(g,vdat,kx,lmode,nmax,msg_len,lmsg)

end if

!call end_log(mggather_end)
  
return
end subroutine mgcollect

! this version of mgcollect uses MPI_allgather and is optmised
! for face decomposition
subroutine mgcollect_face(g,vdat,kx,lmode,nmax,msg_len,lmsg)

use cc_mpi

implicit none

include 'mpif.h'

integer, intent(in) :: g, kx, lmode, nmax, msg_len, lmsg
integer i,k,iq_a,iq_b,iq_c,iq_d
integer nrow,ncol,ilen_a,ilen_b
integer xproc,ir,ic,is,ie,js,je,jj
integer(kind=4) :: ierr, ltype, ilen, lcomm
real, dimension(:,:), intent(inout) :: vdat
real, dimension(lmsg) :: tdat
real, dimension(lmsg*nmax) :: tdat_g

#ifdef r8i8
ltype=MPI_DOUBLE_PRECISION
#else
ltype=MPI_REAL
#endif  

! prep data for sending around the merge
nrow=mg(g)%ipan/mg(g)%merge_row  ! number of points along a row per processor
ncol=msg_len/nrow              ! number of points along a col per processor

! Use MPI to optimise this gather

do k=1,kx
  iq_a=1+(k-1)*msg_len
  iq_b=k*msg_len
  tdat(iq_a:iq_b)=vdat(1:msg_len,k)
end do

ilen=lmsg
lcomm=mg(g)%comm
if (lmode==1) then
  lcomm=mg(g)%comm_mlo
end if
call MPI_AllGather(tdat,ilen,ltype,tdat_g,ilen,ltype,lcomm,ierr)

! add blanks for missing processors
if (lmode==1) then
  ilen_a=0
  ilen_b=0
  tdat=tdat_g
  do i=1,nmax
    if (bnds(mg(g)%merge_list(i,1))%mlomsk==1) then
      tdat_g(ilen_a+1:ilen_a+lmsg)=tdat(ilen_b+1:ilen_b+lmsg)
      ilen_b=ilen_b+lmsg
    else
      tdat_g(ilen_a+1:ilen_a+lmsg)=0.
    end if
    ilen_a=ilen_a+lmsg
  end do
end if

do xproc=1,nmax
  ir=mod(xproc-1,mg(g)%merge_row)+1   ! index for proc row
  ic=(xproc-1)/mg(g)%merge_row+1      ! index for proc col
  is=(ir-1)*nrow+1
  ie=ir*nrow
  js=(ic-1)*ncol+1
  je=ic*ncol
  do k=1,kx
    do jj=js,je
      iq_a=is+(jj-1)*mg(g)%ipan
      iq_b=iq_a+nrow-1
      iq_c=1+(jj-js)*nrow+(k-1)*msg_len+(xproc-1)*lmsg
      iq_d=iq_c+nrow-1
      vdat(iq_a:iq_b,k)=tdat_g(iq_c:iq_d)
    end do
  end do
end do
  
return
end subroutine mgcollect_face

! this version of mgcollect uses point-to-point communications
! and is optimised for uniform decomposition
subroutine mgcollect_uniform(g,vdat,kx,lmode,nmax,msg_len,lmsg)

use cc_mpi

implicit none

include 'mpif.h'
include 'newmpar.h'

integer, intent(in) :: g, kx, lmode, nmax, msg_len, lmsg
integer n,k,np,nm,nx,i,nrow,ncol
integer xproc,yproc,msg_off
integer ir,ic,ida,iq_a,iq_b,iq_c,iq_d
integer is,js,je,jj
integer(kind=4) :: rproc, sproc
integer(kind=4) :: sreq, rreq, ierr, itag=0, ltype, ilen, lcomm
integer(kind=4), dimension(MPI_STATUS_SIZE,2*nproc) :: status
integer, dimension(-npan:npan*(nmax-1)) :: rarry, sarry, roff, soff
integer, dimension(0:nproc-1) :: rlist, slist, pr, ps
integer, dimension(nproc) :: rp, sp
real, dimension(:,:), intent(inout) :: vdat
real, dimension(lmsg*npan,0:(nmax-1)*npan) :: rrtn
real, dimension(lmsg*npan,(nmax-1)*npan) :: sdep
logical ftest

#ifdef r8i8
ltype=MPI_DOUBLE_PRECISION
#else
ltype=MPI_REAL
#endif  

! prep data for sending around the merge
nrow=mg(g)%ipan/mg(g)%merge_row  ! number of points along a row per processor
ncol=msg_len/nrow              ! number of points along a col per processor

! Send them all and let MPI sort them out...
! This approach allows us to send data for different panels without using multiple MPI_AllGathers
! for each panel

! intialise arrays
rreq =0
sreq =0
rlist=0
slist=0
rp   =-1
sp   =-1
pr   =0
ps   =0
rarry=0
sarry=0
roff =0
soff =0

! loop over panels
do n=1,npan
  nx=mg(g)%merge_pos(n)  ! the location of this processor in the merge
  msg_off=(n-1)*msg_len ! offset for input array
  
  ! loop over merge members
  do xproc=1,nmax-1   
    ida=n+npan*(xproc-1) ! index for packing arrays
    
    ! Recv processor
    np=modulo(nx+xproc,nmax)
    if (np==0) np=nmax 
    rproc=mg(g)%merge_list(np,n)
    ftest=lmode==0.or.bnds(rproc)%mlomsk==1
    if (ftest) then
      if (rlist(rproc)==0) then
        rreq=rreq+1     ! number of MPI Recv requests
        rp(rreq)=rproc  ! Processor associated with this request
        pr(rproc)=rreq  ! Request associated with this processor
      end if
      rlist(rproc)=rlist(rproc)+1 ! Number of sub-messages Recv from this processor
      rarry(ida)=pr(rproc)        ! Request number as a function of packing index
      roff(ida)=rlist(rproc)      ! Offset number as a function of packing index
    end if
    
    ! Send processor
    nm=modulo(nx-xproc,nmax) 
    if (nm==0) nm=nmax 
    sproc=mg(g)%merge_list(nm,n)
    ftest=lmode==0.or.bnds(sproc)%mlomsk==1
    if (ftest) then
      if (slist(sproc)==0) then
        sreq=sreq+1     ! number of MPI Send requests
        sp(sreq)=sproc  ! Processor associated with this request
        ps(sproc)=sreq  ! Request associated with this processor
      end if
      slist(sproc)=slist(sproc)+1 ! Number of sub-messages Send from this processor
      sarry(ida)=ps(sproc)        ! Request number as a function of packing index
      soff(ida)=slist(sproc)      ! Offset number as a function of packing index
     
      ! Pack data into send arrays
      do k=1,kx
        sdep(1+(k-1)*msg_len+(soff(ida)-1)*lmsg:k*msg_len+(soff(ida)-1)*lmsg,sarry(ida))=vdat(1+msg_off:msg_len+msg_off,k)
      end do
    end if
    
  end do

  ! data for myid  
  !xproc=0
  ida=n-npan                ! index for packing arrays
  rlist(myid)=rlist(myid)+1 ! Number of sub-messages Recv from this processor
  rarry(ida)=0              ! Request number as a function of packing index
  roff(ida)=rlist(myid)     ! Offset number as a function of packing index
  do k=1,kx
    rrtn(1+(k-1)*msg_len+(roff(ida)-1)*lmsg:k*msg_len+(roff(ida)-1)*lmsg,rarry(ida))=vdat(1+msg_off:msg_len+msg_off,k)
  end do

end do

! MPI Recv
nreq=0
do i=1,rreq
  rproc=rp(i)
  ilen=lmsg*rlist(rproc)
  nreq=nreq+1
  call MPI_IRecv( rrtn(:,i), ilen, ltype, rproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
end do
  
! MPI Send
do i=1,sreq
  sproc=sp(i)
  ilen=lmsg*slist(sproc)
  nreq=nreq+1
  call MPI_ISend( sdep(:,i), ilen, ltype, sproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
end do

if (nreq>0) then
  call start_log(mpiwait_begin)
  call MPI_Waitall(nreq,ireq,status,ierr)
  call end_log(mpiwait_end)
end if
  
! unpack buffers
vdat=0.
do n=1,npan
  nx=mg(g)%merge_pos(n)  ! the location of this processor in the merge
  do yproc=1,nmax
    xproc=modulo(yproc-nx,nmax)         ! processor data to unpack
    ida=n+npan*(xproc-1)                ! index for packing arrays
    ir=mod(yproc-1,mg(g)%merge_row)+1   ! index for proc row
    ic=(yproc-1)/mg(g)%merge_row+1      ! index for proc col
    
    is=(ir-1)*nrow+1
    js=(ic-1)*ncol+1
    je=ic*ncol
    do k=1,kx
      do jj=js,je
        iq_a=is+(jj-1)*mg(g)%ipan+(n-1)*msg_len*nmax
        iq_b=iq_a+nrow-1
        iq_c=1+(jj-js)*nrow+(k-1)*msg_len+(roff(ida)-1)*lmsg
        iq_d=iq_c+nrow-1
        vdat(iq_a:iq_b,k)=rrtn(iq_c:iq_d,rarry(ida))
      end do
    end do
  end do
end do
  
return
end subroutine mgcollect_uniform

! Initialise multi-grid arrays
subroutine mgsor_init

use cc_mpi

implicit none

include 'newmpar.h'

integer g,np,iq,iqq,iproc,xlen,nn,ii,jj
integer mipan,mjpan,hipan,hjpan,mil_g,ia,ja
integer i,j,n,jextra,mg_npan,mxpr,mypr
integer cid,ix,jx,colour,rank,ncol,nrow
integer drow,dcol,npanx,na,nx,ny,sii,eii,sjj,ejj
logical lglob

if (.not.sorfirst) return

! Begin full initialisation
mg_maxsize=ifull+iextra ! first guess
lglob=(nproc==1)        ! Global gather flag
mipan=ipan              ! local number of rows
mjpan=jpan              ! local number of columns
mil_g=il_g              ! global grid size
mxpr=il_g/ipan          ! number of processors over rows
mypr=il_g/jpan          ! number of processors over columns


! calculate number of levels
! Let coarsest grid have 4 (2x2) points on each face. If grid size
! is not a power of 2, just make it as small as possible.
mg_maxlevel=1
g=1
do while (mod(il_g,2**g)==0)
  g=g+1
end do
mg_maxlevel=g

if (myid==0) then
  write(6,*) "Initialising multi-grid arrays"
  write(6,*) "il_g,mg_maxlevel ",il_g,mg_maxlevel
end if

allocate(mg(mg_maxlevel))
allocate(mg_bnds(mg_maxlevel,0:nproc-1))
do iproc=0,nproc-1
  do g=1,mg_maxlevel
    mg_bnds(g,iproc)%len=0
  end do
end do


! calculate fine grid for finest grid level
g=1
mg(1)%globgath=.false.
mg_npan=npan
mg(1)%merge_len=1
mg(1)%merge_row=1

allocate(mg(1)%fproc(mil_g,mil_g,0:npanels))
mg(1)%fproc(:,:,:)=fproc(:,:,:)

! check if coarse grid needs mgcollect
if (mod(mipan,2)/=0.or.mod(mjpan,2)/=0) then
  if (mod(mxpr,2)==0.and.mod(mypr,2)==0) then
    if (myid==0) then
      write(6,*) "Multi-grid gather4 at level ",g,mipan,mjpan
    end if
    
   ! This case occurs when there are multiple processors on a panel.
   ! Consequently, npan should be 1 or 6.
   if (npan>1.and.npan<6) then
     write(6,*) "ERROR: Invalid gather4"
     call ccmpi_abort(-1)
   end if
    
    mg(1)%merge_len=4
    mg(1)%merge_row=2
    mxpr=mxpr/2
    mypr=mypr/2
    hipan=mipan
    hjpan=mjpan
    mipan=mipan*2
    mjpan=mjpan*2

    allocate(mg(1)%merge_list(4,npan),mg(1)%merge_pos(npan))

    ! find my processor and surrounding members of the gather
    do n=1,npan
      nn=n-noff

      ix=-1
      jx=-1
      do jj=1,mil_g,hjpan
        do ii=1,mil_g,hipan
          if (mg(1)%fproc(ii,jj,nn)==myid) then
            ix=ii
            jx=jj
            exit
          end if
        end do
        if (ix>0) exit
      end do
      if (ix<1) then
        write(6,*) "ERROR: Cannot locate processor in gather4"
        stop
      end if
      ii=ix
      jj=jx
      if (mod((ii-1)/hipan,2)/=0) ii=ii-hipan
      if (mod((jj-1)/hjpan,2)/=0) jj=jj-hjpan
      ix=ix-ii ! offset for myid from ii
      jx=jx-jj ! offset for myid from jj

       
      mg(1)%merge_list(1,n)=mg(1)%fproc(ii,      jj      ,nn)
      mg(1)%merge_list(2,n)=mg(1)%fproc(ii+hipan,jj      ,nn)
      mg(1)%merge_list(3,n)=mg(1)%fproc(ii,      jj+hjpan,nn)
      mg(1)%merge_list(4,n)=mg(1)%fproc(ii+hipan,jj+hjpan,nn)
       
      do j=1,mil_g,mjpan
        do i=1,mil_g,mipan
          cid=mg(1)%fproc(i+ix,j+jx,nn) ! processor in same merge position as myid
                                        ! we will maintain communications with this processor
          do ja=1,mjpan
            do ia=1,mipan
              ! update fproc map with processor that owns this data
              mg(1)%fproc(i+ia-1,j+ja-1,nn)=cid
            end do
          end do
        end do
      end do
        
      do j=1,4
        if (mg(1)%merge_list(j,n)==myid) then
          mg(1)%merge_pos(n)=j
          exit
        end if
      end do
      
    end do

    ! fix any remaining panels
    if (npan==1) then
      do nn=1,npanels
        do j=1,mil_g,mjpan
          do i=1,mil_g,mipan
            cid=mg(g)%fproc(i+ix,j+jx,nn) ! processor in same merge position as myid
                                          ! we will maintain communications with this processor
            do ja=1,mjpan
              do ia=1,mipan
                ! update fproc map with processor that owns this data
                mg(g)%fproc(i+ia-1,j+ja-1,nn)=cid
              end do
            end do
          end do
        end do
      end do

      ! define local comm
      colour=mg(1)%merge_list(1,1)
      rank=mg(1)%merge_pos(1)-1
      call ccmpi_commsplit(mg(1)%comm,comm_world,colour,rank)
    end if
      
  else
    write(6,*) "ERROR: Grid m=1 requires gatherall for multi-grid solver"
    call ccmpi_abort(-1)
  end if
end if

mg(1)%ifull=mipan*mjpan*mg_npan
mg(1)%ifull_fine=mg(1)%ifull/4
jextra=2*(mipan+mjpan+2)*mg_npan
if (lglob) then
  jextra=0
end if
mg_maxsize=max(mg_maxsize,mg(1)%ifull+jextra)

allocate(mg(1)%fine(mg(1)%ifull_fine))
iqq=0
do n=1,mg_npan
  do jj=1,mjpan,2
    do ii=1,mipan,2
      iq=indx(ii,jj,n-1,mipan,mjpan)
      iqq=iqq+1
      mg(1)%fine(iqq)=iq
    end do
  end do
end do

np=mg(1)%ifull
allocate(mg(1)%in(np),mg(1)%ie(np),mg(1)%is(np),mg(1)%iw(np))
allocate(mg(1)%ine(np),mg(1)%inw(np),mg(1)%ise(np),mg(1)%isw(np))

call mg_index(1,mil_g,mipan,mjpan)


! loop over upscaled grids
do g=2,mg_maxlevel

  mipan=mipan/2
  mjpan=mjpan/2
  mil_g=mil_g/2

  ! assign processors to each grid point
  ! Below are default values for fproc which are modified
  ! when we collect data over processors
  allocate(mg(g)%fproc(mil_g,mil_g,0:npanels))

  ! Calculate size of grid at this level
  do nn=0,npanels
    do jj=1,2*mil_g,2
      ja=(jj-1)/2+1
      do ii=1,2*mil_g,2
        ia=(ii-1)/2+1
        mg(g)%fproc(ia,ja,nn)=mg(g-1)%fproc(ii,jj,nn)
      end do
    end do
  end do


  ! first assume no gather for upscaled grid
  mg(g)%merge_len=1
  mg(g)%merge_row=1
  mg(g)%globgath=.false.
  
  ! check for multi-grid gather
  if (mod(mipan,2)/=0.or.mod(mjpan,2)/=0) then ! grid cannot be subdivided on current processor
  
    if (mod(mxpr,2)==0.and.mod(mypr,2)==0) then ! collect data over adjacent processors
      if (myid==0) then
        write(6,*) "Multi-grid gather4 at level           ",g,mipan,mjpan
      end if

      ! This case occurs when there are multiple processors on a panel.
      ! Consequently, npan should be 1 or 6.
      if (npan>1.and.npan<6) then
        write(6,*) "ERROR: Invalid gather4"
        call ccmpi_abort(-1)
      end if

      mg(g)%merge_len=4
      mg(g)%merge_row=2
      mxpr=mxpr/2
      mypr=mypr/2
      hipan=mipan ! size of grid before mg collect
      hjpan=mjpan
      mipan=mipan*2
      mjpan=mjpan*2

      allocate(mg(g)%merge_list(4,npan),mg(g)%merge_pos(npan))
      
      do n=1,npan
        nn=n-noff

        ! find my processor and surrounding members of the gather
        ix=-1
        jx=-1
        do jj=1,mil_g,hjpan
          do ii=1,mil_g,hipan
            if (mg(g)%fproc(ii,jj,nn)==myid) then
              ix=ii
              jx=jj
              exit
            end if
          end do
          if (ix>0) exit
        end do
        if (ix<1) then
          write(6,*) "ERROR: Cannot locate processor in gather4"
          call ccmpi_abort(-1)
        end if
        ii=ix
        jj=jx
        if (mod((ii-1)/hipan,2)/=0) ii=ii-hipan ! left corner of merge
        if (mod((jj-1)/hjpan,2)/=0) jj=jj-hjpan ! bottom corner of merge
        ix=ix-ii ! offset for myid from ii
        jx=jx-jj ! offset for myid from jj
       
        mg(g)%merge_list(1,n)=mg(g)%fproc(ii,      jj      ,nn)
        mg(g)%merge_list(2,n)=mg(g)%fproc(ii+hipan,jj      ,nn)
        mg(g)%merge_list(3,n)=mg(g)%fproc(ii,      jj+hjpan,nn)
        mg(g)%merge_list(4,n)=mg(g)%fproc(ii+hipan,jj+hjpan,nn)

        do j=1,mil_g,mjpan
          do i=1,mil_g,mipan
            cid=mg(g)%fproc(i+ix,j+jx,nn) ! processor in same merge position as myid
                                          ! we will maintain communications with this processor
            do ja=1,mjpan
              do ia=1,mipan
                ! update fproc map with processor that owns this data
                mg(g)%fproc(i+ia-1,j+ja-1,nn)=cid
              end do
            end do
          end do
        end do

        do j=1,4
          if (mg(g)%merge_list(j,n)==myid) then
            mg(g)%merge_pos(n)=j
            exit
          end if
        end do

      end do
      
      ! fix any remaining panels
      if (npan==1) then
        do nn=0,npanels
          do j=1,mil_g,mjpan
            do i=1,mil_g,mipan
              cid=mg(g)%fproc(i+ix,j+jx,nn) ! processor in same merge position as myid
                                            ! we will maintain communications with this processor
              do ja=1,mjpan
                do ia=1,mipan
                  ! update fproc map with processor that owns this data
                  mg(g)%fproc(i+ia-1,j+ja-1,nn)=cid
                end do
              end do
            end do
          end do
        end do
      end if
    
    else if (.not.lglob) then ! collect all data to one processor
      if (myid==0) then
        write(6,*) "Multi-grid gatherall at level         ",g,mipan,mjpan
      end if
      lglob=.true.
#ifdef uniform_decomp
      mg(g)%merge_len=mxpr*mypr
#else
      mg(g)%merge_len=min(6*mxpr*mypr,nproc)
      mg(g)%globgath=.true.
      mg_npan=6
#endif

      hipan=mipan ! size of grid before collect
      hjpan=mjpan

      mg(g)%merge_row=mxpr
      mipan=mipan*mxpr
      mjpan=mjpan*mypr
      mxpr=1
      mypr=1
      
      ! find gather members
#ifdef uniform_decomp
      allocate(mg(g)%merge_list(mg(g)%merge_len,npan),mg(g)%merge_pos(npan))
      do n=1,npan
        nn=n-noff
        iqq=0
        do jj=1,mil_g,hjpan
          do ii=1,mil_g,hipan
            iqq=iqq+1
            mg(g)%merge_list(iqq,n)=mg(g)%fproc(ii,jj,nn)
          end do
        end do
        if (iqq/=mg(g)%merge_len) then
          write(6,*) "ERROR: merge_len mismatch ",iqq,mg(g)%merge_len,g
          stop
        end if
        do i=1,mg(g)%merge_len
          if (mg(g)%merge_list(i,n)==myid) then
            mg(g)%merge_pos(n)=i
            exit
          end if
        end do
      end do
#else
      allocate(mg(g)%merge_list(mg(g)%merge_len,1),mg(g)%merge_pos(1))      
      iqq=0
      do n=1,6/npan
        nn=(n-1)*npan
        do jj=1,mil_g,hjpan
          do ii=1,mil_g,hipan
            iqq=iqq+1
            mg(g)%merge_list(iqq,1)=mg(g)%fproc(ii,jj,nn)
          end do
        end do
      end do
      if (iqq/=mg(g)%merge_len) then
        write(6,*) "ERROR: merge_len mismatch ",iqq,mg(g)%merge_len,g
        stop
      end if
      do i=1,mg(g)%merge_len
        if (mg(g)%merge_list(i,1)==myid) then
          mg(g)%merge_pos(1)=i
          exit
        end if
      end do
#endif

      ! modify fproc for remaining processor
      mg(g)%fproc(:,:,:)=myid

    else ! all data is already on one processor
      if (g/=mg_maxlevel) then
        write(6,*) "ERROR: g.ne.mg_maxlevel ",g,mg_maxlevel
        stop
      end if
      if (myid==0) then
        write(6,*) "Multi-grid toplevel                   ",g,mipan,mjpan
      end if
      allocate(mg(g)%merge_pos(1))
      mg(g)%merge_pos(1)=1
    end if

    ! define local comm
    if (mg(g)%merge_len>1) then
      if (npan==1.or.mg(g)%globgath) then
        colour=mg(g)%merge_list(1,1)
        rank=mg(g)%merge_pos(1)-1
        call ccmpi_commsplit(mg(g)%comm,comm_world,colour,rank)
      end if
    end if
  
  else
    if (myid==0) then
      write(6,*) "Multi-grid local subdivision at level ",g,mipan,mjpan
    end if
    ! no messages sent, but we allocate this arrays for the
    ! coarse calculation below
    allocate(mg(g)%merge_pos(1))
    mg(g)%merge_pos(1)=1
  end if
  

  ! total number of points over all panels on this processor
  mg(g)%ipan=mipan
  np=mipan*mjpan*mg_npan
  jextra=2*(mipan+mjpan+2)*mg_npan
  nrow=mg(g)%merge_row
  ncol=mg(g)%merge_len/nrow
  drow=mipan/nrow
  dcol=mjpan/ncol
  npanx=mg_npan

  ! this trick handles a global interpolation
  if (lglob) then
    jextra=0
    npanx=1
    dcol=6*mjpan/ncol
  end if

  ! ifine is the index of the SW point on the next finer grid.
  ! Because of the staggering this means (1,1) maps to (1,1).
  ! Using each grid's own indices, the neighbour of (i,j) is (2i-1,2j-1).

  mg(g)%ifull=np
  mg(g)%ifull_fine=np/4
  mg_maxsize=max(mg_maxsize,np+jextra)
    
  if (g<mg_maxlevel) then
    np=mg(g)%ifull_fine
    allocate(mg(g)%fine(np))
    ! mipan and mjpan should always be an even number here
    iqq=0
    do n=1,mg_npan
      do jj=1,mjpan,2
        do ii=1,mipan,2
          iq=indx(ii,jj,n-1,mipan,mjpan)
          iqq=iqq+1
          mg(g)%fine(iqq)=iq
        end do
      end do
    end do
  end if

  np=mg(g)%ifull
  allocate(mg(g)%in(np),mg(g)%ie(np),mg(g)%is(np),mg(g)%iw(np))
  allocate(mg(g)%ine(np),mg(g)%inw(np),mg(g)%ise(np),mg(g)%isw(np))

  call mg_index(g,mil_g,mipan,mjpan)

  ! Set up pointers for neighbours on the next coarser grid.
  ! ic1 is the nearest neigbour, ic2 and ic3 are equally distant and
  ! ic4 is the most distant.
  ! These are best defined separately on each face. The code here
  ! requires an even number of points to ensure the coarser grid
  ! always sits inside the fine grid.
  ! For odd i (i = 2*ic - 1), the nearest coarse neighbour is to the
  ! E at ic.
  ! For even i (i=2*ic), the nearest coarse neighbour is to the W
  ! at ic.
  ! For odd j, neighbor is to N, even j it is to S.

  mg(g)%ifull_coarse=mg(g-1)%ifull
  np=mg(g)%ifull_coarse
  allocate(mg(g)%coarse_a(np),mg(g)%coarse_b(np),mg(g)%coarse_c(np),mg(g)%coarse_d(np))

  do n=1,npanx
  
    nn=min(n,mg(g)%merge_len)
    na=mg(g)%merge_pos(nn)
    nx=mod(na-1,nrow)+1
    ny=(na-1)/nrow+1
    sii=(nx-1)*drow+1
    eii=nx*drow
    sjj=(ny-1)*dcol+1
    ejj=ny*dcol
    
    do jj=sjj,ejj
      ja=2*(jj-sjj)+1
      do ii=sii,eii
        ia=2*(ii-sii)+1

        iqq=indx(ii,jj,n-1,mipan,mjpan)   ! coarse grid

        ! odd, odd          
        iq =indx(ia,ja,n-1,2*drow,2*dcol)     ! fine grid
        mg(g)%coarse_a(iq)=          iqq
        mg(g)%coarse_b(iq)= mg(g)%is(iqq)
        mg(g)%coarse_c(iq)= mg(g)%iw(iqq)
        mg(g)%coarse_d(iq)=mg(g)%isw(iqq)
          
        ! odd, even
        iq =indx(ia,ja+1,n-1,2*drow,2*dcol)   ! fine grid
        mg(g)%coarse_a(iq)=          iqq
        mg(g)%coarse_b(iq)= mg(g)%in(iqq)
        mg(g)%coarse_c(iq)= mg(g)%iw(iqq)
        mg(g)%coarse_d(iq)=mg(g)%inw(iqq)
          
        ! even, odd
        iq =indx(ia+1,ja,n-1,2*drow,2*dcol)   ! fine grid 
        mg(g)%coarse_a(iq)=          iqq
        mg(g)%coarse_b(iq)= mg(g)%is(iqq)
        mg(g)%coarse_c(iq)= mg(g)%ie(iqq)
        mg(g)%coarse_d(iq)=mg(g)%ise(iqq)
          
        ! even, even
        iq =indx(ia+1,ja+1,n-1,2*drow,2*dcol) ! fine grid
        mg(g)%coarse_a(iq)=          iqq
        mg(g)%coarse_b(iq)= mg(g)%ie(iqq)
        mg(g)%coarse_c(iq)= mg(g)%in(iqq)
        mg(g)%coarse_d(iq)=mg(g)%ine(iqq)
      
      end do
    end do
  end do

end do

! free some memory
do g=1,mg_maxlevel
  deallocate(mg(g)%fproc)
end do

sorfirst=.false.
if (myid==0) then
  write(6,*) "Finished initialising multi-grid arrays"
end if

return
end subroutine mgsor_init

subroutine mgzz_init(zz,zzn,zze,zzw,zzs)

use cc_mpi

implicit none

include 'newmpar.h'

integer g,np,iqq,iq
real, dimension(ifull), intent(in) :: zz,zzn,zze,zzw,zzs
real, dimension(mg_maxsize,5) :: dum
real, parameter :: dfac=0.25 ! adjustment for grid spacing

  ! Later equations use zz/delta**2. It's more efficient to put
  ! this scaling in here. delta = delta_0 * 2**(m-1) so increases
  ! by a factor of 2 each grid. Therefore using a factor of 1/4
  ! in the averaging calculation will compensate.

if (myid==0) then
  write(6,*) "Initialising atmosphere multi-grid coupling arrays"
end if

if (zzfirst) then
  do g=1,mg_maxlevel
    np=mg(g)%ifull
    allocate(mg(g)%zzn(np),mg(g)%zze(np),mg(g)%zzs(np),mg(g)%zzw(np),mg(g)%zz(np))
  end do
  zzfirst=.false.
end if

! no need for bounds call as all indices lie within this processor
g=1
np=mg(1)%ifull

dum(1:ifull,1)=zzn
dum(1:ifull,2)=zzs
dum(1:ifull,3)=zze
dum(1:ifull,4)=zzw
dum(1:ifull,5)=zz
if (mg(1)%merge_len>1) then
  call mgcollect(1,dum)
end if
mg(1)%zzn(:)=dum(1:np,1)
mg(1)%zzs(:)=dum(1:np,2)
mg(1)%zze(:)=dum(1:np,3)
mg(1)%zzw(:)=dum(1:np,4)
mg(1)%zz(:)=dum(1:np,5)

do g=2,mg_maxlevel
  np=mg(g)%ifull

  do iq=1,mg(g-1)%ifull_fine
    iqq=mg(g-1)%fine(iq)
    dum(iq,1)=0.25*dfac*(mg(g-1)%zzn(iqq)+mg(g-1)%zzn(mg(g-1)%in(iqq))+mg(g-1)%zzn(mg(g-1)%ine(iqq))+mg(g-1)%zzn(mg(g-1)%ie(iqq)))
    dum(iq,2)=0.25*dfac*(mg(g-1)%zzs(iqq)+mg(g-1)%zzs(mg(g-1)%in(iqq))+mg(g-1)%zzs(mg(g-1)%ine(iqq))+mg(g-1)%zzs(mg(g-1)%ie(iqq)))
    dum(iq,3)=0.25*dfac*(mg(g-1)%zze(iqq)+mg(g-1)%zze(mg(g-1)%in(iqq))+mg(g-1)%zze(mg(g-1)%ine(iqq))+mg(g-1)%zze(mg(g-1)%ie(iqq)))
    dum(iq,4)=0.25*dfac*(mg(g-1)%zzw(iqq)+mg(g-1)%zzw(mg(g-1)%in(iqq))+mg(g-1)%zzw(mg(g-1)%ine(iqq))+mg(g-1)%zzw(mg(g-1)%ie(iqq)))
    dum(iq,5)=0.25*dfac*(mg(g-1)%zz(iqq) +mg(g-1)%zz(mg(g-1)%in(iqq)) +mg(g-1)%zz(mg(g-1)%ine(iqq)) +mg(g-1)%zz(mg(g-1)%ie(iqq)))
  end do
  if (mg(g)%merge_len>1) then
    call mgcollect(g,dum)
  end if
  mg(g)%zzn=dum(1:np,1)
  mg(g)%zzs=dum(1:np,2)
  mg(g)%zze=dum(1:np,3)
  mg(g)%zzw=dum(1:np,4)
  mg(g)%zz=dum(1:np,5)
end do

if (myid==0) then
  write(6,*) "Finished initialising atmosphere multi-grid coupling arrays"
end if

return
end subroutine mgzz_init

! Set up the indices required for the multigrid scheme.
subroutine mg_index(g,mil_g,mipan,mjpan)

use cc_mpi

implicit none

include 'mpif.h'
include 'newmpar.h'

integer, intent(in) :: g, mil_g, mipan, mjpan
integer, dimension(6*mil_g*mil_g) :: mg_qproc, mg_colourmask
integer, dimension(6*mil_g*mil_g) :: jn_g, je_g, js_g, jw_g, jne_g, jse_g, jsw_g, jnw_g
integer, parameter, dimension(0:5) :: npann=(/ 1, 103, 3, 105, 5, 101 /)
integer, parameter, dimension(0:5) :: npane=(/ 102, 2, 104, 4, 100, 0 /)
integer, parameter, dimension(0:5) :: npanw=(/ 5, 105, 1, 101, 3, 103 /)
integer, parameter, dimension(0:5) :: npans=(/ 104, 0, 100, 2, 102, 4 /)
integer(kind=4), dimension(MPI_STATUS_SIZE,2*nproc) :: status
integer, dimension(npanels+1) :: mioff, mjoff
integer, dimension(2*(mipan+mjpan+2)*(npanels+1)) :: dum
integer i, j, n, iq, iqq, iqg, ii, mfull, mfull_g, ncount
integer mg_colour_np, iloc, jloc, nloc, jextra
integer iext, iql, iproc, xlen, jx, nc, xlev
integer(kind=4) :: itag=0, rproc, sproc, ierr
logical lflag, lglob

! size of this grid
mfull_g = 6*mil_g*mil_g


! calculate processor map in iq coordinates
ncount=0
lglob=.true.
do n=0,npanels
  lflag=.true.
  do j=1,mil_g
    do i=1,mil_g
      iq=indx(i,j,n,mil_g,mil_g)
      mg_qproc(iq)=mg(g)%fproc(i,j,n)
      if (mg_qproc(iq)/=myid) lglob=.false.
      ! ncount>=npan usually indicates that a global gather has occured
      if (lflag.and.mg_qproc(iq)==myid) then
        ncount=ncount+1
        mioff(ncount)=i-1
        mjoff(ncount)=j-1
        lflag=.false.
      end if
    end do
  end do
end do

if (ncount==0) then
  write(6,*) "ERROR: Cannot find myid in mg_proc"
  write(6,*) "myid,g ",myid,g
  write(6,*) "mg_proc ",maxval(mg_qproc),minval(mg_qproc),count(mg_qproc==myid)
  call ccmpi_abort(-1)
end if


! calculate global indices
do iq = 1, mfull_g
  jn_g(iq) = iq + mil_g
  js_g(iq) = iq - mil_g
  je_g(iq) = iq + 1
  jw_g(iq) = iq - 1
end do

do n = 0, npanels
  if (npann(n) < 100) then
    do ii = 1, mil_g
      jn_g(indx(ii,mil_g,n,mil_g,mil_g)) = indx(ii,1,npann(n),mil_g,mil_g)
    end do
  else
    do ii = 1, mil_g
      jn_g(indx(ii,mil_g,n,mil_g,mil_g)) = indx(1,mil_g+1-ii,npann(n)-100,mil_g,mil_g)
    end do
  endif
  if (npane(n) < 100) then
    do ii = 1, mil_g
      je_g(indx(mil_g,ii,n,mil_g,mil_g)) = indx(1,ii,npane(n),mil_g,mil_g)
    end do
  else
    do ii = 1, mil_g
      je_g(indx(mil_g,ii,n,mil_g,mil_g)) = indx(mil_g+1-ii,1,npane(n)-100,mil_g,mil_g)
    end do
  endif
  if (npanw(n) < 100) then
    do ii = 1, mil_g
      jw_g(indx(1,ii,n,mil_g,mil_g)) = indx(mil_g,ii,npanw(n),mil_g,mil_g)
    end do
  else
    do ii = 1, mil_g
      jw_g(indx(1,ii,n,mil_g,mil_g)) = indx(mil_g+1-ii,mil_g,npanw(n)-100,mil_g,mil_g)
    end do
  endif
  if (npans(n) < 100) then
    do ii = 1, mil_g
      js_g(indx(ii,1,n,mil_g,mil_g)) = indx(ii,mil_g,npans(n),mil_g,mil_g)
    end do
  else
    do ii = 1, mil_g
      js_g(indx(ii,1,n,mil_g,mil_g)) = indx(mil_g,mil_g+1-ii,npans(n)-100,mil_g,mil_g)
    end do
  endif
end do ! n loop

jnw_g = jn_g(jw_g)
jne_g = jn_g(je_g)
jse_g = js_g(je_g)
jsw_g = js_g(jw_g)

do n = 0, npanels
  ! Following treats unusual panel boundaries
  if (npanw(n) >= 100) then
    do j = 1, mil_g
      iq = indx(1,j,n,mil_g,mil_g)
      jnw_g(iq) = jw_g(jw_g(iq))
      jsw_g(iq) = je_g(jw_g(iq))
    end do
  endif
  if (npane(n) >= 100) then
    do j = 1, mil_g
      iq = indx(mil_g,j,n,mil_g,mil_g)
      jne_g(iq) = jw_g(je_g(iq))
      jse_g(iq) = je_g(je_g(iq))
    end do
  endif
end do

! Calculate local indices on this processor
if (lglob) then
  mg(g)%in=jn_g
  mg(g)%is=js_g
  mg(g)%ie=je_g
  mg(g)%iw=jw_g
  mg(g)%ine=jne_g
  mg(g)%inw=jnw_g
  mg(g)%ise=jse_g
  mg(g)%isw=jsw_g
else
  ! This only occurs with grids prior to globgath.  So npan and noff are still valid.
  do n=1,npan
    do j=1,mjpan
      do i=1,mipan
        iq = indx(i,j,n-1,mipan,mjpan) ! Local
        iqg = indx(i+mioff(n),j+mjoff(n),n-noff,mil_g,mil_g) ! Global

        iqq = jn_g(iqg)    ! Global neighbour index
        rproc = mg_qproc(iqq) ! Processor that has this point
        if ( rproc == myid ) then ! Just copy the value
          ! Convert global iqq to local value
          call indv_mpix(iqq,iloc,jloc,nloc,mil_g,mioff,mjoff,noff)
          mg(g)%in(iq) = indx(iloc,jloc,nloc-1,mipan,mjpan)
        end if

        iqq = js_g(iqg)    ! Global neighbour index
        rproc = mg_qproc(iqq) ! Processor that has this point
        if ( rproc == myid ) then ! Just copy the value
          ! Convert global iqq to local value
          call indv_mpix(iqq,iloc,jloc,nloc,mil_g,mioff,mjoff,noff)
          mg(g)%is(iq) = indx(iloc,jloc,nloc-1,mipan,mjpan)
        end if

        iqq = je_g(iqg)    ! Global neighbour index
        rproc = mg_qproc(iqq) ! Processor that has this point
        if ( rproc == myid ) then ! Just copy the value
          ! Convert global iqq to local value
          call indv_mpix(iqq,iloc,jloc,nloc,mil_g,mioff,mjoff,noff)
          mg(g)%ie(iq) = indx(iloc,jloc,nloc-1,mipan,mjpan)
        end if

        iqq = jw_g(iqg)    ! Global neighbour index
        rproc = mg_qproc(iqq) ! Processor that has this point
        if ( rproc == myid ) then ! Just copy the value
          ! Convert global iqq to local value
          call indv_mpix(iqq,iloc,jloc,nloc,mil_g,mioff,mjoff,noff)
          mg(g)%iw(iq) = indx(iloc,jloc,nloc-1,mipan,mjpan)
        end if

        ! Note that the model only needs a limited set of the diagonal
        ! index arrays
        iqq = jne_g(iqg)    ! Global neighbour index
        rproc = mg_qproc(iqq) ! Processor that has this point
        if ( rproc == myid ) then ! Just copy the value
          ! Convert global iqq to local value
          call indv_mpix(iqq,iloc,jloc,nloc,mil_g,mioff,mjoff,noff)
          mg(g)%ine(iq) = indx(iloc,jloc,nloc-1,mipan,mjpan)
        end if

        iqq = jse_g(iqg)    ! Global neighbour index
        rproc = mg_qproc(iqq) ! Processor that has this point
        if ( rproc == myid ) then ! Just copy the value
          ! Convert global iqq to local value
          call indv_mpix(iqq,iloc,jloc,nloc,mil_g,mioff,mjoff,noff)
          mg(g)%ise(iq) = indx(iloc,jloc,nloc-1,mipan,mjpan)
        end if

        iqq = jnw_g(iqg)    ! Global neighbour index
        rproc = mg_qproc(iqq) ! Processor that has this point
        if ( rproc == myid ) then ! Just copy the value
          ! Convert global iqq to local value
          call indv_mpix(iqq,iloc,jloc,nloc,mil_g,mioff,mjoff,noff)
          mg(g)%inw(iq) = indx(iloc,jloc,nloc-1,mipan,mjpan)
        end if

        iqq = jsw_g(iqg)    ! Global neighbour index
        rproc = mg_qproc(iqq) ! Processor that has this point
        if ( rproc == myid ) then ! Just copy the value
          ! Convert global iqq to local value
          call indv_mpix(iqq,iloc,jloc,nloc,mil_g,mioff,mjoff,noff)
          mg(g)%isw(iq) = indx(iloc,jloc,nloc-1,mipan,mjpan)
        end if

      end do
    end do
  end do


  ! Calculate local indices in halo
  jextra=2*(mipan+mjpan+2)*npan
  iext=0
  mg_bnds(g,:)%rlen=0
  mg_bnds(g,:)%slen=0
  mg_bnds(g,:)%rlenx=0
  mg_bnds(g,:)%slenx=0
  do n=1,npan

    !     Start with N edge
    j=mjpan
    do i=1,mipan
      iq = indx(i+mioff(n),j+mjoff(n),n-noff,mil_g,mil_g)
      iqq = jn_g(iq)
      ! Which processor has this point
      rproc = mg_qproc(iqq)
      if ( rproc == myid ) cycle ! Don't add points already on this proc.
      iql = indx(i,j,n-1,mipan,mjpan)  !  Local index
      ! Add this point to request list
      call mgcheck_bnds_alloc(g, rproc, jextra, iext)
      mg_bnds(g,rproc)%rlen = mg_bnds(g,rproc)%rlen + 1
      mg_bnds(g,rproc)%request_list(mg_bnds(g,rproc)%rlen) = iqq
      ! Increment extended region index
      iext = iext + 1
      mg_bnds(g,rproc)%unpack_list(mg_bnds(g,rproc)%rlen) = iext
      mg(g)%in(iql) = mg(g)%ifull+iext
    end do

    !     E edge
    i = mipan
    do j=1,mjpan
      iq = indx(i+mioff(n),j+mjoff(n),n-noff,mil_g,mil_g)
      iqq = je_g(iq)
      ! Which processor has this point
      rproc = mg_qproc(iqq)
      if ( rproc == myid ) cycle ! Don't add points already on this proc.
      iql = indx(i,j,n-1,mipan,mjpan)  !  Local index
      ! Add this point to request list
      call mgcheck_bnds_alloc(g, rproc, jextra, iext)
      mg_bnds(g,rproc)%rlen = mg_bnds(g,rproc)%rlen + 1
      mg_bnds(g,rproc)%request_list(mg_bnds(g,rproc)%rlen) = iqq
      ! Increment extended region index
      iext = iext + 1
      mg_bnds(g,rproc)%unpack_list(mg_bnds(g,rproc)%rlen) = iext
      mg(g)%ie(iql) = mg(g)%ifull+iext
    end do

    !     W edge
    i = 1
    do j=1,mjpan
      iq = indx(i+mioff(n),j+mjoff(n),n-noff,mil_g,mil_g)
      iqq = jw_g(iq)
      ! Which processor has this point
      rproc = mg_qproc(iqq)
      if ( rproc == myid ) cycle ! Don't add points already on this proc.
      iql = indx(i,j,n-1,mipan,mjpan)  !  Local index
      ! Add this point to request list
      call mgcheck_bnds_alloc(g, rproc, jextra, iext)
      mg_bnds(g,rproc)%rlen = mg_bnds(g,rproc)%rlen + 1
      mg_bnds(g,rproc)%request_list(mg_bnds(g,rproc)%rlen) = iqq
      ! Increment extended region index
      iext = iext + 1
      mg_bnds(g,rproc)%unpack_list(mg_bnds(g,rproc)%rlen) = iext
      mg(g)%iw(iql) = mg(g)%ifull+iext
    end do

    !     S edge
    j=1
    do i=1,mipan
      iq = indx(i+mioff(n),j+mjoff(n),n-noff,mil_g,mil_g)
      iqq = js_g(iq)
      ! Which processor has this point
      rproc = mg_qproc(iqq)
      if ( rproc == myid ) cycle ! Don't add points already on this proc.
      iql = indx(i,j,n-1,mipan,mjpan)  !  Local index
      ! Add this point to request list
      call mgcheck_bnds_alloc(g, rproc, jextra, iext)
      mg_bnds(g,rproc)%rlen = mg_bnds(g,rproc)%rlen + 1
      mg_bnds(g,rproc)%request_list(mg_bnds(g,rproc)%rlen) = iqq
      ! Increment extended region index
      iext = iext + 1
      mg_bnds(g,rproc)%unpack_list(mg_bnds(g,rproc)%rlen) = iext
      mg(g)%is(iql) = mg(g)%ifull+iext
    end do
  end do ! n=1,npan

  mg_bnds(g,:)%rlenx = mg_bnds(g,:)%rlen  ! so that they're appended.
      
! Now handle the special corner values that need to be remapped
! This adds to rlen, so needs to come before the _XX stuff.
  do n=1,npan
    ! NE
    iq = indx(mipan,mjpan,n-1,mipan,mjpan)
    iqg = indx(mipan+mioff(n),mjpan+mjoff(n),n-noff,mil_g,mil_g)
    iqq = jne_g(iqg)
    ! Which processor has this point
    rproc = mg_qproc(iqq)
    if ( rproc /= myid ) then ! Add to list
      call mgcheck_bnds_alloc(g, rproc, jextra, iext)
      mg_bnds(g,rproc)%rlenx = mg_bnds(g,rproc)%rlenx + 1
      mg_bnds(g,rproc)%request_list(mg_bnds(g,rproc)%rlenx) = iqq
      ! Increment extended region index
      iext = iext + 1
      mg_bnds(g,rproc)%unpack_list(mg_bnds(g,rproc)%rlenx) = iext
      mg(g)%ine(iq) = mg(g)%ifull+iext
    end if

    ! SE
    iq = indx(mipan,1,n-1,mipan,mjpan)
    iqg = indx(mipan+mioff(n),1+mjoff(n),n-noff,mil_g,mil_g)
    iqq = jse_g(iqg)
    ! Which processor has this point
    rproc = mg_qproc(iqq)
    if ( rproc /= myid ) then ! Add to list
      call mgcheck_bnds_alloc(g, rproc, jextra, iext)
      mg_bnds(g,rproc)%rlenx = mg_bnds(g,rproc)%rlenx + 1
      mg_bnds(g,rproc)%request_list(mg_bnds(g,rproc)%rlenx) = iqq
      ! Increment extended region index
      iext = iext + 1
      mg_bnds(g,rproc)%unpack_list(mg_bnds(g,rproc)%rlenx) = iext
      mg(g)%ise(iq) = mg(g)%ifull+iext
    end if

    ! WN
    iq = indx(1,mjpan,n-1,mipan,mjpan)
    iqg = indx(1+mioff(n),mjpan+mjoff(n),n-noff,mil_g,mil_g)
    iqq = jnw_g(iqg)
    ! Which processor has this point
    rproc = mg_qproc(iqq)
    if ( rproc /= myid ) then ! Add to list
      call mgcheck_bnds_alloc(g, rproc, jextra, iext)
      mg_bnds(g,rproc)%rlenx = mg_bnds(g,rproc)%rlenx + 1
      mg_bnds(g,rproc)%request_list(mg_bnds(g,rproc)%rlenx) = iqq
      ! Increment extended region index
      iext = iext + 1
      mg_bnds(g,rproc)%unpack_list(mg_bnds(g,rproc)%rlenx) = iext
      mg(g)%inw(iq) = mg(g)%ifull+iext
    end if

    ! SW
    iq = indx(1,1,n-1,mipan,mjpan)
    iqg = indx(1+mioff(n),1+mjoff(n),n-noff,mil_g,mil_g)
    iqq = jsw_g(iqg)
    ! Which processor has this point
    rproc = mg_qproc(iqq)
    if ( rproc /= myid ) then ! Add to list
      call mgcheck_bnds_alloc(g, rproc, jextra, iext)
      mg_bnds(g,rproc)%rlenx = mg_bnds(g,rproc)%rlenx + 1
      mg_bnds(g,rproc)%request_list(mg_bnds(g,rproc)%rlenx) = iqq
      ! Increment extended region index
      iext = iext + 1
      mg_bnds(g,rproc)%unpack_list(mg_bnds(g,rproc)%rlenx) = iext
      mg(g)%isw(iq) = mg(g)%ifull+iext
    end if

  end do

! Set up the diagonal index arrays. Most of the points here will have
! already been added to copy lists above. The corners are handled
! separately here. This means some points may be copied twice but it's
! a small overhead.
  do n=1,npan
    do j=1,mjpan
      do i=1,mipan
        iq = indx(i,j,n-1,mipan,mjpan)   ! Local
        ! Except at corners, ien = ine etc.
        if ( i > 1 ) then
          mg(g)%inw(iq) = mg(g)%in(mg(g)%iw(iq))
          mg(g)%isw(iq) = mg(g)%is(mg(g)%iw(iq))
        else
          if ( j < mjpan ) mg(g)%inw(iq) = mg(g)%iw(mg(g)%in(iq))
          if ( j > 1 )     mg(g)%isw(iq) = mg(g)%iw(mg(g)%is(iq))
        end if
        if ( i < mipan ) then
          ! ie will be defined
          mg(g)%ine(iq) = mg(g)%in(mg(g)%ie(iq))
          mg(g)%ise(iq) = mg(g)%is(mg(g)%ie(iq))
        else
          ! i = ipan, ie will have been remapped
          if ( j > 1 )     mg(g)%ise(iq) = mg(g)%ie(mg(g)%is(iq))
          if ( j < mjpan ) mg(g)%ine(iq) = mg(g)%ie(mg(g)%in(iq))
        end if
        ! if ( j > 1 ) then
        !   mg(g)%ies(iq) = mg(g)%ie(mg(g)%is(iq))
        !   mg(g)%iws(iq) = mg(g)%iw(mg(g)%is(iq))
        ! else
        !   if ( i < mipan ) mg(g)%ies(iq)=mg(g)%is(mg(g)%ie(iq))
        !   if ( i > 1 )     mg(g)%iws(iq)=mg(g)%is(mg(g)%iw(iq))
        ! end if
        ! if ( j < mjpan ) then
        !   mg(g)%ien(iq) = mg(g)%ie(mg(g)%in(iq))
        !   mg(g)%iwn(iq) = mg(g)%iw(mg(g)%in(iq))
        ! else
        !   if ( i < mipan) mg(g)%ien(iq) = mg(g)%in(mg(g)%ie(iq))
        !   if ( i > 1 )    mg(g)%iwn(iq) = mg(g)%in(mg(g)%iw(iq))
        ! end if
      end do
    end do
  end do


! Now, for each processor send the list of points I want.
! The state of being a neighbour is reflexive so only expect to
! recv from those processors I send to (are there grid arrangements for
! which this would not be true?)
! Get the complete request lists by using rlen2
  nreq = 0
  do iproc = 1,nproc-1  !
    sproc = modulo(myid+iproc,nproc)  ! Send to
    if (mg_bnds(g,sproc)%rlenx > 0 ) then
      nreq = nreq + 1
      ! Use the maximum size in the recv call.
      call MPI_IRecv( mg_bnds(g,sproc)%send_list(1), mg_bnds(g,sproc)%len, &
                 MPI_INTEGER, sproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
    end if
  end do
  do iproc = 1,nproc-1  !
    sproc = modulo(myid+iproc,nproc)  ! Send to
    if (mg_bnds(g,sproc)%rlenx > 0 ) then
      ! Send list of requests
      nreq = nreq + 1
      call MPI_ISend( mg_bnds(g,sproc)%request_list(1), mg_bnds(g,sproc)%rlenx, &
                 MPI_INTEGER, sproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
    end if
  end do      
  if ( nreq > 0 ) then
    call MPI_Waitall(nreq,ireq,status,ierr)
  end if

! Now get the actual sizes from the status
  nreq = 0
  do iproc = 1,nproc-1  !
    sproc = modulo(myid+iproc,nproc)  ! Recv from
    if (mg_bnds(g,sproc)%rlenx > 0 ) then
      ! First half of nreq are recv
      nreq = nreq + 1
      call MPI_Get_count(status(1,nreq), MPI_INTEGER, ncount, ierr)
      ! This the number of points I have to send to rproc.
      mg_bnds(g,sproc)%slenx = ncount
    end if
  end do

! For rlen, just communicate the lengths. The indices have 
! already been taken care of.
  nreq = 0
  do iproc = 1,nproc-1  !
    ! Send and recv from same proc
    sproc = modulo(myid+iproc,nproc)  ! Send to
    if (mg_bnds(g,sproc)%rlenx > 0 ) then
      nreq = nreq + 1
      call MPI_IRecv( mg_bnds(g,sproc)%slen, 1, &
                 MPI_INTEGER, sproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
    end if
  end do
  do iproc = 1,nproc-1  !
    ! Send and recv from same proc
    sproc = modulo(myid+iproc,nproc)  ! Send to
    if (mg_bnds(g,sproc)%rlenx > 0 ) then
      nreq = nreq + 1
      call MPI_ISend( mg_bnds(g,sproc)%rlen, 1, &
                 MPI_INTEGER, sproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
    end if
  end do
  if ( nreq > 0 ) then
    call MPI_Waitall(nreq,ireq,status,ierr)
  end if

  ! At the moment send_lists use global indices. Convert these to local.
  do iproc = 1,nproc-1  !
    sproc = modulo(myid+iproc,nproc)  ! Send to
    do iq=1,mg_bnds(g,sproc)%slenx
      ! send_list(iq) is global point index, i, j, n are local
      call indv_mpix(mg_bnds(g,sproc)%send_list(iq),i,j,n,mil_g,mioff,mjoff,noff)
      mg_bnds(g,sproc)%send_list(iq) = indx(i,j,n-1,mipan,mjpan)
    end do
  end do
  do iq=1,mg_bnds(g,myid)%rlenx
    call indv_mpix(mg_bnds(g,myid)%request_list(iq),i,j,n,mil_g,mioff,mjoff,noff)
    mg_bnds(g,myid)%request_list(iq) = indx(i,j,n-1,mipan,mjpan)
  end do

  ! reduce array size where possible
  do iproc=0,nproc-1
    xlen=max(mg_bnds(g,iproc)%rlenx,mg_bnds(g,iproc)%slenx)
    if (mg_bnds(g,iproc)%len>xlen) then
      dum(1:xlen)=mg_bnds(g,iproc)%request_list(1:xlen)
      deallocate(mg_bnds(g,iproc)%request_list)
      allocate(mg_bnds(g,iproc)%request_list(xlen))
      mg_bnds(g,iproc)%request_list(1:xlen)=dum(1:xlen)
      dum(1:xlen)=mg_bnds(g,iproc)%unpack_list(1:xlen)
      deallocate(mg_bnds(g,iproc)%unpack_list)
      allocate(mg_bnds(g,iproc)%unpack_list(xlen))
      mg_bnds(g,iproc)%unpack_list(1:xlen)=dum(1:xlen)
      dum(1:xlen)=mg_bnds(g,iproc)%send_list(1:xlen)
      deallocate(mg_bnds(g,iproc)%send_list)
      allocate(mg_bnds(g,iproc)%send_list(xlen))
      mg_bnds(g,iproc)%send_list(1:xlen)=dum(1:xlen)
      mg_bnds(g,iproc)%len=xlen
    end if

    ! set-up buffers
    xlev=max(kl,ol)
    if (mg_bnds(g,iproc)%rlenx>0) then
      allocate(mg_bnds(g,iproc)%rbuf(xlev*mg_bnds(g,iproc)%rlenx))
    end if
    if (mg_bnds(g,iproc)%slenx>0) then
      allocate(mg_bnds(g,iproc)%sbuf(xlev*mg_bnds(g,iproc)%slenx))
    end if
  end do

end if


! calculate colours
if (g==mg_maxlevel) then
  
  ! always a three colour mask for coarse grid
  do n=0,npanels
    do j=1,mil_g
      do i=1,mil_g
        iq=indx(i,j,n,mil_g,mil_g)

        jx = mod(i+j+n*mil_g,2)
        select case( n+jx*(npanels+1) )
          case(0,1,3,4)
            mg_colourmask(iq)=1
          case(2,5,6,9)
            mg_colourmask(iq)=2
          case(7,8,10,11)
            mg_colourmask(iq)=3
        end select
      end do
    end do
  end do
  
  mg_colour_np=max(count(mg_colourmask==1),count(mg_colourmask==2),count(mg_colourmask==3))
  allocate(col_iq(mg_colour_np,3),col_iqn(mg_colour_np,3),col_iqe(mg_colour_np,3),col_iqs(mg_colour_np,3),col_iqw(mg_colour_np,3))
  
  mg_ifullc=0
  col_iq=0
  col_iqn=0
  col_iqe=0
  col_iqs=0
  col_iqw=0
  do iq=1,mg(g)%ifull
    nc=mg_colourmask(iq)
    mg_ifullc(nc)=mg_ifullc(nc)+1
    iqq=mg_ifullc(nc)
    col_iq(iqq,nc)=iq
    col_iqn(iqq,nc)=mg(g)%in(iq)
    col_iqe(iqq,nc)=mg(g)%ie(iq)
    col_iqs(iqq,nc)=mg(g)%is(iq)
    col_iqw(iqq,nc)=mg(g)%iw(iq)
  end do

end if

return
end subroutine mg_index

subroutine mgcheck_bnds_alloc(g,iproc,jextra,iext)

use cc_mpi

implicit none

integer, intent(in) :: g,iproc,jextra,iext

if (mg_bnds(g,iproc)%len<=0) then
  allocate(mg_bnds(g,iproc)%request_list(jextra))
  allocate(mg_bnds(g,iproc)%unpack_list(jextra))
  allocate(mg_bnds(g,iproc)%send_list(jextra))
  mg_bnds(g,iproc)%len=jextra
else
  if (iext>jextra) then
    write(6,*) "ERROR: MG grid undersized in mgcheck_bnds_alloc"
    write(6,*) "iext,jextra,g,iproc,myid ",iext,jextra,g,iproc,myid
    stop
  end if
end if

return
end subroutine mgcheck_bnds_alloc

subroutine indv_mpix(iq, i, j, n, mil_g, mioff, mjoff, mnoff)

implicit none

include 'newmpar.h'

integer , intent(in) :: iq, mil_g, mnoff
integer, dimension(npanels+1), intent(in) :: mioff, mjoff
integer , intent(out) :: i
integer , intent(out) :: j
integer , intent(out) :: n
integer :: ierr,ierr2

! Calculate local i, j, n from global iq

! Global i, j, n
n = (iq - 1)/(mil_g*mil_g)
j = 1 + (iq - n*mil_g*mil_g - 1)/mil_g
i = iq - (j - 1)*mil_g - n*mil_g*mil_g

! Reduced to values on my processor
n = n + mnoff  
j = j - mjoff(n)
i = i - mioff(n)

return
end subroutine indv_mpix

function indx(i,j,n,il,jl) result(iq)

implicit none

integer, intent(in) :: i, j, n, il, jl
integer iq

iq = i+(j-1)*il+n*il*jl

return
end function indx


end module mgsolve