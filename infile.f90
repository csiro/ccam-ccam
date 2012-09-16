module infile
      
! This module contains routines for reading netcdf files,
! vertical interpolation and some calendar functions
      
! This version supports parallel localhist input files
            
private
public histrd1,histrd4s,vertint,datefix,getzinp,ncmsg
public histopen,histclose
      
integer, parameter :: ncidmax=10000 ! maximum number of parallel input files
integer, dimension(:), allocatable, save :: pncid
integer, dimension(:), allocatable, save :: pncidold
integer, dimension(:), allocatable, save :: pnoff
integer, dimension(:,:), allocatable, save :: pioff,pjoff
integer, save :: mynproc,fnproc
integer, save :: pil_g,pjl_g,pil,pjl,pnpan
logical, save :: ptest
      
contains

!--------------------------------------------------------------
! Interface for reading 2D+time fields
subroutine histrd1(ncid,iarchi,ier,name,ik,jk,var,ifull)
      
use cc_mpi
      
implicit none
      
include 'mpif.h'
include 'parm.h'
      
integer ncid,iarchi,ier,ik,jk,ifull
character(len=*) name
real, dimension(:), intent(inout) :: var ! may be dummy argument from myid/=0

if (ifull/=ik*jk.and.ptest) then
  ! read local arrays without gather and distribute
       
  if (size(var)/=ifull) then
    write(6,*) "ERROR: Incorrect use of dummy var in histrd1"
    call MPI_Abort(MPI_COMM_WORLD,-1,ier)
  end if

  call hr1p(iarchi,ier,name,.true.,var)

  if(ier==0.and.mod(ktau,nmaxpr)==0.and.myid==0)then
    write(6,'("done histrd1 ",a8,i4,i3)') name,ier,iarchi
  end if

else        
      
  ! split up processors to save memory.  No need to allocate global
  ! arrays on myid/=0.
  if (myid==0) then
    if (size(var)/=ifull) then
      write(6,*) "ERROR: Incorrect use of dummy var in histrd1"
      call MPI_Abort(MPI_COMM_WORLD,-1,ier)
    end if

    call hr1a(ncid,iarchi,ier,name,ik,jk,var,ifull)
  else
    if(ifull/=ik*jk)then
      ! read local arrays with gather and distribute

      if (size(var)/=ifull) then
        write(6,*) "ERROR: Incorrect use of dummy var in histrd"
        call MPI_Abort(MPI_COMM_WORLD,-1,ier)
      end if

      call hr1p(iarchi,ier,name,.false.)
      call ccmpi_distribute(var)

    else
      ! read global arrays for myid==0

      call hr1p(iarchi,ier,name,.false.)

    end if
  end if
end if
      
return
end subroutine histrd1    

!--------------------------------------------------------------------
! Gather 2D+time fields on myid==0
subroutine hr1a(ncid,iarchi,ier,name,ik,jk,var,ifull)
      
use cc_mpi
      
implicit none
      
include 'parm.h'
      
integer, intent(in) :: ncid,iarchi,ik,jk,ifull
integer, intent(out) :: ier
character(len=*), intent(in) :: name
real, dimension(ifull), intent(inout) :: var
real, dimension(ik*jk) :: globvar
real vmax,vmin

call hr1p(iarchi,ier,name,.false.,globvar)

if(ier==0.and.mod(ktau,nmaxpr)==0)then
  vmax = maxval(globvar)
  vmin = minval(globvar)
  write(6,'("done histrd1 ",a8,i4,i3,3e14.6)') name,ier,iarchi,vmin,vmax,globvar(id+(jd-1)*ik)
end if
      
if (ifull==ik*jk) then
  ! read global arrays for myid==0
  if (ier==0) then
    var(:)=globvar(:)
  end if
else
  ! read local arrays with gather and distribute
  call ccmpi_distribute(var,globvar)
endif

return
end subroutine hr1a  

!--------------------------------------------------------------
! This subroutine reads 2D+time input files
      
! when qtest=.true. the input grid decomposition should
! match the current processor decomposition.  We can then
! skip the MPI gather and distribute steps.
subroutine hr1p(iarchi,ier,name,qtest,var)
      
use cc_mpi
      
implicit none

include 'newmpar.h'
include 'mpif.h'
include 'netcdf.inc'

integer, intent(in) :: iarchi
integer, intent(out) :: ier
integer(kind=2), dimension(pil*pjl*pnpan) :: ivar
integer, dimension(3) :: start,count
integer, dimension(MPI_STATUS_SIZE) :: status
integer idv,ipf,ip,jpf,jpmax,jpmod,n
integer nctype,ierb,ca,cb
integer no,i,j,iq,iqi
integer :: itag=0
real, dimension(:), intent(inout), optional :: var
real, dimension(pil*pjl*pnpan) :: rvar
real addoff,sf
logical, intent(in) :: qtest
character(len=*), intent(in) :: name

if (qtest) then
  if (mynproc>1) then
    write(6,*) "ERROR: Invalid use of qtest"
    call MPI_Abort(MPI_COMM_WORLD,-1,ierb)
  end if
  if (.not.present(var)) then
    write(6,*) "ERROR: Missing var in hr1p with qtest"
    call MPI_Abort(MPI_COMM_WORLD,-1,ierb)
  end if
  if (size(var)/=pil*pjl*pnpan) then
    write(6,*) "ERROR: Invalid var size in hr1p with qtest"
    write(6,*) "size(var),pil,pjl,pnpan ",size(var),pil,pjl,pnpan
    call MPI_Abort(MPI_COMM_WORLD,-1,ierb)
  end if
else
  if (myid==0) then
    if (.not.present(var)) then
      write(6,*) "ERROR: Missing var in hr1p"
      call MPI_Abort(MPI_COMM_WORLD,-1,ierb)
    end if
    if (size(var)/=pil_g*pjl_g) then
      write(6,*) "ERROR: Invalid var size in hr1p"
      write(6,*) "size(var),pil_g,pjl_g ",size(var),pil_g,pjl_g
      call MPI_Abort(MPI_COMM_WORLD,-1,ierb)
    end if
  else
    if (present(var)) then
      write(6,*) "ERROR: Invalid use of var in hr1p"
      call MPI_Abort(MPI_COMM_WORLD,-1,ierb)
    end if
  end if
end if

start = (/ 1, 1, iarchi /)
count = (/ pil, pjl*pnpan, 1 /)
ier=0
      
do ipf=0,mynproc-1
  rvar=0.

  ! get variable idv
  idv=ncvid(pncid(ipf),name,ier)
  if(ier/=0)then
    if (myid==0.and.ipf==0) then
      write(6,*) '***absent field for ncid,name,idv,ier: ',pncid(0),name,idv,ier
    end if
  else
    ! obtain scaling factors and offsets from attributes
    call ncagt(pncid(ipf),idv,'add_offset',addoff,ier)
    if (ier/=0) addoff=0.
    call ncagt(pncid(ipf),idv,'scale_factor',sf,ier)
    if (ier/=0) sf=1.
    ! read in all data
    ier=nf_inq_vartype(pncid(ipf),idv,nctype)
    call ncmsg("nctype",ier)
    select case(nctype)
      case(nf_float)
        call ncvgt(pncid(ipf),idv,start,count,rvar,ier)
      case(nf_short)
        call ncvgt(pncid(ipf),idv,start,count,ivar,ier)
        rvar(:)=real(ivar(:))
      case DEFAULT
        write(6,*) "ERROR: Unknown NetCDF format"
        call MPI_Abort(MPI_COMM_WORLD,-1,ierb)
    end select
    call ncmsg("hr1p",ier)
    ! unpack data
    rvar=rvar*sf+addoff
  end if ! ier
      
  if (qtest) then
    ! restart file
    var=rvar
  else
    ! mesonest file
     
    ! recompose grid
    if (myid==0) then
      ! recompose own grid
      ip=ipf*nproc
      do n=0,pnpan-1
        no=n-pnoff(ip)+1
        ca=pioff(ip,no)
        cb=pjoff(ip,no)+no*pil_g
        do j=1,pjl
          do i=1,pil
            iq=i+ca+(j+cb-1)*pil_g
            iqi=i+(j-1)*pil+n*pil*pjl
            var(iq)=rvar(iqi)
          end do
        end do
      end do
      ! recompose other processors
      jpmax=nproc
      jpmod=mod(fnproc,nproc)
      if (ipf==mynproc-1.and.jpmod/=0) then
        jpmax=jpmod
      end if
      do jpf=1,jpmax-1
        call MPI_Recv(rvar,pil*pjl*pnpan,MPI_REAL,jpf,itag,MPI_COMM_WORLD,status,ierb)
        ip=ipf*nproc+jpf
        do n=0,pnpan-1
          no=n-pnoff(ip)+1
          ca=pioff(ip,no)
          cb=pjoff(ip,no)+no*pil_g
          do j=1,pjl
            do i=1,pil
              iq=i+ca+(j+cb-1)*pil_g
              iqi=i+(j-1)*pil+n*pil*pjl
              var(iq)=rvar(iqi)
            end do
          end do
        end do
      end do
    else
      call MPI_SSend(rvar,pil*pjl*pnpan,MPI_REAL,0,itag,MPI_COMM_WORLD,ierb)
    end if
  end if ! qtest

end do ! ipf

return
end subroutine hr1p

!--------------------------------------------------------------   
! Interface for reading 3D+time fields
subroutine histrd4s(ncid,iarchi,ier,name,ik,jk,kk,var,ifull)
      
use cc_mpi
      
implicit none
      
include 'mpif.h'
include 'parm.h'
      
integer, intent(in) :: ncid,iarchi,ik,jk,kk,ifull
integer, intent(out) :: ier
character(len=*), intent(in) :: name
real, dimension(:), intent(inout) :: var ! may be dummy argument from myid/=0
      
if (ifull/=ik*jk.and.ptest) then
  ! read local arrays without gather and distribute
       
  if (size(var)/=ifull) then
    write(6,*) "ERROR: Incorrect use of dummy var in histrd4s"
    call MPI_Abort(MPI_COMM_WORLD,-1,ier)
  end if

  call hr4p(iarchi,ier,name,kk,.true.,var)

  if(ier==0.and.mod(ktau,nmaxpr)==0.and.myid==0)then
    write(6,'("done histrd4 ",a8,i3,i4,i3)') name,kk,ier,iarchi
  endif

else        
      
  ! split up processors to save memory.  No need to allocate global
  ! arrays on myid/=0.
  if (myid==0) then
    if (size(var)/=ifull) then
      write(6,*) "ERROR: Incorrect use of dummy var in histrd4s"
      call MPI_Abort(MPI_COMM_WORLD,-1,ier)
    end if
 
    call hr4sa(ncid,iarchi,ier,name,ik,jk,kk,var,ifull)
  else
    if(ifull/=ik*jk)then
      ! read local arrays with gather and distribute

      if (size(var)/=ifull) then
        write(6,*) "ERROR: Incorrect use of dummy var in histrd"
        call MPI_Abort(MPI_COMM_WORLD,-1,ier)
      end if
 
      call hr4p(iarchi,ier,name,kk,.false.)
      call ccmpi_distribute(var)

    else
      ! read global arrays for myid==0

      call hr4p(iarchi,ier,name,kk,.false.)

    end if
  end if
end if

return
end subroutine histrd4s      

!--------------------------------------------------------------------
! Gather 3D+time fields on myid==0

subroutine hr4sa(ncid,iarchi,ier,name,ik,jk,kk,var,ifull)
      
use cc_mpi
      
implicit none
      
include 'parm.h'
      
integer, intent(in) :: ncid,iarchi,ik,jk,kk,ifull
integer, intent(out) :: ier
character(len=*), intent(in) :: name
real, dimension(ifull) :: var
real, dimension(ik*jk) :: globvar
real vmax,vmin
      
call hr4p(iarchi,ier,name,kk,.false.,globvar)     

if(ier==0.and.mod(ktau,nmaxpr)==0)then
  vmax = maxval(globvar)
  vmin = minval(globvar)
  write(6,'("done histrd4s ",a6,i3,i4,i3,3f12.4)') name,kk,ier,iarchi,vmin,vmax,globvar(id+(jd-1)*ik)
end if

! Have to return correct value of ier on all processes because it's 
! used for initialisation in calling routine
if(ifull==ik*jk)then
  ! read global arrays for myid==0
  if (ier==0) then
    var(:)=globvar(:)
  end if
else
  ! read local arrays with gather and distribute
  call ccmpi_distribute(var,globvar)
endif

return
end subroutine hr4sa      

!--------------------------------------------------------------
! This subroutine reads 3D+time input files

! when qtest=.true. the input grid decomposition should
! match the current processor decomposition.  We can then
! skip the MPI gather and distribute steps.
subroutine hr4p(iarchi,ier,name,kk,qtest,var)
      
use cc_mpi
      
implicit none

include 'newmpar.h'
include 'mpif.h'
include 'netcdf.inc'

integer, intent(in) :: iarchi,kk
integer, intent(out) :: ier
integer(kind=2), dimension(pil*pjl*pnpan) :: ivar
integer, dimension(4) :: start,count
integer, dimension(MPI_STATUS_SIZE) :: status
integer idv,ipf,ip,jpf,jpmax,jpmod,n
integer nctype,ierb,ca,cb
integer no,i,j,iq,iqi
integer :: itag=0
real, dimension(:), intent(inout), optional :: var
real, dimension(pil*pjl*pnpan) :: rvar
real addoff,sf
logical, intent(in) :: qtest
character(len=*), intent(in) :: name

if (qtest) then
  if (mynproc>1) then
     write(6,*) "ERROR: Invalid use of qtest"
     call MPI_Abort(MPI_COMM_WORLD,-1,ierb)
  end if
  if (.not.present(var)) then
    write(6,*) "ERROR: Missing var in hr4p with qtest"
    call MPI_Abort(MPI_COMM_WORLD,-1,ierb)
  end if
  if (size(var)/=pil*pjl*pnpan) then
    write(6,*) "ERROR: Invalid var size in hr4p with qtest"
    write(6,*) "size(var),pil,pjl,pnpan ",size(var),pil,pjl,pnpan
    call MPI_Abort(MPI_COMM_WORLD,-1,ierb)
  end if
else
  if (myid==0) then
    if (.not.present(var)) then
      write(6,*) "ERROR: Missing var in hr4p"
      call MPI_Abort(MPI_COMM_WORLD,-1,ierb)
    end if
    if (size(var)/=pil_g*pjl_g) then
      write(6,*) "ERROR: Invalid var size in hr4p"
      write(6,*) "size(var),pil_g,pjl_g ",size(var),pil_g,pjl_g
      call MPI_Abort(MPI_COMM_WORLD,-1,ierb)
    end if
  else
    if (present(var)) then
      write(6,*) "ERROR: Invalid use of var in hr4p"
      call MPI_Abort(MPI_COMM_WORLD,-1,ierb)
    end if
  end if
end if

start = (/ 1, 1, kk, iarchi /)
count = (/ pil, pjl*pnpan, 1, 1 /)
ier=0
      
do ipf=0,mynproc-1
  rvar=0.

  ! get variable idv
  idv=ncvid(pncid(ipf),name,ier)
  if(ier/=0)then
    if (myid==0.and.ipf==0) then
      write(6,*) '***absent field for ncid,name,idv,ier: ',pncid(0),name,idv,ier
    end if
  else
    ! obtain scaling factors and offsets from attributes
    call ncagt(pncid(ipf),idv,'add_offset',addoff,ier)
    if (ier/=0) addoff=0.
    call ncagt(pncid(ipf),idv,'scale_factor',sf,ier)
    if (ier/=0) sf=1.
    ! read in all data
    ier=nf_inq_vartype(pncid(ipf),idv,nctype)
    call ncmsg("nctype",ier)
    select case(nctype)
      case(nf_float)
        call ncvgt(pncid(ipf),idv,start,count,rvar,ier)
      case(nf_short)
        call ncvgt(pncid(ipf),idv,start,count,ivar,ier)
        rvar(:)=real(ivar(:))
      case DEFAULT
        write(6,*) "ERROR: Unknown NetCDF format"
        call MPI_Abort(MPI_COMM_WORLD,-1,ierb)
    end select
    call ncmsg("hr4p",ier)
    ! unpack data
    rvar=rvar*sf+addoff
  end if ! ier

  if (qtest) then
    ! restart file
    var=rvar
  else
    ! mesonest file
      
    ! recompose grid
    if (myid==0) then
      ! recompose own grid
      ip=ipf*nproc
      do n=0,pnpan-1
        no=n-pnoff(ip)+1
        ca=pioff(ip,no)
        cb=pjoff(ip,no)+no*pil_g
        do j=1,pjl
          do i=1,pil
            iq=i+ca+(j+cb-1)*pil_g
            iqi=i+(j-1)*pil+n*pil*pjl
            var(iq)=rvar(iqi)
          end do
        end do
      end do
      ! recompose grids from other processors
      jpmax=nproc
      jpmod=mod(fnproc,nproc)
      if (ipf==mynproc-1.and.jpmod/=0) then
        jpmax=jpmod
      end if
      do jpf=1,jpmax-1
        call MPI_Recv(rvar,pil*pjl*pnpan,MPI_REAL,jpf,itag,MPI_COMM_WORLD,status,ierb)
        ip=ipf*nproc+jpf
        do n=0,pnpan-1
          no=n-pnoff(ip)+1
          ca=pioff(ip,no)
          cb=pjoff(ip,no)+no*pil_g
          do j=1,pjl
            do i=1,pil
              iq=i+ca+(j+cb-1)*pil_g
              iqi=i+(j-1)*pil+n*pil*pjl
              var(iq)=rvar(iqi)
            end do
          end do
        end do
      end do
    else
      call MPI_SSend(rvar,pil*pjl*pnpan,MPI_REAL,0,itag,MPI_COMM_WORLD,ierb)
    end if
  end if ! qtest

end do ! ipf

return
end subroutine hr4p

!--------------------------------------------------------------
! Trap netcdf error messages
subroutine ncmsg(txt,ncstatus)

use cc_mpi

implicit none

include 'mpif.h'
include 'netcdf.inc'

integer, intent(in) :: ncstatus
integer ierr
character(len=*), intent(in) :: txt

if (ncstatus/=0) then
  write(6,*) txt," ",nf_strerror(ncstatus)
  call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
end if

return
end subroutine ncmsg

!--------------------------------------------------------------
! This subroutine opens parallel input files
subroutine histopen(ncid,ifile)
      
use cc_mpi
      
implicit none
      
include 'newmpar.h'
include 'mpif.h'
include 'netcdf.inc'
      
integer, parameter :: nihead   = 54
      
integer, dimension(nihead) :: ahead
integer, dimension(0:5) :: duma,dumb
integer, dimension(5) :: idum
integer, intent(out) :: ncid
integer ierr,resid,is,ipf,dmode
integer ipin,nxpr,nypr,dumr
character(len=*), intent(in) :: ifile
character(len=170) :: pfile
character(len=7) :: fdecomp
      
call ncpopt(0) ! turn off fatal errors on all processors
if (myid==0) then
  ! attempt to open single file with myid==0
  ierr=nf_open(ifile,nf_nowrite,ncid)
  fnproc=1
  dmode=0
      
  ! attempt to open parallel files
  if (ierr/=0) then
    write(pfile,"(a,'.',i4.4)") trim(ifile), 0
    ierr=nf_open(pfile,nf_nowrite,ncid)
    if (ierr/=0) then
      write(6,*) "ERROR: Cannot open ",pfile
      call ncmsg("pfile open",ierr)
    end if
    write(6,*) "Found parallel input file ",ifile
    ierr=nf_get_att_int(ncid,nf_global,"nproc",fnproc)
    call ncmsg("nproc",ierr)
    fdecomp=''
    ierr=nf_get_att_text(ncid,nf_global,"decomp",fdecomp)
    call ncmsg("decomp",ierr)
    select case(fdecomp)
      case('face')
        dmode=1
      case('uniform')
        dmode=2
      case default
        write(6,*) "ERROR: Unknown decomposition ",fdecomp
        call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
    end select
    call ncmsg("nproc",ierr)
    if (fnproc>ncidmax) then
      write(6,*) "ERROR: Exceeded maximum number of input files"
      call MPI_Abort(MPI_COMM_WORLD,-1,ierr)    
    end if
  else
    ierr=nf_get_att_int(ncid,nf_global,"nproc",dumr)
    if (ierr==0) then
      write(6,*) "ERROR: Incorrect base filename"
      call MPI_Abort(MPI_COMM_WORLD,-1,ierr)    
    end if
    write(6,*) "Found single input file ",ifile
  end if

  ierr=nf_get_att_int(ncid,nf_global,"int_header",ahead)
  call ncmsg("int_header",ierr)
  pil_g=ahead(1)
  pjl_g=ahead(2)
        
  if (allocated(pioff)) then
    deallocate(pioff,pjoff,pnoff)
  end if
  allocate(pioff(0:fnproc-1,0:5),pjoff(0:fnproc-1,0:5))
  allocate(pnoff(0:fnproc-1))
        
  select case(dmode)
    case(0) ! no decomposition
      pnpan=1
      pnoff=1
      pioff=0
      pjoff=0
      pil=pil_g
      pjl=pjl_g
    case(1) ! face decomposition
      pnpan=max(1,6/fnproc)
      do ipf=0,fnproc-1
        call face_set(pil,pjl,pnoff(ipf),duma,dumb,pnpan,pil_g,ipf,fnproc,nxpr,nypr)
        pioff(ipf,:)=duma
        pjoff(ipf,:)=dumb
      end do
    case(2) ! uniform decomposition
      pnpan=6
      do ipf=0,fnproc-1
        call uniform_set(pil,pjl,pnoff(ipf),duma,dumb,pnpan,pil_g,ipf,fnproc,nxpr,nypr)
        pioff(ipf,:)=duma
        pjoff(ipf,:)=dumb
      end do
  end select

  ptest=.false.
#ifdef uniform_decomp
  if (dmode==2) then
    if (nproc==fnproc) then
      if (pil_g==il_g.and.pjl_g==jl_g) then
        ptest=.true.
      end if
    end if
  end if
#else
  if (dmode==1) then
    if (nproc==fnproc) then
      if (pil_g==il_g.and.pjl_g==jl_g) then
        ptest=.true.
      end if
    end if
  end if
#endif

  write(6,*) "Found pil_g,pjl_g,fnproc,dmode,ptest ",pil_g,pjl_g,fnproc,dmode,ptest

  idum(1)=fnproc
  idum(2)=pil
  idum(3)=pjl
  idum(4)=pnpan
  idum(5)=0
  if (ptest) idum(5)=1
end if

! Broadcast file metadata
call MPI_Bcast(idum,5,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
fnproc=idum(1)
pil=idum(2)
pjl=idum(3)
pnpan=idum(4)
ptest=(idum(5)==1)

! calculate number of files to be read on this processor
mynproc=fnproc/nproc
resid=mod(fnproc,nproc)
if (myid<resid) mynproc=mynproc+1
      
if (allocated(pncid)) then
  deallocate(pncid)
end if
if (mynproc>0) then
  allocate(pncid(0:mynproc-1))
end if
      
is=0
if (myid==0) then 
  is=1
  pncid(0)=ncid
end if
      
do ipf=is,mynproc-1
  ipin=ipf*nproc+myid
  write(pfile,"(a,'.',i4.4)") trim(ifile), ipin
  ierr=nf_open(pfile,nf_nowrite,pncid(ipf))
  if (ierr/=0) then
    write(6,*) "ERROR: Cannot open ",pfile
    call ncmsg("pfile open",ierr)
  end if
end do
     
return
end subroutine histopen
      
!--------------------------------------------------------------
! This subroutine closes parallel input files
subroutine histclose
      
implicit none
    
include 'netcdf.inc'
     
integer ipf,ipin,plen,ierr
      
if (allocated(pncidold)) then
  plen=size(pncidold)
  do ipf=0,plen-1
    ierr=nf_close(pncidold(ipf))
  end do
  deallocate(pncidold)
end if
      
if (allocated(pncid)) then
  plen=size(pncid)
  allocate(pncidold(0:plen-1))
  pncidold=pncid
end if

return
end subroutine histclose

!--------------------------------------------------------------
! transforms 3d array from dimension kk in vertical to kl   jlm

! This version of vertint can interpolate from host models with
! a greater number of vertical levels than the nested model.
subroutine vertint(told,t,n,kk,sigin)

use cc_mpi, only : myid
use sigs_m

implicit none

include 'newmpar.h'
include 'parm.h'

integer, intent(in) :: kk,n
integer klapse,k,kin,iq
real, dimension(ifull,kl), intent(out) :: t
real, dimension(ifull,kl) :: told
real, dimension(kk), intent(in) :: sigin
real, dimension(kl) :: ka,kb,wta,wtb
      
if (kk==kl) then
  if (all(abs(sig-sigin)<0.0001)) then
     t=told
     return
   end if
end if
      
klapse=0
do k=1,kl
  if(sig(k)>=sigin(1))then
    ka(k)=2
    kb(k)=1
    wta(k)=0.
    wtb(k)=1.
    klapse=k   ! i.e. T lapse correction for k<=klapse
  elseif(sig(k)<=sigin(kk))then   ! at top
    ka(k)=kk
    kb(k)=kk-1
    wta(k)=1.
    wtb(k)=0.
  else
    do kin=2,kk-1
      if(sig(k)>sigin(kin)) exit
    enddo     ! kin loop
    ka(k)=kin
    kb(k)=kin-1
    wta(k)=(sigin(kin-1)-sig(k))/(sigin(kin-1)-sigin(kin))
    wtb(k)=(sig(k)-sigin(kin)  )/(sigin(kin-1)-sigin(kin))
  endif  !  (sig(k)>=sigin(1)) ... ...
enddo   !  k loop
if (myid==0) then
  write(6,*) 'in vertint kk,kl ',kk,kl
  write(6,"('sigin',10f7.4)") (sigin(k),k=1,kk)
  write(6,"('sig  ',10f7.4)") sig
  write(6,*) 'ka ',ka
  write(6,*) 'kb ',kb
  write(6,"('wta',10f7.4)") wta
  write(6,"('wtb',10f7.4)") wtb
endif   !  (myid==0)

do k=1,kl
  do iq=1,ifull
    ! N.B. "a" denotes "above", "b" denotes "below"
    t(iq,k)=wta(k)*told(iq,ka(k))+wtb(k)*told(iq,kb(k))
  enddo   ! iq loop
enddo    ! k loop

if(n==1.and.klapse/=0)then  ! for T lapse correction
  do k=1,klapse
    do iq=1,ifull
      ! assume 6.5 deg/km, with dsig=.1 corresponding to 1 km
      t(iq,k)=t(iq,k)+(sig(k)-sigin(1))*6.5/.1
    enddo   ! iq loop
  enddo    ! k loop
else if (n==2)then  ! for qg do a -ve fix
  t(:,:)=max(t(:,:),1.e-6)
else if (n==5)then  ! for qfg, qlg do a -ve fix
  t(:,:)=max(t(:,:),0.)
endif
      
return
end subroutine vertint

!--------------------------------------------------------------------
! This subroutine advances input date by the amount of time defined by mtimer_r
subroutine datefix(kdate_r,ktime_r,mtimer_r)

implicit none

include 'newmpar.h'
include 'parm.h'

integer leap
common/leap_yr/leap  ! 1 to allow leap years

integer, intent(inout) :: kdate_r,ktime_r,mtimer_r
integer, dimension(12) :: mdays
integer iyr,imo,iday,ihr,imins
integer mtimerh,mtimerm,mtimer
integer minsday,minsyr

data mdays/31,28,31,30,31,30,31,31,30,31,30,31/
data minsday/1440/,minsyr/525600/

if(kdate_r>=00600000.and.kdate_r<=00991231)then   ! old 1960-1999
  kdate_r=kdate_r+19000000
  write(6,*) 'For Y2K kdate_r altered to: ',kdate_r
endif
iyr=kdate_r/10000
imo=(kdate_r-10000*iyr)/100
iday=kdate_r-10000*iyr-100*imo
ihr=ktime_r/100
imins=ktime_r-100*ihr
write(6,*) 'entering datefix iyr,imo,iday,ihr,imins,mtimer_r: ', &
                             iyr,imo,iday,ihr,imins,mtimer_r

mdays(2)=28
if (leap==1) then
  if(mod(iyr,4)==0)mdays(2)=29
  if(mod(iyr,100)==0)mdays(2)=28
  if(mod(iyr,400)==0)mdays(2)=29
end if
do while (mtimer_r>minsday*mdays(imo))
  mtimer_r=mtimer_r-minsday*mdays(imo)
  imo=imo+1
  if(imo>12)then
    imo=1
    iyr=iyr+1
    if (leap==1) then
      mdays(2)=28      
      if(mod(iyr,4)==0)mdays(2)=29
      if(mod(iyr,100)==0)mdays(2)=28
      if(mod(iyr,400)==0)mdays(2)=29
    end if
  endif
enddo
if(diag)write(6,*)'b datefix iyr,imo,iday,ihr,imins,mtimer_r: ', &
                             iyr,imo,iday,ihr,imins,mtimer_r
do while (mtimer_r>minsday)
  mtimer_r=mtimer_r-minsday
  iday=iday+1
enddo
if(diag)write(6,*)'c datefix iyr,imo,iday,ihr,imins,mtimer_r: ', &
                            iyr,imo,iday,ihr,imins,mtimer_r

! at this point mtimer_r has been reduced to fraction of a day
mtimerh=mtimer_r/60
mtimerm=mtimer_r-mtimerh*60  ! minutes left over
ihr=ihr+mtimerh
imins=imins+mtimerm
if(imins==58.or.imins==59)then
  ! allow for roundoff for real timer from old runs
  write(6,*)'*** imins increased to 60 from imins = ',imins
  imins=60
endif
if(imins>59)then
  imins=imins-60
  ihr=ihr+1
endif
if(diag)write(6,*)'d datefix iyr,imo,iday,ihr,imins,mtimer_r: ', &
                             iyr,imo,iday,ihr,imins,mtimer_r
if(ihr>23)then
  ihr=ihr-24
  iday=iday+1
endif
if(diag)write(6,*)'e datefix iyr,imo,iday,ihr,imins,mtimer_r: ', &
                             iyr,imo,iday,ihr,imins,mtimer_r

if(iday>mdays(imo))then
  iday=iday-mdays(imo)
  imo=imo+1
  if(imo>12)then
    imo=imo-12
    iyr=iyr+1
  endif
endif

kdate_r=iday+100*(imo+100*iyr)
ktime_r=ihr*100+imins
mtimer=0
if(diag)write(6,*)'end datefix iyr,imo,iday,ihr,imins,mtimer_r: ', &
                               iyr,imo,iday,ihr,imins,mtimer_r
write(6,*)'leaving datefix kdate_r,ktime_r: ',kdate_r,ktime_r

return
end subroutine datefix

!--------------------------------------------------------------------
! Set up number of minutes from beginning of year
subroutine getzinp(fjd,jyear,jmonth,jday,jhour,jmin,mins)

use cc_mpi

implicit none

include 'dates.h'
include 'mpif.h'
include 'parm.h'

integer leap
common/leap_yr/leap  ! 1 to allow leap years

integer, intent(out) :: jyear,jmonth,jday,jhour,jmin ! start date of run
integer, intent(out) :: mins                         ! elapsed time from start of year
integer mstart,elp,ierr
integer, dimension(12) :: ndoy
integer, dimension(12), parameter :: odoy=(/0,31,59,90,120,151,181,212,243,273,304,334/)      ! days from beginning of year (1st Jan is 0)
real, intent(out) :: fjd

jyear =kdate/10000
jmonth=(kdate-jyear*10000)/100
jday  =kdate-jyear*10000-jmonth*100
jhour =ktime/100
jmin  =ktime-jhour*100
      
if (jmonth<1.or.jmonth>12) then
  write(6,*) "ERROR: Invalid month ",jmonth
  call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
end if 

ndoy=odoy
if (leap==1) then
  if (mod(jyear,4)  ==0) ndoy(3:12)=odoy(3:12)+1
  if (mod(jyear,100)==0) ndoy(3:12)=odoy(3:12)
  if (mod(jyear,400)==0) ndoy(3:12)=odoy(3:12)+1
end if

mstart=1440*(ndoy(jmonth)+jday-1) + 60*jhour + jmin ! mins from start of year
! mtimer contains number of minutes since the start of the run.
mins = mtimer + mstart

if(nhstest<0)then  ! aquaplanet test
  fjd = 79.+float(mod(mins,1440))/1440.  ! set to 21 March +frac of day
  mins=nint(fjd*1440.)
else
  fjd = float(mod(mins,(ndoy(12)+31)*1440))/1440.    ! 525600 = 1440*365
endif

return
end subroutine getzinp
      
end module infile
