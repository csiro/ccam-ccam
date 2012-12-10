module infile
      
! This module contains routines for reading netcdf files,
! vertical interpolation and some calendar functions.
! This module also contains the interface to the netcdf
! library.
      
! This version supports parallel localhist input files
! Multiple processors read as many input files as
! supplied in parallel and then send this data to the
! head processor (i.e., for interpolation and redistribution)
! The code can also identify restart files, in which case
! no additional message passing is required.
            
private
public vertint,datefix,getzinp,ncmsg
public histopen,histclose,histrd1,histrd4s
public attrib,histwrt3,histwrt4,freqwrite
public ccnf_open,ccnf_create,ccnf_close,ccnf_sync,ccnf_enddef,ccnf_redef,ccnf_nofill,ccnf_inq_varid,ccnf_inq_dimid
public ccnf_inq_dimlen,ccnf_def_dim,ccnf_def_dimu,ccnf_def_var,ccnf_def_var0,ccnf_get_var,ccnf_get_varg
public ccnf_get_vara,ccnf_get_var1,ccnf_get_att,ccnf_get_attg,ccnf_read,ccnf_put_var,ccnf_put_var1
public ccnf_put_vara,ccnf_put_att,ccnf_put_attg

interface ccnf_get_att
  module procedure ccnf_get_att_text, ccnf_get_att_real
end interface ccnf_get_att
interface ccnf_get_attg
  module procedure ccnf_get_att_intg1i, ccnf_get_att_intg2i
  module procedure ccnf_get_att_realg1r, ccnf_get_att_realg2r
end interface ccnf_get_attg
interface ccnf_get_vara
  module procedure ccnf_get_vara_real2r, ccnf_get_vara_real3r, ccnf_get_vara_real4r 
  module procedure ccnf_get_vara_int2i 
  module procedure ccnf_get_vara_double4d
end interface ccnf_get_vara
interface ccnf_get_var1
  module procedure ccnf_get_var1_real, ccnf_get_var1_int
end interface ccnf_get_var1
interface ccnf_get_var
  module procedure ccnf_get_var_real, ccnf_get_var_int
end interface ccnf_get_var
interface ccnf_put_att
  module procedure ccnf_put_att_text
end interface ccnf_put_att
interface ccnf_put_attg
  module procedure ccnf_put_att_realg1, ccnf_put_att_realg2
  module procedure ccnf_put_att_intg1, ccnf_put_att_intg2
  module procedure ccnf_put_att_textg
end interface ccnf_put_attg
interface ccnf_put_vara
  module procedure ccnf_put_vara_real2r, ccnf_put_vara_real3r
  module procedure ccnf_put_vara_int2i
  module procedure ccnf_put_vara_double2r
end interface ccnf_put_vara
interface ccnf_put_var
  module procedure ccnf_put_var_text2r
  module procedure ccnf_put_var_int2i, ccnf_put_var_int3i
end interface ccnf_put_var
interface ccnf_put_var1
  module procedure ccnf_put_var1_real, ccnf_put_var1_int, ccnf_put_var1_double
end interface ccnf_put_var1

integer, dimension(:), allocatable, save :: pnoff
integer, dimension(:,:), allocatable, save :: pioff,pjoff
integer, dimension(:), allocatable, save :: pncid
integer, dimension(:), allocatable, save :: pncidold
integer, save :: mynproc,fnproc
integer, save :: pil_g,pjl_g,pil,pjl,pnpan
integer, save :: comm_ip,comm_ipold
logical, save :: ptest

integer(kind=2), parameter :: minv = -32500
integer(kind=2), parameter :: maxv = 32500
integer(kind=2), parameter :: missval = -32501
      
contains

!--------------------------------------------------------------
! Interface for reading 2D+time fields
subroutine histrd1(ncid,iarchi,ier,name,ik,jk,var,ifull)
      
use cc_mpi
      
implicit none
      
include 'parm.h'
      
integer ncid,iarchi,ier,ik,jk,ifull
character(len=*) name
real, dimension(:), intent(inout) :: var ! may be dummy argument from myid/=0

if (ifull/=ik*jk.and.ptest) then
  ! read local arrays without gather and distribute

#ifdef debug       
  if (size(var)/=ifull) then
    write(6,*) "ERROR: Incorrect use of dummy var in histrd1"
    call ccmpi_abort(-1)
  end if
#endif

  call hr1p(iarchi,ier,name,.true.,var)

  if(ier==0.and.mod(ktau,nmaxpr)==0.and.myid==0)then
    write(6,'("done histrd1 ",a8,i4,i3)') name,ier,iarchi
  end if

else        
      
  ! split up processors to save memory.  No need to allocate global
  ! arrays on myid/=0.
  if (myid==0) then
#ifdef debug
    if (size(var)/=ifull) then
      write(6,*) "ERROR: Incorrect use of dummy var in histrd1"
      call ccmpi_abort(-1)
    end if
#endif

    call hr1a(ncid,iarchi,ier,name,ik,jk,var,ifull)
  else
    if(ifull/=ik*jk)then
      ! read local arrays with gather and distribute

#ifdef debug
      if (size(var)/=ifull) then
        write(6,*) "ERROR: Incorrect use of dummy var in histrd"
        call ccmpi_abort(-1)
      end if
#endif

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
real, dimension(:), intent(inout) :: var
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
  var(:)=globvar(:)
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
#ifndef usenc3
use netcdf
#endif
      
implicit none

include 'newmpar.h'
#ifdef usenc3
include 'netcdf.inc'
#endif

integer, intent(in) :: iarchi
integer, intent(out) :: ier
integer, dimension(3) :: start,count
integer ipf,jpmax,iptst2,lcomm
integer ler,idv
real, dimension(:), intent(inout), optional :: var
real, dimension(pil*pjl*pnpan) :: rvar
real addoff,sf
logical, intent(in) :: qtest
character(len=*), intent(in) :: name

#ifdef debug
if (qtest) then
  if (mynproc>1) then
    write(6,*) "ERROR: Invalid use of qtest"
    write(6,*) "Expecting 1 file, but required to load ",mynproc
    call ccmpi_abort(-1)
  end if
  if (.not.present(var)) then
    write(6,*) "ERROR: Missing var in hr1p with qtest"
    call ccmpi_abort(-1)
  end if
  if (size(var)/=pil*pjl*pnpan) then
    write(6,*) "ERROR: Invalid var size in hr1p with qtest"
    write(6,*) "size(var),pil,pjl,pnpan ",size(var),pil,pjl,pnpan
    call ccmpi_abort(-1)
  end if
else
  if (myid==0) then
    if (.not.present(var)) then
      write(6,*) "ERROR: Missing var in hr1p"
      call ccmpi_abort(-1)
    end if
    if (size(var)/=pil_g*pjl_g) then
      write(6,*) "ERROR: Invalid var size in hr1p"
      write(6,*) "size(var),pil_g,pjl_g ",size(var),pil_g,pjl_g
      call ccmpi_abort(-1)
    end if
  else
    if (present(var)) then
      write(6,*) "ERROR: Invalid use of var in hr1p"
      call ccmpi_abort(-1)
    end if
  end if
end if
#endif

start = (/ 1, 1, iarchi /)
count = (/ pil, pjl*pnpan, 1 /)
ier=0

iptst2=mod(fnproc,nproc)
      
do ipf=0,mynproc-1
  rvar=0.

  jpmax=nproc
  if (ipf==mynproc-1.and.myid<iptst2) then
    jpmax=iptst2
  end if
  
  ! get variable idv
#ifdef usenc3
  ier=nf_inq_varid(pncid(ipf),name,idv)
  if(ier/=nf_noerr)then
#else
  ier=nf90_inq_varid(pncid(ipf),name,idv)
  if(ier/=nf90_noerr)then
#endif
    if (myid==0.and.ipf==0) then
      write(6,*) '***absent field for ncid,name,idv,ier: ',pncid(0),name,idv,ier
    end if
  else
    ! obtain scaling factors and offsets from attributes
#ifdef usenc3
    ler=nf_get_att_real(pncid(ipf),idv,'add_offset',addoff)
    if (ler/=nf_noerr) addoff=0.
    ler=nf_get_att_real(pncid(ipf),idv,'scale_factor',sf)
    if (ler/=nf_noerr) sf=1.
    ier=nf_get_vara_real(pncid(ipf),idv,start,count,rvar)
    call ncmsg(name,ier)
#else
    ler=nf90_get_att(pncid(ipf),idv,'add_offset',addoff)
    if (ler/=nf90_noerr) addoff=0.
    ler=nf90_get_att(pncid(ipf),idv,'scale_factor',sf)
    if (ler/=nf90_noerr) sf=1.
    ier=nf90_get_var(pncid(ipf),idv,rvar,start=start,count=count)
    call ncmsg(name,ier)
#endif
    ! unpack compressed data
    rvar=rvar*sf+addoff
  end if ! ier
      
  if (qtest) then
    ! e.g., restart file
    var=rvar
  else
    ! e.g., mesonest file
    
    lcomm=comm_world
    if (jpmax<nproc) lcomm=comm_ip
    
    if (myid==0) then
      call host_hr1p(lcomm,jpmax,ipf,rvar,var)
    else
      call proc_hr1p(lcomm,rvar)
    end if
  
  end if ! qtest

end do ! ipf

return
end subroutine hr1p

subroutine host_hr1p(lcomm,jpmax,ipf,rvar,var)

use cc_mpi

implicit none

include 'newmpar.h'

integer, intent(in) :: lcomm,jpmax,ipf
integer jpf,ip,n,no,ca,cb,j,iq,iqi
real, dimension(:), intent(out) :: var
real, dimension(pil*pjl*pnpan), intent(in) :: rvar
real, dimension(pil*pjl*pnpan*nproc) :: gvar 

call ccmpi_gatherx(gvar,rvar,0,lcomm)
do jpf=0,jpmax-1
  ip=ipf*nproc+jpf
  do n=0,pnpan-1
    no=n-pnoff(ip)+1
    ca=pioff(ip,no)
    cb=pjoff(ip,no)+no*pil_g
    do j=1,pjl
      iq=ca+(j+cb-1)*pil_g
      iqi=(j-1)*pil+n*pil*pjl+jpf*pil*pjl*pnpan
      var(iq+1:iq+pil)=gvar(iqi+1:iqi+pil)
    end do
  end do
end do

return
end subroutine host_hr1p

subroutine proc_hr1p(lcomm,rvar)

use cc_mpi

implicit none

integer, intent(in) :: lcomm
real, dimension(pil*pjl*pnpan), intent(in) :: rvar
real, dimension(0) :: gvar 

call ccmpi_gatherx(gvar,rvar,0,lcomm)

return
end subroutine proc_hr1p

!--------------------------------------------------------------   
! Interface for reading 3D+time fields
subroutine histrd4s(ncid,iarchi,ier,name,ik,jk,kk,var,ifull)
      
use cc_mpi
      
implicit none
      
include 'parm.h'
      
integer, intent(in) :: ncid,iarchi,ik,jk,kk,ifull
integer, intent(out) :: ier
character(len=*), intent(in) :: name
real, dimension(:), intent(inout) :: var ! may be dummy argument from myid/=0
      
if (ifull/=ik*jk.and.ptest) then
  ! read local arrays without gather and distribute

#ifdef debug       
  if (size(var)/=ifull) then
    write(6,*) "ERROR: Incorrect use of dummy var in histrd4s"
    call ccmpi_abort(-1)
  end if
#endif

  call hr4p(iarchi,ier,name,kk,.true.,var)

  if(ier==0.and.mod(ktau,nmaxpr)==0.and.myid==0)then
    write(6,'("done histrd4 ",a8,i3,i4,i3)') name,kk,ier,iarchi
  endif

else        
      
  ! split up processors to save memory.  No need to allocate global
  ! arrays on myid/=0.
  if (myid==0) then
#ifdef debug
    if (size(var)/=ifull) then
      write(6,*) "ERROR: Incorrect use of dummy var in histrd4s"
      call ccmpi_abort(-1)
    end if
#endif
 
    call hr4sa(ncid,iarchi,ier,name,ik,jk,kk,var,ifull)
  else
    if(ifull/=ik*jk)then
      ! read local arrays with gather and distribute

#ifdef debug
      if (size(var)/=ifull) then
        write(6,*) "ERROR: Incorrect use of dummy var in histrd"
        call ccmpi_abort(-1)
      end if
#endif
 
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
real, dimension(:) :: var
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
  var(:)=globvar(:)
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
#ifndef usenc3
use netcdf
#endif
      
implicit none

include 'newmpar.h'
#ifdef usenc3
include 'netcdf.inc'
#endif

integer, intent(in) :: iarchi,kk
integer, intent(out) :: ier
integer, dimension(4) :: start,ncount
integer ipf,jpmax,iptst2,lcomm
integer ler,idv
real, dimension(:), intent(inout), optional :: var
real, dimension(pil*pjl*pnpan) :: rvar
real, dimension(pil*pjl*pnpan*nproc) :: gvar
real addoff,sf
logical, intent(in) :: qtest
character(len=*), intent(in) :: name

#ifdef debug
if (qtest) then
  if (mynproc>1) then
     write(6,*) "ERROR: Invalid use of qtest"
     call ccmpi_abort(-1)
  end if
  if (.not.present(var)) then
    write(6,*) "ERROR: Missing var in hr4p with qtest"
    call ccmpi_abort(-1)
  end if
  if (size(var)/=pil*pjl*pnpan) then
    write(6,*) "ERROR: Invalid var size in hr4p with qtest"
    write(6,*) "size(var),pil,pjl,pnpan ",size(var),pil,pjl,pnpan
    call ccmpi_abort(-1)
  end if
else
  if (myid==0) then
    if (.not.present(var)) then
      write(6,*) "ERROR: Missing var in hr4p"
      call ccmpi_abort(-1)
    end if
    if (size(var)/=pil_g*pjl_g) then
      write(6,*) "ERROR: Invalid var size in hr4p"
      write(6,*) "size(var),pil_g,pjl_g ",size(var),pil_g,pjl_g
      call ccmpi_abort(-1)
    end if
  else
    if (present(var)) then
      write(6,*) "ERROR: Invalid use of var in hr4p"
      call ccmpi_abort(-1)
    end if
  end if
end if
#endif

start = (/ 1, 1, kk, iarchi /)
ncount = (/ pil, pjl*pnpan, 1, 1 /)
ier=0

iptst2=mod(fnproc,nproc)
      
do ipf=0,mynproc-1
  rvar=0.

  jpmax=nproc
  if (ipf==mynproc-1.and.myid<iptst2) then
    jpmax=iptst2
  end if

  ! get variable idv
#ifdef usenc3
  ier=nf_inq_varid(pncid(ipf),name,idv)
  if(ier/=nf_noerr)then
#else
  ier=nf90_inq_varid(pncid(ipf),name,idv)
  if(ier/=nf90_noerr)then
#endif
    if (myid==0.and.ipf==0.and.kk==1) then
      write(6,*) '***absent field for ncid,name,idv,ier: ',pncid(0),name,idv,ier
    end if
  else
    ! obtain scaling factors and offsets from attributes
#ifdef usenc3
    ler=nf_get_att_real(pncid(ipf),idv,'add_offset',addoff)
    if (ler/=nf_noerr) addoff=0.
    ler=nf_get_att_real(pncid(ipf),idv,'scale_factor',sf)
    if (ler/=nf_noerr) sf=1.
    ier=nf_get_vara_real(pncid(ipf),idv,start,ncount,rvar)
    call ncmsg(name,ier)
#else
    ler=nf90_get_att(pncid(ipf),idv,'add_offset',addoff)
    if (ler/=nf90_noerr) addoff=0.
    ler=nf90_get_att(pncid(ipf),idv,'scale_factor',sf)
    if (ler/=nf90_noerr) sf=1.
    ier=nf90_get_var(pncid(ipf),idv,rvar,start=start,count=ncount)
    call ncmsg(name,ier)
#endif
    ! unpack data
    rvar=rvar*sf+addoff
  end if ! ier

  if (qtest) then
    ! expected restart file
    var=rvar
  else
    ! expected mesonest file

    lcomm=comm_world
    if (jpmax<nproc) lcomm=comm_ip

    if (myid==0) then
      call host_hr1p(lcomm,jpmax,ipf,rvar,var)
    else
      call proc_hr1p(lcomm,rvar)
    end if

  end if ! qtest

end do ! ipf

return
end subroutine hr4p

!--------------------------------------------------------------
! Trap netcdf error messages
subroutine ncmsg(txt,ierr)

use cc_mpi
#ifndef usenc3
use netcdf
#endif

implicit none

#ifdef usenc3
include 'netcdf.inc'
#endif

integer, intent(in) :: ierr
character(len=*), intent(in) :: txt

#ifdef usenc3
if (ierr/=nf_noerr) then
  write(6,*) txt," ",nf_strerror(ierr)
#else
if (ierr/=nf90_noerr) then
  write(6,*) txt," ",nf90_strerror(ierr)
#endif
  call ccmpi_abort(-1)
end if

return
end subroutine ncmsg

!--------------------------------------------------------------
! This subroutine opens parallel input files
subroutine histopen(ncid,ifile,ier)
      
use cc_mpi
#ifndef usenc3
use netcdf
#endif
      
implicit none
      
include 'newmpar.h'
#ifdef usenc3
include 'netcdf.inc'
#endif
      
integer, parameter :: nihead   = 54
      
integer, dimension(nihead) :: ahead
integer, dimension(0:5) :: duma,dumb
integer, dimension(6) :: idum
integer, intent(out) :: ncid,ier
integer ler,dumr
integer resid,is,ipf,dmode
integer ipin,nxpr,nypr
integer ltst
character(len=*), intent(in) :: ifile
character(len=170) :: pfile
character(len=7) :: fdecomp
logical omode

if (myid==0) then
  ! attempt to open single file with myid==0
#ifdef usenc3
  ier=nf_open(ifile,nf_nowrite,ncid)
#else
  ier=nf90_open(ifile,nf90_nowrite,ncid)
#endif
  fnproc=1
  dmode=0
  pil=0
  pjl=0
  pnpan=0
  ptest=.false.
      
  ! attempt to open parallel files
#ifdef usenc3
  if (ier/=nf_noerr) then
    write(pfile,"(a,'.',i6.6)") trim(ifile), 0
    ier=nf_open(pfile,nf_nowrite,ncid)
    if (ier/=nf_noerr) then
#else
  if (ier/=nf90_noerr) then
    write(pfile,"(a,'.',i6.6)") trim(ifile), 0
    ier=nf90_open(pfile,nf90_nowrite,ncid)
    if (ier/=nf90_noerr) then
#endif
      write(6,*) "WARN: Cannot open ",pfile
    else  
      write(6,*) "Found parallel input file ",ifile
      fdecomp=''
#ifdef usenc3
      ler=nf_get_att_int(ncid,nf_global,"nproc",fnproc)
      call ncmsg("nproc",ler)
      ler=nf_get_att_text(ncid,nf_global,"decomp",fdecomp)
      call ncmsg("decomp",ler)
#else
      ler=nf90_get_att(ncid,nf90_global,"nproc",fnproc)
      call ncmsg("nproc",ler)
      ler=nf90_get_att(ncid,nf90_global,"decomp",fdecomp)
      call ncmsg("decomp",ler)
#endif
      select case(fdecomp)
        case('face')
          dmode=1
        case('uniform')
          dmode=2
        case default
          write(6,*) "ERROR: Unknown decomposition ",fdecomp
          call ccmpi_abort(-1)
      end select
    end if
  else
    ! nproc should only exist in multi-file input
#ifdef usenc3
    ler=nf_get_att_int(ncid,nf_global,"nproc",dumr)
    if (ler==nf_noerr) then
#else
    ler=nf90_get_att(ncid,nf90_global,"nproc",dumr)
    if (ler==nf90_noerr) then
#endif
      write(6,*) "ERROR: Incorrect base filename"
      call ccmpi_abort(-1)
    end if
    write(6,*) "Found single input file ",ifile
  end if

#ifdef usenc3
  if (ier==nf_noerr) then
    ler=nf_get_att_int(ncid,nf_global,"int_header",ahead)
    call ncmsg("int_header",ler)
#else
  if (ier==nf90_noerr) then
    ler=nf90_get_att(ncid,nf90_global,"int_header",ahead)
    call ncmsg("int_header",ler)
#endif
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
    
  end if

  idum(1)=fnproc
  idum(2)=pil
  idum(3)=pjl
  idum(4)=pnpan
  idum(5)=0
  if (ptest) idum(5)=1
  idum(6)=ier  
end if

! Broadcast file metadata
call ccmpi_bcast(idum(1:5),0,comm_world)
fnproc=idum(1)
pil=idum(2)
pjl=idum(3)
pnpan=idum(4)
ptest=(idum(5)==1)
ier=idum(6)

#ifdef usenc3
if (ier/=nf_noerr) return
#else
if (ier/=nf90_noerr) return
#endif

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
  write(pfile,"(a,'.',i6.6)") trim(ifile), ipin
#ifdef usenc3
  ler=nf_open(pfile,nf_nowrite,pncid(ipf))
  if (ler/=nf_noerr) then
    write(6,*) "ERROR: Cannot open ",pfile
    write(6,*) nf_strerror(ler)
#else
  ler=nf90_open(pfile,nf90_nowrite,pncid(ipf))
  if (ler/=nf90_noerr) then
    write(6,*) "ERROR: Cannot open ",pfile
    write(6,*) nf90_strerror(ler)
#endif
    call ccmpi_abort(-1)
  end if
end do

ltst=0
if (myid<resid) ltst=1
call ccmpi_commsplit(comm_ip,comm_world,ltst,myid)
    
return
end subroutine histopen
      
!--------------------------------------------------------------
! This subroutine closes parallel input files
subroutine histclose

use cc_mpi
#ifndef usenc3
use netcdf
#endif
      
implicit none

#ifdef usenc3
include 'netcdf.inc'
#endif
     
integer ipf,ipin,plen
integer ierr
      
if (allocated(pncidold)) then
  plen=size(pncidold)
  do ipf=0,plen-1
#ifdef usenc3
    ierr=nf_close(pncidold(ipf))
#else
    ierr=nf90_close(pncidold(ipf))
#endif
  end do
  call ccmpi_commfree(comm_ipold)
  deallocate(pncidold)
end if
      
if (allocated(pncid)) then
  plen=size(pncid)
  allocate(pncidold(0:plen-1))
  pncidold=pncid
  comm_ipold=comm_ip
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
real, dimension(ifull,kk) :: told
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
  call ccmpi_abort(-1)
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

!--------------------------------------------------------------------
! DEFINE ATTRIBUTES
subroutine attrib(cdfid,dim,ndim,name,lname,units,xmin,xmax,daily,itype)

use cc_mpi
#ifndef usenc3
use netcdf
#endif

implicit none

#ifdef usenc3
include 'netcdf.inc'
#endif

integer, intent(in) :: cdfid, itype, ndim
integer, intent(in) :: daily
integer, dimension(ndim), intent(in) :: dim
integer vtype, ier, idv, lsize
real, intent(in) :: xmin, xmax
real scalef, addoff
character(len=*), intent(in) :: name
character(len=*), intent(in) :: lname
character(len=*), intent(in) :: units

#ifdef usenc3
#ifdef r8i8
  vtype = nf_double
#else
  vtype = nf_float
#endif
if (itype==1) vtype = nf_short
ier = nf_def_var(cdfid, name, vtype, ndim, dim, idv)
call ncmsg("def_var",ier)
lsize=len_trim(lname)
ier = nf_put_att_text(cdfid,idv,'long_name',lsize,lname)
call ncmsg("long_name",ier)
lsize=len_trim(units)
if(lsize>0)then
  ier = nf_put_att_text(cdfid,idv,'units',lsize,units)
  call ncmsg("units",ier)
endif
if (vtype == nf_short) then
  ier = nf_put_att_int2(cdfid,idv,'valid_min',nf_short,1,minv)
  call ncmsg("valid_min",ier)
  ier = nf_put_att_int2(cdfid,idv,'valid_max',nf_short,1,maxv)
  call ncmsg("valid_max",ier)
  ier = nf_put_att_int2(cdfid,idv,'missing_value',nf_short,1,missval)
  call ncmsg("missing_value",ier)
  scalef=(xmax-xmin)/(real(maxv)-real(minv))
  addoff=xmin-scalef*minv
  ier = nf_put_att_real(cdfid,idv,'add_offset',nf_float,1,addoff)
  call ncmsg("add_offset",ier)
  ier = nf_put_att_real(cdfid,idv,'scale_factor',nf_float,1,scalef)
  call ncmsg("scale_factor",ier)
else
  ier = nf_put_att_real(cdfid,idv,'missing_value',nf_float,1,nf_fill_float)
  call ncmsg("missing_value",ier)
endif
ier = nf_put_att_text(cdfid,idv,'FORTRAN_format',5,'G11.4')
call ncmsg("FORTRAN_format",ier)
if(daily>0)then
  ier = nf_put_att_text(cdfid,idv,'valid_time',5,'daily')
  call ncmsg("valid_time",ier)
endif
#else
#ifdef r8i8
  vtype = nf90_double
#else
  vtype = nf90_float
#endif
if (itype==1) vtype = nf90_short
ier = nf90_def_var(cdfid, name, vtype, dim, idv, deflate_level=1)
call ncmsg("def_var",ier)
lsize=len_trim(lname)
ier = nf90_put_att(cdfid,idv,'long_name',lname)
call ncmsg("long_name",ier)
lsize=len_trim(units)
if(lsize>0)then
  ier = nf90_put_att(cdfid,idv,'units',units)
  call ncmsg("units",ier)
endif
if (vtype == nf90_short) then
  ier = nf90_put_att(cdfid,idv,'valid_min',minv)
  call ncmsg("valid_min",ier)
  ier = nf90_put_att(cdfid,idv,'valid_max',maxv)
  call ncmsg("valid_max",ier)
  ier = nf90_put_att(cdfid,idv,'missing_value',missval)
  call ncmsg("missing_value",ier)
  scalef=(xmax-xmin)/(real(maxv)-real(minv))
  addoff=xmin-scalef*minv
  ier = nf90_put_att(cdfid,idv,'add_offset',addoff)
  call ncmsg("add_offset",ier)
  ier = nf90_put_att(cdfid,idv,'scale_factor',scalef)
  call ncmsg("scale_factor",ier)
else
  ier = nf90_put_att(cdfid,idv,'missing_value',nf90_fill_float)
  call ncmsg("missing_value",ier)
endif
ier = nf90_put_att(cdfid,idv,'FORTRAN_format','G11.4')
call ncmsg("FORTRAN_format",ier)
if(daily>0)then
  ier = nf90_put_att(cdfid,idv,'valid_time','daily')
  call ncmsg("valid_time",ier)
endif
#endif
      
return
end subroutine attrib

!--------------------------------------------------------------------
! 3D NETCDF WRITE ARRAY ROUTINES
subroutine histwrt3(var,sname,idnc,iarch,local,lwrite)

use cc_mpi              ! CC MPI routines
#ifndef usenc3
use netcdf              ! Netcdf parameters
#endif

implicit none

include 'newmpar.h'     ! Grid parameters
#ifdef usenc3
include 'netcdf.inc'    ! Netcdf parameters
#endif

integer, intent(in) :: idnc, iarch
real, dimension(ifull), intent(in) :: var
real, dimension(ifull,1) :: wvar
character(len=*), intent(in) :: sname
logical, intent(in) :: local, lwrite

wvar(:,1)=var
if (.not.lwrite) then
#ifdef usenc3
  wvar=real(nf_fill_float)
#else
  wvar=real(nf90_fill_float)
#endif
end if

if (local) then
  call fw3l(wvar,sname,idnc,iarch,1)
else if (myid==0) then
  call fw3a(wvar,sname,idnc,iarch,1)
else
  call ccmpi_gather(wvar)
end if

return
end subroutine histwrt3

subroutine freqwrite(fncid,cname,fiarch,istep,local,datain)

use cc_mpi               ! CC MPI routines
      
implicit none
      
include 'newmpar.h'      ! Grid parameters
      
integer, intent(in) :: fncid,fiarch,istep
real, dimension(ifull,istep), intent(in) :: datain
logical, intent(in) :: local
character(len=*), intent(in) :: cname
      
if (local) then
  call fw3l(datain,cname,fncid,fiarch,istep)
elseif (myid==0) then
  call fw3a(datain,cname,fncid,fiarch,istep)
else
  call ccmpi_gather(datain)
endif
     
return
end subroutine freqwrite

subroutine fw3l(var,sname,idnc,iarch,istep)

use cc_mpi               ! CC MPI routines
#ifndef usenc3
use netcdf               ! Netcdf parameters
#endif
      
implicit none
      
include 'newmpar.h'      ! Grid parameters
#ifdef usenc3
include 'netcdf.inc'     ! Netcdf parameters
#endif
include 'parm.h'         ! Model configuration
      
integer, intent(in) :: idnc,iarch,istep
integer mid, vtype, ier
integer, dimension(3) :: start, ncount
integer(kind=2), dimension(ifull,istep) :: ipack
real, dimension(ifull,istep), intent(in) :: var
real addoff, scale_f
character(len=*), intent(in) :: sname

start = (/ 1, 1, iarch /)
ncount = (/ il, jl, istep /)

#ifdef usenc3
ier=nf_inq_varid(idnc,sname,mid)
call ncmsg(sname,ier)
ier=nf_inq_vartype(idnc,mid,vtype)
if(vtype == nf_short)then
#else
ier=nf90_inq_varid(idnc,sname,mid)
call ncmsg(sname,ier)
ier=nf90_inquire_variable(idnc,mid,xtype=vtype)
if(vtype == nf90_short)then
#endif
  if (all(var>9.8E36)) then
    ipack=missval
  else
#ifdef usenc3  
    ier=nf_get_att_real(idnc,mid,'add_offset',addoff)
    ier=nf_get_att_real(idnc,mid,'scale_factor',scale_f)
#else
    ier=nf90_get_att(idnc,mid,'add_offset',addoff)
    ier=nf90_get_att(idnc,mid,'scale_factor',scale_f)
#endif
    ipack=max(min(nint((var-addoff)/scale_f),maxv),minv)
  end if
#ifdef usenc3
  ier=nf_put_vara_int2(idnc,mid,start,ncount,ipack)
else
#ifdef r8i8
  ier=nf_put_vara_double(idnc,mid,start,ncount,var)
#else
  ier=nf_put_vara_real(idnc,mid,start,ncount,var)
#endif
#else
  ier=nf90_put_var(idnc,mid,ipack,start=start,count=ncount)
else
ier=nf90_put_var(idnc,mid,var,start=start,count=ncount)
#endif
end if
call ncmsg(sname,ier)

if (mod(ktau,nmaxpr)==0.and.myid==0) then
#ifdef usenc3
  if (any(var==nf_fill_float)) then
#else
  if (any(var==nf90_fill_float)) then
#endif
    write(6,'("histwrt3 ",a7,i4,a7)') sname,iarch,"missing"
  else
    write(6,'("histwrt3 ",a7,i4)') sname,iarch
  end if
end if

return
end subroutine fw3l
      
subroutine fw3a(var,sname,idnc,iarch,istep)

use cc_mpi               ! CC MPI routines
#ifndef usenc3
use netcdf               ! Netcdf parameters
#endif

implicit none

include 'newmpar.h'      ! Grid parameters
#ifdef usenc3
include 'netcdf.inc'     ! Netcdf parameters
#endif
include 'parm.h'         ! Model configuration

integer, intent(in) :: idnc, iarch, istep
integer mid, vtype, ier
integer imn, imx, jmn, jmx, iq
integer, dimension(3) :: start, ncount
integer(kind=2), dimension(ifull_g,istep) :: ipack
real, dimension(ifull,istep), intent(in) :: var
real, dimension(ifull_g,istep) :: globvar
real addoff, scale_f
real varn, varx
character(len=*), intent(in) :: sname
      
call ccmpi_gather(var, globvar)

start = (/ 1, 1, iarch /)
ncount = (/ il_g, jl_g, istep /)

!     find variable index
#ifdef usenc3
ier=nf_inq_varid(idnc,sname,mid)
call ncmsg(sname,ier)
ier=nf_inq_vartype(idnc,mid,vtype)
if (vtype == nf_short) then
#else
ier=nf90_inq_varid(idnc,sname,mid)
call ncmsg(sname,ier)
ier=nf90_inquire_variable(idnc,mid,xtype=vtype)
if (vtype == nf90_short) then
#endif
  if (all(globvar>9.8e36)) then
    ipack=missval
  else
#ifdef usenc3
    ier=nf_get_att_real(idnc,mid,'add_offset',addoff)
    ier=nf_get_att_real(idnc,mid,'scale_factor',scale_f)
#else
    ier=nf90_get_att(idnc,mid,'add_offset',addoff)
    ier=nf90_get_att(idnc,mid,'scale_factor',scale_f)
#endif
    ipack=max(min(nint((globvar-addoff)/scale_f),maxv),minv)
  endif
#ifdef usenc3
  ier=nf_put_vara_int2(idnc,mid,start,ncount,ipack)
else
#ifdef r8i8
  ier=nf_put_vara_double(idnc,mid,start,ncount,globvar)
#else
  ier=nf_put_vara_real(idnc,mid,start,ncount,globvar)
#endif
#else
  ier=nf90_put_var(idnc,mid,ipack,start=start,count=ncount)
else
  ier=nf90_put_var(idnc,mid,globvar,start=start,count=ncount)
#endif
endif
call ncmsg(sname,ier)

if(mod(ktau,nmaxpr)==0)then
#ifdef usenc3
  if (any(globvar==nf_fill_float)) then
#else
  if (any(globvar==nf90_fill_float)) then
#endif
    write(6,'("histwrt3 ",a7,i4,a7)') sname,iarch,"missing"
  else
    varn = minval(globvar(:,1))
    varx = maxval(globvar(:,1))
    ! This should work ???? but sum trick is more portable???
    ! iq = minloc(globvar,dim=1)
    iq = sum(minloc(globvar(:,1)))
    ! Convert this 1D index to 2D
    imn = 1 + modulo(iq-1,il_g)
    jmn = 1 + (iq-1)/il_g
    iq = sum(maxloc(globvar(:,1)))
    ! Convert this 1D index to 2D
    imx = 1 + modulo(iq-1,il_g)
    jmx = 1 + (iq-1)/il_g
    write(6,'("histwrt3 ",a7,i4,f12.4,2i4,f12.4,2i4,f12.4)') &
                   sname,iarch,varn,imn,jmn,varx,imx,jmx,    &
                   globvar(id+(jd-1)*il_g,1)
  end if
endif

return
end subroutine fw3a

!--------------------------------------------------------------------
! 4D NETCDF WRITE ARRAY ROUTINES
subroutine histwrt4(var,sname,idnc,iarch,local,lwrite)

use cc_mpi              ! CC MPI routines
#ifndef usenc3
use netcdf              ! Netcdf parameters
#endif

implicit none

include 'newmpar.h'     ! Grid parameters
#ifdef usenc3
include 'netcdf.inc'    ! Netcdf parameters
#endif

integer, intent(in) :: idnc, iarch
real, dimension(ifull,kl), intent(in) :: var
real, dimension(ifull,kl) :: wvar
character(len=*), intent(in) :: sname
logical, intent(in) :: local,lwrite

wvar=var
if (.not.lwrite) then
#ifdef usenc3
  wvar=real(nf_fill_float)
#else
  wvar=real(nf90_fill_float)
#endif
endif

if (local) then
  call hw4l(wvar,sname,idnc,iarch)
elseif (myid==0) then
  call hw4a(wvar,sname,idnc,iarch)
else
  call ccmpi_gather(wvar)
endif

return
end subroutine histwrt4
      
subroutine hw4l(var,sname,idnc,iarch)

use cc_mpi               ! CC MPI routines
#ifndef usenc3
use netcdf               ! Netcdf parameters
#endif

implicit none

include 'newmpar.h'      ! Grid parameters
#ifdef usenc3
include 'netcdf.inc'     ! Netcdf parameters
#endif
include 'parm.h'         ! Model configuration

integer, intent(in) :: idnc, iarch
integer mid, vtype, ier
integer, dimension(4) :: start, ncount
integer(kind=2), dimension(ifull,kl) :: ipack
real, dimension(ifull,kl), intent(in) :: var
real addoff, scale_f
character(len=*), intent(in) :: sname

start = (/ 1, 1, 1, iarch /)
ncount = (/ il, jl, kl, 1 /)

#ifdef usenc3
ier=nf_inq_varid(idnc,sname,mid)
call ncmsg(sname,ier)
ier = nf_inq_vartype(idnc, mid, vtype)
if(vtype == nf_short)then
#else
ier=nf90_inq_varid(idnc,sname,mid)
call ncmsg(sname,ier)
ier = nf90_inquire_variable(idnc, mid, xtype=vtype)
if(vtype == nf90_short)then
#endif
  if(all(var>9.8e36))then
    ipack=missval
  else
#ifdef usenc3
    ier=nf_get_att_real(idnc,mid,'add_offset',addoff)
    ier=nf_get_att_real(idnc,mid,'scale_factor',scale_f)
#else
    ier=nf90_get_att(idnc,mid,'add_offset',addoff)
    ier=nf90_get_att(idnc,mid,'scale_factor',scale_f)
#endif
    ipack=max(min(nint((var-addoff)/scale_f),maxv),minv)
  endif
#ifdef usenc3
  ier=nf_put_vara_int2(idnc,mid,start,ncount,ipack)
else
#ifdef r8i8
  ier=nf_put_vara_double(idnc,mid,start,ncount,var)
#else
  ier=nf_put_vara_real(idnc,mid,start,ncount,var)
#endif
#else
  ier=nf90_put_var(idnc,mid,ipack,start=start,count=ncount)
else
  ier=nf90_put_var(idnc,mid,var,start=start,count=ncount)
#endif
endif
call ncmsg(sname,ier)

if(mod(ktau,nmaxpr)==0.and.myid==0)then
#ifdef usenc3
  if (any(var==nf_fill_float)) then
#else
  if (any(var==nf90_fill_float)) then
#endif
    write(6,'("histwrt4 ",a7,i4,a7)') sname,iarch,"missing"
  else
    write(6,'("histwrt4 ",a7,i4)') sname,iarch
  end if
endif

return
end subroutine hw4l      

subroutine hw4a(var,sname,idnc,iarch)

use cc_mpi              ! CC MPI routines
#ifndef usenc3
use netcdf              ! Netcdf parameters
#endif

implicit none

include 'newmpar.h'     ! Grid parameters
#ifdef usenc3
include 'netcdf.inc'    ! Netcdf parameters
#endif
include 'parm.h'        ! Model configuration

integer, intent(in) :: idnc, iarch
integer mid, vtype, ier
integer imx, jmx, kmx, iq
integer, dimension(4) :: start, ncount
integer, dimension(2) :: max_result
integer(kind=2), dimension(ifull_g,kl) :: ipack
real addoff, scale_f
real varn, varx
real, dimension(ifull,kl), intent(in) :: var
real, dimension(ifull_g,kl) :: globvar
character(len=*), intent(in) :: sname
      
call ccmpi_gather(var, globvar)
start(1) = 1
start(2) = 1
start(3) = 1
start(4) = iarch
ncount(1) = il_g
ncount(2) = jl_g
ncount(3) = kl
ncount(4) = 1

!     find variable index
#ifdef usenc3
ier=nf_inq_varid(idnc,sname,mid)
call ncmsg(sname,ier)
ier = nf_inq_vartype(idnc, mid, vtype)
if(vtype == nf_short)then
#else
ier=nf90_inq_varid(idnc,sname,mid)
call ncmsg(sname,ier)
ier = nf90_inquire_variable(idnc, mid, xtype=vtype)
if(vtype == nf90_short)then
#endif
  if(all(globvar>9.8e36))then
    ipack=missval
  else
#ifdef usenc3
    ier=nf_get_att_real(idnc,mid,'add_offset',addoff)
    ier=nf_get_att_real(idnc,mid,'scale_factor',scale_f)
#else
    ier=nf90_get_att(idnc,mid,'add_offset',addoff)
    ier=nf90_get_att(idnc,mid,'scale_factor',scale_f)
#endif
    ipack=max(min(nint((globvar-addoff)/scale_f),maxv),minv)
  endif
#ifdef usenc3
  ier=nf_put_vara_int2(idnc,mid,start,ncount,ipack)
else
#ifdef r8i8
  ier=nf_put_vara_double(idnc,mid,start,ncount,globvar)
#else
  ier=nf_put_vara_real(idnc,mid,start,ncount,globvar)
#endif
#else
  ier=nf90_put_var(idnc,mid,ipack,start=start,count=ncount)
else
  ier=nf90_put_var(idnc,mid,globvar,start=start,count=ncount)
#endif
endif
call ncmsg(sname,ier)

if(mod(ktau,nmaxpr)==0)then
#ifdef usenc3
  if (any(globvar==nf_fill_float)) then
#else
  if (any(globvar==nf90_fill_float)) then
#endif
    write(6,'("histwrt4 ",a7,i4,a7)') sname,iarch,"missing"
  else
    varn = minval(globvar)
    varx = maxval(globvar)
    max_result = maxloc(globvar)
    kmx = max_result(2)
    iq = max_result(1)
    ! Convert this 1D index to 2D
    imx = 1 + modulo(iq-1,il_g)
    jmx = 1 + (iq-1)/il_g
    write(6,'("histwrt4 ",a7,i4,2f12.4,3i4,f12.4)') sname,iarch,varn,varx,imx,jmx,kmx,globvar(id+(jd-1)*il_g,nlv)
  end if
end if

return
end subroutine hw4a      

!--------------------------------------------------------------------

subroutine ccnf_open(fname,ncid,status)

use cc_mpi
#ifndef usenc3
use netcdf
#endif

implicit none

#ifdef usenc3
include 'netcdf.inc'
#endif

integer, intent(out) :: ncid
integer, intent(out), optional :: status
integer ncstatus
character(len=*), intent(in) :: fname

#ifdef usenc3
ncstatus = nf_open(fname,nf_nowrite,ncid)
#else
ncstatus = nf90_open(fname,nf90_nowrite,ncid)
#endif

if (present(status)) then
  status=ncstatus
else
  call ncmsg("open",ncstatus)
endif

return
end subroutine ccnf_open

subroutine ccnf_create(fname,ncid)

use cc_mpi
#ifndef usenc3
use netcdf
#endif

implicit none

#ifdef usenc3
include 'netcdf.inc'
#endif

integer, intent(out) :: ncid
integer ncstatus
character(len=*), intent(in) :: fname

#ifdef usenc3
  ncstatus = nf_create(fname,nf_clobber,ncid)
#else
  ncstatus = nf90_create(fname,nf90_netcdf4,ncid)
#endif
call ncmsg("create",ncstatus)

return
end subroutine ccnf_create

subroutine ccnf_close(ncid)

use cc_mpi
#ifndef usenc3
use netcdf
#endif

implicit none

#ifdef usenc3
include 'netcdf.inc'
#endif

integer, intent(in) :: ncid
integer ncstatus

#if usenc3
ncstatus = nf_close(ncid)
#else
ncstatus = nf90_close(ncid)
#endif
call ncmsg("close",ncstatus)

return
end subroutine ccnf_close

subroutine ccnf_nofill(ncid)

use cc_mpi
#ifndef usenc3
use netcdf
#endif

implicit none

#ifdef usenc3
include 'netcdf.inc'
#endif

integer, intent(in) :: ncid
integer ncstatus,lomode

#ifdef usenc3
ncstatus=nf_set_fill(ncid,nf_nofill,lomode)
#else
ncstatus=nf90_set_fill(ncid,nf90_nofill,lomode)
#endif
call ncmsg("nofill",ncstatus)

return
end subroutine ccnf_nofill

subroutine ccnf_inq_dimid(ncid,dname,did,tst)

use cc_mpi
#ifndef usenc3
use netcdf
#endif

implicit none

#ifdef usenc3
include 'netcdf.inc'
#endif

integer, intent(in) :: ncid
integer, intent(out) :: did
integer ncstatus
logical ltst
logical, intent(out), optional :: tst
character(len=*), intent(in) :: dname

#ifdef usenc3
ncstatus=nf_inq_dimid(ncid,dname,did)
ltst=(ncstatus/=nf_noerr)
#else
ncstatus=nf90_inq_dimid(ncid,dname,did)
ltst=(ncstatus/=nf90_noerr)
#endif

if (present(tst)) then
  tst=ltst
else
  call ncmsg("dimid",ncstatus)
end if

return
end subroutine ccnf_inq_dimid

subroutine ccnf_inq_dimlen(ncid,dname,dlen,failok)

use cc_mpi
#ifndef usenc3
use netcdf
#endif

implicit none

#ifdef usenc3
include 'netcdf.inc'
#endif

integer, intent(in) :: ncid
integer, intent(inout) :: dlen
integer ncstatus,ldid,ldlen
character(len=*), intent(in) :: dname
logical, intent(in), optional :: failok
logical ftest

ftest=.false.
if (present(failok)) ftest=failok

ldlen=dlen
#ifdef usenc3
ncstatus=nf_inq_dimid(ncid,dname,ldid)
if (ncstatus/=nf_noerr) then
  if (ftest) return
  write(6,*) nf_strerror(ncstatus)
  call ccmpi_abort(-1)
end if
ncstatus=nf_inq_dimlen(ncid,ldid,ldlen)
if (ncstatus/=nf_noerr) then
  if (ftest) return
  write(6,*) nf_strerror(ncstatus)
  call ccmpi_abort(-1)
end if
#else
ncstatus=nf90_inq_dimid(ncid,dname,ldid)
if (ncstatus/=nf90_noerr) then
  if (ftest) return
  write(6,*) nf90_strerror(ncstatus)
  call ccmpi_abort(-1)
end if
ncstatus=nf90_inquire_dimension(ncid,ldid,len=ldlen)
if (ncstatus/=nf90_noerr) then
  if (ftest) return
  write(6,*) nf90_strerror(ncstatus)
  call ccmpi_abort(-1)
end if
#endif
dlen=ldlen

return
end subroutine ccnf_inq_dimlen

subroutine ccnf_inq_varid(ncid,vname,vid,tst)

use cc_mpi
#ifndef usenc3
use netcdf
#endif

implicit none

#ifdef usenc3
include 'netcdf.inc'
#endif

integer, intent(in) :: ncid
integer, intent(out) :: vid
integer ncstatus
character(len=*), intent(in) :: vname
logical, intent(out), optional :: tst
logical ltst

#ifdef usenc3
ncstatus = nf_inq_varid(ncid,vname,vid)
ltst=(ncstatus/=nf_noerr)
#else
ncstatus = nf90_inq_varid(ncid,vname,vid)
ltst=(ncstatus/=nf90_noerr)
#endif

if (present(tst)) then
  tst=ltst
else
  call ncmsg("varid",ncstatus)
end if

return
end subroutine ccnf_inq_varid

subroutine ccnf_def_dim(ncid,dname,nsize,did)

use cc_mpi
#ifndef usenc3
use netcdf
#endif

implicit none

#ifdef usenc3
include 'netcdf.inc'
#endif

integer, intent(in) :: ncid,nsize
integer, intent(out) :: did
integer ncstatus
character(len=*), intent(in) :: dname

#ifdef usenc3
ncstatus=nf_def_dim(ncid,dname,nsize,did)
#else
ncstatus=nf90_def_dim(ncid,dname,nsize,did)
#endif
call ncmsg("def_dim",ncstatus)

return
end subroutine ccnf_def_dim

subroutine ccnf_def_dimu(ncid,dname,did)

use cc_mpi
#ifndef usenc3
use netcdf
#endif

implicit none

#ifdef usenc3
include 'netcdf.inc'
#endif

integer, intent(in) :: ncid
integer, intent(out) :: did
integer ncstatus
character(len=*), intent(in) :: dname

#ifdef usenc3
ncstatus=nf_def_dim(ncid,dname,nf_unlimited,did)
#else
ncstatus=nf90_def_dim(ncid,dname,nf90_unlimited,did)
#endif
call ncmsg("def_dimu",ncstatus)

return
end subroutine ccnf_def_dimu

subroutine ccnf_def_var(ncid,vname,vtype,vndim,dims,vid,deflate)

use cc_mpi
#ifndef usenc3
use netcdf
#endif

implicit none

#ifdef usenc3
include 'netcdf.inc'
#endif

integer, intent(in) :: ncid,vndim
integer, intent(in), optional :: deflate
integer, intent(out) :: vid
integer, dimension(vndim), intent(in) :: dims
integer ncstatus,lvid,ltype,ldef
character(len=*), intent(in) :: vname
character(len=*), intent(in) :: vtype

select case(vtype)
#ifdef usenc3
  case('int')
    ltype=nf_int
  case('float')
    ltype=nf_float
  case('double')
    ltype=nf_double
  case('char')
    ltype=nf_char
#else
  case('int')
    ltype=nf90_int
  case('float')
    ltype=nf90_float
  case('double')
    ltype=nf90_double
  case('char')
    ltype=nf90_char
#endif
  case default
    write(6,*) "ERROR: Unknown option for ccnf_def_var ",vtype
    call ccmpi_abort(-1)
end select

#ifdef usenc3
ncstatus = nf_def_var(ncid,vname,ltype,vndim,dims,vid)
#else
ldef=1
if (present(deflate)) ldef=deflate
ncstatus = nf90_def_var(ncid,vname,ltype,dims,vid,deflate_level=ldef)
#endif
call ncmsg("def_var",ncstatus)

return
end subroutine ccnf_def_var

subroutine ccnf_def_var0(ncid,vname,vtype,vid)

use cc_mpi
#ifndef usenc3
use netcdf
#endif

implicit none

#ifdef usenc3
include 'netcdf.inc'
#endif

integer, intent(in) :: ncid
integer, intent(out) :: vid
integer ncstatus,ltype
character(len=*), intent(in) :: vname
character(len=*), intent(in) :: vtype

select case(vtype)
#ifdef usenc3
  case('int')
    ltype=nf_int
  case('float')
    ltype=nf_float
  case('double')
    ltype=nf_double
  case('char')
    ltype=nf_char
#else
  case('int')
    ltype=nf90_int
  case('float')
    ltype=nf90_float
  case('double')
    ltype=nf90_double
  case('char')
    ltype=nf90_char
#endif
  case default
    write(6,*) "ERROR: Unknown option for ccnf_def_var ",vtype
    call ccmpi_abort(-1)
end select

#ifdef usenc3
ncstatus = nf_def_var(ncid,vname,ltype,0,1,vid)
#else
ncstatus = nf90_def_var(ncid,vname,ltype,1,vid)
#endif
call ncmsg("def_var0",ncstatus)

return
end subroutine ccnf_def_var0

subroutine ccnf_enddef(ncid)

use cc_mpi
#ifndef usenc3
use netcdf
#endif

implicit none

#ifdef usenc3
include 'netcdf.inc'
#endif

integer, intent(in) :: ncid
integer ncstatus

#ifdef usenc3
ncstatus=nf_enddef(ncid)
#else
ncstatus=nf90_enddef(ncid)
#endif
call ncmsg("enddef",ncstatus)

return
end subroutine ccnf_enddef

subroutine ccnf_redef(ncid)

use cc_mpi
#ifndef usenc3
use netcdf
#endif

implicit none

#ifdef usenc3
include 'netcdf.inc'
#endif

integer, intent(in) :: ncid
integer ncstatus

#ifdef usenc3
ncstatus=nf_redef(ncid)
#else
ncstatus=nf90_redef(ncid)
#endif
call ncmsg("redef",ncstatus)

return
end subroutine ccnf_redef

subroutine ccnf_sync(ncid)

use cc_mpi
#ifndef usenc3
use netcdf
#endif

implicit none

#ifdef usenc3
include 'netcdf.inc'
#endif

integer, intent(in) :: ncid
integer ncstatus

#ifdef outsync
#ifdef usenc3
ncstatus = nf_sync(ncid)
#else
ncstatus = nf90_sync(ncid)
#endif
call ncmsg("sync",ncstatus)
#endif

return
end subroutine ccnf_sync

subroutine ccnf_get_var_real(ncid,vname,vdat)

use cc_mpi
#ifndef usenc3
use netcdf
#endif

implicit none

#ifdef usenc3
include 'netcdf.inc'
#endif

integer, intent(in) :: ncid
integer lvid,ncstatus
real, dimension(:), intent(out) :: vdat
character(len=*), intent(in) :: vname

#ifdef usenc3
ncstatus = nf_inq_varid(ncid,vname,lvid)
call ncmsg("get_var_varid",ncstatus)
ncstatus = nf_get_var_real(ncid,lvid,vdat)
call ncmsg("get_var",ncstatus)
#else
ncstatus = nf90_inq_varid(ncid,vname,lvid)
call ncmsg("get_var_varid",ncstatus)
ncstatus = nf90_get_var(ncid,lvid,vdat)
call ncmsg("get_var",ncstatus)
#endif

return
end subroutine ccnf_get_var_real

subroutine ccnf_get_var_int(ncid,vname,vdat)

use cc_mpi
#ifndef usenc3
use netcdf
#endif

implicit none

#ifdef usenc3
include 'netcdf.inc'
#endif

integer, intent(in) :: ncid
integer lvid,ncstatus
integer, dimension(:), intent(out) :: vdat
character(len=*), intent(in) :: vname

#ifdef usenc3
ncstatus = nf_inq_varid(ncid,vname,lvid)
call ncmsg("get_var_varid",ncstatus)
ncstatus = nf_get_var_int(ncid,lvid,vdat)
call ncmsg("get_var",ncstatus)
#else
ncstatus = nf90_inq_varid(ncid,vname,lvid)
call ncmsg("get_var_varid",ncstatus)
ncstatus = nf90_get_var(ncid,lvid,vdat)
call ncmsg("get_var",ncstatus)
#endif

return
end subroutine ccnf_get_var_int

subroutine ccnf_get_var1_real(ncid,vid,start,vdat)

use cc_mpi
#ifndef usenc3
use netcdf
#endif

implicit none

#ifdef usenc3
include 'netcdf.inc'
#endif

integer, intent(in) :: ncid,vid,start
integer ncstatus
#ifndef usenc3
integer, dimension(1) :: lstart,lcount
real, dimension(1) :: ldat
#endif
real, intent(out) :: vdat

#ifdef usenc3
ncstatus = nf_get_var1_real(ncid,vid,start,vdat)
#else
lstart=start
lcount=1
ncstatus = nf90_get_var(ncid,vid,ldat,start=lstart,count=lcount)
vdat=ldat(1)
#endif
call ncmsg("get_var1",ncstatus)

return
end subroutine ccnf_get_var1_real

subroutine ccnf_get_var1_int(ncid,vid,start,vdat)

use cc_mpi
#ifndef usenc3
use netcdf
#endif

implicit none

#ifdef usenc3
include 'netcdf.inc'
#endif

integer, intent(in) :: ncid,vid,start
integer ncstatus
#ifndef usenc3
integer, dimension(1) :: lstart,lcount
integer, dimension(1) :: ldat
#endif
integer, intent(out) :: vdat

#ifdef usenc3
ncstatus = nf_get_var1_int(ncid,vid,start,vdat)
#else
lstart=start
lcount=1
ncstatus = nf90_get_var(ncid,vid,ldat,start=lstart,count=lcount)
vdat=ldat(1)
#endif
call ncmsg("get_var1",ncstatus)

return
end subroutine ccnf_get_var1_int

subroutine ccnf_get_vara_real2r(ncid,vid,start,ncount,vdat)

use cc_mpi
#ifndef usenc3
use netcdf
#endif

implicit none

#ifdef usenc3
include 'netcdf.inc'
#endif

integer, intent(in) :: ncid,vid
integer, dimension(:) :: start,ncount
integer ncstatus
real, dimension(:), intent(out) :: vdat

#ifdef usenc3
ncstatus=nf_get_vara_real(ncid,vid,start,ncount,vdat)
#else
ncstatus=nf90_get_var(ncid,vid,vdat,start=start,count=ncount)
#endif
call ncmsg("get_vara",ncstatus)

return
end subroutine ccnf_get_vara_real2r 

subroutine ccnf_get_vara_real3r(ncid,vid,start,ncount,vdat)

use cc_mpi
#ifndef usenc3
use netcdf
#endif

implicit none

#ifdef usenc3
include 'netcdf.inc'
#endif

integer, intent(in) :: ncid,vid
integer, dimension(:) :: start,ncount
integer ncstatus
real, dimension(:,:), intent(out) :: vdat

#ifdef usenc3
ncstatus=nf_get_vara_real(ncid,vid,start,ncount,vdat)
#else
ncstatus=nf90_get_var(ncid,vid,vdat,start=start,count=ncount)
#endif
call ncmsg("get_vara",ncstatus)

return
end subroutine ccnf_get_vara_real3r

subroutine ccnf_get_vara_real4r(ncid,vid,start,ncount,vdat)

use cc_mpi
#ifndef usenc3
use netcdf
#endif

implicit none

#ifdef usenc3
include 'netcdf.inc'
#endif

integer, intent(in) :: ncid,vid
integer, dimension(:) :: start,ncount
integer ncstatus
real, dimension(:,:,:), intent(out) :: vdat

#ifdef usenc3
ncstatus=nf_get_vara_real(ncid,vid,start,ncount,vdat)
#else
ncstatus=nf90_get_var(ncid,vid,vdat,start=start,count=ncount)
#endif
call ncmsg("get_vara",ncstatus)

return
end subroutine ccnf_get_vara_real4r

subroutine ccnf_get_vara_int2i(ncid,vid,start,ncount,vdat)

use cc_mpi
#ifndef usenc3
use netcdf
#endif

implicit none

#ifdef usenc3
include 'netcdf.inc'
#endif

integer, intent(in) :: ncid,vid
integer, dimension(:) :: start,ncount
integer ncstatus
integer, dimension(:), intent(out) :: vdat

#ifdef usenc3
ncstatus=nf_get_vara_int(ncid,vid,start,ncount,vdat)
#else
ncstatus=nf90_get_var(ncid,vid,vdat,start=start,count=ncount)
#endif
call ncmsg("get_vara",ncstatus)

return
end subroutine ccnf_get_vara_int2i

subroutine ccnf_get_vara_double4d(ncid,vid,start,ncount,vdat)

use cc_mpi
#ifndef usenc3
use netcdf
#endif

implicit none

#ifdef usenc3
include 'netcdf.inc'
#endif

integer, intent(in) :: ncid,vid
integer, dimension(:) :: start,ncount
integer ncstatus
double precision, dimension(:,:,:), intent(out) :: vdat

#ifdef usenc3
ncstatus=nf_get_vara_double(ncid,vid,start,ncount,vdat)
#else
ncstatus=nf90_get_var(ncid,vid,vdat,start=start,count=ncount)
#endif
call ncmsg("get_vara",ncstatus)

return
end subroutine ccnf_get_vara_double4d

subroutine ccnf_get_att_text(ncid,vid,aname,atext)

use cc_mpi
#ifndef usenc3
use netcdf
#endif

implicit none

#ifdef usenc3
include 'netcdf.inc'
#endif

integer, intent(in) :: ncid,vid
integer ncstatus
character(len=*), intent(in) :: aname
character(len=*), intent(out) :: atext

atext=''
#ifdef usenc3
ncstatus = nf_get_att_text(ncid,vid,aname,atext)
#else
ncstatus = nf90_get_att(ncid,vid,aname,atext)
#endif
call ncmsg("get_att",ncstatus)

return
end subroutine ccnf_get_att_text

subroutine ccnf_get_att_real(ncid,vid,aname,vdat)

use cc_mpi
#ifndef usenc3
use netcdf
#endif

implicit none

#ifdef usenc3
include 'netcdf.inc'
#endif

integer, intent(in) :: ncid,vid
integer ncstatus
character(len=*), intent(in) :: aname
real, intent(out) :: vdat

#ifdef usenc3
ncstatus = nf_get_att_real(ncid,vid,aname,vdat)
#else
ncstatus = nf90_get_att(ncid,vid,aname,vdat)
#endif
call ncmsg("get_att",ncstatus)


return
end subroutine ccnf_get_att_real

subroutine ccnf_get_att_realg1r(ncid,aname,vdat)

use cc_mpi
#ifndef usenc3
use netcdf
#endif

implicit none

#ifdef usenc3
include 'netcdf.inc'
#endif

integer, intent(in) :: ncid
integer ncstatus
character(len=*), intent(in) :: aname
real, intent(out) :: vdat

#ifdef usenc3
ncstatus = nf_get_att_real(ncid,nf_global,aname,vdat)
#else
ncstatus = nf90_get_att(ncid,nf90_global,aname,vdat)
#endif
call ncmsg("get_attg",ncstatus)

return
end subroutine ccnf_get_att_realg1r

subroutine ccnf_get_att_realg2r(ncid,aname,vdat)

use cc_mpi
#ifndef usenc3
use netcdf
#endif

implicit none

#ifdef usenc3
include 'netcdf.inc'
#endif

integer, intent(in) :: ncid
integer ncstatus
character(len=*), intent(in) :: aname
real, dimension(:), intent(out) :: vdat

#ifdef usenc3
ncstatus = nf_get_att_real(ncid,nf_global,aname,vdat)
#else
ncstatus = nf90_get_att(ncid,nf90_global,aname,vdat)
#endif
call ncmsg("get_attg",ncstatus)

return
end subroutine ccnf_get_att_realg2r

subroutine ccnf_get_att_intg1i(ncid,aname,vdat)

use cc_mpi
#ifndef usenc3
use netcdf
#endif

implicit none

#ifdef usenc3
include 'netcdf.inc'
#endif

integer, intent(in) :: ncid
integer ncstatus
character(len=*), intent(in) :: aname
integer, intent(out) :: vdat

#ifdef usenc3
ncstatus = nf_get_att_int(ncid,nf_global,aname,vdat)
#else
ncstatus = nf90_get_att(ncid,nf90_global,aname,vdat)
#endif
call ncmsg("get_attg",ncstatus)

return
end subroutine ccnf_get_att_intg1i

subroutine ccnf_get_att_intg2i(ncid,aname,vdat)

use cc_mpi
#ifndef usenc3
use netcdf
#endif

implicit none

#ifdef usenc3
include 'netcdf.inc'
#endif

integer, intent(in) :: ncid
integer ncstatus
character(len=*), intent(in) :: aname
integer, dimension(:), intent(out) :: vdat

#ifdef usenc3
ncstatus = nf_get_att_int(ncid,nf_global,aname,vdat)
#else
ncstatus = nf90_get_att(ncid,nf90_global,aname,vdat)
#endif
call ncmsg("get_attg",ncstatus)

return
end subroutine ccnf_get_att_intg2i

subroutine ccnf_read(fname,vname,vdat)

use cc_mpi
#ifndef usenc3
use netcdf
#endif

implicit none

include 'newmpar.h'
#ifdef usenc3
include 'netcdf.inc'
#endif

integer ncstatus,lncid,lvid
real, dimension(ifull), intent(out) :: vdat
real, dimension(ifull_g) :: vdat_g
character(len=*), intent(in) :: fname
character(len=*), intent(in) :: vname

if (myid==0) then
#ifdef usenc3
  ncstatus = nf_open(fname,nf_nowrite,lncid)
  call ncmsg(fname,ncstatus)
  ncstatus = nf_inq_varid(lncid,vname,lvid)
  call ncmsg(fname,ncstatus)
  ncstatus = nf_get_var_real(lncid,lvid,vdat_g)
  call ncmsg(fname,ncstatus)
  ncstatus = nf_close(lncid)
  call ncmsg(fname,ncstatus)
#else
  ncstatus = nf90_open(fname,nf90_nowrite,lncid)
  call ncmsg(fname,ncstatus)
  ncstatus = nf90_inq_varid(lncid,vname,lvid)
  call ncmsg(fname,ncstatus)
  ncstatus = nf90_get_var(lncid,lvid,vdat_g)
  call ncmsg(fname,ncstatus)
  ncstatus = nf90_close(lncid)
  call ncmsg(fname,ncstatus)
#endif
  call ccmpi_distribute(vdat,vdat_g)
else
  call ccmpi_distribute(vdat)
end if

return
end subroutine ccnf_read

subroutine ccnf_put_var_text2r(ncid,vid,vtxt)

use cc_mpi
#ifndef usenc3
use netcdf
#endif

implicit none

#ifdef usenc3
include 'netcdf.inc'
#endif

integer, intent(in) :: ncid,vid
integer ncstatus
character(len=*), dimension(:), intent(in) :: vtxt

#ifdef usenc3
ncstatus = nf_put_var_text(ncid,vid,vtxt)
#else
ncstatus = nf90_put_var(ncid,vid,vtxt)
#endif
call ncmsg("put_var",ncstatus)

return
end subroutine ccnf_put_var_text2r

subroutine ccnf_put_var_int2i(ncid,vid,vdat)

use cc_mpi
#ifndef usenc3
use netcdf
#endif

implicit none

#ifdef usenc3
include 'netcdf.inc'
#endif

integer, intent(in) :: ncid,vid
integer ncstatus
integer, dimension(:), intent(in) :: vdat

#ifdef usenc3
ncstatus = nf_put_var_int(ncid,vid,vdat)
#else
ncstatus = nf90_put_var(ncid,vid,vdat)
#endif
call ncmsg("put_var",ncstatus)

return
end subroutine ccnf_put_var_int2i

subroutine ccnf_put_var_int3i(ncid,vid,vdat)

use cc_mpi
#ifndef usenc3
use netcdf
#endif

implicit none

#ifdef usenc3
include 'netcdf.inc'
#endif

integer, intent(in) :: ncid,vid
integer ncstatus
integer, dimension(:,:), intent(in) :: vdat

#ifdef usenc3
ncstatus = nf_put_var_int(ncid,vid,vdat)
#else
ncstatus = nf90_put_var(ncid,vid,vdat)
#endif
call ncmsg("put_var",ncstatus)

return
end subroutine ccnf_put_var_int3i

subroutine ccnf_put_var1_int(ncid,vid,start,vdat)

use cc_mpi
#ifndef usenc3
use netcdf
#endif

implicit none

#ifdef usenc3
include 'netcdf.inc'
#endif

integer, intent(in) :: ncid,vid,start
integer ncstatus
#ifndef usenc3
integer, dimension(1) :: lstart,lcount
integer, dimension(1) :: ldat
#endif
integer, intent(in) :: vdat

#ifdef usenc3
ncstatus=nf_put_var1_int(ncid,vid,start,vdat)
#else
lstart=start
lcount=1
ldat(1)=vdat
ncstatus=nf90_put_var(ncid,vid,ldat,start=lstart,count=lcount)
#endif
call ncmsg("put_var1",ncstatus)

return
end subroutine ccnf_put_var1_int

subroutine ccnf_put_var1_real(ncid,vid,start,vdat)

use cc_mpi
#ifndef usenc3
use netcdf
#endif

implicit none

#ifdef usenc3
include 'netcdf.inc'
#endif

integer, intent(in) :: ncid,vid,start
integer ncstatus
#ifndef usenc3
integer, dimension(1) :: lstart,lcount
real, dimension(1) :: ldat
#endif
real, intent(in) :: vdat

#ifdef usenc3
#ifdef i8r8
ncstatus=nf_put_var1_double(ncid,vid,start,vdat)
#else
ncstatus=nf_put_var1_real(ncid,vid,start,vdat)
#endif
#else
lstart=start
lcount=1
ldat(1)=vdat
ncstatus=nf90_put_var(ncid,vid,ldat,start=lstart,count=lcount)
#endif
call ncmsg("put_var1",ncstatus)

return
end subroutine ccnf_put_var1_real

subroutine ccnf_put_var1_double(ncid,vid,start,vdat)

use cc_mpi
#ifndef usenc3
use netcdf
#endif

implicit none

#ifdef usenc3
include 'netcdf.inc'
#endif

integer, intent(in) :: ncid,vid,start
integer ncstatus
#ifndef usenc3
integer, dimension(1) :: lstart,lcount
double precision, dimension(1) :: ldat
#endif
double precision, intent(in) :: vdat

#ifdef usenc3
ncstatus=nf_put_var1_double(ncid,vid,start,vdat)
#else
lstart=start
lcount=1
ldat(1)=vdat
ncstatus=nf90_put_var(ncid,vid,ldat,start=lstart,count=lcount)
#endif
call ncmsg("put_var1",ncstatus)

return
end subroutine ccnf_put_var1_double

subroutine ccnf_put_vara_real2r(ncid,vid,start,ncount,vdat)

use cc_mpi
#ifndef usenc3
use netcdf
#endif

implicit none

#ifdef usenc3
include 'netcdf.inc'
#endif

integer, intent(in) :: ncid,vid
integer ncstatus
integer, dimension(:), intent(in) :: start,ncount
real, dimension(:), intent(in) :: vdat

#ifdef usenc3
#ifdef i8r8
ncstatus=nf_put_vara_double(ncid,vid,start,ncount,vdat)
#else
ncstatus=nf_put_vara_real(ncid,vid,start,ncount,vdat)
#endif
#else
ncstatus=nf90_put_var(ncid,vid,vdat,start=start,count=ncount)
#endif
call ncmsg("put_vara",ncstatus)

return
end subroutine ccnf_put_vara_real2r

subroutine ccnf_put_vara_real3r(ncid,vid,start,ncount,vdat)

use cc_mpi
#ifndef usenc3
use netcdf
#endif

implicit none

#ifdef usenc3
include 'netcdf.inc'
#endif

integer, intent(in) :: ncid,vid
integer ncstatus
integer, dimension(:), intent(in) :: start,ncount
real, dimension(:,:), intent(in) :: vdat

#ifdef usenc3
#ifdef i8r8
ncstatus=nf_put_vara_double(ncid,vid,start,ncount,vdat)
#else
ncstatus=nf_put_vara_real(ncid,vid,start,ncount,vdat)
#endif
#else
ncstatus=nf90_put_var(ncid,vid,vdat,start=start,count=ncount)
#endif
call ncmsg("put_vara",ncstatus)

return
end subroutine ccnf_put_vara_real3r

subroutine ccnf_put_vara_double2r(ncid,vid,start,ncount,vdat)

use cc_mpi
#ifndef usenc3
use netcdf
#endif

implicit none

#ifdef usenc3
include 'netcdf.inc'
#endif

integer, intent(in) :: ncid,vid
integer ncstatus
integer, dimension(:), intent(in) :: start,ncount
double precision, dimension(:), intent(in) :: vdat

#ifdef usenc3
ncstatus=nf_put_vara_double(ncid,vid,start,ncount,vdat)
#else
ncstatus=nf90_put_var(ncid,vid,vdat,start=start,count=ncount)
#endif
call ncmsg("put_vara",ncstatus)

return
end subroutine ccnf_put_vara_double2r

subroutine ccnf_put_vara_int2i(ncid,vid,start,ncount,vdat)

use cc_mpi
#ifndef usenc3
use netcdf
#endif

implicit none

#ifdef usenc3
include 'netcdf.inc'
#endif

integer, intent(in) :: ncid,vid
integer ncstatus
integer, dimension(:), intent(in) :: start,ncount
integer, dimension(:), intent(in) :: vdat

#ifdef usenc3
ncstatus=nf_put_vara_int(ncid,vid,start,ncount,vdat)
#else
ncstatus=nf90_put_var(ncid,vid,vdat,start=start,count=ncount)
#endif
call ncmsg("put_vara",ncstatus)

return
end subroutine ccnf_put_vara_int2i

subroutine ccnf_put_att_text(ncid,vid,aname,asize,atext)

use cc_mpi
#ifndef usenc3
use netcdf
#endif

implicit none

#ifdef usenc3
include 'netcdf.inc'
#endif

integer, intent(in) :: ncid,vid,asize
integer ncstatus
character(len=*), intent(in) :: aname
character(len=asize), intent(in) :: atext

#ifdef usenc3
ncstatus=nf_put_att_text(ncid,vid,aname,asize,atext)
#else
ncstatus=nf90_put_att(ncid,vid,aname,atext)
#endif
call ncmsg("put_att",ncstatus)

return
end subroutine ccnf_put_att_text

subroutine ccnf_put_att_textg(ncid,aname,atext)

use cc_mpi
#ifndef usenc3
use netcdf
#endif

implicit none

#ifdef usenc3
include 'netcdf.inc'
#endif

integer, intent(in) :: ncid
integer(kind=4) ncstatus,lsize
character(len=*), intent(in) :: aname
character(len=*), intent(in) :: atext

lsize=len_trim(atext)
#ifdef usenc3
ncstatus=nf_put_att_text(ncid,nf_global,aname,lsize,atext)
#else
ncstatus=nf90_put_att(ncid,nf90_global,aname,atext)
#endif
call ncmsg("put_attg",ncstatus)

return
end subroutine ccnf_put_att_textg

subroutine ccnf_put_att_intg1(ncid,aname,vdat)

use cc_mpi
#ifndef usenc3
use netcdf
#endif

implicit none

#ifdef usenc3
include 'netcdf.inc'
#endif

integer, intent(in) :: ncid
integer ncstatus
integer, intent(in) :: vdat
integer, dimension(1) :: ldat
character(len=*), intent(in) :: aname

ldat(1)=vdat
#ifdef usenc3
ncstatus=nf_put_att_int(ncid,nf_global,aname,nf_int,1,ldat)
#else
ncstatus=nf90_put_att(ncid,nf90_global,aname,ldat)
#endif
call ncmsg("put_attg",ncstatus)

return
end subroutine ccnf_put_att_intg1

subroutine ccnf_put_att_intg2(ncid,aname,vsize,vdat)

use cc_mpi
#ifndef usenc3
use netcdf
#endif

implicit none

#ifdef usenc3
include 'netcdf.inc'
#endif

integer, intent(in) :: ncid,vsize
integer ncstatus
integer, dimension(:), intent(in) :: vdat
character(len=*), intent(in) :: aname

#ifdef usenc3
ncstatus=nf_put_att_int(ncid,nf_global,aname,nf_int,vsize,vdat)
#else
ncstatus=nf90_put_att(ncid,nf90_global,aname,vdat)
#endif
call ncmsg("put_att",ncstatus)

return
end subroutine ccnf_put_att_intg2

subroutine ccnf_put_att_realg1(ncid,aname,vdat)

use cc_mpi
#ifndef usenc3
use netcdf
#endif

implicit none

#ifdef usenc3
include 'netcdf.inc'
#endif

integer, intent(in) :: ncid
integer ncstatus
real, intent(in) :: vdat
real, dimension(1) :: ldat
character(len=*), intent(in) :: aname

ldat(1)=vdat
#ifdef usenc3
#ifdef i8r8
ncstatus=nf_put_att_double(ncid,nf_global,aname,nf_float,1,ldat)
#else
ncstatus=nf_put_att_real(ncid,nf_global,aname,nf_float,1,ldat)
#endif
#else
ncstatus=nf90_put_att(ncid,nf90_global,aname,ldat)
#endif
call ncmsg("put_attg",ncstatus)

return
end subroutine ccnf_put_att_realg1

subroutine ccnf_put_att_realg2(ncid,aname,vsize,vdat)

use cc_mpi
#ifndef usenc3
use netcdf
#endif

implicit none

#ifdef usenc3
include 'netcdf.inc'
#endif

integer, intent(in) :: ncid,vsize
integer ncstatus
real, dimension(:), intent(in) :: vdat
character(len=*), intent(in) :: aname

#ifdef usenc3
#ifdef i8r8
ncstatus=nf_put_att_double(ncid,nf_global,aname,nf_float,vsize,vdat)
#else
ncstatus=nf_put_att_real(ncid,nf_global,aname,nf_float,vsize,vdat)
#endif
#else
ncstatus=nf90_put_att(ncid,nf90_global,aname,vdat)
#endif
call ncmsg("put_attg",ncstatus)

return
end subroutine ccnf_put_att_realg2

end module infile
