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
public histrd1,histrd4s,vertint,datefix,getzinp,ncmsg
public histopen,histclose
public attrib,histwrt3,histwrt4,freqwrite
public ccnf_open,ccnf_create,ccnf_close,ccnf_sync,ccnf_enddef,ccnf_redef,ccnf_nofill,ccnf_inq_varid,ccnf_inq_dimid
public ccnf_inq_dimlen,ccnf_def_dim,ccnf_def_dimu,ccnf_def_var,ccnf_def_var0,ccnf_get_var_real,ccnf_get_var_realg
public ccnf_get_var_int,ccnf_get_att_intg,ccnf_get_vara_real,ccnf_get_vara_int,ccnf_get_vara_double
public ccnf_get_var1_int,ccnf_get_var1_real,ccnf_get_att_text,ccnf_get_att_real,ccnf_get_att_realg
public ccnf_read,ccnf_put_var_text,ccnf_put_var_int,ccnf_put_var1_int,ccnf_put_var1_real,ccnf_put_var1_double
public ccnf_put_vara_int,ccnf_put_vara_real,ccnf_put_vara_double,ccnf_put_att_text,ccnf_put_att_textg
public ccnf_put_att_intg,ccnf_put_att_realg

interface ccnf_get_att_intg
  module procedure ccnf_get_att_intg1i, ccnf_get_att_intg2i
end interface ccnf_get_att_intg
interface ccnf_get_att_realg
  module procedure ccnf_get_att_realg1r, ccnf_get_att_realg2r
end interface ccnf_get_att_realg
interface ccnf_get_vara_real
  module procedure ccnf_get_vara_real2r, ccnf_get_vara_real3r, ccnf_get_vara_real4r 
end interface ccnf_get_vara_real
interface ccnf_get_vara_int
  module procedure ccnf_get_vara_int2i 
end interface ccnf_get_vara_int
interface ccnf_get_vara_double
  module procedure ccnf_get_vara_double4d
end interface ccnf_get_vara_double
interface ccnf_put_vara_real
  module procedure ccnf_put_vara_real2r, ccnf_put_vara_real3r
end interface ccnf_put_vara_real
interface ccnf_put_vara_double
  module procedure ccnf_put_vara_double2r
end interface ccnf_put_vara_double
interface ccnf_put_vara_int
  module procedure ccnf_put_vara_int2i
end interface ccnf_put_vara_int
interface ccnf_put_var_text
  module procedure ccnf_put_var_text2r
end interface ccnf_put_var_text
interface ccnf_put_var_int
  module procedure ccnf_put_var_int2i, ccnf_put_var_int3i
end interface ccnf_put_var_int

integer, dimension(:), allocatable, save :: pnoff
integer, dimension(:,:), allocatable, save :: pioff,pjoff
integer(kind=4), dimension(:), allocatable, save :: pncid
integer(kind=4), dimension(:), allocatable, save :: pncidold
integer(kind=2), parameter :: minv = -32500
integer(kind=2), parameter :: maxv = 32500
integer(kind=2), parameter :: missval = -32501
integer, save :: mynproc,fnproc
integer, save :: pil_g,pjl_g,pil,pjl,pnpan
integer, save :: comm_ip,comm_ipold
logical, save :: ptest
      
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
      
implicit none

include 'newmpar.h'
include 'netcdf.inc'

integer, intent(in) :: iarchi
integer, intent(out) :: ier
integer(kind=2), dimension(pil*pjl*pnpan) :: ivar
integer(kind=4), dimension(3) :: start,count
integer ipf,ip,jpf,jpmax,iptst2,n
integer ierb,ca,cb
integer no,i,j,iq,iqi,lcomm
integer(kind=4) :: ler,idv,nctype
real, dimension(:), intent(inout), optional :: var
real(kind=4), dimension(pil*pjl*pnpan) :: lvar
real, dimension(pil*pjl*pnpan) :: rvar
real, dimension(pil*pjl*pnpan*nproc) :: gvar 
real(kind=4) laddoff,lsf
real addoff,sf
double precision, dimension(pil*pjl*pnpan) :: dvar
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
  ivar=0

  jpmax=nproc
  if (ipf==mynproc-1.and.myid<iptst2) then
    jpmax=iptst2
  end if
  
  ! get variable idv
  ler=nf_inq_varid(pncid(ipf),name,idv)
  ier=ler
  if(ier/=0)then
    if (myid==0.and.ipf==0) then
      write(6,*) '***absent field for ncid,name,idv,ier: ',pncid(0),name,idv,ier
    end if
  else
    ! obtain scaling factors and offsets from attributes
    addoff=0.
    sf=1.
    call ncagt(pncid(ipf),idv,'add_offset',laddoff,ler)
    if (ler==0) addoff=laddoff
    call ncagt(pncid(ipf),idv,'scale_factor',lsf,ler)
    if (ler==0) sf=lsf
    ! read in all data
    ler=nf_inq_vartype(pncid(ipf),idv,nctype)
    if (ler/=0) then
      write(6,*) nf_strerror(ler)
      call ccmpi_abort(-1)
    end if
    select case(nctype)
      case(nf_double)
        call ncvgt(pncid(ipf),idv,start,count,dvar,ler)
        rvar=dvar(:)
      case(nf_float)
        call ncvgt(pncid(ipf),idv,start,count,lvar,ler)
        rvar(:)=lvar(:)
      case(nf_short)
        call ncvgt(pncid(ipf),idv,start,count,ivar,ler)
        rvar(:)=real(ivar(:))
      case DEFAULT
        write(6,*) "ERROR: Unknown NetCDF format"
        call ccmpi_abort(-1)
    end select
    if (ler/=0) then
      write(6,*) nf_strerror(ler)
      call ccmpi_abort(-1)
    end if
    ier=ler
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
  
    call ccmpi_gatherx(gvar,rvar,0,lcomm)
    if (myid==0) then
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
      
implicit none

include 'newmpar.h'
include 'netcdf.inc'

integer, intent(in) :: iarchi,kk
integer, intent(out) :: ier
integer(kind=2), dimension(pil*pjl*pnpan) :: ivar
integer, dimension(4) :: start,count
integer ipf,ip,jpf,jpmax,iptst2,n
integer ierb,ca,cb
integer no,i,j,iq,iqi,lcomm
integer(kind=4) idv,ler,nctype
real, dimension(:), intent(inout), optional :: var
real(kind=4), dimension(pil*pjl*pnpan) :: lvar
real, dimension(pil*pjl*pnpan) :: rvar
real, dimension(pil*pjl*pnpan*nproc) :: gvar
real addoff,sf
real(kind=4) laddoff,lsf
double precision, dimension(pil*pjl*pnpan) :: dvar
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
count = (/ pil, pjl*pnpan, 1, 1 /)
ier=0

iptst2=mod(fnproc,nproc)
      
do ipf=0,mynproc-1
  rvar=0.
  ivar=0

  jpmax=nproc
  if (ipf==mynproc-1.and.myid<iptst2) then
    jpmax=iptst2
  end if

  ! get variable idv
  ler=nf_inq_varid(pncid(ipf),name,idv)
  ier=ler
  if(ier/=0)then
    if (myid==0.and.ipf==0) then
      write(6,*) '***absent field for ncid,name,idv,ier: ',pncid(0),name,idv,ier
    end if
  else
    ! obtain scaling factors and offsets from attributes
    addoff=0.
    sf=1.
    call ncagt(pncid(ipf),idv,'add_offset',laddoff,ler)
    if (ler==0) addoff=laddoff
    call ncagt(pncid(ipf),idv,'scale_factor',lsf,ler)
    if (ler==0) sf=lsf
    ! read in all data
    ler=nf_inq_vartype(pncid(ipf),idv,nctype)
    if (ler/=0) then
      write(6,*) nf_strerror(ler)
      call ccmpi_abort(-1)
    end if
    select case(nctype)
      case(nf_double)
        call ncvgt(pncid(ipf),idv,start,count,dvar,ler)
        rvar=dvar(:)
      case(nf_float)
        call ncvgt(pncid(ipf),idv,start,count,lvar,ler)
        rvar(:)=lvar(:)
      case(nf_short)
        call ncvgt(pncid(ipf),idv,start,count,ivar,ler)
        rvar(:)=real(ivar(:))
      case DEFAULT
        write(6,*) "ERROR: Unknown NetCDF format"
        call ccmpi_abort(-1)
    end select
    if (ler/=0) then
      write(6,*) nf_strerror(ler)
      call ccmpi_abort(-1)
    end if
    ier=ler
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
  
    call ccmpi_gatherx(gvar,rvar,0,lcomm)
    if (myid==0) then
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

include 'netcdf.inc'

integer, intent(in) :: ncstatus
integer(kind=4) ierr
character(len=*), intent(in) :: txt

ierr=ncstatus
if (ierr/=nf_noerr) then
  write(6,*) txt," ",nf_strerror(ierr)
  call ccmpi_abort(-1)
end if

return
end subroutine ncmsg

!--------------------------------------------------------------
! This subroutine opens parallel input files
subroutine histopen(ncid,ifile,ier)
      
use cc_mpi
      
implicit none
      
include 'newmpar.h'
include 'netcdf.inc'
      
integer, parameter :: nihead   = 54
      
integer(kind=4), dimension(nihead) :: ahead
integer, dimension(0:5) :: duma,dumb
integer, dimension(6) :: idum
integer, intent(out) :: ncid,ier
integer(kind=4) lncid,ler,lfnproc,dumr
integer resid,is,ipf,dmode
integer ipin,nxpr,nypr
integer ltst
character(len=*), intent(in) :: ifile
character(len=170) :: pfile
character(len=7) :: fdecomp
      
call ncpopt(0) ! turn off fatal errors on all processors
if (myid==0) then
  ! attempt to open single file with myid==0
  ler=nf_open(ifile,nf_nowrite,lncid)
  ncid=lncid
  ier=ler
  fnproc=1
  dmode=0
  pil=0
  pjl=0
  pnpan=0
  ptest=.false.
      
  ! attempt to open parallel files
  if (ier/=0) then
    write(pfile,"(a,'.',i6.6)") trim(ifile), 0
    ler=nf_open(pfile,nf_nowrite,lncid)
    ncid=lncid
    ier=ler
    if (ier/=0) then
      write(6,*) "WARN: Cannot open ",pfile
    else  
      write(6,*) "Found parallel input file ",ifile
      ler=nf_get_att_int(lncid,nf_global,"nproc",lfnproc)
      if (ler/=0) then
        write(6,*) nf_strerror(ler)
        call ccmpi_abort(-1)
      end if
      fnproc=lfnproc
      fdecomp=''
      ler=nf_get_att_text(lncid,nf_global,"decomp",fdecomp)
      if (ler/=0) then
        write(6,*) nf_strerror(ler)
        call ccmpi_abort(-1)
      end if
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
    ler=nf_get_att_int(lncid,nf_global,"nproc",dumr)
    ! nproc should only exist in multi-file input
    if (ler==0) then
      write(6,*) "ERROR: Incorrect base filename"
      call ccmpi_abort(-1)
    end if
    write(6,*) "Found single input file ",ifile
  end if

  if (ier==0) then
    ler=nf_get_att_int(lncid,nf_global,"int_header",ahead)
    if (ler/=0) then
      write(6,*) nf_strerror(ler)
      call ccmpi_abort(-1)
    end if
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

if (ier/=0) return

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
  ler=nf_open(pfile,nf_nowrite,pncid(ipf))
  if (ler/=0) then
    write(6,*) "ERROR: Cannot open ",pfile
    write(6,*) nf_strerror(ler)
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
      
implicit none
    
include 'netcdf.inc'
     
integer ipf,ipin,plen
integer(kind=4) ierr
      
if (allocated(pncidold)) then
  plen=size(pncidold)
  do ipf=0,plen-1
    ierr=nf_close(pncidold(ipf))
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

implicit none

include 'netcdf.inc'

integer, intent(in) :: cdfid, itype, ndim
integer, intent(in) :: daily
integer, dimension(ndim), intent(in) :: dim
integer(kind=4), dimension(ndim) :: ldim
integer(kind=4) :: lncid, vtype, ier, idv, lsize, lndim
real, intent(in) :: xmin, xmax
real(kind=4) scalef, addoff
character(len=*), intent(in) :: name
character(len=*), intent(in) :: lname
character(len=*), intent(in) :: units
      
if (itype==1) then
  vtype = nf_short
else
#ifdef r8i8
  vtype = nf_double
#else
  vtype = nf_float
#endif
end if

lncid=cdfid
lndim=ndim
ldim=dim
ier = nf_def_var(lncid, name, vtype, lndim, ldim, idv)
if (ier/=0) then
  write(6,*) "ERROR: Cannot define ",trim(name)
  write(6,*) nf_strerror(ier)
  call ccmpi_abort(-1)
end if
lsize=len_trim(lname)
ier = nf_put_att_text(lncid,idv,'long_name',lsize,lname)
if (ier/=0) then
  write(6,*) "ERROR: Cannot define long_name for ",trim(name)
  write(6,*) nf_strerror(ier)
  call ccmpi_abort(-1)
end if
lsize=len_trim(units)
if(lsize>0)then
  ier = nf_put_att_text(lncid,idv,'units',lsize,units)
  if (ier/=0) then
    write(6,*) "ERROR: Cannot define units for ",trim(name)
    write(6,*) nf_strerror(ier)
    call ccmpi_abort(-1)
  end if
endif
if (vtype == nf_short) then
  ier = nf_put_att_int2(lncid,idv,'valid_min',nf_short,1,minv)
  if (ier/=0) then
    write(6,*) "ERROR: Cannot define valid_min for ",trim(name)
    write(6,*) nf_strerror(ier)
    call ccmpi_abort(-1)
  end if
  ier = nf_put_att_int2(lncid,idv,'valid_max',nf_short,1,maxv)
  if (ier/=0) then
    write(6,*) "ERROR: Cannot define valid_max for ",trim(name)
    write(6,*) nf_strerror(ier)
    call ccmpi_abort(-1)
  end if
  ier = nf_put_att_int2(lncid,idv,'missing_value',nf_short,1,missval)
  if (ier/=0) then
    write(6,*) "ERROR: Cannot define missing_value for ",trim(name)
    write(6,*) nf_strerror(ier)
    call ccmpi_abort(-1)
  end if
  scalef=(xmax-xmin)/(real(maxv)-real(minv))
  addoff=xmin-scalef*minv
  ier = nf_put_att_real(lncid,idv,'add_offset',nf_float,1,addoff)
  if (ier/=0) then
    write(6,*) "ERROR: Cannot define add_offset for ",trim(name)
    write(6,*) nf_strerror(ier)
    call ccmpi_abort(-1)
  end if
  ier = nf_put_att_real(lncid,idv,'scale_factor',nf_float,1,scalef)
  if (ier/=0) then
    write(6,*) "ERROR: Cannot define scale_factor for ",trim(name)
    write(6,*) nf_strerror(ier)
    call ccmpi_abort(-1)
  end if
else
  ier = nf_put_att_real(lncid,idv,'missing_value',nf_float,1,nf_fill_float)
  if (ier/=0) then
    write(6,*) "ERROR: Cannot define missing_value for ",trim(name)
    write(6,*) nf_strerror(ier)
    call ccmpi_abort(-1)
  end if
endif
ier = nf_put_att_text(lncid,idv,'FORTRAN_format',5,'G11.4')
  if (ier/=0) then
    write(6,*) "ERROR: Cannot define FORTRAN_format for ",trim(name)
    write(6,*) nf_strerror(ier)
    call ccmpi_abort(-1)
  end if
if(daily>0)then
  ier = nf_put_att_text(lncid,idv,'valid_time',5,'daily')
  if (ier/=0) then
    write(6,*) "ERROR: Cannot define valid_time for ",trim(name)
    write(6,*) nf_strerror(ier)
    call ccmpi_abort(-1)
  end if
endif
      
return
end subroutine attrib

!--------------------------------------------------------------------
! 3D NETCDF WRITE ARRAY ROUTINES
subroutine histwrt3(var,sname,idnc,iarch,local,lwrite)

use cc_mpi              ! CC MPI routines

implicit none

include 'newmpar.h'     ! Grid parameters
include 'netcdf.inc'    ! Netcdf parameters

integer, intent(in) :: idnc, iarch
real, dimension(ifull), intent(in) :: var
real, dimension(ifull,1) :: wvar
character(len=*), intent(in) :: sname
logical, intent(in) :: local, lwrite

wvar(:,1)=var
if (.not.lwrite) then
  wvar=real(nf_fill_float)
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
      
implicit none
      
include 'newmpar.h'      ! Grid parameters
include 'netcdf.inc'     ! Netcdf parameters
include 'parm.h'         ! Model configuration
      
integer, intent(in) :: idnc,iarch,istep
integer(kind=4) mid, vtype, ier, lncid
integer(kind=4), dimension(3) :: start, ncount
integer(kind=2), dimension(ifull,istep) :: ipack
real, dimension(ifull,istep), intent(in) :: var
real(kind=4) addoff, scale_f
character(len=*), intent(in) :: sname

start = (/ 1, 1, iarch /)
ncount = (/ il, jl, istep /)

lncid=idnc
ier=nf_inq_varid(lncid,sname,mid)
if (ier/=nf_noerr) then
  write(6,*) "ERROR: Cannot find varid ",trim(sname)
  write(6,*) nf_strerror(ier)
  call ccmpi_abort(-1)
end if

! Check variable type
ier = nf_inq_vartype(lncid, mid, vtype)
if(vtype == ncshort)then
  if (all(var>9.8E36)) then
    ipack=missval
  else
    ier=nf_get_att_real(lncid,mid,'add_offset',addoff)
    ier=nf_get_att_real(lncid,mid,'scale_factor',scale_f)
    ipack=max(min(nint((var-addoff)/scale_f),maxv),minv)
  end if
  ier=nf_put_vara_int2(lncid,mid,start,ncount,ipack)
else
#ifdef r8i8
  ier=nf_put_vara_double(lncid,mid,start,ncount,var)
#else
  ier=nf_put_vara_real(lncid,mid,start,ncount,var)
#endif
end if
if (ier/=0) then
  write(6,*) "ERROR: Cannot write ",sname
  write(6,*) nf_strerror(ier)
  call ccmpi_abort(-1)
end if

if(mod(ktau,nmaxpr)==0.and.myid==0)then
  if (any(var==nf_fill_float)) then
    write(6,'("histwrt3 ",a7,i4,a7)') sname,iarch,"missing"
  else
    write(6,'("histwrt3 ",a7,i4)') sname,iarch
  end if
end if

return
end subroutine fw3l
      
subroutine fw3a(var,sname,idnc,iarch,istep)

use cc_mpi               ! CC MPI routines

implicit none

include 'newmpar.h'      ! Grid parameters
include 'netcdf.inc'     ! Netcdf parameters
include 'parm.h'         ! Model configuration

integer, intent(in) :: idnc, iarch, istep
integer(kind=4) mid, vtype, ier, lncid
integer imn, imx, jmn, jmx, iq
integer(kind=4), dimension(3) :: start, ncount
integer(kind=2), dimension(ifull_g,istep) :: ipack
real, dimension(ifull,istep), intent(in) :: var
real, dimension(ifull_g,istep) :: globvar
real(kind=4) addoff, scale_f
real varn, varx
character(len=*), intent(in) :: sname
      
call ccmpi_gather(var, globvar)
start(1) = 1
start(2) = 1
start(3) = iarch
ncount(1) = il_g
ncount(2) = jl_g
ncount(3) = istep

!     find variable index
lncid=idnc
ier=nf_inq_varid(lncid,sname,mid)
if (ier/=nf_noerr) then
  write(6,*) "ERROR: Cannot find varid ",trim(sname)
  write(6,*) nf_strerror(ier)
  call ccmpi_abort(-1)
end if

!     Check variable type
ier = nf_inq_vartype(lncid, mid, vtype)
if (vtype == nf_short) then
  if (all(globvar>9.8e36)) then
    ipack=missval
  else
    ier=nf_get_att_real(lncid,mid,'add_offset',addoff)
    ier=nf_get_att_real(lncid,mid,'scale_factor',scale_f)
    ipack=max(min(nint((globvar-addoff)/scale_f),maxv),minv)
  endif
  ier=nf_put_vara_int2(lncid,mid,start,ncount,ipack)
else
#ifdef r8i8
  ier=nf_put_vara_double(lncid,mid,start,ncount,var)
#else
  ier=nf_put_vara_real(lncid,mid,start,ncount,var)
#endif
endif
if (ier/=0) then
  write(6,*) "ERROR: Cannot write ",sname
  write(6,*) nf_strerror(ier)
  call ccmpi_abort(-1)
end if

if(mod(ktau,nmaxpr)==0)then
  if (any(globvar==nf_fill_float)) then
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

implicit none

include 'newmpar.h'     ! Grid parameters
include 'netcdf.inc'    ! Netcdf parameters

integer, intent(in) :: idnc, iarch
real, dimension(ifull,kl), intent(in) :: var
real, dimension(ifull,kl) :: wvar
character(len=*), intent(in) :: sname
logical, intent(in) :: local,lwrite

if (lwrite) then
  wvar=var
else
  wvar=real(nf_fill_float)
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

implicit none

include 'newmpar.h'      ! Grid parameters
include 'netcdf.inc'     ! Netcdf parameters
include 'parm.h'         ! Model configuration

integer, intent(in) :: idnc, iarch
integer(kind=4) lncid, mid, vtype, ier
integer(kind=4), dimension(4) :: start, ncount
integer(kind=2), dimension(ifull,kl) :: ipack
real, dimension(ifull,kl), intent(in) :: var
real(kind=4) addoff, scale_f
character(len=*), intent(in) :: sname

start = (/ 1, 1, 1, iarch /)
ncount = (/ il, jl, kl, 1 /)

lncid=idnc
ier=nf_inq_varid(lncid,sname,mid)
if (ier/=nf_noerr) then
  write(6,*) "ERROR: Cannot find varid ",trim(sname)
  write(6,*) nf_strerror(ier)
  call ccmpi_abort(-1)
end if

!     Check variable type
ier = nf_inq_vartype(lncid, mid, vtype)
if(vtype == nf_short)then
  if(all(var>9.8e36))then
    ipack=missval
  else
    ier=nf_get_att_real(lncid,mid,'add_offset',addoff)
    ier=nf_get_att_real(lncid,mid,'scale_factor',scale_f)
    ipack=max(min(nint((var-addoff)/scale_f),maxv),minv)
  endif
  ier=nf_put_vara_int2(lncid,mid,start,ncount,ipack)
else
#ifdef r8i8
  ier=nf_put_vara_double(lncid,mid,start,ncount,var)
#else
  ier=nf_put_vara_real(lncid,mid,start,ncount,var)
#endif
endif
if (ier/=0) then
  write(6,*) "ERROR: Cannot write ",sname
  write(6,*) nf_strerror(ier)
  call ccmpi_abort(-1)
end if

if(mod(ktau,nmaxpr)==0.and.myid==0)then
  if (any(var==nf_fill_float)) then
    write(6,'("histwrt4 ",a7,i4,a7)') sname,iarch,"missing"
  else
    write(6,'("histwrt4 ",a7,i4)') sname,iarch
  end if
endif

return
end subroutine hw4l      

subroutine hw4a(var,sname,idnc,iarch)

use cc_mpi              ! CC MPI routines

implicit none

include 'newmpar.h'     ! Grid parameters
include 'netcdf.inc'    ! Netcdf parameters
include 'parm.h'        ! Model configuration

integer, intent(in) :: idnc, iarch
integer(kind=4) lncid, mid, vtype, ier
integer imx, jmx, kmx, iq
integer(kind=4), dimension(4) :: start, ncount
integer, dimension(2) :: max_result
integer(kind=2), dimension(ifull_g,kl) :: ipack
real(kind=4) addoff, scale_f
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
lncid=idnc
ier=nf_inq_varid(lncid,sname,mid)
if (ier/=nf_noerr) then
  write(6,*) "ERROR: Cannot find varid ",trim(sname)
  write(6,*) nf_strerror(ier)
  call ccmpi_abort(-1)
end if

!     Check variable type
ier = nf_inq_vartype(lncid, mid, vtype)
if(vtype == nf_short)then
  if(all(globvar>9.8e36))then
    ipack=missval
  else
    ier=nf_get_att_real(lncid,mid,'add_offset',addoff)
    ier=nf_get_att_real(lncid,mid,'scale_factor',scale_f)
    ipack=max(min(nint((globvar-addoff)/scale_f),maxv),minv)
  endif
  ier=nf_put_vara_int2(lncid,mid,start,ncount,ipack)
else
#ifdef r8i8
  ier=nf_put_vara_double(lncid,mid,start,ncount,var)
#else
  ier=nf_put_vara_real(lncid,mid,start,ncount,var)
#endif
endif
if (ier/=0) then
  write(6,*) "ERROR: Cannot write ",sname
  write(6,*) nf_strerror(ier)
  call ccmpi_abort(-1)
end if

if(mod(ktau,nmaxpr)==0)then
  if (any(globvar==nf_fill_float)) then
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

implicit none

include 'netcdf.inc'

integer, intent(out) :: ncid,status
integer(kind=4) ncstatus, lncid
character(len=*), intent(in) :: fname

ncstatus = nf_open(fname,nf_nowrite,lncid)
ncid=lncid
status=ncstatus

return
end subroutine ccnf_open

subroutine ccnf_create(fname,ncid)

use cc_mpi

implicit none

include 'netcdf.inc'

integer, intent(out) :: ncid
integer(kind=4) ncstatus, lncid
character(len=*), intent(in) :: fname

#ifdef usenc3
  ncstatus = nf_create(fname,nf_clobber,lncid)
#else
  ncstatus = nf_create(fname,nf_netcdf4,lncid)
#endif
if (ncstatus/=nf_noerr) then
  write(6,*) nf_strerror(ncstatus)
  call ccmpi_abort(-1)
end if
ncid=lncid

return
end subroutine ccnf_create

subroutine ccnf_close(ncid)

use cc_mpi

implicit none

include 'netcdf.inc'

integer, intent(in) :: ncid
integer(kind=4) lncid,ncstatus

lncid=ncid
ncstatus = nf_close(lncid)
if (ncstatus/=nf_noerr) then
  write(6,*) nf_strerror(ncstatus)
  call ccmpi_abort(-1)
end if

return
end subroutine ccnf_close

subroutine ccnf_nofill(ncid)

use cc_mpi

implicit none

include 'netcdf.inc'

integer, intent(in) :: ncid
integer(kind=4) ncstatus,lncid,lomode

lncid=ncid
ncstatus=nf_set_fill(lncid,nf_nofill,lomode)
if (ncstatus/=nf_noerr) then
  write(6,*) "ERROR: nofill: ",nf_strerror(ncstatus)
  call ccmpi_abort(-1)
end if

return
end subroutine ccnf_nofill

subroutine ccnf_inq_dimid(ncid,dname,did,tst)

use cc_mpi

implicit none

include 'netcdf.inc'

integer, intent(in) :: ncid
integer, intent(out) :: did
integer(kind=4) lncid, ldid, ncstatus
logical, intent(out) :: tst
character(len=*), intent(in) :: dname

lncid=ncid
ncstatus=nf_inq_dimid(lncid,dname,ldid)
did=ldid
tst=(ncstatus/=nf_noerr)

return
end subroutine ccnf_inq_dimid

subroutine ccnf_inq_dimlen(ncid,dname,dlen,failok)

use cc_mpi

implicit none

include 'netcdf.inc'

integer, intent(in) :: ncid
integer, intent(inout) :: dlen
integer(kind=4) lncid,ldlen,ncstatus,ldid
character(len=*), intent(in) :: dname
logical, intent(in), optional :: failok
logical ftest

ftest=.false.
if (present(failok)) ftest=failok

lncid=ncid
ncstatus=nf_inq_dimid(lncid,dname,ldid)
if (ncstatus/=nf_noerr) then
  if (ftest) return
  write(6,*) nf_strerror(ncstatus)
  call ccmpi_abort(-1)
end if
ldlen=dlen
ncstatus=nf_inq_dimlen(lncid,ldid,ldlen)
if (ncstatus/=nf_noerr) then
  if (ftest) return
  write(6,*) nf_strerror(ncstatus)
  call ccmpi_abort(-1)
end if
dlen=ldlen

return
end subroutine ccnf_inq_dimlen

subroutine ccnf_inq_varid(ncid,vname,vid,tst)

use cc_mpi

implicit none

include 'netcdf.inc'

integer, intent(in) :: ncid
integer, intent(out) :: vid
integer(kind=4) ncstatus,lncid,lvid
character(len=*), intent(in) :: vname
logical, intent(out) :: tst

lncid=ncid
ncstatus = nf_inq_varid(lncid,vname,lvid)
tst=(ncstatus/=nf_noerr)
vid=lvid

return
end subroutine ccnf_inq_varid

subroutine ccnf_def_dim(ncid,dname,nsize,did)

use cc_mpi

implicit none

include 'netcdf.inc'

integer, intent(in) :: ncid,nsize
integer, intent(out) :: did
integer(kind=4) ncstatus,lncid,lsize,ldid
character(len=*), intent(in) :: dname

lncid=ncid
lsize=nsize
ncstatus=nf_def_dim(lncid,dname,lsize,ldid)
if (ncstatus/=nf_noerr) then
  write(6,*) nf_strerror(ncstatus)
  call ccmpi_abort(-1)
end if
did=ldid

return
end subroutine ccnf_def_dim

subroutine ccnf_def_dimu(ncid,dname,did)

use cc_mpi

implicit none

include 'netcdf.inc'

integer, intent(in) :: ncid
integer, intent(out) :: did
integer(kind=4) ncstatus,lncid,lsize,ldid
character(len=*), intent(in) :: dname

lncid=ncid
lsize=nf_unlimited
ncstatus=nf_def_dim(lncid,dname,lsize,ldid)
if (ncstatus/=nf_noerr) then
  write(6,*) nf_strerror(ncstatus)
  call ccmpi_abort(-1)
end if
did=ldid

return
end subroutine ccnf_def_dimu

subroutine ccnf_def_var(ncid,vname,vtype,vndim,dims,vid)

use cc_mpi

implicit none

include 'netcdf.inc'

integer, intent(in) :: ncid,vndim
integer, intent(out) :: vid
integer, dimension(vndim), intent(in) :: dims
integer(kind=4) ncstatus,lncid,lndim,lvid,ltype
integer(kind=4), dimension(vndim) :: ldims
character(len=*), intent(in) :: vname
character(len=*), intent(in) :: vtype

lncid=ncid
lndim=vndim
ldims=dims
select case(vtype)
  case('int')
    ltype=nf_int
  case('float')
    ltype=nf_float
  case('double')
    ltype=nf_double
  case('char')
    ltype=nf_char
  case default
    write(6,*) "ERROR: Unknown option for ccnf_def_var ",vtype
    call ccmpi_abort(-1)
end select

ncstatus = nf_def_var(lncid,vname,ltype,lndim,ldims,lvid)
if (ncstatus/=nf_noerr) then
  write(6,*) nf_strerror(ncstatus)
  call ccmpi_abort(-1)
end if
vid=lvid

return
end subroutine ccnf_def_var

subroutine ccnf_def_var0(ncid,vname,vtype,vid)

use cc_mpi

implicit none

include 'netcdf.inc'

integer, intent(in) :: ncid
integer, intent(out) :: vid
integer(kind=4) ncstatus,lncid,lvid,ltype
character(len=*), intent(in) :: vname
character(len=*), intent(in) :: vtype

lncid=ncid
select case(vtype)
  case('int')
    ltype=nf_int
  case('float')
    ltype=nf_float
  case('double')
    ltype=nf_double
  case('char')
    ltype=nf_char
  case default
    write(6,*) "ERROR: Unknown option for ccnf_def_var ",vtype
    call ccmpi_abort(-1)
end select

ncstatus = nf_def_var(lncid,vname,ltype,0,1,lvid)
if (ncstatus/=nf_noerr) then
  write(6,*) nf_strerror(ncstatus)
  call ccmpi_abort(-1)
end if
vid=lvid

return
end subroutine ccnf_def_var0

subroutine ccnf_enddef(ncid)

use cc_mpi

implicit none

include 'netcdf.inc'

integer, intent(in) :: ncid
integer(kind=4) ncstatus,lncid

lncid=ncid
ncstatus=nf_enddef(lncid)
if (ncstatus/=nf_noerr) then
  write(6,*) nf_strerror(ncstatus)
  call ccmpi_abort(-1)
end if

return
end subroutine ccnf_enddef

subroutine ccnf_redef(ncid)

use cc_mpi

implicit none

include 'netcdf.inc'

integer, intent(in) :: ncid
integer(kind=4) ncstatus,lncid

lncid=ncid
ncstatus=nf_redef(ncid)
if (ncstatus/=nf_noerr) then
  write(6,*) nf_strerror(ncstatus)
  call ccmpi_abort(-1)
end if

return
end subroutine ccnf_redef

subroutine ccnf_sync(ncid)

use cc_mpi

implicit none

include 'netcdf.inc'

integer, intent(in) :: ncid
integer(kind=4) ncstatus,lncid

#ifdef outsync
lncid=ncid
ncstatus = nf_sync(lncid)
if (ncstatus/=nf_noerr) then
  write(6,*) nf_strerror(ncstatus)
  call ccmpi_abort(-1)
end if
#endif

return
end subroutine ccnf_sync

subroutine ccnf_get_var_real(ncid,vname,vdat)

use cc_mpi

implicit none

include 'netcdf.inc'

integer, intent(in) :: ncid
integer(kind=4) lncid,lvid,ncstatus
real, dimension(:), intent(out) :: vdat
real(kind=4), dimension(size(vdat)) :: ldat
character(len=*), intent(in) :: vname

lncid=ncid
ncstatus = nf_inq_varid(lncid,vname,lvid)
if (ncstatus/=nf_noerr) then
  write(6,*) nf_strerror(ncstatus)
  call ccmpi_abort(-1)
end if
ncstatus = nf_get_var_real(lncid,lvid,ldat)
if (ncstatus/=nf_noerr) then
  write(6,*) nf_strerror(ncstatus)
  call ccmpi_abort(-1)
end if
vdat=ldat

return
end subroutine ccnf_get_var_real

subroutine ccnf_get_var_int(ncid,vname,vdat)

use cc_mpi

implicit none

include 'netcdf.inc'

integer, intent(in) :: ncid
integer(kind=4) lncid,lvid,ncstatus
integer, dimension(:), intent(out) :: vdat
integer(kind=4), dimension(size(vdat)) :: ldat
character(len=*), intent(in) :: vname

lncid=ncid
ncstatus = nf_inq_varid(lncid,vname,lvid)
if (ncstatus/=nf_noerr) then
  write(6,*) nf_strerror(ncstatus)
  call ccmpi_abort(-1)
end if
ncstatus = nf_get_var_int(lncid,lvid,ldat)
if (ncstatus/=nf_noerr) then
  write(6,*) nf_strerror(ncstatus)
  call ccmpi_abort(-1)
end if
vdat=ldat

return
end subroutine ccnf_get_var_int

subroutine ccnf_get_var1_real(ncid,vid,start,vdat)

use cc_mpi

implicit none

include 'netcdf.inc'

integer, intent(in) :: ncid,vid,start
integer(kind=4) lncid,lvid,ncstatus,lstart
real, intent(out) :: vdat
real(kind=4) ldat

lncid=ncid
lvid=vid
lstart=start
ncstatus = nf_get_var1_real(lncid,lvid,lstart,ldat)
if (ncstatus/=nf_noerr) then
  write(6,*) nf_strerror(ncstatus)
  call ccmpi_abort(-1)
end if
vdat=ldat

return
end subroutine ccnf_get_var1_real

subroutine ccnf_get_var1_int(ncid,vid,start,vdat)

use cc_mpi

implicit none

include 'netcdf.inc'

integer, intent(in) :: ncid,vid,start
integer(kind=4) lncid,lvid,ncstatus,lstart
integer, intent(out) :: vdat
integer(kind=4) ldat

lncid=ncid
lvid=vid
lstart=start
ncstatus = nf_get_var1_int(lncid,lvid,lstart,ldat)
if (ncstatus/=nf_noerr) then
  write(6,*) nf_strerror(ncstatus)
  call ccmpi_abort(-1)
end if
vdat=ldat

return
end subroutine ccnf_get_var1_int

subroutine ccnf_get_vara_real2r(ncid,vid,start,count,vdat)

use cc_mpi

implicit none

include 'netcdf.inc'

integer, intent(in) :: ncid,vid
integer, dimension(:) :: start,count
integer(kind=4) ncstatus,lncid,lvid
integer(kind=4), dimension(size(start)) :: lstart,lcount
real, dimension(:), intent(out) :: vdat
real(kind=4), dimension(size(vdat)) :: ldat

lncid=ncid
lvid=vid
lstart=start
lcount=count
ncstatus=nf_get_vara_real(lncid,lvid,lstart,lcount,ldat)
if (ncstatus/=nf_noerr) then
  write(6,*) nf_strerror(ncstatus)
  call ccmpi_abort(-1)
end if
vdat=ldat

return
end subroutine ccnf_get_vara_real2r 

subroutine ccnf_get_vara_real3r(ncid,vid,start,count,vdat)

use cc_mpi

implicit none

include 'netcdf.inc'

integer, intent(in) :: ncid,vid
integer, dimension(:) :: start,count
integer(kind=4) ncstatus,lncid,lvid
integer(kind=4), dimension(size(start)) :: lstart,lcount
real, dimension(:,:), intent(out) :: vdat
real(kind=4), dimension(size(vdat,1),size(vdat,2)) :: ldat

lncid=ncid
lvid=vid
lstart=start
lcount=count
ncstatus=nf_get_vara_real(lncid,lvid,lstart,lcount,ldat)
if (ncstatus/=nf_noerr) then
  write(6,*) nf_strerror(ncstatus)
  call ccmpi_abort(-1)
end if
vdat=ldat

return
end subroutine ccnf_get_vara_real3r

subroutine ccnf_get_vara_real4r(ncid,vid,start,count,vdat)

use cc_mpi

implicit none

include 'netcdf.inc'

integer, intent(in) :: ncid,vid
integer, dimension(:) :: start,count
integer(kind=4) ncstatus,lncid,lvid
integer(kind=4), dimension(size(start)) :: lstart,lcount
real, dimension(:,:,:), intent(out) :: vdat
real(kind=4), dimension(size(vdat,1),size(vdat,2),size(vdat,3)) :: ldat

lncid=ncid
lvid=vid
lstart=start
lcount=count
ncstatus=nf_get_vara_real(lncid,lvid,lstart,lcount,ldat)
if (ncstatus/=nf_noerr) then
  write(6,*) nf_strerror(ncstatus)
  call ccmpi_abort(-1)
end if
vdat=ldat

return
end subroutine ccnf_get_vara_real4r

subroutine ccnf_get_vara_int2i(ncid,vid,start,count,vdat)

use cc_mpi

implicit none

include 'netcdf.inc'

integer, intent(in) :: ncid,vid
integer, dimension(:) :: start,count
integer(kind=4) ncstatus,lncid,lvid
integer(kind=4), dimension(size(start)) :: lstart,lcount
integer, dimension(:), intent(out) :: vdat
integer(kind=4), dimension(size(vdat)) :: ldat

lncid=ncid
lvid=vid
lstart=start
lcount=count
ncstatus=nf_get_vara_int(lncid,lvid,lstart,lcount,ldat)
if (ncstatus/=nf_noerr) then
  write(6,*) nf_strerror(ncstatus)
  call ccmpi_abort(-1)
end if
vdat=ldat

return
end subroutine ccnf_get_vara_int2i

subroutine ccnf_get_vara_double4d(ncid,vid,start,count,vdat)

use cc_mpi

implicit none

include 'netcdf.inc'

integer, intent(in) :: ncid,vid
integer, dimension(:) :: start,count
integer(kind=4) ncstatus,lncid,lvid
integer(kind=4), dimension(size(start)) :: lstart,lcount
double precision, dimension(:,:,:), intent(out) :: vdat

lncid=ncid
lvid=vid
lstart=start
lcount=count
ncstatus=nf_get_vara_double(lncid,lvid,lstart,lcount,vdat)
if (ncstatus/=nf_noerr) then
  write(6,*) nf_strerror(ncstatus)
  call ccmpi_abort(-1)
end if

return
end subroutine ccnf_get_vara_double4d

subroutine ccnf_get_att_text(ncid,vid,aname,atext)

use cc_mpi

implicit none

include 'netcdf.inc'

integer, intent(in) :: ncid,vid
integer(kind=4) ncstatus,lncid,lvid
character(len=*), intent(in) :: aname
character(len=*), intent(out) :: atext

lncid=ncid
lvid=vid
atext=''
ncstatus = nf_get_att_text(lncid,lvid,aname,atext)
if (ncstatus/=nf_noerr) then
  write(6,*) nf_strerror(ncstatus)
  call ccmpi_abort(-1)
end if

return
end subroutine ccnf_get_att_text

subroutine ccnf_get_att_real(ncid,vid,aname,vdat)

use cc_mpi

implicit none

include 'netcdf.inc'

integer, intent(in) :: ncid,vid
integer(kind=4) ncstatus,lncid,lvid
character(len=*), intent(in) :: aname
real, intent(out) :: vdat
real(kind=4) ldat

lncid=ncid
lvid=vid
ncstatus = nf_get_att_real(lncid,lvid,aname,ldat)
if (ncstatus/=nf_noerr) then
  write(6,*) nf_strerror(ncstatus)
  call ccmpi_abort(-1)
end if
vdat=ldat

return
end subroutine ccnf_get_att_real

subroutine ccnf_get_att_realg1r(ncid,aname,vdat)

use cc_mpi

implicit none

include 'netcdf.inc'

integer, intent(in) :: ncid
integer(kind=4) ncstatus,lncid
character(len=*), intent(in) :: aname
real, intent(out) :: vdat
real(kind=4) ldat

lncid=ncid
ncstatus = nf_get_att_real(lncid,nf_global,aname,ldat)
if (ncstatus/=nf_noerr) then
  write(6,*) nf_strerror(ncstatus)
  call ccmpi_abort(-1)
end if
vdat=ldat

return
end subroutine ccnf_get_att_realg1r

subroutine ccnf_get_att_realg2r(ncid,aname,vdat)

use cc_mpi

implicit none

include 'netcdf.inc'

integer, intent(in) :: ncid
integer(kind=4) ncstatus,lncid
character(len=*), intent(in) :: aname
real, dimension(:), intent(out) :: vdat
real(kind=4), dimension(size(vdat)) :: ldat

lncid=ncid
ncstatus = nf_get_att_real(lncid,nf_global,aname,ldat)
if (ncstatus/=nf_noerr) then
  write(6,*) nf_strerror(ncstatus)
  call ccmpi_abort(-1)
end if
vdat=ldat

return
end subroutine ccnf_get_att_realg2r

subroutine ccnf_get_att_intg1i(ncid,aname,vdat)

use cc_mpi

implicit none

include 'netcdf.inc'

integer, intent(in) :: ncid
integer(kind=4) ncstatus,lncid
character(len=*), intent(in) :: aname
integer, intent(out) :: vdat
integer(kind=4) ldat

lncid=ncid
ncstatus = nf_get_att_int(lncid,nf_global,aname,ldat)
if (ncstatus/=nf_noerr) then
  write(6,*) nf_strerror(ncstatus)
  call ccmpi_abort(-1)
end if
vdat=ldat

return
end subroutine ccnf_get_att_intg1i

subroutine ccnf_get_att_intg2i(ncid,aname,vdat)

use cc_mpi

implicit none

include 'netcdf.inc'

integer, intent(in) :: ncid
integer(kind=4) ncstatus,lncid
character(len=*), intent(in) :: aname
integer, dimension(:), intent(out) :: vdat
integer(kind=4), dimension(size(vdat)) :: ldat

lncid=ncid
ncstatus = nf_get_att_int(lncid,nf_global,aname,ldat)
if (ncstatus/=nf_noerr) then
  write(6,*) nf_strerror(ncstatus)
  call ccmpi_abort(-1)
end if
vdat=ldat

return
end subroutine ccnf_get_att_intg2i

subroutine ccnf_read(fname,vname,vdat)

use cc_mpi

implicit none

include 'newmpar.h'
include 'netcdf.inc'

integer(kind=4) ncstatus,lncid,lvid
real, dimension(ifull), intent(out) :: vdat
real, dimension(ifull_g) :: vdat_g
real(kind=4), dimension(ifull_g) :: ldat_g
character(len=*), intent(in) :: fname
character(len=*), intent(in) :: vname

if (myid==0) then
  ncstatus = nf_open(fname,nf_nowrite,lncid)
  if (ncstatus/=nf_noerr) then
    write(6,*) nf_strerror(ncstatus)
    call ccmpi_abort(-1)
  end if
  ncstatus = nf_inq_varid(lncid,vname,lvid)
  if (ncstatus/=nf_noerr) then
    write(6,*) nf_strerror(ncstatus)
    call ccmpi_abort(-1)
  end if
  ncstatus = nf_get_var_real(lncid,lvid,ldat_g)
  if (ncstatus/=nf_noerr) then
    write(6,*) nf_strerror(ncstatus)
    call ccmpi_abort(-1)
  end if
  ncstatus = nf_close(lncid)
  if (ncstatus/=nf_noerr) then
    write(6,*) nf_strerror(ncstatus)
    call ccmpi_abort(-1)
  end if
  vdat_g=ldat_g
  call ccmpi_distribute(vdat,vdat_g)
else
  call ccmpi_distribute(vdat)
end if

return
end subroutine ccnf_read

subroutine ccnf_put_var_text2r(ncid,vid,vtxt)

use cc_mpi

implicit none

include 'netcdf.inc'

integer, intent(in) :: ncid,vid
integer(kind=4) ncstatus,lncid,lvid
character(len=*), dimension(:), intent(in) :: vtxt

lncid=ncid
lvid=vid
ncstatus = nf_put_var_text(lncid,lvid,vtxt)
if (ncstatus/=nf_noerr) then
  write(6,*) nf_strerror(ncstatus)
  call ccmpi_abort(-1)
end if

return
end subroutine ccnf_put_var_text2r

subroutine ccnf_put_var_int2i(ncid,vid,vdat)

use cc_mpi

implicit none

include 'netcdf.inc'

integer, intent(in) :: ncid,vid
integer(kind=4) ncstatus,lncid,lvid
integer, dimension(:), intent(in) :: vdat
integer(kind=4), dimension(size(vdat)) :: ldat

lncid=ncid
lvid=vid
ldat=vdat
ncstatus = nf_put_var_int(lncid,lvid,ldat)
if (ncstatus/=nf_noerr) then
  write(6,*) nf_strerror(ncstatus)
  call ccmpi_abort(-1)
end if

return
end subroutine ccnf_put_var_int2i

subroutine ccnf_put_var_int3i(ncid,vid,vdat)

use cc_mpi

implicit none

include 'netcdf.inc'

integer, intent(in) :: ncid,vid
integer(kind=4) ncstatus,lncid,lvid
integer, dimension(:,:), intent(in) :: vdat
integer(kind=4), dimension(size(vdat,1),size(vdat,2)) :: ldat

lncid=ncid
lvid=vid
ldat=vdat
ncstatus = nf_put_var_int(lncid,lvid,ldat)
if (ncstatus/=nf_noerr) then
  write(6,*) nf_strerror(ncstatus)
  call ccmpi_abort(-1)
end if

return
end subroutine ccnf_put_var_int3i

subroutine ccnf_put_var1_int(ncid,vid,start,vdat)

use cc_mpi

implicit none

include 'netcdf.inc'

integer, intent(in) :: ncid,vid,start
integer(kind=4) ncstatus,lncid,lvid,lstart
integer, intent(in) :: vdat
integer(kind=4) :: ldat

lncid=ncid
lvid=vid
lstart=start
ldat=vdat
ncstatus=nf_put_var1_int(lncid,lvid,lstart,ldat)
if (ncstatus/=nf_noerr) then
  write(6,*) nf_strerror(ncstatus)
  call ccmpi_abort(-1)
end if

return
end subroutine ccnf_put_var1_int

subroutine ccnf_put_var1_real(ncid,vid,start,vdat)

use cc_mpi

implicit none

include 'netcdf.inc'

integer, intent(in) :: ncid,vid,start
integer(kind=4) ncstatus,lncid,lvid,lstart
real, intent(in) :: vdat

lncid=ncid
lvid=vid
lstart=start
#ifdef i8r8
ncstatus=nf_put_var1_double(lncid,lvid,lstart,vdat)
#else
ncstatus=nf_put_var1_real(lncid,lvid,lstart,vdat)
#endif
if (ncstatus/=nf_noerr) then
  write(6,*) nf_strerror(ncstatus)
  call ccmpi_abort(-1)
end if

return
end subroutine ccnf_put_var1_real

subroutine ccnf_put_var1_double(ncid,vid,start,vdat)

use cc_mpi

implicit none

include 'netcdf.inc'

integer, intent(in) :: ncid,vid,start
integer(kind=4) ncstatus,lncid,lvid,lstart
double precision, intent(in) :: vdat

lncid=ncid
lvid=vid
lstart=start
ncstatus=nf_put_var1_double(lncid,lvid,lstart,vdat)
if (ncstatus/=nf_noerr) then
  write(6,*) nf_strerror(ncstatus)
  call ccmpi_abort(-1)
end if

return
end subroutine ccnf_put_var1_double

subroutine ccnf_put_vara_real2r(ncid,vid,start,count,vdat)

use cc_mpi

implicit none

include 'netcdf.inc'

integer, intent(in) :: ncid,vid
integer(kind=4) ncstatus,lncid,lvid
integer, dimension(:), intent(in) :: start,count
integer(kind=4), dimension(size(start)) :: lstart,lcount
real, dimension(:), intent(in) :: vdat

lncid=ncid
lvid=vid
lstart=start
lcount=count
#ifdef i8r8
ncstatus=nf_put_vara_double(lncid,lvid,lstart,lcount,vdat)
#else
ncstatus=nf_put_vara_real(lncid,lvid,lstart,lcount,vdat)
#endif
if (ncstatus/=nf_noerr) then
  write(6,*) nf_strerror(ncstatus)
  call ccmpi_abort(-1)
end if

return
end subroutine ccnf_put_vara_real2r

subroutine ccnf_put_vara_real3r(ncid,vid,start,count,vdat)

use cc_mpi

implicit none

include 'netcdf.inc'

integer, intent(in) :: ncid,vid
integer(kind=4) ncstatus,lncid,lvid
integer, dimension(:), intent(in) :: start,count
integer(kind=4), dimension(size(start)) :: lstart,lcount
real, dimension(:,:), intent(in) :: vdat

lncid=ncid
lvid=vid
lstart=start
lcount=count
#ifdef i8r8
ncstatus=nf_put_vara_double(lncid,lvid,lstart,lcount,vdat)
#else
ncstatus=nf_put_vara_real(lncid,lvid,lstart,lcount,vdat)
#endif
if (ncstatus/=nf_noerr) then
  write(6,*) nf_strerror(ncstatus)
  call ccmpi_abort(-1)
end if

return
end subroutine ccnf_put_vara_real3r

subroutine ccnf_put_vara_double2r(ncid,vid,start,count,vdat)

use cc_mpi

implicit none

include 'netcdf.inc'

integer, intent(in) :: ncid,vid
integer(kind=4) ncstatus,lncid,lvid
integer, dimension(:), intent(in) :: start,count
integer(kind=4), dimension(size(start)) :: lstart,lcount
double precision, dimension(:), intent(in) :: vdat

lncid=ncid
lvid=vid
lstart=start
lcount=count
ncstatus=nf_put_vara_double(lncid,lvid,lstart,lcount,vdat)
if (ncstatus/=nf_noerr) then
  write(6,*) nf_strerror(ncstatus)
  call ccmpi_abort(-1)
end if

return
end subroutine ccnf_put_vara_double2r

subroutine ccnf_put_vara_int2i(ncid,vid,start,count,vdat)

use cc_mpi

implicit none

include 'netcdf.inc'

integer, intent(in) :: ncid,vid
integer(kind=4) ncstatus,lncid,lvid
integer, dimension(:), intent(in) :: start,count
integer(kind=4), dimension(size(start)) :: lstart,lcount
integer, dimension(:), intent(in) :: vdat
integer(kind=4), dimension(size(vdat)) :: ldat

lncid=ncid
lvid=vid
lstart=start
lcount=count
ldat=vdat
ncstatus=nf_put_vara_int(lncid,lvid,lstart,lcount,ldat)
if (ncstatus/=nf_noerr) then
  write(6,*) nf_strerror(ncstatus)
  call ccmpi_abort(-1)
end if

return
end subroutine ccnf_put_vara_int2i

subroutine ccnf_put_att_text(ncid,vid,aname,asize,atext)

use cc_mpi

implicit none

include 'netcdf.inc'

integer, intent(in) :: ncid,vid,asize
integer(kind=4) ncstatus,lncid,lvid,lsize
character(len=*), intent(in) :: aname
character(len=asize), intent(in) :: atext

lncid=ncid
lvid=vid
lsize=asize
ncstatus=nf_put_att_text(lncid,lvid,aname,lsize,atext)
if (ncstatus/=nf_noerr) then
  write(6,*) nf_strerror(ncstatus)
  call ccmpi_abort(-1)
end if

return
end subroutine ccnf_put_att_text

subroutine ccnf_put_att_textg(ncid,aname,asize,atext)

use cc_mpi

implicit none

include 'netcdf.inc'

integer, intent(in) :: ncid,asize
integer(kind=4) ncstatus,lncid,lsize
character(len=*), intent(in) :: aname
character(len=asize), intent(in) :: atext

lncid=ncid
lsize=asize
ncstatus=nf_put_att_text(lncid,nf_global,aname,lsize,atext)
if (ncstatus/=nf_noerr) then
  write(6,*) nf_strerror(ncstatus)
  call ccmpi_abort(-1)
end if

return
end subroutine ccnf_put_att_textg

subroutine ccnf_put_att_intg(ncid,aname,vsize,vdat)

use cc_mpi

implicit none

include 'netcdf.inc'

integer, intent(in) :: ncid,vsize
integer(kind=4) ncstatus,lncid,lsize
integer, dimension(:), intent(in) :: vdat
integer(kind=4), dimension(size(vdat)) :: ldat
character(len=*), intent(in) :: aname

lncid=ncid
lsize=vsize
ldat=vdat
ncstatus=nf_put_att_int(lncid,nf_global,aname,nf_int,lsize,ldat)
if (ncstatus/=nf_noerr) then
  write(6,*) nf_strerror(ncstatus)
  call ccmpi_abort(-1)
end if

return
end subroutine ccnf_put_att_intg

subroutine ccnf_put_att_realg(ncid,aname,vsize,vdat)

use cc_mpi

implicit none

include 'netcdf.inc'

integer, intent(in) :: ncid,vsize
integer(kind=4) ncstatus,lncid,lsize
real, dimension(:), intent(in) :: vdat
character(len=*), intent(in) :: aname

lncid=ncid
lsize=vsize
#ifdef i8r8
ncstatus=nf_put_att_double(lncid,nf_global,aname,nf_float,lsize,vdat)
#else
ncstatus=nf_put_att_real(lncid,nf_global,aname,nf_float,lsize,vdat)
#endif
if (ncstatus/=nf_noerr) then
  write(6,*) nf_strerror(ncstatus)
  call ccmpi_abort(-1)
end if

return
end subroutine ccnf_put_att_realg

end module infile
