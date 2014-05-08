module infile
      
! This module contains routines for reading netcdf files,
! vertical interpolation and some calendar functions.
! This module also contains the interface to the netcdf
! library.
      
! This version of infile.f90 supports parallel localhist input
! files.  Multiple processors read as many input files as
! supplied in parallel and then send this data to all processors
! for interpolation. The code can also identify restart files,
! in which case no additional message passing is required.

#ifdef usenc3
use netcdf_m
#else
use netcdf
#endif

implicit none
            
private
public vertint,datefix,getzinp,ncmsg
public histopen,histclose,histrd1,histrd4s,pfall
public attrib,histwrt3,histwrt4,freqwrite,surfread
public ccnf_open,ccnf_create,ccnf_close,ccnf_sync,ccnf_enddef,ccnf_redef,ccnf_nofill,ccnf_inq_varid
public ccnf_inq_dimid,ccnf_inq_dimlen,ccnf_def_dim,ccnf_def_dimu,ccnf_def_var,ccnf_def_var0
public ccnf_get_var,ccnf_get_vara,ccnf_get_var1,ccnf_get_att,ccnf_get_attg,ccnf_read,ccnf_put_var
public ccnf_put_var1,ccnf_put_vara,ccnf_put_att,ccnf_put_attg

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
#ifndef i8r8
  module procedure ccnf_get_vara_double4d
#endif
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
#ifndef i8r8
  module procedure ccnf_put_vara_double2r
#endif
end interface ccnf_put_vara
interface ccnf_put_var
  module procedure ccnf_put_var_text2r
  module procedure ccnf_put_var_int2i, ccnf_put_var_int3i
end interface ccnf_put_var
interface ccnf_put_var1
  module procedure ccnf_put_var1_real, ccnf_put_var1_int
#ifndef i8r8
  module procedure ccnf_put_var1_double
#endif
end interface ccnf_put_var1

integer, dimension(:), allocatable, save :: pnoff
integer, dimension(:,:), allocatable, save :: pioff,pjoff
integer(kind=4), dimension(:), allocatable, save :: pncid
integer(kind=4), dimension(:), allocatable, save :: pncidold
integer, save :: mynproc,fnproc
integer, save :: pil_g,pjl_g,pil,pjl,pnpan
integer, save :: comm_ip,comm_ipold
logical, save :: ptest,pfall

integer(kind=2), parameter :: minv = -32500
integer(kind=2), parameter :: maxv =  32500
integer(kind=2), parameter :: missval = -32501
      
contains

!--------------------------------------------------------------
! Interface for reading 2D+time fields
subroutine histrd1(iarchi,ier,name,ik,var,ifull)
      
use cc_mpi

implicit none

include 'parm.h'
      
integer iarchi,ier,ik,ifull
character(len=*) name
real, dimension(:), intent(inout) :: var ! may be dummy argument from myid/=0

if ( ifull/=6*ik*ik .and. ptest ) then
  ! read local arrays without gather and distribute
  call hr1p(iarchi,ier,name,.true.,var)
  if ( ier==0 .and. mod(ktau,nmaxpr)==0 .and. myid==0 ) then
    write(6,'("done histrd1 ",a8,i4,i3)') name,ier,iarchi
  end if

else if ( myid==0 ) then
  ! split up processors to save memory.  No need to allocate global
  ! arrays on myid/=0.
  call hr1a(iarchi,ier,name,ik,var,ifull)

else if ( ifull/=6*ik*ik ) then
  ! read local arrays with gather and distribute
  call hr1p(iarchi,ier,name,.false.)
  call ccmpi_distribute(var)

else
  ! read global arrays for myid==0
  call hr1p(iarchi,ier,name,.false.)

end if
      
return
end subroutine histrd1    

!--------------------------------------------------------------------
! Gather 2D+time fields on myid==0
subroutine hr1a(iarchi,ier,name,ik,var,ifull)
      
use cc_mpi
      
implicit none
      
include 'parm.h'
      
integer, intent(in) :: iarchi, ik, ifull
integer, intent(out) :: ier
character(len=*), intent(in) :: name
real, dimension(:), intent(inout) :: var
real, dimension(6*ik*ik) :: globvar
real vmax, vmin

call hr1p(iarchi,ier,name,.false.,globvar)

if ( ier==0 .and. mod(ktau,nmaxpr)==0 ) then
  vmax = maxval(globvar)
  vmin = minval(globvar)
  write(6,'("done histrd1 ",a8,i4,i3,3e14.6)') name,ier,iarchi,vmin,vmax,globvar(id+(jd-1)*ik)
end if

if ( ifull==6*ik*ik ) then
  ! read global arrays for myid==0
  var(1:ifull)=globvar(:)
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

integer, intent(in) :: iarchi
integer, intent(out) :: ier
integer(kind=4), dimension(3) :: start,ncount
integer ipf,jpmax,iptst2,lcomm
integer(kind=4) :: idv,lier,ler
real, dimension(:), intent(inout), optional :: var
real, dimension(pil*pjl*pnpan) :: rvar
real addoff,sf
real(kind=4) :: laddoff,lsf
logical, intent(in) :: qtest
character(len=*), intent(in) :: name

start  = (/ 1, 1, iarchi /)
ncount = (/ pil, pjl*pnpan, 1 /)
ier=0

iptst2=mod(fnproc,nproc)
      
do ipf=0,mynproc-1

  rvar=0. ! default for missing field
  
  ! get variable idv
#ifdef usenc3
  lier=nf_inq_varid(pncid(ipf),name,idv)
  ier=lier
  if(lier/=nf_noerr)then
    if (myid==0.and.ipf==0) then
      write(6,*) '***absent field for ncid,name,idv,ier: ',pncid(0),name,idv,ier
    end if
  else
    ! obtain scaling factors and offsets from attributes
    ler=nf_get_att_real(pncid(ipf),idv,'add_offset',laddoff)
    addoff=laddoff
    if (ler/=nf_noerr) addoff=0.
    ler=nf_get_att_real(pncid(ipf),idv,'scale_factor',lsf)
    sf=lsf
    if (ler/=nf_noerr) sf=1.
#ifdef i8r8
    lier=nf_get_vara_double(pncid(ipf),idv,start,ncount,rvar)
#else
    lier=nf_get_vara_real(pncid(ipf),idv,start,ncount,rvar)
#endif
    ier=lier
    call ncmsg(name,ier)
    ! unpack compressed data
    rvar=rvar*sf+addoff
  end if ! ier
#else
  lier=nf90_inq_varid(pncid(ipf),name,idv)
  ier=lier
  if(lier/=nf90_noerr)then
    if (myid==0.and.ipf==0) then
      write(6,*) '***absent field for ncid,name,idv,ier: ',pncid(0),name,idv,ier
    end if
  else
    ! obtain scaling factors and offsets from attributes
    ler=nf90_get_att(pncid(ipf),idv,'add_offset',addoff)
    if (ler/=nf90_noerr) addoff=0.
    ler=nf90_get_att(pncid(ipf),idv,'scale_factor',sf)
    if (ler/=nf90_noerr) sf=1.
    lier=nf90_get_var(pncid(ipf),idv,rvar,start=start,count=ncount)
    ier=lier
    call ncmsg(name,ier)
    ! unpack compressed data
    rvar=rvar*sf+addoff
  end if ! ier
#endif
      
  if (qtest) then
    ! e.g., restart file
    var(1:pil*pjl*pnpan)=rvar
  else
    ! e.g., mesonest file
    if (ipf==mynproc-1.and.myid<iptst2) then
      jpmax=iptst2
      lcomm=comm_ip
    else
      jpmax=nproc    
      lcomm=comm_world    
    end if
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
subroutine histrd4s(iarchi,ier,name,ik,kk,var,ifull)
      
use cc_mpi
      
implicit none
      
include 'parm.h'
      
integer, intent(in) :: iarchi,ik,kk,ifull
integer, intent(out) :: ier
character(len=*), intent(in) :: name
real, dimension(:), intent(inout) :: var ! may be dummy argument from myid/=0
      
if (ifull/=6*ik*ik.and.ptest) then
  ! read local arrays without gather and distribute
  call hr4p(iarchi,ier,name,kk,.true.,var)
  if(ier==0.and.mod(ktau,nmaxpr)==0.and.myid==0)then
    write(6,'("done histrd4 ",a8,i3,i4,i3)') name,kk,ier,iarchi
  endif

else if (myid==0) then  
  ! split up processors to save memory.  No need to allocate global
  ! arrays on myid/=0.
  call hr4sa(iarchi,ier,name,ik,kk,var,ifull)

else if(ifull/=6*ik*ik)then
  ! read local arrays with gather and distribute
  call hr4p(iarchi,ier,name,kk,.false.)
  call ccmpi_distribute(var)

else
  ! read global arrays for myid==0
  call hr4p(iarchi,ier,name,kk,.false.)

end if

return
end subroutine histrd4s      

!--------------------------------------------------------------------
! Gather 3D+time fields on myid==0

subroutine hr4sa(iarchi,ier,name,ik,kk,var,ifull)
      
use cc_mpi
      
implicit none
      
include 'parm.h'
      
integer, intent(in) :: iarchi,ik,kk,ifull
integer, intent(out) :: ier
character(len=*), intent(in) :: name
real, dimension(:) :: var
real, dimension(6*ik*ik) :: globvar
real vmax,vmin
      
call hr4p(iarchi,ier,name,kk,.false.,globvar)     

if(ier==0.and.mod(ktau,nmaxpr)==0)then
  vmax = maxval(globvar)
  vmin = minval(globvar)
  write(6,'("done histrd4s ",a6,i3,i4,i3,3f12.4)') name,kk,ier,iarchi,vmin,vmax,globvar(id+(jd-1)*ik)
end if

! Have to return correct value of ier on all processes because it's 
! used for initialisation in calling routine
if(ifull==6*ik*ik)then
  ! read global arrays for myid==0
  var(1:ifull)=globvar(:)
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

integer, intent(in) :: iarchi,kk
integer, intent(out) :: ier
integer(kind=4), dimension(4) :: start,ncount
integer ipf,jpmax,iptst2,lcomm
integer(kind=4) :: ler,idv,lier
real, dimension(:), intent(inout), optional :: var
real, dimension(pil*pjl*pnpan) :: rvar
real, dimension(pil*pjl*pnpan*nproc) :: gvar
real addoff,sf
real(kind=4) :: laddoff,lsf
logical, intent(in) :: qtest
character(len=*), intent(in) :: name

start = (/ 1, 1, kk, iarchi /)
ncount = (/ pil, pjl*pnpan, 1, 1 /)
ier=0

iptst2=mod(fnproc,nproc)
      
do ipf=0,mynproc-1

  rvar=0. ! default value for missing field

  ! get variable idv
#ifdef usenc3
  lier=nf_inq_varid(pncid(ipf),name,idv)
  ier=lier
  if(lier/=nf_noerr)then
    if (myid==0.and.ipf==0.and.kk==1) then
      write(6,*) '***absent field for ncid,name,idv,ier: ',pncid(0),name,idv,ier
    end if
  else
    ! obtain scaling factors and offsets from attributes
    ler=nf_get_att_real(pncid(ipf),idv,'add_offset',laddoff)
    addoff=laddoff
    if (ler/=nf_noerr) addoff=0.
    ler=nf_get_att_real(pncid(ipf),idv,'scale_factor',lsf)
    sf=lsf
    if (ler/=nf_noerr) sf=1.
#ifdef i8r8
    lier=nf_get_vara_double(pncid(ipf),idv,start,ncount,rvar)
#else
    lier=nf_get_vara_real(pncid(ipf),idv,start,ncount,rvar)
#endif
    ier=lier
    call ncmsg(name,ier)
    ! unpack data
    rvar=rvar*sf+addoff
  end if ! ier
#else
  lier=nf90_inq_varid(pncid(ipf),name,idv)
  ier=lier
  if(lier/=nf90_noerr)then
    if (myid==0.and.ipf==0.and.kk==1) then
      write(6,*) '***absent field for ncid,name,idv,ier: ',pncid(0),name,idv,ier
    end if
  else
    ! obtain scaling factors and offsets from attributes
    ler=nf90_get_att(pncid(ipf),idv,'add_offset',laddoff)
    addoff=laddoff
    if (ler/=nf90_noerr) addoff=0.
    ler=nf90_get_att(pncid(ipf),idv,'scale_factor',lsf)
    sf=lsf
    if (ler/=nf90_noerr) sf=1.
    lier=nf90_get_var(pncid(ipf),idv,rvar,start=start,count=ncount)
    ier=lier
    call ncmsg(name,ier)
    ! unpack data
    rvar=rvar*sf+addoff
  end if ! ier
#endif

  if (qtest) then
    ! expected restart file
    var(1:pil*pjl*pnpan)=rvar
  else
    ! expected mesonest file
    if (ipf==mynproc-1.and.myid<iptst2) then
      jpmax=iptst2
      lcomm=comm_ip
    else
      jpmax=nproc
      lcomm=comm_world
    end if
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

implicit none

integer, intent(in) :: ierr
integer(kind=4) :: lierr
character(len=*), intent(in) :: txt

lierr=ierr
#ifdef usenc3
if (lierr/=nf_noerr) then
  write(6,*) txt," ",nf_strerror(lierr)
#else
if (lierr/=nf90_noerr) then
  write(6,*) txt," ",nf90_strerror(lierr)
#endif
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
      
integer, parameter :: nihead   = 54
      
integer, dimension(nihead) :: ahead
integer, dimension(0:5) :: duma,dumb
integer, dimension(6) :: idum
integer, intent(out) :: ncid,ier
integer resid,is,ipf,dmode
integer ipin,nxpr,nypr
integer ltst,der
integer(kind=4), dimension(nihead) :: lahead
integer(kind=4) :: lncid,lier,ler,lidum
character(len=*), intent(in) :: ifile
character(len=170) :: pfile
character(len=8) :: fdecomp
logical omode

if (myid==0) then
  ! attempt to open single file with myid==0
#ifdef usenc3
  lier=nf_open(ifile,nf_nowrite,lncid)
#else
  lier=nf90_open(ifile,nf90_nowrite,lncid)
#endif
  ncid=lncid
  ier=lier
  fnproc=1      ! number of files to be read over all processors
  dmode=0       ! Single file (dmode=0), Face decomposition (dmode=1), Depreciated (dmode=2) or Uniform decomposition (dmode=3)
  pil=0         ! Number of X grid points within a file panel
  pjl=0         ! Number of Y grid points within a file panel
  pnpan=0       ! Number of panels in file
  ptest=.false. ! Files match current processor (e.g., Restart file), allowing MPI gather/scatter to be avoided
  pfall=.false. ! Every processor has been assigned at least one file, no need to Bcast metadata data
      
  ! attempt to open parallel files
#ifdef usenc3
  if (lier/=nf_noerr) then
    write(pfile,"(a,'.',i6.6)") trim(ifile), 0
    lier=nf_open(pfile,nf_nowrite,lncid)
    ncid=lncid
    ier=lier
    if (lier/=nf_noerr) then
      write(6,*) "WARN: Cannot open ",pfile
      write(6,*) "WARN: Cannot open ",ifile
    else  
      write(6,*) "Found parallel input file ",ifile
      fdecomp=''
      ler=nf_get_att_int(lncid,nf_global,"nproc",lidum)
      fnproc=lidum
      der=ler
      call ncmsg("nproc",der)
      ler=nf_get_att_text(lncid,nf_global,"decomp",fdecomp)
      der=ler
      call ncmsg("decomp",der)
#else
  if (lier/=nf90_noerr) then
    write(pfile,"(a,'.',i6.6)") trim(ifile), 0
    lier=nf90_open(pfile,nf90_nowrite,lncid)
    ncid=lncid
    ier=lier
    if (lier/=nf90_noerr) then
      write(6,*) "WARN: Cannot open ",pfile
      write(6,*) "WARN: Cannot open ",ifile
    else  
      write(6,*) "Found parallel input file ",ifile
      fdecomp=''
      ler=nf90_get_att(lncid,nf90_global,"nproc",lidum)
      fnproc=lidum
      der=ler
      call ncmsg("nproc",der)
      ler=nf90_get_att(lncid,nf90_global,"decomp",fdecomp)
      der=ler
      call ncmsg("decomp",der)
#endif
      select case(fdecomp)
        case('face')
          dmode=1
        case('uniform')  ! old uniform
          dmode=2
        case('uniform1') ! new uniform (Dix style)
          dmode=3
        case default
          write(6,*) "ERROR: Unknown decomposition ",fdecomp
          call ccmpi_abort(-1)
      end select
    end if
  else
    ! nproc should only exist in multi-file input
#ifdef usenc3
    ler=nf_get_att_int(lncid,nf_global,"nproc",lidum)
    if (ler==nf_noerr) then
#else
    ler=nf90_get_att(lncid,nf90_global,"nproc",lidum)
    if (ler==nf90_noerr) then
#endif
      write(6,*) "ERROR: Incorrect base filename"
      call ccmpi_abort(-1)
    end if
    write(6,*) "Found single input file ",ifile
  end if

#ifdef usenc3
  if (lier==nf_noerr) then
    ler=nf_get_att_int(lncid,nf_global,"int_header",lahead)
    ahead=lahead
    der=ler
    call ncmsg("int_header",der)
#else
  if (lier==nf90_noerr) then
    ler=nf90_get_att(lncid,nf90_global,"int_header",lahead)
    ahead=lahead
    der=ler
    call ncmsg("int_header",der)
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
      case(2) ! old uniform decomposition
        pnpan=6
        do ipf=0,fnproc-1
          call uniform_set(pil,pjl,pnoff(ipf),duma,dumb,pnpan,pil_g,ipf,fnproc,nxpr,nypr)
          pioff(ipf,:)=duma
          pjoff(ipf,:)=dumb
        end do
      case(3) ! new uniform decomposition
        pnpan=6
        do ipf=0,fnproc-1
          call dix_set(pil,pjl,pnoff(ipf),duma,dumb,pnpan,pil_g,ipf,fnproc,nxpr,nypr)
          pioff(ipf,:)=duma
          pjoff(ipf,:)=dumb
        end do
    end select

    ptest=.false.
#ifdef uniform_decomp
    if (dmode==3) then
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
  if (ptest) then
    idum(5)=1
  else
    idum(5)=0      
  end if
  idum(6)=ier  
end if

! Broadcast file metadata
call ccmpi_bcast(idum(1:6),0,comm_world)
fnproc=idum(1)
pil   =idum(2)
pjl   =idum(3)
pnpan =idum(4)
ptest =(idum(5)==1)
ier   =idum(6)

lier=ier
#ifdef usenc3
if (lier/=nf_noerr) return
#else
if (lier/=nf90_noerr) return
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
      
if (myid==0) then 
  is=1
  pncid(0)=ncid
else
  is=0
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

if (myid<resid) then
  ltst=1
else
  ltst=0
end if
call ccmpi_commsplit(comm_ip,comm_world,ltst,myid)

pfall=(fnproc>=nproc)
if (mynproc>0) then
  ncid=pncid(0)
end if

return
end subroutine histopen
      
!--------------------------------------------------------------
! This subroutine closes parallel input files
subroutine histclose

use cc_mpi
      
implicit none

integer ipf,ipin,plen
integer(kind=4) :: ierr
      
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
integer, dimension(kl) :: ka,kb
real, dimension(:,:), intent(out) :: t
real, dimension(ifull,kk), intent(in) :: told
real, dimension(kk), intent(in) :: sigin
real, dimension(kl) :: wta,wtb
      
if (kk==kl) then
  if (all(abs(sig-sigin)<0.0001)) then
     t(1:ifull,1:kl)=told(:,:)
     return
   end if
end if
      
klapse=0
do k=1,kl
  if ( sig(k)>=sigin(1) ) then
    ka(k)=2
    kb(k)=1
    wta(k)=0.
    wtb(k)=1.
    klapse=k   ! i.e. T lapse correction for k<=klapse
  else if ( sig(k)<=sigin(kk) ) then   ! at top
    ka(k)=kk
    kb(k)=kk-1
    wta(k)=1.
    wtb(k)=0.
  else
    do kin=2,kk-1
      if ( sig(k)>sigin(kin) ) exit
    enddo     ! kin loop
    ka(k)=kin
    kb(k)=kin-1
    wta(k)=(sigin(kin-1)-sig(k))/(sigin(kin-1)-sigin(kin))
    wtb(k)=(sig(k)-sigin(kin)  )/(sigin(kin-1)-sigin(kin))
  endif  !  (sig(k)>=sigin(1)) ... ...
enddo   !  k loop
#ifdef debug
if ( myid==0 ) then
  write(6,*) 'in vertint kk,kl ',kk,kl
  write(6,"('sigin',10f7.4)") (sigin(k),k=1,kk)
  write(6,"('sig  ',10f7.4)") sig
  write(6,*) 'ka ',ka
  write(6,*) 'kb ',kb
  write(6,"('wta',10f7.4)") wta
  write(6,"('wtb',10f7.4)") wtb
endif   !  (myid==0)
#endif

do k=1,kl
  ! N.B. "a" denotes "above", "b" denotes "below"
  t(1:ifull,k)=wta(k)*told(:,ka(k))+wtb(k)*told(:,kb(k))
enddo    ! k loop

if ( n==1 .and. klapse/=0 ) then  ! for T lapse correction
  do k=1,klapse
    ! assume 6.5 deg/km, with dsig=.1 corresponding to 1 km
    t(1:ifull,k)=t(1:ifull,k)+(sig(k)-sigin(1))*6.5/.1
  enddo    ! k loop
else if ( n==2 ) then  ! for qg do a -ve fix
  t(1:ifull,:)=max(t(1:ifull,:),1.e-6)
else if ( n==5 ) then  ! for qfg, qlg do a -ve fix
  t(1:ifull,:)=max(t(1:ifull,:),0.)
endif
      
return
end subroutine vertint

!--------------------------------------------------------------------
! This subroutine advances input date by the amount of time defined by mtimer_r
subroutine datefix(kdate_r,ktime_r,mtimer_r)

use cc_mpi

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
  if (myid==0) then
    write(6,*) 'For Y2K kdate_r altered to: ',kdate_r
  end if
endif
iyr=kdate_r/10000
imo=(kdate_r-10000*iyr)/100
iday=kdate_r-10000*iyr-100*imo
ihr=ktime_r/100
imins=ktime_r-100*ihr
if (myid==0) then
  write(6,*) 'entering datefix iyr,imo,iday,ihr,imins,mtimer_r: ', &
                               iyr,imo,iday,ihr,imins,mtimer_r
end if

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
  if (myid==0) then
    write(6,*)'*** imins increased to 60 from imins = ',imins
  end if
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
if (myid==0) then
  write(6,*)'leaving datefix kdate_r,ktime_r: ',kdate_r,ktime_r
end if

return
end subroutine datefix

!--------------------------------------------------------------------
! Set up number of minutes from beginning of year
subroutine getzinp(fjd,jyear,jmonth,jday,jhour,jmin,mins,allleap)

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
! days from beginning of year (1st Jan is 0)
integer, dimension(12), parameter :: odoy=(/0,31,59,90,120,151,181,212,243,273,304,334/) 
real, intent(out) :: fjd
logical, intent(in), optional :: allleap
logical lleap

lleap=.false.
if (present(allleap)) lleap=allleap

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
if (leap==1.or.lleap) then
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

integer, intent(in) :: cdfid, itype, ndim
integer, intent(in) :: daily
integer, dimension(ndim), intent(in) :: dim
integer ier
integer(kind=4) :: vtype, idv, lcdfid, lndim, lier, lsize, lnum
integer(kind=4), dimension(ndim) :: ldim
real, intent(in) :: xmin, xmax
real scalef, addoff
real(kind=4) :: lscalef,laddoff
character(len=*), intent(in) :: name
character(len=*), intent(in) :: lname
character(len=*), intent(in) :: units

#ifdef usenc3
#ifdef i8r8
  vtype = nf_double
#else
  vtype = nf_float
#endif
if (itype==1) vtype = nf_short
lcdfid=cdfid
lndim=ndim
ldim=dim
lier = nf_def_var(lcdfid, name, vtype, lndim, ldim, idv)
ier=lier
call ncmsg("def_var",ier)
lsize=len_trim(lname)
lier = nf_put_att_text(lcdfid,idv,'long_name',lsize,lname)
ier=lier
call ncmsg("long_name",ier)
lsize=len_trim(units)
if(lsize>0)then
  lier = nf_put_att_text(lcdfid,idv,'units',lsize,units)
  ier=lier
  call ncmsg("units",ier)
endif
if (vtype == nf_short) then
  lnum=1
  lier = nf_put_att_int2(lcdfid,idv,'valid_min',nf_short,lnum,minv)
  ier=lier
  call ncmsg("valid_min",ier)
  lnum=1
  lier = nf_put_att_int2(lcdfid,idv,'valid_max',nf_short,lnum,maxv)
  ier=lier
  call ncmsg("valid_max",ier)
  lnum=1
  lier = nf_put_att_int2(lcdfid,idv,'missing_value',nf_short,lnum,missval)
  ier=lier
  call ncmsg("missing_value",ier)
  scalef=(xmax-xmin)/(real(maxv)-real(minv))
  addoff=xmin-scalef*minv
  lscalef=scalef
  laddoff=addoff
  lnum=1
  lier = nf_put_att_real(lcdfid,idv,'add_offset',nf_float,lnum,laddoff)
  ier=lier
  call ncmsg("add_offset",ier)
  lnum=1
  lier = nf_put_att_real(lcdfid,idv,'scale_factor',nf_float,lnum,lscalef)
  ier=lier
  call ncmsg("scale_factor",ier)
else
  lnum=1
  lier = nf_put_att_real(lcdfid,idv,'missing_value',nf_float,lnum,nf_fill_float)
  ier=lier
  call ncmsg("missing_value",ier)
endif
lnum=5
lier = nf_put_att_text(lcdfid,idv,'FORTRAN_format',lnum,'G11.4')
ier=lier
call ncmsg("FORTRAN_format",ier)
if(daily>0)then
  lnum=5
  lier = nf_put_att_text(lcdfid,idv,'valid_time',lnum,'daily')
  ier=lier
  call ncmsg("valid_time",ier)
endif
#else
#ifdef i8r8
  vtype = nf90_double
#else
  vtype = nf90_float
#endif
if (itype==1) vtype = nf90_short
lcdfid=cdfid
ldim=dim
lnum=1
lier = nf90_def_var(lcdfid, name, vtype, ldim, idv, deflate_level=lnum)
ier=lier
call ncmsg("def_var",ier)
lsize=len_trim(lname)
lier = nf90_put_att(lcdfid,idv,'long_name',lname)
ier=lier
call ncmsg("long_name",ier)
lsize=len_trim(units)
if(lsize>0)then
  lier = nf90_put_att(lcdfid,idv,'units',units)
  ier=lier
  call ncmsg("units",ier)
endif
if (vtype == nf90_short) then
  lier = nf90_put_att(lcdfid,idv,'valid_min',minv)
  ier=lier
  call ncmsg("valid_min",ier)
  lier = nf90_put_att(lcdfid,idv,'valid_max',maxv)
  ier=lier
  call ncmsg("valid_max",ier)
  lier = nf90_put_att(lcdfid,idv,'missing_value',missval)
  ier=lier
  call ncmsg("missing_value",ier)
  scalef=(xmax-xmin)/(real(maxv)-real(minv))
  addoff=xmin-scalef*minv
  lscalef=scalef
  laddoff=addoff
  lier = nf90_put_att(lcdfid,idv,'add_offset',laddoff)
  ier=lier
  call ncmsg("add_offset",ier)
  lier = nf90_put_att(lcdfid,idv,'scale_factor',lscalef)
  ier=lier
  call ncmsg("scale_factor",ier)
else
  lier = nf90_put_att(lcdfid,idv,'missing_value',nf90_fill_float)
  ier=lier
  call ncmsg("missing_value",ier)
endif
lier = nf90_put_att(lcdfid,idv,'FORTRAN_format','G11.4')
ier=lier
call ncmsg("FORTRAN_format",ier)
if(daily>0)then
  lier = nf90_put_att(lcdfid,idv,'valid_time','daily')
  ier=lier
  call ncmsg("valid_time",ier)
endif
#endif
      
return
end subroutine attrib

!--------------------------------------------------------------------
! 3D NETCDF WRITE ARRAY ROUTINES
subroutine histwrt3(var,sname,idnc,iarch,local,lwrite)

use cc_mpi              ! CC MPI routines

implicit none

include 'newmpar.h'     ! Grid parameters

integer, intent(in) :: idnc, iarch
real, dimension(ifull), intent(in) :: var
real, dimension(ifull,1) :: wvar
character(len=*), intent(in) :: sname
logical, intent(in) :: local, lwrite

wvar(:,1)=var(:)
if (.not.lwrite) then
#ifdef usenc3
  wvar(:,1)=real(nf_fill_float)
#else
  wvar(:,1)=real(nf90_fill_float)
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
      
implicit none
      
include 'newmpar.h'      ! Grid parameters
include 'parm.h'         ! Model configuration
      
integer, intent(in) :: idnc,iarch,istep
integer ier, iq, i
integer(kind=4) :: lidnc, lier, mid, vtype
integer(kind=4), dimension(3) :: start, ncount
integer(kind=2), dimension(ifull,istep) :: ipack
real, dimension(ifull,istep), intent(in) :: var
real addoff, scale_f
real(kind=4) :: laddoff,lscale_f
character(len=*), intent(in) :: sname

start = (/ 1, 1, iarch /)
ncount = (/ il, jl, istep /)

lidnc=idnc
#ifdef usenc3
lier=nf_inq_varid(lidnc,sname,mid)
ier=lier
call ncmsg(sname,ier)
lier=nf_inq_vartype(lidnc,mid,vtype)
if(vtype == nf_short)then
  if (all(var>9.8E36)) then
    ipack=missval
  else
    lier=nf_get_att_real(lidnc,mid,'add_offset',laddoff)
    lier=nf_get_att_real(lidnc,mid,'scale_factor',lscale_f)
    addoff=real(laddoff)
    scale_f=real(lscale_f)
    do i=1,istep
      do iq=1,ifull
        ipack(iq,i)=nint(max(min((var(iq,i)-addoff)/scale_f,real(maxv)),real(minv)),2)
      end do
    end do
  end if
  lier=nf_put_vara_int2(lidnc,mid,start,ncount,ipack)
else
#ifdef i8r8
  lier=nf_put_vara_double(lidnc,mid,start,ncount,var)
#else
  lier=nf_put_vara_real(lidnc,mid,start,ncount,var)
#endif
end if
ier=lier
call ncmsg(sname,ier)
if (mod(ktau,nmaxpr)==0.and.myid==0) then
  if (any(var==real(nf_fill_float))) then
    write(6,'("histwrt3 ",a7,i4,a7)') sname,iarch,"missing"
  else
    write(6,'("histwrt3 ",a7,i4)') sname,iarch
  end if
end if
#else
lier=nf90_inq_varid(lidnc,sname,mid)
ier=lier
call ncmsg(sname,ier)
lier=nf90_inquire_variable(lidnc,mid,xtype=vtype)
if(vtype == nf90_short)then
  if (all(var>9.8E36)) then
    ipack=missval
  else
    lier=nf90_get_att(lidnc,mid,'add_offset',laddoff)
    lier=nf90_get_att(lidnc,mid,'scale_factor',lscale_f)
    addoff=real(laddoff)
    scale_f=real(lscale_f)
    do i=1,istep
      do iq=1,ifull
        ipack(iq,i)=nint(max(min((var(iq,i)-addoff)/scale_f,real(maxv)),real(minv)),2)
      end do
    end do
  end if
  lier=nf90_put_var(lidnc,mid,ipack,start=start,count=ncount)
else
  lier=nf90_put_var(lidnc,mid,var,start=start,count=ncount)
end if
ier=lier
call ncmsg(sname,ier)
if (mod(ktau,nmaxpr)==0.and.myid==0) then
  if (any(var==real(nf90_fill_float))) then
    write(6,'("histwrt3 ",a7,i4,a7)') sname,iarch,"missing"
  else
    write(6,'("histwrt3 ",a7,i4)') sname,iarch
  end if
end if
#endif

return
end subroutine fw3l
      
subroutine fw3a(var,sname,idnc,iarch,istep)

use cc_mpi               ! CC MPI routines

implicit none

include 'newmpar.h'      ! Grid parameters
include 'parm.h'         ! Model configuration

integer, intent(in) :: idnc, iarch, istep
integer ier
integer imn, imx, jmn, jmx, iq, i
integer(kind=4) :: lidnc, mid, vtype, lier
integer(kind=4), dimension(3) :: start, ncount
integer(kind=2), dimension(ifull_g,istep) :: ipack
real, dimension(ifull,istep), intent(in) :: var
real, dimension(ifull_g,istep) :: globvar
real addoff, scale_f
real varn, varx
real(kind=4) :: laddoff, lscale_f
character(len=*), intent(in) :: sname
      
call ccmpi_gather(var, globvar)

start = (/ 1, 1, iarch /)
ncount = (/ il_g, jl_g, istep /)

!     find variable index
lidnc=idnc
#ifdef usenc3
lier=nf_inq_varid(lidnc,sname,mid)
ier=lier
call ncmsg(sname,ier)
lier=nf_inq_vartype(lidnc,mid,vtype)
if (vtype == nf_short) then
#else
lier=nf90_inq_varid(lidnc,sname,mid)
ier=lier
call ncmsg(sname,ier)
lier=nf90_inquire_variable(lidnc,mid,xtype=vtype)
if (vtype == nf90_short) then
#endif
  if (all(globvar>9.8e36)) then
    ipack=missval
  else
#ifdef usenc3
    lier=nf_get_att_real(lidnc,mid,'add_offset',laddoff)
    lier=nf_get_att_real(lidnc,mid,'scale_factor',lscale_f)
#else
    lier=nf90_get_att(lidnc,mid,'add_offset',laddoff)
    lier=nf90_get_att(lidnc,mid,'scale_factor',lscale_f)
#endif
    addoff=laddoff
    scale_f=lscale_f
    do i=1,istep
      do iq=1,ifull_g
        ipack(iq,i)=nint(max(min((globvar(iq,i)-addoff)/scale_f,real(maxv)),real(minv)),2)
      end do
    end do
  endif
#ifdef usenc3
  lier=nf_put_vara_int2(lidnc,mid,start,ncount,ipack)
else
#ifdef i8r8
  lier=nf_put_vara_double(lidnc,mid,start,ncount,globvar)
#else
  lier=nf_put_vara_real(lidnc,mid,start,ncount,globvar)
#endif
#else
  lier=nf90_put_var(lidnc,mid,ipack,start=start,count=ncount)
else
  lier=nf90_put_var(lidnc,mid,globvar,start=start,count=ncount)
#endif
endif
ier=lier
call ncmsg(sname,ier)

if(mod(ktau,nmaxpr)==0)then
#ifdef usenc3
  if (any(globvar==real(nf_fill_float))) then
#else
  if (any(globvar==real(nf90_fill_float))) then
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

implicit none

include 'newmpar.h'     ! Grid parameters

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

implicit none

include 'newmpar.h'      ! Grid parameters
include 'parm.h'         ! Model configuration

integer, intent(in) :: idnc, iarch
integer ier, iq, k
integer(kind=4) :: mid, vtype, lier, lidnc
integer(kind=4), dimension(4) :: start, ncount
integer(kind=2), dimension(ifull,kl) :: ipack
real, dimension(ifull,kl), intent(in) :: var
real addoff, scale_f
real(kind=4) :: laddoff, lscale_f
character(len=*), intent(in) :: sname

start = (/ 1, 1, 1, iarch /)
ncount = (/ il, jl, kl, 1 /)

lidnc=idnc
#ifdef usenc3
lier=nf_inq_varid(lidnc,sname,mid)
ier=lier
call ncmsg(sname,ier)
lier = nf_inq_vartype(lidnc, mid, vtype)
if(vtype == nf_short)then
#else
lier=nf90_inq_varid(lidnc,sname,mid)
ier=lier
call ncmsg(sname,ier)
lier = nf90_inquire_variable(lidnc, mid, xtype=vtype)
if(vtype == nf90_short)then
#endif
  if(all(var>9.8e36))then
    ipack=missval
  else
#ifdef usenc3
    lier=nf_get_att_real(lidnc,mid,'add_offset',laddoff)
    lier=nf_get_att_real(lidnc,mid,'scale_factor',lscale_f)
#else
    lier=nf90_get_att(lidnc,mid,'add_offset',laddoff)
    lier=nf90_get_att(lidnc,mid,'scale_factor',lscale_f)
#endif
    addoff=laddoff
    scale_f=lscale_f
    do k=1,kl
      do iq=1,ifull
        ipack(iq,k)=nint(max(min((var(iq,k)-addoff)/scale_f,real(maxv)),real(minv)),2)
      end do
    end do
  endif
#ifdef usenc3
  lier=nf_put_vara_int2(lidnc,mid,start,ncount,ipack)
else
#ifdef i8r8
  lier=nf_put_vara_double(lidnc,mid,start,ncount,var)
#else
  lier=nf_put_vara_real(lidnc,mid,start,ncount,var)
#endif
#else
  lier=nf90_put_var(lidnc,mid,ipack,start=start,count=ncount)
else
  lier=nf90_put_var(lidnc,mid,var,start=start,count=ncount)
#endif
endif
ier=lier
call ncmsg(sname,ier)

if(mod(ktau,nmaxpr)==0.and.myid==0)then
#ifdef usenc3
  if (any(var==real(nf_fill_float))) then
#else
  if (any(var==real(nf90_fill_float))) then
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

implicit none

include 'newmpar.h'     ! Grid parameters
include 'parm.h'        ! Model configuration

integer, intent(in) :: idnc, iarch
integer ier
integer imx, jmx, kmx, iq, k
integer(kind=4) :: mid, vtype, lier, lidnc
integer(kind=4), dimension(4) :: start, ncount
integer, dimension(2) :: max_result
integer(kind=2), dimension(ifull_g,kl) :: ipack
real addoff, scale_f
real varn, varx
real(kind=4) :: laddoff, lscale_f
real, dimension(ifull,kl), intent(in) :: var
real, dimension(ifull_g,kl) :: globvar
character(len=*), intent(in) :: sname
      
call ccmpi_gather(var, globvar)
start = (/ 1, 1, 1, iarch /)
ncount = (/ il_g, jl_g, kl, 1 /)

!     find variable index
lidnc=idnc
#ifdef usenc3
lier=nf_inq_varid(lidnc,sname,mid)
ier=lier
call ncmsg(sname,ier)
lier = nf_inq_vartype(lidnc, mid, vtype)
if(vtype == nf_short)then
#else
lier=nf90_inq_varid(lidnc,sname,mid)
ier=lier
call ncmsg(sname,ier)
lier = nf90_inquire_variable(lidnc, mid, xtype=vtype)
if(vtype == nf90_short)then
#endif
  if(all(globvar>9.8e36))then
    ipack=missval
  else
#ifdef usenc3
    lier=nf_get_att_real(lidnc,mid,'add_offset',laddoff)
    lier=nf_get_att_real(lidnc,mid,'scale_factor',lscale_f)
#else
    lier=nf90_get_att(lidnc,mid,'add_offset',laddoff)
    lier=nf90_get_att(lidnc,mid,'scale_factor',lscale_f)
#endif
    addoff=laddoff
    scale_f=lscale_f
    do k=1,kl
      do iq=1,ifull_g
        ipack(iq,k)=nint(max(min((globvar(iq,k)-addoff)/scale_f,real(maxv)),real(minv)),2)
      end do
    end do
  endif
#ifdef usenc3
  lier=nf_put_vara_int2(lidnc,mid,start,ncount,ipack)
else
#ifdef i8r8
  lier=nf_put_vara_double(lidnc,mid,start,ncount,globvar)
#else
  lier=nf_put_vara_real(lidnc,mid,start,ncount,globvar)
#endif
#else
  lier=nf90_put_var(lidnc,mid,ipack,start=start,count=ncount)
else
  lier=nf90_put_var(lidnc,mid,globvar,start=start,count=ncount)
#endif
endif
ier=lier
call ncmsg(sname,ier)

if(mod(ktau,nmaxpr)==0)then
#ifdef usenc3
  if (any(globvar==real(nf_fill_float))) then
#else
  if (any(globvar==real(nf90_fill_float))) then
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

implicit none

integer, intent(out) :: ncid
integer, intent(out), optional :: status
integer ncstatus
integer(kind=4) :: lncid,lncstatus
character(len=*), intent(in) :: fname

#ifdef usenc3
lncstatus = nf_open(fname,nf_nowrite,lncid)
#else
lncstatus = nf90_open(fname,nf90_nowrite,lncid)
#endif
ncid=lncid
ncstatus=lncstatus

if (present(status)) then
  status=ncstatus
else
  call ncmsg("open",ncstatus)
endif

return
end subroutine ccnf_open

subroutine ccnf_create(fname,ncid)

use cc_mpi

implicit none

integer, intent(out) :: ncid
integer ncstatus
integer(kind=4) :: lncid, lncstatus
character(len=*), intent(in) :: fname

#ifdef usenc3
  lncstatus = nf_create(fname,nf_clobber,lncid)
#else
  lncstatus = nf90_create(fname,nf90_netcdf4,lncid)
#endif
ncid=lncid
ncstatus=lncstatus

call ncmsg("create",ncstatus)

return
end subroutine ccnf_create

subroutine ccnf_close(ncid)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer ncstatus
integer(kind=4) :: lncid,lncstatus

lncid=ncid
#if usenc3
lncstatus = nf_close(lncid)
#else
lncstatus = nf90_close(lncid)
#endif
ncstatus=lncstatus
call ncmsg("close",ncstatus)

return
end subroutine ccnf_close

subroutine ccnf_nofill(ncid)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer ncstatus
integer(kind=4) :: lncid,lncstatus,lomode

lncid=ncid
#ifdef usenc3
lncstatus=nf_set_fill(lncid,nf_nofill,lomode)
#else
lncstatus=nf90_set_fill(lncid,nf90_nofill,lomode)
#endif
ncstatus=lncstatus
call ncmsg("nofill",ncstatus)

return
end subroutine ccnf_nofill

subroutine ccnf_inq_dimid(ncid,dname,did,tst)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer, intent(out) :: did
integer ncstatus
integer(kind=4) :: lncid,ldid,lncstatus
logical ltst
logical, intent(out), optional :: tst
character(len=*), intent(in) :: dname

lncid=ncid
#ifdef usenc3
lncstatus=nf_inq_dimid(lncid,dname,ldid)
ltst=(lncstatus/=nf_noerr)
#else
lncstatus=nf90_inq_dimid(lncid,dname,ldid)
ltst=(lncstatus/=nf90_noerr)
#endif
did=ldid
ncstatus=lncstatus

if (present(tst)) then
  tst=ltst
else
  call ncmsg("dimid",ncstatus)
end if

return
end subroutine ccnf_inq_dimid

subroutine ccnf_inq_dimlen(ncid,dname,dlen,failok)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer, intent(inout) :: dlen
integer(kind=4) :: lncid,lncstatus,ldid,ldlen
character(len=*), intent(in) :: dname
logical, intent(in), optional :: failok
logical ftest

ftest=.false.
if (present(failok)) ftest=failok

lncid=ncid
ldlen=dlen
#ifdef usenc3
lncstatus=nf_inq_dimid(lncid,dname,ldid)
if (lncstatus/=nf_noerr) then
  if (ftest) return
  write(6,*) nf_strerror(lncstatus)
  call ccmpi_abort(-1)
end if
lncstatus=nf_inq_dimlen(lncid,ldid,ldlen)
if (lncstatus/=nf_noerr) then
  if (ftest) return
  write(6,*) nf_strerror(lncstatus)
  call ccmpi_abort(-1)
end if
#else
lncstatus=nf90_inq_dimid(lncid,dname,ldid)
if (lncstatus/=nf90_noerr) then
  if (ftest) return
  write(6,*) nf90_strerror(lncstatus)
  call ccmpi_abort(-1)
end if
lncstatus=nf90_inquire_dimension(lncid,ldid,len=ldlen)
if (lncstatus/=nf90_noerr) then
  if (ftest) return
  write(6,*) nf90_strerror(lncstatus)
  call ccmpi_abort(-1)
end if
#endif
dlen=ldlen

return
end subroutine ccnf_inq_dimlen

subroutine ccnf_inq_varid(ncid,vname,vid,tst)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer, intent(out) :: vid
integer ncstatus
integer(kind=4) :: lncid,lvid,lncstatus
character(len=*), intent(in) :: vname
logical, intent(out), optional :: tst
logical ltst

lncid=ncid
#ifdef usenc3
lncstatus = nf_inq_varid(lncid,vname,lvid)
ltst=(lncstatus/=nf_noerr)
#else
lncstatus = nf90_inq_varid(lncid,vname,lvid)
ltst=(lncstatus/=nf90_noerr)
#endif
vid=lvid
ncstatus=lncstatus

if (present(tst)) then
  tst=ltst
else
  call ncmsg(vname,ncstatus)
end if

return
end subroutine ccnf_inq_varid

subroutine ccnf_def_dim(ncid,dname,nsize,did)

use cc_mpi

implicit none

integer, intent(in) :: ncid,nsize
integer, intent(out) :: did
integer ncstatus
integer(kind=4) :: lncid,lnsize,ldid,lncstatus
character(len=*), intent(in) :: dname

lncid=ncid
lnsize=nsize
#ifdef usenc3
lncstatus=nf_def_dim(lncid,dname,lnsize,ldid)
#else
lncstatus=nf90_def_dim(lncid,dname,lnsize,ldid)
#endif
did=ldid
ncstatus=lncstatus
call ncmsg("def_dim",ncstatus)

return
end subroutine ccnf_def_dim

subroutine ccnf_def_dimu(ncid,dname,did)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer, intent(out) :: did
integer ncstatus
integer(kind=4) :: lncid,ldid,lncstatus
character(len=*), intent(in) :: dname

lncid=ncid
#ifdef usenc3
lncstatus=nf_def_dim(lncid,dname,nf_unlimited,ldid)
#else
lncstatus=nf90_def_dim(lncid,dname,nf90_unlimited,ldid)
#endif
did=ldid
ncstatus=lncstatus
call ncmsg("def_dimu",ncstatus)

return
end subroutine ccnf_def_dimu

subroutine ccnf_def_var(ncid,vname,vtype,vndim,dims,vid,deflate)

use cc_mpi

implicit none

integer, intent(in) :: ncid,vndim
integer, intent(in), optional :: deflate
integer, intent(out) :: vid
integer, dimension(vndim), intent(in) :: dims
integer ncstatus
integer(kind=4) :: lncid,lvndim,ldef,ltype,lvid,lncstatus
integer(kind=4), dimension(vndim) :: ldims
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

lncid=ncid
ldims=dims
#ifdef usenc3
lvndim=vndim
lncstatus = nf_def_var(lncid,vname,ltype,lvndim,ldims,lvid)
#else
ldef=1
if (present(deflate)) ldef=deflate
lncstatus = nf90_def_var(lncid,vname,ltype,ldims,lvid,deflate_level=ldef)
#endif
vid=lvid
ncstatus=lncstatus
call ncmsg("def_var",ncstatus)

return
end subroutine ccnf_def_var

subroutine ccnf_def_var0(ncid,vname,vtype,vid)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer, intent(out) :: vid
integer ncstatus
integer(kind=4) :: lncid,lvid,lncstatus,ltype,lone,lzero
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

lncid=ncid
#ifdef usenc3
lone=1
lzero=0
lncstatus = nf_def_var(lncid,vname,ltype,lzero,lone,lvid)
#else
lncstatus = nf90_def_var(lncid,vname,ltype,lvid)
#endif
vid=lvid
ncstatus=lncstatus
call ncmsg("def_var0",ncstatus)

return
end subroutine ccnf_def_var0

subroutine ccnf_enddef(ncid)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer ncstatus
integer(kind=4) :: lncid,lncstatus

lncid=ncid
#ifdef usenc3
lncstatus=nf_enddef(lncid)
#else
lncstatus=nf90_enddef(lncid)
#endif
ncstatus=lncstatus
call ncmsg("enddef",ncstatus)

return
end subroutine ccnf_enddef

subroutine ccnf_redef(ncid)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer ncstatus
integer(kind=4) :: lncid,lncstatus

lncid=ncid
#ifdef usenc3
lncstatus=nf_redef(lncid)
#else
lncstatus=nf90_redef(lncid)
#endif
ncstatus=lncstatus
call ncmsg("redef",ncstatus)

return
end subroutine ccnf_redef

subroutine ccnf_sync(ncid)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer ncstatus
integer(kind=4) :: lncid,lncstatus

#ifdef outsync
lncid=ncid
#ifdef usenc3
lncstatus = nf_sync(lncid)
#else
lncstatus = nf90_sync(lncid)
#endif
ncstatus=lncstatus
call ncmsg("sync",ncstatus)
#endif

return
end subroutine ccnf_sync

subroutine ccnf_get_var_real(ncid,vname,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer ncstatus
integer(kind=4) :: lvid,lncstatus,lncid
real, dimension(:), intent(out) :: vdat
character(len=*), intent(in) :: vname

lncid=ncid
#ifdef usenc3
lncstatus = nf_inq_varid(lncid,vname,lvid)
ncstatus=lncstatus
call ncmsg("get_var_varid",ncstatus)
#ifdef i8r8
lncstatus = nf_get_var_double(lncid,lvid,vdat)
#else
lncstatus = nf_get_var_real(lncid,lvid,vdat)
#endif
ncstatus=lncstatus
call ncmsg("get_var",ncstatus)
#else
lncstatus = nf90_inq_varid(lncid,vname,lvid)
ncstatus=lncstatus
call ncmsg("get_var_varid",ncstatus)
lncstatus = nf90_get_var(lncid,lvid,vdat)
ncstatus=lncstatus
call ncmsg("get_var",ncstatus)
#endif

return
end subroutine ccnf_get_var_real

subroutine ccnf_get_var_int(ncid,vname,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer ncstatus
integer, dimension(:), intent(out) :: vdat
integer(kind=4) :: lncid,lncstatus,lvid
integer(kind=4), dimension(size(vdat)) :: lvdat
character(len=*), intent(in) :: vname

lncid=ncid
#ifdef usenc3
lncstatus = nf_inq_varid(lncid,vname,lvid)
ncstatus=lncstatus
call ncmsg("get_var_varid",ncstatus)
lncstatus = nf_get_var_int(lncid,lvid,lvdat)
vdat=lvdat
ncstatus=lncstatus
call ncmsg("get_var",ncstatus)
#else
lncstatus = nf90_inq_varid(lncid,vname,lvid)
ncstatus=lncstatus
call ncmsg("get_var_varid",ncstatus)
lncstatus = nf90_get_var(lncid,lvid,vdat)
ncstatus=lncstatus
call ncmsg("get_var",ncstatus)
#endif

return
end subroutine ccnf_get_var_int

subroutine ccnf_get_var1_real(ncid,vid,start,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid,vid,start
integer ncstatus
integer(kind=4) :: lncid,lvid,lncstatus
#ifdef usenc3
integer(kind=4) :: lstart
#else
integer(kind=4), dimension(1) :: lstart,lcount
real, dimension(1) :: ldat
#endif
real, intent(out) :: vdat

lncid=ncid
lvid=vid
#ifdef usenc3
lstart=start
#ifdef i8r8
lncstatus = nf_get_var1_double(lncid,lvid,lstart,vdat)
#else
lncstatus = nf_get_var1_real(lncid,lvid,lstart,vdat)
#endif
#else
lstart=start
lcount=1
lncstatus = nf90_get_var(lncid,lvid,ldat,start=lstart,count=lcount)
vdat=ldat(1)
#endif
ncstatus=lncstatus
call ncmsg("get_var1",ncstatus)

return
end subroutine ccnf_get_var1_real

subroutine ccnf_get_var1_int(ncid,vid,start,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid,vid,start
integer ncstatus
integer(kind=4) :: lncid,lvid,lncstatus
#ifdef usenc3
integer(kind=4) :: lstart,ldat
#else
integer(kind=4), dimension(1) :: lstart,lcount
integer, dimension(1) :: ldat
#endif
integer, intent(out) :: vdat

lncid=ncid
lvid=vid
#ifdef usenc3
lstart=start
#ifdef i8r8
lncstatus = nf_get_var1_int(lncid,lvid,lstart,ldat)
#else
lncstatus = nf_get_var1_int(lncid,lvid,lstart,ldat)
#endif
vdat=ldat
#else
lstart=start
lcount=1
lncstatus = nf90_get_var(lncid,lvid,ldat,start=lstart,count=lcount)
vdat=ldat(1)
#endif
ncstatus=lncstatus
call ncmsg("get_var1",ncstatus)

return
end subroutine ccnf_get_var1_int

subroutine ccnf_get_vara_real2r(ncid,vid,start,ncount,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid,vid
integer, dimension(:), intent(in) :: start,ncount
integer ncstatus
integer(kind=4) :: lncid,lvid,lncstatus
integer(kind=4), dimension(size(start)) :: lstart
integer(kind=4), dimension(size(ncount)) :: lncount
real, dimension(:), intent(out) :: vdat

lncid=ncid
lvid=vid
lstart=start
lncount=ncount
#ifdef usenc3
#ifdef i8r8
lncstatus=nf_get_vara_double(lncid,lvid,lstart,lncount,vdat)
#else
lncstatus=nf_get_vara_real(lncid,lvid,lstart,lncount,vdat)
#endif
#else
lncstatus=nf90_get_var(lncid,lvid,vdat,start=lstart,count=lncount)
#endif
ncstatus=lncstatus
call ncmsg("get_vara",ncstatus)

return
end subroutine ccnf_get_vara_real2r 

subroutine ccnf_get_vara_real3r(ncid,vid,start,ncount,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid,vid
integer, dimension(:) :: start,ncount
integer ncstatus
integer(kind=4) :: lncid,lvid,lncstatus
integer(kind=4), dimension(size(start)) :: lstart
integer(kind=4), dimension(size(ncount)) :: lncount
real, dimension(:,:), intent(out) :: vdat

lncid=ncid
lvid=vid
lstart=start
lncount=ncount
#ifdef usenc3
#ifdef i8r8
lncstatus=nf_get_vara_double(lncid,lvid,lstart,lncount,vdat)
#else
lncstatus=nf_get_vara_real(lncid,lvid,lstart,lncount,vdat)
#endif
#else
lncstatus=nf90_get_var(lncid,lvid,vdat,start=lstart,count=lncount)
#endif
ncstatus=lncstatus
call ncmsg("get_vara",ncstatus)

return
end subroutine ccnf_get_vara_real3r

subroutine ccnf_get_vara_real4r(ncid,vid,start,ncount,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid,vid
integer, dimension(:) :: start,ncount
integer ncstatus
integer(kind=4) :: lncid,lvid,lncstatus
integer(kind=4), dimension(size(start)) :: lstart
integer(kind=4), dimension(size(ncount)) :: lncount
real, dimension(:,:,:), intent(out) :: vdat

lncid=ncid
lvid=vid
lstart=start
lncount=ncount
#ifdef usenc3
#ifdef i8r8
lncstatus=nf_get_vara_double(lncid,lvid,lstart,lncount,vdat)
#else
lncstatus=nf_get_vara_real(lncid,lvid,lstart,lncount,vdat)
#endif
#else
lncstatus=nf90_get_var(lncid,lvid,vdat,start=lstart,count=lncount)
#endif
ncstatus=lncstatus
call ncmsg("get_vara",ncstatus)

return
end subroutine ccnf_get_vara_real4r

subroutine ccnf_get_vara_int2i(ncid,vid,start,ncount,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid,vid
integer, dimension(:) :: start,ncount
integer ncstatus
integer, dimension(:), intent(out) :: vdat
integer(kind=4) :: lncid,lvid,lncstatus
integer(kind=4), dimension(size(start)) :: lstart
integer(kind=4), dimension(size(ncount)) :: lncount
integer(kind=4), dimension(size(vdat)) :: lvdat

lncid=ncid
lvid=vid
lstart=start
lncount=ncount
#ifdef usenc3
lncstatus=nf_get_vara_int(lncid,lvid,lstart,lncount,lvdat)
vdat=lvdat
#else
lncstatus=nf90_get_var(lncid,lvid,vdat,start=lstart,count=lncount)
#endif
ncstatus=lncstatus
call ncmsg("get_vara",ncstatus)

return
end subroutine ccnf_get_vara_int2i

#ifndef i8r8
subroutine ccnf_get_vara_double4d(ncid,vid,start,ncount,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid,vid
integer, dimension(:) :: start,ncount
integer ncstatus
integer(kind=4) :: lncid,lvid,lncstatus
integer(kind=4), dimension(size(start)) :: lstart
integer(kind=4), dimension(size(ncount)) :: lncount
real(kind=8), dimension(:,:,:), intent(out) :: vdat

lncid=ncid
lvid=vid
lstart=start
lncount=ncount
#ifdef usenc3
lncstatus=nf_get_vara_double(lncid,lvid,lstart,lncount,vdat)
#else
lncstatus=nf90_get_var(lncid,lvid,vdat,start=lstart,count=lncount)
#endif
ncstatus=lncstatus
call ncmsg("get_vara",ncstatus)

return
end subroutine ccnf_get_vara_double4d
#endif

subroutine ccnf_get_att_text(ncid,vid,aname,atext,ierr)

use cc_mpi

implicit none

integer, intent(in) :: ncid,vid
integer, intent(out), optional :: ierr
integer ncstatus
integer(kind=4) :: lncid,lvid,lncstatus
character(len=*), intent(in) :: aname
character(len=*), intent(out) :: atext

atext=''
lncid=ncid
lvid=vid
#ifdef usenc3
lncstatus = nf_get_att_text(lncid,lvid,aname,atext)
#else
lncstatus = nf90_get_att(lncid,lvid,aname,atext)
#endif
ncstatus=lncstatus
if (present(ierr)) then
  ierr=ncstatus
else
  call ncmsg(aname,ncstatus)
end if

return
end subroutine ccnf_get_att_text

subroutine ccnf_get_att_real(ncid,vid,aname,vdat,ierr)

use cc_mpi

implicit none

integer, intent(in) :: ncid,vid
integer, intent(out), optional :: ierr
integer ncstatus
integer(kind=4) :: lncid,lvid,lncstatus
character(len=*), intent(in) :: aname
real, intent(out) :: vdat

lncid=ncid
lvid=vid
#ifdef usenc3
#ifdef i8r8
lncstatus = nf_get_att_double(lncid,lvid,aname,vdat)
#else
lncstatus = nf_get_att_real(lncid,lvid,aname,vdat)
#endif
#else
lncstatus = nf90_get_att(lncid,lvid,aname,vdat)
#endif
ncstatus=lncstatus
if (present(ierr)) then
  ierr=ncstatus
else
  call ncmsg(aname,ncstatus)
end if


return
end subroutine ccnf_get_att_real

subroutine ccnf_get_att_realg1r(ncid,aname,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer ncstatus
integer(kind=4) :: lncid,lncstatus
character(len=*), intent(in) :: aname
real, intent(out) :: vdat

lncid=ncid
#ifdef usenc3
#ifdef i8r8
lncstatus = nf_get_att_double(lncid,nf_global,aname,vdat)
#else
lncstatus = nf_get_att_real(lncid,nf_global,aname,vdat)
#endif
#else
lncstatus = nf90_get_att(lncid,nf90_global,aname,vdat)
#endif
ncstatus=lncstatus
call ncmsg("get_attg",ncstatus)

return
end subroutine ccnf_get_att_realg1r

subroutine ccnf_get_att_realg2r(ncid,aname,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer ncstatus
integer(kind=4) :: lncid,lncstatus
character(len=*), intent(in) :: aname
real, dimension(:), intent(out) :: vdat

lncid=ncid
#ifdef usenc3
#ifdef i8r8
lncstatus = nf_get_att_double(lncid,nf_global,aname,vdat)
#else
lncstatus = nf_get_att_real(lncid,nf_global,aname,vdat)
#endif
#else
lncstatus = nf90_get_att(lncid,nf90_global,aname,vdat)
#endif
ncstatus=lncstatus
call ncmsg("get_attg",ncstatus)

return
end subroutine ccnf_get_att_realg2r

subroutine ccnf_get_att_intg1i(ncid,aname,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer, intent(out) :: vdat
integer ncstatus
integer(kind=4) :: lncid,lncstatus,lvdat
character(len=*), intent(in) :: aname

lncid=ncid
#ifdef usenc3
lncstatus = nf_get_att_int(lncid,nf_global,aname,lvdat)
#else
lncstatus = nf90_get_att(lncid,nf90_global,aname,lvdat)
#endif
vdat=lvdat
ncstatus=lncstatus
call ncmsg("get_attg",ncstatus)

return
end subroutine ccnf_get_att_intg1i

subroutine ccnf_get_att_intg2i(ncid,aname,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer, dimension(:), intent(out) :: vdat
integer ncstatus
integer(kind=4) :: lncid,lncstatus
integer(kind=4), dimension(size(vdat)) :: lvdat
character(len=*), intent(in) :: aname

lncid=ncid
#ifdef usenc3
lncstatus = nf_get_att_int(lncid,nf_global,aname,lvdat)
#else
lncstatus = nf90_get_att(lncid,nf90_global,aname,lvdat)
#endif
vdat=lvdat
ncstatus=lncstatus
call ncmsg("get_attg",ncstatus)

return
end subroutine ccnf_get_att_intg2i

subroutine ccnf_read(fname,vname,vdat)

use cc_mpi

implicit none

include 'newmpar.h'

integer ncstatus
integer(kind=4) :: lncid,lvid,lncstatus
real, dimension(ifull), intent(out) :: vdat
real, dimension(ifull_g) :: vdat_g
character(len=*), intent(in) :: fname
character(len=*), intent(in) :: vname

if (myid==0) then
#ifdef usenc3
  lncstatus = nf_open(fname,nf_nowrite,lncid)
  ncstatus=lncstatus
  call ncmsg(fname,ncstatus)
  lncstatus = nf_inq_varid(lncid,vname,lvid)
  ncstatus=lncstatus
  call ncmsg(fname,ncstatus)
#ifdef i8r8
  lncstatus = nf_get_var_double(lncid,lvid,vdat_g)
#else
  lncstatus = nf_get_var_real(lncid,lvid,vdat_g)
#endif
  ncstatus=lncstatus
  call ncmsg(fname,ncstatus)
  lncstatus = nf_close(lncid)
  ncstatus=lncstatus
  call ncmsg(fname,ncstatus)
#else
  lncstatus = nf90_open(fname,nf90_nowrite,lncid)
  ncstatus=lncstatus
  call ncmsg(fname,ncstatus)
  lncstatus = nf90_inq_varid(lncid,vname,lvid)
  ncstatus=lncstatus
  call ncmsg(fname,ncstatus)
  lncstatus = nf90_get_var(lncid,lvid,vdat_g)
  ncstatus=lncstatus
  call ncmsg(fname,ncstatus)
  lncstatus = nf90_close(lncid)
  ncstatus=lncstatus
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

implicit none

integer, intent(in) :: ncid,vid
integer ncstatus
integer(kind=4) :: lncid,lvid,lncstatus
character(len=*), dimension(:), intent(in) :: vtxt

lncid=ncid
lvid=vid
#ifdef usenc3
lncstatus = nf_put_var_text(lncid,lvid,vtxt)
#else
lncstatus = nf90_put_var(lncid,lvid,vtxt)
#endif
ncstatus=lncstatus
call ncmsg("put_var",ncstatus)

return
end subroutine ccnf_put_var_text2r

subroutine ccnf_put_var_int2i(ncid,vid,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid,vid
integer ncstatus
integer, dimension(:), intent(in) :: vdat
integer(kind=4) :: lncid,lvid,lncstatus
integer(kind=4), dimension(size(vdat)) :: lvdat

lncid=ncid
lvid=vid
#ifdef usenc3
lvdat=vdat
lncstatus = nf_put_var_int(lncid,lvid,lvdat)
#else
lncstatus = nf90_put_var(lncid,lvid,vdat)
#endif
ncstatus=lncstatus
call ncmsg("put_var",ncstatus)

return
end subroutine ccnf_put_var_int2i

subroutine ccnf_put_var_int3i(ncid,vid,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid,vid
integer ncstatus
integer, dimension(:,:), intent(in) :: vdat
integer(kind=4) :: lncid,lvid,lncstatus
integer(kind=4), dimension(size(vdat,1),size(vdat,2)) :: lvdat

lncid=ncid
lvid=vid
#ifdef usenc3
lvdat=vdat
lncstatus = nf_put_var_int(lncid,lvid,lvdat)
#else
lncstatus = nf90_put_var(lncid,lvid,vdat)
#endif
ncstatus=lncstatus
call ncmsg("put_var",ncstatus)

return
end subroutine ccnf_put_var_int3i

subroutine ccnf_put_var1_int(ncid,vid,start,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid,vid,start
integer, intent(in) :: vdat
integer ncstatus
integer(kind=4) :: lncid,lvid,lncstatus
#ifdef usenc3
integer(kind=4) :: lstart,ldat
#else
integer(kind=4), dimension(1) :: lstart,lcount
integer, dimension(1) :: ldat
#endif

lncid=ncid
lvid=vid
#ifdef usenc3
lstart=start
ldat=vdat
lncstatus=nf_put_var1_int(lncid,lvid,lstart,ldat)
#else
lstart=start
lcount=1
ldat(1)=vdat
lncstatus=nf90_put_var(lncid,lvid,ldat,start=lstart,count=lcount)
#endif
ncstatus=lncstatus
call ncmsg("put_var1",ncstatus)

return
end subroutine ccnf_put_var1_int

subroutine ccnf_put_var1_real(ncid,vid,start,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid,vid,start
integer ncstatus
integer(kind=4) :: lncid,lvid,lncstatus
#ifdef usenc3
integer(kind=4) :: lstart
#else
integer(kind=4), dimension(1) :: lstart,lcount
real, dimension(1) :: ldat
#endif
real, intent(in) :: vdat

lncid=ncid
lvid=vid
#ifdef usenc3
lstart=start
#ifdef i8r8
lncstatus=nf_put_var1_double(lncid,lvid,lstart,vdat)
#else
lncstatus=nf_put_var1_real(lncid,lvid,lstart,vdat)
#endif
#else
lstart=start
lcount=1
ldat(1)=vdat
lncstatus=nf90_put_var(lncid,lvid,ldat,start=lstart,count=lcount)
#endif
ncstatus=lncstatus
call ncmsg("put_var1",ncstatus)

return
end subroutine ccnf_put_var1_real

#ifndef i8r8
subroutine ccnf_put_var1_double(ncid,vid,start,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid,vid,start
integer ncstatus
integer(kind=4) :: lncid,lvid,lncstatus
#ifdef usenc3
integer(kind=4) :: lstart
#else
integer(kind=4), dimension(1) :: lstart,lcount
real(kind=8), dimension(1) :: ldat
#endif
real(kind=8), intent(in) :: vdat

lncid=ncid
lvid=vid
#ifdef usenc3
lstart=start
lncstatus=nf_put_var1_double(lncid,lvid,lstart,vdat)
#else
lstart=start
lcount=1
ldat(1)=vdat
lncstatus=nf90_put_var(lncid,lvid,ldat,start=lstart,count=lcount)
#endif
ncstatus=lncstatus
call ncmsg("put_var1",ncstatus)

return
end subroutine ccnf_put_var1_double
#endif

subroutine ccnf_put_vara_real2r(ncid,vid,start,ncount,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid,vid
integer ncstatus
integer, dimension(:), intent(in) :: start,ncount
integer(kind=4) :: lncid,lvid,lncstatus
integer(kind=4), dimension(size(start)) :: lstart
integer(kind=4), dimension(size(ncount)) :: lncount
real, dimension(:), intent(in) :: vdat

lncid=ncid
lvid=vid
lstart=start
lncount=ncount
#ifdef usenc3
#ifdef i8r8
lncstatus=nf_put_vara_double(lncid,lvid,lstart,lncount,vdat)
#else
lncstatus=nf_put_vara_real(lncid,lvid,lstart,lncount,vdat)
#endif
#else
lncstatus=nf90_put_var(lncid,lvid,vdat,start=lstart,count=lncount)
#endif
ncstatus=lncstatus
call ncmsg("put_vara_real2r",ncstatus)

return
end subroutine ccnf_put_vara_real2r

subroutine ccnf_put_vara_real3r(ncid,vid,start,ncount,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid,vid
integer ncstatus
integer, dimension(:), intent(in) :: start,ncount
integer(kind=4) :: lncid,lvid,lncstatus
integer(kind=4), dimension(size(start)) :: lstart
integer(kind=4), dimension(size(ncount)) :: lncount
real, dimension(:,:), intent(in) :: vdat

lncid=ncid
lvid=vid
lstart=start
lncount=ncount
#ifdef usenc3
#ifdef i8r8
lncstatus=nf_put_vara_double(lncid,lvid,lstart,lncount,vdat)
#else
lncstatus=nf_put_vara_real(lncid,lvid,lstart,lncount,vdat)
#endif
#else
lncstatus=nf90_put_var(lncid,lvid,vdat,start=lstart,count=lncount)
#endif
ncstatus=lncstatus
call ncmsg("put_vara",ncstatus)

return
end subroutine ccnf_put_vara_real3r

#ifndef i8r8
subroutine ccnf_put_vara_double2r(ncid,vid,start,ncount,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid,vid
integer ncstatus
integer, dimension(:), intent(in) :: start,ncount
integer(kind=4) :: lncid,lvid,lncstatus
integer(kind=4), dimension(size(start)) :: lstart
integer(kind=4), dimension(size(ncount)) :: lncount
real(kind=8), dimension(:), intent(in) :: vdat

lncid=ncid
lvid=vid
lstart=start
lncount=ncount
#ifdef usenc3
lncstatus=nf_put_vara_double(lncid,lvid,lstart,lncount,vdat)
#else
lncstatus=nf90_put_var(lncid,lvid,vdat,start=lstart,count=lncount)
#endif
ncstatus=lncstatus
call ncmsg("put_vara",ncstatus)

return
end subroutine ccnf_put_vara_double2r
#endif

subroutine ccnf_put_vara_int2i(ncid,vid,start,ncount,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid,vid
integer ncstatus
integer, dimension(:), intent(in) :: start,ncount
integer, dimension(:), intent(in) :: vdat
integer(kind=4) :: lncid,lvid,lncstatus
integer(kind=4), dimension(size(start)) :: lstart
integer(kind=4), dimension(size(ncount)) :: lncount
integer(kind=4), dimension(size(vdat)) :: lvdat

lncid=ncid
lvid=vid
lstart=start
lncount=ncount
#ifdef usenc3
lvdat=vdat
lncstatus=nf_put_vara_int(lncid,lvid,lstart,lncount,lvdat)
#else
lncstatus=nf90_put_var(lncid,lvid,vdat,start=lstart,count=lncount)
#endif
ncstatus=lncstatus
call ncmsg("put_vara",ncstatus)

return
end subroutine ccnf_put_vara_int2i

subroutine ccnf_put_att_text(ncid,vid,aname,asize,atext)

use cc_mpi

implicit none

integer, intent(in) :: ncid,vid,asize
integer ncstatus
integer(kind=4) :: lncid,lvid,lncstatus,lasize
character(len=*), intent(in) :: aname
character(len=asize), intent(in) :: atext

lncid=ncid
lvid=vid
#ifdef usenc3
lasize=asize
lncstatus=nf_put_att_text(lncid,lvid,aname,lasize,atext)
#else
lncstatus=nf90_put_att(lncid,lvid,aname,atext)
#endif
ncstatus=lncstatus
call ncmsg("put_att",ncstatus)

return
end subroutine ccnf_put_att_text

subroutine ccnf_put_att_textg(ncid,aname,atext)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer ncstatus
integer(kind=4) :: lsize,lncid,lncstatus
character(len=*), intent(in) :: aname
character(len=*), intent(in) :: atext

lncid=ncid
#ifdef usenc3
lsize=len_trim(atext)
lncstatus=nf_put_att_text(lncid,nf_global,aname,lsize,atext)
#else
lncstatus=nf90_put_att(lncid,nf90_global,aname,atext)
#endif
ncstatus=lncstatus
call ncmsg("put_attg",ncstatus)

return
end subroutine ccnf_put_att_textg

subroutine ccnf_put_att_intg1(ncid,aname,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer ncstatus
integer, intent(in) :: vdat
integer(kind=4) :: lncid,lncstatus,lone
integer(kind=4), dimension(1) :: ldat
character(len=*), intent(in) :: aname

lncid=ncid
ldat(1)=vdat
#ifdef usenc3
lone=1
lncstatus=nf_put_att_int(lncid,nf_global,aname,nf_int,lone,ldat)
#else
lncstatus=nf90_put_att(lncid,nf90_global,aname,ldat)
#endif
ncstatus=lncstatus
call ncmsg("put_attg",ncstatus)

return
end subroutine ccnf_put_att_intg1

subroutine ccnf_put_att_intg2(ncid,aname,vsize,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid,vsize
integer ncstatus
integer, dimension(:), intent(in) :: vdat
integer(kind=4) :: lncid,lvsize,lncstatus
integer(kind=4), dimension(size(vdat)) :: lvdat
character(len=*), intent(in) :: aname

lncid=ncid
lvdat=vdat
#ifdef usenc3
lvsize=vsize
lncstatus=nf_put_att_int(lncid,nf_global,aname,nf_int,lvsize,lvdat)
#else
lncstatus=nf90_put_att(lncid,nf90_global,aname,lvdat)
#endif
ncstatus=lncstatus
call ncmsg("put_att",ncstatus)

return
end subroutine ccnf_put_att_intg2

subroutine ccnf_put_att_realg1(ncid,aname,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid
integer ncstatus
integer(kind=4) :: lncid,lncstatus
real, intent(in) :: vdat
real, dimension(1) :: ldat
character(len=*), intent(in) :: aname

lncid=ncid
ldat(1)=vdat
#ifdef usenc3
#ifdef i8r8
lncstatus=nf_put_att_double(lncid,nf_global,aname,nf_float,1,ldat)
#else
lncstatus=nf_put_att_real(lncid,nf_global,aname,nf_float,1,ldat)
#endif
#else
lncstatus=nf90_put_att(lncid,nf90_global,aname,ldat)
#endif
ncstatus=lncstatus
call ncmsg("put_attg",ncstatus)

return
end subroutine ccnf_put_att_realg1

subroutine ccnf_put_att_realg2(ncid,aname,vsize,vdat)

use cc_mpi

implicit none

integer, intent(in) :: ncid,vsize
integer ncstatus
integer(kind=4) :: lncid,lvsize,lncstatus
real, dimension(:), intent(in) :: vdat
character(len=*), intent(in) :: aname

lncid=ncid
#ifdef usenc3
lvsize=vsize
#ifdef i8r8
lncstatus=nf_put_att_double(lncid,nf_global,aname,nf_float,lvsize,vdat)
#else
lncstatus=nf_put_att_real(lncid,nf_global,aname,nf_float,lvsize,vdat)
#endif
#else
lncstatus=nf90_put_att(lncid,nf90_global,aname,vdat)
#endif
ncstatus=lncstatus
call ncmsg("put_attg",ncstatus)

return
end subroutine ccnf_put_att_realg2

subroutine surfread(dat,varname,netcdfid,filename)

use cc_mpi

implicit none

include 'newmpar.h'

integer, intent(in), optional :: netcdfid
integer ifully
character(len=*), intent(in), optional :: filename
character(len=*), intent(in) :: varname
real, dimension(:), intent(out) :: dat

ifully=size(dat)
if (myid==0) then
  if (present(filename)) then
    call surfreadglob(dat,varname,filename=filename)  
  else if (present(netcdfid)) then
    call surfreadglob(dat,varname,netcdfid=netcdfid)  
  else
    write(6,*) 'Failed to specify input file'
    call ccmpi_abort(-1)
  end if
else
  if (ifully==ifull) then
    call ccmpi_distribute(dat)
  end if
end if

return
end subroutine surfread

! Read surface data and distribute over processors
! This version suports both netcdf and text file formats
subroutine surfreadglob(dat,vname,netcdfid,filename)

use cc_mpi

implicit none

include 'newmpar.h'   ! Grid parameters
include 'parm.h'      ! Model configuration
include 'parmgeom.h'  ! Coordinate data

integer, intent(in), optional :: netcdfid
integer, dimension(3) :: spos,npos
integer ifully,ncidx,iernc,varid,ierr
integer ilx,jlx
character(len=*), intent(in), optional :: filename
character(len=*), intent(in) :: vname
character(len=47) header
real, dimension(:), intent(out) :: dat
real, dimension(ifull_g) :: glob2d
real rlong0x,rlat0x,schmidtx,dsx

ifully=size(dat)

if (present(filename)) then
  call ccnf_open(filename,ncidx,iernc)
else if (present(netcdfid)) then
  ncidx=netcdfid
  iernc=0
end if

if (iernc==0) then ! Netcdf file

  call ccnf_inq_dimlen(ncidx,'longitude',ilx)
  call ccnf_inq_dimlen(ncidx,'latitude',jlx)
  call ccnf_get_attg(ncidx,'lon0',rlong0x)
  call ccnf_get_attg(ncidx,'lat0',rlat0x)
  call ccnf_get_attg(ncidx,'schmidt',schmidtx)
  if(ilx/=il_g.or.jlx/=jl_g.or.rlong0x/=rlong0.or.rlat0x/=rlat0.or.schmidtx/=schmidt) then
    write(6,*) 'wrong data file supplied'
    call ccmpi_abort(-1)
  end if

  spos(1:3)=1
  npos(1)=il_g
  npos(2)=6*il_g
  npos(3)=1
  call ccnf_inq_varid(ncidx,vname,varid)
  call ccnf_get_vara(ncidx,varid,spos,npos,glob2d)
  if (present(filename)) then
    call ccnf_close(ncidx)
  end if
  
else ! ASCII file

  open(87,file=filename,status='old')
  read(87,*,iostat=ierr) ilx,jlx,rlong0x,rlat0x,schmidtx,dsx,header
  if ( ierr == 0 ) then
    write(6,*) ilx,jlx,rlong0x,rlat0x,schmidtx,dsx,header
    if(ilx/=il_g.or.jlx/=jl_g.or.rlong0x/=rlong0.or.rlat0x/=rlat0.or.schmidtx/=schmidt) then
      write(6,*) 'wrong data file supplied'
      call ccmpi_abort(-1)
    end if
    read(87,*) glob2d
    close(87)
  else if ( ierr < 0 ) then ! Error, so really unformatted file
    close(87)
    write(6,*) 'now doing unformatted read'
    open(87,file=filename,status='old',form='unformatted')
    read(87) glob2d
    close(87)
  else
    write(6,*) "error in surfreadglob",trim(filename),ierr
    call ccmpi_abort(-1)
  end if
  
end if

! distrubte data over processors
if (ifully==ifull) then
  call ccmpi_distribute(dat, glob2d)
else if (ifully==ifull_g) then
  dat=glob2d
else
  write(6,*) "ERROR: Invalid array size in surfreadglob"
  call ccmpi_abort(-1)
end if
  
return
end subroutine surfreadglob


end module infile
