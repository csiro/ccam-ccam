module onthefly_m
    
! Main netcdf input routines.  Host grid is automatically
! interpolated to nested model grid.  Three options are
!   nested=0  Initial conditions
!   nested=1  Nudging fields
!   nested=2  Surface data recycling
      
! This version supports the parallel file routines contained
! in infile.f90.  Hence, restart files do not require any
! gathers and scatters.

! In the case where the grid needs to be interpolated, a copy
! of the input data is sent to all processors and each
! processor performs its own interpolation.

implicit none
    
private
public onthefly, retopo
    
integer, parameter :: nord = 3        ! 1 for bilinear, 3 for bicubic interpolation
integer, save :: ik, jk, kk, ok, nsibx
integer dk
integer, dimension(:,:), allocatable, save :: nface4
real, save :: rlong0x, rlat0x, schmidtx
real, dimension(3,3), save :: rotpoles, rotpole
real, dimension(:,:), allocatable, save :: xg4, yg4
real, dimension(:), allocatable, save :: axs_a, ays_a, azs_a
real, dimension(:), allocatable, save :: bxs_a, bys_a, bzs_a
real, dimension(:), allocatable, save :: sigin
logical iotest, newfile
    
contains

! *****************************************************************************
! Main interface for input data that reads grid metadata
    
subroutine onthefly(nested,kdate_r,ktime_r,psl,zss,tss,sicedep,fracice,t,u,v,qg,tgg,wb,wbice,snowd,qfg, &
                    qlg,qrg,tggsn,smass,ssdn,ssdnn,snage,isflag,mlodwn,ocndwn,xtgdwn)

use aerosolldr       ! LDR prognostic aerosols
use cc_mpi           ! CC MPI routines
use infile           ! Input file routines
use mlo              ! Ocean physics and prognostic arrays
use soil_m           ! Soil and surface data

implicit none

include 'newmpar.h'  ! Grid parameters
include 'darcdf.h'   ! Netcdf data
include 'parm.h'     ! Model configuration
include 'stime.h'    ! File date data

integer, parameter :: nihead = 54
integer, parameter :: nrhead = 14

integer, intent(in) :: nested
integer, intent(out) :: kdate_r, ktime_r
integer, save :: maxarchi
integer ier, mtimer, k, ierx, idvkd, idvkt, idvmt
integer, dimension(nihead) :: nahead
integer, dimension(ifull), intent(out) :: isflag
real timer
real, dimension(ifull,wlev,4), intent(out) :: mlodwn
real, dimension(ifull,kl,naero), intent(out) :: xtgdwn
real, dimension(ifull,ms), intent(out) :: wb, wbice, tgg
real, dimension(ifull,3), intent(out) :: tggsn, smass, ssdn
real, dimension(ifull,2), intent(out) :: ocndwn
real, dimension(:,:), intent(out) :: t, u, v, qg, qfg, qlg, qrg
real, dimension(ifull), intent(out) :: psl, zss, tss, fracice, snowd
real, dimension(ifull), intent(out) :: sicedep, ssdnn, snage
real, dimension(nrhead) :: ahead
real, dimension(14) :: rdum
logical ltest, tst

call START_LOG(onthefly_begin)
!--------------------------------------------------------------------
! pfall indicates all processors have a parallel input file and there
! is no need to broadcast metadata (see infile.f90).  Otherwise read
! metadata on myid=0 and broadcast that data to all processors.
if ( myid==0 .or. pfall ) then
  if ( myid==0 ) write(6,*) 'Entering onthefly for nested,ktau = ',nested,ktau
  
  ! Locate new file and read grid metadata --------------------------
  if ( ncid/=ncidold ) then
    if ( myid==0 ) write(6,*) 'Reading new file metadata'
    iarchi=1   ! initial time index for input file
    maxarchi=0 ! initial number of timesteps in input file
    ok=0       ! initial number of ocean levels
    call ccnf_get_attg(ncid,'int_header',nahead)
    call ccnf_get_attg(ncid,'real_header',ahead)
    ik      =nahead(1)  ! grid size
    jk      =nahead(2)  ! grid size
    kk      =nahead(3)  ! vertical levels
    nsibx   =nahead(44) ! land-surface parameterisation
    rlong0x =ahead(5)   ! longitude
    rlat0x  =ahead(6)   ! latitude
    schmidtx=ahead(7)   ! schmidt factor
    if ( schmidtx<=0. .or. schmidtx>1. ) then
      ! backwards compatibility option
      rlong0x =ahead(6)
      rlat0x  =ahead(7)
      schmidtx=ahead(8)
    endif  ! (schmidtx<=0..or.schmidtx>1.)        
    call ccnf_inq_dimlen(ncid,'time',maxarchi)
    call ccnf_inq_dimlen(ncid,'olev',ok,failok=.true.)
    if ( myid==0 ) then
      write(6,*) "Found ik,jk,kk,ok ",ik,jk,kk,ok
      write(6,*) "      maxarchi ",maxarchi
      write(6,*) "      rlong0x,rlat0x,schmidtx ",rlong0x,rlat0x,schmidtx
    end if
  end if
  
  ! search for required date ----------------------------------------
  if ( myid==0 ) write(6,*)'Search for kdate_s,ktime_s >= ',kdate_s,ktime_s
  ltest=.true.     ! flag indicates that the date is not yet found
  iarchi=iarchi-1  ! move time index back one step to check current position in file
  ierx=0           ! indicates normal mtimer format or backwards compatibility mode
  call ccnf_inq_varid(ncid,'kdate',idvkd,tst)
  call ccnf_inq_varid(ncid,'ktime',idvkt,tst)
  call ccnf_inq_varid(ncid,'mtimer',idvmt,tst)
  if ( tst ) then
    ! backwards compatability option
    ierx=1
    call ccnf_inq_varid(ncid,'timer',idvmt,tst)
  end if
  ! start search for required date/time
  do while( ltest .and. iarchi<maxarchi )
    ! could read this as one array, but we only usually need to advance 1 step
    iarchi=iarchi+1
    call ccnf_get_var1(ncid,idvkd,iarchi,kdate_r)
    call ccnf_get_var1(ncid,idvkt,iarchi,ktime_r)
    if ( ierx==0 ) then
      call ccnf_get_var1(ncid,idvmt,iarchi,mtimer)
      timer=mtimer/60.
    else
      timer=0.
      call ccnf_get_var1(ncid,idvmt,iarchi,timer)
      mtimer=nint(timer*60.)
    endif
    if ( mtimer>0 ) then
      ! calculate date if mtimer>0
      call datefix(kdate_r,ktime_r,mtimer)
    end if
    ! ltest = .false. when correct date is found
    ltest=2400*(kdate_r-kdate_s)-1200*nsemble+(ktime_r-ktime_s)<0
  end do
  if ( nsemble/=0 ) then
    kdate_r=kdate_s
    ktime_r=ktime_s
  end if
  if ( ltest ) then
    ! ran out of file before correct date was located
    ktime_r=-1
  end if
  if ( myid==0 ) then
    write(6,*) 'After search ltest,iarchi =',ltest,iarchi
    write(6,*) '             kdate_r,ktime_r =',kdate_r,ktime_r
  end if

endif  ! ( myid==0 .or. pfall )

! if metadata is not read by all processors, then broadcast ---------
if ( .not.pfall ) then
  rdum(1)=rlong0x
  rdum(2)=rlat0x
  rdum(3)=schmidtx
  ! kdate_r is too large to represent as a single real, so
  ! we split kdate_r into year, month and day
  rdum(4)=real(kdate_r/10000)
  rdum(5)=real(kdate_r/100-nint(rdum(4))*100)
  rdum(6)=real(kdate_r-nint(rdum(4))*10000-nint(rdum(5))*100)
  rdum(7)=real(ktime_r)
  if ( ncid/=ncidold ) then
    rdum(8)=1.
  else
    rdum(8)=0.
  end if
  rdum(9) =real(ik)
  rdum(10)=real(jk)
  rdum(11)=real(kk)
  rdum(12)=real(ok)
  rdum(13)=real(iarchi)
  rdum(14)=real(nsibx)
  call ccmpi_bcast(rdum(1:14),0,comm_world)
  rlong0x =rdum(1)
  rlat0x  =rdum(2)
  schmidtx=rdum(3)
  kdate_r =nint(rdum(4))*10000+nint(rdum(5))*100+nint(rdum(6))
  ktime_r =nint(rdum(7))
  newfile =(nint(rdum(8))==1)
  ik      =nint(rdum(9))
  jk      =nint(rdum(10))
  kk      =nint(rdum(11))
  ok      =nint(rdum(12))
  iarchi  =nint(rdum(13))
  nsibx   =nint(rdum(14))
else
  newfile =(ncid/=ncidold)
end if

! mark current file as read for metadata
if ( newfile ) ncidold = ncid

! trap error if correct date/time is not located --------------------
if ( ktime_r<0 ) then
  if ( nested==2 ) then
    if ( myid==0 ) then
      write(6,*) "WARN: Cannot locate date/time in input file"
    end if
    return
  else
    write(6,*) "ERROR: Cannot locate date/time in input file"
    call ccmpi_abort(-1)
  end if
end if
!--------------------------------------------------------------------
      
! Here we call ontheflyx with different automatic array
! sizes.  This means the arrays are correct for interpolation
! and file i/o on myid==0, as well as the arrays are smaller
! on myid/=0 when they are not needed.  This way we avoid
! having to maintain multiple ontheflyx subroutines.
      
! Note that if histrd fails to find a variable, it returns
! zero in the output array
      
if ( myid==0 ) then
  dk=ik ! non-zero automatic array size in onthefly_work
else
  dk=0  ! zero automatic array size in onthefly_work
end if
call onthefly_work(nested,kdate_r,ktime_r,psl,zss,tss,sicedep,fracice,t,u,v,qg,tgg,wb,wbice, &
                   snowd,qfg,qlg,qrg,tggsn,smass,ssdn,ssdnn,snage,isflag,mlodwn,ocndwn,      &
                   xtgdwn)
if ( myid==0 ) write(6,*) "Leaving onthefly"

call END_LOG(onthefly_end)

return
end subroutine onthefly


! *****************************************************************************
! Read data from netcdf file
      
! arrays are typically read as global and then distributed to
! processor local arrays.  This allows for more flexibility
! with diagnosed fields.  However if there is one file per process
! (e.g., for restart files), then there is no need for message
! passing.  Data is usually read in as 2D fields which avoids
! memory problems when the host grid size is significantly
! larger than the regional grid size.
subroutine onthefly_work(nested,kdate_r,ktime_r,psl,zss,tss,sicedep,fracice,t,u,v,qg,tgg,wb,wbice, &
                         snowd,qfg,qlg,qrg,tggsn,smass,ssdn,ssdnn,snage,isflag,mlodwn,ocndwn,      &
                         xtgdwn)
      
use aerosolldr, only : ssn,naero               ! LDR aerosol scheme
use ateb, only : atebdwn                       ! Urban
use cable_def_types_mod, only : ncs, ncp       ! CABLE dimensions
use casadimension, only : mplant,mlitter,msoil ! CASA dimensions
use carbpools_m                                ! Carbon pools
use cc_mpi                                     ! CC MPI routines
use cfrac_m                                    ! Cloud fraction
use cloudmod                                   ! Prognostic strat cloud
use extraout_m                                 ! Additional diagnostics      
use infile                                     ! Input file routines
use latlong_m                                  ! Lat/lon coordinates
use mlo, only : wlev,micdwn,mloregrid          ! Ocean physics and prognostic arrays
use mlodynamics                                ! Ocean dynamics
use morepbl_m                                  ! Additional boundary layer diagnostics
use nharrs_m, only : phi_nh,lrestart           ! Non-hydrostatic atmosphere arrays
use nsibd_m, only : isoilm                     ! Land-surface arrays
use river                                      ! River routing
use savuvt_m                                   ! Saved dynamic arrays
use savuv1_m                                   ! Saved dynamic arrays
use screen_m                                   ! Screen level diagnostics
use sigs_m                                     ! Atmosphere sigma levels
use soil_m                                     ! Soil and surface data
use tkeeps, only : tke,eps,zidry               ! TKE-EPS boundary layer
use tracers_m                                  ! Tracer data
use utilities                                  ! Grid utilities
use vecsuv_m                                   ! Map to cartesian coordinates
use vvel_m, only : dpsldt,sdot                 ! Additional vertical velocity
use xarrs_m, only : pslx                       ! Saved dynamic arrays
use workglob_m                                 ! Additional grid interpolation
use work2_m                                    ! Diagnostic arrays

implicit none

include 'newmpar.h'                            ! Grid parameters
include 'const_phys.h'                         ! Physical constants
include 'darcdf.h'                             ! Netcdf data
include 'kuocom.h'                             ! Convection parameters
include 'parm.h'                               ! Model configuration
include 'parmdyn.h'                            ! Dynamics parmaters
include 'parmgeom.h'                           ! Coordinate data
include 'soilv.h'                              ! Soil parameters
include 'stime.h'                              ! File date data

real, parameter :: iotol=1.E-5      ! tolarance for iotest grid matching
      
integer, intent(in) :: kdate_r, ktime_r, nested
integer idv, isoil, nud_test
integer lev, levkk, ier, ierr, igas
integer nemi, id2, jd2, idjd2
integer i, j, k, mm, iq, ii, jj, np, numneg
integer, dimension(dk*dk*6) :: isoilm_a
integer, dimension(ifull), intent(out) :: isflag
integer, dimension(7+3*ms) :: ierc
integer, dimension(3), save :: iers
integer, dimension(1) :: str, cnt
real(kind=8), dimension(:,:), allocatable, save :: xx4, yy4 ! large common arrays
real(kind=8), dimension(dk*dk*6):: z_a, x_a, y_a
real, dimension(ifull,wlev,4), intent(out) :: mlodwn
real, dimension(ifull,kl,naero), intent(out) :: xtgdwn
real, dimension(ifull,ok,2) :: mloin
real, dimension(ifull,kk) :: u_k, v_k
real, dimension(ifull,2), intent(out) :: ocndwn
real, dimension(ifull,ms), intent(out) :: wb, wbice, tgg
real, dimension(ifull,3), intent(out) :: tggsn, smass, ssdn
real, dimension(:,:), intent(out) :: t, u, v, qg, qfg, qlg, qrg
real, dimension(ifull), intent(out) :: psl, zss, tss, fracice
real, dimension(ifull), intent(out) :: snowd, sicedep, ssdnn, snage
real, dimension(ifull) :: dum6, tss_l, tss_s, pmsl
real, dimension(dk*dk*6) :: ucc, vcc
real, dimension(dk*dk*6) :: fracice_a, sicedep_a
real, dimension(dk*dk*6) :: tss_l_a, tss_s_a
real, dimension(dk*dk*6) :: t_a_lev, psl_a, tss_a
real, dimension(dk*dk*6) :: wts_a  ! not used here or defined in call setxyz
real, dimension(:), allocatable, save :: zss_a, ocndep_l
real, dimension(kk+3) :: dumr
real rlongd, rlatd
character(len=8) vname
character(len=3) trnum
logical tsstest, tst
logical, dimension(:), allocatable, save :: land_a, sea_a

! land-sea mask method (nemi=3 use soilt, nemi=2 use tgg, nemi=1 use zs)
nemi = 3
      
! test if retopo fields are required
if ( nud_p==0 .and. nud_t==0 .and. nud_q==0 ) then
  nud_test = 0
else
  nud_test = 1
end if
      
! Determine if interpolation is required
iotest=6*ik*ik==ifull_g .and. abs(rlong0x-rlong0)<iotol .and. abs(rlat0x-rlat0)<iotol .and. &
       abs(schmidtx-schmidt)<iotol .and. nsib==nsibx
if ( iotest ) then
  io_in = 1   ! no interpolation
else
  io_in = -1  ! interpolation
end if
if ( myid==0 ) write(6,*) "Interpolation iotest,io_in =",iotest,io_in

!--------------------------------------------------------------------
! Allocate interpolation, vertical level and mask arrays
! dk is only non-zero on myid==0
if ( .not.allocated(nface4) ) then
  allocate(nface4(ifull,4),xg4(ifull,4),yg4(ifull,4))
end if
if ( newfile ) then
  if ( allocated(sigin) ) then
    deallocate(sigin,land_a,sea_a)
    deallocate(axs_a,ays_a,azs_a)
    deallocate(bxs_a,bys_a,bzs_a)          
  end if
  allocate(sigin(kk),land_a(dk*dk*6),sea_a(dk*dk*6))
  allocate(axs_a(dk*dk*6),ays_a(dk*dk*6),azs_a(dk*dk*6))
  allocate(bxs_a(dk*dk*6),bys_a(dk*dk*6),bzs_a(dk*dk*6))
end if
      
!--------------------------------------------------------------------
! Determine input grid coordinates and interpolation arrays
if ( newfile .and. .not.iotest ) then
  ! xx4 and yy4 could be replaced with sharded arrays in MPI-3
  allocate(xx4(1+4*ik,1+4*ik),yy4(1+4*ik,1+4*ik))

  if ( m_fly==1 ) then
    rlong4_l(:,1) = rlongg(:)*180./pi
    rlat4_l(:,1)  = rlatt(:)*180./pi
  end if
          
  if ( myid==0 ) then
    write(6,*) "Defining input file grid"
!   following setxyz call is for source data geom    ****   
    do iq = 1,dk*dk*6
      axs_a(iq) = iq
      ays_a(iq) = iq
      azs_a(iq) = iq
    enddo      
    call setxyz(ik,rlong0x,rlat0x,-schmidtx,x_a,y_a,z_a,wts_a,axs_a,ays_a,azs_a,bxs_a,bys_a,bzs_a,xx4,yy4)
  end if ! (myid==0)

  call ccmpi_bcastr8(xx4,0,comm_world)
  call ccmpi_bcastr8(yy4,0,comm_world)
  
  ! calculate the rotated coords for host and model grid
  rotpoles = calc_rotpole(rlong0x,rlat0x)
  rotpole  = calc_rotpole(rlong0,rlat0)
  if ( ktau<3 .and. myid==0 ) then
    write(6,*)'m_fly,nord ',m_fly,nord
    write(6,*)'kdate_r,ktime_r,ktau,ds',kdate_r,ktime_r,ktau,ds
    write(6,*)'rotpoles:'
    do i=1,3
      write(6,'(3x,2i1,5x,2i1,5x,2i1,5x,3f8.4)') (i,j,j=1,3),(rotpoles(i,j),j=1,3)
    enddo
  endif                  ! (ktau<3.and.myid==0)
  if ( nmaxpr==1 .and. myid==0 ) then
    write(6,*)'in onthefly rotpole:'
    do i=1,3
      write(6,'(3x,2i1,5x,2i1,5x,2i1,5x,3f8.4)') (i,j,j=1,3),(rotpole(i,j),j=1,3)
    enddo
    write(6,*)'xx4,yy4 ',xx4(id,jd),yy4(id,jd)
    write(6,*)'before latltoij for id,jd: ',id,jd
    write(6,*)'rlong0x,rlat0x,schmidtx ',rlong0x,rlat0x,schmidtx
  endif                  ! (nmaxpr==1.and.myid==0)

  ! setup interpolation arrays
  do mm=1,m_fly  !  was 4, now may be set to 1 in namelist
    do iq=1,ifull
      call latltoij(rlong4_l(iq,mm),rlat4_l(iq,mm),       & !input
                    rlong0x,rlat0x,schmidtx,              & !input
                    xg4(iq,mm),yg4(iq,mm),nface4(iq,mm),  & !output (source)
                    xx4,yy4,ik)
    enddo
  enddo
  deallocate(xx4,yy4)
       
end if ! newfile .and. .not.iotest
      
! -------------------------------------------------------------------
! read time invariant data when file is first opened
! need global zss_a for (potentially) landsea mask and psl interpolation
! need global isoilm_a for (potentially) landsea mask
if ( newfile ) then

  ! read vertical levels and missing data checks
  if ( myid==0 .or. pfall ) then
    if ( myid==0 ) write(6,*) "Reading time invariant fields"
    call ccnf_inq_varid(ncid,'lev',idv,tst)
    if ( tst ) call ccnf_inq_varid(ncid,'layer',idv,tst)
    if ( tst ) call ccnf_inq_varid(ncid,'sigma',idv,tst)
    str(1)=1
    cnt(1)=kk
    call ccnf_get_vara(ncid,idv,str(1:1),cnt(1:1),sigin)
    if ( myid==0 ) write(6,'("sigin=",(9f7.4))') (sigin(k),k=1,kk)
    ! check for missing data
    iers(1:3)=0
    call ccnf_inq_varid(ncid,'mixr',idv,tst)
    if ( tst ) iers(1)=-1
    call ccnf_inq_varid(ncid,'siced',idv,tst)
    if ( tst ) iers(2)=-1
    call ccnf_inq_varid(ncid,'fracice',idv,tst)
    if ( tst ) iers(3)=-1
  end if
  
  ! bcast data to all processors unless all processes are reading input files
  if ( .not.pfall ) then
    dumr(1:kk)     =sigin
    dumr(kk+1:kk+3)=real(iers(1:3))
    call ccmpi_bcast(dumr(1:kk+3),0,comm_world)
    sigin    =dumr(1:kk)
    iers(1:3)=nint(dumr(kk+1:kk+3))
  end if

  ! determine whether surface temperature needs to be interpolated (tsstest=.false.)
  tsstest=(iers(2)==0.and.iers(3)==0.and.iotest)
  if ( myid==0 ) write(6,*) "tsstest,iers ",tsstest,iers(1:3)
  if ( allocated(zss_a) ) deallocate(zss_a)
  if ( tsstest ) then
    ! load local surface temperature
    allocate(zss_a(ifull))
    call histrd1(iarchi,ier,'zht',ik,zss_a,ifull)
  else
    ! load global surface temperature
    allocate(zss_a(6*dk*dk))
    call histrd1(iarchi,ier,'zht',  ik,zss_a,6*ik*ik)
    call histrd1(iarchi,ier,'soilt',ik,ucc,  6*ik*ik)
    if ( myid==0 ) then
      isoilm_a=nint(ucc)
      if ( all(isoilm_a==0) ) isoilm_a=-1 ! missing value flag
    end if
  end if
  
  ! read host ocean bathymetry data
  if ( nmlo/=0 .and. abs(nmlo)<=9 ) then
    if ( .not.allocated(ocndep_l) ) allocate(ocndep_l(ifull))
    call gethist1('ocndepth',ocndep_l)
  end if
  if ( myid==0 ) write(6,*) "Finished reading fixed fields"
  
else
  ! use saved metadata  
  tsstest=(iers(2)==0.and.iers(3)==0.and.iotest)        
endif ! newfile ..else..

! -------------------------------------------------------------------
! detemine the reference level below sig=0.9 (used to calculate psl)
lev=0
do while( sig(lev+1)>0.9 ) ! nested grid
  lev=lev+1
end do
levkk=0
do while( sigin(levkk+1)>0.9 ) ! host grid
  levkk=levkk+1
end do
if ( levkk == 0 ) then
  write(6,*) "ERROR: Invalid sigma levels"
  write(6,*) "sigin = ",sigin
  call ccmpi_abort(-1)
end if

!--------------------------------------------------------------------
! Read surface pressure
! psf read when nested=0 or nested=1.and.nud_p/=0
psl_a=0.
psl=0.
if ( nested==0 .or. ( nested==1 .and. nud_test/=0 ) ) then
  if ( iotest ) then
    call histrd1(iarchi,ier,'psf',ik,psl,ifull)
  else
    call histrd1(iarchi,ier,'psf',ik,psl_a,6*ik*ik)
  end if
endif

! -------------------------------------------------------------------
! Read surface temperature 
! read global tss to diagnose sea-ice or land-sea mask
if ( tsstest ) then
  call histrd1(iarchi,ier,'tsu',ik,tss,ifull)
  zss=zss_a ! used saved zss arrays
else
  call histrd1(iarchi,ier,'tsu',ik,tss_a,6*ik*ik)
      
  ! set up land-sea mask from either soilt, tss or zss
  if ( newfile .and. myid==0 ) then
    if ( nemi==3 ) then 
      land_a(:)=isoilm_a(:)>0
      numneg=count(.not.land_a)
      if ( any(isoilm_a(:)<0) ) nemi=2
    end if ! (nemi==3)
    if ( nemi==2 ) then
      numneg=0
      do iq=1,dk*dk*6
        if ( tss_a(iq)>0. ) then ! over land
          land_a(iq)=.true.
        else                     ! over sea
          land_a(iq)=.false.
          numneg=numneg+1
        endif               ! (tss(iq)>0) .. else ..
      enddo
      if ( numneg==0 ) nemi=1  ! should be using zss in that case
    endif !  (nemi==2)
    tss_a=abs(tss_a)
    if ( nemi==1 ) then
      land_a(:)=zss_a(:)>0.
      numneg=count(.not.land_a)
    endif ! (nemi==1)
    write(6,*)'Land-sea mask using nemi = ',nemi
    sea_a=.not.land_a
  end if ! (newfile.and.myid==0)
end if ! (tsstest) ..else..

      
!--------------------------------------------------------------
! Read ocean data for nudging (sea-ice is read below)
! read when nested=0 or nested==1.and.nud/=0 or nested=2
if ( nmlo/=0 .and. abs(nmlo)<=9 ) then
  ! fixed ocean depth
  ocndwn(:,1)=ocndep_l
  ! ocean potential temperature
  ! ocean temperature and soil temperature use the same arrays
  ! as no fractional land or sea cover is allowed in CCAM
  if ( ( nested/=1 .or. nud_sst/=0 ) .and. ok>0 ) then
    do k=1,ok
      if (k<=ms) then
        write(vname,'("tgg",I1.1)') k
      else
        write(vname,'("tgg",I2.2)') k
      end if
      call fillhist1(vname,mloin(:,k,1),land_a)
    end do
    call mloregrid(ok,ocndwn(:,1),mloin(:,:,1),mlodwn(:,:,1),0)
    if ( all(mlodwn(:,:,1)==0.) ) mlodwn(:,:,1)=293.
  else
    mlodwn(:,:,1)=293.
  end if ! (nestesd/=1.or.nud_sst/=0) ..else..
  ! ocean salinity
  if ( ( nested/=1 .or. nud_sss/=0 ) .and. ok>0 ) then
    do k=1,ok
      write(vname,'("sal",I2.2)') k
      call fillhist1(vname,mloin(:,k,1),land_a)
    end do
    call mloregrid(ok,ocndwn(:,1),mloin(:,:,1),mlodwn(:,:,2),1)
    mlodwn(:,:,2)=max(mlodwn(:,:,2),0.)
  else
    mlodwn(:,:,2)=34.72   
  end if ! (nestesd/=1.or.nud_sss/=0) ..else..
  ! ocean currents
  if ( ( nested/=1 .or. nud_ouv/=0 ) .and. ok>0 ) then
    if ( iotest ) then
      do k=1,ok              
        write(vname,'("uoc",I2.2)') k
        call histrd1(iarchi,ier,vname,ik,mloin(:,k,1),ifull)
        write(vname,'("voc",I2.2)') k
        call histrd1(iarchi,ier,vname,ik,mloin(:,k,2),ifull)
      end do
    else
      do k=1,ok
        write(vname,'("uoc",I2.2)') k
        call histrd1(iarchi,ier,vname,ik,ucc,6*ik*ik)
        write(vname,'("voc",I2.2)') k
        call histrd1(iarchi,ier,vname,ik,vcc,6*ik*ik)
        call fill_cc(ucc,land_a)
        call fill_cc(vcc,land_a)
        call interpwind(mloin(:,k,1),mloin(:,k,2),ucc,vcc)
      end do
    end if ! iotest
    call mloregrid(ok,ocndwn(:,1),mloin(:,:,1),mlodwn(:,:,3),2)
    call mloregrid(ok,ocndwn(:,1),mloin(:,:,2),mlodwn(:,:,4),3)
  else
    mlodwn(:,:,3:4)=0.               
  end if ! (nestesd/=1.or.nud_ouv/=0) ..else..
  ! water surface height
  if ( nested/=1 .or. nud_sfh/=0 ) then
    call fillhist1('ocheight',ocndwn(:,2),land_a)
  else
    ocndwn(:,2)=0.
  end if ! (nested/=1.or.nud_sfh/=0) ..else..
end if
!--------------------------------------------------------------


!--------------------------------------------------------------
! read sea ice here for prescribed SSTs configuration and for
! mixed-layer-ocean
if ( tsstest ) then
  call histrd1(iarchi,ier,'siced',  ik,sicedep,ifull)
  call histrd1(iarchi,ier,'fracice',ik,fracice,ifull)
else
  call histrd1(iarchi,ier,'siced',  ik,sicedep_a,6*ik*ik)
  call histrd1(iarchi,ier,'fracice',ik,fracice_a,6*ik*ik)
        
  ! diagnose sea-ice if required
  if ( myid==0 ) then
    if ( iers(2)==0 ) then  ! i.e. sicedep read in 
      if (iers(3)/=0 ) then ! i.e. sicedep read in; fracice not read in
        where ( sicedep_a>0. )
          fracice_a=1.
        endwhere
      endif  ! (ierr/=0)  fracice
    else     ! sicedep not read in
      if ( iers(3)/=0 ) then  ! neither sicedep nor fracice read in
        sicedep_a(:)=0.  ! Oct 08
        fracice_a(:)=0.
        write(6,*)'pre-setting siced in onthefly from tss'
        where ( abs(tss_a) <= 271.6 ) ! for ERA-Interim
          sicedep_a=1.  ! Oct 08   ! previously 271.2
          fracice_a=1.
        endwhere
      else  ! i.e. only fracice read in;  done in indata, nestin
            ! but needed here for onthefly (different dims) 28/8/08        
        where ( fracice_a>.01 )
          sicedep_a=2.
        elsewhere
          sicedep_a=0.
          fracice_a=0.
        endwhere
      endif  ! (iers(3)/=0)
    endif    ! (iers(2)/=0) .. else ..    for sicedep
         
    ! fill surface temperature and sea-ice
    tss_l_a=abs(tss_a)
    tss_s_a=abs(tss_a)
    call fill_cc(tss_l_a,sea_a)
    call fill_cc(tss_s_a,land_a)
    call fill_cc(sicedep_a,land_a)
    call fill_cc(fracice_a,land_a)
  end if ! myid==0

  if ( iotest ) then
    ! This case occurs for missing sea-ice data
    if ( myid==0 ) then
      call ccmpi_distribute(zss,zss_a)
      call ccmpi_distribute(tss_l,tss_l_a)
      call ccmpi_distribute(tss_s,tss_s_a)
      call ccmpi_distribute(sicedep,sicedep_a)
      call ccmpi_distribute(fracice,fracice_a)
    else
      call ccmpi_distribute(zss)
      call ccmpi_distribute(tss_l)
      call ccmpi_distribute(tss_s)
      call ccmpi_distribute(sicedep)
      call ccmpi_distribute(fracice)
    end if
!   incorporate other target land mask effects
    where ( land(1:ifull) )
      sicedep=0.
      fracice=0.
      tss=tss_l
    elsewhere
      tss=tss_s
    end where
  else
!   The routine doints4 does the gather, calls ints4 and redistributes
    call doints4(zss_a,     zss)
    call doints4(tss_l_a,   tss_l)
    call doints4(tss_s_a,   tss_s)
    call doints4(fracice_a, fracice)
    call doints4(sicedep_a, sicedep)
!   incorporate other target land mask effects
    where ( land(1:ifull) )
      tss(1:ifull)=tss_l
    elsewhere
      tss(1:ifull)=tss_s   ! no sign switch in CCAM
    end where
    where ( land(1:ifull) .or. sicedep<0.05 ) ! for sflux
      sicedep=0.
      fracice=0.
    end where
  end if ! iotest
end if ! (tsstest) ..else..

! to be depeciated !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!if (nspecial==44.or.nspecial==46) then
!  do iq=1,ifull
!    rlongd=rlongg(iq)*180./pi
!    rlatd=rlatt(iq)*180./pi
!    if (rlatd>=-43..and.rlatd<=-30.) then
!      if (rlongd>=155..and.rlongd<=170.) then
!        tss(iq)=tss(iq)+1.
!      end if
!    end if
!  end do
!end if
!if (nspecial==45.or.nspecial==46) then
!  do iq=1,ifull
!    rlongd=rlongg(iq)*180./pi
!    rlatd=rlatt(iq)*180./pi
!    if (rlatd>=-15..and.rlatd<=-5.) then
!      if (rlongd>=150..and.rlongd<=170.) then
!        tss(iq)=tss(iq)+1.
!      end if
!    end if
!  end do
!end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! -------------------------------------------------------------------
! read atmospheric fields for nested=0 or nested=1.and.nud/=0

! air temperature
! read for nested=0 or nested=1.and.(nud_t/=0.or.nud_p/=0)
if ( nested==0 .or. ( nested==1 .and. nud_test/=0 ) ) then
  if ( iotest ) then
    do k=1,kk
      call histrd4s(iarchi,ier,'temp',ik,k,u_k(:,k),ifull)        !     temperature
    end do
  else
    do k=1,kk
      call histrd4s(iarchi,ier,'temp',ik,k,ucc,6*ik*ik)           !     temperature
      if ( k==levkk .and. myid==0 ) t_a_lev=ucc ! store for psl calculation below
      call doints4(ucc,u_k(:,k))                ! ints4 on source grid
    end do
  end if ! iotest
  call vertint(u_k,t,1,kk,sigin)
else
  t(1:ifull,:)=300.    
end if ! (nested==0.or.(nested==1.and.(nud_t/=0.or.nud_p/=0)))

! winds
! read for nested=0 or nested=1.and.nud_uv/=0
if ( nested==0 .or. ( nested==1 .and. nud_uv/=0 ) ) then
  ! to reduce memory footprint, we now have to alternatively read
  ! u and v.  This is a bit inefficent for disk accessing,
  ! but makes it possible to downscale large grids
  if ( iotest ) then
    do k=1,kk            
      call histrd4s(iarchi,ier,'u',ik,k,u_k(:,k),ifull)           !     u wind component
      call histrd4s(iarchi,ier,'v',ik,k,v_k(:,k),ifull)           !     v wind component
    end do
  else
    do k=1,kk
      call histrd4s(iarchi,ier,'u',ik,k,ucc,6*ik*ik)              !     u wind component
      call histrd4s(iarchi,ier,'v',ik,k,vcc,6*ik*ik)              !     v wind component
      call interpwind(u_k(:,k),v_k(:,k),ucc,vcc)
    end do
!   interpolate all required arrays to new C-C positions
!   don't need to do map factors and Coriolis on target grid
  end if ! iotest
  call vertint(u_k ,u, 3,kk,sigin)
  call vertint(v_k ,v, 4,kk,sigin)
else
  u(1:ifull,:)=0.
  v(1:ifull,:)=0.
end if ! (nested==0.or.(nested==1.and.nud_uv/=0))

! mixing ratio
! read for nested=0 or nested=1.and.nud_q/=0
if ( nested==0 .or. ( nested==1 .and. nud_q/=0 ) ) then
  if ( iers(1)==0 ) then
    call gethist4('mixr',qg,2)                                   !     mixing ratio
  else
    call gethist4('q',qg,2)                                      !     mixing ratio
  end if
else
  qg(1:ifull,:)=qgmin
end if ! (nested==0.or.(nested==1.and.nud_q/=0))

! re-grid surface pressure by mapping to MSLP, interpolating and then map to surface pressure
! requires psl_a, zss, zss_a, t and t_a_lev
if ( nested==0 .or. ( nested==1 .and. nud_test/=0 ) ) then
  if ( .not.iotest ) then
    if ( myid==0 ) then
      ! ucc holds pmsl_a
      call mslpx(ucc,psl_a,zss_a,t_a_lev,sigin(levkk))  ! needs pmsl (preferred)
    end if
    call doints4(ucc,pmsl)
!   invert pmsl to get psl
    call to_pslx(pmsl,psl,zss,t(:,lev),lev)  ! on target grid
  end if ! .not.iotest
end if


!------------------------------------------------------------
! Aerosol data
if ( abs(iaero)>=2 .and. ( nested/=1 .or. nud_aero/=0 ) ) then
  do i=1,naero
    select case(i)
      case(1)
        vname='dms'
      case(2)
        vname='so2'
      case(3)
        vname='so4'
      case(4)
        vname='bco'
      case(5)
        vname='bci'
      case(6)
        vname='oco'
      case(7)
        vname='oci'
      case(8)
        vname='dust1'
      case(9)
        vname='dust2'
      case(10)
        vname='dust3'
      case(11)
        vname='dust4'
      case default
        write(6,*) "ERROR: Unknown aerosol type ",i
        call ccmpi_abort(-1)
    end select
    call gethist4(vname,xtgdwn(:,:,i),5)
  end do
end if


!**************************************************************
! This is the end of reading the nudging arrays
!**************************************************************


!--------------------------------------------------------------
! The following data is only read for initial conditions
if ( nested/=1 ) then

  !------------------------------------------------------------------
  ! check soil variables
  if ( myid==0 .or. pfall ) then
    ierc(7:7+3*ms)=0
    if ( ccycle==0 ) then
      !call ccnf_inq_varid(ncid,'cplant1',idv,tst)
      !if ( tst ) ierc(7)=-1
      ierc(7)=-1
    else
      call ccnf_inq_varid(ncid,'glai',idv,tst)
      if ( tst ) ierc(7)=-1
    end if
    do k=1,ms
      write(vname,'("tgg",I1.1)') k
      call ccnf_inq_varid(ncid,vname,idv,tst)
      if ( tst ) ierc(7+k)=-1
      write(vname,'("wetfrac",I1.1)') k
      call ccnf_inq_varid(ncid,vname,idv,tst)
      if ( tst ) ierc(7+ms+k)=-1
      write(vname,'("wb",I1.1)') k
      call ccnf_inq_varid(ncid,vname,idv,tst)
      if ( tst ) ierc(7+2*ms+k)=-1
    end do
  end if
  
  ! -----------------------------------------------------------------
  ! verify if input is a restart file
  if ( nested==0 ) then
    if ( myid==0 .or. pfall ) then
      if ( kk==kl .and. iotest ) then
        lrestart=.true.
        call ccnf_inq_varid(ncid,'omega',idv,tst)
        if ( tst ) lrestart=.false.
        call ccnf_inq_varid(ncid,'zgnhs',idv,tst)
        if ( tst ) lrestart=.false.
        call ccnf_inq_varid(ncid,'sdot',idv,tst)
        if ( tst ) lrestart=.false.
        call ccnf_inq_varid(ncid,'pslx',idv,tst)
        if ( tst ) lrestart=.false.
        call ccnf_inq_varid(ncid,'savu',idv,tst)
        if ( tst ) lrestart=.false.
        call ccnf_inq_varid(ncid,'savv',idv,tst)
        if ( tst ) lrestart=.false.
        call ccnf_inq_varid(ncid,'savu1',idv,tst)
        if ( tst ) lrestart=.false.
        call ccnf_inq_varid(ncid,'savv1',idv,tst)
        if ( tst ) lrestart=.false.
        call ccnf_inq_varid(ncid,'savu2',idv,tst)
        if ( tst ) lrestart=.false.
        call ccnf_inq_varid(ncid,'savv2',idv,tst)
        if ( tst ) lrestart=.false.
        call ccnf_inq_varid(ncid,'nstag',idv,tst)
        if ( tst ) then
          lrestart=.false.
        else 
          call ccnf_get_var1(ncid,idv,iarchi,ierc(3))
        end if
        call ccnf_inq_varid(ncid,'nstagu',idv,tst)
        if ( tst ) then
          lrestart=.false.
        else 
          call ccnf_get_var1(ncid,idv,iarchi,ierc(4))
        end if
        call ccnf_inq_varid(ncid,'nstagoff',idv,tst)
        if ( tst ) then
          lrestart=.false.
        else 
          call ccnf_get_var1(ncid,idv,iarchi,ierc(5))
        end if
        if ( abs(nmlo)>=3 .and. abs(nmlo)<=9 ) then
          if ( ok==wlev ) then
            call ccnf_inq_varid(ncid,'oldu101',idv,tst)
            if ( tst ) lrestart=.false.
            call ccnf_inq_varid(ncid,'oldv101',idv,tst)
            if ( tst ) lrestart=.false.
            call ccnf_inq_varid(ncid,'oldu201',idv,tst)
            if ( tst ) lrestart=.false.
            call ccnf_inq_varid(ncid,'oldv201',idv,tst)
            if ( tst ) lrestart=.false.
            call ccnf_inq_varid(ncid,'ipice',idv,tst)
            if ( tst ) lrestart=.false.
            call ccnf_inq_varid(ncid,'nstagoffmlo',idv,tst)
            if ( tst ) then
              lrestart=.false.
            else
              call ccnf_get_var1(ncid,idv,iarchi,ierc(6))
            end if
          else
            lrestart=.false.
          end if
        end if
      else
        lrestart=.false.
      end if
      ierc(1:2)=0
      if ( lrestart ) ierc(1)=1
      call ccnf_inq_varid(ncid,'u10',idv,tst)
      if ( tst ) ierc(2)=-1
    end if
    if ( .not.pfall ) then
      call ccmpi_bcast(ierc(1:7+3*ms),0,comm_world)
    end if
    lrestart=(ierc(1)==1)
    if ( lrestart ) then
      nstag      =ierc(3)
      nstagu     =ierc(4)
      nstagoff   =ierc(5)
      nstagoffmlo=ierc(6)
      if ( myid==0 ) then
        write(6,*) "Continue stagging from"
        write(6,*) "nstag,nstagu,nstagoff ",nstag,nstagu,nstagoff
        if ( abs(nmlo)>=3 .and. abs(nmlo)<=9 ) then
          write(6,*) "nstagoffmlo ",nstagoffmlo
        end if
      end if
    end if
  else
    if ( .not.pfall ) then
      call ccmpi_bcast(ierc(7:7+3*ms),0,comm_world)
    end if            
  end if ! nested==0 ..else..
        
  !------------------------------------------------------------------
  ! Read snow and soil tempertaure
  call gethist1('snd',snowd)
  where ( .not.land(1:ifull) .and. (sicedep==0. .or. nmlo==0) )
    snowd=0.
  end where
  do k=1,ms 
    if ( ierc(7+k)==0 ) then
      write(vname,'("tgg",I1.1)') k
    else if ( k<=3 .and. ierc(7+2)==0 ) then
      vname="tgg2"
    else if ( k<=3 ) then
      vname="tb3"
    else if ( ierc(7+6)==0 ) then
      vname="tgg6"
    else
      vname="tb2"
    end if
    if ( iotest ) then
      if ( k==1 .and. ierc(7+1)/=0 ) then
        tgg(:,k)=tss
      else
        call histrd1(iarchi,ier,vname,ik,tgg(:,k),ifull)
      end if
    else
      if ( k==1 .and. ierc(7+1)/=0 ) then
        ucc(1:dk*dk*6)=tss_a(1:dk*dk*6)
      else
        call histrd1(iarchi,ier,vname,ik,ucc,6*ik*ik)
      end if
      call fill_cc(ucc,sea_a)
      call doints4(ucc,tgg(:,k))
    end if
  end do
  if ( .not.iotest ) then
    where ( snowd>0. .and. land(1:ifull) )
      tgg(:,1)=min( tgg(:,1), 270.1 )
    endwhere
  end if

  !--------------------------------------------------
  ! Read MLO sea-ice data
  if ( nmlo/=0 .and. abs(nmlo)<=9 ) then
    if ( .not.allocated(micdwn) ) allocate(micdwn(ifull,11))
    do k=1,7
      select case(k)
        case(1,2,3,4)
          write(vname,'("tggsn",I1.1)') k
          call fillhist1(vname,micdwn(:,k),land_a)
          if ( all(micdwn(:,k)==0.) ) micdwn(:,k)=280.
        case(5)
          micdwn(:,k)=fracice ! read above with nudging arrays
        case(6)
          micdwn(:,k)=sicedep ! read above with nudging arrays
        case(7)
          micdwn(:,k)=snowd*1.E-3
      end select
    end do
    call fillhist1('sto',micdwn(:,8),land_a)
    if ( iotest ) then
      call histrd1(iarchi,ier,'uic',ik,micdwn(:,9),ifull)
      call histrd1(iarchi,ier,'vic',ik,micdwn(:,10),ifull)
    else
      call histrd1(iarchi,ier,'uic',ik,ucc,6*ik*ik)
      call histrd1(iarchi,ier,'vic',ik,vcc,6*ik*ik)
      call fill_cc(ucc,land_a)
      call fill_cc(vcc,land_a)
      call interpwind(micdwn(:,9),micdwn(:,10),ucc,vcc)
    end if ! iotest
    call fillhist1('icesal',micdwn(:,11),land_a)
    if ( abs(nmlo)>=2 ) then
      call gethist1('swater',watbdy)
    end if
  end if

  !------------------------------------------------------------------
  ! Read soil moisture
  wb=20.5
  do k=1,ms
    if ( ierc(7+ms+k)==0 ) then
      write(vname,'("wetfrac",I1.1)') k
    else if ( ierc(7+2*ms+k)==0 ) then
      write(vname,'("wb",I1.1)') k
    else if ( k<2 .and. ierc(7+2*ms+2)==0 ) then
      vname="wb2"
    else if ( k<2 ) then
      vname="wfg"
    else if ( ierc(7+2*ms+6)==0 ) then
      vname="wb6"
    else
      vname="wfb"
    end if
    if ( iotest ) then
      call histrd1(iarchi,ier,vname,ik,wb(:,k),ifull)
      if ( ierc(7+ms+k)==0 ) then
        wb(:,k)=wb(:,k)+20. ! flag for fraction of field capacity
      end if
    else
      call histrd1(iarchi,ier,vname,ik,ucc,6*ik*ik)
      if ( ierc(7+ms+k)==0 ) then
        ucc=ucc+20.
      end if
      call fill_cc(ucc,sea_a)
      call doints4(ucc,wb(:,k))
    end if ! iotest
  end do
  !unpack field capacity into volumetric soil moisture
  if ( any(wb(:,:)>10.) ) then
    if ( mydiag ) write(6,*) "Unpacking wetfrac to wb",wb(idjd,1)
    wb(:,:)=wb(:,:)-20.
    do iq=1,ifull
      isoil=isoilm(iq)
      wb(iq,:)=(1.-wb(iq,:))*swilt(isoil)+wb(iq,:)*sfc(isoil)
    end do
    if ( mydiag ) write(6,*) "giving wb",wb(idjd,1)
  end if
  call fillhist1('wetfac',wetfac,sea_a)
  where ( .not.land )
    wetfac=1.
  end where

  !------------------------------------------------------------------
  ! Read 10m wind speeds for special sea roughness length calculations
  if ( nested==0 ) then
    if ( ierc(2)==0 ) then
      call gethist1('u10',u10)
    else
      u10=sqrt(u(1:ifull,1)**2+v(1:ifull,1)**2)*log(10./0.001)/log(zmin/0.001)
    end if
  end if

  !------------------------------------------------------------------
  ! Read sensible heat flux for convection
  call gethist1('fg',fg)
  
  !------------------------------------------------------------------
  ! Read boundary layer height for TKE-eps mixing and aerosols
  call gethist1('pblh',pblh)
  pblh=max(pblh,1.)
  if ( nvmix==6 .and. nested==0 ) then
    if ( iotest ) then
      call histrd1(iarchi,ier,'dpblh',ik,zidry,ifull)
    else
      zidry=pblh 
    end if ! iotest
    zidry=max(zidry,1.)
  end if

  !------------------------------------------------------------------
  ! Read CABLE/CASA aggregate C+N+P pools
  if ( nsib>=6 ) then
    if ( ccycle==0 ) then
      !if ( ierc(7)==0 ) then
      !  do k=1,ncp
      !    write(vname,'("cplant",I1.1)') k
      !    call fillhist1(vname,cplant(:,k),sea_a)
      !  end do
      !  do k=1,ncs
      !    write(vname,'("csoil",I1.1)') k
      !    call fillhist1(vname,csoil(:,k),sea_a)
      !  end do
      !end if
    else
      if ( ierc(7)==0 ) then
        do k=1,mplant
          write(vname,'("cplant",I1.1)') k
          call fillhist1(vname,cplant(:,k),sea_a)
          write(vname,'("nplant",I1.1)') k
          call fillhist1(vname,niplant(:,k),sea_a)
          write(vname,'("pplant",I1.1)') k
          call fillhist1(vname,pplant(:,k),sea_a)
        end do
        do k=1,mlitter
          write(vname,'("clitter",I1.1)') k
          call fillhist1(vname,clitter(:,k),sea_a)
          write(vname,'("nlitter",I1.1)') k
          call fillhist1(vname,nilitter(:,k),sea_a)
          write(vname,'("plitter",I1.1)') k
          call fillhist1(vname,plitter(:,k),sea_a)
        end do         
        do k=1,msoil
          write(vname,'("csoil",I1.1)') k
          call fillhist1(vname,csoil(:,k),sea_a)
          write(vname,'("nsoil",I1.1)') k
          call fillhist1(vname,nisoil(:,k),sea_a)
          write(vname,'("psoil",I1.1)') k
          call fillhist1(vname,psoil(:,k),sea_a)
        end do
        call fillhist1('glai',glai,sea_a)
      end if ! ierc(7)==0
    end if ! ccycle==0 ..else..
  end if ! if nsib==6.or.nsib==7

  !------------------------------------------------------------------
  ! Read urban data
  if ( nurban/=0 ) then
    if ( .not.allocated(atebdwn) ) allocate(atebdwn(ifull,28))
    do k=1,28
      select case(k)
        case(1)
          vname='rooftgg1'
        case(2)
          vname='rooftgg2'
        case(3)
          vname='rooftgg3'
        case(4)
          vname='rooftgg4'
        case(5)
          vname='waletgg1'
        case(6)
          vname='waletgg2'
        case(7)
          vname='waletgg3'
        case(8)
          vname='waletgg4'
        case(9)
          vname='walwtgg1'
        case(10)
          vname='walwtgg2'
        case(11)
          vname='walwtgg3'
        case(12)
          vname='walwtgg4'
        case(13)
          vname='roadtgg1'
        case(14)
          vname='roadtgg2'
        case(15)
          vname='roadtgg3'
        case(16)
          vname='roadtgg4'
        case(17)
          vname='urbnsmc'
        case(18)
          vname='urbnsmr'
        case(19)
          vname='roofwtr'
        case(20)
          vname='roadwtr'
        case(21)
          vname='urbwtrc'
        case(22)
          vname='urbwtrr'
        case(23)
          vname='roofsnd'
        case(24)
          vname='roadsnd'
        case(25)
          vname='roofden'
        case(26)
          vname='roadden'
        case(27)
          vname='roofsna'
        case(28)
          vname='roadsna'
      end select
      if ( iotest ) then
        call histrd1(iarchi,ier,vname,ik,atebdwn(:,k),ifull)
      else
        call histrd1(iarchi,ier,vname,ik,ucc,6*ik*ik)
        where ( ucc>=399. )
          ucc=999.
        end where
        call fill_cc(ucc,sea_a)
        call doints4(ucc,atebdwn(:,k))
      end if ! iotest
      if ( all(atebdwn(:,k)==0.) ) then
        select case(k)
          case(1:16)
            atebdwn(:,k)=300.
          case(25:26)
            atebdwn(:,k)=100.
          case(27:28)
            atebdwn(:,k)=0.85
        end select
      end if
    end do
  end if

  ! -----------------------------------------------------------------
  ! Read omega for restart
  if ( nested==0 ) then
    ! only for restart - no interpolation
    dpsldt=-999.
    if ( lrestart ) then
      do k=1,kk 
       call histrd4s(iarchi,ier,'omega',ik,k,dpsldt(:,k),ifull)
       dpsldt(:,k)=dpsldt(:,k)/(1.e5*exp(psl(1:ifull)))
      enddo  ! k loop
    end if
  end if

  ! -----------------------------------------------------------------
  ! Read cloud fields
  if ( nested==0 ) then
    call gethist4('qfg',qfg,5)               ! CLOUD FROZEN WATER
    call gethist4('qlg',qlg,5)               ! CLOUD LIQUID WATER
    call gethist4('qrg',qrg,5)               ! RAIN
    call gethist4('cfrac',cfrac,5)           ! CLOUD FRACTION
    call gethist4('rfrac',rfrac,5)           ! RAIN FRACTION
    if ( ncloud>=4 ) then
      call gethist4('stratcf',stratcloud,5)  ! STRAT CLOUD FRACTION
      call gethist4('strat_nt',nettend,5)    ! STRAT NET TENDENCY
    end if ! (ncloud>=4)
  end if   ! (nested==0)

  !------------------------------------------------------------------
  ! TKE-eps data
  if ( nvmix==6 .and. nested==0 ) then
    call gethist4('tke',tke,5)
    if ( all(tke(1:ifull,:)==0.) ) tke(1:ifull,:)=1.5E-4
    call gethist4('eps',eps,5)
    if  (all(eps(1:ifull,:)==0.) ) eps(1:ifull,:)=1.E-7
  end if

  !------------------------------------------------------------------
  ! Tracer data
  if ( ngas>0 ) then              
    do igas=1,ngas              
      write(trnum,'(i3.3)') igas
      call gethist4('tr'//trnum,tr(:,:,igas),7)
    enddo                       
  endif                         

  !------------------------------------------------------------------
  ! Aerosol data ( non-nudged or diagnostic )
  if ( abs(iaero)>=2 ) then
    do i=naero+1,naero+2
      select case(i)
        case(12)
          vname='seasalt1'
        case(13)
          vname='seasalt2'
        case default
          write(6,*) "ERROR: Unknown aerosol type ",i
          call ccmpi_abort(-1)
      end select
      call gethist4(vname,ssn(:,:,i-naero),5)
    end do
    ! Factor 1.e3 to convert to g/m2, x 3 to get sulfate from sulfur
    so4t(:)=0.
    do k=1,kl
      so4t(:)=so4t(:)+3.e3*xtgdwn(:,k,3)*(-1.e5*exp(psl(1:ifull))*dsig(k))/grav
    enddo
  end if
  
  ! -----------------------------------------------------------------
  ! Restart fields
  if ( nested==0 ) then
    ! GEOPOTENTIAL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    phi_nh=0.
    if ( lrestart ) then
      do k=1,kk 
        call histrd4s(iarchi,ier,'zgnhs',ik,k,phi_nh(:,k),ifull)
      enddo  ! k loop
    end if

    ! SDOT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    sdot=-999.
    if ( lrestart ) then
      sdot(:,1)=0.
      do k=1,kk 
        call histrd4s(iarchi,ier,'sdot',ik,k,sdot(:,k+1),ifull)
      enddo  ! k loop
    end if

    ! PSLX !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    pslx=-999.
    if ( lrestart ) then
      do k=1,kk 
       call histrd4s(iarchi,ier,'pslx',ik,k,pslx(:,k),ifull)
      enddo  ! k loop
    end if
          
    ! SAVU !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    savu=-999.
    if ( lrestart ) then
      do k=1,kk 
        call histrd4s(iarchi,ier,'savu',ik,k,savu(:,k),ifull)
      enddo  ! k loop
    end if
          
    ! SAVV !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    savv=-999.
    if ( lrestart ) then
      do k=1,kk 
        call histrd4s(iarchi,ier,'savv',ik,k,savv(:,k),ifull)
      enddo  ! k loop
    end if

    ! SAVU1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    savu1=-999.
    if ( lrestart ) then
      do k=1,kk 
        call histrd4s(iarchi,ier,'savu1',ik,k,savu1(:,k),ifull)
      enddo  ! k loop
    end if
          
    ! SAVV1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    savv1=-999.
    if ( lrestart ) then
      do k=1,kk 
        call histrd4s(iarchi,ier,'savv1',ik,k,savv1(:,k),ifull)
      enddo  ! k loop
    end if

    ! SAVU2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    savu2=-999.
    if ( lrestart ) then
      do k=1,kk 
        call histrd4s(iarchi,ier,'savu2',ik,k,savu2(:,k),ifull)
      enddo  ! k loop
    end if
          
    ! SAVV2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    savv2=-999.
    if ( lrestart ) then
      do k=1,kk 
        call histrd4s(iarchi,ier,'savv2',ik,k,savv2(:,k),ifull)
      enddo  ! k loop
    end if

    ! OCEAN DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if ( abs(nmlo)>=3 .and. abs(nmlo)<=9 ) then
      oldu1=0.
      oldu2=0.
      oldv1=0.
      oldv2=0.
      ipice=0.
      if ( lrestart ) then
        do k=1,ok
          write(vname,'("oldu1",I2.2)') k
          call histrd1(iarchi,ier,vname,ik,oldu1(:,k),ifull)
          write(vname,'("oldv1",I2.2)') k
          call histrd1(iarchi,ier,vname,ik,oldv1(:,k),ifull)
          write(vname,'("oldu2",I2.2)') k
          call histrd1(iarchi,ier,vname,ik,oldu2(:,k),ifull)
          write(vname,'("oldv2",I2.2)') k
          call histrd1(iarchi,ier,vname,ik,oldv2(:,k),ifull)
        end do
        call histrd1(iarchi,ier,'ipice',ik,ipice,ifull)
      end if
    end if
       
  end if ! (nested==0)

  ! -----------------------------------------------------------------
  ! soil ice and snow data
  do k=1,ms
    write(vname,'("wbice",I1.1)') k
    call gethist1(vname,wbice(:,k))    ! SOIL ICE
  end do
  if ( nmlo==0 .or. abs(nmlo)>9 ) then ! otherwise already read above
    do k=1,3
      write(vname,'("tggsn",I1.1)') k
      call gethist1(vname,tggsn(:,k))
      if ( all(tggsn(:,k)==0.) ) tggsn(:,k)=280.
    end do
  end if
  do k=1,3
    write(vname,'("smass",I1.1)') k
    call gethist1(vname,smass(:,k))
  end do
  do k=1,3
    write(vname,'("ssdn",I1.1)') k
    call gethist1(vname,ssdn(:,k))
    if ( all(ssdn(:,k)==0.) ) then
      where ( snowd>100. )
        ssdn(:,k)=240.
      elsewhere
        ssdn(:,k)=140.
      end where
    end if
  end do
  ssdnn=ssdn(:,1)
  call gethist1('snage',snage)
  call gethist1('sflag',dum6)
  isflag=nint(dum6)

  ! -----------------------------------------------------------------
  ! Misc fields
  ! sgsave is needed for convection
  call gethist1('sgsave',sgsave)
        
endif    ! (nested/=1)

!**************************************************************
! This is the end of reading the initial arrays
!**************************************************************         

! -------------------------------------------------------------------
! tgg holds file surface temperature when no MLO
if ( nmlo==0 .or. abs(nmlo)>9 ) then
  where ( .not.land )
    tgg(:,1)=tss
  end where
end if

! -------------------------------------------------------------------
! set-up for next read of file
iarchi=iarchi+1
kdate_s=kdate_r
ktime_s=ktime_r+1

if ( myid==0 .and. nested==0 ) then
  write(6,*) "Final lrestart ",lrestart
end if

return
end subroutine onthefly_work


! *****************************************************************************
! INTERPOLATION ROUTINES                         

! Main interface
! Note that sx is a global array for all processors

! MJT notes - We only need to send sx data for panels contained in
! nface.  Hence we could use RMA to only send the faces that are
! needed for the interpolation.
                         
subroutine doints4(s,sout)
      
use cc_mpi                 ! CC MPI routines
      
implicit none
     
include 'newmpar.h'        ! Grid parameters
include 'parm.h'           ! Model configuration
      
integer :: iq, mm, n, i
integer :: n_n, n_e, n_w, n_s, np1, nm1
integer :: ik2
real, dimension(6*dk*dk), intent(in) :: s
real, dimension(ifull), intent(inout) :: sout
real, dimension(ifull,m_fly) :: wrk
real, dimension(-1:ik+2,-1:ik+2,0:npanels) :: sx ! large common array

call START_LOG(otf_ints_begin)

if ( dk>0 ) then
  ik2 = ik*ik
  !     first extend s arrays into sx - this one -1:il+2 & -1:il+2
  sx(1:ik,1:ik,0:npanels) = reshape(s(1:ik2*6), (/ik,ik,6/))
  do n=0,npanels,2
    n_w=mod(n+5,6)*ik2
    n_e=mod(n+2,6)*ik2
    n_n=mod(n+1,6)*ik2
    n_s=mod(n+4,6)*ik2
    np1=(n+1)*ik2
    do i=1,ik
      sx(0,i,n)   =s(i*ik+n_w)
      sx(-1,i,n)  =s(i*ik-1+n_w)
      sx(ik+1,i,n)=s(ik+1-i+n_e)
      sx(ik+2,i,n)=s(2*ik+1-i+n_e)
      sx(i,ik+1,n)=s(i+np1)
      sx(i,ik+2,n)=s(i+ik+np1)
      sx(i,0,n)   =s((1-i)*ik+ik2+n_s)
      sx(i,-1,n)  =s((1-i)*ik-1+ik2+n_s)
    enddo
    sx(-1,0,n)     =s(2*ik+n_w)         ! wws
    sx(0,-1,n)     =s(ik2-ik+n_s)       ! wss
    sx(0,0,n)      =s(ik+n_w)           ! ws
    sx(ik+1,0,n)   =s(ik+n_e)           ! es  
    sx(ik+2,0,n)   =s(ik-1+n_e)         ! ees 
    sx(-1,ik+1,n)  =s(ik2-ik+n_w)       ! wwn
    sx(0,ik+2,n)   =s(ik2-1+n_w)        ! wnn
    sx(ik+2,ik+1,n)=s(2+n_e)            ! een  
    sx(ik+1,ik+2,n)=s(1+ik+n_e)         ! enn  
    sx(0,ik+1,n)   =s(ik2+n_w)          ! wn  
    sx(ik+1,ik+1,n)=s(1+n_e)            ! en  
    sx(ik+1,-1,n)  =s(2*ik+n_e)         ! ess  
  enddo  ! n loop
  do n=1,npanels,2
    n_w=mod(n+4,6)*ik2
    n_e=mod(n+1,6)*ik2
    n_n=mod(n+2,6)*ik2
    n_s=mod(n+5,6)*ik2
    nm1=(n-1)*ik2
    do i=1,ik
      sx(0,i,n)   =s(1-i+ik2+n_w)
      sx(-1,i,n)  =s(1-i-ik+ik2+n_w)
      sx(ik+1,i,n)=s(1+(i-1)*ik+n_e)
      sx(ik+2,i,n)=s(2+(i-1)*ik+n_e)
      sx(i,ik+1,n)=s(1-i*ik+ik2+n_n)
      sx(i,ik+2,n)=s(2-i*ik+ik2+n_n)
      sx(i,0,n)   =s(i-ik+ik2+nm1)
      sx(i,-1,n)  =s(i-2*ik+ik2+nm1)
    enddo
    sx(-1,0,n)     =s(ik2-1+n_w)       ! wws
    sx(0,-1,n)     =s(2-ik+ik2+n_s)    ! wss
    sx(0,0,n)      =s(ik2+n_w)         ! ws
    sx(ik+1,0,n)   =s(1+n_e)           ! es
    sx(ik+2,0,n)   =s(1+ik+n_e)        ! ees
    sx(-1,ik+1,n)  =s(2-ik+ik2+n_w)    ! wwn   
    sx(0,ik+2,n)   =s(1-2*ik+ik2+n_w)  ! wnn  
    sx(ik+2,ik+1,n)=s(1-2*ik+ik2+n_e)  ! een  
    sx(ik+1,ik+2,n)=s(2-ik+ik2+n_e)    ! enn  
    sx(0,ik+1,n)   =s(1-ik+ik2+n_w)    ! wn  
    sx(ik+1,ik+1,n)=s(1-ik+ik2+n_e)    ! en  
    sx(ik+1,-1,n)  =s(2+n_e)           ! ess  
  enddo  ! n loop
  !     for ew interpolation, sometimes need (different from ns):
  !          (-1,0),   (0,0),   (0,-1)   (-1,il+1),   (0,il+1),   (0,il+2)
  !         (il+1,0),(il+2,0),(il+1,-1) (il+1,il+1),(il+2,il+1),(il+1,il+2)
end if

call ccmpi_bcast(sx,0,comm_world)
  
if ( nord==1 ) then ! bilinear
  do mm=1,m_fly     !  was 4, now may be 1
    call ints_blb(sx,wrk(:,mm),nface4(:,mm),xg4(:,mm),yg4(:,mm))
  enddo
else                ! bicubic
  do mm=1,m_fly  !  was 4, now may be 1
    call intsb(sx,wrk(:,mm),nface4(:,mm),xg4(:,mm),yg4(:,mm))
  enddo
endif   ! (nord==1)  .. else ..

sout=sum(wrk,2)/real(m_fly)

call END_LOG(otf_ints_end)

return
end subroutine doints4

subroutine intsb(sx,sout,nface_l,xg_l,yg_l)   ! N.B. sout here
      
!     same as subr ints, but with sout passed back and no B-S      
!     s is input; sout is output array
!     later may wish to save idel etc between array calls
!     this one does linear interp in x on outer y sides
!     doing x-interpolation before y-interpolation
!     This is a global routine 

implicit none
      
include 'newmpar.h'  ! Grid parameters
include 'parm.h'     ! Model configuration

integer, dimension(ifull), intent(in) :: nface_l
integer :: idel, jdel, nn
integer :: n, iq
real, dimension(ifull), intent(inout) :: sout
real aaa, c1, c2, c3, c4, xxg, yyg
real, intent(in), dimension(ifull) :: xg_l, yg_l
real, dimension(-1:ik+2,-1:ik+2,0:npanels), intent(in) :: sx
real, dimension(4) :: r

!     this is intsb           EW interps done first

do iq=1,ifull   ! runs through list of target points
  n=nface_l(iq)
  idel=int(xg_l(iq))
  xxg=xg_l(iq)-idel
! yg here goes from .5 to il +.5
  jdel=int(yg_l(iq))
  yyg=yg_l(iq)-jdel
  do nn=2,3       ! N.B.
    c1=sx(idel-1,jdel+nn-2,n)
    c2=sx(idel  ,jdel+nn-2,n)
    c3=sx(idel+1,jdel+nn-2,n)
    c4=sx(idel+2,jdel+nn-2,n)
    r(nn)=((1.-xxg)*((2.-xxg)*((1.+xxg)*c2-xxg*c1/3.)-xxg*(1.+xxg)*c4/3.)+xxg*(1.+xxg)*(2.-xxg)*c3)/2.
  enddo    ! nn loop
!       r       ={(1-x     )*{(2-x     )*[(1+x     )*c2-x     *c1/3]
!         -x     *(1+x     )*c4/3}
!         +x    *(1+x     )*(2-x     )*c3}/2
  do nn=1,4,3       ! N.B.
    c2=sx(idel  ,jdel+nn-2,n)
    c3=sx(idel+1,jdel+nn-2,n)
    r(nn)=(1.-xxg)*c2 +xxg*c3
  enddo    ! nn loop
!       array(iq)=((1.-yyg)*((2.-yyg)*((1.+yyg)*r(2)-yyg*r(1)/3.)
!    .             -yyg*(1.+yyg)*r(4)/3.)
!    .             +yyg*(1.+yyg)*(2.-yyg)*r(3))/2.
!      following does Bermejo Staniforth
  aaa=((1.-yyg)*((2.-yyg)*((1.+yyg)*r(2)-yyg*r(1)/3.)-yyg*(1.+yyg)*r(4)/3.)+yyg*(1.+yyg)*(2.-yyg)*r(3))/2.
  aaa=min( aaa , max( sx(idel,jdel,n),sx(idel+1,jdel,n),sx(idel,jdel+1,n),sx(idel+1,jdel+1,n) ) )
  sout(iq)=max( aaa , min( sx(idel,jdel,n),sx(idel+1,jdel,n),sx(idel,jdel+1,n),sx(idel+1,jdel+1,n) ) )
enddo    ! iq loop

return
end subroutine intsb

subroutine ints_blb(sx,sout,nface_l,xg_l,yg_l) 
      
!     this one does bi-linear interpolation only

implicit none
      
include 'newmpar.h'  ! Grid parameters
include 'parm.h'     ! Model configuration

integer :: n, iq, idel, jdel
integer, intent(in), dimension(ifull) :: nface_l
real, dimension(ifull), intent(inout) :: sout
real, intent(in), dimension(ifull) :: xg_l, yg_l
real, dimension(-1:ik+2,-1:ik+2,0:npanels), intent(in) :: sx
real :: xxg, yyg

do iq=1,ifull  ! runs through list of target points
  n=nface_l(iq)
  idel=int(xg_l(iq))
  xxg=xg_l(iq)-idel
  jdel=int(yg_l(iq))
  yyg=yg_l(iq)-jdel
  sout(iq)=yyg*(xxg*sx(idel+1,jdel+1,n)+(1.-xxg)*sx(idel,jdel+1,n))+(1.-yyg)*(xxg*sx(idel+1,jdel,n)+(1.-xxg)*sx(idel,jdel,n))
enddo    ! iq loop

return
end subroutine ints_blb

! *****************************************************************************
! FILL ROUTINES

subroutine fill_cc(a_io,land_a)
      
! routine fills in interior of an array which has undefined points

use cc_mpi          ! CC MPI routines

implicit none

integer :: nrem, i, ii, iq, j, n, neighb
integer :: iminb, imaxb, jminb, jmaxb
integer, dimension(0:5) :: imin,imax,jmin,jmax
real, parameter :: value=999.       ! missing value flag
real, dimension(6*dk*dk) :: a_io, a
logical, dimension(6*dk*dk) :: land_a
logical, dimension(4) :: mask

! only perform fill on myid==0
if ( dk==0 ) return

where ( land_a )
  a_io=value
end where
if ( all(abs(a_io-value)<1.E-6) ) return

call START_LOG(otf_fill_begin)

imin=1
imax=dk
jmin=1
jmax=dk
          
nrem = 1    ! Just for first iteration
do while ( nrem > 0)
  nrem=0
  do iq=1,dk*dk*6
    a(iq)=a_io(iq)
  enddo
  ! MJT restricted fill
  do n=0,5
    iminb=dk
    imaxb=1
    jminb=dk
    jmaxb=1
    do j=jmin(n),jmax(n)
      do i=imin(n),imax(n)
        iq=indx(i,j,n,dk,dk)
        if ( a(iq)==value ) then
          ! ic(iq)[1] is in_g for host grid
          ! ic(iq)[2] is is_g for host grid
          ! ic(iq)[3] is ie_g for host grid
          ! ic(iq)[4] is iw_g for host grid
          mask=a(ic(iq))/=value
          neighb=count(mask)
          if ( neighb>0 ) then
            a_io(iq)=sum(a(ic(iq)),mask)/real(neighb)
          else
            iminb=min(i,iminb)
            imaxb=max(i,imaxb)
            jminb=min(j,jminb)
            jmaxb=max(j,jmaxb)
            nrem=nrem+1   ! current number of points without a neighbour
          endif
        endif
      end do
    end do
    imin(n)=iminb
    imax(n)=imaxb
    jmin(n)=jminb
    jmax(n)=jmaxb
  end do
end do
      
call END_LOG(otf_fill_end)
      
return
end subroutine fill_cc

function ic(iq) result (iqq)

implicit none

integer, intent(in) :: iq
integer, dimension(4) :: iqq
integer i, j, n
integer, parameter, dimension(0:5) :: npann=(/1,103,3,105,5,101/)
integer, parameter, dimension(0:5) :: npane=(/102,2,104,4,100,0/)
integer, parameter, dimension(0:5) :: npanw=(/5,105,1,101,3,103/)
integer, parameter, dimension(0:5) :: npans=(/104,0,100,2,102,4/)

n=(iq-1)/(dk*dk)
j=(iq-1-n*dk*dk)/dk+1
i=iq-(j-1)*dk-n*dk*dk
! North
if (j==dk) then
  if (npann(n)<100) then
    iqq(1)=i+npann(n)*dk*dk
  else
    iqq(1)=1+(dk-i)*dk+(npann(n)-100)*dk*dk
  end if
else
  iqq(1)=iq+dk
end if
! South
if (j==1) then
  if (npans(n)<100) then
    iqq(2)=i+(dk-1)*dk+npans(n)*dk*dk
  else
    iqq(2)=dk+(dk-i)*dk+(npans(n)-100)*dk*dk
  end if
else
  iqq(2)=iq-dk
end if
! East
if (i==dk) then
  if (npane(n)<100) then
    iqq(3)=1+(j-1)*dk+npane(n)*dk*dk
  else
    iqq(3)=dk+1-j+(npane(n)-100)*dk*dk
  end if
else
  iqq(3)=iq+1
end if
! West
if (i==1) then
  if (npanw(n)<100) then
    iqq(4)=dk+(j-1)*dk+npanw(n)*dk*dk
  else
    iqq(4)=dk+1-j+(dk-1)*dk+(npanw(n)-100)*dk*dk
  end if
else
  iqq(4)=iq-1
end if

end function ic

! *****************************************************************************
! OROGRAPHIC ADJUSTMENT ROUTINES

subroutine mslpx(pmsl,psl,zs,t,siglev)
      
use cc_mpi, only : mydiag ! CC MPI routines
use sigs_m                ! Atmosphere sigma levels
      
!     this one will ignore negative zs (i.e. over the ocean)
implicit none
      
include 'newmpar.h'       ! Grid parameters
include 'const_phys.h'    ! Physical constants
include 'parm.h'          ! Model configuration

integer iq
real siglev, c, con, conr
real, dimension(6*dk*dk) :: pmsl, psl, zs, t
real, dimension(6*dk*dk) :: dlnps, phi1, tav, tsurf

c=grav/stdlapse
conr=c/rdry
con=siglev**(rdry/c)/c
phi1(:)=t(:)*rdry*(1.-siglev)/siglev ! phi of sig(lev) above sfce
tsurf(:)=t(:)+phi1(:)*stdlapse/grav
tav(:)=tsurf(:)+max(0.,zs(:))*.5*stdlapse/grav
dlnps(:)=max(0.,zs(:))/(rdry*tav(:))
pmsl(:)=1.e5*exp(psl(:)+dlnps(:))

return
end subroutine mslpx
      
subroutine to_pslx(pmsl,psl,zs,t,lev)

use cc_mpi, only : mydiag  ! CC MPI routines
use sigs_m                 ! Atmosphere sigma levels
      
implicit none
      
include 'newmpar.h'        ! Grid parameters
include 'const_phys.h'     ! Physical constants
include 'parm.h'           ! Model configuration
      
integer iq, lev
real, dimension(ifull) :: pmsl, psl, zs, t
real, dimension(ifull) :: dlnps, phi1, tav, tsurf

phi1(:)=t(:)*rdry*(1.-sig(lev))/sig(lev) ! phi of sig(lev) above sfce
tsurf(:)=t(:)+phi1(:)*stdlapse/grav
tav(:)=tsurf(:)+max(0.,zs(:))*.5*stdlapse/grav
dlnps(:)=max(0.,zs(:))/(rdry*tav(:))
psl(:)=log(1.e-5*pmsl(:)) -dlnps(:)

#ifdef debug
if ( nmaxpr==1 .and. mydiag ) then
  write(6,*)'to_psl lev,sig(lev) ',lev,sig(lev)
  write(6,*)'zs,t_lev,psl,pmsl ',zs(idjd),t(idjd),psl(idjd),pmsl(idjd)
endif
#endif

return
end subroutine to_pslx

subroutine retopo(psl,zsold,zs,t,qg)
!     in Jan 07, renamed recalc of ps from here, to reduce confusion     
!     this version (Aug 2003) allows -ve zsold (from spectral model),
!     but assumes new zs is positive for atmospheric purposes
!     this routine redefines psl, t to compensate for zsold going to zs
!     (but does not overwrite zs, ps themselves here)
!     called by indata and nestin for newtop>=1
!     nowadays just for ps and atmospheric fields Mon  08-23-1999
use cc_mpi, only : mydiag
use diag_m
use sigs_m

implicit none

include 'newmpar.h'
include 'const_phys.h'
include 'parm.h'

real, dimension(ifull), intent(inout) :: psl
real, dimension(ifull), intent(in) :: zsold, zs
real, dimension(:,:), intent(inout) :: t, qg
real, dimension(ifull) :: psnew, psold, pslold
real, dimension(kl) :: told, qgold
real sig2
integer iq, k, kkk

pslold(1:ifull)=psl(1:ifull)
psold(1:ifull)=1.e5*exp(psl(1:ifull))
psl(1:ifull)=psl(1:ifull)+(zsold(1:ifull)-zs(1:ifull))/(rdry*t(1:ifull,1))
psnew(1:ifull)=1.e5*exp(psl(1:ifull))

!     now alter temperatures to compensate for new topography
if(ktau<100.and.mydiag)then
  write(6,*) 'retopo: zsold,zs,psold,psnew ',zsold(idjd),zs(idjd),psold(idjd),psnew(idjd)
  write(6,*) 'retopo: old t ',(t(idjd,k),k=1,kl)
  write(6,*) 'retopo: old qg ',(qg(idjd,k),k=1,kl)
endif  ! (ktau.lt.100)
do iq=1,ifull
  do k=1,kl
    qgold(k)=qg(iq,k)
    told(k)=t(iq,k)
  enddo  ! k loop
  do k=1,kl-1
    sig2=sig(k)*psnew(iq)/psold(iq)
    if(sig2>=sig(1))then
!     assume 6.5 deg/km, with dsig=.1 corresponding to 1 km
      t(iq,k)=told(1)+(sig2-sig(1))*6.5/.1  
    else
      do kkk=2,kl
        if(sig2>sig(kkk)) exit
      end do
      t(iq,k)=(told(kkk)*(sig(kkk-1)-sig2)+told(kkk-1)*(sig2-sig(kkk)))/(sig(kkk-1)-sig(kkk))
      qg(iq,k)=(qgold(kkk)*(sig(kkk-1)-sig2)+qgold(kkk-1)*(sig2-sig(kkk)))/(sig(kkk-1)-sig(kkk))
    endif
  enddo  ! k loop
enddo   ! iq loop
if(ktau<100.and.mydiag)then
  write(6,*) 'retopo: new t ',(t(idjd,k),k=1,kl)
  write(6,*) 'retopo: new qg ',(qg(idjd,k),k=1,kl)
endif  ! (ktau.lt.100)

return
end subroutine retopo

! *****************************************************************************
! VECTOR INTERPOLATION ROUTINES

subroutine interpwind(uct,vct,ucc,vcc)
      
use cc_mpi           ! CC MPI routines
use vecsuv_m         ! Map to cartesian coordinates
      
implicit none
      
include 'newmpar.h'  ! Grid parameters
      
integer iq
real, dimension(6*dk*dk), intent(inout) :: ucc, vcc
real, dimension(6*dk*dk) :: wcc
real, dimension(ifull), intent(out) :: uct, vct
real, dimension(ifull) :: wct
real uc, vc, wc, newu, newv, neww

! dk is only non-zero on myid==0
do iq=1,6*dk*dk
  ! first set up winds in Cartesian "source" coords            
  uc=axs_a(iq)*ucc(iq) + bxs_a(iq)*vcc(iq)
  vc=ays_a(iq)*ucc(iq) + bys_a(iq)*vcc(iq)
  wc=azs_a(iq)*ucc(iq) + bzs_a(iq)*vcc(iq)
  ! now convert to winds in "absolute" Cartesian components
  ucc(iq)=uc*rotpoles(1,1)+vc*rotpoles(1,2)+wc*rotpoles(1,3)
  vcc(iq)=uc*rotpoles(2,1)+vc*rotpoles(2,2)+wc*rotpoles(2,3)
  wcc(iq)=uc*rotpoles(3,1)+vc*rotpoles(3,2)+wc*rotpoles(3,3)
end do  ! iq loop
! interpolate all required arrays to new C-C positions
! do not need to do map factors and Coriolis on target grid
call doints4(ucc, uct)
call doints4(vcc, vct)
call doints4(wcc, wct)
do iq=1,ifull
  ! now convert to "target" Cartesian components (transpose used)
  newu=uct(iq)*rotpole(1,1)+vct(iq)*rotpole(2,1)+wct(iq)*rotpole(3,1)
  newv=uct(iq)*rotpole(1,2)+vct(iq)*rotpole(2,2)+wct(iq)*rotpole(3,2)
  neww=uct(iq)*rotpole(1,3)+vct(iq)*rotpole(2,3)+wct(iq)*rotpole(3,3)
  ! then finally to "target" local x-y components
  uct(iq) = ax(iq)*newu + ay(iq)*newv + az(iq)*neww
  vct(iq) = bx(iq)*newu + by(iq)*newv + bz(iq)*neww
end do  ! iq loop
      
return
end subroutine interpwind

! *****************************************************************************
! FILE IO ROUTINES

! This version reads and interpolates a surface field
subroutine gethist1(vname,varout)
      
use infile             ! Input file routines
      
implicit none
      
include 'newmpar.h'    ! Grid parameters
include 'darcdf.h'     ! Netcdf data

integer ier
real, dimension(:), intent(out) :: varout
real, dimension(6*dk*dk) :: ucc
character(len=*), intent(in) :: vname
      
if ( iotest ) then
  call histrd1(iarchi,ier,vname,ik,varout,ifull)
else
  call histrd1(iarchi,ier,vname,ik,ucc,6*ik*ik)
  call doints4(ucc,varout)
end if ! iotest

return
end subroutine gethist1

! This version reads, fills and interpolates a surface field
subroutine fillhist1(vname,varout,mask_a)
      
use infile             ! Input file routines
      
implicit none
      
include 'newmpar.h'    ! Grid parameters
include 'darcdf.h'     ! Netcdf data
      
integer ier
real, dimension(:), intent(out) :: varout
real, dimension(6*dk*dk) :: ucc
logical, dimension(6*dk*dk), intent(in) :: mask_a
character(len=*), intent(in) :: vname
      
if ( iotest ) then
  call histrd1(iarchi,ier,vname,ik,varout,ifull)
else
  call histrd1(iarchi,ier,vname,ik,ucc,6*ik*ik)
  call fill_cc(ucc,mask_a)
  call doints4(ucc,varout)
end if ! iotest
      
return
end subroutine fillhist1

! This version reads and interpolates 3D atmospheric fields
subroutine gethist4(vname,varout,vmode)
      
use infile             ! Input file routines
      
implicit none
      
include 'newmpar.h'    ! Grid parameters
include 'darcdf.h'     ! Netcdf data
      
integer, intent(in) :: vmode
integer k, ier
real, dimension(:,:), intent(out) :: varout
real, dimension(6*dk*dk) :: ucc
real, dimension(ifull,kk) :: u_k
character(len=*), intent(in) :: vname
      
if ( iotest ) then
  do k=1,kk          
    call histrd4s(iarchi,ier,vname,ik,k,u_k(:,k),ifull)
  end do
else
  do k=1,kk
    call histrd4s(iarchi,ier,vname,ik,k,ucc,6*ik*ik)
    call doints4(ucc,u_k(:,k))
  end do
end if ! iotest
call vertint(u_k,varout,vmode,kk,sigin)
      
return
end subroutine gethist4   

end module onthefly_m
