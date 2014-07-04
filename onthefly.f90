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

integer, parameter :: nord = 3        ! 1 for bilinear, 3 for bicubic interpolation
integer, save :: ik, jk, kk, ok, nsibx
integer, dimension(:,:), allocatable, save :: nface4
real, save :: rlong0x, rlat0x, schmidtx
real, dimension(3,3), save :: rotpoles, rotpole
real, dimension(:,:), allocatable, save :: xg4, yg4
real, dimension(:), allocatable, save :: axs_a, ays_a, azs_a
real, dimension(:), allocatable, save :: bxs_a, bys_a, bzs_a
real, dimension(:), allocatable, save :: sigin
logical iotest, newfile
    
contains
    
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
integer, save :: ncidold = -1
integer ier,mtimer,k,ierx,idvkd,idvkt,idvmt
integer, dimension(nihead) :: nahead
integer, dimension(ifull), intent(out) :: isflag
real timer
real, dimension(ifull,wlev,4), intent(out) :: mlodwn
real, dimension(ifull,kl,naero), intent(out) :: xtgdwn
real, dimension(ifull,ms), intent(out) :: wb,wbice,tgg
real, dimension(ifull,3), intent(out) :: tggsn,smass,ssdn
real, dimension(ifull,2), intent(out) :: ocndwn
real, dimension(:,:), intent(out) :: t,u,v,qg,qfg,qlg,qrg
real, dimension(ifull), intent(out) :: psl,zss,tss,fracice,snowd,sicedep
real, dimension(ifull), intent(out) :: ssdnn,snage
real, dimension(nrhead) :: ahead
real, dimension(14) :: rdum
logical ltest,tst

call START_LOG(onthefly_begin)
!--------------------------------------------------------------
! pfall indicates all processors have an input file and there
! is no need to broadcast metadata (see infile.f90)
if ( myid==0 .or. pfall ) then
  if ( myid==0 ) write(6,*) 'Entering onthefly for nested,ktau = ',nested,ktau
  if ( ncid/=ncidold ) then
    if ( myid==0 ) write(6,*) 'Reading new file metadata'
    iarchi=1
    call ccnf_get_attg(ncid,'int_header',nahead)
    call ccnf_get_attg(ncid,'real_header',ahead)
    ik      =nahead(1)
    jk      =nahead(2)
    kk      =nahead(3)
    nsibx   =nahead(44)
    rlong0x =ahead(5)
    rlat0x  =ahead(6)
    schmidtx=ahead(7)
    if ( schmidtx<=0. .or. schmidtx>1. ) then
      ! backwards compatibility option
      rlong0x =ahead(6)
      rlat0x  =ahead(7)
      schmidtx=ahead(8)
    endif  ! (schmidtx<=0..or.schmidtx>1.)        
    maxarchi=0
    call ccnf_inq_dimlen(ncid,'time',maxarchi)
    ok=0
    call ccnf_inq_dimlen(ncid,'olev',ok,failok=.true.)
    if ( myid==0 ) then
      write(6,*) "Found ik,jk,kk,ok ",ik,jk,kk,ok
      write(6,*) "      maxarchi ",maxarchi
      write(6,*) "      rlong0x,rlat0x,schmidtx ",rlong0x,rlat0x,schmidtx
    end if
  end if
  if ( myid==0 ) write(6,*)'Search for kdate_s,ktime_s >= ',kdate_s,ktime_s
  ltest=.true.
  iarchi=iarchi-1
  call ccnf_inq_varid(ncid,'kdate',idvkd,tst)
  call ccnf_inq_varid(ncid,'ktime',idvkt,tst)
  call ccnf_inq_varid(ncid,'mtimer',idvmt,tst)
  if ( tst ) then
    ierx=1
    call ccnf_inq_varid(ncid,'timer',idvmt,tst)
  else
    ierx=0
  end if
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
  ! store metadata for possible bcast
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
endif  ! ( myid==0 .or. pfall )

! if metadata is not read by all processors, then broadcast
if ( .not.pfall ) then
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
  newfile=(ncid/=ncidold)            
end if
      
! close old file if a new file is opened
if ( newfile ) then
  if ( ncidold/=-1 ) then
    if ( myid==0 ) then
      write(6,*) 'Closing old input file'
    end if
    call histclose
  end if
  ncidold=ncid
end if

! trap error if correct date/time is not located      
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
!--------------------------------------------------------------
      
! Here we call ontheflyx with different automatic array
! sizes.  This means the arrays are correct for interpolation
! and file i/o on myid==0, as well as the arrays are smaller
! on myid/=0 when they are not needed.  This way we avoid
! having to maintain multiple ontheflyx subroutines.
      
! Note that if histrd fails to find a variable, it returns
! zero in the output array
      
if ( myid==0 ) then
  call ontheflyx(nested,kdate_r,ktime_r,psl,zss,tss,sicedep,fracice,t,u,v,qg,tgg,wb,wbice,snowd,qfg,qlg, &
                 qrg,tggsn,smass,ssdn,ssdnn,snage,isflag,ik,mlodwn,ocndwn,xtgdwn)
  write(6,*) "Leaving onthefly"
else
  call ontheflyx(nested,kdate_r,ktime_r,psl,zss,tss,sicedep,fracice,t,u,v,qg,tgg,wb,wbice,snowd,qfg,qlg, &
                 qrg,tggsn,smass,ssdn,ssdnn,snage,isflag,0,mlodwn,ocndwn,xtgdwn)
end if

call END_LOG(onthefly_end)

return
end subroutine onthefly
      
! Read data from netcdf file
      
! arrays are typically read as global and then distributed to
! processor local arrays.  This allows for more flexibility
! with diagnosed fields.  Data is usually read in as 2D
! fields which avoids memory problems when the host grid
! size is significantly larger than the regional grid size.
subroutine ontheflyx(nested,kdate_r,ktime_r,psl,zss,tss,sicedep,fracice,t,u,v,qg,tgg,wb,wbice,snowd,qfg,qlg,  &
                     qrg,tggsn,smass,ssdn,ssdnn,snage,isflag,dk,mlodwn,ocndwn,xtgdwn)
      
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
      
integer, intent(in) :: kdate_r, ktime_r, nested, dk
integer idv, isoil, nud_test
integer lev, levkk, ier, ierr, igas
integer nemi, id2, jd2, idjd2
integer i, j, k, mm, iq, ii, jj, np, numneg
integer, dimension(:), allocatable, save :: isoilm_a
integer, dimension(ifull), intent(out) :: isflag
integer, dimension(7+3*ms) :: ierc
integer, dimension(3), save :: iers
integer, dimension(1) :: str, cnt
real(kind=8), dimension(:,:), allocatable, save :: xx4,yy4
real(kind=8), dimension(dk*dk*6):: z_a,x_a,y_a
real, dimension(ifull,wlev,4), intent(out) :: mlodwn
real, dimension(ifull,kl,naero), intent(out) :: xtgdwn
real, dimension(ifull,ok,2) :: mloin
real, dimension(ifull,2), intent(out) :: ocndwn
real, dimension(ifull,ms), intent(out) :: wb,wbice,tgg
real, dimension(ifull,3), intent(out) :: tggsn,smass,ssdn
real, dimension(:,:), intent(out) :: t,u,v,qg,qfg,qlg,qrg
real, dimension(ifull,kk) :: u_k,v_k
real, dimension(ifull), intent(out) :: psl,zss,tss,fracice
real, dimension(ifull), intent(out) :: snowd,sicedep,ssdnn,snage
real, dimension(ifull) :: dum6
real, dimension(ifull) :: tss_l, tss_s, pmsl
real, dimension(ik*ik*6) :: ucc,vcc
real, dimension(ik*ik*6) :: fracice_a,sicedep_a
real, dimension(ik*ik*6) :: tss_l_a,tss_s_a
real, dimension(dk*dk*6) :: t_a_lev,psl_a,tss_a
real, dimension(dk*dk*6) :: wts_a  ! not used here or defined in call setxyz
real, dimension(:), allocatable, save :: zss_a,ocndep_l
real, dimension(kk+3) :: dumr
real rlongd, rlatd
character(len=8) vname
character(len=3) trnum
logical, dimension(:), allocatable, save :: land_a,sea_a
logical tsstest,tst

! land-sea mask method (nemi=3 use soilt, nemi=2 use tgg, nemi=1 use zs)
nemi=3
      
! test if retopo fields are required
if ( nud_p==0 .and. nud_t==0 .and. nud_q==0 ) then
  nud_test=0
else
  nud_test=1
end if
      
! Determine if interpolation is required
iotest=6*ik*ik==ifull_g .and. abs(rlong0x-rlong0)<iotol .and. abs(rlat0x-rlat0)<iotol .and. &
       abs(schmidtx-schmidt)<iotol .and. nsib==nsibx
if ( iotest ) then
  io_in=1   ! no interpolation
else
  io_in=-1  ! interpolation
end if
if ( myid==0 ) write(6,*) "Interpolation iotest,io_in =",iotest,io_in

!--------------------------------------------------------------
! Allocate interpolation, vertical level and mask arrays
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
      
!--------------------------------------------------------------
! Determine input grid coordinates and interpolation arrays
if ( newfile .and. .not.iotest ) then
  ! xx4 and yy4 could be replaced with sharded arrays in MPI-3
  allocate(xx4(1+4*ik,1+4*ik),yy4(1+4*ik,1+4*ik))

  if ( m_fly==1 ) then
    rlong4_l(:,1)=rlongg(:)*180./pi
    rlat4_l(:,1)=rlatt(:)*180./pi
  end if    
          
  if ( myid==0 ) then
    write(6,*) "Defining input file grid"
!   N.B. -ve ik in call setxyz preserves TARGET rlat4, rlong4     
!   following setxyz call is for source data geom    ****   
    do iq=1,ik*ik*6
      axs_a(iq)=iq
      ays_a(iq)=iq
      azs_a(iq)=iq
    enddo      
    call setxyz(ik,rlong0x,rlat0x,-schmidtx,x_a,y_a,z_a,wts_a,axs_a,ays_a,azs_a,bxs_a,bys_a,bzs_a,xx4,yy4,myid)
  end if ! (myid==0)

  call ccmpi_bcastr8(xx4,0,comm_world)
  call ccmpi_bcastr8(yy4,0,comm_world)
     
  rotpoles = calc_rotpole(rlong0x,rlat0x)
  if ( ktau<3 .and. myid==0 ) then
    write(6,*)'m_fly,nord ',m_fly,nord
    write(6,*)'kdate_r,ktime_r,ktau,ds',kdate_r,ktime_r,ktau,ds
    write(6,*)'rotpoles:'
    do i=1,3
      write(6,'(3x,2i1,5x,2i1,5x,2i1,5x,3f8.4)') (i,j,j=1,3),(rotpoles(i,j),j=1,3)
    enddo
  endif                  ! (ktau<3.and.myid==0)

! rotpole(1,) is x-axis of rotated coords in terms of orig Cartesian
! rotpole(2,) is y-axis of rotated coords in terms of orig Cartesian
! rotpole(3,) is z-axis of rotated coords in terms of orig Cartesian
  rotpole = calc_rotpole(rlong0,rlat0)
  if ( nmaxpr==1 .and. myid==0 ) then
    write(6,*)'in onthefly rotpole:'
    do i=1,3
      write(6,'(3x,2i1,5x,2i1,5x,2i1,5x,3f8.4)') (i,j,j=1,3),(rotpole(i,j),j=1,3)
    enddo
    write(6,*)'xx4,yy4 ',xx4(id,jd),yy4(id,jd)
    write(6,*)'before latltoij for id,jd: ',id,jd
    write(6,*)'rlong0x,rlat0x,schmidtx ',rlong0x,rlat0x,schmidtx
  endif                  ! (nmaxpr==1.and.myid==0)

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
      
       
! read time invariant data when file is first opened
! need global zss_a for (potentially) landsea mask and psl interpolation
! need global isoilm_a for (potentially) landsea mask
if ( newfile ) then
  if ( myid==0 .or. pfall ) then
    if ( myid==0 ) write(6,*) "Reading time invariant fields"
    call ccnf_inq_varid(ncid,'lev',idv,tst)
    if ( tst ) then
      call ccnf_inq_varid(ncid,'layer',idv,tst)
    end if
    if ( tst ) then
      call ccnf_get_attg(ncid,'sigma',sigin)
    else
      str(1)=1
      cnt(1)=kk
      call ccnf_get_vara(ncid,idv,str(1:1),cnt(1:1),sigin)
    end if
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
  if ( allocated(isoilm_a) ) deallocate(isoilm_a)
  if ( tsstest ) then
    ! load local surface temperature
    allocate(zss_a(ifull))
    call histrd1(iarchi,ier,'zht',ik,zss_a,ifull)
  else
    ! load global surface temperature   
    allocate(zss_a(6*ik*ik)) ! allocate for all processors to avoid memory copies
                             ! can replace with allocate on myid==0 and then copy
                             ! to ucc for doints4 if memory becomes a problem
    if ( myid==0 ) allocate(isoilm_a(6*ik*ik))
    call histrd1(iarchi,ier,'zht',  ik,zss_a,6*ik*ik)
    call histrd1(iarchi,ier,'soilt',ik,ucc,  6*ik*ik)
    if ( myid==0 ) then
      isoilm_a=nint(ucc)
      if ( all(isoilm_a==0) ) isoilm_a=-1 ! missing value flag
    end if
  end if
  if ( nmlo/=0 .and. abs(nmlo)<=9 ) then
    if ( .not.allocated(ocndep_l) ) allocate(ocndep_l(ifull))
    call gethist1('ocndepth',ocndep_l)
  end if
  if ( myid==0 ) write(6,*) "Finished reading fixed fields"
else
  if ( myid==0 ) write(6,*) "Using saved fixed fields"
  tsstest=(iers(2)==0.and.iers(3)==0.and.iotest)        
endif ! newfile ..else..
      
! detemine the level below sig=0.9 (used to calculate psl)
lev=0
do while( sig(lev+1)>0.9 ) ! nested grid
  lev=lev+1
end do
levkk=0
do while( sigin(levkk+1)>0.9 ) ! host grid
  levkk=levkk+1
end do      
if (myid==0) write(6,*) "Ref height lev,levkk =",lev,levkk

!--------------------------------------------------------------
! Begin reading host data for current time step
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
      do iq=1,ik*ik*6
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
    deallocate(isoilm_a)
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
      call filhist1(vname,dk,mloin(:,k,1),land_a)
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
      call filhist1(vname,dk,mloin(:,k,1),land_a)
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
        call fill_cc(ucc,dk,0,land_a)
        call fill_cc(vcc,dk,0,land_a)
        call interpwind(mloin(:,k,1),mloin(:,k,2),ucc,vcc,dk)
      end do
    end if ! iotest
    call mloregrid(ok,ocndwn(:,1),mloin(:,:,1),mlodwn(:,:,3),2)
    call mloregrid(ok,ocndwn(:,1),mloin(:,:,2),mlodwn(:,:,4),3)
  else
    mlodwn(:,:,3:4)=0.               
  end if ! (nestesd/=1.or.nud_ouv/=0) ..else..
  ! water surface height
  if ( nested/=1 .or. nud_sfh/=0 ) then
    call filhist1('ocheight',dk,ocndwn(:,2),land_a)
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
    call fill_cc(tss_l_a,ik,0,sea_a)
    call fill_cc(tss_s_a,ik,0,land_a)
    call fill_cc(sicedep_a,ik,0,land_a)
    call fill_cc(fracice_a,ik,0,land_a)
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
    call doints4(zss_a,    zss)
    call doints4(tss_l_a,  tss_l)
    call doints4(tss_s_a,  tss_s)
    call doints4(fracice_a,fracice)
    call doints4(sicedep_a,sicedep)
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
if (nspecial==44.or.nspecial==46) then
  do iq=1,ifull
    rlongd=rlongg(iq)*180./pi
    rlatd=rlatt(iq)*180./pi
    if (rlatd>=-43..and.rlatd<=-30.) then
      if (rlongd>=155..and.rlongd<=170.) then
        tss(iq)=tss(iq)+1.
      end if
    end if
  end do
end if
if (nspecial==45.or.nspecial==46) then
  do iq=1,ifull
    rlongd=rlongg(iq)*180./pi
    rlatd=rlatt(iq)*180./pi
    if (rlatd>=-15..and.rlatd<=-5.) then
      if (rlongd>=150..and.rlongd<=170.) then
        tss(iq)=tss(iq)+1.
      end if
    end if
  end do
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
  call vertint(u_k ,t, 1,kk,sigin)
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
      call interpwind(u_k(:,k),v_k(:,k),ucc,vcc,dk)
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
    call gethist4('mixr',qg,2)                                  !     mixing ratio
  else
    call gethist4('q',qg,2)                                     !     mixing ratio
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

  !--------------------------------------------------
  ! check soil variables
  if ( myid==0 .or. pfall ) then
    ierc(7:7+3*ms)=0
    if ( ccycle==0 ) then
      call ccnf_inq_varid(ncid,'cplant1',idv,tst)
      if ( tst ) ierc(7)=-1
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
        
  !--------------------------------------------------

  ! SNOW !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call filhist1('snd',dk,snowd,sea_a)

  ! SOIL TEMPERATURE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
      call fill_cc(ucc,dk,0,sea_a)
      call doints4(ucc,tgg(:,k))
    end if
  end do
  if ( .not.iotest ) then
    where ( snowd>0. )
      tgg(:,1)=min(tgg(:,1),270.1)
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
          call filhist1(vname,dk,micdwn(:,k),land_a)
          if ( all(micdwn(:,k)==0.) ) micdwn(:,k)=280.
        case(5)
          micdwn(:,k)=fracice ! read above with nudging arrays
        case(6)
          micdwn(:,k)=sicedep ! read above with nudging arrays
        case(7)
          micdwn(:,k)=snowd*1.E-3
      end select
    end do
    call filhist1('sto',dk,micdwn(:,8),land_a)
    if ( iotest ) then
      call histrd1(iarchi,ier,'uic',ik,micdwn(:,9),ifull)
      call histrd1(iarchi,ier,'vic',ik,micdwn(:,10),ifull)
    else
      call histrd1(iarchi,ier,'uic',ik,ucc,6*ik*ik)
      call histrd1(iarchi,ier,'vic',ik,vcc,6*ik*ik)
      call fill_cc(ucc,dk,0,land_a)
      call fill_cc(vcc,dk,0,land_a)
      call interpwind(micdwn(:,9),micdwn(:,10),ucc,vcc,dk)
    end if ! iotest
    call filhist1('icesal',dk,micdwn(:,11),land_a)
    if ( abs(nmlo)>=2 ) then
      call gethist1('swater',watbdy)
      call gethist1('ssalin',salbdy)
    end if
  end if
  !--------------------------------------------------

  ! SOIL MOISTURE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
      call fill_cc(ucc,dk,0,sea_a)
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
  call filhist1('wetfac',dk,wetfac,sea_a)
  where ( .not.land )
    wetfac=1.
  end where

  !--------------------------------------------------
  ! Read 10m wind speeds for special sea roughness length calculations
  if ( nested==0 ) then
    if ( ierc(2)==0 ) then
      call gethist1('u10',u10)
    else
      u10=sqrt(u(1:ifull,1)**2+v(1:ifull,1)**2)*log(10./0.001)/log(zmin/0.001)
    end if
  end if

  !--------------------------------------------------
  ! Read boundary layer height for TKE-eps mixing
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

  !--------------------------------------------------
  ! Read CABLE/CASA aggregate C+N+P pools
  if ( nsib>=6 ) then
    if ( ccycle==0 ) then
      if ( ierc(7)==0 ) then
        do k=1,ncp
          write(vname,'("cplant",I1.1)') k
          call filhist1(vname,dk,cplant(:,k),sea_a)
        end do
        do k=1,ncs
          write(vname,'("csoil",I1.1)') k
          call filhist1(vname,dk,csoil(:,k),sea_a)
        end do
      end if
    else
      if ( ierc(7)==0 ) then
        do k=1,mplant
          write(vname,'("cplant",I1.1)') k
          call filhist1(vname,dk,cplant(:,k),sea_a)
          write(vname,'("nplant",I1.1)') k
          call filhist1(vname,dk,niplant(:,k),sea_a)
          write(vname,'("pplant",I1.1)') k
          call filhist1(vname,dk,pplant(:,k),sea_a)
        end do
        do k=1,mlitter
          write(vname,'("clitter",I1.1)') k
          call filhist1(vname,dk,clitter(:,k),sea_a)
          write(vname,'("nlitter",I1.1)') k
          call filhist1(vname,dk,nilitter(:,k),sea_a)
          write(vname,'("plitter",I1.1)') k
          call filhist1(vname,dk,plitter(:,k),sea_a)
        end do         
        do k=1,msoil
          write(vname,'("csoil",I1.1)') k
          call filhist1(vname,dk,csoil(:,k),sea_a)
          write(vname,'("nsoil",I1.1)') k
          call filhist1(vname,dk,nisoil(:,k),sea_a)
          write(vname,'("psoil",I1.1)') k
          call filhist1(vname,dk,psoil(:,k),sea_a)
        end do
        call filhist1('glai',dk,glai,sea_a)
      end if ! ierc(7)==0
    end if ! ccycle==0 ..else..
  end if ! if nsib==6.or.nsib==7

  !--------------------------------------------------
  ! Read urban data
  if ( nurban/=0 ) then
    if ( .not.allocated(atebdwn) ) allocate(atebdwn(ifull,24))
    do k=1,24
      select case(k)
        case(1)
          vname='rooftgg1'
        case(2)
          vname='rooftgg2'
        case(3)
          vname='rooftgg3'
        case(4)
          vname='waletgg1'
        case(5)
          vname='waletgg2'
        case(6)
          vname='waletgg3'
        case(7)
          vname='walwtgg1'
        case(8)
          vname='walwtgg2'
        case(9)
          vname='walwtgg3'
        case(10)
          vname='roadtgg1'
        case(11)
          vname='roadtgg2'
        case(12)
          vname='roadtgg3'
        case(13)
          vname='urbnsmc'
        case(14)
          vname='urbnsmr'
        case(15)
          vname='roofwtr'
        case(16)
          vname='roadwtr'
        case(17)
          vname='urbwtrc'
        case(18)
          vname='urbwtrr'
        case(19)
          vname='roofsnd'
        case(20)
          vname='roadsnd'
        case(21)
          vname='roofden'
        case(22)
          vname='roadden'
        case(23)
          vname='roofsna'
        case(24)
          vname='roadsna'
      end select
      if ( iotest ) then
        call histrd1(iarchi,ier,vname,ik,atebdwn(:,k),ifull)
      else
        call histrd1(iarchi,ier,vname,ik,ucc,6*ik*ik)
        where ( ucc>=399. )
          ucc=999.
        end where
        call fill_cc(ucc,dk,0,sea_a)
        call doints4(ucc,atebdwn(:,k))
      end if ! iotest
      if ( all(atebdwn(:,k)==0.) ) then
        select case(k)
          case(1:12)
            atebdwn(:,k)=300.
          case(21:22)
            atebdwn(:,k)=100.
          case(23:24)
            atebdwn(:,k)=0.85
        end select
      end if
    end do
  end if
  !--------------------------------------------------
        
  if ( nested==0 ) then
    ! OMEGA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! only for restart - no interpolation
    dpsldt=-999.
    if ( lrestart ) then
      do k=1,kk 
       call histrd4s(iarchi,ier,'omega',ik,k,dpsldt(:,k),ifull)
       dpsldt(:,k)=dpsldt(:,k)/(1.e5*exp(psl(1:ifull)))
      enddo  ! k loop
    end if

    ! CLOUD FROZEN WATER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call gethist4('qfg',qfg,5)
    ! CLOUD LIQUID WATER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call gethist4('qlg',qlg,5)
    ! RAIN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call gethist4('qrg',qrg,5)
    ! CLOUD FRACTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call gethist4('cfrac',cfrac,5)
    ! RAIN FRACTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call gethist4('cfrain',cffall,5)
    if ( ncloud>=3 ) then
      ! STRAT CLOUD FRACTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call gethist4('stratcf',stratcloud,5)
      ! STRAT NET TENDENCY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call gethist4('strat_nt',nettend,5)
    end if ! (ncloud>=3)
  end if   ! (nested==0)

  !--------------------------------------------------
  ! TKE-eps data
  if ( nvmix==6 .and. nested==0 ) then
    call gethist4('tke',tke,5)
    if ( all(tke(1:ifull,:)==0.) ) tke(1:ifull,:)=1.5E-4
    call gethist4('eps',eps,5)
    if  (all(eps(1:ifull,:)==0.) ) eps(1:ifull,:)=1.E-7
  end if

  !------------------------------------------------------------
  ! Tracer data
  if ( ngas>0 ) then              
    do igas=1,ngas              
      write(trnum,'(i3.3)') igas
      call gethist4('tr'//trnum,tr(:,:,igas),7)
    enddo                       
  endif                         

  !------------------------------------------------------------
  ! Aerosol data
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
  
  if ( nested==0 ) then
    ! GEOPOTENTIAL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! only for restart - no interpolation
    phi_nh=0.
    if ( lrestart ) then
      do k=1,kk 
        call histrd4s(iarchi,ier,'zgnhs',ik,k,phi_nh(:,k),ifull)
      enddo  ! k loop
    end if

    ! SDOT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! only for restart - no interpolation
    sdot=-999.
    if ( lrestart ) then
      sdot(:,1)=0.
      do k=1,kk 
        call histrd4s(iarchi,ier,'sdot',ik,k,sdot(:,k+1),ifull)
      enddo  ! k loop
    end if

    ! PSLX !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! only for restart - no interpolation
    pslx=-999.
    if ( lrestart ) then
      do k=1,kk 
       call histrd4s(iarchi,ier,'pslx',ik,k,pslx(:,k),ifull)
      enddo  ! k loop
    end if
          
    ! SAVU !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! only for restart - no interpolation
    savu=-999.
    if ( lrestart ) then
      do k=1,kk 
        call histrd4s(iarchi,ier,'savu',ik,k,savu(:,k),ifull)
      enddo  ! k loop
    end if
          
    ! SAVV !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! only for restart - no interpolation
    savv=-999.
    if ( lrestart ) then
      do k=1,kk 
        call histrd4s(iarchi,ier,'savv',ik,k,savv(:,k),ifull)
      enddo  ! k loop
    end if

    ! SAVU1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! only for restart - no interpolation
    savu1=-999.
    if ( lrestart ) then
      do k=1,kk 
        call histrd4s(iarchi,ier,'savu1',ik,k,savu1(:,k),ifull)
      enddo  ! k loop
    end if
          
    ! SAVV1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! only for restart - no interpolation
    savv1=-999.
    if ( lrestart ) then
      do k=1,kk 
        call histrd4s(iarchi,ier,'savv1',ik,k,savv1(:,k),ifull)
      enddo  ! k loop
    end if

    ! SAVU2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! only for restart - no interpolation
    savu2=-999.
    if ( lrestart ) then
      do k=1,kk 
        call histrd4s(iarchi,ier,'savu2',ik,k,savu2(:,k),ifull)
      enddo  ! k loop
    end if
          
    ! SAVV2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! only for restart - no interpolation
    savv2=-999.
    if ( lrestart ) then
      do k=1,kk 
        call histrd4s(iarchi,ier,'savv2',ik,k,savv2(:,k),ifull)
      enddo  ! k loop
    end if
          
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

  ! SOIL ICE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do k=1,ms
    write(vname,'("wbice",I1.1)') k
    call filhist1(vname,dk,wbice(:,k),sea_a)
  end do

  ! SNOW !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if ( nmlo==0 .or. abs(nmlo)>9 ) then ! otherwise already read above
    do k=1,3
      write(vname,'("tggsn",I1.1)') k
      call filhist1(vname,dk,tggsn(:,k),sea_a)
      if ( all(tggsn(:,k)==0.) ) tggsn(:,k)=280.
      where( .not.land )
        tggsn(:,k)=280.
      end where
    end do
  end if
  do k=1,3
    write(vname,'("smass",I1.1)') k
    call filhist1(vname,dk,smass(:,k),sea_a)
  end do
  do k=1,3
    write(vname,'("ssdn",I1.1)') k
    call filhist1(vname,dk,ssdn(:,k),sea_a)
    if ( all(ssdn(:,k)==0.) ) then
      where ( snowd>100. )
        ssdn(:,k)=240.
      elsewhere
        ssdn(:,k)=140.
      end where
    end if
  end do
  ssdnn=ssdn(:,1)
  call filhist1('snage',dk,snage,sea_a)

  call filhist1('sflag',dk,dum6,sea_a)
  isflag=nint(dum6)

  ! sgsave is needed for convection
  call gethist1('sgsave',sgsave)
        
endif    ! (nested/=1)

!**************************************************************
! This is the end of reading the initial arrays
!**************************************************************         
          
! tgg holds file surface temperature when no MLO
if ( nmlo==0 .or. abs(nmlo)>9 ) then
  where ( .not.land )
    tgg(:,1)=tss
  end where
end if

! set-up for next read of file
iarchi=iarchi+1
kdate_s=kdate_r
ktime_s=ktime_r+1

if ( myid==0 .and. nested==0 ) then
  write(6,*) "Final lrestart ",lrestart
end if

return
end subroutine ontheflyx

subroutine doints4(s,sout)  ! does calls to intsb
      
use cc_mpi                 ! CC MPI routines
      
implicit none
     
include 'newmpar.h'        ! Grid parameters
include 'parm.h'           ! Model configuration
      
integer, parameter :: ntest=0
integer :: iq, mm
real, dimension(ik*ik*6), intent(inout) :: s
real, dimension(ifull), intent(inout) :: sout
real, dimension(ifull,4) :: wrk

call ccmpi_bcast(s,0,comm_world)

if ( nord==1 ) then
  do mm=1,m_fly  !  was 4, now may be 1
    call ints_blb(s,wrk(:,mm),nface4(:,mm),xg4(:,mm),yg4(:,mm))
  enddo
else
  do mm=1,m_fly  !  was 4, now may be 1
    call intsb(s,wrk(:,mm),nface4(:,mm),xg4(:,mm),yg4(:,mm))
  enddo
endif   ! (nord==1)  .. else ..

#ifdef debug
if(ntest>0.and.mydiag)then
   write(6,*)'in ints4 for id,jd,nord: ',id,jd,nord
   write(6,*)'nface4(1-4) ',(nface4(idjd,mm),mm=1,4)
   write(6,*)'xg4(1-4) ',(xg4(idjd,mm),mm=1,4)
   write(6,*)'yg4(1-4) ',(yg4(idjd,mm),mm=1,4)
   write(6,*)'wrk(1-4) ',(wrk(idjd,mm),mm=1,4)
endif
#endif

if ( m_fly==1 ) then
  sout=wrk(:,1)
else
! average 4 m values to central value
  sout(:)=.25*sum(wrk,2)
endif    ! (m_fly==1)

return
end subroutine doints4

subroutine intsb(s,sout,nface_l,xg_l,yg_l)   ! N.B. sout here
      
!     same as subr ints, but with sout passed back and no B-S      
!     s is input; sout is output array
!     later may wish to save idel etc between array calls
!     this one does linear interp in x on outer y sides
!     doing x-interpolation before y-interpolation
!     This is a global routine 

use cc_mpi, only : indx

implicit none
      
include 'newmpar.h'  ! Grid parameters
include 'parm.h'     ! Model configuration

integer, dimension(ifull), intent(in) :: nface_l
integer :: idel, jdel, nn
integer :: i, j, n, iq, n_n, n_e, n_w, n_s
real, dimension(ik*ik*6), intent(in) :: s
real, dimension(ifull), intent(inout) :: sout
real aaa, c1, c2, c3, c4, xxg, yyg
real, intent(in), dimension(ifull) :: xg_l, yg_l
real, dimension(-1:ik+2,-1:ik+2,0:npanels) :: sx
real, dimension(4) :: r

!     this is intsb           EW interps done first
!     first extend s arrays into sx - this one -1:il+2 & -1:il+2
sx(1:ik,1:ik,0:npanels) = reshape(s(1:ik*ik*6), (/ik,ik,6/))
do n=0,npanels,2
  n_w=mod(n+5,6)
  n_e=mod(n+2,6)
  n_n=mod(n+1,6)
  n_s=mod(n+4,6)
  do j=1,ik
    sx(0,j,n)   =s(indx(ik,j,n_w,ik,ik))
    sx(-1,j,n)  =s(indx(ik-1,j,n_w,ik,ik))
    sx(ik+1,j,n)=s(indx(ik+1-j,1,n_e,ik,ik))
    sx(ik+2,j,n)=s(indx(ik+1-j,2,n_e,ik,ik))
  enddo
  do i=1,ik
    sx(i,ik+1,n)=s(indx(i,1,n+1,ik,ik))
    sx(i,ik+2,n)=s(indx(i,2,n+1,ik,ik))
    sx(i,0,n)   =s(indx(ik,ik+1-i,n_s,ik,ik))
    sx(i,-1,n)  =s(indx(ik-1,ik+1-i,n_s,ik,ik))
  enddo
  sx(-1,0,n)     =s(indx(ik,2,n_w,ik,ik))    ! wws
  sx(0,-1,n)     =s(indx(ik,ik-1,n_s,ik,ik)) ! wss
  sx(0,0,n)      =s(indx(ik,1,n_w,ik,ik))    ! ws
  sx(ik+1,0,n)   =s(indx(ik,1,n_e,ik,ik))    ! es  
  sx(ik+2,0,n)   =s(indx(ik-1,1,n_e,ik,ik))  ! ees 
  sx(-1,ik+1,n)  =s(indx(ik,ik-1,n_w,ik,ik)) ! wwn
  sx(0,ik+2,n)   =s(indx(ik-1,ik,n_w,ik,ik)) ! wnn
  sx(ik+2,ik+1,n)=s(indx(2,1,n_e,ik,ik))     ! een  
  sx(ik+1,ik+2,n)=s(indx(1,2,n_e,ik,ik))     ! enn  
  sx(0,ik+1,n)   =s(indx(ik,ik,n_w,ik,ik))   ! wn  
  sx(ik+1,ik+1,n)=s(indx(1,1,n_e,ik,ik))     ! en  
  sx(ik+1,-1,n)  =s(indx(ik,2,n_e,ik,ik))    ! ess  
enddo  ! n loop
do n=1,npanels,2
  n_w=mod(n+4,6)
  n_e=mod(n+1,6)
  n_n=mod(n+2,6)
  n_s=mod(n+5,6)
  do j=1,ik
    sx(0,j,n)   =s(indx(ik+1-j,ik,n_w,ik,ik))
    sx(-1,j,n)  =s(indx(ik+1-j,ik-1,n_w,ik,ik))
    sx(ik+1,j,n)=s(indx(1,j,n_e,ik,ik))
    sx(ik+2,j,n)=s(indx(2,j,n_e,ik,ik))
  enddo
  do i=1,ik
    sx(i,ik+1,n)=s(indx(1,ik+1-i,n_n,ik,ik))
    sx(i,ik+2,n)=s(indx(2,ik+1-i,n_n,ik,ik))
    sx(i,0,n)   =s(indx(i,ik,n-1,ik,ik))
    sx(i,-1,n)  =s(indx(i,ik-1,n-1,ik,ik))
  enddo
  sx(-1,0,n)     =s(indx(ik-1,ik,n_w,ik,ik)) ! wws
  sx(0,-1,n)     =s(indx(2,ik,n_s,ik,ik))    ! wss
  sx(0,0,n)      =s(indx(ik,ik,n_w,ik,ik))   ! ws
  sx(ik+1,0,n)   =s(indx(1,1,n_e,ik,ik))     ! es
  sx(ik+2,0,n)   =s(indx(1,2,n_e,ik,ik))     ! ees
  sx(-1,ik+1,n)  =s(indx(2,ik,n_w,ik,ik))    ! wwn   
  sx(0,ik+2,n)   =s(indx(1,ik-1,n_w,ik,ik))  ! wnn  
  sx(ik+2,ik+1,n)=s(indx(1,ik-1,n_e,ik,ik))  ! een  
  sx(ik+1,ik+2,n)=s(indx(2,ik,n_e,ik,ik))    ! enn  
  sx(0,ik+1,n)   =s(indx(1,ik,n_w,ik,ik))    ! wn  
  sx(ik+1,ik+1,n)=s(indx(1,ik,n_e,ik,ik))    ! en  
  sx(ik+1,-1,n)  =s(indx(2,1,n_e,ik,ik))     ! ess  
enddo  ! n loop
!     for ew interpolation, sometimes need (different from ns):
!          (-1,0),   (0,0),   (0,-1)   (-1,il+1),   (0,il+1),   (0,il+2)
!         (il+1,0),(il+2,0),(il+1,-1) (il+1,il+1),(il+2,il+1),(il+1,il+2)

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

subroutine ints_blb(s,sout,nface_l,xg_l,yg_l) 
      
!     this one does bi-linear interpolation only

use cc_mpi, only : indx

implicit none
      
include 'newmpar.h'  ! Grid parameters
include 'parm.h'     ! Model configuration

integer :: i, j, n, iq, idel, jdel
integer :: n_n, n_e, n_w, n_s
integer, intent(in), dimension(ifull) :: nface_l
real, dimension(ik*ik*6), intent(inout) :: s
real, dimension(ifull), intent(inout) :: sout
real, intent(in), dimension(ifull) :: xg_l, yg_l
real, dimension(-1:ik+2,-1:ik+2,0:npanels) :: sx
real :: xxg, yyg

!     first extend s arrays into sx - this one -1:il+2 & -1:il+2
!                    but for bi-linear only need 0:il+1 &  0:il+1
sx(1:ik,1:ik,0:npanels) = reshape(s(1:ik*ik*6), (/ik,ik,6/))
do n=0,npanels,2
  n_w=mod(n+5,6)
  n_e=mod(n+2,6)
  n_n=mod(n+1,6)
  n_s=mod(n+4,6)
  do j=1,ik
    sx(0,j,n)   =s(indx(ik,j,n_w,ik,ik))
    sx(ik+1,j,n)=s(indx(ik+1-j,1,n_e,ik,ik))
  enddo
  do i=1,ik
    sx(i,ik+1,n)=s(indx(i,1,n+1,ik,ik))
    sx(i,0,n)   =s(indx(ik,ik+1-i,n_s,ik,ik))
  enddo
  sx(0,0,n)      =s(indx(ik,1,n_w,ik,ik))  ! ws
  sx(ik+1,0,n)   =s(indx(ik,1,n_e,ik,ik))  ! es  
  sx(0,ik+1,n)   =s(indx(ik,ik,n_w,ik,ik)) ! wn  
  sx(ik+1,ik+1,n)=s(indx(1,1,n_e,ik,ik))   ! en  
enddo  ! n loop
do n=1,npanels,2
  n_w=mod(n+4,6)
  n_e=mod(n+1,6)
  n_n=mod(n+2,6)
  n_s=mod(n+5,6)
  do j=1,ik
    sx(0,j,n)   =s(indx(ik+1-j,ik,n_w,ik,ik))
    sx(ik+1,j,n)=s(indx(1,j,n_e,ik,ik))
  enddo
  do i=1,ik
    sx(i,ik+1,n)=s(indx(1,ik+1-i,n_n,ik,ik))
    sx(i,0,n)   =s(indx(i,ik,n-1,ik,ik))
  enddo
  sx(0,0,n)      =s(indx(ik,ik,n_w,ik,ik))  ! ws
  sx(ik+1,0,n)   =s(indx(1,1,n_e,ik,ik))    ! es
  sx(0,ik+1,n)   =s(indx(1,ik,n_w,ik,ik))   ! wn  
  sx(ik+1,ik+1,n)=s(indx(1,ik,n_e,ik,ik))   ! en  
enddo  ! n loop
      

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

subroutine fill_cc(a_io,dk,ndiag,land_a)
      
!     this version holds whole array in memory      
!     routine fills in interior of an array which has undefined points

use cc_mpi          ! CC MPI routines
      
implicit none
      
include 'newmpar.h' ! Grid parameters

integer :: nrem, i, ii, dk, iq, j, n, neighb, ndiag
integer :: iminb, imaxb, jminb, jmaxb
integer, save :: olddk = 0
integer, dimension(:,:), allocatable, save :: ic
integer, dimension(0:5) :: imin,imax,jmin,jmax
integer, dimension(0:5) :: npann,npane,npanw,npans
real, parameter :: value=999.       ! missing value flag
real, dimension(6*dk*dk) :: a_io, a
real av     
logical, dimension(6*dk*dk) :: land_a
logical, dimension(4) :: mask

data npann/1,103,3,105,5,101/,npane/102,2,104,4,100,0/
data npanw/5,105,1,101,3,103/,npans/104,0,100,2,102,4/

if ( myid/=0 ) return

where ( land_a )
  a_io=value
end where
if ( all(abs(a_io-value)<1.E-6) ) return

call START_LOG(otf_fill_begin)
      
if ( dk/=olddk ) then
  olddk=dk
  if ( allocated(ic) ) deallocate(ic)
  allocate(ic(4,dk*dk*6))
  do iq=1,dk*dk*6
    ic(1,iq)=iq+dk
    ic(2,iq)=iq-dk
    ic(3,iq)=iq+1
    ic(4,iq)=iq-1
  enddo   ! iq loop
  do n=0,npanels
    if ( npann(n)<100 ) then
      do ii=1,dk
        ic(1,indx(ii,dk,n,dk,dk))=indx(ii,1,npann(n),dk,dk)
      enddo    ! ii loop
    else
      do ii=1,dk
        ic(1,indx(ii,dk,n,dk,dk))=indx(1,dk+1-ii,npann(n)-100,dk,dk)
      enddo    ! ii loop
    endif      ! (npann(n)<100)
    if ( npane(n)<100 ) then
      do ii=1,dk
        ic(3,indx(dk,ii,n,dk,dk))=indx(1,ii,npane(n),dk,dk)
      enddo    ! ii loop
    else
      do ii=1,dk
        ic(3,indx(dk,ii,n,dk,dk))=indx(dk+1-ii,1,npane(n)-100,dk,dk)
      enddo    ! ii loop
    endif      ! (npane(n)<100)
    if ( npanw(n)<100 ) then
      do ii=1,dk
        ic(4,indx(1,ii,n,dk,dk))=indx(dk,ii,npanw(n),dk,dk)
      enddo    ! ii loop
    else
      do ii=1,dk
        ic(4,indx(1,ii,n,dk,dk))=indx(dk+1-ii,dk,npanw(n)-100,dk,dk)
      enddo    ! ii loop
    endif      ! (npanw(n)<100)
    if ( npans(n)<100 ) then
      do ii=1,dk
        ic(2,indx(ii,1,n,dk,dk))=indx(ii,dk,npans(n),dk,dk)
      enddo    ! ii loop
    else
      do ii=1,dk
        ic(2,indx(ii,1,n,dk,dk))=indx(dk,dk+1-ii,npans(n)-100,dk,dk)
      enddo    ! ii loop
    endif      ! (npans(n)<100)
  enddo      ! n loop
end if ! olddk/=dk

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
          mask=a(ic(:,iq))/=value
          neighb=count(mask)
          if ( neighb>0 ) then
            av=sum(a(ic(:,iq)),mask)
            a_io(iq)=av/real(neighb)
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

subroutine mslpx(pmsl,psl,zs,t,siglev)
      
!     generalized from ifull to allow new onthefly    0808
!     can replace usual mslp sometime

use cc_mpi, only : mydiag ! CC MPI routines
use sigs_m                ! Atmosphere sigma levels
      
!     this one will ignore negative zs (i.e. over the ocean)
implicit none
      
include 'newmpar.h'       ! Grid parameters
include 'const_phys.h'    ! Physical constants
include 'parm.h'          ! Model configuration

integer iq
real siglev
real, dimension(6*ik*ik) :: pmsl,psl,zs,t
real c, con, conr
real, dimension(6*ik*ik) :: dlnps, phi1, tav, tsurf

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
!     generalized from ifull to allow new onthefly    0808
!     can replace usual mslp sometime

use cc_mpi, only : mydiag  ! CC MPI routines
use sigs_m                 ! Atmosphere sigma levels
      
!     this one will ignore negative zs (i.e. over the ocean)
implicit none
      
include 'newmpar.h'        ! Grid parameters
include 'const_phys.h'     ! Physical constants
include 'parm.h'           ! Model configuration
      
integer iq,lev
real, dimension(ifull) :: pmsl,psl,zs,t
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
      
subroutine interpwind(uct,vct,ucc,vcc,dk)
      
use cc_mpi           ! CC MPI routines
use vecsuv_m         ! Map to cartesian coordinates
      
implicit none
      
include 'newmpar.h'  ! Grid parameters
      
integer, intent(in) :: dk
integer iq
real, dimension(6*dk*dk) :: uc,vc,wc
real, dimension(6*ik*ik), intent(inout) :: ucc,vcc
real, dimension(6*ik*ik) :: wcc
real, dimension(ifull), intent(out) :: uct,vct
real, dimension(ifull) :: wct
      
do iq=1,dk*dk*6
  ! first set up winds in Cartesian "source" coords            
  uc(iq)=axs_a(iq)*ucc(iq) + bxs_a(iq)*vcc(iq)
  vc(iq)=ays_a(iq)*ucc(iq) + bys_a(iq)*vcc(iq)
  wc(iq)=azs_a(iq)*ucc(iq) + bzs_a(iq)*vcc(iq)
  ! now convert to winds in "absolute" Cartesian components
  ucc(iq)=uc(iq)*rotpoles(1,1)+vc(iq)*rotpoles(1,2)+wc(iq)*rotpoles(1,3)
  vcc(iq)=uc(iq)*rotpoles(2,1)+vc(iq)*rotpoles(2,2)+wc(iq)*rotpoles(2,3)
  wcc(iq)=uc(iq)*rotpoles(3,1)+vc(iq)*rotpoles(3,2)+wc(iq)*rotpoles(3,3)
end do
! interpolate all required arrays to new C-C positions
! don't need to do map factors and Coriolis on target grid
call doints4(ucc,  uct)
call doints4(vcc,  vct)
call doints4(wcc,  wct)
do iq=1,ifull
  ! now convert to "target" Cartesian components (transpose used)
  ucc(iq)=uct(iq)*rotpole(1,1)+vct(iq)*rotpole(2,1)+wct(iq)*rotpole(3,1)
  vcc(iq)=uct(iq)*rotpole(1,2)+vct(iq)*rotpole(2,2)+wct(iq)*rotpole(3,2)
  wcc(iq)=uct(iq)*rotpole(1,3)+vct(iq)*rotpole(2,3)+wct(iq)*rotpole(3,3)
  ! then finally to "target" local x-y components
  uct(iq) = ax(iq)*ucc(iq) + ay(iq)*vcc(iq) +  az(iq)*wcc(iq)
  vct(iq) = bx(iq)*ucc(iq) + by(iq)*vcc(iq) +  bz(iq)*wcc(iq)
enddo               ! iq loop
      
return
end subroutine interpwind

subroutine gethist1(vname,varout)
      
use infile             ! Input file routines
      
implicit none
      
include 'newmpar.h'    ! Grid parameters
include 'darcdf.h'     ! Netcdf data
      
integer ier
real, dimension(:), intent(out) :: varout
real, dimension(6*ik*ik) :: ucc
character(len=*), intent(in) :: vname
      
if ( iotest ) then
  call histrd1(iarchi,ier,vname,ik,varout,ifull)
else
  call histrd1(iarchi,ier,vname,ik,ucc,6*ik*ik)
  call doints4(ucc,varout)
end if ! iotest

return
end subroutine gethist1
     
subroutine filhist1(vname,dk,varout,mask_a)
      
use infile             ! Input file routines
      
implicit none
      
include 'newmpar.h'    ! Grid parameters
include 'darcdf.h'     ! Netcdf data
      
integer, intent(in) :: dk
integer ier
real, dimension(:), intent(out) :: varout
real, dimension(6*ik*ik) :: ucc
logical, dimension(6*dk*dk), intent(in) :: mask_a
character(len=*), intent(in) :: vname
      
if ( iotest ) then
  call histrd1(iarchi,ier,vname,ik,varout,ifull)
else
  call histrd1(iarchi,ier,vname,ik,ucc,6*ik*ik)
  call fill_cc(ucc,dk,0,mask_a)
  call doints4(ucc,varout)
end if ! iotest
      
return
end subroutine filhist1
     
subroutine gethist4(vname,varout,vmode)
      
use infile             ! Input file routines
      
implicit none
      
include 'newmpar.h'    ! Grid parameters
include 'darcdf.h'     ! Netcdf data
      
integer, intent(in) :: vmode
integer k, ier
real, dimension(:,:), intent(out) :: varout
real, dimension(6*ik*ik) :: ucc
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
