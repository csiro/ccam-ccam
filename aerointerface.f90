
! This module interfaces aerosol schemes with CCAM

! Currently, only LDR aerosols are supported.  However, the xtg arrays should be moved to
! here in the future, so that there is a common interface for advection and outcdf.

module aerointerface

implicit none

private
public load_aerosolldr,aerocalc,aerodrop
public ppfprec,ppfmelt,ppfsnow,ppfconv,ppfevap,ppfsubl,pplambs,ppmrate, &
       ppmaccr,ppfstay,ppqfsed,pprscav

real, dimension(:,:,:,:), allocatable, save :: oxidantprev
real, dimension(:,:,:,:), allocatable, save :: oxidantnow
real, dimension(:,:,:,:), allocatable, save :: oxidantnext
real, dimension(:), allocatable, save :: rlat,rlev
real, dimension(:,:), allocatable, save :: ppfprec,ppfmelt,ppfsnow,ppfconv  ! data saved from LDR cloud scheme
real, dimension(:,:), allocatable, save :: ppfevap,ppfsubl,pplambs,ppmrate  ! data saved from LDR cloud scheme
real, dimension(:,:), allocatable, save :: ppmaccr,ppfstay,ppqfsed,pprscav  ! data saved from LDR cloud scheme
integer, save :: ilon,ilat,ilev

real, parameter :: wlc      = 0.2e-3         ! LWC of deep conv cloud (kg/m**3)

contains

! Load aerosols emissions from netcdf
subroutine load_aerosolldr(aerofile,oxidantfile)
      
use aerosolldr
use cc_mpi
      
implicit none

include 'dates.h'      
include 'netcdf.inc'
include 'newmpar.h'
include 'mpif.h'
include 'parmgeom.h'
      
integer ncstatus,ncid,i,j,varid,tilg
integer ierr,jyear,jmonth
integer, dimension(2) :: spos,npos
integer, dimension(4) :: sposs,nposs
real, dimension(:), allocatable :: dumg
real, dimension(ifull) :: duma
real tlat,tlon,tschmidt
character(len=*), intent(in) :: aerofile,oxidantfile

allocate(ppfprec(ifull,kl),ppfmelt(ifull,kl))
allocate(ppfsnow(ifull,kl),ppfconv(ifull,kl))
allocate(ppfevap(ifull,kl),ppfsubl(ifull,kl))
allocate(pplambs(ifull,kl),ppmrate(ifull,kl))
allocate(ppmaccr(ifull,kl),ppfstay(ifull,kl))
allocate(ppqfsed(ifull,kl),pprscav(ifull,kl))
ppfprec=0.
ppfmelt=0.
ppfsnow=0.
ppfconv=0.
ppfevap=0.
ppfsubl=0.
pplambs=0.
ppmrate=0.
ppmaccr=0.
ppfstay=0.
ppqfsed=0.
pprscav=0.

call aldrinit(ifull,iextra,kl)

if (myid==0) then
  allocate(dumg(ifull_g))
  ncstatus=nf_open(aerofile,nf_nowrite,ncid)
  if (ncstatus.ne.nf_noerr) then
    write(6,*) "ERROR: Cannot open ",trim(aerofile)
    stop
  end if
  write(6,*) "Reading ",trim(aerofile)
  ! check dimensions and location
  ncstatus=nf_get_att_real(ncid,nf_global,'lat0',tlat)
  ncstatus=nf_get_att_real(ncid,nf_global,'lon0',tlon)
  ncstatus=nf_get_att_real(ncid,nf_global,'schmidt0',tschmidt)
  if (rlong0.ne.tlon.or.rlat0.ne.tlat.or.schmidt.ne.tschmidt) then
    write(6,*) "ERROR: Grid mismatch for ",trim(aerofile)
    write(6,*) "rlong0,rlat0,schmidt ",rlong0,rlat0,schmidt
    write(6,*) "tlon,tlat,tschmidt   ",tlon,tlat,tschmidt
    stop
  end if
  ncstatus = nf_inq_dimid(ncid,'longitude',varid)
  ncstatus = nf_inq_dimlen(ncid,varid,tilg)
  if (tilg.ne.il_g) then
    write (6,*) "ERROR: Grid mismatch for ",trim(aerofile)
    write (6,*) "il_g,tilg ",il_g,tilg
    stop
  end if
  ! load emission fields
  spos=1
  npos(1)=il_g
  npos(2)=il_g*6
  ncstatus = nf_inq_varid(ncid,'so2a1',varid)
  ncstatus = nf_get_vara_real(ncid,varid,spos,npos,dumg)
  call ccmpi_distribute(duma,dumg)
  call aldrloademiss(1,duma)
  ncstatus = nf_inq_varid(ncid,'so2a2',varid)
  ncstatus = nf_get_vara_real(ncid,varid,spos,npos,dumg)
  call ccmpi_distribute(duma,dumg)
  call aldrloademiss(2,duma)
  ncstatus = nf_inq_varid(ncid,'bca1',varid)
  ncstatus = nf_get_vara_real(ncid,varid,spos,npos,dumg)
  call ccmpi_distribute(duma,dumg)
  call aldrloademiss(3,duma)
  ncstatus = nf_inq_varid(ncid,'bca2',varid)
  ncstatus = nf_get_vara_real(ncid,varid,spos,npos,dumg)
  call ccmpi_distribute(duma,dumg)
  call aldrloademiss(4,duma)
  ncstatus = nf_inq_varid(ncid,'oca1',varid)
  ncstatus = nf_get_vara_real(ncid,varid,spos,npos,dumg)
  call ccmpi_distribute(duma,dumg)
  call aldrloademiss(5,duma)
  ncstatus = nf_inq_varid(ncid,'oca2',varid)
  ncstatus = nf_get_vara_real(ncid,varid,spos,npos,dumg)
  call ccmpi_distribute(duma,dumg)
  call aldrloademiss(6,duma)
  ncstatus = nf_inq_varid(ncid,'so2b1',varid)
  ncstatus = nf_get_vara_real(ncid,varid,spos,npos,dumg)
  call ccmpi_distribute(duma,dumg)
  call aldrloademiss(7,duma)
  ncstatus = nf_inq_varid(ncid,'so2b2',varid)
  ncstatus = nf_get_vara_real(ncid,varid,spos,npos,dumg)
  call ccmpi_distribute(duma,dumg)
  call aldrloademiss(8,duma)
  ncstatus = nf_inq_varid(ncid,'bcb1',varid)
  ncstatus = nf_get_vara_real(ncid,varid,spos,npos,dumg)
  call ccmpi_distribute(duma,dumg)
  call aldrloademiss(9,duma)
  ncstatus = nf_inq_varid(ncid,'bcb2',varid)
  ncstatus = nf_get_vara_real(ncid,varid,spos,npos,dumg)
  call ccmpi_distribute(duma,dumg)
  call aldrloademiss(10,duma)
  ncstatus = nf_inq_varid(ncid,'ocb1',varid)
  ncstatus = nf_get_vara_real(ncid,varid,spos,npos,dumg)
  call ccmpi_distribute(duma,dumg)
  call aldrloademiss(11,duma)
  ncstatus = nf_inq_varid(ncid,'ocb2',varid)
  ncstatus = nf_get_vara_real(ncid,varid,spos,npos,dumg)
  call ccmpi_distribute(duma,dumg)
  call aldrloademiss(12,duma)
  ncstatus = nf_inq_varid(ncid,'dmso',varid)
  ncstatus = nf_get_vara_real(ncid,varid,spos,npos,dumg)
  call ccmpi_distribute(duma,dumg)
  call aldrloademiss(13,duma)
  ncstatus = nf_inq_varid(ncid,'dmst',varid)
  ncstatus = nf_get_vara_real(ncid,varid,spos,npos,dumg)
  call ccmpi_distribute(duma,dumg)
  call aldrloademiss(14,duma)
  ncstatus = nf_inq_varid(ncid,'ocna',varid)
  ncstatus = nf_get_vara_real(ncid,varid,spos,npos,dumg)
  call ccmpi_distribute(duma,dumg)
  call aldrloademiss(15,duma)
  ncstatus = nf_inq_varid(ncid,'vso2',varid)
  ncstatus = nf_get_vara_real(ncid,varid,spos,npos,dumg)
  call ccmpi_distribute(duma,dumg)
  call aldrloademiss(16,duma)
  ! load dust fields
  ncstatus = nf_inq_varid(ncid,'sandem',varid)
  ncstatus = nf_get_vara_real(ncid,varid,spos,npos,dumg)
  call ccmpi_distribute(duma,dumg)
  call aldrloaderod(1,1,duma)
  ncstatus = nf_inq_varid(ncid,'siltem',varid)
  ncstatus = nf_get_vara_real(ncid,varid,spos,npos,dumg)
  call ccmpi_distribute(duma,dumg)
  call aldrloaderod(2,1,duma)
  ncstatus = nf_inq_varid(ncid,'clayem',varid)
  ncstatus = nf_get_vara_real(ncid,varid,spos,npos,dumg)
  call ccmpi_distribute(duma,dumg)
  call aldrloaderod(3,1,duma)
  ncstatus=nf_close(ncid)
  deallocate(dumg)
  ! load oxidant fields
  ncstatus=nf_open(oxidantfile,nf_nowrite,ncid)
  if (ncstatus.ne.nf_noerr) then
    write(6,*) "ERROR: Cannot open ",trim(oxidantfile)
    stop
  end if
  write(6,*) "Reading ",trim(oxidantfile)
  ! check dimensions and location
  ncstatus = nf_inq_dimid(ncid,'lon',varid)
  ncstatus = nf_inq_dimlen(ncid,varid,ilon)
  ncstatus = nf_inq_dimid(ncid,'lat',varid)
  ncstatus = nf_inq_dimlen(ncid,varid,ilat)
  ncstatus = nf_inq_dimid(ncid,'lev',varid)
  ncstatus = nf_inq_dimlen(ncid,varid,ilev)
  call MPI_Bcast(ilon,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(ilat,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(ilev,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  allocate(oxidantprev(ilon,ilat,ilev,4))
  allocate(oxidantnow(ilon,ilat,ilev,4))
  allocate(oxidantnext(ilon,ilat,ilev,4))
  allocate(rlat(ilat),rlev(ilev))
  sposs=1
  nposs(1)=ilon
  nposs(2)=ilat
  nposs(3)=ilev
  nposs(4)=1
  jyear=kdate/10000
  jmonth=(kdate-jyear*10000)/100
  ncstatus = nf_inq_varid(ncid,'lat',varid)
  ncstatus = nf_get_vara_real(ncid,varid,sposs(2),nposs(2),rlat)
  ncstatus = nf_inq_varid(ncid,'lev',varid)
  ncstatus = nf_get_vara_real(ncid,varid,sposs(3),nposs(3),rlev)
  call MPI_Bcast(rlat,ilat,MPI_REAL,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(rlev,ilev,MPI_REAL,0,MPI_COMM_WORLD,ierr)
  do i=1,3
    select case(i)
      case(1)
        sposs(4)=jmonth-1
      case(2)
        sposs(4)=jmonth
      case(3)
        sposs(4)=jmonth+1
    end select
    if (sposs(4).lt.1) sposs(4)=12
    if (sposs(4).gt.12) sposs(4)=1
    do j=1,4
      select case(j)
        case(1)
          ncstatus = nf_inq_varid(ncid,'OH',varid)
        case(2)
          ncstatus = nf_inq_varid(ncid,'H2O2',varid)
        case(3)
          ncstatus = nf_inq_varid(ncid,'O3',varid)
        case(4)
          ncstatus = nf_inq_varid(ncid,'NO2',varid)
      end select
      select case(i)
        case(1)
          ncstatus = nf_get_vara_real(ncid,varid,sposs,nposs,oxidantprev(:,:,:,j))
          call MPI_Bcast(oxidantprev(:,:,:,j),ilon*ilat*ilev,MPI_REAL,0,MPI_COMM_WORLD,ierr)
        case(2)
          ncstatus = nf_get_vara_real(ncid,varid,sposs,nposs,oxidantnow(:,:,:,j))
          call MPI_Bcast(oxidantnow(:,:,:,j),ilon*ilat*ilev,MPI_REAL,0,MPI_COMM_WORLD,ierr)
        case(3)
          ncstatus = nf_get_vara_real(ncid,varid,sposs,nposs,oxidantnext(:,:,:,j))
          call MPI_Bcast(oxidantnext(:,:,:,j),ilon*ilat*ilev,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      end select
    end do
  end do
  ncstatus=nf_close(ncid)
else
  ! load emission fields
  do i=1,16
    call ccmpi_distribute(duma)
    call aldrloademiss(i,duma)
  end do
  ! load dust fields
  do i=1,3
    call ccmpi_distribute(duma)
    call aldrloaderod(i,1,duma)
  end do
  ! load oxidant fields
  call MPI_Bcast(ilon,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(ilat,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(ilev,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  allocate(oxidantprev(ilon,ilat,ilev,4))
  allocate(oxidantnow(ilon,ilat,ilev,4))
  allocate(oxidantnext(ilon,ilat,ilev,4))
  allocate(rlat(ilat),rlev(ilev))
  call MPI_Bcast(rlat,ilat,MPI_REAL,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(rlev,ilev,MPI_REAL,0,MPI_COMM_WORLD,ierr)
  do i=1,3
    do j=1,4
      select case(i)
        case(1)
          call MPI_Bcast(oxidantprev(:,:,:,j),ilon*ilat*ilev,MPI_REAL,0,MPI_COMM_WORLD,ierr)
        case(2)
          call MPI_Bcast(oxidantnow(:,:,:,j),ilon*ilat*ilev,MPI_REAL,0,MPI_COMM_WORLD,ierr)
        case(3)
          call MPI_Bcast(oxidantnext(:,:,:,j),ilon*ilat*ilev,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      end select
    end do
  end do
end if

return
end subroutine load_aerosolldr

subroutine aerocalc(odcalc)

use aerosolldr
use arrays_m
use cfrac_m
use extraout_m
use kuocomb_m
use latlong_m
use liqwpar_m
use map_m
use morepbl_m
use nsibd_m
use ozoneread
use pbl_m
use screen_m
use sigs_m
use soil_m
use soilsnow_m
use vegpar_m
use work2_m

implicit none

include 'newmpar.h'
include 'const_phys.h'
include 'dates.h'
include 'kuocom.h'
include 'parm.h'
include 'soilv.h'

integer jyear,jmonth,jday,jhour,jmin,mstart,mins
integer j,k
integer, dimension(12) :: ndoy
real, dimension(ilon) :: rlon
real, dimension(ifull,kl) :: oxout,zg,clcon,pccw,rhoa
real, dimension(ifull) :: blon,blat,dxy,mcmax,cldcon,wg
real, dimension(kl+1) :: sigh
logical, intent(in) :: odcalc
data ndoy/0,31,59,90,120,151,181,212,243,273,304,334/

! update oxidant fields
jyear=kdate/10000
jmonth=(kdate-jyear*10000)/100
jday=kdate-jyear*10000-jmonth*100
jhour=ktime/100
jmin=ktime-jhour*100
mstart=1440*(ndoy(jmonth)+jday-1) + 60*jhour + jmin
mins = mtimer + mstart

! update oxidant fields when radiation is updated (typically once per hour)
if (odcalc) then
  do j=1,ilon
    rlon(j)=real(j-1)*360./real(ilon)
  end do
  blon=rlongg*180./pi
  blat=rlatt*180./pi

  do j=1,4      
    call fieldinterpolate(oxout(:,:),blon,blat,oxidantprev(:,:,:,j),oxidantnow(:,:,:,j), &
                          oxidantnext(:,:,:,j),rlon,rlat,rlev,ifull,kl,ilon,ilat,ilev,mins,sig,ps)
    do k=1,kl
      call aldrloadoxidant(k+(j-1)*kl,oxout(:,k))
    end do
  end do
end if

! set-up half levels ------------------------------------------------
sigh(1:kl) = sigmh(1:kl)
sigh(kl+1) = 0.
zg(:,1)=bet(1)*t(1:ifull,1)/grav
do k=2,kl
  zg(:,k)=zg(:,k-1)+(bet(k)*t(1:ifull,k)+betm(k)*t(1:ifull,k-1))/grav
end do
dxy=ds*ds/(em*em)
mcmax=0.1*max(vlai,0.1)
do k=1,kl
  rhoa(:,k)=ps*sig(k)/(rdry*t(:,k)) !density of air
end do
where (land)
  wg=(wb(:,1)-swilt(isoilm))/(sfc(isoilm)-swilt(isoilm))
elsewhere
  wg=1.
end where

! estimate convective cloud fraction
! from leoncld.f
where (ktsav<kl-1)
  cldcon=min(acon+bcon*log(1.+condc*86400./dt),0.8) !NCAR
elsewhere
  cldcon=0.
end where
clcon=0.
pccw=0.
if (nmr.ge.1) then
  do k=1,kl
    where (k.le.ktsav.and.k.ge.kbsav+1)
      clcon(:,k)=cldcon ! maximum overlap
      pccw(:,kl+1-k)=wlc/rhoa(:,k)
    end where
  end do
else
  do k=1,kl
    where (k.le.ktsav.and.k.ge.kbsav+1)
      clcon(:,k)=1.-(1.-cldcon)**(1./real(ktsav-kbsav+2)) !Random overlap
      pccw(:,kl+1-k)=wlc/rhoa(:,k)      
    end where
  end do
end if

! update prognostic aerosols
call aldrcalc(dt,sig,sigh,dsig,zg,cansto,mcmax,wg,pblh,ps,  &
              tss,t(1:ifull,:),condx,condc,snowd,sgsave,fg, &
              eg,u10,ustar,zo,land,fracice,sigmf,           &
              qlg(1:ifull,:),qfg(1:ifull,:),cfrac,clcon,    &
              pccw,dxy,rhoa,ppfprec,ppfmelt,ppfsnow,        &
              ppfconv,ppfevap,ppfsubl,pplambs,ppmrate,      &
              ppmaccr,ppfstay,ppqfsed,pprscav)

! Factor 1.e3 to convert to g/m2, x 3 to get sulfate from sulfur
so4t(:)=0.
do k=1,kl
  so4t(:)=so4t(:)+3.e3*xtg(:,k,3)*(-ps(:)*dsig(k))/grav
enddo

return
end subroutine aerocalc

subroutine aerodrop(iaero,istart,imax,kl,cdn,rhoa,land,rlatt)

use aerosolldr

implicit none

include 'const_phys.h'
include 'cparams.h'

integer, intent(in) :: iaero,istart,imax,kl
integer k
real, dimension(imax,kl), intent(out) :: cdn
real, dimension(imax,kl), intent(in) :: rhoa
real, dimension(imax), intent(in) :: rlatt
logical, dimension(imax), intent(in) :: land

select case(abs(iaero))
  case(2)
    call cldrop(istart,imax,cdn,rhoa)
  case default
    where (land(:).and.rlatt(:)>0.)
      cdn(:,1)=cdropl_nh
    else where (land(:))
      cdn(:,1)=cdropl_sh
    else where (rlatt(:)>0.)
      cdn(:,1)=cdrops_nh
    else where
      cdn(:,1)=cdrops_sh
    end where
    do k=2,kl
      cdn(:,k)=cdn(:,1)
    enddo
end select

return
end subroutine aerodrop

end module aerointerface