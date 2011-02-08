
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
real, dimension(ifull) :: vt
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
  where (blon.lt.0.)
    blon=blon+360.
  end where
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
  wg=max(min((wb(:,1)-swilt(isoilm))/(sfc(isoilm)-swilt(isoilm)),1.),0.)
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

! calculate transfer velocity
call getvt(vt,zg(:,1))

! update prognostic aerosols
call aldrcalc(dt,sig,sigh,dsig,zg,cansto,mcmax,wg,pblh,ps,  &
              tss,t(1:ifull,:),condx,condc,snowd,sgsave,fg, &
              eg,u10,ustar,zo,land,fracice,sigmf,           &
              qlg(1:ifull,:),qfg(1:ifull,:),cfrac,clcon,    &
              pccw,dxy,rhoa,vt,ppfprec,ppfmelt,ppfsnow,     &
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

! This subroutine computes the transfer velocity
subroutine getvt(vt,zmin)
      
use arrays_m   ! for ps,t,qg,u,v
use pbl_m      ! for tss
use sigs_m     ! for sig
use work2_m    ! for wetfac,zo
      
implicit none

include 'const_phys.h'      
include 'newmpar.h'
include 'parm.h'
      
integer ic,iq
real, dimension(ifull), intent(out) :: vt
real, dimension(ifull), intent(in) :: zmin
real, dimension(ifull) :: smixr,thetav,thetavstar,tvs
real, dimension(ifull) :: integralh,lzom,lzoh,qsats
real, dimension(ifull) :: z_on_l,zt_on_l,ph0,ph1
real, dimension(ifull) :: integralm,z0_on_l,pm0,pm1,ustar,umag
real scrp
integer, parameter ::  nc     = 5
real, parameter    ::  vkar   = 0.4
real, parameter    ::  a_1    = 1.
real, parameter    ::  b_1    = 2./3.
real, parameter    ::  c_1    = 5.
real, parameter    ::  d_1    = 0.35
real, parameter    ::  lna    = 2.3
real, parameter    ::  z0     = 1.5
real, parameter    ::  z10    = 10.

! calculate qsat at surface
call getqsat(qsats,tss,ps(1:ifull))
smixr=wetfac*qsats+(1.-wetfac)*min(qsats,qg(:,1))

! calculate thetav at surface and 1st model level
scrp=sig(1)**(rdry/cp)
thetav=t(1:ifull,1)*(1.+0.61*qg(1:ifull,1))/scrp
tvs=tss*(1.+0.61*smixr)
      
! calculate wind speed at first model level
umag=sqrt(u(1:ifull,1)*u(1:ifull,1)+v(1:ifull,1)*v(1:ifull,1))
umag=max(umag,vmodmin)

! Roughness length for momentum and heat (i.e., log(z0m) and log(z0h))
lzom=log(zmin/zo)
lzoh=lna+lzom

! Dyer and Hicks approach 
thetavstar=vkar*(thetav-tvs)/lzoh ! 1st guess
ustar     =vkar*umag/lzom         ! 1st guess
do ic=1,nc
  z_on_l=vkar*zmin*grav*thetavstar/(thetav*ustar**2)
  z_on_l=min(z_on_l,10.)
  z0_on_l  = z_on_l*exp(-lzom)
  zt_on_l  = z_on_l*exp(-lzoh)
  where (z_on_l.lt.0.)
    pm0     = (1.-16.*z0_on_l)**(-0.25)
    pm1     = (1.-16.*z_on_l)**(-0.25)
    ph0     = (1.-16.*zt_on_l)**(-0.5)
    ph1     = (1.-16.*z_on_l)**(-0.5)
    integralm = lzom-2.*log((1.+1./pm1)/(1.+1./pm0))      &
                -log((1.+1./pm1**2)/(1.+1./pm0**2))       &
                +2.*(atan(1./pm1)-atan(1./pm0))
    integralh = lzoh-2.*log((1.+1./ph1)/(1.+1./ph0))
  elsewhere
    !--------------Beljaars and Holtslag (1991) momentum & heat            
    pm0 = -(a_1*z0_on_l+b_1*(z0_on_l-(c_1/d_1))*exp(-d_1*z0_on_l)+b_1*c_1/d_1)
    pm1 = -(a_1*z_on_l+b_1*(z_on_l-(c_1/d_1))*exp(-d_1*z_on_l)+b_1*c_1/d_1)
    ph0 = -((1.+(2./3.)*a_1*zt_on_l)**1.5+b_1*(zt_on_l-(c_1/d_1))*exp(-d_1*zt_on_l)+b_1*c_1/d_1-1.)
    ph1 = -((1.+(2./3.)*a_1*z_on_l)**1.5+b_1*(z_on_l-(c_1/d_1))*exp(-d_1*z_on_l)+b_1*c_1/d_1-1.)
    integralm = lzom-(pm1-pm0)    
    integralh = lzoh-(ph1-ph0)         
  endwhere
  thetavstar=vkar*(thetav-tvs)/integralh
  ustar     =vkar*umag/integralm
end do
      
vt=vkar*ustar/integralh
      
return
end subroutine getvt

subroutine getqsat(qsat,temp,ps)

implicit none

include 'newmpar.h'

real, dimension(ifull), intent(in) :: temp,ps
real, dimension(ifull), intent(out) :: qsat
real, dimension(0:220) :: table
real, dimension(ifull) :: esatf,tdiff,rx
integer, dimension(ifull) :: ix

table(0:4)=    (/ 1.e-9, 1.e-9, 2.e-9, 3.e-9, 4.e-9 /)                                !-146C
table(5:9)=    (/ 6.e-9, 9.e-9, 13.e-9, 18.e-9, 26.e-9 /)                             !-141C
table(10:14)=  (/ 36.e-9, 51.e-9, 71.e-9, 99.e-9, 136.e-9 /)                          !-136C
table(15:19)=  (/ 0.000000188, 0.000000258, 0.000000352, 0.000000479, 0.000000648 /)  !-131C
table(20:24)=  (/ 0.000000874, 0.000001173, 0.000001569, 0.000002090, 0.000002774 /)  !-126C
table(25:29)=  (/ 0.000003667, 0.000004831, 0.000006340, 0.000008292, 0.00001081 /)   !-121C
table(30:34)=  (/ 0.00001404, 0.00001817, 0.00002345, 0.00003016, 0.00003866 /)       !-116C
table(35:39)=  (/ 0.00004942, 0.00006297, 0.00008001, 0.0001014, 0.0001280 /)         !-111C
table(40:44)=  (/ 0.0001613, 0.0002026, 0.0002538, 0.0003170, 0.0003951 /)            !-106C
table(45:49)=  (/ 0.0004910, 0.0006087, 0.0007528, 0.0009287, 0.001143 /)             !-101C
table(50:55)=  (/ .001403, .001719, .002101, .002561, .003117, .003784 /)             !-95C
table(56:63)=  (/ .004584, .005542, .006685, .008049, .009672,.01160,.01388,.01658 /) !-87C
table(64:72)=  (/ .01977, .02353, .02796,.03316,.03925,.04638,.05472,.06444,.07577 /) !-78C
table(73:81)=  (/ .08894, .1042, .1220, .1425, .1662, .1936, .2252, .2615, .3032 /)   !-69C
table(82:90)=  (/ .3511, .4060, .4688, .5406, .6225, .7159, .8223, .9432, 1.080 /)    !-60C
table(91:99)=  (/ 1.236, 1.413, 1.612, 1.838, 2.092, 2.380, 2.703, 3.067, 3.476 /)    !-51C
table(100:107)=(/ 3.935,4.449, 5.026, 5.671, 6.393, 7.198, 8.097, 9.098 /)            !-43C
table(108:116)=(/ 10.21, 11.45, 12.83, 14.36, 16.06, 17.94, 20.02, 22.33, 24.88 /)    !-34C
table(117:126)=(/ 27.69, 30.79, 34.21, 37.98, 42.13, 46.69,51.70,57.20,63.23,69.85 /) !-24C 
table(127:134)=(/ 77.09, 85.02, 93.70, 103.20, 114.66, 127.20, 140.81, 155.67 /)      !-16C
table(135:142)=(/ 171.69, 189.03, 207.76, 227.96 , 249.67, 272.98, 298.00, 324.78 /)  !-8C
table(143:150)=(/ 353.41, 383.98, 416.48, 451.05, 487.69, 526.51, 567.52, 610.78 /)   !0C
table(151:158)=(/ 656.62, 705.47, 757.53, 812.94, 871.92, 934.65, 1001.3, 1072.2 /)   !8C
table(159:166)=(/ 1147.4, 1227.2, 1311.9, 1401.7, 1496.9, 1597.7, 1704.4, 1817.3 /)   !16C
table(167:174)=(/ 1936.7, 2063.0, 2196.4, 2337.3, 2486.1, 2643.0, 2808.6, 2983.1 /)   !24C
table(175:182)=(/ 3167.1, 3360.8, 3564.9, 3779.6, 4005.5, 4243.0, 4492.7, 4755.1 /)   !32C
table(183:190)=(/ 5030.7, 5320.0, 5623.6, 5942.2, 6276.2, 6626.4, 6993.4, 7377.7 /)   !40C
table(191:197)=(/ 7780.2, 8201.5, 8642.3, 9103.4, 9585.5, 10089.0, 10616.0 /)         !47C
table(198:204)=(/ 11166.0, 11740.0, 12340.0, 12965.0, 13617.0, 14298.0, 15007.0 /)    !54C
table(205:211)=(/ 15746.0, 16516.0, 17318.0, 18153.0, 19022.0, 19926.0, 20867.0 /)    !61C
table(212:218)=(/ 21845.0, 22861.0, 23918.0, 25016.0, 26156.0, 27340.0, 28570.0 /)    !68C
table(219:220)=(/ 29845.0, 31169.0 /)

tdiff=min(max( temp-123.16, 0.), 219.)
rx=tdiff-aint(tdiff)
ix=int(tdiff)
esatf=(1.-rx)*table(ix)+ rx*table(ix+1)
qsat=0.622*esatf/(ps-esatf)

return
end subroutine getqsat

end module aerointerface