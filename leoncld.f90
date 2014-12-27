module leoncld_mod
    
private
public leoncld

contains
    
! This subroutine is the interface for the LDR cloud microphysics
      
subroutine leoncld(cfrac,rfrac)
      
use aerointerface
use arrays_m
use cc_mpi, only : mydiag, myid
use cloudmod
use diag_m
use estab
use kuocomb_m
use latlong_m
use liqwpar_m  ! ifullw
use morepbl_m
use nharrs_m
use nlin_m
use prec_m
use sigs_m
use soil_m     ! land
use tracers_m  ! ngas, nllp, ntrac
use vvel_m
use work3f_m
      
implicit none
      
include 'newmpar.h'
include 'const_phys.h' !Input physical constants
include 'cparams.h'    !Input cloud scheme parameters
include 'kuocom.h'     ! acon,bcon,Rcm
include 'parm.h'
include 'params.h'
      
! Local variables
integer iq,k,ncl
real ccw !Stuff for convective cloud fraction
real fl(ifull)

real prf(ifull,kl)     !Pressure on full levels (hPa)
real dprf(ifull,kl)    !Pressure thickness (hPa)
real rhoa(ifull,kl)    !Air density (kg/m3)
real dz(ifull,kl)      !Layer thickness (m)
real cdso4(ifull,kl)   !Cloud droplet conc (#/m3)
real cfrac(ifull,kl)   !Cloud fraction (passed back to globpe)
real rfrac(ifull+iextra,kl)  !Rain fraction (passed back to globpe)
real ccov(ifull,kl)    !Cloud cover (may differ from cloud frac if vertically subgrid)
real cfa(ifull,kl)     !Cloud fraction in which autoconv occurs (option in newrain.f)
real qca(ifull,kl)     !Cloud water mixing ratio in cfa(:,:)    (  "    "     "     )
real fluxc(ifull,kl)   !Flux of convective rainfall in timestep (kg/m**2)
real ccrain(ifull,kl)  !Convective raining cloud cover
real precs(ifull)      !Amount of stratiform precipitation in timestep (mm)
real preci(ifull)      !Amount of stratiform snowfall in timestep (mm)
real wcon(ifull)       !Convective cloud water content (in-cloud, prescribed)
real clcon(ifull,kl)   !Convective cloud fraction in layer 
real qsg(ifull,kl)     !Saturation mixing ratio
real qcl(ifull,kl)     !Vapour mixing ratio inside convective cloud
real qenv(ifull,kl)    !Vapour mixing ratio outside convective cloud
real tenv(ifull,kl)    !Temperature outside convective cloud
real tnhs(ifull,kl)    !Non-hydrostatic temperature adjusement
real tv(ifull,kl)

! These outputs are not used in this model at present
real qevap(ifull,kl)
real qsubl(ifull,kl)
real qauto(ifull,kl)
real qcoll(ifull,kl)
real qaccr(ifull,kl)
real fluxr(ifull,kl)
real fluxi(ifull,kl)
real fluxm(ifull,kl)
real pqfsed(ifull,kl)
real pfstayice(ifull,kl)
real pfstayliq(ifull,kl)
real slopes(ifull,kl)
real prscav(ifull,kl)

integer kbase(ifull),ktop(ifull) !Bottom and top of convective cloud 

! set-up params.h
ln2=ifull
nl=kl
nlp=nl+1
nlm=nl-1

do k=1,kl   
  prf(:,k)=0.01*ps(1:ifull)*sig(k)    !ps is SI units
  dprf(:,k)=-0.01*ps(1:ifull)*dsig(k) !dsig is -ve
  tv(:,k)=t(1:ifull,k)*(1.+0.61*qg(1:ifull,k)-qlg(1:ifull,k)-qfg(1:ifull,k))
  rhoa(:,k)=100.*prf(:,k)/(rdry*tv(1:ifull,k))
  do iq=1,ifull
    qsg(iq,k)=qsat(100.*prf(iq,k),t(iq,k))
  enddo
enddo
      
! Non-hydrostatic terms
tnhs(:,1)=phi_nh(:,1)/bet(1)
do k=2,kl
  ! representing non-hydrostatic term as a correction to air temperature
  tnhs(:,k)=(phi_nh(:,k)-phi_nh(:,k-1)-betm(k)*tnhs(:,k-1))/bet(k)
end do
      
! Calculate droplet concentration from aerosols (for non-convective faction of grid-box)
call aerodrop(1,ifull,kl,cdso4,rhoa,land,rlatt,outconv=.true.)

kbase(:)=0  ! default
ktop(:) =0  ! default
dz(:,:)=100.*dprf/(rhoa*grav)*(1.+tnhs(1:ifull,:)/tv)
fluxc(:,:)=0.
ccrain(:,:)=0.1  !Assume this for now
precs(:)=0.
preci(:)=0.

!     Set up convective cloud column
call convectivecloudfrac(clcon)
where ( ktsav<kl-1 )
  ktop   = ktsav
  kbase  = kbsav+1
  wcon   = wlc
elsewhere
  wcon   = 0.
end where

if(nmaxpr==1.and.mydiag)then
  if(ktau==1) then
    write(6,*)'in leoncloud acon,bcon,Rcm ',acon,bcon,Rcm
  end if
  write(6,*) 'entering leoncld'
  write (6,"('qg  ',9f8.3/4x,9f8.3)") qg(idjd,:)
  write (6,"('qf  ',9f8.3/4x,9f8.3)") qfg(idjd,:)
  write (6,"('ql  ',9f8.3/4x,9f8.3)") qlg(idjd,:)
endif

!     Calculate convective cloud fraction and adjust moisture variables 
!     before calling newcloud
if ( ncloud<=3 ) then
  ! diagnose cloud fraction
  if ( nmr>0 ) then ! Max/Rnd cloud overlap
    do k=1,kl
      where ( clcon(:,k)>0. )
        !ccw=wcon(:)/rhoa(:,k)  !In-cloud l.w. mixing ratio
        qccon(:,k) = clcon(:,k)*wcon(:)/rhoa(:,k)
        qcl(:,k)   = max(qsg(:,k),qg(1:ifull,k))  ! jlm
        qenv(1:ifull,k) = max(1.e-8,qg(1:ifull,k)-clcon(:,k)*qcl(:,k))/(1.-clcon(:,k))
        qcl(:,k) = (qg(1:ifull,k)-(1.-clcon(:,k))*qenv(1:ifull,k))/clcon(:,k)
        qlg(1:ifull,k) = qlg(1:ifull,k)/(1.-clcon(:,k))
        qfg(1:ifull,k) = qfg(1:ifull,k)/(1.-clcon(:,k))
      elsewhere
        clcon(:,k)      = 0.
        qccon(:,k)      = 0.
        qcl(1:ifull,k)  = 0.
        qenv(1:ifull,k) = qg(1:ifull,k)
      end where
    end do
        
  else ! usual random cloud overlap
    do k=1,kl
      do iq=1,ifull
        if(clcon(iq,k)>0.)then
          ccw=wcon(iq)/rhoa(iq,k)  !In-cloud l.w. mixing ratio
          qccon(iq,k)=clcon(iq,k)*ccw
          qcl(iq,k)=max(qsg(iq,k),qg(iq,k))  ! jlm
          qenv(iq,k)=max(1.e-8,qg(iq,k)-clcon(iq,k)*qcl(iq,k))/(1.-clcon(iq,k))
          qcl(iq,k)=(qg(iq,k)-(1.-clcon(iq,k))*qenv(iq,k))/clcon(iq,k)
          qlg(iq,k)=qlg(iq,k)/(1.-clcon(iq,k))
          qfg(iq,k)=qfg(iq,k)/(1.-clcon(iq,k))
        else
          clcon(iq,k)=0.
          qccon(iq,k)=0.
          qcl(iq,k)=0.
          qenv(iq,k)=qg(iq,k)
        endif
      enddo
    enddo
      
  end if
else
  ! prognostic cloud fraction
  clcon(:,:)      = 0.
  qccon(:,:)      = 0.
  qcl(1:ifull,:)  = 0.
  qenv(1:ifull,:) = qg(1:ifull,:)
end if
      
tenv(:,:)=t(1:ifull,:) !Assume T is the same in and out of convective cloud

if ( nmaxpr==1 .and. mydiag ) then
  write(6,*) 'before newcloud',ktau
  write (6,"('t   ',9f8.2/4x,9f8.2)") t(idjd,:)
  write (6,"('qg  ',9f8.3/4x,9f8.3)") qg(idjd,:)
  write (6,"('qf  ',9f8.3/4x,9f8.3)") qfg(idjd,:)
  write (6,"('ql  ',9f8.3/4x,9f8.3)") qlg(idjd,:)
  write (6,"('qnv ',9f8.3/4x,9f8.3)") qenv(idjd,:)
  write (6,"('qsg ',9f8.3/4x,9f8.3)") qsg(idjd,:)
  write (6,"('qcl ',9f8.3/4x,9f8.3)") qcl(idjd,:)
  write (6,"('clc ',9f8.3/4x,9f8.3)") clcon(idjd,:)
  write(6,*) 'kbase,ktop ',kbase(idjd),ktop(idjd)
endif

!     Calculate cloud fraction and cloud water mixing ratios
call newcloud(dt,land,prf,kbase,ktop,rhoa,cdso4,tenv,qenv,qlg,qfg,cfrac,ccov,cfa,qca)

if(nmaxpr==1.and.mydiag)then
  write(6,*) 'after newcloud',ktau
  write (6,"('tnv ',9f8.2/4x,9f8.2)") tenv(idjd,:)
  write (6,"('qg0 ',9f8.3/4x,9f8.3)") qg(idjd,:)
  write (6,"('qf  ',9f8.3/4x,9f8.3)") qfg(idjd,:)
  write (6,"('ql  ',9f8.3/4x,9f8.3)") qlg(idjd,:)
  write (6,"('qnv ',9f8.3/4x,9f8.3)") qenv(idjd,:) ! really new qg
endif

!     Weight output variables according to non-convective fraction of grid-box            
do k=1,kl
  t(1:ifull,k)=clcon(:,k)*t(1:ifull,k)+(1.-clcon(:,k))*tenv(:,k)
  qg(1:ifull,k)=clcon(:,k)*qcl(:,k)+(1.-clcon(:,k))*qenv(:,k)
  where ( k>=kbase .and. k<=ktop )
    cfrac(:,k)=cfrac(:,k)*(1.-clcon(:,k))
    ccov(:,k)=ccov(:,k)*(1.-clcon(:,k))              
    qlg(1:ifull,k)=qlg(1:ifull,k)*(1.-clcon(:,k))
    qfg(1:ifull,k)=qfg(1:ifull,k)*(1.-clcon(:,k))
    cfa(:,k)=cfa(:,k)*(1.-clcon(:,k))
    qca(:,k)=qca(:,k)*(1.-clcon(:,k))              
  end where
enddo

if(nmaxpr==1.and.mydiag)then
  write(6,*) 'before newrain',ktau
  write (6,"('t   ',9f8.2/4x,9f8.2)") t(idjd,:)
  write (6,"('qg  ',9f8.3/4x,9f8.3)") qg(idjd,:)
  write (6,"('qf  ',9f8.3/4x,9f8.3)") qfg(idjd,:)
  write (6,"('ql  ',9f8.3/4x,9f8.3)") qlg(idjd,:)
endif
if(diag)then
  call maxmin(t,' t',ktau,1.,kl)
  call maxmin(qg,'qg',ktau,1.e3,kl)
  call maxmin(qfg,'qf',ktau,1.e3,kl)
  call maxmin(qlg,'ql',ktau,1.e3,kl)
endif

!     Add convective cloud water into fields for radiation
!     cfrad replaced by updating cfrac Oct 2005
!     Moved up 16/1/06 (and ccov,cfrac NOT UPDATED in newrain)
!     done because sometimes newrain drops out all qlg, ending up with 
!     zero cloud (although it will be rediagnosed as 1 next timestep)
do k=1,kl
  fl=max(0.,min(1.,(t(1:ifull,k)-ticon)/(273.15-ticon)))
  qlrad(:,k)=qlg(:,k)+fl*qccon(:,k)
  qfrad(:,k)=qfg(:,k)+(1.-fl)*qccon(:,k)
enddo

!     Calculate precipitation and related processes
call newrain(land,dt,fluxc,rhoa,dz,ccrain,prf,cdso4,cfa,qca,t,qlg,qfg,qrg, &
          precs,qg,cfrac,rfrac,ccov,preci,qevap,qsubl,qauto,qcoll,qaccr,   &
          fluxr,fluxi,fluxm,pfstayice,pfstayliq,pqfsed,slopes,prscav)

if(nmaxpr==1.and.mydiag)then
  write(6,*) 'after newrain',ktau
  write (6,"('t   ',9f8.2/4x,9f8.2)") t(idjd,:)
  write (6,"('qg  ',9f8.3/4x,9f8.3)") qg(idjd,:)
  write (6,"('qf  ',9f8.3/4x,9f8.3)") qfg(idjd,:)
  write (6,"('ql  ',9f8.3/4x,9f8.3)") qlg(idjd,:)
end if
if(diag)then
  call maxmin(t,' t',ktau,1.,kl)
  call maxmin(qg,'qg',ktau,1.e3,kl)
  call maxmin(qfg,'qf',ktau,1.e3,kl)
  call maxmin(qlg,'ql',ktau,1.e3,kl)
endif

!--------------------------------------------------------------
! Store data needed by prognostic aerosol scheme
if (abs(iaero)>=2) then
  ppfprec(:,1)=0. !At TOA
  ppfmelt(:,1)=0. !At TOA
  ppfsnow(:,1)=0. !At TOA
  do k=1,kl-1
    ppfprec(:,kl+1-k)=(fluxr(:,k+1)+fluxm(:,k))/dt !flux *entering* layer k
    ppfmelt(:,kl+1-k)=fluxm(:,k)/dt                !flux melting in layer k
    ppfsnow(:,kl+1-k)=(fluxi(:,k+1)-fluxm(:,k))/dt !flux *entering* layer k
  enddo
  do k=1,kl
    ppfevap(:,kl+1-k)=qevap(:,k)*rhoa(:,k)*dz(:,k)/dt
    ppfsubl(:,kl+1-k)=qsubl(:,k)*rhoa(:,k)*dz(:,k)/dt !flux sublimating or staying in k
    pplambs(:,kl+1-k)=slopes(:,k)
    ppmrate(:,kl+1-k)=(qauto(:,k)+qcoll(:,k))/dt
    ppmaccr(:,kl+1-k)=qaccr(:,k)/dt
    ppfstayice(:,kl+1-k)=pfstayice(:,k)
    ppfstayliq(:,kl+1-k)=pfstayliq(:,k)
    ppqfsed(:,kl+1-k)=pqfsed(:,k)
    pprscav(:,kl+1-k)=prscav(:,k)
  enddo
end if
!--------------------------------------------------------------

!     Add convective cloud water into fields for radiation
!     cfrad replaced by updating cfrac Oct 2005
!     Moved up 16/1/06 (and ccov,cfrac NOT UPDATED in newrain)
!     done because sometimes newrain drops out all qlg, ending up with 
!     zero cloud (although it will be rediagnosed as 1 next timestep)
do k=1,kl
  cfrac(:,k)=min(1.,ccov(:,k)+clcon(:,k))
enddo

!========================= Jack's diag stuff =========================
!if(ncfrp==1)then  ! from here to near end; Jack's diag stuff
!  do iq=1,icfrp
!    tautot(iq)=0.
!    cldmax(iq)=0.
!    ctoptmp(iq)=0.
!    ctoppre(iq)=0.
!    do k=1,kl
!      fice(iq,k)=0.
!    enddo
!    kcldfmax(iq)=0.
!  enddo
!!      cfrp data
!  do k=1,kl-1
!    do iq=1,icfrp
!      taul(iq,k)=0.
!      taui(iq,k)=0.
!      Reffl=0.
!      if(cfrac(iq,k)>0.)then
!        tau_sfac=1.
!        fice(iq,k) = qfrad(iq,k)/(qfrad(iq,k)+qlrad(iq,k)) ! 16/1/06
!!            Liquid water clouds
!        if(qlg(iq,k)>1.0e-8)then
!          Wliq=rhoa(iq,k)*qlg(iq,k)/(cfrac(iq,k)*(1-fice(iq,k))) !kg/m^3
!          if(.not.land(iq))then !sea
!            rk=0.8
!          else            !land
!            rk=0.67
!          endif
!! Reffl is the effective radius at the top of the cloud (calculated following
!! Martin etal 1994, JAS 51, 1823-1842) due to the extra factor of 2 in the
!! formula for reffl. Use mid cloud value of Reff for emissivity.
!          Reffl=(3*2*Wliq/(4*rhow*pi*rk*cdso4(iq,k)))**(1./3)
!          qlpath=Wliq*dz(iq,k)
!          taul(iq,k)=tau_sfac*1.5*qlpath/(rhow*Reffl)
!        endif ! qlg
!! Ice clouds
!        if(qfg(iq,k)>1.0e-8)then
!          Wice=rhoa(iq,k)*qfg(iq,k)/(cfrac(iq,k)*fice(iq,k)) !kg/m**3
!          sigmai = aice*Wice**bice !visible ext. coeff. for ice
!          taui(iq,k)=sigmai*dz(iq,k) !visible opt. depth for ice
!          taui(iq,k)=tau_sfac*taui(iq,k)
!        endif ! qfg
!      endif !cfrac
!    enddo ! iq
!  enddo ! k
!! Code to get vertically integrated value...
!! top down to get highest level with cfrac=cldmax (kcldfmax)
!  do k=kl-1,1,-1
!    do iq=1,icfrp
!      tautot(iq)=tautot(iq)+cfrac(iq,k)*(fice(iq,k)*taui(iq,k)+(1.-fice(iq,k))*taul(iq,k))
!      if(cfrac(iq,k)>cldmax(iq)) kcldfmax(iq)=k
!      cldmax(iq)=max(cldmax(iq),cfrac(iq,k))
!    enddo ! iq
!  enddo ! k
!
!  do iq=1,icfrp
!    if(cldmax(iq)>1.e-10) then
!      tautot(iq)=tautot(iq)/cldmax(iq)
!
!      cfd=0.
!      do k=kl,kcldfmax(iq),-1
!        fcf = max(0.,cfrac(iq,k)-cfd) ! cld frac. from above
!        ctoptmp(iq)=ctoptmp(iq)+fcf*t(iq,k)/cldmax(iq)
!        ctoppre(iq)=ctoppre(iq)+fcf*prf(iq,k)/cldmax(iq)
!        cfd=max(cfrac(iq,k),cfd)
!      enddo ! k=kl,kcldfmax(iq),-1
!
!    endif ! (cldmax(iq).gt.1.e-10) then
!  enddo   ! iq
!endif    ! ncfrp.eq.1
!========================= end of Jack's diag stuff ======================

!     factor of 2 is because used .5 in newrain.f (24-mar-2000, jjk)
do iq=1,ifullw        
  condx(iq)=condx(iq)+precs(iq)*2.
  conds(iq)=conds(iq)+preci(iq)*2.
  precip(iq)=precip(iq)+precs(iq)*2.
enddo

return
end subroutine leoncld


! from argumen
!      ttg - temperature (K)
!      qtg - water vapour mixing ratio (kg/kg) - called qenv in leoncld
!      qlg - cloud liquid water mixing ratio (kg/kg)
!      qfg - cloud ice mixing ratio (kg/kg)
!
! Output:
!
! from arguments
!      cfrac - cloudy fraction of grid box
!      ccov - cloud cover looking from above (currently = cloud fraction)
! 
!******************************************************************************

 subroutine newcloud(tdt,land,prf,kbase,ktop,rhoa,cdrop,ttg,qtg,qlg,qfg,cfrac,ccov,cfa,qca)

! This routine is part of the prognostic cloud water scheme

use cc_mpi, only : mydiag
use cloudmod
use diag_m      
use estab, only : esdiffx, qsati
use map_m
use sigs_m

implicit none

! Global parameters
include 'newmpar.h'
include 'const_phys.h' ! Input physical constants
include 'cparams.h'    ! Input cloud scheme parameters
include 'kuocom.h'     ! Input cloud scheme parameters rcrit_l & rcrit_s
include 'params.h'     ! Input model grid dimensions (modified params.h for CCAM)
include 'parm.h'

! Argument list
real tdt,den1
logical land(ifull)
integer kbase(ifull)
integer ktop(ifull)
real prf(ifull,nl)
real rhoa(ifull,nl)
real ttg(ifull,nl)
real qtg(ifull,nl)
real qlg(ifull+iextra,nl)
real qfg(ifull+iextra,nl)
real cfrac(ifull,nl)
real ccov(ifull,nl)
real cfa(ifull,nl)
real qca(ifull,nl)
real qsl(ifull,nl)
real qsw(ifull,nl)

! Local work arrays and variables
real qcg(ifull,nl)
real qtot(ifull,nl),tliq(ifull,nl)
real fice(ifull,nl)
real qcold(ifull,nl)
real cdrop(ifull,nl)
real rcrit(ifull,nl)
real qsi(ifull,nl)
real qfnew(ifull,nl)

integer k
integer mg

real al
real alf
real aprpr
real bprpr
real cice
real cm0
real crate
real decayfac
real deles
real delq
real dqsdt
real es
real fd
real fl
real hlrvap
real pk
real qc
real qfdep
real qi0
real qs
real rhoic
real root6
real root6i
real tk
real qcic,qcrit,qc2,qto,wliq,r3c,r6c,eps,beta6

! Start code : ----------------------------------------------------------

root6=sqrt(6.)
root6i=1./root6

if(diag.and.mydiag)then
  write(6,*) 'entering newcloud'
  write(6,'(a,30f10.3)') 'prf ',(prf(idjd,k),k=1,nl)
  write(6,'(a,30f10.3)') 'ttg ',(ttg(idjd,k),k=1,nl)
  write(6,*) 'qtg ',(qtg(idjd,k),k=1,nl)
  write(6,*) 'qlg ',(qlg(idjd,k),k=1,nl)
  write(6,*) 'qfg ',(qfg(idjd,k),k=1,nl)
endif

! First melt cloud ice or freeze cloud water to give correct ice fraction fice.
! Then calculate the cloud conserved variables qtot and tliq.
! Note that qcg is the total cloud water (liquid+frozen)

where(ttg>=tfrz)
  fice=0.
elsewhere(ttg>=tice.and.qfg(1:ifull,:)>1.e-12)
  fice=min(qfg(1:ifull,:)/(qfg(1:ifull,:)+qlg(1:ifull,:)),1.)
elsewhere(ttg>=tice)
  fice=0.
elsewhere
  fice=1.
end where
qcg(:,:)=qlg(1:ifull,:)+qfg(1:ifull,:)
qcold(:,:)=qcg(:,:)
qfnew=fice(:,:)*qcg(:,:)
ttg(:,:)=ttg(:,:)+hlfcp*(qfnew-qfg(1:ifull,:)) !Release L.H. of fusion
qfg(1:ifull,:)=fice(:,:)*qcg(:,:)
qlg(1:ifull,:)=max(0.,qcg(:,:)-qfg(1:ifull,:))

if(diag.and.mydiag)then
  write(6,*) 'within newcloud'
  write(6,*) 'ttg ',ttg(idjd,:)
  write(6,*) 'qcold ',qcold(idjd,:)
  write(6,*) 'qcg ',qcg(idjd,:)
  write(6,*) 'qlg ',qlg(idjd,:)
  write(6,*) 'qfg ',qfg(idjd,:)
  write(6,*) 'fice ',fice(idjd,:)
endif
      
if (ncloud<3) then ! usual diagnostic cloud
      
! Precompute the array of critical relative humidities 

  if(nclddia==-3)then
    do k=1,nl
      where(land(1:ifull))
        rcrit(:,k)=max( rcrit_l, (1.-16.*(1.-sig(k))**3) )
      elsewhere
        rcrit(:,k)=max( rcrit_s, (1.-16.*(1.-sig(k))**3) )
      end where
    enddo
  elseif(nclddia<0)then
    do k=1,nl
      where(land(1:ifull))
        rcrit(:,k)=max( rcrit_l, (1.-4.*(1.-sig(k))**2) )
      elsewhere
        rcrit(:,k)=max( rcrit_s, (1.-4.*(1.-sig(k))**2) )
      end where
    enddo
  elseif(nclddia==1)then
    do k=1,nl
      where(land(1:ifull))
        rcrit(:,k)=max( rcrit_l, sig(k)**3 )
      elsewhere
        rcrit(:,k)=max( rcrit_s, sig(k)**3 )
      end where
    enddo
  elseif(nclddia==2)then
    do k=1,nl
      where(land(1:ifull))
        rcrit(:,k)=rcrit_l
      elsewhere
        rcrit(:,k)=rcrit_s
      end where
    enddo
  elseif(nclddia==3)then
    do k=1,nl
      where(land(1:ifull))
        rcrit(:,k)=max( rcrit_l, sig(k)**2 )          ! .75 for R21 Mk2
      elsewhere
        rcrit(:,k)=max( rcrit_s, sig(k)**2 )          ! .85 for R21 Mk2
      end where
    enddo
  elseif(nclddia==4)then
    do k=1,nl
      where(land(1:ifull))
        rcrit(:,k)=max( rcrit_l, sig(k) )             ! .75 for test Mk2/3
      elsewhere
        rcrit(:,k)=max( rcrit_s, sig(k) )             ! .9  for test Mk2/3
      end where
    enddo
  elseif(nclddia==5)then  ! default till May 08
    do k=1,nl
      where(land(1:ifull))
        rcrit(:,k)=max( rcrit_l, min(.99,sig(k)) )    ! .75 for same as T63
      elsewhere
        rcrit(:,k)=max( rcrit_s, min(.99,sig(k)) )    ! .85 for same as T63
      end where
    enddo
  elseif(nclddia==6)then
    do k=1,nl
      where(land(1:ifull))
        rcrit(:,k)=max(rcrit_l*(1.-.15*sig(k)),sig(k)**4)
      elsewhere
        rcrit(:,k)=max(rcrit_s*(1.-.15*sig(k)),sig(k)**4)
      end where
    enddo
  elseif(nclddia==7)then
    do k=1,nl
      where(land(1:ifull))
        rcrit(:,k)=max(rcrit_l*(1.-.2*sig(k)),sig(k)**4)
      elsewhere
        rcrit(:,k)=max(rcrit_s*(1.-.2*sig(k)),sig(k)**4)
      end where
    enddo
  elseif(nclddia>7)then  ! e.g. 12    JLM
    do k=1,nl  ! typically set rcrit_l=.75,  rcrit_s=.85
      do mg=1,ifull
        tk=em(mg)*208498./ds
        fl=(1+nclddia)*tk/(1.+nclddia*tk)
!         for rcit_l=.75 & nclddia=12 get rcrit=(.797,.9, .9375, .971, .985) for (50,10, 5, 2, 1) km
        if(land(mg))then
          rcrit(mg,k)=max(1.-fl*(1.-rcrit_l),sig(k)**3)        
        else
          rcrit(mg,k)=max(1.-fl*(1.-rcrit_s),sig(k)**3)         
        endif
      enddo
    enddo
  endif  ! (nclddia<0)  .. else ..

! Calculate cloudy fraction of grid box (cfrac) and gridbox-mean cloud water
! using the triangular PDF of Smith (1990)

  do k=1,nl
    do mg=1,ifull
      hlrvap=(hl+fice(mg,k)*hlf)/rvap
      qtot(mg,k)=qtg(mg,k)+qcg(mg,k)
      tliq(mg,k)=ttg(mg,k)-hlcp*qcg(mg,k)-hlfcp*qfg(mg,k)

! Calculate qs and gam=(L/cp)*dqsdt,  at temperature tliq
      pk=100.0*prf(mg,k)
      qsi(mg,k)=qsati(pk,tliq(mg,k))           !Ice value
      deles=esdiffx(tliq(mg,k))                              ! MJT suggestion
      qsl(mg,k)=qsi(mg,k)+epsil*deles/pk       !qs over liquid
      qsw(mg,k)=fice(mg,k)*qsi(mg,k)+(1.-fice(mg,k))*qsl(mg,k) !Weighted qs at temperature Tliq
      qs=qsw(mg,k)
      dqsdt=qs*hlrvap/tliq(mg,k)**2
!     qvc(mg,k)=qs !Vapour mixing ratio in cloud

      al=1./(1.+(hlcp+fice(mg,k)*hlfcp)*dqsdt)    !Smith's notation
      qc=qtot(mg,k)-qs

      delq=(1.-rcrit(mg,k))*qs      !UKMO style (equivalent to above)
      cfrac(mg,k)=1.
      qcg(mg,k)=al*qc
      if(qc<delq)then
        cfrac(mg,k)=max(1.e-6 , 1.-.5*((qc-delq)/delq)**2)     ! for roundoff
        qcg(mg,k)=max(1.e-8,al*(qc-(qc-delq)**3/(6.*delq**2))) ! for roundoff
      endif
      if(qc<=0.)then
        cfrac(mg,k)=max(1.e-6 , .5*((qc+delq)/delq)**2)    ! for roundoff
        qcg(mg,k)=max(1.e-8, al*(qc+delq)**3/(6.*delq**2)) ! for roundoff
      endif
      if(qc<=-delq)then
        cfrac(mg,k)=0.
        qcg(mg,k)=0.
      endif

! Calculate the cloud fraction (cfa) in which ql exceeds qcrit, and
! the corresponding gridbox-mean cloud water mixing ratio qca. 
! This (qca) is the cloud-water mixing ratio inside cfa times cfa.
! The new variable qc2 is like qc above, but is used for integration limits
! only, not the integrand

      if(cfrac(mg,k)>0.)then
        qcic=qcg(mg,k)/cfrac(mg,k) !Mean in cloud value

! Following few lines are for Yangang Liu's new scheme (2004: JAS, GRL)
! Need to do first-order estimate of qcrit using mean in-cloud qc (qcic)

        Wliq = max( 1.e-10, 1000. * qcic * rhoa(mg,k)) !g/m3
        R6c = 4.09e-4 * ( 1.15e23*1.e-6*cdrop(mg,k) / Wliq**2 ) ** (1./6.)
        eps = 1. - 0.7 * exp(-0.003e-6*cdrop(mg,k)) !mid range
        beta6 = ((1.+3.*eps**2)*(1.+4.*eps**2)*(1.+5.*eps**2) / ((1.+eps**2)*(1.+2.*eps**2)) )**(1./6.)
        R3c = 1.e-6*R6c/beta6 !in metres
        qcrit=(4.*pi/3.)*rhow*R3c**3*cdrop(mg,k)/rhoa(mg,k) !New qcrit

        qc2=qtot(mg,k)-qs-qcrit/al
        cfa(mg,k)=1.
        qca(mg,k)=al*qc
        if(qc2<delq)then
          cfa(mg,k)=1.-0.5*((qc2-delq)/delq)**2
          qto=(qtot(mg,k)-delq+2.*(qs+qcrit/al))/3.
          qca(mg,k)=al*(qtot(mg,k) - qto + cfa(mg,k)*(qto-qs))
        endif
        if(qc2<=0.)then
          cfa(mg,k)=0.5*((qc2+delq)/delq)**2
          qca(mg,k)=cfa(mg,k)*(al/3.)*(2.*qcrit/al + qc+delq)
        endif
        if(qc2<=-delq)then
          cfa(mg,k)=0.
          qca(mg,k)=0.
        endif
      else
        cfa(mg,k)=0.
        qca(mg,k)=0.
      endif

    enddo
  enddo

  if(diag.and.mydiag)then
    write(6,*) 'rcrit ',rcrit(idjd,:)
    write(6,*) 'qtot ',qtot(idjd,:)
    write(6,*) 'qsi',qsi(mg,k)
    write(6,*) 'tliq',tliq(idjd,:)
    write(6,*) 'qsl ',qsl(idjd,:)
    write(6,*) 'qsw ',qsw(idjd,:)
    write(6,*) 'cfrac ',cfrac(idjd,:)
    write(6,*) 'qc  ',  qtot(idjd,:)-qsw(idjd,:)
    write(6,*) 'qcg ',qcg(idjd,:)
    write(6,*) 'delq ', (1.-rcrit(idjd,:))*qsw(idjd,:)
  endif

! Assume condensation or evaporation retains ice fraction fice.
! Introduce a time-decay factor for cirrus (as suggested by results of Khvorostyanov & Sassen,
! JAS, 55, 1822-1845, 1998). Their suggested range for the time constant is 0.5 to 2 hours.
! The grid-box-mean values of qtg and ttg are adjusted later on (below).

  decayfac = exp ( -tdt/7200. )              ! Try 2 hrs
  !decayfac = 0.                             ! Instant adjustment (old scheme)
  where(ttg>=Tice)
    qfg(1:ifull,:) = fice*qcg
    qlg(1:ifull,:) = qcg - qfg(1:ifull,:)
  elsewhere                                    ! Cirrus T range
    qfg(1:ifull,:) = qcold*decayfac + qcg*(1.-decayfac)
    qlg(1:ifull,:) = 0.
    qcg(1:ifull,:) = qfg(1:ifull,:)
  end where
  
else ! prognostic cloud
      
  ! Tiedtke prognostic cloud model
  ! MJT notes - we use ttg instead of tliq
  qtot(:,:)=qtg(:,:)+qcg(:,:)
  tliq(:,:)=ttg(:,:)-hlcp*qcg(:,:)-hlfcp*qfg(1:ifull,:)
  do k=1,nl
    do mg=1,ifull
      pk=100.0*prf(mg,k)
      qsi(mg,k)=qsati(pk,ttg(mg,k))      ! Ice value
      deles=esdiffx(ttg(mg,k))
      qsl(mg,k)=qsi(mg,k)+epsil*deles/pk ! Liquid value
    end do
  end do
  qsw(:,:)=fice*qsi+(1.-fice)*qsl        ! Weighted qs at temperature Tliq
  call progcloud(cfrac,qcg,qtot,prf,rhoa,fice,qsw,ttg)
        
  ! Use 'old' autoconversion with prognostic cloud
  cfa(:,:)=0.
  qca(:,:)=0.

  where(ttg>=Tice)
    qfg(1:ifull,:) = fice*qcg
    qlg(1:ifull,:) = qcg - qfg(1:ifull,:)
  elsewhere
    qfg(1:ifull,:) = qcg
    qlg(1:ifull,:) = 0.
    qcg(1:ifull,:) = qfg(1:ifull,:)
  end where
  
end if ! ncloud<3 ..else..


! Do the vapour deposition calculation in mixed-phase clouds:
! Calculate deposition on cloud ice, assuming es(T) is the weighted value of the 
! liquid and ice values.

do k=1,nl   ! was nl-1 until Sept '04
  do mg=1,ifull
    if(cfrac(mg,k)>0.)then
      Tk=tliq(mg,k)+hlcp*(qlg(mg,k)+qfg(mg,k))/cfrac(mg,k) !T in liq cloud
      fl=qlg(mg,k)/max(qfg(mg,k)+qlg(mg,k),1.e-30)
      if(Tk<tfrz.and.qlg(mg,k)>1.e-8)then
        pk=100*prf(mg,k)
        qs=qsati(pk,Tk)
        es=qs*pk/0.622 !ice value
        Aprpr=hl/(rKa*Tk)*(hls/(rvap*Tk)-1.)
        Bprpr=rvap*Tk/((Dva/pk)*es)
        deles=(1.-fice(mg,k))*esdiffx(Tk)
        Cice=1.e3*exp(12.96*deles/es - 0.639) !Meyers et al 1992

        cm0=1.e-12 !Initial crystal mass
        qi0=cm0*Cice/rhoa(mg,k) !Initial ice mixing ratio

! Next 2 lines are for assumption of fully mixed ql and qf (also a line further down).

        qi0=max(qi0, qfg(mg,k)/cfrac(mg,k)) !Assume all qf and ql are mixed
        fd=1.   !Fraction of cloud in which deposition occurs
!       fd=fl   !Or, use option of adjacent ql,qf

        alf=1./3.
        rhoic=700.
        Crate=7.8*((Cice/rhoa(mg,k))**2/rhoic)**(1./3.)*deles/((Aprpr+Bprpr)*es)
        qfdep=fd*cfrac(mg,k)*sqrt(((2./3.)*Crate*tdt+qi0**(2./3.))**3)

! Also need this line for fully-mixed option...
        qfdep = qfdep - qfg(mg,k)

        qfdep=min(qfdep,qlg(mg,k))
        qlg(mg,k)=qlg(mg,k)-qfdep
        qfg(mg,k)=qfg(mg,k)+qfdep
      endif

      fice(mg,k)=qfg(mg,k)/max(qfg(mg,k)+qlg(mg,k),1.e-30)

    endif
  enddo
enddo  

! Calculate new values of vapour mixing ratio and temperature

qtg(:,:)=qtot(:,:)-qcg(:,:)
ttg(:,:)=tliq(:,:)+hlcp*qcg(:,:)+hlfcp*qfg(1:ifull,:)
ccov(:,:)=cfrac(:,:)      !Do this for now

! Vertically sub-grid cloud

do k=2,nl-1
  where (cfrac(1:ifull,k)>1.0e-2.and.cfrac(1:ifull,k+1)==0..and.cfrac(1:ifull,k-1)==0.)
    ccov(1:ifull,k)=sqrt(cfrac(1:ifull,k))
  end where
enddo
     
if(diag.and.mydiag)then
   write(6,*) 'at end of newcloud'
   write(6,*) 'ttg ',ttg(idjd,:)
   write(6,*) 'qcg ',qcg(idjd,:)
   write(6,*) 'qlg ',qlg(idjd,:)
   write(6,*) 'qfg ',qfg(idjd,:)
   write(6,*) 'qtg ',qtg(idjd,:)
endif

return
end subroutine newcloud

! ** Based on GCM revision 1.54 **
! This routine is part of the prognostic cloud scheme. It calculates rainfall
! and the evaporation of rain, and also calls icefall, which does the frozen 
! precipitation. It is called by progcld.
!
! INPUT/OUTPUT
!
! Input:
!
! parameters from include file params.h
!      nl - number of vertical levels
!
! see also include files physparams.h (model physical constants)
!                        cparams.h    (cloud scheme parameters)
!
! from arguments
!      land - logical variable for surface type ( = T for land points)
!      tdt - leapfrog timestep (seconds)
!      fluxc - flux of convective rain in timestep (kg/m**2)
!      rhoa - air density (kg/m**3)
!      dz - layer thicknes (m)
!      ccrain - raining convective cloud fraction at level k
!      prf - pressure at full levels (in hPa. NB: not SI units)
!
! In/Out:
!
! from arguments
!      ttg - temperature (K)
!      qlg - cloud liquid water mixing ratio (kg/kg)
!      qfg - cloud ice mixing ratio (kg/kg)
!      precs - amount of stratiform precipitation in timestep (mm)
!      qtg - water vapour mixing ratio (kg/kg) - called qg in C-CAM
!      cfrac - stratiform cloud fraction
!      ccov - stratiform cloud *cover* looking from above (currently = cfrac)
!
! Output:
!
! from arguments
!      preci - amount of stratiform snowfall in timestep (mm)
!      qevap - evaporation of rainfall (kg/kg)
!      qsubl - sublimation of snowfall (kg/kg)
!      qauto - autoconversion of cloud liquid water (kg/kg)
!      qcoll - collection by rain of cloud liquid water (kg/kg)
!      qaccr - accretion by snow of cloud liquid water (kg/kg)
!
!**************************************************************************

subroutine newrain(land,tdt,fluxc,rhoa,dz,ccrain,prf,cdrop,cfa,qca,ttg,qlg,qfg,qrg,precs,qtg,cfrac,cffall,ccov, &
                   preci,qevap,qsubl,qauto,qcoll,qaccr,fluxr,fluxi,fluxm,pfstayice,pfstayliq,pqfsed,slopes,     &
                   prscav)

use cc_mpi, only : mydiag
use estab, only : esdiffx, qsati, pow75
use kuocomb_m
use morepbl_m  !condx        

implicit none

! Global parameters
include 'newmpar.h'
include 'const_phys.h' !Input physical constants
include 'cparams.h'    !Input cloud scheme parameters
include 'kuocom.h'     !acon,bcon,Rcm,ktsav,nevapls
include 'params.h'     !Input model grid dimensions (modified PARAMS.f for CCAM)
include 'parm.h'

! Argument list
logical land(ifull)
real tdt
real fluxc(ifull,nl)
real rhoa(ifull,nl)
real dz(ifull,nl)
real ccrain(ifull,nl)
real prf(ifull,nl)
real cdrop(ifull,nl)
real ttg(ifull+iextra,nl)
real qlg(ifull+iextra,nl)
real qfg(ifull+iextra,nl)
real qrg(ifull+iextra,nl)
real precs(ifull)
real qtg(ifull+iextra,nl)
real cfrac(ifull,nl)
real cffall(ifull+iextra,nl)
real ccov(ifull,nl)
real preci(ifull)
real qevap(ifull,nl)
real qsubl(ifull,nl)
real qauto(ifull,nl)
real qcoll(ifull,nl)
real qaccr(ifull,nl)
real cfa(ifull,nl)
real qca(ifull,nl)
real pqfsed(ifull,nl)
real pfstayice(ifull,nl)
real pfstayliq(ifull,nl)
real slopes(ifull,nl)
real prscav(ifull,nl)
real fluxr(ifull,nl)
real fluxi(ifull,nl)
real fluxm(ifull,nl)

! Local work arrays and variables
real clfr(ifull,nl)
real qsg(ifull,nl)
real cifr(ifull,nl)
real frclr(ifull,nl)
real qsl(ifull,nl)
real fluxa(ifull,nl)
real clfra(ifull)
real ccra(ifull)
real cfrain(ifull,nl)
real cfmelt(ifull,nl)
real fluxrain(ifull)
real fracr(ifull,nl)
real mxclfr(ifull)
real rdclfr(ifull)
real vr(ifull,nl)
real rhor(ifull,nl-1)
real fout(ifull,nl-1)
real fthru(ifull,nl-1)

integer k,mb,mg,ns,nt

real apr
real bpr
real bl
real cdt
real cev
real clrevap
real coll
real crate
real delt
real dql
real dqsdt
real es
real evap
real fcol
real fr
real frb
real frc
real pk
real qcic
real qcrit
real ql
real ql1
real ql2
real qpf
real R6c,R3c,beta6,eps
real rhodz
real satevap
real selfcoll
real tk
real Wliq
real cfla,dqla,qla
real alph
real rhorin,rhorout
real cffluxin,cffluxout
real mxovr,rdovr


! Start code : ----------------------------------------------------------


delt=tdt

do k=1,nl
  do mg=1,ifull
    fracr(mg,k)=0.
    fluxr(mg,k)=0.
    frclr(mg,k)=0.
    fluxm(mg,k)=0.
    fluxi(mg,k)=0.
    fluxa(mg,k)=0.
    qevap(mg,k)=0.
    qauto(mg,k)=0.
    qcoll(mg,k)=0.
    cfrain(mg,k)=0.
    pk=100.0*prf(mg,k)
    qsg(mg,k)=qsati(pk,ttg(mg,k))
    cifr(mg,k)=cfrac(mg,k)*qfg(mg,k)/max(qlg(mg,k)+qfg(mg,k),1.E-30)
    clfr(mg,k)=cfrac(mg,k)*qlg(mg,k)/max(qlg(mg,k)+qfg(mg,k),1.E-30)
  enddo
enddo

!**************** Cut here to insert new auto scheme ********************            
if (ncloud>0.and.ncloud<3) then
! Using new (subgrid) autoconv scheme... 
  do k=nl-1,1,-1
    do mg=1,ifull
      cfrain(mg,k)=0.0
      rhodz=rhoa(mg,k)*dz(mg,k)

      if(clfr(mg,k)>0.)then

        ql=qlg(mg,k)
        cfla=0.
        dqla=0.
        if(cfa(mg,k)>0.)then
          cfla=cfa(mg,k)*clfr(mg,k)/(clfr(mg,k)+cifr(mg,k))
          qla=qca(mg,k)/cfa(mg,k)

! Following few lines are for Yangang Liu's new scheme (2004: JAS, GRL)

          Wliq = max(1.e-10, 1000. * qla * rhoa(mg,k)) !g/m3
          R6c = 4.09e-4 * ( 1.15e23*1.e-6*cdrop(mg,k) / Wliq**2 )**(1./6.)
          eps = 1. - 0.7 * exp(-0.003e-6*cdrop(mg,k)) !mid range
          beta6 = ((1.+3.*eps**2)*(1.+4.*eps**2)*(1.+5.*eps**2) / ((1.+eps**2)*(1.+2.*eps**2)) )**(1./6.)
          R3c = 1.e-6*R6c/beta6 !in metres
          qcrit=(4.*pi/3.)*rhow*R3c**3*Cdrop(mg,k)/rhoa(mg,k) !New qcrit

          if(qla<=qcrit)then
            ql2=qla
          else

! Following is Liu & Daum (JAS, 2004)
            Crate=1.9e17*(0.75*rhoa(mg,k)/(pi*rhow))**2*beta6**6/cdrop(mg,k)
            ql1=qla/sqrt(1.+2.*crate*qla**2*delt)

            ql1=max(ql1, qcrit) !Intermediate qlg after auto
            Frb=dz(mg,k)*rhoa(mg,k)*(qla-ql1)/delt
            cdt=delt*0.5*Ecol*0.24*pow75(Frb)
            selfcoll=min(ql1,ql1*cdt)
            ql2=ql1-selfcoll
            cfrain(mg,k)=cfla
          endif
          dqla=cfla*(qla-ql2)
          ql=max(1.e-10,qlg(mg,k)-dqla)
        endif
        dql=qlg(mg,k)-ql

        qauto(mg,k)=qauto(mg,k)+dql
        qlg(mg,k)=qlg(mg,k)-dql
        fluxa(mg,k)=dql*rhodz
      endif
    enddo
  enddo

! Or, using old autoconv scheme... also used by prognostic cloud scheme
else
  do k=nl-1,1,-1
    do mg=1,ifull
      cfrain(mg,k)=0.0
      rhodz=rhoa(mg,k)*dz(mg,k)

      if(clfr(mg,k)>0.)then
      
        qcrit=(4.*pi/3.)*rhow*Rcm**3*cdrop(mg,k)/rhoa(mg,k)
        qcic=qlg(mg,k)/clfr(mg,k) !In cloud value

        if(qcic<qcrit)then
          ql=qlg(mg,k)
        else
          Crate=Aurate*rhoa(mg,k)*(rhoa(mg,k)/(cdrop(mg,k)*rhow))**(1./3.)
          ql1=1./pow75(qcic**(-4./3.)+(4./3.)*Crate*delt)
          ql1=max(ql1, qcrit) !Intermediate qlg after auto
          Frb=dz(mg,k)*rhoa(mg,k)*(qcic-ql1)/delt
          cdt=delt*0.5*Ecol*0.24*pow75(Frb) !old
          selfcoll=min(ql1,ql1*cdt)
          ql2=ql1-selfcoll
          ql=clfr(mg,k)*ql2
          cfrain(mg,k)=clfr(mg,k)
        endif
        dql=qlg(mg,k)-ql

        qauto(mg,k)=qauto(mg,k)+dql
        qlg(mg,k)=qlg(mg,k)-dql
        fluxa(mg,k)=dql*rhodz
      endif
    enddo
  enddo
endif ! (ncloud>0)

! Call frozen precipitation routine

call icefall(tdt,rhoa,dz,prf,ttg,qsg,qlg,qfg,qtg,cfrac,cfmelt,fluxi,fluxm,clfr,cifr,qsubl,qaccr,pfstayice,pqfsed,slopes)

! Set up prognostic rain - MJT
! The following has been modified according to LDR's flux divergence calculation
! (see icefall.f).  LDR's original scheme can be recovered by setting fout=1 and
! fthru=1 or using ncloud<=1.

do k=nl-1,1,-1
  do mg=1,ifull
    ! combine autoconversion and prognostic rain
    rhodz=rhoa(mg,k)*dz(mg,k)
    qauto(mg,k)=qauto(mg,k)+qrg(mg,k)
    rhor(mg,k)=qrg(mg,k)*rhoa(mg,k)
    cfrain(mg,k)=max(cfrain(mg,k),cffall(mg,k)) ! max overlap autoconversion and rain from previous time step
  end do
end do
      
if (ncloud<=1) then
  fout=1.
  fthru=1.
end if

! The following has been modified to track the random overlap rain fraction (rdclfr)
! and the max overlap rain fraction (mxclfr) so than both random overlaped and
! maximum/random overlaped clouds are supported - MJT

clfra(1:ifull)=1.e-6
ccra(1:ifull)=0.
fluxrain(1:ifull)=0.
prscav(1:ifull,nl)=0.
pfstayliq(1:ifull,1)=0.
mxclfr(1:ifull)=0. ! max overlap rain fraction
rdclfr(1:ifull)=0. ! rnd overlap rain fraction
        
! Now work down through the levels...
        
do k=nl-1,1,-1
  do mg=1,ifull
    rhodz=rhoa(mg,k)*dz(mg,k)
    evap=0.

    ! The following flag detects maximum/random overlap clouds
    ! that are separated by a clear layer
    if ((clfr(mg,k)<1.e-10.and.cfrain(mg,k)<1.e-10).or.nmr==0) then
      ! combine max overlap from above cloud with net random overlap
      rdclfr(mg)=rdclfr(mg)+mxclfr(mg)-rdclfr(mg)*mxclfr(mg)
      mxclfr(mg)=0.
    end if

! Add flux of melted snow to fluxrain

    fluxrain(mg)=fluxrain(mg)+fluxm(mg,k)+fluxa(mg,k)
            
! Evaporation of rain

    qpf=fluxrain(mg)/rhodz !Mix ratio of rain which falls into layer  ! MJT suggestion
    clrevap=(1.-clfr(mg,k))*qpf                                       ! MJT suggestion
    if(fluxrain(mg)>0.)then
      pk=100.0*prf(mg,k)
      qsg(mg,k)=qsati(pk,ttg(mg,k))
      if(ttg(mg,k)<tfrz.and.ttg(mg,k)>=tice)then
        qsl(mg,k)=qsg(mg,k)+epsil*esdiffx(ttg(mg,k))/(100.0*prf(mg,k))            ! MJT suggestion
      else
        qsl(mg,k)=qsg(mg,k)
      endif             !qsl is qs value over liquid surface
      Tk=ttg(mg,k)
      es=qsl(mg,k)*pk/epsil 
      Apr=(hl/(rKa*Tk))*(hl/(rvap*Tk)-1.)
      Bpr=rvap*Tk/((Dva/pk)*es)
      Fr=fluxrain(mg)/delt/clfra(mg)
      Cev=clfra(mg)*3.8e2*sqrt(Fr/rhoa(mg,k))/(qsl(mg,k)*(Apr+Bpr))
      dqsdt=hl*qsl(mg,k)/(rvap*ttg(mg,k)**2)
      bl=1.+0.5*Cev*delt*(1.+hlcp*dqsdt)
      evap=delt*(Cev/bl)*(qsl(mg,k)-qtg(mg,k))
      satevap=(qsl(mg,k)-qtg(mg,k))/(1.+hlcp*dqsdt) !Evap to saturate
!     Vr=11.3*Fr**(1./9.)/sqrt(rhoa(mg,k))  !Actual fall speed
!     Vr=5./sqrt(rhoa(mg,k))                !Nominal fall speed

      evap=max(0., min(evap,satevap,qpf,clrevap))
      if(nevapls==-1)evap=0.
      if(nevapls==-2.and.k<=ktsav(mg).and.condx(mg)>0.)evap=0.
      if(nevapls==-3)evap=.5*evap
      if(nevapls==-4.and.k<=ktsav(mg).and.condx(mg)>0.)evap=.5*evap ! usual
      qevap(mg,k)=qevap(mg,k)+evap
      qtg(mg,k)=qtg(mg,k)+evap
      ttg(mg,k)=ttg(mg,k)-hlcp*evap
    endif
    frclr(mg,k)=rhodz*(clrevap-evap) !over delt ! MJT suggestion
            

! Now do the collection term.

    if(fluxrain(mg)>0.)then
      Fr=fluxrain(mg)/clfra(mg)/delt
      mxovr=min(mxclfr(mg),clfr(mg,k))                ! max overlap
      mxovr=max(cfrain(mg,k),mxovr)
      rdovr=rdclfr(mg)*clfr(mg,k)                     ! rnd overlap
      cfrain(mg,k)=mxovr+rdovr-mxovr*rdovr            ! combine collection
    else
      Fr=0.
    endif

    if(fluxc(mg,k+1)>0.)then
      Frc=max(0.,fluxc(mg,k+1)/max(ccra(mg),0.01)/tdt) ! over tdt
      rdovr=clfr(mg,k)*ccra(mg)                                ! rnd overlap
      cfrain(mg,k)=cfrain(mg,k)+rdovr-cfrain(mg,k)*rdovr  
      !cfrain(mg,k)=max(cfrain(mg,k),min(clfr(mg,k),ccra(mg))) ! max overlap
    else
      Frc=0.
    endif

! The collection term comprises collection by stratiform rain falling from
! above (Fr), stratiform rain released in this grid box (Frb), and
! convective rain (Frc).
! Frb term now done above.

    fcol=min(1.,mxclfr(mg)/(1.e-20+clfr(mg,k))) !max overlap
    fcol=fcol+rdclfr(mg)-fcol*rdclfr(mg)        !rnd overlap
    cdt=delt*Ecol*0.24*(fcol*pow75(Fr)+ccra(mg)*pow75(Frc))
    prscav(mg,k)=delt*0.24*fcol*pow75(Fr) !Strat only

    coll=min(qlg(mg,k),qlg(mg,k)*cdt/(1.+0.5*cdt))
    qcoll(mg,k)=qcoll(mg,k)+coll
    qlg(mg,k)=qlg(mg,k)-coll
    fluxrain(mg)=fluxrain(mg)+coll*rhodz

! subtract evaporated rain

    fluxrain(mg)=fluxrain(mg)-rhodz*evap
    fluxrain(mg)=max(fluxrain(mg),0.) !To avoid roundoff -ve's

! Calculate the raining cloud cover down to this level, for stratiform (clfra)
! and convective (ccra).

    cfrain(mg,k)=min(1.,cfrain(mg,k)+cfmelt(mg,k)-cfrain(mg,k)*cfmelt(mg,k)) ! rnd overlap
    if (frclr(mg,k)<1.e-15) then
      rdclfr(mg)=0.
      mxclfr(mg)=0.
    end if
    mxclfr(mg)=max(mxclfr(mg),cfrain(mg,k)) !max overlap
    clfra(mg)=max(1.e-15,rdclfr(mg)+mxclfr(mg)-rdclfr(mg)*mxclfr(mg)) !rnd overlap the mx and rd rain fractions
    ccra(mg)=max(ccra(mg),ccrain(mg,k)) !always max overlap for convective rainfall - MJT
    fracr(mg,k)=clfra(mg)

    ! Calculate rain fall speed (MJT)
    if (ncloud>1) then
      Fr=fluxrain(mg)/delt/clfra(mg)
      vr(mg,k)=11.3*Fr**(1./9.)/sqrt(rhoa(mg,k))  !Actual fall speed
      vr(mg,k)=max(vr(mg,k),0.1)
      alph=delt*vr(mg,k)/dz(mg,k)
      fout(mg,k)=1.-exp(-alph)       !analytical
      fthru(mg,k)=1.-fout(mg,k)/alph !analytical
    end if
            
    pfstayliq(mg,nlp-k)=fluxrain(mg)*(1.-fthru(mg,k))/tdt

! Compute fluxes into the box
    cffluxin=clfra(mg)-cfrain(mg,k)
    rhorin=fluxrain(mg)/dz(mg,k)

! Compute the fluxes of rain leaving the box
! Use the flux-divergent form as in Rotstayn (QJRMS, 1997)
    cffluxout=cfrain(mg,k)*fout(mg,k)
    rhorout=rhor(mg,k)*fout(mg,k)
            
! Update the rhor and cffall fields
    cffall(mg,k)=cfrain(mg,k)-cffluxout+cffluxin*(1.-fthru(mg,k))
    rhor(mg,k)=rhor(mg,k)-rhorout+rhorin*(1.-fthru(mg,k))
    fluxrain(mg)=rhorout*dz(mg,k)+fluxrain(mg)*fthru(mg,k) 
! Now fluxrain is flux leaving layer k
    fluxr(mg,k)=fluxr(mg,k)+fluxrain(mg)
            
  enddo
enddo

! Re-create qrg field

qrg(1:ifull,nl)=0.
do k=1,nl-1
  qrg(1:ifull,k)=rhor(1:ifull,k)/rhoa(1:ifull,k)
enddo

! Factor 0.5 here accounts for leapfrog scheme

precs(1:ifull)=precs(1:ifull)+0.5*(fluxr(1:ifull,1)+fluxi(1:ifull,1))
preci(1:ifull)=preci(1:ifull)+0.5*fluxi(1:ifull,1)

! Remove small amounts of cloud

do k=1,nl
  where (qlg(1:ifull,k)<1e-10.or.clfr(1:ifull,k)<1e-5)
    qtg(1:ifull,k)=qtg(1:ifull,k)+qlg(1:ifull,k)
    ttg(1:ifull,k)=ttg(1:ifull,k)-hlcp*qlg(1:ifull,k)
    qlg(1:ifull,k)=0.
    clfr(1:ifull,k)=0.
  end where
  where (qfg(1:ifull,k)<1e-10.or.cifr(1:ifull,k)<1e-5)
    qtg(1:ifull,k)=qtg(1:ifull,k)+qfg(1:ifull,k)
    ttg(1:ifull,k)=ttg(1:ifull,k)-hlscp*qfg(1:ifull,k)
    qfg(1:ifull,k)=0.
    cifr(1:ifull,k)=0.
  end where
  where (qrg(1:ifull,k)<1e-10.or.cffall(1:ifull,k)<1e-5)
    qtg(1:ifull,k)=qtg(1:ifull,k)+qrg(1:ifull,k)
    ttg(1:ifull,k)=ttg(1:ifull,k)-hlcp*qrg(1:ifull,k)
    qrg(1:ifull,k)=0.
    cffall(1:ifull,k)=0.
  end where
enddo

!      Adjust cloud fraction (and cloud cover) after precipitation
if(nmaxpr==1.and.mydiag)then
  write(6,*) 'diags from newrain for idjd ',idjd
  write (6,"('cfrac ',9f8.3/6x,9f8.3)") cfrac(idjd,:)
  write (6,"('cftemp',9f8.3/6x,9f8.3)") cifr(idjd,:)+clfr(idjd,:)
  write (6,"('ccov_in',9f8.3/6x,9f8.3)") ccov(idjd,:)
endif

return
end subroutine newrain


! This routine is part of the prognostic cloud scheme. It calculates the frozen 
! precipitation. It is called by newrain.
!
! INPUT/OUTPUT
!
! Input:
!
! parameters from include file params.h
!      nl - number of vertical levels
!
! see also include files physparams.h (model physical constants)
!                        cparams.h    (cloud scheme parameters)
!
! from arguments
!
!      lg - latitude index (ranges from 1 at poles to LAT at equator)
!      tdt - leapfrog timestep (seconds)
!      rhoa - air density (kg/m**3)
!      dz - layer thicknes (m)
!      prf - pressure at full levels (in hPa. NB: not SI units)
!
! In/Out:
!
! from arguments
!
!      ttg - temperature (K)
!      qsg - saturation mixing ratio (kg/kg)
!      qlg - cloud liquid water mixing ratio (kg/kg)
!      qfg - cloud ice mixing ratio (kg/kg)
!      qtg - water vapour mixing ratio (kg/kg) - called qg in C-CAM
!      cfrac - stratiform cloud fraction
!      cfmelt - fraction of grid box occupied by falling ice that melts
!
! Output:
!
! from arguments
!
!      fluxi - flux of falling ice in timestep (kg/m**2)
!      fluxm - flux of falling ice that melts in timestep (kg/m**2)
!      clfr - liquid cloud fraction
!      qsubl - sublimation of snowfall (kg/kg)
!      qaccr - accretion by snow of cloud liquid water (kg/kg)
!      pqfsed - (dqf/qf) due to ice falling out of layer (fraction)
!      pfstay - incoming flux of ice staying in layer  (kg/m**2/s)
!
!******************************************************************************

subroutine icefall(tdt,rhoa,dz,prf,ttg,qsg,qlg,qfg,qtg,cfrac,cfmelt,fluxi,fluxm,clfr,cifr,qsubl,qaccr,pfstay,pqfsed,slopes)

use cc_mpi, only : mydiag
use kuocomb_m
use morepbl_m  !condx  

implicit none

! Global parameters
include 'newmpar.h'
include 'const_phys.h' !Input physical constants
include 'cparams.h'    !Input cloud scheme parameters
include 'kuocom.h'     !ktsav,ldr,nevapls
include 'parm.h'       !just for nmaxpr
include 'params.h'     !Input model grid dimensions (modified PARAMS.f for CCAM)

! Argument list
real tdt
real rhoa(ifull,nl)
real dz(ifull,nl)
real prf(ifull,nl)
real ttg(ifull+iextra,nl)
real qsg(ifull,nl)
real qlg(ifull+iextra,nl)
real qfg(ifull+iextra,nl)
real qtg(ifull+iextra,nl)
real cfrac(ifull,nl)
real cfmelt(ifull,nl)
real fluxi(ifull,nl)
real fluxm(ifull,nl)
real clfr(ifull,nl)
real cifr(ifull,nl)
real qsubl(ifull,nl)
real qaccr(ifull,nl)
real pqfsed(ifull,nl)
real pfstay(ifull,nl)

! Local work arrays and variables
real Csbsav(ifull,nl-1)
real cifra(ifull)
real fluxice(ifull)
real fout(ifull,nl-1)
real fthru(ifull,nl-1)
real gam(ifull,nl-1)
real qfdiv(ifull,nl)
real rhoi(ifull,nl)
real rica(ifull)
real slopes(ifull,nl)
real vi2(ifull,nl)
real mxclfr(ifull)
real rdclfr(ifull)

integer k,mg

real alph
real alphaf
real aprpr
real bf
real bprpr
real caccr
real cdt
real cffluxin
real cffluxout
real csb
real curly
real delt
real dqf
real dqs
real dttg
real es
real fsclr
real pk
real qif
real ql
real rhoiin
real rhoiout
real sublflux
real tc
real tk
real viin
 
! Start code : ----------------------------------------------------------

! Set up timestep for ice processes

      delt=tdt

! Convert from mixing ratio to density of ice, and work out ice cloud fraction

do k=1,nl
  fluxi(1:ifull,k)=0.
  ! N.B. qfg >= 0 from dynamics, but occasionally it seems to become -ve,
  ! apparently from newrain.f, or else lower down in this subroutine
  qfg(1:ifull,k)=max( 0.,qfg(1:ifull,k) )  
  rhoi(1:ifull,k)=rhoa(1:ifull,k)*qfg(1:ifull,k) 
  cifr(1:ifull,k)=cfrac(1:ifull,k)*qfg(1:ifull,k)/max(qlg(1:ifull,k)+qfg(1:ifull,k),1.e-30)
  clfr(1:ifull,k)=max(cfrac(1:ifull,k)-cifr(1:ifull,k),0.)
  qsubl(1:ifull,k)=0.
  qaccr(1:ifull,k)=0.
  qfdiv(1:ifull,k)=0.
  cfmelt(1:ifull,k)=0.
end do

if(diag.and.mydiag)then  ! JLM
  write(6,*) 'entering icefall.'
  write(6,*) 'cfrac ',cfrac(idjd,:)
  write(6,*) 'cifr ',cifr(idjd,:)
  write(6,*) 'clfr ',clfr(idjd,:)
  write(6,*) 'qlg ',qlg(idjd,:)
  write(6,*) 'qfg ',qfg(idjd,:)
endif  ! (diag.and.mydiag)

! Set up ice fall speed field and other arrays

vi2(1:ifull,nl)=0.1
slopes(1:ifull,nl)=0.

select case(abs(ldr))
  case(1)  ! 1 for R21 runs, like prev lw=22
    do k=nl-1,1,-1
      do mg=1,ifull
        if(cifr(mg,k)>=1.e-10)then
          vi2(mg,k)=3.23*(rhoi(mg,k)/cifr(mg,k))**0.17
        else
          vi2(mg,k)=vi2(mg,k+1)
        endif
      enddo
    enddo

  case(2)
    do k=nl-1,1,-1
      do mg=1,ifull
        vi2(mg,k)=vi2(mg,k+1)
        if(cifr(mg,k)>=1.e-10)then
          vi2(mg,k)=0.9*3.23*(rhoi(mg,k)/cifr(mg,k))**0.17
        endif
      enddo
    enddo

  case(3)
    do k=nl-1,1,-1
      do mg=1,ifull
        vi2(mg,k)=vi2(mg,k+1)
        if(cifr(mg,k)>=1.e-10)then
          vi2(mg,k)=max(0.1,2.05+0.35*log10(qfg(mg,k)/cifr(mg,k)))
        endif
      enddo
    enddo

  case(4)
    do k=nl-1,1,-1
      do mg=1,ifull
        vi2(mg,k)=vi2(mg,k+1)
        if(cifr(mg,k)>=1.e-10)then
          vi2(mg,k)=1.4*3.23*(rhoi(mg,k)/cifr(mg,k))**0.17
        endif
      enddo
    enddo

!     following are alternative slightly-different versions of above
!     used for I runs from 29/4/05 till 30/8/05
!     for given qfg, large cifr implies small ice crystals, 
!     with a small fall speed. 
!     Note that for very small qfg, cifr is small.
!     But rhoi is like qfg, so ratio should also be small and OK.
  case(11) ! 1 for R21 runs, like prev lw=22
    do k=nl-1,1,-1
      do mg=1,ifull
        vi2(mg,k)=max( vi2(mg,k+1),3.23*(rhoi(mg,k)/max(cifr(mg,k),1.e-30))**0.17 )
      enddo
    enddo

  case(22)
    do k=nl-1,1,-1
      do mg=1,ifull
        vi2(mg,k)=max( vi2(mg,k+1),.9*3.23*(rhoi(mg,k)/max(cifr(mg,k),1.e-30))**0.17 )
      enddo
    enddo

  case(33)
    do k=nl-1,1,-1
      do mg=1,ifull
!          following max gives vi2=.1 for qfg=cifr=0
        vi2(mg,k)=max( vi2(mg,k+1),2.05 +0.35*log10(max(qfg(mg,k),2.68e-36)/max(cifr(mg,k),1.e-30)) )
      enddo
    enddo
end select

do k=nl-1,1,-1
  do mg=1,ifull
    tc=ttg(mg,k)-tfrz
    slopes(mg,k)=1.6e3*10**(-0.023*tc)
    alphaf=hls*qsg(mg,k)/(rvap*ttg(mg,k)**2)
    gam(mg,k)=hlscp*alphaf !(L/cp)*dqsdt (HBG notation)

! Set up the Rate constant for snow sublimation
          
    Tk=ttg(mg,k)
    pk=100.*prf(mg,k)
    es=qsg(mg,k)*pk/epsil
    Aprpr=(hls/(rKa*Tk))*(hls/(rvap*Tk)-1.)
    Bprpr=rvap*Tk/((Dva/pk)*es)
    curly=0.65*slopes(mg,k)**2+0.493*slopes(mg,k)*sqrt(slopes(mg,k)*vi2(mg,k+1)*rhoa(mg,k)/um)
    if(nevapls==-1)curly=0.
    if(nevapls==-2.and.condx(mg)>0..and.k<=ktsav(mg))curly=0.

! Define the rate constant for sublimation of snow, omitting factor rhoi

    Csbsav(mg,k)=4*curly/(rhoa(mg,k)*qsg(mg,k)*(Aprpr+Bprpr)*pi*vi2(mg,k+1)*rhosno)

! Set up the parameters for the flux-divergence calculation

    alph=delt*vi2(mg,k)/dz(mg,k)
    fout(mg,k)=1.-exp(-alph) !analytical
    fthru(mg,k)=1.-fout(mg,k)/alph !analytical
  enddo
enddo

! Save sedimentation rate for aerosol scheme

do k=1,nl-1
  pqfsed(1:ifull,k)=fout(1:ifull,k)
enddo
pqfsed(1:ifull,nl)=0.

! The following has been modified to track the random overlap rain fraction (rdclfr)
! and the max overlap rain fraction (mxclfr) so than both random overlaped and
! max/random overlaped clouds are supported - MJT

! Assume no cloud at top level

fluxice(1:ifull)=0.
cifra(1:ifull)=0.
rica(1:ifull)=0.
pfstay(1:ifull,nl)=0.
mxclfr(1:ifull)=0. ! max overlap rain fraction
rdclfr(1:ifull)=0. ! rnd overlap rain fraction

! Now work down through the levels...

do k=nl-1,1,-1
  do mg=1,ifull

    sublflux=0.
    fsclr=0.
    caccr=0.
    dqf=0.

    ! The following flag detects max/random overlap clouds
    ! that are separated by a clear layer
    if (cifr(mg,k)<1.e-10.or.nmr==0) then
      ! combine max overlap from last cloud with net random overlap
      rdclfr(mg)=rdclfr(mg)+mxclfr(mg)-rdclfr(mg)*mxclfr(mg)
      mxclfr(mg)=0.
    end if

! Melt falling ice if > 0 deg C

    if(ttg(mg,k)>tfrz.and.fluxice(mg)>0.)then
      qif=fluxice(mg)/(dz(mg,k)*rhoa(mg,k))      !Mixing ratio of ice
      dttg=-hlfcp*qif
      ttg(mg,k)=ttg(mg,k)+dttg
      qsg(mg,k)=qsg(mg,k)+gam(mg,k)*dttg/hlscp
      fluxm(mg,k)=fluxm(mg,k)+fluxice(mg)
      cfmelt(mg,k)=cifra(mg)
      fluxice(mg)=0.
      cifra(mg)=0.
      rdclfr(mg)=0.
      mxclfr(mg)=0.
    endif


! Compute the sublimation of ice falling from level k+1 into level k

    fsclr=(1.-cifr(mg,k)-clfr(mg,k))*fluxice(mg)
    if(fluxice(mg)>0..and.qtg(mg,k)<qsg(mg,k))then ! sublime snow
      Csb=Csbsav(mg,k)*fluxice(mg)/delt !LDR
      bf=1.+0.5*Csb*delt*(1.+gam(mg,k))
      dqs=max(0.,delt*(Csb/bf)*(qsg(mg,k)-qtg(mg,k)))
      dqs=min(dqs,(qsg(mg,k)-qtg(mg,k))/(1.+gam(mg,k))) !Don't supersat.
      sublflux=min(dqs*rhoa(mg,k)*dz(mg,k),fsclr)
      fluxice(mg)=fluxice(mg)-sublflux
      fsclr=fsclr-sublflux
      dqs=sublflux/(rhoa(mg,k)*dz(mg,k))
      qsubl(mg,k)=qsubl(mg,k)+dqs
      qtg(mg,k)=qtg(mg,k)+dqs
      dttg=-hlscp*dqs
      ttg(mg,k)=ttg(mg,k)+dttg
      qsg(mg,k)=qsg(mg,k)+gam(mg,k)*dttg/hlscp
    endif

! Save this for the wet deposition scheme.

    pfstay(mg,k)=fluxice(mg)*(1.-fthru(mg,k))/tdt !Flux staying in layer k

! Accretion of cloud water by falling ice
! This calculation uses the incoming fluxice without subtracting sublimation
! (since subl occurs only outside cloud), so add sublflux back to fluxice.
            
    if(fluxice(mg)+sublflux>0.and.qlg(mg,k)>1.e-10)then
      ql=qlg(mg,k)
      cdt=Eac*slopes(mg,k)*(fluxice(mg)+sublflux)/(2.*rhosno)
      dqf=min(ql,cifra(mg)*ql,ql*cdt/(1.+0.5*cdt))

      clfr(mg,k)=clfr(mg,k)*(1.-dqf/qlg(mg,k))
      caccr=clfr(mg,k)*dqf/qlg(mg,k)
      qlg(mg,k)=qlg(mg,k)-dqf
      qaccr(mg,k)=qaccr(mg,k)+dqf
      fluxice(mg)=fluxice(mg)+rhoa(mg,k)*dz(mg,k)*dqf
      dttg=hlfcp*dqf
      ttg(mg,k)=ttg(mg,k)+dttg
      qsg(mg,k)=qsg(mg,k)+gam(mg,k)*dttg/hlscp
    endif


    if(fsclr<1.e-15)then
      rdclfr(mg)=0.
      mxclfr(mg)=0.
    end if
    mxclfr(mg)=max(mxclfr(mg),cifr(mg,k)+caccr) !max overlap
    cifra(mg)=max(0.01,mxclfr(mg)+rdclfr(mg)-mxclfr(mg)*rdclfr(mg)) !rnd overlap the mx and rd ice fractions

! Compute fluxes into the box
            
    if(fluxice(mg)>0.)then
      rhoiin=fluxice(mg)/dz(mg,k)
      cffluxin=min(1.,rhoiin/rica(mg))
    else
      rhoiin=0.
      cffluxin=0.
    endif


! Compute the fluxes of ice and cloud amount leaving the box
            
    if(cifr(mg,k)>=1.e-10)then

! Use the flux-divergent form as in Rotstayn (QJRMS, 1997)
      rhoiout=rhoi(mg,k)*fout(mg,k)
      cffluxout=cifr(mg,k)*fout(mg,k)

      rica(mg)=rhoi(mg,k)/cifr(mg,k) !in cloud rhoi above
    else !Keep value of rica from above
      rhoiout=0.
      cffluxout=0.
    endif
            
! Update the rhoi and cifr fields

    cifr(mg,k)=min(1.-clfr(mg,k),(cifr(mg,k)-cffluxout)+cffluxin*(1.-fthru(mg,k)))
    rhoi(mg,k)=(rhoi(mg,k)-rhoiout)+rhoiin*(1.-fthru(mg,k))
    fluxice(mg)=rhoiout*dz(mg,k)+fluxice(mg)*fthru(mg,k) 
! Now fluxice is flux leaving layer k
    fluxi(mg,k)=fluxi(mg,k)+fluxice(mg)

  enddo
enddo

! End of small timestep loop

! Re-create qfg field

do k=1,nl
  qfg(1:ifull,k)=rhoi(1:ifull,k)/rhoa(1:ifull,k)
enddo

! Diagnostics for debugging
if(diag.and.mydiag)then  ! JLM
  write(6,*) 'near end of icefall.'
  write(6,*) 'vi2',vi2(idjd,:)
  write(6,*) 'cfraci ',cfrac(idjd,:)
  write(6,*) 'cifr',cifr(idjd,:)
  write(6,*) 'clfr',clfr(idjd,:)
  write(6,*) 'ttg',ttg(idjd,:)
  write(6,*) 'qsg',qsg(idjd,:)         
  write(6,*) 'qlg',qlg(idjd,:)
  write(6,*) 'qfg',qfg(idjd,:)
  write(6,*) 'qsubl',qsubl(idjd,:)
  write(6,*) 'rhoa',rhoa(idjd,:)
  write(6,*) 'rhoi',rhoi(idjd,:)
  write(6,*) 'fluxi ',fluxi(idjd,:)
  write(6,*) 'gam',gam(idjd,:)
  write(6,*) 'fout',fout(idjd,:)
  write(6,*) 'fthru',fthru(idjd,:)
  write(6,*) 'pqfsed',pqfsed(idjd,:)
  write(6,*) 'cfmelt',cfmelt(idjd,:)
  write(6,*) 'fluxm',fluxm(idjd,:)
  write(6,*) 'cifra,fluxice,rica',cifra(idjd),fluxice(idjd),rica(idjd)
endif  ! (diag.and.mydiag)

return
end subroutine icefall

end module leoncld_mod