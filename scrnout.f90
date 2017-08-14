! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2016 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
! This file is part of the Conformal Cubic Atmospheric Model (CCAM)
!
! CCAM is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! CCAM is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with CCAM.  If not, see <http://www.gnu.org/licenses/>.

!------------------------------------------------------------------------------
    
subroutine scrnout(zo,ustar,factch,wetfac,qsttg,qgscrn,tscrn,uscrn,u10,rhscrn,af,aft,ri,vmod, &
                   bprm,cms,chs,chnsea,nalpha)

use arrays_m
use cc_mpi, only : mydiag,myid
use const_phys
use diag_m
use estab
use liqwpar_m  ! ifullw,qfg,qlg  just for diags
use morepbl_m  ! condx,fg,eg
use newmpar_m
use nsibd_m    ! rsmin,ivegt,sigmf,ssdn,rmc
use parm_m
use pbl_m
use permsurf_m
use prec_m     ! just for diags
use sigs_m
use soil_m
use soilsnow_m

implicit none

! from Feb '05 scales uscrn, u10 to remove vmod>2 effect
! allows for zo>zscr from 30/7/04
! has fixer at bottom to ensure tscrn and qgscrn bounded (& non-neg)
integer, parameter :: nits=2  ! nits=2 for 2 iterations
integer, parameter :: ntest=0 ! ntest= 0 for diags off; ntest>= 1 for diags 
integer iq,nalpha,numq,numt,numu
real, parameter :: vkar=.4, zscr=1.8
real alf,chnscr,chnsea,fact,qtgnet
real rho,srcp,vfact,zscronzt,zscrr,ztv,z10
real aft2ice,aft10ice,aft38ice,zt,zscrt,z10t
real bprm,cms,chs,denha,denma,ri2x,root,vmag,zlog,es
real deltaq,deltat
real, dimension(ifull) :: ri,vmod,af,aft
real, dimension(ifull) :: qgscrn,tscrn,uscrn,u10,rhscrn
real, dimension(ifull) :: wetfac,factch,qsttg,ustar,zo,theta
real, dimension(ifull) :: qsurf,tsurf,af2,af10,afroot2,afroot10
real, dimension(ifull) :: aft2,aft10,rich2,rich10,qstarx,tstarx
real, dimension(ifull) :: fh2,fh10,fh38,fm2,fm10,fm38

srcp =sig(1)**(rdry/cp)
ztv=exp(vkar/sqrt(chn10))/10.  ! proper inverse of ztsea
zscronzt=zscr*ztv
chnscr=(vkar/log(zscronzt))**2

if(diag.or.ntest>0)then
  call maxmin(qg,'qg',ktau,1.e3,kl)
  call maxmin(qfg,'qf',ktau,1.e3,kl)
  call maxmin(qlg,'ql',ktau,1.e3,kl)
  call maxmin(t,' t',ktau,1.,kl)
  call maxmin(wb,'wb',ktau,1.,ms)
  call maxmin(tggsn,'tggsn',ktau,1.,3)
  call maxmin(tgg,'tgg',ktau,1.,ms)
  call maxmin(tss,'ts',ktau,1.,1)
  call maxmin(precip,'pr',ktau,1.,1)
  call maxmin(eg,'eg',ktau,1.,1)
  call maxmin(fg,'fg',ktau,1.,1)
  call maxmin(ps,'ps',ktau,.01,1)
  if(mydiag)then
    write(6,*) 'entering scrnout; wetfac,tss,qsttg: ',wetfac(idjd),tss(idjd),qsttg(idjd)
  endif
endif
         
! calculate "surface" specific humidity
! N.B. qsttg is coming in as before-calling-sib3 value of qs_tss
if(nalpha==1) then  ! usual
  do iq=1,ifull
    qsurf(iq) = wetfac(iq)*qsttg(iq) + (1.-wetfac(iq))*min(qg(iq,1),qsttg(iq))  ! for v. cold surfaces
  enddo   ! iq=1,ifull
else
  do iq=1,ifull
    qtgnet=qsttg(iq)*wetfac(iq) -qg(iq,1)
    if(qtgnet>0.) then
      qsurf(iq) = qsttg(iq)*wetfac(iq)
    else
      qsurf(iq) = .1*qsttg(iq)*wetfac(iq) + .9*qg(iq,1)
    endif
  enddo   ! iq loop
endif   ! (nalpha==1) then .. else ..

do iq=1,ifull
  zscrr=max(zscr,zo(iq)+1.)
  z10=max(10.,zo(iq)+2.)
  zlog=log(zo(iq))
  afroot2(iq)=vkar/(log(zscrr)-zlog)
  af2(iq)=afroot2(iq)*afroot2(iq)
  ! constants for calculating 2m, 10m winds
  afroot10(iq) = vkar/(log(z10)-zlog)
  af10(iq) = afroot10(iq)*afroot10(iq)
  zt=zo(iq)/factch(iq)**2 ! MJT urban
  zscrt=max(zscr,zt+1.)     ! Aug '05
  z10t=max(z10,zt+2.)       ! Aug '05
  aft2(iq)=vkar*afroot2(iq)/log(zscrt/zt)     ! land only
  aft10(iq)=vkar*afroot10(iq)/log(z10t/zt)    ! land only
  tsurf(iq)=tss(iq)

  aft2ice=(vkar/log(zscrr/.001))**2
  aft10ice=(vkar/log(z10/.001))**2
  aft38ice=(vkar/log(zmin/.001))**2 ! needs recalc, as aft is ice only
  if(.not.land(iq))then 
    tsurf(iq)=fracice(iq)*tggsn(iq,1)+(1.-fracice(iq))*tpan(iq) ! Dec05 ! MJT seaice
    aft2(iq)=fracice(iq)*aft2ice+(1.-fracice(iq))*chnscr
    aft10(iq)=fracice(iq)*aft10ice+(1.-fracice(iq))*chn10
    aft(iq)=fracice(iq)*aft38ice+(1.-fracice(iq))*chnsea
  endif   ! (.not.land(iq)) 
enddo

! the following assumes a suitable ri(level 1) has been found
do iq=1,ifull              ! all valid points 
  zscrr=max(zscr,zo(iq)+1.)
  z10=max(10.,zo(iq)+2.)
  theta(iq)=t(iq,1)/srcp   ! can use same value from sflux
  deltat=tsurf(iq)-theta(iq)
  deltaq=qsurf(iq)-qg(iq,1)
  rho=ps(iq)/(rdry*tsurf(iq))
  ! section for solving for screen height values of various parameters
  if (ri(iq) > 0.) then
    ! stable case     
    fm38(iq)=1./(1.+bprm*ri(iq))**2
    fh38(iq)=fm38(iq)
    ! rich2(iq)=ri(iq)*zscrr/zmin        ! usually a good first guess
    fact=af2(iq)/(af(iq)*fm38(iq))
    alf=ri(iq)*zscrr*fact*sqrt(fact)*aft(iq)*fh38(iq)*sqrt(fh38(iq))/(zmin*aft2(iq))
    rich2(iq)=( -1. + sqrt(1.+4.*bprm*alf) ) /(2.*bprm)
    fm2(iq)=1./(1.+bprm*rich2(iq))**2
    fh2(iq)=fm2(iq)
    fact=af10(iq)/(af(iq)*fm38(iq))
    alf=ri(iq)*z10*fact*sqrt(fact)*aft(iq)*fh38(iq)*sqrt(fh38(iq))/(zmin*aft10(iq))
    rich10(iq)=( -1. + sqrt(1.+4.*bprm*alf) ) /(2.*bprm)
    fm10(iq)=1./(1.+bprm*rich10(iq))**2
    fh10(iq)=fm10(iq)
    ! alternative formulae - not checked out	  
  else
    ! unstable case
    ! Starting value of the Richardson number.
    root=sqrt(-ri(iq)*zmin/zo(iq))        
    denma=1.+cms*2.*bprm*af(iq)*root                         
    fm38(iq)=1.-2.*bprm *ri(iq)/denma                           
    denha=1.+chs*2.*bprm*factch(iq)*aft(iq)*root             
    fh38(iq)=1.-2.*bprm *ri(iq)/denha                        
    rich2(iq)=ri(iq)*zscrr/zmin        ! usually a good first guess
    root=sqrt(-rich2(iq)*zscrr/zo(iq))        
    denma=1.+cms*2.*bprm*af2(iq)*root                         
    fm2(iq)=1.-2.*bprm *rich2(iq)/denma                           
    denha=1.+chs*2.*bprm*factch(iq)*aft2(iq)*root             
    fh2(iq)=1.-2.*bprm *rich2(iq)/denha          
    rich10(iq)=ri(iq)*z10/zmin        ! usually a good first guess
    root=sqrt(-rich10(iq)*z10/zo(iq))        
    denma=1.+cms*2.*bprm*af10(iq)*root      
    fm10(iq)=1.-2.*bprm *rich10(iq)/denma                           
    denha=1.+chs*2.*bprm*factch(iq)*aft10(iq)*root             
    fh10(iq)=1.-2.*bprm *rich10(iq)/denha                        
    if(nits==2)then
      fact=af2(iq)*fm2(iq)/(af(iq)*fm38(iq))
      rich2(iq)=ri(iq)*zscrr*fact*sqrt(fact)*aft(iq)*fh38(iq)/(zmin*aft2(iq)*fh2(iq))
      root=sqrt(-rich2(iq)*zscrr/zo(iq))        
      denma=1.+cms*2.*bprm*af2(iq)*root                         
      fm2(iq)=1.-2.*bprm *rich2(iq)/denma                           
      denha=1.+chs*2.*bprm*factch(iq)*aft2(iq)*root             
      fh2(iq)=1.-2.*bprm *rich2(iq)/denha          
      fact=af10(iq)*fm10(iq)/(af(iq)*fm38(iq))
      rich10(iq)=ri(iq)*z10*fact*sqrt(fact)*aft(iq)*fh38(iq)/(zmin*aft10(iq)*fh10(iq))
      root=sqrt(-rich10(iq)*z10/zo(iq))        
      denma=1.+cms*2.*bprm*af10(iq)*root      
      fm10(iq)=1.-2.*bprm *rich10(iq)/denma                           
      denha=1.+chs*2.*bprm*factch(iq)*aft10(iq)*root             
      fh10(iq)=1.-2.*bprm *rich10(iq)/denha                        
    endif  ! (nits==2)
  endif

  ! in calculating tstar & qstar, use derived ri38, NOT eg and fg	
  tstarx(iq)=aft(iq)*fh38(iq)*deltat/sqrt(af(iq)*fm38(iq)) ! MJT cable
  qstarx(iq)=aft(iq)*fh38(iq)*deltaq/sqrt(af(iq)*fm38(iq)) ! MJT cable
  fact=sqrt(af2(iq)*fm2(iq))/(aft2(iq)*fh2(iq))
  tscrn(iq)=tsurf(iq)-fact*tstarx(iq)
       
  ! N.B. over unstable sea points, may sometimes get supersat qgscrn
  ! also over stable snow points, as tscrn & qgscrn derived independently in
  ! location with large vertical gradients
  qgscrn(iq)=qsurf(iq)-fact*qstarx(iq)
     
  ! screen wind speeds
  vfact=sqrt(u(iq,1)**2+v(iq,1)**2)/vmod(iq)
  uscrn(iq) = vfact*vmod(iq)*sqrt(af(iq)*fm38(iq))/(afroot2(iq)*sqrt(fm2(iq))) ! MJT cable
  u10(iq) = vfact*vmod(iq)*sqrt(af(iq)*fm38(iq))/(afroot10(iq)*sqrt(fm10(iq))) ! MJT cable
  if(ntest==-1)then
    vmag=sqrt(u(iq,1)**2+v(iq,1)**2)
    if(u10(iq)>vmag)then
      write(6,*) 'strange iq,vmod,vmag,u10/vmod,u10/vmag ',iq,vmod(iq),vmag,u10(iq)/vmod(iq),u10(iq)/vmag
      write(6,*) 'ri,rich10,ustar,fm10 ',ri(iq),rich10(iq),ustar(iq),fm10(iq)
    endif
  endif
  ! use tscrn for calc rhscrn  23/12/05	
  es=establ(tscrn(iq))
  qsttg(iq)= .622*es/(ps(iq)-es)  ! recalc     
  rhscrn(iq) = 100.*qgscrn(iq)/qsttg(iq)
enddo

if(diag.and.mydiag)then
  iq=idjd
  fact=sqrt(af2(iq)*fm2(iq))/(aft2(iq)*fh2(iq))
  write(6,*) 'in scrnout qsurf,qstarx,qgscrn,qg1 ',qsurf(iq),qstarx(iq),qgscrn(iq),qg(iq,1)
  write(6,*) 'tsurf,tstarx,tscrn,t1 ',tsurf(iq),tstarx(iq),tscrn(iq),t(iq,1)
  write(6,*) 'tpan,tgg2,ri,fact ',tpan(iq),tgg(iq,2),ri(iq),fact
  write(6,*) 'in scrnout; epot,eg,fg: ',epot(idjd),eg(idjd),fg(idjd)
  write(6,*) 'aft,vmod,fh38,ustar ',aft(iq),vmod(iq),fh38(iq),ustar(iq)
  write(6,*) 'af2,fh2,fm2 ',af2(iq),fh2(iq),fm2(iq) 
  write(6,*) 'new qsttg: ',qsttg(iq)
  write(6,*) 'uscrn,u10,rhscrn ',uscrn(iq),u10(iq),rhscrn(iq)
endif

if(ntest>=1.and.myid==0)then
  numt=0
  numq=0
  numu=0
  do iq=1,ifull
    zscrr=max(zscr,zo(iq)+1.)
    z10=max(10.,zo(iq)+2.)
    vmag=sqrt(u(iq,1)**2+v(iq,1)**2)
    if(uscrn(iq)>vmag.or.u10(iq)>vmag)then
      numu=numu+1
      write(6,"('strange u iq,land,uscrn,u10,vmag,vmod,tstarx ',i5,l2,6f7.3)") &
           iq,land(iq),uscrn(iq),u10(iq),vmag,vmod(iq),tstarx(iq)
      write(6,"('more u isoil,iveg,snowd,fg',2i3,6f8.2)") isoilm(iq),ivegt(iq),snowd(iq),fg(iq)
      write(6,*) 'tgg1,tggsn1,tss,theta ',tgg(iq,1),tggsn(iq,1),tss(iq),theta(iq)
    endif
    if(tscrn(iq)<min(theta(iq),tsurf(iq)).or.tscrn(iq)>max(theta(iq),tsurf(iq)).or.iq==idjd )then
      numt=numt+1
      write(6,"('checking T iq,land,sicedep,tsurf,tscrn,theta,tstarx ',i5,l2,4f7.2,f7.3)") iq,land(iq),sicedep(iq),  &
             tsurf(iq),tscrn(iq),theta(iq),tstarx(iq)
      write(6,"('tss,tgg1,tpan,tgg3,t1,fracice ',6f7.2)") tss(iq),tgg(iq,1),tpan,tgg(iq,3),t(iq,1),fracice(iq)
      write(6,*) 'zo ',zo(iq)
      write(6,*) 'ri2,ri10,ri38 ',rich2(iq),rich10(iq),ri(iq)
      write(6,*) 'cm2,ch2,cm10,ch10 ',cms*2.*bprm*af2(iq)*sqrt(zscrr/zo(iq)),                                   &
                     cms*2.*bprm*af10(iq)*sqrt(z10/zo(iq)),chs*2.*bprm*aft2(iq)*factch(iq)*sqrt(zscrr/zo(iq)),  &
                     chs*2.*bprm*aft10(iq)*factch(iq)*sqrt(z10/zo(iq))  
      write(6,*) 'af2,af10,af38 ',af2(iq),af10(iq),af(iq)
      write(6,*) 'aft2,aft10,aft38 ',aft2(iq),aft10(iq),aft(iq)
      write(6,*) 'rat2,rat10,rat38 ',sqrt(af2(iq))/aft2(iq),sqrt(af10(iq))/aft10(iq),sqrt(af(iq))/aft(iq)    
      write(6,*) 'fh2,fh10,fh38 ',fh2(iq),fh10(iq),fh38(iq)
      write(6,*) 'fm2,fm10,fm38 ',fm2(iq),fm10(iq),fm38(iq)
      write(6,*) 'vmod,fact2,fact10,fact38 ',vmod(iq),sqrt(af2(iq)*fm2(iq))/(aft2(iq)*fh2(iq)), &
                    sqrt(af10(iq)*fm10(iq))/(aft10(iq)*fh10(iq)),sqrt(af(iq)*fm38(iq))/(aft(iq)*fh38(iq))
      if(sicedep(iq)>0.)write(6,*) 'tpan,tgg3,fracice ',tpan(iq),tgg(iq,3),fracice(iq)       
    endif
    if(qgscrn(iq)<min(qg(iq,1),qsurf(iq)).or.qgscrn(iq)>max(qg(iq,1),qsurf(iq)))then
      numq=numq+1
      write(6,*) 'strange qg iq,land,sicedep,qstarx,t,ri ',iq,land(iq),sicedep(iq),qstarx(iq),t(iq,1),ri(iq),numq
      write(6,*) 'qsurf,qgscrn,qg1 ',qsurf(iq),qgscrn(iq),qg(iq,1)
      write(6,*) 'tsurf,tscrn,theta ',tsurf(iq),tscrn(iq),theta(iq)
    endif
  enddo
  write(6,*) 'numq,numt,numu',numq,numt,numu
  if(ntest==2)then
    do iq=1,ifull,50
      ! prints first guess and final values
      ri2x=rich2(iq)
      fact=af2(iq)*fm2(iq)/(af(iq)*fm38(iq))
      if(ri(iq)<0.)ri2x=ri(iq)*zscrr*fact*sqrt(fact)*aft(iq)*fh38(iq)/(zmin*aft2(iq)*fh2(iq))
      write (6,"('iq,ri2a,ri2,ri2x,ri10a,ri10 ',i5,5f9.5)") iq,zscrr*ri(iq)/zmin,rich2(iq),ri2x,z10*ri(iq)/zmin,rich10(iq)
    enddo
  endif
endif   ! (ntest>=1)      

! consistency fixer for tscrn and qgscrn (odd pts)  March 2004
! not needed from March 2005 with new version of scrnout
! N.B. T & qg problems are for land and sea-ice points

tscrn(:)=tscrn(:)-.018  ! apply 1.8 m approx. adiab. corrn.
if(ntest==-1.and.myid==0)then
  do iq=1,ifull
    if(.not.land(iq).and.sicedep(iq)<1.e-20)then 
      write(23,'(f6.2,2f9.6,f6.2,2i6,3f7.4)') u10(iq),cduv(iq)/vmod(iq),zo(iq),vmod(iq),ktau,iq,ri(iq),fh10(iq),ustar(iq)
    endif
  enddo        
endif

return
end subroutine scrnout
      
subroutine screencalc(ifull,qscrn,rhscrn,tscrn,uscrn,u10,ustar,tstar,qstar,thetavstar,zo,zoh,zoq,stemp, &
                      temp,smixr,mixr,umag,ps,zmin,sig)
 
use const_phys
use estab
use parm_m

implicit none

integer, intent(in) :: ifull
integer ic
real, dimension(ifull), intent(out) :: qscrn,rhscrn,tscrn,uscrn,u10
real, dimension(ifull), intent(out) :: ustar,tstar,qstar,thetavstar
real, dimension(ifull), intent(in) :: zo,zoh,zoq,stemp,temp,umag
real, dimension(ifull), intent(in) :: smixr,mixr,ps,zmin
real, dimension(ifull) :: lzom,lzoh,lzoq,thetav,sthetav
real, dimension(ifull) :: z_on_l,z0_on_l,zt_on_l,z10_on_l
real, dimension(ifull) :: pm0,ph0,pm1,ph1,integralm,integralh
real, dimension(ifull) :: neutral,neutral10,pm10
real, dimension(ifull) :: integralm10,zq_on_l,integralq
real, dimension(ifull) :: esatb,qsatb,umagn
real, dimension(ifull) :: pq0,pq1
real, intent(in) :: sig
real scrp
integer, parameter ::  nc     = 5
real, parameter    ::  vkar   = 0.4
real, parameter    ::  a_1    = 1.
real, parameter    ::  b_1    = 2./3.
real, parameter    ::  c_1    = 5.
real, parameter    ::  d_1    = 0.35
real, parameter    ::  z0     = 1.5
real, parameter    ::  z10    = 10.

scrp             = (sig)**(rdry/cp)
thetav(1:ifull)  = temp(1:ifull)*(1.+0.61*mixr(1:ifull))/scrp
sthetav(1:ifull) = stemp(1:ifull)*(1.+0.61*smixr(1:ifull))
umagn(1:ifull)   = max(umag(1:ifull),vmodmin)

! Roughness length for heat
lzom(1:ifull) = log(zmin(1:ifull)/zo(1:ifull))
lzoh(1:ifull) = log(zmin(1:ifull)/zoh(1:ifull))
lzoq(1:ifull) = log(zmin(1:ifull)/zoq(1:ifull))

! Dyer and Hicks approach 
thetavstar(1:ifull) = vkar*(thetav(1:ifull)-sthetav(1:ifull))/lzoh(1:ifull)
ustar(1:ifull)      = vkar*umagn(1:ifull)/lzom(1:ifull)
do ic = 1,nc
  z_on_l(1:ifull)  = vkar*zmin(1:ifull)*grav*thetavstar(1:ifull)/(thetav(1:ifull)*ustar(1:ifull)**2)
  z_on_l(1:ifull)  = min(z_on_l(1:ifull),10.)
  z0_on_l(1:ifull) = z_on_l(1:ifull)*zo(1:ifull)/zmin(1:ifull)
  zt_on_l(1:ifull) = z_on_l(1:ifull)*zoh(1:ifull)/zmin(1:ifull)
  zq_on_l(1:ifull) = z_on_l(1:ifull)*zoq(1:ifull)/zmin(1:ifull)
  where ( z_on_l(1:ifull)<0. )
    pm0(1:ifull) = (1.-16.*z0_on_l(1:ifull))**(-0.25)
    ph0(1:ifull) = (1.-16.*zt_on_l(1:ifull))**(-0.5)
    pq0(1:ifull) = (1.-16.*zq_on_l(1:ifull))**(-0.5)
    pm1(1:ifull) = (1.-16.*z_on_l(1:ifull))**(-0.25)
    ph1(1:ifull) = (1.-16.*z_on_l(1:ifull))**(-0.5)
    pq1(1:ifull) = ph1(1:ifull)
    integralm(1:ifull) = lzom(1:ifull)-2.*log((1.+1./pm1(1:ifull))/(1.+1./pm0(1:ifull))) &
                      -log((1.+1./pm1(1:ifull)**2)/(1.+1./pm0(1:ifull)**2))              &
                      +2.*(atan(1./pm1(1:ifull))-atan(1./pm0(1:ifull)))
    integralh(1:ifull) = lzoh(1:ifull)-2.*log((1.+1./ph1(1:ifull))/(1.+1./ph0(1:ifull)))
    integralq(1:ifull) = lzoq(1:ifull)-2.*log((1.+1./pq1(1:ifull))/(1.+1./pq0(1:ifull)))
  elsewhere
    !--------------Beljaars and Holtslag (1991) momentum & heat            
    pm0(1:ifull) = -(a_1*z0_on_l(1:ifull)+b_1*(z0_on_l(1:ifull)-(c_1/d_1))*exp(-d_1*z0_on_l(1:ifull)) &
                   +b_1*c_1/d_1)
    pm1(1:ifull) = -(a_1*z_on_l(1:ifull)+b_1*(z_on_l(1:ifull)-(c_1/d_1))*exp(-d_1*z_on_l(1:ifull)) &
                   +b_1*c_1/d_1)
    ph0(1:ifull) = -((1.+(2./3.)*a_1*zt_on_l(1:ifull))**1.5+b_1*(zt_on_l(1:ifull)-(c_1/d_1))*exp(-d_1*zt_on_l(1:ifull)) &
               +b_1*c_1/d_1-1.)
    ph1(1:ifull) = -((1.+(2./3.)*a_1*z_on_l(1:ifull))**1.5+b_1*(z_on_l(1:ifull)-(c_1/d_1))*exp(-d_1*z_on_l(1:ifull))    &
                +b_1*c_1/d_1-1.)
    pq0(1:ifull) = -((1.+(2./3.)*a_1*zq_on_l(1:ifull))**1.5+b_1*(zq_on_l(1:ifull)-(c_1/d_1))*exp(-d_1*zq_on_l(1:ifull)) &
                +b_1*c_1/d_1-1.)
    pq1(1:ifull) = ph1(1:ifull)
    integralm(1:ifull) = lzom(1:ifull)-(pm1(1:ifull)-pm0(1:ifull))    
    integralh(1:ifull) = lzoh(1:ifull)-(ph1(1:ifull)-ph0(1:ifull))
    integralq(1:ifull) = lzoq(1:ifull)-(pq1(1:ifull)-pq0(1:ifull))
  endwhere
  integralm(1:ifull)  = max(integralm(1:ifull),1.e-10)
  integralh(1:ifull)  = max(integralh(1:ifull),1.e-10)
  integralq(1:ifull)  = max(integralq(1:ifull),1.e-10)
  thetavstar(1:ifull) = vkar*(thetav(1:ifull)-sthetav(1:ifull))/integralh(1:ifull)
  ustar(1:ifull)      = vkar*umagn(1:ifull)/integralm(1:ifull)
end do
tstar(1:ifull) = vkar*(temp(1:ifull)-stemp(1:ifull))/integralh(1:ifull)
qstar(1:ifull) = vkar*(mixr(1:ifull)-smixr(1:ifull))/integralq(1:ifull)
      
! estimate screen diagnostics
z_on_l(1:ifull)    = vkar*zmin(1:ifull)*grav*thetavstar(1:ifull)/(thetav(1:ifull)*ustar(1:ifull)**2)
z_on_l(1:ifull)    = min(z_on_l(1:ifull),10.)
z0_on_l(1:ifull)   = z0*z_on_l(1:ifull)/zmin(1:ifull)
z10_on_l(1:ifull)  = z10*z_on_l(1:ifull)/zmin(1:ifull)
z0_on_l(1:ifull)   = min(z0_on_l(1:ifull),10.)
z10_on_l(1:ifull)  = min(z10_on_l(1:ifull),10.)
neutral(1:ifull)   = log(zmin(1:ifull)/z0)
neutral10(1:ifull) = log(zmin(1:ifull)/z10)
where ( z_on_l(1:ifull)<0. )
  ph0(1:ifull)  = (1.-16.*z0_on_l(1:ifull))**(-0.50)
  ph1(1:ifull)  = (1.-16.*z_on_l(1:ifull))**(-0.50)
  pm0(1:ifull)  = (1.-16.*z0_on_l(1:ifull))**(-0.25)
  pm1(1:ifull)  = (1.-16.*z_on_l(1:ifull))**(-0.25)
  pm10(1:ifull) = (1.-16.*z10_on_l(1:ifull))**(-0.25)  
  integralh(1:ifull)   = neutral(1:ifull)-2.*log((1.+1./ph1(1:ifull))/(1.+1./ph0(1:ifull)))
  integralq(1:ifull)   = integralh(1:ifull)
  integralm(1:ifull)   = neutral(1:ifull)-2.*log((1.+1./pm1(1:ifull))/(1.+1./pm0(1:ifull))) &
                     -log((1.+1./pm1(1:ifull)**2)/(1.+1./pm0(1:ifull)**2))                  &
                     +2.*(atan(1./pm1(1:ifull))-atan(1./pm0(1:ifull)))
  integralm10(1:ifull) = neutral10(1:ifull)-2.*log((1.+1./pm1(1:ifull))/(1.+1./pm10(1:ifull))) &
                     -log((1.+1./pm1(1:ifull)**2)/(1.+1./pm10(1:ifull)**2))                    &
                     +2.*(atan(1./pm1(1:ifull))-atan(1./pm10(1:ifull)))     
elsewhere
!-------Beljaars and Holtslag (1991) heat function
  ph0(1:ifull)  = -((1.+(2./3.)*a_1*z0_on_l(1:ifull))**1.5+b_1*(z0_on_l(1:ifull)-(c_1/d_1)) &
              *exp(-d_1*z0_on_l(1:ifull))+b_1*c_1/d_1-1.)
  ph1(1:ifull)  = -((1.+(2./3.)*a_1*z_on_l(1:ifull))**1.5+b_1*(z_on_l(1:ifull)-(c_1/d_1))   &
              *exp(-d_1*z_on_l(1:ifull))+b_1*c_1/d_1-1.)
  pm0(1:ifull)  = -(a_1*z0_on_l(1:ifull)+b_1*(z0_on_l(1:ifull)-(c_1/d_1))*exp(-d_1*z0_on_l(1:ifull))+b_1*c_1/d_1)
  pm10(1:ifull) = -(a_1*z10_on_l(1:ifull)+b_1*(z10_on_l(1:ifull)-(c_1/d_1))*exp(-d_1*z10_on_l(1:ifull))+b_1*c_1/d_1)
  pm1(1:ifull)  = -(a_1*z_on_l(1:ifull)+b_1*(z_on_l(1:ifull)-(c_1/d_1))*exp(-d_1*z_on_l(1:ifull))+b_1*c_1/d_1)
  integralh(1:ifull)   = neutral(1:ifull)-(ph1(1:ifull)-ph0(1:ifull))
  integralq(1:ifull)   = integralh(1:ifull)
  integralm(1:ifull)   = neutral(1:ifull)-(pm1(1:ifull)-pm0(1:ifull))
  integralm10(1:ifull) = neutral10(1:ifull)-(pm1(1:ifull)-pm10(1:ifull))
endwhere
integralh(1:ifull)  = max(integralh(1:ifull), 1.e-10)
integralq(1:ifull)  = max(integralq(1:ifull), 1.e-10)
integralm(1:ifull)  = max(integralm(1:ifull), 1.e-10)
integralm10(1:ifull)  = max(integralm10(1:ifull), 1.e-10)

tscrn(1:ifull)  = temp(1:ifull) - tstar(1:ifull)*integralh(1:ifull)/vkar
esatb(1:ifull)  = establ(tscrn(1:ifull))
qscrn(1:ifull)  = mixr(1:ifull) - qstar(1:ifull)*integralq(1:ifull)/vkar
qscrn(1:ifull)  = max(qscrn(1:ifull), qgmin)
qsatb(1:ifull)  = 0.622*esatb(1:ifull)/(ps(1:ifull)-esatb(1:ifull))
rhscrn(1:ifull) = 100.*min(qscrn(1:ifull)/qsatb(1:ifull), 1.)
uscrn(1:ifull)  = max(umagn(1:ifull)-ustar(1:ifull)*integralm(1:ifull)/vkar, 0.)
u10(1:ifull)    = max(umagn(1:ifull)-ustar(1:ifull)*integralm10(1:ifull)/vkar, 0.)
      
return
end subroutine screencalc
      
subroutine autoscrn
      
use arrays_m
use const_phys
use estab
use extraout_m
use liqwpar_m
use mlo
use morepbl_m, only : urban_tas, urban_zom, urban_zoh, urban_zoq, &
                      urban_ts, urban_wetfac
use newmpar_m
use nharrs_m
use nsibd_m, only : sigmu
use parm_m
use pbl_m
use permsurf_m
use screen_m
use sigs_m
use soil_m
use soilsnow_m
use work2_m
      
implicit none
      
integer iq
real, dimension(ifull) :: umag, zminx, smixr
real, dimension(ifull) :: ou, ov, atu, atv, iu, iv
real, dimension(ifull) :: au, av, es, rho
real, dimension(ifull) :: u_qgscrn, u_rhscrn, u_tscrn, u_uscrn, u_u10
real, dimension(ifull) :: u_ustar, u_tstar, u_qstar, u_thetavstar
real, dimension(ifull) :: u_zo, u_zoh, u_zoq, u_tss, u_smixr
      
zminx(:) = (bet(1)*t(1:ifull,1)+phi_nh(:,1))/grav
zminx(:) = max( zminx(:), 5. )
if ( nmlo/=0 ) then
  iu(:) = 0.
  iv(:) = 0.
  ou(:) = 0.
  ov(:) = 0.
  call mloexport(2,ou,1,0)
  call mloexport(3,ov,1,0)
  call mloexpice(iu,9,0)
  call mloexpice(iv,10,0)
  ou(:) = (1.-fracice(:))*ou(:) + fracice(:)*iu(:)
  ov(:) = (1.-fracice(:))*ov(:) + fracice(:)*iv(:)
else
  ou(:) = 0.
  ov(:) = 0.
end if
au(:)   = u(1:ifull,1) - ou(:)
av(:)   = v(1:ifull,1) - ov(:)
umag(:) = max( sqrt(au(:)*au(:)+av(:)*av(:)), vmodmin )
es(1:ifull) = establ(tss(1:ifull))
qsttg(1:ifull) = 0.622*es(:)/(ps(1:ifull)-es(:))
smixr(:) = wetfac(:)*qsttg(:) + (1.-wetfac(:))*min( qsttg(:), qg(1:ifull,1) )
      
call screencalc(ifull,qgscrn,rhscrn,tscrn,uscrn,u10,ustar,tstar,qstar,thetavstar, &
                zo,zoh,zoq,tss,t(1:ifull,1),smixr,qg(1:ifull,1),umag,ps(1:ifull), &
                zminx,sig(1))

rho(:) = ps(1:ifull)/(rdry*tss(:))
cduv(:) = ustar(:)**2/umag(:)
taux(:) = rho(:)*cduv(:)*au(:)
tauy(:) = rho(:)*cduv(:)*av(:)

atu(:)   = au(:)*uscrn(:)/umag(:) + ou(:)
atv(:)   = av(:)*uscrn(:)/umag(:) + ov(:)
uscrn(:) = sqrt(atu(:)*atu(:)+atv(:)*atv(:))
atu(:)   = au(:)*u10(:)/umag(:) + ou(:)
atv(:)   = av(:)*u10(:)/umag(:) + ov(:)      
u10(:)   = sqrt(atu(:)*atu(:)+atv(:)*atv(:))


! urban tile
where ( sigmu>0. )
  u_zo    = urban_zom
  u_zoh   = urban_zoh
  u_zoq   = urban_zoq
  u_tss   = urban_ts
  u_smixr = urban_wetfac*qsttg + (1.-urban_wetfac)*min( qsttg, qg(1:ifull,1) )
elsewhere
  u_zo    = zo
  u_zoh   = zoh
  u_zoq   = zoq
  u_tss   = tss
  u_smixr = smixr
end where
call screencalc(ifull,u_qgscrn,u_rhscrn,u_tscrn,u_uscrn,u_u10,u_ustar,u_tstar,u_qstar,u_thetavstar, &
                u_zo,u_zoh,u_zoq,u_tss,t(1:ifull,1),u_smixr,qg(1:ifull,1),umag,ps(1:ifull),zminx,   &
                sig(1))
urban_tas = u_tscrn

return
end subroutine autoscrn
