! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2024 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
subroutine scrnout(zo,ustar,zoh,wetfac,qsttg,qgscrn,tscrn,uscrn,u10,rhscrn,af,aft,ri,vmod, &
                   bprm,cms,chs,chnsea,is,ie)

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
integer, intent(in) :: is, ie
integer iq,numq,numt,numu
real, parameter :: vkar=.4, zscr=1.8
real alf,chnscr,chnsea,fact,qtgnet
real rho,srcp,vfact,zscronzt,zscrr,ztv,z10
real aft2ice,aft10ice,aft38ice,zt,zscrt,z10t
real bprm,cms,chs,denha,denma,ri2x,root,vmag,zlog,es
real deltaq,deltat
real, dimension(is:ie) :: ri,vmod,af,aft
real, dimension(is:ie) :: qgscrn,tscrn,uscrn,u10,rhscrn
real, dimension(is:ie) :: wetfac,zoh,qsttg,ustar,zo,theta
real, dimension(is:ie) :: qsurf,tsurf,af2,af10,afroot2,afroot10
real, dimension(is:ie) :: aft2,aft10,rich2,rich10,qstarx,tstarx
real, dimension(is:ie) :: fh2,fh10,fh38,fm2,fm10,fm38
real, dimension(is:ie) :: factch

factch(is:ie) = sqrt(zo(is:ie)/zoh(is:ie))
srcp =sig(1)**(rdry/cp)
ztv=exp(vkar/sqrt(chn10))/10.  ! proper inverse of ztsea
zscronzt=zscr*ztv
chnscr=(vkar/log(zscronzt))**2

if(diag.or.ntest>0.and.ntiles==1)then
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
  do iq=is,ie
    qsurf(iq) = wetfac(iq)*qsttg(iq) + (1.-wetfac(iq))*min(qg(iq,1),qsttg(iq))  ! for v. cold surfaces
  enddo   ! iq=is,ie
else
  do iq=is,ie
    qtgnet=qsttg(iq)*wetfac(iq) -qg(iq,1)
    if(qtgnet>0.) then
      qsurf(iq) = qsttg(iq)*wetfac(iq)
    else
      qsurf(iq) = .1*qsttg(iq)*wetfac(iq) + .9*qg(iq,1)
    endif
  enddo   ! iq loop
endif   ! (nalpha==1) then .. else ..

do iq=is,ie
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
do iq=is,ie              ! all valid points 
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

if(diag.and.mydiag.and.ntiles==1)then
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

if(ntest>=1.and.myid==0.and.ntiles==1)then
  numt=0
  numq=0
  numu=0
  do iq=is,ie
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
    do iq=is,ie,50
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

tscrn(is:ie)=tscrn(is:ie)-.018  ! apply 1.8 m approx. adiab. corrn.
if(ntest==-1.and.myid==0.and.ntiles==1)then
  do iq=is,ie
    if(.not.land(iq).and.sicedep(iq)<1.e-20)then 
      write(23,'(f6.2,2f9.6,f6.2,2i6,3f7.4)') u10(iq),cduv(iq)/vmod(iq),zo(iq),vmod(iq),ktau,iq,ri(iq),fh10(iq),ustar(iq)
    endif
  enddo        
endif

return
end subroutine scrnout
      
subroutine screencalc(imax,qscrn,rhscrn,tscrn,uscrn,u10,ustar,tstar,qstar,thetavstar, &
                      zo,zoh,zoq,stemp,temp,smixr,mixr,umag,ps,zmin,sig)
 
use const_phys
use estab
use parm_m

implicit none

integer, intent(in) :: imax
integer ic
real, dimension(imax), intent(out) :: qscrn,rhscrn,tscrn,uscrn,u10
real, dimension(imax), intent(out) :: ustar,tstar,qstar,thetavstar
real, dimension(imax), intent(in) :: zo,zoh,zoq,stemp,temp,umag
real, dimension(imax), intent(in) :: smixr,mixr,ps,zmin
real, dimension(imax) :: lzom,lzoh,lzoq,thetav,sthetav
real, dimension(imax) :: z_on_l,z0_on_l,zt_on_l,z10_on_l
real, dimension(imax) :: pm0,ph0,pm1,ph1,integralm,integralh
real, dimension(imax) :: neutral,neutral10,pm10
real, dimension(imax) :: integralm10,zq_on_l,integralq
real, dimension(imax) :: esatb,qsatb,umagn
real, dimension(imax) :: pq0,pq1
real, intent(in) :: sig
real scrp
integer, parameter ::  nc     = 5 ! number of iterations
real, parameter    ::  vkar   = 0.4
real, parameter    ::  a_1    = 1.
real, parameter    ::  b_1    = 2./3.
real, parameter    ::  c_1    = 5.
real, parameter    ::  d_1    = 0.35
real, parameter    ::  z0     = 1.5
real, parameter    ::  z10    = 10.

scrp            = (sig)**(rdry/cp)
thetav(1:imax)  = temp(1:imax)*(1.+0.61*mixr(1:imax))/scrp
sthetav(1:imax) = stemp(1:imax)*(1.+0.61*smixr(1:imax))
umagn(1:imax)   = max(umag(1:imax),vmodmin)

! Roughness length for heat
lzom(1:imax) = log(zmin(1:imax)/max(zo(1:imax),1.e-10))
lzoh(1:imax) = log(zmin(1:imax)/max(zoh(1:imax),1.e-10))
lzoq(1:imax) = log(zmin(1:imax)/max(zoq(1:imax),1.e-10))

! Dyer and Hicks approach 
thetavstar(1:imax) = vkar*(thetav(1:imax)-sthetav(1:imax))/lzoh(1:imax)
ustar(1:imax)      = vkar*umagn(1:imax)/lzom(1:imax)
do ic = 1,nc
  z_on_l(1:imax)  = vkar*zmin(1:imax)*grav*thetavstar(1:imax)/(thetav(1:imax)*ustar(1:imax)**2)
  z_on_l(1:imax)  = min(z_on_l(1:imax),10.)
  z0_on_l(1:imax) = z_on_l(1:imax)*zo(1:imax)/zmin(1:imax)
  zt_on_l(1:imax) = z_on_l(1:imax)*zoh(1:imax)/zmin(1:imax)
  zq_on_l(1:imax) = z_on_l(1:imax)*zoq(1:imax)/zmin(1:imax)
  where ( z_on_l(1:imax)<0. )
    pm0(1:imax) = (1.-16.*z0_on_l(1:imax))**(-0.25)
    ph0(1:imax) = (1.-16.*zt_on_l(1:imax))**(-0.5)
    pq0(1:imax) = (1.-16.*zq_on_l(1:imax))**(-0.5)
    pm1(1:imax) = (1.-16.*z_on_l(1:imax))**(-0.25)
    ph1(1:imax) = (1.-16.*z_on_l(1:imax))**(-0.5)
    pq1(1:imax) = ph1(1:imax)
    integralm(1:imax) = lzom(1:imax)-2.*log((1.+1./pm1(1:imax))/(1.+1./pm0(1:imax))) &
                      -log((1.+1./pm1(1:imax)**2)/(1.+1./pm0(1:imax)**2))            &
                      +2.*(atan(1./pm1(1:imax))-atan(1./pm0(1:imax)))
    integralh(1:imax) = lzoh(1:imax)-2.*log((1.+1./ph1(1:imax))/(1.+1./ph0(1:imax)))
    integralq(1:imax) = lzoq(1:imax)-2.*log((1.+1./pq1(1:imax))/(1.+1./pq0(1:imax)))
  elsewhere
    !--------------Beljaars and Holtslag (1991) momentum & heat            
    pm0(1:imax) = -(a_1*z0_on_l(1:imax)+b_1*(z0_on_l(1:imax)-(c_1/d_1))*exp(-d_1*z0_on_l(1:imax)) &
                   +b_1*c_1/d_1)
    pm1(1:imax) = -(a_1*z_on_l(1:imax)+b_1*(z_on_l(1:imax)-(c_1/d_1))*exp(-d_1*z_on_l(1:imax))    &
                   +b_1*c_1/d_1)
    ph0(1:imax) = -((1.+(2./3.)*a_1*zt_on_l(1:imax))**1.5+b_1*(zt_on_l(1:imax)-(c_1/d_1))*exp(-d_1*zt_on_l(1:imax)) &
               +b_1*c_1/d_1-1.)
    ph1(1:imax) = -((1.+(2./3.)*a_1*z_on_l(1:imax))**1.5+b_1*(z_on_l(1:imax)-(c_1/d_1))*exp(-d_1*z_on_l(1:imax))    &
                +b_1*c_1/d_1-1.)
    pq0(1:imax) = -((1.+(2./3.)*a_1*zq_on_l(1:imax))**1.5+b_1*(zq_on_l(1:imax)-(c_1/d_1))*exp(-d_1*zq_on_l(1:imax)) &
                +b_1*c_1/d_1-1.)
    pq1(1:imax) = ph1(1:imax)
    integralm(1:imax) = lzom(1:imax)-(pm1(1:imax)-pm0(1:imax))    
    integralh(1:imax) = lzoh(1:imax)-(ph1(1:imax)-ph0(1:imax))
    integralq(1:imax) = lzoq(1:imax)-(pq1(1:imax)-pq0(1:imax))
  endwhere
  integralm(1:imax)  = max(integralm(1:imax),1.e-10)
  integralh(1:imax)  = max(integralh(1:imax),1.e-10)
  integralq(1:imax)  = max(integralq(1:imax),1.e-10)
  thetavstar(1:imax) = vkar*(thetav(1:imax)-sthetav(1:imax))/integralh(1:imax)
  ustar(1:imax)      = vkar*umagn(1:imax)/integralm(1:imax)
end do
tstar(1:imax) = vkar*(temp(1:imax)-stemp(1:imax))/integralh(1:imax)
qstar(1:imax) = vkar*(mixr(1:imax)-smixr(1:imax))/integralq(1:imax)
      
! estimate screen diagnostics
z_on_l(1:imax)    = vkar*zmin(1:imax)*grav*thetavstar(1:imax)/(thetav(1:imax)*ustar(1:imax)**2)
z_on_l(1:imax)    = min(z_on_l(1:imax),10.)
z0_on_l(1:imax)   = z0*z_on_l(1:imax)/zmin(1:imax)
z10_on_l(1:imax)  = z10*z_on_l(1:imax)/zmin(1:imax)
z0_on_l(1:imax)   = min(z0_on_l(1:imax),10.)
z10_on_l(1:imax)  = min(z10_on_l(1:imax),10.)
neutral(1:imax)   = log(zmin(1:imax)/z0)
neutral10(1:imax) = log(zmin(1:imax)/z10)
where ( z_on_l(1:imax)<0. )
  ph0(1:imax)  = (1.-16.*z0_on_l(1:imax))**(-0.50)
  ph1(1:imax)  = (1.-16.*z_on_l(1:imax))**(-0.50)
  pm0(1:imax)  = (1.-16.*z0_on_l(1:imax))**(-0.25)
  pm1(1:imax)  = (1.-16.*z_on_l(1:imax))**(-0.25)
  pm10(1:imax) = (1.-16.*z10_on_l(1:imax))**(-0.25)  
  integralh(1:imax)   = neutral(1:imax)-2.*log((1.+1./ph1(1:imax))/(1.+1./ph0(1:imax)))
  integralq(1:imax)   = integralh(1:imax)
  integralm(1:imax)   = neutral(1:imax)-2.*log((1.+1./pm1(1:imax))/(1.+1./pm0(1:imax))) &
                     -log((1.+1./pm1(1:imax)**2)/(1.+1./pm0(1:imax)**2))                  &
                     +2.*(atan(1./pm1(1:imax))-atan(1./pm0(1:imax)))
  integralm10(1:imax) = neutral10(1:imax)-2.*log((1.+1./pm1(1:imax))/(1.+1./pm10(1:imax))) &
                     -log((1.+1./pm1(1:imax)**2)/(1.+1./pm10(1:imax)**2))                    &
                     +2.*(atan(1./pm1(1:imax))-atan(1./pm10(1:imax)))     
elsewhere
!-------Beljaars and Holtslag (1991) heat function
  ph0(1:imax)  = -((1.+(2./3.)*a_1*z0_on_l(1:imax))**1.5+b_1*(z0_on_l(1:imax)-(c_1/d_1)) &
              *exp(-d_1*z0_on_l(1:imax))+b_1*c_1/d_1-1.)
  ph1(1:imax)  = -((1.+(2./3.)*a_1*z_on_l(1:imax))**1.5+b_1*(z_on_l(1:imax)-(c_1/d_1))   &
              *exp(-d_1*z_on_l(1:imax))+b_1*c_1/d_1-1.)
  pm0(1:imax)  = -(a_1*z0_on_l(1:imax)+b_1*(z0_on_l(1:imax)-(c_1/d_1))*exp(-d_1*z0_on_l(1:imax))+b_1*c_1/d_1)
  pm10(1:imax) = -(a_1*z10_on_l(1:imax)+b_1*(z10_on_l(1:imax)-(c_1/d_1))*exp(-d_1*z10_on_l(1:imax))+b_1*c_1/d_1)
  pm1(1:imax)  = -(a_1*z_on_l(1:imax)+b_1*(z_on_l(1:imax)-(c_1/d_1))*exp(-d_1*z_on_l(1:imax))+b_1*c_1/d_1)
  integralh(1:imax)   = neutral(1:imax)-(ph1(1:imax)-ph0(1:imax))
  integralq(1:imax)   = integralh(1:imax)
  integralm(1:imax)   = neutral(1:imax)-(pm1(1:imax)-pm0(1:imax))
  integralm10(1:imax) = neutral10(1:imax)-(pm1(1:imax)-pm10(1:imax))
endwhere
integralh(1:imax)   = max(integralh(1:imax), 1.e-10)
integralq(1:imax)   = max(integralq(1:imax), 1.e-10)
integralm(1:imax)   = max(integralm(1:imax), 1.e-10)
integralm10(1:imax) = max(integralm10(1:imax), 1.e-10)

! 2m air temperature
tscrn(1:imax)  = temp(1:imax) - tstar(1:imax)*integralh(1:imax)/vkar
! 2m mixing ratio
esatb(1:imax)  = establ(tscrn(1:imax))
qscrn(1:imax)  = mixr(1:imax) - qstar(1:imax)*integralq(1:imax)/vkar
qscrn(1:imax)  = max(qscrn(1:imax), qgmin)
! 2m relative humidity
qsatb(1:imax)  = 0.622*esatb(1:imax)/(ps(1:imax)-esatb(1:imax))
rhscrn(1:imax) = 100.*min(qscrn(1:imax)/qsatb(1:imax), 1.)
! 2m wind speed
uscrn(1:imax)  = max(umagn(1:imax)-ustar(1:imax)*integralm(1:imax)/vkar, 0.)
! 10m wind speed
u10(1:imax)    = max(umagn(1:imax)-ustar(1:imax)*integralm10(1:imax)/vkar, 0.)

return
end subroutine screencalc
      
subroutine autoscrn(is,ie)
      
use arrays_m
use const_phys
use estab
use extraout_m
use liqwpar_m
use mlo_ctrl
use morepbl_m, only : urban_tas, urban_zom, urban_zoh, urban_zoq, &
                      urban_ts, urban_wetfac, wsgs, fg, eg
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
use tkeeps, only : tke
use work2_m
      
implicit none
      
integer, intent(in) :: is, ie
integer :: tile, iq, k, k850, k950
real :: rsig, wspd1, wspd2, u850, u950, x850, x950
real :: thetav, rhos, wt0, wq0, wtv0, z_on_l, fl,pil, th
real, dimension(is:ie) :: umag, zminx, smixr
real, dimension(is:ie) :: ou, ov, atu, atv, iu, iv
real, dimension(is:ie) :: au, av, es, rho
real, dimension(is:ie) :: u_qgscrn, u_rhscrn, u_uscrn, u_u10
real, dimension(is:ie) :: u_ustar, u_tstar, u_qstar, u_thetavstar
real, dimension(is:ie) :: new_zo, new_zoh, new_zoq, new_tss, new_smixr

tile = ie/imax

zminx(is:ie) = bet(1)*t(is:ie,1)/grav
zminx(is:ie) = max( zminx(is:ie), 5. )
iu(is:ie) = 0.
iv(is:ie) = 0.
ou(is:ie) = 0.
ov(is:ie) = 0.
if ( nmlo/=0 ) then
  call mloexport("u",ou(is:ie),1,0,water_g(tile),depth_g(tile))
  call mloexport("v",ov(is:ie),1,0,water_g(tile),depth_g(tile))
  call mloexpice("u",iu(is:ie),0,ice_g(tile),depth_g(tile))
  call mloexpice("v",iv(is:ie),0,ice_g(tile),depth_g(tile))
  ou(is:ie) = (1.-fracice(is:ie))*ou(is:ie) + fracice(is:ie)*iu(is:ie)
  ov(is:ie) = (1.-fracice(is:ie))*ov(is:ie) + fracice(is:ie)*iv(is:ie)
end if
au(is:ie)   = u(is:ie,1) - ou(is:ie)
av(is:ie)   = v(is:ie,1) - ov(is:ie)
umag(is:ie) = max( sqrt(au(is:ie)*au(is:ie)+av(is:ie)*av(is:ie)), vmodmin )
es(is:ie) = establ(tss(is:ie))
qsttg(is:ie) = 0.622*es(is:ie)/(ps(is:ie)-es(is:ie))
smixr(is:ie) = wetfac(is:ie)*qsttg(is:ie) + (1.-wetfac(is:ie))*min( qsttg(is:ie), qg(is:ie,1) )

new_zo(is:ie) = zo(is:ie)
new_zoh(is:ie) = zoh(is:ie)
new_zoq(is:ie) = zoq(is:ie)
new_tss(is:ie) = tss(is:ie)
new_smixr(is:ie) = smixr(is:ie)
if ( zo_clearing>0. ) then
  where ( land(is:ie) )  
    new_zo(is:ie)  = min( zo_clearing, new_zo(is:ie) )
    new_zoh(is:ie) = min( 0.1*zo_clearing, new_zoh(is:ie) )
    new_zoq(is:ie) = min( 0.1*zo_clearing, new_zoq(is:ie) )
  end where  
end if
if ( rescrn==2 ) then
  where ( land(is:ie) )
    new_zo(is:ie) = max( new_zo(is:ie), zobgin )
    new_zoh(is:ie) = max( new_zoh(is:ie), zobgin/7.4 )
    new_zoq(is:ie) = max( new_zoq(is:ie), zobgin/7.4 )
  end where
end if

call screencalc(ie-is+1,qgscrn(is:ie),rhscrn(is:ie),tscrn(is:ie),uscrn(is:ie),u10(is:ie), &
                ustar(is:ie),tstar(is:ie),qstar(is:ie),thetavstar(is:ie),new_zo(is:ie),   &
                new_zoh(is:ie),new_zoq(is:ie),new_tss(is:ie),t(is:ie,1),new_smixr(is:ie), &
                qg(is:ie,1),umag(is:ie),ps(is:ie),zminx(is:ie),sig(1))

rho(is:ie) = ps(is:ie)/(rdry*tss(is:ie))
cduv(is:ie) = ustar(is:ie)**2/max(umag(is:ie),vmodmin)
taux(is:ie) = rho(is:ie)*cduv(is:ie)*au(is:ie)
tauy(is:ie) = rho(is:ie)*cduv(is:ie)*av(is:ie)

atu(is:ie)   = au(is:ie)*uscrn(is:ie)/umag(is:ie) + ou(is:ie)
atv(is:ie)   = av(is:ie)*uscrn(is:ie)/umag(is:ie) + ov(is:ie)
uscrn(is:ie) = sqrt(atu(is:ie)*atu(is:ie)+atv(is:ie)*atv(is:ie))
atu(is:ie)   = au(is:ie)*u10(is:ie)/umag(is:ie) + ou(is:ie)
atv(is:ie)   = av(is:ie)*u10(is:ie)/umag(is:ie) + ov(is:ie)      
u10(is:ie)   = sqrt(atu(is:ie)*atu(is:ie)+atv(is:ie)*atv(is:ie))


! Calculate wind gusts
select case(ugs_meth)
  case(0) ! Schreur et al (2008) "Theory of a TKE based parameterisation of wind gusts" HIRLAM newsletter 54.
    wsgs(is:ie) = wgcoeff*2.185*ustar(is:ie) + u10(is:ie)
  case(1) ! Schreur et al (2008) "Theory of a TKE based parameterisation of wind gusts" HIRLAM newsletter 54.
    if ( .not.(nvmix==6.or.nvmix==9) ) then
      write(6,*) "ERROR: ugs_meth=1 requires nvmix=6 or nvmix=9"
      stop -1
    end if
    do iq = is,ie
      rsig = (1.-0.069*exp(-2.3*u10(iq)*wg_tau/10.))*exp(-0.116*(u10(iq)*wg_tau/10.)**0.555)
      wsgs(iq) = wgcoeff*rsig*sqrt(2.*tke(iq,1)) + u10(iq)
    end do
  case(2)
    if ( .not.(nvmix==6.or.nvmix==9) ) then
      write(6,*) "ERROR: ugs_meth=2 requires nvmix=6 or nvmix=9"
      stop -1
    end if
    do iq = is,ie
      wsgs(iq) = 3.5*sqrt(tke(iq,1)) + u10(iq)
    end do
  case(3)
    wsgs(is:ie) = wgcoeff*2.185*0.4/log(10./0.03)*u10(is:ie) + u10(is:ie)  
  case(4) ! ECMWF IFS DOCUMENTATION � Cy33r1
    ! https://www.ecmwf.int/sites/default/files/elibrary/2009/9227-part-iv-physical-processes.pdf
    do k = 1,kl-1
      if ( sig(k)>=0.95 .and. sig(k+1)<=0.95 ) then
        k950 = k
        x950 = (0.95-sig(k))/(sig(k+1)-sig(k))
      end if
      if ( sig(k)>=0.85 .and. sig(k+1)<=0.85 ) then
        k850 = k
        x850 = (0.85-sig(k))/(sig(k+1)-sig(k))
        exit
      end if
    end do  
    do iq = is,ie
      wspd1 = sqrt(u(iq,k850)**2+v(iq,k850)**2)
      wspd2 = sqrt(u(iq,k850+1)**2+v(iq,k850+1)**2)
      u850 = wspd1*(1.-x850) + wspd2*x850
      wspd1 = sqrt(u(iq,k950)**2+v(iq,k950)**2)
      wspd2 = sqrt(u(iq,k950+1)**2+v(iq,k950+1)**2)
      u950 = wspd1*(1.-x950) + wspd2*x950
      !pil = (ps(iq)/1.e5)**(rdry/cp)
      !th = t(iq,1)/pil
      !thetav = th*(1.+0.61*qg(iq,1)-qlg(iq,1)-qfg(iq,1))
      !rhos = ps(iq)/(rdry*tss(iq))
      !wt0 = fg(iq)/(rhos*cp)
      !wq0 = eg(iq)/(rhos*hl)
      !wtv0 = wt0 + th*0.61*wq0
      !z_on_l = -0.4*1000.*grav*wtv0/(thetav*max(ustar(iq)**3,1.E-10))
      !if ( z_on_l < 0. ) then
      !  fl = (1.-z_on_l*0.5/12.)**3
      !else
      !  fl = 1.
      !end if
      !wsgs(iq) = u10(iq) + 7.2*fl*ustar(iq) + 0.3*max(0.,u850-u950)
      wsgs(iq) = u10(iq) + 7.2*ustar(iq) + 0.3*max(0.,u850-u950)
    end do  
  case(10) 
    ! Schreur et al (2008) "Theory of a TKE based parameterisation of wind gusts" HIRLAM newsletter 54.
    ! + convective wind gust  
    do k = 1,kl-1
      if ( sig(k)>=0.95 .and. sig(k+1)<=0.95 ) then
        k950 = k
        x950 = (0.95-sig(k))/(sig(k+1)-sig(k))
      end if
      if ( sig(k)>=0.85 .and. sig(k+1)<=0.85 ) then
        k850 = k
        x850 = (0.85-sig(k))/(sig(k+1)-sig(k))
        exit
      end if
    end do  
    do iq = is,ie
      wspd1 = sqrt(u(iq,k850)**2+v(iq,k850)**2)
      wspd2 = sqrt(u(iq,k850+1)**2+v(iq,k850+1)**2)
      u850 = wspd1*(1.-x850) + wspd2*x850
      wspd1 = sqrt(u(iq,k950)**2+v(iq,k950)**2)
      wspd2 = sqrt(u(iq,k950+1)**2+v(iq,k950+1)**2)
      u950 = wspd1*(1.-x950) + wspd2*x950      
      wsgs(iq) = wgcoeff*2.185*ustar(iq) + u10(iq) + 0.3*max(0.,u850-u950)
    end do  
  case(11)
    ! Schreur et al (2008) "Theory of a TKE based parameterisation of wind gusts" HIRLAM newsletter 54.
    ! + convective wind gust  
    if ( .not.(nvmix==6.or.nvmix==9) ) then
      write(6,*) "ERROR: ugs_meth=11 requires nvmix=6 or nvmix=9"
      stop -1
    end if
    do k = 1,kl-1
      if ( sig(k)>=0.95 .and. sig(k+1)<=0.95 ) then
        k950 = k
        x950 = (0.95-sig(k))/(sig(k+1)-sig(k))
      end if
      if ( sig(k)>=0.85 .and. sig(k+1)<=0.85 ) then
        k850 = k
        x850 = (0.85-sig(k))/(sig(k+1)-sig(k))
        exit
      end if
    end do  
    do iq = is,ie
      wspd1 = sqrt(u(iq,k850)**2+v(iq,k850)**2)
      wspd2 = sqrt(u(iq,k850+1)**2+v(iq,k850+1)**2)
      u850 = wspd1*(1.-x850) + wspd2*x850
      wspd1 = sqrt(u(iq,k950)**2+v(iq,k950)**2)
      wspd2 = sqrt(u(iq,k950+1)**2+v(iq,k950+1)**2)
      u950 = wspd1*(1.-x950) + wspd2*x950      
      rsig = (1.-0.069*exp(-2.3*u10(iq)*wg_tau/10.))*exp(-0.116*(u10(iq)*wg_tau/10.)**0.555)
      wsgs(iq) = wgcoeff*rsig*sqrt(2.*tke(iq,1)) + u10(iq) + 0.3*max(0.,u850-u950)
    end do
  case(21) ! Schreur et al (2008) "Theory of a TKE based parameterisation of wind gusts" HIRLAM newsletter 54.
    if ( .not.(nvmix==6.or.nvmix==9) ) then
      write(6,*) "ERROR: ugs_meth=21 requires nvmix=6 or nvmix=9"
      stop -1
    end if
    do iq = is,ie
      wsgs(iq) = wgcoeff*sqrt(2.*tke(iq,1)) + u10(iq)
    end do    
  case default
    write(6,*) "ERROR: Unknown method ugs_meth = ",ugs_meth
    stop -1
end select


! urban tile
new_zo(is:ie) = zo(is:ie)
new_zoh(is:ie) = zoh(is:ie)
new_zoq(is:ie) = zoq(is:ie)
new_tss(is:ie) = tss(is:ie)
new_smixr(is:ie) = smixr(is:ie)
where ( sigmu(is:ie)>0. )
  new_tss(is:ie)   = urban_ts(is:ie)
  new_smixr(is:ie) = urban_wetfac(is:ie)*qsttg(is:ie) + (1.-urban_wetfac(is:ie))*min( qsttg(is:ie), qg(is:ie,1) )
  new_zo(is:ie)  = urban_zom(is:ie)
  new_zoh(is:ie) = urban_zoh(is:ie)
  new_zoq(is:ie) = urban_zoq(is:ie)
end where
call screencalc(ie-is+1,u_qgscrn(is:ie),u_rhscrn(is:ie),urban_tas(is:ie),u_uscrn(is:ie),u_u10(is:ie),   &
                u_ustar(is:ie),u_tstar(is:ie),u_qstar(is:ie),u_thetavstar(is:ie),new_zo(is:ie),         &
                new_zoh(is:ie),new_zoq(is:ie),new_tss(is:ie),t(is:ie,1),new_smixr(is:ie),               &
                qg(is:ie,1),umag(is:ie),ps(is:ie),zminx(is:ie),sig(1))

return
end subroutine autoscrn
    
! This version provides a local update for use with physical parameterisations    
subroutine update_u10m(is,ie,u10m)
      
use arrays_m, only : t, u, v, qg, ps
use const_phys, only : grav
use estab
use mlo_ctrl
use newmpar_m
use parm_m
use pbl_m, only : tss
use sigs_m, only : bet, sig
use soilsnow_m, only : fracice
use work2_m, only : zo, zoh, zoq, wetfac
      
implicit none
      
integer, intent(in) :: is, ie
integer :: tile
real, dimension(imax), intent(out) :: u10m
real, dimension(is:ie) :: umag, zminx, smixr
real, dimension(is:ie) :: ou, ov, atu, atv, iu, iv
real, dimension(is:ie) :: au, av, es, rho
real, dimension(is:ie) :: cduv, qsttg, qgscrn, tscrn
real, dimension(is:ie) :: rhscrn, uscrn, u10, u10gs_var
real, dimension(is:ie) :: ustar, tstar, qstar, thetavstar
    
tile = ie/imax

zminx(is:ie) = bet(1)*t(is:ie,1)/grav
zminx(is:ie) = max( zminx(is:ie), 5. )
if ( nmlo/=0 ) then
  iu(is:ie) = 0.
  iv(is:ie) = 0.
  ou(is:ie) = 0.
  ov(is:ie) = 0.
  call mloexport("u",ou(is:ie),1,0,water_g(tile),depth_g(tile))
  call mloexport("v",ov(is:ie),1,0,water_g(tile),depth_g(tile))
  call mloexpice("u",iu(is:ie),0,ice_g(tile),depth_g(tile))
  call mloexpice("v",iv(is:ie),0,ice_g(tile),depth_g(tile))
  ou(is:ie) = (1.-fracice(is:ie))*ou(is:ie) + fracice(is:ie)*iu(is:ie)
  ov(is:ie) = (1.-fracice(is:ie))*ov(is:ie) + fracice(is:ie)*iv(is:ie)
else
  ou(is:ie) = 0.
  ov(is:ie) = 0.
end if
au(is:ie)   = u(is:ie,1) - ou(is:ie)
av(is:ie)   = v(is:ie,1) - ov(is:ie)
umag(is:ie) = max( sqrt(au(is:ie)*au(is:ie)+av(is:ie)*av(is:ie)), vmodmin )
es(is:ie) = establ(tss(is:ie))
qsttg(is:ie)  = 0.622*es(is:ie)/(ps(is:ie)-es(is:ie))
smixr(is:ie) = wetfac(is:ie)*qsttg(is:ie) + (1.-wetfac(is:ie))*min( qsttg(is:ie), qg(is:ie,1) )

call screencalc(ie-is+1,qgscrn(is:ie),rhscrn(is:ie),tscrn(is:ie),uscrn(is:ie),u10(is:ie),   &
                ustar(is:ie),tstar(is:ie),qstar(is:ie),thetavstar(is:ie),  &
                zo(is:ie),zoh(is:ie),zoq(is:ie),tss(is:ie),t(is:ie,1),smixr(is:ie),         &
                qg(is:ie,1),umag(is:ie),ps(is:ie),zminx(is:ie),sig(1))

atu(is:ie) = au(is:ie)*u10(is:ie)/umag(is:ie) + ou(is:ie)
atv(is:ie) = av(is:ie)*u10(is:ie)/umag(is:ie) + ov(is:ie)      
u10m(:)    = sqrt(atu(is:ie)*atu(is:ie)+atv(is:ie)*atv(is:ie))

return
end subroutine update_u10m    


! Calculate CAPE and CIN
! This differs from the value used in convection, but is comparable with
! other models.  Hence we refer to this as diagnosed CAPE and CIN.
subroutine capecalc

use arrays_m
use const_phys
use estab
use newmpar_m
use prec_m
use sigs_m

implicit none

integer js, je, tile
integer k, n, iq, icount, nloop, ktop
integer, parameter :: kmax = 1 ! default for source parcel at surface
real, dimension(imax,kl) :: pl, tl, pil, th, thv
real, dimension(imax) :: th2, pl2, tl2, thv2, qv2, b2
real, dimension(imax) :: narea, ql2, qi2, qt, capel, cinl
real, dimension(imax) :: qs
real b1, dp, pl1, tl1, th1, qv1, ql1, qi1, thv1, thlast, pil2, fliq, fice
real tbarl, qvbar, qlbar, qibar, lhv, lhs, lhf, rm, cpm, qsat_save
real dz, frac, parea
!real, parameter :: pinc = 100. ! Pressure increment (Pa) - smaller is more accurate
real, parameter :: pinc = 1000.
real, parameter :: lv1 = 2501000. + (4190.-cpv)*273.15
real, parameter :: lv2 = 4190. - cpv
real, parameter :: ls1 = 2836017. + (2118.636-cpv)*273.15
real, parameter :: ls2 = 2188.636 - cpv
logical not_converged

! Following code is based on Bryan (NCAR) citing
! Bolton 1980 MWR p1046 and Bryan and Fritsch 2004 MWR p2421

#ifdef GPU
!$acc parallel loop copyin(ps,t,ntiles,imax,kl) copyout(cape_d,cin_d) present(sig)       &
!$acc   private(pl,tl,pil,th,thv,th2,pl2,tl2,thv2,qv2,ql2,qi2,qt,b2,narea,capel,cinl,qs) &
!$acc   private(js,je,k,iq,tile,b1,dp,nloop,pl1,tl1,th1,qv1,ql1,qi1,thv1,thlast,pil2)    &
!$acc   private(icount,not_converged,fliq,fice,tbarl,qvbar,qlbar,qibar,rm,cpm,qsat_save) &
!$acc   private(dz,frac,parea,n,ktop,lhv,lhs,lhf)
#else
!$omp do schedule(static) private(pl,tl,pil,th,thv,th2,pl2,tl2,thv2,qv2,ql2,qi2,qt,b2) &
!$omp   private(narea,capel,cinl,qs,js,je,k,iq,tile,b1,dp,nloop,pl1,tl1,th1,qv1,ql1)   &
!$omp   private(qi1,thv1,thlast,pil2,icount,not_converged,fliq,fice,tbarl,qvbar,qlbar) &
!$omp   private(qibar,rm,cpm,qsat_save,dz,frac,parea,n,ktop,lhv,lhs,lhf)
#endif
do tile = 1,ntiles
  js = (tile-1)*imax + 1
  je = tile*imax    

  ! limit height to 10hPa
  ktop = 1
  do k = 1,kl
    if ( 1.e5*sig(k)>1.e4 ) ktop = k
  end do
  
  ! create local variables
  do k = 1,ktop
    pl(:,k) = ps(js:je)*sig(k)
    tl(:,k) = t(js:je,k)
    pil(:,k) = (pl(:,k)/1.e5)**(rdry/cp)
    qs(:) = qsat(pl(:,k),tl(:,k)) 
    th(:,k) = tl(:,k)/pil(:,k)
    thv(:,k) = th(:,k)*(1.+1.61*qs(:))/(1.+qs(:))
  end do  

  ! define initial parcel properties
  th2(:) = th(:,kmax)
  pl2(:) = pl(:,kmax)
  tl2(:) = tl(:,kmax)
  thv2(:) = thv(:,kmax)
  qv2(:) = qsat(pl(:,kmax),tl(:,kmax)) 
  ql2(:) = 0.
  qi2(:) = 0.
  qt(:) = qv2(:)
  b2(:) = 0.
  narea(:) = 0.
  capel(:) = 0.
  cinl(:) = 0.
  !cape_lvl(:,:) = 0.
  !cin_lvl(:,:) = 0.

  ! start ascent of parcel
  do k = kmax+1,ktop
    do iq = 1,imax  
      
      b1 = b2(iq)
      nloop = 1 + int( 1.e5*(sig(k-1)-sig(k))/pinc )
      dp = (pl(iq,k-1)-pl(iq,k))/real(nloop)  
  
      do n = 1,nloop
            
        pl1 = pl2(iq)
        tl1 = tl2(iq)
        th1 = th2(iq)
        qv1 = qv2(iq)
        ql1 = ql2(iq)
        qi1 = qi2(iq)
        thv1 = thv2(iq)
     
        pl2(iq) = pl2(iq) - dp
        pil2 = (pl2(iq)/1.e5)**(rdry/cp)
        thlast = th1
      
        not_converged = .true.
        do icount = 1,50
          if ( not_converged ) then

            tl2(iq) = thlast*pil2
            fliq = max(min((tl2(iq)-233.15)/(273.15-233.15),1.),0.)
            fice = 1. - fliq

            qsat_save = qsat(pl2(iq),tl2(iq))

            qv2(iq) = min( qt(iq), qsat_save )
            qi2(iq) = max( fice*(qt(iq)-qv2(iq)), 0. )
            ql2(iq) = max( qt(iq)-qv2(iq)-qi2(iq), 0. )

            tbarl = 0.5*(tl1+tl2(iq))
            qvbar = 0.5*(qv1+qv2(iq))
            qlbar = 0.5*(ql1+ql2(iq))
            qibar = 0.5*(qi1+qi2(iq))

            lhv = lv1 - lv2*tbarl
            lhs = ls1 - ls2*tbarl
            lhf = lhs - lhv

            rm = rdry + rvap*qvbar
            cpm = cp + cpv*qvbar + 4190.*qlbar + 2118.636*qibar
            th2(iq) = th1*exp( lhv*(ql2(iq)-ql1)/(cpm*tbarl)    &
                   +lhs*(qi2(iq)-qi1)/(cpm*tbarl)               &
                   +(rm/cpm-rdry/cp)*log(pl2(iq)/pl1) )

            if ( abs(th2(iq)-thlast)>0.0002 ) then
              thlast = thlast + 0.2*(th2(iq)-thlast)
            else
              not_converged = .false.
            end if
            
          end if  ! not_converged  
        end do    ! ncount

        ! pseudoadiabat
        qt(iq)  = qv2(iq)
        ql2(iq) = 0.
        qi2(iq) = 0.

      end do     ! n loop  

      thv2(iq) = th2(iq)*(1.+1.61*qv2(iq))/(1.+qv2(iq)+ql2(iq)+qi2(iq))
      b2(iq) = grav*( thv2(iq)-thv(iq,k) )/thv(iq,k)
      dz = -(cp/grav)*0.5*(thv(iq,k)+thv(iq,k-1))*(pil(iq,k)-pil(iq,k-1))

      ! calculate contributions to CAPE and CIN
      if ( b2(iq)>=0. .and. b1<0. ) then
        ! first time entering positive region
        frac = b2(iq)/(b2(iq)-b1)
        parea = 0.5*b2(iq)*dz*frac
        narea(iq) = narea(iq)-0.5*b1*dz*(1.-frac)
        cinl(iq) = cinl(iq) + narea(iq)
        capel(iq) = capel(iq) + max(0.,parea)
        narea(iq) = 0.
      else if ( b2(iq)<0. .and. b1>0. ) then  
        ! first time entering negative region  
        frac = b1/(b1-b2(iq))
        parea = 0.5*b1*dz*frac
        narea(iq) = -0.5*b2(iq)*dz*(1.-frac)
        capel(iq) = capel(iq) + max(0.,parea)
      else if ( b2(iq)<0. ) then 
        ! continue negative buoyancy region
        parea = 0.
        narea(iq) = narea(iq)-0.5*dz*(b1+b2(iq))
      else
        ! continue positive buoyancy region  
        parea = 0.5*dz*(b1+b2(iq))
        narea(iq) = 0.
        capel(iq) = capel(iq) + max(0.,parea)
      end if
  
      !cape_lvl(iq,k) = parea
      !cin_lvl(iq,k) = narea(iq)
      
    end do ! iq loop  
  end do   ! k loop
    
  cape_d(js:je) = capel(:)
  cin_d(js:je) = cinl(:)

end do
#ifdef GPU
!$acc end parallel loop
#else
!$omp end do nowait
#endif

return
end subroutine capecalc
