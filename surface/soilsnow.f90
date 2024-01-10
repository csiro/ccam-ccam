! Conformal Cubic Atmospheric Model
    
! Copyright 2015 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
!    N.B. check settings of  parameter (ncondxpr & nmeth). snmin to parm.h
!    This batch of code has only the vector version soilsnowv
!************************* soilsnowv follows  ****some to be vectorized*****
subroutine soilsnowv

use cc_mpi, only : mydiag
use const_phys
use diag_m
use morepbl_m  ! need runoff
use newmpar_m
use nsibd_m    ! soilm
use parm_m
use permsurf_m
use soil_m     ! land
use soilsnow_m
use soilv_m
use work3_m
use work3b_m

implicit none

integer, parameter :: ntest=0   ! 3: forces 3-layer snow, 1: for snow diag prints
!        for snow diag prints set ntest to 1 throughout
!        or, usefully can edit 'ntest>0' to 'ktau>nnn'
!----------------------------------------------------------------------
! Inputs:
!     iq   - current grid point
!     isoil - soil type
!     dt    - time step
!     timi  - start time step
!     ga    - ground heat flux W/m^2
!     condxpr - precip (liquid and solid) 
!     fes   - soil evaporation (W/m2)
! Output
!     runoff - total runoff
!----------------------------------------------------------------------

integer k,iq,ip,isoil
real, dimension(ifull) :: ggflux
real, dimension(ifull,ms) :: gammzz
real, dimension(3) :: etac
real csice,cswat,rhowat,cgsnow,rhosnow
real z1snow,ccoef,tggd,tr11,sdd
real ossdn2,tr1,xx,pr,osm1,excm
real excd,osm2,osm3,sd1,sm1,sfl
real sicefreeze,sicemelt

data csice /2.100e3/, cswat /4.218e3/, rhowat /1000./  ! for calgammv
data cgsnow/2090./,rhosnow/200./                       ! for calgammv

!     update land points.
if(ktau==1)then
  call work3b_init(ifull,ms)
  if(ntest==3)snmin=.11   ! to force 3-layer snow for testing
  ! N.B. snmin should exceed sum of layer depths, i.e. .11 m	 
  do k=1,ms
    do iq=1,ifull
      wblf(iq,k)   = -99.  ! just to put something in array
    enddo
  enddo
  do iq=1,ifull
    if(isflag(iq)==0)then
      ssdnn(iq) = ssdn(iq,1)
    else
      sdepth(iq,1)=smass(iq,1)/ssdn(iq,1)
      sdepth(iq,2)=smass(iq,2)/ssdn(iq,2)
      sdepth(iq,3)=smass(iq,3)/ssdn(iq,3)
      z1snow = sdepth(iq,1)+sdepth(iq,2)+sdepth(iq,3) ! real snow depth in m
      ssdnn(iq)  = snowd(iq)/z1snow
    endif
  enddo
endif   ! (ktau==1)

snowflx(:)=0.   ! 23/12/05
if ( mydiag ) sdepth(idjd,1:3)=0.  ! default to avoid *** in globpe.f display
do ip=1,ipland         ! all land points in this nsib>=1 loop
  iq=iperm(ip)
  isoil = isoilm(iq)
  do k=1,ms
    wblf(iq,k)   = (wb(iq,k)-wbice(iq,k))/ssat(isoil)  ! for stempv
    wbfice(iq,k) = wbice(iq,k)/ssat(isoil)
  enddo
  if(snowd(iq)<1.e-7) then
    isflag(iq) = 0
    ssdn(iq,1) = 140.
    ssdnn(iq)  = 140.
    snowd(iq)=0.     ! to reduce any risk of underflow
  elseif(snowd(iq)<snmin*ssdnn(iq)) then
    ccoef = 0.
    if(ssdn(iq,1) >=150.) ccoef = .046
    tggd        = min(tfrz,tgg(iq,1))
    if(isflag(iq)==1)then
      tggd = min(tfrz,tggsn(iq,1))
      ssdn(iq,1) = ssdnn(iq)
    endif
    ssdn(iq,1)  = max( 140. , ssdn(iq,1)+dt*ssdn(iq,1)*2.8e-6* exp(-.03*(273.1-tggd)-ccoef*(ssdn(iq,1)-150.)) )
    tr11 =  max(0.,1.-osnowd(iq)/snowd(iq))
    ssdn(iq,1) = tr11*140.+ (1.-tr11)*ssdn(iq,1)
    ssdnn(iq)   = ssdn(iq,1)
    isflag(iq)  = 0
    sdepth(iq,1)    = snowd(iq)/ssdn(iq,1)

  else  ! sufficient snow now
    if(isflag(iq)==0) then
      tggsn(iq,1) = tgg(iq,1)
      tggsn(iq,2) = tgg(iq,1)
      tggsn(iq,3) = tgg(iq,1)
      ssdn(iq,2)  = ssdn(iq,1)
      ssdn(iq,3)  = ssdn(iq,1)
      sdepth(iq,1)= .07
      sdd=snowd(iq)/ssdn(iq,1)-.07
      if(snowd(iq)>20.)then
        sdepth(iq,2)=max(.02 , .3*sdd)
        sdepth(iq,3)=max(.02 , .7*sdd)
      else
        sdepth(iq,2)=max(.02 , .45*sdd)
        sdepth(iq,3)=max(.02 , .55*sdd)
      endif
      smass(iq,1) = .07*ssdn(iq,1)
      smass(iq,2) = sdepth(iq,2)*ssdn(iq,2)
      smass(iq,3) = sdepth(iq,3)*ssdn(iq,3)
    endif

    ossdn2 = ssdn(iq,2)
    if(ntest>0.and.iq==idjd.and.mydiag)then
      write(6,*) 'soilsnow  snowd,osnowd ',snowd(iq),osnowd(iq)
      write(6,*) 'ssdn a ',(ssdn(iq,k),k=1,3)
      write(6,*) 'smass a ',(smass(iq,k),k=1,3)
      write(6,*) 'sdepth a ',(sdepth(iq,k),k=1,3)
    endif
    do k=1,3
      ccoef  = 0.
      if(ssdn(iq,k) >=150.) ccoef=4.6e-2
      tggd   = min(tfrz,tggsn(iq,k))
      ssdn(iq,k) = ssdn(iq,k)+dt*ssdn(iq,k)*3.1e-6*exp(-.03*(273.1-tggd)-ccoef*(ssdn(iq,k)-150.))
      etac(k)=3.e7*exp( .021*ssdn(iq,k)+.081*(273.1-tggd) )  ! same as:
    enddo
    ssdn(iq,1)=ssdn(iq,1)+dt*grav*.5*.07*ssdn(iq,1)*ssdn(iq,1)/etac(1)
    ssdn(iq,2)=ssdn(iq,2)+dt*grav*ssdn(iq,2)*(.07*ssdn(iq,1)+.5*smass(iq,2))/etac(2)
    ssdn(iq,3)=ssdn(iq,3)+dt*grav*ssdn(iq,3)*(.07*ssdn(iq,1)+smass(iq,2)+.5*smass(iq,3))/etac(3)
    tr1  =  snowd(iq)-osnowd(iq)
    xx=max(0. , .07-smass(iq,1)/ssdn(iq,1))
    pr=min(smass(iq,2)/(smass(iq,3)+smass(iq,2)) , .9)
    if(ntest>0.and.iq==idjd.and.mydiag)then
      write(6,*) 'ssdn b ',(ssdn(iq,k),k=1,3)
    endif
    if(tr1>=0.) then
      ssdn(iq,1)=max((smass(iq,1)+tr1)/(smass(iq,1)/ssdn(iq,1)+tr1/140.) , 140.)
      osm1        = smass(iq,1)
      smass(iq,1) = .07*ssdn(iq,1)
      sdepth(iq,1)= .07
      excm        = osm1+tr1-smass(iq,1)
      excd        = excm/ssdn(iq,1)
      osm2        = smass(iq,2)
      smass(iq,2) = max(.01,smass(iq,2)+.4*excm)
      ssdn(iq,2)=max(140. , min(500.,smass(iq,2)/(osm2/ssdn(iq,2)+.4*excd)))
      sdepth(iq,2)=max(.02,smass(iq,2)/ssdn(iq,2))
      osm3        = smass(iq,3)
      smass(iq,3) = max(.01 , snowd(iq)-smass(iq,1)-smass(iq,2))
      sdepth(iq,3)= max(.02 , osm3/ssdn(iq,3)+.6*excm/ssdn(iq,2))
      ssdn(iq,3)  = max(140., min(500. , smass(iq,3)/sdepth(iq,3)))
      if(ssdn(iq,3)<ssdn(iq,2)) then
        ssdn(iq,3) = ssdn(iq,2)
        sdepth(iq,3) =max(.02,smass(iq,3)/ssdn(iq,3))
      endif

    else
      ! snow melting
      sdepth(iq,1)   = .07
      sd1         = max(.005,smass(iq,1)/ssdn(iq,1)) !current depth of
      ! the 1st layer
      sm1         = max(.01,smass(iq,1))   !current mass of the 1st layer
      excd        = .07-sd1
      smass(iq,1)=max(140.*.07,min(500. , sd1*ssdn(iq,1)+excd*ssdn(iq,2) ) )
      ssdn(iq,1)  = smass(iq,1)/.07
      excm        = smass(iq,1)-sm1
      excd        = excm/ssdn(iq,1)
      osm2        = smass(iq,2)
      smass(iq,2) = max(.01,smass(iq,2)-pr*excm)
      sdepth(iq,2)= max(.02,osm2/ssdn(iq,2)-pr*excd)
      ssdn(iq,2)  = max(140.,min(500.,smass(iq,2)/sdepth(iq,2)))
 
      if( ssdn(iq,2) < ossdn2 ) then
        ssdn(iq,2)=ossdn2
        smass(iq,2)=.45*(snowd(iq)-smass(iq,1))
        sdepth(iq,2)=max(.02,smass(iq,2)/ssdn(iq,2))
      endif
      smass(iq,3)  = max(.01 , snowd(iq)-smass(iq,1)-smass(iq,2))
      sdepth(iq,3) = max(.02 , smass(iq,3)/ssdn(iq,3))
    endif   ! (tr1>=0.) .. else ..
    if(ntest>0.and.iq==idjd.and.mydiag)then
      write(6,*) 'ssdn c ',(ssdn(iq,k),k=1,3)
      write(6,*) 'smass c ',(smass(iq,k),k=1,3)
      write(6,*) 'sdepth c ',(sdepth(iq,k),k=1,3)
    endif

    ! --------------------end of snowprv-------------------------------
    isflag(iq) = 1
    z1snow = sdepth(iq,1)+sdepth(iq,2)+sdepth(iq,3) ! real snow depth in m
    ssdnn(iq)  = (ssdn(iq,1)*sdepth(iq,1)+ssdn(iq,2)*sdepth(iq,2) + ssdn(iq,3)*sdepth(iq,3))/z1snow
  endif
enddo   ! land points

if((ntest>0.or.diag).and.mydiag) then
  if (land(idjd))then  !MJT bugfix
    write(6,*) 'in soilsnowv before stempv,  ktau= ',ktau
    write(6,*) 'ga,dt,ssdn ',ga(idjd),dt,(ssdn(idjd,k),k=1,3)
    write(6,*) 'osnowd,snowd,isflag',osnowd(idjd),snowd(idjd),isflag(idjd)
    write(6,*) 'tggsn ',(tggsn(idjd,k),k=1,3)
    write(6,*) 'tgg ',(tgg(idjd,k),k=1,ms)
    write(6,*) 'wb ',(wb(idjd,k),k=1,ms)
    write(6,*) 'wbice ',(wbice(idjd,k),k=1,ms)
    write(6,*) 'wblf ',(wblf(idjd,k),k=1,ms)
    write(6,*) 'wbfice ',(wbfice(idjd,k),k=1,ms)
    write(6,*) 'sdepth c2 ',(sdepth(idjd,k),k=1,3)
  endif
endif

call stempv(gammzz) 

do ip=1,ipland  ! all land points in this nsib=1 or 3 loop
  iq=iperm(ip)
  isoil = isoilm(iq)
  if(isflag(iq)==1 )then
    ggflux(iq)=sgflux(iq)
  else
    ggflux(iq)=gflux(iq)    ! as calculated in stemp
  endif
  do k=1,ms
    if(tgg(iq,k)<tfrz.and. .99*wb(iq,k)-wbice(iq,k)>.001) then
      sfl=(tfrz-tgg(iq,k))*gammzz(iq,k)
      sicefreeze=min(max(0.,(.99*wb(iq,k)-wbice(iq,k)))*zse(k)*1000. ,sfl/hlf)
      wbice(iq,k)=min(wbice(iq,k)+sicefreeze/(zse(k)*1000.),.99*wb(iq,k))
      gammzz(iq,k)=max( (1.-ssat(isoil))*css(isoil)*rhos(isoil)+(wb(iq,k)-wbice(iq,k))*cswat*rhowat+wbice(iq,k)*csice*rhowat*.9, &
                        css(isoil)*rhos(isoil)  ) * zse(k)
      if(k==1.and.isflag(iq)==0)gammzz(iq,k)=gammzz(iq,k) + cgsnow*snowd(iq)  ! changed back 21/5/01
      tgg(iq,k)=tgg(iq,k)+sicefreeze*hlf/gammzz(iq,k)
      if(ntest>0.and.iq==idjd.and.mydiag)then
        write(6,*) 'D k,tgg,sicefreeze,wbice,gammzz ',k,tgg(iq,k),sicefreeze,wbice(iq,k),gammzz(iq,k)
      endif

    elseif( tgg(iq,k)>tfrz.and.wbice(iq,k)>0.) then
      sfl=(tgg(iq,k)-tfrz)*gammzz(iq,k)
      sicemelt=min(wbice(iq,k)*zse(k)*1000.,sfl/hlf)
      wbice(iq,k)=max(0.,wbice(iq,k)-sicemelt/(zse(k)*1000.))
      gammzz(iq,k)=max( (1.-ssat(isoil))*css(isoil)*rhos(isoil)+(wb(iq,k)-wbice(iq,k))*cswat*rhowat+wbice(iq,k)*csice*rhowat*.9, &
                        css(isoil)*rhos(isoil)  ) * zse(k)
      if(k==1.and.isflag(iq)==0)gammzz(iq,k)=gammzz(iq,k) + cgsnow*snowd(iq)  ! changed back 21/5/01
      tgg(iq,k)=tgg(iq,k)-sicemelt*hlf/gammzz(iq,k)
      if(ntest>0.and.iq==idjd.and.mydiag)write(6,*) 'E tgg k',k,tgg(iq,k)
    endif
  enddo   ! k=1,ms
  if(isoil==9)then
    do k=2,ms
      tgg(iq,k)=min(tgg(iq,k),273.1)  ! jlm 7/6/00
    enddo
  endif

enddo   ! ip loop for land points

call surfbv(gammzz)

if((ntest>0.or.diag).and.mydiag) then
  if(land(idjd))then ! MJT bugfix
    write(6,*) 'after surfbv,isflag ',isflag(idjd)
    write(6,*) 'tgg ',(tgg(idjd,k),k=1,ms)
    write(6,*) 'wb ',(wb(idjd,k),k=1,ms)
    write(6,*) 'wblf ',(wblf(idjd,k),k=1,ms)
    write(6,*) 'wbfice ',(wbfice(idjd,k),k=1,ms)
  endif
endif

return 
end subroutine soilsnowv

!***********************************************************************

subroutine surfbv(gammzz)

use arrays_m
use cc_mpi, only : mydiag
use const_phys
use morepbl_m  ! need runoff
use newmpar_m
use nsibd_m
use parm_m
use permsurf_m
use sigs_m
use soil_m     ! land,sice,sicedep,alb
use soilsnow_m
use soilv_m
use work3_m
use work3b_m

implicit none

integer, parameter :: ntest=0    ! 3: forces 3-layer snow, 1: for snow diag prints
integer, parameter :: ncondxpr=1 ! 0: old sfce scheme, 1: jlm mid-level suggestion
integer, parameter :: newsmelt=1 ! 0: old, 1: new from Aug 2003
!integer, parameter :: nglacier=1 ! 0 original, 1 off, 2 new from Eva; to parm.h

integer k,iq,ip,isoil
real smelt,sgamm,segg,evapsn
real snowflx0,snowflx1,snowflxr,snowflx2
real snowflx3,totwet,dtotw,weting
real sinfil,dwb,rnof5,smasstot
real, dimension(ifull) :: rnof1,rnof2
real, dimension(ifull,ms) :: gammzz
real, dimension(3) :: smelt1
real, dimension(9) :: c3
data c3/1.255, .334, .138, .521, .231, .199, .375, .623, .334/

if(ktau==1.and.mydiag)then
  write(6,*) 'ncondxpr,nglacier ',ncondxpr,nglacier
endif
if((ntest>0.or.diag).and.mydiag)then
  if(land(idjd))then !MJT bugfix
    write(6,*) 'entering surfbv  condxpr',condxpr(idjd)
    write(6,*) 'osnowd,snowd,isflag',osnowd(idjd),snowd(idjd),isflag(idjd)
    write(6,*) 'tggsn ',(tggsn(idjd,k),k=1,3)
    write(6,*) 'tgg ',(tgg(idjd,k),k=1,ms)
    write(6,*) 't ',(t(idjd,k),k=1,kl)
    write(6,*) 'wb ',(wb(idjd,k),k=1,ms)
    write(6,*) 'wbice ',(wbice(idjd,k),k=1,ms)
    write(6,*) 'gammzz ',(gammzz(idjd,k),k=1,ms)
  endif
endif
      
do ip=1,ipland  ! all land points in this nsib=1/3 loop
  iq=iperm(ip)
  isoil = isoilm(iq)
  smelt=0.
  osnowd(iq)=snowd(iq)
  if(condspr(iq)>0.)then  ! just using ncondxpr=1 treatment now
    if(isflag(iq)==0)then
      snowd(iq)=max(snowd(iq) + condspr(iq), 0.)
      sno(iq)=sno(iq)+condspr(iq)  ! snow precip accum in mm
      tgg(iq,1)=tgg(iq,1)+condspr(iq)*hlf/gammzz(iq,1)
      condxpr(iq)=condxpr(iq)-condspr(iq)
      condspr(iq)=0.
    else  ! i.e. isflag(iq)=1
      snowd(iq)=max(snowd(iq) + condspr(iq), 0.)
      sno(iq)=sno(iq)+condspr(iq)  ! snow precip accum in mm
      do k=1,3
        sgamm  = ssdn(iq,k)*2105. * sdepth(iq,k)
        tggsn(iq,k)=tggsn(iq,k)+condspr(iq)*hlf*smass(iq,k)/(sgamm*osnowd(iq))
      enddo
      condxpr(iq)=condxpr(iq)-condspr(iq)
      condspr(iq)=0.
    endif     ! (isflag(iq)==0) ... else ...
  endif       ! (condxpr(iq)>0.)

  ! snow evaporation and melting
  segg=fes(iq)/hl
  if(snowd(iq)>.1)then
    evapsn=min(snowd(iq),dt*fes(iq)/(hl+hlf))
    snowd(iq)=snowd(iq)-evapsn
    segg=0.
  endif

  if(snowd(iq)>0.)then
    if(isflag(iq)==0) then
!          snow covered land
!          following done in sflux  via  ga= ... +cls*egg + ...
!          tgg(iq,1)=tgg(iq,1)-evapsn*(hlf+hl)/gammzz(iq,1) 
      if(tgg(iq,1)>=tfrz)then
!           land,snow,melting; N.B.  divide snowflx by dt to get flux
        snowflx0=(tgg(iq,1)-tfrz)*gammzz(iq,1)
!           prevent snow depth going negative
        smelt= min(snowflx0/hlf ,snowd(iq))
        snowd(iq)=snowd(iq)-smelt
        tgg(iq,1)= tgg(iq,1)-smelt*hlf/gammzz(iq,1)
        snowflx(iq)=smelt*hlf/dt
        if(ntest>0.and.iq==idjd.and.mydiag)then
          write(6,*) 'in surfbv b'
          write(6,*) 'tgg ',(tgg(idjd,k),k=1,ms)
          write(6,*) 'snowflx0,snowflx,gammzz,smelt,snowd ',snowflx0,snowflx(iq),gammzz(iq,1),smelt,snowd(iq) 
        endif
      endif  ! (tgg(iq,1)>=tfrz)then
    else     ! 3-layer scheme,  isflag=1
      k=1
      smelt1(k)=0.
      if(tggsn(iq,k)>tfrz)then
        sgamm   = ssdn(iq,k)*2105. * sdepth(iq,k)
        snowflx1=(tggsn(iq,k)-tfrz)*sgamm
        smelt1(k)= min(snowflx1/hlf ,0.9*smass(iq,k))
        smass(iq,k) = smass(iq,k) - smelt1(k)
        tggsn(iq,k)= min(tggsn(iq,k)-smelt1(k)*hlf/sgamm, tfrz)
        if(newsmelt==1.and.isoilm(iq)==9)then !snow melt refreezing   
          sgamm=ssdn(iq,k+1)*2105.*sdepth(iq,k+1)
          snowflxr=smelt1(k)*hlf/dt
          tggsn(iq,k+1)=tggsn(iq,k+1)+snowflxr*dt/sgamm
          ssdn(iq,k+1)=min((smass(iq,k+1)+smelt1(k))/(smass(iq,k+1)/ssdn(iq,k+1)+smelt1(k)/1000.),450.)
          smass(iq,k+1) = smass(iq,k+1) + smelt1(k)
          sdepth(iq,k+1)=smass(iq,k+1)/ssdn(iq,k+1) ! to give new sgamm
          smelt1(k)=0.
        endif  ! (newsmelt==1.and.isoilm(iq)==9)
      endif   ! (tggsn(iq,k)>tfrz)
      k=2
      smelt1(k)=0.
      if(tggsn(iq,k)>tfrz)then
        sgamm   = ssdn(iq,k)*2105. * sdepth(iq,k)
        snowflx2=(tggsn(iq,k)-tfrz)*sgamm
        smelt1(k)= min(snowflx2/hlf ,0.9*smass(iq,k))
        smass(iq,k) = smass(iq,k) - smelt1(k)
        tggsn(iq,k)= min(tggsn(iq,k)-smelt1(k)*hlf/sgamm, tfrz)
        if(newsmelt==1.and.isoilm(iq)==9)then !snow melt refreezing   
          sgamm=ssdn(iq,k+1)*2105.*sdepth(iq,k+1)
          snowflxr=smelt1(k)*hlf/dt
          tggsn(iq,k+1)=tggsn(iq,k+1)+snowflxr*dt/sgamm
          ssdn(iq,k+1)=min((smass(iq,k+1)+smelt1(k))/(smass(iq,k+1)/ssdn(iq,k+1)+smelt1(k)/1000.),450.)
          smass(iq,k+1) = smass(iq,k+1) + smelt1(k)
          sdepth(iq,k+1)=smass(iq,k+1)/ssdn(iq,k+1) ! to give new sgamm
          smelt1(k)=0.
        endif  ! (newsmelt==1.and.isoilm(iq)==9)
      endif   ! (tggsn(iq,k)>tfrz)
      k=3
      smelt1(k)=0.
      if( tggsn(iq,k)>tfrz ) then
        sgamm   = ssdn(iq,k)*2105. * sdepth(iq,k)
        snowflx3=(tggsn(iq,k)-tfrz)*sgamm
        smelt1(k)= min(snowflx3/hlf ,0.9*smass(iq,k))
        smass(iq,k) = smass(iq,k) - smelt1(k)
        tggsn(iq,k)= min(tggsn(iq,k)-smelt1(k)*hlf/sgamm, tfrz)
      endif   ! (tggsn(iq,k)>tfrz)
      smelt=smelt1(1)+smelt1(2)+smelt1(3)
      snowd(iq)=snowd(iq) - smelt
      snowflx(iq)=smelt*hlf/dt
    endif  ! (isflag(iq)==0) .. else ..
  endif    !  (snowd(iq)>0.)t

  totwet=condxpr(iq)+smelt
  snowmelt(iq) = snowmelt(iq) + smelt
  dtotw=totwet*86400./dt
  rnof1(iq)=max(0. ,dtotw-150.)*(dt/86400.)   ! presumably in mm
  weting=totwet-rnof1(iq)
  sinfil=.8*min((ssat(isoil)-wb(iq,1))*zse(1)*1000.,weting)
  rnof1(iq)=rnof1(iq)+max(0. , weting-sinfil)
  weting=totwet-rnof1(iq)
  fwtop(iq)=weting/dt-segg
enddo               ! ip loop

if((ntest>0.or.diag).and.mydiag)then
  if(land(idjd))then !MJT bugfix
    write(6,*) 'in surfbv before smoisturev  condxpr',condxpr(idjd)
    write(6,*) 'osnowd,snowd,isflag',osnowd(idjd),snowd(idjd),isflag(idjd)
    write(6,*) 'tggsn_c ',(tggsn(idjd,k),k=1,3)
    write(6,*) 'tgg ',(tgg(idjd,k),k=1,ms)
  endif
endif

call smoisturev

if((ntest>0.or.diag).and.mydiag)then
  if(land(idjd))then !MJT bugfix
    write(6,*) 'in surfbv after smoisturev '
    write(6,*) 'osnowd,snowd,isflag,ssat,runoff',osnowd(idjd),snowd(idjd),isflag(idjd),ssat(isoilm(idjd)),runoff(idjd)
    write(6,*) 'tggsn_d ',(tggsn(idjd,k),k=1,3)
  endif
endif

do ip=1,ipland  ! all land points in this nsib=1  or 3 loop
  iq=iperm(ip)
  isoil = isoilm(iq)
  do k=1,ms
    rnof1(iq)=rnof1(iq)+max(wb(iq,k)-ssat(isoil),0.)*1000.*zse(k)
    wb(iq,k)=min(wb(iq,k),ssat(isoil))
  enddo
  dwb=max((wb(iq,ms)-sfc(isoil))*c3(isoil)/86400. , 0.)     ! 23/4/99
  ! for deep runoff use wb-sfc, but this value not to exceed .99*wb-wbice
  dwb=max(min(wb(iq,ms)-sfc(isoil),.99*wb(iq,ms)-wbice(iq,ms))*c3(isoil)/86400. , 0.)     ! 1/9/00
  rnof2(iq)=zse(ms)*1000.*dwb*dt                              ! MJT cable
  wb(iq,ms)=wb(iq,ms)-dwb*dt
  runoff_surface(iq)=runoff_surface(iq)+rnof1(iq)
  runoff(iq)=runoff(iq)+rnof1(iq)+rnof2(iq)  ! accumulated mm ! MJT cable
enddo               ! ip loop

!---  glacier formation
if(nglacier==0)then  ! crashes with tggsn1 going v cold
  do iq=1,ifull
    if(snowd(iq)>400.)then
      rnof5=snowd(iq)-400.
      runoff_surface(iq)=runoff_surface(iq)+rnof5
      runoff(iq)=runoff(iq)+rnof5
!----      change local tg to account for energy - clearly not best method
      sgamm   = ssdn(iq,1)*2105. * sdepth(iq,1)
      tggsn(iq,1)=tggsn(iq,1)-rnof5*hlf/sgamm
      snowd(iq)=400.
    endif ! (snowd(iq)>400.)
  enddo  ! iq loop
endif    ! (nglacier==0)
if(nglacier==2)then  ! new from Eva 3/10/02           
  do iq=1,ifull
    if(snowd(iq)>400.)then
      rnof5=snowd(iq)-400.
      runoff_surface(iq)=runoff_surface(iq)+rnof5
      runoff(iq)=runoff(iq)+rnof5
!----      change local tg to account for energy - clearly not best method
      if(isflag(iq)==0)then
        tgg(iq,1)= tgg(iq,1)-rnof5*hlf/gammzz(iq,1)
        snowd(iq)=400.
      else
        smasstot=smass(iq,1)+smass(iq,2)+smass(iq,3)
        do k=1,3
          sgamm   = ssdn(iq,k)*2105. * sdepth(iq,k)
          smelt1(k)= min(rnof5*smass(iq,k)/smasstot,0.9*smass(iq,k))
          smass(iq,k) = smass(iq,k) - smelt1(k)
          snowd(iq)=snowd(iq)-smelt1(k)
          tggsn(iq,k)= tggsn(iq,k)-smelt1(k)*hlf/sgamm
          snowmelt(iq) = snowmelt(iq) + smelt1(k)
        enddo
      endif ! (isflag(iq)==0) ... else ...
    endif   ! (snowd(iq)>400.)
  enddo    ! iq loop
endif     ! (nglacier==2)

if((ntest>0.or.diag).and.mydiag)then
  if(land(idjd))then !MJT bugfix
    iq=idjd
    write(6,*) 'end surfbv  rnof1,runoff ',rnof1(idjd),runoff(idjd)
    sgamm   = ssdn(iq,1)*2105. * sdepth(iq,1)
    write(6,*) 'snowd,isflag,sgamm ',snowd(idjd),isflag(idjd),sgamm
    write(6,*) 'tggsn_d ',(tggsn(idjd,k),k=1,3)
    write(6,*) 'tgg ',(tgg(idjd,k),k=1,ms)
    write(6,*) 'wb ',(wb(idjd,k),k=1,ms)
  endif
endif
return
end subroutine surfbv

!***********************************************************************

subroutine smoisturev
 
use cc_mpi, only : myid,mydiag
use newmpar_m
use nsibd_m
use parm_m
use permsurf_m
use soil_m           ! land
use soilsnow_m
use soilv_m
use work3_m
use work3b_m
 
implicit none
 
integer, parameter :: ntest=0  ! 2 for funny pre-set for idjd
integer, parameter :: nmeth=-1 ! 1 for full implicit, 2 for simpler implicit
!                            3 for simple implicit D, explicit K jlm pref
!                            4 for simple implicit D, implicit K  
!                            0 for simple implicit D, new jlm TVD K  
!                           -1 for simple implicit D, new jlm TVD K constrained 
      
!
!     solves implicit soil moisture equation
!
!     fwtop  - water flux into the surface (precip-evap)
!     dt   - time step
!     isoil - soil type
!
integer num,k,iq,ip,isoil,iqmx,iqmn,ml 
real, dimension(ifull,ms+3) :: at,bt,ct
real, dimension(ms+1) :: wbh,z1,z2,z3
real, dimension(ifull,0:ms) :: fluxh,delt
real, dimension(ifull,ms) :: dtt
real, dimension(mxst), save :: pwb_min
real, dimension(ms) :: ssatcurr
real, dimension(ms+1) :: z1mult
real totwba,wblfmx,wblfmn
real wbl_kp,wh,wbl_k
real hydss,speed_k,rat,phi
real fluxhi,fluxlo,ssatcurr_k
real wbh_k,fact,wbicefrac
real hsbhh,pwb_wbh,z3_k,totwbb
real totwblb,wbficemx,fact2
real pwb,totwbc,totwblc
real zsetot
real rhowat
data rhowat /1000./

isoil = 0 ! for cray compiler

if(ktau==1)then
  num=1
  do isoil=1,mxst
    pwb_min(isoil)=(swilt(isoil)/ssat(isoil))**ibp2(isoil)
  enddo
  if (myid==0) then
    write(6,*) 'in smoisturev; nmeth,ntest = ',nmeth,ntest  
  end if
endif  ! (ktau==1)
if((ntest>0.or.diag).and.mydiag)then
  if(land(idjd))then !MJT bugfix
    isoil=isoilm(idjd)
    write(6,*) 'entering smoisturev i2bp3,swilt,sfc,ssat: ',i2bp3(isoil),swilt(isoil),sfc(isoil),ssat(isoil)
    if(ntest==2)then   ! just to test conservation
      if(ktau==1)wb(idjd,ms)=swilt(isoil)
      fwtop(idjd)=0.
    endif 
    write (6,"('wb   ',6f8.3)") (wb(idjd,k),k=1,ms)
    write (6,"('wbice',6f8.3)") (wbice(idjd,k),k=1,ms)
    totwba=0.
    do k=1,ms
      totwba=totwba+zse(k)*wb(idjd,k)      ! diagnostic
    enddo
  endif
endif

do k=1,ms ! preset to allow for non-land & snow points in trimb
  at(:,k)=0.
  bt(:,k)=1.
  ct(:,k)=0.
enddo
z1mult(1)=0.      ! corresponds to 2b+3
z1mult(ms+1)=0.   ! corresponds to 2b+3
z1(1)=0.      !  i.e. K(.5),    value at surface
z1(ms+1)=0.   !  i.e. K(ms+.5), value at bottom

wblfmx=0.
wblfmn=1.

if(nmeth<=0)then    ! ip loop split March '03
!      jlm split TVD version
  do ip=1,ipland  ! all land points 
    iq=iperm(ip)
    delt(iq,0)=0.
    fluxh(iq,0)=0.
    fluxh(iq,ms)=0.
  enddo   ! ip loop
  do k=1,ms-1
    do ip=1,ipland  ! all land points 
      iq=iperm(ip)
      isoil = isoilm(iq)
      wbl_k=wb(iq,k)-wbice(iq,k)      ! for calc. speed etc
      wbl_kp=wb(iq,k+1)-wbice(iq,k+1)    
      delt(iq,k)=wbl_kp-wbl_k         
!     wh=(zse(k+1)*wbl(k)+zse(k)*wbl(k+1))/(zse(k)+zse(k+1))
!     especially to allow for isolated frozen layers, use min speed
      wh=min(wbl_k,wbl_kp)
!     with 50% wbice, reduce hyds by 1.e-5
      hydss=hyds(isoil)*(1.-min(2.*wbice(iq,k)/wb(iq,k),.99999))
      speed_k=hydss*(wh/ssat(isoil))**(i2bp3(isoil)-1)
!     update wb by TVD method
      rat=delt(iq,k-1)/(delt(iq,k)+sign(1.e-20,delt(iq,k)))
      phi=max(0.,min(1.,2.*rat),min(2.,rat)) ! 0 for -ve rat
      fluxhi=wh
      fluxlo=wbl_k
      if(ntest>0.and.iq==idjd.and.mydiag)then
        write(6,*) 'in TVD for k= ',k
        write(6,*) 'wbl,wh,hydss ',wbl_k,wh,hydss
        write(6,*) 'speeda,speedb,fluxhi,fluxlo,delt,rat,phi ',speed_k,.5*zse(k)/dt,fluxhi,fluxlo,delt(iq,k),rat,phi
      endif
!        scale speed to grid lengths per dt & limit speed for stability
      speed_k=min(speed_k,.5*zse(k)/dt)  !  1. OK too for stability
      fluxh(iq,k)=speed_k*(fluxlo+phi*(fluxhi-fluxlo))
    enddo    ! ip loop
  enddo  ! k loop

!      update wb by TVD method
  do k=ms,1,-1
    do ip=1,ipland  ! all land points 
      iq=iperm(ip)
      isoil = isoilm(iq)
      if(nmeth==-1)then  ! each new wb constrained by ssat
        fluxh(iq,k-1)=min(fluxh(iq,k-1),(ssat(isoil)-wb(iq,k))*zse(k)/dt +fluxh(iq,k))
      endif   ! (nmeth==-1)
      wb(iq,k)=wb(iq,k)+dt*(fluxh(iq,k-1)-fluxh(iq,k))/zse(k)
!        re-calculate wblf
      ssatcurr_k=ssat(isoil)-wbice(iq,k)
      dtt(iq,k)=dt/(zse(k)*ssatcurr_k)
!        this defn of wblf has different meaning from previous one in surfbv
!        N.B. are imposing wbice<wb, so wblf <1
      wblf(iq,k)=(wb(iq,k)-wbice(iq,k))/ssatcurr_k
    enddo   ! ip loop
  enddo    ! k=ms,1,-1 loop

  do k=2,ms ! wbh_k represents wblf(k-.5)
    do ip=1,ipland  ! all land points in this nsib=1 or 3 loop
      iq=iperm(ip)
      isoil = isoilm(iq)
      ssatcurr_k=ssat(isoil)-wbice(iq,k)
      wbh_k=(zse(k)*wblf(iq,k-1)+zse(k-1)*wblf(iq,k))/(zse(k)+zse(k-1))
!        wbh_k=min(wblf(iq,k-1),wblf(iq,k)) ! jlm to avoid wbice problems
      fact=wbh_k**(ibp2(isoil)-1)   ! i.e. wbh**(bch+1)
!        with 50% wbice, reduce hbsh by 1.e-5
      wbicefrac=max(wbice(iq,k-1)/wb(iq,k-1),wbice(iq,k)/wb(iq,k))
      hsbhh=hsbh(isoil)*(1.-min(2.*wbicefrac,.99999))
      pwb_wbh=hsbhh*max( pwb_min(isoil),wbh_k*fact )
!        moisture diffusivity (D) is  wbh*pwb; hsbh includes b
      z3_k=pwb_wbh/zshh(k)            !  i.e. D(k-.5)/zshh(k)
      at(iq,k) = -dtt(iq,k)*z3_k  ! where dtt=dt/(zse(k)*ssatcurr_k)
      ct(iq,k-1) = -dtt(iq,k-1)*z3_k
    enddo   ! ip loop
  enddo   ! k loop
      
  do k=1,ms
    do iq=1,ifull  ! can do all points
      bt(iq,k)=1.-at(iq,k)-ct(iq,k)
    enddo   ! iq loop
  enddo    ! k loop

  if((ntest>0.or.diag).and.mydiag)then
    if(land(idjd))then !MJT bugfix
      write(6,*) 'midway through nmeth<=0'
      write(6,*) 'fluxh ',(fluxh(iq,k),k=1,ms)
      write (6,"('wb   ',6f8.3)") (wb(idjd,k),k=1,ms)
      write (6,"('wblf ',6f8.3)") (wblf(idjd,k),k=1,ms)    
      totwbb=0.
      totwblb=0.
      do k=1,ms
        totwbb=totwbb+zse(k)*wb(idjd,k)         ! diagnostic
        totwblb=totwblb+zse(k)*wblf(idjd,k)     ! diagnostic
      enddo
      write(6,*) 'nmeth, b+2, 2b+3: ',nmeth,ibp2(isoil),i2bp3(isoil)
      write (6,"('wb   ',6f8.3)") (wb(idjd,k),k=1,ms)
      write (6,"('wbice',6f8.3)") (wbice(idjd,k),k=1,ms)
      write (6,"('wblf ',6f8.3)") (wblf(idjd,k),k=1,ms)
      write(6,*) 'zse ',zse
      write(6,*) 'at ',(at(idjd,k),k=1,ms)
      write(6,*) 'bt ',(bt(idjd,k),k=1,ms)
      write(6,*) 'ct ',(ct(idjd,k),k=1,ms)
    endif
  endif  ! (ntest>0)

  do ip=1,ipland  ! all land points 
    iq=iperm(ip)
    isoil = isoilm(iq)
    wblf(iq,1)=wblf(iq,1)+dtt(iq,1)*fwtop(iq)/rhowat
  enddo   ! ip loop
endif    ! (nmeth<=0)

if(nmeth>0)then    ! ip loop split March '03
  do ip=1,ipland  ! all land points in this nsib=1 or 3 loop
    iq=iperm(ip)
    isoil = isoilm(iq)
    wbficemx=0.
    do k=1,ms
      ssatcurr(k)=ssat(isoil)-wbice(iq,k)
!        this defn of wblf has different meaning from previous one in surfbv
!        N.B. are imposing wbice<wb, so wblf <1
      wblf(iq,k)=(wb(iq,k)-wbice(iq,k))/ssatcurr(k)
      wbfice(iq,k)=wbice(iq,k)/ssat(isoil)
      wbficemx=max(wbficemx,wbfice(iq,k))
      dtt(iq,k)=dt/(zse(k)*ssatcurr(k))
    enddo

    if(nmeth==1)then  ! full implicit method
      do k=2,ms
        wbh(k)=(zse(k)*wblf(iq,k-1)+zse(k-1)*wblf(iq,k))/(zse(k)+zse(k-1))
        fact=wbh(k)**(ibp2(isoil)-1)   ! i.e. wbh**(bch+1)
        fact2=fact*fact
        pwb = hsbh(isoil)*fact
!          moisture diffusivity (D) is  wbh*pwb
!          other term (K) is wbh*hyds(isoil)*fact2
        z1(k)=wbh(k)*( (i2bp3(isoil)-1)*hyds(isoil)*fact2-ibp2(isoil)*pwb*(wblf(iq,k)-wblf(iq,k-1))/zshh(k) )
        z2(k)=-i2bp3(isoil)*hyds(isoil)*fact2 + ibp2(isoil)*pwb*(wblf(iq,k)-wblf(iq,k-1))/zshh(k)
        z3(k)=pwb*wbh(k)/zshh(k)
        at(iq,k) = dtt(iq,k)*( z2(k)*.5*zse(k)/zshh(k) -z3(k) )
      enddo
      do k=1,ms-1
        ml = max (k-1,1)
        ct(iq,k)=dtt(iq,k)*( -z2(k+1)*.5*zse(k)/zshh(k+1) -z3(k+1))
        bt(iq,k)=1.+dtt(iq,k)*( -z2(k+1)*.5*zse(k+1)/zshh(k+1) + z2(k)*.5*zse(ml)/zshh(k)+z3(k+1)+z3(k) )
      enddo
      bt(iq,ms)=1.+dtt(iq,ms)*( z2(ms)*.5*zse(ms)/zshh(ms) +z3(ms) )
      do k=1,ms
        wblf(iq,k) = wblf(iq,k)+dtt(iq,k)*( z1(k+1) - z1(k) )
      enddo
    endif   ! (nmeth==1)  ! full implicit method

    if(nmeth>=2)then  ! part implicit method
      do k=2,ms
        z1mult(k)=i2bp3(isoil)   ! corresponds to 2b+3
      enddo
      do k=2,ms ! wbh(k) represents wblf(k-.5)
        wbh(k)=(zse(k)*wblf(iq,k-1)+zse(k-1)*wblf(iq,k))/(zse(k)+zse(k-1))
        fact=wbh(k)**(ibp2(isoil)-1)   ! i.e. wbh**(bch+1)
        if(nmeth==2)pwb_wbh=hsbh(isoil)*wbh(k)*fact
        if(nmeth>=3)pwb_wbh=hsbh(isoil)*max( pwb_min(isoil),wbh(k)*fact )
        fact2=fact*fact
!          moisture diffusivity (D) is  wbh*pwb
!          other term (K) is wbh*hyds(isoil)*fact2
        z1(k)=hyds(isoil)*fact2     !  i.e. K(k-.5)/wbh(k)
        z3(k)=pwb_wbh/zshh(k)            !  i.e. D(k-.5)/zshh(k)
        at(iq,k) = -dtt(iq,k)*z3(k)
        ct(iq,k-1) = -dtt(iq,k-1)*z3(k)
      enddo
      do k=1,ms
        bt(iq,k)=1.-at(iq,k)-ct(iq,k)
      enddo
      if(nmeth==4)then   ! for simple implicit D, implicit K
        bt(iq,1)=bt(iq,1)+dtt(iq,1)*z1mult(1+1)*z1(1+1)*zse(1+1)/(zse(1)+zse(1+1))
        do k=2,ms
          at(iq,k)=at(iq,k)-dtt(iq,k)*z1mult(k)*z1(k)*zse(k)/(zse(k)+zse(k-1))
          ct(iq,k-1)=ct(iq,k-1)+dtt(iq,k-1)*z1mult(k)*z1(k)*zse(k-1)/(zse(k)+zse(k-1))
          bt(iq,k)=bt(iq,k)-dtt(iq,k)*z1mult(k)*z1(k)*zse(k-1)/(zse(k)+zse(k-1)) &
                  +dtt(iq,k)*z1mult(k+1)*z1(k+1)*zse(k+1)/(zse(k)+zse(k+1))
        enddo
      endif  ! (nmeth==4)
      do k=2,ms
        z1(k)=wbh(k)*z1(k)     !  i.e. now K(k-.5)
      enddo
!         the following top & bottom b.c.'s will preserve a uniform column
      z1(1) =min(z1(2),z1(ms))   !  N.B. z1 are here +ve
      z1(ms+1)=z1(1)
      if(wbficemx<.75)then ! no gravit. term if too much ice 11/12/00
        do k=1,ms
          if(nmeth==4)then
            wblf(iq,k)=wblf(iq,k)+dtt(iq,k)*((z1mult(k+1)-1.)*z1(k+1) - (z1mult(k)-1.)*z1(k) )
          else
            wblf(iq,k) = wblf(iq,k)+dtt(iq,k)*( z1(k) - z1(k+1) )
          endif  ! (nmeth==4) .. else ..
        enddo
      endif  ! (wbficemx<.75)
    endif    ! (nmeth>=2)

    if(ntest>0)then
      do k=1,ms
        if(wblf(iq,k)>wblfmx)then
          wblfmx=wblf(iq,k)
          iqmx=iq
        endif
        if(wblf(iq,k)<wblfmn)then
          wblfmn=wblf(iq,k)
          iqmn=iq
        endif
      enddo
    endif
    if(ntest>0.and.iq==idjd.and.mydiag)then
      totwbb=0.
      totwblb=0.
      do k=1,ms
        totwbb=totwbb+zse(k)*wb(iq,k)         ! diagnostic
        totwblb=totwblb+zse(k)*wblf(iq,k)     ! diagnostic
      enddo
      write(6,*) 'nmeth, b+2, 2b+3: ',nmeth,ibp2(isoil),i2bp3(isoil)
      write (6,"('wb   ',6f8.3)") (wb(idjd,k),k=1,ms)
      write (6,"('wbice',6f8.3)") (wbice(idjd,k),k=1,ms)
      write (6,"('wblf ',6f8.3)") (wblf(idjd,k),k=1,ms)
      write (6,"('wbh  ',7f8.3)") wbh
      write (6,"('ssatcurr',6f8.3)") ssatcurr
      write(6,*) 'pwb_wbh,pwb_min* for ms ',pwb_wbh,hsbh(isoil)*pwb_min(isoil)
      write(6,*) 'wblfmx,wblfmn,iqmx,iqmn ',wblfmx,wblfmn,iqmx,iqmn
      write(6,*) 'zse ',zse
      write(6,*) 'zshh ',zshh
      write(6,*) 'at ',(at(iq,k),k=1,ms)
      write(6,*) 'bt ',(bt(iq,k),k=1,ms)
      write(6,*) 'ct ',(ct(iq,k),k=1,ms)
    endif  ! (ntest>0.and.iq==idjd)

    if(nmeth==3)then
!     artificial fix applied here for safety (explicit nmeth only)
      do k=1,ms
        wblf(iq,k)=max(0.,min(wblf(iq,k),1.))
      enddo
    endif   ! (nmeth==3)

    wblf(iq,1)=wblf(iq,1)+dtt(iq,1)*fwtop(iq)/rhowat
  enddo   ! ip loop
endif   ! (nmeth>0)

call trimb(at,bt,ct,wblf,ms)               ! B

!     *** following loop needed sopt on SX5 212 compiler!!
do ip=1,ipland  ! all land points in this nsib=1  or 3 loop
  iq=iperm(ip)
  isoil = isoilm(iq)
  do k=1,ms
    ssatcurr(k)=ssat(isoil)-wbice(iq,k)
    wb(iq,k)=wblf(iq,k)*ssatcurr(k)+wbice(iq,k)
    wbice(iq,k)=min(wbice(iq,k),.99*wb(iq,k))
  enddo
enddo
if((ntest>0.or.diag).and.mydiag)then
  if(land(idjd))then !MJT bugfix
    write(6,*) 'at end of smoisturev,fwtop ',fwtop(idjd)
    write(6,*) 'tgg ',(tgg(idjd,k),k=1,ms)
    write (6,"('wb   ',6f8.3)") (wb(idjd,k),k=1,ms)
    write (6,"('wbice',6f8.3)") (wbice(idjd,k),k=1,ms)
    write (6,"('wblf ',6f8.3)") (wblf(idjd,k),k=1,ms)
    totwbc=0.
    totwblc=0.
    zsetot=0.
    do k=1,ms
      totwbc=totwbc+zse(k)*wb(idjd,k)       ! diagnostic
      totwblc=totwblc+zse(k)*wblf(idjd,k)   ! diagnostic
      zsetot=zsetot+zse(k)
    enddo
    write(6,*) 'totwba,totwbb,totwbc ',totwba,totwbb,totwbc
    write(6,*) 'totwblb,totwblc ',totwblb,totwblc
    write(6,*) 'with totwbc/zsetot: ',totwbc/zsetot
  endif
endif

return
end subroutine smoisturev

!***********************************************************************

subroutine stempv(gammzz)

use cc_mpi, only : mydiag
use newmpar_m
use nsibd_m
use parm_m
use permsurf_m
use soil_m     ! land
use soilsnow_m
use soilv_m
use work2_m
use work3_m
use work3b_m

implicit none

integer, parameter :: ntest=0

!     calculates temperatures of the soil 
!     tgsoil - new soil/ice temperature
!     ga - heat flux from the atmosphere (ground heat flux)
!     ccnsw - soil conductivity
!     dt  - time step 
      
integer k,ip,iq,isoil
real, dimension(ifull,-2:ms) :: at,bt,ct
real, dimension(ifull,-2:ms) :: tggdm
real, dimension(ifull) :: coefa,coefb
real, dimension(ifull,ms) :: gammzz
real, dimension(ifull,ms) :: ccnsw
real, dimension(-2:ms) :: rhs
real, dimension(-2:ms+1) :: coeff
real, dimension(3) :: sconds
real csice,cswat,rhowat,cgsnow,rhosnow
real ew,ccf,ei,scondss,xx,dtg,sgamm
real xy
      
data csice /2.100e3/, cswat /4.218e3/, rhowat /1000./  ! for calgammv
data cgsnow/2090./,rhosnow/200./                       ! for calgammv

do k=-2,ms ! preset to allow for non-land & snow points in trimb
  at(:,k)=0.
  bt(:,k)=1.
  ct(:,k)=0.
enddo

do k=1,ms
  do ip=1,ipland  ! all land points in this nsib=1 or 3 loop
    iq=iperm(ip)
    isoil = isoilm(iq)
    if(isoil==9)then
      ccnsw(iq,k)=2.5
    else
      ew=wblf(iq,k)*ssat(isoil)
      ccf=sqrt(.5/max(.25,min(wblf(iq,k),.5)))  ! jlm: equiv to above
      ei=wbfice(iq,k)*ssat(isoil)
      ccnsw(iq,k)= min(cnsd(isoil)*exp(ew*log(60.)+ei*log(250.)), 2.2)*ccf
    endif  ! (isoil==9) ... else ...
  enddo
enddo

coeff(1)=0.
coeff(ms+1)=0.
do ip=1,ipland  ! all land points in this nsib=1 or 3 loop
  iq=iperm(ip)
  if( isflag(iq)==0) then
    isoil = isoilm(iq)
    scondss = max(0.2,min(2.576e-6*ssdn(iq,1)*ssdn(iq,1)+.074,1.))   ! or should it be 0.8?
    xx=max(0.,snowd(iq)/ssdnn(iq))
    xy=zse(1)/(zse(1)+xx)
    if( xx >0.) ccnsw(iq,1)=ccnsw(iq,1)*xy + (1.-xy)*scondss
    coeff(2)=2./((zse(1)+xx)/ccnsw(iq,1)+zse(2)/ccnsw(iq,2))
    coefa(iq)=0.         ! jlm for B
    coefb(iq)=coeff(2)   ! jlm for B
    do k=3,ms
      coeff(k)=2./(zse(k-1)/ccnsw(iq,k-1)+zse(k)/ccnsw(iq,k))
    enddo

    k=1
    gammzz(iq,k)=max( (1.-ssat(isoil))*css(isoil)*rhos(isoil)+ssat(isoil)*(wblf(iq,k)*cswat*rhowat+wbfice(iq,k)*csice*rhowat*.9), &
                      css(isoil)*rhos(isoil)  ) * zse(k)
    gammzz(iq,k)=gammzz(iq,k) +  cgsnow*snowd(iq) 
    dtg=dt/gammzz(iq,k)
    at(iq,k)= -dtg*coeff(k)
    ct(iq,k)= -dtg*coeff(k+1)     ! c3(ms)=0 & not really used
    bt(iq,k)= 1.-at(iq,k)-ct(iq,k)

    do k=2,ms
      gammzz(iq,k)=max((1.-ssat(isoil))*css(isoil)*rhos(isoil)+ssat(isoil)*(wblf(iq,k)*cswat*rhowat+wbfice(iq,k)*csice*rhowat*.9), &
                       css(isoil)*rhos(isoil)  ) * zse(k)
      dtg=dt/gammzz(iq,k)
      at(iq,k)= -dtg*coeff(k)
      ct(iq,k)= -dtg*coeff(k+1)     ! c3(ms)=0 & not really used
      bt(iq,k)= 1.-at(iq,k)-ct(iq,k)
    enddo   ! k=1,ms loop

    bt(iq,1)=bt(iq,1)-dgdtg(iq)*dt/gammzz(iq,1)             ! 9/3/99
    tgg(iq,1)=tgg(iq,1)+(ga(iq)-tgg(iq,1)*dgdtg(iq))*dt/gammzz(iq,1)  ! 9/3/99
  endif  ! ( isflag(iq)== 0 )
enddo   ! ip=1,ipland           land points
if(ntest>0.and.mydiag)then
  iq=idjd
  isoil=isoilm(iq)
  write(6,*) 'ga,gammzz ',ga(iq),gammzz(iq,1)
  write(6,*) 'rhs_tgg ',(tgg(iq,k),k=1,ms)
  write(6,*) 'dgdtg ',dgdtg(iq)
  write(6,*) 'ssat,css,rhos,cswat,rhowat,csice ',ssat(isoil),css(isoil),rhos(isoil),cswat,rhowat,csice
  write(6,*) 'wblf1,wbfice1,zse1,cgsnow ',wblf(iq,1),wbfice(iq,1),zse(1),cgsnow
  write(6,*) 'at ',(at(iq,k),k=1,ms)
  write(6,*) 'bt ',(bt(iq,k),k=1,ms)
  write(6,*) 'ct ',(ct(iq,k),k=1,ms)
endif  ! (ntest>0)

coeff(1-3)=0.
! ****** next cdir nodep does not work on SX5 ****************      

do ip=1,ipland  ! all land points in this nsib=1 or 3 loop
!      N.B. vectorized version assumes tggsn precedes tgg in memory
  iq=iperm(ip)
  if( isflag(iq)/=0) then   ! 3-layer snow points done here
    isoil = isoilm(iq)
    do k=1,3
      sconds(k)=max(.2,min(2.576e-6*ssdn(iq,k)*ssdn(iq,k)+.074,1.))
    enddo
    coeff(1)=2./(sdepth(iq,3)/sconds(3) +zse(1)/ccnsw(iq,1))
    do k=2,ms
      coeff(k)=2./(zse(k-1)/ccnsw(iq,k-1)+zse(k)/ccnsw(iq,k))
    enddo
    k=3
    coeff(k-3)=2./(sdepth(iq,k-1)/sconds(k-1)+sdepth(iq,k)/sconds(k))
    sgamm   = ssdn(iq,k)*2105. * sdepth(iq,k)
    dtg=dt/sgamm
    at(iq,k-3) = -dtg*coeff(k-3)
    ct(iq,k-3) = -dtg*coeff(k-2)
    bt(iq,k-3)= 1.-at(iq,k-3)-ct(iq,k-3)
    k=2  
    coeff(k-3)=2./(sdepth(iq,k-1)/sconds(k-1)+sdepth(iq,k)/sconds(k))
    sgamm   = ssdn(iq,k)*2105. * sdepth(iq,k)
    dtg=dt/sgamm
    at(iq,k-3) = -dtg*coeff(k-3)
    ct(iq,k-3) = -dtg*coeff(k-2)
    bt(iq,k-3)= 1.-at(iq,k-3)-ct(iq,k-3)
    k=1 
    sgamm   = ssdn(iq,k)*2105. * sdepth(iq,k)
    dtg=dt/sgamm
    at(iq,k-3) = -dtg*coeff(k-3)
    ct(iq,k-3) = -dtg*coeff(k-2)
    bt(iq,k-3)= 1.-at(iq,k-3)-ct(iq,k-3)
     
    coefa(iq)=coeff(2-3)     ! jlm B
    coefb(iq)=coeff(4-3)     ! jlm B

    do k=1,ms
      gammzz(iq,k)=max((1.-ssat(isoil))*css(isoil)*rhos(isoil)+ssat(isoil)*(wblf(iq,k)*cswat*rhowat+wbfice(iq,k)*csice*rhowat*.9), &
                       css(isoil)*rhos(isoil)  ) * zse(k)
      dtg=dt/gammzz(iq,k)
      at(iq,k)= -dtg*coeff(k)
      ct(iq,k) = -dtg*coeff(k+1)      ! c3(ms)=0 & not really used
      bt(iq,k)= 1.-at(iq,k)-ct(iq,k)
    enddo

    sgamm   = ssdn(iq,1)*2105. * sdepth(iq,1)
    bt(iq,-2)=bt(iq,-2)-dgdtg(iq)*dt/sgamm             ! 9/5/02
    tggsn(iq,1)=tggsn(iq,1)+(ga(iq)-tggsn(iq,1)*dgdtg(iq))*dt/sgamm ! 9/5/-2

    rhs(1-3)=tggsn(iq,1)    ! A
    if(ntest>0.and.iq==idjd.and.mydiag)then
      write(6,*) 'in stempv 3-layer snow code '
      write(6,*) 'ccnsw ',(ccnsw(iq,k),k=1,ms)
      write(6,*) 'sdepth d ',(sdepth(iq,k),k=1,3)
      write(6,*) 'sconds ',sconds
      write(6,*) 'coeff ',coeff
      write(6,*) 'at ',(at(iq,k),k=-2,ms)
      write(6,*) 'bt ',(bt(iq,k),k=-2,ms)
      write(6,*) 'ct ',(ct(iq,k),k=-2,ms)
      write(6,*) 'rhs:tggsn,tgg ',(tggsn(iq,k),k=1,3),(tgg(iq,k),k=1,ms)
    endif  ! (ntest>0.and.iq==idjd)
  endif  ! ( isflag(iq).ne. 0 )
enddo   ! ip=1,ipland           land points

!     note in the following that tgg and tggsn are stacked together
tggdm(:,-2:0)=tggsn ! MJT bugfix
tggdm(:,1:ms)=tgg   ! MJT bugfix
call trimb(at(1,-2),bt(1,-2),ct(1,-2),tggdm,ms+3)           ! B
tggsn=tggdm(:,-2:0) ! MJT bugfix
tgg=tggdm(:,1:ms)   ! MJT bugfix

do k=2,ms
  do iq=1,ifull
    isoil = isoilm(iq)
    if(isoil==9)then
      tgg(iq,k)=min(tgg(iq,k),273.1)  ! jlm 7/6/00
    endif
  enddo   
enddo

do ip=1,ipland  ! all land points in this nsib=1 or 3 loop
  iq=iperm(ip)
  sgflux(iq)=coefa(iq)*(tggsn(iq,1)-tggsn(iq,2))
  gflux(iq) =coefb(iq)*(  tgg(iq,1)-  tgg(iq,2))  ! +ve downwards
enddo   ! ip=1,ipland           land points
if((ntest>0.or.diag).and.mydiag)then
  if(land(idjd))then !MJT bugfix
    write(6,*) 'at end of stempv '
    write (6,"('wb   ',6f8.3)") (wb(idjd,k),k=1,ms)
    write (6,"('wbice',6f8.3)") (wbice(idjd,k),k=1,ms)
    write (6,"('wblf ',6f8.3)") (wblf(idjd,k),k=1,ms)
    write(6,*) 'tggsn ',(tggsn(idjd,k),k=1,3)
    write(6,*) 'tgg ',(tgg(idjd,k),k=1,ms)
  endif
endif  ! (ntest>0)

return
end subroutine stempv

!***********************************************************************

subroutine snowprv(iq)    ! N.B. this one is not vectorized

use const_phys
use newmpar_m
use nsibd_m
use parm_m
use soil_m  ! land
use soilsnow_m
use soilv_m
use work3b_m

implicit none

integer, parameter :: newsmlt=1 ! 0: old, 1: new from Aug 2003

integer, intent(in) :: iq
integer k
real, dimension(3) :: etac
real sdd,ossdn1,ossdn2,ossdn3
real ccoef,tggd,tr1,xx,pr
real osm1,excm,excd,osm2,osm3
real sd1,sm1,otggsn1

if( isflag(iq)==0) then
  tggsn(iq,1) = tgg(iq,1)
  tggsn(iq,2) = tgg(iq,1)
  tggsn(iq,3) = tgg(iq,1)
  ssdn(iq,2)  = ssdn(iq,1)
  ssdn(iq,3)  = ssdn(iq,1)
  sdepth(iq,1)= .07
  sdd=(snowd(iq)-.07*ssdn(iq,1))/ssdn(iq,1)
  sdepth(iq,2)=max(.02,0.45*sdd)
  sdepth(iq,3)=max(.02,0.55*sdd)
  if(snowd(iq)>20.)sdepth(iq,2)=max(.02,0.3*sdd)
  if(snowd(iq)>20.)sdepth(iq,3)=max(.02,0.7*sdd)
  smass(iq,1) = .07*ssdn(iq,1)
  smass(iq,2) = sdepth(iq,2)*ssdn(iq,2)
  smass(iq,3) = sdepth(iq,3)*ssdn(iq,3)
endif

ossdn1 = ssdn(iq,1)
ossdn2 = ssdn(iq,2)
ossdn3 = ssdn(iq,3)
do k=1,3
  ccoef  = 0.
  if(ssdn(iq,k) >=150.) ccoef=4.6e-2
  tggd   = min(tfrz,tggsn(iq,k))
  ssdn(iq,k) = ssdn(iq,k)+dt*ssdn(iq,k)*3.1e-6*exp(-.03*(273.1-tggd)-ccoef*(ssdn(iq,k)-150.))
  etac(k)=3.e7*exp( .021*ssdn(iq,k)+.081*(273.1-tggd) )  ! same as:
enddo
ssdn(iq,1)=ssdn(iq,1)+dt*grav*.5*.07*ssdn(iq,1)*ssdn(iq,1)/etac(1)
ssdn(iq,2)=ssdn(iq,2)+dt*grav*ssdn(iq,2)*(.07*ssdn(iq,1)+.5*smass(iq,2))/etac(2)
ssdn(iq,3)=ssdn(iq,3)+dt*grav*ssdn(iq,3)*(.07*ssdn(iq,1)+smass(iq,2)+.5*smass(iq,3))/etac(3)
 
tr1  =  snowd(iq)-osnowd(iq)
xx=max(0.,.07-smass(iq,1)/ssdn(iq,1))
pr=min(smass(iq,2)/(smass(iq,3)+smass(iq,2)),.9)

if( tr1>=0.) then
  tr1        =  tr1/140.
  ossdn1     = ssdn(iq,1)
  ssdn(iq,1)=max((smass(iq,1)+tr1*140.)/(smass(iq,1)/ossdn1+tr1),140.)
  osm1        = smass(iq,1)
  smass(iq,1) = .07*ssdn(iq,1)
  sdepth(iq,1)= .07
  excm        = osm1+tr1*140.-smass(iq,1)
  excd        = excm/ssdn(iq,1)

  osm2        = smass(iq,2)
  ossdn2      = ssdn(iq,2)
  smass(iq,2) = max(.01,smass(iq,2)+0.4*excm)
  ssdn(iq,2)=max(140.,min(500.,smass(iq,2)/(osm2/ossdn2+.4*excd)))
  sdepth(iq,2)=max(.02,smass(iq,2)/ssdn(iq,2))

  osm3        = smass(iq,3)
  smass(iq,3) = max(.01,snowd(iq)-smass(iq,1)-smass(iq,2))
  sdepth(iq,3) = max(.02,osm3/ssdn(iq,3)+0.6*excm/ssdn(iq,2))
  ssdn(iq,3)  = max(140.,min(500.,smass(iq,3)/sdepth(iq,3)))
  if(ssdn(iq,3)<ssdn(iq,2)) then
    ssdn(iq,3) = ssdn(iq,2)
    sdepth(iq,3) =max(.02,smass(iq,3)/ssdn(iq,3))
  endif

else
  ! snow melting
  sdepth(iq,1)   = .07
  sd1= max(.005,smass(iq,1)/ssdn(iq,1)) !current depth of 1st layer
  sm1= max(.01,smass(iq,1))             !current mass of 1st layer
  excd        = .07-sd1
  smass(iq,1)=max(140.*.07,min(500. , sd1*ssdn(iq,1)+excd*ssdn(iq,2) ) )
  ssdn(iq,1)=smass(iq,1)/.07
  if(newsmlt==0)then  ! old way
    excm        = smass(iq,1)-sm1
    excd        = excm/ssdn(iq,1)
    osm2        = smass(iq,2)
    smass(iq,2) = max(.01,smass(iq,2)-pr*excm)
    sdepth(iq,2)= max(.02,osm2/ssdn(iq,2)-pr*excd)
    ssdn(iq,2)  = max(140.,min(500.,smass(iq,2)/sdepth(iq,2)))
    if( ssdn(iq,2) < ossdn2 ) then
      ssdn(iq,2)=ossdn2
      smass(iq,2)=.45*(snowd(iq)-smass(iq,1))
      sdepth(iq,2)=max(.02,smass(iq,2)/ssdn(iq,2))
    endif
    smass(iq,3)  = max(.01 , snowd(iq)-smass(iq,1)-smass(iq,2))
    sdepth(iq,3) = max(.02 , smass(iq,3)/ssdn(iq,3))
  endif  ! (newsmlt==0)
  if(newsmlt==1)then  ! new way
    excm       = max(0.,smass(iq,1)-sm1)
    otggsn1    = tggsn(iq,1)
    tggsn(iq,1)= tggsn(iq,1)*sm1/smass(iq,1) + (1.-sm1/smass(iq,1))*tggsn(iq,2)
    smass(iq,2) = max(.01,smass(iq,2)-pr*excm)
    sdepth(iq,2)= max(.02,smass(iq,2)/ssdn(iq,2))
    smass(iq,3) = max(.01 , snowd(iq)-smass(iq,1)-smass(iq,2))
    sdepth(iq,3)= max(.02 , smass(iq,3)/ssdn(iq,3))
    if(sdepth(iq,3)<sdepth(iq,2)) then
      smass(iq,2)=.45*(snowd(iq)-smass(iq,1))
      sdepth(iq,2)=max(.02,smass(iq,2)/ssdn(iq,2))
      smass(iq,3) = max(.01 , snowd(iq)-smass(iq,1)-smass(iq,2))
      sdepth(iq,3)= max(.02 , smass(iq,3)/ssdn(iq,3))
    endif  ! (sdepth(iq,3)<sdepth(iq,2))
  endif    ! (newsmlt==1)
endif      ! ( tr1>=0.) ... else ...

return
end subroutine snowprv

!***********************************************************************

subroutine trimb(a,b,c,rhs,kmax)
!     like trim, but work arrays in work3f
!     rhs initially contains rhs; leaves with answer (jlm)
!     n.b. this one does not assume b = 1-a-c

use newmpar_m

implicit none

integer, intent(in) :: kmax
integer iq,k
real, dimension(ifull,kl) :: a,b,c
real, dimension(ifull,kl) :: e,g,temp,rhs

!     this routine solves the system
!       a(k)*u(k-1)+b(k)*u(k)+c(k)*u(k+1)=rhs(k)    for k=2,kmax-1
!       with  b(k)*u(k)+c(k)*u(k+1)=rhs(k)          for k=1
!       and   a(k)*u(k-1)+b(k)*u(k)=rhs(k)          for k=kmax

!     the Thomas algorithm is used

do iq=1,ifull
  e(iq,1)=c(iq,1)/b(iq,1)
enddo
do k=2,kmax-1
  do iq=1,ifull
    temp(iq,k)= 1./(b(iq,k)-a(iq,k)*e(iq,k-1))
    e(iq,k)=c(iq,k)*temp(iq,k)
  enddo ! iq loop
enddo  ! k loop

do iq=1,ifull
  g(iq,1)=rhs(iq,1)/b(iq,1)
enddo
do k=2,kmax-1
  do iq=1,ifull
    g(iq,k)=(rhs(iq,k)-a(iq,k)*g(iq,k-1))*temp(iq,k)
  enddo ! iq loop
enddo  ! k loop

!     do back substitution to give answer now
do iq=1,ifull
  rhs(iq,kmax)=(rhs(iq,kmax)-a(iq,kmax)*g(iq,kmax-1))/(b(iq,kmax)-a(iq,kmax)*e(iq,kmax-1))
enddo
do k=kmax-1,1,-1
  do iq=1,ifull
    rhs(iq,k)=g(iq,k)-e(iq,k)*rhs(iq,k+1)
  enddo ! iq loop
enddo  ! k loop

return
end subroutine trimb
