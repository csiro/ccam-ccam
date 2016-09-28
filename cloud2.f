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

c This is the interface between the Fels-Schwarzkopf radiation scheme and
c LDR's prognostic cloud water scheme. It is called by radrive.
c
c INPUT/OUTPUT:
c Input:
c
c      dprf - pressure thickness at each sigma level
c      prf  - pressure at full levels
c      coszro - zenith angle at grid pt. (used in SW scheme)
c      cldoff - flag set to T for clear sky radiation calculations
c      lg     - latitude index
c      ttg    - temperature
c      qlg    - cloud liquid water mixing ratio (kg/kg)
c      qfg    - cloud ice mixing ratio (kg/kg)
c      cfrac  - total cloud fraction (stratiform + convective)
c      qccon  - cloud water mixing ratio of convective clouds (kg/kg)
c
c Output:
c
c in common/radisw in RADISW.f (passed back to radiation scheme)
c      camt   - cloud amounts (locations specified by ktop/kbtm indices) 
c      cirab  - absorptivity of clouds in the near IR band (used in SW scheme)
c      cirrf  - reflectivity of clouds in the near IR band (used in SW scheme)
c      cuvrf - reflectivity of clouds in the visible band (used in SW scheme)
c      emcld - cloud emissivity (used in LW scheme)
c      kbtm  - index of (data level) pressure of cloud bottom (used in LW)
c      kbtmsw- index of (flux level) pressure of cloud bottom (used in SW)
c      ktop  - index of (data level) pressure of cloud top (used in LW)
c      ktopsw- index of (flux level) pressure of cloud top (used in SW)
c      nclds - no. clouds at each grid point
c
c in arguments
c      clat - cloud amount diagnostic (upside-down version of cfrac array)
c      clh - high level cloud diagnostic
c      cll - low level cloud diagnostic
c      clm - mid level cloud diagnostic
c
c******************************************************************************
 
      subroutine cloud2(cldoff,lg,ttg,qlg,qfg,cfrac,qccon,
     &                  cdrop,land,sigh,prf,dprf,cosz,     !Inputs
     &                  cll,clm,clh)                       !Outputs

      use diag_m
      use cc_mpi, only : mydiag,myid
      use const_phys
      use leoncld_mod
      use newmpar_m
      use parm_m
      use radisw_m  !Output various things (see above) and input coszro
      use sigs_m
      implicit none
C Global parameters
      include 'kuocom.h'     ! ldr
      include 'rdparm.h'     !Input radiation scheme parameters
      include 'hcon.h'       !Input radiation physical constants

C Argument list
      logical cldoff
      integer lg
      real ttg(imax,l)
      real qlg(imax,l)
      real qfg(imax,l)
      real cfrac(imax,l)
      real qccon(imax,l)
      logical land(imax)
      real sigh(kl+1)
      real prf(imax,l)
      real dprf(imax,l)
      real cosz(imax)
      real cll(imax)
      real clm(imax)
      real clh(imax)
      real mx(imax)

c      real Refflm(imax)
c      real cldliq(imax)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      integer naerosol_i(2)
      data naerosol_i / 2*0 /

      integer jclb,nlow,nmid,ktied,kbgel,klcmc,klowrn,kmidrn
      real aftsea
      common/levdata/jclb,nlow,nmid,aftsea,ktied,kbgel
     &,klcmc,klowrn,kmidrn

C Local work arrays and variables
      real taul(imax,l), taui(imax,l) !Visible optical depths
      real Reffl(imax,l), Reffi(imax,l) !Effective radii (um)
      real Emw(imax,l), Abw(imax,l), Rew1(imax,l), Rew2(imax,l) !Water clouds
      real Emi(imax,l), Abi(imax,l), Rei1(imax,l), Rei2(imax,l) !Ice clouds
      real qlptot(imax),taultot(imax)
      real fice(imax,l),tau_sfac(imax,l)
      real rk(imax,kl),cdrop(imax,kl)

      integer k
      integer mg
      integer nc

      integer kb,km,kk,kt,kp,pos(1)
      integer, dimension(imax,l) :: cldtop,cldbtm

      real ab, cfl, cldht, deltai, diffk, dz, em, eps, f1, f2, fcon
      real onem, qlpath, refac, refac1, refac2, re1, re2, rhoa
      real siglow, sigmid
      real sigmai, tciwc, tclwc, temp_correction, tmid
      real trani, tranw, wice, wliq

      real csum,ctemp,ttgsum,qlsum,qfsum,qconsum
      real refflsum,reffisum,dzsum,qlpathsum,qipathsum
      real refflhold,reffihold,ficehold

      integer, save :: istart
      data istart/0/
C Start code : ----------------------------------------------------------

c**** Set up levels for low/mid/high cloud counting (for cloud2.f)

      if(istart==0)then

c..   The FULL level is 1 below the half level indicator
c     N.B. for ISCCP, low is at 680 mb and mid at 440 mb
c.... "low" cloud up to half level closest to p/Ps=.68 (~800mbs)
c.... "mid" cloud up to half level closest to p/Ps=.44 (~400mbs)
c       siglow=.8
c       sigmid=.4
        siglow=.68
        sigmid=.44
        f1=1.
        f2=1.
        do k=1,kl-1
         if(abs(sigh(k+1)-siglow)<f1)then
           f1=abs(sigh(k+1)-siglow)
           nlow=k
         endif
         if(abs(sigh(k+1)-sigmid)<f2)then
           f2=abs(sigh(k+1)-sigmid)
           nmid=k
         endif
        enddo
        if(mydiag)then
         print *,'in cloud2, nlow,nmid,sigh_vals: ',
     &                       nlow,nmid,sigh(nlow+1),sigh(nmid+1)
        endif
        istart=1
      endif

      if(ndi<0.and.nmaxpr==1.and.mydiag)then
        print *,'Start cloud2 imax,idjd,myid ',imax,idjd,myid,cldoff
        if(idjd<=imax)then
          write(6,"('ttg ',18f7.2)")(ttg(idjd,k),k=1,kl)
          write(6,"('qlg   ',9g10.3)")(qlg(idjd,k),k=1,kl)
          write(6,"('qfg   ',9g10.3)")(qfg(idjd,k),k=1,kl)
          write(6,"('cfrac ',9g10.3)")(cfrac(idjd,k),k=1,kl)
          write(6,"('prf ',18f7.1)")(prf(idjd,k),k=1,kl)
          write(6,"('dprf ',18f7.1)")(dprf(idjd,k),k=1,kl)
          write(6,"('qccon ',9g10.3)")(qccon(idjd,k),k=1,kl)
          write(6,"('cdrop ',9g10.3)")(cdrop(idjd,k),k=1,kl)
          write(6,"('qo3   ',9g10.3)")(qo3(idjd,k),k=1,kl)  ! used by lwr88
          write(6,"('cosz ',9g10.3)")cosz(idjd)
       endif
      endif

c***INITIALIZE THE CLOUD AND CLOUD INDEX FIELDS
c     Except for the ground layer (nc=1) the assumption is that
c     no cloud exists. also, for this purpose, the cloud index is set
c     at one (p=0)
c     Don't set cirrf, cuvrf, cirab at the surface because these are set
c     the albedo by radfs.
         camt(:,1)=zero
         emcld(:,1)=one
         ktop(:,1)=1
         kbtm(:,1)=1
         ktopsw(:,1)=1
         kbtmsw(:,1)=1
      do k=2,lp1
         camt(:,k)=zero
         emcld(:,k)=one
         cirrf(:,k)=0.
         cuvrf(:,k)=0.
         cirab(:,k)=0.
         ktop(:,k)=1
         kbtm(:,k)=1
         ktopsw(:,k)=1
         kbtmsw(:,k)=1
      enddo
c***NOW SET CLOUD AND CLOUD INDEX FIELDS DEPENDING ON THE NO. OF CLOUDS
      nc=1
c---FIRST, THE ground layer (nc=1)
      emcld(:,nc)=one
      camt(:,nc)=one
      ktop(:,nc)=lp1
      kbtm(:,nc)=lp1
      ktopsw(:,nc)=lp1
      kbtmsw(:,nc)=lp1

      cll(:)=0.
      clm(:)=0.
      clh(:)=0.
      nclds(:)=0
      If (.not.cldoff) Then
c        cldliq(:)=0.
         qlptot(:)=0.
         taultot(:)=0.      
        
c Diagnose low, middle and high clouds; nlow,nmid are set up in initax.f
        
        if (nmr>0) then ! max-rnd overlap
        
          mx=0.
          do k=1,nlow
            mx=max(mx,cfrac(:,k))
            where (cfrac(:,k)==0.)
              cll=cll+mx*(1.-cll)
              mx=0.
            end where
          end do
          mx=0.
          do k=nlow+1,nmid
            mx=max(mx,cfrac(:,k))
            where (cfrac(:,k)==0.)
              clm=clm+mx*(1.-clm)
              mx=0.
            end where
          end do
          mx=0.
          do k=nmid+1,kl-1
            mx=max(mx,cfrac(:,k))
            where (cfrac(:,k)==0.)
              clh=clh+mx*(1.-clh)
              mx=0.
            end where
          end do          
          
        else ! usual random overlap
        
        do k=1,nlow
            cll(:)=cll(:)+cfrac(:,k)-cll(:)*cfrac(:,k)
        enddo
c        Note that cll, clm, clh are just diagnostic, so no effect from foll.          
c        if(nclddia<0)then
c          clm(:)=cfrac(:,1)
c          do k=2,nlow  !  find max value of low clouds
c           clm(:)=max(clm(:),cfrac(:,k))
c          enddo
c!         either use max_value for cll
c          if(nclddia==-1)cll(:)=clm(:)
c!         or take average of max_value, and random_overlap value          
c          if(nclddia==-2)cll(:)=.5*(cll(:)+clm(:))
c          clm(:)=0.
c        endif  ! (nclddia<0)
        do k=nlow+1,nmid
            clm(:)=clm(:)+cfrac(:,k)-clm(:)*cfrac(:,k)
        enddo
        do k=nmid+1,kl-1
            clh(:)=clh(:)+cfrac(:,k)-clh(:)*cfrac(:,k)
        enddo
        
        end if ! (nmr>0) ..else..

c Set up rk and cdrop (now as cdso4 from radriv90.f)

c This is the Liu and Daum scheme for relative dispersion (Nature, 419, 580-581 and pers. comm.)

        do k=1,kl-1
          do mg=1,imax
c            eps = 1. - 0.7 * exp(-0.008e-6*cdrop(mg,k)) !Upper bound
c            eps = 1. - 0.7 * exp(-0.001e-6*cdrop(mg,k)) !Lower bound
            eps = 1. - 0.7 * exp(-0.003e-6*cdrop(mg,k)) !mid range
!           beta = ((1+2*eps**2)**2/(1+eps**2))**(1./3)
!           rk(mg,k)=1./(beta**3)
            rk(mg,k)= (1.+eps**2)/(1.+2.*eps**2)**2
!           beta=1./rk(mg,k)**(1./3)   ! calc for diag only

C***            write(26,'(2f12.3)')1.e-6*cdrop(mg,k),beta
C***            eps = 1. - 0.7 * exp(-0.008e-6*cdrop(mg,k)) !Upper bound
C***            beta = ((1+2*eps**2)**2/(1+eps**2))**(1./3)
C***            write(27,'(2f12.3)')1.e-6*cdrop(mg,k),beta
          enddo
        enddo

c        do k=1,kl
c         do mg=1,imax
c          if(land(mg))then 
c            rk(mg,k)=0.67 ! land
c	   else
c            rk(mg,k)=0.8  !sea
c          endif
c         enddo
c        enddo

       
        !--------------------------------------------------------------------------------
        ! MJT CHANGE - mr
        ! Attempt to remove vertical dependence using GFDL AM2 approach
        ! Here we determine the top and bottom of the cloud using the largest overlap
        ! Cloud fraction is the largest cloud fraction
        ! Bulk cloud properties are calculated from average Reff and LWP.
        if (nmr==1) then

        taul=0.
        taui=0.
        Reffl=0.
        Reffi=0.
        Emw=0.
        Emi=0.
        cldtop=0
        cldbtm=0
        
c Locate clouds
          do mg=1,imax
            k=1
            do while (k.le.kl-1)
              if (cfrac(mg,k)>0.) then

                ! find cloud levels from k to kk
                ! break large changes in cloud properties into seperate clouds
                kk=k
                do while (cfrac(mg,kk+1)>0..and.kk<kl-1)
                  kk=kk+1
                end do
                ! find maximum cloud fraction level
                pos=maxloc(cfrac(mg,k:kk))
                km=pos(1)+k-1

                ! find cloud bottom
                csum=cfrac(mg,k)
                kb=k
                do kp=k+1,km
                  ctemp=max(0.,cfrac(mg,kp)-cfrac(mg,kp-1))
                  if (ctemp>csum) then
                    csum=ctemp
                    kb=kp
                  end if
                end do
                
                ! find cloud top
                csum=cfrac(mg,kk)
                kt=kk
                do kp=kk-1,km,-1
                  ctemp=max(0.,cfrac(mg,kp)-cfrac(mg,kp+1))
                  if (ctemp>csum) then
                    csum=ctemp
                    kt=kp
                  end if
                end do

                ! Store cloud data
                nclds(mg)=nclds(mg)+1
                nc=nclds(mg)+1
                camt(mg,nc)=cfrac(mg,km) ! max cloud fraction
                tau_sfac(mg,nc)=1.
                kbtm(mg,nc)=kl+1-kb
                kbtmsw(mg,nc)=kl+1-kb+1
                ktop(mg,nc)=kl+1-kt
                ktopsw(mg,nc)=kl+1-kt
                cldtop(mg,nc)=kk
                cldbtm(mg,nc)=k
                
                k=kk
              end if
              k=k+1
            end do
          end do

c Liquid water clouds
c Ice clouds : Choose scheme according to resolution
       IF(ldr.gt.0)THEN  ! 1,2,3  corresponds to previous lw=22 option
        
        refac1=0.85
        refac2=0.95
        do mg=1,imax
          do kp=1,nclds(mg)
            nc=kp+1
            dzsum=0.
            refflsum=0.
            reffisum=0.
            qlpathsum=0.
            qipathsum=0.
            qfsum=0.
            qlsum=0.
            do k=cldbtm(mg,nc),cldtop(mg,nc)
              ficehold=qfg(mg,k)/(qfg(mg,k)+qlg(mg,k))
              rhoa=100.*prf(mg,k)/(rdry*ttg(mg,k))
              dz=(dprf(mg,k)/prf(mg,k))*rdry*ttg(mg,k)/grav
              cfl=cfrac(mg,k)*(1.-ficehold)
              if (qlg(mg,k).gt.1.E-8) then
                Wliq=rhoa*qlg(mg,k)/cfl     !kg/m^3
c Reffl is the effective radius at the top of the cloud (calculated following
c Martin etal 1994, JAS 51, 1823-1842) due to the extra factor of 2 in the
c formula for reffl. Use mid cloud value of Reff for emissivity.
                Refflhold=
     &           (3.*2.*Wliq/(4.*rhow*pi*rk(mg,k)*Cdrop(mg,k)))**(1./3.)
                refflsum=refflsum+rhoa*qlg(mg,k)*dz*Refflhold
                qlpathsum=qlpathsum+rhoa*qlg(mg,k)*dz
                qlsum=qlsum+qlg(mg,k)*dz
              end if

              if (qfg(mg,k).gt.1.E-8) then
                Wice=rhoa*qfg(mg,k)/(cfrac(mg,k)*ficehold) !kg/m**3
                Reffihold=min(150.e-6,3.73e-4*Wice**0.216) !Lohmann et al.(1999)
                reffisum=reffisum+rhoa*qfg(mg,k)*dz*Reffihold
                qipathsum=qipathsum+rhoa*qfg(mg,k)*dz
                qfsum=qfsum+qfg(mg,k)*dz
              end if

              dzsum=dzsum+dz
            end do
            if (qfsum+qlsum.gt.0.) then
              fice(mg,nc)=qfsum/(qfsum+qlsum)
            else
              fice(mg,nc)=0.
            end if
            if (qlpathsum.gt.0.) then
              Reffl(mg,nc)=refflsum/qlpathsum
              qlpathsum=qlpathsum/(camt(mg,nc)*(1.-fice(mg,nc)))
              taul(mg,nc)=tau_sfac(mg,nc)*1.5*qlpathsum
     &                 /(rhow*Reffl(mg,nc))
              qlptot(mg)=qlptot(mg)+qlpathsum
              taultot(mg)=taultot(mg)+taul(mg,nc)*camt(mg,nc)
              cldht=dzsum/1000.           ! in km
              tclwc=qlpathsum*1000./dzsum ! in g/m**3
              tranw=exp(-1.66*cldht*50.885*tclwc**0.769917)
              Emw(mg,nc)=1.-tranw
            end if

            if (qipathsum.gt.0.) then
              Reffi(mg,nc)=reffisum/qipathsum
              qipathsum=qipathsum/(camt(mg,nc)*fice(mg,nc))
              taui(mg,nc)=1.5*qipathsum/(rhoice*Reffi(mg,nc))
              deltai=min(0.5*taui(mg,nc),45.) !IR optical depth for ice.
              taui(mg,nc)=tau_sfac(mg,nc)*taui(mg,nc)
c Ice-cloud emissivity following Platt
              if(taui(mg,nc).gt.0.4)then
                diffk=1.6
              else
                diffk=1.8
              endif
              Emi(mg,nc)=1.-exp(-diffk*deltai)
            end if
          end do
        end do
        
       ELSE  ! i.e. for ldr = -1,-2,-3

        refac1=0.90
        refac2=1.00
        do mg=1,imax
          do kp=1,nclds(mg)
            nc=kp+1
            dzsum=0.
            refflsum=0.
            ttgsum=0.
            qlpathsum=0.
            qipathsum=0.
            qlsum=0.
            qfsum=0.
            do k=cldbtm(mg,nc),cldtop(mg,nc)
              ficehold=qfg(mg,k)/(qfg(mg,k)+qlg(mg,k))
              rhoa=100.*prf(mg,k)/(rdry*ttg(mg,k))
              dz=(dprf(mg,k)/prf(mg,k))*rdry*ttg(mg,k)/grav
              cfl=cfrac(mg,k)*(1-ficehold)
              if (qlg(mg,k).gt.1.E-8) then
                Wliq=rhoa*qlg(mg,k)/cfl     !kg/m^3
c Reffl is the effective radius at the top of the cloud (calculated following
c Martin etal 1994, JAS 51, 1823-1842) due to the extra factor of 2 in the
c formula for reffl. Use mid cloud value of Reff for emissivity.
                Refflhold=
     &           (3.*2.*Wliq/(4.*rhow*pi*rk(mg,k)*Cdrop(mg,k)))**(1./3.)
                refflsum=refflsum+rhoa*qlg(mg,k)*dz*Refflhold
                qlpathsum=qlpathsum+rhoa*qlg(mg,k)*dz
                qlsum=qlsum+qlg(mg,k)*dz
              end if            

              if (qfg(mg,k).gt.1.E-8) then
                Wice=rhoa*qfg(mg,k)/(cfrac(mg,k)*ficehold) !kg/m**3
                qipathsum=qipathsum+rhoa*qfg(mg,k)*dz                
                qfsum=qfsum+qfg(mg,k)*dz
              end if

              dzsum=dzsum+dz
              ttgsum=ttgsum+dz*ttg(mg,k)
            end do
            if (qfsum+qlsum.gt.0.) then
              fice(mg,nc)=qfsum/(qfsum+qlsum)
            else
              fice(mg,nc)=0.
            end if
            if (qlpathsum.gt.0.) then
              Reffl(mg,nc)=refflsum/qlpathsum
              qlpathsum=qlpathsum/(camt(mg,nc)*(1.-fice(mg,nc)))
              taul(mg,nc)=tau_sfac(mg,nc)*1.5*qlpathsum
     &                   /(rhow*Reffl(mg,nc))
              qlptot(mg)=qlptot(mg)+qlpathsum
              taultot(mg)=taultot(mg)+taul(mg,nc)*camt(mg,nc)
              cldht=dzsum/1000.           ! in km
              tclwc=qlpathsum*1000./dzsum ! in g/m**3
              tranw=exp(-1.66*cldht*50.885*tclwc**0.769917)
              Emw(mg,nc)=1.-tranw
            end if
            
            if (qipathsum.gt.0.) then
              qipathsum=qipathsum/(camt(mg,nc)*fice(mg,nc))
              sigmai=aice*(qipathsum/dzsum)**bice !visible ext. coeff. for ice
              taui(mg,nc)=sigmai*dzsum !visible opt. depth for ice
              Reffi(mg,nc)=1.5*qipathsum/(rhoice*taui(mg,nc))
              taui(mg,nc)=tau_sfac(mg,nc)*taui(mg,nc)
c Ice-cloud emissivity following the Sunshine scheme
              tmid=ttgsum/dzsum-tfrz      !in celsius
              cldht=dzsum/1000.           !in km
              tciwc=qipathsum*1000./dzsum !in g/m**3

              temp_correction=1.047E+00+tmid
     &        *(-9.13e-05+tmid
     &        *(2.026e-04-1.056e-05*tmid))
              temp_correction=max(1.0, temp_correction)
              trani = real(exp(-1.66*temp_correction * tciwc*cldht
     &              / (0.630689d-01+0.265874*tciwc)))
c-- Limit ice cloud emissivities
c             trani=min(0.70,trani)
              kk=(2*(kl+1)-kbtm(mg,nc)-ktop(mg,nc))/2
              if(kk.gt.nlow) trani=min(0.70,trani)
              Emi(mg,nc) = 1.0 - trani    ! em is (1 - transmittance)
            end if
          end do
        end do
        
        ENDIF ! (ldr.gt.0) .. ELSE ..
          
c Calculate the SW cloud radiative properties for liquid water and ice clouds
c respectively, following Tony Slingo's (1989) Delta-Eddington scheme.

        call slingo (Reffl, taul, cosz, !inputs
     &       Rew1, Rew2, Abw )  !outputs
        
        call slingi (Reffi, taui, cosz, !inputs
     &       Rei1, Rei2, Abi )  !outputs

        onem=1.-1.e-6   ! to avoid possible later 0/0
        do mg=1,imax        
          do kp=1,nclds(mg)
            nc=kp+1
            qfsum=0.
            qlsum=0.
            qconsum=0.
            do k=cldbtm(mg,nc),cldtop(mg,nc)
              dz=(dprf(mg,k)/prf(mg,k))*rdry*ttg(mg,k)/grav
              qconsum=qconsum+qccon(mg,k)*dz
              qlsum=qlsum+qlg(mg,k)*dz
              qfsum=qfsum+qfg(mg,k)*dz
            end do
            if (qlsum+qfsum.gt.0.) then
              fcon=min(1.,qconsum/(qlsum+qfsum))
            else
              fcon=0.
            end if
c Mk3 with direct aerosol effect :
              refac=refac1*fcon+refac2*(1.-fcon)

              Rei1(mg,nc)=min(refac*Rei1(mg,nc),onem)
              Rei2(mg,nc)=min(refac*Rei2(mg,nc),onem-2.*Abi(mg,nc))
              Rew1(mg,nc)=min(refac*Rew1(mg,nc),onem)
              Rew2(mg,nc)=min(refac*Rew2(mg,nc),onem-2.*Abw(mg,nc))

              Re1=fice(mg,nc)*Rei1(mg,nc)+(1.-fice(mg,nc))*Rew1(mg,nc)
              Re2=fice(mg,nc)*Rei2(mg,nc)+(1.-fice(mg,nc))*Rew2(mg,nc)
              Em=fice(mg,nc)*Emi(mg,nc)+(1.-fice(mg,nc))*Emw(mg,nc)
              Ab=fice(mg,nc)*Abi(mg,nc)+(1.-fice(mg,nc))*Abw(mg,nc)
              emcld(mg,nc)=Em
              cuvrf(mg,nc)=Re1
              cirrf(mg,nc)=Re2
              cirab(mg,nc)=2.*Ab
              
         enddo
        enddo


        else ! usual random overlap (nmr=0) and M/R without bulk cloud properties (nmr=2), see swr99.f and clo89.f


c Define the emissivity (Em), and the SW properties (Re, Ab) for liquid (w)
c and ice (i) clouds respectively.
        
        do k=1,kl-1
          do mg=1,imax
            taul(mg,k)=0.
            taui(mg,k)=0.
            Reffl(mg,k)=0.
            Reffi(mg,k)=0.
            Emw(mg,k)=0.
            Emi(mg,k)=0.
            if(cfrac(mg,k).gt.0)then
              tau_sfac(mg,k)=1.
              fice(mg,k) = qfg(mg,k)/(qfg(mg,k)+qlg(mg,k))
            endif         !cfrac
          enddo
        enddo

c Liquid water clouds
        do k=1,kl-1
          do mg=1,imax
            if((cfrac(mg,k).gt.0).and.(qlg(mg,k).gt.1.0e-8))then
              rhoa=100*prf(mg,k)/(rdry*ttg(mg,k))
              dz=(dprf(mg,k)/prf(mg,k))*rdry*ttg(mg,k)/grav
              cfl=cfrac(mg,k)*(1-fice(mg,k))
              Wliq=rhoa*qlg(mg,k)/cfl     !kg/m^3
c              cldliq(mg)=cldliq(mg)+cfl-cldliq(mg)*cfl !Liquid cloud cover
                
c Reffl is the effective radius at the top of the cloud (calculated following
c Martin etal 1994, JAS 51, 1823-1842) due to the extra factor of 2 in the
c formula for reffl. Use mid cloud value of Reff for emissivity.
                
              Reffl(mg,k)=
     &           (3.*2.*Wliq/(4.*rhow*pi*rk(mg,k)*Cdrop(mg,k)))**(1./3.)
c              Reffl(mg,k)=
c     &             (3*Wliq/(4*rhow*pi*rk(mg,k)*Cdrop(mg,k)))**(1./3)
              qlpath=Wliq*dz
              taul(mg,k)=tau_sfac(mg,k)*1.5*qlpath/(rhow*Reffl(mg,k))
              qlptot(mg)=qlptot(mg)+qlpath*cfl
              taultot(mg)=taultot(mg)+taul(mg,k)*cfl

c Water cloud emissivity according to Martin Platt

C***          deltvl=taul(mg,k)*1.26 !Mult by 2^(1/3) so using mid cloud Reff
C***          deltal=min(0.4*deltvl, 45.) !IR optical depth for liq.
C***          if(deltvl.gt.0.4)then
C***            diffk=1.6
C***          else
C***            diffk=1.8
C***          endif
C***          Emw(mg,k) = 1.0 - exp(-diffk*deltal) !em of strat water clouds

c Or, water-cloud emissivity following the Sunshine scheme

              cldht=dz/1000.       !in km
              tclwc=Wliq*1000.     !in g/m**3
              tranw = exp( -1.66 * cldht * 50.885 * tclwc ** 0.769917)
              Emw(mg,k) = 1.0 - tranw    ! em is (1 - transmittance)
            endif         !cfrac
          enddo
        enddo
              
c Ice clouds : Choose scheme according to resolution
!      IF(lw.eq.22)THEN
       IF(ldr.gt.0)THEN  ! 1,2,3  corresponds to previous lw=22 option
        
        refac1=0.85
        refac2=0.95
        do k=1,kl-1
          do mg=1,imax
            if((cfrac(mg,k).gt.0).and.(qfg(mg,k).gt.1.0e-8))then
              rhoa=100*prf(mg,k)/(rdry*ttg(mg,k))
              dz=(dprf(mg,k)/prf(mg,k))*rdry*ttg(mg,k)/grav
              Wice=rhoa*qfg(mg,k)/(cfrac(mg,k)*fice(mg,k)) !kg/m**3
!             Reffi(mg,k)=3.73e-4*Wice**0.216 !Lohmann et al. (1999)
              Reffi(mg,k)=min(150.e-6,3.73e-4*Wice**0.216) !Lohmann et al.(1999)
              taui(mg,k)=1.5*Wice*dz/(rhoice*Reffi(mg,k))
              deltai=min(0.5*taui(mg,k), 45.) !IR optical depth for ice.
              taui(mg,k)=tau_sfac(mg,k)*taui(mg,k)

c Ice-cloud emissivity following Platt

              if(taui(mg,k).gt.0.4)then
                diffk=1.6
              else
                diffk=1.8
              endif
              Emi(mg,k) = 1.0 - exp(-diffk*deltai)
            endif         !cfrac
          enddo
        enddo

       ELSE  ! i.e. for ldr = -1,-2,-3

        refac1=0.90
        refac2=1.00
        do k=1,kl-1
          do mg=1,imax
            if((cfrac(mg,k).gt.0).and.(qfg(mg,k).gt.1.0e-8))then
              rhoa=100*prf(mg,k)/(rdry*ttg(mg,k))
              dz=(dprf(mg,k)/prf(mg,k))*rdry*ttg(mg,k)/grav
              Wice=rhoa*qfg(mg,k)/(cfrac(mg,k)*fice(mg,k)) !kg/m**3
              sigmai = aice*Wice**bice !visible ext. coeff. for ice
              taui(mg,k)=sigmai*dz !visible opt. depth for ice
              Reffi(mg,k)=1.5*Wice*dz/(rhoice*taui(mg,k))
              taui(mg,k)=tau_sfac(mg,k)*taui(mg,k)

c Ice-cloud emissivity following the Sunshine scheme

              tmid=ttg(mg,k)-tfrz  !in celsius
              cldht=dz/1000.       !in km
              tciwc=Wice*1000.     !in g/m**3

              temp_correction=1.047E+00+tmid
     &        *(-9.13e-05+tmid
     &        *(2.026e-04-1.056e-05*tmid))
              temp_correction=max(1.0, temp_correction)
              trani = real(exp(-1.66*temp_correction * tciwc*cldht
     &              / (0.630689d-01+0.265874*tciwc)))
c-- Limit ice cloud emissivities
c             trani=min(0.70,trani)
              if(k.gt.nlow) trani=min(0.70,trani)
              Emi(mg,k) = 1.0 - trani    ! em is (1 - transmittance)
            endif         !cfrac
          enddo
        enddo

       ENDIF  ! (ldr.gt.0) .. ELSE ..

c Calculate the effective radius of liquid water clouds seen from above

C***        do mg=1,imax
C***          if(taultot(mg).gt.0.)then
C***            Refflm(mg)=1.5*qlptot(mg)/(rhow*taultot(mg))
C***          else
C***            Refflm(mg)=0.
C***          endif
C***        enddo
          
c Calculate the SW cloud radiative properties for liquid water and ice clouds
c respectively, following Tony Slingo's (1989) Delta-Eddington scheme.

!       do k=1,kl
!         do mg=1,imax
!           if(cfrac(mg,k).gt.0.)then
C***              write(26,'(5f12.4)')cfrac(mg,k),taul(mg,k),taui(mg,k),
C***     &             coszro(mg),fice(mg,k)
C***              write(27,'(3f12.4)')
C***     &             cfrac(mg,k),1.e6*reffl(mg,k),1.e6*reffi(mg,k)
!           endif
!        enddo
!      enddo

        call slingo (Reffl, taul, cosz, !inputs
     &       Rew1, Rew2, Abw )  !outputs
        
        call slingi (Reffi, taui, cosz, !inputs
     &       Rei1, Rei2, Abi )  !outputs

        onem=1.-1.e-6   ! to avoid possible later 0/0
        do k=1,kl-1
          do mg=1,imax
            if(cfrac(mg,k).gt.0.)then
              fcon=min(1.,qccon(mg,k)/(qlg(mg,k)+qfg(mg,k)))
c Original refac :
c             refac=0.7*fcon+0.9*(1-fcon)
c Mk3 with no direct aerosol effect :
c             refac=0.7*fcon+0.85*(1-fcon)
c Mk3 with direct aerosol effect :
              refac=refac1*fcon+refac2*(1.-fcon)

              Rei1(mg,k)=min(refac*Rei1(mg,k),onem)
              Rei2(mg,k)=min(refac*Rei2(mg,k),onem-2.*Abi(mg,k))
              Rew1(mg,k)=min(refac*Rew1(mg,k),onem)
              Rew2(mg,k)=min(refac*Rew2(mg,k),onem-2.*Abw(mg,k))
            endif
          enddo
        enddo

c Weight cloud properties by liquid/ice fraction
        
        do k=1,kl-1
          do mg=1,imax
            if(cfrac(mg,k).gt.0.)then

              Re1 = fice(mg,k)*Rei1(mg,k) + (1.-fice(mg,k))*Rew1(mg,k)
              Re2 = fice(mg,k)*Rei2(mg,k) + (1.-fice(mg,k))*Rew2(mg,k)
              Em = fice(mg,k)*Emi(mg,k) + (1.-fice(mg,k))*Emw(mg,k)
!             if(qlg(mg,k).gt.0)write(26,'(2g12.3)')qlg(mg,k),rei1(mg,k)
!             if(qfg(mg,k).gt.0)write(27,'(2g12.3)')qfg(mg,k),emi(mg,k)

c             if(prf(mg,k).gt.800.) Em = 1.
              Ab = fice(mg,k)*Abi(mg,k) + (1-fice(mg,k))*Abw(mg,k)

!             write(30,*)'mg, nclds(mg) ',mg, nclds(mg)
              nclds(mg)=nclds(mg)+1
              nc=nclds(mg)+1
              camt(mg,nc)=cfrac(mg,k)
              ktop(mg,nc)=kl+1-k
              kbtm(mg,nc)=kl+1-k
              emcld(mg,nc)=Em
              ktopsw(mg,nc)=kl+1-k
              kbtmsw(mg,nc)=kl+1-k+1
              cuvrf(mg,nc)=Re1
              cirrf(mg,nc)=Re2
              cirab(mg,nc)=2*Ab

            endif
          enddo
        enddo

        
        end if ! (nmr.eq.1) ..else..
        !--------------------------------------------------------------------------------


        if(ndi<0.and.nmaxpr==1.and.idjd<=imax.and.mydiag)then
         print *,'After cloud2 myid',myid
         write(6,"('nclds ',i3)") nclds(idjd)
         write(6,"('ktop  ',18i3)")(ktop(idjd,k),k=1,kl)
         write(6,"('kbtm  ',18i3)")(kbtm(idjd,k),k=1,kl)
         write(6,"('ktopsw',18i3)")(ktopsw(idjd,k),k=1,kl)
         write(6,"('kbtmsw',18i3)")(kbtmsw(idjd,k),k=1,kl)
         write(6,"('Rew ',18f7.4)")
     &           ((Rew1(idjd,k)+Rew2(idjd,k))/2,k=1,kl)
         write(6,"('Rei ',18f7.4)")
     &           ((Rei1(idjd,k)+Rei2(idjd,k))/2,k=1,kl)
         write(6,"('Abw ',18f7.4)")(Abw(idjd,k),k=1,kl)
         write(6,"('Abi ',18f7.4)")(Abi(idjd,k),k=1,kl)
         write(6,"('Emw ',18f7.4)")(Emw(idjd,k),k=1,kl)
         write(6,"('Emi ',18f7.4)")(Emi(idjd,k),k=1,kl)
         write(6,"('camt ',18f7.4)")(camt(idjd,k),k=1,kl)
         write(6,"('emcld ',18f7.4)")(emcld(idjd,k),k=1,kl)
         write(6,"('cuvrf ',18f7.4)")(cuvrf(idjd,k),k=1,kl)
         write(6,"('cirrf ',18f7.4)")(cirrf(idjd,k),k=1,kl)
         write(6,"('cirab ',18f7.4)")(cirab(idjd,k),k=1,kl)
         write(6,"('Reffl ',9g10.3)")(Reffl(idjd,k),k=1,kl)
         write(6,"('Reffi ',9g10.3)")(Reffi(idjd,k),k=1,kl)
        endif
      endIf   !.not.cldoff      

      return
      end

c******************************************************************************

c Calculate SW radiative properties for water clouds using delta-Eddington
c scheme as given by Slingo (1989) JAS 46, 1419-1427.
c Coefficients di and fi modified to use Reff in SI units.

      subroutine slingo(reff, tau, mu0,       !inputs
     &                  refl1, refl2, abso )  !outputs

      use newmpar_m
      implicit none
C Global parameters
      include 'rdparm.h'
      integer nbands
      parameter (nbands=4)

C Argument list
      real reff(imax,kl)
      real tau(imax,kl)
      real mu0(imax)
      real refl1(imax,kl)
      real refl2(imax,kl)
      real abso(imax,kl)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks

C Local work arrays and variables
      real(kind=8) exparg,denom,epsilon,omega,f,omwf

      integer i
      integer k
      integer mg

      real absband
      real alpha1
      real alpha2
      real alpha3
      real alpha4
      real beta
      real beta0
      real e
      real gam1
      real gam2
      real gi
      real rdif
      real rm
      real rdir
      real tdb
      real tdif
      real tdir
      real ttot
      real u1
      real u2

C Local data, functions etc
      real ci(nbands)
      data ci / -5.62e-8, -6.94e-6, 4.64e-4,  2.01e-1  /

      real di(nbands)
c      data di / 1.63e-7,  2.35e-5,  1.24e-3,  7.56e-3  /
      data di / 1.63e-1,  2.35e+1,  1.24e+3,  7.56e+3  / !Si units

      real ei(nbands)
      data ei / 0.829,    0.794,    0.754,    0.826    /

      real fi(nbands)
c      data fi / 2.482e-3, 4.226e-3, 6.560e-3, 4.353e-3 /
      data fi / 2.482e+3, 4.226e+3, 6.560e+3, 4.353e+3 / !SI units

      real wi(nbands)
      data wi / 0.459760, 0.326158, 0.180608, 0.033474 /
 
      real wi24
      data wi24 / 0.540240 / !Sum of wi, i=2 to 4

C Start code : ----------------------------------------------------------

      do k=1,kl
        do mg=1,imax
          refl1(mg,k)=0.
          refl2(mg,k)=0.
          abso(mg,k)=0.
        enddo
      enddo

      do i=1,nbands
        do k=1,kl-1
          do mg=1,imax
            if(tau(mg,k).gt.0..and.mu0(mg).gt.0.)then
               reff(mg,k)=min(20.e-6,max(4.e-6,reff(mg,k)))
               omega=1-(ci(i)+di(i)*reff(mg,k))     
               gi=ei(i)+fi(i)*reff(mg,k)
               beta0=(3./7.)*(1-gi)
               beta=0.5-0.75*mu0(mg)*gi/(1+gi)
               f=gi**2
               U1=7./4.
               U2=real((7./4)*(1.-(1-omega)/(7*omega*beta0)))
               alpha1=real(U1*(1.-omega*(1-beta0)))
               alpha2=real(U2*omega*beta0)
               alpha3=real((1-f)*omega*beta)
               alpha4=real((1-f)*omega*(1-beta))
               epsilon=sqrt(alpha1**2-alpha2**2)
               rM=real(alpha2/(alpha1+epsilon))
               E=real(exp(-epsilon*tau(mg,k)))
               omwf=1-omega*f
               denom=omwf**2-epsilon**2*mu0(mg)**2
               gam1=real((omwf*alpha3-mu0(mg)*
     &               (alpha1*alpha3+alpha2*alpha4))/denom) 
               gam2=real((-omwf*alpha4-mu0(mg)*
     &               (alpha1*alpha4+alpha2*alpha3))/denom)
               exparg=min(real(70.0,8),omwf*tau(mg,k)/mu0(mg))
               Tdb=real(exp(-exparg))
               Rdif=rM*(1-E**2)/(1-(E*rM)**2)
               Tdif=E*(1-rM**2)/(1-(E*rM)**2)
               Rdir=-gam2*Rdif-gam1*Tdb*Tdif+gam1
               Tdir=-gam2*Tdif-gam1*Tdb*Rdif+gam2*Tdb
               Ttot=Tdb+Tdir
               Absband=1-Rdir-Ttot
               Absband=max(0., Absband) !Needed for 32 bit
c               if(absband.gt.1..or.absband.lt.0)then
c                 print*,'Warning slingo: band, abs =',i,absband
c               endif
               abso(mg,k)=abso(mg,k)+Absband*wi(i)
               if(i.eq.1)then
                 refl1(mg,k)=Rdir
               else
                 refl2(mg,k)=refl2(mg,k)+Rdir*wi(i)/wi24
               endif
             endif
           enddo
         enddo
       enddo

       return
       end

c******************************************************************************

c Slingo type scheme for ice cloud SW properties (similar to water clouds)
c Single scattering albedo follows Francis etal (1994) QJRMS 120, 809--848.
c Lambdas modified to use Reff in SI units.

      subroutine slingi(reff, tau, mu0,       !inputs
     &                  refl1, refl2, abso )  !outputs

      use const_phys
      use newmpar_m
      implicit none
C Global parameters
      include 'rdparm.h'
      integer nbands
      parameter (nbands=4)

C Argument list
      real reff(imax,kl)
      real tau(imax,kl)
      real mu0(imax)
      real refl1(imax,kl)
      real refl2(imax,kl)
      real abso(imax,kl)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks

C Local work arrays and variables
      real(kind=8) vm !For 32 bit
      real(kind=8) exparg,denom,epsilon,omega,f,omwf

      integer i
      integer k
      integer mg

      real absband
      real alpha1
      real alpha2
      real alpha3
      real alpha4
      real beta
      real beta0
      real e
      real gam1
      real gam2
      real gi
      real rdif
      real rm
      real rdir
      real tdb
      real tdif
      real tdir
      real ttot
      real u1
      real u2

C Local data, functions etc
      real lambda(nbands)
      data lambda / 0.470e-6, 0.940e-6, 1.785e-6, 3.190e-6 /

      real nprime(nbands)
      data nprime / 5.31e-9,  8.74e-7,  3.13e-4,  1.07e-1  /

      real wi(nbands)
      data wi     / 0.459760, 0.326158, 0.180608, 0.033474 /

      real wi24
      data wi24   / 0.540240 / !Sum of wi, i=2 to 4

C Start code : ----------------------------------------------------------

      do k=1,kl
        do mg=1,imax
          refl1(mg,k)=0.
          refl2(mg,k)=0.
          abso(mg,k)=0.
        enddo
      enddo

      do i=1,nbands
        do k=1,kl-1
          do mg=1,imax
            if(tau(mg,k).gt.0..and.mu0(mg).gt.0.)then
               vm=2*pi*reff(mg,k)*nprime(i)/lambda(i)
               omega=0.5+(vm+3)/(6*(vm+1)**3)
               gi=0.8  !Constant asymmetry parameter
               beta0=(3./7.)*(1-gi)
               beta=0.5-0.75*mu0(mg)*gi/(1+gi)
               f=gi**2
               U1=7./4.
               U2=real((7./4)*(1.-(1-omega)/(7*omega*beta0)))
               alpha1=real(U1*(1.-omega*(1-beta0)))
               alpha2=real(U2*omega*beta0)
               alpha3=real((1-f)*omega*beta)
               alpha4=real((1-f)*omega*(1-beta))
               epsilon=sqrt(alpha1**2-alpha2**2)
               rM=real(alpha2/(alpha1+epsilon))
               E=real(exp(-epsilon*tau(mg,k)))
               omwf=1-omega*f
               denom=omwf**2-epsilon**2*mu0(mg)**2
               gam1=real((omwf*alpha3-mu0(mg)*
     &               (alpha1*alpha3+alpha2*alpha4))/denom)
               gam2=real((-omwf*alpha4-mu0(mg)*
     &               (alpha1*alpha4+alpha2*alpha3))/denom)
               exparg=min(real(70.0,8),omwf*tau(mg,k)/mu0(mg))
               Tdb=real(exp(-exparg))
               Rdif=rM*(1-E**2)/(1-(E*rM)**2)
               Tdif=E*(1-rM**2)/(1-(E*rM)**2)
               Rdir=-gam2*Rdif-gam1*Tdb*Tdif+gam1
               Tdir=-gam2*Tdif-gam1*Tdb*Rdif+gam2*Tdb
               Ttot=Tdb+Tdir
               Absband=1-Rdir-Ttot
               Absband=max(0., Absband) !Needed for 32 bit
c               if(absband.gt.1..or.absband.lt.0)then
c                 print*,'Warning slingi: band, abs =',i,absband
c               endif
               abso(mg,k)=abso(mg,k)+Absband*wi(i)
               if(i.eq.1)then
                 refl1(mg,k)=Rdir
               else
                 refl2(mg,k)=refl2(mg,k)+Rdir*wi(i)/wi24
               endif
             endif
           enddo
         enddo
       enddo

       return
       end
