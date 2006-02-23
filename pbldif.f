      subroutine pbldif(theta,rkh,rkm,uav,vav)
!     vectorized version      
      use cc_mpi, only : mydiag, myid
      use diag_m
      implicit none
      integer ntest,nrkmin,npblmin,kmax
      parameter (ntest=0)
      parameter (nrkmin=2)  ! 1 original (& from 0510); 2 new; 3 newer
      parameter (npblmin=4) ! 1 original (best for Oz); 2 new ; 3,4 newer
      include 'newmpar.h'
      parameter (kmax=3*kl/4)  ! changed from kl/2 on 30/1/06
C------------------------------------------------------------------------
C 
C Atmospheric boundary layer computation.
C
C Nonlocal scheme that determines eddy diffusivities based on a
C diagnosed boundary layer height and a turbulent velocity scale;
C also, countergradient effects for heat and moisture, and constituents
C are included, along with temperature and humidity perturbations which 
C measure the strength of convective thermals in the lower part of the 
C atmospheric boundary layer.
C
C For more information, see Holtslag, A.A.M., and B.A. Boville, 1993:
C Local versus Nonlocal Boundary-Layer Diffusion in a Global Climate
C Model. J. Clim., vol. 6., p. 1825--1842.
c
c Updated by Holtslag and Hack to exclude the surface layer from the
c definition of the boundary layer Richardson number. Ri is now defined
c across the outer layer of the pbl (between the top of the surface
c layer and the pbl top) instead of the full pbl (between the surface and
c the pbl top). For simplicity, the surface layer is assumed to be the
c region below the first model level (otherwise the boundary layer depth 
c determination would require iteration).
C
C NOTE that all calculation in this module is at temperature points (DARLAM)
C------------------------------Code history--------------------------------
C
C Original version:  B. Boville
C Standardized:      J. Rosinski, June 1992
C Reviewed:          B. Boville, P. Rasch, August 1992
C Reviewed:          B. Boville, P. Rasch, April 1996
C
C Modified for boundary layer height diagnosis: Bert Holtslag, june 1994
C >>>>>>>>>  (Use ricr = 0.3 in this formulation)
C
C-----------------------------------------------------------------------
      include 'arrays.h'       !t
      include 'const_phys.h'
      include 'extraout.h'     !ustar
      include 'kuocom.h'
      include 'morepbl.h'      !fg,eg
      include 'parm.h'         !dtin
      include 'sigs.h'         !sig,sigmh
      include 'map.h'          !f (coriolis)
C------------------------------Arguments--------------------------------
C
C Input arguments:u,v,fg,eg,theta,ustar,uav,vav
      
C Input & Output arguments
      real rkm(ifull,kl)           ! eddy diffusivity for momentum [m2/s]
      real rkh(ifull,kl)           ! eddy diffusivity for heat [m2/s]
      real theta(ifull,kl)         ! potential temperature [K]
C     also qg                      ! mixing ratio [kg/kg}

      real cfrac
      common/cfrac/cfrac(ifull,kl)     ! globpe,radriv90,vertmix,convjlm
      common/work3/cgh(ifull,kl), ! counter-gradient term for heat [K/m]
     .             cgq(ifull,kl), ! counter-gradient term for constituents
     .             rino(ifull,kl),dum3(ifull,2*kl)
      common/work3d/zg(ifull,kl)
      integer iq,iflag
      real c1,zg,cgh,cgq,ztodtgor,delsig,tmp1,sigotbk,sigotbkm1
      real dum2,dum3
      real cgs                     ! counter-gradient star (cg/flux)
!     real pblh(ifull)             ! boundary-layer height [m] - in morepbl.h
C     Not used at the moment
C     real tpert(il)               ! convective temperature excess
C     real qpert(il)               ! convective humidity excess
C
C---------------------------Local parameters----------------------------
C
      real tiny                    ! lower bound for wind magnitude
      parameter (tiny=1.e-36)
C
C---------------------------Local workspace-----------------------------
C
      integer k                 ! level index

      real heatv                ! surface virtual heat flux
      real thvref               ! reference level virtual temperature
      real tkv                  ! model level potential temperature
      real therm                ! thermal virtual temperature excess
      real pmid                 ! midpoint pressures
      real phiminv              ! inverse phi function for momentum
      real phihinv              ! inverse phi function for heat 
      real wm                   ! turbulent velocity scale for momentum
      real vvk                  ! velocity magnitude squared
      real rkhfs                ! surface kinematic heat flux [mK/s]
      real rkqfs                ! sfc kinematic constituent flux [m/s]
      real zmzp                 ! level height halfway between zm and zp
      real rino                 ! bulk Richardson no. from level to ref lev
      real tlv                  ! ref. level pot tmp + tmp excess
      real fak1                 ! k*ustar*pblh
      real fak2                 ! k*wm*pblh
      real fak3                 ! fakn*wstr/wm 
      real pblk                 ! level eddy diffusivity for momentum
      real pr                   ! Prandtl number for eddy diffusivities
      real zm                   ! current level height
      real zp                   ! current level height + one level up
      real zl                   ! zmzp / Obukhov length
      real zh                   ! zmzp / pblh      at half levels
      real zzh                  ! (1-(zmzp/pblh))**2      at half levels
      real wstr                 ! w*, convective velocity scale
      real rrho                 ! 1./bottom level density (temporary)
      real obklen               ! Obukhov length
      real term                 ! intermediate calculation
      real fac                  ! interpolation factor
c     real pblmin               ! min pbl height due to mechanical mixing

c     logical unstbl(il)        ! pts w/unstbl pbl (positive virtual ht flx)
C------------------------------Commons----------------------------------
      real uav(ifull,kl),vav(ifull,kl)
      common/work2/heatv(ifull),rkhfs(ifull),rkqfs(ifull),thvref(ifull),
     .       phihinv(ifull),phiminv(ifull),tlv(ifull),wm(ifull),
     .       wstr(ifull),obklen(ifull),iflag(ifull),dum2(7*ifull)

      real betam   ! Constant in wind gradient expression
      real betas   ! Constant in surface layer gradient expression
      real betah   ! Constant in temperature gradient expression
      real fak     ! Constant in surface temperature excess
      real fakn    ! Constant in turbulent prandtl number
      real ricr    ! Critical richardson number
      real sffrac  ! Surface layer fraction of boundary layer
      real vk      ! Von Karman's constant
      real ccon    ! fak * sffrac * vk
      real binm    ! betam * sffrac
      real binh    ! betah * sffrac

      real zkmin   ! Minimum kneutral*f(ri)
      real rkmin   ! minimum eddy coeffs based on Hourdin et. al. (2001)

C Initialize COMCON

      data betam, betas, betah, sffrac / 15.0, 5.0, 15.0, 0.1 /
      data fakn,fak,vk,zkmin,c1/7.2, 8.5, 0.4, 0.01, 0.61/
      data ricr /0.25/
      
      fac=100.
      if(nlocal==6)fac=10.
      binh   = betah*sffrac
      binm   = betam*sffrac
      ccon   = fak*sffrac*vk
C****************************************************************
      do iq=1,ifull
       zg(iq,1)=bet(1)*t(iq,1)/grav
      enddo      ! iq loop
      do k=2,kl
       do iq=1,ifull
        zg(iq,k)=zg(iq,k-1)+
     .           (bet(k)*t(iq,k)+betm(k)*t(iq,k-1))/grav
       enddo        ! iq loop
      enddo         ! k  loop
      cgh(:,:)=0.   ! 3D
      cgq(:,:)=0.   ! 3D
      if(ktau==1)print *,'in pbldif nrkmin,npblmin: ',nrkmin,npblmin 

C Compute kinematic surface fluxes
         do iq=1,ifull
           pmid=ps(iq)*sigmh(1) 
           rrho = rdry*t(iq,1)/pmid
           ustar(iq) = max(ustar(iq),0.01)
           rkhfs(iq) = fg(iq)*rrho/cp           !khfs=w'theta'
           rkqfs(iq) = eg(iq)*rrho/hl              !kqfs=w'q'

C Compute various arrays for use later:

           thvref(iq) = theta(iq,1)*(1.0 + 0.61*qg(iq,1))
           heatv(iq)  = rkhfs(iq) + 0.61*theta(iq,1)*rkqfs(iq)
           wm(iq)     = 0.
c          obklen at t point
           obklen(iq) = -thvref(iq)*ustar(iq)**3/
     $             (grav*vk*(heatv(iq) + sign(1.e-10,heatv(iq))))
         
C >>>> Define first a new factor fac=100 for use in Richarson number
C      Calculate virtual potential temperature first level
C      and initialize pbl height to z1 i.e  1st full level

           pblh(iq) = zg(iq,1)    
           rino(iq,1) = 0.
          enddo

C PBL height calculation:
C Search for level of pbl. Scan upward until the Richardson number between
C the first level and the current level exceeds the "critical" value.
C Richardson no. is computed using eq. (4.d.18) NCAR technical report, CCM3)
C
          iflag(:)=0
          do k=2,kmax
           do iq=1,ifull
            vvk = (uav(iq,k) - uav(iq,1))**2 
     $          + (vav(iq,k) - vav(iq,1))**2
     $          + fac*ustar(iq)**2
c           vvk = max(vvk,tiny)
            tkv = theta(iq,k)*(1. + 0.61*qg(iq,k))
            rino(iq,k) = grav*(tkv - thvref(iq))*(zg(iq,k)-zg(iq,1))
     $                    /max(thvref(iq)*vvk,tiny)
            if(rino(iq,k)>=ricr.and.iflag(iq)==0)then
              pblh(iq) = zg(iq,k-1) + (ricr - rino(iq,k-1))/
     $                                (rino(iq,k) - rino(iq,k-1))
     $                               *(zg(iq,k) - zg(iq,k-1))
              iflag(iq)=1
            endif  ! (rino(iq,k)>=ricr.and.iflag(iq)==0)
           enddo  ! iq loop
          enddo   ! k loop
          if(nmaxpr==1.and.mydiag)then
            write (6,"('zg',9f8.1)") zg(idjd,1:kmax)
            write (6,"('rino_pa',9f8.3)") rino(idjd,1:kmax)
          endif

C Set pbl height to maximum value where computation exceeds number of
C layers allowed
 
C Improve estimate of pbl height for the unstable points.
C Find unstable points (virtual heat flux is positive):
C
         do iq=1,ifull
          if(heatv(iq)>0.)then  ! unstable case
            phiminv(iq) =     (1. - binm*pblh(iq)/obklen(iq))**(1./3.)
            wm(iq)= ustar(iq)*phiminv(iq)
C           therm: 2nd term in eq. (4.d.19):
C           temperature excess due to convective thermal
            therm = heatv(iq)*fak/wm(iq)
C           eq. (4.d.19) : tlv then used in eq. (4.d.18) to improve pblh
            tlv(iq) = thvref(iq) + therm
          end if  
        end do

C Improve pblh estimate for unstable conditions using the
C convective temperature excess:

        do k=2,kmax
         do iq=1,ifull
          vvk = (uav(iq,k) - uav(iq,1))**2 
     $        + (vav(iq,k) - vav(iq,1))**2
     $        + fac*ustar(iq)**2
          vvk = max(vvk,tiny)
          tkv = theta(iq,k)*(1. + 0.61*qg(iq,k))
          rino(iq,k) = grav*(tkv - tlv(iq))*(zg(iq,k)-zg(iq,1))
     $                  /max(thvref(iq)*vvk,tiny)     ! (see (4.d.18)
         enddo  !  i loop
        enddo   !  k loop

        iflag(:)=0
        do k=2,kmax
         do iq=1,ifull
          if(heatv(iq)>0..and.iflag(iq)==0)then  ! unstable case
            pblh(iq) = zg(iq,kl)    ! large default for unstable case
            if(rino(iq,k)>=ricr)then
              pblh(iq) = zg(iq,k-1) + (ricr - rino(iq,k-1))/
     $                                (rino(iq,k) - rino(iq,k-1))
     $                                *(zg(iq,k) - zg(iq,k-1))
              iflag(iq)=1  ! i.e. found it
            endif  ! (rino(iq,k)>=ricr)
          endif    ! (heatv(iq)>0..and.iflag(iq)==0)
         enddo     ! i loop
        enddo      ! k loop

C Points for which pblh exceeds number of pbl layers allowed;
C set to maximum
 
C PBL height must be greater than some minimum mechanical mixing depth
C Several investigators have proposed minimum mechanical mixing depth
C relationships as a function of the local friction velocity, u*.  We 
C make use of a linear relationship of the form h = c u* where c=700.
C The scaling arguments that give rise to this relationship most often 
C represent the coefficient c as some constant over the local coriolis
C parameter.  Here we make use of the experimental results of Koracin 
C and Berkowicz (1988) [BLM, Vol 43] for which they recommend 0.07/f
C where f was evaluated at 39.5 N and 52 N.  Thus we use a typical mid
C latitude value for f so that c = 0.07/f = 700.
 
      if(npblmin==1)pblh(:) = max(pblh(:),min(200.,700.*ustar(:)))
      if(npblmin==2)pblh(:) =
     &              max(pblh(:),.07*ustar(:)/max(.5e-4,abs(f(1:ifull))))
      if(npblmin==3)pblh(:) =                      ! to ~agree 39.5N
     &              max(pblh(:),.07*ustar(:)/max(1.e-4,abs(f(1:ifull))))
      if(npblmin==4)pblh(:) = max(pblh(:),50.)
c     if(npblmin==4)then  !  older Zilit., stable only
c       do iq=1,ifull
c        pblmin=.5*sqrt(ustar(iq)*max(0.,obklen(iq))/
c    .                            max(.5e-4,abs(f(iq))))
c        pblh(iq) = max(pblh(iq),pblmin)
c      end do
c     endif  ! (npblmin==4)

C pblh is now available; do preparation for diffusivity calculation:

C Do additional preparation for unstable cases only, set temperature
C and moisture perturbations depending on stability.

         do iq=1,ifull
          if(heatv(iq)>0.)then  ! unstable case
            phiminv(iq) =     (1. - binm*pblh(iq)/obklen(iq))**(1./3.)
            phihinv(iq) = sqrt(1. - binh*pblh(iq)/obklen(iq))
            wm(iq)      = ustar(iq)*phiminv(iq)
            wstr(iq)    = (heatv(iq)*grav*pblh(iq)/thvref(iq))**(1./3.)
          end if
        end do

C Main level loop to compute the diffusivities and 
C counter-gradient terms:

         if(nlocal==3)then
!          suppress nonlocal scheme over the sea   jlm
           do iq=1,ifull
            if(.not.land(iq))pblh(iq)=0.
           enddo
         endif  !  (nlocal==3)

         if(nlocal==4)then
!          suppress nonlocal scheme for column if cloudy layer in pbl   jlm
           do k=1,kl/2
            do iq=1,ifull
             if(zg(iq,k)<pblh(iq).and.cfrac(iq,k)>0.)pblh(iq)=0.
            enddo
           enddo
         endif  !  (nlocal==4)

         if(nlocal==5)then
!          suppress nonlocal scheme for column if cloudy layer in pbl   jlm
!          restores pblh at the bottom to have it available in convjlm/vertmix
           do k=1,kl/2
            do iq=1,ifull
              if(zg(iq,k)<pblh(iq).and.cfrac(iq,k)>0.)
     &         pblh(iq)=-pblh(iq)  
            enddo
           enddo
         endif  !  (nlocal==5)

        do k=1,kmax-1
         if(nlocal==2)then
!          suppress nonlocal scheme if cloudy layers in pbl   jlm
!          note this allows layers below to be done as nonlocal
           do iq=1,ifull
            if(zg(iq,k)<pblh(iq).and.cfrac(iq,k)>0.)pblh(iq)=0.
           enddo
         endif  !  (nlocal==2)

C Find levels within boundary layer:
C This is where Kh is at half model levels 
C zmzp = 0.5*(zm + zp)

           do iq=1,ifull
            fak1 = ustar(iq)*pblh(iq)*vk
            zm = zg(iq,k)
            zp = zg(iq,k+1)
!           if ( zkmin==0.0 .and. zp>pblh(iq)) zp = pblh(iq) ! zkmin=.01
            if (zm < pblh(iq)) then
              zmzp = 0.5*(zm + zp)
              zh = zmzp/pblh(iq)
              zl = zmzp/obklen(iq)
              zzh= 0.
              if (zh<=1.0) zzh = (1. - zh)**2

C stblev for points zm < plbh and stable and neutral
C unslev for points zm < plbh and unstable
              if(heatv(iq)>0.)then  ! unstable case
                fak2   = wm(iq)*pblh(iq)*vk
C unssrf, unstable within surface layer of pbl
C Unstable for surface layer; counter-gradient terms zero
                if (zh<sffrac) then
!                 term = cbrt(1. - betam*zl)
                  term =     (1. - betam*zl)**(1./3.)
                  pblk = fak1*zh*zzh*term
                  pr = term/sqrt(1. - betah*zl)
                else
C unsout, unstable within outer   layer of pbl
C Unstable for outer layer; counter-gradient terms non-zero:
                  pblk = fak2*zh*zzh
                  fak3 = fakn*wstr(iq)/wm(iq)
                  cgs     = fak3/(pblh(iq)*wm(iq))
                  cgh(iq,k) = rkhfs(iq)*cgs                 !eq. (4.d.17)
                  cgq(iq,k) = rkqfs(iq)*cgs                 !eq. (4.d.17)
                  pr = phiminv(iq)/phihinv(iq) + ccon*fak3/fak
                end if
                rkm(iq,k) = max(pblk,rkm(iq,k))
                rkh(iq,k) = max(pblk/pr,rkh(iq,k))
              elseif(nlocal>0)then    ! following are stable or neutral
C Stable and neutral points; set diffusivities; counter-gradient
C terms zero for stable case:
C term: pblk is Kc in eq. (4.d.16)
c but reverts to Louis stable treatment for nlocal=-1
                pblk = fak1*zh*zzh/(betas + zl)
                if (zl<=1.) then   ! 0 < z/L < 1.
                  pblk = fak1*zh*zzh/(1. + betas*zl)
                endif
                if(nrkmin==2)rkmin=vk*ustar(iq)*zmzp*zzh
                if(nrkmin==3)rkmin=max(rkh(iq,k),vk*ustar(iq)*zmzp*zzh)
                if(nrkmin==1.or.nlocal==6)rkmin=rkh(iq,k)
                if(ntest==1)then
                  if(iq==idjd)then
                    print *,'in pbldif k,ustar,zmzp,zh,zl,zzh ',
     .                                 k,ustar(iq),zmzp,zh,zl,zzh
                    print *,'rkh_L,rkmin,pblk,fak1,pblh ',
     .                       rkh(iq,k),rkmin,pblk,fak1,pblh(iq)
                  endif  ! (iq==idjd)
                endif    ! (ntest==1)
                rkm(iq,k) = max(pblk,rkmin)        
                rkh(iq,k) = rkm(iq,k)
              endif      ! (heatv(iq)>0.)    unstbl(i)
            endif        ! zm < pblh(iq)
           enddo         ! iq=1,ifull
       enddo             !end of k loop
       if(diag.and.mydiag)then
         if(heatv(idjd)>0.)   ! not meaningful or used otherwise
     &      write (6,"('rino_pb',9f8.3)") rino(idjd,1:kmax)
         print *,'ricr,obklen,heatv,pblh ',
     &            ricr,obklen(idjd),heatv(idjd),pblh(idjd)
         write (6,"('rkh_p',9f9.3/5x,9f9.3)") rkh(idjd,1:kl-2)
       endif

      ztodtgor = dtin*grav/rdry
C     update theta and qtg due to counter gradient
      do k=2,kmax-1
         do iq=1,ifull
            delsig = (sigmh(k+1)-sigmh(k))
            tmp1 = ztodtgor/delsig
            sigotbk=sigmh(k+1)/(0.5*(t(iq,k+1) + t(iq,k)))
            sigotbkm1=sigmh(k)/(0.5*(t(iq,k-1) + t(iq,k)))
            theta(iq,k) = theta(iq,k) + tmp1*
     $        (sigotbk*rkh(iq,k)*cgh(iq,k) -
     $         sigotbkm1*rkh(iq,k-1)*cgh(iq,k-1))
            qg(iq,k) = qg(iq,k) + tmp1*
     $        (sigotbk*rkh(iq,k)*cgq(iq,k) -
     $         sigotbkm1*rkh(iq,k-1)*cgq(iq,k-1))
         end do

      end do
      k=1
      do iq=1,ifull
         delsig = (sigmh(k+1)-sigmh(k))
         tmp1 = ztodtgor/delsig
         sigotbk=sigmh(k+1)/(0.5*(t(iq,k+1) + t(iq,k)))
         theta(iq,k) = theta(iq,k) + tmp1*
     $        sigotbk*rkh(iq,k)*cgh(iq,k)
         qg(iq,k) = qg(iq,k) + tmp1*
     $        sigotbk*rkh(iq,k)*cgq(iq,k)
      end do

      if (ntest>0) then
         print *,'pbldif'
         print *,'rkh= ',(rkh(idjd,k),k=1,kl)
         print *,'theta= ',(theta(idjd,k),k=1,kl)
         print *,'qg= ',(qg(idjd,k),k=1,kl)
         print *,'cgh= ',(cgh(idjd,k),k=1,kl)
         print *,'cgq= ',(cgq(idjd,k),k=1,kl)
       endif
C
C    Check for neg qtg's and put the original vertical
C    profile back if a neg value is found. A neg value implies that the
C    quasi-equilibrium conditions assumed for the countergradient term are
C    strongly violated.
C    Original code rewritten by Rosinski 7/8/91 to vectorize in longitude.

!     simpler alernative
c      qg=max(qg,qgmin)   ! 3D guided by McCormick et al 1993

      if(nlocal==5)then
!       restoring pblh to have it available in convjlm/vertmix jlm
        pblh(:)=abs(pblh(:))
      endif  !  (nlocal==5)

      return
      end
