      subroutine leoncld(cfrac)

      implicit none
      include 'newmpar.h'
      include 'liqwpar.h' ! ifullw
      include 'const_phys.h' !Input physical constants
      include 'cparams.h'    !Input cloud scheme parameters

      include 'arrays.h'
c     include 'constant.h'
      include 'dava.h'    ! davt
      include 'kuocom.h'  ! kbsav,ktsav,convfact,convpsav,ndavconv
      include 'map.h'     ! land
      include 'morepbl.h'
      include 'nlin.h'
      include 'parm.h'
      include 'prec.h'
      include 'sigs.h'
      include 'soil.h'
      include 'tracers.h'  ! ngas, nllp, ntrac
      include 'vvel.h'
      include 'xarrs.h'

      real delq,dels,rnrt3d
      common/work3f/delq(ifull,kl),dels(ifull,kl),rnrt3d(ifull,kl) !Pass in rnrt3d from convjlm 

c Local variables
      integer iq,k
      real tdiff,tm,establ,table

      real prf(ifullw,kl)     !Pressure on full levels (hPa)
      real dprf(ifullw,kl)    !Pressure thickness (hPa)
      real rhoa(ifullw,kl)    !Air density (kg/m3)
      real dz(ifullw,kl)      !Layer thickness (m)
      real cdso4(ifullw,kl)   !Cloud droplet conc (#/m3)
      real cfrac(ifullw,kl)   !Cloud fraction
      real ccov(ifullw,kl)    !Cloud cover (may differ from cloud frac if vertically subgrid)
      real cfa(ifullw,kl)     !Cloud fraction in which autoconv occurs (option in newrain.f)
      real qca(ifullw,kl)     !Cloud water mixing ratio in cfa(:,:)    (  "    "     "     )
      real fluxc(ifullw,kl)   !Flux of convective rainfall in timestep (kg/m**2)
      real ccrain(ifullw,kl)  !Convective raining cloud cover
      real precs(ifullw)      !Amount of stratiform precipitation in timestep (mm)
      real preci(ifullw)      !Amount of stratiform snowfall in timestep (mm)
      real dum3d(ifullw,kl)   !Dummy 3d outputs from newrain (not needed here)

      integer kbase(ifullw),ktop(ifullw) !Bottom and top of convective cloud 
      include 'establ.h'

      do k=1,kl   
         do iq=1,ifull
C***           if(t(iq,k).gt.343.16)then
C***             print*,'Start leoncld: Warning, iq,k,t = ',iq,k,t(iq,k)
C***           endif

          prf(iq,k)=0.01*ps(iq)*sig(k) !Looks like ps is SI units
          dprf(iq,k)=-0.01*ps(iq)*dsig(k) !dsig is -ve
          rhoa(iq,k)=100.*prf(iq,k)/(rdry*t(iq,k))
          if(land(iq))then
            cdso4(iq,k)=cdropl
          else
            cdso4(iq,k)=cdrops
          endif
        enddo
      enddo

      kbase(:)=0  !Not used for now
      ktop(:)=0
      dz(:,:)=100.*dprf(:,:)/(rhoa(:,:)*grav)
c      fluxc(:,:)=rnrt3d(:,:)*1.e-3*dt ! kg/m2 (should be same level as rnrt3d)
      fluxc(:,:)=0. !For now... above line may be wrong
      ccrain(:,:)=0.1  !Assume this for now
      precs(:)=0.

c Calculate cloud fraction and cloud water mixing ratios

      call newcloud(dt,1,land,prf,kbase,ktop,rhoa,cdso4, !Inputs
     &     t,qg,qlg,qfg,   !In and out
     &     cfrac,ccov,cfa,qca)   !Outputs


c Calculate precipitation and related processes
      
      call newrain(land,1,dt,fluxc,rhoa,dz,ccrain,prf,cdso4,  !Inputs
     &    cfa,qca,                                            !Inputs
     &    t,qlg,qfg,precs,qg,cfrac,ccov,                        !In and Out
     &    preci,dum3d,dum3d,dum3d,dum3d,dum3d,dum3d,dum3d,dum3d,!Outputs
     &    dum3d,dum3d,dum3d,dum3d)                            !Outputs

C***      do k=1,kl   
C***         do iq=1,ifull
C***           if(t(iq,k).gt.343.16)then
C***             print*,'End leoncld: Warning, iq,k,t = ',iq,k,t(iq,k)
C***           endif
C***         enddo
C***       enddo

      return
      end



