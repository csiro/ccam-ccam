      subroutine hordifgt(iaero)    !  globpea version    N.B. k loop in here
!     usual scheme
      use aerosolldr
      use arrays_m
      use cc_mpi
      use cfrac_m
      use diag_m
      use dpsdt_m
      use indices_m
      use liqwpar_m
      use map_m
      use morepbl_m
      use nlin_m
      use parmhdff_m
      use sigs_m
      use tkeeps, only : tke,eps,shear,tkeww => ww,tkewx => dwdx,
     &                   tkewy => dwdy ! MJT tke
      use vecsuv_m
      use vvel_m
      implicit none
!      integer, parameter :: nhorjlm=1 ! 1 for jlm 3D deformation rather than Smagorinsky
c     called from globpe (now not tendencies),
c     called for -ve nhor
c     for +ve nhor see hordifg 
c     It has -ve nhorps option:
c        nhorps=0  does everything
c        nhorps=-1 does only T & qg horiz diff.
c        nhorps=-2 does only u &  v horiz diff.
c        nhorps=-3 does only qg     horiz diff.
c     and u,v have same options as T (e.g.nhor=-157)
c     this one has got map factors
c     N.B. no trace_gases yet
c     has jlm nhorx option as last digit of nhor, e.g. -157
      include 'newmpar.h'
      include 'const_phys.h'
      include 'parm.h'
 
      integer nhor,nhorps,khor,khdif,nhorjlm
      real hdifmax
      common/parmhdff/nhor,nhorps,khor,khdif,hdifmax,nhorjlm
     
      real, dimension(ifull+iextra,kl) :: uc, vc, wc, ee, ff, xfact,
     &                                    yfact, t_kh
      real, dimension(ifull) :: ptemp, tx_fact, ty_fact
      real, dimension(ifull,kl) :: zg,gg,bb,dd,ppb     ! MJT smag
      real, dimension(ifull,kl) :: theta,thetav,qs     ! MJT smag
      real, dimension(ifull,kl) :: dudx,dudy,dvdx,dvdy ! MJT smag
      real, dimension(ifull,kl) :: dwdx,dwdy,dwdz      ! MJT smag
      real, dimension(ifull,kl) :: dudz,dvdz           ! MJT smag
      real, dimension(ifull+iextra,kl) :: ww           ! MJT smag
      real es                                          ! MJT emag
      integer, parameter :: nf=2
!     Local variables
      integer iq, k, nhora, nhorx, iaero, l
      real cc, delphi, emi, hdif, ucc, vcc, wcc
      integer, save :: kmax=-1 ! MJT smag
      !integer i, j, n, ind
      !ind(i,j,n)=i+(j-1)*il+n*il*il  ! *** for n=0,5
      include 'establ.h'

c     nhorx used in hordif  ! previous code effectively has nhorx=0
c           = 1 u, v, T, qg  diffusion reduced near mountains
c           = 4              diffusion not reduced near mountains    
c           = 7 u, v, T, qg  diffusion reduced near mountains (bottom 2/3 only)

c     khor (set up in darlam.f): if non-zero increases the diffusion in the
c     upper model levels   (never used these days)
c         khor > 0 progessively doubles for khor upper levels
c         khor < 0 linearly increases above sig = 0.2 (up to factor of 1-khor)

c     this horizontal diffusion routine allows compensation for
c     being on sigma surfaces through reductions of kh over orography
c     this is called for -nhor.ge.50; the value of nhor (in metres) gives the
c     scaling value for deltaz
c     in namelist khdif is fujio's a**2, e.g. 4.
!     if(nhorps.gt.1)stop 'nhorps > 1 not permitted in hordifgt'

c     set up topography reduction factors for each type of location
c     expect power nf to be about 1 or 2 (see data statement)
      delphi=1.e6  ! turns off reduction (can also use nhorx=4)
      if(abs(nhor).ge.50)then
         nhora=10*(abs(nhor)/10)    ! e.g. 150  for nhor=-157
         nhorx=abs(nhor)-nhora      ! e.g.   7  for nhor=-157
         delphi=nhora*grav
      endif

      do iq=1,ifull
       ptemp(iq)=ps(iq)**.286
       tx_fact(iq)=1./(1.+(abs(zs(ie(iq))-zs(iq))/delphi)**nf)
       ty_fact(iq)=1./(1.+(abs(zs(in(iq))-zs(iq))/delphi)**nf)
      enddo   !  iq loop
c     above code independent of k

      if(diag.and.mydiag)then
         print *,'hordifgt u ',(u(idjd,k),k=1,kl)
         print *,'hordifgt v ',(v(idjd,k),k=1,kl)
         print *,'ax,ay,az ',ax(idjd),ay(idjd),az(idjd)
         print *,'bx,by,bz ',bx(idjd),by(idjd),bz(idjd)
      endif

      !--------------------------------------------------------------
      ! MJT smag
      ! replace kmax=2*kl/3 with sig(kmax)=0.2, no impact for 18 level
      ! version
      if (kmax.lt.0) then
        do k=1,kl
          if (sig(k).ge.0.2) kmax=k
        end do
      end if
      !--------------------------------------------------------------

      !--------------------------------------------------------------
      ! MJT tke ! MJT smag
      ! Calculate du/dx,dv/dx,du/dy,dv/dy but use cartesian vectors
      ! so as to avoid changes in vector direction across panel boundaries.
      ! Should be compatible with gnomic grid.
      if (nhorjlm==0.or.nhorjlm==4.or.nvmix==6) then
        ! neglect terrain following component for now
        ! u=ax*uc+ay*vc+az*wc
        ! dudx=u(ie)-u(iw)=ax*(uc(ie)-uc(iw))+ay*(vc(ie)-vc(iw))+az*(wc(ie)-wc(iw))
        ! dudy=u(in)-u(is)=ax*(uc(in)-uc(is))+ay*(vc(in)-vc(is))+az*(wc(in)-wc(is))
        do k=1,kl
          uc(1:ifull,k)=ax(1:ifull)*u(1:ifull,k)
          vc(1:ifull,k)=ay(1:ifull)*u(1:ifull,k)
          wc(1:ifull,k)=az(1:ifull)*u(1:ifull,k)
        end do
        call bounds(uc)
        call bounds(vc)
        call bounds(wc)
        do k=1,kl
          dudx(:,k)=(ax(1:ifull)*(uc(ie,k)-uc(iw,k))
     &              +ay(1:ifull)*(vc(ie,k)-vc(iw,k))
     &              +az(1:ifull)*(wc(ie,k)-wc(iw,k)))
     &              *0.5*em(1:ifull)/ds
          dudy(:,k)=(ax(1:ifull)*(uc(in,k)-uc(is,k))
     &              +ay(1:ifull)*(vc(in,k)-vc(is,k))
     &              +az(1:ifull)*(wc(in,k)-wc(is,k)))
     &              *0.5*em(1:ifull)/ds
        end do
        ! v=bx*uc+by*vc+bz*wc
        ! dvdx=v(ie)-v(iw)=bx*(uc(ie)-uc(iw))+by*(vc(ie)-vc(iw))+bz*(wc(ie)-wc(iw))
        ! dvdy=v(in)-v(is)=bx*(uc(in)-uc(is))+by*(vc(in)-vc(is))+bz*(wc(in)-wc(is))
        do k=1,kl
          uc(1:ifull,k)=bx(1:ifull)*v(1:ifull,k)
          vc(1:ifull,k)=by(1:ifull)*v(1:ifull,k)
          wc(1:ifull,k)=bz(1:ifull)*v(1:ifull,k)
        end do
        call bounds(uc)
        call bounds(vc)
        call bounds(wc)
        do k=1,kl
          dvdx(:,k)=(bx(1:ifull)*(uc(ie,k)-uc(iw,k))
     &              +by(1:ifull)*(vc(ie,k)-vc(iw,k))
     &              +bz(1:ifull)*(wc(ie,k)-wc(iw,k)))
     &              *0.5*em(1:ifull)/ds
          dvdy(:,k)=(bx(1:ifull)*(uc(in,k)-uc(is,k))
     &              +by(1:ifull)*(vc(in,k)-vc(is,k))
     &              +bz(1:ifull)*(wc(in,k)-wc(is,k)))
     &              *0.5*em(1:ifull)/ds
        end do
  
        ! calculate model level heights in meters (zg)
        ! calculate vertical velocity in m/s (ww)
        zg(:,1)=bet(1)*t(1:ifull,1)/grav
        do k=2,kl
          zg(:,k)=zg(:,k-1)+(bet(k)*t(1:ifull,k)
     &                      +betm(k)*t(1:ifull,k-1))/grav
        end do
        do k=1,kl        
          ! omega=ps*dpsldt
          ww(1:ifull,k)=(dpsldt(:,k)/sig(k)-dpsdt/(860.*ps(1:ifull)))
     &           *(-rdry/grav)*t(1:ifull,k)*(1.+0.61*qg(1:ifull,k))
        end do
        call bounds(ww)
        do k=1,kl
          dwdx(:,k)=(ww(ie,k)-ww(iw,k))*0.5*em(1:ifull)/ds
          dwdy(:,k)=(ww(in,k)-ww(is,k))*0.5*em(1:ifull)/ds
        end do
        dwdz(:,1)=(ww(1:ifull,2)-ww(1:ifull,1))
     &           /(zg(:,2)-zg(:,1))
        dudz(:,1)=(u(1:ifull,2)-u(1:ifull,1))
     &           /(zg(:,2)-zg(:,1))
        dvdz(:,1)=(v(1:ifull,2)-v(1:ifull,1))
     &           /(zg(:,2)-zg(:,1))
        do k=2,kl-1
          dwdz(:,k)=(ww(1:ifull,k+1)-ww(1:ifull,k-1))
     &             /(zg(:,k+1)-zg(:,k-1))
          dudz(:,k)=(u(1:ifull,k+1)-u(1:ifull,k-1))
     &             /(zg(:,k+1)-zg(:,k-1))
          dvdz(:,k)=(v(1:ifull,k+1)-v(1:ifull,k-1))
     &             /(zg(:,k+1)-zg(:,k-1))
        end do
        dwdz(:,kl)=(ww(1:ifull,kl)-ww(1:ifull,kl-1))
     &            /(zg(:,kl)-zg(:,kl-1))
        dudz(:,kl)=(u(1:ifull,kl)-u(1:ifull,kl-1))
     &            /(zg(:,kl)-zg(:,kl-1))
        dvdz(:,kl)=(v(1:ifull,kl)-v(1:ifull,kl-1))
     &            /(zg(:,kl)-zg(:,kl-1))
        
        ! Calculate bouyancy = kh * ppb (ppb=-N^2)
        ! Based on Durran and Klemp JAS 1982, where
        ! bb is for dry air and dd is for cloudy/saturated
        ! air
        do k=1,kl
          do iq=1,ifull
            theta(iq,k)=t(iq,k)*sig(k)**(-roncp)
            thetav(iq,k)=theta(iq,k)*(1.+0.61*qg(iq,k))
            es = establ(t(iq,k))
            qs(iq,k)=.622*es/(ps(iq)*sig(k)-es)
            gg(iq,k)=(1.+hl*qs(iq,k)/(rdry*t(iq,k)))
     &      /(1.+hl*hl*qs(iq,k)/(cp*rvap*t(iq,k)*t(iq,k)))
          end do
        end do
        bb(:,1)=-grav/(zg(:,2)-zg(:,1))*gg(:,1)
     &    *((theta(:,2)-theta(:,1))/theta(:,1)
     &    +hl*qs(:,1)/(cp*t(1:ifull,1))*(qs(:,2)-qs(:,1)))
     &    +grav/(zg(:,2)-zg(:,1))*(qs(:,2)+qlg(1:ifull,2)
     &    +qfg(1:ifull,2)-qs(:,1)-qlg(1:ifull,1)-qfg(1:ifull,1))
        dd(:,1)=-grav/(zg(:,2)-zg(:,1))
     &    *(thetav(:,2)-thetav(:,1))/thetav(:,1)
     &    +grav/(zg(:,2)-zg(:,1))*(qlg(1:ifull,2)+qfg(1:ifull,2)
     &    -qlg(1:ifull,1)-qfg(1:ifull,1))
        do k=2,kl-1
          bb(:,k)=-grav/(zg(:,k+1)-zg(:,k-1))*gg(:,k)
     &      *((theta(:,k+1)-theta(:,k-1))/theta(:,k)
     &      +hl*qs(:,k)/(cp*t(1:ifull,k))*(qs(:,k+1)-qs(:,k-1)))
     &      +grav/(zg(:,k+1)-zg(:,k-1))*(qs(:,k+1)+qlg(1:ifull,k+1)
     &      +qfg(1:ifull,k+1)-qs(:,k-1)-qlg(1:ifull,k-1)
     &      -qfg(1:ifull,k-1))
          dd(:,k)=-grav/(zg(:,k+1)-zg(:,k-1))
     &      *(thetav(:,k+1)-thetav(:,k-1))/thetav(:,k)
     &      +grav/(zg(:,k+1)-zg(:,k-1))*(qlg(1:ifull,k+1)
     &      +qfg(1:ifull,k+1)-qlg(1:ifull,k-1)-qfg(1:ifull,k-1))
        end do
        bb(:,kl)=-grav/(zg(:,kl)-zg(:,kl-1))*gg(:,kl)
     &    *((theta(:,kl)-theta(:,kl-1))/theta(:,kl)
     &    +hl*qs(:,kl)/(cp*t(1:ifull,kl))*(qs(:,kl)-qs(:,kl-1)))
     &    +grav/(zg(:,kl)-zg(:,kl-1))*(qs(:,kl)+qlg(1:ifull,kl)
     &    +qfg(1:ifull,kl)-qs(:,kl-1)-qlg(1:ifull,kl-1)
     &    -qfg(1:ifull,kl-1))
        dd(:,kl)=-grav/(zg(:,kl)-zg(:,kl-1))
     &    *(thetav(:,kl)-thetav(:,kl-1))/thetav(:,kl)
     &    +grav/(zg(:,kl)-zg(:,kl-1))*(qlg(1:ifull,kl)+qfg(1:ifull,kl)
     &    -qlg(1:ifull,kl-1)-qfg(1:ifull,kl-1))    
        !if (buoymeth.eq.1) then
        !  where (cfrac.gt.0.5)
        !    ppb=bb(:,2:kl)
        !  elsewhere
        !    ppb=cc(:,2:kl)
        !  end where
        !else
          ppb=cfrac*bb+(1.-cfrac)*dd
        !end if
      end if
      if (nhorjlm==1.or.nhorjlm==2.or.
     &    nhorps==0.or.nhorps==-2) then ! usual deformation for nhorjlm=1 or nhorjlm=2
        
        do k=1,kl
!        in hordifgt, need to calculate Cartesian components 
         do iq=1,ifull
            uc(iq,k) = ax(iq)*u(iq,k) + bx(iq)*v(iq,k)
            vc(iq,k) = ay(iq)*u(iq,k) + by(iq)*v(iq,k)
            wc(iq,k) = az(iq)*u(iq,k) + bz(iq)*v(iq,k)
         enddo
        end do
        call bounds(uc)
        call bounds(vc)
        call bounds(wc)
      
      end if
!      !--------------------------------------------------------------

      select case(nhorjlm)
       case(0,4) ! MJT smag
       !-------------------------------------------------------------
       ! MJT smag
       ! Smagorinsky REDUX
       ! This is based on 3D Smagorinsky closure (see WRF)
       ! This seems to help when the grid spacing is less than
       ! the boundary layer height.
         do k=1,kl
           hdif=dt*hdiff(k) ! N.B.  hdiff(k)=khdif*.1
           do iq=1,ifull
             cc=2.*dudx(iq,k)**2+2.*dvdy(iq,k)**2
             cc=cc+2.*dwdz(iq,k)**2
             cc=cc+(dvdx(iq,k)+dudy(iq,k))**2
             cc=cc+(dwdx(iq,k)+dudz(iq,k))**2
             cc=cc+(dwdy(iq,k)+dvdz(iq,k))**2
             cc=cc+3.*ppb(iq,k) ! 3*ppb=-N^2/Pr
             cc=max(cc,1.E-10)
             t_kh(iq,k)=sqrt(cc)*hdif/(em(iq)*em(iq))  ! this one with em in D terms
           end do
         end do
         if (nhorjlm==4) then
           do k=1,kl
             t_kh(1:ifull,k)=t_kh(1:ifull,k)*max(1.,
     &         pblh*pblh*em(1:ifull)*em(1:ifull)/(ds*ds))
           end do
         end if

       case(1)
c      jlm deformation scheme using 3D uc, vc, wc
         do k=1,kl
            hdif=dt*hdiff(k)/ds ! N.B.  hdiff(k)=khdif*.1  ! MJT bug fix
            do iq=1,ifull
               cc = (uc(ie(iq),k)-uc(iw(iq),k))**2 +
     &              (uc(in(iq),k)-uc(is(iq),k))**2 +
     &              (vc(ie(iq),k)-vc(iw(iq),k))**2 +
     &              (vc(in(iq),k)-vc(is(iq),k))**2 +
     &              (wc(ie(iq),k)-wc(iw(iq),k))**2 +
     &              (wc(in(iq),k)-wc(is(iq),k))**2
!              N.B. using double grid length
               t_kh(iq,k)= .5*sqrt(cc)*hdif/em(iq) ! this one without em in D terms
            enddo               !  iq loop
         enddo

       case(2)
c      jlm deformation scheme using 3D uc, vc, wc and omega (1st rough scheme)
         do k=1,kl
            hdif=dt*hdiff(k)/ds ! N.B.  hdiff(k)=khdif*.1  ! MJT bug fix
            do iq=1,ifull
               cc = (uc(ie(iq),k)-uc(iw(iq),k))**2 +
     &              (uc(in(iq),k)-uc(is(iq),k))**2 +
     &              (vc(ie(iq),k)-vc(iw(iq),k))**2 +
     &              (vc(in(iq),k)-vc(is(iq),k))**2 +
     &              (wc(ie(iq),k)-wc(iw(iq),k))**2 +
     &              (wc(in(iq),k)-wc(is(iq),k))**2 +
     & .01*(dpsldt(ie(iq),k)*ps(ie(iq))-dpsldt(iw(iq),k)*ps(iw(iq)))**2+
     & .01*(dpsldt(in(iq),k)*ps(in(iq))-dpsldt(is(iq),k)*ps(is(iq)))**2 
!         approx 1 Pa/s = .1 m/s     
!              N.B. using double grid length
               t_kh(iq,k)= .5*sqrt(cc)*hdif/em(iq) ! this one without em in D terms
            enddo               !  iq loop
         enddo

       case(3)
         t_kh=0. ! no diffusion (i.e., for pure nvmix.eq.6)
                 ! Probably works best for grid scales that
                 ! are less than 500 m, since that is the
                 ! maximum length scale allowed by the
                 ! prognostic eddy diffusivity - MJT

        case DEFAULT                                          ! MJT smag
         write(6,*) "ERROR: Unknown option nhorjlm=",nhorjlm  ! MJT smag
         stop                                                 ! MJT smag
       end select
       
       ! Calculate horizontal diffusion based on prognostic TKE
       ! This can be combined with the filters above operate
       ! over a large range of grid length scales
       if (nvmix.eq.6) then                                       ! MJT tke
         tke=max(tke,1.5E-8)                                      ! MJT tke
         eps=min(eps,(0.09**0.75)*(tke**1.5)/5.)                  ! MJT tke
         eps=max(eps,(0.09**0.75)*(tke**1.5)/500.)                ! MJT tke
         eps=max(eps,1.E-10)                                      ! MJT tke
         hdif=dt*0.09/(ds*ds)                                     ! MJT tke
         do k=1,kl                                                ! MJT tke
           t_kh(1:ifull,k)= max(max(tke(1:ifull,k)*tke(1:ifull,k) ! MJT tke
     &     /eps(1:ifull,k),1.E-7)*hdif,t_kh(1:ifull,k))           ! MJT tke
         end do                                                   ! MJT tke
         do k=1,kl                                                ! MJT tke
           ! vertical component included in tkeeps.f90 and        ! MJT tke
           ! additional terms for changing orography are          ! MJT tke
           ! neglected for now                                    ! MJT tke
           shear(:,k)=t_kh(1:ifull,k)*ds*ds/dt*(                  ! MJT tke
     &       (2.*dudx(:,k))**2+(2.*dvdy(:,k))**2                  ! MJT tke
     &      +(dudy(:,k)+dvdx(:,k))**2)                            ! MJT tke
         end do                                                   ! MJT tke
         tkewx=dwdx ! save dwdx for tke                           ! MJT tke
         tkewy=dwdy ! save dwdy for tke                           ! MJT tke
         tkeww=ww   ! save ww for tke                             ! MJT tke
       end if                                                     ! MJT tke

      call bounds(t_kh)
      do k=1,kl
         do iq=1,ifull
            xfact(iq,k) = (t_kh(ie(iq),k)+t_kh(iq,k))*.5
            yfact(iq,k) = (t_kh(in(iq),k)+t_kh(iq,k))*.5
         enddo
         !if((nhorx.ge.7.and.k.le.2*kl/3).or.nhorx.eq.1)then
         if((nhorx.ge.7.and.k.le.kmax).or.nhorx.eq.1)then ! MJT smag
            do iq=1,ifull
               xfact(iq,k) = xfact(iq,k)*tx_fact(iq)
               yfact(iq,k) = yfact(iq,k)*ty_fact(iq)
            enddo               !  iq loop
         endif                  ! (nhorx.ge.7.and.k.le.2*kl/3).or.nhorx.eq.1
      end do
      call boundsuv(xfact,yfact)

      if(nhorps.eq.0.or.nhorps.eq.-2)then ! for nhorps=-1,-3 don't diffuse u,v
         do k=1,kl
            do iq=1,ifull
               emi=1./em(iq)**2
               ucc = ( uc(iq,k)*emi +
     &                 xfact(iq,k)*uc(ie(iq),k) +
     &                 xfact(iwu(iq),k)*uc(iw(iq),k) +
     &                 yfact(iq,k)*uc(in(iq),k) +
     &                 yfact(isv(iq),k)*uc(is(iq),k) ) /
     &              ( emi + xfact(iq,k) + xfact(iwu(iq),k) +
     &                yfact(iq,k)+yfact(isv(iq),k) )
               vcc = ( vc(iq,k)*emi +
     &                 xfact(iq,k)*vc(ie(iq),k) +
     &                 xfact(iwu(iq),k)*vc(iw(iq),k) +
     &                 yfact(iq,k)*vc(in(iq),k) +
     &                 yfact(isv(iq),k)*vc(is(iq),k) ) /
     &               ( emi + xfact(iq,k) + xfact(iwu(iq),k) +
     &                 yfact(iq,k)+yfact(isv(iq),k) )
               wcc = ( wc(iq,k)*emi +
     &                 xfact(iq,k)*wc(ie(iq),k) +
     &                 xfact(iwu(iq),k)*wc(iw(iq),k) +
     &                 yfact(iq,k)*wc(in(iq),k) +
     &                 yfact(isv(iq),k)*wc(is(iq),k) ) /
     &               ( emi + xfact(iq,k) + xfact(iwu(iq),k) +
     &                 yfact(iq,k) + yfact(isv(iq),k) )
               u(iq,k) = ax(iq)*ucc + ay(iq)*vcc + az(iq)*wcc
               v(iq,k) = bx(iq)*ucc + by(iq)*vcc + bz(iq)*wcc
            enddo   !  iq loop
         end do
       endif   ! nhorps.ge.0

       if(diag.and.mydiag)then
          do k=1,kl
             print *,'k,id,jd,idjd ',k,id,jd,idjd
             print *,'k, xfact, xfactw ',k,xfact(idjd,k),
     &                                   xfact(iwu(idjd),k)
             print *,'k, yfact, yfacts ',k,yfact(idjd,k),
     &                                   yfact(isv(idjd),k)
             print *,'k, uc,uce,ucw,ucn,ucs '
     &        ,k,uc(idjd,k),uc(ie(idjd),k),uc(iw(idjd),k)
     &        ,uc(in(idjd),k),uc(is(idjd),k)
             print *,'k,ee,ff,u,v ',
     &            k,ee(idjd,k),ff(idjd,k),u(idjd,k),v(idjd,k)
          end do
       endif

       ! MJT - apply horizontal diffusion to TKE and EPS terms
       if (nvmix.eq.6) then
         ee(1:ifull,:)=tke(1:ifull,:)
         call bounds(ee)
         do k=1,kl
           do iq=1,ifull
             emi=1./em(iq)**2
             tke(iq,k) = ( ee(iq,k)*emi +
     &                     xfact(iq,k)*ee(ie(iq),k) +
     &                     xfact(iwu(iq),k)*ee(iw(iq),k) +
     &                     yfact(iq,k)*ee(in(iq),k) +
     &                     yfact(isv(iq),k)*ee(is(iq),k) ) /
     &                   ( emi + xfact(iq,k) + xfact(iwu(iq),k) +
     &                     yfact(iq,k)+yfact(isv(iq),k))
           enddo           !  iq loop
         end do
         ee(1:ifull,:)=eps(1:ifull,:)
         call bounds(ee)
         do k=1,kl
           do iq=1,ifull
             emi=1./em(iq)**2
             eps(iq,k) = ( ee(iq,k)*emi +
     &                     xfact(iq,k)*ee(ie(iq),k) +
     &                     xfact(iwu(iq),k)*ee(iw(iq),k) +
     &                     yfact(iq,k)*ee(in(iq),k) +
     &                     yfact(isv(iq),k)*ee(is(iq),k) ) /
     &                   ( emi + xfact(iq,k) + xfact(iwu(iq),k) +
     &                     yfact(iq,k)+yfact(isv(iq),k))
           end do           !  iq loop
         end do
       end if
       
       if (nhorjlm==0.or.nhorjlm==4) then
         ! increase by 1/Prandtl number for scalars
         xfact=xfact*3.
         yfact=yfact*3.
       end if

       if(nhorps.ne.-2)then   ! for nhorps=-2 don't diffuse T, qg
c        do t diffusion based on potential temperature ff
          do k=1,kl
             do iq=1,ifull
                ee(iq,k)=qg(iq,k)
                ff(iq,k)=t(iq,k)/ptemp(iq) ! watch out for Chen!
             enddo              !  iq loop
          end do
          call bounds(ee)
          call bounds(ff)
          if(nhorps.ne.-3)then  ! for nhorps=-3 don't diffuse T; only qg
             do k=1,kl
                do iq=1,ifull
                   emi=1./em(iq)**2
                   t(iq,k)= ptemp(iq) *
     &                      ( ff(iq,k)*emi +
     &                        xfact(iq,k)*ff(ie(iq),k) +
     &                        xfact(iwu(iq),k)*ff(iw(iq),k) +
     &                        yfact(iq,k)*ff(in(iq),k) +
     &                        yfact(isv(iq),k)*ff(is(iq),k) ) /
     &                      ( emi + xfact(iq,k) + xfact(iwu(iq),k) +
     &                        yfact(iq,k) + yfact(isv(iq),k) )
                enddo           !  iq loop
             end do
          endif                 ! (nhorps.ne.-3)
          do k=1,kl
             do iq=1,ifull
                emi=1./em(iq)**2
                qg(iq,k) = ( ee(iq,k)*emi +
     &                       xfact(iq,k)*ee(ie(iq),k) +
     &                       xfact(iwu(iq),k)*ee(iw(iq),k) +
     &                       yfact(iq,k)*ee(in(iq),k) +
     &                       yfact(isv(iq),k)*ee(is(iq),k) ) /
     &                     ( emi + xfact(iq,k) + xfact(iwu(iq),k) +
     &                       yfact(iq,k)+yfact(isv(iq),k))
             end do              !  iq loop
          end do
       endif                    ! (nhorps.ge.-1)
       
       ! MJT aerosols
       if (abs(iaero).eq.2) then
         do l=1,naero
           ee(1:ifull,:)=xtg(1:ifull,:,l)
           call bounds(ee)
           do k=1,kl
             do iq=1,ifull
               emi=1./em(iq)**2
               xtg(iq,k,l) = ( ee(iq,k)*emi +
     &                     xfact(iq,k)*ee(ie(iq),k) +
     &                     xfact(iwu(iq),k)*ee(iw(iq),k) +
     &                     yfact(iq,k)*ee(in(iq),k) +
     &                     yfact(isv(iq),k)*ee(is(iq),k) ) /
     &                   ( emi + xfact(iq,k) + xfact(iwu(iq),k) +
     &                     yfact(iq,k)+yfact(isv(iq),k))
             enddo           !  iq loop
           end do
         end do
       end if


      return
      end
