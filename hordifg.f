      subroutine hordifgt   !  globpea version    N.B. k loop in here
!     usual scheme
      use aerosolldr
      use arrays_m
      use cc_mpi
      use cfrac_m
      use cloudmod
      use diag_m
      use dpsdt_m
      use indices_m
      use liqwpar_m
      use map_m
      use morepbl_m
      use nharrs_m
      use nlin_m
      use parmhdff_m
      use savuvt_m
      use savuv1_m
      use sigs_m
      use tkeeps, only : tke,eps,shear,mintke,mineps,cm0,minl,maxl
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
      include 'kuocom.h'
      include 'parm.h'
      include 'parmdyn.h'
 
      real, dimension(ifull+iextra,kl) :: uc, vc, wc, xfact,
     &                                    yfact, t_kh
      real, dimension(ifull+iextra,kl,nagg) :: ff
      real, dimension(ifull,kl) :: dudx,dudy,dudz
      real, dimension(ifull,kl) :: dvdx,dvdy,dvdz
      real, dimension(ifull,kl) :: dwdx,dwdy,dwdz
      real, dimension(ifull,kl) :: base
      real, dimension(ifull) :: ptemp, tx_fact, ty_fact
      real, dimension(ifull) :: sx_fact, sy_fact
      real, dimension(ifull) :: emi,r1,r2
      real, dimension(ifull,2) :: zgh      
      real, dimension(ifull+iextra,kl) :: ww,uav,vav
      real, dimension(ifull,kl) :: zg,tnhs
      integer, parameter :: nf=2
!     Local variables
      integer iq, k, nhora, nhorx, ntr
      real cc, delphi, hdif, ucc, vcc, wcc
      integer, save :: kmax=-1
      !integer i, j, n, ind
      !ind(i,j,n)=i+(j-1)*il+n*il*il  ! *** for n=0,5

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

      if (kmax<0) then
        kmax=1
        do while (sig(kmax)>=0.25.and.kmax<kl)
          kmax=kmax+1
        end do
      end if

c     set up topography reduction factors for each type of location
c     expect power nf to be about 1 or 2 (see data statement)

      delphi=1.e6  ! turns off reduction (can also use nhorx=4)
      if (abs(nhor)>=50) then
        nhora=10*(abs(nhor)/10)    ! e.g. 150  for nhor=-157
        nhorx=abs(nhor)-nhora      ! e.g.   7  for nhor=-157
        delphi=nhora*grav
      endif

      do iq=1,ifull
        emi(iq)=1./(em(iq)*em(iq))
        ptemp(iq)=ps(iq)**.286
        tx_fact(iq)=1./(1.+(abs(zs(ie(iq))-zs(iq))/delphi)**nf)
        ty_fact(iq)=1./(1.+(abs(zs(in(iq))-zs(iq))/delphi)**nf)
        sx_fact(iq)=1./(1.+(0.5*abs(zs(ie(iq))-zs(iw(iq)))/delphi)**nf)
        sy_fact(iq)=1./(1.+(0.5*abs(zs(in(iq))-zs(is(iq)))/delphi)**nf)
      enddo   !  iq loop
c     above code independent of k

#ifdef debug
      if(diag.and.mydiag)then
        write(6,*) 'hordifgt u ',(u(idjd,k),k=1,kl)
        write(6,*) 'hordifgt v ',(v(idjd,k),k=1,kl)
        write(6,*) 'ax,ay,az ',ax(idjd),ay(idjd),az(idjd)
        write(6,*) 'bx,by,bz ',bx(idjd),by(idjd),bz(idjd)
      endif
#endif

      !--------------------------------------------------------------
      ! This option is for prognostic TKE
      if (nhorjlm==0.or.nvmix==6) then
        ! Calculate du/dx,dv/dx,du/dy,dv/dy, etc 

        ! calculate height on full levels and non-hydrostatic temp correction
        tnhs(:,1)=phi_nh(:,1)/bet(1)
        zg(:,1)=(zs(1:ifull)+bet(1)*(t(1:ifull,1)+tnhs(:,1)))/grav
        do k=2,kl
         tnhs(:,k)=(phi_nh(:,k)-phi_nh(:,k-1)-betm(k)*tnhs(:,k-1))
     &             /bet(k)
         zg(:,k)=zg(:,k-1)+(bet(k)*(t(1:ifull,k)+tnhs(:,k))
     &                    +betm(k)*(t(1:ifull,k-1)+tnhs(:,k-1)))/grav
        end do ! k  loop

        do k=1,kl        
          ! weighted horizontal velocities
          uav(1:ifull,k)=av_vmod*u(1:ifull,k)
     &                 +(1.-av_vmod)*savu(1:ifull,k)
          vav(1:ifull,k)=av_vmod*v(1:ifull,k)
     &                 +(1.-av_vmod)*savv(1:ifull,k)
        
          ! calculate vertical velocity in m/s
          ! omega=ps*dpsldt
          ! ww = -R/g * (T+Tnhs) * dpsldt/sig
          ww(1:ifull,k)=(dpsldt(:,k)/sig(k))
     &        *(-rdry/grav)*(t(1:ifull,k)+tnhs(:,k))
     &        *(1.+0.61*qg(1:ifull,k)-qlg(1:ifull,k)-qfg(1:ifull,k))
        end do
        
        call boundsuv_allvec(uav,vav)
        call bounds(ww)

        do k=1,kl
          dudx(:,k)=0.5*(uav(ieu,k)-uav(iwu,k))*em(1:ifull)/ds
          dudy(:,k)=0.5*(uav(inu,k)-uav(isu,k))*em(1:ifull)/ds
          dvdx(:,k)=0.5*(vav(iev,k)-vav(iwv,k))*em(1:ifull)/ds
          dvdy(:,k)=0.5*(vav(inv,k)-vav(isv,k))*em(1:ifull)/ds
          dwdx(:,k)=0.5*(ww(ie,k)-ww(iw,k))*em(1:ifull)/ds
          dwdy(:,k)=0.5*(ww(in,k)-ww(is,k))*em(1:ifull)/ds
        end do
        
        ! calculate vertical gradients
        zgh(:,2)=ratha(1)*zg(:,2)+rathb(1)*zg(:,1) ! upper half level
        r1=uav(1:ifull,1)
        r2=ratha(1)*uav(1:ifull,2)+rathb(1)*uav(1:ifull,1)          
        dudz(1:ifull,1)=(r2-r1)/(zgh(1:ifull,2)-zg(1:ifull,1))
        r1=vav(1:ifull,1)
        r2=ratha(1)*vav(1:ifull,2)+rathb(1)*vav(1:ifull,1)
        dvdz(1:ifull,1)=(r2-r1)/(zgh(1:ifull,2)-zg(1:ifull,1))
        r1=ww(1:ifull,1)
        r2=ratha(1)*ww(1:ifull,2)+rathb(1)*ww(1:ifull,1)          
        dwdz(1:ifull,1)=(r2-r1)/(zgh(1:ifull,2)-zg(1:ifull,1))
        do k=2,kl-1
          zgh(:,1)=zgh(:,2) ! lower half level
          zgh(:,2)=ratha(k)*zg(:,k+1)+rathb(k)*zg(:,k) ! upper half level
          r1=ratha(k-1)*uav(1:ifull,k)+rathb(k-1)*uav(1:ifull,k-1)
          r2=ratha(k)*uav(1:ifull,k+1)+rathb(k)*uav(1:ifull,k)          
          dudz(1:ifull,k)=(r2-r1)/(zgh(1:ifull,2)-zgh(1:ifull,1))
          r1=ratha(k-1)*vav(1:ifull,k)+rathb(k-1)*vav(1:ifull,k-1)
          r2=ratha(k)*vav(1:ifull,k+1)+rathb(k)*vav(1:ifull,k)          
          dvdz(1:ifull,k)=(r2-r1)/(zgh(1:ifull,2)-zgh(1:ifull,1))
          r1=ratha(k-1)*ww(1:ifull,k)+rathb(k-1)*ww(1:ifull,k-1)
          r2=ratha(k)*ww(1:ifull,k+1)+rathb(k)*ww(1:ifull,k)          
          dwdz(1:ifull,k)=(r2-r1)/(zgh(1:ifull,2)-zgh(1:ifull,1))
        end do
        zgh(:,1)=zgh(:,2) ! lower half level
        r1=ratha(kl-1)*uav(1:ifull,kl)+rathb(kl-1)*uav(1:ifull,kl-1)
        r2=uav(1:ifull,kl)          
        dudz(1:ifull,kl)=(r2-r1)/(zg(1:ifull,kl)-zgh(1:ifull,1))
        r1=ratha(kl-1)*vav(1:ifull,kl)+rathb(kl-1)*vav(1:ifull,kl-1)
        r2=vav(1:ifull,kl)          
        dvdz(1:ifull,kl)=(r2-r1)/(zg(1:ifull,kl)-zgh(1:ifull,1))
        r1=ratha(kl-1)*ww(1:ifull,kl)+rathb(kl-1)*ww(1:ifull,kl-1)
        r2=ww(1:ifull,kl)          
        dwdz(1:ifull,kl)=(r2-r1)/(zg(1:ifull,kl)-zgh(1:ifull,1))
        
      end if   ! nhorjlm==0.or.nvmix==6
      
      if (nhorjlm==1.or.nhorjlm==2.or.
     &    nhorps==0.or.nhorps==-2) then ! usual deformation for nhorjlm=1 or nhorjlm=2
        
        do k=1,kl
!        in hordifgt, need to calculate Cartesian components 
         do iq=1,ifull
            ff(iq,k,1) = ax(iq)*u(iq,k) + bx(iq)*v(iq,k)
            ff(iq,k,2) = ay(iq)*u(iq,k) + by(iq)*v(iq,k)
            ff(iq,k,3) = az(iq)*u(iq,k) + bz(iq)*v(iq,k)
         enddo
        end do
        call bounds(ff(:,:,1:3))
        uc(:,:)=ff(:,:,1)
        vc(:,:)=ff(:,:,2)
        wc(:,:)=ff(:,:,3)
      
      end if
!      !--------------------------------------------------------------

      select case(nhorjlm)
       case(0)
         ! This is based on 2D Smagorinsky closure
         do k=1,kl
           hdif=dt*hdiff(k) ! N.B.  hdiff(k)=khdif*.1
           r1=(dudx(:,k)-dvdy(:,k))**2
     &       +(dvdx(:,k)+dudy(:,k))**2
           t_kh(1:ifull,k)= sqrt(r1)*hdif*emi(:)
         end do
         call bounds(t_kh,nehalf=.true.)
         do k=1,kl
            do iq=1,ifull
               xfact(iq,k) = (t_kh(ie(iq),k)+t_kh(iq,k))*.5
               yfact(iq,k) = (t_kh(in(iq),k)+t_kh(iq,k))*.5
            end do
         end do
             
       case(1)
c      jlm deformation scheme using 3D uc, vc, wc
         do k=1,kl
            hdif=dt*hdiff(k)/ds ! N.B.  hdiff(k)=khdif*.1
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
         call bounds(t_kh,nehalf=.true.)
         do k=1,kl
            do iq=1,ifull
               xfact(iq,k) = (t_kh(ie(iq),k)+t_kh(iq,k))*.5
               yfact(iq,k) = (t_kh(in(iq),k)+t_kh(iq,k))*.5
            end do
         end do

       case(2)
c      jlm deformation scheme using 3D uc, vc, wc and omega (1st rough scheme)
         do k=1,kl
            hdif=dt*hdiff(k)/ds ! N.B.  hdiff(k)=khdif*.1
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
         call bounds(t_kh,nehalf=.true.)
         do k=1,kl
            do iq=1,ifull
               xfact(iq,k) = (t_kh(ie(iq),k)+t_kh(iq,k))*.5
               yfact(iq,k) = (t_kh(in(iq),k)+t_kh(iq,k))*.5
            end do
         end do

       case(3)
        if (nvmix==6) then
         hdif=dt*cm0
         do k=1,kl
           tke(1:ifull,k)=max(tke(1:ifull,k),mintke)
           r1=(cm0**0.75)*tke(1:ifull,k)
     &                 *sqrt(tke(1:ifull,k))
           eps(1:ifull,k)=min(eps(1:ifull,k),
     &                 r1/minl)
           eps(1:ifull,k)=max(eps(1:ifull,k),
     &                 r1/maxl,mineps)
           t_kh(1:ifull,k)=tke(1:ifull,k)*tke(1:ifull,k)
     &     /eps(1:ifull,k)*hdif*emi(1:ifull)
         end do
         call bounds(t_kh,nehalf=.true.)
         do k=1,kl
          do iq=1,ifull
           xfact(iq,k) = (t_kh(ie(iq),k)+t_kh(iq,k))*.5
           yfact(iq,k) = (t_kh(in(iq),k)+t_kh(iq,k))*.5
          end do
         end do
        else
         t_kh=0. ! no diffusion (i.e., for pure nvmix.eq.6)
                 ! Probably works better for grid scales that
                 ! are less than 500 m
        end if

       case DEFAULT
         write(6,*) "ERROR: Unknown option nhorjlm=",nhorjlm
         stop
      end select
       
      ! Calculate horizontal diffusion based on prognostic TKE
      ! This can be combined with the diffusion coefficents above
      ! so as to operate over a large range of grid length scales
      if (nvmix==6) then
        if (nhorx==1) then
          do k=1,kl
            shear(:,k)=2.*((dudx(:,k)*sx_fact)**2
     &                    +(dvdy(:,k)*sy_fact)**2+dwdz(:,k)**2)
     &              +(dudy(:,k)*sy_fact
     &               +dvdx(:,k)*sx_fact)**2
     &              +(dudz(:,k)+dwdx(:,k)*sx_fact)**2
     &              +(dvdz(:,k)+dwdy(:,k)*sy_fact)**2
          end do
        else if (nhorx>=7) then
          do k=1,kmax
            shear(:,k)=2.*((dudx(:,k)*sx_fact)**2
     &                    +(dvdy(:,k)*sy_fact)**2+dwdz(:,k)**2)
     &              +(dudy(:,k)*sy_fact
     &               +dvdx(:,k)*sx_fact)**2
     &              +(dudz(:,k)+dwdx(:,k)*sx_fact)**2
     &              +(dvdz(:,k)+dwdy(:,k)*sy_fact)**2
          end do
          do k=kmax+1,kl
            shear(:,k)=2.*(dudx(:,k)**2+dvdy(:,k)**2
     &                    +dwdz(:,k)**2)
     &              +(dudy(:,k)+dvdx(:,k))**2
     &              +(dudz(:,k)+dwdx(:,k))**2
     &              +(dvdz(:,k)+dwdy(:,k))**2
          end do
        else
          shear(:,:)=2.*(dudx(:,:)**2+dvdy(:,:)**2
     &            +dwdz(:,:)**2)
     &            +(dudy(:,:)+dvdx(:,:))**2
     &            +(dudz(:,:)+dwdx(:,:))**2
     &            +(dvdz(:,:)+dwdy(:,:))**2
        end if
      end if

      if (nhorx==1) then
         do k=1,kl
            do iq=1,ifull
               xfact(iq,k) = xfact(iq,k)*tx_fact(iq)
               yfact(iq,k) = yfact(iq,k)*ty_fact(iq)
            end do               !  iq loop
         end do    
      else if (nhorx>=7) then
         do k=1,kmax
            do iq=1,ifull
               xfact(iq,k) = xfact(iq,k)*tx_fact(iq)
               yfact(iq,k) = yfact(iq,k)*ty_fact(iq)
            end do               !  iq loop
         end do
      end if
      call boundsuv(xfact,yfact,stag=-9) ! MJT - can use stag=-9 option that will
                                         ! only update iwu and isv values

      do k=1,kl
        base(:,k)=emi+xfact(1:ifull,k)+xfact(iwu,k)
     &               +yfact(1:ifull,k)+yfact(isv,k)
      end do

      if(nhorps==0.or.nhorps==-2)then ! for nhorps=-1,-3 don't diffuse u,v
         do k=1,kl
            do iq=1,ifull
               ucc = ( uc(iq,k)*emi(iq) +
     &                 xfact(iq,k)*uc(ie(iq),k) +
     &                 xfact(iwu(iq),k)*uc(iw(iq),k) +
     &                 yfact(iq,k)*uc(in(iq),k) +
     &                 yfact(isv(iq),k)*uc(is(iq),k) ) /
     &              base(iq,k)
               vcc = ( vc(iq,k)*emi(iq) +
     &                 xfact(iq,k)*vc(ie(iq),k) +
     &                 xfact(iwu(iq),k)*vc(iw(iq),k) +
     &                 yfact(iq,k)*vc(in(iq),k) +
     &                 yfact(isv(iq),k)*vc(is(iq),k) ) /
     &               base(iq,k)
               wcc = ( wc(iq,k)*emi(iq) +
     &                 xfact(iq,k)*wc(ie(iq),k) +
     &                 xfact(iwu(iq),k)*wc(iw(iq),k) +
     &                 yfact(iq,k)*wc(in(iq),k) +
     &                 yfact(isv(iq),k)*wc(is(iq),k) ) /
     &               base(iq,k)
               u(iq,k) = ax(iq)*ucc + ay(iq)*vcc + az(iq)*wcc
               v(iq,k) = bx(iq)*ucc + by(iq)*vcc + bz(iq)*wcc
            enddo   !  iq loop
         end do
      endif   ! nhorps.ge.0

#ifdef debug
      if(diag.and.mydiag)then
          do k=1,kl
            write(6,*) 'k,id,jd,idjd ',k,id,jd,idjd
            write(6,*) 'k, xfact, xfactw ',k,xfact(idjd,k),
     &                                   xfact(iwu(idjd),k)
            write(6,*) 'k, yfact, yfacts ',k,yfact(idjd,k),
     &                                   yfact(isv(idjd),k)
            write(6,*) 'k, uc,uce,ucw,ucn,ucs '
     &        ,k,uc(idjd,k),uc(ie(idjd),k),uc(iw(idjd),k)
     &        ,uc(in(idjd),k),uc(is(idjd),k)
            write(6,*) 'k,u,v ',
     &            k,u(idjd,k),v(idjd,k)
          end do
      endif
#endif

      if(nhorps/=-2)then   ! for nhorps=-2 don't diffuse T, qg
c       do t diffusion based on potential temperature ff
        if(nhorps/=-3)then  ! for nhorps=-3 don't diffuse T; only qg
          do k=1,kl
            do iq=1,ifull
              ff(iq,k,1)=t(iq,k)/ptemp(iq) ! watch out for Chen!
              ff(iq,k,2)=qg(iq,k)
            enddo              !  iq loop
          end do
          call bounds(ff(:,:,1:2))
          do k=1,kl
            do iq=1,ifull
              t(iq,k)= ptemp(iq) *
     &                ( ff(iq,k,1)*emi(iq) +
     &                  xfact(iq,k)*ff(ie(iq),k,1) +
     &                  xfact(iwu(iq),k)*ff(iw(iq),k,1) +
     &                  yfact(iq,k)*ff(in(iq),k,1) +
     &                  yfact(isv(iq),k)*ff(is(iq),k,1) ) /
     &                base(iq,k)
              qg(iq,k) = ( ff(iq,k,2)*emi(iq) +
     &                  xfact(iq,k)*ff(ie(iq),k,2) +
     &                  xfact(iwu(iq),k)*ff(iw(iq),k,2) +
     &                  yfact(iq,k)*ff(in(iq),k,2) +
     &                  yfact(isv(iq),k)*ff(is(iq),k,2) ) /
     &                base(iq,k)
            end do              !  iq loop
          end do
        else
          do k=1,kl
            do iq=1,ifull
              ff(iq,k,1)=qg(iq,k)
            enddo              !  iq loop
          end do
          call bounds(ff(:,:,1:1))
          do k=1,kl
            do iq=1,ifull
              qg(iq,k) = ( ff(iq,k,1)*emi(iq) +
     &                     xfact(iq,k)*ff(ie(iq),k,1) +
     &                     xfact(iwu(iq),k)*ff(iw(iq),k,1) +
     &                     yfact(iq,k)*ff(in(iq),k,1) +
     &                     yfact(isv(iq),k)*ff(is(iq),k,1) ) /
     &                   base(iq,k)
             end do              !  iq loop
          end do
        endif                 ! (nhorps.ne.-3)

        ! cloud microphysics
        if (ldr/=0.and.ncloud/=0) then
          ff(1:ifull,:,1)=qlg(1:ifull,:)
          ff(1:ifull,:,2)=qfg(1:ifull,:)
          ff(1:ifull,:,3)=qrg(1:ifull,:)
          ff(1:ifull,:,4)=cffall(1:ifull,:)
          if (ncloud>=3) then
            ff(1:ifull,:,5)=stratcloud(1:ifull,:)
            call bounds(ff(:,:,1:5))
          else
            call bounds(ff(:,:,1:4))
          end if
          do k=1,kl
           do iq=1,ifull
            qlg(iq,k) = ( ff(iq,k,1)*emi(iq) +
     &        xfact(iq,k)*ff(ie(iq),k,1) +
     &        xfact(iwu(iq),k)*ff(iw(iq),k,1) +
     &        yfact(iq,k)*ff(in(iq),k,1) +
     &        yfact(isv(iq),k)*ff(is(iq),k,1) ) /
     &        base(iq,k)
            qfg(iq,k) = ( ff(iq,k,2)*emi(iq) +
     &        xfact(iq,k)*ff(ie(iq),k,2) +
     &        xfact(iwu(iq),k)*ff(iw(iq),k,2) +
     &        yfact(iq,k)*ff(in(iq),k,2) +
     &        yfact(isv(iq),k)*ff(is(iq),k,2) ) /
     &        base(iq,k)
            qrg(iq,k) = ( ff(iq,k,3)*emi(iq) +
     &        xfact(iq,k)*ff(ie(iq),k,3) +
     &        xfact(iwu(iq),k)*ff(iw(iq),k,3) +
     &        yfact(iq,k)*ff(in(iq),k,3) +
     &        yfact(isv(iq),k)*ff(is(iq),k,3) ) /
     &        base(iq,k)
            cffall(iq,k) = ( ff(iq,k,4)*emi(iq) +
     &        xfact(iq,k)*ff(ie(iq),k,4) +
     &        xfact(iwu(iq),k)*ff(iw(iq),k,4) +
     &        yfact(iq,k)*ff(in(iq),k,4) +
     &        yfact(isv(iq),k)*ff(is(iq),k,4) ) /
     &        base(iq,k)
           end do
          end do
          if (ncloud>=3) then
            do k=1,kl
              do iq=1,ifull
                stratcloud(iq,k) = ( ff(iq,k,5)*emi(iq) +
     &            xfact(iq,k)*ff(ie(iq),k,5) +
     &            xfact(iwu(iq),k)*ff(iw(iq),k,5) +
     &            yfact(iq,k)*ff(in(iq),k,5) +
     &            yfact(isv(iq),k)*ff(is(iq),k,5) ) /
     &           base(iq,k)
              end do
            end do
          end if
        end if                 ! (ldr.ne.0)
      endif                    ! (nhorps.ne.-2)

      ! apply horizontal diffusion to TKE and EPS terms
      if (nvmix==6) then
        ff(1:ifull,:,1)=tke(1:ifull,:)
        ff(1:ifull,:,2)=eps(1:ifull,:)
        call bounds(ff(:,:,1:2))
        do k=1,kl
          do iq=1,ifull
           tke(iq,k) = ( ff(iq,k,1)*emi(iq) +
     &                     xfact(iq,k)*ff(ie(iq),k,1) +
     &                     xfact(iwu(iq),k)*ff(iw(iq),k,1) +
     &                     yfact(iq,k)*ff(in(iq),k,1) +
     &                     yfact(isv(iq),k)*ff(is(iq),k,1) ) /
     &                   base(iq,k)
           eps(iq,k) = ( ff(iq,k,2)*emi(iq) +
     &                     xfact(iq,k)*ff(ie(iq),k,2) +
     &                     xfact(iwu(iq),k)*ff(iw(iq),k,2) +
     &                     yfact(iq,k)*ff(in(iq),k,2) +
     &                     yfact(isv(iq),k)*ff(is(iq),k,2) ) /
     &                   base(iq,k)
          end do
        end do
      end if
       
      ! prgnostic aerosols
      if (abs(iaero)==2) then
        ff(1:ifull,:,1:naero)=xtg(1:ifull,:,:)
        call bounds(ff(:,:,1:naero))
        do ntr=1,naero
          do k=1,kl
            do iq=1,ifull
              xtg(iq,k,ntr) =
     &          ( ff(iq,k,ntr)*emi(iq) +
     &          xfact(iq,k)*ff(ie(iq),k,ntr) +
     &          xfact(iwu(iq),k)*ff(iw(iq),k,ntr) +
     &          yfact(iq,k)*ff(in(iq),k,ntr) +
     &          yfact(isv(iq),k)*ff(is(iq),k,ntr) ) /
     &          base(iq,k)
            end do
          end do
        end do
      end if

      return
      end
