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
      use nharrs_m
      use nlin_m
      use parmhdff_m
      use sigs_m
      use tkeeps, only : tke,eps,shear,mintke,cm0,minl,maxl
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
 
      real, dimension(ifull+iextra,kl) :: uc, vc, wc, ee, ff, xfact,
     &                                    yfact, t_kh
      real, dimension(ifull) :: ptemp, tx_fact, ty_fact
      real, dimension(ifull) :: r1,r2,r3,r4
      real, dimension(ifull) :: z1,z2,z3,z4
      real, dimension(ifull+iextra,kl) :: ww
      real, dimension(ifull+iextra,0:kl-1) :: zgh
      real, dimension(ifull,kl) :: dudx,dudy,dvdx,dvdy
      real, dimension(ifull,kl) :: dudz,dvdz
      real, dimension(ifull,kl) :: dwdx,dwdy,dwdz
      integer, parameter :: nf=2
!     Local variables
      integer iq, k, nhora, nhorx, iaero, l
      real cc, delphi, emi, hdif, ucc, vcc, wcc
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

c     set up topography reduction factors for each type of location
c     expect power nf to be about 1 or 2 (see data statement)

      if (kmax.lt.0) then
        kmax=1
        do k=2,kl
          if (sig(k).ge.0.2) then
            kmax=k
          else
            exit
          end if
        end do
      end if

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
      ! This option is for prognostic TKE
      if (nhorjlm==0.or.nvmix==6) then
        ! Calculate du/dx,dv/dx,du/dy,dv/dy, etc 

        ! estimate height from geopotential at half levels
        zgh(1:ifull,0)=zs(1:ifull)/grav
        do k=1,kl-1
          zgh(1:ifull,k)=(ratha(k)*phi(:,k+1)+rathb(k)*phi(:,k))/grav
        end do
        call bounds(zgh)
       
        ! calculate horizontal gradients using POM stencil (see Shchepetkin and McWilliams 2003)
        call boundsuv(u,v,allvec=.true.)
        r1=0.
        r2=ratha(1)*u(ieu,2)+rathb(1)*u(ieu,1)
        r3=0.
        r4=ratha(1)*u(iwu,2)+rathb(1)*u(iwu,1)
        z1=zgh(ie,0)
        z2=zgh(ie,1)
        z3=zgh(iw,0)
        z4=zgh(iw,1)
        dudx(:,1)=0.25*((r1+r2-r3-r4)-(r2-r1+r4-r3)*(z1+z2-z3-z4)
     &            /max(z2-z1+z4-z3,1.E-8))*em(1:ifull)/ds      
        r2=ratha(1)*v(iev,2)+rathb(1)*v(iev,1)
        r4=ratha(1)*v(iwv,2)+rathb(1)*v(iwv,1)
        dvdx(:,1)=0.25*((r1+r2-r3-r4)-(r2-r1+r4-r3)*(z1+z2-z3-z4)
     &            /max(z2-z1+z4-z3,1.E-8))*em(1:ifull)/ds 
        r2=ratha(1)*u(inu,2)+rathb(1)*u(inu,1)
        r4=ratha(1)*u(isu,2)+rathb(1)*u(isu,1)
        z1=zgh(in,0)
        z2=zgh(in,1)
        z3=zgh(is,0)
        z4=zgh(is,1)
        dudy(:,1)=0.25*((r1+r2-r3-r4)-(r2-r1+r4-r3)*(z1+z2-z3-z4)
     &            /max(z2-z1+z4-z3,1.E-8))*em(1:ifull)/ds 
        r2=ratha(1)*v(inv,2)+rathb(1)*v(inv,1)
        r4=ratha(1)*v(isv,2)+rathb(1)*v(isv,1)
        dvdy(:,1)=0.25*((r1+r2-r3-r4)-(r2-r1+r4-r3)*(z1+z2-z3-z4)
     &            /max(z2-z1+z4-z3,1.E-8))*em(1:ifull)/ds 
        do k=2,kl-1
          r1=ratha(k-1)*u(ieu,k)+rathb(k-1)*u(ieu,k-1)
          r2=ratha(k)*u(ieu,k+1)+rathb(k)*u(ieu,k)
          r3=ratha(k-1)*u(iwu,k)+rathb(k-1)*u(iwu,k-1)
          r4=ratha(k)*u(iwu,k+1)+rathb(k)*u(iwu,k)
          z1=zgh(ie,k-1)
          z2=zgh(ie,k)
          z3=zgh(iw,k-1)
          z4=zgh(iw,k)
          dudx(:,k)=0.25*((r1+r2-r3-r4)-(r2-r1+r4-r3)*(z1+z2-z3-z4)
     &              /max(z2-z1+z4-z3,1.E-8))*em(1:ifull)/ds
          r1=ratha(k-1)*v(iev,k)+rathb(k-1)*v(iev,k-1)
          r2=ratha(k)*v(iev,k+1)+rathb(k)*v(iev,k)
          r3=ratha(k-1)*v(iwv,k)+rathb(k-1)*v(iwv,k-1)
          r4=ratha(k)*v(iwv,k+1)+rathb(k)*v(iwv,k)
          dvdx(:,k)=0.25*((r1+r2-r3-r4)-(r2-r1+r4-r3)*(z1+z2-z3-z4)
     &              /max(z2-z1+z4-z3,1.E-8))*em(1:ifull)/ds
          r1=ratha(k-1)*u(inu,k)+rathb(k-1)*u(inu,k-1)
          r2=ratha(k)*u(inu,k+1)+rathb(k)*u(inu,k)
          r3=ratha(k-1)*u(isu,k)+rathb(k-1)*u(isu,k-1)
          r4=ratha(k)*u(isu,k+1)+rathb(k)*u(isu,k)
          z1=zgh(in,k-1)
          z2=zgh(in,k)
          z3=zgh(is,k-1)
          z4=zgh(is,k)
          dudy(:,k)=0.25*((r1+r2-r3-r4)-(r2-r1+r4-r3)*(z1+z2-z3-z4)
     &              /max(z2-z1+z4-z3,1.E-8))*em(1:ifull)/ds
          r1=ratha(k-1)*v(inv,k)+rathb(k-1)*v(inv,k-1)
          r2=ratha(k)*v(inv,k+1)+rathb(k)*v(inv,k)
          r3=ratha(k-1)*v(isv,k)+rathb(k-1)*v(isv,k-1)
          r4=ratha(k)*v(isv,k+1)+rathb(k)*v(isv,k)
          dvdy(:,k)=0.25*((r1+r2-r3-r4)-(r2-r1+r4-r3)*(z1+z2-z3-z4)
     &              /max(z2-z1+z4-z3,1.E-8))*em(1:ifull)/ds
        end do
        dudx(:,kl)=dudx(:,kl-1) ! just a simple approximiation
        dvdx(:,kl)=dvdx(:,kl-1)
        dudy(:,kl)=dudy(:,kl-1)
        dvdy(:,kl)=dvdy(:,kl-1)

        ! calculate vertical velocity in m/s
        do k=1,kl        
          ! omega=ps*dpsldt
          ww(1:ifull,k)=(dpsldt(:,k)/sig(k)-dpsdt/(860.*ps(1:ifull)))
     &        *(-rdry/grav)*t(1:ifull,k)*(1.+0.61*qg(1:ifull,k)
     &        -qlg(1:ifull,k)-qfg(1:ifull,k))
        end do
        call bounds(ww)
        r1=0.
        r2=ratha(1)*ww(ie,2)+rathb(1)*ww(ie,1)
        r3=0.
        r4=ratha(1)*ww(iw,2)+rathb(1)*ww(iw,1)
        z1=zgh(ie,0)
        z2=zgh(ie,1)
        z3=zgh(iw,0)
        z4=zgh(iw,1)
        dwdx(:,1)=0.25*((r1+r2-r3-r4)-(r2-r1+r4-r3)*(z1+z2-z3-z4)
     &            /max(z2-z1+z4-z3,1.E-8))*em(1:ifull)/ds
        r2=ratha(1)*ww(in,2)+rathb(1)*ww(in,1)
        r4=ratha(1)*ww(is,2)+rathb(1)*ww(is,1)
        z1=zgh(in,0)
        z2=zgh(in,1)
        z3=zgh(is,0)
        z4=zgh(is,1)
        dwdy(:,1)=0.25*((r1+r2-r3-r4)-(r2-r1+r4-r3)*(z1+z2-z3-z4)
     &            /max(z2-z1+z4-z3,1.E-8))*em(1:ifull)/ds
        do k=2,kl-1
          r1=ratha(k-1)*ww(ie,k)+rathb(k-1)*ww(ie,k-1)
          r2=ratha(k)*ww(ie,k+1)+rathb(k)*ww(ie,k)
          r3=ratha(k-1)*ww(iw,k)+rathb(k-1)*ww(iw,k-1)
          r4=ratha(k)*ww(iw,k+1)+rathb(k)*ww(iw,k)
          z1=zgh(ie,k-1)
          z2=zgh(ie,k)
          z3=zgh(iw,k-1)
          z4=zgh(iw,k)
          dwdx(:,k)=0.25*((r1+r2-r3-r4)-(r2-r1+r4-r3)*(z1+z2-z3-z4)
     &              /max(z2-z1+z4-z3,1.E-8))*em(1:ifull)/ds
        
          r1=ratha(k-1)*ww(in,k)+rathb(k-1)*ww(in,k-1)
          r2=ratha(k)*ww(in,k+1)+rathb(k)*ww(in,k)
          r3=ratha(k-1)*ww(is,k)+rathb(k-1)*ww(is,k-1)
          r4=ratha(k)*ww(is,k+1)+rathb(k)*ww(is,k)
          z1=zgh(in,k-1)
          z2=zgh(in,k)
          z3=zgh(is,k-1)
          z4=zgh(is,k)
          dwdy(:,k)=0.25*((r1+r2-r3-r4)-(r2-r1+r4-r3)*(z1+z2-z3-z4)
     &              /max(z2-z1+z4-z3,1.E-8))*em(1:ifull)/ds
        end do
        dwdx(:,kl)=dwdx(:,kl-1)
        dwdy(:,kl)=dwdy(:,kl-1)
        
        ! calculate vertical gradients
        r1=0.
        r2=ratha(1)*u(1:ifull,2)+rathb(1)*u(1:ifull,1)          
        dudz(:,1)=(r2-r1)/(zgh(1:ifull,1)-zgh(1:ifull,0))
        r2=ratha(1)*v(1:ifull,2)+rathb(1)*v(1:ifull,1)          
        dvdz(:,1)=(r2-r1)/(zgh(1:ifull,1)-zgh(1:ifull,0))          
        r2=ratha(1)*ww(1:ifull,2)+rathb(1)*ww(1:ifull,1)          
        dwdz(:,1)=(r2-r1)/(zgh(1:ifull,1)-zgh(1:ifull,0))
        do k=2,kl-1
          r1=ratha(k-1)*u(1:ifull,k)+rathb(k-1)*u(1:ifull,k-1)
          r2=ratha(k)*u(1:ifull,k+1)+rathb(k)*u(1:ifull,k)          
          dudz(:,k)=(r2-r1)/(zgh(1:ifull,k)-zgh(1:ifull,k-1))
          r1=ratha(k-1)*v(1:ifull,k)+rathb(k-1)*v(1:ifull,k-1)
          r2=ratha(k)*v(1:ifull,k+1)+rathb(k)*v(1:ifull,k)          
          dvdz(:,k)=(r2-r1)/(zgh(1:ifull,k)-zgh(1:ifull,k-1))          
          r1=ratha(k-1)*ww(1:ifull,k)+rathb(k-1)*ww(1:ifull,k-1)
          r2=ratha(k)*ww(1:ifull,k+1)+rathb(k)*ww(1:ifull,k)          
          dwdz(:,k)=(r2-r1)/(zgh(1:ifull,k)-zgh(1:ifull,k-1))
        end do
        dudz(:,kl)=0.
        dvdz(:,kl)=0.
        dudz(:,kl)=0.

      end if   ! nhorjlm==0.or.nvmix==6
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
       case(0)
         ! This is based on 2D Smagorinsky closure
         do k=1,kl
           hdif=dt*hdiff(k) ! N.B.  hdiff(k)=khdif*.1
           do iq=1,ifull
             cc=(dudx(iq,k)-dvdy(iq,k))**2
             cc=cc+(dudy(iq,k)+dvdx(iq,k))**2
             t_kh(iq,k)=sqrt(cc)*hdif/(em(iq)*em(iq))  ! this one with em in D terms
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

       case(3)
         t_kh=0. ! no diffusion (i.e., for pure nvmix.eq.6)
                 ! Probably works better for grid scales that
                 ! are less than 500 m

       case DEFAULT
         write(6,*) "ERROR: Unknown option nhorjlm=",nhorjlm
         stop
      end select
       
      ! Calculate horizontal diffusion based on prognostic TKE
      ! This can be combined with the diffusion coefficents above
      ! so as to operate over a large range of grid length scales
      if (nvmix.eq.6) then
        tke(1:ifull,:)=max(tke(1:ifull,:),mintke)
        eps(1:ifull,:)=min(eps(1:ifull,:),
     &                 (cm0**0.75)*tke(1:ifull,:)
     &                 *sqrt(tke(1:ifull,:))/minl)
        eps(1:ifull,:)=max(eps(1:ifull,:),
     &                 (cm0**0.75)*tke(1:ifull,:)
     &                 *sqrt(tke(1:ifull,:))/maxl)
        hdif=dt*cm0/(ds*ds)
        do k=1,kl
          t_kh(1:ifull,k)=max(tke(1:ifull,k)*tke(1:ifull,k)
     &    /eps(1:ifull,k)*hdif,t_kh(1:ifull,k))

          shear(:,k)=2.*(dudx(:,k)**2+dvdy(:,k)**2+dwdz(:,k)**2)
     &              +(dudy(:,k)+dvdx(:,k))**2
     &              +(dudz(:,k)+dwdx(:,k))**2
     &              +(dvdz(:,k)+dwdy(:,k))**2
        end do
      end if

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
       
      if(nhorps.ne.-2)then   ! for nhorps=-2 don't diffuse T, qg
c       do t diffusion based on potential temperature ff
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
     &                ( ff(iq,k)*emi +
     &                  xfact(iq,k)*ff(ie(iq),k) +
     &                  xfact(iwu(iq),k)*ff(iw(iq),k) +
     &                  yfact(iq,k)*ff(in(iq),k) +
     &                  yfact(isv(iq),k)*ff(is(iq),k) ) /
     &                ( emi + xfact(iq,k) + xfact(iwu(iq),k) +
     &                  yfact(iq,k) + yfact(isv(iq),k) )
            enddo           !  iq loop
          end do
        endif                 ! (nhorps.ne.-3)
        do k=1,kl
          do iq=1,ifull
            emi=1./em(iq)**2
            qg(iq,k) = ( ee(iq,k)*emi +
     &                   xfact(iq,k)*ee(ie(iq),k) +
     &                   xfact(iwu(iq),k)*ee(iw(iq),k) +
     &                   yfact(iq,k)*ee(in(iq),k) +
     &                   yfact(isv(iq),k)*ee(is(iq),k) ) /
     &                 ( emi + xfact(iq,k) + xfact(iwu(iq),k) +
     &                   yfact(iq,k)+yfact(isv(iq),k))
           end do              !  iq loop
        end do
        ! cloud microphysics
        ee=qlg(1:ifull,:)
        call bounds(ee)
        do k=1,kl
          qlg(1:ifull,k) = ( ee(1:ifull,k)/em(1:ifull)**2 +
     &      xfact(1:ifull,k)*ee(ie,k) +
     &      xfact(iwu,k)*ee(iw,k) +
     &      yfact(1:ifull,k)*ee(in,k) +
     &      yfact(isv,k)*ee(is,k) ) /
     &      ( 1./em(1:ifull)**2 +
     &        xfact(1:ifull,k) + xfact(iwu,k) +
     &        yfact(1:ifull,k) + yfact(isv,k) )
        end do
        ee=qfg(1:ifull,:)
        call bounds(ee)
        do k=1,kl
          qfg(1:ifull,k) = ( ee(1:ifull,k)/em(1:ifull)**2 +
     &      xfact(1:ifull,k)*ee(ie,k) +
     &      xfact(iwu,k)*ee(iw,k) +
     &      yfact(1:ifull,k)*ee(in,k) +
     &      yfact(isv,k)*ee(is,k) ) /
     &      ( 1./em(1:ifull)**2 +
     &        xfact(1:ifull,k) + xfact(iwu,k) +
     &        yfact(1:ifull,k) + yfact(isv,k) )
        end do
        ee=qrg(1:ifull,:)
        call bounds(ee)
        do k=1,kl
          qrg(1:ifull,k) = ( ee(1:ifull,k)/em(1:ifull)**2 +
     &      xfact(1:ifull,k)*ee(ie,k) +
     &      xfact(iwu,k)*ee(iw,k) +
     &      yfact(1:ifull,k)*ee(in,k) +
     &      yfact(isv,k)*ee(is,k) ) /
     &      ( 1./em(1:ifull)**2 +
     &        xfact(1:ifull,k) + xfact(iwu,k) +
     &        yfact(1:ifull,k) + yfact(isv,k) )
        end do
      endif                    ! (nhorps.ne.-2)
       
      ! MJT aerosols
      if (abs(iaero).eq.2) then
        do l=1,naero
          ee(1:ifull,:)=xtg(1:ifull,:,l)
          call bounds(ee)
          do k=1,kl
            do iq=1,ifull
              emi=1./em(iq)**2
              xtg(iq,k,l) = ( ee(iq,k)*emi +
     &                    xfact(iq,k)*ee(ie(iq),k) +
     &                    xfact(iwu(iq),k)*ee(iw(iq),k) +
     &                    yfact(iq,k)*ee(in(iq),k) +
     &                    yfact(isv(iq),k)*ee(is(iq),k) ) /
     &                  ( emi + xfact(iq,k) + xfact(iwu(iq),k) +
     &                    yfact(iq,k)+yfact(isv(iq),k))
            enddo           !  iq loop
          end do
        end do
      end if

      return
      end
