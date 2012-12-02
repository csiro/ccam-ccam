      subroutine mslp(pmsl,psl,zs,t)
      use cc_mpi, only : mydiag
      use sigs_m
      implicit none
!     this one will ignore negative zs (i.e. over the ocean)
      integer, parameter :: meth=1 ! 0 for original, 1 for other jlm - always now
      include 'newmpar.h'
      include 'const_phys.h'
      include 'parm.h'
      real, dimension(ifull), intent(out) :: pmsl
      real, dimension(ifull), intent(in) :: psl,zs
      real, dimension(ifull,kl), intent(in) :: t
      integer, save :: lev = -1
      real c,conr,con
      real, dimension(ifull) :: phi1,tsurf,tav,dlnps
      
      c=grav/stdlapse
      conr=c/rdry
      if (lev<0) then
        lev=1
        do while (sig(lev+1)<=0.9)
          lev=lev+1
        end do
      end if
      con=sig(lev)**(rdry/c)/c
      
c     if(meth.eq.0)then
c       do iq=1,ifull
c        pmsl(iq)=ps(iq)*(1.+con*zs(iq)/t(iq,lev))**conr
c       enddo
c     endif  ! (meth.eq.0)

      if(meth==1)then
         phi1(:)=t(:,lev)*rdry*(1.-sig(lev))/sig(lev) ! phi of sig(lev) above sfce
         tsurf(:)=t(:,lev)+phi1(:)*stdlapse/grav
         !tav(:)=tsurf(:)+max(0.,zs(1:ifull))*.5*stdlapse/grav
         !dlnps(:)=max(0.,zs(1:ifull))/(rdry*tav)
         ! MJT suggestion
         tav(:)=tsurf(:)+zs(1:ifull)*.5*stdlapse/grav
         dlnps(:)=zs(1:ifull)/(rdry*tav(:))
         pmsl(:)=1.e5*exp(psl(:)+dlnps(:))
      endif  ! (meth.eq.1)
      
      if(nmaxpr==1.and.mydiag)then
        print *,'meth,lev,sig(lev) ',meth,lev,sig(lev)
        print *,'zs,t_lev,psl,pmsl ',
     .           zs(idjd),t(idjd,lev),psl(idjd),pmsl(idjd)
      endif
      
      return
      end
      
      subroutine to_psl(pmsl,psl,zs,t)
      use cc_mpi, only : mydiag
      use sigs_m
      implicit none
      include 'newmpar.h'
      include 'const_phys.h'
      include 'parm.h'
      real, dimension(ifull), intent(out) :: psl
      real, dimension(ifull), intent(in) :: pmsl,zs
      real, dimension(ifull,kl), intent(in) :: t
      integer, save :: lev = -1
      real c,conr,con
      real, dimension(ifull) :: dlnps,phi1,tsurf,tav

      c=grav/stdlapse
      conr=c/rdry
      if (lev<0) then
        lev=1
        do while (sig(lev+1)<=0.9)
          lev=lev+1
        end do
      end if
      con=sig(lev)**(rdry/c)/c

      phi1=t(:,lev)*rdry*(1.-sig(lev))/sig(lev) ! phi of sig(lev) above sfce
      tsurf=t(:,lev)+phi1*stdlapse/grav
      !tav=tsurf+max(0.,zs(1:ifull))*.5*stdlapse/grav
      !dlnps=max(0.,zs(1:ifull))/(rdry*tav)
      tav=tsurf+zs(1:ifull)*.5*stdlapse/grav
      dlnps=zs(1:ifull)/(rdry*tav)
      psl=log(1.e-5*pmsl) -dlnps

      if(nmaxpr==1.and.mydiag)then
        print *,'to_psl lev,sig(lev) ',lev,sig(lev)
        print *,'zs,t_lev,psl,pmsl ',
     .           zs(idjd),t(idjd,lev),psl(idjd),pmsl(idjd)
      endif

      return
      end
