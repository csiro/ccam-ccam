      subroutine staguv(u,v,uout,vout)
!     stripped down version just with nstag=nstagu=-3, mstagpt=-3      
!     staguv    may be called from adjust5, upglobal
!     unstaguv  may be called from adjust5,  nonlin
c     nstag now in parm.h  ! for nstag=0   staguv uses cubic interpolation
c                          ! for nstag=1   staguv goes via cartesian components
c                          ! for nstag=2   staguv linear with em terms
c                          ! for nstag=3-5 staguv: jlm reversible 2/1/98
c                          ! for nstag=7   staguv: jlm reversible Akima
c                          ! for nstag=10  staguv: jlm rev_stag
c                          !    only -3 used nowadays
c     nstagu now in parm.h ! same but for unstaguv
c     N.B. staguv & unstaguv previously required 2D arrays as input
c     - no longer, as transferred here to uin and vin
c     assume k level already given
      implicit none
      include 'newmpar.h'
      include 'indices.h' ! in,is,iw,ie,inn,iss,iww,iee
      include 'map.h'
      include 'parm.h'
      include 'parmdyn.h'
      include 'vecsuv.h'   ! vecsuv info
      include 'vecsuva.h'  ! vecsuva info
      real u(ifull,kl),v(ifull,kl),uout(ifull,kl),vout(ifull,kl)
      real ua(ifull,2*kl),va(ifull,kl),ud(ifull,2*kl),vd(ifull,kl) ! work arrays
      real ub(ifull,2*kl),vb(ifull,kl)                       ! work arrays
      equivalence (ua(1,kl+1),va),(ud(1,kl+1),vd)    ! to ensure contiguous
      equivalence (ub(1,kl+1),vb)                     ! to ensure contiguous
      real uin(ifull,2*kl),vin(ifull,kl)                     ! work arrays
      equivalence (uin(1,kl+1),vin)                   ! to ensure contiguous
      integer, parameter :: itnmax=3
      integer :: i, iq, itn, j, k

c     unstaggered u & v as input; staggered as output
      do k=1,kl
         do iq=1,ifull
            uin(iq,k) = u(iq,k)
            vin(iq,k) = v(iq,k)
         end do
      end do

!     if(nstag.eq.-3)then
!       call reva6(uin,vin,uout,vout,ua,va,ub,vb,ud,vd)
!     endif    !  (nstag.eq.-3)
!cdir nodep
c       do iq=1,ifull   ! precalculate rhs terms
c         ud(iwu(iq))= uin(iwu(iq))/2.+uin(iq)+uin(ieu(iq))/10.
c         vd(isv(iq))= vin(isv(iq))/2.+vin(iq)+vin(inv(iq))/10.
c       enddo
      do k=1,kl
         do iq=2,ifull-1
            ud(iq-1,k) = uin(iq-1,k)/2. + uin(iq,k) + uin(iq+1,k)/10.
         end do
         do iq=il+1,ifull-il
          vd(iq-il,k) = vin(iq-il,k)/2. + vin(iq,k) + vin(iq+il,k)/10.
         end do
         do j=1,jl
            iq=1+(j-1)*il
            ud(iwu(iq),k) = uin(iwu(iq),k)/2. + uin(iq,k) +
     &                       uin(ieu(iq),k)/10.
            iq=il+(j-1)*il
            ud(iwu(iq),k) = uin(iwu(iq),k)/2. + uin(iq,k) +
     &                       uin(ieu(iq),k)/10.
         end do
         do j=1,jl,2*il
            do i=1,il
               iq=i+(j-1)*il
               vd(isv(iq),k) = vin(isv(iq),k)/2. + vin(iq,k) +
     &                          vin(inv(iq),k)/10.
            end do
         end do
         do j=2*il,jl,2*il
            do i=1,il
               iq=i+(j-1)*il
               vd(isv(iq),k) = vin(isv(iq),k)/2. + vin(iq,k) +
     &                          vin(inv(iq),k)/10.
            end do
         end do
      end do ! k
c        do iq=1,ifull
c          ua(iq)=ud(iq)-ud(ieu(iq))/2.   ! 1st guess
c          va(iq)=vd(iq)-vd(inv(iq))/2.   ! 1st guess
c        enddo
!cdir nodep

      do k=1,kl
         do iq=1,ifull-1
            ua(iq,k) = ud(iq,k) - ud(iq+1,k)/2.
         end do
         do iq=1,ifull-il
            va(iq,k) = vd(iq,k) - vd(iq+il,k)/2.
         end do
         do j=1,jl
            iq=il+(j-1)*il
            ua(iq,k) = ud(iq,k) - ud(ieu(iq),k)/2. ! 1st guess
         end do
         do j=2*il,jl,2*il
            do i=1,il
               iq=i+(j-1)*il
               va(iq,k) = vd(iq,k) - vd(inv(iq),k)/2. ! 1st guess
            end do
         end do
      end do  ! k loop

      do itn=1,itnmax           ! each loop is a double iteration
!cdir nodep
c         do iq=1,ifull
c          ub(iq)=ua(ieu(iq))
c          vb(iq)=va(inv(iq))
c         end do
         do k=1,kl
            do iq=1,ifull-1
               ub(iq,k) = ua(iq+1,k)
            end do
            do iq=1,ifull-il
               vb(iq,k) = va(iq+il,k)
            end do
            do j=1,jl
               iq=il+(j-1)*il
               ub(iq,k) = ua(ieu(iq),k)
            end do
            do j=2*il,jl,2*il
               do i=1,il
                  iq=i+(j-1)*il
                  vb(iq,k) = va(inv(iq),k)
               end do
            end do
         end do
!cdir nodep
c       do iq=1,ifull
c         uin(iq)=(ud(iq)-.5*ud(ieu(iq))
c    &                 -ua(iwu(iq))/10. +ub(ieu(iq))/4.)/.95
c         vin(iq)=(vd(iq)-.5*vd(inv(iq))
c    &                 -va(isv(iq))/10. +vb(inv(iq))/4.)/.95
c        end do
         do k=1,kl
            do iq=2,ifull-1
               uin(iq,k) = (ud(iq,k)-.5*ud(iq+1,k)
     &                 -ua(iq-1,k)/10. +ub(iq+1,k)/4.)/.95
            end do
            do iq=il+1,ifull-il
               vin(iq,k) = (vd(iq,k)-.5*vd(iq+il,k)
     &                -va(iq-il,k)/10. +vb(iq+il,k)/4.)/.95
            end do
            do j=1,jl
               iq=1+(j-1)*il
               uin(iq,k) = (ud(iq,k)-.5*ud(ieu(iq),k)
     &                 -ua(iwu(iq),k)/10. +ub(ieu(iq),k)/4.)/.95
               iq=il+(j-1)*il
               uin(iq,k) = (ud(iq,k)-.5*ud(ieu(iq),k)
     &                 -ua(iwu(iq),k)/10. +ub(ieu(iq),k)/4.)/.95
            end do
            do j=1,jl,2*il
               do i=1,il
                  iq=i+(j-1)*il
                  vin(iq,k) = (vd(iq,k)-.5*vd(inv(iq),k)
     &                   -va(isv(iq),k)/10. +vb(inv(iq),k)/4.)/.95
               end do
            end do
            do j=2*il,jl,2*il
               do i=1,il
                  iq=i+(j-1)*il
                  vin(iq,k) = (vd(iq,k)-.5*vd(inv(iq),k)
     &                   -va(isv(iq),k)/10. +vb(inv(iq),k)/4.)/.95
               end do
            end do
         end do                 ! k

!cdir nodep
c         do iq=1,ifull
c          ub(iq)=uin(ieu(iq))
c          vb(iq)=vin(inv(iq))
c         end do
         do k=1,kl
            do iq=1,ifull-1
               ub(iq,k) = uin(iq+1,k)
            end do
            do iq=1,ifull-il
               vb(iq,k) = vin(iq+il,k)
            end do
            do j=1,jl
               iq=il+(j-1)*il
               ub(iq,k) = uin(ieu(iq),k)
            end do
            do j=2*il,jl,2*il
               do i=1,il
                  iq=i+(j-1)*il
                  vb(iq,k) = vin(inv(iq),k)
               end do
            end do
         end do ! k
!cdir nodep
c        do iq=1,ifull
c         ua(iq)=(ud(iq)-.5*ud(ieu(iq))
c    &                 -uin(iwu(iq))/10. +ub(ieu(iq))/4.)/.95
c         va(iq)=(vd(iq)-.5*vd(inv(iq))
c    &                 -vin(isv(iq))/10. +vb(inv(iq))/4.)/.95
c        end do
         do k=1,kl
            do iq=2,ifull-1
               ua(iq,k) = (ud(iq,k)-.5*ud(iq+1,k)
     &                 -uin(iq-1,k)/10. +ub(iq+1,k)/4.)/.95
            end do
            do iq=il+1,ifull-il
               va(iq,k) = (vd(iq,k)-.5*vd(iq+il,k)
     &                 -vin(iq-il,k)/10. +vb(iq+il,k)/4.)/.95
            end do
            do j=1,jl
               iq=1+(j-1)*il
               ua(iq,k) = (ud(iq,k)-.5*ud(ieu(iq),k)
     &                 -uin(iwu(iq),k)/10. +ub(ieu(iq),k)/4.)/.95
               iq=il+(j-1)*il
               ua(iq,k) = (ud(iq,k)-.5*ud(ieu(iq),k)
     &                 -uin(iwu(iq),k)/10. +ub(ieu(iq),k)/4.)/.95
            end do
            do j=1,jl,2*il
               do i=1,il
                  iq=i+(j-1)*il
                  va(iq,k) = (vd(iq,k)-.5*vd(inv(iq),k)
     &                 -vin(isv(iq),k)/10. +vb(inv(iq),k)/4.)/.95
               end do
            end do
            do j=2*il,jl,2*il
               do i=1,il
                  iq=i+(j-1)*il
                  va(iq,k) = (vd(iq,k)-.5*vd(inv(iq),k)
     &                 -vin(isv(iq),k)/10. +vb(inv(iq),k)/4.)/.95
               end do
            end do
         end do ! k
      end do                    ! itn=1,itnmax
	 
      do k=1,kl
         do iq=1,ifull          ! final values for output
            uout(iq,k) = ua(iq,k)
            vout(iq,k) = va(iq,k)
         end do
      end do

      return

      end


      subroutine unstaguv(u,v,uout,vout)
c     staggered u & v as input; unstaggered as output
      implicit none
      include 'newmpar.h'
      include 'indices.h' ! in,is,iw,ie,inn,iss,iww,iee
      include 'map.h'
      include 'parm.h'
      include 'parmdyn.h'
      include 'vecsuv.h'   ! vecsuv info
      include 'vecsuva.h'  ! vecsuva info
      real u(ifull,kl),v(ifull,kl),uout(ifull,kl),vout(ifull,kl)
      real ua(ifull,2*kl),va(ifull,kl),ud(ifull,2*kl),vd(ifull,kl) ! work arrays
      real ub(ifull,2*kl),vb(ifull,kl)                       ! work arrays
      equivalence (ua(1,kl+1),va),(ud(1,kl+1),vd)    ! to ensure contiguous
      equivalence (ub(1,kl+1),vb)                     ! to ensure contiguous
      real uin(ifull,2*kl),vin(ifull,kl)                     ! work arrays
      equivalence (uin(1,kl+1),vin)                   ! to ensure contiguous
      integer, parameter :: itnmax=3
      integer :: i, iq, itn, j, k

      do k=1,kl
         do iq=1,ifull
            uin(iq,k) = u(iq,k)
            vin(iq,k) = v(iq,k)
         end do
      end do

!     if(nstagu.eq.-3)then
!       call revb6(uin,vin,uout,vout,ua,va,ub,vb,ud,vd)
!     endif    !  (nstagu.eq.-3)
 
!cdir nodep
c       do iq=1,ifull   ! precalculate rhs terms
c         ud(ieu(iq))= uin(ieu(iq))/2.+uin(iq)+uin(iwu(iq))/10.
c         vd(inv(iq))= vin(inv(iq))/2.+vin(iq)+vin(isv(iq))/10.
c       end do
      do k=1,kl
         do iq=2,ifull-1
            ud(iq+1,k) = uin(iq+1,k)/2. + uin(iq,k) + uin(iq-1,k)/10.
         end do
         do iq=il+1,ifull-il
            vd(iq+il,k) = vin(iq+il,k)/2. + vin(iq,k) + vin(iq-il,k)/10.
         end do
         do j=1,jl
            iq=1+(j-1)*il
            ud(ieu(iq),k) = uin(ieu(iq),k)/2. + uin(iq,k) +
     &                      uin(iwu(iq),k)/10.
            iq=il+(j-1)*il
            ud(ieu(iq),k) = uin(ieu(iq),k)/2. + uin(iq,k) +
     &                      uin(iwu(iq),k)/10.
         end do
         do j=1,jl,2*il
            do i=1,il
               iq=i+(j-1)*il
               vd(inv(iq),k) = vin(inv(iq),k)/2. + vin(iq,k) +
     &                         vin(isv(iq),k)/10.
            end do
         end do
         do j=2*il,jl,2*il
            do i=1,il
               iq=i+(j-1)*il
               vd(inv(iq),k) = vin(inv(iq),k)/2. + vin(iq,k) +
     &                         vin(isv(iq),k)/10.
            end do
         end do
      end do ! k
c        do iq=1,ifull
c          ua(iq)=ud(iq)-ud(iwu(iq))/2.   ! 1st guess
c          va(iq)=vd(iq)-vd(isv(iq))/2.   ! 1st guess
c        end do
!cdir nodep
      do k=1,kl
         do iq=2,ifull
            ua(iq,k) = ud(iq,k)-ud(iq-1,k)/2.
         end do
         do iq=il+1,ifull
            va(iq,k) = vd(iq,k)-vd(iq-il,k)/2.
         end do
         do j=1,jl
            iq=1+(j-1)*il
            ua(iq,k) = ud(iq,k)-ud(iwu(iq),k)/2. ! 1st guess
         end do
         do j=1,jl,2*il
            do i=1,il
               iq=i+(j-1)*il
               va(iq,k) = vd(iq,k)-vd(isv(iq),k)/2. ! 1st guess
            end do
         end do
      end do ! k

      do itn=1,itnmax           ! each loop is a double iteration
!cdir nodep
c         do iq=1,ifull
c          ub(iq)=ua(iwu(iq))
c          vb(iq)=va(isv(iq))
c         end do
!cdir nodep
         do k=1,kl
            do iq=2,ifull
               ub(iq,k) = ua(iq-1,k)
            end do
            do iq=il+1,ifull
               vb(iq,k) = va(iq-il,k)
            end do
            do j=1,jl
               iq=1+(j-1)*il
               ub(iq,k) = ua(iwu(iq),k)
            end do
            do j=1,jl,2*il
               do i=1,il
                  iq=i+(j-1)*il
                  vb(iq,k) = va(isv(iq),k)
               end do
            end do
         end do ! k
!cdir nodep
c        do iq=1,ifull
c         uin(iq)=(ud(iq)-.5*ud(iwu(iq))
c    &                 -ua(ieu(iq))/10. +ub(iwu(iq))/4.)/.95
c         vin(iq)=(vd(iq)-.5*vd(isv(iq))
c    &                 -va(inv(iq))/10. +vb(isv(iq))/4.)/.95
c        end do
         do k=1,kl
            do iq=2,ifull-1
               uin(iq,k) = (ud(iq,k)-.5*ud(iq-1,k)
     &                 -ua(iq+1,k)/10. +ub(iq-1,k)/4.)/.95
            end do
            do iq=il+1,ifull-il
               vin(iq,k) = (vd(iq,k)-.5*vd(iq-il,k)
     &                 -va(iq+il,k)/10. +vb(iq-il,k)/4.)/.95
            end do
            do j=1,jl
               iq=1+(j-1)*il
               uin(iq,k) = (ud(iq,k)-.5*ud(iwu(iq),k)
     &                 -ua(ieu(iq),k)/10. +ub(iwu(iq),k)/4.)/.95
               iq=il+(j-1)*il
               uin(iq,k) = (ud(iq,k)-.5*ud(iwu(iq),k)
     &                 -ua(ieu(iq),k)/10. +ub(iwu(iq),k)/4.)/.95
            end do
            do j=1,jl,2*il
               do i=1,il
                  iq=i+(j-1)*il
                  vin(iq,k) = (vd(iq,k)-.5*vd(isv(iq),k)
     &                 -va(inv(iq),k)/10. +vb(isv(iq),k)/4.)/.95
               end do
            end do
            do j=2*il,jl,2*il
               do i=1,il
                  iq=i+(j-1)*il
                  vin(iq,k) = (vd(iq,k)-.5*vd(isv(iq),k)
     &                 -va(inv(iq),k)/10. +vb(isv(iq),k)/4.)/.95
               end do
            end do
         end do ! k

!cdir nodep
c         do iq=1,ifull
c          ub(iq)=uin(iwu(iq))
c          vb(iq)=vin(isv(iq))
c         end do
!cdir nodep
         do k=1,kl
            do iq=2,ifull
               ub(iq,k) = uin(iq-1,k)
            end do
            do iq=il+1,ifull
               vb(iq,k) = vin(iq-il,k)
            end do
            do j=1,jl
               iq=1+(j-1)*il
               ub(iq,k) = uin(iwu(iq),k)
            end do
            do j=1,jl,2*il
               do i=1,il
                  iq=i+(j-1)*il
                  vb(iq,k) = vin(isv(iq),k)
               end do
            end do
         end do ! k
!cdir nodep
c        do iq=1,ifull
c         ua(iq)=(ud(iq)-.5*ud(iwu(iq))
c    &                 -uin(ieu(iq))/10. +ub(iwu(iq))/4.)/.95
c         va(iq)=(vd(iq)-.5*vd(isv(iq))
c    &                 -vin(inv(iq))/10. +vb(isv(iq))/4.)/.95
c        end do
         do k=1,kl
            do iq=2,ifull-1
               ua(iq,k) = (ud(iq,k)-.5*ud(iq-1,k)
     &                 -uin(iq+1,k)/10. +ub(iq-1,k)/4.)/.95
            end do
            do iq=il+1,ifull-il
               va(iq,k) = (vd(iq,k)-.5*vd(iq-il,k)
     &              -vin(iq+il,k)/10. +vb(iq-il,k)/4.)/.95
            end do
            do j=1,jl
               iq=1+(j-1)*il
               ua(iq,k) = (ud(iq,k)-.5*ud(iwu(iq),k)
     &                 -uin(ieu(iq),k)/10. +ub(iwu(iq),k)/4.)/.95
               iq=il+(j-1)*il
               ua(iq,k) = (ud(iq,k)-.5*ud(iwu(iq),k)
     &                 -uin(ieu(iq),k)/10. +ub(iwu(iq),k)/4.)/.95
            end do
            do j=1,jl,2*il
               do i=1,il
                  iq=i+(j-1)*il
                  va(iq,k) = (vd(iq,k)-.5*vd(isv(iq),k)
     &                 -vin(inv(iq),k)/10. +vb(isv(iq),k)/4.)/.95
               end do
            end do
            do j=2*il,jl,2*il
               do i=1,il
                  iq=i+(j-1)*il
                  va(iq,k) = (vd(iq,k)-.5*vd(isv(iq),k)
     &                 -vin(inv(iq),k)/10. +vb(isv(iq),k)/4.)/.95
               end do
            end do
         end do ! k
      end do                    ! itn=1,itnmax
      
      do k=1,kl
         do iq=1,ifull          ! final values for output
            uout(iq,k) = ua(iq,k)
            vout(iq,k) = va(iq,k)
         end do
      end do
      return
      end        
