      subroutine helmsol(accel,helm,s,rhs)

!     Solve Helmholtz equation.
!     For conformal-cubic this requires a 3 phase scheme
!     while for conformal-octagon it requires a 4 phase scheme
!     rather than simple red-black

      include 'newmpar.h'
      include 'indices.h' ! in,is,iw,ie,inn,iss,iww,iee
      include 'parm.h'
      include 'parmdyn.h'
!     itmax is maximum number of iterations allowed
      parameter(itmax=100,mev1=il+1-il/2*2,mev2=3-mev1)
c     mev1 = 1 for il even (2 for il odd)
c     mev2 = 2 for il even (1 for il odd)
!     Arguments
      real accel            ! SOR acceleration factor 
      real helm(ifull)      ! Helmholtz coefficients
      real s(ifull)         ! Solution
      real rhs(ifull)       ! RHS
      common/work2/zz(ifull),zzn(ifull),zze(ifull),zzw(ifull),
     . zzs(ifull),dum(il,jl,13)
!     real z(ifull,5)       ! Point coefficients, approx 1,1,1,1,-4.
      logical close_enough
      integer ip(2*il*il,3)! Cube only
      logical first
      save ip, first
      real, dimension(ifull) :: accfactor

      integer nface6(4,3)   ! Faces to use in each phase    (c-cub)
      integer ioff6(4,3)    ! Starting offset for each face (c-cub)
      integer nface14(8,4)  ! Faces to use in each phase    (c-oct)
      integer ioff14(8,4)   ! Starting offset for each face (c-oct)
      data nface6 / 0, 1, 3, 4,   0, 2, 3, 5,   1, 2, 4, 5 /
!     data ioff6 / 1, 1, 1, 1,   2, 1, 2, 1,   2, 2, 2, 2 /    ! up till 27/4/97
      data ioff6 / 1,mev1,mev2,2,  2,1,mev1,mev2,  mev2,2,1,mev1 / ! jlm general
!     data ioff6 / 1, 1, 2, 2,   2, 1, 1, 2,   2, 2, 1, 1 /        ! jlm even il
!     data ioff6 / 1, 2, 1, 2,   2, 1, 2, 1,   1, 2, 1, 2 /        ! jlm odd il

      data nface14/ 2, 4, 5,10,11,12,-1,-1,   2, 7, 8,10, 0, 1,-1,-1,
     .       1, 3, 5, 6, 7, 9,11,13,   4, 6, 8, 9,12,13, 0, 3/
      data ioff14 /1,mev1,1,2,mev2,2,-1,-1,   2,1,mev1,1,2,mev2,-1,-1, ! general
     .       mev1,1,2,mev2,2,1,mev1,mev2, mev2,mev1,mev2,2,1,mev1,1,2/ ! general
!     data ioff14 / 1, 1, 1, 2, 2, 2,-1,-1,   2, 1, 1, 1, 2, 2,-1,-1,  ! even il
!    .       1, 1, 2, 2, 2, 1, 1, 2,      2, 1, 2, 2, 1, 1, 1, 2/      ! even il
!     data ioff14 / 1, 2, 1, 2, 1, 2,-1,-1,   2, 1, 2, 1, 2, 1,-1,-1,  ! odd il
!    .       2, 1, 2, 1, 2, 1, 2, 1,      1, 2, 1, 2, 1, 2, 1, 2/      ! odd il
      data first /.true./

      ind(i,j,n)=i+(j-1)*il+n*il*il  ! *** for n=0,npanels

      if ( first ) then
         do iphase = 1,3         !  mrd code to run fast on NEC
            np = 0
            do iface=1,4
               if = nface6(iface,iphase)
               istart = ioff6(iface,iphase) ! 1 or 2
               do j=1,il
                  istart = 3 - istart ! 1 or 2 alternately
                  do i=istart,il,2
                     iq=ind(i,j,if)
                     np = np + 1
                     ip(np,iphase) = iq
                  end do
               end do
            end do
         end do
         first = .false.
      end if

      close_enough = .false.
      iter = 0

      accfactor(:) = accel/(helm(:)-zz(:))
      if(npanels.eq.5)then
        do while ( iter.lt.itmax .and. .not.close_enough )
           dsolmax=0.
           smax=0.
           close_enough=.true.
           ! Using acceleration of 1.0 on the first iteration improves 
           ! convergence when there's a good initial guess.
           if ( iter == 0 ) then
              accfactor(:) = 1.0/(helm(:)-zz(:))
           else
              accfactor(:) = accel/(helm(:)-zz(:))
           end if
           do iphase = 1,3
*cdir nodep
              do i=1,2*il*il
                 iq=ip(i,iphase)
                 dsol= ( zzn(iq)*s(in(iq)) + zzw(iq)*s(iw(iq)) +
     &                         zze(iq)*s(ie(iq)) + zzs(iq)*s(is(iq)) +
     &                      ( zz(iq)-helm(iq) )*s(iq) - rhs(iq) )
     &                      * accfactor(iq)
                 s(iq) = s(iq) + dsol
                 dsolmax=max(dsolmax,abs(dsol))
                 smax=max(smax,abs(s(iq)))
              end do
           end do
           if (dsolmax.gt.restol*smax)close_enough=.false.
           iter = iter + 1
        end do
      elseif(npanels.eq.13)then
        do while ( iter.lt.itmax .and. .not.close_enough )
           dsolmax=0.
           smax=0.
           close_enough=.true.
           do iphase = 1,4
              do iface=1,6+(iphase-1)/2*2   ! to 6 for 1 & 2; to 8 for 3 & 4
                 if = nface14(iface,iphase)
                 istart = ioff14(iface,iphase)  ! 1 or 2
                 do j=1,il
                    istart = 3 - istart ! 1 or 2 alternately
                    do i=istart,il,2
                       iq=ind(i,j,if)
                       dsol= ( zzn(iq)*s(in(iq)) + zzw(iq)*s(iw(iq)) +
     &                         zze(iq)*s(ie(iq)) + zzs(iq)*s(is(iq)) +
     &                      ( zz(iq)-helm(iq) )*s(iq) - rhs(iq) )
     &                      * accfactor(iq)
                       s(iq) = s(iq) + dsol
                       dsolmax=max(dsolmax,abs(dsol))
                       smax=max(smax,abs(s(iq)))
c                      if ( abs(dsol) .gt. abs(restol*s(iq)) )
c    &                      close_enough=.false.
                    end do
                 end do
              end do
           end do
           if (dsolmax.gt.restol*smax)close_enough=.false.
           iter = iter + 1
        end do
      endif    !  (npanels.eq.5)elseif(npanels.eq.13)

      if (diag.or.ktau.lt.6)print*,'helmsol acc, iterations ',
     .                        accel, iter
      return
      end
