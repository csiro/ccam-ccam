      subroutine helmsol(helm,s,rhs)

!     Solve Helmholtz equation using simple conjugate gradient method.

      implicit none
      include 'newmpar.h'
      include 'indices.h' ! in,is,iw,ie,inn,iss,iww,iee
      include 'parm.h'
      include 'parmdyn.h'
      integer, parameter :: itmax=100 ! maximum number of iterations allowed
!     Arguments
      real helm(ifull)      ! Helmholtz coefficients
      real s(ifull)         ! Solution
      real rhs(ifull)       ! RHS
      real zz(ifull),zzn(ifull),zze(ifull),zzw(ifull),
     . zzs(ifull),dum(il,jl,13)
      common/work2/zz,zzn,zze,zzw,zzs,dum
!     real z(ifull,5)       ! Point coefficients, approx 1,1,1,1,-4.
      real, dimension(ifull) :: fac, r, d, h

      real delta_0, delta_1, tau, alpha, beta, smag
      integer iq, iter
      real :: dsol, dsolmax, smax

      fac(:) = -1.0/(helm(:)-zz(:))

      delta_0 = 0.
      do iq=1,ifull
         r(iq) = ( zze(iq)*s(ie(iq)) + zzw(iq)*s(iw(iq)) +
     &             zzn(iq)*s(in(iq)) + zzs(iq)*s(is(iq)) -
     &             rhs(iq) ) * fac(iq) + s(iq)

         delta_0 = delta_0 + r(iq)*r(iq)
      end do

      d(:) = -r(:)
      do iter = 1, itmax

         alpha = 0.
         do iq=1,ifull
            h(iq) = ( zze(iq)*d(ie(iq)) + zzw(iq)*d(iw(iq)) +
     &                zzn(iq)*d(in(iq)) + zzs(iq)*d(is(iq)) ) * fac(iq)
     &              + d(iq)
            alpha = alpha + d(iq)*h(iq)
         end do

         tau = delta_0 / alpha
         delta_1 = 0.
         smag = 0.
         do iq=1,ifull
            s(iq) = s(iq) + tau * d(iq)
            r(iq) = r(iq) + tau * h(iq)
            delta_1 = delta_1 + r(iq)*r(iq) ! Magnitude of residual
            smag = smag + s(iq)*s(iq)
         end do
         if ( sqrt(delta_1) < restol*sqrt(smag) ) then
            exit
         end if
         beta = delta_1 / delta_0
         delta_0 = delta_1
         d(:) = -r(:) + beta * d(:)
      end do

      if (diag .or. ktau<6) then
         dsolmax = 0.
         smax = 0.
         do iq=1,ifull
            dsol= ( zzn(iq)*s(in(iq)) + zzw(iq)*s(iw(iq)) +
     &                         zze(iq)*s(ie(iq)) + zzs(iq)*s(is(iq)) +
     &                      ( zz(iq)-helm(iq) )*s(iq) - rhs(iq) )
            dsolmax=max(dsolmax,abs(dsol))
            smax=max(smax,abs(s(iq)))
         end do
         print*,'helmsol iterations ', iter, "final error", dsolmax/smax
      end if

      return
      end
