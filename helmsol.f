      subroutine helmsol(helm,s,rhs)

!     Solve Helmholtz equation using simple conjugate gradient method.
!     Each mode is solved separately. Highest numbered modes 
!     converge fastest.

      implicit none
      include 'newmpar.h'
      include 'indices.h' ! in,is,iw,ie,inn,iss,iww,iee
      include 'parm.h'
      include 'parmdyn.h'
      integer, parameter :: itmax=100 ! maximum number of iterations allowed
!     Arguments
      real helm(ifull,kl)      ! Helmholtz coefficients
      real s(ifull,kl)         ! Solution
      real rhs(ifull,kl)       ! RHS
      real zz(ifull),zzn(ifull),zze(ifull),zzw(ifull),
     . zzs(ifull),dum(il,jl,13)
      common/work2/zz,zzn,zze,zzw,zzs,dum
!     real z(ifull,5)       ! Point coefficients, approx 1,1,1,1,-4.
      real, dimension(ifull,kl) :: fac, r, d, h

      real, dimension(kl) :: delta_0, delta_1, tau, alpha, beta, smag
      integer iq, iter, k, klim
      real, dimension(kl) :: dsolmax, smax
      real :: dsol

      do k=1,kl
         fac(:,k) = -1.0/(helm(:,k)-zz(:))
      end do

      delta_0 = 0.
      do k=1,kl
         do iq=1,ifull
            r(iq,k) = ( zze(iq)*s(ie(iq),k) + zzw(iq)*s(iw(iq),k) +
     &                  zzn(iq)*s(in(iq),k) + zzs(iq)*s(is(iq),k) -
     &                  rhs(iq,k) ) * fac(iq,k) + s(iq,k)

            delta_0(k) = delta_0(k) + r(iq,k)*r(iq,k)
         end do
      end do

      d(:,:) = -r(:,:)
      klim = kl ! All modes at first
      do iter = 1, itmax

         alpha = 0.
         do k=1,klim
            do iq=1,ifull
               h(iq,k) = ( zze(iq)*d(ie(iq),k) + zzw(iq)*d(iw(iq),k) +
     &                     zzn(iq)*d(in(iq),k) + zzs(iq)*d(is(iq),k) ) *
     &                   fac(iq,k) + d(iq,k)
               alpha(k) = alpha(k) + d(iq,k)*h(iq,k)
            end do
         end do

         tau(1:klim) = delta_0(1:klim) / alpha(1:klim)
         delta_1 = 0.
         smag = 0.
         do k=1,klim
            do iq=1,ifull
               s(iq,k) = s(iq,k) + tau(k) * d(iq,k)
               r(iq,k) = r(iq,k) + tau(k) * h(iq,k)
               delta_1(k) = delta_1(k) + r(iq,k)*r(iq,k) ! Magnitude of residual
               smag(k) = smag(k) + s(iq,k)*s(iq,k)
            end do
         end do
!        Check which modes have converged
         do k=klim,1,-1
            if ( sqrt(delta_1(k)) > restol*sqrt(smag(k)) ) then
               ! This mode hasn't converged yet
               exit
            end if
         end do
!        Now k is the lowest mode yet to converge
         klim = k
         if ( klim == 0 ) exit
         beta(1:klim) = delta_1(1:klim) / delta_0(1:klim)
         delta_0(1:klim) = delta_1(1:klim)
         do k=1,klim
            d(:,k) = -r(:,k) + beta(k) * d(:,k)
         end do
         if (diag .or. ktau<6) then
            print*, "Iterations", iter, klim, delta_1(1:klim)
         end if
      end do

      if (diag .or. ktau<6) then
         dsolmax = 0.
         smax = 0.
         do k=1,kl
            do iq=1,ifull
               dsol = ( zzn(iq)*s(in(iq),k) + zzw(iq)*s(iw(iq),k) +
     &                  zze(iq)*s(ie(iq),k) + zzs(iq)*s(is(iq),k) +
     &                      ( zz(iq)-helm(iq,k) )*s(iq,k) - rhs(iq,k) )
               dsolmax(k) = max(dsolmax(k),abs(dsol))
               smax(k) = max(smax(k),abs(s(iq,k)))
            end do
         end do
         print*,'helmsol iterations ', iter, "final error", dsolmax/smax
      end if

      return
      end
