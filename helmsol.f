      subroutine helmsol(zz,zzn,zze,zzw,zzs,helm,s,rhs)

!     Solve Helmholtz equation using simple conjugate gradient method.
!     Each mode is solved separately. Highest numbered modes 
!     converge fastest.

      use cc_mpi
      implicit none
      include 'newmpar.h'
      include 'indices.h' ! in,is,iw,ie,inn,iss,iww,iee
      include 'parm.h'
      include 'parmdyn.h'
      include 'mpif.h'
      integer, parameter :: itmax=100 ! maximum number of iterations allowed
!     Arguments
      real, intent(in), dimension(ifull) :: zz,zzn,zze,zzw,zzs
      real helm(ifull+iextra,kl)      ! Helmholtz coefficients
      real s(ifull+iextra,kl)         ! Solution
      real rhs(ifull+iextra,kl)       ! RHS
!     real z(ifull,5)       ! Point coefficients, approx 1,1,1,1,-4.
      real, dimension(ifull,kl) :: fac, r, h
      real, dimension(ifull+iextra,kl) :: d

      real, dimension(kl) :: delta_0, delta_1, tau, alpha, beta, smag
      real, dimension(kl) :: gdelta_0, gdelta_1, gsmag, galpha
      real, dimension(2*kl) :: arr, garr
      integer iq, iter, k, klim, ierr

      call start_log(helm_begin)

      ! Assume bounds(s) has been done before call.
      do k=1,kl
         fac(1:ifull,k) = -1.0/(helm(1:ifull,k)-zz(1:ifull))
      end do

      delta_0 = 0.
      do k=1,kl
         do iq=1,ifull
            r(iq,k) = ( zze(iq)*s(ie(iq),k) + zzw(iq)*s(iw(iq),k) +
     &                  zzn(iq)*s(in(iq),k) + zzs(iq)*s(is(iq),k) -
     &                  rhs(iq,k) ) * fac(iq,k) + s(iq,k)

!            if ( k == 1 ) then
!               print*, "HELMSOL R", myid, iq, r(iq,k)
!            end if
            delta_0(k) = delta_0(k) + r(iq,k)*r(iq,k)
         end do
      end do
      call MPI_ALLREDUCE ( delta_0, gdelta_0, kl, MPI_REAL, MPI_SUM,
     &                     MPI_COMM_WORLD, ierr )
      if ( diag .and. mydiag ) then
         print*, "HELMSOL HELM", helm(idjd,:)
         print*, "HELMSOL S", s(idjd,:)
         print*, "HELMSOL RHS", rhs(idjd,:)
         print*, "HELMSOL R", r(idjd,:)
         print*, "GDELTA0", gdelta_0
      end if
      d(1:ifull,:) = -r(:,:)
      klim = kl ! All modes at first
      do iter = 1, itmax

         call bounds(d, klim=klim)
         alpha = 0.
         do k=1,klim
            do iq=1,ifull
               h(iq,k) = ( zze(iq)*d(ie(iq),k) + zzw(iq)*d(iw(iq),k) +
     &                     zzn(iq)*d(in(iq),k) + zzs(iq)*d(is(iq),k) ) *
     &                   fac(iq,k) + d(iq,k)
               alpha(k) = alpha(k) + d(iq,k)*h(iq,k)
            end do
         end do
         call MPI_ALLREDUCE ( alpha, galpha, klim, MPI_REAL, MPI_SUM,
     &                        MPI_COMM_WORLD, ierr )
         tau(1:klim) = gdelta_0(1:klim) / galpha(1:klim)
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
         ! Combine delta1 and smag into single array for reduction.
         arr(1:klim) = delta_1(1:klim)
         arr(klim+1:2*klim) = smag(1:klim)
         call MPI_ALLREDUCE ( arr, garr, 2*klim, MPI_REAL, MPI_SUM,
     &                        MPI_COMM_WORLD, ierr )
         gdelta_1(1:klim) = garr(1:klim)
         gsmag(1:klim) = garr(klim+1:2*klim)

!        Check which modes have converged
         do k=klim,1,-1
            if ( sqrt(gdelta_1(k)) > restol*sqrt(gsmag(k)) ) then
               ! This mode hasn't converged yet
               exit
            end if
         end do
!        Now k is the lowest mode yet to converge
         klim = k
         if ( klim == 0 ) exit
         beta(1:klim) = gdelta_1(1:klim) / gdelta_0(1:klim)
         gdelta_0(1:klim) = gdelta_1(1:klim)
         do k=1,klim
            d(1:ifull,k) = -r(:,k) + beta(k) * d(1:ifull,k)
         end do
         if ( (diag .or. ktau<6) .and. myid == 0 ) then
            print*, "Iterations", iter, klim, gdelta_1(1:klim)
         end if
      end do

      call end_log(helm_end)
      return
      end
