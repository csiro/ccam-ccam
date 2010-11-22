   subroutine helmsol(zz,zzn,zze,zzw,zzs,helm,s,rhs)

!     Solve Helmholtz equation using preconditioned conjugate gradient method.
!     Each mode is solved separately. Highest numbered modes 
!     converge fastest.
!     Note that the ILU preconditioner is expensive to set up. This code 
!     works best with ntbar=0 so the matrix doesn't change in time.

      use cc_mpi
      use ilu_m
      use indices_m
      use sumdd_m
      implicit none
      include 'newmpar.h'
      include 'parm.h'
      include 'parmdyn.h'
      include 'mpif.h'
      integer, parameter :: itmax=100 ! maximum number of iterations allowed
!     Arguments
      real, intent(in), dimension(ifull) :: zz,zzn,zze,zzw,zzs
!     WHY are helm and rhs ifull+iextra?????????
!     Not just for printa call ?????
      real, intent(in) :: helm(ifull+iextra,kl)      ! Helmholtz coefficients
      real, intent(inout) :: s(ifull+iextra,kl)         ! Solution
      real, intent(in) :: rhs(ifull+iextra,kl)       ! RHS
      real, dimension(ifull,kl) :: fac, invfac, v, sx
      real, dimension(ifull+iextra,kl) :: d, r, h

      real, dimension(kl) :: delta_0, delta_1, tau, alpha, beta, smag, &
                             gamma_1, sigma, delta
      real, dimension(kl) :: gdelta_0, gdelta_1, gsmag, galpha, ggamma_0, &
                             ggamma_1, gsigma, gdelta
      real, dimension(3*kl) :: arr, garr
      integer :: iq, iter, k, klim, ierr
      logical, save :: ilustart = .true.
      real, save :: factest
#ifdef sumdd
      complex, dimension(3*kl) :: local_sum, global_sum
!     Temporary array for the drpdr_local function
      real, dimension(ifull) :: tmparr, tmparr2 
#endif

      call start_log(helm_begin)

      klim = kl ! All modes at first

      do k=1,klim
         fac(1:ifull,k) = zz(1:ifull) - helm(1:ifull,k)
         invfac(1:ifull,k) = -1.0/(helm(1:ifull,k)-zz(1:ifull))
      end do
      if ( precon /= 0 ) then
         if ( ilustart) then
            call iludecomp(precon,fac(:,1:precon),zzn,zze,zzs,zzw)
            ilustart = .false.
!           Save this as a check on whether array has changed.
!           This will catch both changes in the time step and the reference
!           temperature profile
            factest = fac(1,1) 
         else
            if ( factest /= fac(1,1) ) then
               call iludecomp(precon,fac(:,1:precon),zzn,zze,zzs,zzw)
               factest = fac(1,1) 
            end if
         end if
      end if

      ! Use D'Azevedo method
      smag = 0.
#ifdef sumdd
      local_sum = (0.,0.)
#endif
      do k=1,klim
         if ( k <= precon) then
            do iq=1,ifull
               r(iq,k) = - ( zze(iq)*s(ie(iq),k) + zzw(iq)*s(iw(iq),k) + &
                             zzn(iq)*s(in(iq),k) + zzs(iq)*s(is(iq),k) - &
                             rhs(iq,k) ) - s(iq,k)*fac(iq,k)
#ifdef sumdd
               tmparr(iq) = s(iq,k)*s(iq,k)
#else
               smag(k) = smag(k) + s(iq,k)*s(iq,k)
#endif
            end do
            call ilusolve(h,r,k)
            d(1:ifull,k) = h(1:ifull,k)
         else
            do iq=1,ifull
               r(iq,k) = - ( zze(iq)*s(ie(iq),k) + zzw(iq)*s(iw(iq),k) + &
                 	     zzn(iq)*s(in(iq),k) + zzs(iq)*s(is(iq),k) - &
                             rhs(iq,k) ) * invfac(iq,k) - s(iq,k)
               d(iq,k) = r(iq,k)
#ifdef sumdd
               tmparr(iq) = s(iq,k)*s(iq,k)
#else
               smag(k) = smag(k) + s(iq,k)*s(iq,k)
#endif
            end do
         end if
#ifdef sumdd
         call drpdr_local(tmparr, local_sum(2*klim+k))
#endif
      end do

      call bounds(d, klim=klim)
      sigma = 0.
      gamma_1 = 0.
      do k=1,klim
         if ( k <= precon ) then
            do iq=1,ifull
               v(iq,k) = zze(iq)*d(ie(iq),k) + zzw(iq)*d(iw(iq),k) + &
                    zzn(iq)*d(in(iq),k) + zzs(iq)*d(is(iq),k) + &
                    d(iq,k)*fac(iq,k) 
#ifdef sumdd
               tmparr(iq) = d(iq,k)*v(iq,k)
               tmparr2(iq) = r(iq,k)*h(iq,k)
#else
               sigma(k) = sigma(k) + d(iq,k)*v(iq,k)
               gamma_1(k) = gamma_1(k) + r(iq,k)*h(iq,k)
#endif
            end do
         else
            do iq=1,ifull
               v(iq,k) = ( zze(iq)*d(ie(iq),k) + zzw(iq)*d(iw(iq),k) + &
                           zzn(iq)*d(in(iq),k) + zzs(iq)*d(is(iq),k) ) *  &
                           invfac(iq,k) + d(iq,k)
#ifdef sumdd
               tmparr(iq) = d(iq,k)*v(iq,k)
               tmparr2(iq) = r(iq,k)*r(iq,k)
#else
               sigma(k) = sigma(k) + d(iq,k)*v(iq,k)
               gamma_1(k) = gamma_1(k) + r(iq,k)*r(iq,k)
#endif
            end do
         end if
#ifdef sumdd
         call drpdr_local(tmparr, local_sum(k))
         call drpdr_local(tmparr2, local_sum(klim+k))
#endif
      end do
      call start_log(reduce_begin)
#ifdef sumdd
      call MPI_Allreduce (local_sum, global_sum, 2*klim, MPI_COMPLEX, &
                          MPI_SUMDR, MPI_COMM_WORLD, ierr) 
      gsigma(1:klim) = real(global_sum(1:klim))
      ggamma_1(1:klim) = real(global_sum(klim+1:2*klim))
      local_sum(1:2*klim) = (0.,0.) ! Still need the later part.
#else
      arr(1:klim) = sigma(1:klim)
      arr(klim+1:2*klim) = gamma_1(1:klim)
      call MPI_Allreduce ( arr, garr, 2*klim, MPI_REAL, MPI_SUM,       &
                           MPI_COMM_WORLD, ierr )
      gsigma(1:klim) = garr(1:klim)
      ggamma_1(1:klim) = garr(klim+1:2*klim)
#endif
      call end_log(reduce_end)
      alpha(1:klim) = ggamma_1(1:klim) / gsigma(1:klim)
      do k=1,klim
         s(1:ifull,k) = s(1:ifull,k) + alpha(k) * d(1:ifull,k)
         r(1:ifull,k) = r(1:ifull,k) - alpha(k) * v(1:ifull,k)
         if ( k <= precon) then
            call ilusolve(h,r,k)
         else
            ! Need this so that a single bounds call does both h and r
            h(1:ifull,k) = r(1:ifull,k)
         end if
      end do
      ggamma_0(1:klim) = ggamma_1(1:klim)
         
      do iter=2,itmax
         call bounds(h, klim=klim)
         gamma_1(1:klim) = 0.0
         delta(1:klim) = 0.0
         do k=1,klim
            if ( k <= precon ) then
               do iq=1,ifull
                  sx(iq,k) = zze(iq)*h(ie(iq),k) + zzw(iq)*h(iw(iq),k) + &
                             zzn(iq)*h(in(iq),k) + zzs(iq)*h(is(iq),k) + &
                             h(iq,k)*fac(iq,k)
#ifdef sumdd
                  tmparr(iq) = r(iq,k)*h(iq,k)
                  tmparr2(iq) = sx(iq,k)*h(iq,k)
#else
                  gamma_1(k) = gamma_1(k) + r(iq,k)*h(iq,k)
                  delta(k) = delta(k) + sx(iq,k)*h(iq,k)
#endif
               end do
            else
               do iq=1,ifull
                  ! Here h is a copy of r
                  sx(iq,k) = ( zze(iq)*h(ie(iq),k) + zzw(iq)*h(iw(iq),k) + &
                               zzn(iq)*h(in(iq),k) + zzs(iq)*h(is(iq),k) ) *&
                               invfac(iq,k) + h(iq,k)
#ifdef sumdd
                  tmparr(iq) = h(iq,k)*h(iq,k)
                  tmparr2(iq) = sx(iq,k)*h(iq,k)
#else
                  gamma_1(k) = gamma_1(k) + h(iq,k)*h(iq,k)
                  delta(k) = delta(k) + sx(iq,k)*h(iq,k)
#endif
               end do
            end if
#ifdef sumdd
            call drpdr_local(tmparr, local_sum(k))
            call drpdr_local(tmparr2, local_sum(klim+k))
#endif
         end do
            
         call start_log(reduce_begin)
#ifdef sumdd
         call MPI_Allreduce (local_sum, global_sum, 3*klim, MPI_COMPLEX, &
                             MPI_SUMDR, MPI_COMM_WORLD, ierr) 
         ggamma_1(1:klim) = real(global_sum(1:klim))
         gdelta(1:klim) = real(global_sum(klim+1:2*klim))
         gsmag(1:klim) = real(global_sum(2*klim+1:3*klim))
         local_sum = (0.,0.)
#else
         arr(1:klim) = gamma_1(1:klim)
         arr(klim+1:2*klim) = delta(1:klim)
         arr(2*klim+1:3*klim) = smag(1:klim)
         call MPI_Allreduce ( arr, garr, 3*klim, MPI_REAL, MPI_SUM,     &
                              MPI_COMM_WORLD, ierr )
         ggamma_1(1:klim) = garr(1:klim)
         gdelta(1:klim) = garr(klim+1:2*klim)
         gsmag(1:klim) = garr(2*klim+1:3*klim)
#endif
         call end_log(reduce_end)
         if ( (diag .or. ktau<6) .and. myid == 0 ) then
            write(6,'("Iterations",i4,i3,6g13.6/(10x,6g13.6))')   &
     &                 iter, klim, sqrt(abs(ggamma_1(1:klim)))
         end if
         !  Check which modes have converged
         do k=klim,1,-1
            if ( sqrt(abs(ggamma_1(k))) >= restol*sqrt(gsmag(k)) ) then
               exit
            end if
         end do
!        Now k is the lowest mode yet to converge
         klim = k
         if ( klim == 0 ) exit

         beta(1:klim) = ggamma_1(1:klim) / ggamma_0(1:klim)
         ggamma_0(1:klim) = ggamma_1(1:klim)
         gsigma(1:klim) = gdelta(1:klim) - beta(1:klim)**2*gsigma(1:klim)
         alpha(1:klim) = ggamma_1(1:klim) / gsigma(1:klim)
         smag(1:klim) = 0.
         do k=1,klim
            if ( k <= precon) then
               do iq=1,ifull
                  d(iq,k) = h(iq,k) + beta(k) * d(iq,k)
                  v(iq,k) = sx(iq,k) + beta(k) * v(iq,k)
                  s(iq,k) = s(iq,k) + alpha(k) * d(iq,k)
                  r(iq,k) = r(iq,k) - alpha(k) * v(iq,k)
#ifdef sumdd
                  tmparr(iq) = s(iq,k)*s(iq,k)
#else
                  smag(k) = smag(k) + s(iq,k)*s(iq,k)
#endif
               end do
               call ilusolve(h,r,k)
            else
               do iq=1,ifull
                  ! Use h in place of r
                  d(iq,k) = h(iq,k) + beta(k) * d(iq,k)
                  v(iq,k) = sx(iq,k) + beta(k) * v(iq,k)
                  s(iq,k) = s(iq,k) + alpha(k) * d(iq,k)
                  h(iq,k) = h(iq,k) - alpha(k) * v(iq,k)
#ifdef sumdd
                  tmparr(iq) = s(iq,k)*s(iq,k)
#else
                  smag(k) = smag(k) + s(iq,k)*s(iq,k)
#endif
               end do
            end if
#ifdef sumdd
            call drpdr_local(tmparr, local_sum(2*klim+k))
#endif
         end do

      end do

      call end_log(helm_end)

   end subroutine helmsol
