module ilu_m
   ! Routines for working with the sparse form of the incomplete LU
   ! decomposition.
   ! Values outside the processors local region are just ignored (no attempt
   ! to rebalance sum)
   implicit none
   private
   include 'newmpar.h'
   public :: iludecomp, ilusolve! , ilumul
   integer, dimension(ifull,4), private, save :: ileft
   integer, dimension(ifull), private, save   :: nleft, nright
   real, dimension(ifull,4,kl), private, save :: pleft
   real, dimension(ifull,kl), private, save   :: ppinv
contains
   subroutine iludecomp(ilumax,fac,zzn,zze,zzs,zzw)
      use cc_mpi
      use indices_m
      integer, intent(in) :: ilumax
      real, dimension(ifull), intent(in)  :: zzn, zze, zzs, zzw
      real, dimension(ifull,ilumax), intent(in)  :: fac
      real, dimension(ifull,4,ilumax) :: pp4
      real, dimension(ifull,ilumax) :: pp
      integer, dimension(ifull,4) :: ii4
      integer :: i, j, k, ir, itmp, jtmp, idir, jdir, ni, nj
      real, dimension(ilumax) :: e, dval
      integer, dimension(4) :: rlist, list2

      
      ! Need thse arrays so can index over the various cases
      ii4(:,1) = in
      ii4(:,2) = ie
      ii4(:,3) = is
      ii4(:,4) = iw
      pp = fac
      do k=1,ilumax
         pp4(:,1,k) = zzn
         pp4(:,2,k) = zze
         pp4(:,3,k) = zzs
         pp4(:,4,k) = zzw
      end do

      do ir=1,ifull
         dval = 1.0/pp(ir,:)
         ! i = r+1 .. n
         ! Look at all the neighbours. If i_x(r) > r, then this is a point
         ! below too by matrix symmetry.
         ! Use <= ifull check to get only points in this processor's region
         ni = count(ii4(ir,:)>ir .and. ii4(ir,:)<=ifull)
         rlist(1:ni) = pack(ii4(ir,:), ii4(ir,:)>ir .and. ii4(ir,:)<=ifull)
         do itmp=1,ni
            i = rlist(itmp)
            ! By the symmetry a(i,r) is non zero. To get value need to know
            ! for which direction x is i_x(i) = r. Direction index idir
            ! Use sum trick to get a scalar.
            idir = sum(pack( (/1, 2, 3, 4 /), ii4(i,:) == ir ))
            e = pp4(i,idir,:)*dval
            pp4(i,idir,:) = e
            
            ! j = r+1 .. n

            ! Point on diagonal certainly meets a[i,j] != 0, but is there
            ! another point in the same row?
            j = i
            if ( any ( ii4(ir,:) == j ) ) then
               jdir = sum(pack( (/1, 2, 3, 4 /), ii4(ir,:) == j ))
               pp(j,:) = pp(j,:) - e*pp4(ir,jdir,:)
            end if

            ! Loop checking other points for which j>r
            nj = count(ii4(i,:)>ir .and. ii4(ir,:)<=ifull)
            list2(1:nj) = pack(ii4(i,:), ii4(i,:)>ir .and. ii4(ir,:)<=ifull)
            do jtmp=1,nj ! points for which a(i,j) is non zero
               j = list2(jtmp)
               ! Check if a(r,j) is non-zero
               if ( any ( ii4(ir,:) == j ) ) then
                  jdir = sum(pack( (/1, 2, 3, 4 /), ii4(ir,:) == j ))
                  pp(j,:) = pp(j,:) - e*pp4(ir,jdir,:)
               end if
            end do
         end do
      end do

      nleft = 0
      do i=1,ifull
         do idir=1,4
            ! Pack indices and values that are to the left of the diagonal
            if ( ii4(i,idir) < i ) then
               nleft(i) = nleft(i) + 1
               ileft(i,nleft(i)) = ii4(i,idir)
               pleft(i,nleft(i),1:ilumax) = pp4(i,idir,:)
            end if
         end do
         nright(i) = nleft(i)
         do idir=1,4
            ! Pack indices and values that are to the right of the diagonal
            if ( ii4(i,idir) > i .and. ii4(i,idir)<= ifull ) then
               nright(i) = nright(i) + 1
               ileft(i,nright(i)) = ii4(i,idir)
               pleft(i,nright(i),1:ilumax) = pp4(i,idir,:)
            end if
         end do
      end do
      ppinv(:,1:ilumax) = 1./pp
      
   end subroutine iludecomp
   
   subroutine ilusolve(x,rhs,k)
      use cc_mpi
      use indices_m
      real, dimension(:,:), intent(in)  :: rhs
      real, dimension(:,:), intent(out) :: x
      integer, intent(in) :: k
      real, dimension(ifull) :: y
      integer :: i, itmp
      real :: tmpsum


      call start_log(precon_begin)
      ! Solve Cx=rhs, where C is in sparse LU form

!      ifull = size(rhs)
      ! First solve L y = rhs
      y(1) = rhs(1,k)
      do i=2,ifull
         ! Sum over points to left of diagonal
         tmpsum = 0.0
         do itmp=1,nleft(i)
            tmpsum = tmpsum+pleft(i,itmp,k)*y(ileft(i,itmp))
         end do
         y(i) = rhs(i,k) - tmpsum
      end do
      
      ! Now solve Ux = y by backsubstitution
      x(ifull,k) = y(ifull)*ppinv(ifull,k)
      do i=ifull-1,1,-1
         tmpsum = 0.0
         do itmp=nleft(i)+1,nright(i)
            tmpsum = tmpsum+pleft(i,itmp,k)*x(ileft(i,itmp),k)
         end do
         x(i,k) = (y(i) -tmpsum) * ppinv(i,k)
      end do

      call end_log(precon_end)
   end subroutine ilusolve

end module ilu_m
