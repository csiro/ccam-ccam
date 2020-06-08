! Conformal Cubic Atmospheric Model
    
! Copyright 2015 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
! This file is part of the Conformal Cubic Atmospheric Model (CCAM)
!
! CCAM is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! CCAM is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with CCAM.  If not, see <http://www.gnu.org/licenses/>.

!------------------------------------------------------------------------------
    
module diag_m

#ifndef scm
   implicit none
   public :: printa, maxmin, average, diagvals
   private :: maxmin1, maxmin2, printa1, printa2
   interface maxmin
      module procedure maxmin1, maxmin2
   end interface
   interface printa
      module procedure printa1, printa2
   end interface
   interface diagvals
      module procedure diagvals_r, diagvals_i, diagvals_l
   end interface
contains
   subroutine printa2(name,a,ktau,level,i1,i2,j1,j2,bias,facti)
      ! printa has to work with arrays dimension both ifull and ifull+iextra
      ! printb has printj entry, and automatic choice of fact if facti=0.
      ! Have both 1D and multi-level versions.
      use cc_mpi
      use newmpar_m
      character(len=*), intent(in) :: name
      real, dimension(:,:), intent(in) :: a
      integer, intent(in) :: ktau, level, i1, i2, j1, j2
      real, intent(in) :: bias, facti

      call printa1(name,a(:,level),ktau,level,i1,i2,j1,j2,bias,facti)
   end subroutine printa2

   subroutine printa1(name,a,ktau,level,i1,i2,j1,j2,bias,facti)
      ! printa has to work with arrays dimension both ifull and ifull+iextra
      ! printb has printj entry, and automatic choice of fact if facti=0.
      ! Have both 1D and multi-level versions.
      use cc_mpi
      use newmpar_m
      character(len=*), intent(in) :: name
      real, dimension(:), intent(in) :: a
      integer, intent(in) :: ktau, level, i1, i2, j1, j2
      real, intent(in) :: bias, facti
      integer i, j, ja, jb, n, n2, ilocal, jlocal, nlocal
      real fact, atmp

      ! The indices i1, i2, j1, j2 are global
      n = (j1-1)/il_g ! Global n
      j = j1 - n*il_g ! Value on face
      if ( myid == fproc(i1,j,n) ) then
         ! Check if whole region is on this processor
         n2 = (j2-1)/il_g
         if ( fproc(i2, j2-n2*il_g, n2) /= myid ) then
           write(6,*)"Warning, printa region covers more than one processor"
           return    
         end if
         ja=j1
         jb=min(ja+24,j2)
         fact=facti
         ! Calculate factor from the middle of the face 
         ! This will be consistent for 1-6 processors at least
         nlocal = n + noff
         if(abs(facti)<1.e-20) then ! facti==0.
            atmp = abs(a(indp(ipan/2,jpan/2,nlocal)))
            if ( atmp > 0 ) then
               fact = 10./atmp
            else 
               fact = 1.0
            end if
         end if
         write(6,9) name,ktau,level,bias,fact
 9       format(/1x,a4,' ktau =',i7,'  level =',i3,'  addon =',g8.2,'  has been mult by',1pe8.1)
         write(6,91) (j,j=j1,j2)
91       format(4x,25i11)
         do i=i1,i2
            write(unit=*,fmt="(i5)", advance="no") i
            do j=ja,jb
               n = (j-1)/il_g ! Global n
               nlocal = n + noff
               ilocal = i - ioff
               jlocal = j - n*il_g - joff
               write(unit=*,fmt="(f11.6)", advance="no") (a(indp(ilocal,jlocal,nlocal))-bias)*fact
            end do
            write(*,*)
         end do
      end if

   end subroutine printa1

! has more general format & scaling factor  with subroutine average at bottom
   subroutine maxmin2(u,char,ktau,fact,kup)
      use cc_mpi
      use newmpar_m
      character(len=2), intent(in) :: char
      integer, intent(in) :: ktau, kup
      real, intent(in) :: fact
      real, dimension(:,:), intent(in) :: u
      real, dimension(2,kup) :: umin, umax
      integer, dimension(2,kup) :: ijumax,ijumin
      integer, dimension(1) :: imax,imin
      integer iqg
      integer i, j, k
      ! gumax(1,:) is true maximum, gumax(2,:) used for the location
      real, dimension(2,kup) :: gumax, gumin
      real, dimension(kup) :: gout
      real gmax, gmin
      
      gumax = 0.
      gumin = 0.

      do k=1,kup
         umax(1,k) = maxval(u(1:ifull,k))
         umin(1,k) = minval(u(1:ifull,k))
         ! Simpler to use real to hold the integer location. 
         ! No rounding problem for practical numbers of points
         ! Convert this to a global index
         !-----------------------------------------------------------
         imax=maxloc(u(1:ifull,k))
         imin=minloc(u(1:ifull,k))
         umax(2,k)=real(iq2iqg(imax(1)))
         umin(2,k)=real(iq2iqg(imin(1)))         
         !-----------------------------------------------------------  
      end do

      call ccmpi_reduce(umax,gumax,"maxloc",0,comm_world)
      call ccmpi_reduce(umin,gumin,"minloc",0,comm_world)

      if ( myid==0 ) then
          
         do k=1,kup
            iqg = nint(gumax(2,k))
            ! Convert to global i, j indices
            j = 1 + (iqg-1)/il_g
            i = iqg - (j-1)*il_g
            ijumax(:,k) = (/ i, j /)
            iqg = nint(gumin(2,k))
            j = 1 + (iqg-1)/il_g
            i = iqg - (j-1)*il_g
            ijumin(:,k) = (/ i, j /)
         end do
         
        if( abs(fact)<1.e-20 )then  !fact==0.
          gmax=-1.e20
          gmin=1.e20
          do k=1,kup
           gmax=max(gmax,gumax(1,k))
           gmin=min(gmin,gumin(1,k))
          enddo
          write(6,*) 'Gmax,Gmin ',gmax,gmin
          if(abs(gmax)<1.e-20)gmax=1. ! gmax==0.
          if(abs(gmin)<1.e-20)gmin=1. ! gmin==0.
          if ( kup>10 ) then
          write(6,981) ktau,char,gumax(1,1:10)/abs(gmax),     &
                          char,gumax(1,11:kup)/abs(gmax)
          end if
          write(6,977) ktau,ijumax
          write(6,982) ktau,char,gumin(1,:)/abs(gmin)
          write(6,977) ktau,ijumin
          return
        endif  ! (fact==0.)
        
        gumax(1,:)=gumax(1,:)*fact
        gumin(1,:)=gumin(1,:)*fact

       if(kup==1)then
        write(6,970) ktau,char,gumax(1,1),char,gumin(1,1)
970     format(i7,1x,a2,'max ',f8.3,3x,a2,'min ',f8.3)
        write(6,9705) ktau,ijumax(:,1),ijumin(:,1)
9705    format(i7,'  posij',i4,i4,10x,i3,i4)
        return
       endif   !  (kup.eq.1)

       if(gumax(1,1)>=1000.)then   ! for radon
        gout(:) = gumax(1,:)
        write(6,961) ktau,char,gout
!961     format(i7,1x,a2,'max ',10f7.1/(14x,10f7.1)/(14x,10f7.1))
961     format(i7,a3,'max ',10f7.1/(14x,10f7.1)/(14x,10f7.1))
        write(6,977) ktau,ijumax
        gout(:) = gumin(1,:)
        write(6,962) ktau,char,gout
!962     format(i7,1x,a2,'min ',10f7.1/(14x,10f7.1)/(14x,10f7.1))
962     format(i7,a3,'min ',10f7.1/(14x,10f7.1)/(14x,10f7.1))
        write(6,977) ktau,ijumin
       elseif(kup<=10)then  ! format for tggsn
        write(6,971) ktau,char,(gumax(1,k),k=1,kup)
        write(6,977) ktau,ijumax
        write(6,972) ktau,char,(gumin(1,k),k=1,kup)
        write(6,977) ktau,ijumin
       elseif(gumax(1,kup)>30.)then  ! format for T, & usually u,v
        gout(1:kup) = gumax(1,1:kup)
        write(6,971) ktau,char,gout(1:10),char,gout(11:kup)
!971     format(i7,1x,a2,'max ',10f7.2/(a10,'maX ',10f7.2)/(14x,10f7.2))
971     format(i7,a3,'max ',10f7.2/(a10,'maX ',10f7.2)/(14x,10f7.2))
        write(6,977) ktau,ijumax
        gout(:) = gumin(1,:)
        write(6,972) ktau,char,gout
!972     format(i7,1x,a2,'min ',10f7.2/(14x,10f7.2)/(14x,10f7.2))
972     format(i7,a3,'min ',10f7.2/(14x,10f7.2)/(14x,10f7.2))
        write(6,977) ktau,ijumin
977     format(i7,'  posij',10(i3,i4)/(14x,10(i3,i4))/(14x,10(i3,i4)))
       else  ! for qg & sd
        gout(1:kup) = gumax(1,1:kup)
        write(6,981) ktau,char,gout(1:10),char,gout(11:kup)
!981     format(i7,1x,a2,'max ',10f7.3/(a10,'maX ',10f7.3)/(14x,10f7.3))
981     format(i7,a3,'max ',10f7.3/(a10,'maX ',10f7.3)/(14x,10f7.3))        
        write(6,977) ktau,ijumax
        gout(1:kup) = gumin(1,1:kup)
        write(6,982) ktau,char,gout
!982     format(i7,1x,a2,'min ',10f7.3/(14x,10f7.3)/(14x,10f7.3))
982     format(i7,a3,'min ',10f7.3/(14x,10f7.3)/(14x,10f7.3))        
        write(6,977) ktau,ijumin
       endif
      endif ! myid == 0
      return
   end subroutine maxmin2

   subroutine maxmin1(u,char,ktau,fact,kup)
      use cc_mpi
      use newmpar_m
      character(len=2), intent(in) :: char
      integer, intent(in) :: ktau, kup
      real, intent(in) :: fact
      real, dimension(:), intent(in) :: u
      real, dimension(2) :: umin, umax
      integer, dimension(2) :: ijumax,ijumin
      integer, dimension(1) :: imax,imin
      integer i, j
      integer iqg
      real, dimension(2) :: gumax, gumin


      umax(1) = maxval(u(1:ifull))*fact
      umin(1) = minval(u(1:ifull))*fact
      !-----------------------------------------------------------
      ! MJT bug fix
      imax=maxloc(u(1:ifull))
      imin=minloc(u(1:ifull))
      umax(2)=real(iq2iqg(imax(1)))
      umin(2)=real(iq2iqg(imin(1)))
      !-----------------------------------------------------------
      call ccmpi_reduce(umax,gumax,"maxloc",0,comm_world)
      call ccmpi_reduce(umin,gumin,"minloc",0,comm_world)

      if ( myid == 0 ) then
        iqg = nint(gumax(2))
        ! Convert to global i, j indices
        j = 1 + (iqg-1)/il_g
        i = iqg - (j-1)*il_g
        ijumax(:) = (/ i, j /)
        iqg = nint(gumin(2))
        j = 1 + (iqg-1)/il_g
        i = iqg - (j-1)*il_g
        ijumin(:) = (/ i, j /)

        write(6,970) ktau,char,gumax(1),char,gumin(1)
970     format(i7,1x,a2,'max ',f10.3,3x,a2,'min ',f10.3)
        write(6,9705) ktau,ijumax,ijumin
9705    format(i7,'  posij',i5,i5,10x,i4,i4)
      end if ! myid == 0
      return
   end subroutine maxmin1

   subroutine average(speed,spmean_g,spavge_g)
      use cc_mpi
      use newmpar_m
      use sigs_m
      use sumdd_m
      use xyzinfo_m
      real, dimension(:,:), intent(in) :: speed
      real, dimension(:), intent(out) :: spmean_g
      real, intent(out) :: spavge_g
      real, dimension(ifull) :: tmpb
      complex, dimension(kl) :: tmpc,tmpc_g
      integer k
      
      tmpc=(0.,0.)
      do k=1,kl
         tmpb=speed(1:ifull,k)*wts(1:ifull)
         call drpdr_local(tmpb,tmpc(k))
      end do
      tmpc_g=(0.,0.)
      call ccmpi_reduce(tmpc,tmpc_g,"sumdr",0,comm_world)
      spmean_g=real(tmpc_g)
      spavge_g = 0.0
      do k=1,kl
         spavge_g = spavge_g-dsig(k)*spmean_g(k) ! dsig is -ve
      end do

   end subroutine average

   function diagvals_r(a) result (res)
      use cc_mpi
      use newmpar_m
      use parm_m
      real, intent(in), dimension(:) :: a
      real, dimension(9) :: res
      integer :: i, j, n, jf, ilocal, jlocal, nloc, iq

!     Return the equivalent of arr(ii+(jj-1)*il),ii=id-1,id+1),jj=jd-1,jd+1)
!     Note that this doesn't get off-processor values correctly!
!     Restrict range so that it still works if id=1 etc
      iq = 0
      res = 0. ! As a sort of missing value
      do j=min(jd-1,jl_g),max(jd+1,1)
         do i=max(id-1,1),min(id+1,il_g)
            iq = iq + 1
            n = (j-1)/il_g  ! Global n
            jf = j - n*il_g ! Value on face
            if ( fproc(i, jf, n) == myid ) then
               nloc = n + noff
               ilocal = i - ioff
               jlocal = j - n*il_g - joff
               res(iq) = a(min(max(indp(ilocal,jlocal,nloc),1),size(a)))
            end if
         end do
      end do
   end function diagvals_r

   function diagvals_i(a) result (res)
      use cc_mpi
      use newmpar_m
      use parm_m
      integer, intent(in), dimension(:) :: a
      integer, dimension(9) :: res
      integer :: i, j, n, jf, ilocal, jlocal, nloc, iq

!     Return the equivalent of arr(ii+(jj-1)*il),ii=id-1,id+1),jj=jd-1,jd+1)
      iq = 0
      res = 0
      do j=max(jd-1,1),min(jd+1,jl_g)
         do i=max(id-1,1),min(id+1,il_g)
            iq = iq + 1
            n = (j-1)/il_g  ! Global n
            jf = j - n*il_g ! Value on face
            if ( fproc(i, jf, n) == myid ) then
               nloc = n + noff
               ilocal = i - ioff
               jlocal = j - n*il_g - joff
               res(iq) = a(min(max(indp(ilocal,jlocal,nloc),1),size(a)))
            end if
         end do
      end do
   end function diagvals_i

   function diagvals_l(a) result (res)
      use cc_mpi
      use newmpar_m
      use parm_m
      logical, intent(in), dimension(:) :: a
      logical, dimension(9) :: res
      integer :: i, j, n, jf, ilocal, jlocal, nloc, iq

!     Return the equivalent of arr(ii+(jj-1)*il),ii=id-1,id+1),jj=jd-1,jd+1)
      iq = 0
      res = .false.
      do j=max(jd-1,1),min(jd+1,jl_g)
         do i=max(id-1,1),min(id+1,il_g)
            iq = iq + 1
            n = (j-1)/il_g  ! Global n
            jf = j - n*il_g ! Value on face
            if ( fproc(i, jf, n) == myid ) then
               nloc = n + noff
               ilocal = i - ioff
               jlocal = j - n*il_g - joff
               res(iq) = a(min(max(indp(ilocal,jlocal,nloc),1),size(a)))
            end if
         end do
      end do
   end function diagvals_l
#endif

#ifdef scm
   implicit none

   public :: printa, maxmin
   private :: maxmin1, maxmin2, printa1, printa2
   interface maxmin
      module procedure maxmin1, maxmin2
   end interface
   interface printa
      module procedure printa1, printa2
   end interface
contains
    
   subroutine maxmin2(u,char,ktau,fact,kup)
      character(len=2), intent(in) :: char
      integer, intent(in) :: ktau, kup
      real, intent(in) :: fact
      real, dimension(:,:), intent(in) :: u
   end subroutine maxmin2

   subroutine maxmin1(u,char,ktau,fact,kup)
      character(len=2), intent(in) :: char
      integer, intent(in) :: ktau, kup
      real, intent(in) :: fact
      real, dimension(:), intent(in) :: u
   end subroutine maxmin1

   subroutine printa2(name,a,ktau,level,i1,i2,j1,j2,bias,facti)
      character(len=*), intent(in) :: name
      real, dimension(:,:), intent(in) :: a
      integer, intent(in) :: ktau, level, i1, i2, j1, j2
      real, intent(in) :: bias, facti
   end subroutine printa2

   subroutine printa1(name,a,ktau,level,i1,i2,j1,j2,bias,facti)
      character(len=*), intent(in) :: name
      real, dimension(:), intent(in) :: a
      integer, intent(in) :: ktau, level, i1, i2, j1, j2
      real, intent(in) :: bias, facti
   end subroutine printa1

   
#endif
   
end module diag_m
