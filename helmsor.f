      subroutine helmsor(zz,zzn,zze,zzw,zzs,helm,s,rhs)
!     Solve Helmholtz equation - experimental jlm version
!     e.g. use  -2540, -2335, -2635, -2425 (in decr. order)

      use cc_mpi
      use diag_m
      use ilu_m
      use sumdd_m
      implicit none
      include 'newmpar.h'
      include 'indices.h' ! in,is,iw,ie,inn,iss,iww,iee
      include 'parm.h'
      include 'parmdyn.h'
      include 'vecs.h'
      include 'mpif.h'
!     integer, parameter :: meth=3      ! 3b, 5c1 good
      integer, parameter :: ntest=0 
      integer, parameter :: itmax=300 ! maximum number of iterations allowed
!     Arguments
      real, intent(in), dimension(ifull) :: zz,zzn,zze,zzw,zzs
!     WHY are helm and rhs ifull+iextra?????????
!     Not just for printa call ?????
      real, intent(in) :: helm(ifull+iextra,kl)      ! Helmholtz coefficients
      real, intent(inout) :: s(ifull+iextra,kl)      ! Solution
      real, intent(in) :: rhs(ifull+iextra,kl)       ! RHS
      real, dimension(ifull,kl) :: sb, sa, snew, dsol
      integer, save, dimension(ifull) :: mask
      integer, dimension(kl) :: iters
      real, dimension(kl) ::  dsolmax, dsolmax_g, smax, smax_g
      real, dimension(kl) ::  smin, smin_g
      real aa, bb, cc, axel, accel(kl)
      integer iq, iter, k, nx, j, jx, i, klim, ierr, meth, nx_max
      logical first
      save  first, meth, nx_max, axel, accel
      data first /.true./

      if(first)then
        if(precon==-1)precon=-2325  ! i.e. 2, 3, .25
        nx_max=abs(precon)/1000
        meth=abs(precon)/100-10*nx_max
        axel=-.01*real(precon)-10*nx_max-meth
        if(myid==0)print *,'in helmsor nx_max,meth,axel: ',
     &                                 nx_max,meth,axel
        mask(:)=1
        do j=1,jl
         do i=1,il,2
          jx=mod(i+j,2)    ! 0 or 1  starting (1,1) on panel
          iq=i+jx+(j-1)*il
          if(nx_max==2)then
            mask(iq)=2
          else
            if(j>il.and.j<=3*il)mask(iq)=3
            if(j>3*il.and.j<=5*il)mask(iq)=2
            jx=mod(i+j+1,2)  ! 0 or 1  starting (2,1) on panel
            iq=i+jx+(j-1)*il
            if(j<=2*il)mask(iq)=2
            if(j>4*il)mask(iq)=3
          endif 
         enddo
c       print *,'j ',j,(mask(iq),iq=1+(j-1)*il,6+(j-1)*il)
       enddo
       first=.false.
      endif  ! (first)
      if(ktau==1)then
       do k=1,kl
        call optmx(il,schmidt,dt,bam(k),accel(k))
        if(myid==0)print *,'k,accel ',k,accel(k)
       enddo
      endif
 
      if(precon<-2899)go to 5  ! e.g. -2900 or -3900
      klim=kl
      iter = 1
      do while ( iter<itmax .and. klim>1)
       call bounds(s, klim=klim,nrows=2)
       do k=1,klim        
        do nx=1,nx_max
         do iq=1,ifull
          if(mask(iq)==nx)then
             dsol(iq,k)=( zzn(iq)*s(in(iq),k) + zzw(iq)*s(iw(iq),k)    
     &          +zze(iq)*s(ie(iq),k) + zzs(iq)*s(is(iq),k)
     &          +(zz(iq)-helm(iq,k))*s(iq,k) - rhs(iq,k))      
     &               /(helm(iq,k)-zz(iq))
             snew(iq,k) = s(iq,k) + dsol(iq,k)
c            snew(iq,k) = s(iq,k) + dsol(iq,k)*accel(k)  ! no help
          endif
         enddo

!        following are jlm methods for improving guess
         if(iter>=3.and.meth==3)then   ! qls
          do iq=1,ifull
           if(mask(iq)==nx)then
            aa=(sb(iq,k)-3*sa(iq,k)+3*s(iq,k)+19*snew(iq,k))/20.  
            bb=(9*sb(iq,k)-17*sa(iq,k)-13*s(iq,k)+21*snew(iq,k))/20.  
c           snew(iq,k)=aa+.25*bb+.0625*cc !3c
            snew(iq,k)=aa+axel*bb          
           endif
          enddo
         endif  ! meth=3
         if(iter>=3.and.meth==4)then   ! oscill
          do iq=1,ifull
           if(mask(iq)==nx)then
            aa=(7.*snew(iq,k)+3.*s(iq,k)-3.*sa(iq,k)+sb(iq,k))/8. ! oscill
            bb=snew(iq,k)-.5*s(iq,k)-sa(iq,k)+.5*sb(iq,k)         ! oscill
c           cc=.25*(snew(iq,k)-s(iq,k)-sa(iq,k)+sb(iq,k))         !             aa=(sb(iq,k)-3*sa(iq,k)+3*s(iq,k)+19*snew(iq,k))/20.  
            snew(iq,k)=aa+axel*bb          
           endif
          enddo
         endif  ! meth=4
         if(iter>=3.and.meth==5)then   ! wqls
          do iq=1,ifull
           if(mask(iq)==nx)then
            aa=(2*sb(iq,k)-6*sa(iq,k)+6*s(iq,k)+68*snew(iq,k))/70.  
            bb=(5*sb(iq,k)-8*sa(iq,k)-13*s(iq,k)+16*snew(iq,k))/14.  
c           cc=(3*sb(iq,k)-2*sa(iq,k)-5*s(iq,k)+4*snew(iq,k))/14.  
            snew(iq,k)=aa+axel*bb          
          endif
          enddo
         endif  ! meth=5
         if(iter>=3.and.meth==6)then   ! wqls again
          do iq=1,ifull
           if(mask(iq)==nx)then
            aa=(2*sb(iq,k)-6*sa(iq,k)+6*s(iq,k)+68*snew(iq,k))/70.  
            bb=(5*sb(iq,k)-8*sa(iq,k)-13*s(iq,k)+16*snew(iq,k))/14.  
            cc=(3*sb(iq,k)-2*sa(iq,k)-5*s(iq,k)+4*snew(iq,k))/14.  
            snew(iq,k)=aa+axel*(bb+axel*cc)          
          endif
          enddo
         endif  ! meth=5

         do iq=1,ifull
          if(mask(iq)==nx)then
c           rotate s files
            sb(iq,k)=sa(iq,k)
            sa(iq,k)=s(iq,k)
            s(iq,k)=snew(iq,k)
          endif
         enddo
        enddo  ! nx loop
        iters(k)=iter
      enddo ! k loop   
      if(ntest>0.and.meth>=1.and.mydiag)
     &  write (6,"('Iter ,s',i4,4f14.5)") iter,(s(iq,1),iq=1,4)
      do k=1,klim
       if(iter==1)then
         smax(k) = maxval(s(1:ifull,k))
         smin(k) = minval(s(1:ifull,k))
       endif
       dsolmax(k) = maxval(abs(dsol(1:ifull,k)))
      enddo
      if(iter==1)then
        call MPI_Reduce( smax, smax_g, klim, MPI_REAL, MPI_MAX, 0,
     &                    MPI_COMM_WORLD, ierr )
        call MPI_Reduce( smin, smin_g, klim, MPI_REAL, MPI_MIN, 0,
     &                    MPI_COMM_WORLD, ierr )
      endif
      call MPI_Reduce( dsolmax, dsolmax_g, klim, MPI_REAL, MPI_MAX, 0,
     &                    MPI_COMM_WORLD, ierr )
      if(myid==0.and.ntest>0)then
        print *,'smin_g ',smin_g(:)
        print *,'smax_g ',smax_g(:)
        print *,'dsolmax_g ',dsolmax_g(:)
      endif  ! (myid==0)
      do k=klim,1,-1
       if(dsolmax_g(k)<restol*(smax_g(k)-smin_g(k)))then
c        print *,'k,klim,iter,restol ',k,klim,iter,restol
         klim=k
       endif
      enddo
      call MPI_Bcast(klim,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr) 
      iter = iter + 1
      enddo   ! while( iter<itmax .and. klim>1)

      if(myid==0.and.(diag.or.ktau<6))then
        do k=1,kl
         print*,'helmjlm ktau,k,Iterations ',ktau,k,iters(k)
        enddo
      endif
      return

5     continue  
      klim=kl
      iter = 1
      do while ( iter<itmax .and. klim>1)
       call bounds(s, klim=klim,nrows=2)
       do k=1,klim        
        do nx=1,nx_max
         do iq=1,ifull
          if(mask(iq)==nx)then
             dsol(iq,k)=( zzn(iq)*s(in(iq),k) + zzw(iq)*s(iw(iq),k)    
     &          +zze(iq)*s(ie(iq),k) + zzs(iq)*s(is(iq),k)
     &          +(zz(iq)-helm(iq,k))*s(iq,k) - rhs(iq,k))      
     &               /(helm(iq,k)-zz(iq))
             snew(iq,k) = s(iq,k) + accel(k)*dsol(iq,k)
          endif
         enddo

         do iq=1,ifull
          if(mask(iq)==nx)then
            s(iq,k)=snew(iq,k)
          endif
         enddo
        enddo  ! nx loop
        iters(k)=iter
      enddo ! k loop   
      if(ntest>0.and.mydiag)
     &  write (6,"('Iter ,s',i4,4f14.5)") iter,(s(iq,1),iq=1,4)
      do k=1,klim
       if(iter==1)then
         smax(k) = maxval(s(1:ifull,k))
         smin(k) = minval(s(1:ifull,k))
         if(ntest>0.and.mydiag)print *,'k,smax,smin ',k,smax(k),smin(k)
       endif
c      write (6,"('iter,k ,s',2i4,4f14.5)") iter,k,(s(iq,k),iq=1,4)
       dsolmax(k) = maxval(abs(dsol(1:ifull,k)))
c       smin_g(k)=1.e20
c       smax_g(k)=-1.e20
c       do iq=1,ifull
c        smin_g(k)=min(smin_g(k),abs(dsol(iq,k)))
c        smax_g(k)=min(smax_g(k),abs(dsol(iq,k)))
c       enddo
      enddo
      if(iter==1)then
        call MPI_Reduce( smax, smax_g, klim, MPI_REAL, MPI_MAX, 0,
     &                    MPI_COMM_WORLD, ierr )
        call MPI_Reduce( smin, smin_g, klim, MPI_REAL, MPI_MIN, 0,
     &                    MPI_COMM_WORLD, ierr )
      endif
      call MPI_Reduce( dsolmax, dsolmax_g, klim, MPI_REAL, MPI_MAX, 0,
     &                    MPI_COMM_WORLD, ierr )
      if(myid==0.and.ntest>0)then
        print *,'smin_g ',smin_g(:)
        print *,'smax_g ',smax_g(:)
        print *,'dsolmax_g ',dsolmax_g(:)
      endif  ! (myid==0)
      do k=klim,1,-1
       if(dsolmax_g(k)<restol*(smax_g(k)-smin_g(k)))then
c        print *,'k,klim,iter,restol ',k,klim,iter,restol
         klim=k
       endif
      enddo
      call MPI_Bcast(klim,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr) 
      iter = iter + 1
      enddo   ! while( iter<itmax .and. klim>1)

      if(myid==0.and.(diag.or.ktau<6))then
        do k=1,kl
         print*,'helmjlm ktau,k,Iterations ',ktau,k,iters(k)
        enddo
      endif
      return

      end
