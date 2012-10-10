      subroutine helmsor(zz,zzn,zze,zzw,zzs,helm,s,rhs)
!     Solve Helmholtz equation - experimental jlm version
!     e.g. use  -2540, -2335, -2635, -2425 (in decr. order)

      use cc_mpi
      use diag_m
      use ilu_m
      use indices_m
      use sumdd_m
      use vecs_m
      implicit none
      include 'newmpar.h'
      include 'parm.h'
      include 'parmdyn.h'
      include 'parmgeom.h'  ! rlong0,rlat0,schmidt  - briefly
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
      integer, allocatable, save, dimension(:) :: mask
      integer, dimension(kl) :: iters
      integer, dimension(3),save :: ifullx    ! MJT pack
      integer, dimension(:,:), allocatable, save :: iqx,iqn,iqe ! MJT pack
      integer, dimension(:,:), allocatable, save :: iqw,iqs     ! MJT pack
      real, dimension(kl) ::  dsolmax, dsolmax_g, smax, smax_g
      real, dimension(kl) ::  smin, smin_g
      real, dimension(:), allocatable, save :: accel ! MJT pack
      real aa(ifull), bb(ifull), cc(ifull), axel
      real ctst1,ctst2,mm,ww,yy
      integer iq, iter, k, nx, j, jx, i,klim,klimnew, ierr, meth, nx_max
      integer ifx,itstest,its1,its2
      save  meth, nx_max, axel

      call start_log(helm_begin) ! MJT
      
      itstest=1 ! just a default
      ctst1=9.E9
      ctst2=8.E9
      its1=-1
      its2=0
      
      if (.not.allocated(mask)) then
        allocate(mask(ifull))
        allocate(iqx(ifull,3),iqn(ifull,3),iqe(ifull,3))
        allocate(iqw(ifull,3),iqs(ifull,3))
        allocate(accel(kl))

        if(precon==-1)precon=-2325  ! i.e. 2, 3, .25
        nx_max=abs(precon)/1000
        meth=abs(precon)/100-10*nx_max
        axel=-.01*real(precon)-10*nx_max-meth
        if(myid==0)print *,'in helmsor nx_max,meth,axel: ',
     &                                 nx_max,meth,axel
        !print *,'kind precon, axel ',kind(precon),kind(axel)
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
       ifullx=0                                ! MJT pack
       iqx=0                                   ! MJT pack
       iqn=0                                   ! MJT pack
       iqe=0                                   ! MJT pack
       iqw=0                                   ! MJT pack
       iqs=0                                   ! MJT pack
       do iq=1,ifull                           ! MJT pack
         ifullx(mask(iq))=ifullx(mask(iq))+1   ! MJT pack
         iqx(ifullx(mask(iq)),mask(iq))=iq     ! MJT pack
         iqn(ifullx(mask(iq)),mask(iq))=in(iq) ! MJT pack
         iqe(ifullx(mask(iq)),mask(iq))=ie(iq) ! MJT pack
         iqw(ifullx(mask(iq)),mask(iq))=iw(iq) ! MJT pack
         iqs(ifullx(mask(iq)),mask(iq))=is(iq) ! MJT pack
       end do                                  ! MJT pack

       do k=1,kl
        call optmx(il_g,schmidt,dt,bam(k),accel(k))
	if(il_g>il)accel(k)=1.+.55*(accel(k)-1.)  ! for uniform-dec 22/4/08
c       if(il_g==il)accel(k)=1.+.55*(accel(k)-1.) ! just a test
        if(myid==0)print *,'k,accel ',k,accel(k)
       enddo
      endif
c$      print *,'myid,ktau,precon ',myid,ktau,precon
 
      if(precon<-2899)go to 5  ! e.g. -2900 or -3900
      klim=kl
      iter = 1
      do while ( iter<itmax .and. klim>1)
       call bounds(s, klim=klim)
       do k=1,klim        
        do nx=1,nx_max
!         do iq=1,ifull
!          if(mask(iq)==nx)then
           ifx=ifullx(nx)
           dsol(iqx(1:ifx,nx),k)=
     &        ( zzn(iqx(1:ifx,nx))*s(iqn(1:ifx,nx),k)
     &        + zzw(iqx(1:ifx,nx))*s(iqw(1:ifx,nx),k)
     &        + zze(iqx(1:ifx,nx))*s(iqe(1:ifx,nx),k)
     &        + zzs(iqx(1:ifx,nx))*s(iqs(1:ifx,nx),k)
     &        +(zz(iqx(1:ifx,nx))
     &        -helm(iqx(1:ifx,nx),k))*s(iqx(1:ifx,nx),k)
     &        - rhs(iqx(1:ifx,nx),k))      
     &        /(helm(iqx(1:ifx,nx),k)-zz(iqx(1:ifx,nx)))
           snew(iqx(1:ifx,nx),k) = s(iqx(1:ifx,nx),k)
     &        + dsol(iqx(1:ifx,nx),k)
c            snew(iq,k) = s(iq,k) + dsol(iq,k)*accel(k)  ! no help
!          endif
!         enddo

!        following are jlm methods for improving guess
         if(iter>=3)then
         select case(meth)
          case(3)
!         if(iter>=3.and.meth==3)then   ! qls
!          do iq=1,ifull
!           if(mask(iq)==nx)then
            aa(iqx(1:ifx,nx))=
     &       (sb(iqx(1:ifx,nx),k)
     &       -3.*sa(iqx(1:ifx,nx),k)
     &       +3.*s(iqx(1:ifx,nx),k)
     &       +19.*snew(iqx(1:ifx,nx),k))/20.  
            bb(iqx(1:ifx,nx))=(9.*sb(iqx(1:ifx,nx),k)
     &       -17.*sa(iqx(1:ifx,nx),k)
     &       -13.*s(iqx(1:ifx,nx),k)
     &       +21.*snew(iqx(1:ifx,nx),k))/20.  
c           snew(iq,k)=aa+.25*bb+.0625*cc !3c
            snew(iqx(1:ifx,nx),k)=aa(iqx(1:ifx,nx))
     &       +axel*bb(iqx(1:ifx,nx))
!           endif
!          enddo
!         endif  ! meth=3
          case(4)
!         if(iter>=3.and.meth==4)then   ! oscill
!          do iq=1,ifull
!           if(mask(iq)==nx)then
            aa(iqx(1:ifx,nx))=(7.*snew(iqx(1:ifx,nx),k)
     &       +3.*s(iqx(1:ifx,nx),k)
     &       -3.*sa(iqx(1:ifx,nx),k)
     &       +sb(iqx(1:ifx,nx),k))/8. ! oscill
            bb(iqx(1:ifx,nx))=snew(iqx(1:ifx,nx),k)
     &       -.5*s(iqx(1:ifx,nx),k)
     &       -sa(iqx(1:ifx,nx),k)
     &       +.5*sb(iqx(1:ifx,nx),k)         ! oscill
c           cc=.25*(snew(iq,k)-s(iq,k)-sa(iq,k)+sb(iq,k))         !             aa=(sb(iq,k)-3*sa(iq,k)+3*s(iq,k)+19*snew(iq,k))/20.  
            snew(iqx(1:ifx,nx),k)=aa(iqx(1:ifx,nx))
     &       +axel*bb(iqx(1:ifx,nx))
!           endif
!          enddo
!         endif  ! meth=4
          case(5)
!         if(iter>=3.and.meth==5)then   ! wqls
!          do iq=1,ifull
!           if(mask(iq)==nx)then
            aa(iqx(1:ifx,nx))=(2.*sb(iqx(1:ifx,nx),k)
     &       -6.*sa(iqx(1:ifx,nx),k)+6*s(iqx(1:ifx,nx),k)
     &       +68.*snew(iqx(1:ifx,nx),k))/70.  
            bb(iqx(1:ifx,nx))=(5.*sb(iqx(1:ifx,nx),k)
     &       -8.*sa(iqx(1:ifx,nx),k)
     &       -13.*s(iqx(1:ifx,nx),k)
     &       +16.*snew(iqx(1:ifx,nx),k))/14.  
c           cc=(3*sb(iq,k)-2*sa(iq,k)-5*s(iq,k)+4*snew(iq,k))/14.  
            snew(iqx(1:ifx,nx),k)=aa(iqx(1:ifx,nx))
     &       +axel*bb(iqx(1:ifx,nx))          
!          endif
!          enddo
!         endif  ! meth=5
          case(6)
!         if(iter>=3.and.meth==6)then   ! wqls again
!          do iq=1,ifull
!           if(mask(iq)==nx)then
            aa(iqx(1:ifx,nx))=(2.*sb(iqx(1:ifx,nx),k)
     &       -6.*sa(iqx(1:ifx,nx),k)+6*s(iqx(1:ifx,nx),k)
     &       +68.*snew(iqx(1:ifx,nx),k))/70.  
            bb(iqx(1:ifx,nx))=(5.*sb(iqx(1:ifx,nx),k)
     &       -8.*sa(iqx(1:ifx,nx),k)-13*s(iqx(1:ifx,nx),k)
     &       +16.*snew(iqx(1:ifx,nx),k))/14.  
            cc(iqx(1:ifx,nx))=(3.*sb(iqx(1:ifx,nx),k)
     &       -2.*sa(iqx(1:ifx,nx),k)-5*s(iqx(1:ifx,nx),k)
     &       +4.*snew(iqx(1:ifx,nx),k))/14.  
            snew(iqx(1:ifx,nx),k)=aa(iqx(1:ifx,nx))
     &       +axel*(bb(iqx(1:ifx,nx))
     &       +axel*cc(iqx(1:ifx,nx)))          
!          endif
!          enddo
!         endif  ! meth=5
         end select
         end if

!         do iq=1,ifull
!          if(mask(iq)==nx)then
c           rotate s files
            sb(iqx(1:ifx,nx),k)=sa(iqx(1:ifx,nx),k)
            sa(iqx(1:ifx,nx),k)=s(iqx(1:ifx,nx),k)
            s(iqx(1:ifx,nx),k)=snew(iqx(1:ifx,nx),k)
!          endif
!         enddo
        enddo  ! nx loop
        iters(k)=iter
      enddo ! k loop   
      if(ntest>0.and.meth>=1.and.mydiag)
     &  write (6,"('Iter ,s',i4,4f14.5)") iter,(s(iq,1),iq=1,4)
      if(iter==1)then
       do k=1,klim
         smax(k) = maxval(s(1:ifull,k))
         smin(k) = minval(s(1:ifull,k))
       end do
        call MPI_AllReduce( smax, smax_g, klim, MPI_REAL, MPI_MAX,
     &                    MPI_COMM_WORLD, ierr )
        call MPI_AllReduce( smin, smin_g, klim, MPI_REAL, MPI_MIN,
     &                    MPI_COMM_WORLD, ierr )
      endif
      if (iter>=itstest) then
     
      do k=1,klim
       dsolmax(k) = maxval(abs(dsol(1:ifull,k)))
      enddo
      call MPI_AllReduce( dsolmax, dsolmax_g, klim, MPI_REAL, MPI_MAX, 
     &                    MPI_COMM_WORLD, ierr )
      if(myid==0.and.ntest>0)then
        print *,'smin_g ',smin_g(:)
        print *,'smax_g ',smax_g(:)
        print *,'dsolmax_g ',dsolmax_g(:)
      endif  ! (myid==0)
      klimnew=klim
      do k=klim,1,-1
       if(dsolmax_g(k)<restol*(smax_g(k)-smin_g(k)))then
c        print *,'k,klim,iter,restol ',k,klim,iter,restol
         klimnew=k
       endif
      enddo
      klim=klimnew
      !call MPI_Bcast(klim,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      
      its1=its2
      ctst1=ctst2
      its2=iter
      ctst2=dsolmax_g(1)
      mm=(ctst2-ctst1)/real(its2-its1)
      ww=ctst1-mm*real(its1)
      yy=restol*(smax_g(1)-smin_g(1))
      itstest=nint((yy-ww)/mm)
      itstest=max(itstest,iter+1)
     
      end if ! iter>=itstest
      iter = iter + 1
      enddo   ! while( iter<itmax .and. klim>1)

      if(myid==0.and.(diag.or.ktau<6))then
        do k=1,kl
         print*,'helmjlm ktau,k,Iterations ',ktau,k,iters(k)
        enddo
      endif
      call end_log(helm_end) ! MJT
      return

5     continue  
      klim=kl
      iter = 1
      do while ( iter<itmax .and. klim>1)
       if(ntest==1.and.diag)print *,'myid,iter a ',myid,iter
       call bounds(s, klim=klim)
       if(ntest==1.and.diag)print *,'myid,iter b ',myid,iter
       do k=1,klim        
        do nx=1,nx_max
        !------------------------------------------------------------
        ! MJT pack
          ifx=ifullx(nx)
          dsol(iqx(1:ifx,nx),k)=
     &       ( zzn(iqx(1:ifx,nx))*s(iqn(1:ifx,nx),k)
     &       + zzw(iqx(1:ifx,nx))*s(iqw(1:ifx,nx),k)
     &       + zze(iqx(1:ifx,nx))*s(iqe(1:ifx,nx),k)
     &       + zzs(iqx(1:ifx,nx))*s(iqs(1:ifx,nx),k)
     &       +(zz(iqx(1:ifx,nx))
     &       -helm(iqx(1:ifx,nx),k))*s(iqx(1:ifx,nx),k)
     &       - rhs(iqx(1:ifx,nx),k))      
     &       /(helm(iqx(1:ifx,nx),k)-zz(iqx(1:ifx,nx)))
          snew(iqx(1:ifx,nx),k) = s(iqx(1:ifx,nx),k)
     &       + accel(k)*dsol(iqx(1:ifx,nx),k)
          s(iqx(1:ifx,nx),k)=snew(iqx(1:ifx,nx),k)
         !------------------------------------------------------------
        enddo  ! nx loop
        iters(k)=iter
      enddo ! k loop   
      if(ntest>0.and.diag)
     &  write (6,"('myid,Iter ,s',i4,4f14.5)")myid,iter,(s(iq,1),iq=1,4)
      if(iter==1)then
        do k=1,klim
         smax(k) = maxval(s(1:ifull,k))
         smin(k) = minval(s(1:ifull,k))
!         if(ntest>0.and.diag)print *,'myid,k,smax,smin ',
!     &                                myid,k,smax(k),smin(k)
        enddo
        if(ntest>0.and.diag)print *,' before smax call myid ', myid
        call MPI_AllReduce( smax, smax_g, klim, MPI_REAL, MPI_MAX,
     &                      MPI_COMM_WORLD, ierr )
        if(ntest>0.and.diag)print *,' before smin call myid ', myid
        call MPI_AllReduce( smin, smin_g, klim, MPI_REAL, MPI_MIN,
     &                      MPI_COMM_WORLD, ierr )
        if(ntest>0.and.myid==0)then
          print *,'ktau,myid,smin_g ',ktau,myid,smin_g(:)
          print *,'ktau,myid,smax_g ',ktau,myid,smax_g(:)
        endif  ! (myid==0)
      end if
      if (iter>=itstest) then
      
      do k=1,klim
c      write (6,"('iter,k ,s',2i4,4f14.5)") iter,k,(s(iq,k),iq=1,4)
       dsolmax(k) = maxval(abs(dsol(1:ifull,k)))
      enddo  ! k loop
      call MPI_AllReduce( dsolmax, dsolmax_g, klim, MPI_REAL, MPI_MAX, 
     &                    MPI_COMM_WORLD, ierr )
      if(ntest>0)then
        print *,'ktau,myid,iter,dsolmax ',ktau,myid,iter,dsolmax(:)
      endif  ! (myid==0)
      klimnew=klim
      do k=klim,1,-1
       if(dsolmax_g(k)<restol*(smax_g(k)-smin_g(k)))then
         klimnew=k
       endif
      enddo
      if(ntest>0)print *,'ktau,myid,klim,klimnew ',
     &                    ktau,myid,klim,klimnew
      klim=klimnew
    !  call MPI_AllReduce( klimnew, klim, 1, MPI_INTEGER, MPI_MAX,
    ! &                    MPI_COMM_WORLD, ierr )
     
      its1=its2
      ctst1=ctst2
      its2=iter
      ctst2=dsolmax_g(1)
      mm=(ctst2-ctst1)/real(its2-its1)
      ww=ctst1-mm*real(its1)
      yy=restol*(smax_g(1)-smin_g(1))
      itstest=nint((yy-ww)/mm)
      itstest=max(itstest,iter+1)
    
      end if ! iter>itstest
      iter = iter + 1
      enddo   ! while( iter<itmax .and. klim>1)

      if(myid==0)then
        if(nmaxpr==1)print*,'helmjlm ktau,k,Iterations ',ktau,1,iters(1)
        if(diag.or.ktau<6)then
         do k=2,kl
          print*,'helmjlm ktau,k,Iterations ',ktau,k,iters(k)
         enddo
        endif
      endif
      
      call end_log(helm_end) ! MJT
      return
      end
