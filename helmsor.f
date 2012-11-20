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
!     integer, parameter :: meth=3      ! 3b, 5c1 good
      integer, parameter :: ntest=0 
      integer, parameter :: itmax=400 ! maximum number of iterations allowed
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
      integer, dimension(3), save :: ifullx
      integer, dimension(:,:), allocatable, save :: iqx,iqn,iqe
      integer, dimension(:,:), allocatable, save :: iqw,iqs
      real, dimension(kl) ::  dsolmax, dsolmax_g, smax, smax_g
      real, dimension(kl) ::  smin, smin_g
      real, dimension(:), allocatable, save :: accel
      real aa(ifull), bb(ifull), cc(ifull), axel
      real gd, ci, itserr1, itserr2
      integer iq, iter, k, nx, j, jx, i, klim, klimnew, ierr, meth
      integer ifx, nx_max, iqg, ig, jg, ng, tg, n
      integer itstest, itc, itsave1, itsave2
      integer, dimension(1) :: idum
      save  meth, nx_max, axel

      call start_log(helm_begin)
      
      if (.not.allocated(mask)) then
        allocate(mask(ifull))
        allocate(iqx(ifull,3),iqn(ifull,3),iqe(ifull,3))
        allocate(iqw(ifull,3),iqs(ifull,3),accel(kl))

        if(precon==-1)precon=-2325  ! i.e. 2, 3, .25
        nx_max=abs(precon)/1000
        meth=abs(precon)/100-10*nx_max
        axel=-.01*real(precon)-10*nx_max-meth
        if(myid==0)write(6,*)'in helmsor nx_max,meth,axel: ',
     &                                   nx_max,meth,axel
        
        ! global index method
        if (nx_max==3) then ! 3 colour

        mask=0
        do n=1,npan
         do j=1,jpan
          do i=1,ipan
           iq = indp(i,j,n)   ! Local
           iqg = indg(i,j,n)  ! Global
           
           tg=iqg-1
           ng=tg/(il_g*il_g)
           tg=tg-ng*il_g*il_g
           jg=tg/il_g
           tg=tg-jg*il_g
           ig=tg
           ig=ig+1
           jg=jg+1
           
           jx=mod(ig+jg+ng*il_g,2)
           select case(ng)
            case(0,3)
             if (jx==0) then
              mask(iq)=1
             else
              mask(iq)=2
             end if
            case(1,4)
             if (jx==0) then
              mask(iq)=1
             else
              mask(iq)=3
             end if
            case(2,5)
             if (jx==0) then
              mask(iq)=2
             else
              mask(iq)=3
             end if
           end select
          end do
         end do
        end do
      
       else ! 2 colour

        mask=0
        do n=1,npan
         do j=1,jpan
          do i=1,ipan
           iq = indp(i,j,n)   ! Local
           iqg = indg(i,j,n)  ! Global
           
           tg=iqg-1
           ng=tg/(il_g*il_g)
           tg=tg-ng*il_g*il_g
           jg=tg/il_g
           tg=tg-jg*il_g
           ig=tg
           ig=ig+1
           jg=jg+1
           
           jx=mod(ig+jg+ng*il_g,2)
           if (jx==0) then
            mask(iq)=1
           else
            mask(iq)=2
           end if
          end do
         end do
        end do

       end if ! if (nx_max==3) ..else..

     
!        ! old local indices
!        mask(:)=1
!        do j=1,jl
!         do i=1,il,2
!          jx=mod(i+j,2)    ! 0 or 1  starting (1,1) on panel
!          iq=i+jx+(j-1)*il
!          if(nx_max==2)then
!            mask(iq)=2
!          else
!            if(j>il.and.j<=3*il)mask(iq)=3
!            if(j>3*il.and.j<=5*il)mask(iq)=2
!            jx=mod(i+j+1,2)  ! 0 or 1  starting (2,1) on panel
!            iq=i+jx+(j-1)*il
!            if(j<=2*il)mask(iq)=2
!            if(j>4*il)mask(iq)=3
!          endif 
!         enddo
!       enddo
       
       ! Pack colour indices
       ifullx=0
       iqx=0
       iqn=0
       iqe=0
       iqw=0
       iqs=0
       do iq=1,ifull
         ifullx(mask(iq))=ifullx(mask(iq))+1
         iqx(ifullx(mask(iq)),mask(iq))=iq
         iqn(ifullx(mask(iq)),mask(iq))=in(iq)
         iqe(ifullx(mask(iq)),mask(iq))=ie(iq)
         iqw(ifullx(mask(iq)),mask(iq))=iw(iq)
         iqs(ifullx(mask(iq)),mask(iq))=is(iq)
       end do

       do k=1,kl
        !MJT - issues with convergence.  Use Gauss-Seidel
        ! until accel can be re-calculated
        !call optmx(il_g,schmidt,dt,bam(k),accel(k))
        accel(k)=1. ! gauss-seidel
        
        ! MJT - not sure about the following line
!       if(il_g>il)accel(k)=1.+.55*(accel(k)-1.)  ! for uniform-dec 22/4/08
c       if(il_g==il)accel(k)=1.+.55*(accel(k)-1.) ! just a test
        if(myid==0)write(6,*)'k,accel ',k,accel(k)
       enddo
      endif

      if(precon>=-2899) then  ! e.g. not -2900 or -3900
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
       call ccmpi_reduce(smax(1:klim),smax_g(1:klim),"max",0,
     &                   comm_world)
       call ccmpi_reduce(smin(1:klim),smin_g(1:klim),"min",0,
     &                   comm_world)
      endif
      do k=1,klim
       dsolmax(k) = maxval(abs(dsol(1:ifull,k)))
      enddo
      call ccmpi_reduce(dsolmax(1:klim),dsolmax_g(1:klim),"max",0,
     &                  comm_world)
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
      idum(1)=klim
      call ccmpi_bcast(idum(1:1),0,comm_world)
      klim=idum(1)
      iter = iter + 1
      enddo   ! while( iter<itmax .and. klim>1)

      if(myid==0.and.(diag.or.ktau<6))then
        do k=1,kl
         print*,'helmjlm ktau,k,Iterations ',ktau,k,iters(k)
        enddo
      endif

      else ! e.g. -2900 or -3900

      klim=kl
      iter = 1
      itsave2=0
      itserr2=9.E9
      itstest=1
      itc=0
      call bounds(s, klim=klim)      
      do while ( iter<itmax .and. klim>0)
      do nx=1,nx_max
        ifx=ifullx(nx)
        do k=1,klim
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
        enddo ! k loop
        call bounds(s, klim=klim)  ! Need this to work for different colours
      enddo  ! nx loop  
      do k=1,klim
        iters(k)=iter
      end do
       if((ntest>0.or.nmaxpr==1).and.diag)
     &  write (6,"('myid,Iter ,s',i4,4f14.5)")myid,iter,(s(iq,1),iq=1,4)
       if(iter==1)then
        do k=1,klim
         smax(k) = maxval(s(1:ifull,k))
         smin(k) = minval(s(1:ifull,k))
!         if(ntest>0.and.diag)print *,'myid,k,smax,smin ',
!     &                                myid,k,smax(k),smin(k)
        enddo
        if(ntest>0.and.diag)write(6,*)' before smax call myid ',myid
        call ccmpi_allreduce(smax(1:klim),smax_g(1:klim),"max",
     &                       comm_world)
        if(ntest>0.and.diag)write(6,*)' before smin call myid ',myid
        call ccmpi_allreduce(smin(1:klim),smin_g(1:klim),"min",
     &                       comm_world)
        if((ntest>0.or.nmaxpr==1).and.myid==0)then
          write(6,*)'ktau,smin_g ',ktau,smin_g(:)
          write(6,*)'ktau,smax_g ',ktau,smax_g(:)
        endif  ! (myid==0)
       end if
       if (iter>=itstest) then
        do k=1,klim
c        write (6,"('iter,k ,s',2i4,4f14.5)") iter,k,(s(iq,k),iq=1,4)
         dsolmax(k) = maxval(abs(dsol(1:ifull,k)))
        enddo  ! k loop
        if(ntest>0)then
         write(6,*)'ktau,myid,iter,dsolmax ',ktau,myid,iter,dsolmax(:)
        endif  ! (myid==0)
        klimnew=klim
        call ccmpi_allreduce(dsolmax(1:klim),dsolmax_g(1:klim),"max",
     &                       comm_world)
        do k=klim,1,-1
         if(dsolmax_g(k)<restol*(smax_g(k)-smin_g(k)))then
           klimnew=k-1
         endif
        enddo
        klim=klimnew

        itsave1=itsave2
        itsave2=iter
        itserr1=itserr2
        itserr2=log10(dsolmax_g(1))
       
        gd=(itserr2-itserr1)/real(itsave2-itsave1)
        ci=itserr2-gd*real(itsave2)
        if (gd/=0.) then
         itstest=nint((log10(restol*(smax_g(1)-smin_g(1)))-ci)/gd)
         itstest=max(itstest,iter+1)
        else
         itstest=iter+1
        end if
        if (myid==0.and.nmaxpr==1) then
          write(6,*) "iter,itstest ",iter,itstest
        end if
        itc=itc+1
       end if ! iter>=itstest
       iter = iter + 1
      enddo   ! while( iter<itmax .and. klim>1)

      if(myid==0)then
        if(nmaxpr==1) then
          write(6,*)'helmjlm ktau,k,Iterations ',ktau,1,iters(1)
        end if
        if(diag.or.ktau<6)then
         do k=2,kl
          write(6,*)'helmjlm ktau,k,Iterations ',ktau,k,iters(k)
         enddo
         write(6,*) "itc ",itc
        endif
      endif
      
      end if ! precon>=-2899 ..else..
      
      call end_log(helm_end)
      return
      end
