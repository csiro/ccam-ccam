      subroutine helmsor(zz,zzn,zze,zzw,zzs,helm,s,irhs)
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
      real, dimension(ifull), intent(in) :: zz,zzn,zze,zzw,zzs
      real, dimension(ifull,kl), intent(in) :: helm         ! Helmholtz coefficients
      real, dimension(ifull+iextra,kl), intent(inout) :: s  ! Solution
      real, dimension(ifull,kl), intent(in) :: irhs         ! RHS
      real, dimension(ifull,kl) :: rhs
      real, dimension(ifullx,kl,maxcolour) :: helmc,rhsc
      real, dimension(ifull,kl) :: sb, sa, snew
      real, dimension(ifull,kl) :: dsol
      real, dimension(ifullx,maxcolour) :: zznc,zzec
      real, dimension(ifullx,maxcolour) :: zzwc,zzsc
      real, dimension(kl) ::  dsolmax, dsolmax_g, smax, smax_g
      real, dimension(kl) ::  smin, smin_g, savg
      real, dimension(:), allocatable, save :: accel
      real, dimension(ifull) :: aa, bb, cc
      real axel
      real, save :: dtsave = 0.
      integer iq, iter, k, nx, j, jx, i, klim, klimnew, ierr, meth
      integer nx_max, iqg, ig, jg, ng, tg, n
      integer, dimension(kl) :: iters
      integer, dimension(1) :: idum
      save  meth, nx_max, axel

#include "log.h"

      START_LOG(helm)
      
      rhs=irhs ! allows subroutine to modify rhs
      
      if (dt/=dtsave) then
        dtsave=dt
        if (.not.allocated(accel)) then
          allocate(accel(kl))
        end if

       if(precon==-1)precon=-2325  ! i.e. 2, 3, .25
       nx_max=abs(precon)/1000
       meth=abs(precon)/100-10*nx_max
       axel=-.01*real(precon)-10*nx_max-meth
       if(myid==0)write(6,*)'in helmsor nx_max,meth,axel: ',
     &                                   nx_max,meth,axel
 
       if (nx_max/=maxcolour) then
        if (myid==0) then
         write(6,*) "WARN: mismatched number of colours"
         write(6,*) "changing nx_max ",nx_max,maxcolour
        end if
        nx_max=maxcolour
       end if
  
       if (il_g<=200) then
         ! usual
         do k=1,kl
          call optmx(il_g,schmidt,dt,bam(k),accel(k))
          ! MJT - not sure about the following line
          accel(k)=1.+.55*(accel(k)-1.) ! MJT suggestion
!         if(il_g>il)accel(k)=1.+.55*(accel(k)-1.)  ! for uniform-dec 22/4/08
c         if(il_g==il)accel(k)=1.+.55*(accel(k)-1.) ! just a test
          if(myid==0)write(6,*)'k,accel ',k,accel(k)
         enddo
       else
         ! large grid
         accel=1. ! gauss-seidel
       end if
      endif

      if(precon>=-2899) then  ! e.g. not -2900 or -3900
      klim=kl
      iter = 1
      do while ( iter<itmax .and. klim>1)
       call bounds(s, klim=klim)
       do k=1,klim        
        do nx=1,nx_max
           dsol(iqx(:,nx),k)=
     &        ( zzn(iqx(:,nx))*s(iqn(:,nx),k)
     &        + zzw(iqx(:,nx))*s(iqw(:,nx),k)
     &        + zze(iqx(:,nx))*s(iqe(:,nx),k)
     &        + zzs(iqx(:,nx))*s(iqs(:,nx),k)
     &        +(zz(iqx(:,nx))
     &        -helm(iqx(:,nx),k))*s(iqx(:,nx),k)
     &        - rhs(iqx(:,nx),k))      
     &        /(helm(iqx(:,nx),k)-zz(iqx(:,nx)))
           snew(iqx(:,nx),k) = s(iqx(:,nx),k)
     &        + dsol(iqx(:,nx),k)

!        following are jlm methods for improving guess
         if(iter>=3)then
         select case(meth)
          case(3)
!         if(iter>=3.and.meth==3)then   ! qls
            aa(iqx(:,nx))=
     &       (sb(iqx(:,nx),k)
     &       -3.*sa(iqx(:,nx),k)
     &       +3.*s(iqx(:,nx),k)
     &       +19.*snew(iqx(:,nx),k))/20.  
            bb(iqx(:,nx))=(9.*sb(iqx(:,nx),k)
     &       -17.*sa(iqx(:,nx),k)
     &       -13.*s(iqx(:,nx),k)
     &       +21.*snew(iqx(:,nx),k))/20.  
c           snew(iq,k)=aa+.25*bb+.0625*cc !3c
            snew(iqx(:,nx),k)=aa(iqx(:,nx))
     &       +axel*bb(iqx(:,nx))
!         endif  ! meth=3
          case(4)
!         if(iter>=3.and.meth==4)then   ! oscill
            aa(iqx(:,nx))=(7.*snew(iqx(:,nx),k)
     &       +3.*s(iqx(:,nx),k)
     &       -3.*sa(iqx(:,nx),k)
     &       +sb(iqx(:,nx),k))/8. ! oscill
            bb(iqx(:,nx))=snew(iqx(:,nx),k)
     &       -.5*s(iqx(:,nx),k)
     &       -sa(iqx(:,nx),k)
     &       +.5*sb(iqx(:,nx),k)         ! oscill
c           cc=.25*(snew(iq,k)-s(iq,k)-sa(iq,k)+sb(iq,k))         !             aa=(sb(iq,k)-3*sa(iq,k)+3*s(iq,k)+19*snew(iq,k))/20.  
            snew(iqx(:,nx),k)=aa(iqx(:,nx))
     &       +axel*bb(iqx(:,nx))
!         endif  ! meth=4
          case(5)
!         if(iter>=3.and.meth==5)then   ! wqls
            aa(iqx(:,nx))=(2.*sb(iqx(:,nx),k)
     &       -6.*sa(iqx(:,nx),k)+6*s(iqx(:,nx),k)
     &       +68.*snew(iqx(:,nx),k))/70.  
            bb(iqx(:,nx))=(5.*sb(iqx(:,nx),k)
     &       -8.*sa(iqx(:,nx),k)
     &       -13.*s(iqx(:,nx),k)
     &       +16.*snew(iqx(:,nx),k))/14.  
c           cc=(3*sb(iq,k)-2*sa(iq,k)-5*s(iq,k)+4*snew(iq,k))/14.  
            snew(iqx(:,nx),k)=aa(iqx(:,nx))
     &       +axel*bb(iqx(:,nx))          
!         endif  ! meth=5
          case(6)
!         if(iter>=3.and.meth==6)then   ! wqls again
            aa(iqx(:,nx))=(2.*sb(iqx(:,nx),k)
     &       -6.*sa(iqx(:,nx),k)+6*s(iqx(:,nx),k)
     &       +68.*snew(iqx(:,nx),k))/70.  
            bb(iqx(:,nx))=(5.*sb(iqx(:,nx),k)
     &       -8.*sa(iqx(:,nx),k)-13*s(iqx(:,nx),k)
     &       +16.*snew(iqx(:,nx),k))/14.  
            cc(iqx(:,nx))=(3.*sb(iqx(:,nx),k)
     &       -2.*sa(iqx(:,nx),k)-5*s(iqx(:,nx),k)
     &       +4.*snew(iqx(:,nx),k))/14.  
            snew(iqx(:,nx),k)=aa(iqx(:,nx))
     &       +axel*(bb(iqx(:,nx))
     &       +axel*cc(iqx(:,nx)))          
!         endif  ! meth=5
         end select
         end if

            sb(iqx(:,nx),k)=sa(iqx(:,nx),k)
            sa(iqx(:,nx),k)=s(iqx(:,nx),k)
            s(iqx(:,nx),k)=snew(iqx(:,nx),k)
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
      do k=1,klim
       smax(k) = maxval(s(1:ifull,k))
       smin(k) = minval(s(1:ifull,k))
!       if(ntest>0.and.diag)print *,'myid,k,smax,smin ',
!     &                              myid,k,smax(k),smin(k)
      enddo
      if(ntest>0.and.diag)write(6,*)' before smax call myid ',myid
      call ccmpi_allreduce(smax(1:klim),smax_g(1:klim),"max",
     &                     comm_world)
      if(ntest>0.and.diag)write(6,*)' before smin call myid ',myid
      call ccmpi_allreduce(smin(1:klim),smin_g(1:klim),"min",
     &                     comm_world)
      if((ntest>0.or.nmaxpr==1).and.myid==0)then
        write(6,*)'ktau,smin_g ',ktau,smin_g(:)
        write(6,*)'ktau,smax_g ',ktau,smax_g(:)
      endif  ! (myid==0)

      ! JLM suggestion
      do k=1,kl
        savg(k)=0.5*(smax_g(k)+smin_g(k))
        s(1:ifull,k)=s(1:ifull,k)-savg(k)
        rhs(:,k)=rhs(:,k)+(helm(:,k)-zz-zzn-zzs-zze-zzw)*savg(k)
      end do

      do nx=1,nx_max
        zznc(:,nx) =zzn(iqx(:,nx))
        zzwc(:,nx) =zzw(iqx(:,nx))
        zzec(:,nx) =zze(iqx(:,nx))
        zzsc(:,nx) =zzs(iqx(:,nx))
        do k=1,kl
          helmc(:,k,nx)=helm(iqx(:,nx),k)-zz(iqx(:,nx))
          rhsc(:,k,nx) =rhs(iqx(:,nx),k)
        end do
      end do
      
      call bounds(s)     
      do while ( iter<itmax .and. klim>0)
       do nx=1,nx_max
        do k=1,klim
          dsol(iqx(:,nx),k)=
     &       ( zznc(:,nx)*s(iqn(:,nx),k)
     &       + zzwc(:,nx)*s(iqw(:,nx),k)
     &       + zzec(:,nx)*s(iqe(:,nx),k)
     &       + zzsc(:,nx)*s(iqs(:,nx),k)
     &       -helmc(:,k,nx)*s(iqx(:,nx),k)
     &       -rhsc(:,k,nx) )/helmc(:,k,nx)
          s(iqx(:,nx),k) = s(iqx(:,nx),k)
     &       + accel(k)*dsol(iqx(:,nx),k)
        enddo ! k loop
        call bounds_colour(s, nx, klim=klim)
       enddo  ! nx loop  
       do k=1,klim
        iters(k)=iter
       end do
!       if((ntest>0.or.nmaxpr==1).and.diag)
!     &  write (6,"('myid,Iter ,s',i4,4f14.5)")
!     &  myid,iter,(s(iq,1),iq=1,4)
       
       
       do k=1,klim
c       write (6,"('iter,k ,s',2i4,4f14.5)") iter,k,(s(iq,k),iq=1,4)
        dsolmax(k) = maxval(abs(dsol(1:ifull,k)))
       enddo  ! k loop
!       if(ntest>0)then
!        write(6,*)'ktau,myid,iter,dsolmax ',ktau,myid,iter,dsolmax(:)
!       endif  ! (myid==0)
       klimnew=klim
       call ccmpi_allreduce(dsolmax(1:klim),dsolmax_g(1:klim),"max",
     &                      comm_world)
       do k=klim,1,-1
        if(dsolmax_g(k)<restol*(smax_g(k)-smin_g(k)))then
          klimnew=k-1
        endif
       enddo
       klim=klimnew

       iter = iter + 1
      enddo   ! while( iter<itmax .and. klim>1)
      
      do k=1,kl
        s(1:ifull,k)=s(1:ifull,k)+savg(k)
      end do

      if(myid==0)then
        if(nmaxpr==1) then
          write(6,*)'helmjlm ktau,k,Iterations ',ktau,1,iters(1)
        end if
        if(diag.or.ktau<6.or.iters(1)>itmax-5)then
         do k=1,kl
          write(6,*)'helmjlm ktau,k,Iterations ',ktau,k,iters(k)
         enddo
         !write(6,*) "itc ",itc
        endif
      endif
      
      end if ! precon>=-2899 ..else..
      
      END_LOG(helm)
      return
      end
