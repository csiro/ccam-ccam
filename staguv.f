      subroutine staguv(u,v,uout,vout)
!     stripped down version just with nstag=nstagu=-3, mstagpt=-3      
!     staguv    may be called from adjust5, upglobal
!     unstaguv  may be called from adjust5,  nonlin
c     nstag now in parm.h  ! for nstag=0   staguv uses cubic interpolation
c                          ! for nstag=1   staguv goes via cartesian components
c                          ! for nstag=2   staguv linear with em terms
c                          ! for nstag=3-5 staguv: jlm reversible 2/1/98
c                          ! for nstag=7   staguv: jlm reversible Akima
c                          ! for nstag=10  staguv: jlm rev_stag
c                          !    only -3 used nowadays
c     nstagu now in parm.h ! same but for unstaguv
c     N.B. staguv & unstaguv previously required 2D arrays as input
c     - no longer, as transferred here to uin and vin
c     assume k level already given
      use cc_mpi
      implicit none
      include 'newmpar.h'
      include 'indices.h' ! in,is,iw,ie,inn,iss,iww,iee
      include 'map.h'
      include 'parm.h'
      include 'parmdyn.h'
      include 'vecsuv.h'   ! vecsuv info
      include 'vecsuva.h'  ! vecsuva info
      real, dimension(ifull,kl), intent(in)  :: u, v
      real, dimension(ifull,kl), intent(out) :: uout, vout
      real, dimension(ifull+iextra,kl) :: ua, va, ub, vb, ud, vd,
     &                                    uin, vin
      integer, parameter :: itnmax=3
      integer :: i, iq, itn, j, k

      call start_log(stag_begin)
c     unstaggered u & v as input; staggered as output

      ! Copying could be avoided if input arrays were dimensioned ifull+iextra
      do k=1,kl
         do iq=1,ifull
            uin(iq,k) = u(iq,k)
            vin(iq,k) = v(iq,k)
         end do
      end do
      call boundsuv(uin,vin)

      ! This is the nstag=3, mstagpt=-3 code.
      ! The MPI version can't have index arrays on the LHS.

      do k=1,kl
!cdir nodep
         do iq=1,ifull          ! precalculate rhs terms
            ud(iq,k)= uin(iwu(iq),k)/10.+uin(iq,k)+uin(ieu(iq),k)/2.
            vd(iq,k)= vin(isv(iq),k)/10.+vin(iq,k)+vin(inv(iq),k)/2.
         enddo
      enddo

      call boundsuv(ud,vd)

      do k=1,kl
!cdir nodep
         do iq=1,ifull
            ua(iq,k)=ud(iq,k)-ud(iwu(iq),k)/2. ! 1st guess
            va(iq,k)=vd(iq,k)-vd(isv(iq),k)/2. ! 1st guess
         enddo
      enddo

      do itn=1,itnmax           ! each loop is a double iteration
         call boundsuv(ua,va)
         do k=1,kl
!cdir nodep
            do iq=1,ifull
               ub(iq,k)=ua(iwu(iq),k)
               vb(iq,k)=va(isv(iq),k)
            enddo
         enddo
         call boundsuv(ub,vb)
         do k=1,kl
!cdir nodep
            do iq=1,ifull
               uin(iq,k)=(ud(iq,k)-.5*ud(iwu(iq),k)
     &                 -ua(ieu(iq),k)/10. +ub(iwu(iq),k)/4.)/.95
               vin(iq,k)=(vd(iq,k)-.5*vd(isv(iq),k)
     &                 -va(inv(iq),k)/10. +vb(isv(iq),k)/4.)/.95
            enddo
         enddo
         call boundsuv(uin,vin)
         do k=1,kl
!cdir nodep
            do iq=1,ifull
               ub(iq,k)=uin(iwu(iq),k)
               vb(iq,k)=vin(isv(iq),k)
            enddo
         end do
         call boundsuv(ub,vb)
!cdir nodep
         do k=1,kl
            do iq=1,ifull
               ua(iq,k)=(ud(iq,k)-.5*ud(iwu(iq),k)
     &              -uin(ieu(iq),k)/10. +ub(iwu(iq),k)/4.)/.95
               va(iq,k)=(vd(iq,k)-.5*vd(isv(iq),k)
     &              -vin(inv(iq),k)/10. +vb(isv(iq),k)/4.)/.95
            enddo
         end do
      end do                     ! itn=1,itnmax
	 
      do k=1,kl
         do iq=1,ifull          ! final values for output
            uout(iq,k) = ua(iq,k)
            vout(iq,k) = va(iq,k)
         end do
      end do

      call end_log(stag_end)
      return

      end


      subroutine unstaguv(u,v,uout,vout)
c     staggered u & v as input; unstaggered as output
      use cc_mpi
      implicit none
      include 'newmpar.h'
      include 'indices.h' ! in,is,iw,ie,inn,iss,iww,iee
      include 'map.h'
      include 'parm.h'
      include 'parmdyn.h'
      include 'vecsuv.h'   ! vecsuv info
      include 'vecsuva.h'  ! vecsuva info
      real, dimension(ifull,kl), intent(in)  :: u, v
      real, dimension(ifull,kl), intent(out) :: uout, vout
      real, dimension(ifull+iextra,kl) :: ua, va, ub, vb, ud, vd,
     &                                    uin, vin
      integer, parameter :: itnmax=3
      integer :: i, iq, itn, j, k

      call start_log(stag_begin)
      do k=1,kl
         do iq=1,ifull
            uin(iq,k) = u(iq,k)
            vin(iq,k) = v(iq,k)
         end do
      end do
      call boundsuv(uin,vin)

      do k=1,kl
!cdir nodep
         do iq=1,ifull          ! precalculate rhs terms
            ud(iq,k)= uin(ieu(iq),k)/10.+uin(iq,k)+uin(iwu(iq),k)/2.
            vd(iq,k)= vin(inv(iq),k)/10.+vin(iq,k)+vin(isv(iq),k)/2.
         enddo
      enddo
      call boundsuv(ud,vd)
      do k=1,kl
!cdir nodep
         do iq=1,ifull
            ua(iq,k)=ud(iq,k)-ud(ieu(iq),k)/2. ! 1st guess
            va(iq,k)=vd(iq,k)-vd(inv(iq),k)/2. ! 1st guess
         enddo
      enddo

      do itn=1,itnmax           ! each loop is a double iteration
         call boundsuv(ua,va)
         do k=1,kl
!cdir nodep
            do iq=1,ifull
               ub(iq,k)=ua(ieu(iq),k)
               vb(iq,k)=va(inv(iq),k)
            enddo
         enddo
         call boundsuv(ub,vb)
         do k=1,kl
!cdir nodep
            do iq=1,ifull
               uin(iq,k)=(ud(iq,k)-.5*ud(ieu(iq),k)
     &              -ua(iwu(iq),k)/10. +ub(ieu(iq),k)/4.)/.95
               vin(iq,k)=(vd(iq,k)-.5*vd(inv(iq),k)
     &              -va(isv(iq),k)/10. +vb(inv(iq),k)/4.)/.95
            enddo
         enddo
         call boundsuv(uin,vin)
         do k=1,kl
!cdir nodep
            do iq=1,ifull
               ub(iq,k)=uin(ieu(iq),k)
               vb(iq,k)=vin(inv(iq),k)
            enddo
         enddo
         call boundsuv(ub,vb)
         do k=1,kl
!cdir nodep
            do iq=1,ifull
               ua(iq,k)=(ud(iq,k)-.5*ud(ieu(iq),k)
     &              -uin(iwu(iq),k)/10. +ub(ieu(iq),k)/4.)/.95
               va(iq,k)=(vd(iq,k)-.5*vd(inv(iq),k)
     &              -vin(isv(iq),k)/10. +vb(inv(iq),k)/4.)/.95
            enddo
         enddo
      enddo                     ! itn=1,itnmax
      
      do k=1,kl
         do iq=1,ifull          ! final values for output
            uout(iq,k) = ua(iq,k)
            vout(iq,k) = va(iq,k)
         end do
      end do
      call end_log(stag_end)
      return
      end        
