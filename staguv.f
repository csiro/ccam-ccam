      subroutine staguv(u,v,uout,vout)
!     stripped down version just with nstag=nstagu=3, mstagpt=-3      
!     also now includes nstag=nstagu=4 & 5
!     staguv    may be called from adjust5, upglobal
!     unstaguv  may be called from adjust5,  nonlin
c     nstag now in parm.h  ! for nstag=3-5 staguv: jlm reversible 
c                          ! -ve switches every abs(nstag) time steps
c     nstagu now in parm.h ! same but for unstaguv
c     N.B. staguv & unstaguv previously required 2D arrays as input
c     - no longer, as transferred here to uin and vin
c     unstaggered u & v as input; staggered as output
      use cc_mpi
c     use diag_m             ! for calls to maxmin
      implicit none
      include 'newmpar.h'
      include 'indices.h' ! in,is,iw,ie,inn,iss,iww,iee
      include 'map.h'
      include 'parm.h'
      include 'parmdyn.h'
      include 'vecsuv.h'   ! vecsuv info
      real, dimension(ifull,kl), intent(in)  :: u, v
      real, dimension(ifull,kl), intent(out) :: uout, vout
      real, dimension(ifull+iextra,kl) :: ua, va, ud, vd,
     &                                    uin, vin
      integer, parameter :: ntest=0    ! usually 0, 1 for test prints
      integer, parameter :: itnmax=3
      integer :: iq, itn, k, i, j

      call start_log(stag_begin)
      if(nmaxpr==1.and.mydiag)then
        print *,'  stag_ktau,nstag,nstagu',ktau,nstag,nstagu
      endif

      ! Copying could be avoided if input arrays were dimensioned ifull+iextra
      do k=1,kl
         do iq=1,ifull
            uin(iq,k) = u(iq,k)
            vin(iq,k) = v(iq,k)
         end do
      end do
      
      if(abs(nstag)<3)then
        call boundsuv(uin,vin,nrows=2)
        do k=1,kl
         do iq=1,ifull
          uout(iq,k)=(9.*(uin(ieu(iq),k)+uin(iq,k))
     &                   -uin(iwu(iq),k)-uin(ieeu(iq),k))/16.
          vout(iq,k)=(9.*(vin(inv(iq),k)+vin(iq,k))
     &                   -vin(isv(iq),k)-vin(innv(iq),k))/16.
c          uout(iq,k)=.5*(uin(ieu(iq),k)+uin(iq,k))
c          vout(iq,k)=.5*(vin(inv(iq),k)+vin(iq,k))
         enddo   ! iq loop
        enddo
        return
      endif  ! (nstag==0)

      if ( nstag==3 ) then
         call boundsuv(uin,vin,nrows=2)
         if(ntest==1)then
           print *,'staguv diags'
           write (6,"(2x,4i8,6x,4i8)") (i,i=1,4),(i,i=1,4)
           do j=93,96
            write (6,"(i4,4f8.3,6x,4f8.3)") j,(u(i+(j-1)*il,4),i=1,4),
     &                                        (v(i+(j-1)*il,4),i=1,4)
           enddo          
           write (6,"(2x,4i8,6x,4i8)") (i,i=1,4),(i,i=1,4)
           do j=189,192
            write (6,"(i4,4f8.3,6x,4f8.3)") j,(u(i+(j-1)*il,4),i=1,4),
     &                                        (v(i+(j-1)*il,4),i=1,4)
           enddo          
           write (6,"(2x,4i8,6x,4i8)") (i,i=1,4),(i,i=1,4)
           do j=285,288
            write (6,"(i4,4f8.3,6x,4f8.3)") j,(u(i+(j-1)*il,4),i=1,4),
     &                                        (v(i+(j-1)*il,4),i=1,4)
           enddo          
           do j=95,288,96
            do i=1,2
             iq=i+(j-1)*il
             write (6,"('i,j,uin(ieu),uin(ieeu) ',2i4,2f8.3)")
     &                   i,j,uin(ieu(iq),4),uin(ieeu(iq),4) 
             write (6,"('i,j,uin(iwu),uin(iwwu) ',2i4,2f8.3)")
     &                   i,j,uin(iwu(iq),4),uin(iwwu(iq),4) 
             write (6,"('i,j,vin(inv),vin(innv) ',2i4,2f8.3)")
     &                   i,j,vin(inv(iq),4),vin(innv(iq),4) 
             write (6,"('i,j,vin(isv),vin(issv) ',2i4,2f8.3)")
     &                   i,j,vin(isv(iq),4),vin(issv(iq),4) 
            enddo
           enddo
         endif  ! (ntest==1)
         
         do k=1,kl
!cdir nodep
            do iq=1,ifull       ! precalculate rhs terms with iwwu2 & issv2
               ud(iq,k)=uin(iq,k)/2.+uin(ieu(iq),k)+uin(ieeu(iq),k)/10.
               vd(iq,k)=vin(iq,k)/2.+vin(inv(iq),k)+vin(innv(iq),k)/10.
            end do
         end do

         call boundsuv(ud,vd)
         do k=1,kl
!cdir nodep
            do iq=1,ifull
               ua(iq,k)=ud(iq,k)-ud(ieu(iq),k)/2. ! 1st guess
               va(iq,k)=vd(iq,k)-vd(inv(iq),k)/2. ! 1st guess
            end do
         end do

         do itn=1,itnmax        ! each loop is a double iteration

            call boundsuv(ua,va,nrows=2)
            do k=1,kl
!cdir nodep
               do iq=1,ifull
                  uin(iq,k)=(ud(iq,k)-.5*ud(ieu(iq),k)
     &                 -ua(iwu(iq),k)/10. +ua(ieeu(iq),k)/4.)/.95
                  vin(iq,k)=(vd(iq,k)-.5*vd(inv(iq),k)
     &                 -va(isv(iq),k)/10. +va(innv(iq),k)/4.)/.95
               end do
            end do

            call boundsuv(uin,vin,nrows=2)
            do k=1,kl
!cdir nodep
               do iq=1,ifull
                  ua(iq,k)=(ud(iq,k)-.5*ud(ieu(iq),k)
     &                 -uin(iwu(iq),k)/10. +uin(ieeu(iq),k)/4.)/.95
                  va(iq,k)=(vd(iq,k)-.5*vd(inv(iq),k)
     &                 -vin(isv(iq),k)/10. +vin(innv(iq),k)/4.)/.95
               end do
            end do

         end do                  ! itn=1,itnmax

      else if ( nstag==4 ) then
         call boundsuv(uin,vin)
         do k=1,kl
!cdir nodep
            do iq=1,ifull       ! precalculate rhs terms
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

         do itn=1,itnmax        ! each loop is a double iteration
            call boundsuv(ua,va,nrows=2)

            do k=1,kl
!cdir nodep
               do iq=1,ifull
                  uin(iq,k)=(ud(iq,k)-.5*ud(iwu(iq),k)
     &                 -ua(ieu(iq),k)/10. +ua(iwwu(iq),k)/4.)/.95
                  vin(iq,k)=(vd(iq,k)-.5*vd(isv(iq),k)
     &                 -va(inv(iq),k)/10. +va(issv(iq),k)/4.)/.95
               enddo
            enddo

            call boundsuv(uin,vin,nrows=2)
!cdir nodep
            do k=1,kl
               do iq=1,ifull
                  ua(iq,k)=(ud(iq,k)-.5*ud(iwu(iq),k)
     &                 -uin(ieu(iq),k)/10. +uin(iwwu(iq),k)/4.)/.95
                  va(iq,k)=(vd(iq,k)-.5*vd(isv(iq),k)
     &                 -vin(inv(iq),k)/10. +vin(issv(iq),k)/4.)/.95
               enddo
            end do
         end do                 ! itn=1,itnmax
	 
      else
         print*, "Error, unsupported nstag option:", nstag
         stop
      end if

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
      real, dimension(ifull,kl), intent(in)  :: u, v
      real, dimension(ifull,kl), intent(out) :: uout, vout
      real, dimension(ifull+iextra,kl) :: ua, va, ud, vd,
     &                                    uin, vin
      integer, parameter :: itnmax=3
      integer :: iq, itn, k, num
      save num
      data num/0/

      call start_log(stag_begin)
      if(nmaxpr==1.and.mydiag)then
        print *,'unstag_ktau,nstag,nstagu',ktau,nstag,nstagu
      endif
      do k=1,kl
         do iq=1,ifull
            uin(iq,k) = u(iq,k)
            vin(iq,k) = v(iq,k)
         end do
      end do
      
      if(abs(nstagu)<3)then
        call boundsuv(uin,vin,nrows=2)
        do k=1,kl
         do iq=1,ifull
          uout(iq,k)=(9.*(uin(iwu(iq),k)+uin(iq,k))
     &                   -uin(iwwu(iq),k)-uin(ieu(iq),k))/16.
          vout(iq,k)=(9.*(vin(isv(iq),k)+vin(iq,k))
     &                   -vin(issv(iq),k)-vin(inv(iq),k))/16.
           uout(iq,k)=.5*(uin(iwu(iq),k)+uin(iq,k))  ! for linear tests
           vout(iq,k)=.5*(vin(isv(iq),k)+vin(iq,k))  ! for linear tests
         enddo   ! iq loop
        enddo
        return
      endif  ! (nstagu==0)

      if ( nstagu==3 ) then
         call boundsuv(uin,vin,nrows=2)
         do k=1,kl
!cdir nodep
            do iq=1,ifull       ! precalculate rhs terms with iwwu2 & issv2
               ud(iq,k)=uin(iq,k)/2.+uin(iwu(iq),k)+uin(iwwu(iq),k)/10.
               vd(iq,k)=vin(iq,k)/2.+vin(isv(iq),k)+vin(issv(iq),k)/10.
            end do
         end do

         call boundsuv(ud,vd)
         do k=1,kl
!cdir nodep
            do iq=1,ifull
               ua(iq,k)=ud(iq,k)-ud(iwu(iq),k)/2. ! 1st guess
               va(iq,k)=vd(iq,k)-vd(isv(iq),k)/2. ! 1st guess
            end do
         end do

         do itn=1,itnmax        ! each loop is a double iteration

            call boundsuv(ua,va,nrows=2)
            do k=1,kl
!cdir nodep
               do iq=1,ifull
                  uin(iq,k)=(ud(iq,k)-.5*ud(iwu(iq),k)
     &                 -ua(ieu(iq),k)/10. +ua(iwwu(iq),k)/4.)/.95
                  vin(iq,k)=(vd(iq,k)-.5*vd(isv(iq),k)
     &                 -va(inv(iq),k)/10. +va(issv(iq),k)/4.)/.95
               end do
            end do

            call boundsuv(uin,vin,nrows=2)
            do k=1,kl
!cdir nodep
               do iq=1,ifull
                  ua(iq,k)=(ud(iq,k)-.5*ud(iwu(iq),k)
     &                 -uin(ieu(iq),k)/10. +uin(iwwu(iq),k)/4.)/.95
                  va(iq,k)=(vd(iq,k)-.5*vd(isv(iq),k)
     &                 -vin(inv(iq),k)/10. +vin(issv(iq),k)/4.)/.95
               end do
            end do

         end do                 ! itn=1,itnmax

      else if ( nstagu==4 ) then
         call boundsuv(uin,vin)
         do k=1,kl
!cdir nodep
            do iq=1,ifull       ! precalculate rhs terms
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

         do itn=1,itnmax        ! each loop is a double iteration
            call boundsuv(ua,va,nrows=2)
            do k=1,kl
!cdir nodep
               do iq=1,ifull
                  uin(iq,k)=(ud(iq,k)-.5*ud(ieu(iq),k)
     &                 -ua(iwu(iq),k)/10. +ua(ieeu(iq),k)/4.)/.95
                  vin(iq,k)=(vd(iq,k)-.5*vd(inv(iq),k)
     &                 -va(isv(iq),k)/10. +va(innv(iq),k)/4.)/.95
               enddo
            enddo
            call boundsuv(uin,vin,nrows=2)
            do k=1,kl
!cdir nodep
               do iq=1,ifull
                  ua(iq,k)=(ud(iq,k)-.5*ud(ieu(iq),k)
     &                 -uin(iwu(iq),k)/10. +uin(ieeu(iq),k)/4.)/.95
                  va(iq,k)=(vd(iq,k)-.5*vd(inv(iq),k)
     &                 -vin(isv(iq),k)/10. +vin(innv(iq),k)/4.)/.95
               enddo
            enddo
         enddo                  ! itn=1,itnmax
      
      else
         print*, "Error, unsupported nstagu option:", nstagu
         stop
      end if

      do k=1,kl
         do iq=1,ifull          ! final values for output
            uout(iq,k) = ua(iq,k)
            vout(iq,k) = va(iq,k)
         end do
      end do
      call end_log(stag_end)
      return
      end        
