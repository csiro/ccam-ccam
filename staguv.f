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
      include 'vecsuva.h'  ! vecsuva info
      real, dimension(ifull,kl), intent(in)  :: u, v
      real, dimension(ifull,kl), intent(out) :: uout, vout
      real, dimension(ifull+iextra,kl) :: ua, va, ud, vd,
     &                                    uin, vin
      common/ktau_stag/ktau_stag
      integer, parameter :: itnmax=3
      integer :: iq, itn, k, nstagin, num, ktau_stag
      save nstagin,num
      data num/0/

      call start_log(stag_begin)
      if(num==0)then
        num=1
	 if(nstag==5)then   ! to be backward compatible with pre-Oct '04
	   nstag=-1
	   nstagu=4
	 endif
        nstagin=nstag
	 nstag=nstagu	 
      endif
      ktau_stag=ktau
!     N.B. swapping only done in unstaguv, during calls from nonlin      
c     print *,'ktau,nstag,nstagu,mod,mod2 ',
c    .         ktau,nstag,nstagu,mod(ktau,abs(nstagin)),mod(ktau,2)

      ! Copying could be avoided if input arrays were dimensioned ifull+iextra
      do k=1,kl
         do iq=1,ifull
            uin(iq,k) = u(iq,k)
            vin(iq,k) = v(iq,k)
         end do
      end do

c     if ( nstag==3 .or. (nstagin==5 .and. modulo(ktau,2)==0) ) then
      if ( nstag==3 ) then
c         print *,'doing nstag3 for ktau = ',ktau

         call boundsuv(uin,vin,nrows=2)
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

c     else if ( nstag==4 .or. (nstagin==5 .and. modulo(ktau,2)==1)) then
      else if ( nstag==4 ) then
c         print *,'doing nstag4 for ktau = ',ktau

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
      include 'vecsuva.h'  ! vecsuva info
      real, dimension(ifull,kl), intent(in)  :: u, v
      real, dimension(ifull,kl), intent(out) :: uout, vout
      real, dimension(ifull+iextra,kl) :: ua, va, ud, vd,
     &                                    uin, vin
      common/ktau_stag/ktau_stag
      integer, parameter :: itnmax=3
      integer :: iq, itn, k, nstagin, num, ktau_stag
      save nstagin,num
      data num/0/

      call start_log(stag_begin)
      if(num==0)then
        num=1
	 if(nstag==5)then   ! to be backward compatible with pre-Oct '04
	   nstag=-1
	   nstagu=4
	 endif
	 if(nstag==6)then  
	   nstagu=4
	 endif
        nstagin=nstag
	 nstag=nstagu
	 ktau_stag=0   ! only set to ktau in call to staguv
      endif
!     N.B. swapping only done in unstaguv, during calls from nonlin      
      if(num<ktau)then  ! following only for very first time each ktau
        if(nstagin<0.and.mod(ktau,abs(nstagin))==0)then
          nstag=7-nstagu   ! swap between 3 & 4
	   nstagu=nstag
          num=ktau
        endif
        if(nstagin==6.and.ktau==ktau_stag)then
!         this swapping only done in unstaguv, during calls from adjust5     
          nstag=7-nstagu   ! swap between 3 & 4
	   nstagu=nstag
          num=ktau
        endif
      endif  !  (num<ktau)
      if(diag.and.mydiag)then
        print *,'uns ktau,nstag,nstagu,mod,mod2 ',
     &           ktau,nstag,nstagu,mod(ktau,abs(nstagin)),mod(ktau,2)
        print *,'nstagin,ktau_stag ',nstagin,ktau_stag
      endif
      do k=1,kl
         do iq=1,ifull
            uin(iq,k) = u(iq,k)
            vin(iq,k) = v(iq,k)
         end do
      end do

c     if ( nstagu==3 .or. (nstagu==5 .and. modulo(ktau,2)==0)) then
      if ( nstagu==3 ) then
c         print *,'doing unstag3 for ktau = ',ktau

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

c     else if ( nstagu==4 .or. (nstagu==5 .and. modulo(ktau,2)==1)) then
      else if ( nstagu==4 ) then
c         print *,'doing unstag4 for ktau = ',ktau

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
