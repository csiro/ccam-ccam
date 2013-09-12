      subroutine ints(ntr,s,intsch,nface,xg,yg,nfield)
c     parameter (mhint):   0 for simple; 2 for Bessel (was called mbess)
!     jlm finds Bessel (mhint=2) gives bigger overshoots near discontinuities
c     this one includes Bermejo & Staniforth option
c     parameter (mh_bs): B&S on/off depending on value of nfield
c     can put sx into work array (only 2d)
c     later may wish to save idel etc between array calls
c     this one does linear interp in x on outer y sides
c     doing x-interpolation before y-interpolation
!     nfield: 1 (psl), 2 (u, v), 3 (T), 4 (gases)
      use cc_mpi
      use indices_m
      implicit none
      integer, parameter :: ntest=0
      include 'newmpar.h'
      include 'parm.h'
      include 'parmhor.h'    ! has mh_bs
      integer, intent(in) :: ntr
      integer intsch, nfield
      integer nface(ifull,kl)
      real xg(ifull,kl),yg(ifull,kl)      ! now passed through call
      ! Need work common for this
      real sx(ntr,-1:ipan+2,-1:jpan+2,1:npan,kl)
      real s(ifull+iextra,kl,ntr)
      integer idel, iq, jdel, nn
      real xxg, yyg
      real, dimension(ntr) :: c1, c2, c3, c4
      real, dimension(ntr) :: a3, a4, sss
      real, dimension(ntr) :: cmax, cmin
      real, dimension(ntr,4) :: r
      integer i, j, k, n, ind, ip, jp, iproc
      integer ii
      ! This is really indp, just repeated here to get inlining to work
      ind(i,j,n)=i+(j-1)*ipan+(n-1)*ipan*jpan  ! *** for n=1,npan

      call start_log(ints_begin)
      call bounds(s,nrows=2)

!======================== start of intsch=1 section ====================
      if(intsch==1)then

        do nn=1,ntr
          do k=1,kl              ! A single main k loop uses cache better

            do n=1,npan         ! first simple copy into larger array
              do j=1,jpan
                do i=1,ipan
                  sx(nn,i,j,n,k) = s(ind(i,j,n),k,nn)
                enddo         ! i loop
              enddo           ! j loop
            
c       this is intsb           EW interps done first
c       first extend s arrays into sx - this one -1:il+2 & -1:il+2
              do j=1,jpan
                sx(nn,0,j,n,k)      = s(iw(ind(1,j,n)),k,nn)
                sx(nn,-1,j,n,k)     = s(iww(ind(1,j,n)),k,nn)
                sx(nn,ipan+1,j,n,k) = s(ie(ind(ipan,j,n)),k,nn)
                sx(nn,ipan+2,j,n,k) = s(iee(ind(ipan,j,n)),k,nn)
              enddo            ! j loop
              do i=1,ipan
                sx(nn,i,0,n,k)      = s(is(ind(i,1,n)),k,nn)
                sx(nn,i,-1,n,k)     = s(iss(ind(i,1,n)),k,nn)
                sx(nn,i,jpan+1,n,k) = s(in(ind(i,jpan,n)),k,nn)
                sx(nn,i,jpan+2,n,k) = s(inn(ind(i,jpan,n)),k,nn)
              enddo            ! i loop
c        for ew interpolation, sometimes need (different from ns):
c            (-1,0),   (0,0),   (0,-1)   (-1,il+1),   (0,il+1),   (0,il+2)
c           (il+1,0),(il+2,0),(il+1,-1) (il+1,il+1),(il+2,il+1),(il+1,il+2)
!        sx(nn,-1,0,n)=sx(nn,0,2,n)            !  s(iwws6(1,1,n))
!        sx(nn,0,0,n)=sx(nn,0,1,n)             !  s(iws6(1,1,n))
!        sx(nn,0,-1,n)=sx(nn,-1,1,n)           !  s(iwss6(1,1,n)) or s(isws6(1,1,n))
!        sx(nn,il+1,0,n)=sx(nn,il+1,1,n)       !  s(ies6(il,1,n))
!        sx(nn,il+2,0,n)=sx(nn,il+1,2,n)       !  s(iees6(il,1,n))
!        sx(nn,il+1,-1,n)=sx(nn,il+2,1,n)      !  s(iess6(il,1,n))
!        sx(nn,-1,il+1,n)=sx(nn,0,il-1,n)      !  s(iwwn6(1,il,n))
!        sx(nn,0,il+2,n)=sx(nn,-1,il,n)        !  s(iwnn6(1,il,n))
!        sx(nn,il+2,il+1,n)=sx(nn,il+1,il-1,n) !  s(ieen6(il,il,n))
!        sx(nn,il+1,il+2,n)=sx(nn,il+2,il,n)   !  s(ienn6(il,il,n))
!        sx(nn,0,il+1,n)=sx(nn,0,il,n)         !  s(iwn6(1,il,n))   **
!        sx(nn,il+1,il+1,n)=sx(nn,il+1,il,n)   !  s(ien6(il,il,n))  **

              sx(nn,-1,0,n,k)          = s(lwws(n),k,nn)
              sx(nn,0,0,n,k)           = s(iws(ind(1,1,n)),k,nn)
              sx(nn,0,-1,n,k)          = s(lwss(n),k,nn)
              sx(nn,ipan+1,0,n,k)      = s(ies(ind(ipan,1,n)),k,nn)
              sx(nn,ipan+2,0,n,k)      = s(lees(n),k,nn)
              sx(nn,ipan+1,-1,n,k)     = s(less(n),k,nn)
              sx(nn,-1,jpan+1,n,k)     = s(lwwn(n),k,nn)
              sx(nn,0,jpan+2,n,k)      = s(lwnn(n),k,nn)
              sx(nn,ipan+2,jpan+1,n,k) = s(leen(n),k,nn)
              sx(nn,ipan+1,jpan+2,n,k) = s(lenn(n),k,nn)
              sx(nn,0,jpan+1,n,k)      = s(iwn(ind(1,jpan,n)),k,nn)
              sx(nn,ipan+1,jpan+1,n,k) = s(ien(ind(ipan,jpan,n)),k,nn)
            enddo               ! n loop
          end do                ! k loop
        end do                  ! nn loop

! Loop over points that need to be calculated for other processes
        if(nfield<mh_bs)then
          do ii=1,neighnum
            iproc=neighlistsend(ii)
            do iq=1,drlen(iproc)
              !  Convert face index from 0:npanels to array indices
              ip = min(il_g,max(1,nint(dpoints(iproc)%a(2,iq))))
              jp = min(il_g,max(1,nint(dpoints(iproc)%a(3,iq))))
              n = nint(dpoints(iproc)%a(1,iq)) + noff ! Local index
              !  Need global face index in fproc call
#ifdef debug
              if ( fproc(ip,jp,n-noff) /= myid ) then
                print*, "Error in ints", myid, n, iq, iproc,
     &                  dpoints(iproc)%a(:,iq)
                call ccmpi_abort(-1)
              end if
#endif
              idel = int(dpoints(iproc)%a(2,iq))
              xxg = dpoints(iproc)%a(2,iq) - idel
              jdel = int(dpoints(iproc)%a(3,iq))
              yyg = dpoints(iproc)%a(3,iq) - jdel
              k = nint(dpoints(iproc)%a(4,iq))
              idel = idel - ioff(n-noff)
              jdel = jdel - joff(n-noff)
              
              
              c1 = sx(:,idel-1,jdel,n,k) ! manually unrolled loop
              c2 = sx(:,idel  ,jdel,n,k)
              c3 = sx(:,idel+1,jdel,n,k)
              c4 = sx(:,idel+2,jdel,n,k)
              if(mhint==2)then ! Bessel interp
                a4 = c4-c1+3.*(c2-c3)
                a3 = c1-2.*c2+c3-a4
                r(:,2) = c2+.5*xxg*(c3-c1 +xxg*(a3+xxg*a4))
              else
                r(:,2) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*c2
     &                    -xxg*c1/3.)
     &                    -xxg*(1.+xxg)*c4/3.)
     &                    +xxg*(1.+xxg)*(2.-xxg)*c3)/2.
              endif         !  (mhint==2)
              c1 = sx(:,idel-1,jdel+1,n,k)
              c2 = sx(:,idel  ,jdel+1,n,k)
              c3 = sx(:,idel+1,jdel+1,n,k)
              c4 = sx(:,idel+2,jdel+1,n,k)
              if(mhint==2)then ! Bessel interp
                a4 = c4-c1+3.*(c2-c3)
                a3 = c1-2.*c2+c3-a4
                r(:,3) = c2+.5*xxg*(c3-c1 +xxg*(a3+xxg*a4))
              else
                r(:,3) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*c2
     &                    -xxg*c1/3.)
     &                    -xxg*(1.+xxg)*c4/3.)
     &                    +xxg*(1.+xxg)*(2.-xxg)*c3)/2.
              endif         !  (mhint==2)
c           r = {(1-x)*{(2-x     )*[(1+x     )*c2-x     *c1/3]
c                -x*(1+x)*c4/3}
c                +x*(1+x)*(2-x)*c3}/2
              do nn=1,4,3   ! N.B.
                c2 = sx(:,idel  ,jdel+nn-2,n,k)
                c3 = sx(:,idel+1,jdel+nn-2,n,k)
                r(:,nn) = (1.-xxg)*c2 +xxg*c3
              enddo         ! nn loop
              if(mhint==2)then ! Bessel interp
                a4 = r(:,4)-r(:,1)+3.*(r(:,2)-r(:,3))
                a3 = r(:,1)-2.*r(:,2)+r(:,3)-a4
                do nn=1,ntr
                  sextra(iproc)%a(nn+(iq-1)*ntr) = r(nn,2) +
     &                    0.5*yyg*(r(nn,3)-r(nn,1)
     &                    +yyg*(a3(nn)+yyg*a4(nn)))
                end do
              else
                do nn=1,ntr
                  sextra(iproc)%a(nn+(iq-1)*ntr) =
     &                    ((1.-yyg)*((2.-yyg)*
     &                    ((1.+yyg)*r(nn,2)-yyg*r(nn,1)/3.)
     &                    -yyg*(1.+yyg)*r(nn,4)/3.)
     &                    +yyg*(1.+yyg)*(2.-yyg)*r(nn,3))/2.
                end do
              endif         !  (mhint==2)
            enddo            ! iq loop
          end do              ! iproc loop
        else                   ! (nfield<mh_bs)
          do ii=1,neighnum
            iproc=neighlistsend(ii)
            do iq=1,drlen(iproc)
              !  Convert face index from 0:npanels to array indices
              ip = min(il_g,max(1,nint(dpoints(iproc)%a(2,iq))))
              jp = min(il_g,max(1,nint(dpoints(iproc)%a(3,iq))))
              n = nint(dpoints(iproc)%a(1,iq)) + noff ! Local index
              !  Need global face index in fproc call
#ifdef debug
              if ( fproc(ip,jp,n-noff) /= myid ) then
                print*, "Error in ints", myid, n, iq, iproc,
     &                   dpoints(iproc)%a(:,iq)
                call ccmpi_abort(-1)
              end if
#endif
              idel = int(dpoints(iproc)%a(2,iq))
              xxg = dpoints(iproc)%a(2,iq) - idel
              jdel = int(dpoints(iproc)%a(3,iq))
              yyg = dpoints(iproc)%a(3,iq) - jdel
              k = nint(dpoints(iproc)%a(4,iq))
              idel = idel - ioff(n-noff)
              jdel = jdel - joff(n-noff)
              c1 = sx(:,idel-1,jdel,n,k) ! manually unrolled loop
              c2 = sx(:,idel  ,jdel,n,k)
              c3 = sx(:,idel+1,jdel,n,k)
              c4 = sx(:,idel+2,jdel,n,k)
              cmin = min( 1.e20,c2,c3) ! Bermejo & Staniforth
              cmax = max(-1.e20,c2,c3) ! Bermejo & Staniforth
              if(mhint==2)then ! Bessel interp
                a4 = c4-c1+3.*(c2-c3)
                a3 = c1-2.*c2+c3-a4
                r(:,2) = c2+.5*xxg*(c3-c1 +xxg*(a3+xxg*a4))
              else
                r(:,2) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*c2
     &                    -xxg*c1/3.)
     &                    -xxg*(1.+xxg)*c4/3.)
     &                    +xxg*(1.+xxg)*(2.-xxg)*c3)/2.
              endif         !  (mhint==2)
              c1 = sx(:,idel-1,jdel+1,n,k)
              c2 = sx(:,idel  ,jdel+1,n,k)
              c3 = sx(:,idel+1,jdel+1,n,k)
              c4 = sx(:,idel+2,jdel+1,n,k)
              cmin = min(cmin,c2,c3) ! Bermejo & Staniforth
              cmax = max(cmax,c2,c3) ! Bermejo & Staniforth
              if(mhint==2)then ! Bessel interp
                a4 = c4-c1+3.*(c2-c3)
                a3 = c1-2.*c2+c3-a4
                r(:,3) = c2+.5*xxg*(c3-c1 +xxg*(a3+xxg*a4))
              else
                r(:,3) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*c2
     &                    -xxg*c1/3.)
     &                    -xxg*(1.+xxg)*c4/3.)
     &                    +xxg*(1.+xxg)*(2.-xxg)*c3)/2.
              endif         !  (mhint==2)
c           r = {(1-x)*{(2-x     )*[(1+x     )*c2-x     *c1/3]
c                -x*(1+x)*c4/3}
c                +x*(1+x)*(2-x)*c3}/2
              do nn=1,4,3   ! N.B.
                c2 = sx(:,idel  ,jdel+nn-2,n,k)
                c3 = sx(:,idel+1,jdel+nn-2,n,k)
                r(:,nn) = (1.-xxg)*c2 +xxg*c3
              enddo         ! nn loop
              if(mhint==2)then ! Bessel interp
                a4 = r(:,4)-r(:,1)+3.*(r(:,2)-r(:,3))
                a3 = r(:,1)-2.*r(:,2)+r(:,3)-a4
                sss = r(:,2)+.5*yyg*(r(:,3)-r(:,1)
     &                +yyg*(a3+yyg*a4))
              else
                sss = ((1.-yyg)*((2.-yyg)*
     &                ((1.+yyg)*r(:,2)-yyg*r(:,1)/3.)
     &                -yyg*(1.+yyg)*r(:,4)/3.)
     &                +yyg*(1.+yyg)*(2.-yyg)*r(:,3))/2.
              endif         !  (mhint==2)
              do nn=1,ntr
                sextra(iproc)%a(nn+(iq-1)*ntr) = 
     &             min(max(cmin(nn),sss(nn)),cmax(nn)) ! Bermejo & Staniforth
              end do
            enddo            ! iq loop
          end do              ! iproc loop
        endif                  ! (nfield<mh_bs)  .. else ..

        call intssync_send(ntr)

        if(nfield<mh_bs)then
          do k=1,kl
            do iq=1,ifull    ! non Berm-Stan option
!             Convert face index from 0:npanels to array indices
              idel=int(xg(iq,k))
              xxg=xg(iq,k)-idel
              jdel=int(yg(iq,k))
              yyg=yg(iq,k)-jdel
              ! Now make them proper indices in this processor's region
              idel = idel - ioff(nface(iq,k))
              jdel = jdel - joff(nface(iq,k))
              n = nface(iq,k) + noff ! Make this a local index

              if ( idel < 0 .or. idel > ipan .or. jdel < 0 .or.
     &             jdel > jpan .or. n < 1 .or. n > npan ) then
                 cycle      ! Will be calculated on another processor
              end if

              c1 = sx(:,idel-1,jdel,n,k) ! manually unrolled loop
              c2 = sx(:,idel  ,jdel,n,k)
              c3 = sx(:,idel+1,jdel,n,k)
              c4 = sx(:,idel+2,jdel,n,k)
              if(mhint==2)then ! Bessel interp
                a4 = c4-c1+3.*(c2-c3)
                a3 = c1-2.*c2+c3-a4
                r(:,2) = c2+.5*xxg*(c3-c1 +xxg*(a3+xxg*a4))
              else
                r(:,2) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*c2
     &                    -xxg*c1/3.)
     &                    -xxg*(1.+xxg)*c4/3.)
     &                    +xxg*(1.+xxg)*(2.-xxg)*c3)/2.
              endif         !  (mhint==2)
              c1 = sx(:,idel-1,jdel+1,n,k)
              c2 = sx(:,idel  ,jdel+1,n,k)
              c3 = sx(:,idel+1,jdel+1,n,k)
              c4 = sx(:,idel+2,jdel+1,n,k)
              if(mhint==2)then ! Bessel interp
                a4 = c4-c1+3.*(c2-c3)
                a3 = c1-2.*c2+c3-a4
                r(:,3) = c2+.5*xxg*(c3-c1 +xxg*(a3+xxg*a4))
              else
                r(:,3) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*c2
     &                    -xxg*c1/3.)
     &                    -xxg*(1.+xxg)*c4/3.)
     &                    +xxg*(1.+xxg)*(2.-xxg)*c3)/2.
              endif         !  (mhint==2)
c           r = {(1-x)*{(2-x     )*[(1+x     )*c2-x     *c1/3]
c                -x*(1+x)*c4/3}
c                +x*(1+x)*(2-x)*c3}/2
              do nn=1,4,3   ! N.B.
                c2 = sx(:,idel  ,jdel+nn-2,n,k)
                c3 = sx(:,idel+1,jdel+nn-2,n,k)
                r(:,nn) = (1.-xxg)*c2 +xxg*c3
              enddo         ! nn loop
              if(mhint==2)then ! Bessel interp
                a4 = r(:,4)-r(:,1)+3.*(r(:,2)-r(:,3))
                a3 = r(:,1)-2.*r(:,2)+r(:,3)-a4
                s(iq,k,:) = r(:,2)+.5*yyg*(r(:,3)-r(:,1)
     &                     +yyg*(a3+yyg*a4))
              else
                s(iq,k,:) = ((1.-yyg)*((2.-yyg)*
     &                    ((1.+yyg)*r(:,2)-yyg*r(:,1)/3.)
     &                    -yyg*(1.+yyg)*r(:,4)/3.)
     &                    +yyg*(1.+yyg)*(2.-yyg)*r(:,3))/2.
              endif         !  (mhint==2)
            enddo            ! iq loop
          enddo               ! k loop
        else                ! (nfield<mh_bs)
          do k=1,kl
            do iq=1,ifull    ! Berm-Stan option here e.g. qg & gases
              idel=int(xg(iq,k))
              xxg=xg(iq,k)-idel
              jdel=int(yg(iq,k))
              yyg=yg(iq,k)-jdel
              ! Now make them proper indices in this processor's region
              idel = idel - ioff(nface(iq,k))
              jdel = jdel - joff(nface(iq,k))
              n = nface(iq,k) + noff ! Make this a local index

              if ( idel < 0 .or. idel > ipan .or. jdel < 0 .or.
     &             jdel > jpan .or. n < 1 .or. n > npan ) then
                cycle      ! Will be calculated on another processor
              end if

              c1 = sx(:,idel-1,jdel,n,k) ! manually unrolled loop
              c2 = sx(:,idel  ,jdel,n,k)
              c3 = sx(:,idel+1,jdel,n,k)
              c4 = sx(:,idel+2,jdel,n,k)
              cmin = min( 1.e20,c2,c3) ! Bermejo & Staniforth
              cmax = max(-1.e20,c2,c3) ! Bermejo & Staniforth
              if(mhint==2)then ! Bessel interp
                a4 = c4-c1+3.*(c2-c3)
                a3 = c1-2.*c2+c3-a4
                r(:,2) = c2+.5*xxg*(c3-c1 +xxg*(a3+xxg*a4))
              else
                r(:,2) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*c2
     &                    -xxg*c1/3.)
     &                    -xxg*(1.+xxg)*c4/3.)
     &                    +xxg*(1.+xxg)*(2.-xxg)*c3)/2.
              endif         !  (mhint==2)
              c1 = sx(:,idel-1,jdel+1,n,k)
              c2 = sx(:,idel  ,jdel+1,n,k)
              c3 = sx(:,idel+1,jdel+1,n,k)
              c4 = sx(:,idel+2,jdel+1,n,k)
              cmin = min(cmin,c2,c3) ! Bermejo & Staniforth
              cmax = max(cmax,c2,c3) ! Bermejo & Staniforth
              if(mhint==2)then ! Bessel interp
                a4 = c4-c1+3.*(c2-c3)
                a3 = c1-2.*c2+c3-a4
                r(:,3) = c2+.5*xxg*(c3-c1 +xxg*(a3+xxg*a4))
              else
                r(:,3) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*c2
     &                    -xxg*c1/3.)
     &                    -xxg*(1.+xxg)*c4/3.)
     &                    +xxg*(1.+xxg)*(2.-xxg)*c3)/2.
              endif         !  (mhint==2)
c           r = {(1-x)*{(2-x     )*[(1+x     )*c2-x     *c1/3]
c                -x*(1+x)*c4/3}
c                +x*(1+x)*(2-x)*c3}/2
              do nn=1,4,3   ! N.B.
                c2 = sx(:,idel  ,jdel+nn-2,n,k)
                c3 = sx(:,idel+1,jdel+nn-2,n,k)
                r(:,nn) = (1.-xxg)*c2 +xxg*c3
              enddo         ! nn loop
              if(mhint==2)then ! Bessel interp
                a4 = r(:,4)-r(:,1)+3.*(r(:,2)-r(:,3))
                a3 = r(:,1)-2.*r(:,2)+r(:,3)-a4
                sss = r(:,2)+.5*yyg*(r(:,3)-r(:,1)
     &                +yyg*(a3+yyg*a4))
              else
                sss = ((1.-yyg)*((2.-yyg)*
     &                ((1.+yyg)*r(:,2)-yyg*r(:,1)/3.)
     &                -yyg*(1.+yyg)*r(:,4)/3.)
     &                +yyg*(1.+yyg)*(2.-yyg)*r(:,3))/2.
              endif         !  (mhint==2)
              s(iq,k,:) = min(max(cmin,sss),cmax) ! Bermejo & Staniforth
            enddo            ! iq loop
          enddo               ! k loop
        endif               ! (nfield<mh_bs)  .. else ..
            
!========================   end of intsch=1 section ====================
      else     ! if(intsch==1)then
!======================== start of intsch=2 section ====================
c       this is intsc           NS interps done first
c       first extend s arrays into sx - this one -1:il+2 & -1:il+2
        do nn=1,ntr
          do k=1,kl              ! A single main k loop uses cache better
            do n=1,npan         ! first simple copy into larger array
              do j=1,jpan
                do i=1,ipan
                  sx(nn,i,j,n,k)=s(ind(i,j,n),k,nn)
                enddo         ! i loop
              enddo           ! j loop

              do j=1,jpan
                sx(nn,0,j,n,k) = s(iw(ind(1,j,n)),k,nn)
                sx(nn,-1,j,n,k) = s(iww(ind(1,j,n)),k,nn)
                sx(nn,ipan+1,j,n,k) = s(ie(ind(ipan,j,n)),k,nn)
                sx(nn,ipan+2,j,n,k) = s(iee(ind(ipan,j,n)),k,nn)
              enddo            ! j loop
              do i=1,ipan
                sx(nn,i,0,n,k) = s(is(ind(i,1,n)),k,nn)
                sx(nn,i,-1,n,k) = s(iss(ind(i,1,n)),k,nn)
                sx(nn,i,jpan+1,n,k) = s(in(ind(i,jpan,n)),k,nn)
                sx(nn,i,jpan+2,n,k) = s(inn(ind(i,jpan,n)),k,nn)
              enddo            ! i loop
c        for ns interpolation, sometimes need (different from ew):
c            (-1,0),   (0,0),   (0,-1)   (-1,il+1),   (0,il+1),   (0,il+2)
c          (il+1,0),(il+2,0),(il+1,-1) (il+1,il+1),(il+2,il+1),(il+1,il+2)
!        sx(nn,-1,0,n,k)=sx(nn,1,-1,n,k)           !  s(isww6(1,1,n))
!        sx(nn,0,0,n,k)=sx(nn,1,0,n,k)             !  s(isw6(1,1,n))
!        sx(nn,0,-1,n,k)=sx(nn,2,0,n,k)            !  s(issw6(1,1,n))
!        sx(nn,il+2,0,n,k)=sx(nn,il,-1,n,k)        !  s(isee6(il,1,n))
!        sx(nn,il+1,-1,n,k)=sx(nn,il-1,0,n,k)      !  s(isse6(il,1,n))
!        sx(nn,-1,il+1,n,k)=sx(nn,1,il+2,n,k)      !  s(inww6(1,il,n))
!        sx(nn,0,il+1,n,k)=sx(nn,1,il+1,n,k)       !  s(inw6(1,il,n))
!        sx(nn,0,il+2,n,k)=sx(nn,2,il+1,n,k)       !  s(innw6(1,il,n))
!        sx(nn,il+2,il+1,n,k)=sx(nn,il,il+2,n,k)   !  s(inee6(il,il,n))
!        sx(nn,il+1,il+2,n,k)=sx(nn,il-1,il+1,n,k) !  s(inne6(il,il,n))
!        sx(nn,il+1,0,n,k)=sx(nn,il,0,n,k)         !  s(ise6(il,1,n))    **
!        sx(nn,il+1,il+1,n,k)=sx(nn,il,il+1,n,k)   !  s(ine6(il,il,n))   **

              sx(nn,-1,0,n,k)=s(lsww(n),k,nn)
              sx(nn,0,0,n,k) = s(isw(ind(1,1,n)),k,nn)
              sx(nn,0,-1,n,k) = s(lssw(n),k,nn)
              sx(nn,ipan+2,0,n,k) = s(lsee(n),k,nn)
              sx(nn,ipan+1,-1,n,k) = s(lsse(n),k,nn)
              sx(nn,-1,jpan+1,n,k) = s(lnww(n),k,nn)
              sx(nn,0,jpan+1,n,k) = s(inw(ind(1,jpan,n)),k,nn)
              sx(nn,0,jpan+2,n,k) = s(lnnw(n),k,nn)
              sx(nn,ipan+2,jpan+1,n,k) = s(lnee(n),k,nn)
              sx(nn,ipan+1,jpan+2,n,k) = s(lnne(n),k,nn)
              sx(nn,ipan+1,0,n,k)    = s(ise(ind(ipan,1,n)),k,nn)
              sx(nn,ipan+1,jpan+1,n,k) = s(ine(ind(ipan,jpan,n)),k,nn)
            enddo               ! n loop
          end do                ! k loop
        end do                  ! nn loop

!       For other processes
        if(nfield<mh_bs)then
          do ii=1,neighnum
            iproc=neighlistsend(ii)
            do iq=1,drlen(iproc)
              !  Convert face index from 0:npanels to array indices
              ip = min(il_g,max(1,nint(dpoints(iproc)%a(2,iq))))
              jp = min(il_g,max(1,nint(dpoints(iproc)%a(3,iq))))
              n = nint(dpoints(iproc)%a(1,iq)) + noff ! Local index
              !  Need global face index in fproc call
#ifdef debug
              if ( fproc(ip,jp,n-noff) /= myid ) then
                 print*, "Error in ints_bl", myid, n, iq, iproc,
     &             dpoints(iproc)%a(:,iq)
                 call ccmpi_abort(-1)
              end if
#endif
              idel = int(dpoints(iproc)%a(2,iq))
              xxg = dpoints(iproc)%a(2,iq) - idel
              jdel = int(dpoints(iproc)%a(3,iq))
              yyg = dpoints(iproc)%a(3,iq) - jdel
              k = nint(dpoints(iproc)%a(4,iq))
              idel = idel - ioff(n-noff)
              jdel = jdel - joff(n-noff)
              c1 = sx(:,idel,jdel-1,n,k) ! manually unrolled loop
              c2 = sx(:,idel,jdel  ,n,k)
              c3 = sx(:,idel,jdel+1,n,k)
              c4 = sx(:,idel,jdel+2,n,k)
              if(mhint==2)then ! Bessel interp
                a4 = c4-c1+3.*(c2-c3)
                a3 = c1-2.*c2+c3-a4
                r(:,2) = c2+.5*yyg*(c3-c1 +yyg*(a3+yyg*a4))
              else
                r(:,2) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*c2
     &                    -yyg*c1/3.)
     &                    -yyg*(1.+yyg)*c4/3.)
     &                    +yyg*(1.+yyg)*(2.-yyg)*c3)/2.
              endif         !  (mhint==2)
              c1 = sx(:,idel+1,jdel-1,n,k)
              c2 = sx(:,idel+1,jdel  ,n,k)
              c3 = sx(:,idel+1,jdel+1,n,k)
              c4 = sx(:,idel+1,jdel+2,n,k)
              if(mhint==2)then ! Bessel interp
                a4 = c4-c1+3.*(c2-c3)
                a3 = c1-2.*c2+c3-a4
                r(:,3) = c2+.5*yyg*(c3-c1 +yyg*(a3+yyg*a4))
              else
                r(:,3) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*c2
     &                    -yyg*c1/3.)
     &                    -yyg*(1.+yyg)*c4/3.)
     &                    +yyg*(1.+yyg)*(2.-yyg)*c3)/2.
              endif         !  (mhint==2)
c          r = {(1-y)*{(2-y     )*[(1+y     )*c2-y     *c1/3]
c               -y*(1+y)*c4/3}
c               +y*(1+y)*(2-y)*c3}/2
              do nn=1,4,3   ! N.B.
                c2 = sx(:,idel+nn-2,jdel  ,n,k)
                c3 = sx(:,idel+nn-2,jdel+1,n,k)
                r(:,nn) = (1.-yyg)*c2 +yyg*c3
              enddo         ! nn loop
              if(mhint==2)then ! Bessel interp
                a4 = r(:,4)-r(:,1)+3.*(r(:,2)-r(:,3))
                a3 = r(:,1)-2.*r(:,2)+r(:,3)-a4
                do nn=1,ntr
                  sextra(iproc)%a(nn+(iq-1)*ntr) = r(nn,2)+
     &                   0.5*xxg*(r(nn,3)-r(nn,1) +xxg*(a3(nn)
     &                   +xxg*a4(nn)))
                end do
              else
                do nn=1,ntr
                  sextra(iproc)%a(nn+(iq-1)*ntr) =
     &                   ((1.-xxg)*((2.-xxg)*
     &                   ((1.+xxg)*r(nn,2)-xxg*r(nn,1)/3.)
     &                   -xxg*(1.+xxg)*r(nn,4)/3.)
     &                   +xxg*(1.+xxg)*(2.-xxg)*r(nn,3))/2.
                end do
              endif         !  (mhint==2)
            enddo            ! iq loop
          end do              ! iproc
        else                   ! (nfield<mh_bs)
          do ii=1,neighnum
            iproc=neighlistsend(ii)
            do iq=1,drlen(iproc)
              !  Convert face index from 0:npanels to array indices
              ip = min(il_g,max(1,nint(dpoints(iproc)%a(2,iq))))
              jp = min(il_g,max(1,nint(dpoints(iproc)%a(3,iq))))
              n = nint(dpoints(iproc)%a(1,iq)) + noff ! Local index
              !  Need global face index in fproc call
#ifdef debug
              if ( fproc(ip,jp,n-noff) /= myid ) then
                 print*, "Error in ints_bl", myid, n, iq, iproc,
     &             dpoints(iproc)%a(:,iq)
                 call ccmpi_abort(-1)
              end if
#endif
              idel = int(dpoints(iproc)%a(2,iq))
              xxg = dpoints(iproc)%a(2,iq) - idel
              jdel = int(dpoints(iproc)%a(3,iq))
              yyg = dpoints(iproc)%a(3,iq) - jdel
              k = nint(dpoints(iproc)%a(4,iq))
              idel = idel - ioff(n-noff)
              jdel = jdel - joff(n-noff)
              c1 = sx(:,idel,jdel-1,n,k) ! manually unrolled loop
              c2 = sx(:,idel,jdel  ,n,k)
              c3 = sx(:,idel,jdel+1,n,k)
              c4 = sx(:,idel,jdel+2,n,k)
              cmin = min( 1.e20,c2,c3) ! Bermejo & Staniforth
              cmax = max(-1.e20,c2,c3) ! Bermejo & Staniforth
              if(mhint==2)then ! Bessel interp
                a4 = c4-c1+3.*(c2-c3)
                a3 = c1-2.*c2+c3-a4
                r(:,2) = c2+.5*yyg*(c3-c1 +yyg*(a3+yyg*a4))
              else
                r(:,2) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*c2
     &                    -yyg*c1/3.)
     &                    -yyg*(1.+yyg)*c4/3.)
     &                    +yyg*(1.+yyg)*(2.-yyg)*c3)/2.
              endif         !  (mhint==2)
              c1 = sx(:,idel+1,jdel-1,n,k)
              c2 = sx(:,idel+1,jdel  ,n,k)
              c3 = sx(:,idel+1,jdel+1,n,k)
              c4 = sx(:,idel+1,jdel+2,n,k)
              cmin = min(cmin,c2,c3) ! Bermejo & Staniforth
              cmax = max(cmax,c2,c3) ! Bermejo & Staniforth
              if(mhint==2)then ! Bessel interp
                a4 = c4-c1+3.*(c2-c3)
                a3 = c1-2.*c2+c3-a4
                r(:,3) = c2+.5*yyg*(c3-c1 +yyg*(a3+yyg*a4))
              else
                r(:,3) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*c2
     &                    -yyg*c1/3.)
     &                    -yyg*(1.+yyg)*c4/3.)
     &                    +yyg*(1.+yyg)*(2.-yyg)*c3)/2.
              endif         !  (mhint==2)
c           r = {(1-y)*{(2-y     )*[(1+y     )*c2-y     *c1/3]
c                -y*(1+y)*c4/3}
c                +y*(1+y)*(2-y)*c3}/2
              do nn=1,4,3   ! N.B.
                c2 = sx(:,idel+nn-2,jdel  ,n,k)
                c3 = sx(:,idel+nn-2,jdel+1,n,k)
                r(:,nn) = (1.-yyg)*c2 +yyg*c3
              enddo         ! nn loop
              if(mhint==2)then ! Bessel interp
                a4 = r(:,4)-r(:,1)+3.*(r(:,2)-r(:,3))
                a3 = r(:,1)-2.*r(:,2)+r(:,3)-a4
                sss = r(:,2)+.5*xxg*(r(:,3)-r(:,1)
     &               +xxg*(a3+xxg*a4))
              else
                sss = ((1.-xxg)*((2.-xxg)*
     &               ((1.+xxg)*r(:,2)-xxg*r(:,1)/3.)
     &               -xxg*(1.+xxg)*r(:,4)/3.)
     &               +xxg*(1.+xxg)*(2.-xxg)*r(:,3))/2.
              endif         !  (mhint==2)
              do nn=1,ntr
                sextra(iproc)%a(nn+(iq-1)*ntr) =
     &             min(max(cmin(nn),sss(nn)),cmax(nn)) ! Bermejo & Staniforth
              end do
            enddo            ! iq loop
          end do              ! iproc
        endif                  ! (nfield<mh_bs)  .. else ..

        call intssync_send(ntr)

        if(nfield<mh_bs)then
          do k=1,kl
            do iq=1,ifull    ! non Berm-Stan option
!             Convert face index from 0:npanels to array indices
              idel=int(xg(iq,k))
              xxg=xg(iq,k)-idel
              jdel=int(yg(iq,k))
              yyg=yg(iq,k)-jdel
              ! Now make them proper indices in this processor's region
              idel = idel - ioff(nface(iq,k))
              jdel = jdel - joff(nface(iq,k))
              n = nface(iq,k) + noff ! Make this a local index

              if ( idel < 0 .or. idel > ipan .or. jdel < 0 .or.
     &             jdel > jpan .or. n < 1 .or. n > npan ) then
                 cycle      ! Will be calculated on another processor
              end if

              c1 = sx(:,idel,jdel-1,n,k) ! manually unrolled loop
              c2 = sx(:,idel,jdel  ,n,k)
              c3 = sx(:,idel,jdel+1,n,k)
              c4 = sx(:,idel,jdel+2,n,k)
              if(mhint==2)then ! Bessel interp
                a4 = c4-c1+3.*(c2-c3)
                a3 = c1-2.*c2+c3-a4
                r(:,2) = c2+.5*yyg*(c3-c1 +yyg*(a3+yyg*a4))
              else
                r(:,2) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*c2
     &                    -yyg*c1/3.)
     &                    -yyg*(1.+yyg)*c4/3.)
     &                    +yyg*(1.+yyg)*(2.-yyg)*c3)/2.
              endif         !  (mhint==2)
              c1 = sx(:,idel+1,jdel-1,n,k)
              c2 = sx(:,idel+1,jdel  ,n,k)
              c3 = sx(:,idel+1,jdel+1,n,k)
              c4 = sx(:,idel+1,jdel+2,n,k)
              if(mhint==2)then ! Bessel interp
                a4 = c4-c1+3.*(c2-c3)
                a3 = c1-2.*c2+c3-a4
                r(:,3) = c2+.5*yyg*(c3-c1 +yyg*(a3+yyg*a4))
              else
                r(:,3) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*c2
     &                    -yyg*c1/3.)
     &                    -yyg*(1.+yyg)*c4/3.)
     &                    +yyg*(1.+yyg)*(2.-yyg)*c3)/2.
              endif         !  (mhint==2)
c          r = {(1-y)*{(2-y     )*[(1+y     )*c2-y     *c1/3]
c               -y*(1+y)*c4/3}
c               +y*(1+y)*(2-y)*c3}/2
              do nn=1,4,3   ! N.B.
                c2 = sx(:,idel+nn-2,jdel  ,n,k)
                c3 = sx(:,idel+nn-2,jdel+1,n,k)
                r(:,nn) = (1.-yyg)*c2 +yyg*c3
              enddo         ! nn loop
              if(mhint==2)then ! Bessel interp
                a4 = r(:,4)-r(:,1)+3.*(r(:,2)-r(:,3))
                a3 = r(:,1)-2.*r(:,2)+r(:,3)-a4
                s(iq,k,:) = r(:,2)+.5*xxg*(r(:,3)-r(:,1)
     &                   +xxg*(a3+xxg*a4))
              else
                    s(iq,k,:) = ((1.-xxg)*((2.-xxg)*
     &                   ((1.+xxg)*r(:,2)-xxg*r(:,1)/3.)
     &                   -xxg*(1.+xxg)*r(:,4)/3.)
     &                   +xxg*(1.+xxg)*(2.-xxg)*r(:,3))/2.
              endif         !  (mhint==2)
            enddo            ! iq loop
          enddo               ! k loop
        else                ! (nfield<mh_bs)
          do k=1,kl
            do iq=1,ifull    ! Berm-Stan option here e.g. qg & gases
              idel=int(xg(iq,k))
              xxg=xg(iq,k)-idel
              jdel=int(yg(iq,k))
              yyg=yg(iq,k)-jdel
              ! Now make them proper indices in this processor's region
              idel = idel - ioff(nface(iq,k))
              jdel = jdel - joff(nface(iq,k))
              n = nface(iq,k) + noff ! Make this a local index

              if ( idel < 0 .or. idel > ipan .or. jdel < 0 .or.
     &             jdel > jpan .or. n < 1 .or. n > npan ) then
                cycle      ! Will be calculated on another processor
              end if

              c1 = sx(:,idel,jdel-1,n,k) ! manually unrolled loop
              c2 = sx(:,idel,jdel  ,n,k)
              c3 = sx(:,idel,jdel+1,n,k)
              c4 = sx(:,idel,jdel+2,n,k)
              cmin = min( 1.e20,c2,c3) ! Bermejo & Staniforth
              cmax = max(-1.e20,c2,c3) ! Bermejo & Staniforth
              if(mhint==2)then ! Bessel interp
                a4 = c4-c1+3.*(c2-c3)
                a3 = c1-2.*c2+c3-a4
                r(:,2) = c2+.5*yyg*(c3-c1 +yyg*(a3+yyg*a4))
              else
                r(:,2) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*c2
     &                   -yyg*c1/3.)
     &                   -yyg*(1.+yyg)*c4/3.)
     &                   +yyg*(1.+yyg)*(2.-yyg)*c3)/2.
              endif         !  (mhint==2)
              c1 = sx(:,idel+1,jdel-1,n,k)
              c2 = sx(:,idel+1,jdel  ,n,k)
              c3 = sx(:,idel+1,jdel+1,n,k)
              c4 = sx(:,idel+1,jdel+2,n,k)
              cmin = min(cmin,c2,c3) ! Bermejo & Staniforth
              cmax = max(cmax,c2,c3) ! Bermejo & Staniforth
              if(mhint==2)then ! Bessel interp
                a4 = c4-c1+3.*(c2-c3)
                a3 = c1-2.*c2+c3-a4
                r(:,3) = c2+.5*yyg*(c3-c1 +yyg*(a3+yyg*a4))
              else
                r(:,3) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*c2
     &                   -yyg*c1/3.)
     &                   -yyg*(1.+yyg)*c4/3.)
     &                   +yyg*(1.+yyg)*(2.-yyg)*c3)/2.
              endif         !  (mhint==2)
c           r = {(1-y)*{(2-y     )*[(1+y     )*c2-y     *c1/3]
c                -y*(1+y)*c4/3}
c                +y*(1+y)*(2-y)*c3}/2
              do nn=1,4,3   ! N.B.
                c2 = sx(:,idel+nn-2,jdel  ,n,k)
                c3 = sx(:,idel+nn-2,jdel+1,n,k)
                r(:,nn) = (1.-yyg)*c2 +yyg*c3
              enddo         ! nn loop
              if(mhint==2)then ! Bessel interp
                a4 = r(:,4)-r(:,1)+3.*(r(:,2)-r(:,3))
                a3 = r(:,1)-2.*r(:,2)+r(:,3)-a4
                sss = r(:,2)+.5*xxg*(r(:,3)-r(:,1)
     &                    +xxg*(a3+xxg*a4))
              else
                sss = ((1.-xxg)*((2.-xxg)*
     &                   ((1.+xxg)*r(:,2)-xxg*r(:,1)/3.)
     &                   -xxg*(1.+xxg)*r(:,4)/3.)
     &                   +xxg*(1.+xxg)*(2.-xxg)*r(:,3))/2.
              endif         !  (mhint==2)
              s(iq,k,:) = min(max(cmin,sss),cmax) ! Bermejo & Staniforth
            enddo            ! iq loop
          enddo               ! k loop
        endif               ! (nfield<mh_bs)  .. else ..

      endif                     ! (intsch==1) .. else ..
!========================   end of intsch=1 section ====================

      call intssync_recv(s)
      
      call end_log(ints_end)
      return
      end

      subroutine ints_bl(s,intsch,nface,xg,yg)  ! not usually called

      use cc_mpi
      use indices_m
      implicit none
      integer, parameter :: ntest=0
      include 'newmpar.h'
      include 'parm.h'
      include 'parmhor.h'    ! has mh_bs
      integer intsch
      integer nface(ifull,kl)
      real xg(ifull,kl),yg(ifull,kl)      ! now passed through call
      ! Need work common for this
      real sx(0:ipan+1,0:jpan+1,1:npan,kl)
      real s(ifull+iextra,kl,1)
      real duma(ifull+iextra,kl)
      integer idel, iq, jdel, nn
      real xxg, yyg
      integer i, j, k, n, ind, ip, jp, iproc
      integer ii
      ! This is really indp, just repeated here to get inlining to work
      ind(i,j,n)=i+(j-1)*ipan+(n-1)*ipan*jpan  ! *** for n=1,npan

c     this one does bi-linear interpolation only
c     first extend s arrays into sx - this one -1:il+2 & -1:il+2
c                    but for bi-linear only need 0:il+1 &  0:il+1
      call start_log(ints_begin)
      call bounds(s,corner=.true.)
      do k=1,kl
         do n=1,npan
            do j=1,jpan
               do i=1,ipan
                  sx(i,j,n,k) = s(ind(i,j,n),k,1)
               enddo            ! i loop
            enddo               ! j loop
            do j=1,jpan
               sx(0,j,n,k) = s(iw(ind(1,j,n)),k,1)
               sx(ipan+1,j,n,k) = s(ie(ind(ipan,j,n)),k,1)
            enddo               ! j loop
            do i=1,ipan
               sx(i,0,n,k) = s(is(ind(i,1,n)),k,1)
               sx(i,jpan+1,n,k) = s(in(ind(i,jpan,n)),k,1)
            enddo               ! i loop

            sx(0,0,n,k) = s(iws(ind(1,1,n)),k,1)
            sx(ipan+1,0,n,k) = s(ies(ind(ipan,1,n)),k,1)
            sx(0,jpan+1,n,k)    = s(iwn(ind(1,jpan,n)),k,1)
            sx(ipan+1,jpan+1,n,k) = s(ien(ind(ipan,jpan,n)),k,1)
         enddo                  ! n loop
      enddo                     ! k loop

! Loop over points that need to be calculated for other processes
      do ii=1,neighnum
         iproc=neighlistsend(ii)
         do iq=1,drlen(iproc)
           !  Convert face index from 0:npanels to array indices
            ip = min(il_g,max(1,nint(dpoints(iproc)%a(2,iq))))
            jp = min(il_g,max(1,nint(dpoints(iproc)%a(3,iq))))
            n = nint(dpoints(iproc)%a(1,iq)) + noff ! Local index
         !  Need global face index in fproc call
#ifdef debug
            if ( fproc(ip,jp,n-noff) /= myid ) then
               print*, "Error in ints_bl", myid, n, iq, iproc,
     &              dpoints(iproc)%a(:,iq)
               call ccmpi_abort(-1)
            end if
#endif
            idel = int(dpoints(iproc)%a(2,iq))
            xxg = dpoints(iproc)%a(2,iq) - idel
            jdel = int(dpoints(iproc)%a(3,iq))
            yyg = dpoints(iproc)%a(3,iq) - jdel
            k = nint(dpoints(iproc)%a(4,iq))
            idel = idel - ioff(n-noff)
            jdel = jdel - joff(n-noff)
            sextra(iproc)%a(iq) = yyg*( xxg*sx(idel+1,jdel+1,n,k)
     &                              +(1.-xxg)*sx(idel,jdel+1,n,k))
     &                   +(1.-yyg)*(      xxg*sx(idel+1,jdel,n,k)
     &                               +(1.-xxg)*sx(idel,jdel,n,k))
         end do
      end do

      call intssync_send(1)

      do k=1,kl
         do iq=1,ifull
!           Convert face index from 0:npanels to array indices
            idel=int(xg(iq,k))
            xxg=xg(iq,k)-idel
            jdel=int(yg(iq,k))
            yyg=yg(iq,k)-jdel
            ! Now make them proper indices in this processor's region
            idel = idel - ioff(nface(iq,k))
            jdel = jdel - joff(nface(iq,k))
            n = nface(iq,k) + noff ! Make this a local index

            if ( idel < 0 .or. idel > ipan .or. jdel < 0 .or.
     &           jdel > jpan .or. n < 1 .or. n > npan ) then
               cycle            ! Will be calculated on another processor
            end if

            s(iq,k,1) =      yyg*(      xxg*sx(idel+1,jdel+1,n,k)
     &                            +(1.-xxg)*sx(idel,jdel+1,n,k))
     &                 +(1.-yyg)*(      xxg*sx(idel+1,jdel,n,k)
     &                            +(1.-xxg)*sx(idel,jdel,n,k))
         enddo                  ! iq loop
      end do                    ! k

      call intssync_recv(s)

      call end_log(ints_end)
      return
      end
