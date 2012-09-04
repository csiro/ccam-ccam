      subroutine ints(s,intsch,nface,xg,yg,nfield)    ! with ints_bl entry
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
      include 'mpif.h'
      integer intsch, nfield
      integer nface(ifull,kl)
      real xg(ifull,kl),yg(ifull,kl)      ! now passed through call
      ! Need work common for this
      real sx(-1:ipan+2,-1:jpan+2,1:npan,kl)
      real s(ifull+iextra,kl),r(4)
      integer idel, iq, jdel, nn
      real a3, a4, c1, c2, c3, c4, cmax, cmin, sss, xxg, yyg
      integer i, j, k, n, ind, ip, jp, iproc, ierr
      ! This is really indp, just repeated here to get inlining to work
      ind(i,j,n)=i+(j-1)*ipan+(n-1)*ipan*jpan  ! *** for n=1,npan

      call start_log(ints_begin)
      call bounds(s,nrows=2)

!======================== start of intsch=1 section ====================
      if(intsch==1)then

         do k=1,kl              ! A single main k loop uses cache better

            do n=1,npan         ! first simple copy into larger array
               do j=1,jpan
                  do i=1,ipan
                     sx(i,j,n,k) = s(ind(i,j,n),k)
                  enddo         ! i loop
               enddo            ! j loop
            enddo               ! n loop
            
c       this is intsb           EW interps done first
c       first extend s arrays into sx - this one -1:il+2 & -1:il+2
            do n=1,npan
               do j=1,jpan
                  sx(0,j,n,k) = s(iw(ind(1,j,n)),k)
                  sx(-1,j,n,k) = s(iww(ind(1,j,n)),k)
                  sx(ipan+1,j,n,k) = s(ie(ind(ipan,j,n)),k)
                  sx(ipan+2,j,n,k) = s(iee(ind(ipan,j,n)),k)
               enddo            ! j loop
               do i=1,ipan
                  sx(i,0,n,k) = s(is(ind(i,1,n)),k)
                  sx(i,-1,n,k) = s(iss(ind(i,1,n)),k)
                  sx(i,jpan+1,n,k) = s(in(ind(i,jpan,n)),k)
                  sx(i,jpan+2,n,k) = s(inn(ind(i,jpan,n)),k)
               enddo            ! i loop
c        for ew interpolation, sometimes need (different from ns):
c            (-1,0),   (0,0),   (0,-1)   (-1,il+1),   (0,il+1),   (0,il+2)
c           (il+1,0),(il+2,0),(il+1,-1) (il+1,il+1),(il+2,il+1),(il+1,il+2)
!        sx(-1,0,n)=sx(0,2,n)            !  s(iwws6(1,1,n))
!        sx(0,0,n)=sx(0,1,n)             !  s(iws6(1,1,n))
!        sx(0,-1,n)=sx(-1,1,n)           !  s(iwss6(1,1,n)) or s(isws6(1,1,n))
!        sx(il+1,0,n)=sx(il+1,1,n)       !  s(ies6(il,1,n))
!        sx(il+2,0,n)=sx(il+1,2,n)       !  s(iees6(il,1,n))
!        sx(il+1,-1,n)=sx(il+2,1,n)      !  s(iess6(il,1,n))
!        sx(-1,il+1,n)=sx(0,il-1,n)      !  s(iwwn6(1,il,n))
!        sx(0,il+2,n)=sx(-1,il,n)        !  s(iwnn6(1,il,n))
!        sx(il+2,il+1,n)=sx(il+1,il-1,n) !  s(ieen6(il,il,n))
!        sx(il+1,il+2,n)=sx(il+2,il,n)   !  s(ienn6(il,il,n))
!        sx(0,il+1,n)=sx(0,il,n)         !  s(iwn6(1,il,n))   **
!        sx(il+1,il+1,n)=sx(il+1,il,n)   !  s(ien6(il,il,n))  **

               sx(-1,0,n,k) = s(lwws(n),k)
               sx(0,0,n,k) = s(iws(ind(1,1,n)),k)
               sx(0,-1,n,k) = s(lwss(n),k)
               sx(ipan+1,0,n,k) = s(ies(ind(ipan,1,n)),k)
               sx(ipan+2,0,n,k) = s(lees(n),k)
               sx(ipan+1,-1,n,k) = s(less(n),k)
               sx(-1,jpan+1,n,k) = s(lwwn(n),k)
               sx(0,jpan+2,n,k) = s(lwnn(n),k)
               sx(ipan+2,jpan+1,n,k) = s(leen(n),k)
               sx(ipan+1,jpan+2,n,k) = s(lenn(n),k)
               sx(0,jpan+1,n,k)    = s(iwn(ind(1,jpan,n)),k)
               sx(ipan+1,jpan+1,n,k) = s(ien(ind(ipan,jpan,n)),k)
            enddo               ! n loop

            if(ntest==1.and.mydiag.and.nfield==4.and.k==nlv)then !just moisture
!             just set up for single processor
              if(nproc==1)then
                n=1+jd/il
                jdel=jd+il-n*il
                print *,'qg,qlg,qfg for id,jd,n,jdel ',id,jd,n,jdel
                print *,'ipan,jpan,npan ',ipan,jpan,npan
                iq=idjd
                print *,'xg,yg,nface ',xg(iq,k),yg(iq,k),nface(iq,k)
                print *,'ioff,joff,noff ',ioff(:),joff(:),noff
                do j=jdel+2,jdel-2,-1
                 write (6,"('qg#5 ',5f8.3)") 
     &            (1000.*sx(i,j,n,k),i=id-2,id+2)
                enddo
              endif
            endif

            if(nfield<mh_bs)then
               do iq=1,ifull    ! non Berm-Stan option
!                 Convert face index from 0:npanels to array indices
                  idel=int(xg(iq,k))
                  xxg=xg(iq,k)-idel
                  jdel=int(yg(iq,k))
                  yyg=yg(iq,k)-jdel
                  ! Now make them proper indices in this processor's region
                  idel = idel - ioff(nface(iq,k))
                  jdel = jdel - joff(nface(iq,k))
                  n = nface(iq,k) + noff ! Make this a local index
#ifdef SX
                  ! For better vectorisation calculate all values
                  idel = max(0,min(idel,ipan))
                  jdel = max(0,min(jdel,jpan))
                  n = max(1,min(n,npan))
#else
                  if ( idel < 0 .or. idel > ipan .or. jdel < 0 .or.
     &                 jdel > jpan .or. n < 1 .or. n > npan ) then
                     cycle      ! Will be calculated on another processor
                  end if
#endif
c          if(ntest==1)then
c            if(xxg<0..or.xxg>1..or.yyg<0..or.yyg>1.)
c    &        print *,'intsb problem iq xxg,yyg ',iq,xxg,yyg
c          endif
                  c1 = sx(idel-1,jdel,n,k) ! manually unrolled loop
                  c2 = sx(idel  ,jdel,n,k)
                  c3 = sx(idel+1,jdel,n,k)
                  c4 = sx(idel+2,jdel,n,k)
                  if(mhint==2)then ! Bessel interp
                     a4 = c4-c1+3.*(c2-c3)
                     a3 = c1-2.*c2+c3-a4
                     r(2) = c2+.5*xxg*(c3-c1 +xxg*(a3+xxg*a4))
                  else
                     r(2) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*c2-xxg*c1/3.)
     &                    -xxg*(1.+xxg)*c4/3.)
     &                    +xxg*(1.+xxg)*(2.-xxg)*c3)/2.
                  endif         !  (mhint==2)
                  c1 = sx(idel-1,jdel+1,n,k)
                  c2 = sx(idel  ,jdel+1,n,k)
                  c3 = sx(idel+1,jdel+1,n,k)
                  c4 = sx(idel+2,jdel+1,n,k)
                  if(mhint==2)then ! Bessel interp
                     a4 = c4-c1+3.*(c2-c3)
                     a3 = c1-2.*c2+c3-a4
                     r(3) = c2+.5*xxg*(c3-c1 +xxg*(a3+xxg*a4))
                  else
                     r(3) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*c2-xxg*c1/3.)
     &                    -xxg*(1.+xxg)*c4/3.)
     &                    +xxg*(1.+xxg)*(2.-xxg)*c3)/2.
                  endif         !  (mhint==2)
c           r = {(1-x)*{(2-x     )*[(1+x     )*c2-x     *c1/3]
c                -x*(1+x)*c4/3}
c                +x*(1+x)*(2-x)*c3}/2
                  do nn=1,4,3   ! N.B.
                     c2 = sx(idel  ,jdel+nn-2,n,k)
                     c3 = sx(idel+1,jdel+nn-2,n,k)
                     r(nn) = (1.-xxg)*c2 +xxg*c3
                  enddo         ! nn loop
                  if(mhint==2)then ! Bessel interp
                     a4 = r(4)-r(1)+3.*(r(2)-r(3))
                     a3 = r(1)-2.*r(2)+r(3)-a4
                     s(iq,k) = r(2)+.5*yyg*(r(3)-r(1) +yyg*(a3+yyg*a4))
                  else
                     s(iq,k) = ((1.-yyg)*((2.-yyg)*
     &                    ((1.+yyg)*r(2)-yyg*r(1)/3.)
     &                    -yyg*(1.+yyg)*r(4)/3.)
     &                    +yyg*(1.+yyg)*(2.-yyg)*r(3))/2.
                  endif         !  (mhint==2)
               enddo            ! iq loop
            else                ! (nfield<mh_bs)
               do iq=1,ifull    ! Berm-Stan option here e.g. qg & gases
                  idel=int(xg(iq,k))
                  xxg=xg(iq,k)-idel
                  jdel=int(yg(iq,k))
                  yyg=yg(iq,k)-jdel
                  ! Now make them proper indices in this processor's region
                  idel = idel - ioff(nface(iq,k))
                  jdel = jdel - joff(nface(iq,k))
                  n = nface(iq,k) + noff ! Make this a local index
#ifdef SX
                  ! For better vectorisation calculate all values
                  idel = max(0,min(idel,ipan))
                  jdel = max(0,min(jdel,jpan))
                  n = max(1,min(n,npan))
#else
                  if ( idel < 0 .or. idel > ipan .or. jdel < 0 .or.
     &                 jdel > jpan .or. n < 1 .or. n > npan ) then
                     cycle      ! Will be calculated on another processor
                  end if
#endif
                  c1 = sx(idel-1,jdel,n,k) ! manually unrolled loop
                  c2 = sx(idel  ,jdel,n,k)
                  c3 = sx(idel+1,jdel,n,k)
                  c4 = sx(idel+2,jdel,n,k)
                  cmin = min( 1.e20,c2,c3) ! Bermejo & Staniforth
                  cmax = max(-1.e20,c2,c3) ! Bermejo & Staniforth
                  if(mhint==2)then ! Bessel interp
                     a4 = c4-c1+3.*(c2-c3)
                     a3 = c1-2.*c2+c3-a4
                     r(2) = c2+.5*xxg*(c3-c1 +xxg*(a3+xxg*a4))
                  else
                     r(2) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*c2-xxg*c1/3.)
     &                    -xxg*(1.+xxg)*c4/3.)
     &                    +xxg*(1.+xxg)*(2.-xxg)*c3)/2.
                  endif         !  (mhint==2)
                  c1 = sx(idel-1,jdel+1,n,k)
                  c2 = sx(idel  ,jdel+1,n,k)
                  c3 = sx(idel+1,jdel+1,n,k)
                  c4 = sx(idel+2,jdel+1,n,k)
                  cmin = min(cmin,c2,c3) ! Bermejo & Staniforth
                  cmax = max(cmax,c2,c3) ! Bermejo & Staniforth
                  if(mhint==2)then ! Bessel interp
                     a4 = c4-c1+3.*(c2-c3)
                     a3 = c1-2.*c2+c3-a4
                     r(3) = c2+.5*xxg*(c3-c1 +xxg*(a3+xxg*a4))
                  else
                     r(3) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*c2-xxg*c1/3.)
     &                    -xxg*(1.+xxg)*c4/3.)
     &                    +xxg*(1.+xxg)*(2.-xxg)*c3)/2.
                  endif         !  (mhint==2)
c           r = {(1-x)*{(2-x     )*[(1+x     )*c2-x     *c1/3]
c                -x*(1+x)*c4/3}
c                +x*(1+x)*(2-x)*c3}/2
                  do nn=1,4,3   ! N.B.
                     c2 = sx(idel  ,jdel+nn-2,n,k)
                     c3 = sx(idel+1,jdel+nn-2,n,k)
                     r(nn) = (1.-xxg)*c2 +xxg*c3
                  enddo         ! nn loop
                  if(mhint==2)then ! Bessel interp
                     a4 = r(4)-r(1)+3.*(r(2)-r(3))
                     a3 = r(1)-2.*r(2)+r(3)-a4
                     sss = r(2)+.5*yyg*(r(3)-r(1) +yyg*(a3+yyg*a4))
                  else
                     sss = ((1.-yyg)*((2.-yyg)*
     &                    ((1.+yyg)*r(2)-yyg*r(1)/3.)
     &                    -yyg*(1.+yyg)*r(4)/3.)
     &                    +yyg*(1.+yyg)*(2.-yyg)*r(3))/2.
                  endif         !  (mhint==2)
                  s(iq,k) = min(max(cmin,sss),cmax) ! Bermejo & Staniforth
               enddo            ! iq loop
            endif               ! (nfield<mh_bs)  .. else ..
            
         end do                 ! k

! Loop over points that need to be calculated for other processes

         if(nfield<mh_bs)then
            do iproc=0,nproc-1
               if ( iproc == myid ) then
                  cycle
               end if
!!!         print*, "INTS OTHER", myid, iproc, drlen(iproc)
               do iq=1,drlen(iproc)
                  !  Convert face index from 0:npanels to array indices
                  ip = min(il_g,max(1,nint(dpoints(iproc)%a(2,iq))))    ! MJT memory
                  jp = min(il_g,max(1,nint(dpoints(iproc)%a(3,iq))))    ! MJT memory
                  n = nint(dpoints(iproc)%a(1,iq)) + noff ! Local index ! MJT memory
                  !  Need global face index in fproc call
#ifdef debug
                  if ( fproc(ip,jp,n-noff) /= myid ) then
                     print*, "Error in ints", myid, n, iq, iproc,
     &                     dpoints(iproc)%a(:,iq)                       ! MJT memory
                     call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
                  end if
#endif
                  idel = int(dpoints(iproc)%a(2,iq))                    ! MJT memory
                  xxg = dpoints(iproc)%a(2,iq) - idel                   ! MJT memory
                  jdel = int(dpoints(iproc)%a(3,iq))                    ! MJT memory
                  yyg = dpoints(iproc)%a(3,iq) - jdel                   ! MJT memory
                  k = nint(dpoints(iproc)%a(4,iq))                      ! MJT memory
                  idel = idel - ioff(n-noff)
                  jdel = jdel - joff(n-noff)
                  c1 = sx(idel-1,jdel,n,k) ! manually unrolled loop
                  c2 = sx(idel  ,jdel,n,k)
                  c3 = sx(idel+1,jdel,n,k)
                  c4 = sx(idel+2,jdel,n,k)
                  if(mhint==2)then ! Bessel interp
                     a4 = c4-c1+3.*(c2-c3)
                     a3 = c1-2.*c2+c3-a4
                     r(2) = c2+.5*xxg*(c3-c1 +xxg*(a3+xxg*a4))
                  else
                     r(2) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*c2-xxg*c1/3.)
     &                    -xxg*(1.+xxg)*c4/3.)
     &                    +xxg*(1.+xxg)*(2.-xxg)*c3)/2.
                  endif         !  (mhint==2)
                  c1 = sx(idel-1,jdel+1,n,k)
                  c2 = sx(idel  ,jdel+1,n,k)
                  c3 = sx(idel+1,jdel+1,n,k)
                  c4 = sx(idel+2,jdel+1,n,k)
                  if(mhint==2)then ! Bessel interp
                     a4 = c4-c1+3.*(c2-c3)
                     a3 = c1-2.*c2+c3-a4
                     r(3) = c2+.5*xxg*(c3-c1 +xxg*(a3+xxg*a4))
                  else
                     r(3) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*c2-xxg*c1/3.)
     &                    -xxg*(1.+xxg)*c4/3.)
     &                    +xxg*(1.+xxg)*(2.-xxg)*c3)/2.
                  endif         !  (mhint==2)
c           r = {(1-x)*{(2-x     )*[(1+x     )*c2-x     *c1/3]
c                -x*(1+x)*c4/3}
c                +x*(1+x)*(2-x)*c3}/2
                  do nn=1,4,3   ! N.B.
                     c2 = sx(idel  ,jdel+nn-2,n,k)
                     c3 = sx(idel+1,jdel+nn-2,n,k)
                     r(nn) = (1.-xxg)*c2 +xxg*c3
                  enddo         ! nn loop
                  if(mhint==2)then ! Bessel interp
                     a4 = r(4)-r(1)+3.*(r(2)-r(3))
                     a3 = r(1)-2.*r(2)+r(3)-a4
                     sextra(iproc)%a(iq) = r(2) +                       ! MJT memory
     &                    0.5*yyg*(r(3)-r(1) +yyg*(a3+yyg*a4))
                  else
                     sextra(iproc)%a(iq) = ((1.-yyg)*((2.-yyg)*         ! MJT memory
     &                    ((1.+yyg)*r(2)-yyg*r(1)/3.)
     &                    -yyg*(1.+yyg)*r(4)/3.)
     &                    +yyg*(1.+yyg)*(2.-yyg)*r(3))/2.
                  endif         !  (mhint==2)
               enddo            ! iq loop
            end do              ! iproc loop
         else                   ! (nfield<mh_bs)
            do iproc=0,nproc-1
               if ( iproc == myid ) then
                  cycle
               end if
!!!         print*, "INTS OTHER", myid, iproc, drlen(iproc)
               do iq=1,drlen(iproc)
                  !  Convert face index from 0:npanels to array indices
                  ip = min(il_g,max(1,nint(dpoints(iproc)%a(2,iq))))    ! MJT memory
                  jp = min(il_g,max(1,nint(dpoints(iproc)%a(3,iq))))    ! MJT memory
                  n = nint(dpoints(iproc)%a(1,iq)) + noff ! Local index ! MJT memory
                  !  Need global face index in fproc call
#ifdef debug
                  if ( fproc(ip,jp,n-noff) /= myid ) then
                     print*, "Error in ints", myid, n, iq, iproc,
     &                     dpoints(iproc)%a(:,iq)                       ! MJT memory
                     call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
                  end if
#endif
                  idel = int(dpoints(iproc)%a(2,iq))                    ! MJT memory
                  xxg = dpoints(iproc)%a(2,iq) - idel                   ! MJT memory
                  jdel = int(dpoints(iproc)%a(3,iq))                    ! MJT memory
                  yyg = dpoints(iproc)%a(3,iq) - jdel                   ! MJT memory
                  k = nint(dpoints(iproc)%a(4,iq))                      ! MJT memory
                  idel = idel - ioff(n-noff)
                  jdel = jdel - joff(n-noff)
                  c1 = sx(idel-1,jdel,n,k) ! manually unrolled loop
                  c2 = sx(idel  ,jdel,n,k)
                  c3 = sx(idel+1,jdel,n,k)
                  c4 = sx(idel+2,jdel,n,k)
                  cmin = min( 1.e20,c2,c3) ! Bermejo & Staniforth
                  cmax = max(-1.e20,c2,c3) ! Bermejo & Staniforth
                  if(mhint==2)then ! Bessel interp
                     a4 = c4-c1+3.*(c2-c3)
                     a3 = c1-2.*c2+c3-a4
                     r(2) = c2+.5*xxg*(c3-c1 +xxg*(a3+xxg*a4))
                  else
                     r(2) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*c2-xxg*c1/3.)
     &                    -xxg*(1.+xxg)*c4/3.)
     &                    +xxg*(1.+xxg)*(2.-xxg)*c3)/2.
                  endif         !  (mhint==2)
                  c1 = sx(idel-1,jdel+1,n,k)
                  c2 = sx(idel  ,jdel+1,n,k)
                  c3 = sx(idel+1,jdel+1,n,k)
                  c4 = sx(idel+2,jdel+1,n,k)
                  cmin = min(cmin,c2,c3) ! Bermejo & Staniforth
                  cmax = max(cmax,c2,c3) ! Bermejo & Staniforth
                  if(mhint==2)then ! Bessel interp
                     a4 = c4-c1+3.*(c2-c3)
                     a3 = c1-2.*c2+c3-a4
                     r(3) = c2+.5*xxg*(c3-c1 +xxg*(a3+xxg*a4))
                  else
                     r(3) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*c2-xxg*c1/3.)
     &                    -xxg*(1.+xxg)*c4/3.)
     &                    +xxg*(1.+xxg)*(2.-xxg)*c3)/2.
                  endif         !  (mhint==2)
c           r = {(1-x)*{(2-x     )*[(1+x     )*c2-x     *c1/3]
c                -x*(1+x)*c4/3}
c                +x*(1+x)*(2-x)*c3}/2
                  do nn=1,4,3   ! N.B.
                     c2 = sx(idel  ,jdel+nn-2,n,k)
                     c3 = sx(idel+1,jdel+nn-2,n,k)
                     r(nn) = (1.-xxg)*c2 +xxg*c3
                  enddo         ! nn loop
                  if(mhint==2)then ! Bessel interp
                     a4 = r(4)-r(1)+3.*(r(2)-r(3))
                     a3 = r(1)-2.*r(2)+r(3)-a4
                     sss = r(2)+.5*yyg*(r(3)-r(1) +yyg*(a3+yyg*a4))
                  else
                     sss = ((1.-yyg)*((2.-yyg)*
     &                    ((1.+yyg)*r(2)-yyg*r(1)/3.)
     &                    -yyg*(1.+yyg)*r(4)/3.)
     &                    +yyg*(1.+yyg)*(2.-yyg)*r(3))/2.
                  endif         !  (mhint==2)
                  sextra(iproc)%a(iq) = min(max(cmin,sss),cmax) ! Bermejo & Staniforth ! MJT memory
               enddo            ! iq loop
            end do              ! iproc loop
         endif                  ! (nfield<mh_bs)  .. else ..
            
!========================   end of intsch=1 section ====================
      else     ! if(intsch==1)then
!======================== start of intsch=2 section ====================
c       this is intsc           NS interps done first
c       first extend s arrays into sx - this one -1:il+2 & -1:il+2
         do k=1,kl              ! A single main k loop uses cache better
            do n=1,npan         ! first simple copy into larger array
               do j=1,jpan
                  do i=1,ipan
                     sx(i,j,n,k)=s(ind(i,j,n),k)
                  enddo         ! i loop
               enddo            ! j loop
            enddo               ! n loop

            do n=1,npan
               do j=1,jpan
                  sx(0,j,n,k) = s(iw(ind(1,j,n)),k)
                  sx(-1,j,n,k) = s(iww(ind(1,j,n)),k)
                  sx(ipan+1,j,n,k) = s(ie(ind(ipan,j,n)),k)
                  sx(ipan+2,j,n,k) = s(iee(ind(ipan,j,n)),k)
               enddo            ! j loop
               do i=1,ipan
                  sx(i,0,n,k) = s(is(ind(i,1,n)),k)
                  sx(i,-1,n,k) = s(iss(ind(i,1,n)),k)
                  sx(i,jpan+1,n,k) = s(in(ind(i,jpan,n)),k)
                  sx(i,jpan+2,n,k) = s(inn(ind(i,jpan,n)),k)
               enddo            ! i loop
c        for ns interpolation, sometimes need (different from ew):
c            (-1,0),   (0,0),   (0,-1)   (-1,il+1),   (0,il+1),   (0,il+2)
c          (il+1,0),(il+2,0),(il+1,-1) (il+1,il+1),(il+2,il+1),(il+1,il+2)
!        sx(-1,0,n)=sx(1,-1,n)           !  s(isww6(1,1,n))
!        sx(0,0,n)=sx(1,0,n)             !  s(isw6(1,1,n))
!        sx(0,-1,n)=sx(2,0,n)            !  s(issw6(1,1,n))
!        sx(il+2,0,n)=sx(il,-1,n)        !  s(isee6(il,1,n))
!        sx(il+1,-1,n)=sx(il-1,0,n)      !  s(isse6(il,1,n))
!        sx(-1,il+1,n)=sx(1,il+2,n)      !  s(inww6(1,il,n))
!        sx(0,il+1,n)=sx(1,il+1,n)       !  s(inw6(1,il,n))
!        sx(0,il+2,n)=sx(2,il+1,n)       !  s(innw6(1,il,n))
!        sx(il+2,il+1,n)=sx(il,il+2,n)   !  s(inee6(il,il,n))
!        sx(il+1,il+2,n)=sx(il-1,il+1,n) !  s(inne6(il,il,n))
!        sx(il+1,0,n)=sx(il,0,n)         !  s(ise6(il,1,n))    **
!        sx(il+1,il+1,n)=sx(il,il+1,n)   !  s(ine6(il,il,n))   **

               sx(-1,0,n,k)=s(lsww(n),k)
               sx(0,0,n,k) = s(isw(ind(1,1,n)),k)
               sx(0,-1,n,k) = s(lssw(n),k)
               sx(ipan+2,0,n,k) = s(lsee(n),k)
               sx(ipan+1,-1,n,k) = s(lsse(n),k)
               sx(-1,jpan+1,n,k) = s(lnww(n),k)
               sx(0,jpan+1,n,k) = s(inw(ind(1,jpan,n)),k)
               sx(0,jpan+2,n,k) = s(lnnw(n),k)
               sx(ipan+2,jpan+1,n,k) = s(lnee(n),k)
               sx(ipan+1,jpan+2,n,k) = s(lnne(n),k)
               sx(ipan+1,0,n,k)    = s(ise(ind(ipan,1,n)),k)
               sx(ipan+1,jpan+1,n,k) = s(ine(ind(ipan,jpan,n)),k)
            enddo               ! n loop

            if(nfield<mh_bs)then
               do iq=1,ifull    ! non Berm-Stan option
!                 Convert face index from 0:npanels to array indices
                  idel=int(xg(iq,k))
                  xxg=xg(iq,k)-idel
                  jdel=int(yg(iq,k))
                  yyg=yg(iq,k)-jdel
                  ! Now make them proper indices in this processor's region
                  idel = idel - ioff(nface(iq,k))
                  jdel = jdel - joff(nface(iq,k))
                  n = nface(iq,k) + noff ! Make this a local index
#ifdef SX
                  ! For better vectorisation calculate all values
                  idel = max(0,min(idel,ipan))
                  jdel = max(0,min(jdel,jpan))
                  n = max(1,min(n,npan))
#else
                  if ( idel < 0 .or. idel > ipan .or. jdel < 0 .or.
     &                 jdel > jpan .or. n < 1 .or. n > npan ) then
                     cycle      ! Will be calculated on another processor
                  end if
#endif
                  c1 = sx(idel,jdel-1,n,k) ! manually unrolled loop
                  c2 = sx(idel,jdel  ,n,k)
                  c3 = sx(idel,jdel+1,n,k)
                  c4 = sx(idel,jdel+2,n,k)
                  if(mhint==2)then ! Bessel interp
                     a4 = c4-c1+3.*(c2-c3)
                     a3 = c1-2.*c2+c3-a4
                     r(2) = c2+.5*yyg*(c3-c1 +yyg*(a3+yyg*a4))
                  else
                     r(2) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*c2-yyg*c1/3.)
     &                    -yyg*(1.+yyg)*c4/3.)
     &                    +yyg*(1.+yyg)*(2.-yyg)*c3)/2.
                  endif         !  (mhint==2)
                  c1 = sx(idel+1,jdel-1,n,k)
                  c2 = sx(idel+1,jdel  ,n,k)
                  c3 = sx(idel+1,jdel+1,n,k)
                  c4 = sx(idel+1,jdel+2,n,k)
                  if(mhint==2)then ! Bessel interp
                     a4 = c4-c1+3.*(c2-c3)
                     a3 = c1-2.*c2+c3-a4
                     r(3) = c2+.5*yyg*(c3-c1 +yyg*(a3+yyg*a4))
                  else
                     r(3) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*c2-yyg*c1/3.)
     &                    -yyg*(1.+yyg)*c4/3.)
     &                    +yyg*(1.+yyg)*(2.-yyg)*c3)/2.
                  endif         !  (mhint==2)
c          r = {(1-y)*{(2-y     )*[(1+y     )*c2-y     *c1/3]
c               -y*(1+y)*c4/3}
c               +y*(1+y)*(2-y)*c3}/2
                  do nn=1,4,3   ! N.B.
                     c2 = sx(idel+nn-2,jdel  ,n,k)
                     c3 = sx(idel+nn-2,jdel+1,n,k)
                     r(nn) = (1.-yyg)*c2 +yyg*c3
                  enddo         ! nn loop
                  if(mhint==2)then ! Bessel interp
                     a4 = r(4)-r(1)+3.*(r(2)-r(3))
                     a3 = r(1)-2.*r(2)+r(3)-a4
                     s(iq,k) = r(2)+.5*xxg*(r(3)-r(1) +xxg*(a3+xxg*a4))
                  else
                     s(iq,k) = ((1.-xxg)*((2.-xxg)*
     &                    ((1.+xxg)*r(2)-xxg*r(1)/3.)
     &                    -xxg*(1.+xxg)*r(4)/3.)
     &                    +xxg*(1.+xxg)*(2.-xxg)*r(3))/2.
                  endif         !  (mhint==2)
               enddo            ! iq loop
            else                ! (nfield<mh_bs)
               do iq=1,ifull    ! Berm-Stan option here e.g. qg & gases
                  idel=int(xg(iq,k))
                  xxg=xg(iq,k)-idel
                  jdel=int(yg(iq,k))
                  yyg=yg(iq,k)-jdel
                  ! Now make them proper indices in this processor's region
                  idel = idel - ioff(nface(iq,k))
                  jdel = jdel - joff(nface(iq,k))
                  n = nface(iq,k) + noff ! Make this a local index
#ifdef SX
                  ! For better vectorisation calculate all values
                  idel = max(0,min(idel,ipan))
                  jdel = max(0,min(jdel,jpan))
                  n = max(1,min(n,npan))
#else
                  if ( idel < 0 .or. idel > ipan .or. jdel < 0 .or.
     &                 jdel > jpan .or. n < 1 .or. n > npan ) then
                     cycle      ! Will be calculated on another processor
                  end if
#endif
                  c1 = sx(idel,jdel-1,n,k) ! manually unrolled loop
                  c2 = sx(idel,jdel  ,n,k)
                  c3 = sx(idel,jdel+1,n,k)
                  c4 = sx(idel,jdel+2,n,k)
                  cmin = min( 1.e20,c2,c3) ! Bermejo & Staniforth
                  cmax = max(-1.e20,c2,c3) ! Bermejo & Staniforth
                  if(mhint==2)then ! Bessel interp
                     a4 = c4-c1+3.*(c2-c3)
                     a3 = c1-2.*c2+c3-a4
                     r(2) = c2+.5*yyg*(c3-c1 +yyg*(a3+yyg*a4))
                  else
                     r(2) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*c2-yyg*c1/3.)
     &                    -yyg*(1.+yyg)*c4/3.)
     &                    +yyg*(1.+yyg)*(2.-yyg)*c3)/2.
                  endif         !  (mhint==2)
                  c1 = sx(idel+1,jdel-1,n,k)
                  c2 = sx(idel+1,jdel  ,n,k)
                  c3 = sx(idel+1,jdel+1,n,k)
                  c4 = sx(idel+1,jdel+2,n,k)
                  cmin = min(cmin,c2,c3) ! Bermejo & Staniforth
                  cmax = max(cmax,c2,c3) ! Bermejo & Staniforth
                  if(mhint==2)then ! Bessel interp
                     a4 = c4-c1+3.*(c2-c3)
                     a3 = c1-2.*c2+c3-a4
                     r(3) = c2+.5*yyg*(c3-c1 +yyg*(a3+yyg*a4))
                  else
                     r(3) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*c2-yyg*c1/3.)
     &                    -yyg*(1.+yyg)*c4/3.)
     &                    +yyg*(1.+yyg)*(2.-yyg)*c3)/2.
                  endif         !  (mhint==2)
c           r = {(1-y)*{(2-y     )*[(1+y     )*c2-y     *c1/3]
c                -y*(1+y)*c4/3}
c                +y*(1+y)*(2-y)*c3}/2
                  do nn=1,4,3   ! N.B.
                     c2 = sx(idel+nn-2,jdel  ,n,k)
                     c3 = sx(idel+nn-2,jdel+1,n,k)
                     r(nn) = (1.-yyg)*c2 +yyg*c3
                  enddo         ! nn loop
                  if(mhint==2)then ! Bessel interp
                     a4 = r(4)-r(1)+3.*(r(2)-r(3))
                     a3 = r(1)-2.*r(2)+r(3)-a4
                     sss = r(2)+.5*xxg*(r(3)-r(1) +xxg*(a3+xxg*a4))
                  else
                     sss = ((1.-xxg)*((2.-xxg)*
     &                    ((1.+xxg)*r(2)-xxg*r(1)/3.)
     &                    -xxg*(1.+xxg)*r(4)/3.)
     &                    +xxg*(1.+xxg)*(2.-xxg)*r(3))/2.
                  endif         !  (mhint==2)
                  s(iq,k) = min(max(cmin,sss),cmax) ! Bermejo & Staniforth
               enddo            ! iq loop
            endif               ! (nfield<mh_bs)  .. else ..
         end do                 ! k

!        For other processes
         if(nfield<mh_bs)then
            do iproc=0,nproc-1
               if ( iproc == myid ) then
                  cycle
               end if
               do iq=1,drlen(iproc)
                  !  Convert face index from 0:npanels to array indices
                  ip = min(il_g,max(1,nint(dpoints(iproc)%a(2,iq))))    ! MJT memory
                  jp = min(il_g,max(1,nint(dpoints(iproc)%a(3,iq))))    ! MJT memory
                  n = nint(dpoints(iproc)%a(1,iq)) + noff ! Local index ! MJT memory
                  !  Need global face index in fproc call
#ifdef debug
                  if ( fproc(ip,jp,n-noff) /= myid ) then
                     print*, "Error in ints_bl", myid, n, iq, iproc,
     &                 dpoints(iproc)%a(:,iq)                           ! MJT memory
                     call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
                  end if
#endif
                  idel = int(dpoints(iproc)%a(2,iq))                    ! MJT memory
                  xxg = dpoints(iproc)%a(2,iq) - idel                   ! MJT memory
                  jdel = int(dpoints(iproc)%a(3,iq))                    ! MJT memory
                  yyg = dpoints(iproc)%a(3,iq) - jdel                   ! MJT memory
                  k = nint(dpoints(iproc)%a(4,iq))                      ! MJT memory
                  idel = idel - ioff(n-noff)
                  jdel = jdel - joff(n-noff)
                  c1 = sx(idel,jdel-1,n,k) ! manually unrolled loop
                  c2 = sx(idel,jdel  ,n,k)
                  c3 = sx(idel,jdel+1,n,k)
                  c4 = sx(idel,jdel+2,n,k)
                  if(mhint==2)then ! Bessel interp
                     a4 = c4-c1+3.*(c2-c3)
                     a3 = c1-2.*c2+c3-a4
                     r(2) = c2+.5*yyg*(c3-c1 +yyg*(a3+yyg*a4))
                  else
                     r(2) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*c2-yyg*c1/3.)
     &                    -yyg*(1.+yyg)*c4/3.)
     &                    +yyg*(1.+yyg)*(2.-yyg)*c3)/2.
                  endif         !  (mhint==2)
                  c1 = sx(idel+1,jdel-1,n,k)
                  c2 = sx(idel+1,jdel  ,n,k)
                  c3 = sx(idel+1,jdel+1,n,k)
                  c4 = sx(idel+1,jdel+2,n,k)
                  if(mhint==2)then ! Bessel interp
                     a4 = c4-c1+3.*(c2-c3)
                     a3 = c1-2.*c2+c3-a4
                     r(3) = c2+.5*yyg*(c3-c1 +yyg*(a3+yyg*a4))
                  else
                     r(3) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*c2-yyg*c1/3.)
     &                    -yyg*(1.+yyg)*c4/3.)
     &                    +yyg*(1.+yyg)*(2.-yyg)*c3)/2.
                  endif         !  (mhint==2)
c          r = {(1-y)*{(2-y     )*[(1+y     )*c2-y     *c1/3]
c               -y*(1+y)*c4/3}
c               +y*(1+y)*(2-y)*c3}/2
                  do nn=1,4,3   ! N.B.
                     c2 = sx(idel+nn-2,jdel  ,n,k)
                     c3 = sx(idel+nn-2,jdel+1,n,k)
                     r(nn) = (1.-yyg)*c2 +yyg*c3
                  enddo         ! nn loop
                  if(mhint==2)then ! Bessel interp
                     a4 = r(4)-r(1)+3.*(r(2)-r(3))
                     a3 = r(1)-2.*r(2)+r(3)-a4
                     sextra(iproc)%a(iq) = r(2)+                        ! MJT memory
     &                    0.5*xxg*(r(3)-r(1) +xxg*(a3+xxg*a4))
                  else
                     sextra(iproc)%a(iq) = ((1.-xxg)*((2.-xxg)*         ! MJT memory
     &                    ((1.+xxg)*r(2)-xxg*r(1)/3.)
     &                    -xxg*(1.+xxg)*r(4)/3.)
     &                    +xxg*(1.+xxg)*(2.-xxg)*r(3))/2.
                  endif         !  (mhint==2)
               enddo            ! iq loop
            end do              ! iproc
         else                   ! (nfield<mh_bs)
            do iproc=0,nproc-1
               if ( iproc == myid ) then
                  cycle
               end if
               do iq=1,drlen(iproc)
                  !  Convert face index from 0:npanels to array indices
                  ip = min(il_g,max(1,nint(dpoints(iproc)%a(2,iq))))    ! MJT memory
                  jp = min(il_g,max(1,nint(dpoints(iproc)%a(3,iq))))    ! MJT memory
                  n = nint(dpoints(iproc)%a(1,iq)) + noff ! Local index ! MJT memory
                  !  Need global face index in fproc call
#ifdef debug
                  if ( fproc(ip,jp,n-noff) /= myid ) then
                     print*, "Error in ints_bl", myid, n, iq, iproc,
     &                 dpoints(iproc)%a(:,iq)                           ! MJT memory
                     call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
                  end if
#endif
                  idel = int(dpoints(iproc)%a(2,iq))                    ! MJT memory
                  xxg = dpoints(iproc)%a(2,iq) - idel                   ! MJT memory
                  jdel = int(dpoints(iproc)%a(3,iq))                    ! MJT memory
                  yyg = dpoints(iproc)%a(3,iq) - jdel                   ! MJT memory
                  k = nint(dpoints(iproc)%a(4,iq))                      ! MJT memory
                  idel = idel - ioff(n-noff)
                  jdel = jdel - joff(n-noff)
                  c1 = sx(idel,jdel-1,n,k) ! manually unrolled loop
                  c2 = sx(idel,jdel  ,n,k)
                  c3 = sx(idel,jdel+1,n,k)
                  c4 = sx(idel,jdel+2,n,k)
                  cmin = min( 1.e20,c2,c3) ! Bermejo & Staniforth
                  cmax = max(-1.e20,c2,c3) ! Bermejo & Staniforth
                  if(mhint==2)then ! Bessel interp
                     a4 = c4-c1+3.*(c2-c3)
                     a3 = c1-2.*c2+c3-a4
                     r(2) = c2+.5*yyg*(c3-c1 +yyg*(a3+yyg*a4))
                  else
                     r(2) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*c2-yyg*c1/3.)
     &                    -yyg*(1.+yyg)*c4/3.)
     &                    +yyg*(1.+yyg)*(2.-yyg)*c3)/2.
                  endif         !  (mhint==2)
                  c1 = sx(idel+1,jdel-1,n,k)
                  c2 = sx(idel+1,jdel  ,n,k)
                  c3 = sx(idel+1,jdel+1,n,k)
                  c4 = sx(idel+1,jdel+2,n,k)
                  cmin = min(cmin,c2,c3) ! Bermejo & Staniforth
                  cmax = max(cmax,c2,c3) ! Bermejo & Staniforth
                  if(mhint==2)then ! Bessel interp
                     a4 = c4-c1+3.*(c2-c3)
                     a3 = c1-2.*c2+c3-a4
                     r(3) = c2+.5*yyg*(c3-c1 +yyg*(a3+yyg*a4))
                  else
                     r(3) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*c2-yyg*c1/3.)
     &                    -yyg*(1.+yyg)*c4/3.)
     &                    +yyg*(1.+yyg)*(2.-yyg)*c3)/2.
                  endif         !  (mhint==2)
c           r = {(1-y)*{(2-y     )*[(1+y     )*c2-y     *c1/3]
c                -y*(1+y)*c4/3}
c                +y*(1+y)*(2-y)*c3}/2
                  do nn=1,4,3   ! N.B.
                     c2 = sx(idel+nn-2,jdel  ,n,k)
                     c3 = sx(idel+nn-2,jdel+1,n,k)
                     r(nn) = (1.-yyg)*c2 +yyg*c3
                  enddo         ! nn loop
                  if(mhint==2)then ! Bessel interp
                     a4 = r(4)-r(1)+3.*(r(2)-r(3))
                     a3 = r(1)-2.*r(2)+r(3)-a4
                     sss = r(2)+.5*xxg*(r(3)-r(1) +xxg*(a3+xxg*a4))
                  else
                     sss = ((1.-xxg)*((2.-xxg)*
     &                    ((1.+xxg)*r(2)-xxg*r(1)/3.)
     &                    -xxg*(1.+xxg)*r(4)/3.)
     &                    +xxg*(1.+xxg)*(2.-xxg)*r(3))/2.
                  endif         !  (mhint==2)
                  sextra(iproc)%a(iq) = min(max(cmin,sss),cmax) ! Bermejo & Staniforth ! MJT memory
               enddo            ! iq loop
            end do              ! iproc
         endif                  ! (nfield<mh_bs)  .. else ..

      endif                     ! (intsch==1) .. else ..
!========================   end of intsch=1 section ====================

      call intssync(s)
      call end_log(ints_end)
      return

      entry ints_bl(s,intsch,nface,xg,yg)  ! not usually called

c     this one does bi-linear interpolation only
c     first extend s arrays into sx - this one -1:il+2 & -1:il+2
c                    but for bi-linear only need 0:il+1 &  0:il+1
      call start_log(ints_begin)
      call bounds(s,corner=.true.)
      do k=1,kl
         do n=1,npan
            do j=1,jpan
               do i=1,ipan
                  sx(i,j,n,k) = s(ind(i,j,n),k)
               enddo            ! i loop
            enddo               ! j loop
            do j=1,jpan
               sx(0,j,n,k) = s(iw(ind(1,j,n)),k)
               sx(ipan+1,j,n,k) = s(ie(ind(ipan,j,n)),k)
            enddo               ! j loop
            do i=1,ipan
               sx(i,0,n,k) = s(is(ind(i,1,n)),k)
               sx(i,jpan+1,n,k) = s(in(ind(i,jpan,n)),k)
            enddo               ! i loop

            sx(0,0,n,k) = s(iws(ind(1,1,n)),k)
            sx(ipan+1,0,n,k) = s(ies(ind(ipan,1,n)),k)
            sx(0,jpan+1,n,k)    = s(iwn(ind(1,jpan,n)),k)
            sx(ipan+1,jpan+1,n,k) = s(ien(ind(ipan,jpan,n)),k)
         enddo                  ! n loop

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
#ifdef SX
            ! For better vectorisation calculate all values
            idel = max(0,min(idel,ipan))
            jdel = max(0,min(jdel,jpan))
            n = max(1,min(n,npan))
#else
            if ( idel < 0 .or. idel > ipan .or. jdel < 0 .or.
     &           jdel > jpan .or. n < 1 .or. n > npan ) then
               cycle            ! Will be calculated on another processor
            end if
#endif
            s(iq,k) =      yyg*(      xxg*sx(idel+1,jdel+1,n,k)
     &                 +(1.-xxg)*sx(idel,jdel+1,n,k))
     &      +(1.-yyg)*(      xxg*sx(idel+1,jdel,n,k)
     &                 +(1.-xxg)*sx(idel,jdel,n,k))
         enddo                  ! iq loop
      end do                    ! k

! Loop over points that need to be calculated for other processes
      do iproc=0,nproc-1
         if ( iproc == myid ) then
            cycle
         end if
         do iq=1,drlen(iproc)
           !  Convert face index from 0:npanels to array indices
            ip = min(il_g,max(1,nint(dpoints(iproc)%a(2,iq))))          ! MJT memory
            jp = min(il_g,max(1,nint(dpoints(iproc)%a(3,iq))))          ! MJT memory
            n = nint(dpoints(iproc)%a(1,iq)) + noff ! Local index       ! MJT memory
         !  Need global face index in fproc call
#ifdef debug
            if ( fproc(ip,jp,n-noff) /= myid ) then
               print*, "Error in ints_bl", myid, n, iq, iproc,
     &              dpoints(iproc)%a(:,iq)                              ! MJT memory
               call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
            end if
#endif
            idel = int(dpoints(iproc)%a(2,iq))                          ! MJT memory
            xxg = dpoints(iproc)%a(2,iq) - idel                         ! MJT memory
            jdel = int(dpoints(iproc)%a(3,iq))                          ! MJT memory
            yyg = dpoints(iproc)%a(3,iq) - jdel                         ! MJT memory
            k = nint(dpoints(iproc)%a(4,iq))                            ! MJT memory
            idel = idel - ioff(n-noff)
            jdel = jdel - joff(n-noff)
            sextra(iproc)%a(iq) = yyg*( xxg*sx(idel+1,jdel+1,n,k)       ! MJT memory
     &                              +(1.-xxg)*sx(idel,jdel+1,n,k))
     &                   +(1.-yyg)*(      xxg*sx(idel+1,jdel,n,k)
     &                               +(1.-xxg)*sx(idel,jdel,n,k))
         end do
      end do

      call intssync(s)

      call end_log(ints_end)
      return
      end
