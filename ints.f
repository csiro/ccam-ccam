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
      implicit none
      integer, parameter :: ntest=0
      include 'newmpar.h'
      include 'parm.h'
      include 'parmhor.h'    ! has mh_bs
      integer intsch, nfield
      integer nface(ifull,kl)
      real xg(ifull,kl),yg(ifull,kl)      ! now passed through call
      ! Need work common for this
      real sx(-1:il+2,-1:il+2,0:npanels,kl)
      real s(ifull,kl),r(4)
      include 'indices.h' ! in,is,iw,ie,inn,iss,iww,iee
      integer   in6(il,il,0:npanels),is6(il,il,0:npanels)
     &         ,iw6(il,il,0:npanels),ie6(il,il,0:npanels)
     &         ,inn6(il,il,0:npanels),iss6(il,il,0:npanels)
     &         ,iww6(il,il,0:npanels),iee6(il,il,0:npanels)
     &         ,ine6(il,il,0:npanels),ise6(il,il,0:npanels)
     &         ,ien6(il,il,0:npanels),iwn6(il,il,0:npanels)
      equivalence (in,in6),(is,is6),(iw,iw6),(ie,ie6)
      equivalence (inn,inn6),(iss,iss6),(iww,iww6),(iee,iee6)
      equivalence (ine,ine6),(ise,ise6),(ien,ien6),(iwn,iwn6)
      integer idel, iq, ix, jdel, nn
      real a3, a4, c1, c2, c3, c4, cmax, cmin, sss, xxg, yyg
      integer i, j, k, n, ind
      ind(i,j,n)=i+(j-1)*il+n*il*il  ! *** for n=0,npanels

      do k=1,kl
         do n=0,npanels         ! first simple copy into larger array
            do j=1,il
               do i=1,il
                  sx(i,j,n,k)=s(ind(i,j,n),k)
               enddo            ! i loop
            enddo               ! j loop
         enddo                  ! n loop
      end do

      if(intsch.eq.1)then
c       this is intsb           EW interps done first
c       first extend s arrays into sx - this one -1:il+2 & -1:il+2
         do k=1,kl
            do n=0,npanels
               do j=1,il
                  sx(0,j,n,k)=s(iw6(1,j,n),k)
                  sx(-1,j,n,k)=s(iww6(1,j,n),k)
                  sx(il+1,j,n,k)=s(ie6(il,j,n),k)
                  sx(il+2,j,n,k)=s(iee6(il,j,n),k)
               enddo            ! j loop
               do i=1,il
                  sx(i,0,n,k)=s(is6(i,1,n),k)
                  sx(i,-1,n,k)=s(iss6(i,1,n),k)
                  sx(i,il+1,n,k)=s(in6(i,il,n),k)
                  sx(i,il+2,n,k)=s(inn6(i,il,n),k)
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

               sx(-1,0,n,k)=s(lwws(n),k)
               sx(0,0,n,k)=s(lws(n),k)
               sx(0,-1,n,k)=s(lwss(n),k)
               sx(il+1,0,n,k)=s(les(n),k)
               sx(il+2,0,n,k)=s(lees(n),k)
               sx(il+1,-1,n,k)=s(less(n),k)
               sx(-1,il+1,n,k)=s(lwwn(n),k)
               sx(0,il+2,n,k)=s(lwnn(n),k)
               sx(il+2,il+1,n,k)=s(leen(n),k)
               sx(il+1,il+2,n,k)=s(lenn(n),k)
               sx(0,il+1,n,k)   =s(iwn6(1,il,n),k)
               sx(il+1,il+1,n,k)=s(ien6(il,il,n),k)
            enddo               ! n loop

        if(nfield.lt.mh_bs)then
          do iq=1,ifull           ! non Berm-Stan option
           n=nface(iq,k)
           idel=int(xg(iq,k))
           xxg=xg(iq,k)-idel
           jdel=int(yg(iq,k))
           yyg=yg(iq,k)-jdel
c          if(ntest.eq.1)then
c            if(xxg.lt.0..or.xxg.gt.1..or.yyg.lt.0..or.yyg.gt.1.)
c    &        print *,'intsb problem iq xxg,yyg ',iq,xxg,yyg
c          endif
            c1=sx(idel-1,jdel,n,k)   ! manually unrolled loop
            c2=sx(idel  ,jdel,n,k)
            c3=sx(idel+1,jdel,n,k)
            c4=sx(idel+2,jdel,n,k)
            if(mhint.eq.2)then   ! Bessel interp
              a4=c4-c1+3.*(c2-c3)
              a3=c1-2.*c2+c3-a4
              r(2)=c2+.5*xxg*(c3-c1 +xxg*(a3+xxg*a4))
            else
              r(2)=((1.-xxg)*((2.-xxg)*((1.+xxg)*c2-xxg*c1/3.)
     &          -xxg*(1.+xxg)*c4/3.)
     &          +xxg*(1.+xxg)*(2.-xxg)*c3)/2.
            endif   !  (mhint.eq.2)
            c1=sx(idel-1,jdel+1,n,k)
            c2=sx(idel  ,jdel+1,n,k)
            c3=sx(idel+1,jdel+1,n,k)
            c4=sx(idel+2,jdel+1,n,k)
            if(mhint.eq.2)then   ! Bessel interp
              a4=c4-c1+3.*(c2-c3)
              a3=c1-2.*c2+c3-a4
              r(3)=c2+.5*xxg*(c3-c1 +xxg*(a3+xxg*a4))
            else
              r(3)=((1.-xxg)*((2.-xxg)*((1.+xxg)*c2-xxg*c1/3.)
     &          -xxg*(1.+xxg)*c4/3.)
     &          +xxg*(1.+xxg)*(2.-xxg)*c3)/2.
            endif   !  (mhint.eq.2)
c           r = {(1-x)*{(2-x     )*[(1+x     )*c2-x     *c1/3]
c                -x*(1+x)*c4/3}
c                +x*(1+x)*(2-x)*c3}/2
           do nn=1,4,3       ! N.B.
            c2=sx(idel  ,jdel+nn-2,n,k)
            c3=sx(idel+1,jdel+nn-2,n,k)
            r(nn)=(1.-xxg)*c2 +xxg*c3
           enddo    ! nn loop
           if(mhint.eq.2)then   ! Bessel interp
             a4=r(4)-r(1)+3.*(r(2)-r(3))
             a3=r(1)-2.*r(2)+r(3)-a4
             s(iq,k)=r(2)+.5*yyg*(r(3)-r(1) +yyg*(a3+yyg*a4))
           else
             s(iq,k)=((1.-yyg)*((2.-yyg)*((1.+yyg)*r(2)-yyg*r(1)/3.)
     &        -yyg*(1.+yyg)*r(4)/3.)
     &        +yyg*(1.+yyg)*(2.-yyg)*r(3))/2.
            endif   !  (mhint.eq.2)
          enddo    ! iq loop
        else                      ! (nfield.lt.mh_bs)
          do iq=1,ifull           ! Berm-Stan option here e.g. qg & gases
           n=nface(iq,k)
           idel=int(xg(iq,k))
           xxg=xg(iq,k)-idel
           jdel=int(yg(iq,k))
           yyg=yg(iq,k)-jdel
            c1=sx(idel-1,jdel,n,k)     ! manually unrolled loop
            c2=sx(idel  ,jdel,n,k)
            c3=sx(idel+1,jdel,n,k)
            c4=sx(idel+2,jdel,n,k)
            cmin=min( 1.e20,c2,c3)   ! Bermejo & Staniforth
            cmax=max(-1.e20,c2,c3)   ! Bermejo & Staniforth
            if(mhint.eq.2)then   ! Bessel interp
              a4=c4-c1+3.*(c2-c3)
              a3=c1-2.*c2+c3-a4
              r(2)=c2+.5*xxg*(c3-c1 +xxg*(a3+xxg*a4))
            else
              r(2)=((1.-xxg)*((2.-xxg)*((1.+xxg)*c2-xxg*c1/3.)
     &          -xxg*(1.+xxg)*c4/3.)
     &          +xxg*(1.+xxg)*(2.-xxg)*c3)/2.
            endif   !  (mhint.eq.2)
            c1=sx(idel-1,jdel+1,n,k)
            c2=sx(idel  ,jdel+1,n,k)
            c3=sx(idel+1,jdel+1,n,k)
            c4=sx(idel+2,jdel+1,n,k)
            cmin=min(cmin,c2,c3)   ! Bermejo & Staniforth
            cmax=max(cmax,c2,c3)   ! Bermejo & Staniforth
            if(mhint.eq.2)then   ! Bessel interp
              a4=c4-c1+3.*(c2-c3)
              a3=c1-2.*c2+c3-a4
              r(3)=c2+.5*xxg*(c3-c1 +xxg*(a3+xxg*a4))
            else
              r(3)=((1.-xxg)*((2.-xxg)*((1.+xxg)*c2-xxg*c1/3.)
     &          -xxg*(1.+xxg)*c4/3.)
     &          +xxg*(1.+xxg)*(2.-xxg)*c3)/2.
            endif   !  (mhint.eq.2)
c           r = {(1-x)*{(2-x     )*[(1+x     )*c2-x     *c1/3]
c                -x*(1+x)*c4/3}
c                +x*(1+x)*(2-x)*c3}/2
           do nn=1,4,3       ! N.B.
            c2=sx(idel  ,jdel+nn-2,n,k)
            c3=sx(idel+1,jdel+nn-2,n,k)
            r(nn)=(1.-xxg)*c2 +xxg*c3
           enddo    ! nn loop
           if(mhint.eq.2)then   ! Bessel interp
             a4=r(4)-r(1)+3.*(r(2)-r(3))
             a3=r(1)-2.*r(2)+r(3)-a4
             sss=r(2)+.5*yyg*(r(3)-r(1) +yyg*(a3+yyg*a4))
           else
             sss=((1.-yyg)*((2.-yyg)*((1.+yyg)*r(2)-yyg*r(1)/3.)
     &        -yyg*(1.+yyg)*r(4)/3.)
     &        +yyg*(1.+yyg)*(2.-yyg)*r(3))/2.
          endif   !  (mhint.eq.2)
          s(iq,k)=min(max(cmin,sss),cmax)   ! Bermejo & Staniforth
          enddo    ! iq loop
        endif      ! (nfield.lt.mh_bs)  .. else ..

      end do ! k

      else     ! if(intsch.eq.1)then
c       this is intsc           NS interps done first
c       first extend s arrays into sx - this one -1:il+2 & -1:il+2
         do k=1,kl
            do n=0,npanels
               do j=1,il
                  sx(0,j,n,k)=s(iw6(1,j,n),k)
                  sx(-1,j,n,k)=s(iww6(1,j,n),k)
                  sx(il+1,j,n,k)=s(ie6(il,j,n),k)
                  sx(il+2,j,n,k)=s(iee6(il,j,n),k)
               enddo            ! j loop
               do i=1,il
                  sx(i,0,n,k)=s(is6(i,1,n),k)
                  sx(i,-1,n,k)=s(iss6(i,1,n),k)
                  sx(i,il+1,n,k)=s(in6(i,il,n),k)
                  sx(i,il+2,n,k)=s(inn6(i,il,n),k)
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
               sx(0,0,n,k)=s(lsw(n),k)
               sx(0,-1,n,k)=s(lssw(n),k)
               sx(il+2,0,n,k)=s(lsee(n),k)
               sx(il+1,-1,n,k)=s(lsse(n),k)
               sx(-1,il+1,n,k)=s(lnww(n),k)
               sx(0,il+1,n,k)=s(lnw(n),k)
               sx(0,il+2,n,k)=s(lnnw(n),k)
               sx(il+2,il+1,n,k)=s(lnee(n),k)
               sx(il+1,il+2,n,k)=s(lnne(n),k)
               sx(il+1,0,n,k)   =s(ise6(il,1,n),k)
               sx(il+1,il+1,n,k)=s(ine6(il,il,n),k)
            enddo               ! n loop

        if(nfield.lt.mh_bs)then
          do iq=1,ifull           ! non Berm-Stan option
           n=nface(iq,k)
           idel=int(xg(iq,k))
           xxg=xg(iq,k)-idel
           jdel=int(yg(iq,k))
           yyg=yg(iq,k)-jdel
c          if(ntest.eq.1)then
c            if(xxg.lt.0..or.xxg.gt.1..or.yyg.lt.0..or.yyg.gt.1.)
c    &        print *,'intsc problem iq xxg,yyg ',iq,xxg,yyg
c          endif
            c1=sx(idel,jdel-1,n,k)       ! manually unrolled loop
            c2=sx(idel,jdel  ,n,k)
            c3=sx(idel,jdel+1,n,k)
            c4=sx(idel,jdel+2,n,k)
            if(mhint.eq.2)then   ! Bessel interp
              a4=c4-c1+3.*(c2-c3)
              a3=c1-2.*c2+c3-a4
              r(2)=c2+.5*yyg*(c3-c1 +yyg*(a3+yyg*a4))
            else
              r(2)=((1.-yyg)*((2.-yyg)*((1.+yyg)*c2-yyg*c1/3.)
     &          -yyg*(1.+yyg)*c4/3.)
     &          +yyg*(1.+yyg)*(2.-yyg)*c3)/2.
            endif   !  (mhint.eq.2)
            c1=sx(idel+1,jdel-1,n,k)
            c2=sx(idel+1,jdel  ,n,k)
            c3=sx(idel+1,jdel+1,n,k)
            c4=sx(idel+1,jdel+2,n,k)
            if(mhint.eq.2)then   ! Bessel interp
              a4=c4-c1+3.*(c2-c3)
              a3=c1-2.*c2+c3-a4
              r(3)=c2+.5*yyg*(c3-c1 +yyg*(a3+yyg*a4))
            else
              r(3)=((1.-yyg)*((2.-yyg)*((1.+yyg)*c2-yyg*c1/3.)
     &          -yyg*(1.+yyg)*c4/3.)
     &          +yyg*(1.+yyg)*(2.-yyg)*c3)/2.
            endif   !  (mhint.eq.2)
c          r = {(1-y)*{(2-y     )*[(1+y     )*c2-y     *c1/3]
c               -y*(1+y)*c4/3}
c               +y*(1+y)*(2-y)*c3}/2
           do nn=1,4,3       ! N.B.
            c2=sx(idel+nn-2,jdel  ,n,k)
            c3=sx(idel+nn-2,jdel+1,n,k)
            r(nn)=(1.-yyg)*c2 +yyg*c3
           enddo    ! nn loop
           if(mhint.eq.2)then   ! Bessel interp
             a4=r(4)-r(1)+3.*(r(2)-r(3))
             a3=r(1)-2.*r(2)+r(3)-a4
             s(iq,k)=r(2)+.5*xxg*(r(3)-r(1) +xxg*(a3+xxg*a4))
           else
             s(iq,k)=((1.-xxg)*((2.-xxg)*((1.+xxg)*r(2)-xxg*r(1)/3.)
     &        -xxg*(1.+xxg)*r(4)/3.)
     &        +xxg*(1.+xxg)*(2.-xxg)*r(3))/2.
           endif   !  (mhint.eq.2)
          enddo    ! iq loop
        else                      ! (nfield.lt.mh_bs)
          do iq=1,ifull           ! Berm-Stan option here e.g. qg & gases
           n=nface(iq,k)
           idel=int(xg(iq,k))
           xxg=xg(iq,k)-idel
           jdel=int(yg(iq,k))
           yyg=yg(iq,k)-jdel
            c1=sx(idel,jdel-1,n,k)       ! manually unrolled loop
            c2=sx(idel,jdel  ,n,k)
            c3=sx(idel,jdel+1,n,k)
            c4=sx(idel,jdel+2,n,k)
            cmin=min( 1.e20,c2,c3)   ! Bermejo & Staniforth
            cmax=max(-1.e20,c2,c3)   ! Bermejo & Staniforth
            if(mhint.eq.2)then   ! Bessel interp
              a4=c4-c1+3.*(c2-c3)
              a3=c1-2.*c2+c3-a4
              r(2)=c2+.5*yyg*(c3-c1 +yyg*(a3+yyg*a4))
            else
              r(2)=((1.-yyg)*((2.-yyg)*((1.+yyg)*c2-yyg*c1/3.)
     &          -yyg*(1.+yyg)*c4/3.)
     &          +yyg*(1.+yyg)*(2.-yyg)*c3)/2.
            endif   !  (mhint.eq.2)
            c1=sx(idel+1,jdel-1,n,k)
            c2=sx(idel+1,jdel  ,n,k)
            c3=sx(idel+1,jdel+1,n,k)
            c4=sx(idel+1,jdel+2,n,k)
            cmin=min(cmin,c2,c3)   ! Bermejo & Staniforth
            cmax=max(cmax,c2,c3)   ! Bermejo & Staniforth
            if(mhint.eq.2)then   ! Bessel interp
              a4=c4-c1+3.*(c2-c3)
              a3=c1-2.*c2+c3-a4
              r(3)=c2+.5*yyg*(c3-c1 +yyg*(a3+yyg*a4))
            else
              r(3)=((1.-yyg)*((2.-yyg)*((1.+yyg)*c2-yyg*c1/3.)
     &          -yyg*(1.+yyg)*c4/3.)
     &          +yyg*(1.+yyg)*(2.-yyg)*c3)/2.
            endif   !  (mhint.eq.2)
c           r = {(1-y)*{(2-y     )*[(1+y     )*c2-y     *c1/3]
c                -y*(1+y)*c4/3}
c                +y*(1+y)*(2-y)*c3}/2
           do nn=1,4,3       ! N.B.
            c2=sx(idel+nn-2,jdel  ,n,k)
            c3=sx(idel+nn-2,jdel+1,n,k)
            r(nn)=(1.-yyg)*c2 +yyg*c3
           enddo    ! nn loop
           if(mhint.eq.2)then   ! Bessel interp
             a4=r(4)-r(1)+3.*(r(2)-r(3))
             a3=r(1)-2.*r(2)+r(3)-a4
             sss=r(2)+.5*xxg*(r(3)-r(1) +xxg*(a3+xxg*a4))
           else
             sss=((1.-xxg)*((2.-xxg)*((1.+xxg)*r(2)-xxg*r(1)/3.)
     &        -xxg*(1.+xxg)*r(4)/3.)
     &        +xxg*(1.+xxg)*(2.-xxg)*r(3))/2.
           endif   !  (mhint.eq.2)
           s(iq,k)=min(max(cmin,sss),cmax)   ! Bermejo & Staniforth
          enddo    ! iq loop
        endif      ! (nfield.lt.mh_bs)  .. else ..
      end do ! k
      endif    ! (intsch.eq.1) .. else ..

      return

      entry ints_bl(s,intsch,nface,xg,yg)  ! not usually called

c     this one does bi-linear interpolation only
c     first extend s arrays into sx - this one -1:il+2 & -1:il+2
c                    but for bi-linear only need 0:il+1 &  0:il+1
      do k=1,kl
         do n=0,npanels
            do j=1,il
               do i=1,il
                  sx(i,j,n,k)=s(ind(i,j,n),k)
               enddo            ! i loop
            enddo               ! j loop
            do j=1,il
               sx(0,j,n,k)=s(iw6(1,j,n),k)
c       sx(-1,j,n)=s(iww6(1,j,n))
               sx(il+1,j,n,k)=s(ie6(il,j,n),k)
c       sx(il+2,j,n)=s(iee6(il,j,n))
            enddo               ! j loop
            do i=1,il
               sx(i,0,n,k)=s(is6(i,1,n),k)
c       sx(i,-1,n)=s(iss6(i,1,n))
               sx(i,il+1,n,k)=s(in6(i,il,n),k)
c       sx(i,il+2,n)=s(inn6(i,il,n))
            enddo               ! i loop

c      sx(-1,0,n)=s(lwws(n))
            sx(0,0,n,k)=s(lws(n),k)
c      sx(0,-1,n)=s(lwss(n))
            sx(il+1,0,n,k)=s(les(n),k)
c      sx(il+2,0,n)=s(lees(n))
c      sx(il+1,-1,n)=s(less(n))
c      sx(-1,il+1,n)=s(lwwn(n))
c      sx(0,il+2,n)=s(lwnn(n))
c      sx(il+2,il+1,n)=s(leen(n))
c      sx(il+1,il+2,n)=s(lenn(n))
            sx(0,il+1,n,k)   =s(iwn6(1,il,n),k)
            sx(il+1,il+1,n,k)=s(ien6(il,il,n),k)
         enddo                  ! n loop

      do iq=1,ifull
       n=nface(iq,k)
       idel=int(xg(iq,k))
       xxg=xg(iq,k)-idel
       jdel=int(yg(iq,k))
       yyg=yg(iq,k)-jdel
       s(iq,k)=     yyg*(      xxg*sx(idel+1,jdel+1,n,k)
     &                 +(1.-xxg)*sx(idel,jdel+1,n,k))
     &      +(1.-yyg)*(      xxg*sx(idel+1,jdel,n,k)
     &                 +(1.-xxg)*sx(idel,jdel,n,k))
      enddo    ! iq loop
      end do ! k
      return
      end
