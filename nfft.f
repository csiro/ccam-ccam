c    nfft.f contains: rfft, hfft, xrfft, xhfft, cfft, dfft, rumble, xcfft,
c                     xdfft, ycfft, ydfft, blockdata dat235
c                                               *****************
c                                               *   nfft2.for   *
c                                               *  purser 1994  *
c                                               *****************
c
c------------------------------------------------------------------------------c
c   r.j.purser, national meteorological center, washington d.c.  1994          c
c                   subroutines rfft, hfft                                     c
c                                                                              c
c   fourier analyze (rfft) or synthesize (hfft) a line of real data            c
c <-> r       data and transform                                               c
c --> n       number of complex data of series (product of 2's,3's,5's)        c
c --> period  period of one cycle                                              c
c <-> w       coefficient array of n reals (re)computed only when              c
c             jumble(0) differs from n                                         c
c <-> jumble  coefficient array of n integers (re)computed only when           c
c             jumble(0) set to zero on input to this routine                   c
c --- t       work array of size n                                             c
c------------------------------------------------------------------------------c
      subroutine rfft(r,n,period,w,jumble,t)
      dimension r(0:*),t(0:*),w(0:*),jumble(0:*)
      nh=n/2
      nm=n-1
      do i=0,nm
       t(i)=0.
      enddo
      call cfft(r,t,n,period,w,jumble)
      do i=nh+1,nm
       r(i)=t(i)
      enddo
      return

      entry     hfft(r,n,period,w,jumble,t)
      nh=n/2
      periodi=n/period
      t(0)=0
      t(nh)=0
      do i=nh+1,n-1
       t(n-i)=r(i)
       t(i)=-r(i)
       r(i)=r(n-i)
      enddo
      call cfft(r,t,n,periodi,w,jumble)
      return
      end

c------------------------------------------------------------------------------c
c   r.j.purser, national meteorological center, washington d.c.  1994          c
c                   subroutine xrfft                                           c
c                                                                              c
c   fourier analyze several lines of real data                                 c
c <-> rb      data and transform                                               c
c --> n       number of data in each cycle (must be product of 2's,3's,5's)    c
c --> ny      number of lines of data in the array                             c
c --> ncb     first fortran dimension of array of data (at least as big as n)  c
c --> period  period of one cycle                                              c
c <-> w       coefficient array of n reals (re)computed only when              c
c             jumble(0) differs from n                                         c
c <-> jumble  coefficient array of n integers (re)computed only when           c
c             jumble(0) set to zero on input to this routine                   c
c --- t1      work array. of size n when ny is odd, not used when ny even      c
c------------------------------------------------------------------------------c
      subroutine xrfft(rb,n,ny,ncb,period,w,jumble,t1)
      logical even
      dimension rb(0:ncb-1,*),t1(0:*),w(0:*),jumble(0:*)
      ncb2=ncb*2
      nyh=ny/2
      nyh2=nyh*2
      nm=n-1
      nh=n/2
      even=nyh2.eq.ny
      periodh=period/2
      call xcfft(rb(0,1),rb(0,2),n,nyh,ncb2,periodh,w,jumble)
      do iy=1,nyh2
        rb(0,iy)=rb(0,iy)*2
      enddo
      if(nh*2.eq.n)then
      do iy=1,nyh2
        rb(nh,iy)=rb(nh,iy)*2
      enddo
      endif
      do j=nh+1,nm
       i=n-j
       do iy=1,nyh2,2
        t         =rb(i,iy)  -rb(j,iy)
        rb(i,iy)  =rb(i,iy)  +rb(j,iy)
        rb(j,iy) =-rb(i,iy+1)+rb(j,iy+1)
        rb(i,iy+1)=rb(i,iy+1)+rb(j,iy+1)
        rb(j,iy+1)=t
       enddo
      enddo
      if(.not.even)then
       do i=0,nm
        t1(i)=0.
       enddo
       call cfft(rb(0,ny),t1,n,period,w,jumble)
       do i=nh+1,nm
        rb(i,ny)=t1(i)
       enddo
      endif
      return
      end

c------------------------------------------------------------------------------c
c   r.j.purser, national meteorological center, washington d.c.  1994          c
c                   subroutine xhfft                                           c
c                                                                              c
c   fourier synthesize several lines of real data (inverse of xrfft)           c
c <-> rb      data and transform                                               c
c --> n       number of data in each cycle (must be product of 2's,3's,5's)    c
c --> ny      number of lines of data in the array                             c
c --> ncb     first fortran dimension of array of data (at least as big as n)  c
c --> period  period of one cycle                                              c
c <-> w       coefficient array of n reals (re)computed only when              c
c             jumble(0) differs from n                                         c
c <-> jumble  coefficient array of n integers (re)computed only when           c
c             jumble(0) set to zero on input to this routine                   c
c --- t1      work array. of size n when ny is odd, not used when ny even      c
c------------------------------------------------------------------------------c
      subroutine xhfft(rb,n,ny,ncb,period,w,jumble,t1)
      logical even
      dimension rb(0:ncb-1,*),t1(0:*),w(0:*),jumble(0:*)
      ncb2=ncb*2
      nyh=ny/2
      nyh2=nyh*2
      nm=n-1
      nh=n/2
      even=nyh2.eq.ny
      periodi=n/period
      do j=nh+1,nm
       i=n-j
       do iy=1,nyh2,2
        t         =rb(j,iy)
        rb(j,iy)  =rb(i,iy)  +rb(j,iy+1)
        rb(i,iy)  =rb(i,iy)  -rb(j,iy+1)
        rb(j,iy+1)=rb(i,iy+1)-t
        rb(i,iy+1)=rb(i,iy+1)+t
       enddo
      enddo
      call xcfft(rb(0,1),rb(0,2),n,nyh,ncb2,periodi,w,jumble)
      if(.not.even)then
       t1(0)=0.
       t1(nh)=0.
       do j=nh+1,nm
        t1(n-j)=rb(j,ny)
        t1(j)  =-rb(j,ny)
        rb(j,ny)=rb(n-j,ny)
       enddo
       call cfft(rb(0,ny),t1,n,periodi,w,jumble)
      endif
      return
      end
c                                               *****************
c                                               *   nfft1.for   *
c                                               *  purser 1994  *
c                                               *****************
c
      blockdata dat235
      common/fft235/ ln2,ln3,ln5
      save/fft235/
      data ln2/0/,ln3/0/,ln5/0/
      end

c------------------------------------------------------------------------------c
c   r.j.purser, national meteorological center, washington d.c.  1994          c
c                   subroutines cfft, dfft                                     c
c                                                                              c
c   fourier analyze (cfft) or synthesize (dfft) a line of complex data         c
c <-> rb      real part of data and transform                                  c
c <-> qb      imaginary part of data and transform                             c
c --> n       number of complex data of series (product of 2's,3's,5's)        c
c --> period  period of one cycle                                              c
c <-> w       coefficient array of n reals (re)computed only when              c
c             jumble(0)=0                                                      c
c <-> jumble  coefficient array of n integers (re)computed only when           c
c             jumble(0) set to zero on input to this routine                   c
c                                                                              c
c   for related routines, see subroutine rumble (in which the                  c
c   coefficient arrays are initialized), xcfft, ycfft (matrix versions         c
c   of cfft), xdfft, ydfft (matrix versions of dfft), rfft (fourier            c
c   analysis of real data), hfft (fourier synthesis of real data), and         c
c   matrix counterparts, xrfft, yrfft, xhfft, yhfft.                           c
c------------------------------------------------------------------------------c
      subroutine cfft(rb,qb,n,period,w,jumble)
      logical get235
      common/fft235/ ln2,ln3,ln5
      save/fft235/
      dimension rb(0:*),qb(0:*),w(0:*),jumble(0:*)
      rfac=period/n
      goto 300
      entry      dfft(rb,qb,n,period,w,jumble)
      rfac=1./period

c  for fourier synthesis, reverse the order of wavenumbers:
      do j=n/2+1,n-1
       i=n-j
       t    =rb(j)
       rb(j)=rb(i)
       rb(i)=t
       t    =qb(j)
       qb(j)=qb(i)
       qb(i)=t
      enddo

300   if(n.ne.2**ln2*3**ln3*5**ln5)then
      if(.not.get235(n))stop'prime factors of n are not only 2, 3, 5'
      endif
      if(n.ne.jumble(0))call rumble(jumble,w)
      nm=n-1
      nh=n/2

c  scale and permute the data:
      rb(0)=rb(0)*rfac
      qb(0)=qb(0)*rfac
      do i=1,nm
       j=jumble(i)
       if(j.gt.i)then
        t1=rb(i)
        rb(i)=rb(j)*rfac
        rb(j)=t1
        t1=qb(i)
        qb(i)=qb(j)*rfac
        qb(j)=t1
       else
        rb(i)=rb(i)*rfac
        qb(i)=qb(i)*rfac
       endif
      enddo

c  transform the data:
      ma=1
      mb=n
c  radix 4
      ls=0
      do l=2,ln2,2
      mb=mb/4
      ls=l
      ma4=ma*4
      mb2=mb*2
      do j=0,ma-1
      jmb=j*mb
      jmb2=j*mb2
      rf1=w(jmb)
      qf1=w(nh+jmb)
      rf2=w(jmb2)
      qf2=w(nh+jmb2)
      rf3=rf1*rf2-qf1*qf2
      qf3=rf1*qf2+qf1*rf2
      do i=0,nm,ma4
      k0=i+j
      k1=k0+ma
      k2=k1+ma
      k3=k2+ma
      r0=rb(k0)
      r1=rb(k1)
      r2=rb(k2)
      r3=rb(k3)
      q0=qb(k0)
      q1=qb(k1)
      q2=qb(k2)
      q3=qb(k3)
      t1=r3*rf3-q3*qf3  ! q13
      q3=r3*qf3+q3*rf3  ! r13
      r3=r2*rf1-q2*qf1  ! r12
      r2=r2*qf1+q2*rf1  ! q12
      q2=q3-r2          ! r23
      r2=q3+r2          ! q22
      q3=r3+t1          ! r22
      t1=r3-t1          ! q23
      r3=r1*rf2-q1*qf2  ! r11
      q1=r1*qf2+q1*rf2  ! q11
      r1=r0-r3          ! r21
      r0=r0+r3          ! r20
      rb(k3)=r1-q2      ! r3
      rb(k1)=r1+q2      ! r1
      q2=q0+q1          ! q20
      q1=q0-q1          ! q21
      qb(k0)=q2+r2      ! q0
      qb(k2)=q2-r2      ! q2
      rb(k2)=r0-q3      ! r2
      rb(k0)=r0+q3      ! r0
      qb(k3)=q1-t1      ! q3
      qb(k1)=q1+t1      ! q1
      enddo
      enddo
      ma=ma4
      enddo
      if(ls.ne.ln2)then
c  radix 2
      mb=mb/2
      ma2=ma*2
      do j=0,ma-1
      jmb=j*mb
      do i=0,nm,ma2
      k0=j+i
      k1=k0+ma
      rf1=w(jmb)
      qf1=w(nh+jmb)
      r0=rb(k0)
      q0=qb(k0)
      r1=rb(k1)
      q1=qb(k1)
      t1=r1*qf1+q1*rf1 ! q11
      q1=r1*rf1-q1*qf1 ! r11
      rb(k1)=r0-q1     ! r1
      rb(k0)=r0+q1     ! r0
      qb(k1)=q0-t1     ! q1
      qb(k0)=q0+t1     ! q0
      enddo
      enddo
      ma=ma2
      endif
c  radix 3
      rep=-.5
      rec=1.5
      qep=.5*sqrt(3.)
      do l=1,ln3
      mb=mb/3
      ma3=ma*3
      do j=0,ma-1
      jmb=j*mb
      rf1=w(jmb)
      qf1=w(nh+jmb)
      rf2=rf1*rf1-qf1*qf1
      qf2=2*rf1*qf1
      do i=0,nm,ma3
      k0=i+j
      k1=k0+ma
      k2=k1+ma
      r1=rb(k1)
      q1=qb(k1)
      r2=rb(k2)
      q2=qb(k2)
      t1=r2*qf2+q2*rf2  ! r12
      q2=r2*rf2-q2*qf2  ! q12
      r2=r1*qf1+q1*rf1  ! q11
      r1=r1*rf1-q1*qf1  ! r11
      q1=r2+t1          ! q21
      r2=(r2-t1)*qep    ! r22
      t1=r1+q2          ! r21
      r1=(r1-q2)*qep    ! q22
      rb(k0)=rb(k0)+t1  ! r0
      qb(k0)=qb(k0)+q1  ! q0
      t1=rb(k0)-t1*rec  ! r21
      q1=qb(k0)-q1*rec  ! q21
      qb(k2)=q1-r1      ! q2
      qb(k1)=q1+r1      ! q1
      rb(k1)=t1-r2      ! r1
      rb(k2)=t1+r2      ! r2
      enddo
      enddo
      ma=ma3
      enddo

      if(ln5.gt.0)then
c  radix 5
      nze=n/5
      rze=w(nze)
      qze=w(nh+nze)
      rzc=1.-rze
      ret=rze*rze-qze*qze
      qet=2*rze*qze
      rec=1.-ret
      do l=1,ln5
      mb=mb/5
      ma5=ma*5
      do j=0,ma-1
      jmb=j*mb
      jmb2=jmb*2
      rf1=w(jmb)
      qf1=w(nh+jmb)
      rf2=w(jmb2)
      qf2=w(nh+jmb2)
      rf3=rf1*rf2-qf1*qf2
      qf3=rf1*qf2+qf1*rf2
      rf4=rf2*rf2-qf2*qf2
      qf4=2*rf2*qf2
      do i=0,nm,ma5
      k0=i+j
      k1=k0+ma
      k2=k1+ma
      k3=k2+ma
      k4=k3+ma
      r1=rb(k1)
      r2=rb(k2)
      r3=rb(k3)
      r4=rb(k4)
      q1=qb(k1)
      q2=qb(k2)
      q3=qb(k3)
      q4=qb(k4)
      t1=r1*qf1+q1*rf1        ! q11
      r1=r1*rf1-q1*qf1        ! r11
      q1=r4*rf4-q4*qf4        ! q14
      r4=r4*qf4+q4*rf4        ! r14
      q4=r1-q1                ! q24
      r1=r1+q1                ! r21
      q1=t1+r4                ! q21
      r4=t1-r4                ! r24
      t1=r3*rf3-q3*qf3        ! q13
      r3=r3*qf3+q3*rf3        ! r13
      q3=r2*qf2+q2*rf2        ! q12
      r2=r2*rf2-q2*qf2        ! r12
      q2=q3+r3                ! q22
      r3=q3-r3                ! r23
      q3=r2-t1                ! q23
      r2=r2+t1                ! r22
      rb(k0)=rb(k0)+r1+r2     ! r0
      qb(k0)=qb(k0)+q1+q2     ! q0
      t1=       r4*qze+r3*qet ! r34
      r3=       r3*qze-r4*qet ! r33
      r4=rb(k0)-r2*rzc-r1*rec ! r32
      r1=rb(k0)-r1*rzc-r2*rec ! r31
      rb(k2)=r4+r3            ! r2
      rb(k3)=r4-r3            ! r3
      rb(k4)=r1+t1            ! r4
      rb(k1)=r1-t1            ! r1
      t1=qb(k0)-q1*rzc-q2*rec ! q31
      q2=qb(k0)-q2*rzc-q1*rec ! q32
      q1=       q3*qze-q4*qet ! q33
      q4=       q4*qze+q3*qet ! q34
      qb(k3)=q2+q1            ! q3
      qb(k2)=q2-q1            ! q2
      qb(k1)=t1+q4            ! q1
      qb(k4)=t1-q4            ! q4
      enddo
      enddo
      ma=ma5
      enddo
      endif
      return
      end

c------------------------------------------------------------------------------c
c   r.j.purser, national meteorological center, washington d.c.  1994          c
c                   function get235                                            c
c                                                                              c
c   analyze the factorization of n in terms of 2's, 3's, and 5's and verify    c
c   that these are the only prime factors. set ln2, ln3, ln5 respectively to   c
c   the number of powers of 2, 3, 5.                                           c
c --> n       number of data along the line of fourier transformation          c
c------------------------------------------------------------------------------c
      function get235(n)
      logical get235
      common/fft235/ ln2,ln3,ln5
      get235=.true.
      nkt=n
      nm=n-1
      nh=n/2
      ln2=-1
400   nk=nkt
      nkt=nk/2
      ln2=ln2+1
      if(2*nkt.eq.nk)goto 400
      nkt=nk
      ln3=-1
401   nk=nkt
      nkt=nk/3
      ln3=ln3+1
      if(3*nkt.eq.nk)goto 401
      nkt=nk
      ln5=-1
402   nk=nkt
      nkt=nk/5
      ln5=ln5+1
      if(5*nkt.eq.nk)goto 402
      if(nk.ne.1)get235=.false.
      return
      end

c------------------------------------------------------------------------------c
c   r.j.purser, national meteorological center, washington d.c.  1994          c
c                   subroutine rumble                                          c
c                                                                              c
c   initialize coefficient arrays jumble and tumble for use with the           c
c   fast-fourier-transform routines when the number n of data has the prime-   c
c   factorization 2**ln2*3**ln3*5**ln5                                         c
c                                                                              c
c <-- jumble  permutation of n data encoded as a sequence of transpositions    c
c <-- tumble  trigonometric coefficients for use in fft. the first half are    c
c             the cosines, the second half the sines, of uniformly increasing  c
c             relevant angles                                                  c
c------------------------------------------------------------------------------c
      subroutine rumble(jumble,tumble)
      parameter(ml=20)
      common/fft235/ ln2,ln3,ln5
      dimension nd(ml),md(ml),jumble(0:*),tumble(0:*)
      ln=ln2+ln3+ln5
      n=2**ln2*3**ln3*5**ln5
      nm=n-1
      nh=n/2
      pi2on=8.*atan(1.)/n
      do i=0,nh-1
      ang=pi2on*i
      tumble(i)=cos(ang)
      tumble(i+nh)=sin(ang)
      enddo
      id=1
      is=0
      do i=1,ln5
       is=is+1
       md(is)=id
       id=id*5
      enddo
      do i=1,ln3
       is=is+1
       md(is)=id
       id=id*3
      enddo
      do i=1,ln2
       is=is+1
       md(is)=id
       id=id*2
      enddo
      id=1
      do i=1,ln2
       nd(is)=id
       id=id*2
       is=is-1
      enddo
      do i=1,ln3
       nd(is)=id
       id=id*3
       is=is-1
      enddo
      do i=1,ln5
       nd(is)=id
       id=id*5
       is=is-1
      enddo
      jumble(0)=n
      do i=1,nm
       ir=i
       j=0
       do l=1,ln
        kd=ir/nd(l)
        ir=ir-kd*nd(l)
        j=j+kd*md(l)
       enddo
       jumble(i)=j
      enddo
      do i=1,nm
       j=jumble(i)
       if(j.lt.i)then
400        j=jumble(j)
           if(j.lt.i)goto 400
        jumble(i)=j
       endif
      enddo
      return
      end

      subroutine xcfft(rb,qb,n,ny,ncb,period,w,jumble)
      logical get235
      common/fft235/ ln2,ln3,ln5
      dimension rb(0:ncb-1,*),qb(0:ncb-1,*),w(0:*),jumble(0:*)
      rfac=period/n
300   if(n.ne.2**ln2*3**ln3*5**ln5)then
      if(.not.get235(n))stop'prime factors of n are not only 2, 3, 5'
      endif
      if(n.ne.jumble(0))call rumble(jumble,w)
      nm=n-1
      nh=n/2

c  scale and permute the data:
      do iy=1,ny
       rb(0,iy)=rb(0,iy)*rfac
       qb(0,iy)=qb(0,iy)*rfac
      enddo
      do i=1,nm
       j=jumble(i)
       if(j.gt.i)then
        do iy=1,ny
         t=rb(i,iy)
         rb(i,iy)=rb(j,iy)*rfac
         rb(j,iy)=t
         t=qb(i,iy)
         qb(i,iy)=qb(j,iy)*rfac
         qb(j,iy)=t
        enddo
       else
        do iy=1,ny
         rb(i,iy)=rb(i,iy)*rfac
         qb(i,iy)=qb(i,iy)*rfac
        enddo
       endif
      enddo
c  transform the data:
      ma=1
      mb=n
c  radix 4
      ls=0
      do l=2,ln2,2
       mb=mb/4
       ls=l
       ma4=ma*4
       mb2=mb*2
       do j=0,ma-1
        jmb=j*mb
        jmb2=j*mb2
        rf1=w(jmb)
        qf1=w(nh+jmb)
        rf2=w(jmb2)
        qf2=w(nh+jmb2)
        rf3=rf1*rf2-qf1*qf2
        qf3=rf1*qf2+qf1*rf2
        do i=0,nm,ma4
         k0=i+j
         k1=k0+ma
         k2=k1+ma
         k3=k2+ma
         do iy=1,ny
          t        =rb(k3,iy)*rf3-qb(k3,iy)*qf3    ! q13
          qb(k3,iy)=rb(k3,iy)*qf3+qb(k3,iy)*rf3    ! r13
          rb(k3,iy)=rb(k2,iy)*rf1-qb(k2,iy)*qf1    ! r12
          rb(k2,iy)=rb(k2,iy)*qf1+qb(k2,iy)*rf1    ! q12
          qb(k2,iy)=qb(k3,iy)-rb(k2,iy)            ! r23
          rb(k2,iy)=qb(k3,iy)+rb(k2,iy)            ! q22
          qb(k3,iy)=rb(k3,iy)+t                    ! r22
          t        =rb(k3,iy)-t                    ! q23
          rb(k3,iy)=rb(k1,iy)*rf2-qb(k1,iy)*qf2    ! r11
          qb(k1,iy)=rb(k1,iy)*qf2+qb(k1,iy)*rf2    ! q11
          rb(k1,iy)=rb(k0,iy)-rb(k3,iy)            ! r21
          rb(k0,iy)=rb(k0,iy)+rb(k3,iy)            ! r20
          rb(k3,iy)=rb(k1,iy)-qb(k2,iy)            ! r3
          rb(k1,iy)=rb(k1,iy)+qb(k2,iy)            ! r1
          qb(k2,iy)=qb(k0,iy)+qb(k1,iy)            ! q20
          qb(k1,iy)=qb(k0,iy)-qb(k1,iy)            ! q21
          qb(k0,iy)=qb(k2,iy)+rb(k2,iy)            ! q0
          qb(k2,iy)=qb(k2,iy)-rb(k2,iy)            ! q2
          rb(k2,iy)=rb(k0,iy)-qb(k3,iy)            ! r2
          rb(k0,iy)=rb(k0,iy)+qb(k3,iy)            ! r0
          qb(k3,iy)=qb(k1,iy)-t                    ! q3
          qb(k1,iy)=qb(k1,iy)+t                    ! q1
         enddo
        enddo
       enddo
       ma=ma4
      enddo

      if(ls.ne.ln2)then
c  radix 2
      mb=mb/2
      ma2=ma*2
      do j=0,ma-1
       jmb=j*mb
       do i=0,nm,ma2
        k0=j+i
        k1=k0+ma
        rf1=w(jmb)
        qf1=w(nh+jmb)
        do iy=1,ny
         t        =rb(k1,iy)*qf1+qb(k1,iy)*rf1 ! q11
         qb(k1,iy)=rb(k1,iy)*rf1-qb(k1,iy)*qf1 ! r11
         rb(k1,iy)=rb(k0,iy)-qb(k1,iy)         ! r1
         rb(k0,iy)=rb(k0,iy)+qb(k1,iy)         ! r0
         qb(k1,iy)=qb(k0,iy)-t                 ! q1
         qb(k0,iy)=qb(k0,iy)+t                 ! q0
        enddo
       enddo
      enddo
      ma=ma2
      endif
c  radix 3
      rep=-.5
      rec=1.5
      qep=.5*sqrt(3.)
      do l=1,ln3
       mb=mb/3
       ma3=ma*3
       do j=0,ma-1
        jmb=j*mb
        rf1=w(jmb)
        qf1=w(nh+jmb)
        rf2=rf1*rf1-qf1*qf1
        qf2=2*rf1*qf1
        do i=0,nm,ma3
         k0=i+j
         k1=k0+ma
         k2=k1+ma
         do iy=1,ny
          t        =rb(k2,iy)*qf2+qb(k2,iy)*rf2     ! r12
          qb(k2,iy)=rb(k2,iy)*rf2-qb(k2,iy)*qf2     ! q12
          rb(k2,iy)=rb(k1,iy)*qf1+qb(k1,iy)*rf1     ! q11
          rb(k1,iy)=rb(k1,iy)*rf1-qb(k1,iy)*qf1     ! r11
          qb(k1,iy)=rb(k2,iy)+t                     ! q21
          rb(k2,iy)=(rb(k2,iy)-t)*qep               ! r22
          t        =rb(k1,iy)+qb(k2,iy)             ! r21
          rb(k1,iy)=(rb(k1,iy)-qb(k2,iy))*qep       ! q22
          rb(k0,iy)=rb(k0,iy)+t                     ! r0
          qb(k0,iy)=qb(k0,iy)+qb(k1,iy)             ! q0
          t        =rb(k0,iy)-t*rec                 ! r21
          qb(k1,iy)=qb(k0,iy)-qb(k1,iy)*rec         ! q21
          qb(k2,iy)=qb(k1,iy)-rb(k1,iy)             ! q2
          qb(k1,iy)=qb(k1,iy)+rb(k1,iy)             ! q1
          rb(k1,iy)=t        -rb(k2,iy)             ! r1
          rb(k2,iy)=t        +rb(k2,iy)             ! r2
         enddo
        enddo
       enddo
       ma=ma3
      enddo

      if(ln5.gt.0)then
c  radix 5
      nze=n/5
      rze=w(nze)
      qze=w(nh+nze)
      rzc=1.-rze
      ret=rze*rze-qze*qze
      qet=2*rze*qze
      rec=1.-ret
      do l=1,ln5
       mb=mb/5
       ma5=ma*5
       do j=0,ma-1
        jmb=j*mb
        jmb2=jmb*2
        rf1=w(jmb)
        qf1=w(nh+jmb)
        rf2=w(jmb2)
        qf2=w(nh+jmb2)
        rf3=rf1*rf2-qf1*qf2
        qf3=rf1*qf2+qf1*rf2
        rf4=rf2*rf2-qf2*qf2
        qf4=2*rf2*qf2
        do i=0,nm,ma5
         k0=i+j
         k1=k0+ma
         k2=k1+ma
         k3=k2+ma
         k4=k3+ma
         do iy=1,ny
          t        =rb(k1,iy)*qf1+qb(k1,iy)*rf1           ! q11
          rb(k1,iy)=rb(k1,iy)*rf1-qb(k1,iy)*qf1           ! r11
          qb(k1,iy)=rb(k4,iy)*rf4-qb(k4,iy)*qf4           ! q14
          rb(k4,iy)=rb(k4,iy)*qf4+qb(k4,iy)*rf4           ! r14
          qb(k4,iy)=rb(k1,iy)-qb(k1,iy)                   ! q24
          rb(k1,iy)=rb(k1,iy)+qb(k1,iy)                   ! r21
          qb(k1,iy)=t        +rb(k4,iy)                   ! q21
          rb(k4,iy)=t        -rb(k4,iy)                   ! r24
          t        =rb(k3,iy)*rf3-qb(k3,iy)*qf3           ! q13
          rb(k3,iy)=rb(k3,iy)*qf3+qb(k3,iy)*rf3           ! r13
          qb(k3,iy)=rb(k2,iy)*qf2+qb(k2,iy)*rf2           ! q12
          rb(k2,iy)=rb(k2,iy)*rf2-qb(k2,iy)*qf2           ! r12
          qb(k2,iy)=qb(k3,iy)+rb(k3,iy)                   ! q22
          rb(k3,iy)=qb(k3,iy)-rb(k3,iy)                   ! r23
          qb(k3,iy)=rb(k2,iy)-t                           ! q23
          rb(k2,iy)=rb(k2,iy)+t                           ! r22
          rb(k0,iy)=rb(k0,iy)+rb(k1,iy)+rb(k2,iy)         ! r0
          qb(k0,iy)=qb(k0,iy)+qb(k1,iy)+qb(k2,iy)         ! q0
          t        =rb(k4,iy)*qze+rb(k3,iy)*qet           ! r34
          rb(k3,iy)=rb(k3,iy)*qze-rb(k4,iy)*qet           ! r33
          rb(k4,iy)=rb(k0,iy)-rb(k2,iy)*rzc-rb(k1,iy)*rec ! r32
          rb(k1,iy)=rb(k0,iy)-rb(k1,iy)*rzc-rb(k2,iy)*rec ! r31
          rb(k2,iy)=rb(k4,iy)+rb(k3,iy)                   ! r2
          rb(k3,iy)=rb(k4,iy)-rb(k3,iy)                   ! r3
          rb(k4,iy)=rb(k1,iy)+t                           ! r4
          rb(k1,iy)=rb(k1,iy)-t                           ! r1
          t        =qb(k0,iy)-qb(k1,iy)*rzc-qb(k2,iy)*rec ! q31
          qb(k2,iy)=qb(k0,iy)-qb(k2,iy)*rzc-qb(k1,iy)*rec ! q32
          qb(k1,iy)=qb(k3,iy)*qze-qb(k4,iy)*qet           ! q33
          qb(k4,iy)=qb(k4,iy)*qze+qb(k3,iy)*qet           ! q34
          qb(k3,iy)=qb(k2,iy)+qb(k1,iy)                   ! q3
          qb(k2,iy)=qb(k2,iy)-qb(k1,iy)                   ! q2
          qb(k1,iy)=t        +qb(k4,iy)                   ! q1
          qb(k4,iy)=t        -qb(k4,iy)                   ! q4
         enddo
        enddo
       enddo
       ma=ma5
      enddo
      endif
      return
      end

      subroutine xdfft(rb,qb,n,ny,ncb,period,w,jumble)
      logical get235
      common/fft235/ ln2,ln3,ln5
      dimension rb(0:ncb-1,*),qb(0:ncb-1,*),w(0:*),jumble(0:*)
      rfac=1./period

c  for fourier synthesis, reverse the order of wavenumbers:
      do j=n/2+1,n-1
       i=n-j
       do iy=1,ny
        t       =rb(j,iy)
        rb(j,iy)=rb(i,iy)
        rb(i,iy)=t
        t       =qb(j,iy)
        qb(j,iy)=qb(i,iy)
        qb(i,iy)=t
       enddo
      enddo

300   if(n.ne.2**ln2*3**ln3*5**ln5)then
      if(.not.get235(n))stop'prime factors of n are not only 2, 3, 5'
      endif
      if(n.ne.jumble(0))call rumble(jumble,w)
      nm=n-1
      nh=n/2

c  scale and permute the data:
      do iy=1,ny
       rb(0,iy)=rb(0,iy)*rfac
       qb(0,iy)=qb(0,iy)*rfac
      enddo
      do i=1,nm
       j=jumble(i)
       if(j.gt.i)then
        do iy=1,ny
         t=rb(i,iy)
         rb(i,iy)=rb(j,iy)*rfac
         rb(j,iy)=t
         t=qb(i,iy)
         qb(i,iy)=qb(j,iy)*rfac
         qb(j,iy)=t
        enddo
       else
        do iy=1,ny
         rb(i,iy)=rb(i,iy)*rfac
         qb(i,iy)=qb(i,iy)*rfac
        enddo
       endif
      enddo
c  transform the data:
      ma=1
      mb=n
c  radix 4
      ls=0
      do l=2,ln2,2
       mb=mb/4
       ls=l
       ma4=ma*4
       mb2=mb*2
       do j=0,ma-1
        jmb=j*mb
        jmb2=j*mb2
        rf1=w(jmb)
        qf1=w(nh+jmb)
        rf2=w(jmb2)
        qf2=w(nh+jmb2)
        rf3=rf1*rf2-qf1*qf2
        qf3=rf1*qf2+qf1*rf2
        do i=0,nm,ma4
         k0=i+j
         k1=k0+ma
         k2=k1+ma
         k3=k2+ma
         do iy=1,ny
          t        =rb(k3,iy)*rf3-qb(k3,iy)*qf3    ! q13
          qb(k3,iy)=rb(k3,iy)*qf3+qb(k3,iy)*rf3    ! r13
          rb(k3,iy)=rb(k2,iy)*rf1-qb(k2,iy)*qf1    ! r12
          rb(k2,iy)=rb(k2,iy)*qf1+qb(k2,iy)*rf1    ! q12
          qb(k2,iy)=qb(k3,iy)-rb(k2,iy)            ! r23
          rb(k2,iy)=qb(k3,iy)+rb(k2,iy)            ! q22
          qb(k3,iy)=rb(k3,iy)+t                    ! r22
          t        =rb(k3,iy)-t                    ! q23
          rb(k3,iy)=rb(k1,iy)*rf2-qb(k1,iy)*qf2    ! r11
          qb(k1,iy)=rb(k1,iy)*qf2+qb(k1,iy)*rf2    ! q11
          rb(k1,iy)=rb(k0,iy)-rb(k3,iy)            ! r21
          rb(k0,iy)=rb(k0,iy)+rb(k3,iy)            ! r20
          rb(k3,iy)=rb(k1,iy)-qb(k2,iy)            ! r3
          rb(k1,iy)=rb(k1,iy)+qb(k2,iy)            ! r1
          qb(k2,iy)=qb(k0,iy)+qb(k1,iy)            ! q20
          qb(k1,iy)=qb(k0,iy)-qb(k1,iy)            ! q21
          qb(k0,iy)=qb(k2,iy)+rb(k2,iy)            ! q0
          qb(k2,iy)=qb(k2,iy)-rb(k2,iy)            ! q2
          rb(k2,iy)=rb(k0,iy)-qb(k3,iy)            ! r2
          rb(k0,iy)=rb(k0,iy)+qb(k3,iy)            ! r0
          qb(k3,iy)=qb(k1,iy)-t                    ! q3
          qb(k1,iy)=qb(k1,iy)+t                    ! q1
         enddo
        enddo
       enddo
       ma=ma4
      enddo

      if(ls.ne.ln2)then
c  radix 2
      mb=mb/2
      ma2=ma*2
      do j=0,ma-1
       jmb=j*mb
       do i=0,nm,ma2
        k0=j+i
        k1=k0+ma
        rf1=w(jmb)
        qf1=w(nh+jmb)
        do iy=1,ny
         t        =rb(k1,iy)*qf1+qb(k1,iy)*rf1 ! q11
         qb(k1,iy)=rb(k1,iy)*rf1-qb(k1,iy)*qf1 ! r11
         rb(k1,iy)=rb(k0,iy)-qb(k1,iy)         ! r1
         rb(k0,iy)=rb(k0,iy)+qb(k1,iy)         ! r0
         qb(k1,iy)=qb(k0,iy)-t                 ! q1
         qb(k0,iy)=qb(k0,iy)+t                 ! q0
        enddo
       enddo
      enddo
      ma=ma2
      endif
c  radix 3
      rep=-.5
      rec=1.5
      qep=.5*sqrt(3.)
      do l=1,ln3
       mb=mb/3
       ma3=ma*3
       do j=0,ma-1
        jmb=j*mb
        rf1=w(jmb)
        qf1=w(nh+jmb)
        rf2=rf1*rf1-qf1*qf1
        qf2=2*rf1*qf1
        do i=0,nm,ma3
         k0=i+j
         k1=k0+ma
         k2=k1+ma
         do iy=1,ny
          t        =rb(k2,iy)*qf2+qb(k2,iy)*rf2     ! r12
          qb(k2,iy)=rb(k2,iy)*rf2-qb(k2,iy)*qf2     ! q12
          rb(k2,iy)=rb(k1,iy)*qf1+qb(k1,iy)*rf1     ! q11
          rb(k1,iy)=rb(k1,iy)*rf1-qb(k1,iy)*qf1     ! r11
          qb(k1,iy)=rb(k2,iy)+t                     ! q21
          rb(k2,iy)=(rb(k2,iy)-t)*qep               ! r22
          t        =rb(k1,iy)+qb(k2,iy)             ! r21
          rb(k1,iy)=(rb(k1,iy)-qb(k2,iy))*qep       ! q22
          rb(k0,iy)=rb(k0,iy)+t                     ! r0
          qb(k0,iy)=qb(k0,iy)+qb(k1,iy)             ! q0
          t        =rb(k0,iy)-t*rec                 ! r21
          qb(k1,iy)=qb(k0,iy)-qb(k1,iy)*rec         ! q21
          qb(k2,iy)=qb(k1,iy)-rb(k1,iy)             ! q2
          qb(k1,iy)=qb(k1,iy)+rb(k1,iy)             ! q1
          rb(k1,iy)=t        -rb(k2,iy)             ! r1
          rb(k2,iy)=t        +rb(k2,iy)             ! r2
         enddo
        enddo
       enddo
       ma=ma3
      enddo

      if(ln5.gt.0)then
c  radix 5
      nze=n/5
      rze=w(nze)
      qze=w(nh+nze)
      rzc=1.-rze
      ret=rze*rze-qze*qze
      qet=2*rze*qze
      rec=1.-ret
      do l=1,ln5
       mb=mb/5
       ma5=ma*5
       do j=0,ma-1
        jmb=j*mb
        jmb2=jmb*2
        rf1=w(jmb)
        qf1=w(nh+jmb)
        rf2=w(jmb2)
        qf2=w(nh+jmb2)
        rf3=rf1*rf2-qf1*qf2
        qf3=rf1*qf2+qf1*rf2
        rf4=rf2*rf2-qf2*qf2
        qf4=2*rf2*qf2
        do i=0,nm,ma5
         k0=i+j
         k1=k0+ma
         k2=k1+ma
         k3=k2+ma
         k4=k3+ma
         do iy=1,ny
          t        =rb(k1,iy)*qf1+qb(k1,iy)*rf1           ! q11
          rb(k1,iy)=rb(k1,iy)*rf1-qb(k1,iy)*qf1           ! r11
          qb(k1,iy)=rb(k4,iy)*rf4-qb(k4,iy)*qf4           ! q14
          rb(k4,iy)=rb(k4,iy)*qf4+qb(k4,iy)*rf4           ! r14
          qb(k4,iy)=rb(k1,iy)-qb(k1,iy)                   ! q24
          rb(k1,iy)=rb(k1,iy)+qb(k1,iy)                   ! r21
          qb(k1,iy)=t        +rb(k4,iy)                   ! q21
          rb(k4,iy)=t        -rb(k4,iy)                   ! r24
          t        =rb(k3,iy)*rf3-qb(k3,iy)*qf3           ! q13
          rb(k3,iy)=rb(k3,iy)*qf3+qb(k3,iy)*rf3           ! r13
          qb(k3,iy)=rb(k2,iy)*qf2+qb(k2,iy)*rf2           ! q12
          rb(k2,iy)=rb(k2,iy)*rf2-qb(k2,iy)*qf2           ! r12
          qb(k2,iy)=qb(k3,iy)+rb(k3,iy)                   ! q22
          rb(k3,iy)=qb(k3,iy)-rb(k3,iy)                   ! r23
          qb(k3,iy)=rb(k2,iy)-t                           ! q23
          rb(k2,iy)=rb(k2,iy)+t                           ! r22
          rb(k0,iy)=rb(k0,iy)+rb(k1,iy)+rb(k2,iy)         ! r0
          qb(k0,iy)=qb(k0,iy)+qb(k1,iy)+qb(k2,iy)         ! q0
          t        =rb(k4,iy)*qze+rb(k3,iy)*qet           ! r34
          rb(k3,iy)=rb(k3,iy)*qze-rb(k4,iy)*qet           ! r33
          rb(k4,iy)=rb(k0,iy)-rb(k2,iy)*rzc-rb(k1,iy)*rec ! r32
          rb(k1,iy)=rb(k0,iy)-rb(k1,iy)*rzc-rb(k2,iy)*rec ! r31
          rb(k2,iy)=rb(k4,iy)+rb(k3,iy)                   ! r2
          rb(k3,iy)=rb(k4,iy)-rb(k3,iy)                   ! r3
          rb(k4,iy)=rb(k1,iy)+t                           ! r4
          rb(k1,iy)=rb(k1,iy)-t                           ! r1
          t        =qb(k0,iy)-qb(k1,iy)*rzc-qb(k2,iy)*rec ! q31
          qb(k2,iy)=qb(k0,iy)-qb(k2,iy)*rzc-qb(k1,iy)*rec ! q32
          qb(k1,iy)=qb(k3,iy)*qze-qb(k4,iy)*qet           ! q33
          qb(k4,iy)=qb(k4,iy)*qze+qb(k3,iy)*qet           ! q34
          qb(k3,iy)=qb(k2,iy)+qb(k1,iy)                   ! q3
          qb(k2,iy)=qb(k2,iy)-qb(k1,iy)                   ! q2
          qb(k1,iy)=t        +qb(k4,iy)                   ! q1
          qb(k4,iy)=t        -qb(k4,iy)                   ! q4
         enddo
        enddo
       enddo
       ma=ma5
      enddo
      endif
      return
      end

      subroutine ycfft(rb,qb,n,nx,ndx,ncb,period,w,jumble)
      logical get235
      common/fft235/ ln2,ln3,ln5
      dimension rb(ncb,0:*),qb(ncb,0:*),w(0:*),jumble(0:*)
      rfac=period/n
      goto 300

      entry      ydfft(rb,qb,n,nx,ndx,ncb,period,w,jumble)
      rfac=1./period

c  for fourier synthesis, reverse the order of wavenumbers:
      do j=n/2+1,n-1
       i=n-j
       do ix=1,nx,ndx
        t       =rb(ix,j)
        rb(ix,j)=rb(ix,i)
        rb(ix,i)=t
        t       =qb(ix,j)
        qb(ix,j)=qb(ix,i)
        qb(ix,i)=t
       enddo
      enddo

300   if(n.ne.2**ln2*3**ln3*5**ln5)then
      if(.not.get235(n))stop'prime factors of n are not only 2, 3, 5'
      endif
      if(n.ne.jumble(0))call rumble(jumble,w)
      nm=n-1
      nh=n/2

c  scale and permute the data:
      do ix=1,nx,ndx
       rb(ix,0)=rb(ix,0)*rfac
       qb(ix,0)=qb(ix,0)*rfac
      enddo
      do i=1,nm
       j=jumble(i)
       if(j.gt.i)then
        do ix=1,nx,ndx
         t=rb(ix,i)
         rb(ix,i)=rb(ix,j)*rfac
         rb(ix,j)=t
         t=qb(ix,i)
         qb(ix,i)=qb(ix,j)*rfac
         qb(ix,j)=t
        enddo
       else
        do ix=1,nx,ndx
         rb(ix,i)=rb(ix,i)*rfac
         qb(ix,i)=qb(ix,i)*rfac
        enddo
       endif
      enddo

c  transform the data:
      ma=1
      mb=n
c  radix 4
      ls=0
      do l=2,ln2,2
       mb=mb/4
       ls=l
       ma4=ma*4
       mb2=mb*2
       do j=0,ma-1
        jmb=j*mb
        jmb2=j*mb2
        rf1=w(jmb)
        qf1=w(nh+jmb)
        rf2=w(jmb2)
        qf2=w(nh+jmb2)
        rf3=rf1*rf2-qf1*qf2
        qf3=rf1*qf2+qf1*rf2
        do i=0,nm,ma4
         k0=i+j
         k1=k0+ma
         k2=k1+ma
         k3=k2+ma
         do ix=1,nx,ndx
          t        =rb(ix,k3)*rf3-qb(ix,k3)*qf3     ! q13
          qb(ix,k3)=rb(ix,k3)*qf3+qb(ix,k3)*rf3     ! r13
          rb(ix,k3)=rb(ix,k2)*rf1-qb(ix,k2)*qf1     ! r12
          rb(ix,k2)=rb(ix,k2)*qf1+qb(ix,k2)*rf1     ! q12
          qb(ix,k2)=qb(ix,k3)-rb(ix,k2)             ! r23
          rb(ix,k2)=qb(ix,k3)+rb(ix,k2)             ! q22
          qb(ix,k3)=rb(ix,k3)+t                     ! r22
          t        =rb(ix,k3)-t                     ! q23
          rb(ix,k3)=rb(ix,k1)*rf2-qb(ix,k1)*qf2     ! r11
          qb(ix,k1)=rb(ix,k1)*qf2+qb(ix,k1)*rf2     ! q11
          rb(ix,k1)=rb(ix,k0)-rb(ix,k3)             ! r21
          rb(ix,k0)=rb(ix,k0)+rb(ix,k3)             ! r20
          rb(ix,k3)=rb(ix,k1)-qb(ix,k2)             ! r3
          rb(ix,k1)=rb(ix,k1)+qb(ix,k2)             ! r1
          qb(ix,k2)=qb(ix,k0)+qb(ix,k1)             ! q20
          qb(ix,k1)=qb(ix,k0)-qb(ix,k1)             ! q21
          qb(ix,k0)=qb(ix,k2)+rb(ix,k2)             ! q0
          qb(ix,k2)=qb(ix,k2)-rb(ix,k2)             ! q2
          rb(ix,k2)=rb(ix,k0)-qb(ix,k3)             ! r2
          rb(ix,k0)=rb(ix,k0)+qb(ix,k3)             ! r0
          qb(ix,k3)=qb(ix,k1)-t                     ! q3
          qb(ix,k1)=qb(ix,k1)+t                     ! q1
         enddo
        enddo
       enddo
       ma=ma4
      enddo
      if(ls.ne.ln2)then

c  radix 2
      mb=mb/2
      ma2=ma*2
      do j=0,ma-1
       jmb=j*mb
       do i=0,nm,ma2
        k0=j+i
        k1=k0+ma
        rf1=w(jmb)
        qf1=w(nh+jmb)
        do ix=1,nx,ndx
         t        =rb(ix,k1)*qf1+qb(ix,k1)*rf1 ! q11
         qb(ix,k1)=rb(ix,k1)*rf1-qb(ix,k1)*qf1 ! r11
         rb(ix,k1)=rb(ix,k0)-qb(ix,k1)         ! r1
         rb(ix,k0)=rb(ix,k0)+qb(ix,k1)         ! r0
         qb(ix,k1)=qb(ix,k0)-t                 ! q1
         qb(ix,k0)=qb(ix,k0)+t                 ! q0
        enddo
       enddo
      enddo
      ma=ma2
      endif

c  radix 3
      rep=-.5
      rec=1.5
      qep=.5*sqrt(3.)
      do l=1,ln3
       mb=mb/3
       ma3=ma*3
       do j=0,ma-1
        jmb=j*mb
        rf1=w(jmb)
        qf1=w(nh+jmb)
        rf2=rf1*rf1-qf1*qf1
        qf2=2*rf1*qf1
        do i=0,nm,ma3
         k0=i+j
         k1=k0+ma
         k2=k1+ma
         t        =rb(ix,k2)*qf2+qb(ix,k2)*rf2     ! r12
         qb(ix,k2)=rb(ix,k2)*rf2-qb(ix,k2)*qf2     ! q12
         rb(ix,k2)=rb(ix,k1)*qf1+qb(ix,k1)*rf1     ! q11
         rb(ix,k1)=rb(ix,k1)*rf1-qb(ix,k1)*qf1     ! r11
         qb(ix,k1)=rb(ix,k2)+t                     ! q21
         rb(ix,k2)=(rb(ix,k2)-t        )*qep       ! r22
         t        =rb(ix,k1)+qb(ix,k2)             ! r21
         rb(ix,k1)=(rb(ix,k1)-qb(ix,k2))*qep       ! q22
         rb(ix,k0)=rb(ix,k0)+t                     ! r0
         qb(ix,k0)=qb(ix,k0)+qb(ix,k1)             ! q0
         t        =rb(ix,k0)-t        *rec         ! r21
         qb(ix,k1)=qb(ix,k0)-qb(ix,k1)*rec         ! q21
         qb(ix,k2)=qb(ix,k1)-rb(ix,k1)             ! q2
         qb(ix,k1)=qb(ix,k1)+rb(ix,k1)             ! q1
         rb(ix,k1)=t        -rb(ix,k2)             ! r1
         rb(ix,k2)=t        +rb(ix,k2)             ! r2
        enddo
       enddo
       ma=ma3
      enddo
      if(ln5.gt.0)then

c  radix 5
      nze=n/5
      rze=w(nze)
      qze=w(nh+nze)
      rzc=1.-rze
      ret=rze*rze-qze*qze
      qet=2*rze*qze
      rec=1.-ret
      do l=1,ln5
       mb=mb/5
       ma5=ma*5
       do j=0,ma-1
        jmb=j*mb
        jmb2=jmb*2
        rf1=w(jmb)
        qf1=w(nh+jmb)
        rf2=w(jmb2)
        qf2=w(nh+jmb2)
        rf3=rf1*rf2-qf1*qf2
        qf3=rf1*qf2+qf1*rf2
        rf4=rf2*rf2-qf2*qf2
        qf4=2*rf2*qf2
        do i=0,nm,ma5
         k0=i+j
         k1=k0+ma
         k2=k1+ma
         k3=k2+ma
         k4=k3+ma
         t        =rb(ix,k1)*qf1+qb(ix,k1)*rf1           ! q11
         rb(ix,k1)=rb(ix,k1)*rf1-qb(ix,k1)*qf1           ! r11
         qb(ix,k1)=rb(ix,k4)*rf4-qb(ix,k4)*qf4           ! q14
         rb(ix,k4)=rb(ix,k4)*qf4+qb(ix,k4)*rf4           ! r14
         qb(ix,k4)=rb(ix,k1)-qb(ix,k1)                   ! q24
         rb(ix,k1)=rb(ix,k1)+qb(ix,k1)                   ! r21
         qb(ix,k1)=t        +rb(ix,k4)                   ! q21
         rb(ix,k4)=t        -rb(ix,k4)                   ! r24
         t        =rb(ix,k3)*rf3-qb(ix,k3)*qf3           ! q13
         rb(ix,k3)=rb(ix,k3)*qf3+qb(ix,k3)*rf3           ! r13
         qb(ix,k3)=rb(ix,k2)*qf2+qb(ix,k2)*rf2           ! q12
         rb(ix,k2)=rb(ix,k2)*rf2-qb(ix,k2)*qf2           ! r12
         qb(ix,k2)=qb(ix,k3)+rb(ix,k3)                   ! q22
         rb(ix,k3)=qb(ix,k3)-rb(ix,k3)                   ! r23
         qb(ix,k3)=rb(ix,k2)-t                           ! q23
         rb(ix,k2)=rb(ix,k2)+t                           ! r22
         rb(ix,k0)=rb(ix,k0)+rb(ix,k1)+rb(ix,k2)         ! r0
         qb(ix,k0)=qb(ix,k0)+qb(ix,k1)+qb(ix,k2)         ! q0
         t        =rb(ix,k4)*qze+rb(ix,k3)*qet           ! r34
         rb(ix,k3)=rb(ix,k3)*qze-rb(ix,k4)*qet           ! r33
         rb(ix,k4)=rb(ix,k0)-rb(ix,k2)*rzc-rb(ix,k1)*rec ! r32
         rb(ix,k1)=rb(ix,k0)-rb(ix,k1)*rzc-rb(ix,k2)*rec ! r31
         rb(ix,k2)=rb(ix,k4)+rb(ix,k3)                   ! r2
         rb(ix,k3)=rb(ix,k4)-rb(ix,k3)                   ! r3
         rb(ix,k4)=rb(ix,k1)+t                           ! r4
         rb(ix,k1)=rb(ix,k1)-t                           ! r1
         t        =qb(ix,k0)-qb(ix,k1)*rzc-qb(ix,k2)*rec ! q31
         qb(ix,k2)=qb(ix,k0)-qb(ix,k2)*rzc-qb(ix,k1)*rec ! q32
         qb(ix,k1)=qb(ix,k3)*qze-qb(ix,k4)*qet           ! q33
         qb(ix,k4)=qb(ix,k4)*qze+qb(ix,k3)*qet           ! q34
         qb(ix,k3)=qb(ix,k2)+qb(ix,k1)                   ! q3
         qb(ix,k2)=qb(ix,k2)-qb(ix,k1)                   ! q2
         qb(ix,k1)=t        +qb(ix,k4)                   ! q1
         qb(ix,k4)=t        -qb(ix,k4)                   ! q4
        enddo
       enddo
       ma=ma5
      enddo
      endif
      return
      end
