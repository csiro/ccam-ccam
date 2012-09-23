!  this is eig derived from eignew, but used in-line in C-CAM
      subroutine eig(sigin,sigmhin,tbar,lapsbot,isoth,dtin,eps,nsig,
     &               bet,betm,nh)
      use vecs_m
      implicit none
      include 'newmpar.h'
      integer nh,nsig,lapsbot,isoth
      integer nchng,k
      integer l
      integer, parameter :: neig = 1
      real tem
      real eps,dtin,dt
      real sigin(kl),sigmhin(kl)
      real sig(kl),sigmh(kl+1),tbar(kl)
      real bet(kl),betm(kl)

!     lapsbot=1 gives zero lowest t lapse for phi calc
      write(6,*) 'this run configured with kl = ',kl
      write(6,*) 'entering eig tbar,lapsbot,isoth,dtin,eps,nh: ',
     &                      tbar(1),lapsbot,isoth,dtin,eps,nh
      dt=dtin
      sig(:)=sigin(:)
      sigmh(1:kl)=sigmhin(:)
      sigmh(kl+1)=0.
!     expect data from bottom up
      if(sig(1)<sig(kl))then
        call flip3(sig,1,1,kl,1)
        bam(1:kl)=sigmhin(:)
        call flip3(bam,1,1,kl,1)  
        sigmh(1)=0.
        sigmh(2:kl+1)=bam(1:kl)
      endif  ! (sig(1)<sig(kl))
      write(6,*) 'final sig values: ',sig
      write(6,*) 'final sigmh values: ',sigmh
      open(28,file='eigenv.out')
      call eigs(lapsbot,isoth,tbar,dt,eps,nh,sig,sigmh,
     &          bet,betm)   !------------------------
      write(6,*) 'about to write to 28 '
      write(28,*)kl,lapsbot,isoth,nsig,
     &       '   kl,lapsbot,isoth,nsig'
!     re-order the eigenvectors if necessary
      nchng=1
      do while (nchng/=0)
       nchng=0
       do k=2,kl
        if(bam(k)<bam(k-1)) cycle
        nchng=1
        tem=bam(k)
        bam(k)=bam(k-1)
        bam(k-1)=tem
        do l=1,kl
         tem=emat(l,k)
         emat(l,k)=emat(l,k-1)
         emat(l,k-1)=tem
         tem=einv(k,l)
         einv(k,l)=einv(k-1,l)
         einv(k-1,l)=tem
        end do
       end do
      end do
      write(6,*)'eigenvectors re-ordered'
      write(6,*)'bam',(bam(k),k=1,kl)
      if(neig==1)then
!      write data from bottom up
       if(sig(1)<sig(kl))then
        call flip3( sig,1, 1,kl, 1)
        call flip3( sigmh,1, 1,kl+1, 1)
        call flip3(emat,1, 1,kl,kl)
        call flip3(einv,1,kl,kl, 1)
       endif
      endif
      write(6,*) 'eig '
      do k=1,kl
       write(6,'(i3,15f8.4/3x,15f8.4/3x,15f8.4)')
     &   k,(emat(k,l),l=1,kl)
      enddo
      write(6,*) 'einv'
      do k=1,kl
       write(6,'(i3,15f8.4/3x,15f8.4/3x,15f8.4)')
     &   k,(einv(k,l),l=1,kl)
      enddo
      write(28,'(9e14.6)')
     &  (sig(k),k=1,kl),(tbar(k),k=1,kl),(bam(k),k=1,kl)
     &  ,((emat(k,l),k=1,kl),l=1,kl),((einv(k,l),k=1,kl),l=1,kl),
     &  (bam(k),k=1,kl),((emat(k,l),k=1,kl),l=1,kl) ! just to fill space
      write(28,'(9e14.6)') (sigmh(k),k=1,kl+1)
      end

      subroutine eigs(lapsbot,isoth,tbar,dt,eps,nh,sig,sigmh,
     &                bet,betm)
      use vecs_m
      implicit none
      include 'newmpar.h'
      integer lapsbot,isoth,nh
      integer k,l,irror
      integer indic(kl)
      real dt,eps
      real factg,factr,dp
c     units here are SI, but final output is dimensionless
      real, parameter :: g=9.806
      real, parameter :: cp=1004.64
      real, parameter :: r=287.
c     sets up eigenvectors
      real sig(kl),sigmh(kl+1)
      real dsig(kl)
      real bet(kl),betm(kl),get(kl),getm(kl),gmat(kl,kl)
      real bmat(kl,kl),evimag(kl),veci(kl,kl),sum1(kl)
      real tbar(kl)
      real aa(kl,kl),ab(kl,kl),ac(kl,kl)
      real aaa(kl,kl),cc(kl,kl)
      
      aa=0.
      bmat=0.

      do k=1,kl
       dsig(k)=sigmh(k+1)-sigmh(k)
      enddo
      write(6,*) 'sigmh ',(sigmh(k),k=1,kl)
      write(6,*) 'dsig ',(dsig(k),k=1,kl)

      get(1)=bet(1)/(r*sig(1))
      do k=2,kl
        get(k)=bet(k)/(r*sig(k))
        getm(k)=betm(k)/(r*sig(k-1))
      enddo      
      factg=2./(dt*(1.+eps))      
      factr=factg*r*r*tbar(1)*tbar(1)/(g*g)
      
      bmat(:,:)=0.  ! N.B. bmat includes effect of r/sig weighting
      gmat(:,:)=0.  ! N.B. gmat may include effect of 1/sig**2 weighting
      do k=2,kl
       do l=1,k-1
        bmat(k,l)=bet(l)+betm(l+1)
        gmat(k,l)=factr*(get(l)+getm(l+1))
       enddo ! l loop
      enddo  ! k loop
      do k=1,kl
       bmat(k,k)=bet(k)
       gmat(k,k)=factr*get(k)
      enddo
   
      write(6,*)'dt,eps,factg,tbar,factr ',dt,eps,factg,tbar(1),factr
      write(6,*)'bet ',bet
      write(6,*)'bmat'
      do k=1,kl
       write(6,'(i3,15f8.3/3x,15f8.3/3x,15f8.3)') k,(bmat(k,l),l=1,kl)
      enddo

      write(6,*) 'get ',get
      write(6,*) 'getm ',getm
      write(6,*) 'gmat'
      do k=1,kl
       write(6,'(i3,15f8.0/3x,15f8.0/3x,15f8.0)') k,(gmat(k,l),l=1,kl)
      enddo

!     even newer derivation section
      write(6,*) 'even newer derivation section'

      do k=1,kl
       do l=1,kl
        ab(k,l)=dsig(l)
        ac(k,l)=dsig(l)
       enddo
      enddo
      do k=1,kl
       do l=k,kl
        aa(k,l)=-r*tbar(1)*dsig(l)/(cp*sig(k))
       enddo
      enddo
      do k=1,kl
       aa(k,k)=-r*tbar(1)*(sigmh(k+1)-sig(k))/(cp*sig(k))
       ac(k,k)=ac(k,k)+1.
      enddo
      write(6,*) 'aa'
      do k=1,kl
       write(6,'(i3,15f8.4/3x,15f8.4/3x,15f8.4)') 
     &   k,(aa(k,l),l=1,kl)
      enddo
      write(6,*) 'ac'
      do k=1,kl
       write(6,'(i3,15f8.3/3x,15f8.3/3x,15f8.3)')
     &   k,(ac(k,l),l=1,kl)
      enddo
      if(isoth.eq.1)then  !  extra vadv terms added
        aa(:,:)=aa(:,:)+tbar(1)*ac(:,:)
      endif
      call matm(aaa,bmat,aa)
      cc(:,:)=aaa(:,:)-r*tbar(1)*ab(:,:)
      write(6,*) 'cc'
      do k=1,kl
       write(6,'(i3,15f8.1/3x,15f8.1/3x,15f8.1)')
     &   k,(cc(k,l),l=1,kl)
      enddo

      if(nh>0)then  ! use gmat instead of bmat to derive aaa ! MJT suggestion
        gmat(:,:)=bmat(:,:)*(1.+4.*cp*tbar(1)/
     &           ((g*dt*(1.+eps))**2))
        call matm(aaa,gmat,aa)
        cc(:,:)=aaa(:,:)-r*tbar(1)*ab(:,:)
        write(6,*) 'cc with gmat'
        do k=1,kl
         write(6,'(i3,15f8.1/3x,15f8.1/3x,15f8.1)')
     &     k,(cc(k,l),l=1,kl)
        enddo
      endif  ! (nh==2)
      
      aaa(:,:)=cc(:,:)
      call eigenp(aaa,bam,evimag,emat,veci,indic)
      write(6,*) 'bam',(bam(k),k=1,kl)
      write(6,*) 'eig '
      do k=1,kl
       write(6,'(i3,15f8.4/3x,15f8.4/3x,15f8.4)')
     &   k,(emat(k,l),l=1,kl)
      enddo
      call matinv(cc,sum1,0,dp,irror)
      einv(:,:)=emat(:,:)
      call matinv(einv,sum1,0,dp,irror)
      write(6,*) 'einv'
      do k=1,kl
       write(6,'(i3,15f8.4/3x,15f8.4/3x,15f8.4)')
     &   k,(einv(k,l),l=1,kl)
      enddo

      return
      end

      subroutine flip3(a,il,jl,kl,ll)
      implicit none
      integer, intent(in) :: il,jl,kl,ll
      integer l,j,i,k
      real, dimension(il,jl,kl,ll), intent(inout) :: a
      real tem
      do l=1,ll
       do j=1,jl
        do i=1,il
         do k=1,kl/2
          tem=a(i,j,k,l)
          a(i,j,k,l)=a(i,j,kl+1-k,l)
          a(i,j,kl+1-k,l)=tem
         end do
        end do
       end do
      end do
      return
      end

      subroutine matm(a,b,c)
      implicit none
      include 'newmpar.h'
c     matrix multiplication      a = b * c
      integer k,l,ll,ivec
      real a(kl,kl),b(kl,kl),c(kl,kl)
      do k=1,kl
       do l=1,kl
        a(k,l)=b(k,1)*c(1,l)
        do ll=2,kl
         a(k,l)=a(k,l)+b(k,ll)*c(ll,l)
        end do
       end do
      end do
      return
      end

      subroutine eigenp(a,evr,evi,vecr,veci,indic)
      implicit none
      include 'newmpar.h'
      integer, dimension(kl*kl) :: iwork,local
      integer, dimension(kl) :: indic
      integer i,j,k,l,m
      integer k1,l1,kon,ivec
      real, dimension(kl*kl) :: prfact,subdia,work
      real, dimension(kl,kl) :: a,vecr,veci
      real, dimension(kl) :: evr,evi
      real d1,d2,d3,enorm
      real r,r1,ex,eps
      
c     currently set up for a general number of 10 levels
c a.c.m. algorithm number 343
c revised july, 1970 by n.r.pummeroy, dcr, csiro, canberra
c the following variables were changed
c from single to real in eigenp:
c r,r1,enorm,eps,ex,work,work1,work2,subdia
c see also comments for routines scaler, hesqr, realve, compve
c
c
c this sub.routine finds all the eigenvalues and the
c eigenvectors of a real general matrix of order n.
c
c first in the sub.routine scaler the matrix is scaled so that
c the corresponding rows and columns are approximately
c balanced and then the matrix is normalised so that the
c value of the euclidian norm of the matrix is equal to one.
c
c the eigenvalues are computed by the qr double-step method
c in the sub.routine hesqr.
c the eigenvectors are computed by inverse iteration in
c the sub.routine realve,for the real eigenvalues,or in the
c subroutine compve,for the complex eigenvalues.
c
c the elements of the matrix are to be stored in the first n
c rows and columns of the two dimensional array a. the
c original matrix is destroyed by the sub.routine.
c n is the order of the matrix.
c nm defines the first dimension of the two dimensional
c arrays a,vecr,veci and the dimension of the one
c the real parts of the n computed eigenvalues will be found
c in the first n places of the array evr and the imaginary
c parts in the first n places of the array evi.
c the real components of the normalised eigenvector i
c (i=1,2,...,n) corresponding to the eigenvalue stored in
c evr(i) and evi(i) will be found in the first n places of
c the column i of the two dimensional array vecr and the
c imaginary components in the first n places of the column i
c of the two dimensional array veci.
c
c the real eigenvector is normalised so that the sum of the
c squares of the components is equal to one.
c the complex eigenvector is normalised so that the
c component with the largest value in modulus has its real
c part equal to one and the imaginary part equal to zero.
c
c the array indic indicates the success of the sub.routine
c eigenp as follows
c     value of indic(i)   eigenvalue i   eigenvector i
c            0              not found      not found
c            1              found          not found
c            2              found          found
c
c
      call scaler(a,veci,prfact,enorm)
c the computation of the eigenvalues of the normalised
c matrix.
!  take t=50 significant binary figures.  ex=2**(-t)
!     ex=8.88178418e-16
!  following for 60 binary figures:
      ex=8.674e-19
      call hesqr(a,veci,evr,evi,subdia,indic,eps,ex)
c
c the possible decomposition of the upper-hessenberg matrix
c into the submatrices of lower order is indicated in the
c array local. the decomposition occurs when some
c subdiagonal elements are in modulus less than a small
c positive number eps defined in the sub.routine hesqr . the
c amount of work in the eigenvector problem may be
c minimised in this way.
      j = kl
      i = 1
      local(1) = 1
    2 if(abs(subdia(j-1)).gt.eps)go to 3
      i = i+1
      local(i)=0
    3 j = j-1
      local(i)=local(i)+1
      if(j.ne.1)go to 2
c
c the eigenvector problem.
      k = 1
      kon = 0
      l = local(1)
      m = kl
      do 10 i=1,kl
        ivec = kl-i+1
        if(i.le.l)go to 5
        k = k+1
        m = kl-l
        l = l+local(k)
    5   if(indic(ivec).eq.0)go to 10
        if(evi(ivec).ne.0.0)go to 8
c
c transfer of an upper-hessenberg matrix of the order m from
c the arrays veci and subdia into the array a.
        do 7 k1=1,m
          do 6 l1=k1,m
    6       a(k1,l1) = veci(k1,l1)
          if(k1.eq.1)go to 7
          a(k1,k1-1) = subdia(k1-1)
    7     continue
c
c the computation of the real engenvector ivec of the upper-
c hessenberg matrix corresponding to the real eigenvalue
c evr(ivec).
        call realve(m,ivec,a,vecr,evr,evi,iwork,
     1  work,indic,eps,ex)
        go to 10
c
c the computation of the complex eigenvector ivec of the
c upper-hessenberg matrix corresponding to the complex
c eigenvalue evr(ivec) + i*evi(ivec). if the value of kon is
c not equal to zero then this complex eigenvector has
c already been found from its conjugate.
    8   if(kon.ne.0)go to 9
        kon = 1
      write(6,*) 'attempted call to comove'
      stop
    9   kon = 0
   10   continue
c
c the reconstruction of the matrix used in the reduction of
c matrix a to an upper-hessenberg form by householder method
      do 12 i=1,kl
        do 11 j=i,kl
          a(i,j) = 0.0
   11     a(j,i) = 0.0
   12   a(i,i) = 1.0
      if(kl.le.2)go to 15
      m = kl-2
      do 14 k=1,m
        l = k+1
        do 14 j=2,kl
          d1 = 0.0
          do 13 i=l,kl
            d2 = veci(i,k)
   13       d1 = d1+ d2*a(j,i)
          do 14 i=l,kl
   14       a(j,i) = a(j,i)-veci(i,k)*d1
c
c the computation of the eigenvectors of the original non-
c scaled matrix.
   15 kon = 1
      do 24 i=1,kl
        l = 0
        if(evi(i).eq.0.0)go to 16
        l = 1
        if(kon.eq.0)go to 16
        kon = 0
        go to 24
   16   do 18 j=1,kl
      d1 = 0.0
      d2 = 0.0
          do 17 k=1,kl
            d3 = a(j,k)
            d1 = d1+d3*vecr(k,i)
            if(l.eq.0)go to 17
            d2 = d2+d3*vecr(k,i-1)
   17       continue
          work(j) = d1/prfact(j)
          if(l.eq.0)go to 18
          subdia(j)=d2/prfact(j)
   18     continue
c
c the normalisation of the eigenvectors and the computation
c of the eigenvalues of the original non-normalised matrix.
        if(l.eq.1)go to 21
        d1 = 0.0
        do 19 m=1,kl
   19     d1 = d1+work(m)**2
        d1 = sqrt(d1)
        do 20 m=1,kl
          veci(m,i) = 0.0
   20     vecr(m,i) = work(m)/d1
        evr(i) = evr(i)*enorm
        go to 24
c
   21   kon = 1
        evr(i) = evr(i)*enorm
        evr(i-1) = evr(i)
        evi(i) = evi(i)*enorm
        evi(i-1) =-evi(i)
        r = 0.0
        do 22 j=1,kl
          r1 = work(j)**2 + subdia(j)**2
          if(r.ge.r1)go to 22
          r = r1
          l = j
   22     continue
        d3 = work(l)
        r1 = subdia(l)
        do 23 j=1,kl
          d1 = work(j)
          d2 = subdia(j)
          vecr(j,i) = (d1*d3+d2*r1)/r
          veci(j,i) = (d2*d3-d1*r1)/r
          vecr(j,i-1) = vecr(j,i)
   23     veci(j,i-1) =-veci(j,i)
   24   continue
c
      return
      end
      
      subroutine hesqr(a,h,evr,evi,subdia,indic,eps,ex)
      implicit none
      include 'newmpar.h'
c
c the following real variables were initially single prec.-
c subdia, eps, ex, r, shift
      integer, dimension(kl) :: indic
      integer i,j,k,l,m
      integer m1,maxst,ns
      real, dimension(kl,kl) :: a,h
      real, dimension(kl) :: evr,evi,subdia
      real eps,ex
      real sr,sr2,shift
      real r,s,t,x,y,z
c this sub.routine finds all the eigenvalues of a real
c general matrix. the original matrix a of order n is
c reduced to the upper-hessenberg form h by means of
c similarity transformations(householder method). the matrix
c h is preserved in the upper half of the array h and in the
c array subdia.  the special vectors used in the definition
c of the householder transformation matrices are stored in
c the lower part of the array h.
c nm is the first dimension of the arrays a and h. nm must
c be equal to or greater than n.
c the real parts of the n eigenvalues will be found in the
c first n places of the array evr,and
c the imaginary parts in the first n places of the array evi
c the array indic indicates the success of the routine as
c follows
c     value of indic(i)  eigenvalue i
c            0             not found
c            1               found
c eps is a small positive number that numerically represents
c zero in the program. eps = (euclidian norm of h)*ex ,where
c ex = 2**(-t). t is the number of binary digits in the
c mantissa of a floating point number.
c
c
c
c reduction of the matrix a to an upper-hessenberg form h.
c there are n-2 steps.
      if(kl-2)14,1,2
    1 subdia(1) = a(2,1)
      go to 14
    2 m = kl-2
      do 12 k=1,m
        l = k+1
        s = 0.0
        do 3 i=l,kl
          h(i,k) = a(i,k)
3     s=s+abs(a(i,k))
      if(s.ne.abs(a(k+1,k)))go to 4
        subdia(k) = a(k+1,k)
        h(k+1,k) = 0.0
        go to 12
    4   sr2 = 0.0
        do 5 i=l,kl
          sr = a(i,k)
          sr = sr/s
          a(i,k) = sr
    5     sr2 = sr2+sr*sr
        sr = sqrt(sr2)
        if(a(l,k).lt.0.0)go to 6
        sr = -sr
    6   sr2 = sr2-sr*a(l,k)
        a(l,k) = a(l,k)-sr
        h(l,k) = h(l,k)-sr*s
        subdia(k) = sr*s
        x = s*sqrt(sr2)
        do 7 i=l,kl
          h(i,k) =h(i,k)/x
    7     subdia(i) = a(i,k)/sr2
c premultiplication by the matrix pr.
          do 9 j=l,kl
            sr = 0.0
            do 8 i=l,kl
    8         sr = sr+a(i,k)*a(i,j)
            do 9 i=l,kl
    9         a(i,j) = a(i,j)-subdia(i)*sr
c postmultiplication by the matrix pr.
            do 11 j=1,kl
              sr=0.0
              do 10 i=l,kl
   10           sr = sr+a(j,i)*a(i,k)
              do 11 i=l,kl
   11           a(j,i) = a(j,i)-subdia(i)*sr
   12       continue
      do 13 k=1,m
   13   a(k+1,k) = subdia(k)
c transer of the upper half of the matrix a into the
c array h and the calculation of the small positive number
c eps.
      subdia(kl-1) = a(kl,kl-1)
   14 eps = 0.0
      do 15 k=1,kl
        indic(k) = 0
        if(k.ne.kl)eps = eps+subdia(k)**2
        do 15 i=k,kl
          h(k,i) = a(k,i)
   15     eps = eps + a(k,i)**2
      eps = ex*sqrt(eps)
c
c the qr iterative process. the upper-hessenberg matrix h is
c reduced to the upper-modified triangular form.
c
c determination of the shift of origin for the first step of
c the qr iterative process.
      shift = a(kl,kl-1)
      if(kl.le.2)shift = 0.0
      if(a(kl,kl).ne.0.0)shift = 0.0
      if(a(kl-1,kl).ne.0.0)shift = 0.0
      if(a(kl-1,kl-1).ne.0.0)shift = 0.0
      m = kl
      ns= 0
      maxst = kl*10
c
c testing if the upper half of the matrix is equal to zero.
c if it is equal to zero the qr process is not necessary.
      do 16 i=2,kl
      do 16 k=i,kl
      if(a(i-1,k).ne.0.0)go to 18
   16 continue
      do 17 i=1,kl
      indic(i)=1

      evr(i) = a(i,i)
   17 evi(i) = 0.0
      go to 37
c
c start the main loop of the qr process.
   18 k=m-1
      m1=k
      i = k
c find any decompositions of the matrix.
c jump to 34 if the last submatrix of the decomposition is
c of the order one.
c jump to 35 if the last submatrix of the decomposition is
c of the order two.
      if(k)37,34,19
19    if(abs(a(m,k)).le.eps)go to 34
      if(m-2.eq.0)go to 35
   20 i = i-1
      if(abs(a(k,i)).le.eps)go to 21
      k = i
      if(k.gt.1)go to 20
   21 if(k.eq.m1)go to 35
c transformation of the matrix of the order greater than two
      s = a(m,m)+a(m1,m1)+shift
      sr= a(m,m)*a(m1,m1)-a(m,m1)*a(m1,m)+0.25*shift**2
      a(k+2,k) = 0.0
c calculate x1,y1,z1,for the submatrix obtained by the
c decomposition
      x = a(k,k)*(a(k,k)-s)+a(k,k+1)*a(k+1,k)+sr
      y = a(k+1,k)*(a(k,k)+a(k+1,k+1)-s)
      r = abs(x)+abs(y)
      if(r.eq.0.0)shift = a(m,m-1)
      if(r.eq.0.0)go to 21
      z = a(k+2,k+1)*a(k+1,k)
      shift = 0.0
      ns = ns+1
c
c the loop for one step of the qr process.
      do 33 i=k,m1
      if(i.eq.k)go to 22
c calculate xr,yr,zr.
      x = a(i,i-1)
      y = a(i+1,i-1)
      z = 0.0
      if(i+2.gt.m)go to 22
      z = a(i+2,i-1)
   22 sr2 = abs(x)+abs(y)+abs(z)
      if(sr2.eq.0.0)go to 23
      x = x/sr2
      y = y/sr2
      z = z/sr2
c     print *,'x,y,z,sr2: ',x,y,z,sr2
c     following line for p.c.version
c     if(abs(z).lt.1.e-19)z=0.
   23 s = sqrt(x*x + y*y + z*z)
      if(x.lt.0.0)go to 24
      s = -s
   24 if(i.eq.k)go to 25
      a(i,i-1) = s*sr2
   25 if(sr2.ne.0.0)go to 26
      if(i+3.gt.m)go to 33
      go to 32
   26 sr = 1.0-x/s
      s = x-s
      x = y/s
      y = z/s
c premultiplication by the matrix pr.
      do 28 j=i,m
      s = a(i,j)+a(i+1,j)*x
      if(i+2.gt.m)go to 27
      s = s+a(i+2,j)*y
   27 s = s*sr
      a(i,j) = a(i,j)-s
      a(i+1,j) = a(i+1,j)-s*x
      if(i+2.gt.m)go to 28
      a(i+2,j) = a(i+2,j)-s*y
   28 continue
c postmultiplication by the matrix pr.
      l = i+2
      if(i.lt.m1)go to 29
      l = m
   29 do 31 j=k,l
      s = a(j,i)+a(j,i+1)*x
      if(i+2.gt.m)go to 30
      s = s + a(j,i+2)*y
   30 s = s*sr
      a(j,i) = a(j,i)-s
      a(j,i+1)=a(j,i+1)-s*x
      if(i+2.gt.m)go to 31
      a(j,i+2)=a(j,i+2)-s*y
   31 continue
      if(i+3.gt.m)go to 33
      s = -a(i+3,i+2)*y*sr
   32 a(i+3,i) = s
c     print *,'s,x: ',s,x
      a(i+3,i+1) = s*x
      a(i+3,i+2) = s*y + a(i+3,i+2)
   33 continue
c
      if(ns.gt.maxst)go to 37
      go to 18
c
c compute the last eigenvalue.
   34 evr(m) = a(m,m)
      evi(m) = 0.0
      indic(m) = 1
      m = k
      go to 18
c
c compute the eigenvalues of the last 2x2 matrix obtained by
c the decomposition.
   35 r = 0.5*(a(k,k)+a(m,m))
      s = 0.5*(a(m,m)-a(k,k))
      s = s*s + a(k,m)*a(m,k)
      indic(k) = 1
      indic(m) = 1
      if(s.lt.0.0)go to 36
      t = sqrt(s)
      evr(k) = r-t
      evr(m) = r+t
      evi(k) = 0.0
      evi(m) = 0.0
      m = m-2
      go to 18
   36 t = sqrt(-s)
      evr(k) = r
      evi(k) = t
      evr(m) = r
      evi(m) = -t
      m = m-2
      go to 18
c
   37 return
      end
      
      subroutine matinv(a,b,l,d,irror)
      implicit none
      include 'newmpar.h'
      integer, dimension(kl,2) :: ind
      integer, dimension(kl) :: ipiv
      integer l,irror
      integer i,j,k,m
      integer irow,icol
      real, dimension(kl,kl) :: a,b
      real d
      real amax
c     a is an nxn matrix to be inverted,or containing equation coeffs
c     b is an nxm rhs matrix for equations
c     if l=0,inverse only given.l positive,solutions only.l negative
c      both.   m=abs(l).
c     d contains the determinant of the a matrix on exit
c     a is replaced by the inverse ,b by the solutions.
c     method of gauss-jordon pivotal elimination
      m=iabs(l)
      d=1.0
      do i=1,kl
        ipiv(i)=0
      end do
      do 220 i=1,kl
      amax=0.0
c     search sub-matrix for largest element as pivot
      do  70 j=1,kl
      if(ipiv(j)) 80,30,70
   30 do  60 k=1,kl
      if(ipiv(k)-1) 40,60,80
c     this row column has not been a pivot
   40 if(abs(a(j,k))-amax)60,60,50
   50 irow=j
      icol=k
      amax=abs(a(j,k))
   60 continue
   70 continue
c     pivot found
      ipiv(icol)=ipiv(icol)+1
      if(amax.gt.1.0e-20)go to 90
c     matrix singular,error return
   80 irror=1
      return
   90 if(irow-icol) 95,130,95
c     make pivot a diagonal element by row interchange.
   95 d=-d
      do 100 k=1,kl
      amax=a(irow,k)
      a(irow,k)=a(icol,k)
  100 a(icol,k)=amax
      if(m) 130,130,110
  110 do 120 k=1,m
      amax=b(irow,k)
      b(irow,k)=b(icol,k)
  120 b(icol,k)=amax
  130 ind(i,1)=irow
      ind(i,2)=icol
      amax=a(icol,icol)
      d=d*amax
      a(icol,icol)=1.0
c     divide pivot row by pivot
      do 140 k=1,kl
  140 a(icol,k)=a(icol,k)/amax
      if(m) 170,170,150
  150 do 160 k=1,m
  160 b(icol,k)=b(icol,k)/amax
c     reduce non-pivot rows
  170 do 220 j=1,kl
      if(j-icol) 180,220,180
  180 amax=a(j,icol)
      a(j,icol)=0.0
      do 190  k=1,kl
  190 a(j,k)=a(j,k)-a(icol,k)*amax
      if(m)  220,220,200
  200 do  210 k=1,m
  210 b(j,k)=b(j,k)-b(icol,k)*amax
  220 continue
c     after n pivotal condensations,solutions lie in b matrix
      if(l) 230,230,270
c     for inverse of a, interchange columns
  230 do 260 i=1,kl
      j=kl+1-i
      if(ind(j,1)-ind(j,2)) 240,260,240
  240 irow=ind(j,1)
      icol=ind(j,2)
      do 250  k=1,kl
      amax=a(k,irow)
      a(k,irow)=a(k,icol)
  250 a(k,icol)=amax
  260 continue
  270 irror=0
      return
      end
      
      subroutine realve(m,ivec,a,vecr,evr,evi,
     1 iwork,work,indic,eps,ex)
      implicit none
      include 'newmpar.h'
      integer, dimension(kl) :: iwork,indic
      integer m
      integer i,j,k,l
      integer ivec,iter,ns
      real, dimension(kl,kl) :: a,vecr
      real, dimension(kl) :: evr,evi
      real, dimension(kl) :: work
      real eps,ex
      real r,r1,t,evalue,previs
      real s,sr,bound

c the following real variables were initially single-
c bound,eps,evalue,ex,previs,r,r1,work
c this sub.routine finds the real eigenvector of the real
c upper-hessenberg matrix in the array a,corresponding to
c the real eigenvalue stored in evr(ivec). the inverse
c iteration method is used.
c note the matrix in a is destroyed by the sub.routine.
c n is the order of the upper-hessenberg matrix.
c nm defines the first dimension of the two dimensional
c arrays a and vecr. nm must be equal to or greater than n.
c m is the order of the submatrix obtained by a suitable
c decomposition of the upper-hessenberg matrix if some
c subdiagonal elements are equal to zero. the value of m is
c chosen so that the last n-m components of the eigenvector
c are zero.
c ivec gives the position of the eigenvalue in the array evr
c for which the corresponding eigenvector is computed.
c the array evi would contain the imaginary parts of the n
c eigenvalues if they existed.
c
c the m components of the computed real eigenvector will be
c found in the first m places of the column ivec of the two
c dimensional array vecr.
c
c iwork and work are the working stores used during the
c gaussian elimination and backsubstitution process.
c the array indic indicates the success of the routine as
c follows
c     value of indic(i)   eigenvector i
c            1             not found
c            2               found
c eps is a small positive number that numerically represents
c zero in the program. eps = (euclidian norm of a)*ex,where
c ex = 2**(-t). t is the number of binary digits in the
c mantissa of a floating point number.

      vecr(1,ivec) = 1.0
      if(m.eq.1)go to 24
c small perturbation of equal eigenvalues to obtain a full
c set of eigenvectors.
      evalue = evr(ivec)
      if(ivec.eq.m)go to 2
       k = ivec+1
      r = 0.0
      do 1 i=k,m
      if(evalue.ne.evr(i))go to 1
      if(evi(i).ne.0.0)go to 1
      r = r+3.0
    1 continue
      evalue = evalue+r*ex
    2 do 3 k=1,m
    3 a(k,k) = a(k,k)-evalue
c
c gaussian elimination of the upper-hessenberg matrix a. all
c row interchanges are indicated in the array iwork.all the
c multipliers are stored as the subdiagonal elements of a.
      k = m-1
      do 8 i=1,k
      l = i+1
      iwork(i) = 0
      if(a(i+1,i).ne.0.0)go to 4
      if(a(i,i).ne.0.0)go to 8
      a(i,i) = eps
      go to 8
4     if(abs(a(i,i)).ge.abs(a(i+1,i)))go to 6
      iwork(i) = 1
      do 5 j=i,m
      r = a(i,j)
      a(i,j) = a(i+1,j)
    5 a(i+1,j) = r
    6 r = -a(i+1,i)/a(i,i)
      a(i+1,i) = r
      do 7 j=l,m
    7 a(i+1,j) = a(i+1,j)+r*a(i,j)
    8 continue
      if(a(m,m).ne.0.0)go to 9
      a(m,m) = eps
c
c the vector (1,1,...,1) is stored in the place of the right
c hand side column vector.
    9 do 11 i=1,kl
      if(i.gt.m)go to 10
      work(i) = 1.0
      go to 11
   10 work(i) = 0.0
   11 continue
c
c the inverse iteration is performed on the matrix until the
c infinite norm of the right-hand side vector is greater
c than the bound defined as  0.01(n*ex).
      bound = 0.01/(ex * float(kl))
      ns = 0
      iter = 1
c
c the backsubstitution.
   12 r = 0.0
      do 15 i=1,m
      j = m-i+1
      s = work(j)
      if(j.eq.m)go to 14
      l = j+1
      do 13 k=l,m
      sr = work(k)
   13 s = s - sr*a(j,k)
   14 work(j) = s/a(j,j)
      t = abs(work(j))
      if(r.ge.t)go to 15
      r = t
   15 continue
c
c the computation of the right-hand side vector for the new
c iteration step.
      do 16 i=1,m
   16 work(i) = work(i)/r
c
c the computation of the residuals and comparison of the
c residuals of the two successive steps of the inverse
c iteration. if the infinite norm of the residual vector is
c greater than the infinite norm of the previous residual
c vector the computed eigenvector of the previous step is
c taken as the final eigenvector.
      r1 = 0.0
      do 18 i=1,m
      t = 0.0
      do 17 j=i,m
   17 t = t+a(i,j)*work(j)
      t=abs(t)
      if(r1.ge.t)go to 18
      r1 = t
   18 continue
      if(iter.eq.1)go to 19
      if(previs.le.r1)go to 24
   19 do 20 i=1,m
   20 vecr(i,ivec) = work(i)
      previs = r1
      if(ns.eq.1)go to 24
      if(iter.gt.6)go to 25
      iter = iter+1
      if(r.lt.bound)go to 21
      ns = 1
c
c gaussian elimination of the right-hand side vector.
   21 k = m-1
      do 23 i=1,k
      r = work(i+1)
      if(iwork(i).eq.0)go to 22
      work(i+1)=work(i)+work(i+1)*a(i+1,i)
      work(i) = r
      go to 23
   22 work(i+1)=work(i+1)+work(i)*a(i+1,i)
   23 continue
      go to 12
c
   24 indic(ivec) =2
   25 if(m.eq.kl)go to 27
      j = m+1
      do 26 i=j,kl
   26 vecr(i,ivec) = 0.0
   27 return
      end
      
      subroutine scaler(a,h,prfact,enorm)
      implicit none
      include 'newmpar.h'
c the following real variables were initially single prec.-
c bound1,bound2,enorm
      integer i,j
      integer iter,ncount
      real, dimension(kl,kl) :: a,h
      real, dimension(kl) :: prfact
      real enorm
      real fnorm,column,row
      real bound1,bound2,q,factor
c this sub.routine stores the matrix of the order n from the
c array a into the array h. afterward the matrix in the
c array a is scaled so that the quotient of the absolute sum
c of the off-diagonal elements of column i and the absolute
c sum of the off-diagonal elements of row i lies within the
c values of bound1 and bound2.
c the component i of the eigenvector obtained by using the
c scaled matrix must be divided by the value found in the
c prfact(i) of the array prfact. in this way the eigenvector
c of the non-scaled matrix is obtained.
c
c after the matrix is scaled it is normalised so that the
c value of the euclidian norm is equal to one.
c if the process of scaling was not successful the original
c matrix from the array h would be stored back into a and
c the eigenproblem would be solved by using this matrix.
c nm defines the first dimension of the arrays a and h. nm
c must be greater or equal to n.
c the eigenvalues of the normalised matrix must be
c multiplied by the scalar enorm in order that they become
c the eigenvalues of the non-normalised matrix.
      do i=1,kl
        do j=1,kl
          h(i,j) = a(i,j)
        end do
        prfact(i)= 1.0
      end do
      bound1 = 0.75
      bound2 = 1.33
      iter = 0
    3 ncount = 0
      do 8 i=1,kl
      column = 0.0
      row    = 0.0
        do 4 j=1,kl
          if(i.eq.j)go to 4
      column=column+abs(a(j,i))
      row   =row   +abs(a(i,j))
    4     continue
        if(column.eq.0.0)go to 5
        if(row.eq.0.0)go to 5
        q = column/row
        if(q.lt.bound1)go to 6
        if(q.gt.bound2)go to 6
    5   ncount = ncount + 1
        go to 8
    6   factor = sqrt(q)
        do 7 j=1,kl
          if(i.eq.j)go to 7
          a(i,j) = a(i,j)*factor
          a(j,i) = a(j,i)/factor
    7     continue
        prfact(i) = prfact(i)*factor
    8   continue
      iter = iter+1
      if(iter.gt.30)go to 11
      if(ncount.lt.kl)go to 3
c
      fnorm = 0.0
      do 9 i=1,kl
        do 9 j=1,kl
          q = a(i,j)
    9     fnorm = fnorm+q*q
      fnorm = sqrt(fnorm)
      do 10 i=1,kl
        do 10 j=1,kl
   10     a(i,j)=a(i,j)/fnorm
      enorm = fnorm
      go to 13
c
   11 do 12 i=1,kl
c
c modification suggested by referee in a.c.m.certification
c
        prfact(i)=1.0
        do 12 j=1,kl
   12     a(i,j) = h(i,j)
      enorm = 1.0
c
   13 return
      end
      
      subroutine sigtosigh(sig,sigmh,kl)
      implicit none
      integer, intent(in) :: kl
      integer k
!     these routines are written from top down
      real, dimension(kl), intent(in) :: sig
      real, dimension(kl+1), intent(out) :: sigmh
      sigmh(1)=1.
      sigmh(kl+1)=0.
      do k=1,kl-1
        sigmh(k+1)=.5*(sig(k)+sig(k+1))
      end do
      return
      end
      
      subroutine sightosig(sig,sigmh,kl)
      implicit none
      integer, intent(in) :: kl
      integer k
      real, dimension(kl), intent(out) :: sig
      real, dimension(kl+1), intent(in) :: sigmh
      do k=1,kl
        sig(k)=.5*(sigmh(k+1)+sigmh(k))
      end do
      return
      end
