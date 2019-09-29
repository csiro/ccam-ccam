! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2019 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
! This file is part of the Conformal Cubic Atmospheric Model (CCAM)
!
! CCAM is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! CCAM is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with CCAM.  If not, see <http://www.gnu.org/licenses/>.

!------------------------------------------------------------------------------
    
!  this is eig derived from eignew, but used in-line in C-CAM
subroutine eig(sigin,sigmhin,tbarin,lapsbot,isoth,dtin,epspin,epshin,nsig,betin,betmin,nh)
use cc_mpi, only : myid
use vecs_m, only : emat,einv,bam
use newmpar_m
implicit none
#ifdef pgi
integer, parameter :: r16 = kind(1._8)
#else
integer, parameter :: r16 = kind(1._16)
#endif
integer, intent(in) :: nh,nsig,lapsbot,isoth
integer :: nchng,k,l
integer, parameter :: neig = 1
real, intent(in) :: dtin, epspin, epshin
real, dimension(kl), intent(in) :: sigin,sigmhin
real, dimension(kl), intent(in) :: tbarin,betin,betmin
real(kind=r16) :: tem, dt, epsp, epsh
real(kind=r16), dimension(kl) :: sig,tbar,bet,betm,lbam
real(kind=r16), dimension(kl+1) :: sigmh
real(kind=r16), dimension(kl,kl) :: lemat,leinv

dt = real(dtin,r16)
epsp = real(epspin,r16)
epsh = real(epshin,r16)
sig(:) = real(sigin(:),r16)
sigmh(1:kl) = real(sigmhin(:),r16)
sigmh(kl+1) = 0._r16
tbar = real(tbarin,r16)
bet = real(betin,r16)
betm = real(betmin,r16)
lemat = real(emat,r16)
leinv = real(einv,r16)
lbam = real(bam,r16)
! lapsbot=1 gives zero lowest t lapse for phi calc
if ( myid==0 ) then
  write(6,*) 'this run configured with kl = ',kl
  write(6,*) 'entering eig tbar,lapsbot: ',tbarin(1),lapsbot
  write(6,*) '             isoth,dtin:   ',isoth,dtin
  write(6,*) '             epsp,epsh:    ',epspin,epshin
  write(6,*) '             nh:           ',nh
end if

!     expect data from bottom up
if ( sig(1)<sig(kl) ) then
  call flip3(sig,1,1,kl,1)
  lbam(1:kl) = sigmhin(:)
  call flip3(lbam,1,1,kl,1)  
  sigmh(1) = 0.
  sigmh(2:kl+1) = lbam(1:kl)
endif  ! (sig(1)<sig(kl))
if ( myid==0 ) then
  write(6,*) 'final sig values: ',sig
  write(6,*) 'final sigmh values: ',sigmh
  open(28,file='eigenv.out')
end if
call eigs(isoth,tbar,dt,epsp,epsh,nh,sig,sigmh,bet,betm,lbam,lemat,leinv)
if ( myid==0 ) then
  write(6,*) 'about to write to 28 '
  write(28,*)kl,lapsbot,isoth,nsig,'   kl,lapsbot,isoth,nsig'
end if
!     re-order the eigenvectors if necessary
nchng = 1
do while ( nchng/=0 )
  nchng = 0
  do k = 2,kl
    if ( lbam(k)<lbam(k-1) ) cycle
    nchng = 1
    tem = lbam(k)
    lbam(k) = lbam(k-1)
    lbam(k-1) = tem
    do l = 1,kl
      tem = lemat(l,k)
      lemat(l,k) = lemat(l,k-1)
      lemat(l,k-1) = tem
      tem = leinv(k,l)
      leinv(k,l) = leinv(k-1,l)
      leinv(k-1,l) = tem
    end do
  end do
end do
if ( myid==0 ) then
  write(6,*)'eigenvectors re-ordered'
  write(6,*)'bam',(lbam(k),k=1,kl)
end if
if ( neig==1 ) then
!      write data from bottom up
  if ( sig(1)<sig(kl) ) then
    call flip3( sig,1, 1,kl, 1)
    call flip3( sigmh,1, 1,kl+1, 1)
    call flip3(lemat,1, 1,kl,kl)
    call flip3(leinv,1,kl,kl, 1)
  endif
endif
emat = real(lemat)
einv = real(leinv)
bam = real(lbam)
return
end subroutine eig

subroutine eigs(isoth,tbar,dt,epsp,epsh,nh,sig,sigmh,bet,betm,bam,emat,einv)

use cc_mpi, only : myid
use const_phys
use newmpar_m

implicit none

#ifdef pgi
integer, parameter :: r16 = kind(1._8)
#else
integer, parameter :: r16 = kind(1._16)
#endif
integer isoth, nh
integer k, l, irror
integer, dimension(kl) :: indic
real(kind=r16) dt, epsp, epsh
real(kind=r16) dp !, factg, factr
!     sets up eigenvectors
real(kind=r16), dimension(kl) :: sig, dsig
real(kind=r16), dimension(kl) :: bet, betm, bam
real(kind=r16), dimension(kl) :: evimag, sum1, tbar
real(kind=r16), dimension(kl+1) :: sigmh
real(kind=r16), dimension(kl,kl) :: emat, einv
real(kind=r16), dimension(kl,kl) :: bmat, veci !, gmat
real(kind=r16), dimension(kl,kl) :: aa, ab, ac
real(kind=r16), dimension(kl,kl) :: aaa, cc
      
do k = 1,kl
  dsig(k) = sigmh(k+1) - sigmh(k)
enddo
if ( myid==0 ) then
  write(6,*) 'sigmh ',(sigmh(k),k=1,kl)
  write(6,*) 'dsig ',(dsig(k),k=1,kl)
end if

!get(1) = bet(1)/(rdry*sig(1))
!getm(1) = 0.
!do k = 2,kl
!  get(k) = bet(k)/(rdry*sig(k))
!  getm(k) = betm(k)/(rdry*sig(k-1))
!enddo      
!factg = 2./(dt*(1.+epsp))      
!factr = factg*rdry*rdry*tbar(1)*tbar(1)/(grav*grav)
      
bmat(:,:) = 0._r16  ! N.B. bmat includes effect of r/sig weighting
!gmat(:,:) = 0.  ! N.B. gmat may include effect of 1/sig**2 weighting
do k = 2,kl
  do l = 1,k-1
    bmat(k,l) = bet(l)+betm(l+1)
    !gmat(k,l) = factr*(get(l)+getm(l+1))
  enddo ! l loop
enddo  ! k loop
do k = 1,kl
  bmat(k,k) = bet(k)
  !gmat(k,k) = factr*get(k)
enddo

if ( myid==0 ) then
  write(6,*)'dt,epsp,tbar ',dt,epsp,tbar(1)
  write(6,*)'bet ',bet
  !write(6,*) 'get ',get
  !write(6,*) 'getm ',getm
  !     even newer derivation section
  write(6,*) 'even newer derivation section'
end if

do l = 1,kl
  do k = 1,kl
    ab(k,l) = dsig(l)
    ac(k,l) = dsig(l)
  end do
end do
aa(:,:) = 0._r16
do k = 1,kl
  do l = k,kl
    aa(k,l) = -rdry*tbar(1)*dsig(l)/(cp*sig(k))
  end do
end do
do k = 1,kl
  aa(k,k) = -rdry*tbar(1)*(sigmh(k+1)-sig(k))/(cp*sig(k))
  ac(k,k) = ac(k,k) + 1._r16
end do

if ( isoth==1 ) then  !  extra vadv terms added
  aa(:,:) = aa(:,:) + tbar(1)*ac(:,:)
end if

if ( nh>0 ) then
  ! non-hydrostatic
  bmat(:,:) = bmat(:,:)*real(1.+4.*cp*tbar(1)/((grav*dt)**2*(1.+epsp)*(1.+epsh)),r16)
end if  ! (nh<=0) ..else..
aaa(:,:) = matmul(bmat(:,:), aa(:,:))
cc(:,:) = aaa(:,:) - rdry*tbar(1)*ab(:,:)
      
aaa(:,:) = cc(:,:)
call eigenp(aaa,bam,evimag,emat,veci,indic)
if ( myid==0 ) then
  write(6,*) 'bam',(bam(k),k=1,kl)
end if
call matinv(cc,sum1,0,dp,irror)
einv(:,:) = emat(:,:)
call matinv(einv,sum1,0,dp,irror)

return
end subroutine eigs

subroutine flip3(a,il,jl,kl,ll)
implicit none
#ifdef pgi
integer, parameter :: r16 = kind(1._8)
#else
integer, parameter :: r16 = kind(1._16)
#endif
integer, intent(in) :: il,jl,kl,ll
integer l,j,i,k
real(kind=r16), dimension(il,jl,kl,ll), intent(inout) :: a
real(kind=r16) tem
do l = 1,ll
  do j = 1,jl
    do i = 1,il
      do k = 1,kl/2
        tem = a(i,j,k,l)
        a(i,j,k,l) = a(i,j,kl+1-k,l)
        a(i,j,kl+1-k,l) = tem
      end do
    end do
  end do
end do
return
end subroutine flip3

subroutine eigenp(a,evr,evi,vecr,veci,indic)
use newmpar_m
implicit none
#ifdef pgi
integer, parameter :: r16 = kind(1._8)
#else
integer, parameter :: r16 = kind(1._16)
#endif
integer, dimension(kl*kl) :: iwork,local
integer, dimension(kl) :: indic
integer i,j,k,l,m
integer k1,l1,kon,ivec
real(kind=r16), dimension(:), allocatable :: prfact,subdia,work ! use allocatable for cray compiler
real(kind=r16), dimension(kl,kl) :: a,vecr,veci
real(kind=r16), dimension(kl) :: evr,evi
real(kind=r16) d1,d2,d3,enorm
real(kind=r16) r,r1,ex,eps
      
! a.c.m. algorithm number 343
! revised july, 1970 by n.r.pummeroy, dcr, csiro, canberra
! the following variables were changed
! from single to real in eigenp:
! r,r1,enorm,eps,ex,work,work1,work2,subdia
! see also comments for routines scaler, hesqr, realve, compve
!
!
! this sub.routine finds all the eigenvalues and the
! eigenvectors of a real general matrix of order n.
!
! first in the sub.routine scaler the matrix is scaled so that
! the corresponding rows and columns are approximately
! balanced and then the matrix is normalised so that the
! value of the euclidian norm of the matrix is equal to one.
!
! the eigenvalues are computed by the qr double-step method
! in the sub.routine hesqr.
! the eigenvectors are computed by inverse iteration in
! the sub.routine realve,for the real eigenvalues,or in the
! subroutine compve,for the complex eigenvalues.
!
! the elements of the matrix are to be stored in the first n
! rows and columns of the two dimensional array a. the
! original matrix is destroyed by the sub.routine.
! n is the order of the matrix.
! nm defines the first dimension of the two dimensional
! arrays a,vecr,veci and the dimension of the one
! the real parts of the n computed eigenvalues will be found
! in the first n places of the array evr and the imaginary
! parts in the first n places of the array evi.
! the real components of the normalised eigenvector i
! (i=1,2,...,n) corresponding to the eigenvalue stored in
! evr(i) and evi(i) will be found in the first n places of
! the column i of the two dimensional array vecr and the
! imaginary components in the first n places of the column i
! of the two dimensional array veci.
!
! the real eigenvector is normalised so that the sum of the
! squares of the components is equal to one.
! the complex eigenvector is normalised so that the
! component with the largest value in modulus has its real
! part equal to one and the imaginary part equal to zero.
!
! the array indic indicates the success of the sub.routine
! eigenp as follows
!     value of indic(i)   eigenvalue i   eigenvector i
!            0              not found      not found
!            1              found          not found
!            2              found          found


allocate( prfact(kl), subdia(kl), work(kl) )

call scaler(a,veci,prfact,enorm)
! the computation of the eigenvalues of the normalised
! matrix.
!  take t=50 significant binary figures.  ex=2**(-t)
!     ex=8.88178418e-16
!  following for 60 binary figures:
ex = 8.674e-19_r16
call hesqr(a,veci,evr,evi,subdia,indic,eps,ex)

! the possible decomposition of the upper-hessenberg matrix
! into the submatrices of lower order is indicated in the
! array local. the decomposition occurs when some
! subdiagonal elements are in modulus less than a small
! positive number eps defined in the sub.routine hesqr . the
! amount of work in the eigenvector problem may be
! minimised in this way.
i = 1
local(1) = 1
do j = kl,2,-1
   if ( abs(subdia(j-1))<=eps ) then
     i = i + 1
     local(i) = 0
   end if
  local(i) = local(i) + 1
end do

! the eigenvector problem.
k = 1
kon = 0
l = local(1)
m = kl
do i = 1,kl
  ivec = kl - i + 1
  if ( i>l ) then
    k = k + 1
    m = kl - l
    l = l + local(k)
  end if
  if ( indic(ivec)/=0 ) then
    if ( abs(evi(ivec))<1.e-99_r16 ) then
! transfer of an upper-hessenberg matrix of the order m from
! the arrays veci and subdia into the array a.
      do l1 = 1,m
        a(1,l1) = veci(1,l1)
      end do
      do k1 = 2,m
        do l1 = k1,m
          a(k1,l1) = veci(k1,l1)
        end do
        a(k1,k1-1) = subdia(k1-1)
      end do

! the computation of the real engenvector ivec of the upper-
! hessenberg matrix corresponding to the real eigenvalue
! evr(ivec).
      call realve(m,ivec,a,vecr,evr,evi,iwork,work,indic,eps,ex)

    else

! the computation of the complex eigenvector ivec of the
! upper-hessenberg matrix corresponding to the complex
! eigenvalue evr(ivec) + i*evi(ivec). if the value of kon is
! not equal to zero then this complex eigenvector has
! already been found from its conjugate.
      if(kon==0) then
        kon = 1
        write(6,*) 'attempted call to comove'
        stop
      end if
      kon = 0
    end if
  end if
end do

! the reconstruction of the matrix used in the reduction of
! matrix a to an upper-hessenberg form by householder method
do i = 1,kl
  do j = i,kl
    a(i,j) = 0.0_r16  
    a(j,i) = 0.0_r16
  end do
  a(i,i) = 1.0_r16
end do
m = kl-2
do k = 1,m
  l = k + 1
  do j = 2,kl
    d1 = 0.0_r16
    do i = l,kl
      d2 = veci(i,k)
      d1 = d1+ d2*a(j,i)
    end do
    do i = l,kl
      a(j,i) = a(j,i) - veci(i,k)*d1
    end do
  end do
end do

! the computation of the eigenvectors of the original non-
! scaled matrix.
kon = 1
do i = 1,kl
  l = 0
  if ( abs(evi(i))>1.e-99_r16 ) then
    l = 1
    if ( kon/=0 ) then
      kon = 0
      cycle
    end if
  end if
  do j = 1,kl
    d1 = 0.0_r16
    d2 = 0.0_r16
    do k = 1,kl
      d3 = a(j,k)
      d1 = d1+d3*vecr(k,i)
      if ( l/=0 ) then
        d2 = d2 + d3*vecr(k,i-1)
      end if
    end do
    work(j) = d1/prfact(j)
    if ( l/=0 ) then
      subdia(j) = d2/prfact(j)
    end if
  end do

! the normalisation of the eigenvectors and the computation
! of the eigenvalues of the original non-normalised matrix.
  if(l/=1) then
    d1 = 0.0_r16
    do m = 1,kl
      d1 = d1 + work(m)**2
    end do
    d1 = sqrt(d1)
    do m = 1,kl
      veci(m,i) = 0.0_r16
      vecr(m,i) = work(m)/d1
    end do
    evr(i) = evr(i)*enorm
     
  else

    kon = 1
    evr(i) = evr(i)*enorm
    evr(i-1) = evr(i)
    evi(i) = evi(i)*enorm
    evi(i-1) =-evi(i)
    r = 0.0_r16
    do j = 1,kl
      r1 = work(j)**2 + subdia(j)**2
      if ( r<r1 ) then
        r = r1
        l = j
      end if
    end do
    d3 = work(l)
    r1 = subdia(l)
    do j = 1,kl
      d1 = work(j)
      d2 = subdia(j)
      vecr(j,i) = (d1*d3+d2*r1)/r
      veci(j,i) = (d2*d3-d1*r1)/r
      vecr(j,i-1) = vecr(j,i)
      veci(j,i-1) =-veci(j,i)
    end do
  end if
end do

deallocate( prfact, subdia, work )

return
end subroutine eigenp
      
subroutine hesqr(a,h,evr,evi,subdia,indic,eps,ex)
use newmpar_m
implicit none

#ifdef pgi
integer, parameter :: r16 = kind(1._8)
#else
integer, parameter :: r16 = kind(1._16)
#endif
! the following real variables were initially single prec.-
! subdia, eps, ex, r, shift
integer, dimension(kl) :: indic
integer i,j,k,l,m
integer m1,maxst,ns
real(kind=r16), dimension(kl,kl) :: a,h
real(kind=r16), dimension(kl) :: evr,evi,subdia
real(kind=r16) eps,ex
real(kind=r16) sr,sr2,shift
real(kind=r16) r,s,t,x,y,z
logical fflag
! this sub.routine finds all the eigenvalues of a real
! general matrix. the original matrix a of order n is
! reduced to the upper-hessenberg form h by means of
! similarity transformations(householder method). the matrix
! h is preserved in the upper half of the array h and in the
! array subdia.  the special vectors used in the definition
! of the householder transformation matrices are stored in
! the lower part of the array h.
! nm is the first dimension of the arrays a and h. nm must
! be equal to or greater than n.
! the real parts of the n eigenvalues will be found in the
! first n places of the array evr,and
! the imaginary parts in the first n places of the array evi
! the array indic indicates the success of the routine as
! follows
!     value of indic(i)  eigenvalue i
!            0             not found
!            1               found
! eps is a small positive number that numerically represents
! zero in the program. eps = (euclidian norm of h)*ex ,where
! ex = 2**(-t). t is the number of binary digits in the
! mantissa of a floating point number.
!
!
!
! reduction of the matrix a to an upper-hessenberg form h.
! there are n-2 steps.
m = kl-2
do k=1,m
  l = k+1
  s = 0.0_r16
  do i=l,kl
    h(i,k) = a(i,k)
    s=s+abs(a(i,k))
  end do
  if(abs(s-abs(a(k+1,k)))<1.e-99_8)then
    subdia(k) = a(k+1,k)
    h(k+1,k) = 0.0_r16
  else
    sr2 = 0.0_r16
    do i=l,kl
      sr = a(i,k)
      sr = sr/s
      a(i,k) = sr
      sr2 = sr2+sr*sr
    end do
    sr = sqrt(sr2)
    if(a(l,k)>=0.0_r16)then
      sr = -sr
    end if
    sr2 = sr2-sr*a(l,k)
    a(l,k) = a(l,k)-sr
    h(l,k) = h(l,k)-sr*s
    subdia(k) = sr*s
    x = s*sqrt(sr2)
    do i=l,kl
      h(i,k) =h(i,k)/x
      subdia(i) = a(i,k)/sr2
    end do
! premultiplication by the matrix pr.
    do j=l,kl
      sr = sum(a(l:kl,k)*a(l:kl,j))
      do i=l,kl
        a(i,j) = a(i,j)-subdia(i)*sr
      end do
    end do
! postmultiplication by the matrix pr.
    do j=1,kl
      sr=sum(a(j,l:kl)*a(l:kl,k))
      do i=l,kl
        a(j,i) = a(j,i)-subdia(i)*sr
      end do
    end do
  end if
end do
do k=1,m
  a(k+1,k) = subdia(k)
end do
! transer of the upper half of the matrix a into the
! array h and the calculation of the small positive number
! eps.
subdia(kl-1) = a(kl,kl-1)
eps = 0.0_r16
do k=1,kl
  indic(k) = 0
  if(k/=kl)eps = eps+subdia(k)**2
  do i=k,kl
    h(k,i) = a(k,i)
    eps = eps + a(k,i)**2
  end do
end do
eps = ex*sqrt(eps)

! the qr iterative process. the upper-hessenberg matrix h is
! reduced to the upper-modified triangular form.
!
! determination of the shift of origin for the first step of
! the qr iterative process.
shift = a(kl,kl-1)
if(abs(a(kl,kl))>1.e-99_r16)shift = 0.0_r16
if(abs(a(kl-1,kl))>1.e-99_r16)shift = 0.0_r16
if(abs(a(kl-1,kl-1))>1.e-99_r16)shift = 0.0_r16
m = kl
ns= 0
maxst = kl*10

! testing if the upper half of the matrix is equal to zero.
! if it is equal to zero the qr process is not necessary.
fflag=.true.
do i=2,kl
  do k=i,kl
    if(abs(a(i-1,k))>1.e-99_r16) fflag=.false.
  end do
end do

if (fflag) then
  do i=1,kl
    indic(i)=1
    evr(i) = a(i,i)
    evi(i) = 0.0_r16
  end do
else

! start the main loop of the qr process.
  do while(ns<=maxst)
    k=m-1
    m1=k
    i = k
! find any decompositions of the matrix.
    if(k<0) then
      return
    else if (k==0.or.abs(a(m,max(k,1)))<=eps) then
      ! compute the last eigenvalue.
      evr(m) = a(m,m)
      evi(m) = 0.0_r16
      indic(m) = 1
      m = k
      cycle
    end if
    if(m==2)then
      ! compute the eigenvalues of the last 2x2 matrix obtained by
      ! the decomposition.
      r = 0.5*(a(k,k)+a(m,m))
      s = 0.5*(a(m,m)-a(k,k))
      s = s*s + a(k,m)*a(m,k)
      indic(k) = 1
      indic(m) = 1
      if(s<0.0_r16)then
        t = sqrt(-s)
        evr(k) = r
        evi(k) = t
        evr(m) = r
        evi(m) = -t
      else
        t = sqrt(s)
        evr(k) = r-t
        evr(m) = r+t
        evi(k) = 0.0_r16
        evi(m) = 0.0_r16
      end if
      m = m-2
      cycle
    end if
    do while (k>1)  
      i = i-1
      if(abs(a(k,i))<=eps)exit
      k = i
    end do
    if(k==m1)then
      ! compute the eigenvalues of the last 2x2 matrix obtained by
      ! the decomposition.
      r = 0.5*(a(k,k)+a(m,m))
      s = 0.5*(a(m,m)-a(k,k))
      s = s*s + a(k,m)*a(m,k)
      indic(k) = 1
      indic(m) = 1
      if(s<0.0_r16)then
        t = sqrt(-s)
        evr(k) = r
        evi(k) = t
        evr(m) = r
        evi(m) = -t
      else
        t = sqrt(s)
        evr(k) = r-t
        evr(m) = r+t
        evi(k) = 0.0_r16
        evi(m) = 0.0_r16
      end if
      m = m-2
      cycle
    end if
    r = 0.0_r16
    do while(abs(r)<1.e-99_r16)
! transformation of the matrix of the order greater than two
      s = a(m,m)+a(m1,m1)+shift
      sr= a(m,m)*a(m1,m1)-a(m,m1)*a(m1,m)+0.25*shift**2
      a(k+2,k) = 0.0_r16
! calculate x1,y1,z1,for the submatrix obtained by the
! decomposition
      x = a(k,k)*(a(k,k)-s)+a(k,k+1)*a(k+1,k)+sr
      y = a(k+1,k)*(a(k,k)+a(k+1,k+1)-s)
      r = abs(x)+abs(y)
      if(abs(r)<1.e-99_r16)shift = a(m,m-1)
    end do
    z = a(k+2,k+1)*a(k+1,k)
    shift = 0.0_r16
    ns = ns+1

! the loop for one step of the qr process.
    do i=k,m1
      if(i/=k)then
! calculate xr,yr,zr.
        x = a(i,i-1)
        y = a(i+1,i-1)
        z = 0.0_r16
        if(i+2<=m)then
          z = a(i+2,i-1)
        end if
      end if
      sr2 = abs(x)+abs(y)+abs(z)
      if(abs(sr2)>1.e-99_r16)then
        x = x/sr2
        y = y/sr2
        z = z/sr2
      end if
      s = sqrt(x*x + y*y + z*z)
      if(x>=0.0_r16)then
        s = -s
      end if
      if(i/=k)then
        a(i,i-1) = s*sr2
      end if
      if(abs(sr2)<1.e-99_r16)then
        if(i+3<=m)then
          a(i+3,i) = s
          a(i+3,i+1) = s*x
          a(i+3,i+2) = s*y + a(i+3,i+2)
        end if
        cycle
      end if
      sr = 1.0_r16-x/s
      s = x-s
      x = y/s
      y = z/s
! premultiplication by the matrix pr.
      do j=i,m
        s = a(i,j)+a(i+1,j)*x
        if(i+2<=m)then
          s = s+a(i+2,j)*y
        end if
        s = s*sr
        a(i,j) = a(i,j)-s
        a(i+1,j) = a(i+1,j)-s*x
        if(i+2<=m)then
          a(i+2,j) = a(i+2,j)-s*y
        end if
      end do
! postmultiplication by the matrix pr.
      l = i+2
      if(i>=m1)then
        l = m
      end if
      do j=k,l
        s = a(j,i)+a(j,i+1)*x
        if(i+2<=m)then
          s = s + a(j,i+2)*y
        end if
        s = s*sr
        a(j,i) = a(j,i)-s
        a(j,i+1)=a(j,i+1)-s*x
        if(i+2<=m)then
          a(j,i+2)=a(j,i+2)-s*y
        end if
      end do
      if(i+3<=m)then
        s = -a(i+3,i+2)*y*sr
        a(i+3,i) = s
        a(i+3,i+1) = s*x
        a(i+3,i+2) = s*y + a(i+3,i+2)
      end if
    end do

  end do
end if
     
return
end subroutine hesqr
      
subroutine matinv(a,b,l,d,irror)
use newmpar_m
implicit none
#ifdef pgi
integer, parameter :: r16 = kind(1._8)
#else
integer, parameter :: r16 = kind(1._16)
#endif
integer, dimension(kl,2) :: ind
integer, dimension(kl) :: ipiv
integer, intent(in) :: l
integer, intent(out) :: irror
integer i,j,k,m
integer irow,icol
real(kind=r16), dimension(kl,kl), intent(inout) :: a,b
real(kind=r16), intent(out) :: d
real(kind=r16) :: amax
!     a is an nxn matrix to be inverted,or containing equation coeffs
!     b is an nxm rhs matrix for equations
!     if l=0,inverse only given.l positive,solutions only.l negative
!      both.   m=abs(l).
!     d contains the determinant of the a matrix on exit
!     a is replaced by the inverse ,b by the solutions.
!     method of gauss-jordon pivotal elimination
irow=0
icol=0
m=iabs(l)
d=1.0
do i=1,kl
  ipiv(i)=0
end do
do i=1,kl
  amax=0.0_r16
!       search sub-matrix for largest element as pivot
  do j=1,kl
    if (ipiv(j)<0) then
      irror=1
      return
    else if(ipiv(j)==0) then
      do k=1,kl
        if(ipiv(k)>1) then
          irror=1
          return
        else if(ipiv(k)<1) then
!     this row column has not been a pivot
          if(abs(a(j,k))>amax)then
            irow=j
            icol=k
            amax=abs(a(j,k))
          end if
        end if
      end do
    end if
  end do
!       pivot found
  ipiv(icol)=ipiv(icol)+1
  if(amax<=1.0e-20_r16)then
!         matrix singular,error return
    irror=1
    return
  end if
  if(irow/=icol) then
!     make pivot a diagonal element by row interchange.
    d=-d
    do k=1,kl
      amax=a(irow,k)
      a(irow,k)=a(icol,k)
      a(icol,k)=amax
    end do
    do k=1,m
      amax=b(irow,k)
      b(irow,k)=b(icol,k)
      b(icol,k)=amax
    end do
  end if          
  ind(i,1)=irow
  ind(i,2)=icol
  amax=a(icol,icol)
  d=d*amax
  a(icol,icol)=1.0
!       divide pivot row by pivot
  do k=1,kl
    a(icol,k)=a(icol,k)/amax
  end do
  do k=1,m
    b(icol,k)=b(icol,k)/amax
  end do
!       reduce non-pivot rows
  do j=1,kl
    if(j/=icol)then
      amax=a(j,icol)
      a(j,icol)=0.0
      do k=1,kl
        a(j,k)=a(j,k)-a(icol,k)*amax
      end do
      do k=1,m
        b(j,k)=b(j,k)-b(icol,k)*amax
      end do
    end if
  end do
end do
!     after n pivotal condensations,solutions lie in b matrix
if(l<=0)then
!       for inverse of a, interchange columns
  do i=1,kl
    j=kl+1-i
    if(ind(j,1)/=ind(j,2))then
      irow=ind(j,1)
      icol=ind(j,2)
      do k=1,kl
        amax=a(k,irow)
        a(k,irow)=a(k,icol)
        a(k,icol)=amax
      end do
    end if
  end do
end if
irror=0
return
end subroutine matinv
      
subroutine realve(m,ivec,a,vecr,evr,evi,iwork,work,indic,eps,ex)
use newmpar_m
implicit none
#ifdef pgi
integer, parameter :: r16 = kind(1._8)
#else
integer, parameter :: r16 = kind(1._16)
#endif
integer, dimension(kl) :: iwork,indic
integer m
integer i,j,k,l
integer ivec,iter,ns
real(kind=r16), dimension(kl,kl) :: a,vecr
real(kind=r16), dimension(kl) :: evr,evi
real(kind=r16), dimension(kl) :: work
real(kind=r16) eps,ex
real(kind=r16) r,r1,t,evalue,previs
real(kind=r16) s,sr,bound

! the following real variables were initially single-
! bound,eps,evalue,ex,previs,r,r1,work
! this sub.routine finds the real eigenvector of the real
! upper-hessenberg matrix in the array a,corresponding to
! the real eigenvalue stored in evr(ivec). the inverse
! iteration method is used.
! note the matrix in a is destroyed by the sub.routine.
! n is the order of the upper-hessenberg matrix.
! nm defines the first dimension of the two dimensional
! arrays a and vecr. nm must be equal to or greater than n.
! m is the order of the submatrix obtained by a suitable
! decomposition of the upper-hessenberg matrix if some
! subdiagonal elements are equal to zero. the value of m is
! chosen so that the last n-m components of the eigenvector
! are zero.
! ivec gives the position of the eigenvalue in the array evr
! for which the corresponding eigenvector is computed.
! the array evi would contain the imaginary parts of the n
! eigenvalues if they existed.
!
! the m components of the computed real eigenvector will be
! found in the first m places of the column ivec of the two
! dimensional array vecr.
!
! iwork and work are the working stores used during the
! gaussian elimination and backsubstitution process.
! the array indic indicates the success of the routine as
! follows
!     value of indic(i)   eigenvector i
!            1             not found
!            2               found
! eps is a small positive number that numerically represents
! zero in the program. eps = (euclidian norm of a)*ex,where
! ex = 2**(-t). t is the number of binary digits in the
! mantissa of a floating point number.

previs=0._r16

vecr(1,ivec) = 1.0_r16
! small perturbation of equal eigenvalues to obtain a full
! set of eigenvectors.
evalue = evr(ivec)
if(ivec/=m)then
  k = ivec+1
  r = 0.0_r16
  do i=k,m
    if(abs(evalue-evr(i))<1.e-99_r16.and.abs(evi(i))<1.e-99_r16) then
      r = r+3.0_r16
    end if
  end do
  evalue = evalue+r*ex
end if
do k=1,m
  a(k,k) = a(k,k)-evalue
end do

! gaussian elimination of the upper-hessenberg matrix a. all
! row interchanges are indicated in the array iwork.all the
! multipliers are stored as the subdiagonal elements of a.
k = m-1
do i=1,k
  l = i+1
  iwork(i) = 0
  if (abs(a(i+1,i))<1.e-99_r16.and.abs(a(i,i))<1.e-99_r16) then
    a(i,i)=eps
    cycle
  else if (abs(a(i+1,i))<1.e-99_r16) then
    cycle
  end if
  if(abs(a(i,i))<abs(a(i+1,i)))then
    iwork(i) = 1
    do j=i,m
      r = a(i,j)
      a(i,j) = a(i+1,j)
      a(i+1,j) = r
    end do
  end if
  r = -a(i+1,i)/a(i,i)
  a(i+1,i) = r
  do j=l,m
    a(i+1,j) = a(i+1,j)+r*a(i,j)
  end do
end do
      
if(abs(a(m,m))<1.e-99_r16) then
  a(m,m) = eps
end if

! the vector (1,1,...,1) is stored in the place of the right
! hand side column vector.
do i=1,m
  work(i) = 1.0_r16
end do
do i=m+1,kl
  work(i) = 0.0_r16
end do

! the inverse iteration is performed on the matrix until the
! infinite norm of the right-hand side vector is greater
! than the bound defined as  0.01(n*ex).
bound = 0.01_r16/(ex * float(kl))
ns = 0
iter = 1

! the backsubstitution.
do while(.true.)
  r = 0.0_r16
  do i=1,m
    j = m-i+1
    s = work(j)
    if(j/=m)then
      l = j+1
      do k=l,m
        sr = work(k)
        s = s - sr*a(j,k)
      end do
    end if
    work(j) = s/a(j,j)
    t = abs(work(j))
    r = max(r,t)
  end do

! the computation of the right-hand side vector for the new
! iteration step.
  do i=1,m
    work(i) = work(i)/r
  end do

! the computation of the residuals and comparison of the
! residuals of the two successive steps of the inverse
! iteration. if the infinite norm of the residual vector is
! greater than the infinite norm of the previous residual
! vector the computed eigenvector of the previous step is
! taken as the final eigenvector.
  r1 = 0.0_r16
  do i=1,m
    t = sum(a(i,i:m)*work(i:m))
    t = abs(t)
    r1 = max(r1,t)
  end do
  if(iter>1)then
    if(previs<=r1)exit
  end if
  do i=1,m
    vecr(i,ivec) = work(i)
  end do
  previs = r1
  if(ns==1.or.iter>6)exit
  iter = iter+1
  if(r>=bound)then
    ns = 1
  end if
! gaussian elimination of the right-hand side vector.
  k = m-1
  do i=1,k
    if(iwork(i)==0)then
      work(i+1)=work(i+1)+work(i)*a(i+1,i)
    else
      r = work(i+1)
      work(i+1)=work(i)+work(i+1)*a(i+1,i)
      work(i) = r
    end if
  end do
end do
      
if (iter<=6) then
  indic(ivec) =2
end if
      
if(m/=kl)then
  j = m+1
  do i=j,kl
    vecr(i,ivec) = 0.0_r16
  end do
end if
      
return
end subroutine realve
      
subroutine scaler(a,h,prfact,enorm)
use newmpar_m
implicit none
#ifdef pgi
integer, parameter :: r16 = kind(1._8)
#else
integer, parameter :: r16 = kind(1._16)
#endif
! the following real variables were initially single prec.-
! bound1,bound2,enorm
integer i,j
integer iter,ncount
real(kind=r16), dimension(kl,kl) :: a,h
real(kind=r16), dimension(kl) :: prfact
real(kind=r16) enorm
real(kind=r16) fnorm,column,row
real(kind=r16) bound1,bound2,q,factor
! this sub.routine stores the matrix of the order n from the
! array a into the array h. afterward the matrix in the
! array a is scaled so that the quotient of the absolute sum
! of the off-diagonal elements of column i and the absolute
! sum of the off-diagonal elements of row i lies within the
! values of bound1 and bound2.
! the component i of the eigenvector obtained by using the
! scaled matrix must be divided by the value found in the
! prfact(i) of the array prfact. in this way the eigenvector
! of the non-scaled matrix is obtained.
!
! after the matrix is scaled it is normalised so that the
! value of the euclidian norm is equal to one.
! if the process of scaling was not successful the original
! matrix from the array h would be stored back into a and
! the eigenproblem would be solved by using this matrix.
! nm defines the first dimension of the arrays a and h. nm
! must be greater or equal to n.
! the eigenvalues of the normalised matrix must be
! multiplied by the scalar enorm in order that they become
! the eigenvalues of the non-normalised matrix.
do i=1,kl
  do j=1,kl
    h(i,j) = a(i,j)
  end do
  prfact(i)= 1.0_r16
end do
bound1 = 0.75_r16
bound2 = 1.33_r16
iter = 0
ncount = 0
do while (ncount<kl.and.iter<=2*kl)
  ncount = 0
  do i=1,kl
    column = 0.0_r16
    row    = 0.0_r16
    do j=1,kl
      if(i/=j)then
        column=column+abs(a(j,i))
        row   =row   +abs(a(i,j))
       end if
    end do
    if(column<1.e-99_r16.or.row<1.e-99_r16) then
      ncount = ncount + 1
      cycle
    end if
    q = column/row
    if (q>bound1.or.q<=bound2) then
      ncount = ncount + 1
      cycle
    end if
    factor = sqrt(q)
    do j=1,kl
      if(i/=j)then
        a(i,j) = a(i,j)*factor
        a(j,i) = a(j,i)/factor
      end if
    end do
    prfact(i) = prfact(i)*factor
  end do
  iter = iter+1
end do
      
if (iter<=2*kl) then
  fnorm = 0.0_r16
  do i=1,kl
    do j=1,kl
      q = a(i,j)
      fnorm = fnorm+q*q
    end do
  end do
  fnorm = sqrt(fnorm)
  do i=1,kl
    do j=1,kl
      a(i,j)=a(i,j)/fnorm
    end do
  end do
  enorm = fnorm
      
else
  do i=1,kl
! modification suggested by referee in a.c.m.certification
    prfact(i)=1.0_r16
    do j=1,kl
      a(i,j) = h(i,j)
    end do
  end do
  enorm = 1.0_r16

end if

return
end subroutine scaler
      
subroutine sigtosigh(sig,sigmh,kl)
implicit none
#ifdef pgi
integer, parameter :: r16 = kind(1._8)
#else
integer, parameter :: r16 = kind(1._16)
#endif
integer, intent(in) :: kl
! these routines are written from top down
real(kind=r16), dimension(kl), intent(in) :: sig
real(kind=r16), dimension(kl+1), intent(out) :: sigmh
sigmh(1)=1._r16
sigmh(2:kl)=.5_r16*(sig(1:kl-1)+sig(2:kl))
sigmh(kl+1)=0._r16
return
end subroutine sigtosigh
      
subroutine sightosig(sig,sigmh,kl)
implicit none
#ifdef pgi
integer, parameter :: r16 = kind(1._8)
#else
integer, parameter :: r16 = kind(1._16)
#endif
integer, intent(in) :: kl
real(kind=r16), dimension(kl), intent(out) :: sig
real(kind=r16), dimension(kl+1), intent(in) :: sigmh
sig(1:kl) = .5+r16*(sigmh(2:kl+1)+sigmh(1:kl))
return
end subroutine sightosig
