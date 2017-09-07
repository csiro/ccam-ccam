! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2016 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
module jimcc_m

! This module is based on Martin Dix's jimcc_m.f90.  However, the interface
! is slighly different, hence jimcc.f90

implicit none

private
public jimcc

integer, dimension(0:47), parameter :: kda =                                &
              (/ 0,1,1,2,2,0, 0,1,1,5,5,0, 3,1,1,2,2,3, 3,1,1,5,5,3,        &
                 0,4,4,2,2,0, 0,4,4,5,5,0, 3,4,4,2,2,3, 3,4,4,5,5,3 /)
integer, dimension(0:47), parameter :: kdb =                                &
              (/ 3,10,4,11,5,9, 6,7,1,11,5,0, 3,1,7,8,2,9,  6,4,10,8,2,0,   &
                 0,10,4,2,8,6, 9,7,1,2,8,3,  0,1,7,5,11,6, 9,4,10,5,11,3 /)
integer, dimension(0:47), parameter :: kdc =                                &
              (/ 5,11,3,9,4,10, 5,11,6,0,1,7, 2,8,3,9,7,1,  2,8,6,0,10,4,   &
                 8,2,0,6,4,10, 8,2,9,3,1,7,  11,5,0,6,7,1, 11,5,9,3,10,4 /)
integer, dimension(0:47), parameter :: kna =                                &
              (/ 12,25,26,9,10,17,   18,31,32,3,4,23,                       &
                  0,37,38,21,22,5,   6,43,44,15,16,11,                      &
                 36,1,2,33,34,41,   42,7,8,27,28,47,                        &
                 24,13,14,45,46,29, 30,19,20,39,40,35 /)
!integer, dimension(0:47), parameter :: knb =                                &
!              (/  5,2,1,4,3,0,       11,8,7,10,9,6,                         &
!                 17,14,13,16,15,12, 23,20,19,22,21,18,                      &
!                 29,26,25,28,27,24, 35,32,31,34,33,30,                      &
!                 41,38,37,40,39,36, 47,44,43,46,45,42 /)
integer, dimension(0:47), parameter :: knc =                                &
              (/  1,0,3,2,5,4,       7,6,9,8,11,10,                         &
                 13,12,15,14,17,16, 19,18,21,20,23,22,                      &
                 25,24,27,26,29,28, 31,30,33,32,35,34,                      &
                 37,36,39,38,41,40, 43,42,45,44,47,46 /)
integer, dimension(0:47), save :: ipofig, kgofig
integer, dimension(8,6), save :: igofkg
integer, save :: ig = 0
real, save :: ss
real, parameter :: third = 1./3.
real, parameter :: fourth = 1./4.
real, dimension(30), parameter :: a =                              &
      (/ 1.47713062600964, -0.38183510510174, -0.05573058001191,   &
        -0.00895883606818, -0.00791315785221, -0.00486625437708,   &
        -0.00329251751279, -0.00235481488325, -0.00175870527475,   &
        -0.00135681133278, -0.00107459847699, -0.00086944475948,   &
        -0.00071607115121, -0.00059867100093, -0.00050699063239,   &
        -0.00043415191279, -0.00037541003286, -0.00032741060100,   &
        -0.00028773091482, -0.00025458777519, -0.00022664642371,   &
        -0.00020289261022, -0.00018254510830, -0.00016499474460,   &
        -0.00014976117167, -0.00013646173947, -0.00012478875822,   &
        -0.00011449267279, -0.00010536946150, -0.00009725109376 /)
!real, dimension(30), parameter :: b =                              &
!      (/ 0.67698819751739, 0.11847293456554, 0.05317178134668,     &
!         0.02965810434052, 0.01912447304028, 0.01342565621117,     &
!         0.00998873323180, 0.00774868996406, 0.00620346979888,     &
!         0.00509010874883, 0.00425981184328, 0.00362308956077,     &
!         0.00312341468940, 0.00272360948942, 0.00239838086555,     &
!         0.00213001905118, 0.00190581316131, 0.00171644156404,     &
!         0.00155493768255, 0.00141600715207, 0.00129556597754,     &
!         0.00119042140226, 0.00109804711790, 0.00101642216628,     &
!         0.00094391366522, 0.00087919021224, 0.00082115710311,     &
!         0.00076890728775, 0.00072168382969, 0.00067885087750 /)
real, dimension(2,2,8), parameter :: flip8 =                                &
           reshape ( (/ 1.0,0.0,0.0,1.0,  -1.0,0.0,0.0,1.0,                 &
                        1.0,0.0,0.0,-1.0,  -1.0,0.0,0.0,-1.0,               &
                        0.0,1.0,1.0,0.0,  0.0,-1.0,1.0,0.0,                 &
                        0.0,1.0,-1.0,0.0,  0.0,-1.0,-1.0,0.0 /),            &
                     (/ 2, 2, 8 /) )
real, dimension(3,0:5), save :: f
real, dimension(3,0:11), save :: e
real, dimension(3,3), save :: txe
real, dimension(3,3,0:47), save :: rotg
complex, save :: ci, cip4, cip3oss

contains

subroutine jimcc(em4,ax4,ay4,az4,xx4,yy4,il)
!     like jim6.f but without stretch option
!     hedra1 data is hardwired
!     xx-->xx4, yy-->yy4, fm-->em4, dxa-->ax4, dxb-->ay4, dxc-->az4
implicit none
integer, parameter :: ipanel = 2
integer, parameter :: ngrmax = 1
integer, intent(in) :: il
integer, save :: num = 0
integer np, ngr, i, j
real, dimension(4*il+1,4*il+1), intent(out) :: em4, ax4, ay4, az4
real, dimension(:,:), allocatable :: xa, xb, xc
real(kind=8), dimension(:,:), pointer :: xx4, yy4 ! avoid intent for pointers
!     real dya(np,np),dyb(np,np),dyc(np,np)
!     ngr = 1  at unstaggered positions
!         = 2  at i+.5,j      positions
!         = 3  at i,j+.5      positions
!         = 4  at i+.5,j+.5   positions

if ( size(xx4,1)/=1+4*il .or. size(xx4,2)/=1+4*il ) then
  write(6,*) "ERROR: xx4 argument is invalid in jimcc"
  stop
end if

if ( size(yy4,1)/=1+4*il .or. size(yy4,2)/=1+4*il ) then
  write(6,*) "ERROR: yy4 argument is invalid in jimcc"
  stop
end if

allocate( xa(4*il+1,4*il+1), xb(4*il+1,4*il+1), xc(4*il+1,4*il+1) )

np = 4*il + 1
CALL INROT

do ngr=1,ngrmax
  call rgrid(xa,xb,xc,ax4,ay4,az4,em4,np,ipanel,ngr)
! these are values on the sphere

  if(num==0)then
   do j=1,np,np-1
    do i=1,np,np-1
     write(6,*)'in jimcc xa,xb,xc: ',i,j,xa(i,j),xb(i,j),xc(i,j)
    enddo
   enddo
   do j=np/4,np+1-np/4,np+1-np/4-np/4
    do i=np/4,np+1-np/4,np+1-np/4-np/4
     write(6,*)'in jimcc xa,xb,xc: ',i,j,xa(i,j),xb(i,j),xc(i,j)
    enddo
   enddo
   write(6,*)'in jimcc xa,xb,xc: ',1,5,xa(1,5),xb(1,5),xc(1,5)
   write(6,*)'now imposing 8-fold symmetry (jlm)'
  endif  ! (num==0)

  do i=1,(np+1)/2
    xc(i,1)=max(-1.,xc(i,1)) ! for rounding errors
  enddo
  do j=1,(np+1)/2
   do i=j+1,(np+1/2)  ! rest of bottom LH corner
    xa(j,i)=xa(i,j)
    xb(j,i)=xc(i,j)
    xc(j,i)=xb(i,j)
   enddo
   do i=1,np/2        ! all of bottom RH corner
    xa(np+1-i,j)=xa(i,j)
    xb(np+1-i,j)=-xb(i,j)
    xc(np+1-i,j)=xc(i,j)
   enddo
   do i=1,np          ! all of top half
    xa(i,np+1-j)=xa(i,j)
    xb(i,np+1-j)=xb(i,j)
    xc(i,np+1-j)=-xc(i,j)
   enddo
  enddo

  if(num==0)then
   num=1
   do j=1,np,np-1
    do i=1,np,np-1
     write(6,*)'in jimcc xa,xb,xc: ',i,j,xa(i,j),xb(i,j),xc(i,j)
    enddo
   enddo
   do j=np/4,np+1-np/4,np+1-np/4-np/4
    do i=np/4,np+1-np/4,np+1-np/4-np/4
     write(6,*)'in jimcc xa,xb,xc: ',i,j,xa(i,j),xb(i,j),xc(i,j)
    enddo
   enddo
   write(6,*)'in jimcc xa,xb,xc: ',1,5,xa(1,5),xb(1,5),xc(1,5)
  endif  ! (num==0)

! now convert these to just x,y values on the cube
  do j=1,np
   do i=1,np
    xx4(i,j)=real(xb(i,j),8)/real(xa(i,j),8)
    yy4(i,j)=real(xc(i,j),8)/real(xa(i,j),8)
   enddo
  enddo

 enddo  ! ngr loop
 
 deallocate( xa, xb, xc )

 return
 end subroutine jimcc

!------------------------------------------------------------------------------!
!   R.J.Purser, National Meteorological Center, Washington D.C.  1994          !
!                   SUBROUTINE  INROT                                          !
!                                                                              !
!   Initialize the rotation matrix needed to set the orientation               !
!   of the cube to user's specification, then call general initialization      !
!   routine INHEDRA to make the transformation code ready for use.             !
!                                                                              !
!  ROT is the user-defined orthogonal matrix determining the orientation of    !
!      the mapping-cube with respect to the standard earth-centered            !
!      coordinate basis defined in the comments below.                         !
!  IG1  is an array of 6 group indices. For designated map-panel IPANEL, the   !
!           element IG1(IPANEL) is associated with the octant of this map      !
!           in contact with the map-origin and that abuts the x-axis           !
!   Note: it is up to the user to ensure that ROT is indeed orthogonal and     !
!   that the 6 elements specified in IG1 do indeed belong to the 6 distinct    !
!   faces of the cube in accordance with the numbering convention adopted      !
!   for the triangulation into "group elements".                               !
!------------------------------------------------------------------------------!
SUBROUTINE INROT
implicit none
real, dimension(3,3) :: rot
integer, parameter, dimension(6) :: ig1=(/ 22,13,40,43,9,45 /)
real r2, r3, r6

R2=SQRT(2.)
R3=SQRT(3.)
R6=R2*R3

!  SPECIFY THE ORIENTATION OF THE TRIANGULATED MAPPING-CUBE IN TERMS OF 3
!  ORTHOGONAL UNIT-VECTORS, REFERRED TO HERE AS p, q, r, DEFINED AS FOLLOWS:
!  VECTOR q POINTS TOWARDS THE MID-POINT OF EDGE-3 WHERE TRIANGULAR ELEMENTS
!  24, 29, 41 AND 36 MEET; VECTOR r POINTS TOWARDS VERTEX-0 WHERE ELEMENTS
!  0, 1, 2, 3, 4 AND 5 MEET; VECTOR p POINTS TOWARDS A POSITION ON THE
!  COMMON BOUNDARY OF ELEMENTS 14 AND 15 (IN FACE-3) SUCH THAT IT IS
!  ORTHOGONAL TO BOTH q AND r, OR EQUIVALENTLY, p IS DEFINED AS THE
!  (RIGHT-HANDED) CROSS-PRODUCT, q*r. THE BASIS-VECTORS USED TO EXPRESS
!  p, q, AND r ARE (1): THE UNIT VECTOR POINTING TO LAT,LONG=0,0; (2): THE
!  VECTOR POINTING TO LAT,LONG=0,90E; (3) THE UNIT VECTOR POINTING TO THE
!  NORTH POLE.

!  DEFINE VECTOR p AND MAKE IT THE FIRST COLUMN OF MATRIX ROT:
ROT(1,1)=R6/6.
ROT(2,1)=-R6/6.
ROT(3,1)=-R6/3.

!  DEFINE VECTOR q AND MAKE IT THE SECOND COLUMN OF MATRIX ROT:
ROT(1,2)=R2/2.
ROT(2,2)=R2/2.
ROT(3,2)=0.

!  DEFINE VECTOR r AND MAKE IT THE THIRD COLUMN OF MATRIX ROT:
ROT(1,3)=R3/3.
ROT(2,3)=-R3/3.
ROT(3,3)=R3/3.

!  CUSTOMIZATION OF THE MAPPING TRANSFORMATION IS COMPLETED BY SPECIFYING,
!  FOR EACH NUMBERED MAP-PANEL (FROM 1 TO 6) THE TRIANGULAR ELEMENT THAT
!  CORRESPONDS TO THE 3 RESTRICTIONS IN THE LOCAL COORDINATES x,y:
!          a: x < .5;
!          b: 0. < y;
!          c: y < x.
!  FOR EACH MAP-PANEL, IP, THIS BASIC ELEMENT IS PRESCRIBED IN IG1(IP).
!  THESE, TOGETHER WITH THE ORTHOGONAL MATRIX ROT MADE UP OF p,q,r, ARE
!  PASSED TO THE GENERAL INITIALIZATION ROUTINE INHEDRA, AFTER WHICH, THE
!  MAP-TRANSFORMATION ROUTINES ARE READY FOR USE.
CALL INHEDRA(ROT,IG1)
RETURN
end subroutine inrot

!------------------------------------------------------------------------------!
!   r.j.purser, national meteorological center, washington d.c.  1994          !
!                   subroutine  rgrid                                          !
!                                                                              !
!   set up array of standard earth-centered coordinates of a chosen panel      !
!   of the rancic map. convention used for map coordinates here is that        !
!   each origin is the pole (north or south as appropriate) and the (x,y)      !
!   map coordinates in each panel form a right-handed pair, with x and y both  !
!   between 0 and 1. to choose another panel-corner for origin, or to alter    !
!   chiral convention, just rearrange the table igofkg. for more radical       !
!   change in map-coordinate convention, rewrite this routine!                 !
!                                                                              !
!  <-- xe,ye,ze     earth-centered coordinate of regular map-grid              !
!  --> np           number of points along each grid line (edges included)     !
!  --> ipanel       map-panel index [0,5]                                      !
!------------------------------------------------------------------------------!

subroutine rgrid(xe,ye,ze,dxa,dxb,dxc,em4,np,ipanel,ngr)
!     ngr = 1  at unstaggered positions
!         = 2  at i+.5,j+.5   positions
!         = 3  at i+.5,j      positions
!         = 4  at i,j+.5      positions
!     for staggered grids, actually only use 1,n part of the array
! set up earth-centered coordinates for points of chosen panel of the
!  rancic map
implicit none
integer, intent(in) :: np, ipanel, ngr
integer :: j, jp, i
real, parameter :: stretch=1.
real, parameter :: stretchm=1.-stretch
real, dimension(np,np), intent(out) :: xe, ye, ze
real, dimension(np,np), intent(out) :: dxa,dxb,dxc,em4
real, dimension(np,3,3) :: xc
!real, dimension(np,2) :: dfdx
real, dimension(np) :: xvec, den
real :: d, xadd, yadd, x, y

d=1./real(np-1)
xadd=0.
yadd=0.
if(ngr==2.or.ngr==3)then
  xadd=.5*d
end if
if(ngr==2.or.ngr==4)then
  yadd=.5*d
end if

do j = 0,np-1
 jp = j + 1
 y = real(j)*d  + yadd   ! jlm allows staggered v
 y = .5+(y-.5)*(stretch+stretchm*(2.*y-1.)**2)   !jlm
 do i = 0,np-1
  x = real(i)*d + xadd   ! jlm allows staggered u
  xvec(i+1) = .5+(x-.5)*(stretch+stretchm*(2.*x-1.)**2)   !jlm
 end do
 call vmtoc(xvec,y,ipanel,xe(:,jp),ye(:,jp),ze(:,jp),np)
 call vmtocd(xvec,y,ipanel,xc,em4(:,jp),np)
! return dxa etc as unit vectors
 den(1:np)=sqrt(xc(1:np,1,1)**2+xc(1:np,2,1)**2+xc(1:np,3,1)**2)
 where (den(1:np)<1.e-6)
   den(1:np) = 1.
 end where
 dxa(1:np,jp) = xc(1:np,1,1)/den(1:np)   ! the three components of a vector along dx
 dxb(1:np,jp) = xc(1:np,2,1)/den(1:np)
 dxc(1:np,jp) = xc(1:np,3,1)/den(1:np)
enddo
return
end subroutine rgrid

!------------------------------------------------------------------------------!
!   R.J.Purser, National Meteorological Center, Washington D.C.  1994          !
!                   SUBROUTINE INHEDRA                                         !
!   Initialize variables needed to perform the Rancic-transformation and       !
!   its inverse                                                                !
!------------------------------------------------------------------------------!
SUBROUTINE INHEDRA(ROT,IG1)
implicit none
integer, dimension(6), intent(in) :: ig1
integer, dimension(8) :: igk
integer :: ip, kg, i, j, k, lg
real, dimension(3,3), intent(in) :: rot
real :: r2, r3, r6, r2o2, r3o2, r3o3, r6o3, r3o6, r6o6

!  SET UP CROSS-REFERENCE TABLES CONNECTING GROUP-ELEMENTS IG TO
!  USER-DEFINED NUMBERING AND ORIENTATIONS OF PANELS:
DO IP=1,6
 IGK(1)=IG1(IP)
 IGK(2)=KNA(IGK(1))
 IGK(5)=KNC(IGK(1))
 IGK(6)=KNC(IGK(2))
 IGK(7)=KNA(IGK(5))
 IGK(8)=KNA(IGK(6))
 IGK(3)=KNC(IGK(7))
 IGK(4)=KNC(IGK(8))
 DO KG=1,8
  IG=IGK(KG)
  IGOFKG(KG,IP)=IG
  IPOFIG(IG)=IP
  KGOFIG(IG)=KG
 ENDDO
ENDDO
R2=SQRT(2.)
R3=SQRT(3.)
R6=SQRT(6.)
R2O2=R2/2.
R3O2=R3/2.
R3O3=R3/3.
R6O3=R6/3.
R6O6=R6/6.
R3O6=R3/6.
SS=R2
CI=CMPLX(0.,1.)
CIP4=CI**FOURTH
CIP3OSS=CI**THIRD/SS
F(1,0)=-R6O3
F(2,0)=0.
F(3,0)=R3O3
F(1,1)=R6O6
F(2,1)=-R2O2
F(3,1)=R3O3
F(1,2)=R6O6
F(2,2)=R2O2
F(3,2)=R3O3
E(1,0)=R3O3
E(2,0)=0.
E(3,0)=R6O3
E(1,1)=-R3O6
E(2,1)=.5
E(3,1)=R6O3
E(1,2)=-R3O6
E(2,2)=-.5
E(3,2)=R6O3
E(1,3)=0.
E(2,3)=1.
E(3,3)=0.
E(1,4)=-R3O2
E(2,4)=-.5
E(3,4)=0.
E(1,5)=R3O2
E(2,5)=-.5
E(3,5)=0.
DO J=0,2
 K=J+3
 DO I=1,3
  F(I,K)=-F(I,J)
 ENDDO
ENDDO
DO J=0,5
 K=J+6
 DO I=1,3
  E(I,K)=-E(I,J)
 ENDDO
ENDDO
DO I=1,3
 TXE(1,I)=F(I,0)
 TXE(2,I)=E(I,3)
 TXE(3,I)=E(I,5)
ENDDO
CALL INVMM(TXE)
!  ROTATE THE 6 FACE-VECTORS, F, TO USER-DEFINED ORIENTATION:
f = matmul( rot, f )

!  ROTATE THE 12 EDGE-VECTORS, E, TO USER-DEFINED ORIENTATION:
e = matmul( rot, e )

!  BASED ON THE PRESCRIBED ORIENTATION (DEFINED BY "ROT"),
!  CONSTRUCT THE ROTATION MATRIX ASSOCIATED WITH EACH GROUP ELEMENT LG:

DO LG=0,47
 DO J=1,3
  DO I=1,3
   ROTG(I,J,LG)=F(I,KDA(LG))*TXE(J,1) &
               +E(I,KDB(LG))*TXE(J,2) &
               +E(I,KDC(LG))*TXE(J,3)
  ENDDO
 ENDDO
ENDDO

RETURN
end subroutine inhedra

!----------------------------------------------------------------------------
!   R.J.Purser, National Meteorological Center, Washington D.C.  1994
!                   SUBROUTINE MTOC
!   Transform from map-coordinates to standard earth-centered cartesians
!
!  -->    XX,YY  map-coordinates
!  -->    IPANEL map-panel index
!  <--    XC     standard earth-centered cartesian coordinates
!----------------------------------------------------------------------------
subroutine vmtoc(xx, yy, ipanel, xc, yc, zc, n)
implicit none
integer, intent(in) :: ipanel
integer, intent(in) :: n
integer, dimension(n) :: kg
integer :: i
real, intent(in), dimension(n) :: xx
real, intent(in) :: yy
real, dimension(n), intent(out) :: xc, yc, zc
real, dimension(3) :: xv, tmpc
real, dimension(n) :: x, y, t, xw, yw, h
complex, dimension(n) :: w, z, arg

x = xx
y = yy
kg = 1
where ( x > 0.5 ) 
 kg = kg + 1
 x = 1.0 - x
endwhere
where ( y > 0.5 ) 
 kg = kg + 2
 y = 1.0 - y
endwhere
where (y > x) 
 kg = kg + 4
 t = x
 x = y
 y = t
endwhere

! z=cmplx(x,y)**4
z = cmplx(x,y)*cmplx(x,y)
z = z*z
call vtay (z, a, 30, w)
arg = -ci*w                                ! mrd
where ( abs(arg)<1.e-20 )                  ! mrd
 w = (0.0,0.0)                             ! mrd
elsewhere                                  ! mrd
 w = cip3oss*(-ci*w)**third                ! mrd
endwhere                                   ! mrd
xw = real(w)
yw = aimag(w)
h = 2.0/(1.0 + xw*xw + yw*yw)
do i = 1,n
 xv(1) = xw(i)*h(i)
 xv(2) = yw(i)*h(i)
 xv(3) = h(i) - 1.0
 ig = igofkg(kg(i),ipanel)
 tmpc(:) = matmul ( rotg(:,:,ig), xv )
 xc(i) = tmpc(1)
 yc(i) = tmpc(2)
 zc(i) = tmpc(3)
end do

return
end subroutine vmtoc

!---------------------------------------------------------------------------
!   R.J.Purser, National Meteorological Center, Washington D.C.  1994
!                   SUBROUTINE MTOCD
!   Transform from map-coordinates to standard earth-centered cartesians
!   and simultaneously accumulate the derivative of the transformation in
!   order to provide information on map-scaling factor and relative
!   orientation
!
!  --> XX,YY  map-coordinates (from corner IG, each panel a unit square)
!  --> IPANEL map-panel index
!  <-- XC   augmented jacobian matrix: first two columns represent the
!           derivative with respect to X and Y of the earth-centered
!           cartesian coordinates of the point corresponding to map-image
!           X,Y. These cartesian coordinates themselves are inserted into
!           column-3 of XDC.
!  <-- em4   map-factor at this point
!  <-- DFDX x- and y-derivatives of map-factor here
!---------------------------------------------------------------------------
subroutine vmtocd(xx, yy, ipanel, xc, em4, n)
implicit none
integer, intent(in) :: ipanel
integer, intent(in) :: n
integer, dimension(n) :: kg
integer :: i
real, intent(in), dimension(n) :: xx
real, intent(in) :: yy
real, intent(out), dimension(n) :: em4
real, intent(out), dimension(n,3,3)  :: xc
!real, intent(out), dimension(2,n)    :: dfdx
real, dimension(n,3,3) :: xdc
real, dimension(n,3,2) :: xd
real, dimension(n,2) :: v1
real, dimension(n) :: x, y, t, xw, yw, xwxw, xwyw, ywyw, h, hh, rd, qd, &
                      rdd, qdd, s, dsdx, dsdy, dhdx, dhdy
complex, dimension(n) :: w, z, zu, wu, cd, cdd, arg

x = xx
y = yy
kg = 1
where ( x > 0.5 ) 
 kg = kg + 1
 x = 1.0 - x
endwhere
where ( y > 0.5 ) 
 kg = kg + 2
 y = 1.0 - y
endwhere
where ( y > x ) 
 kg = kg + 4
 t = x
 x = y
 y = t
endwhere
zu = cmplx(x,y)
z = zu**4
call vtaydd (z, a, 30, w, cd, cdd)
arg = -ci*w                                ! mrd
where ( abs(arg)<1.e-20 )                  ! mrd
 wu = (0.0,0.0)                            ! mrd
elsewhere                                  ! mrd
 wu = cip3oss*(-ci*w)**third               ! mrd
endwhere                                   ! mrd
xw = real(wu)
yw = aimag(wu)
xwxw = xw*xw
xwyw = xw*yw
ywyw = yw*yw
h = 2.0/(1.0 + xwxw + ywyw)
hh = h*h
xdc(:,1,3) = xw*h
xdc(:,2,3) = yw*h
xdc(:,3,3) = h - 1.0
xdc(:,1,1) = h - hh*xwxw
xdc(:,2,1) = -hh*xwyw
xdc(:,3,1) = -hh*xw
xdc(:,1,2) = xdc(:,2,1)
xdc(:,2,2) = h - hh*ywyw
xdc(:,3,2) = -hh*yw
where ( abs(z)<1.e-20 ) 
 cd    = 0.0
 cdd   = 0.0
 em4   = 0.0
 v1(:,1) = 0.0
 v1(:,2) = 0.0
elsewhere
 cd = 4.0*wu*cd*z/(3.0*w*zu)
 cdd = 3.0*cd/zu - 2.0*cd*cd/wu + 16.0*wu*z**2*cdd/(3.0*w*zu**2)
 rd = real(cd)
 qd = aimag(cd)
 rdd = real(cdd)
 qdd = aimag(cdd)
 s = sqrt(rd*rd + qd*qd)
 dsdx = (rdd*rd + qdd*qd)/s
 dsdy = (rdd*qd - qdd*rd)/s
 dhdx = -hh*(xw*rd + yw*qd)
 dhdy = -hh*((-xw*qd) + yw*rd)
 em4 = h*s
 v1(:,1) = dhdx*s + h*dsdx
 v1(:,2) = dhdy*s + h*dsdy
endwhere
rd = real(cd)
qd = aimag(cd)
do i=1,3
 xd(:,i,1) = xdc(:,i,1)*rd + xdc(:,i,2)*qd
 xd(:,i,2) = (-xdc(:,i,1)*qd) + xdc(:,i,2)*rd
end do
do i=1,n
 ig = igofkg(kg(i),ipanel)
 !dfdx(i,:) = matmul ( transpose(flip8(:,:,kg(i))), v1(i,:) )
 xdc(i,:,1:2) = matmul ( xd(i,:,:), flip8(:,:,kg(i)) )
 xc(i,:,:) = matmul ( rotg(:,:,ig), xdc(i,:,:) )
end do

return
end subroutine vmtocd


!------------------------------------------------------------------------------!
!   R.J.Purser, National Meteorological Center, Washington D.C.  1994          !
!                   SUBROUTINE  TAY                                            !
!  Evaluate the complex function W of Z whose real                             !
!  Taylor series coefficients are RA.                                          !
!                                                                              !
!  --> Z    function argument (complex)                                        !
!  --> RA   Taylor coefficients (real)                                         !
!  --> N    number of coefficients (starting with the linear term)             !
!  <-- W    Taylor-series approximation of the function (complex)              !
!------------------------------------------------------------------------------!
SUBROUTINE VTAY(Z,RA,N,W)
implicit none
integer, intent(in) :: n
integer i
real, dimension(n), intent(in) :: ra
complex, dimension(:), intent(in) :: z
complex, dimension(:), intent(out) :: w

W=cmplx(0.,0.)
DO I=N,1,-1
 W=(W+RA(I))*Z
ENDDO
RETURN
end subroutine vtay

!------------------------------------------------------------------------------!
!   R.J.Purser, National Meteorological Center, Washington D.C.  1994          !
!                   SUBROUTINE  TAYDD                                          !
!  Evaluate the complex function W of Z whose real                             !
!  Taylor series coefficients are RA, together with its derivative WD.         !
!  and second derivative WDD                                                   !
!                                                                              !
!  --> Z    function argument (complex)                                        !
!  --> RA   Taylor coefficients (real)                                         !
!  --> N    number of coefficients (starting with the linear term)             !
!  <-- W    Taylor-series approximation of the function (complex)              !
!  <-- WD   Taylor series approximation of the derivative of the function W    !
!  <-- WDD  Taylor series approximation of the derivative of the function WD   !
!------------------------------------------------------------------------------!
SUBROUTINE VTAYDD(Z,RA,N,W,WD,WDD)
implicit none
integer, intent(in) :: n
integer :: i
real, dimension(n), intent(in) :: ra
complex, dimension(:), intent(in) :: z
complex, dimension(:), intent(out) :: w, wd, wdd

W=cmplx(0.,0.)
WD=cmplx(0.,0.)
WDD=cmplx(0.,0.)
DO I=N,1,-1
 W=(W+RA(I))*Z
 WD=Z*WD+I*RA(I)
ENDDO
DO I=N,2,-1
 WDD=Z*WDD+I*(I-1)*RA(I)
ENDDO
RETURN
end subroutine vtaydd

!------------------------------------------------------------------------------!
!   R.J.Purser, National Meteorological Center, Washington D.C.  1993          !
!                   SUBROUTINE LUFM                                            !
!  perform l-u decomposition of square matrix a in place with                  !
!  partial pivoting                                                            !
!  For DOUBLE PRECISION version see DLUFM                                      !
!                                                                              !
!  --> a    square matrix to be factorized                                     !
!  <-- ipiv array encoding the pivoting sequence                               !
!  <-- d    indicator for possible sign change of determinant                  !
!  --> m    degree of (active part of) a                                       !
!  --> na   first fortran dimension of a                                       !
!                                                                              !
!------------------------------------------------------------------------------!
subroutine LUFM(a,ipiv,d)
implicit none
integer, dimension(:), intent(out) :: ipiv
integer :: m, j, jp, iquad, i, k, jm
real, dimension(:,:), intent(inout) :: a
real, intent(out) :: d
real :: abig, aa, t, ajj, ajji, aij

m=size(a,1)
d=1.
ipiv(m)=m
do j=1,m-1
 jp=j+1
 abig=abs(a(j,j))
 iquad=j
 do i=jp,m
  aa=abs(a(i,j))
  if (aa <= abig) cycle
  iquad=i
  abig=aa
 end do
 !  swap rows, recording changed sign of determinant
 ipiv(j)=iquad
 if(iquad/=j)then
  d=-d
  do k=1,m
   t=a(j,k)
   a(j,k)=a(iquad,k)
   a(iquad,k)=t
  enddo
 endif
 ajj=a(j,j)
 if(abs(ajj)<1.e-20)then
  jm=j-1
  write(6,*) 'failure in lufact: matrix singular, rank=',jm
  stop
 endif
 ajji=1./ajj
 do i=jp,m
  aij=ajji*a(i,j)
  a(i,j)=aij
  a(i,jp:m)=a(i,jp:m)-aij*a(j,jp:m)
 enddo
enddo
return
end subroutine lufm

!------------------------------------------------------------------------------!
!   R.J.Purser, National Meteorological Center, Washington D.C.  1993          !
!                   SUBROUTINE INVMM                                           !
!  invert matrix, possibly in place (a=b), using the l-u decomposition method  !
!  For DOUBLE PRECISION version see DINVMM                                     !
!                                                                              !
!  --> b    square matrix to be inverted                                       !
!  <-- a    inverse of b                                                       !
!  --> m    degree of (active part of) b and a                                 !
!  --> nb   first fortran dimension of b                                       !
!  --> na   first fortran dimension of a                                       !
!                                                                              !
!   LIMITATION:                                                                !
!    ipiv is an index array, internal to this array, encoding the              !
!    pivoting sequence used. It is given a fortran dimension of NN=500         !
!    in the parameter statement below. If the order of the linear system       !
!    exceeds 500, increase this parameter accordingly                          !
!                                                                              !
!------------------------------------------------------------------------------!
subroutine INVMM(a)
implicit none
real, dimension(:,:), intent(inout) :: a
 integer, dimension(size(a,1)) :: ipiv
 integer :: m, j, i, l
 real :: d, s, t
!-----------------------------------------------

!  Check it's a square matrix
if ( size(a, 1) /= size(a, 2) ) then
 print*, " Can't calculate inverse of non-square matrix "
 print*, " Shape is ", shape(a)
 stop
end if

m = size(a, 1)
call lufm (a, ipiv, d)
!  invert u in place:
do i = 1, m
   a(i,i) = 1.0/a(i,i)
end do
do i = 1, m - 1
   do j = i + 1, m
      s = 0.0
      s = -sum(a(i,i:j-1)*a(i:j-1,j))
      a(i,j) = a(j,j)*s
   end do
end do
!  invert l in place assuming implicitly diagonal elements of unity
do j = 1, m - 1
   do i = j + 1, m
      s = -a(i,j)
      s = s - sum(a(i,j+1:i-1)*a(j+1:i-1,j))
      a(i,j) = s
   end do
end do
!  form the product of u**-1 and l**-1 in place
do j = 1, m - 1
   do i = 1, j
      s = a(i,j)
      s = s + sum(a(i,j+1:m)*a(j+1:m,j))
      a(i,j) = s
   end do
   do i = j + 1, m
      s = 0.0
      s = sum(a(i,i:m)*a(i:m,j))
      a(i,j) = s
   end do
end do
!  permute columns according to ipiv
do j = m - 1, 1, -1
   l = ipiv(j)
   do i = 1, m
      t = a(i,j)
      a(i,j) = a(i,l)
      a(i,l) = t
   end do
end do
return
end subroutine invmm

end module jimcc_m
