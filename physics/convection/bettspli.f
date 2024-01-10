! Conformal Cubic Atmospheric Model
    
! Copyright 2015 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
      
      subroutine bettspli ( nold, xold, yold, y2, nnew,xnew,ynew,p,q) ! was spline
c     this is called by the Betts-Miller parameterization
c
c***********************************************************************
c
c this is a one-dimensional cubic spline fitting routine
c  programmed for a small scalar machine.
c
c  programmed by z.janjic
c
c  nold - number of given values of the function. must be ge 3.
c  xold - locations of the points at which the values of the
c         function are given.
c  yold - the given values of the function at the points xold.
c  y2   - the second derivatives at the points xold.  if natural
c         spline is fitted, y2(1)=0. and y2(nold)=0. must be
c         specified.
c  nnew - number of values of the function to be calculated
c  xnew - locations of the points at which thew values of the
c         function are calculated.  xnew(k) must be ge xold(1)
c         and le xold(nold).
c  ynew - the values of the function calculated
c  p, q - auxilary vectors of the length nold-2.
c
c***********************************************************************
c
      dimension xold(nold), yold(nold), y2(nold), p(nold), q(nold)
      dimension xnew(nnew), ynew(nnew)
      
      ak=0.
      bk=0.
      ck=0.
c
      noldm1 = nold-1
c
      dxl = xold(2) - xold(1)
      dxr = xold(3) - xold(2)
      dydxl = (yold(2)-yold(1))/dxl
      dydxr = (yold(3)-yold(2))/dxr
      rtdxc = .5/(dxl+dxr)
c
      p(1) = rtdxc*(6.*(dydxr-dydxl)-dxl*y2(1))
      q(1) = -rtdxc*dxr
c
      if ( nold.eq.3 ) go to 700
c
      k = 3
c
 100  dxl = dxr
      dydxl = dydxr
      dxr = xold(k+1)-xold(k)
      dydxr = ( yold(k+1) - yold(k) ) / dxr
      dxc = dxl+dxr
      den = 1./(dxl*q(k-2)+dxc+dxc)
c
      p(k-1) = den*(6.*(dydxr-dydxl)-dxl*p(k-2))
      q(k-1) = -den*dxr
c
      k = k+1
      if ( k .lt. nold ) go to 100
c
 700  continue
c
      k = noldm1
c
 200  continue
c
      y2(k) = p(k-1) + q(k-1)*y2(k+1)
c
      k = k-1
c
      if ( k.gt. 1 ) go to 200
c
      k1 = 1
c
 300  continue
c
      xk = xnew(k1)
c
      do 400 k2 = 2 , nold
c
        if ( xold(k2) .le. xk ) go to 400
c
        kold = k2-1
        go to 450
c
 400  continue
c
      ynew(k1) = yold(nold)
      go to 600
c
 450  continue
c
      if ( k1 .eq. 1 ) go to 500
      if ( k .eq. kold ) go to 550
c
 500  continue
c
      k = kold
c
      y2k = y2(k)
      y2kp1 = y2(k+1)
      dx = xold(k+1) - xold(k)
      rdx = 1./dx
c
      ak = .16666667*rdx*(y2kp1-y2k)
      bk = .5*y2k
      ck = rdx*(yold(k+1)-yold(k))-.16666667*dx*(y2kp1+y2k+y2k)
c
 550  continue
c
      x = xk-xold(k)
      xsq = x*x
c
      ynew(k1) = ak*xsq*x+bk*xsq+ck*x+yold(k)
c
 600  continue
c
      k1 = k1+1
c
      if ( k1 .le. nnew ) go to 300
c
      return
c
      end
