      subroutine resetd(x1,x2,x3,x4,n)
      implicit none
      integer n, i
      real x1(n), x2(n), x3(n), x4(n)
      real avg, a1, b1, b2
      do 100 i=1,n
          avg=.25*(x1(i)+x2(i)+x3(i)+x4(i))
          a1=.5*(x2(i)-x4(i))
          b1=.5*(x1(i)-x3(i))
          b2=.25*((x1(i)+x3(i))-(x2(i)+x4(i)))
          x1(i)=avg
          x2(i)=a1
          x3(i)=b1
          x4(i)=b2
 100  continue
      return
      end
