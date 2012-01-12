      subroutine trim(a,c,rhs,it)
c     u initially now contains rhs; leaves with answer u (jlm)
c     n.b. we now always assume b = 1-a-c
      include 'newmpar.h'
!     N.B.  e, g, temp are just work arrays (not passed through at all)     
      real a(ifull,kl),c(ifull,kl),rhs(ifull,kl)
      real e(ifull,kl),g(ifull,kl),temp(ifull,kl)

c     this routine solves the system
c       a(k)*u(k-1)+b(k)*u(k)+c(k)*u(k+1)=rhs(k)    for k=2,kl-1
c       with  b(k)*u(k)+c(k)*u(k+1)=rhs(k)          for k=1
c       and   a(k)*u(k-1)+b(k)*u(k)=rhs(k)          for k=kl

c     the Thomas algorithm is used
c     save - only needed if common/work removed

      if(it.eq.0)then
        do iq=1,ifull	 
         b=1.-a(iq,1)-c(iq,1)
c        print *,'iq,a,b,c ',iq,a(iq,1),b,c(iq,1)
         e(iq,1)=c(iq,1)/b
        enddo
        do k=2,kl-1
         do iq=1,ifull
          b=1.-a(iq,k)-c(iq,k)
          temp(iq,k)= 1./(b-a(iq,k)*e(iq,k-1))
          e(iq,k)=c(iq,k)*temp(iq,k)
         enddo
        enddo
      endif

c     use precomputed values of e array when available
      do iq=1,ifull
       b=1.-a(iq,1)-c(iq,1)
       g(iq,1)=rhs(iq,1)/b
      enddo
      do k=2,kl-1
       do iq=1,ifull
        g(iq,k)=(rhs(iq,k)-a(iq,k)*g(iq,k-1))*temp(iq,k)
       enddo
      enddo

c     do back substitution to give answer now
      do iq=1,ifull
       b=1.-a(iq,kl)-c(iq,kl)
       rhs(iq,kl)=(rhs(iq,kl)-a(iq,kl)*g(iq,kl-1))/
     .            (b-a(iq,kl)*e(iq,kl-1))
      enddo
      do k=kl-1,1,-1
       do iq=1,ifull
        rhs(iq,k)=g(iq,k)-e(iq,k)*rhs(iq,k+1)
       enddo
      enddo
      return
      end
