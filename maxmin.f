c has more general format & scaling factor  with subroutine average at bottom
      subroutine maxmin(u,char,ktau,fact,kup)
      include 'newmpar.h'
      character*2 char
      dimension u(ifull,kl)
     . ,umin(kl),umax(kl),iumax(kl),jumax(kl),iumin(kl),jumin(kl)
      do k=1,kup
       umax(k)=u(1,k)
       umin(k)=u(1,k)
       jumax(k)=1
       iumax(k)=1
       jumin(k)=1
       iumin(k)=1
       do j=1,jl
        do i=1,il
	  iq=i+(j-1)*il
         if(umax(k).lt.u(iq,k))then
           umax(k)=u(iq,k)
           jumax(k)=j
           iumax(k)=i
         endif
         if(umin(k).gt.u(iq,k))then
           umin(k)=u(iq,k)
           jumin(k)=j
           iumin(k)=i
         endif
        enddo   ! i loop
       enddo   ! j loop
       umax(k)=fact*umax(k)
       umin(k)=fact*umin(k)
      enddo   ! k loop
      if(kup.eq.1)then
        print 970,ktau,char,umax(1),char,umin(1)
970     format(i7,1x,a2,'max ',f8.3,3x,a2,'min ',f8.3)
        print 9705,ktau,iumax(1),jumax(1),iumin(1),jumin(1)
9705    format(i7,'  posij',i4,i4,10x,i3,i4)
        return
      endif   !  (kup.eq.1)

      if(umax(1).ge.1000.)then   ! for radon
        print 961,ktau,char,(umax(k),k=1,kup)
961     format(i7,1x,a2,'max ',10f7.1/(14x,10f7.1)/(14x,10f7.1))
        print 977,ktau,(iumax(k),jumax(k),k=1,kup)
        print 962,ktau,char,(umin(k),k=1,kup)
962     format(i7,1x,a2,'min ',10f7.1/(14x,10f7.1)/(14x,10f7.1))
        print 977,ktau,(iumin(k),jumin(k),k=1,kup)
      elseif(umax(kup).gt.30.)then  ! format for T, & usually u,v
c        print 971,ktau,char,(umax(k),k=1,kup)
c971     format(i7,1x,a2,'max ',10f7.2/(14x,10f7.2)/(14x,10f7.2))
        print 971,ktau,char,(umax(k),k=1,10),char,(umax(k),k=11,kup)
971     format(i7,1x,a2,'max ',10f7.2/(a10,'maX ',10f7.2)/(14x,10f7.2))
        print 977,ktau,(iumax(k),jumax(k),k=1,kup)
        print 972,ktau,char,(umin(k),k=1,kup)
972     format(i7,1x,a2,'min ',10f7.2/(14x,10f7.2)/(14x,10f7.2))
        print 977,ktau,(iumin(k),jumin(k),k=1,kup)
977     format(i7,'  posij',10(i3,i4)/(14x,10(i3,i4))/(14x,10(i3,i4)))
      else  ! for qg & sd
        print 981,ktau,char,(umax(k),k=1,kup)
981     format(i7,1x,a2,'max ',10f7.3/(14x,10f7.3)/(14x,10f7.3))
        print 977,ktau,(iumax(k),jumax(k),k=1,kup)
        print 982,ktau,char,(umin(k),k=1,kup)
982     format(i7,1x,a2,'min ',10f7.3/(14x,10f7.3)/(14x,10f7.3))
        print 977,ktau,(iumin(k),jumin(k),k=1,kup)
      endif
      return
      end
      subroutine average(speed,spmean,spavge)
      include 'newmpar.h'
      include 'sigs.h'
      include 'xyzinfo.h'  ! wts
      real spmean(kl),speed(ifull,kl)
      spavge=0.
      do k=1,kl
       spmean(k)=0.
       do iq=1,ifull
        spmean(k)=spmean(k)+speed(iq,k)*wts(iq)
       enddo  !  iq loop
       spavge=spavge-dsig(k)*spmean(k)  ! dsig is -ve
      enddo
      return
      end
