      subroutine printa(name,a,ktau,level,i1,i2,j1,j2,bias,facti)
c     printa has automatic choice of fact if facti=0.
      include 'newmpar.h'
      character*(*) name
      dimension a(il,jl),col(il+jl)
      fact=facti
      if(facti.eq.0.)fact=10./abs(a((i1+i2)/2,(j1+j2)/2))
      print 9 ,name,ktau,level,bias,fact
9     format(/1x,a4,' ktau =',i7,'  level =',i3,'  addon =',g8.2,
     . '  has been mult by',1pe8.1)
      ja=j1
22    jb=min(ja+24,j2)
      print 91,(j,j=ja,jb)
91    format(4x,25i6)
      do i=i1,i2
       do j=ja,jb
        col(j)=(a(i,j)-bias)*fact
       enddo
       print 92,i,(col(j),j=ja,jb)
92     format(i5,25f6.2)
      enddo
      if(jb.eq.j2)return
      ja=jb+1
      go to 22
      end
