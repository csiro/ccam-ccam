c $Log$
c Revision 1.1  2003/08/13 01:24:20  dix043
c Initial revision
c
c Revision 1.1  1996/10/17  05:14:15  mrd
c Initial revision
c
c Revision 1.1  91/07/19  11:40:56  ldr
c Initial revision
c 
      block data extra_s
c     this contains data statments taken from other routines for variabl
c     in common blocks. such declarations must be in block data routines
      parameter (js=9,jsp1=js+1,jsp2=js+2,js2=2*js,jsm1=js-1)
      parameter (js2m2=js2-2)
      parameter (mlocd=4)
      parameter (jf=12,jg2=10,jcf=jf/2,jg=jg2/2,jf2=2*jf)
      parameter (kmax=9,l=kmax,lp1=l+1)
      parameter (lm1=l-1)
      parameter (lp3=l+3)
      parameter (lp1m=360,lp1m1=lp1m-1)
      parameter (lp2=l+2,lm2=l-2,kx=l,jmax=jg2)
      parameter (imax=jg2,jmaxp2=imax+2)
      parameter (kpad=2002,kdm=l)
      parameter (nb=l,nb1=nb+1,kd2=20,kd1p=11,kd2m=kd2-1)
      parameter (nl=l,nl1=nl+1,nl2=nl+2,nl3=nl+3)
      parameter (kd=10)
c     implicit half precision (a-h,o-z)
      common/cabd/cao3sw(js), cah2sw(js), cbsw(js)
      common/swnuin/press(kd),p(kd),r(kd),ro3(kd),cwca(nl2),cwcb(nl2),
     1     coca(nl2),cloud(nl2),kth(nl2),kbh(nl2),pr2(kdm),cosz,tauda,
     2     nc,nc1,nc2,nc3,nc4,solar,g
c     data statements for cao3sw, cah2sw, cbsw taken from imptfh
      data cao3sw /0.,.21,.48,.69,5*0./
      data cah2sw /0.,.21,.48,.69,5*0./
      data cbsw  / 0., .01, .04, .07, 5*0./
c  data statement for g & solar from swnew in extras.
      data g,solar/980.6,2./
      end
