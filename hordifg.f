      subroutine hordifgt    !  globpea version    N.B. k loop in here
!     usual scheme
      parameter (nhorjlm=0)   ! 1 for jlm 3D deformation rather than Smagorinsky
c     called from globpe (now not tendencies),
c     called for -ve nhor
c     for +ve nhor see hordifg 
c     It has -ve nhorps option:
c        nhorps=-1 does only T & qg horiz diff.
c        nhorps=-2 does only u &  v horiz diff.
c        nhorps=-3 does only qg     horiz diff.
c     and u,v have same options as T (e.g.nhor=-157)
c     this one has got map factors
c     N.B. no trace_gases yet
c     has jlm nhorx option as last digit of nhor, e.g. -157
      include 'newmpar.h'
      include 'arrays.h'
      include 'indices.h'
      include 'map.h'
      include 'nlin.h'
      include 'parm.h'
      include 'sigs.h'
      include 'vecsuv.h'   ! vecsuv info
      include 'xarrs.h'
      common/work2/uc(ifull),vc(ifull),wc(ifull),ee(ifull),ff(ifull) 
     .      ,xfact(ifull),yfact(ifull),t_kh(ifull),tx_fact(ifull)
     .      ,ty_fact(ifull),ptemp(ifull),dum(7*ifull)
      data nf/2/
      ind(i,j,n)=i+(j-1)*il+n*il*il  ! *** for n=0,5

c     nhorx used in hordif  ! previous code effectively has nhorx=0
c           = 1 u, v, T, qg  diffusion reduced near mountains
c           = 4              diffusion not reduced near mountains    
c           = 7 u, v, T, qg  diffusion reduced near mountains (bottom 2/3 only)

c     khor (set up in darlam.f): if non-zero increases the diffusion in the
c     upper model levels   (never used these days)
c         khor > 0 progessively doubles for khor upper levels
c         khor < 0 linearly increases above sig = 0.2 (up to factor of 1-khor)

c     this horizontal diffusion routine allows compensation for
c     being on sigma surfaces through reductions of kh over orography
c     this is called for -nhor.ge.50; the value of nhor (in metres) gives the
c     scaling value for deltaz
c     in namelist khdif is fujio's a**2, e.g. 4.
      if(nhorps.gt.0)stop 'nhorps > 0 not permitted in hordifgt'

c     set up topography reduction factors for each type of location
c     expect power nf to be about 1 or 2 (see data statement)
      delphi=1.e6  ! turns off reduction (can also use nhorx=4)
      if(abs(nhor).ge.50)then
         nhora=10*(abs(nhor)/10)    ! e.g. 150  for nhor=-157
         nhorx=abs(nhor)-nhora      ! e.g.   7  for nhor=-157
         delphi=nhora*9.806
      endif

      do iq=1,ifull
       ptemp(iq)=ps(iq)**.286
       tx_fact(iq)=1./(1.+(abs(zs(ie(iq))-zs(iq))/delphi)**nf)
       ty_fact(iq)=1./(1.+(abs(zs(in(iq))-zs(iq))/delphi)**nf)
      enddo   !  iq loop
c     above code independent of k

       if(diag)then
         print *,'hordifgt u ',(u(idjd,k),k=1,kl)
         print *,'hordifgt v ',(v(idjd,k),k=1,kl)
         print *,'ax,ay,az ',ax(idjd),ay(idjd),az(idjd)
         print *,'bx,by,bz ',bx(idjd),by(idjd),bz(idjd)
       endif

      do k=1,kl
       hdif=dt*hdiff(k)/ds  ! N.B.  hdiff(k)=khdif*.1 
!        in hordifgt, need to calculate Cartesian components 
         do iq=1,ifull
          uc(iq)=ax(iq)*u(iq,k) + bx(iq)*v(iq,k)
          vc(iq)=ay(iq)*u(iq,k) + by(iq)*v(iq,k)
          wc(iq)=az(iq)*u(iq,k) + bz(iq)*v(iq,k)
         enddo     ! iq loop

      if(nhorjlm.eq.1)then
c      jlm scheme using 3D uc, vc, wc
       do iq=1,ifull
        cc=(uc(ie(iq))-uc(iw(iq)))**2+(uc(in(iq))-uc(is(iq)))**2
     .    +(vc(ie(iq))-vc(iw(iq)))**2+(vc(in(iq))-vc(is(iq)))**2
     .    +(wc(ie(iq))-wc(iw(iq)))**2+(wc(in(iq))-wc(is(iq)))**2
!       N.B. using double grid length
        t_kh(iq)= .5*sqrt(cc)*hdif/em(iq)  ! this one without em in D terms
       enddo   !  iq loop

      else
c      uses (dv/dx+du/dy)**2 + .5*(du/dx)**2 + .5*(dv/dy)**2
c      following Kikuchi et al. 1981      now Smag. Wed  04-30-1997
!      N.B. original Smag. had m on top (in D formulae) and khdif=3.2
!      More recently (21/9/00) I think original Smag has khdif=0.8
!      Smag's actual diffusion also differentiated Dt and Ds
c      t_kh is kh at t points

c      use ee and ff arrays temporarily for x and y derivs for 2nd deformation
       do iq=1,ifull
        ee(iq)= v(ie(iq),k)-v(iw(iq),k)              ! globpea
        ff(iq)= u(in(iq),k)-u(is(iq),k)              ! globpea
!       some bdy vals changed in next loops
       enddo   !  iq loop

c      some special boundary values switch signs as well as vel components
       do n=0,npanels
c       following treats unusual panel boundaries
        if(npann(n).ge.100)then
          do i=1,il
           iq=ind(i,il,n)
           ff(iq)= -v(in(iq),k)-u(iq-il,k)     ! globpea sign switch too
          enddo  ! i loop
        endif      ! (npann(n).ge.100)
        if(npane(n).ge.100)then
          do j=1,il
           iq=ind(il,j,n)
           ee(iq)= -u(ie(iq),k)-v(iq-1,k)     ! globpea sign switch too
          enddo   ! j loop
        endif      ! (npane(n).ge.100)
        if(npanw(n).ge.100)then
          do j=1,il
           iq=ind(1,j,n)
           ee(iq)= v(iq+1,k)+u(iw(iq),k)      ! globpea sign switch too
          enddo   ! j loop
        endif      ! (npanw(n).ge.100)
        if(npans(n).ge.100)then
          do i=1,il
           iq=ind(i,1,n)
           ff(iq)= u(iq+1,k)+v(is(iq),k)      ! globpea sign switch too
          enddo   ! i loop
        endif      ! (npans(n).ge.100)
       enddo      ! n loop

       do iq=1,ifull
c       better to use ordinary u & v for divergence-type calcs
        aa=.5*(u(ieu(iq),k)-u(iwu(iq),k))    ! globpea code
        bb=.5*(v(inv(iq),k)-v(isv(iq),k))    ! globpea code
        cc=.5*(ee(iq)+ff(iq))   !  .5 because double grid length
!       cc=cc**2+(aa**2+bb**2)*.5   ! this one for Kikuchi
        cc=cc**2+(aa-bb)**2         ! this one for Smagorinsky
        t_kh(iq)= sqrt(cc)*hdif/em(iq)  ! this one without em in D terms
       enddo   !  iq loop
c      ee now finished with (not used till now in globpea)
      endif    !  (nhorjlm.eq.1)

       do iq=1,ifull
        xfact(iq)=(t_kh(ie(iq))+t_kh(iq))*.5
        yfact(iq)=(t_kh(in(iq))+t_kh(iq))*.5
       enddo   !  iq loop
       if((nhorx.ge.7.and.k.le.2*kl/3).or.nhorx.eq.1)then
         do iq=1,ifull
          xfact(iq)=xfact(iq)*tx_fact(iq)
          yfact(iq)=yfact(iq)*ty_fact(iq)
         enddo   !  iq loop
       endif   ! (nhorx.ge.7.and.k.le.2*kl/3).or.nhorx.eq.1

       if(nhorps.eq.0.or.nhorps.eq.-2)then ! for nhorps=-1,-3 don't diffuse u,v
         do iq=1,ifull
          emi=1./em(iq)**2
          ucc=( uc(iq)*emi
     .      +xfact(iq)*uc(ie(iq)) +xfact(iwu2(iq))*uc(iw(iq))
     .      +yfact(iq)*uc(in(iq)) +yfact(isv2(iq))*uc(is(iq)) )
     .      /(emi+xfact(iq)+xfact(iwu2(iq))+yfact(iq)+yfact(isv2(iq)))
          vcc=( vc(iq)*emi
     .      +xfact(iq)*vc(ie(iq)) +xfact(iwu2(iq))*vc(iw(iq))
     .      +yfact(iq)*vc(in(iq)) +yfact(isv2(iq))*vc(is(iq)) )
     .      /(emi+xfact(iq)+xfact(iwu2(iq))+yfact(iq)+yfact(isv2(iq)))
          wcc=( wc(iq)*emi
     .      +xfact(iq)*wc(ie(iq)) +xfact(iwu2(iq))*wc(iw(iq))
     .      +yfact(iq)*wc(in(iq)) +yfact(isv2(iq))*wc(is(iq)) )
     .      /(emi+xfact(iq)+xfact(iwu2(iq))+yfact(iq)+yfact(isv2(iq)))
          u(iq,k)=ax(iq)*ucc +ay(iq)*vcc +az(iq)*wcc
          v(iq,k)=bx(iq)*ucc +by(iq)*vcc +bz(iq)*wcc
         enddo   !  iq loop
       endif   ! nhorps.ge.0

       if(diag)then
         print *,'k,id,jd,idjd ',k,id,jd,idjd
         print *,'k, xfact, xfactw ',k,xfact(idjd),xfact(iwu2(idjd))
         print *,'k, yfact, yfacts ',k,yfact(idjd),yfact(isv2(idjd))
         print *,'k, uc,uce,ucw,ucn,ucs '
     .     ,k,uc(idjd),uc(ie(idjd)),uc(iw(idjd))
     .     ,uc(in(idjd)),uc(is(idjd))
         print *,'k,ee,ff,u,v ',
     .            k,ee(idjd),ff(idjd),u(idjd,k),v(idjd,k)
       endif

       if(nhorps.ne.-2)then   ! for nhorps=-2 don't diffuse T, qg
c        do t diffusion based on potential temperature ff
         do iq=1,ifull
          ee(iq)=qg(iq,k)
          ff(iq)=t(iq,k)/ptemp(iq)  ! watch out for Chen!
         enddo   !  iq loop
         if(nhorps.ne.-3)then   ! for nhorps=-3 don't diffuse T; only qg
           do iq=1,ifull
            emi=1./em(iq)**2
            t(iq,k)= ptemp(iq) * ( ff(iq)*emi
     .        +xfact(iq)*ff(ie(iq)) +xfact(iwu2(iq))*ff(iw(iq))
     .        +yfact(iq)*ff(in(iq)) +yfact(isv2(iq))*ff(is(iq)) )
     .        /(emi+xfact(iq)+xfact(iwu2(iq))+yfact(iq)+yfact(isv2(iq)))
           enddo   !  iq loop
         endif   ! (nhorps.ne.-3)
         do iq=1,ifull
          emi=1./em(iq)**2
          qg(iq,k)=( ee(iq)*emi
     .      +xfact(iq)*ee(ie(iq)) +xfact(iwu2(iq))*ee(iw(iq))
     .      +yfact(iq)*ee(in(iq)) +yfact(isv2(iq))*ee(is(iq)) )
     .      /(emi+xfact(iq)+xfact(iwu2(iq))+yfact(iq)+yfact(isv2(iq)))
         enddo   !  iq loop
       endif    ! (nhorps.ge.-1)

      enddo    !  k loop
      return
      end
