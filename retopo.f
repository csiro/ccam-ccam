      subroutine retopo(psl,zsold,zs,t,qg)
c     this version (Aug 2003) allows -ve zsold (from spectral model),
c     but assumes new zs is positive for atmospheric purposes
c     this routine redefines psl, ps, t to compensate for zsold going to zs
c     (but does not overwrite zs itself here)
c     called by indata and nestin for newtop.ge.1
!     nowadays just for ps and atmospheric fields Mon  08-23-1999
      include 'newmpar.h'
      include 'const_phys.h'
      include 'parm.h'
      include 'sigs.h'
      real psl(ifull),zsold(ifull),zs(ifull)
      real qg(ifull,kl),t(ifull,kl)
c     first dum inserted so as not to conflict with zsb in nestin
      common/work2/dum(ifull),ps(ifull),psold(ifull),dum2(ifull,15)
      real told(kl),qgold(kl)
      do iq=1,ij
       psold(iq)=1.e5*exp(psl(iq))
       psl(iq)=psl(iq)+(zsold(iq)-max(0.,zs(iq)))/(rdry*t(iq,1))
       ps(iq)=1.e5*exp(psl(iq))
      enddo
c     now alter temperatures to compensate for new topography
      if(ktau.lt.100)then
        print *,'entering retopo, old t ',(t(idjd,k),k=1,kl)
        print *,'entering retopo, old qg ',(qg(idjd,k),k=1,kl)
        print *,'zsold,zs,psold,ps ',
     .           zsold(idjd),zs(idjd),psold(idjd),ps(idjd)
      endif  ! (ktau.lt.100)
      do iq=1,ifull
       do k=1,kl
        qgold(k)=qg(iq,k)
        told(k)=t(iq,k)
       enddo  ! k loop
       do k=1,kl-1
        sig2=sig(k)*ps(iq)/psold(iq)
        if(sig2.ge.sig(1))then
c         assume 6.5 deg/km, with dsig=.1 corresponding to 1 km
          t(iq,k)=told(1)+(sig2-sig(1))*6.5/.1  
        else
          do kkk=2,kl
           if(sig2.gt.sig(kkk))go to 526
          enddo
526       t(iq,k)=(told(kkk)*(sig(kkk-1)-sig2)+told(kkk-1)*
     .             (sig2-sig(kkk)))/(sig(kkk-1)-sig(kkk))
          qg(iq,k)=(qgold(kkk)*(sig(kkk-1)-sig2)
     .              +qgold(kkk-1)*(sig2-sig(kkk)))/(sig(kkk-1)-sig(kkk))
        endif
       enddo  ! k loop
      enddo   ! iq loop
      if(ktau.lt.100)then
        print *,'retopo new t ',(t(idjd,k),k=1,kl)
        print *,'retopo new qg ',(qg(idjd,k),k=1,kl)
      endif  ! (ktau.lt.100)

      return
      end
