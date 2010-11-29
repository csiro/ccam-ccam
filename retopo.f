      subroutine retopo(psl,zsold,zs,t,qg)
!     in Jan 07, renamed recalc of ps from here, to reduce confusion     
c     this version (Aug 2003) allows -ve zsold (from spectral model),
c     but assumes new zs is positive for atmospheric purposes
c     this routine redefines psl, t to compensate for zsold going to zs
c     (but does not overwrite zs, ps themselves here)
c     called by indata and nestin for newtop>=1
!     nowadays just for ps and atmospheric fields Mon  08-23-1999
      use cc_mpi, only : mydiag
      use diag_m
      use sigs_m
      implicit none
      include 'newmpar.h'
      include 'const_phys.h'
      include 'parm.h'
      real psl(ifull),zsold(ifull),zs(ifull)
      real qg(ifull,kl),t(ifull,kl)
      real told(kl),qgold(kl),psnew(ifull),psold(ifull),pslold(ifull)
      integer iq,k,kkk
      real sig2
      do iq=1,ifull
       pslold(iq)=psl(iq)
       psold(iq)=1.e5*exp(psl(iq))
       psl(iq)=psl(iq)+(zsold(iq)-zs(iq))/(rdry*t(iq,1))
       psnew(iq)=1.e5*exp(psl(iq))
      enddo
c     now alter temperatures to compensate for new topography
      if(ktau.lt.100.and.mydiag)then
        if(nproc==1)then
          write (6,"('100*psl(b)old#',9f8.2)") 100.*diagvals(pslold)
          write (6,"('100*psl(b)new#',9f8.2)") 100.*diagvals(psl)
	endif
        print *,'retopo: zsold,zs,psold,psnew ',
     .           zsold(idjd),zs(idjd),psold(idjd),psnew(idjd)
        print *,'retopo: old t ',(t(idjd,k),k=1,kl)
        print *,'retopo: old qg ',(qg(idjd,k),k=1,kl)
      endif  ! (ktau.lt.100)
      do iq=1,ifull
       do k=1,kl
        qgold(k)=qg(iq,k)
        told(k)=t(iq,k)
       enddo  ! k loop
       do k=1,kl-1
        sig2=sig(k)*psnew(iq)/psold(iq)
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
      if(ktau.lt.100.and.mydiag)then
        print *,'retopo: new t ',(t(idjd,k),k=1,kl)
        print *,'retopo: new qg ',(qg(idjd,k),k=1,kl)
      endif  ! (ktau.lt.100)

      return
      end
