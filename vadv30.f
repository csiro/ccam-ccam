      subroutine vadv30(tin,tout,uin,uout,vin,vout)   
!     only calls vadvbess (nvad=7) from Aug 2003      
      parameter (ntest=0) !  0: usual   1: for diagnostic prints
      parameter (nvdep=7) !  0  for original;
!                            1  for newer average vels
!                            2  for newest jlm
!                            3  for newest jlm (quadratic right at end)
!                            4+ for new similar to that of depts3d
!                              - only used in vadv30, vadvl_w (not vadvtvd)
      parameter (nvint=2) !  0 simple+1/2; 1 simple;  **set in vadv30in also
!                            2 for Bessel
!                            3 for spline
!                            7 for Akima
c     does t, u,v then qg, [q1, q2,]

      include 'newmpar.h'
      parameter (kl1=kl+1,kl2=kl+2)
      include 'arrays.h'
      include 'indices.h'
      include 'liqwpar.h'  ! ifullw
      include 'map.h'
      include 'parm.h'
      include 'parmdyn.h'
      include 'parmvert.h'
      include 'sigs.h'
      include 'tracers.h'
      include 'vvel.h'
      include 'xarrs.h'
      common/work3b/wrk1(ifull,kl),ders(ifull,kl)   ! just work arrays here
      common/work3c/kdel(ifull,kl)
      common/work3d/sd(ifull,kl)
      common/work3f/st(ifull,kl),anew(ifull,kl),gwrk(ifull,kl)
      real tin(ifull,kl),tout(ifull,kl),uin(ifull,kl),uout(ifull,kl),
     .     vin(ifull,kl),vout(ifull,kl) ! N.B. tin & tout etc. may be equivalent
      real bb(kl),sddk(0:kl+1)
      real dersh(ifull,kl),sdotder(ifull,kl+1)
      equivalence (dersh,gwrk),(sdotder,wrk1)
      data sddk/kl2*0./
      tfact=1./nvadh   ! simpler alternative
      if(ktau.eq.1)then
        print *,'in vadv30 nvdep,nvint: ',nvdep,nvint
      endif

c     note sdot coming through is at level k-.5
c     & was converted in updps/adjust5 to units of grid-steps/timestep  (+ve upwards)

c     Determine departure points st for all i,j,k:
      if(nvdep.eq.1)then
!       this method not so good if trajectory > 1/2 grid length
        do iq=1,ifull
         do k=1,kl   ! tfact gives ~1/2 timestep's worth
          delsd=tfact*(sdot(iq,k+1)-sdot(iq,k))
          st(iq,k)=k- tfact*(sdot(iq,k+1)+sdot(iq,k))/
     .                  (2.+delsd+delsd*delsd/6.)   ! cont. fraction for exp
         enddo   ! k loop
        enddo  ! iq loop
      endif     ! (nvdep.eq.1)
      if(nvdep.eq.0)then  ! original jlm method follows 
        do iq=1,ifull
         do k=1,kl
c         interpolate sdot to full-level sd with cubic polynomials
          sd(iq,k)=(sdot(iq,k)+sdot(iq,k+1))*.5  ! linear at ends
          bb(k)=sdot(iq,k+1)-sdot(iq,k)
         enddo   ! k loop
         do k=2,kl-1
          sd(iq,k)=sd(iq,k)-(bb(k+1)-bb(k-1))/16.   ! interpolated sdot
         enddo   ! k loop
c        build in 1/2 & 1/6 into sdd & sddd
c        OK for sd, sdd to go off ends of array as sdot(1)=0., sdot(kl+1)=0.
         do k=1,kl    ! sdd is kdot*(dkdot/dk)/2  or d(kdot**2)/dk /4
          sddk(k)=(sdot(iq,k+1)**2-sdot(iq,k)**2)/4.
         enddo
         do k=1,kl    ! sddd is kdot*(dsdd/dk)/3
          sddd=(sdot(iq,k)*(sddk(k)-sddk(k-1))
     .         +sdot(iq,k+1)*(sddk(k+1)-sddk(k)))/6.  
          st(iq,k)=k -tfact*(sd(iq,k) -tfact*(sddk(k) -tfact*sddd))
         enddo   ! k loop
        enddo  ! iq loop
      endif   ! (nvdep.eq.0)
       
      if(nvdep.eq.2)then  !  jlm method in 2003
        do iq=1,ifull
         sd(iq,1) =.75*sdot(iq,2) -.125*sdot(iq,3)     ! quadratic at ends
         sd(iq,kl)=.75*sdot(iq,kl)-.125*sdot(iq,kl-1)  ! quadratic at ends
        enddo  ! iq loop
c       interpolate sdot to full-level sd with cubic polynomials
        do k=2,kl-1
         do iq=1,ifull
          sd(iq,k)=(9.*(sdot(iq,k)+sdot(iq,k+1))
     .                 -sdot(iq,k-1)-sdot(iq,k+2))/16.   
         enddo  ! iq loop
        enddo   ! k loop
c       build in 1/2 & 1/6 into sdd & sddd
        do k=1,kl    ! sdd is kdot*(dkdot/dk)/2  or d(kdot**2)/dk /4
         do iq=1,ifull
          sdd=(sdot(iq,k+1)**2-sdot(iq,k)**2)/4.
!	   sddd is kdot*(dsdd/dk)/3 or kdot*(d2(kdot**2)/dk2)/12
          sddd=sd(iq,k)*(sdot(iq,k)**2-2.*sd(iq,k)**2+sdot(iq,k+1)**2)/3.
          st(iq,k)=k -tfact*(sd(iq,k) -tfact*(sdd -tfact*sddd))
         enddo  ! iq loop
        enddo   ! k loop
      endif   ! (nvdep.eq.2) 
       
      if(nvdep.eq.3)then  !  jlm method in 2003 - quadratic right at end
        do iq=1,ifull
         sd(iq,1) =.25*sdot(iq,2)     ! quadratic right at end
         sd(iq,kl)=.25*sdot(iq,kl)    ! quadratic right at end
        enddo  ! iq loop
c       interpolate sdot to full-level sd with cubic polynomials
        do k=2,kl-1
         do iq=1,ifull
          sd(iq,k)=(9.*(sdot(iq,k)+sdot(iq,k+1))
     .                 -sdot(iq,k-1)-sdot(iq,k+2))/16.   
         enddo  ! iq loop
        enddo   ! k loop
c       build in 1/2 & 1/6 into sdd & sddd
        do k=1,kl    ! sdd is kdot*(dkdot/dk)/2  or d(kdot**2)/dk /4
         do iq=1,ifull
          sdd=(sdot(iq,k+1)**2-sdot(iq,k)**2)/4.
!	   sddd is kdot*(dsdd/dk)/3 or kdot*(d2(kdot**2)/dk2)/12
          sddd=sd(iq,k)*(sdot(iq,k)**2-2.*sd(iq,k)**2+sdot(iq,k+1)**2)/3.
          st(iq,k)=k -tfact*(sd(iq,k) -tfact*(sdd -tfact*sddd))
         enddo  ! iq loop
        enddo   ! k loop
      endif   ! (nvdep.eq.3) 
       
      if(nvdep.eq.4)then  !  newer method similar to depts3d - very poor
        do iq=1,ifull
         sd(iq,1) =.25*sdot(iq,2)     ! quadratic right at end
         sd(iq,kl)=.25*sdot(iq,kl)    ! quadratic right at end
        enddo  ! iq loop
c       interpolate sdot to full-level sd with cubic polynomials
        do k=2,kl-1
         do iq=1,ifull
          sd(iq,k)=(9.*(sdot(iq,k)+sdot(iq,k+1))
     .                 -sdot(iq,k-1)-sdot(iq,k+2))/16.   
         enddo    ! iq loop
        enddo     ! k loop
	 ders=-tfact*sd  ! 3D
        do k=1,kl    
         do iq=1,ifull
          st(iq,k)=k+ders(iq,k) 
         enddo  ! iq loop
        enddo   ! k loop
	 do itn=2,3
c	 do itn=2,2  
         do k=2,kl-1
          do iq=1,ifull
           wrk1(iq,k)=-.5*tfact*sd(iq,k)*(ders(iq,k+1)-ders(iq,k-1))/itn
          enddo  ! iq loop
         enddo   ! k loop
         do iq=1,ifull
!         assume ders=0 off ends and varies quadratically to 1st level	
          wrk1(iq,1)=-tfact*sd(iq,1)*4.*ders(iq,1)/itn
          wrk1(iq,kl)=tfact*sd(iq,kl)*4.*ders(iq,kl)/itn
         enddo  ! iq loop
         ders=wrk1   ! 3D
         st=st+ders  ! 3D
	 enddo  ! itn loop
      endif    ! (nvdep.eq.4) 
       
      if(nvdep.eq.5)then  !  newer method similar to depts3d  itn=2
        do iq=1,ifull
         sd(iq,1) =.25*sdot(iq,2)     ! quadratic right at end
         sd(iq,kl)=.25*sdot(iq,kl)    ! quadratic right at end
        enddo  ! iq loop
c       interpolate sdot to full-level sd with cubic polynomials
        do k=2,kl-1
         do iq=1,ifull
          sd(iq,k)=(9.*(sdot(iq,k)+sdot(iq,k+1))
     .                 -sdot(iq,k-1)-sdot(iq,k+2))/16.   
         enddo    ! iq loop
        enddo     ! k loop
	 ders=-tfact*sd  ! 3D
        do k=1,kl    
         do iq=1,ifull
          st(iq,k)=k+ders(iq,k) 
         enddo  ! iq loop
        enddo   ! k loop
	 itn=2   ! spell it out as for nvdep=3
         do k=1,kl
          do iq=1,ifull
           wrk1(iq,k)=.5*tfact*(sdot(iq,k+1)**2-sdot(iq,k)**2)/itn
          enddo  ! iq loop
         enddo   ! k loop
         ders=wrk1   ! 3D
         st=st+ders  ! 3D
 	 itn=3
         do k=2,kl-1
          do iq=1,ifull
           wrk1(iq,k)=-.5*tfact*sd(iq,k)*(ders(iq,k+1)-ders(iq,k-1))/itn
          enddo  ! iq loop
         enddo   ! k loop
         do iq=1,ifull
!         assume ders=0 off ends and varies quadratically to 1st level	
          wrk1(iq,1)=-tfact*sd(iq,1)*4.*ders(iq,1)/itn
          wrk1(iq,kl)=tfact*sd(iq,kl)*4.*ders(iq,kl)/itn
         enddo  ! iq loop
         ders=wrk1   ! 3D
         st=st+ders  ! 3D
      endif    ! (nvdep.eq.5) 
       
      if(nvdep.eq.6)then  !  newer method similar to depts3d  
        do iq=1,ifull
         sd(iq,1) =.25*sdot(iq,2)     ! quadratic right at end
         sd(iq,kl)=.25*sdot(iq,kl)    ! quadratic right at end
        enddo  ! iq loop
c       interpolate sdot to full-level sd with cubic polynomials
        do k=2,kl-1
         do iq=1,ifull
          sd(iq,k)=(9.*(sdot(iq,k)+sdot(iq,k+1))
     .                 -sdot(iq,k-1)-sdot(iq,k+2))/16.   
         enddo    ! iq loop
        enddo     ! k loop
	 ders=-tfact*sd  ! 3D
        do k=1,kl-1    
         do iq=1,ifull
          dersh(iq,k)=-tfact*sdot(iq,k+1) ! N.B. dersh(,0)=dersh(,kl)=0
         enddo  ! iq loop
        enddo   ! k loop
        do k=1,kl    
         do iq=1,ifull
          st(iq,k)=k+ders(iq,k) 
         enddo  ! iq loop
        enddo   ! k loop
	 itn=2   ! spell it out 
        do iq=1,ifull
         wrk1(iq,1)=-tfact*sdot(iq,2)*(ders(iq,2)-ders(iq,1))/itn  ! new dersh
         ders(iq,1)=-tfact*sd(iq,1)*dersh(iq,1)/itn    ! new ders
        enddo  ! iq loop
        do k=2,kl-1
         do iq=1,ifull  ! for new dersh, ders
          wrk1(iq,k)=-tfact*sdot(iq,k+1)*(ders(iq,k+1)-ders(iq,k))/itn  
          ders(iq,k)=-tfact*sd(iq,k)*(dersh(iq,k)-dersh(iq,k-1))/itn  ! new ders
         enddo  ! iq loop
        enddo   ! k loop
        do iq=1,ifull
         ders(iq,kl)=tfact*sd(iq,kl)*dersh(iq,kl-1)/itn  ! new ders
        enddo  ! iq loop
        st=st+ders  ! 3D
  	 itn=3
        do iq=1,ifull
         ders(iq,1)=-tfact*sd(iq,1)*wrk1(iq,k)/itn    ! new ders
        enddo  ! iq loop
        do k=2,kl-1
         do iq=1,ifull
          ders(iq,k)=-tfact*sd(iq,k)*(wrk1(iq,k)-wrk1(iq,k-1))/itn  ! new ders
         enddo  ! iq loop
        enddo   ! k loop
        do iq=1,ifull
         ders(iq,kl)=tfact*sd(iq,kl)*wrk1(iq,kl-1)/itn  ! new ders
        enddo  ! iq loop
        st=st+ders  ! 3D
      endif    ! (nvdep.eq.6) 
       
      if(nvdep.eq.7)then  !  iterative
        do iq=1,ifull
         sd(iq,1) =.25*sdot(iq,2)     ! quadratic right at end
         sd(iq,kl)=.25*sdot(iq,kl)    ! quadratic right at end
         sdotder(iq,1)   =0.
         sdotder(iq,kl+1)=0.
        enddo  ! iq loop
c       interpolate sdot to full-level sd with cubic polynomials
        do k=2,kl-1
         do iq=1,ifull
          sd(iq,k)=(9.*(sdot(iq,k)+sdot(iq,k+1))
     .                 -sdot(iq,k-1)-sdot(iq,k+2))/16.   
         enddo  ! iq loop
        enddo   ! k loop
        do k=2,kl
         do iq=1,ifull
          sdotder(iq,k)=.5*(sdot(iq,k+1)-sdot(iq,k-1))
         enddo  ! iq loop
        enddo   ! k loop
        do k=1,kl    
         do iq=1,ifull
          st(iq,k)=k -tfact*sd(iq,k)        ! 1st guess
         enddo  ! iq loop
        enddo   ! k loop
        if(diag.or.ntest.eq.1)then
          print *,'1st guess:  '
          write (6,"('st  ',9f8.3/4x,9f8.3)") (st(idjd,kk),kk=1,kl)
        endif
        do k=1,kl    
         do iq=1,ifull
!	   following kdel and st apply to sdot interp
          st(iq,k)=max(1.,min(st(iq,k)+.5,kl+1.))
          kdel(iq,k)=max(1,min(int(st(iq,k)),kl))  ! 1 <= kdel <= kl
          st(iq,k)=st(iq,k)-kdel(iq,k)             ! 0 <= st   <= 1.
         enddo  ! iq loop
         do iq=1,ifull  ! calculate upstream values of sdot
          kk=kdel(iq,k)
          gwrk(iq,k)=sdot(iq,kk)+st(iq,k)*(sdotder(iq,kk)
     .                      +st(iq,k)*(3.*(sdot(iq,kk+1)-sdot(iq,kk))
     .                               -2.*sdotder(iq,kk)-sdotder(iq,kk+1)
     .                      +st(iq,k)*(2.*(sdot(iq,kk)-sdot(iq,kk+1))
     .                              +sdotder(iq,kk)+sdotder(iq,kk+1) )))     
         enddo  ! iq loop
        enddo   ! k loop
        do k=1,kl    
         do iq=1,ifull
          st(iq,k)=k -.5*tfact*(sd(iq,k)+gwrk(iq,k))        ! 2nd guess
         enddo  ! iq loop
        enddo   ! k loop
        if(diag.or.ntest.eq.1)then
          print *,'2nd guess:  '
          write (6,"('st  ',9f8.3/4x,9f8.3)") (st(idjd,kk),kk=1,kl)
        endif
        do k=1,kl    
         do iq=1,ifull
!	   following kdel and st apply to sdot interp
          st(iq,k)=max(1.,min(st(iq,k)+.5,kl+1.))
          kdel(iq,k)=max(1,min(int(st(iq,k)),kl))  ! 1 <= kdel <= kl
          st(iq,k)=st(iq,k)-kdel(iq,k)             ! 0 <= st   <= 1.
         enddo  ! iq loop
         do iq=1,ifull  ! calculate upstream values of sdot
          kk=kdel(iq,k)
          gwrk(iq,k)=sdot(iq,kk)+st(iq,k)*(sdotder(iq,kk)
     .                      +st(iq,k)*(3.*(sdot(iq,kk+1)-sdot(iq,kk))
     .                               -2.*sdotder(iq,kk)-sdotder(iq,kk+1)
     .                      +st(iq,k)*(2.*(sdot(iq,kk)-sdot(iq,kk+1))
     .                              +sdotder(iq,kk)+sdotder(iq,kk+1) )))     
         enddo  ! iq loop
        enddo   ! k loop
        do k=1,kl    
         do iq=1,ifull
          st(iq,k)=k -.5*tfact*(sd(iq,k)+gwrk(iq,k))        ! 3rd guess
         enddo  ! iq loop
        enddo   ! k loop
      endif   ! (nvdep.eq.7) 
       
      if(diag.or.ntest.eq.1)then
         print *,'vadv30  ktau,nvint,nvdep:  ',
     .                    ktau,nvint,nvdep
         write (6,"('sdot',9f8.3/4x,9f8.3)") (sdot(idjd,kk),kk=1,kl)
         write (6,"('st  ',9f8.3/4x,9f8.3)") (st(idjd,kk),kk=1,kl)
       endif

!     here transform array st (i.e. nvad=7)
!        this should do Bessell better, with zero gradient top & bottom      
        do k=1,kl
         do iq=1,ifull
             st(iq,k)=min( max(st(iq,k),1.) ,real(kl))  ! 1 <= st   <= kl
             kdel(iq,k)=min(max(int(st(iq,k)),1),kl-1)  ! 1 <= kdel <= kl-1
             st(iq,k)=st(iq,k)-kdel(iq,k)               ! 0 <= st   <= 1.
         enddo  ! iq loop
        enddo   ! k loop
        if(diag.or.ntest.eq.1)then
         write (6,"('Bess st',9f8.3/4x,9f8.3)") (st(idjd,kk),kk=1,kl)
          print *,'new kdel ',(kdel(idjd,k),k=1,kl)
        endif
        call vadvbess(tin,tout,st,kdel,1)                          
        call vadvbess(uin,uout,st,kdel,2)                          
        call vadvbess(vin,vout,st,kdel,2)                          
        call vadvbess(qg,qg,st,kdel,3)                          
        if(ilt.gt.1)then
          do ntr=1,ntrac
           call vadvbess(tr(1,1,ntr),tr(1,1,ntr),st,kdel,3)           ! tr next
          enddo
        endif   ! (ilt.gt.1)
	 return
      end

c     this one also returns t, u, v, q -  but tendencies formed in vadv30
!     watch out for replacing qg in place!
      subroutine vadvbess(t,tout,st,kdel,ifield)
      parameter (ntopp=1)  ! 1 for 1-sided gradient at top & bottom full-levels
c                          ! 2 for zero gradient at top & bottom full-levels
      include 'newmpar.h'
      include 'parm.h'
      include 'parmvert.h'
      include 'sigs.h'
      real t(ifull,kl),tout(ifull,kl),st(ifull,kl)
      dimension kdel(ifull,kl)
      common/work3/tgrad(ifull,kl),toutt(ifull,kl),dum(3*ijk)
c     st() is the sigma displacement array
      if(ktau.eq.1)then
        print *,'in vadvbess with ntopp = ',ntopp
      endif

      if(ntopp.eq.1)then  ! 1-sided
        do iq=1,ifull
         tgrad(iq,1)=t(iq,2)-t(iq,1)
         tgrad(iq,kl)=t(iq,kl)-t(iq,kl-1)
        enddo
      endif  ! (ntopp.eq.2)
      if(ntopp.eq.2)then
        do iq=1,ifull
         tgrad(iq,1)=0.
         tgrad(iq,kl)=0.
        enddo
      endif  ! (ntopp.eq.2)

      do k=2,kl-1
       do iq=1,ifull
        tgrad(iq,k)=.5*(t(iq,k+1)-t(iq,k-1))
       enddo   ! iq loop
      enddo    ! k loop

      do k=1,kl
       do iq=1,ifull
        kk=kdel(iq,k)
        toutt(iq,k)=t(iq,kk)+st(iq,k)*(tgrad(iq,kk)
     .                      +st(iq,k)*(3.*(t(iq,kk+1)-t(iq,kk))
     .                                 -2.*tgrad(iq,kk)-tgrad(iq,kk+1)
     .                      +st(iq,k)*(2.*(t(iq,kk)-t(iq,kk+1))
     .                                 +tgrad(iq,kk)+tgrad(iq,kk+1) )))     
       enddo
      enddo
      if(diag)then
        print *,'t in  ',(t(idjd,k),k=1,kl)
        print *,'tgrad ',(tgrad(idjd,k),k=1,kl)
        print *,'toutt ',(toutt(idjd,k),k=1,kl)
      endif

!     can impose non-negative constraint here
      if(ifield.eq.3)then
        tout=max(toutt,0.)
      else
        tout=toutt
      endif

      if(diag)then
        print *,'tout  ',(tout (idjd,k),k=1,kl)
      endif
      return
      end
