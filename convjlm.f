! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2018 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
      module convjlm_m

      implicit none

      private
      public convjlm, convjlm_init
      
      integer, save :: k500,k600,k700,k900,k980,klon2,komega
!$acc declare create(k500,k600,k700,k900,k980,klon2,komega)
      integer, save :: mcontlnd,mcontsea           
!$acc declare create(mcontlnd,mcontsea)
      real, save :: convt_frac,tied_a,tied_b
!$acc declare create(convt_frac,tied_b)
      real, dimension(:), allocatable, save ::  entrainn
      real, dimension(:), allocatable, save ::  timeconv
      real, dimension(:), allocatable, save ::  alfin
      real, dimension(:,:), allocatable, save :: downex,upin,upin4
      real, dimension(:,:,:), allocatable, save :: detrarr
!$acc declare create(upin,downex,detrarr)

      contains

      subroutine convjlm_init
      
      use cc_mpi, only : myid, ccmpi_abort, mydiag
      use map_m
      use newmpar_m, only : ifull, kl
      use parm_m
      use sigs_m
      use soil_m

      implicit none

      include 'kuocom.h'   ! kbsav,ktsav,convfact,convpsav,ndavconv
      
      integer iq,k,kt,nlayers,ntest,kb
      real frac,summ,sumb
      parameter (ntest=0)      ! 1 or 2 to turn on; -1 for ldr writes
      integer kpos(1)

        kpos=minloc(abs(sig-.98)) ! finds k value closest to sig=.98  level 2 for L18 & L27
        k980=kpos(1)
        kpos=minloc(abs(sig-.9)) ! finds k value closest to sig=.9
        k900=kpos(1)
        kpos=minloc(abs(sig-.7)) ! finds k value closest to sig=.7
        k700=kpos(1)
        kpos=minloc(abs(sig-.6)) ! finds k value closest to sig=.6
        k600=kpos(1)
!$acc update device(k600,k700,k900,k980)
        komega=1
        if(dsig4>.4)then
          kpos=minloc(abs(sig-dsig4)) ! finds k value closest to dsig4
          komega=kpos(1)
        endif
        k500=1
        do while(sig(k500)>0.5)
          k500=k500+1
        enddo
        k500=k500-1    ! level just below .5
!$acc update device(k500,komega)
        if (myid==0) then
          write(6,*) 'k980 ',k980,sig(k980)
          write(6,*) 'k900 ',k900,sig(k900)
          write(6,*) 'k700 ',k700,sig(k700)
          write(6,*) 'k600 ',k600,sig(k600)
          write(6,*) 'k500 ',k500,sig(k500)
          write(6,*) 'komega',komega,sig(komega)
        end if
        allocate(timeconv(ifull))  ! init for ktau=1 (allocate needed for mdelay>0)
        allocate(entrainn(ifull))  ! init for ktau=1 (allocate needed for nevapcc.ne.0)
        allocate(alfin(ifull))     ! init for ktau=1
        allocate(upin(kl,kl))
        allocate(upin4(kl,kl))
        allocate(downex(kl,kl))
        allocate(detrarr(kl,k500,kl))
        detrarr(:,:,:)=1.e20  ! in case someone uses other than methprec=0,4,5,6,7,8
        if(methprec.ne.0.and.(methprec<3.or.methprec>9))then
          write(6,*) "unsupported methprec in convjlm"
          call ccmpi_abort(-1)
        endif
         if(methdetr<0.and.entrain>0.)then
          write(6,*) "not permitted that methdetr<0.and entrain>0."
          call ccmpi_abort(-1)
        endif
       if(methprec==0)detrarr(:,:,:)=0.

      !----------------------------------------------------------------
       entrainn(:)=entrain  ! N.B. for nevapcc.ne.0, entrain handles type of scheme   
       if ((mbase<0.and.mbase.ne.-10).or.mbase>4
     &        .or.alflnd<0.or.alfsea<0)then
         write(6,*) "ERROR: negative alflnd and alfsea or"
          write(6,*) "unsupported mbase convjlm"
          call ccmpi_abort(-1)
       endif
       timeconv(:)=0.
       do k=1,kl
         if(sig(k)> .5)klon2=k
       enddo
!$acc update device(klon2)
       if(myid==0)write(6,*) 'klon2,k500',klon2,k500
!      precalculate detrainment arrays for various methprec
!      methprec gives the vertical distribution of the detrainment
!      methdetr gives the fraction of precip detrained (going into qxcess)
       if(methprec==4)then
       do kb=1,k500
       do kt=1,kl-1
        sumb=0.
        do k=kb+1,kt
         detrarr(k,kb,kt)=(.05+ (sig(kb)-sig(kt))*(sig(kb)-sig(k)))*
     &                      dsig(k)
         sumb=sumb+detrarr(k,kb,kt)
        enddo
        do k=kb+1,kt
         detrarr(k,kb,kt)=detrarr(k,kb,kt)/sumb
        enddo
       enddo
       enddo
       endif  ! (methprec==4)
       if(methprec==5)then  ! rather similar detrainment to old 8, but more general. Now usual
       do kb=1,k500
       do kt=1,kl-1
        sumb=0.
        do k=kb+1,kt
         detrarr(k,kb,kt)=(.1+ (sig(kb)-sig(kt))*(sig(kb)-sig(k)))*
     &                      dsig(k)
         sumb=sumb+detrarr(k,kb,kt)
        enddo
        do k=kb+1,kt
         detrarr(k,kb,kt)=detrarr(k,kb,kt)/sumb
        enddo
       enddo
       enddo
       endif  ! (methprec==5)
       if(methprec==6)then  ! cloud top values larger than for methprec=4, which is larger than methprec=5
       do kb=1,k500
       do kt=1,kl-1
        sumb=0.
        do k=kb+1,kt
         detrarr(k,kb,kt)=(sig(kb)-sig(k))*dsig(k)  ! i.e. (0.+ (sig(kb)-sig(kt))*(sig(kb)-sig(k)))*dsig(k)
         sumb=sumb+detrarr(k,kb,kt)
        enddo
        do k=kb+1,kt
         detrarr(k,kb,kt)=detrarr(k,kb,kt)/sumb
        enddo
       enddo
       enddo
       endif  ! (methprec==6)
       if(methprec==7)then
       do kb=1,k500
       do kt=1,kl-1
         frac=min(1.,(sig(kb)-sig(kt))/.6)
         nlayers=kt-kb
         summ=.5*nlayers*(nlayers+1)
        do k=kb+1,kt
         detrarr(k,kb,kt)=((1.-frac)*(kt-kb)+
     &        (2.*frac-1.)*(k-kb))/
     &       (((1.-frac)*(kt-kb)**2+    
     &        (2.*frac-1.)*summ)) 
        enddo
       enddo
       enddo   
       endif  ! (methprec==7)
       if(methprec==8)then
       do kb=1,k500
       do kt=1,kl-1
         nlayers=kt-kb
         summ=.5*nlayers*(nlayers+1)
        do k=kb+1,kt
         detrarr(k,kb,kt)=(k-kb)/(summ) 
        enddo
       enddo
       enddo   
       endif  ! (methprec==8)
      if(myid==0)then
        write (6,*) "following gives kb=1 values for methprec=",methprec
        do kt=2,kl-1
         write (6,"('detrarr2-kt',17f7.3)") (detrarr(k,1,kt),k=2,kt)  ! just printed for kb=1
        enddo
       endif

!     precalculate below-base downdraft & updraft environmental contrib terms (ktau=1)
       do kb=1,kl-1   ! original proposal, kscsea=0 (small upin near sfce, large downex)
        summ=0.
        sumb=0.
        do k=1,kb
         downex(k,kb)=(sig(k)-sigmh(kb+1))*dsig(k)
         summ=summ+downex(k,kb)
         upin(k,kb)=(sig(k)-1.)*dsig(k)
         sumb=sumb+upin(k,kb)
        enddo
        do k=1,kb
         downex(k,kb)=downex(k,kb)/summ
         upin(k,kb)=upin(k,kb)/sumb
        enddo
       enddo
!     precalculate below-base downdraft & updraft environmental contrib terms
      if(kscsea==-1.or.kscsea==-2)then   ! trying to be backward compatible for upin
!      this one gives linearly increasing mass fluxes at each half level      
       do kb=1,k500
        upin4(1,kb)=(1.-sigmh(1+1))/(1.-sigmh(kb+1))
        upin(1,kb)=upin4(1,kb)
        do k=2,kb
         upin4(k,kb)=(1.-sigmh(k+1))/(1.-sigmh(kb+1))
         upin(k,kb)=upin4(k,kb)-upin4(k-1,kb)
        enddo
        if(ntest==1.and.myid==0)then
          write (6,"('upinx kb ',i3,15f7.3)") kb,(upin(k,kb),k=1,kb)
        endif
       enddo
      endif  ! kscsea==-1.or.kscsea==-2  
      if(kscsea<=-2)then   ! trying to be backward compatible for downex
!      gives relatively larger downdraft fluxes near surface      
       do kb=1,k500
        summ=0.
        do k=1,kb
         downex(k,kb)=dsig(k)
         summ=summ+downex(k,kb)
        enddo
        do k=1,kb
         downex(k,kb)=downex(k,kb)/summ
        enddo
       enddo
      endif  ! kscsea<=-2
!$acc update device(downex,upin,detrarr)
      if(myid==0)then
        do kb=1,kl-1
         write (6,"('upin1-kb ',17f7.3)") (upin(k,kb),k=1,kb)
        if(ntest==1)write (6,"('downex',17f7.3)") (downex(k,kb),k=1,kb)
        enddo
        do kb=1,kl-1
         write (6,"('downex',17f7.3)") (downex(k,kb),k=1,kb)
        enddo
      endif

       if(nint(convtime)==-23)convtime=3030.6
       if(nint(convtime)==-24)convtime=4040.6
        if(convtime>100.)then   ! new general style  May 2014
!         1836.45 is old 36;  4040.6 is old -24; 2020.001 is old .33
          mcontlnd=int(.01*convtime)                 ! in minutes
          mcontsea=int(convtime-100.*real(mcontlnd)) ! in minutes
          convt_frac=convtime-100*mcontlnd-mcontsea  ! changeover sigma value of cloud thickness
          convt_frac=max(convt_frac,1.e-7)  ! allows for zero entry
         elseif(convtime<-100.)then
!          mcontsea=-.01*convtime         ! in minutes
!          mcontlnd=-convtime-100*mcontsea ! fg value in W/m2
!          convt_frac=-convtime-100*mcontsea-mcontlnd  ! changeover sigma value of cloud thickness
!          convt_frac=max(convt_frac,1.e-7)  ! allows for zero entry
          mcontlnd=int(-.01*convtime)                 ! in minutes
          mcontsea=int(-convtime-100.*real(mcontlnd)) ! in minutes
!         convt_frac=100.*(-convtime-100*mcontlnd-mcontsea)  ! fg value in W/m2
          convt_frac=-convtime-100*mcontlnd-mcontsea  ! changeover sigma value of cloud thickness
          convt_frac=max(convt_frac,1.e-7)  ! allows for zero entry
         elseif(convtime>0.)then  ! using old values provided as hours
          mcontlnd=nint(60.*convtime)         ! in minutes
          mcontsea=nint(60.*convtime)         ! in minutes
          convt_frac=.001
          convtime=100*mcontlnd + mcontsea
         else
           write(6,*) 'unsupported convtime value'
           call ccmpi_abort(-1)
        endif  ! (convtime >100) .. else ..
!$acc update device(mcontlnd,mcontsea,convt_frac)
        if(mydiag)write(6,*) 'convtime,mcontlnd,mcontsea,convt_frac',
     &           convtime,mcontlnd,mcontsea,convt_frac

        where(land(1:ifull))
          alfin(:)=alflnd
        elsewhere
          alfin(:)=alfsea
        end where
        
        if(tied_over>0.)then   ! e.g. 2626.   b then a.  2600 to get old -26
          tied_b=int(tied_over/100.)
          tied_a=tied_over-100.*tied_b
        else
          tied_b=abs(tied_over)
          tied_a=0.
        endif
!$acc update device(tied_b)
        if (myid==0) then
          write(6,*) 'ds,tied_over,tied_a,tied_b',
     &                ds,tied_over,tied_a,tied_b
        end if
       if(tied_a>1.)then  ! e.g. 20 or 26  
!         alfin may be reduced over land and sea for finer resolution than 200 km grid            
          do iq=1,ifull
            summ=ds/(em(iq)*208498.)
            alfin(iq)=1.+(alfin(iq)-1.) *
     &      (1.+tied_a)*summ/(1.+tied_a*summ)
            !if(iq<200)print *,'iq,alfin',iq,alfin(iq)
!           tied_over=26 gives factor [1, .964, .900, .794, .529] for ds = [200, 100, 50, 25, 8} km     
!           tied_over=10 gives factor [1, .917, .786, .611, .306] for ds = [200, 100, 50, 25, 8} km     
          enddo
        endif  ! (tied_a>1.)
      !----------------------------------------------------------------

      end subroutine convjlm_init
      
      subroutine convjlm
      
      use arrays_m   
      use aerosolldr
      use cc_mpi, only : mydiag
      use cc_omp
      use cfrac_m
      use extraout_m
      use kuocomb_m
      use liqwpar_m  ! ifullw
      use map_m
      use morepbl_m
      use newmpar_m
      use parm_m
      use prec_m
      use sigs_m
      use soil_m
      use tracers_m  ! ngas, nllp, ntrac
      use vvel_m
      use work2_m   ! for wetfa!    JLM

      implicit none

      include 'kuocom.h'   ! kbsav,ktsav,convfact,convpsav,ndavconv

      integer :: tile, is, ie
      integer :: idjd_t
      integer, dimension(imax)          :: lkbsav, lktsav
      real, dimension(imax,kl,naero)    :: lxtg
      real, dimension(imax,kl,ntrac)    :: ltr
      real, dimension(imax,kl)          :: ldpsldt, lt, lqg
      real, dimension(imax,kl)          :: lfluxtot
      real, dimension(imax,kl)          :: lqlg, lu, lv, lqfg
      real, dimension(imax,kl)          :: lcfrac
      real, dimension(imax,ndust)       :: ldustwd
      real, dimension(imax)             :: lso2wd, lso4wd
      real, dimension(imax)             :: lbcwd, locwd, lsaltwd
      logical :: mydiag_t

!$omp  do schedule(static) private(is,ie),
!$omp& private(ldpsldt,lt,lqg,lfluxtot),
!$omp& private(lxtg,lso2wd,lso4wd,lbcwd,locwd,ldustwd,lsaltwd),
!$omp& private(lqlg,lqfg,lcfrac,lu,lv,ltr,idjd_t,mydiag_t)
!!$acc parallel loop copy(t,qg,qlg,qfg,u,v,xtg,dustwd,so2wd,so4wd),
!!$acc&  copy(bcwd,ocwd,saltwd,tr,precc,precip,timeconv,kbsav,ktsav),
!!$acc&  copyin(dpsldt,cfrac,alfin,ps,pblh,fg,wetfac,land,entrainn),
!!$acc&  copyin(em,sgsave),
!!$acc&  copyout(convpsav,cape,condc,condx,conds,condg)
!!$acc&  private(ldpsldt,lt,lqg,lqlg,lqfg,lcfrac,lu,lv),
!!$acc&  private(lxtg,ldustwd,lso2wd,lso4wd,lbcwd,locwd,lsaltwd),
!!$acc&  private(ltr,lfluxtot)
      do tile=1,ntiles
        is=(tile-1)*imax+1
        ie=tile*imax
        
        idjd_t = mod(idjd-1,imax)+1
        mydiag_t = ((idjd-1)/imax==tile-1).and.mydiag

        ldpsldt   = dpsldt(is:ie,:)
        lt        = t(is:ie,:)
        lqg       = qg(is:ie,:)
        lqlg      = qlg(is:ie,:)
        lqfg      = qfg(is:ie,:)
        lcfrac    = cfrac(is:ie,:)
        lu        = u(is:ie,:)
        lv        = v(is:ie,:) 
        if ( abs(iaero)>=2 ) then
          lxtg    = xtg(is:ie,:,:)
          ldustwd = dustwd(is:ie,:)
          lso2wd  = so2wd(is:ie)
          lso4wd  = so4wd(is:ie)
          lbcwd   = bcwd(is:ie)
          locwd   = ocwd(is:ie)
          lsaltwd = saltwd(is:ie)
        end if
        if ( ngas>0 ) then
          ltr = tr(is:ie,:,:)
        end if

        call convjlm_work(alfin(is:ie),ldpsldt,lt,lqg,
     &       ps(is:ie),lfluxtot,convpsav(is:ie),cape(is:ie),
     &       lxtg,lso2wd,lso4wd,lbcwd,locwd,ldustwd,lsaltwd,
     &       lqlg,condc(is:ie),precc(is:ie),condx(is:ie),
     &       conds(is:ie),condg(is:ie),precip(is:ie),
     &       pblh(is:ie),fg(is:ie),wetfac(is:ie),land(is:ie),
     &       entrainn(is:ie),lu,lv,timeconv(is:ie),em(is:ie),
     &       kbsav(is:ie),ktsav(is:ie),ltr,lqfg,lcfrac,sgsave(is:ie),
     &       idjd_t,mydiag_t,entrain,detrain,mbase,nbase,iterconv,
     &       nuvconv,alfsea,methdetr,methprec,fldown,alflnd,detrainx,
     &       sigkscb,dsig2,sigksct,rhcv,sig_ct,convtime,tied_con,
     &       mdelay,nevapcc,convfact,ncvcloud,ldr,rhmois,imax,kl)     ! jlm convective scheme

        t(is:ie,:)       = lt
        qg(is:ie,:)      = lqg
        qlg(is:ie,:)     = lqlg
        qfg(is:ie,:)     = lqfg
        fluxtot(is:ie,:) = lfluxtot
        u(is:ie,:)       = lu
        v(is:ie,:)       = lv
        if ( abs(iaero)>=2 ) then
          xtg(is:ie,:,:)  = lxtg
          dustwd(is:ie,:) = ldustwd
          so2wd(is:ie) = lso2wd
          so4wd(is:ie) = lso4wd
          bcwd(is:ie) = lbcwd
          ocwd(is:ie) = locwd
          saltwd(is:ie) = lsaltwd
        end if
        if ( ngas>0 ) then
          tr(is:ie,:,:) = ltr
        end if
       
      end do
!!$acc end parallel
!$omp end do nowait
      
      return
      end subroutine convjlm     ! jlm convective scheme

      subroutine convjlm_work(alfin,dpsldt,t,qg,ps,
     &       fluxtot,convpsav,cape,xtg,so2wd,so4wd,bcwd,ocwd,
     &       dustwd,saltwd,qlg,condc,precc,condx,conds,condg,precip,
     &       pblh,fg,wetfac,land,entrainn,u,v,timeconv,em,
     &       kbsav,ktsav,tr,qfg,cfrac,sgsave,
     &       idjd,mydiag,entrain,detrain,mbase,nbase,iterconv,
     &       nuvconv,alfsea,methdetr,methprec,fldown,alflnd,detrainx,
     &       sigkscb,dsig2,sigksct,rhcv,sig_ct,convtime,tied_con,
     &       mdelay,nevapcc,convfact,ncvcloud,ldr,rhmois,imax,kl)     ! jlm convective scheme
!$acc routine vector

!     version 1503e removes various fluxh, and keeps prior defn fluxv
!     version 1503d changes reflux(k) to fluxv()
!     version 1503c changes fluxv(k) to fluxv(k-1)
!     version 1503 replace 1+entracc by fluxv()      
!     version 1502 with preferred mbase = 4, nbase=-2
!     unused switches: detrainx, rhsat, shaltime 
!     unused switches if ksc=0:  kscmom, 
!     unused switches if nkuo=23:  rhsat
!     following will be unused when we remove methprec=3 option: dsig2, sigkscb, sigksct, tied_rh
!     has +ve fldownn depending on delta sigma; (-ve fldown descends from sig=.6))   
!     nevapcc option now affects entrainment
      use aerosolldr, only : itracso2,itracbc,itracoc,itracdu,ndust,
     &                       naero,convscav,itracsa,nsalt
      use cc_mpi, only : ccmpi_abort
!      use cc_omp, only : ntiles
      use const_phys
      use diag_m, only : maxmin
      use estab      
!      use newmpar_m
      use parm_m, only : ktau,dt,nmaxpr,diag,ds,iaero,qgmin
      use parmdyn_m
      use sigs_m
      use tracers_m, only : ngas,ntrac  ! ngas, nllp, ntrac
      
      implicit none
      
      integer itn,iq,k
     .       ,khalfp,kt
     .       ,nlayersp
     .       ,ntest,ntr,nums,nuv
     .       ,ka,kb
      integer idjd
      real convmax,deltaq,delq_av,delt_av
     .    ,den1,den2,den3,dprec
     .    ,facuv,fldownn,fluxup
     .    ,hbase,heatlev
     .    ,pwater,pwater0,qsk,rnrt_k
     .    ,summ,totprec,veldt
     .    ,pk,dz
     .    ,sumb
      real dtsol,detrainn
      logical mydiag
      parameter (ntest=0)      ! 1 or 2 to turn on; -1 for ldr writes
!                               -2,-3 for other detrainn test      
!     parameter (iterconv=3)  ! to kuocom.h
!     parameter (fldown=-.3)  ! to kuocom.h
!     parameter (detrain=.1)  ! to kuocom.h
!     parameter (alflnd=1.20) 
!     parameter (alfsea=1.05) 
!                                3:for RH-enhanced base
!     parameter (ncubase=0)    ! removed; 3 will tend to give kb_sav=1  Feb 2015  
!     parameter (mbase=0)     
!     parameter (methdetr=3)  ! fraction of precip detrained as fn of cloud depth; new scheme for -ve
!     parameter (methprec=5) !  gives the vertical distribution of the detrainment
!     parameter (nuvconv=0)    ! usually 0, >0 or <0 to turn on momentum mixing
!     parameter (nuv=0)        ! usually 0, >0 to turn on new momentum mixing (base layers too)
!     nevapls:  turn off/on ls evap - through parm.h; 0 off, 5 newer UK
      integer, intent(in) :: imax, kl
      integer knet
      real, dimension(imax,kl,naero), intent(inout)    :: xtg
      real, dimension(imax,kl,ntrac), intent(inout)    :: tr
      real, dimension(imax,kl), intent(in)         :: dpsldt
      real, dimension(imax,kl), intent(in)         :: cfrac
      real, dimension(imax,kl), intent(inout)          :: t
      real, dimension(imax,kl), intent(inout)          :: qg
      real, dimension(imax,kl), intent(inout)          :: qlg
      real, dimension(imax,kl), intent(inout)          :: qfg
      real, dimension(imax,kl), intent(inout)          :: u
      real, dimension(imax,kl), intent(inout)          :: v
      real, dimension(imax,kl), intent(out)            :: fluxtot
      real, dimension(imax,ndust), intent(inout)       :: dustwd
      real, dimension(imax), intent(in)                :: alfin
      real, dimension(imax), intent(in)                :: ps
      real, dimension(imax), intent(in)                :: pblh
      real, dimension(imax), intent(in)                :: fg
      real, dimension(imax), intent(in)                :: wetfac
      real, dimension(imax), intent(in)                :: entrainn
      real, dimension(imax), intent(in)                :: em
      real, dimension(imax), intent(in)                :: sgsave
      real, dimension(imax), intent(out)               :: convpsav
      real, dimension(imax), intent(out)               :: cape
      real, dimension(imax), intent(inout)             :: so2wd
      real, dimension(imax), intent(inout)             :: so4wd
      real, dimension(imax), intent(inout)             :: bcwd
      real, dimension(imax), intent(inout)             :: ocwd
      real, dimension(imax), intent(inout)             :: saltwd
      real, dimension(imax), intent(out)               :: condc
      real, dimension(imax), intent(inout)             :: precc
      real, dimension(imax), intent(out)               :: condx
      real, dimension(imax), intent(out)               :: conds
      real, dimension(imax), intent(out)               :: condg
      real, dimension(imax), intent(inout)             :: precip
      real, dimension(imax), intent(inout)             :: timeconv
      integer, dimension(imax), intent(inout)          :: kbsav
      integer, dimension(imax), intent(inout)          :: ktsav
      logical, dimension(imax), intent(in)             :: land
!from kuocom.h common block
      real, intent(in) :: entrain, detrain, alfsea, fldown, alflnd
      real, intent(in) :: detrainx, sigkscb, dsig2, sigksct, rhcv
      real, intent(in) :: sig_ct, convtime, tied_con, convfact, rhmois
      integer, intent(in) :: mbase, nbase, iterconv, nuvconv, methdetr
      integer, intent(in) :: methprec, mdelay, nevapcc, ncvcloud, ldr
      
      integer kbsav_ls(imax),kb_sav(imax),kt_sav(imax)
      integer kkbb(imax),kmin(imax)
      real, dimension(imax,kl) :: qqsav,qliqwsav,xtgscav
      real, dimension(kl) :: fscav,xtgtmp
      real, dimension(kl) :: rho,ttsto,qqsto,qqold,qlsto,qlold
      real, dimension(kl) :: qqnew,qqrem
      real, dimension(imax) :: conrev,alfqarr,omega,omgtst
      real, dimension(imax) :: convtim_deep,aa,bb
      real delq(imax,kl),dels(imax,kl),delu(imax,kl)
      real delv(imax,kl),dqsdt(imax,kl),es(imax,kl) 
      real fldow(imax),fluxq(imax)
      real fluxr(imax),fluxt(imax,kl),hs(imax,kl)  
      real phi(imax,kl),qdown(imax),qliqw(imax,kl),delqliqw(imax,kl)
      real entrsav(imax,kl),detxsav(imax,kl)
      real fluxh(imax,kl),fluxv(imax,0:kl-1)
      real qq(imax,kl),qs(imax,kl),qxcess(imax)
      real qbass(imax,kl-1)
      real rnrt(imax),rnrtc(imax),rnrtcn(imax)
      real s(imax,kl),tdown(imax),tt(imax,kl)
      real dsk(kl),h0(kl),q0(kl),t0(kl)  
      real qplume(imax,kl),splume(imax,kl)
      integer kdown(imax)
      real factr(imax)
      real fluxqs,fluxt_k(kl)
      integer, dimension(imax) :: kb_saved,kt_saved

      kb_saved(:)=kl-1
      kt_saved(:)=kl-1
      facuv=0.
      
      do k=1,kl
       dsk(k)=-dsig(k)    !   dsk = delta sigma (positive)
      enddo     ! k loop
      
      alfqarr(:)=alfin(:)
      omega(1:imax)=dpsldt(1:imax,komega)
       
      tt(:,:)=t(1:imax,:)       
      qq(:,:)=qg(1:imax,:)      
      phi(:,1)=bet(1)*tt(:,1)  ! moved up here May 2012
      do k=2,kl
       do iq=1,imax
        phi(iq,k)=phi(iq,k-1)+bet(k)*tt(iq,k)+betm(k)*tt(iq,k-1)
       enddo     ! iq loop
      enddo      ! k  loop

!      following defines kb_sav (as kkbb), method chosen by nbase
       kkbb(:)=1   ! kkbb is defined as top level within PBL 
       s(1:imax,1)=cp*tt(1:imax,1)+phi(1:imax,1)  ! dry static energy; only level 1 needed here for kkbb calc
       do k=2,k500   
        do iq=1,imax
        if(cp*tt(iq,k)+phi(iq,k)<s(iq,1)+max(cp,.1*cp*abs(nbase)))  ! uses 1 for nbase=-1,-2
     &                    kkbb(iq)=k  ! simple +.1*n deg s check for within PBL  
        enddo    ! iq loop
       enddo     ! k loop
       if(nbase==-12)then   ! retained for MJT
        do k=2,k500
         do iq=1,imax
!          find tentative cloud base ! 
!          (middle of k-1 level, uppermost level below pblh)
          if(phi(iq,k-1)<pblh(iq)*grav)kkbb(iq)=k-1
         enddo    ! iq loop
        enddo     ! k loop	 
       endif  ! (nbase==-12)
       if(nbase==-3.or.nbase==-4)then
         do k=2,k500   
          do iq=1,imax
!          find tentative cloud base ! 
!          (uppermost layer, with approx. bottom of layer below pblh)
           if(.5*(phi(iq,k-1)+phi(iq,k))<pblh(iq)*grav)kkbb(iq)=k
          enddo    ! iq loop
         enddo     ! k loop
       endif  ! (nbase==-3.or.nbase==-4)
          
        if(mbase==-10)then ! fg; qg1   
         do iq=1,imax
          k=kkbb(iq)
           if(fg(iq)>0.)alfqarr(iq)=alfqarr(iq)*                 !  mbase=-10; N.B. qs check done later with qbass
     &               max(qg(iq,1),qg(iq,k980),qg(iq,k))/qg(iq,k) ! MJT suggestion
         enddo 
        endif  ! (mbase-=-10)
          
      if(mbase==0)then ! fg; qg1     was usual from 28/2/14
         do iq=1,imax
          es(iq,1)=establ(tt(iq,1))
          pk=ps(iq)*sig(1)
          qs(iq,1)=.622*es(iq,1)/max(pk-es(iq,1),0.1)  
          k=kkbb(iq)  ! for alfqarr calc
c         N.B. if fg<=0, then alfqarr will keep its input value, e.g. 1.25         
         if(fg(iq)>0.)alfqarr(iq)=alfqarr(iq)*         !  mbase>=0;  N.B. qs check done later with qbass
     &            max(wetfac(iq)*qs(iq,1),qg(iq,k980),qg(iq,k))/qg(iq,k)
         enddo 
        endif  ! (mbase==0)
          
      if(mbase==1)then ! fg; qg1   similar to -10  
         do iq=1,imax
          es(iq,1)=establ(tt(iq,1))
          pk=ps(iq)*sig(1)
          qs(iq,1)=.622*es(iq,1)/max(pk-es(iq,1),0.1)  
          k=kkbb(iq)  ! for alfqarr calc
c         N.B. if fg<=0, then alfqarr will keep its input value, e.g. 1.25         
         if(fg(iq)>0.)alfqarr(iq)=alfqarr(iq)*         !  N.B. qs check done later with qbass
     &               max(qg(iq,1),qg(iq,k980),qg(iq,k))/qg(iq,k)
         enddo 
        endif  ! (mbase==1)
          
      if(mbase==2)then ! fg; qg1    fixed-up mbase=0
         do iq=1,imax
          es(iq,1)=establ(tt(iq,1))
          pk=ps(iq)*sig(1)
          qs(iq,1)=.622*es(iq,1)/max(pk-es(iq,1),0.1)  
          k=kkbb(iq)  ! for alfqarr calc
          if(fg(iq)>0.)alfqarr(iq)=max(wetfac(iq)*qs(iq,1),  !  N.B. qs check done later with qbass
     &           alfqarr(iq)*qg(iq,k980),alfqarr(iq)*qg(iq,k))/qg(iq,k)
         enddo 
        endif  ! (mbase==2)
          
      if(mbase==3)then ! fg; qg1   combined 1 & 2
         do iq=1,imax
          es(iq,1)=establ(tt(iq,1))
          pk=ps(iq)*sig(1)
          qs(iq,1)=.622*es(iq,1)/max(pk-es(iq,1),0.1)  
          k=kkbb(iq)  ! for alfqarr calc
          if(land(iq))then
            if(fg(iq)>0.)alfqarr(iq)=alfqarr(iq)*         !  same as mbase=1
     &            max(qg(iq,1),qg(iq,k980),qg(iq,k))/qg(iq,k)
          else
            if(fg(iq)>0.)alfqarr(iq)=max(qs(iq,1),        !  same as mbase=2
     &           alfqarr(iq)*qg(iq,k980),alfqarr(iq)*qg(iq,k))/qg(iq,k)
          endif
         enddo 
        endif  ! (mbase==3)
          
      if(mbase==4)then ! fg; qg_k   similar to combined 1 & 2 with qs limit - now preferred
         do iq=1,imax
          k=kkbb(iq)  ! for alfqarr calc
          es(iq,k)=establ(tt(iq,k))
          pk=ps(iq)*sig(k)
          qs(iq,k)=.622*es(iq,k)/max(pk-es(iq,k),0.1)  
          if(land(iq))then
            if(fg(iq)>0.)alfqarr(iq)=min(qs(iq,k),alfqarr(iq)*       
     &               max(qg(iq,1),qg(iq,k980),qg(iq,k)))
     &              /max(qg(iq,k),qgmin)
          else
            if(fg(iq)>0.)alfqarr(iq)=max(1.,qs(iq,k)
     &              /max(qg(iq,k),qgmin))
          endif
         enddo 
      endif  ! (mbase==4)

!      if(nproc==1.and.ntiles==1)write(6,*) 'max_alfqarr,alfin:',
!     &                    maxval(alfqarr),maxval(alfin)

#ifndef GPU
      if(ktau==1.and.mydiag)
     &   write(6,"('alfqarr',2f7.3)") alfqarr(idjd)
#endif

!     just does convective; L/S rainfall done later by LDR scheme
      qliqw(:,:)=0.  ! before itn
      kmin(:)=0
      conrev(1:imax)=1000.*ps(1:imax)/(grav*dt) ! factor to convert precip to g/m2/s
      rnrt(:)=0.       ! initialize large-scale rainfall array; before itn
      rnrtc(:)=0.      ! initialize convective  rainfall array; before itn
      kbsav_ls(:)=0    ! for L/S
!!      ktsav(:)=kl-1    ! preset value to show no deep or shallow convection
!!      kbsav(:)=kl-1    ! preset value to show no deep or shallow convection

      nuv=0
      fluxt(:,:)=0.     ! just for some compilers
      fluxtot(:,:)=0.   ! diag for MJT  - total mass flux at level k-.5
   
!__________beginning of convective calculation loop____#################################
      do itn=1,abs(iterconv)
!     calculate geopotential height 
      kb_sav(:)=kl-1   ! preset value for no convection
      kt_sav(:)=kl-1   ! preset value for no convection
      rnrtcn(:)=0.     ! initialize convective rainfall array (pre convpsav)
      delqliqw(:,:)=0. 
      convpsav(:)=0.
      dels(:,:)=1.e-20
      delq(:,:)=0.
      qqsav(:,:)=qq(:,:)       ! for convective scavenging of aerosols  MJT
      qliqwsav(:,:)=qliqw(:,:) ! for convective scavenging of aerosols  MJT
      if(nuvconv.ne.0)then 
        if(nuvconv<0)nuv=abs(nuvconv)  ! Oct 08 nuv=0 for nuvconv>0
        delu(:,:)=0.
        delv(:,:)=0.
        facuv=.1*abs(nuvconv) ! full effect for nuvconv=10
      endif

      do k=1,kl   
       do iq=1,imax
        es(iq,k) = establ(tt(iq,k))
       enddo  ! iq loop
      enddo   ! k loop
      do k=1,kl   
       do iq=1,imax
        pk=ps(iq)*sig(k)
!       qs(iq,k)=max(.622*es(iq,k)/(pk-es(iq,k)),1.5e-6)  
        qs(iq,k)=.622*es(iq,k)/max(pk-es(iq,k),0.1)  
        dqsdt(iq,k)=qs(iq,k)*pk*hlars/(tt(iq,k)**2*max(pk-es(iq,k),1.))
        s(iq,k)=cp*tt(iq,k)+phi(iq,k)  ! dry static energy
!       calculate hs
        hs(iq,k)=s(iq,k)+hl*qs(iq,k)   ! saturated moist static energy
       enddo  ! iq loop
      enddo   ! k loop
      hs(:,:)=max(hs(:,:),s(:,:)+hl*qq(:,:))   ! added 21/7/04
      if(itn==1)then  !--------------------------------------
        do k=1,kl-1
         qbass(:,k)=min(alfqarr(:)*qq(:,k),max(qs(:,k),qq(:,k)))  ! qbass only defined for itn=1
        end do
        if(nbase<=-1.and.nbase>=-3)then  ! i.e. -1, -2 or -3
         do k=1,kl-1
          do iq=1,imax
           ! above PBL; qbass only defined for itn=1   
           if(k>kkbb(iq))qbass(iq,k)=min(alfin(iq)*qq(iq,k),qbass(iq,k))
          enddo
         enddo  ! k loop
        endif  ! (nbase<=-1.and.nbase>=-3)         
c***    by defining qbass just for itn=1, convpsav does not converge as quickly as may be expected,
c***    as qbass may not be reduced by convpsav if already qs-limited    
c***    Also entrain may slow convergence   N.B. qbass only used in next few lines
      endif  ! (itn==1)  !--------------------------------------
      qplume(:,kl)=0. ! just for diag prints
      splume(:,kl)=0. ! just for diag prints
      do k=1,kl-1
       qplume(:,k)=min(alfqarr(:)*qq(:,k),qbass(:,k),   ! from Jan 08
     &         max(qs(:,k),qq(:,k)))   ! to avoid qb increasing with itn
       splume(:,k)=s(:,k)  
      enddo
      if(nbase<=-7)then
       do k=1,kl-1
       do iq=1,imax
        if(k>kkbb(iq))qplume(iq,k)=min(alfsea*qq(iq,k), 
     &         max(qs(iq,k),qq(iq,k)))   ! to avoid qb increasing with itn
        enddo
       enddo
      endif  ! (nbase<=-7)    
      
#ifndef GPU
      if(ktau==1.and.mydiag)then
      write(6,*) 'itn,iterconv,nuv,nuvconv ',itn,iterconv,nuv,nuvconv
      write(6,*) 'ntest,methdetr,methprec,detrain',
     .           ntest,methdetr,methprec,detrain
      write(6,*) 'fldown ',fldown
      write(6,*) 'alflnd,alfsea',alflnd,alfsea
      endif  ! (ktau==1.and.mydiag)
      
      if((ntest>0.or.nmaxpr==1).and.mydiag) then
        iq=idjd
        write (6,"('near beginning of convjlm; ktau',i5,' itn',i1)") 
     &                                         ktau,itn 
        write (6,"('rh   ',12f7.2/(5x,12f7.2))") 
     .             (100.*qq(iq,k)/qs(iq,k),k=1,kl)
        write (6,"('qs   ',12f7.3/(5x,12f7.3))") 
     .             (1000.*qs(iq,k),k=1,kl)
        write (6,"('qq   ',12f7.3/(5x,12f7.3))") 
     .             (1000.*qq(iq,k),k=1,kl)
        if(itn==1)write (6,"('qlg  ',12f7.3/(5x,12f7.3))") 
     .             (1000.*qlg(iq,k),k=1,kl)
        if(itn==1)write (6,"('qfg  ',12f7.3/(5x,12f7.3))") 
     .             (1000.*qfg(iq,k),k=1,kl)
        if(itn==1)write (6,"('qtot ',12f7.3/(5x,12f7.3))")
     &        (1000.*(qq(iq,k)+qlg(iq,k)+qfg(iq,k)),k=1,kl)
        if(itn==1)write (6,"('qtotx',12f7.3/(5x,12f7.3))")
     &        (10000.*dsk(k)*(qq(iq,k)+qlg(iq,k)+qfg(iq,k)),k=1,kl)
        pwater0=0.   ! in mm     
        do k=1,kl
         iq=idjd
         h0(k)=s(iq,k)/cp+hlcp*qq(iq,k)
         q0(k)=qq(iq,k)
         t0(k)=tt(iq,k)
         pwater0=pwater0-dsig(k)*qg(iq,k)*ps(iq)/grav
        enddo
!       following prints are just preliminary values. Others further down  
        write (6,"('qplume',f6.3,11f7.3/(5x,12f7.3))")1000.*qplume(iq,:)
        write (6,"('splume',12f7.2/(5x,12f7.2))")splume(iq,:)/cp 
        write (6,"('hplume',12f7.2/(5x,12f7.2))")
     &            (splume(iq,:)+hl*qplume(iq,:))/cp 
        write (6,"('tt    ',12f7.2/(5x,12f7.2))") tt(iq,:)
        write (6,"('s/cp  ',12f7.2/(5x,12f7.2))") s(iq,:)/cp
        write (6,"('h/cp  ',12f7.2/(5x,12f7.2))") 
     .             s(iq,:)/cp+hlcp*qq(iq,:)
        write (6,"('hb/cp',12f7.2/(5x,12f7.2))") 
     .             s(iq,:)/cp+hlcp*qplume(iq,:)
        write (6,"('hs/cp ',12f7.2/(5x,12f7.2))") hs(iq,:)/cp
        write (6,"('cfracc',12f7.3/6x,12f7.3)") cfrac(iq,:)
        write (6,"('k   ',12i7/(4x,12i7))") (k,k=1,kl)
        write (6,"('es     ',9f7.1/(7x,9f7.1))") es(iq,:)
        write (6,"('dqsdt6p',6p9f7.1/(7x,9f7.1))") dqsdt(iq,:)
        summ=0.
        do k=1,kl
         summ=summ-dsig(k)*(s(iq,k)/cp+hlcp*qq(iq,k))
        enddo
        write(6,*) 'h_sum   ',summ
      endif  ! ((ntest>0.or.nmaxpr==1).and.mydiag)
#endif

!     Following procedure chooses lowest valid cloud bottom (and base)
       kb_sav(:)=kl-1
      if(nbase>=0.or.nbase==-12)then    !  -12 retained for MJT
         do iq=1,imax
!         prescribe allowable cloud base 
           k=kkbb(iq)+1
          if(splume(iq,k-1)+hl*qplume(iq,k-1)>hs(iq,k)
     &       .and.qplume(iq,k-1)>max(qs(iq,k),qq(iq,k)))then   
             kb_sav(iq)=k-1
             kt_sav(iq)=k                   
           endif  
        enddo    ! iq loop
      elseif(nbase==-1.or.nbase==-4)then  ! older option
       kdown(:)=1     ! set tentatively as valid cloud bottom (not needed beyond this k loop)
       do k=klon2,2,-1  ! downwards for nbase=-4 to find cloud base **k loop**
         do iq=1,imax
!        find tentative cloud base, and bottom of below-cloud layer
         if(kdown(iq)>0)then   
          if(splume(iq,k-1)+hl*qplume(iq,k-1)>hs(iq,k)
     &       .and.qplume(iq,k-1)>max(qs(iq,k),qq(iq,k)))then   
             kb_sav(iq)=k-1
             kt_sav(iq)=k
             kdown(iq)=0
          endif  
         endif   ! (kdown(iq)>0)
        enddo    ! iq loop             
       enddo     ! ******************** k loop ****************************
      elseif(nbase==-2.or.nbase==-3)then  ! same as -1 but UPwards; now preferred
       kdown(:)=1     ! set tentatively as valid cloud bottom (not needed beyond this k loop)
       do k=2,k700 ! UPwards for nbase=-2 to find cloud base **k loop**
         do iq=1,imax
!        find tentative cloud base, and bottom of below-cloud layer
         if(kdown(iq)>0.and.k>kkbb(iq).and.
     &       splume(iq,k-1)+hl*qplume(iq,k-1)>hs(iq,k)
     &       .and.qplume(iq,k-1)>max(qs(iq,k),qq(iq,k)))then   
             kdown(iq)=0
             kb_sav(iq)=k-1
             kt_sav(iq)=k
             if(k>kkbb(iq)+1)alfqarr(iq)=alfin(iq)
         endif   ! (kdown(iq)>0...)
c         if(iq==idjd)print *,'k,kb_sav,s+q,hs,q,max,qbass',k,kb_sav(iq),
c     &          (splume(iq,k-1)+hl*qplume(iq,k-1))/cp,hs(iq,k)/cp,
c     & qplume(iq,k-1),max(qs(iq,k),qq(iq,k)),qbass(iq,k-1)
        enddo    ! iq loop             
       enddo     ! ******************** k loop ****************************
       elseif(nbase==-5)then   ! new JLMn
         aa(:)=splume(:,1)+hl*qplume(:,1)
         do k=2,k700+1  
          do iq=1,imax
!          find tentative cloud base as peak of hplume if above PBL
           bb(:)=splume(:,k)+hl*qplume(:,k)
           if(bb(iq)>aa(iq))then
              if(k>kkbb(iq))kkbb(iq)=k
              aa(iq)=bb(iq)
           endif
          enddo    ! iq loop
         enddo     ! k loop
         do iq=1,imax
          k=kkbb(iq)
          if(splume(iq,k)+hl*qplume(iq,k)>hs(iq,k+1)
     &       .and.qplume(iq,k)>max(qs(iq,k+1),qq(iq,k+1)))then     
             kb_sav(iq)=k
             kt_sav(iq)=k+1
          endif  
        enddo    ! iq loop             
       elseif(nbase==-6)then   ! new JLMn
         kdown(:)=1     ! set tentatively as valid cloud bottom (not needed beyond this k loop)
         do k=2,k700+1  
         do iq=1,imax
          if(k>=kkbb(iq).and.splume(iq,k)+hl*qplume(iq,k)>hs(iq,k+1).and
     &  .qplume(iq,k)>max(qs(iq,k+1),qq(iq,k+1)).and.kdown(iq)==1)then
             kb_sav(iq)=k
             kt_sav(iq)=k+1
             kdown(iq)=0
          endif  
        enddo    ! iq loop             
       enddo     ! k loop
      elseif(nbase==-7)then   ! new JLMn   like -6, but with omg & alf
         kdown(:)=1     ! set tentatively as first valid cloud bottom (not needed beyond this k loop)
         do k=2,k700+1  
         do iq=1,imax
          if(k>=kkbb(iq).and.splume(iq,k)+hl*qplume(iq,k)>hs(iq,k+1).and
     &  .qplume(iq,k)>max(qs(iq,k+1),qq(iq,k+1)).and.kdown(iq)==1
     &  .and.dpsldt(iq,k)<0.)then     
             kb_sav(iq)=k
             kt_sav(iq)=k+1
             kdown(iq)=0
          endif  
        enddo    ! iq loop             
       enddo     ! k loop
      elseif(nbase==-8)then   ! new JLMn   like -6 but with alf
         kdown(:)=1     ! set tentatively as first valid cloud bottom (not needed beyond this k loop)
         do k=2,k700+1  
         do iq=1,imax
          if(k>=kkbb(iq).and.splume(iq,k)+hl*qplume(iq,k)>hs(iq,k+1).and
     &  .qplume(iq,k)>max(qs(iq,k+1),qq(iq,k+1)).and.kdown(iq)==1)then
             kb_sav(iq)=k
             kt_sav(iq)=k+1
             kdown(iq)=0
          endif  
        enddo    ! iq loop             
       enddo     ! k loop
       elseif(nbase==-9)then   ! new JLMn  like -5 but with alf
         aa(:)=splume(:,1)+hl*qplume(:,1)
         do k=2,k700+1  
          do iq=1,imax
!          find tentative cloud base as peak of hplume if above PBL
           bb(:)=splume(:,k)+hl*qplume(:,k)
           if(bb(iq)>aa(iq))then
              if(k>kkbb(iq))kkbb(iq)=k
              aa(iq)=bb(iq)
           endif
          enddo    ! iq loop
         enddo     ! k loop
         do iq=1,imax
          k=kkbb(iq)
          if(splume(iq,k)+hl*qplume(iq,k)>hs(iq,k+1)
     &       .and.qplume(iq,k)>max(qs(iq,k+1),qq(iq,k+1)))then     
             kb_sav(iq)=k
             kt_sav(iq)=k+1
          endif  
        enddo    ! iq loop             
      elseif(nbase==-10)then   ! new JLMn   like -6 but with alf (and no q>qs test)
         kdown(:)=1     ! set tentatively as first valid cloud bottom (not needed beyond this k loop)
         do k=2,k700+1  
         do iq=1,imax
          if(k>=kkbb(iq).and.splume(iq,k)+hl*qplume(iq,k)>hs(iq,k+1).
     &          and.kdown(iq)==1)then
            kdown(iq)=0
            if(qplume(iq,k)>max(qs(iq,k+1),qq(iq,k+1)))then
             kb_sav(iq)=k
             kt_sav(iq)=k+1
            endif  
         endif  
        enddo    ! iq loop             
       enddo     ! k loop
      else  ! for all -ve values of nbase except -1,-2,-3,-4, -12
       kdown(:)=1     ! set tentatively as valid cloud bottom (not needed beyond this k loop)
       do k=2,k700+1  ! UPwards to find first suitable cloud base 
         do iq=1,imax
!        N.B. bit noisy if allow cloud base too much above kkbb
         if(kdown(iq)>0.
     &   and.k==kkbb(iq)+1.or.k==max(kkbb(iq),kb_saved(iq))+1)then 
          if(splume(iq,k-1)+hl*qplume(iq,k-1)>hs(iq,k)
     &       .and.qplume(iq,k-1)>max(qs(iq,k),qq(iq,k)))then     
             kb_sav(iq)=k-1
             kt_sav(iq)=k
             kdown(iq)=0
          endif  
         endif   ! (kdown(iq)>0)
        enddo    ! iq loop             
       enddo     ! ******************** k loop ****************************
       endif   !  nbase>=0 .. else ..
#ifndef GPU
       if(ntest>1.and.mydiag)then
        ka=1
        do iq=1,imax
         if(kb_sav(iq)>kkbb(iq).and.kb_sav(iq)<kl-1)
     & write(6,*) 'iq,kkbb,kb_sav,kt_sav',
     &             iq,kkbb(iq),kb_sav(iq),kt_sav(iq)
         if(kb_sav(iq)>ka.and.kb_sav(iq)<kl-1)then
           ka=kb_sav(iq)
           write(6,*) 'itn,iq,running_kb_sav_max,kt_sav',
     &                 itn,iq,ka,kt_sav(iq)
         endif
        enddo
       endif
#endif

      entrsav(:,:)=0.
      detxsav(:,:)=0.
      fluxv(:,:)=0.
      do iq=1,imax
        fluxv(iq,kb_sav(iq))=1.  ! unit reference base mass flux (at level k+.5)
      enddo
#ifndef GPU
      if(nmaxpr==1.and.mydiag)write(6,*) 
     &       'kb_saved,kt_saved,timeconva',
     &       kb_saved(idjd), kt_saved(idjd),timeconv(idjd)
#endif
     
      kdown(:)=1    ! set to show allowed to check for cloud top; removed mdelay stuff
      kt_sav(:)=kl-1  ! added 2/5/15 for safety with supersaturated layers
      do k=2,kl-2   ! upwards to find cloud top  
         do iq=1,imax
          if(k>kb_sav(iq).and.kdown(iq)==1)then 
            entrsav(iq,k)=-entrainn(iq)*dsig(k) ! (linear) entrained mass into plume_k
            if(entrain<0.)entrsav(iq,k)=entrain*fluxv(iq,k-1)*dsig(k) ! (non-linear) entrained mass into plume_k
!           fluxv is net mass into plume_k (assumed well mixed), and
!           subsequently is mass out top of layer for plume_k
!           N.B. qplume is usually > qs at each level
            fluxv(iq,k)=fluxv(iq,k-1)+entrsav(iq,k) ! plume reference mass flux
            qplume(iq,k)=(qplume(iq,k-1)*fluxv(iq,k-1)
     &          +entrsav(iq,k)*qq(iq,k))/fluxv(iq,k) 
            hbase=splume(iq,k-1)+hl*qplume(iq,k-1)
            if(hbase>hs(iq,k).and.   ! added for safety 22/4/15
     &       qplume(iq,k)>max(qq(iq,k),qs(iq,k)))then
             kt_sav(iq)=k
             splume(iq,k)=(splume(iq,k-1)*fluxv(iq,k-1)
     &                      +entrsav(iq,k)*s(iq,k))/fluxv(iq,k)
            else
              kdown(iq)=0
              entrsav(iq,k)=0.  ! not to mess up dels above cloud top
              fluxv(iq,k)=0.    ! not to mess up dels above cloud top
              if(kt_sav(iq)==kl-1)kb_sav(iq)=kl-1
            endif  ! (hbase>hs(iq,k))
#ifndef GPU
            if(ntest>0.and.entrain<0.and.iq==idjd.and.mydiag
     &         )then
              write(6,*) 'k,kb_sav,kt_sav,hbase/cp,hs/cp ',
     &                 k,kb_sav(iq),kt_sav(iq),hbase/cp,hs(iq,k)/cp
            endif
#endif
          endif   ! (k>kb_sav(iq))
         enddo    ! iq loop
      enddo     ! k loop
         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do iq=1,imax
c      next 4 lines ensure no entrainment into top layer (but maybe could allow) Mar 15
       fluxv(iq,kt_sav(iq))=0.
       entrsav(iq,kt_sav(iq))=0.
c      detxsav(iq,kt_sav(iq))=0.
       splume(iq,kt_sav(iq))=splume(iq,kt_sav(iq)-1)
       qplume(iq,kt_sav(iq))=qplume(iq,kt_sav(iq)-1)
      enddo    ! iq loop
      
#ifndef GPU
      if(ntest>0.and.mydiag)then
         print *,'before methdetr<0 part, kb_sav,kt_sav',
     &       kb_sav(idjd),kt_sav(idjd)
         write (6,"('fluxv_UP',15f6.3/(8x,15f6.3))")
     &             fluxv(idjd,1:kt_sav(idjd))
      endif
#endif
      if(methdetr<0)then  
!      calculate new detrainments and modify fluxq and entrsav      
       detxsav(:,:)=0.
       aa(:)=0.  ! this aa is for detxsum
       if(methdetr==-1)then
        do k=2,kl-2  
         do iq=1,imax
!         if(k>kb_sav(iq).and.qplume(iq,k)>max(qs(iq,k),qq(iq,k)))then  ! JLM max added for safety 21/4/15
          if(k>kb_sav(iq))then  
           if(k<=kt_sav(iq))
     &     detxsav(iq,k)=-dsig(k)*(sigmh(kb_sav(iq)+1)-sig(k))  ! (+ve) value at each level k 
           if(k<kt_sav(iq))aa(iq)=aa(iq)*(1.+entrain*dsig(k))- ! aa & entrain are -ve
     &              detxsav(iq,k)  ! flux changes at k+.5 due to modified entrainment & detrainment
c          if(iq==idjd)print *,'k,detxsav,entrdsig,aa',
c     &          k,detxsav(iq,k),entrain*dsig(k),aa(iq)
          endif
         enddo    ! iq loop
        enddo     ! k loop
       else  ! i.e. methdetr=-2 or -3                
        do k=2,kl-2  
         do iq=1,imax
!         if(k>kb_sav(iq).and.qplume(iq,k)>max(qs(iq,k),qq(iq,k)))then  ! JLM max added for safety 21/4/15
          if(k>kb_sav(iq))then  
           if(k<=kt_sav(iq))
     &     detxsav(iq,k)=-dsig(k)*(sigmh(kb_sav(iq)+1)-sig(k))**  ! (+ve) value at each level k 
     &     abs(methdetr)                ! **** use methdetr=-1, -2 or -3
           if(k<kt_sav(iq))aa(iq)=aa(iq)*(1.+entrain*dsig(k))- ! aa & entrain are -ve
     &              detxsav(iq,k)  ! flux changes at k+.5 due to modified entrainment & detrainment
c          if(iq==idjd)print *,'k,detxsav,entrdsig,aa',
c     &          k,detxsav(iq,k),entrain*dsig(k),aa(iq)
          endif
         enddo    ! iq loop
        enddo     ! k loop
       endif   !  (methdetr==-1) .. else ..
       
#ifndef GPU
       if(ntest>0.and.mydiag)then
         print *,'detxsav',detxsav(idjd,1:kt_sav(idjd))
         print *,'entrsav',entrsav(idjd,1:kt_sav(idjd))
         print *,'in methdetr<0 part, aa_a=',aa(idjd)
       endif
#endif
!      and modify fluxv and entrsav  *** assumes entrain<0.  ***     
       do iq=1,imax
         if(kb_sav(iq)<kl-1)aa(iq)=fluxv(iq,kt_sav(iq)-1)/  ! aa now is +ve scaling factor for detxsav
     &         (detxsav(iq,kt_sav(iq))-aa(iq))
       enddo    ! iq loop
       do k=2,kl-2  
        do iq=1,imax
          detxsav(iq,k)=aa(iq)*detxsav(iq,k)
         if(k>kb_sav(iq).and.k<kt_sav(iq))then  ! maybe can omit this "if" later
            entrsav(iq,k)=entrain*fluxv(iq,k-1)*dsig(k) ! (non-linear) entrained mass into plume_k
            fluxv(iq,k)=fluxv(iq,k-1)+entrsav(iq,k)-detxsav(iq,k)
          endif
        enddo    ! iq loop
       enddo     ! k loop
#ifndef GPU
       if(ntest>0.and.mydiag)then
         print *,'in methdetr<0 part, scaled aa=',aa(idjd)
         print *,'detxsav',detxsav(idjd,1:kt_sav(idjd))
         print *,'entrsav',entrsav(idjd,1:kt_sav(idjd))
         write (6,"('fluxv_up',15f6.3/(8x,15f6.3))")
     &             fluxv(idjd,1:kt_sav(idjd))
       endif
#endif
      endif  ! (methdetr<0)
#ifndef GPU
      if(nmaxpr==1.and.mydiag)then
        iq=idjd
        write (6,"('qplume',f6.3,11f7.3/(5x,12f7.3))")1000.*qplume(iq,:)
        write (6,"('splume',12f7.2/(5x,12f7.2))")splume(iq,:)/cp 
        write (6,"('entrsav ',15f6.3/(8x,15f6.3))")
     &            entrsav(iq,1:kt_sav(iq))
        write (6,"('detxsav ',15f6.3/(8x,15f6.3))")
     &            detxsav(iq,1:kt_sav(iq))
      endif
#endif
     
!     apply entrainment and detrainment effects to environment    
      if(methdetr<0)then
       if(methprec==3)then   ! JLMa start
        do iq=1,imax
          dz=sig(kb_sav(iq))-sig(kt_sav(iq))
          fluxr(iq)=dz  ! just for printout
          aa(iq)=min(1.,max(detrain, detrainx
     &  +max(((dz-sigkscb)*detrainx+(dsig2-dz))/(dsig2-sigkscb)
     &       -detrainx,0.)
     &  +min(((dz-dsig2)*detrain+(sigksct-dz)*detrainx)/(sigksct-dsig2)
     &       -detrainx,0.)   ))
        enddo
#ifndef GPU
        if(mydiag.and.nmaxpr==1)
     &   write(6,*) 'kb_sav,kt_sav,depth,aa',
     &   kb_sav(idjd),kt_sav(idjd),fluxr(idjd),aa(idjd)
#endif
       else
        do iq=1,imax
          ! typical detrain is .1 for deep clouds  
          ! as for methdetr=-3 here  
          aa(iq)=min( 1., detrain+(1.-detrain)*( 
     &   (.55-min(sig(kb_sav(iq))-sig(kt_sav(iq)), .55)) /(.6-.14))**3 )
        enddo
       endif
       do k=2,kl-2
        do iq=1,imax
          dels(iq,k)=-entrsav(iq,k)*s(iq,k)   ! entr into updraft
     &      +detxsav(iq,k)*splume(iq,k)
!         following line allows for qq>qs, to avoid drying layer
          qsk=max(qs(iq,k),qq(iq,k))
          delq(iq,k)=-entrsav(iq,k)*qq(iq,k)+detxsav(iq,k)*qsk  ! entr into updraft
c         rnrt_k=detxsav(iq,k)*max(0.,qplume(iq,k)-qsk)     ! not need as such a detxmax will be 0
          rnrt_k=detxsav(iq,k)*(qplume(iq,k)-qsk)
          dels(iq,k)=dels(iq,k)+hl*rnrt_k    ! corresponding precip. heating (as not done separately in plume)
!         part of this as detrained liquid water
!!!       detrainn=min( 1., detrain+(1.-detrain)*(                       ! typical detrain is .15 for deep clouds
!!!  &   (.55-min(sig(kb_sav(iq))-sig(kt_sav(iq)), .55)) /(.6-.14))**3 )  ! as for methdetr=-3 here
          delqliqw(iq,k)=delqliqw(iq,k)+aa(iq)*rnrt_k   ! JLMa  end
          rnrt_k= rnrt_k-delqliqw(iq,k)        
          rnrtcn(iq)=rnrtcn(iq)+rnrt_k
         enddo     ! iq loop
        enddo  ! k loop                        
       else   ! older methdetr>0 method
        do iq=1,imax
         detxsav(iq,kt_sav(iq))=fluxv(iq,kt_sav(iq)-1)  ! also true for methdetr<0 
        enddo
        do k=2,kl-2  
         do iq=1,imax
          dels(iq,k)=-entrsav(iq,k)*s(iq,k)   ! entr into updraft
     &      +detxsav(iq,k)*splume(iq,k)
          delq(iq,k)=-entrsav(iq,k)*qq(iq,k)  ! entr into updraft
!         following line allows for qq>qs, to avoid drying layer
          qsk=max(qs(iq,k),qq(iq,k))  
          rnrt_k=detxsav(iq,k)*max(0.,qplume(iq,k)-qsk)     ! think more about this max some time        
          delq(iq,k)=delq(iq,k)+detxsav(iq,k)*qsk  
          dels(iq,k)=dels(iq,k)+hl*rnrt_k    ! corresponding precip. heating
          rnrtcn(iq)=rnrtcn(iq)+rnrt_k
         enddo     ! iq loop
        enddo  ! k loop                        
      endif   ! (methdetr<0)  .. else
            
!     downdraft calculations
      do iq=1,imax
       ! for downdraft   
       kdown(iq)=min(kl,kb_sav(iq)+nint(.6+.75*(kt_sav(iq)-kb_sav(iq))))
!      following line may give kb_sav=kt_sav=kl, and smaller kdown, but does not matter       
       if(fldown<0.)kdown(iq)=min(kdown(iq),k600) ! makes downdraft below sig=.6
!      qsk is value entering downdraft, qdown is value leaving downdraft
!      tdown is temperature emerging from downdraft (at cloud base)       
       qsk=min(qs(iq,kdown(iq)),qq(iq,kdown(iq))) ! N.B. min for downdraft subs
       tdown(iq)=t(iq,kb_sav(iq))+(s(iq,kdown(iq))-s(iq,kb_sav(iq))
     .      +hl*(qsk-qs(iq,kb_sav(iq))))/(cp+hl*dqsdt(iq,kb_sav(iq)))
!      don't want moistening of cloud base layer (can lead to -ve fluxes) 
!      qdown(iq)=min( qq(iq,kb_sav(iq)) , qs(iq,kb_sav(iq))+
!    .           (tdown(iq)-t(iq,kb_sav(iq)))*dqsdt(iq,kb_sav(iq)) )
       qdown(iq)= qs(iq,kb_sav(iq))+
     .           (tdown(iq)-t(iq,kb_sav(iq)))*dqsdt(iq,kb_sav(iq)) 
       dprec=qdown(iq)-qsk                      ! to be mult by fldownn
!      typically fldown=.3 or -.3      
       fldownn=max(-abs(fldown),
     &                abs(fldown)*(sig(kb_sav(iq))-sig(kt_sav(iq))))
       totprec=rnrtcn(iq)-fldownn*dprec 
       if(tdown(iq)>t(iq,kb_sav(iq)).or.totprec<0.
     &    .or.sig(kb_sav(iq))-sig(kt_sav(iq))<.4)then   ! suppr. downdraft
         fldow(iq)=0.
       else
         fldow(iq)=fldownn
       endif
       rnrtcn(iq)=rnrtcn(iq)-fldow(iq)*dprec        ! already has dsk factor
#ifndef GPU
       if(ntest==1.and.iq==idjd.and.mydiag)then
         write(6,*)'qsk,rnrtcn,totprec',qsk,rnrtcn(iq),totprec 
         write(6,*) 'dprec,rnrtcn ',dprec,rnrtcn(iq)
       endif
#endif
!      add in downdraft contributions at level kdown
       dels(iq,kdown(iq))=dels(iq,kdown(iq))-fldow(iq)*s(iq,kdown(iq))
       delq(iq,kdown(iq))=delq(iq,kdown(iq))-fldow(iq)*qsk
       if(nuv>0)then
         delu(iq,kt_sav(iq))=u(iq,kb_sav(iq))  ! no entrainment effects included for u,v
         delv(iq,kt_sav(iq))=v(iq,kb_sav(iq))   ! so assuming unit upward flux through depth of plume
         delu(iq,kdown(iq))=delu(iq,kdown(iq))-fldow(iq)*u(iq,kdown(iq))
         delv(iq,kdown(iq))=delv(iq,kdown(iq))-fldow(iq)*v(iq,kdown(iq))
       endif  ! (nuv>0)
      enddo  ! iq loop

#ifndef GPU
      if(nmaxpr==1.and.mydiag)then
       iq=idjd
       write (6,"('hplume',12f7.2/(5x,12f7.2))")
     &            (splume(iq,:)+hl*qplume(iq,:))/cp 
       write (6,"('fluxv_up',15f6.3/(8x,15f6.3))")
     &             fluxv(iq,1:kt_sav(iq))
      endif
#endif
!     calculate environment fluxes, fluxh (full levels) and fluxv (now k+.5 levels)
!     N.B. fluxv usually <1 within downdraft layers
      fluxv(:,0)=0.  ! +ve for downwards subsident air calcs below cloud base
      do iq=1,imax
       if(kb_sav(iq)<kl)then
         fluxh(iq,kdown(iq))=fldow(iq)
         do k=1,  kb_sav(iq)   ! +ve into plume ###### BELOW cloud base ###### 
          fluxh(iq,k)=upin(k,kb_sav(iq))-fldow(iq)*downex(k,kb_sav(iq))
          fluxv(iq,k)=fluxv(iq,k-1)+fluxh(iq,k)  ! needed for subsidence calcs below cloud base
!         calculate emergent downdraft properties
          delq(iq,k)=fldow(iq)*downex(k,kb_sav(iq))*qdown(iq)
          dels(iq,k)=fldow(iq)*downex(k,kb_sav(iq))*
     &        (s(iq,kb_sav(iq))+cp*(tdown(iq)-t(iq,kb_sav(iq))))  ! correct
!         subtract contrib from cloud base layers into plume
          delq(iq,k)=delq(iq,k)-upin(k,kb_sav(iq))*qplume(iq,kb_sav(iq))
          dels(iq,k)=dels(iq,k)-upin(k,kb_sav(iq))*splume(iq,kb_sav(iq))
          if(nuv>0)then
             delu(iq,k)=-upin(k,kb_sav(iq))*u(iq,kb_sav(iq))
     &                 +fldow(iq)*downex(k,kb_sav(iq))*u(iq,kdown(iq))
             delv(iq,k)=-upin(k,kb_sav(iq))*v(iq,kb_sav(iq))
     &                 +fldow(iq)*downex(k,kb_sav(iq))*v(iq,kdown(iq))
          endif  ! (nuv>0)
         enddo   ! k loop             ######
         do  k=kb_sav(iq)+1,kdown(iq)-1  ! #### within cloud layers
          fluxv(iq,k)=fluxv(iq,k)-fldow(iq)
         enddo  ! k loop                         #### 
        endif    ! (kb_sav(iq)<kl)
       enddo     ! iq loop

#ifndef GPU
      if(diag.and.mydiag)then
       iq=idjd
       write (6,"('fluxv_dn',15f6.3/(8x,15f6.3))")
     &              fluxv(iq,1:kt_sav(iq))
       write (6,"('rh      ',15f6.3/(8x,15f6.3))")
     &         (qq(iq,k)/qs(iq,k),k=1,kt_sav(iq))
       write (6,"('fluxh   ',15f6.3/(8x,15f6.3))")
     &              fluxh(iq,1:kt_sav(iq))
       write(6,*) 'ktau,itn,fluxv,fldow',
     &             ktau,itn,fluxv(iq,kt_sav(iq)-1),fldow(iq)
       write(6,*) 'ktau,itn,kb_sav,kt_sav,kdown ',
     &          ktau,itn,kb_sav(iq),kt_sav(iq),kdown(iq)
       write(6,*) 'alfqarr,omega7,kb_sav',
     &             alfqarr(iq),1.e7*omega(iq),kb_sav(iq)
       write (6,"('s/cp    b',12f7.2/(5x,12f7.2))")s(iq,1:kt_sav(iq))/cp
       write (6,"('splume/cp',12f7.2/(5x,12f7.2))")
     &             splume(iq,1:kt_sav(iq))/cp
       write(6,*) 'qplume',qplume(iq,1:kt_sav(iq))*1.e3
       write (6,"('hplume/cp',12f7.2/(5x,12f7.2))") 
     &           splume(iq,1:kt_sav(iq))/cp+hlcp*qplume(iq,1:kt_sav(iq))
       write (6,"('hs/cp    ',12f7.2/(5x,12f7.2))") hs(iq,:)/cp
!      following just shows flux into kt layer, and fldown at base and top     
       write (6,"('delsa',9f8.0/(5x,9f8.0))")dels(iq,1:kt_sav(iq))
       write (6,"('delqa',3p9f8.3/(5x,9f8.3))")delq(iq,1:kt_sav(iq))
       write (6,"('dellq',3p9f8.3/(5x,9f8.3))")delqliqw(iq,1:kt_sav(iq))
      endif
#endif

!     subsidence effects
      do k=2,kl-1
       do iq=1,imax
         if(fluxv(iq,k-1)>0.)then  ! downwards
           dels(iq,k-1)=dels(iq,k-1)+fluxv(iq,k-1)*s(iq,k)
           delq(iq,k-1)=delq(iq,k-1)+fluxv(iq,k-1)*qq(iq,k)
           dels(iq,k)=dels(iq,k)-fluxv(iq,k-1)*s(iq,k)
           delq(iq,k)=delq(iq,k)-fluxv(iq,k-1)*qq(iq,k)
           if(nuv>0)then
             delu(iq,k-1)=delu(iq,k-1)+fluxv(iq,k-1)*u(iq,k)
             delv(iq,k-1)=delv(iq,k-1)+fluxv(iq,k-1)*v(iq,k)
             delu(iq,k)=delu(iq,k)-fluxv(iq,k-1)*u(iq,k)
             delv(iq,k)=delv(iq,k)-fluxv(iq,k-1)*v(iq,k)
           endif
         else  ! fluxv may be -ve for kb_sav and below
           dels(iq,k-1)=dels(iq,k-1)+fluxv(iq,k-1)*s(iq,k-1)
           delq(iq,k-1)=delq(iq,k-1)+fluxv(iq,k-1)*qq(iq,k-1)
           dels(iq,k)=dels(iq,k)-fluxv(iq,k-1)*s(iq,k-1)
           delq(iq,k)=delq(iq,k)-fluxv(iq,k-1)*qq(iq,k-1)
           if(nuv>0)then
             delu(iq,k-1)=delu(iq,k-1)+fluxv(iq,k-1)*u(iq,k-1)
             delv(iq,k-1)=delv(iq,k-1)+fluxv(iq,k-1)*v(iq,k-1)
             delu(iq,k)=delu(iq,k)-fluxv(iq,k-1)*u(iq,k-1)
             delv(iq,k)=delv(iq,k)-fluxv(iq,k-1)*v(iq,k-1)
           endif
        endif
       enddo
      enddo
#ifndef GPU
      if(ntest>0.and.mydiag)then
        iq=idjd
        write (6,"('delsb',9f8.0/(5x,9f8.0))")
     &              dels(iq,:)
        write (6,"('delqb',3p9f8.3/(5x,9f8.3))")
     &              delq(iq,:)
        write(6,*) "before diag print of dels,delh "
        iq=idjd
        nlayersp=max(1,nint((kt_sav(iq)-kb_sav(iq)-.1)/methprec)) ! round down
        khalfp=kt_sav(iq)+1-nlayersp
        write(6,*) 'kb_sav,kt_sav',
     .           kb_sav(iq),kt_sav(iq)
        write(6,*) 'khalfp,nlayersp ',khalfp,nlayersp 
        write (6,"('delsd',9f8.0/(5x,9f8.0))")
     &              dels(iq,:)
        write (6,"('delq*hl',9f8.0/(5x,9f8.0))")
     &              delq(iq,:)*2.5e6
        write (6,"('delh ',9f8.0/(5x,9f8.0))")
     &             (dels(iq,k)+hl*delq(iq,k),k=1,kl)
        write (6,"('delhb',9f8.0/(5x,9f8.0))")
     &             (dels(iq,k)+alfqarr(iq)*hl*delq(iq,k),k=1,kl)
        summ=0.
        do k=kb_sav(iq),kt_sav(iq)
         summ=summ+dels(iq,k)+hl*delq(iq,k)
        enddo
        write(6,*) 'qplume,sum_delh ',qplume(iq,kb_sav(iq)),summ
      endif  ! (ntest>0.and.mydiag)
#endif

!     calculate actual delq and dels
      do k=1,kl-1    
       do iq=1,imax
          delq(iq,k)=delq(iq,k)/dsk(k)
          dels(iq,k)=dels(iq,k)/dsk(k)
       enddo     ! iq loop
      enddo      ! k loop

#ifndef GPU
      if(diag.and.mydiag)then   ! JLM
        iq=idjd
        write(6,*) "before convpsav calc, after division by dsk"
        write (6,"('dels ',9f8.0/(5x,9f8.0))")
     &              dels(iq,:)
        write (6,"('delq3p',3p9f8.3/(7x,9f8.3))")
     &              delq(iq,:)
      endif  ! (diag.and.mydiag)   JLM
#endif

!----------------------------------------------------------------
!     calculate base mass flux 
      do iq=1,imax
       fluxq(iq) = 0. ! MJT suggestion   
       if(kb_sav(iq)<kl-1)then   
!          fluxq limiter: alfqarr*new_qq(kb)>=new_qs(kb+1)
!          note satisfied for M=0, so M from following eqn gives cutoff value
!          i.e. alfqarr*[qq+M*delq]_kb=[qs+M*dqsdt*dels/cp]_kb+1
           k=kb_sav(iq)
           fluxq(iq)=max(0.,(alfqarr(iq)*qq(iq,k)-qs(iq,k+1))/ 
     .     (dqsdt(iq,k+1)*dels(iq,k+1)/cp +alfqarr(iq)*abs(delq(iq,k))))
           if(delq(iq,k)>0.)fluxq(iq)=0.    ! delq should be -ve 15/5/08
           convpsav(iq)=fluxq(iq)
       endif    ! (kb_sav(iq)<kl-1)
      enddo     ! iq loop
      
      do k=kl-1,2,-1
       do iq=1,imax
!       want: new_hbas>=new_hs(k), i.e. in limiting case:
!       [h+alfsarr*M*dels+alfqarr*M*hl*delq]_base 
!                                          = [hs+M*dels+M*hlcp*dels*dqsdt]_k
!       Assume alfsarr=1.
!       pre 21/7/04 occasionally got -ve fluxt, for supersaturated layers
        if(k>kb_sav(iq).and.k<=kt_sav(iq))then
          fluxt(iq,k)=(splume(iq,k-1)+hl*qplume(iq,k-1)-hs(iq,k))/
     .       max(1.e-9,dels(iq,k)*(1.+hlcp*dqsdt(iq,k))       ! 0804 for max
     .                -dels(iq,kb_sav(iq))                    ! to avoid zero
     .                -alfqarr(iq)*hl*delq(iq,kb_sav(iq)) )   ! with real*4
          if(sig(k)<sig(kb_sav(iq))                 ! rhcv not RH but for depth-related flux test
     &       -rhcv*(sig(kb_sav(iq))-sig(kt_sav(iq))).and.    !  typically rhcv=.1; 0 for backward-compat
     &       fluxt(iq,k)<convpsav(iq))then  
!         if(fluxt(iq,k)<convpsav(iq))then  ! superseded by preceding line
            convpsav(iq)=fluxt(iq,k)
            kmin(iq)=k   ! level where fluxt is a minimum (diagnostic)
          endif    ! (sig(k)<sig(kb_sav(iq)-....and.)
          if(dels(iq,kb_sav(iq))+alfqarr(iq)*hl*delq(iq,kb_sav(iq))>=0.)
     &             convpsav(iq)=0. 
          if(dels(iq,kt_sav(iq))<=0.)convpsav(iq)=0.    ! JLM 1505 must stabilize
        endif   ! (k>kb_sav(iq).and.k<kt_sav(iq))
       enddo    ! iq loop
      enddo     ! k loop      
#ifndef GPU
      if(diag.and.mydiag)then    ! JLM
       do k=kl-1,2,-1
        iq=idjd
        if(k>kb_sav(iq).and.k<=kt_sav(iq))then
          den1=dels(iq,k)*(1.+hlcp*dqsdt(iq,k))
          den2=dels(iq,kb_sav(iq))
          den3=alfqarr(iq)*hl*delq(iq,kb_sav(iq))
          fluxt_k(k)=(splume(iq,k-1)+hl*qplume(iq,k-1)-hs(iq,k))/
     &                max(1.e-9,den1-den2-den3) 
         write(6,*)'k,dqsdt,den1,den2,den3,fluxt ',
     &              k,dqsdt(iq,k),den1,den2,den3,fluxt_k(k)
         write(6,*)'k,den,num',den1-den2-den3,
     &              k,splume(iq,k-1)+hl*qplume(iq,k-1)-hs(iq,k)     
        endif   ! (k>kb_sav(iq).and.k<kt_sav(iq))
       enddo     ! k loop      
      endif    ! (diag.and.mydiag)   JLM
#endif

#ifndef GPU
      if(ntest==2.and.mydiag)then     !######################
        convmax=0.
        do iq=1,imax
         if(convpsav(iq)>convmax.and.kb_sav(iq)==2)then
           write(6,*) 'ktau,iq,convpsav,fluxt3,kb_sav2,kt_sav ',
     &            ktau,iq,convpsav(iq),fluxt(iq,3),kb_sav(iq),kt_sav(iq)
           convmax=convpsav(iq)
         endif
        enddo
        iq=idjd
        k=kb_sav(iq)
        if(k<kl)then
          fluxqs=(alfqarr(iq)*qq(iq,k)-qs(iq,k+1))/
     &           (dqsdt(iq,k+1)*dels(iq,k+1)/cp -alfqarr(iq)*delq(iq,k))
          write(6,*) 'alfqarr(iq)*qq(iq,k) ',alfqarr(iq)*qq(iq,k)
          write(6,*) 'alfqarr(iq)*delq(iq,k) ',alfqarr(iq)*delq(iq,k)
          write(6,*) 'qs(iq,k+1) ',qs(iq,k+1)
          write(6,*) 'dqsdt(iq,k+1)*dels(iq,k+1)/cp ',
     &             dqsdt(iq,k+1)*dels(iq,k+1)/cp
          write(6,*) 'fluxqs ',fluxqs
        endif
        write(6,"('kmin,fluxq,convpsav',i3,5f9.5)") 
     &    kmin(iq),fluxq(iq),convpsav(iq)      
        write(6,"('delQ*dsk6p',6p9f8.3/(10x,9f8.3))")
     .    (convpsav(iq)*delq(iq,k)*dsk(k),k=1,kl)
        write(6,"('delt*dsk3p',3p9f8.3/(10x,9f8.3))")
     .   (convpsav(iq)*dels(iq,k)*dsk(k)/cp,k=1,kl)
        convmax=0.
        nums=0
        write(6,*) '    ktau   iq nums  kb kt     dsk    flt  flux',
     &         '  rhb  rht   rnrt    rnd'
        do iq=1,imax     
         if(kt_sav(iq)-kb_sav(iq)==1.and.convpsav(iq)>0.)then
           nums=nums+1
           if(convpsav(iq)>convmax)then
             convmax=convpsav(iq)
             write(6,"('bc  ',3i5,2i3,2x,2f7.3,f6.3,2f5.2,2f7.4)") 
     .       ktau,iq,nums,kb_sav(iq),kt_sav(iq),
     .       dsk(kb_sav(iq)),
     .       fluxt(iq,kt_sav(iq)),convpsav(iq),
     .       qq(iq,kb_sav(iq))/qs(iq,kb_sav(iq)),
     .       qq(iq,kt_sav(iq))/qs(iq,kt_sav(iq)),
     &       rnrtcn(iq),rnrtcn(iq)*convpsav(iq)*ps(iq)/grav
           endif
         endif
        enddo  ! iq loop
        convmax=0.
        nums=0
        do iq=1,imax     
         if(kt_sav(iq)-kb_sav(iq)==2.and.convpsav(iq)>0.)then
           nums=nums+1
           if(convpsav(iq)>convmax)then
             convmax=convpsav(iq)
             write(6,"('bcd ',3i5,2i3,2x,2f7.3,f6.3,2f5.2,2f7.4)") 
     .       ktau,iq,nums,kb_sav(iq),kt_sav(iq),
     .       dsk(kb_sav(iq)),
     .       fluxt(iq,kt_sav(iq)),convpsav(iq),
     .       qq(iq,kb_sav(iq))/qs(iq,kb_sav(iq)),
     .       qq(iq,kt_sav(iq))/qs(iq,kt_sav(iq)),
     &       rnrtcn(iq),rnrtcn(iq)*convpsav(iq)*ps(iq)/grav
           endif
         endif
        enddo  ! iq loop
        convmax=0.
        nums=0
        do iq=1,imax     
         if(kt_sav(iq)-kb_sav(iq)==3.and.convpsav(iq)>0.)then
           nums=nums+1
           if(convpsav(iq)>convmax)then
             convmax=convpsav(iq)
             write(6,"('bcde',3i5,2i3,2x,2f7.3,f6.3,2f5.2,2f7.4)") 
     .       ktau,iq,nums,kb_sav(iq),kt_sav(iq),
     .       dsk(kb_sav(iq)),
     .       fluxt(iq,kt_sav(iq)),convpsav(iq),
     .       qq(iq,kb_sav(iq))/qs(iq,kb_sav(iq)),
     .       qq(iq,kt_sav(iq))/qs(iq,kt_sav(iq)),
     &       rnrtcn(iq),rnrtcn(iq)*convpsav(iq)*ps(iq)/grav
           endif
         endif
        enddo  ! iq loop
      endif    ! (ntest==2.and.mydiag)    !######################
#endif
        
#ifndef GPU
      if(nmaxpr==1.and.mydiag)then
        iq=idjd
        write(6,*) 'Total_a delq (g/kg) & delt for this itn'
        write(6,"('delQ',3p12f7.3/(4x,12f7.3))")
     &            convpsav(iq)*delq(iq,:)
        write(6,"('delT',12f7.3/(4x,12f7.3))")
     &            convpsav(iq)*dels(iq,:)/cp
        write (6,"('hb/cpo',12f7.2/(9x,12f7.2))") 
     .            (splume(iq,:)+hl*qplume(iq,:))/cp
        write (6,"('qq_n  ',12f7.3/(9x,12f7.3))") 
     &             (qq(iq,:)+convpsav(iq)*delq(iq,:))*1000.
        write (6,"('new_qbas ',12f7.3/(9x,12f7.3))") 
     &     alfqarr(iq)*(qq(iq,:)+convpsav(iq)*delq(iq,:))*1000.
!       this s and h does not yet include updated phi (or does it?)     
        write (6,"('s/cp_n',12f7.2/(9x,12f7.2))") 
     &             (s(iq,:)+convpsav(iq)*dels(iq,:))/cp 
        write (6,"('h/cp_n',12f7.2/(9x,12f7.2))") 
     &            (s(iq,:)+hl*qq(iq,:)+convpsav(iq)*
     &             (dels(iq,:)+hl*delq(iq,:)))/cp       
        write (6,"('hb/cpn',12f7.2/(9x,12f7.2))") 
     &            (splume(iq,:)+hl*qplume(iq,:)+convpsav(iq)
     &            *(dels(iq,:)+alfqarr(iq)*hl*delq(iq,:)))/cp 
        write (6,"('hs/cpn',12f7.2/(9x,12f7.2))") 
     &             (hs(iq,:)+convpsav(iq)*dels(iq,:)+
     &              (1.+hlcp*dqsdt(iq,:)))/cp  
      endif
#endif

      if(methdetr<0)qxcess(:)=0.
      if(methdetr==0)qxcess(:)=detrain*rnrtcn(:)             ! e.g. .2* gives 20% detrainment
      if(sig_ct<0.)then  ! detrainn for shallow clouds; only Jack was using sig_ct ~ -.8
        do iq=1,imax
         if(sig(kt_sav(iq))>-sig_ct)then  
           qxcess(iq)=rnrtcn(iq)     ! full detrainment
         endif
        enddo  ! iq loop
      endif    !  (sig_ct<0.)    other sig_ct options removed June 2012

!     new detrainn calc producing qxcess to also cope with shallow conv 
!     - usually with methprec=4,5 or 6 too
      if(methdetr==1)then     ! like 3 but uses .55 instead of .6
        do iq=1,imax     
         detrainn=min( 1., detrain+(1.-detrain)*(   ! typical detrainn is .15
     &   (.55-min(sig(kb_sav(iq))-sig(kt_sav(iq)), .55)) /(.6-.14))**3 )
!        diff=0, .14, .2,.3, .4, .6, .7 gives detrainn=1, .71, .52, .29, .22, .18, .15   for .15  
        qxcess(iq)=detrainn*rnrtcn(iq) 
        enddo
      endif   ! (methdetr==1)     
      if(methdetr==2)then     ! like 3 but uses .57 instead of .6
        do iq=1,imax     
         detrainn=min( 1., detrain+(1.-detrain)*(   ! typical detrainn is .15
     &   (.57-min(sig(kb_sav(iq))-sig(kt_sav(iq)), .57)) /(.6-.14))**3 )
!        diff=0, .14, .2,.3, .4, .6, .7 gives detrainn=1, .84, .59, .32, .19, .18, .15   for .15  
        qxcess(iq)=detrainn*rnrtcn(iq) 
        enddo
      endif   ! (methdetr==2)     
      if(methdetr==3)then  ! usual (with .14 from 20/6/12). Used to calculate qxcess
        do iq=1,imax     
         detrainn=min( 1., detrain+(1.-detrain)*(   ! typical detrainn is .15
     &     (.6-min(sig(kb_sav(iq))-sig(kt_sav(iq)), .6)) /(.6-.14))**3 )
!        diff=0, .14, .2, .3, .4, .6, .7 gives detrainn=1, 1, .71, .39, .22, .16, .15   for .15  
!        With .55 instead of .6 would get
!        diff=0, .14, .2,.3, .4, .6, .7 gives detrainn=1, .71, .52, .29, .22, .18, .15   for .15  
!        With .57 instead of .6 would get
!        diff=0, .14, .2,.3, .4, .6, .7 gives detrainn=1, .84, .59, .32, .19, .18, .15   for .15  
!        With detrain=.05 instead of .15 would get
!        diff=0, .14, .2,.3, .4, .6, .7 gives detrainn=1, 1., .68, .31, .13, .05, .05   for .05
!        With detrain=.02 instead of .15 would get
!        diff=0, .14, .2,.3, .4, .6, .7 gives detrainn=1, 1., .66, .29, .10, .02, .02   for .02
!        With detrain=.01 instead of .15 would get
!        diff=0, .14, .2,.3, .4, .6, .7 gives detrainn=1, 1., .66, .28, .09, .01, .01   for .01
        qxcess(iq)=detrainn*rnrtcn(iq) 
        enddo
      endif   ! (methdetr==3)     
      if(methdetr==4)then   ! like 3 but uses .7 instead of .6 - was more usual
        do iq=1,imax     
         detrainn=min( 1., detrain+(1.-detrain)*(   ! typical detrainn is .15
     &     (.7-min(sig(kb_sav(iq))-sig(kt_sav(iq)), .7)) /(.6-.14))**3 )
!        diff=0, .14, .2, .3, .4, .5, .6, .7 gives detrainn=1, 1., 1., .71, .39, .22, .16, .15   for .15  
!        diff=0, .14, .2, .3, .4, .5, .6, .7 gives detrainn=1, 1., 1., .69, .35, .17, .11, .1     for .1  
         qxcess(iq)=detrainn*rnrtcn(iq) 
        enddo
      endif   ! (methdetr==4)     
      if(methdetr==5)then   ! like 4 with .7 but denom is .5
        do iq=1,imax     
         detrainn=min( 1., detrain+(1.-detrain)*(   ! typical detrainn is .15
     &     (.7-min(sig(kb_sav(iq))-sig(kt_sav(iq)), .7)) /.5)**3 )
!        diff=0, .14, .2, .3, .4, .5, .6, .7 gives detrainn=1, 1., 1., .71, .39, .22, .16, .15   for .15 X 
!        diff=0, .14, .2, .3, .4, .5, .6, .7 gives detrainn=1, 1., 1., .69, .35, .17, .11, .1     for .1  X
         qxcess(iq)=detrainn*rnrtcn(iq) 
        enddo
      endif   ! (methdetr==5)     
      if(methdetr==6)then   ! like 4 with .7 but denom is .45
        do iq=1,imax     
         detrainn=min( 1., detrain+(1.-detrain)*(   ! typical detrainn is .15
     &     (.7-min(sig(kb_sav(iq))-sig(kt_sav(iq)), .7)) /.45)**3 )
!        diff=0, .14, .2, .3, .4, .5, .6, .7 gives detrainn=1, 1., 1., .71, .39, .22, .16, .15   for .15 X 
!        diff=0, .14, .2, .3, .4, .5, .6, .7 gives detrainn=1, 1., 1., .69, .35, .17, .11, .1     for .1  X
         qxcess(iq)=detrainn*rnrtcn(iq) 
        enddo
      endif   ! (methdetr==6)     

      if(itn==1)then   ! ************************************************
        cape(:)=0.
        do k=2,kl-1
         do iq=1,imax
          if(k>kb_sav(iq).and.k<=kt_sav(iq))then
            kb=kb_sav(iq)
            cape(iq)=cape(iq)-(splume(iq,kb)+hl*qplume(iq,kb)-hs(iq,k))*
     &              rdry*dsig(k)/(cp*sig(k))            
          endif
         enddo
        enddo  ! k loop
       if(convtime>100.)then  ! land/sea effects, with tied_con irrelevant
           do iq=1,imax    
             if(land(iq))then  ! e.g. for 3060.60 
               convtim_deep(iq)=mcontlnd
             else
               convtim_deep(iq)=mcontsea
             endif
           enddo
       endif  !  (convtime>100.)
       if(convtime<-100.)then  ! newer uses for tied_con (when convtime < -100)
!                                              all with convt_frac treatment
!                       fg or omega<0
!                  2       
!                       sharp-pbl_depth.or.omega_sign
!                 1       
!                       gradual-pbl_depth.and.omega_sign
!                  0       
!                       sharp-pbl_depth.and.omega<0
!                 -1       
!                       magnitude of omega
!        N.B. tied_con now only used convtime < -100    
         if(tied_con>2.)then   ! now omega_900<0 or fg>tied_con test                                         
           do iq=1,imax    
             if(fg(iq)>tied_con.or.dpsldt(iq,k900)<0.)then  ! e.g. for -2099.60 t_c=40
               convtim_deep(iq)=mcontlnd
             else
               convtim_deep(iq)=mcontsea
             endif
           enddo
         elseif(tied_con>1.)then   !  e.g. 1.85       omega.or.pbl test
           do iq=1,imax    
             if(sig(kb_sav(iq))<tied_con-1..or.dpsldt(iq,k900)<0.)then ! e.g. -2040.60 t_c=1.85
               convtim_deep(iq)=mcontlnd
             else
               convtim_deep(iq)=mcontsea
             endif
           enddo
         elseif(tied_con<-.99)then ! for magnitude_of_omega test
           sumb=1.e-8*tied_con  ! a -ve test value
           omgtst(:)=max(min(dpsldt(:,k900)/sumb,1.),0.) ! used for lin interp
           convtim_deep(:)=mcontsea+omgtst(:)*(mcontlnd-mcontsea)
         elseif(tied_con>0.)then   !  e.g. tied_con = .85;  omega.and.pbl test
           do iq=1,imax    
             if(dpsldt(iq,k900)<0.)then 
               convtim_deep(iq)=mcontlnd*
     &          max(1.,(1.-tied_con)/max(.01,1-sig(kb_sav(iq))))
             else
               convtim_deep(iq)=mcontsea*
     &          max(1.,(1.-tied_con)/max(.01,1-sig(kb_sav(iq))))
             endif
           enddo
         elseif(tied_con<0.)then                      !  e.g. .tied_con = -.9
           do iq=1,imax    
             if(sig(kb_sav(iq))<-tied_con.and.dpsldt(iq,k900)<0.)then ! e.g. -2040.60 t_c=-.9
               convtim_deep(iq)=mcontlnd
             else
               convtim_deep(iq)=mcontsea
             endif
           enddo
         else   !   (tied_con==0.)
c           print *,'has tied_con=0'
           do iq=1,imax    
             if(omega(iq)<0.)then 
               convtim_deep(iq)=mcontlnd
             else
               convtim_deep(iq)=mcontsea
             endif
           enddo
         endif  ! (tied_con>1.)  .. else ..
       endif    ! (convtime<-100.)
#ifndef GPU
        if(nmaxpr==1.and.mydiag)then
         iq=idjd
         write(6,*) 'kb_sav,kt_sav',kb_sav(iq),kt_sav(iq)
         k=kt_sav(idjd)
         write(6,*) 'sig_kb_sav,sig_k',sig(kb_sav(iq)),sig(k)
         write(6,*) 'timeconvb,b',timeconv(idjd)
     &    ,max(0.,min(sig(kb_sav(idjd))-sig(k)-.2,.4)*mdelay/.4) 
        endif
#endif
        do iq=1,imax
         if(convpsav(iq)>0.)then ! used by mdelay options
           timeconv(iq)=timeconv(iq)+dt/60. ! now in minutes
         else       ! usual
           timeconv(iq)=0.
         endif
        enddo
      
         do iq=1,imax   
           sumb=min(1.,  (sig(kb_sav(iq))-sig(kt_sav(iq)))/convt_frac)
!         convtim_deep in mins between mcontlnd and mcontsea, linearly
           factr(iq)=min(1.,sumb*dt/(60*convtim_deep(iq)))  ! defined just for itn=1
!           if(iq==idjd)
!     &      write(6,*) 'fg,sumb,tim_deep,thick,dpsldt,factr,kb,kt',
!     &     fg(iq),sumb,convtim_deep(iq),sig(kb_sav(iq))-sig(kt_sav(iq)),
!     &    1.e8*dpsldt(iq,k900),factr(iq),kb_sav(iq),kt_sav(iq)
         enddo
#ifndef GPU
        if(nmaxpr==1.and.mydiag)then
        iq=idjd
        write(6,*)'timeconvc,convpsav,sig_b,sig_t,convt_frac,factr_a',
     &      timeconv(iq),convpsav(iq),sig(kb_sav(iq)),sig(kt_sav(iq)),
     &      convt_frac,factr(iq)
        endif
#endif

        if(tied_b>1.)then ! typical value is 26, used directly with factr
!         tied_b=26 gives factr [1, .964, .900, .794,.5,.216] for ds = [200, 100, 50, 25,8,2} km     
!         tied_b=10 gives factr [1, .917, .786, .611,.3,.100] for ds = [200, 100, 50, 25,8,2} km     
          do iq=1,imax
            summ=ds/(em(iq)*208498.)
            factr(iq)=factr(iq)*(1.+tied_b)*summ/(1.+tied_b*summ)
          enddo
       endif  ! (tied_b>1.)
      do iq=1,imax
         kbsav(iq)=kb_sav(iq)   ! just value at end of itn=1  Apr'15
         ktsav(iq)=kt_sav(iq)  
      enddo   ! iq loop
      endif    ! (itn==1)     ! ************************************************

#ifndef GPU
      if(nmaxpr==1.and.nevapcc.ne.0.and.mydiag)then
       write (6,
     & "('itn,kb_sd,kt_sd,kb,kt,delS,timeconv,entrainn,factr',
     &       5i3,5f7.3)")
     &   itn,kb_saved(idjd),kt_saved(idjd),
     &   kb_sav(idjd),kt_sav(idjd),
     &   sig(kb_sav(idjd))-sig(kt_sav(idjd)),
     &   timeconv(idjd),entrainn(idjd),factr(idjd) 
      endif
#endif

!      do iq=1,imax
!       if(ktsav(iq)==kl-1.and.convpsav(iq)>0.)then
!         kbsav(iq)=kb_sav(iq) 
!         ktsav(iq)=kt_sav(iq)  
!       endif  ! (ktsav(iq)==kl-1.and.convpsav(iq)>0.)
!      enddo   ! iq loop

      if(itn<iterconv.or.iterconv<0)then   ! iterconv<0 calls convjlm twice, and does convfact*
        convpsav(:)=convfact*convpsav(:) ! typically convfact=1.05
!       convpsav(:)=convfact*convpsav(:)*factr(:) ! dangerous tried on 7/3/14 
      endif                              ! (itn<iterconv.or.iterconv<0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      rnrtcn(:)=rnrtcn(:)-qxcess(:)

!     update qq, tt and precip
      rnrtc(:)=rnrtc(:)+convpsav(:)*rnrtcn(:)*conrev(:) ! g/m**2/s
      if(ncvcloud==0)then  ! usual
        do k=1,kl   
         do iq=1,imax
          qq(iq,k)=qq(iq,k)+convpsav(iq)*delq(iq,k)
          tt(iq,k)=tt(iq,k)+convpsav(iq)*dels(iq,k)/cp
         enddo    ! iq loop
        enddo     ! k loop
      else   ! phi correction done for ncvcloud.ne.0
        do k=1,kl   
         do iq=1,imax
          qq(iq,k)=qq(iq,k)+convpsav(iq)*delq(iq,k)
         enddo    ! iq loop
        enddo     ! k loop
!        fluxt here denotes delt, phi denotes del_phi
           do iq=1,imax
            fluxt(iq,1)=convpsav(iq)*dels(iq,1)/(cp+bet(1))
            phi(iq,1)=bet(1)*fluxt(iq,1)  ! phi correction done for ncvcloud.ne.0
            tt(iq,1)=tt(iq,1)+fluxt(iq,1)
           enddo    ! iq loop
          do k=2,kl
           do iq=1,imax
            fluxt(iq,k)=(convpsav(iq)*dels(iq,k)-phi(iq,k-1)
     &                  -betm(k)*fluxt(iq,k-1))/(cp+bet(k))
            phi(iq,k)=phi(iq,k-1)+betm(k)*fluxt(iq,k-1)  ! phi correction done for ncvcloud.ne.0
     &                +bet(k)*fluxt(iq,k)
            tt(iq,k)=tt(iq,k)+fluxt(iq,k)
           enddo    ! iq loop
        enddo     ! k loop
      endif       ! (ncvcloud==0) .. else ..
      
!!!!!!!!!!!!!!! "deep" detrainnment using detrainn !!!!! v3 !!!!!!    
!     N.B. convpsav has been updated here with convfact, but without factr term  
      if(methdetr<0)then
       do k=2,kl-1 
        do iq=1,imax  
         qliqw(iq,k)=qliqw(iq,k)+convpsav(iq)*delqliqw(iq,k)    ! after applying convpsav, accumulating for itns
        enddo  ! iq loop
       enddo   ! k loop
      else
       do k=2,kl-1 
        do iq=1,imax  
         if(k>kb_sav(iq).and.k<=kt_sav(iq))then
           deltaq=convpsav(iq)*qxcess(iq)*
     &                   detrarr(k,kb_sav(iq),kt_sav(iq))/dsk(k)
           qliqw(iq,k)=qliqw(iq,k)+deltaq    ! after applying convpsav, accumulating for itns
         endif  ! (k>kb_sav(iq).and.k<=kt_sav(iq))
        enddo  ! iq loop
       enddo   ! k loop
      endif  ! (methdetr<0)  ..  else .. 

#ifndef GPU
      if(ntest>0.and.mydiag)then
        iq=idjd
        print *,'liqw',convpsav(iq),qxcess(iq),kb_sav(iq),kt_sav(iq),
     &                   detrarr(10,kb_sav(iq),kt_sav(iq)),dsk(10)
!       N.B. convpsav(iq) is already mult by dt jlm: mass flux is convpsav/dt
        write(6,*) "after convection: ktau,itn,kbsav,ktsav ",

     .                 ktau,itn,kb_sav(iq),kt_sav(iq)
        write (6,"('qgc ',12f7.3/(8x,12f7.3))") 
     .             (1000.*qq(iq,k),k=1,kl)
        write (6,"('ttc ',12f7.2/(8x,12f7.2))") 
     .             (tt(iq,k),k=1,kl)
        write(6,*) 'rnrtc,qxcess ',rnrtc(iq),qxcess(iq)
        write(6,*) 'rnrtcn,convpsav ',rnrtcn(iq),convpsav(iq)
        write(6,*) 'ktsav,qplume,qs_ktsav,qq_ktsav ',kt_sav(iq),
     .         qplume(iq,kb_sav(iq)),qs(iq,kt_sav(iq)),qq(iq,kt_sav(iq))
        iq=idjd
        delq_av=0.
        delt_av=0.
        heatlev=0.
        do k=1,kl
         delq_av=delq_av+dsk(k)*convpsav(iq)*delq(iq,k)
         delt_av=delt_av+dsk(k)*convpsav(iq)*dels(iq,k)/cp
         heatlev=heatlev+sig(k)*dsk(k)*convpsav(iq)*dels(iq,k)/cp
        enddo
        write(6,*) 'delq_av,delt_exp,rnd_exp ',
     .         delq_av,-delq_av*hl/cp,-delq_av*conrev(iq)
        if(abs(delt_av)>1.e-20)write(6,*) 
     &        'ktau,itn,kbsav,ktsav,delt_av,heatlev',
     .        ktau,itn,kb_sav(iq),kt_sav(iq),delt_av,heatlev/delt_av
      endif   ! (ntest>0)
#endif
      
!     update u & v using actual delu and delv (i.e. divided by dsk)
      if(nuvconv.ne.0.or.nuv>0)then
#ifndef GPU
        if(ntest>0.and.mydiag)then
          write(6,*) 'u,v before convection'
          write (6,"('u  ',12f7.2/(3x,12f7.2))") u(idjd,:)
          write (6,"('v  ',12f7.2/(3x,12f7.2))") v(idjd,:)
        endif  ! (ntest>0.and.mydiag)
#endif
        do k=1,kl-2   
         do iq=1,imax
          u(iq,k)=u(iq,k)+facuv*factr(iq)*convpsav(iq)*delu(iq,k)/dsk(k)
          v(iq,k)=v(iq,k)+facuv*factr(iq)*convpsav(iq)*delv(iq,k)/dsk(k)
         enddo  ! iq loop
        enddo   ! k loop
#ifndef GPU
        if(ntest>0.and.mydiag)then
          write(6,*) 'u,v after convection'
          write (6,"('u  ',12f7.2/(3x,12f7.2))") u(idjd,:)
          write (6,"('v  ',12f7.2/(3x,12f7.2))") v(idjd,:)
        endif
#endif
      endif     ! (nuvconv.ne.0)

!     section for convective transport of trace gases (jlm 22/2/01)
      if(ngas>0)then
        do ntr=1,ngas
         do k=1,kl-2
          do iq=1,imax
           s(iq,k)=tr(iq,k,ntr)
          enddo    ! iq loop
         enddo     ! k loop
         do iq=1,imax
          if(kt_sav(iq)<kl-1)then
            kb=kb_sav(iq)
            kt=kt_sav(iq)
            veldt=factr(iq)*convpsav(iq)*(1.-fldow(iq)) ! simple treatment
            fluxup=veldt*s(iq,kb)
            fluxup=min(fluxup,tr(iq,kb,ntr)*dsk(kb))
!           remove gas from cloud base layer
            tr(iq,kb,ntr)=tr(iq,kb,ntr)-fluxup/dsk(kb)
!           put flux of gas into top convective layer
            tr(iq,kt,ntr)=tr(iq,kt,ntr)+fluxup/dsk(kt)
            do k=kb+1,kt
             fluxup=s(iq,k)*veldt
             fluxup=min(fluxup,tr(iq,k,ntr)*dsk(k))
             tr(iq,k,ntr)=tr(iq,k,ntr)-fluxup/dsk(k)
             tr(iq,k-1,ntr)=tr(iq,k-1,ntr)+fluxup/dsk(k-1)
            enddo
          endif
         enddo   ! iq loop
        enddo    ! ntr loop    
      endif      ! (ngas>0)

      ! Convective transport of aerosols - MJT
      if ( abs(iaero)==2 ) then
        do ntr = 1,naero
          xtgscav(:,:) = 0.
          s(:,1:kl-2) = xtg(1:imax,1:kl-2,ntr)
!$acc loop private(ttsto,qqsto,qlsto,qqold,qlold,rho,xtgtmp)
          do iq = 1,imax
           if ( kt_sav(iq)<kl-1 ) then
             kb = kb_sav(iq)
             kt = kt_sav(iq)
             knet = kt - kb + 1

             ! convective scavenging of aerosols
             do k = kb,kt
               ttsto(k) = t(iq,k) + factr(iq)*(tt(iq,k)-t(iq,k))
               qqsto(k) = qg(iq,k) + factr(iq)*(qq(iq,k)-qg(iq,k))
               qqold(k) = qg(iq,k) + factr(iq)*(qqsav(iq,k)-qg(iq,k))
               qlsto(k) = factr(iq)*qliqw(iq,k)
               qlold(k) = factr(iq)*qliqwsav(iq,k)
               rho(k) = ps(iq)*sig(k)/(rdry*ttsto(k))
               ! convscav expects change in liquid cloud water, so we assume
               ! change in qg is added to qlg before precipating out
               qqold(k) = qqold(k) - qqsto(k) ! created liquid water
               qqsto(k) = qlsto(k) - qlold(k) ! remaining liquid water
               xtgtmp(k) = xtg(iq,k,3)
             end do
             ! calculate fraction of aerosols that will be scavenged
             call convscav(fscav(kb:kt),qqsto(kb:kt),qqold(kb:kt),
     &                  ttsto(kb:kt),xtgtmp(kb:kt),rho(kb:kt),ntr,
     &                  knet)
               
             veldt = factr(iq)*convpsav(iq)*(1.-fldow(iq)) ! simple treatment
             fluxup = veldt*s(iq,kb)
             ! remove aerosol from lower layer
             fluxup = min( fluxup, xtg(iq,kb,ntr)*dsk(kb) )
             xtg(iq,kb,ntr) = xtg(iq,kb,ntr) - fluxup/dsk(kb)
             ! store aerosol concentration that was scavenged
             xtgscav(iq,kt) = fluxup*fscav(kt)/dsk(kt)
             ! put flux of aerosol into upper layer after scavenging
             xtg(iq,kt,ntr) = xtg(iq,kt,ntr) + fluxup/dsk(kt)
     &           - xtgscav(iq,kt)
             ! continue the calculation along the column
             do k = kb+1,kt
              ! remove aerosol from upper layer  
              fluxup = veldt*s(iq,k)
              fluxup = min( fluxup, xtg(iq,k,ntr)*dsk(k) )
              xtg(iq,k,ntr) = xtg(iq,k,ntr) - fluxup/dsk(k)
              ! store aerosol concentration that was scavenged
              xtgscav(iq,k-1) = fluxup*fscav(k-1)/dsk(k-1)
              ! put flux of aerosol into lower layer after scavenging
              xtg(iq,k-1,ntr) = xtg(iq,k-1,ntr) + 
     &           fluxup/dsk(k-1) - xtgscav(iq,k-1)
             enddo ! k loop
           endif
          enddo    ! iq loop
          ! store wet deposition from scavenging
          if ( ntr==itracso2 ) then
            do k=1,kl
              so2wd=so2wd+xtgscav(:,k)*ps(1:imax)*dsk(k)
     &          /(grav*dt)
            end do
          elseif ( ntr==itracso2+1 ) then
            do k=1,kl
              so4wd=so4wd+xtgscav(:,k)*ps(1:imax)*dsk(k)
     &          /(grav*dt)
            end do
          elseif (ntr==itracbc.or.ntr==itracbc+1) then
            do k=1,kl
              bcwd=bcwd+xtgscav(:,k)*ps(1:imax)*dsk(k)
     &          /(grav*dt)
            end do
          elseif (ntr==itracoc.or.ntr==itracoc+1) then
            do k=1,kl
              ocwd=ocwd+xtgscav(:,k)*ps(1:imax)*dsk(k)
     &          /(grav*dt)
            end do
          elseif (ntr>=itracdu.and.ntr<=itracdu+ndust-1) then
            do k=1,kl
              dustwd(:,ntr-itracdu+1)=dustwd(:,ntr-itracdu+1)
     &          +xtgscav(:,k)*ps(1:imax)*dsk(k)
     &          /(grav*dt)
            end do
          elseif (ntr>=itracsa.and.ntr<=itracsa+nsalt-1) then
            do k=1,kl
              saltwd=saltwd+xtgscav(:,k)*ps(1:imax)*dsk(k)
     &          /(grav*dt)
            end do              
          end if
        end do     ! nt loop
      end if   ! (abs(iaero)==2) 
      
#ifndef GPU
      if(ntest>0.and.mydiag)then
        iq=idjd
        write (6,"('uuc ',12f6.1/(4x,12f6.1))") (u(iq,k),k=1,kl)
        write (6,"('vvc ',12f6.1/(4x,12f6.1))") (v(iq,k),k=1,kl)
      endif  ! (ntest>0.and.mydiag)
#endif

#ifndef GPU
      if(nmaxpr==1.and.mydiag)then
       iq=idjd
       if(ktau==1.and.itn==1)write(6, "(15x,
     &  'ktau itn kb_sav kmin kdown kt_sav cfrac+ entxsav detxsav ',
     &  ' fluxv  fluxvt-1 factr convpsav3 rnrtcn3 cape')")
       if(kb_sav(iq)<kl-1)write (6,"('ktau ... cape ',
     &   i5,i4,']',2i5,2i6,f9.2,5f8.2,3p2f8.3,1pf9.1)")  ktau,itn,
     &   kb_sav(iq),kmin(iq),kdown(iq),kt_sav(iq),
     &   cfrac(iq,kb_sav(iq)+1),entrsav(iq,kdown(iq)),
     &   detxsav(iq,kdown(iq)),fluxv(iq,kdown(iq)),
     &   fluxv(iq,kt_sav(iq)-1),factr(iq),convpsav(iq),
     &   rnrtcn(iq),cape(iq)
       dtsol=.01*sgsave(iq)/(1.+.25*(u(iq,1)**2+v(iq,1)**2))   ! solar heating   ! sea
       write (6,"('pblh,fldow,tdown,qdown,fluxq3',
     &             f8.2,f5.2,2f7.2,3p2f8.3,1pf8.2)")
     &     pblh(iq),fldow(iq),tdown(iq),qdown(iq),fluxq(iq)
       write(6,"('ktau,kkbb,wetfac,dtsol,alfqarr,fg,omega8
     &,omgtst',i5,i3,3f6.3,f8.3,8pf9.3,1pf9.3)") ktau,kkbb(idjd),
     &      wetfac(idjd),dtsol,alfqarr(idjd),fg(idjd),
     &     omega(idjd),omgtst(idjd)
       write(6,"('fluxt3',3p13f10.3)") fluxt(iq,kb_sav(iq)+1:kt_sav(iq))
      endif
#endif

      if(itn==1)then
        kb_saved(:)=kb_sav(:)
        kt_saved(:)=kt_sav(:)
      else
        do iq=1,imax
         if(kt_sav(iq)<kl-1)then
!          itns 2 & 3 only update saved values if convection occurred         
            kb_saved(iq)=kb_sav(iq)
            kt_saved(iq)=kt_sav(iq)
         endif
        enddo
      endif
      do k=1,kl-1
         do iq=1,imax
c         if(fluxv(iq,k)>1.)fluxtot(iq,k)=fluxtot(iq,k)+ 
          fluxtot(iq,k)=fluxtot(iq,k)+   
     &          fluxv(iq,k)*factr(iq)*convpsav(iq)/dt  ! needs div by dt (check)
        enddo
      enddo

      enddo     ! itn=1,abs(iterconv)
!------------------------------------------------ end of iterations ____#################################

#ifndef GPU
      if(ntest>0.and.mydiag)
     &  write (6,"('C k,fluxtot6',i3,f7.2)")
     &                               (k,1.e6*fluxtot(idjd,k),k=1,kl)
      if(nmaxpr==1.and.mydiag)then
        write(6,*) 'convtime,factr,kb_sav,kt_sav',convtime,
     &          factr(idjd),kb_sav(idjd),kt_sav(idjd)
      endif
#endif
      rnrtc(:)=factr(:)*rnrtc(:)    ! N.B. factr applied after itn loop 
      do k=1,kl
      qq(1:imax,k)=qg(1:imax,k)+factr(:)*(qq(1:imax,k)-qg(1:imax,k))
      qliqw(1:imax,k)=factr(:)*qliqw(1:imax,k)      
      tt(1:imax,k)= t(1:imax,k)+factr(:)*(tt(1:imax,k)- t(1:imax,k))
      enddo

!     update qq, tt for evap of qliqw (qliqw arose from moistening detrainment)
      if(ldr.ne.0)then
!       Leon's stuff here, e.g.
        if(abs(rhmois)<1.e-20)then  ! Nov 2012
          do k=1,kl            
!           this is older simpler option, allowing ldr scheme to assign qfg without time complications          
            qlg(1:imax,k)=qlg(1:imax,k)+qliqw(1:imax,k)
          end do
        else
          do k=1,kl              
           do iq=1,imax
            if(tt(iq,k)<253.16)then   ! i.e. -20C
              qfg(iq,k)=qfg(iq,k)+qliqw(iq,k)
              tt(iq,k)=tt(iq,k)+3.35e5*qliqw(iq,k)/cp   ! fusion heating
            else
              qlg(iq,k)=qlg(iq,k)+qliqw(iq,k)
            endif
           enddo
          enddo  ! k loop           
        endif !  (rhmois==0.)
      else      ! for ldr=0
        qq(1:imax,:)=qq(1:imax,:)+qliqw(1:imax,:)         
        tt(1:imax,:)=tt(1:imax,:)-hl*qliqw(1:imax,:)/cp   ! evaporate it
        qliqw(1:imax,:)=0.   ! just for final diags
      endif  ! (ldr.ne.0)
!_______________end of convective calculations__________________
     

#ifndef GPU
      if(nmaxpr==1.and.mydiag)then
        iq=idjd
        write(6,*) 'Total delq (g/kg) & delt after all itns'
        write(6,"('delQ_t',3p12f7.3/(6x,12f7.3))")
     .   (qq(iq,k)+qliqw(iq,k)-qg(iq,k),k=1,kl)
        write(6,"('delT_t',12f7.3/(6x,12f7.3))")
     .   (tt(iq,k)-t(iq,k),k=1,kl)
      endif
#endif
      qg(1:imax,:)=qq(1:imax,:)                   
      condc(1:imax)=.001*dt*rnrtc(1:imax)      ! convective precip for this timestep
      precc(1:imax)=precc(1:imax)+condc(1:imax)       
!    N.B. condx used to update rnd03,...,rnd24. LS added in leoncld.f       
      condx(1:imax)=condc(1:imax)+.001*dt*rnrt(1:imax) ! total precip for this timestep
      conds(:)=0.   ! (:) here for easy use in VCAM
      condg(:)=0.   ! (:) here for easy use in VCAM
      precip(1:imax)=precip(1:imax)+condx(1:imax)
      t(1:imax,:)=tt(1:imax,:)             

#ifndef GPU
      if(ntest>0.or.(ktau<=2.and.nmaxpr==1))then
       if(mydiag)then
        iq=idjd
        write(6,*) 'at end of convjlm: kdown,rnrt,rnrtc',
     &                                     kdown(iq),rnrt(iq),rnrtc(iq)
        phi(iq,1)=bet(1)*tt(iq,1)
        do k=2,kl
         phi(iq,k)=phi(iq,k-1)+bet(k)*tt(iq,k)+betm(k)*tt(iq,k-1)
        enddo      ! k  loop
        do k=1,kl   
         es(iq,k)=establ(tt(iq,k))  ! in diag loop using updated tt
         pk=ps(iq)*sig(k)
         qs(iq,k)=.622*es(iq,k)/max(pk-es(iq,k),0.1)  
         s(iq,k)=cp*tt(iq,k)+phi(iq,k)  ! dry static energy
         hs(iq,k)=s(iq,k)+hl*qs(iq,k)   ! saturated moist static energy
        enddo   ! k loop
        write (6,"('rhx   ',12f7.2/(5x,12f7.2))") 
     .             (100.*qq(iq,k)/qs(iq,k),k=1,kl)
        write (6,"('qsx   ',12f7.3/(5x,12f7.3))") 
     .             (1000.*qs(iq,k),k=1,kl)
        write (6,"('qqx   ',12f7.3/(5x,12f7.3))") 
     .             (1000.*qq(iq,k),k=1,kl)
        write (6,"('ttx   ',12f7.2/(5x,12f7.2))") 
     .             (tt(iq,k),k=1,kl)
        write (6,"('s/cpx ',12f7.2/(5x,12f7.2))") (s(iq,k)/cp,k=1,kl)
        write (6,"('h/cpx ',12f7.2/(5x,12f7.2))") 
     .             (s(iq,k)/cp+hlcp*qq(iq,k),k=1,kl)
        write (6,"('hb/cpx',12f7.2/(5x,12f7.2))") 
     .             (s(iq,k)/cp+hlcp*qplume(iq,k),k=1,kl)
        write (6,"('hs/cpx',12f7.2/(5x,12f7.2))") 
     .             (hs(iq,k)/cp,k=1,kl)
        write (6,"('k  ',12i7/(3x,12i7))") (k,k=1,kl)
        write(6,*) 'following are h,q,t changes during timestep'
        write (6,"('delh ',12f7.3/(5x,12f7.3))") 
     .             (s(iq,k)/cp+hlcp*qq(iq,k)-h0(k),k=1,kl)
        write (6,"('delq ',3p12f7.3/(5x,12f7.3))") 
     .             (qq(iq,k)-q0(k),k=1,kl)
        write (6,"('delt ',12f7.3/(5x,12f7.3))") 
     .             (tt(iq,k)-t0(k),k=1,kl)
        write(6,"('fluxq,convpsav',5f8.5)") ! printed 1st 2 steps
     .           fluxq(iq),convpsav(iq)
        write(6,"('fluxt',9f8.5/(5x,9f8.5))")
     .            (fluxt(iq,k),k=1,kt_sav(iq))
        pwater=0.   ! in mm     
        do k=1,kl
         pwater=pwater-dsig(k)*qg(iq,k)*ps(iq)/grav
        enddo
        write(6,*) 'pwater0,pwater+condx,pwater ',
     .           pwater0,pwater+condx(iq),pwater
        write(6,*) 'D condx ',condx(iq)
        write(6,*) 'precc,precip ',
     .           precc(iq),precip(iq)
       endif   ! (mydiag) needed here for maxmin
       call maxmin(rnrtc,'rc',ktau,1.,1)
      endif
#endif
      
#ifndef GPU
!      if(ntest==-1.and.nproc==1.and.ktau==10.and.ntiles==1)then
!        do k=1,kl
!         do iq=il*il+1,2*il*il
!          den1=qfg(iq,k)+qlg(iq,k)
!          if(land(iq).and.den1>0.)then
!            write(26,'(4f9.4,2i6)')
!     &       t(iq,k),qfg(iq,k)/den1,1000.*qfg(iq,k),1000.*qlg(iq,k),iq,k
!          endif
!         enddo
!        enddo     
!      endif  ! (ntest==-1.and.nproc==1)  
#endif
      return
      end subroutine convjlm_work

      end module convjlm_m
