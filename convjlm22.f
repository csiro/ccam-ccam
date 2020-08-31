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
      module convjlm22_m
      
      implicit none
      
      private
      public convjlm22, convjlm22_init
      
      integer, save :: k500,k600,k700,k900,k950,k970,k980,klon2,komega  
!$acc declare create(k600,k700,k900,klon2,komega)
      integer, save :: mcontlnd,mcontsea           
!$acc declare create(mcontlnd,mcontsea)
      real,save :: convt_frac,tied_a,tied_b
!$acc declare create(convt_frac,tied_b)
      real, dimension(:), allocatable, save :: timeconv, alfin, aug
      real, dimension(:,:), allocatable, save :: downex,upin,upin4
!$acc declare create(upin,downex)
      integer, dimension(:), allocatable, save :: kb_saved
      integer, dimension(:), allocatable, save :: kt_saved

      contains 
      
      subroutine convjlm22_init
      
      use cc_mpi, only : myid, ccmpi_abort, mydiag
      use const_phys
      use map_m
      use newmpar_m, only : ifull, kl
      use parm_m
      use sigs_m
      use soil_m

      implicit none

      include 'kuocom.h'   ! kbsav,ktsav,convfact,convpsav,ndavconv
      
      integer iq,k,ntest,kb
      real summ,sumb
      parameter (ntest=0)      ! 1 or 2 to turn on; -1 for ldr writes
      integer kpos(1)

      kpos=minloc(abs(sig-.98)) ! finds k value closest to sig=.98  level 2 for L18 & L27
      k980=kpos(1)
      kpos=minloc(abs(sig-.97)) ! finds k value closest to sig=.97  level 3 for L18 & L27
      k970=kpos(1)
      kpos=minloc(abs(sig-.95)) ! finds k value closest to sig=.95
      k950=kpos(1)
      kpos=minloc(abs(sig-.9)) ! finds k value closest to sig=.9
      k900=kpos(1)
      kpos=minloc(abs(sig-.7)) ! finds k value closest to sig=.7
      k700=kpos(1)
      kpos=minloc(abs(sig-.6)) ! finds k value closest to sig=.6
      k600=kpos(1)
!$acc update device(k600,k700,k900)
      komega=1  ! JLM
      if(dsig4>.4)then ! JLM
        kpos=minloc(abs(sig-dsig4)) ! finds k value closest to dsig4    JLM
        komega=kpos(1) ! JLM
      endif ! JLM
      k500=1
      do while(sig(k500)>0.5)
        k500=k500+1
      enddo
      k500=k500-1    ! level just below .5
!$acc update device(komega)
      if (myid==0) then
        write(6,*) 'k980 ',k980,sig(k980)
        write(6,*) 'k970 ',k970,sig(k970)
        write(6,*) 'k900 ',k900,sig(k900)
        write(6,*) 'k700 ',k700,sig(k700)
        write(6,*) 'k600 ',k600,sig(k600)
        write(6,*) 'k500 ',k500,sig(k500)
        write(6,*) 'komega',komega,sig(komega)  ! JLM
      end if
      allocate(timeconv(ifull))  ! init for ktau=1 (allocate needed for mdelay>0)
      allocate(alfin(ifull))     ! init for ktau=1
      allocate(aug(ifull))       ! init for ktau=1
      allocate(kb_saved(ifull))  ! init for ktau=1 
      allocate(kt_saved(ifull))  ! init for ktau=1 
      allocate(upin(kl,kl))
      allocate(upin4(kl,kl))
      allocate(downex(kl,kl))
      if ( methdetr>=0 .or. entrain>0 ) then
        write(6,*) "not permitted that methdetr>=0.or.entrain>0."
        call ccmpi_abort(-1)
      endif

      !----------------------------------------------------------------
      if (mbase<0.or.mbase>5.or.nbase.ne.3.or.alflnd<0.or.alfsea<0)then
        write(6,*) "ERROR: negative alflnd and alfsea or"
        write(6,*) "unsupported mbase, nbase in convjlm22"
        call ccmpi_abort(-1)
      endif
      kb_saved(:)=kl-1
      kt_saved(:)=kl-1
      timeconv(:)=0.
      do k=1,kl
        if(sig(k)> .5)klon2=k
      enddo
!$acc update device(klon2)
      if(myid==0)write(6,*) 'klon2,k500',klon2,k500
!     methdetr (for >0) gives the fraction of precip detrained

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
!       this one gives linearly increasing mass fluxes at each half level      
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
      if ( kscsea<=-2 ) then   ! trying to be backward compatible for downex
!       gives relatively larger downdraft fluxes near surface      
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
!$acc update device(downex,upin)
      if (myid==0) then
        do kb=1,kl-1
          write (6,"('upin1-kb ',17f7.3)") (upin(k,kb),k=1,kb)
          if(ntest==1) then
            write(6,"('downex',17f7.3)") (downex(k,kb),k=1,kb)
          end if  
        enddo
        do kb=1,kl-1
          write (6,"('downex',17f7.3)") (downex(k,kb),k=1,kb)
        enddo
      endif

      if(nint(convtime)==-23)convtime=3030.6
      if(nint(convtime)==-24)convtime=4040.6
      if(convtime>100.)then   ! new general style  May 2014
!       1836.45 is old 36;  4040.6 is old -24; 2020.001 is old .33
        mcontlnd=int(.01*convtime)         ! in minutes
        mcontsea=int(convtime-real(100*mcontlnd)) ! in minutes
        convt_frac=convtime-100*mcontlnd-mcontsea  ! changeover sigma value of cloud thickness
        convt_frac=max(convt_frac,1.e-7)  ! allows for zero entry
      elseif(convtime<-100.)then
        mcontlnd=-int(.01*convtime)         ! in minutes
        mcontsea=int(-convtime-real(100*mcontlnd)) ! in minutes
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
        aug(:)=tied_rh*cp  ! 16/05
      elsewhere
        alfin(:)=alfsea
        aug(:)=0.
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
     &              ds,tied_over,tied_a,tied_b
      end if
      if(tied_a>1.)then  ! e.g. 20 or 26  
!       alfin may be reduced over land and sea for finer resolution than 200 km grid            
        do iq=1,ifull
          summ=ds/(em(iq)*208498.)
          alfin(iq)=1.+(alfin(iq)-1.) *
     &      (1.+tied_a)*summ/(1.+tied_a*summ)
          !if(iq<200)print *,'iq,alfin',iq,alfin(iq)
!         tied_over=26 gives factor [1, .964, .900, .794, .529] for ds = [200, 100, 50, 25, 8} km     
!         tied_over=10 gives factor [1, .917, .786, .611, .306] for ds = [200, 100, 50, 25, 8} km     
        enddo
      endif  ! (tied_a>1.)
      !--------------------------------------------------------
         
      end subroutine convjlm22_init
      
      
      subroutine convjlm22
      
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
      real, dimension(imax,kl,naero) :: lxtg
      real, dimension(imax,kl,ntrac) :: ltr
      real, dimension(imax,kl)       :: ldpsldt, lt, lqg
      real, dimension(imax,kl)       :: lfluxtot
      real, dimension(imax,kl)       :: lqlg, lqfg, lcfrac
      real, dimension(imax,kl)       :: lu, lv
      real, dimension(imax,ndust)    :: ldustwd
      real, dimension(imax)          :: lso2wd, lso4wd, lbcwd
      real, dimension(imax)          :: locwd, lsaltwd
      logical :: mydiag_t

!$omp  do schedule(static) private(is,ie),
!$omp& private(ldpsldt,lt,lqg,lqlg,lqfg,lfluxtot,lcfrac),
!$omp& private(lu,lv),
!$omp& private(lxtg,lso2wd,lso4wd,lbcwd,locwd,ldustwd,lsaltwd),
!$omp& private(ltr,idjd_t,mydiag_t)
!$acc parallel loop copy(t,qg,qlg,qfg,u,v,xtg,dustwd,so2wd,so4wd,bcwd),
!$acc&  copy(ocwd,saltwd,tr,precc,precip,timeconv,kbsav,ktsav,cfrac),
!$acc&  copy(cape,condc,condx,conds,condg,kt_saved,kb_saved,aug),
!$acc&  copyin(dpsldt,alfin,ps,pblh,fg,wetfac,land),
!$acc&  copyin(em,sgsave),
!$acc&  copyout(convpsav)
!$acc&  private(ldpsldt,lt,lqg,lqlg,lqfg,lcfrac,lu,lv),
!$acc&  private(lxtg,ldustwd,lso2wd,lso4wd,lbcwd,locwd,lsaltwd),
!$acc&  private(ltr,lfluxtot)
      do tile = 1,ntiles
        is = (tile-1)*imax + 1
        ie = tile*imax
        
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

        ! jlm convective scheme
        call convjlm22_work(alfin(is:ie),ldpsldt,lt,lqg,
     &       ps(is:ie),lfluxtot,convpsav(is:ie),cape(is:ie),
     &       lxtg,lso2wd,lso4wd,lbcwd,locwd,ldustwd,lsaltwd,
     &       lqlg,condc(is:ie),precc(is:ie),condx(is:ie),conds(is:ie),
     &       condg(is:ie),precip(is:ie),
     &       pblh(is:ie),fg(is:ie),wetfac(is:ie),land(is:ie),lu,lv,
     &       timeconv(is:ie),em(is:ie),
     &       kbsav(is:ie),ktsav(is:ie),ltr,lqfg,lcfrac,sgsave(is:ie),
     &       kt_saved(is:ie),kb_saved(is:ie),aug(is:ie),
     &       idjd_t,mydiag_t,entrain,detrain,mbase,iterconv,
     &       nuvconv,alfsea,methdetr,methprec,fldown,alflnd,rhcv,
     &       convtime,nkuo,rhsat,nevapls,
     &       tied_con,mdelay,convfact,ncvcloud,ldr,rhmois,imax,kl)

        t(is:ie,:)       = lt
        qg(is:ie,:)      = lqg
        qlg(is:ie,:)     = lqlg
        qfg(is:ie,:)     = lqfg
        fluxtot(is:ie,:) = lfluxtot
        u(is:ie,:)       = lu
        v(is:ie,:)       = lv
        cfrac(is:ie,:)   = lcfrac
        if ( abs(iaero)>=2 ) then
          xtg(is:ie,:,:)  = lxtg
          dustwd(is:ie,:) = ldustwd
          so2wd(is:ie)    = lso2wd
          so4wd(is:ie)    = lso4wd
          bcwd(is:ie)     = lbcwd
          ocwd(is:ie)     = locwd
          saltwd(is:ie)   = lsaltwd
        end if
        if ( ngas>0 ) then
          tr(is:ie,:,:) = ltr
        end if
        
      end do
!$acc end parallel
!$omp end do nowait

      return
      end subroutine convjlm22     ! jlm convective scheme
       
      
      subroutine convjlm22_work(alfin,dpsldt,t,qg,ps,
     &       fluxtot,convpsav,cape,xtg,so2wd,so4wd,bcwd,ocwd,
     &       dustwd,saltwd,qlg,condc,precc,condx,conds,condg,precip,
     &       pblh,fg,wetfac,land,u,v,timeconv,em,
     &       kbsav,ktsav,tr,qfg,cfrac,sgsave,kt_saved,kb_saved,aug,
     &       idjd,mydiag,entrain,detrain,mbase,iterconv,
     &       nuvconv,alfsea,methdetr,methprec,fldown,alflnd,rhcv,
     &       convtime,nkuo,rhsat,nevapls,
     &       tied_con,mdelay,convfact,ncvcloud,ldr,rhmois,imax,kl)
!$acc routine vector

      !jlm convective scheme - latest and cleaned up
!     unused switches: nevapcc, rhsat, shaltime 
!     unused switches if ksc=0:  kscmom, 
!     unused switches if nkuo=21,22,23:  rhsat
!     unused in convjlm22: dsig2, sigkscb, sigksct, tied_rh
!     has +ve fldownn depending on delta sigma; (-ve fldown descends from sig=.6))   
      use aerosolldr, only : itracso2,itracbc,itracoc,itracdu,ndust,
     &                       naero,convscav,itracsa,nsalt
      use cc_mpi, only : ccmpi_abort
!      use cc_omp, only : imax, ntiles
      use const_phys
      use diag_m, only : maxmin
      use estab      
!      use newmpar_m
      use parm_m, only : ktau,dt,nmaxpr,diag,ds,iaero
      use parmdyn_m
      use sigs_m
      use tracers_m, only : ngas,ntrac
      
      implicit none
      
      integer itn,iq,k,kt,ntest,ntr,nums,nuv,kb
      integer idjd
      real convmax,delq_av,delt_av,den1,den2,den3,dprec
      real facuv,fldownn,fluxup,hbase,heatlev,pwater,pwater0,qsk
      real rnrt_k,summ,totprec,veldt,pk,dz,sumb,bbb,ccc
      real dtsol
      logical mydiag
      parameter (ntest=0)      ! 1 or 2 to turn on; -1 for ldr writes
!     convjlm22 requires methdetr=-1,-2 or -3; entrain -ve; nbase=-10  or +ve
!     parameter (iterconv=3)  ! to kuocom.h
!     parameter (fldown=-.3)  ! to kuocom.h
!     parameter (detrain=.1)  ! to kuocom.h
!     parameter (alflnd=1.20) 
!     parameter (alfsea=1.10) 
!                                3:for RH-enhanced base
!     parameter (methdetr=-1)  ! fraction of precip detrained as fn of cloud depth; new scheme for -ve
!     parameter (nuvconv=0)    ! usually 0, >0 or <0 to turn on momentum mixing
!     parameter (nuv=0)        ! usually 0, >0 to turn on new momentum mixing (base layers too)
!     nevapls:  turn off/on ls evap - through parm.h; 0 off, 5 newer UK

      integer, intent(in) :: imax, kl
      integer knet
      real, dimension(imax,kl,naero), intent(inout)    :: xtg
      real, dimension(imax,kl,ntrac), intent(inout)    :: tr
      real, dimension(imax,kl), intent(in)         :: dpsldt
      real, dimension(imax,kl), intent(inout)      :: cfrac
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
      real, dimension(imax), intent(in)                :: em
      real, dimension(imax), intent(in)                :: sgsave
      real, dimension(imax), intent(inout)             :: cape
      real, dimension(imax), intent(inout)             :: condc
      real, dimension(imax), intent(inout)             :: condx
      real, dimension(imax), intent(inout)             :: conds
      real, dimension(imax), intent(inout)             :: condg
      real, dimension(imax), intent(inout)             :: precc
      real, dimension(imax), intent(inout)             :: precip
      real, dimension(imax), intent(inout)             :: timeconv
      real, dimension(imax), intent(inout)             :: so2wd
      real, dimension(imax), intent(inout)             :: so4wd
      real, dimension(imax), intent(inout)             :: bcwd
      real, dimension(imax), intent(inout)             :: ocwd
      real, dimension(imax), intent(inout)             :: saltwd
      real, dimension(imax), intent(out)               :: convpsav
      real, dimension(imax), intent(inout)             :: aug
      integer, dimension(imax), intent(inout)          :: kbsav
      integer, dimension(imax), intent(inout)          :: ktsav
      integer, dimension(imax), intent(inout)          :: kt_saved
      integer, dimension(imax), intent(inout)          :: kb_saved
      logical, dimension(imax), intent(in)             :: land
!from kuocom.h common block
      real, intent(in) :: entrain, detrain, alfsea, fldown, alflnd
      real, intent(in) :: rhcv, rhsat
      real, intent(in) :: convtime, tied_con, convfact, rhmois
      integer, intent(in) :: mbase, iterconv, nuvconv, methdetr
      integer, intent(in) :: methprec, mdelay, ncvcloud, ldr, nkuo
      integer, intent(in) :: nevapls
      
      integer kbsav_ls(imax),kb_sav(imax),kt_sav(imax)
      integer kkbb(imax),kmin(imax)
      real, dimension(imax,kl) :: qqsav,qliqwsav,xtgscav
      real, dimension(kl) :: fscav,xtgtmp
      real, dimension(kl) :: rho,ttsto,qqsto,qqold,qlsto,qlold
      real, dimension(imax) :: conrev,alfqarr,omega,omgtst
      real, dimension(imax) :: convtim_deep,aa,beta
      real delq(imax,kl),dels(imax,kl),delu(imax,kl)
      real delv(imax,kl),dqsdt(imax,kl),es(imax,kl) 
      real fldow(imax),fluxq(imax)
      real fluxt(imax,kl),hs(imax,kl)  
      real phi(imax,kl),qdown(imax),qliqw(imax,kl),delqliqw(imax,kl)
      real entrsav(imax,kl),detrx(imax,kl)
      real fluxh(imax,kl),fluxv0(imax,0:kl-1),fluxv(imax,0:kl-1)
      real qq(imax,kl),qs(imax,kl)
      real rnrt(imax),rnrtc(imax),rnrtcn(imax)
      real s(imax,kl),tdown(imax),tt(imax,kl)
      real dsk(kl),h0(kl),q0(kl),t0(kl)  
      real qplume(imax,kl),splume(imax,kl)
      integer kdown(imax)
      real factr(imax)
      real fluxqs,fluxt_k(kl)
      real gam, dqrx, rKa, Dva, cdls, cflscon, rhodz
      real qpf, Apr, Bpr, Fr, rhoa, Vr, dtev, qr
      real qgdiff, Cev2, qr2, Cevx, alphal, blx
      real evapls, revq, cfls
      real fluxr(imax)

      do k=1,kl
       dsk(k)=-dsig(k)    !   dsk = delta sigma (positive)
      enddo     ! k loop
     
      alfqarr(:)=alfin(:)
      omega(1:imax)=dpsldt(1:imax,komega)   ! JLM
       
      tt(1:imax,:)=t(1:imax,:)       
      qq(1:imax,:)=qg(1:imax,:)      
      phi(:,1)=bet(1)*tt(:,1)  ! moved up here May 2012
      do k=2,kl
       do iq=1,imax
        phi(iq,k)=phi(iq,k-1)+bet(k)*tt(iq,k)+betm(k)*tt(iq,k-1)
       enddo     ! iq loop
      enddo      ! k  loop
       s(1:imax,1)=cp*tt(1:imax,1)+phi(1:imax,1)  ! dry static energy; only level 1 needed here for kkbb calc
          
!      if(nproc==1.and.ntiles==1)
!     &  write(6,*) 'max_alfqarr,alfin:',maxval(alfqarr),maxval(alfin)

#ifndef GPU
      if(ktau==1.and.mydiag)
     &  write(6,"('alfqarr',2f7.3)") alfqarr(idjd)
#endif

!     just does convective; L/S rainfall done later by LDR scheme
      qliqw(:,:)=0.  ! before itn
      kmin(:)=0
      conrev(:)=1000.*ps(:)/(grav*dt) ! factor to convert precip to g/m2/s
      rnrt(:)=0.       ! initialize large-scale rainfall array; before itn
      rnrtc(:)=0.      ! initialize convective  rainfall array; before itn
      kbsav_ls(:)=0    ! for L/S
!!      ktsav(:)=kl-1    ! preset value to show no deep or shallow convection
!!      kbsav(:)=kl-1    ! preset value to show no deep or shallow convection

      nuv=0
      facuv=.1*abs(nuvconv) ! just initial value for gfortran
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
      beta(:)=0. ! MJT suggestion
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
        es(iq,k)=establ(tt(iq,k))
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
!      following used to define kb_sav (via kkbb), method chosen by nbase
!      if(nbase==3)then  ! allows for cooled surface layer
         kkbb(:)=2
         do k=3,k700+1
           do iq=1,imax
            if(s(iq,k)<max(s(iq,1),s(iq,2))+cp.and.kkbb(iq)==k-1)
     &         kkbb(iq)=k
           enddo    ! iq loop
         enddo      ! k loop   
         do iq=1,imax
          if(kkbb(iq)==2)kkbb(iq)=1  ! i.e. suitable PBL top not found
         enddo    ! iq loop
!      endif  ! (nbase==3)

        !qplume(:,kl)=0. ! 0. just for diag prints  ! JLM
        !splume(:,kl)=0. ! 0. just for diag prints  ! JLM
        qplume(:,:) = 0. ! MJT suggestion
        splume(:,:) = 0. ! MJT suggestion
!      following just for setting up kb_sav values, so k700 OK
       qplume(:,1)=alfin(:)*qq(:,1)
       splume(:,1)=s(:,1)  
       if(mbase==0.or.mbase==2)then  ! simple mbase=0 runs before 29/4/18 not correct          
          do iq=1,imax
           k=kkbb(iq)
           splume(iq,k)=s(iq,k)+aug(iq)
           qplume(iq,k)=alfin(iq)*qq(iq,k)
          enddo  ! iq loop
       endif  ! mbase==0.or.mbase==2)
       if(mbase==1)then  ! finds proper maximum values for qplume & splume
         do k=1,k700+1
          qplume(:,k)=qs(:,k) ! only used over sea
         enddo  ! k loop
         qplume(:,1)=alfin(:)*qq(:,1) ! only relevant over land
         do k=2,k700+1
          do iq=1,imax
           if(k<=kkbb(iq))then
            if(land(iq))then
             qplume(iq,k)=max(qplume(iq,k-1),alfin(iq)*qq(iq,k))
            endif
            splume(iq,k)=max(splume(iq,k-1),s(iq,k))
           endif
          enddo  ! iq loop
         enddo  ! k loop
       endif  ! mbase==1)
       if(mbase==3)then  ! same as mbase=0 but qs over sea; ignores alfsea
          do iq=1,imax
           k=kkbb(iq)
           splume(iq,k)=s(iq,k)+aug(iq)
           qplume(iq,k)=alfin(iq)*qq(iq,k)
           if(.not.land(iq))then
             qplume(iq,k)=qs(iq,k)
           endif
          enddo  ! iq loop
       endif  ! mbase==3)
       if(mbase==4)then  ! finds proper maximum values for qplume & splume
         do k=2,k700+1
          do iq=1,imax
           if(k<=kkbb(iq))then
             qplume(iq,k)=max(qplume(iq,k-1),alfin(iq)*qq(iq,k))
	     
             splume(iq,k)=max(splume(iq,k-1),s(iq,k))
           endif
          enddo  ! iq loop
         enddo  ! k loop
!      this will make the assumption that hbase changes at kb_sav add directly to hplume
!      during the calculation of M	 
       endif  ! mbase==4)
       if(mbase==5)then  ! simple with max for qplume, splume          
          do iq=1,imax
           k=kkbb(iq)
           qplume(iq,k)=alfin(iq)*qq(iq,k)
           splume(iq,k)=max(splume(iq,k-1),s(iq,k))
          enddo  ! iq loop
       endif  ! mbase==5)
!      N.B. only need qplume and splume at level kkbb; values above kkbb calculated later 
!      and values below are not used   	
       do iq=1,imax
         k=kkbb(iq)
         qplume(iq,k)=min(qplume(iq,k),max(qs(iq,k),qq(iq,k))) 
       enddo  ! iq loop

#ifndef GPU
        if(ktau==1.and.mydiag)then
         write(6,*) 'itn,iterconv,nuv,nuvconv ',itn,iterconv,nuv,nuvconv
         write(6,*) 'ntest,methdetr,detrain',
     &               ntest,methdetr,detrain
         write(6,*) 'fldown ',fldown
         write(6,*) 'alflnd,alfsea',alflnd,alfsea
        endif  ! (ktau==1.and.mydiag)
#endif
        
        do iq=1,imax
         alfqarr(iq)=qplume(iq,kkbb(iq))/max(qq(iq,kkbb(iq)),1.e-20) ! MJT suggestion
        enddo    ! iq loop             	 
      endif  ! (itn==1) !--------------------------------------
      
      if(itn>1)then
!      redefine qplume, splume at cloud base
       do iq=1,imax
        k=kkbb(iq)     
        qplume(iq,k)=min(alfqarr(iq)*qq(iq,k)
     &                ,max(qs(iq,k),qq(iq,k)))   
        splume(iq,k)=s(iq,k)+aug(iq)  ! later probably use =max(splume(iq,k),s(iq,k))
       enddo    ! iq loop             	 
      endif  ! (itn>1)

!     Following procedure chooses lowest valid cloud bottom (and base) (moved out of itn=1 loop in 2018)
      kb_sav(:)=kl-1
      kdown(:)=1     
      do iq=1,imax
       k=kkbb(iq)
!      if(k>1.and.splume(iq,k)+hl*qplume(iq,k)>hs(iq,k+1).and. ! do not need K>1
       if(splume(iq,k)+hl*qplume(iq,k)>hs(iq,k+1).and. 
     &   qplume(iq,k)>max(qs(iq,k+1),qq(iq,k+1)))then
         kb_sav(iq)=k
         kt_sav(iq)=k+1
       endif  
      enddo    ! iq loop             

#ifndef GPU
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
         pwater0=pwater0-dsig(k)*qq(iq,k)*ps(iq)/grav
        enddo
!       following prints are just preliminary values. Others further down  
        iq=idjd
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

      entrsav(:,:)=0.
      fluxv0(:,:)=0.  
      do iq=1,imax
        fluxv0(iq,kb_sav(iq))=1.  ! unit reference base mass flux (at level kb+.5)
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
            entrsav(iq,k)=-entrain*dsig(k) ! (linear) entrained mass into plume_k
            if(entrain<0.)entrsav(iq,k)=entrain*fluxv0(iq,k-1)*dsig(k) ! (non-linear) entrained mass into plume_k
!           fluxv is net mass into plume_k (assumed well mixed), and
!           subsequently is mass out top of layer for plume_k
!           N.B. qplume is usually > qs at each level
            fluxv0(iq,k)=fluxv0(iq,k-1)+entrsav(iq,k) ! plume reference mass flux
            qplume(iq,k)=(qplume(iq,k-1)*fluxv0(iq,k-1)
     &          +entrsav(iq,k)*qq(iq,k))/fluxv0(iq,k) 
            hbase=splume(iq,k-1)+hl*qplume(iq,k-1)
            if(hbase>hs(iq,k).and.   ! added for safety 22/4/15
     &       qplume(iq,k)>max(qq(iq,k),qs(iq,k)))then
             kt_sav(iq)=k
             splume(iq,k)=(splume(iq,k-1)*fluxv0(iq,k-1)
     &                      +entrsav(iq,k)*s(iq,k))/fluxv0(iq,k)
            else
              kdown(iq)=0
              entrsav(iq,k)=0.  ! not to mess up dels above cloud top
              fluxv0(iq,k)=0.    ! not to mess up dels above cloud top
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
       fluxv0(iq,kt_sav(iq))=0.
       entrsav(iq,kt_sav(iq))=0.
       splume(iq,kt_sav(iq))=splume(iq,kt_sav(iq)-1)
       qplume(iq,kt_sav(iq))=qplume(iq,kt_sav(iq)-1)
      enddo    ! iq loop
      
#ifndef GPU
      if(ntest>0.and.mydiag)then
         write(6,*) 'before methdetr<0 part, kb_sav,kt_sav',
     &       kb_sav(idjd),kt_sav(idjd)
         write (6,"('fluxv0_UP',15f6.3/(8x,15f6.3))")
     &             fluxv0(idjd,1:kt_sav(idjd))
      endif
#endif
!      calculate new detrainments and modify fluxq and entrsav      
       detrx(:,:)=0.  ! this is actually d0 for next few lines
       aa(:)=0.  ! this aa is sumd, now +ve
       if(methdetr==-1)then  
        do k=2,kl-2  
         do iq=1,imax
          if(k>kb_sav(iq))then  
           if(k<=kt_sav(iq))
     &     detrx(iq,k)=-dsig(k)*(sigmh(kb_sav(iq)+1)-sig(k))  ! +ve d0 value at each level k 
           if(k<kt_sav(iq))aa(iq)=detrx(iq,k)+ ! aa(sumd) & entrain*dsig are +ve
     &                     aa(iq)*(1.+entrain*dsig(k)) 
c          if(iq==idjd) write(6,*) 'k,detrx,entrdsig,aa',
c     &          k,detrx(iq,k),entrain*dsig(k),aa(iq)
          endif
         enddo    ! iq loop
        enddo     ! k loop
       elseif(methdetr==-5)then
        do k=2,kl-2
         do iq=1,imax
          if(k>kb_sav(iq))then
           if(k<=kt_sav(iq))
     &     detrx(iq,k)=-dsig(k)*sqrt(sigmh(kb_sav(iq)+1)-sig(k))  ! +ve d0 value at each level k 
           if(k<kt_sav(iq))aa(iq)=detrx(iq,k)+ ! aa(sumd) & entrain*dsig are +ve
     &                     aa(iq)*(1.+entrain*dsig(k))
          endif
         enddo    ! iq loop
        enddo     ! k loop
       elseif(methdetr==-9)then
        do k=2,kl-2
         do iq=1,imax
          if(k>kb_sav(iq))then
           if(k<=kt_sav(iq))
     &     detrx(iq,k)=-dsig(k)  ! this one constant
           if(k<kt_sav(iq))aa(iq)=detrx(iq,k)+ ! aa(sumd) & entrain*dsig are +ve
     &                     aa(iq)*(1.+entrain*dsig(k))
          endif
         enddo    ! iq loop
        enddo     ! k loop
       else  ! i.e. methdetr=-2 or -3
        do k=2,kl-2  
         do iq=1,imax
!         if(k>kb_sav(iq).and.qplume(iq,k)>max(qs(iq,k),qq(iq,k)))then  ! JLM max added for safety 21/4/15
          if(k>kb_sav(iq))then  
           if(k<=kt_sav(iq))
     &     detrx(iq,k)=-dsig(k)*(sigmh(kb_sav(iq)+1)-sig(k))**  ! +ve d0 value at each level k 
     &     abs(methdetr)                ! **** use methdetr=-1, -2 or -3
           if(k<kt_sav(iq))aa(iq)=detrx(iq,k)+ ! aa(sumd) & entrain*dsig are +ve
     &                     aa(iq)*(1.+entrain*dsig(k)) 
c          if(iq==idjd)write(6,*)'k,detrx,entrdsig,aa',
c     &          k,detrx(iq,k),entrain*dsig(k),aa(iq)
          endif
         enddo    ! iq loop
        enddo     ! k loop
       endif   !  (methdetr==-1) .. else ..

#ifndef GPU
       if(ntest>0.and.mydiag)then
         write(6,*)'detrx',detrx(idjd,1:kt_sav(idjd))
         write(6,*)'entrsav',entrsav(idjd,1:kt_sav(idjd))
         write(6,*)'in methdetr<0 part, aa_a=',aa(idjd)
       endif       
#endif
!      and modify fluxv and entrsav  *** assumes entrain<0.  ***     
       do iq=1,imax
         if(kb_sav(iq)<kl-1)beta(iq)=fluxv0(iq,kt_sav(iq)-1)/  ! beta now is +ve scaling factor for detrx
     &         (detrx(iq,kt_sav(iq))+aa(iq))
       enddo    ! iq loop
       fluxv(:,:)=fluxv0(:,:)  ! fluxv is full flux including dettrainment
       do k=2,kl-2  
        do iq=1,imax
         detrx(iq,k)=beta(iq)*detrx(iq,k)  ! before this line detrx is detr0
         if(k>kb_sav(iq).and.k<kt_sav(iq))then  ! maybe can omit this "if" later
            entrsav(iq,k)=entrain*fluxv(iq,k-1)*dsig(k) ! (non-linear) entrained mass into plume_k
            fluxv(iq,k)=fluxv(iq,k-1)+entrsav(iq,k)-detrx(iq,k)
          endif
        enddo    ! iq loop
       enddo     ! k loop
#ifndef GPU
       if(ntest>0.and.mydiag)then
         write(6,*)'in methdetr<0 part, scaled beta=',beta(idjd)
         write(6,*)'detrx',detrx(idjd,1:kt_sav(idjd))
         write(6,*)'entrsav',entrsav(idjd,1:kt_sav(idjd))
         write (6,"('fluxv_up',15f6.3/(8x,15f6.3))")
     &             fluxv(idjd,1:kt_sav(idjd))
       endif
      
      if(nmaxpr==1.and.mydiag)then
        iq=idjd
        write (6,"('qplume',f6.3,11f7.3/(5x,12f7.3))")1000.*qplume(iq,:)
        write (6,"('splume',12f7.2/(5x,12f7.2))")splume(iq,:)/cp 
        write (6,"('entrsav ',15f6.3/(8x,15f6.3))")
     &            entrsav(iq,1:kt_sav(iq))
        write (6,"('detrx ',15f6.3/(8x,15f6.3))")
     &            detrx(iq,1:kt_sav(iq))
      endif   
#endif
        
!     calculate moistening fraction
      if(methprec==5)then
        do iq=1,imax
          aa(iq)=min( 1., detrain+(1.-detrain)*( ! typical detrain is .1 for deep clouds
     &   (.55-min(sig(kb_sav(iq))-sig(kt_sav(iq)), .55)) /(.55-.09))**3)
!        diff=0, .14, .2,.3, .4, .6, .7 gives detrainn=1, .74, .50, .24, .13, .10, .10   for .1 
!        diff=0, .14, .2,.3, .4, .6, .7 gives detrainn=1, .87, .52, .29, .22, .18, .15   for .15 
        enddo
      elseif(methprec==6)then
        do iq=1,imax
          dz=sig(kb_sav(iq))-sig(kt_sav(iq))-.1
          aa(iq)=min( 1., detrain+(1.-detrain-2.*dz)/(1.+40.*dz*dz))
!        diff=0, .14, .2,.3, .4, .6, .7 gives detrainn=1, .71, .6, .29, .165, ..09, .08   for .15  
        enddo
      elseif(methprec==2)then
        do iq=1,imax
          dz=(sig(kb_sav(iq))-sig(kt_sav(iq)))*ps(iq)*1.e-5
          aa(iq)=min( 1., 1.25/(1.+25.*dz*dz))
!        diff=0, .14, .2,.3, .4, .6, .7 gives detrainn=1, .86, .67, .43, .29, ..15, .11   for 20
!        diff=0, .14, .2,.3, .4, .6, .7 gives detrainn=1, .84, .625, .38, .25, ..125, .094   for 25
        enddo
      elseif(methprec==3)then
        do iq=1,imax
          dz=(sig(kb_sav(iq))-sig(kt_sav(iq)))*ps(iq)*1.e-5
          aa(iq)=min( 1., 1.3/(1.+30.*dz*dz))
!        diff=0, .14, .2,.3, .4, .6, .7 gives detrainn=1, .82, .59, .35, .22, ..11, .08   for 30  
        enddo
      elseif(methprec>0)then  ! e.g. methprec= 20, 25, 30, 60   
       do iq=1,imax
        dz=(sig(kb_sav(iq))-sig(kt_sav(iq)))*ps(iq)*1.e-5
        aa(iq)=max(detrain, min( 1.,(1.+.01*real(methprec))/
     &                 (1.+real(methprec)*dz*dz)) )
!        aa(iq)=min( 1.,(1.+.01*real(methprec))/
!     &                 (1.+real(methprec)*dz*dz))
!       diff=0, .14, .2,.3, .4, .6, .7 gives detrainn=1, .86, .67, .43, .29, ..15, .11   for 20
!       diff=0, .14, .2,.3, .4, .6, .7 gives detrainn=1, .84, .625, .38, .25, ..125, .094   for 25
!       diff=0, .14, .2,.3, .4, .6, .7 gives detrainn=1, .82, .59, .35, .22, ..11, .08   for 30  
! cf prev diff=0, .14, .2,.3, .4, .6, .7 gives detrainn=1, .74, .50, .24, .13, .10, .10   for .1 
       enddo
      endif
       
!     apply entrainment, detrainment and moistening effects to environment    
       do k=2,kl-2
        do iq=1,imax
          dels(iq,k)=-entrsav(iq,k)*s(iq,k)   ! entr into updraft
     &      +detrx(iq,k)*splume(iq,k)
!         following line allows for qq>qs, to avoid drying layer
          qsk=max(qs(iq,k),qq(iq,k))
          delq(iq,k)=-entrsav(iq,k)*qq(iq,k)+detrx(iq,k)*qsk  ! entr into updraft & detr to env
c         rnrt_k=detrx(iq,k)*max(0.,qplume(iq,k)-qsk) ! max not need as such a detxmax will be 0
          rnrt_k=detrx(iq,k)*(qplume(iq,k)-qsk)
          dels(iq,k)=dels(iq,k)+hl*rnrt_k    ! corresponding precip. heating (as not done separately in plume)
!         part of this is detrained liquid water
          delqliqw(iq,k)=delqliqw(iq,k)+aa(iq)*rnrt_k   ! moistening
          rnrt_k= rnrt_k-delqliqw(iq,k)        
          rnrtcn(iq)=rnrtcn(iq)+rnrt_k
         enddo     ! iq loop
        enddo  ! k loop       
       
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
       if(kb_sav(iq)<kl-1)then
         fluxh(iq,kdown(iq))=fldow(iq)
         do k=1,  kb_sav(iq)   ! +ve into plume ###### BELOW cloud base ###### 
          fluxh(iq,k)=upin(k,kb_sav(iq))-fldow(iq)*downex(k,kb_sav(iq))
          fluxv(iq,k)=fluxv(iq,k-1)+fluxh(iq,k)  ! needed for subsidence calcs below cloud base
!         calculate emergent downdraft properties
          delq(iq,k)=fldow(iq)*downex(k,kb_sav(iq))*qdown(iq)
          dels(iq,k)=fldow(iq)*downex(k,kb_sav(iq))*
     &       (s(iq,kb_sav(iq))+cp*(tdown(iq)-t(iq,kb_sav(iq))))  ! correct
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
        write(6,*) 'kb_sav,kt_sav',
     .           kb_sav(iq),kt_sav(iq)
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
!     calculate net base mass flux 
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
      
      if(nkuo==21)then     ! this option incorporates entrainment effects
      aa(:)=0. 
      do k=2,kl-1
       do iq=1,imax
!       want: new_hplume(k)>=new_hs(k)
!       [h+M*aa/|dsig|+alfsarr*M*dels+alfqarr*M*hl*delq]_base 
!                                          = [hs+M*dels+M*hlcp*dels*dqsdt]_k
        if(k>kb_sav(iq).and.k<=kt_sav(iq))then
          ccc=dels(iq,k)*(1.+hlcp*dqsdt(iq,k))
          bbb=dels(iq,kb_sav(iq))+
     &      alfqarr(iq)*hl*delq(iq,kb_sav(iq))
!       aa gives accumulated entrainment contribs from dels & delq, working up from kb_sav+1     
!       N.B. in following line both entrain and dsig are -ve     
        if(k<kt_sav(iq))aa(iq)=aa(iq)+entrain*fluxv0(iq,k-1)*dsig(k)*
     &      (dels(iq,k)+hl*delq(iq,k))
       fluxt(iq,k)=(splume(iq,k-1)+hl*qplume(iq,k-1)-hs(iq,k))/
     &  max(1.e-9,ccc+aa(iq)/dsig(k)-bbb)
          if(sig(k)<sig(kb_sav(iq))                 ! rhcv not RH but for depth-related flux test
     &       -rhcv*(sig(kb_sav(iq))-sig(kt_sav(iq))).and.    !  typically rhcv=.1; 0 for backward-compat
     &       fluxt(iq,k)<convpsav(iq))then  
            convpsav(iq)=fluxt(iq,k)
            kmin(iq)=k   ! level where fluxt is a minimum (diagnostic)
          endif    ! (sig(k)<sig(kb_sav(iq)-....and.)
        endif   ! (k>kb_sav(iq).and.k<kt_sav(iq))
       enddo    ! iq loop
      enddo     ! k loop      
      endif   ! nkuo==21)

      if(nkuo==22)then ! this one like nkuo=23 code (partial entr effects)  
      do k=2,kl-1
       do iq=1,imax
!       want: new_hplume(k)>=new_hs(k) {without entr, hplume(k)=hplume(kbase)  }
!       [h+alfsarr*M*dels+alfqarr*M*hl*delq]_base 
!                                          = [hs+M*dels+M*hlcp*dels*dqsdt]_k
        if(k>kb_sav(iq).and.k<=kt_sav(iq))then
          fluxt(iq,k)=(splume(iq,k-1)+hl*qplume(iq,k-1)-hs(iq,k))/
     .       max(1.e-9,dels(iq,k)*(1.+hlcp*dqsdt(iq,k))       ! 0804 for max
     .                -dels(iq,kb_sav(iq))                    ! to avoid zero
     .                -alfqarr(iq)*hl*delq(iq,kb_sav(iq)) )   ! with real*4
!          hbase=splume(iq,kb_sav(iq))+hl*qplume(iq,kb_sav(iq))
!          den1=splume(iq,k)+hl*qplume(iq,k)-hbase
!          fluxt(iq,k)=(splume(iq,k-1)+hl*qplume(iq,k-1)-hs(iq,k))/
!     .       max(1.e-9,dels(iq,k)*(1.+hlcp*dqsdt(iq,k))       ! 0804 for max
!     .                -dels(iq,kb_sav(iq)) -den1              ! to avoid zero
!     .                -alfqarr(iq)*hl*delq(iq,kb_sav(iq)) )   ! with real*4
          if(sig(k)<sig(kb_sav(iq))                 ! rhcv not RH but for depth-related flux test
     &       -rhcv*(sig(kb_sav(iq))-sig(kt_sav(iq))).and.    !  typically rhcv=.1; 0 for backward-compat
     &       fluxt(iq,k)<convpsav(iq))then  
            convpsav(iq)=fluxt(iq,k)
            kmin(iq)=k   ! level where fluxt is a minimum (diagnostic)
          endif    ! (sig(k)<sig(kb_sav(iq)-....and.)
        endif   ! (k>kb_sav(iq).and.k<kt_sav(iq))
       enddo    ! iq loop
      enddo     ! k loop    

#ifndef GPU
      if(diag.and.mydiag)then    ! JLM
       do k=2,kl-1
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
      endif   ! nkuo==22)

      do iq=1,imax
       if(dels(iq,kb_sav(iq))+alfqarr(iq)*hl*delq(iq,kb_sav(iq))>=0.)
     &          convpsav(iq)=0. 
       if(dels(iq,kt_sav(iq))<=0.)convpsav(iq)=0.    ! JLM 1505 must stabilize
      enddo    ! iq loop
      
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
         else   !   (tied_con==0.  nowadays usual)
c           write(6,*)'has tied_con=0'
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
         write(6,*) 'sig_kb_sav,sig_k',sig(kb_sav(iq)),sig(kt_sav(iq))
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
           sumb=min(1.,  (sig(kb_sav(iq))-sig(kt_sav(iq)))/convt_frac) ! typically /0.6
!  or taub=convtim_deep(iq)* convt_frac/min(convt_frac,sig(kb_sav(iq))-sig(kt_sav(iq))
!         convtim_deep in mins between mcontlnd and mcontsea, linearly
           factr(iq)=min(1.,sumb*dt/(60.*convtim_deep(iq)))  ! defined just for itn=1
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
!         tied_b=26 gives factr [1, .964, .900, .794,.5.216] for ds = [200, 100, 50, 25,8,2} km     
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
       do k=2,kl-1 
        do iq=1,imax  
         qliqw(iq,k)=qliqw(iq,k)+convpsav(iq)*delqliqw(iq,k)    ! after applying convpsav, accumulating for itns
        enddo  ! iq loop
      enddo   ! k loop

#ifndef GPU
      if(ntest>0.and.mydiag)then
        iq=idjd
         write(6,*)'liqw',convpsav(iq),kb_sav(iq),kt_sav(iq)
!       N.B. convpsav(iq) is already mult by dt jlm: mass flux is convpsav/dt
        write(6,*) "after convection: ktau,itn,kbsav,ktsav ",

     .                 ktau,itn,kb_sav(iq),kt_sav(iq)
        write (6,"('qgc ',12f7.3/(8x,12f7.3))") 
     .             (1000.*qq(iq,k),k=1,kl)
        write (6,"('ttc ',12f7.2/(8x,12f7.2))") 
     .             (tt(iq,k),k=1,kl)
        write(6,*) 'rnrtc',rnrtc(iq)
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
     &         delq_av,-delq_av*hl/cp,-delq_av*conrev(iq)
        if(abs(delt_av)>1.e-20)write(6,*) 
     &        'ktau,itn,kbsav,ktsav,delt_av,heatlev',
     &        ktau,itn,kb_sav(iq),kt_sav(iq),delt_av,heatlev/delt_av
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
        do k=1,kl-2   
         do iq=1,imax
          u(iq,k)=u(iq,k)+facuv*factr(iq)*convpsav(iq)*delu(iq,k)/dsk(k)
          v(iq,k)=v(iq,k)+facuv*factr(iq)*convpsav(iq)*delv(iq,k)/dsk(k)
         enddo  ! iq loop
        enddo   ! k loop
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
          xtgscav(1:imax,1:kl) = 0.
          s(1:imax,1:kl) = xtg(1:imax,1:kl,ntr)
!$acc loop private(ttsto,qqsto,qlsto,qqold,qlsto,rho,xtgtmp)
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
     &                  kt-kb+1)
               
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

      if(nmaxpr==1.and.mydiag)then
       iq=idjd
       if(ktau==1.and.itn==1)write (6, "(15x,
     &  'ktau itn kb_sav kmin kdown kt_sav cfrac+ entxsav detrx ',
     &  ' fluxv  fluxvt-1 factr convpsav3 rnrtcn3 cape')")
       if(kb_sav(iq)<kl-1)write (6,"('ktau ... cape ',
     &   i5,i4,']',2i5,2i6,f9.2,5f8.2,3p2f8.3,1pf9.1)")  ktau,itn,
     &   kb_sav(iq),kmin(iq),kdown(iq),kt_sav(iq),
     &   cfrac(iq,kb_sav(iq)+1),entrsav(iq,kdown(iq)),
     &   detrx(iq,kdown(iq)),fluxv(iq,kdown(iq)),
     &   fluxv(iq,kt_sav(iq)-1),factr(iq),convpsav(iq),
     &   rnrtcn(iq),cape(iq)
       dtsol=.01*sgsave(iq)/(1.+.25*(u(iq,1)**2+v(iq,1)**2))   ! solar heating   ! sea
       write (6,"('pblh,fldow,tdown,qdown,fluxq3',
     &             f8.2,f5.2,2f7.2,3p2f8.3,1pf8.2)")
     &     pblh(iq),fldow(iq),tdown(iq),qdown(iq),fluxq(iq)
       write(6,"('ktau,kkbb,wetfac,dtsol,alfqarr,fg,omega,'
     &  'omgtst',i5,i3,3f6.3,f8.3,2f9.3)") ktau,kkbb(idjd),
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
     &   write (6,"('C k,fluxtot6',i3,f7.2)")
     &   (k,1.e6*fluxtot(idjd,k),k=1,kl)
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
        !qliqw(1:imax,:)=0.   ! just for final diags
!       following is old 2014 ldr=0 code for large-scale calculations 
        do k=1,kl   
         do iq=1,imax
          es(iq,k)=establ(tt(iq,k))
         enddo  ! iq loop
         if(sig(k)> .25)klon2=k
        enddo   ! k loop
        do k=klon2,1,-1    ! top down to end up with proper kbsav_ls
         do iq=1,imax
          pk=ps(iq)*sig(k)
          qs(iq,k)=.622*es(iq,k)/max(pk-es(iq,k),1.)  
          if(qq(iq,k)>rhsat*qs(iq,k))then
            kbsav_ls(iq)=k
            gam=hlcp*qs(iq,k)*pk*hlars/(tt(iq,k)**2*max(pk-es(iq,k),1.))
            dqrx=(qq(iq,k)-rhsat*qs(iq,k))/(1.+rhsat*gam)
            tt(iq,k)=tt(iq,k)+hlcp*dqrx
            qq(iq,k)=qq(iq,k)-dqrx
            rnrt(iq)=rnrt(iq)+dqrx*dsk(k)
          endif   ! (qq(iq,k)>rhsat*qs(iq,k))
         enddo    ! iq loop
        enddo     ! k loop
!!!     now do evaporation of L/S precip (ONLY for ldr=0) !!
        rnrt(:)=rnrt(:)*conrev(:)                 
!       here the rainfall rate rnrt has been converted to g/m**2/sec
        fluxr(:)=rnrt(:)*1.e-3*dt ! kg/m2      
#ifndef GPU
        if(nmaxpr==1.and.mydiag)then
          iq=idjd
          write(6,*) 'after large scale rain: kbsav_ls,rnrt,convpsav',
     &                kbsav_ls(iq),rnrt(iq),convpsav(iq)
          write (6,"('qliqw ',12f7.3/(8x,12f7.3))") 
     &              (1000.*qliqw(iq,k),k=1,kl)
          write (6,"('qs ',12f7.3/(8x,12f7.3))") 
     &              (1000.*qs(iq,k),k=1,kl)
          write (6,"('qg ',12f7.3/(8x,12f7.3))") 
     &              (1000.*qq(iq,k),k=1,kl)
          write (6,"('tt ',12f7.2/(8x,12f7.2))") 
     &              (tt(iq,k),k=1,kl)
        endif
#endif
        qliqw(1:imax,:)=0.   ! just for final diags
        if(nevapls>.0)then ! even newer UKMO (just for ldr=0)
         rKa=2.4e-2
         Dva=2.21
         cfls=1. ! cld frac7 large scale
         cflscon=4560.*cfls**.3125
         do k=kl,1,-1  ! JLM
          do iq=1,imax
           if(k<kbsav_ls(iq))then
             rhodz=ps(iq)*dsk(k)/grav
             qpf=fluxr(iq)/rhodz     ! Mix ratio of rain which falls into layer
             pk=ps(iq)*sig(k)
!            es(iq,k)=qs(iq,k)*pk/.622
             Apr=hl*hl/(rKa*tt(iq,k)*(rvap*tt(iq,k))-1.)
             Bpr=rvap*tt(iq,k)*pk/(Dva*es(iq,k))
             Fr=fluxr(iq)/(cfls*dt)
             rhoa=pk/(rdry*tt(iq,k))
             dz=pk/(rhoa*grav)
             Vr=max(.01 , 11.3*Fr**(1./9.)/sqrt(rhoa)) ! Actual fall speed
             dtev=dz/Vr
             qr=fluxr(iq)/(dt*rhoa*Vr)
             qgdiff=qs(iq,k)-qq(iq,k)
             Cev2=cflscon*qgdiff/(qs(iq,k)*(Apr+Bpr))  ! Ignore rhoa**0.12
             qr2=max(0. , qr**.3125 - .3125*Cev2*dtev)**3.2
             Cevx=(qr-qr2)/dtev  ! i.e. Cev*qgdiff
             alphal=hl*qs(iq,k)/(ars*tt(iq,k)**2)
             blx=qgdiff+Cevx*dt*(1.+hlcp*alphal)  ! i.e. bl*qgdiff
             evapls= cfls*dt*Cevx*qgdiff/blx      ! UKMO
             evapls=max(0. , min(evapls,qpf))
             qpf=qpf-evapls
             fluxr(iq)=fluxr(iq)+rhodz*qpf
             revq=evapls
             revq=min(revq , rnrt(iq)/(conrev(iq)*dsk(k)))
!            max needed for roundoff
             rnrt(iq)=max(1.e-10 , rnrt(iq)-revq*dsk(k)*conrev(iq))
             tt(iq,k)=tt(iq,k)-revq*hlcp
             qq(iq,k)=qq(iq,k)+revq
           endif   ! (k<kbsav_ls(iq))
          enddo    ! iq loop
         enddo     ! k loop
        endif      ! (nevapls>0)
        cfrac=.3   ! just a default for radn
      endif  ! (ldr.ne.0)
!_______________end of convective calculations__________________
     
!      if(ldr.ne.0)go to 8
!      obsolete ldr=0 code removed, for large-scale calculations 

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
      precip(1:imax)=precip(1:imax)+condx(1:imax)
      t(1:imax,:)=tt(1:imax,:)             

#ifndef GPU
      if(ntest>0.or.(ktau<=2.and.nmaxpr==1))then  ! diag print near bottom
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
      return
      end subroutine convjlm22_work
      
      end module convjlm22_m
