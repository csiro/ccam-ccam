      module trvmix


      contains
c ***************************************************************************
      subroutine tracervmix(at,ct)
!     next to lines in JLM version of tracervmix - not sure if need
      use cc_mpi, only : mydiag
      use diag_m
      use sigs_m
      use tracermodule, only :tracunit,tracname
      use tracers_m  ! tr
c     this routine does the vertical mixing of tracers
      implicit none
      include 'newmpar.h'
      include 'const_phys.h'      ! grav,fair_molm,fc_molm
      include 'parm.h'          ! dt
      real trsrc(ifull,kl)
      real updtr(ilt*jlt,klt,ngasmax),at(ifull,kl),ct(ifull,kl)
      real trfact,molfact,radfact,co2fact,gasfact
      logical decay,methloss,mcfloss
      integer igas

      trfact = grav * dt / dsig(1)
c rml 25/08/04 factor for units in mol/m2/s
      molfact = 1000.*trfact*fair_molm
      co2fact=1000.*trfact*fair_molm/fc_molm
c     o2fact=1000.*trfact*fair_molm/fo2_molm
c rml 08/11/04 test factor for radon units in Bq/m2/s, conc in Bq/m3
c     density of air at 273K and 1013hPa = 1.293 kg/m3
      radfact = trfact*1.293

      do igas=1,ngas                  

        call trgassflux(igas,trsrc)
c rml 25/08/04 change gasfact to be depend on tracer flux units
        if (trim(tracunit(igas)).eq.'gC/m2/s') then
          gasfact = co2fact
          decay = .false.
        elseif (trim(tracunit(igas)).eq.'mol/m2/s') then
          gasfact = molfact
          decay = .false.
c rml 08/11/04 add radon case
        elseif (trim(tracunit(igas)).eq.'Bq/m2/s') then
          gasfact = radfact
          decay = .true.
        else
c         assume no surface flux so gasfact could be anything but we'll 
c         set it to zero
          gasfact = 0.
        endif
c rml 05/12/05 also set decay for tracer name 'radon' in case not in Bq/m2/s
        if (trim(tracname(igas)).eq.'radon'.or.
     &           tracname(igas)(1:2).eq.'Rn') decay=.true.

! rml 16/2/10 check for methane tracers to set flag to do loss
        methloss=.false.
        if (tracname(igas)(1:7).eq.'methane') methloss=.true.
! rml 30/4/10 check for mcf tracers to set flag to do loss
        mcfloss = .false.
        if (tracname(igas)(1:3).eq.'mcf') mcfloss=.true.
c
c rml 08/11/04 add decay flag
        call gasvmix(updtr(:,:,igas), gasfact, igas, decay,trsrc,
     &               methloss,mcfloss)
      enddo
      call trimt(ngas,at,ct,updtr,0)
      tr(1:ilt*jlt,:,1:ngasmax)=updtr(1:ilt*jlt,:,1:ngasmax)
      return
      end subroutine

c ***************************************************************************
      subroutine trgassflux(igas,trsrc)
      use cable_ccam, only : cbmemiss
      use carbpools_m ! online co2 fluxes
      use define_dimensions, only : ncs, ncp ! Used in carbpool.h
      use nsibd_m     !ivegt (vegetation type)
      use tracermodule, only :co2em,tractype,tracname,tracdaytime
c     this routine put the correct tracer surface flux into trsrc
      implicit none
      include 'newmpar.h'
      include 'dates.h'  ! timeg
!     tml 17/09/07 online tracers by veg type
c     can these common blocks be 'lost' with rewrite of cbm/soilsnow?
      integer igas
      real trsrc(ifull,kl)
!     rml 17/09/07 online tracers by veg type
      integer nchar, mveg

c     initialise (to allow for ocean gridpoints for cbm fluxes)      
!     and non surface layers
      trsrc = 0.

!     rml 2/10/03 allow for online (cbm) tracers or flux from file
!     rml 17/09/07 rewrite online tracer section to include tracers
!     separated by vegetation type (also moved if/elseif to case)
      if (trim(tractype(igas)).eq.'online') then
!       write(6,*) 'adding surface flux for ',trim(tracname(igas))
        if (trim(tracname(igas)(1:3)).eq.'cbm') then
          select case (trim(tracname(igas)))
            case('cbmnep'); trsrc(:,1) = fnee
            case('cbmpn'); trsrc(:,1) = fpn
            case('cbmrp'); trsrc(:,1) = frp
            case('cbmrs'); trsrc(:,1) = frs
            case default ; stop 'unknown online tracer name'
          end select
        else
          nchar = len_trim(tracname(igas))
          read(tracname(igas)(nchar-1:nchar),'(i2)',err=101) mveg
c         write(131,*) 'tracer test: ',mveg
          if (mveg.lt.1.or.mveg.gt.maxval(ivegt)) stop 
     &      'tracer selection: veg type out of range'
          select case (tracname(igas)(1:nchar-2))
            case('gpp')
              call cbmemiss(trsrc(:,1),mveg,1)
            case('plresp')
              call cbmemiss(trsrc(:,1),mveg,2)
            case('slresp')
              call cbmemiss(trsrc(:,1),mveg,3)
            case default 
              stop 'unknown online tracer name'
          end select
        endif
      elseif (trim(tractype(igas)).eq.'daypulseon') then
c       only add flux during day time
        if (tracdaytime(igas,1).lt.tracdaytime(igas,2) .and.
     &          tracdaytime(igas,1).le.timeg .and.
     &          tracdaytime(igas,2).ge.timeg) then
          trsrc(:,1) = co2em(:,igas)
        elseif (tracdaytime(igas,1).gt.tracdaytime(igas,2) .and.
     &          (tracdaytime(igas,1).le.timeg .or.
     &          tracdaytime(igas,2).ge.timeg)) then
          trsrc(:,1) = co2em(:,igas)
        else
          trsrc(:,1) = 0.
        endif
      else
c       emissions from file
        trsrc(:,1) = co2em(:,igas)
      endif

      return
 101  stop 'unknown online tracer name or veg type number error'
      end subroutine
c *****************************************************************
      subroutine gasvmix(temptr,fluxfact,igas,decay,trsrc,methloss,
     &                   mcfloss)
! rml 16/2/10 addition for methane
! rml 30/4/10 addition for mcf
      use arrays_m        ! ps
      use cc_mpi, only : myid
      use sigs_m          ! disg,bet,betm
      use tracermodule, only : oh,strloss,unit_trout,mcfdep,jmcf
      use tracers_m       ! tr
      use xyzinfo_m       ! wts
      implicit none
      include 'newmpar.h'
      include 'const_phys.h' !eradsq,pi,grav
      include 'mpif.h'
      include 'parm.h'    ! dt
      real trsrc(ilt*jlt,kl)
      real temptr(ilt*jlt,klt)
      real loss(ilt*jlt,kl)
! rml 30/4/10 deposition for mcf
      real dep(ilt*jlt),dz(ifull),zg1,zg2
c rml 08/11/04 decay flag to all decay for radon
      logical decay,methloss,mcfloss
      real drate,fluxfact
      integer igas
      real koh,totloss_l,totloss,kohmcf
      parameter(koh=2.45e-12,kohmcf=1.64e-12)
      integer ierr,k,iq

! rml 30/4/10 Since decay, loss and deposition need total tracer field
!  add trback in at start of this subroutine and remove again at end
      tr(1:ilt*jlt,:,igas) = tr(1:ilt*jlt,:,igas)+trback_g(igas)

c rml 08/11/04 decay rate for radon (using units of source, Bq/m2/s, to
c indicate that radon and need decay
      if (decay) then
        drate = exp(-dt*2.11e-6)
      else
        drate = 1.
      endif
c
! rml 16/2/10 methane loss by OH and in stratosphere
      if (methloss) then
        loss(:,:) = tr(1:ilt*jlt,:,igas)*dt*
     &     (koh*exp(-1775./t(1:ilt*jlt,:))*oh(:,:) + strloss(:,:))
!
!       calculate total loss
        totloss_l = 0.
        do k=1,kl
          do iq=1,ilt*jlt
            totloss_l = totloss_l + loss(iq,k)*dsig(k)*ps(iq)*wts(iq)
          enddo
        enddo
        call MPI_Allreduce(totloss_l,totloss,1,MPI_REAL,MPI_SUM,
     &                     MPI_COMM_WORLD,ierr)
!       convert to TgCH4 and write out
        if (myid == 0) then
          totloss = -1.*totloss*4*pi*eradsq*fCH4_MolM/
     &                 (grav*fAIR_MolM*1.e18)
          write(6,*) 'Total loss',ktau,totloss
          write(unit_trout,*) 'Total loss',ktau,totloss
!         accumulate loss over month
          acloss_g(igas) = acloss_g(igas) + totloss
        endif
        dep=1.
      elseif (mcfloss) then
        loss(:,:) = tr(1:ilt*jlt,:,igas)*dt*
     &     (kohmcf*exp(-1520./t(1:ilt*jlt,:))*oh(:,:) + jmcf(:,:))
!       deposition
!       calculate thickness of surface layer
        do iq=1,ifull
!        zg1=bet(1)*t(iq,1)/grav
!        zg2=zg1+(bet(2)*t(iq,2)+betm(2)*t(iq,1))/grav
!        dz(iq) = 0.5*(zg1+zg2) - zs(iq)
! rml 16/7/10 dz formula provided by J McGregor
         dz(iq) = t(iq,1)*(1.-sigmh(2))*rdry/(grav*sig(1))
        enddo      ! iq loop

        dep =  exp(-1.*mcfdep(:)*dt/dz)
      else
        loss = 0.
        dep = 1.
      endif

      temptr(:,1) = tr(1:ilt*jlt,1,igas)*drate *dep(:)
     &              - fluxfact*trsrc(:,1)/ps(1:ilt*jlt)
     &              - loss(:,1)
      temptr(:,2:kl) = tr(1:ilt*jlt,2:kl,igas)*drate 
     &                 - loss(:,2:kl)

!     remove trback from tr and temptr
      tr(1:ilt*jlt,:,igas) = tr(1:ilt*jlt,:,igas)-trback_g(igas)
      temptr(1:ilt*jlt,:) = temptr(1:ilt*jlt,:)-trback_g(igas)
 


      return
      end subroutine
c *********************************************************************
      subroutine trimt(ngas,a,c,rhs,it)
c     This is a copy of trim.f but trying to do all tracers at once.  
c     u initially now contains rhs; leaves with answer u (jlm)
c     n.b. we now always assume b = 1-a-c
      implicit none
      include 'newmpar.h'
!     N.B.  e, g, temp are just work arrays (not passed through at all)     
      real e(ifull,kl),g(ifull,kl,ngas)
      real temp(ifull,kl)
      real a(ifull,kl),c(ifull,kl),rhs(ifull,kl,ngas)
      real b
      integer it,iq,k,ngas


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
       g(iq,1,:)=rhs(iq,1,:)/b
      enddo
      do k=2,kl-1
       do iq=1,ifull
        g(iq,k,:)=(rhs(iq,k,:)-a(iq,k)*g(iq,k-1,:))*temp(iq,k)
       enddo
      enddo

c     do back substitution to give answer now
      do iq=1,ifull
       b=1.-a(iq,kl)-c(iq,kl)
       rhs(iq,kl,:)=(rhs(iq,kl,:)-a(iq,kl)*g(iq,kl-1,:))/
     .            (b-a(iq,kl)*e(iq,kl-1))
      enddo
      do k=kl-1,1,-1
       do iq=1,ifull
        rhs(iq,k,:)=g(iq,k,:)-e(iq,k)*rhs(iq,k+1,:)
       enddo
      enddo
      return
      end subroutine
c ********************************************************************
      end module
