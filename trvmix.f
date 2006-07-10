      module trvmix


      contains
c ***************************************************************************
      subroutine tracervmix(at,ct)
!     next to lines in JLM version of tracervmix - not sure if need
      use cc_mpi, only : mydiag
      use diag_m
      use tracermodule, only :tracunit,tracname
c     this routine does the vertical mixing of tracers
      implicit none
      include 'newmpar.h'
      include 'const_phys.h'      ! grav,fair_molm,fc_molm
      include 'parm.h'          ! dt
      include 'sigs.h'          ! dsig
      include 'tracers.h'       ! tr
      real trsrc(ifull,kl)
      real updtr(ilt*jlt,klt,ngasmax),at(ifull,kl),ct(ifull,kl)
      real trfact,molfact,radfact,co2fact,gasfact
      logical decay
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
        if (trim(tracname(igas)).eq.'radon') decay=.true.
c
c rml 08/11/04 add decay flag
        call gasvmix(updtr(:,:,igas), gasfact, igas, decay,trsrc)
      enddo
      call trimt(ngas,at,ct,updtr,0)
      tr(1:ilt*jlt,:,1:ngasmax)=updtr(1:ilt*jlt,:,1:ngasmax)
      return
      end subroutine

c ***************************************************************************
      subroutine trgassflux(igas,trsrc)
      use tracermodule, only :co2em,tractype,tracname,tracdaytime
c     this routine put the correct tracer surface flux into trsrc
      implicit none
      include 'newmpar.h'
      include 'dates.h'  ! timeg
      include 'carbpools.h' ! online co2 fluxes
c     can these common blocks be 'lost' with rewrite of cbm/soilsnow?
      integer igas
      real trsrc(ifull,kl)

c     initialise (to allow for ocean gridpoints for cbm fluxes)      
!     and non surface layers
      trsrc = 0.

!     rml 2/10/03 allow for online (cbm) tracers or flux from file
      if (trim(tractype(igas)).eq.'online') then
!       write(6,*) 'adding surface flux for ',trim(tracname(igas))
        if (trim(tracname(igas)).eq.'cbmnep') then
          trsrc(:,1) = fnee
        elseif (trim(tracname(igas)).eq.'cbmpn') then
          trsrc(:,1) = fpn 
        elseif (trim(tracname(igas)).eq.'cbmrp') then
          trsrc(:,1) = frp
        elseif (trim(tracname(igas)).eq.'cbmrs') then
          trsrc(:,1) = frs
        else
          stop 'unknown online tracer name'
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
      end subroutine
c *****************************************************************
      subroutine gasvmix(temptr, fluxfact, igas, decay,trsrc)

      implicit none
      include 'newmpar.h'
      include 'arrays.h'  ! ps
      include 'tracers.h' ! tr
      include 'parm.h'    ! dt
      real trsrc(ilt*jlt,kl)
      real temptr(ilt*jlt,klt)
c rml 08/11/04 decay flag to all decay for radon
      logical decay
      real drate,fluxfact
      integer igas

c rml 08/11/04 decay rate for radon (using units of source, Bq/m2/s, to
c indicate that radon and need decay
      if (decay) then
        drate = exp(-dt*2.11e-6)
      else
        drate = 1.
      endif
c
      temptr(:,1) = tr(1:ilt*jlt,1,igas)*drate 
     &              - fluxfact*trsrc(:,1)/ps(1:ilt*jlt)
      temptr(:,2:kl) = tr(1:ilt*jlt,2:kl,igas)*drate 


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
