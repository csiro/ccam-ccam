      subroutine staguv(u,v,uout,vout)   !  unstaguv as an entry
!     stripped down version just with nstag=nstagu=-3, mstagpt=-3      
!     staguv    may be called from adjust5, upglobal
!     unstaguv  may be called from adjust5,  nonlin
c     nstag now in parm.h  ! for nstag=0   staguv uses cubic interpolation
c                          ! for nstag=1   staguv goes via cartesian components
c                          ! for nstag=2   staguv linear with em terms
c                          ! for nstag=3-5 staguv: jlm reversible 2/1/98
c                          ! for nstag=7   staguv: jlm reversible Akima
c                          ! for nstag=10  staguv: jlm rev_stag
c                          !    only -3 used nowadays
c     nstagu now in parm.h ! same but for unstaguv
c     N.B. staguv & unstaguv previously required 2D arrays as input
c     - no longer, as transferred here to uin and vin
c     assume k level already given
      include 'newmpar.h'
      include 'indices.h' ! in,is,iw,ie,inn,iss,iww,iee
      include 'map.h'
      include 'parm.h'
      include 'parmdyn.h'
      include 'vecsuv.h'   ! vecsuv info
      include 'vecsuva.h'  ! vecsuva info
      real u(ifull),v(ifull),uout(ifull),vout(ifull)
      real ua(2*ifull),va(ifull),ud(2*ifull),vd(ifull) ! work arrays
      real ub(2*ifull),vb(ifull)                       ! work arrays
      equivalence (ua(ifull+1),va),(ud(ifull+1),vd)    ! to ensure contiguous
      equivalence (ub(ifull+1),vb)                     ! to ensure contiguous
      real uin(2*ifull),vin(ifull)                     ! work arrays
      equivalence (uin(ifull+1),vin)                   ! to ensure contiguous
      data itnmax/3/

c     unstaggered u & v as input; staggered as output
      do iq=1,ifull
       uin(iq)=u(iq)
       vin(iq)=v(iq)
      enddo   ! iq loop

!     if(nstag.eq.-3)then
!       call reva6(uin,vin,uout,vout,ua,va,ub,vb,ud,vd)
!     endif    !  (nstag.eq.-3)
!cdir nodep
c       do iq=1,ifull   ! precalculate rhs terms
c         ud(iwu2(iq))= uin(iwu2(iq))/2.+uin(iq)+uin(ieu2(iq))/10.
c         vd(isv2(iq))= vin(isv2(iq))/2.+vin(iq)+vin(inv2(iq))/10.
c       enddo
         do iq=2,ifull-1
          ud(iq-1)= uin(iq-1)/2.+uin(iq)+uin(iq+1)/10.
         enddo
         do iq=il+1,ifull-il
          vd(iq-il)= vin(iq-il)/2.+vin(iq)+vin(iq+il)/10.
         enddo
	  do j=1,jl
	   iq=1+(j-1)*il
          ud(iwu2(iq))= uin(iwu2(iq))/2.+uin(iq)+uin(ieu2(iq))/10.
	   iq=il+(j-1)*il
          ud(iwu2(iq))= uin(iwu2(iq))/2.+uin(iq)+uin(ieu2(iq))/10.
	  enddo
	  do j=1,jl,2*il
	   do i=1,il
	    iq=i+(j-1)*il
           vd(isv2(iq))= vin(isv2(iq))/2.+vin(iq)+vin(inv2(iq))/10.
 	   enddo
	  enddo
	  do j=2*il,jl,2*il
	   do i=1,il
	    iq=i+(j-1)*il
           vd(isv2(iq))= vin(isv2(iq))/2.+vin(iq)+vin(inv2(iq))/10.
 	   enddo
	  enddo
c        do iq=1,ifull
c          ua(iq)=ud(iq)-ud(ieu2(iq))/2.   ! 1st guess
c          va(iq)=vd(iq)-vd(inv2(iq))/2.   ! 1st guess
c        enddo
!cdir nodep
         do iq=1,ifull-1
          ua(iq)=ud(iq)-ud(iq+1)/2.
         enddo
         do iq=1,ifull-il
          va(iq)=vd(iq)-vd(iq+il)/2.
         enddo
	  do j=1,jl
	   iq=il+(j-1)*il
          ua(iq)=ud(iq)-ud(ieu2(iq))/2.    ! 1st guess
	  enddo
	  do j=2*il,jl,2*il
	   do i=1,il
	    iq=i+(j-1)*il
           va(iq)=vd(iq)-vd(inv2(iq))/2.   ! 1st guess
 	   enddo
	  enddo

        do itn=1,itnmax  ! each loop is a double iteration
!cdir nodep
c         do iq=1,ifull
c          ub(iq)=ua(ieu2(iq))
c          vb(iq)=va(inv2(iq))
c         enddo
         do iq=1,ifull-1
          ub(iq)=ua(iq+1)
         enddo
         do iq=1,ifull-il
          vb(iq)=va(iq+il)
         enddo
	  do j=1,jl
	   iq=il+(j-1)*il
          ub(iq)=ua(ieu2(iq))
	  enddo
	  do j=2*il,jl,2*il
	   do i=1,il
	    iq=i+(j-1)*il
           vb(iq)=va(inv2(iq))
 	   enddo
	  enddo
!cdir nodep
c       do iq=1,ifull
c         uin(iq)=(ud(iq)-.5*ud(ieu2(iq))
c    .                 -ua(iwu2(iq))/10. +ub(ieu2(iq))/4.)/.95
c         vin(iq)=(vd(iq)-.5*vd(inv2(iq))
c    .                 -va(isv2(iq))/10. +vb(inv2(iq))/4.)/.95
c        enddo
	  do iq=2,ifull-1
          uin(iq)=(ud(iq)-.5*ud(iq+1)
     .                 -ua(iq-1)/10. +ub(iq+1)/4.)/.95
         enddo
         do iq=il+1,ifull-il
          vin(iq)=(vd(iq)-.5*vd(iq+il)
     .                 -va(iq-il)/10. +vb(iq+il)/4.)/.95
         enddo
	  do j=1,jl
	   iq=1+(j-1)*il
          uin(iq)=(ud(iq)-.5*ud(ieu2(iq))
     .                 -ua(iwu2(iq))/10. +ub(ieu2(iq))/4.)/.95
	   iq=il+(j-1)*il
          uin(iq)=(ud(iq)-.5*ud(ieu2(iq))
     .                 -ua(iwu2(iq))/10. +ub(ieu2(iq))/4.)/.95
	  enddo
	  do j=1,jl,2*il
	   do i=1,il
	    iq=i+(j-1)*il
           vin(iq)=(vd(iq)-.5*vd(inv2(iq))
     .                 -va(isv2(iq))/10. +vb(inv2(iq))/4.)/.95
 	   enddo
	  enddo
	  do j=2*il,jl,2*il
	   do i=1,il
	    iq=i+(j-1)*il
           vin(iq)=(vd(iq)-.5*vd(inv2(iq))
     .                 -va(isv2(iq))/10. +vb(inv2(iq))/4.)/.95
 	   enddo
	  enddo

!cdir nodep
c         do iq=1,ifull
c          ub(iq)=uin(ieu2(iq))
c          vb(iq)=vin(inv2(iq))
c         enddo
         do iq=1,ifull-1
          ub(iq)=uin(iq+1)
         enddo
         do iq=1,ifull-il
          vb(iq)=vin(iq+il)
         enddo
	  do j=1,jl
	   iq=il+(j-1)*il
          ub(iq)=uin(ieu2(iq))
	  enddo
	  do j=2*il,jl,2*il
	   do i=1,il
	    iq=i+(j-1)*il
           vb(iq)=vin(inv2(iq))
 	   enddo
	  enddo
!cdir nodep
c        do iq=1,ifull
c         ua(iq)=(ud(iq)-.5*ud(ieu2(iq))
c    .                 -uin(iwu2(iq))/10. +ub(ieu2(iq))/4.)/.95
c         va(iq)=(vd(iq)-.5*vd(inv2(iq))
c    .                 -vin(isv2(iq))/10. +vb(inv2(iq))/4.)/.95
c        enddo
	  do iq=2,ifull-1
          ua(iq)=(ud(iq)-.5*ud(iq+1)
     .                 -uin(iq-1)/10. +ub(iq+1)/4.)/.95
         enddo
         do iq=il+1,ifull-il
          va(iq)=(vd(iq)-.5*vd(iq+il)
     .                 -vin(iq-il)/10. +vb(iq+il)/4.)/.95
         enddo
	  do j=1,jl
	   iq=1+(j-1)*il
          ua(iq)=(ud(iq)-.5*ud(ieu2(iq))
     .                 -uin(iwu2(iq))/10. +ub(ieu2(iq))/4.)/.95
	   iq=il+(j-1)*il
          ua(iq)=(ud(iq)-.5*ud(ieu2(iq))
     .                 -uin(iwu2(iq))/10. +ub(ieu2(iq))/4.)/.95
	  enddo
	  do j=1,jl,2*il
	   do i=1,il
	    iq=i+(j-1)*il
           va(iq)=(vd(iq)-.5*vd(inv2(iq))
     .                 -vin(isv2(iq))/10. +vb(inv2(iq))/4.)/.95
 	   enddo
	  enddo
	  do j=2*il,jl,2*il
	   do i=1,il
	    iq=i+(j-1)*il
           va(iq)=(vd(iq)-.5*vd(inv2(iq))
     .                 -vin(isv2(iq))/10. +vb(inv2(iq))/4.)/.95
 	   enddo
	  enddo
        enddo  ! itn=1,itnmax
	 
        do iq=1,ifull      ! final values for output
         uout(iq)=ua(iq)
         vout(iq)=va(iq)
        enddo   ! iq loop
      return

      entry unstaguv(u,v,uout,vout)
c     staggered u & v as input; unstaggered as output
      do iq=1,ifull
       uin(iq)=u(iq)
       vin(iq)=v(iq)
      enddo   ! iq loop

!     if(nstagu.eq.-3)then
!       call revb6(uin,vin,uout,vout,ua,va,ub,vb,ud,vd)
!     endif    !  (nstagu.eq.-3)
 
!cdir nodep
c       do iq=1,ifull   ! precalculate rhs terms
c         ud(ieu2(iq))= uin(ieu2(iq))/2.+uin(iq)+uin(iwu2(iq))/10.
c         vd(inv2(iq))= vin(inv2(iq))/2.+vin(iq)+vin(isv2(iq))/10.
c       enddo
         do iq=2,ifull-1
          ud(iq+1)= uin(iq+1)/2.+uin(iq)+uin(iq-1)/10.
         enddo
         do iq=il+1,ifull-il
          vd(iq+il)= vin(iq+il)/2.+vin(iq)+vin(iq-il)/10.
         enddo
	  do j=1,jl
	   iq=1+(j-1)*il
          ud(ieu2(iq))= uin(ieu2(iq))/2.+uin(iq)+uin(iwu2(iq))/10.
	   iq=il+(j-1)*il
          ud(ieu2(iq))= uin(ieu2(iq))/2.+uin(iq)+uin(iwu2(iq))/10.
	  enddo
	  do j=1,jl,2*il
	   do i=1,il
	    iq=i+(j-1)*il
           vd(inv2(iq))= vin(inv2(iq))/2.+vin(iq)+vin(isv2(iq))/10.
 	   enddo
	  enddo
	  do j=2*il,jl,2*il
	   do i=1,il
	    iq=i+(j-1)*il
           vd(inv2(iq))= vin(inv2(iq))/2.+vin(iq)+vin(isv2(iq))/10.
 	   enddo
	  enddo
c        do iq=1,ifull
c          ua(iq)=ud(iq)-ud(iwu2(iq))/2.   ! 1st guess
c          va(iq)=vd(iq)-vd(isv2(iq))/2.   ! 1st guess
c        enddo
!cdir nodep
         do iq=2,ifull
          ua(iq)=ud(iq)-ud(iq-1)/2.
         enddo
         do iq=il+1,ifull
          va(iq)=vd(iq)-vd(iq-il)/2.
         enddo
	  do j=1,jl
	   iq=1+(j-1)*il
          ua(iq)=ud(iq)-ud(iwu2(iq))/2.    ! 1st guess
	  enddo
	  do j=1,jl,2*il
	   do i=1,il
	    iq=i+(j-1)*il
           va(iq)=vd(iq)-vd(isv2(iq))/2.   ! 1st guess
 	   enddo
	  enddo

        do itn=1,itnmax  ! each loop is a double iteration
!cdir nodep
c         do iq=1,ifull
c          ub(iq)=ua(iwu2(iq))
c          vb(iq)=va(isv2(iq))
c         enddo
!cdir nodep
         do iq=2,ifull
          ub(iq)=ua(iq-1)
         enddo
         do iq=il+1,ifull
          vb(iq)=va(iq-il)
         enddo
	  do j=1,jl
	   iq=1+(j-1)*il
          ub(iq)=ua(iwu2(iq))
	  enddo
	  do j=1,jl,2*il
	   do i=1,il
	    iq=i+(j-1)*il
           vb(iq)=va(isv2(iq))
 	   enddo
	  enddo
!cdir nodep
c        do iq=1,ifull
c         uin(iq)=(ud(iq)-.5*ud(iwu2(iq))
c    .                 -ua(ieu2(iq))/10. +ub(iwu2(iq))/4.)/.95
c         vin(iq)=(vd(iq)-.5*vd(isv2(iq))
c    .                 -va(inv2(iq))/10. +vb(isv2(iq))/4.)/.95
c        enddo
         do iq=2,ifull-1
          uin(iq)=(ud(iq)-.5*ud(iq-1)
     .                 -ua(iq+1)/10. +ub(iq-1)/4.)/.95
         enddo
         do iq=il+1,ifull-il
          vin(iq)=(vd(iq)-.5*vd(iq-il)
     .                 -va(iq+il)/10. +vb(iq-il)/4.)/.95
         enddo
	  do j=1,jl
	   iq=1+(j-1)*il
          uin(iq)=(ud(iq)-.5*ud(iwu2(iq))
     .                 -ua(ieu2(iq))/10. +ub(iwu2(iq))/4.)/.95
	   iq=il+(j-1)*il
          uin(iq)=(ud(iq)-.5*ud(iwu2(iq))
     .                 -ua(ieu2(iq))/10. +ub(iwu2(iq))/4.)/.95
	  enddo
	  do j=1,jl,2*il
	   do i=1,il
	    iq=i+(j-1)*il
           vin(iq)=(vd(iq)-.5*vd(isv2(iq))
     .                 -va(inv2(iq))/10. +vb(isv2(iq))/4.)/.95
 	   enddo
	  enddo
	  do j=2*il,jl,2*il
	   do i=1,il
	    iq=i+(j-1)*il
           vin(iq)=(vd(iq)-.5*vd(isv2(iq))
     .                 -va(inv2(iq))/10. +vb(isv2(iq))/4.)/.95
 	   enddo
	  enddo

!cdir nodep
c         do iq=1,ifull
c          ub(iq)=uin(iwu2(iq))
c          vb(iq)=vin(isv2(iq))
c         enddo
!cdir nodep
         do iq=2,ifull
          ub(iq)=uin(iq-1)
         enddo
         do iq=il+1,ifull
         vb(iq)=vin(iq-il)
         enddo
	  do j=1,jl
	   iq=1+(j-1)*il
          ub(iq)=uin(iwu2(iq))
	  enddo
	  do j=1,jl,2*il
	   do i=1,il
	    iq=i+(j-1)*il
           vb(iq)=vin(isv2(iq))
 	   enddo
	  enddo
!cdir nodep
c        do iq=1,ifull
c         ua(iq)=(ud(iq)-.5*ud(iwu2(iq))
c    .                 -uin(ieu2(iq))/10. +ub(iwu2(iq))/4.)/.95
c         va(iq)=(vd(iq)-.5*vd(isv2(iq))
c    .                 -vin(inv2(iq))/10. +vb(isv2(iq))/4.)/.95
c        enddo
         do iq=2,ifull-1
          ua(iq)=(ud(iq)-.5*ud(iq-1)
     .                 -uin(iq+1)/10. +ub(iq-1)/4.)/.95
         enddo
         do iq=il+1,ifull-il
          va(iq)=(vd(iq)-.5*vd(iq-il)
     .                 -vin(iq+il)/10. +vb(iq-il)/4.)/.95
         enddo
	  do j=1,jl
	   iq=1+(j-1)*il
          ua(iq)=(ud(iq)-.5*ud(iwu2(iq))
     .                 -uin(ieu2(iq))/10. +ub(iwu2(iq))/4.)/.95
	   iq=il+(j-1)*il
          ua(iq)=(ud(iq)-.5*ud(iwu2(iq))
     .                 -uin(ieu2(iq))/10. +ub(iwu2(iq))/4.)/.95
	  enddo
	  do j=1,jl,2*il
	   do i=1,il
	    iq=i+(j-1)*il
           va(iq)=(vd(iq)-.5*vd(isv2(iq))
     .                 -vin(inv2(iq))/10. +vb(isv2(iq))/4.)/.95
 	   enddo
	  enddo
	  do j=2*il,jl,2*il
	   do i=1,il
	    iq=i+(j-1)*il
           va(iq)=(vd(iq)-.5*vd(isv2(iq))
     .                 -vin(inv2(iq))/10. +vb(isv2(iq))/4.)/.95
 	   enddo
	  enddo
        enddo  ! itn=1,itnmax

        do iq=1,ifull      ! final values for output
         uout(iq)=ua(iq)
         vout(iq)=va(iq)
        enddo   ! iq loop
      return
      end        
