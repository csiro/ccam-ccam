c     co2.f bundles together co2sflux and other gases
c                 as well as co2vmix  and other gases
      subroutine co2sflux
      use cc_mpi, only : mydiag
      include 'newmpar.h'
      include 'aalat.h'     ! along
      include 'arrays.h'
      include 'const_phys.h'
      include 'dates.h'     ! timeg
      include 'extraout.h'
      include 'map.h'       ! zs,land
      include 'nsibd.h'     ! rsmin,ivegt,sigmf,tgf,ssdn,res,rmc,tsigmf
      include 'parm.h'      ! nsib parameter passed
      include 'scamdim.h'
c     include 'soil.h'
      include 'soilsnow.h'  ! tgg,wb,snowd
      include 'tracers.h'
      include 'trcom2.h' ! nstn,slat,slon,istn,jstn, nstn2 etc.
      include 'pbl.h'    ! slwa

      parameter (ftcrit = 293.1, fco2const =3.17097e-7*1.5)
c     parameter (ftcrit = 293.1, fco2const = 3.1536e-6)
      parameter (ntest=0)  ! usually 0, 1 for special diag prints

!!    common/totco2/snetr(ifull,3)
      common/work3/vmixarrs(ifull,kl,3),trsrc(ifull,kl),spare(ifull,kl)

c     global matrices of:
c      real rco2(ifull),                            ! co2 respiration
c     -     phsco2(ifull),                          ! photosynthesis
c     -     pnfco2                                  ! net flux

c     ecosystem characteristics:
      real bmsco2(7),                              ! biomass
     -     nppco2(7),                              ! net primary production
     -     betaco2(7)                              ! turnover time

      data bmsco2/20000.,17000.,3000.,1000.,300.,1000.,0./
      data nppco2/250.e-7,180.e-7,80.e-7,80.e-7,10.e-7,90.e-7,1./
      data betaco2/3.1e-9,3.1e-9,6.3e-9,1.6e-8,1.6e-8,1.6e-8,1./

      integer ivegmap(44)
!     data ivegmap/1,1,2,3,3,5,3,3,3,3,
!    &             3,5,4,4,5,4,5,4,4,4,
!    &             5,4*6,5,5,3,3*5,
      data ivegmap/4*1,11*3,6*4,6*6,3,3*5,
     &             1,2,3,2,3,3,4,3,3,3,5,6,    5/
c      integer ivegmap(13)
c      data ivegmap/1,2,3,2,3,3,4,3,3,3,5,6,7/

!     srcmin1= 1.e6
!     srcmin2= 1.e6
!     srcmax1=-1.e6
!     srcmax2=-1.e6
c     ico2em(19,65)=0.  ! for no sydney emissions c20
      do iq=1,ifull
       pnfco2=0.
!      n.b. c-c model doesn't mind industrial sources over the sea
       co2inem=ico2em(iq)*fco2const
!      timel=timeg+along(iq)/15.        ! local time between 0 and 48 h
!      if(timel.gt.24.and.iq.eq.idjd)then
!        timela=timeg+along(iq)/15.
!        print *,'iq,timeg,timela,timel ',ii,timeg,timela,timel
!        print *,'slwa,sgsave,rgsave',
!     .           slwa(iq),sgsave(iq),rgsave(iq)
!      endif
!      if(timel.gt.24.)timel=timel-24.    ! local time between 0 and 24 h
!      if( timel.gt.6.and.timel.le.18.) then                    ! day
!        if(ico2em(iq).gt.20) then
!          co2inem=ico2em(iq)*fco2const*4./3.
!        else
!          co2inem=ico2em(iq)*fco2const*2.
!        endif
!      else                                                     ! night
!        if(ico2em(iq).gt.20) then
!          co2inem=ico2em(iq)*fco2const*2./3.
!        else
!          co2inem=0.
!        endif
!      endif

        if( land(iq).and.snowd(iq).eq.0.)then
          iveg=ivegt(iq)
          ivegco2=ivegmap(iveg)
          coef=.1
          if(t(iq,1).lt. ftcrit) coef=.05
          rco2= max(betaco2(ivegco2)*bmsco2(ivegco2)*
     -             (1.+(t(iq,1)-ftcrit)*coef) ,0.)
          rco2bg= max(betaco2(5)*bmsco2(5)*
     -               (1.+(t(iq,1)-ftcrit)*coef) ,0.)
          phsco2  =nppco2(ivegco2)*sgsave(iq)/50.
          phsco2bg=nppco2(5)*sgsave(iq)/50.
            pnfco2=tsigmf(iq) *(rco2  -phsco2  ) +
     .             (1.-tsigmf(iq))*(rco2bg-phsco2bg)
          pnfo2= -2.6667*pnfco2

        elseif(.not.land(iq))then      !  sea
          pco2=abs(alat(iq))*4./3.     ! not used in c-c runs
c         if( abs(alat(iq)) .gt. 45. ) pco2=0.
c         if(wspeed(iq).gt.13.) then
c           tk=5.9*wspeed(iq)-49.3
c         else if(wspeed(iq).gt.3.6) then
c           tk=2.85*wspeed(iq)-9.65
c         else
c           tk=.17*wspeed(iq)
c         end if
c         rco2=-2.e-9*tk*pco2
c         pnfco2=rco2
c         pnfo2  = -10.*pnfco2
        endif

        trsrc(iq,1)   = pnfco2+co2inem
c       trsrc(iq,2)   = pnfo2            ! check why ,2
!!      snetr(iq)=snetr(iq)+trsrc(iq)
!       srcmin1=min(srcmin1,trsrc(iq))
!       srcmin2=min(srcmin2,trsrc(iq,2))
!       srcmax1=max(srcmax1,trsrc(iq))
!       srcmax2=max(srcmax2,trsrc(iq,2))

      enddo ! iq=1,ifull

      if(ntest.eq.1.and.mydiag)then
       iq=idjd
       timel=timeg+along(iq)/15.        ! local time between 0 and 48 h
       print *,'ktau,iq,along,timeg,timel ',
     .          ktau,iq,along(iq),timeg,timel
       print *,'slwa,sgsave,rgsave',slwa(iq),sgsave(iq),rgsave(iq)
       iveg=ivegt(iq)
       ivegco2=ivegmap(iveg)
       phsco2  =nppco2(ivegco2)*sgsave(iq)/50.
       print *,'land,iveg,ivegco2,phsco2,tsigmf ',
     .          land(iq),iveg,ivegco2,phsco2,tsigmf(iq)
      endif   ! (ntest.eq.1)
      return
      end
      subroutine o2sflux
      include 'const_phys.h'
      include 'newmpar.h'
      include 'tracers.h'
      common/work3/vmixarrs(il,jl,kl,3),trsrc(il,jl,kl),spare(ifull,kl)

c --- o2 surface flux is computed in co2sflux() and returned in trsrc(i,j,2), so it must
c       be called before the o2 routines
      if( ico2.eq.0) then
         write(*,*)
     . ' o2sflux: co2 must be present for o2 calculations, fatal error'
         stop ' ... terminating in o2sflux'
      end if
      return
      end

      subroutine radonsflux
      include 'const_phys.h'
      include 'newmpar.h'
      include 'aalat.h'
      include 'arrays.h'
      include 'extraout.h'
      include 'map.h'
      include 'nsibd.h'    ! rsmin,ivegt,sigmf,tgf,ssdn,res,rmc,tsigmf
      include 'scamdim.h'
      include 'sigs.h'     ! sigmh
c     include 'soil.h'
      include 'soilsnow.h'  ! tgg,wb,snowd
      include 'tracers.h'
      include 'trcom2.h'   ! itrace,nstn,slat,slon,istn,jstn
      include 'vvel.h'
      common/work3/vmixarrs(ifull,kl,3),trsrc(ifull,kl),spare(ifull,kl)

c     common/totco2/snetr(il,jl,3)

      logical firstcall
      save    firstcall
      data    firstcall/.true./

      real radonm(9,7)
c      data radonm/ .8  ,.8  ,.95 ,.25 ,.67 ,.87 ,.6  ,.4  ,.4
      data radonm/ .8  ,.8  ,.95 ,.25 ,.67 ,.87 ,1.0  ,1.0  ,1.0
     -            ,.8  ,.8  ,.95 ,.25 ,.67 ,.87 ,.64 ,.38 ,.18
     -            ,.7  ,1.8 ,4.4 ,.43 ,.67 ,.9  ,1.5 ,2.4 ,1.5
     -            ,.59 ,1.5 ,2.  ,.5  ,.7  ,1.5 ,1.5 ,1.5 ,1.5
     -            ,.93 ,.94 ,.8  ,.8  ,.76 ,1.84,4.1 ,1.5 ,1.5
     -            ,.93 ,.94 ,.5  ,.84 ,.39 ,1.5 ,2.5 ,1.5 ,1.5
     -            ,.0  ,.0  ,.0  ,4.6 ,2.  ,1.5 ,1.5 ,1.5 ,1.5/

      save x, y

      if(npanels.gt.0)then   ! globpe
        do iq=1,ifull
         if(snowd(iq).gt.0.)then
           trsrc(iq,1)=0.
         else
           trsrc(iq,1)=radonem(iq)
         endif
        enddo
      else   ! i.e. for eva's 125 km darlam runs
*       initialise conversion factors and coefficients
        if( firstcall ) then
!         for darlam iauste etc defined in indata
          iaustw=il/4
          iauste=iaustw+8
          jausts=jl/4
          jaustn=-999
          x=(iauste-iaustw)/9.
          y=(jaustn-jausts)/7.
          firstcall = .not.firstcall
          print *,'radon iaustw iauste jausts jaustn',iaustw,iauste,
     &    jausts,jaustn
        endif
        do j=1,jl
         do i=1,il
	   iq=i+(j-1)*il
          if( land(iq) ) then
           trsrc(iq,1)=1.
           if( j.ge.jausts.and.j.le.jaustn) then
             xi=(i - iaustw)/x
             ii=xi+ 1
             ii=min(max(1,ii),9)
             yj=(j - jausts)/y
             jj=yj+ 1
             jj=min(max(1,jj),7)
             trsrc(iq,1)=radonm(ii,jj)
           endif    ! ( j.ge.jausts.and.j.le.jaustn)
          else
           trsrc(iq,1) = 0.
          endif
         enddo ! i=1,il
        enddo  ! j=1,jl
        stop 'radonsflux not yet ready for darlam'
      endif  ! (npanels.gt.0)
      return
      end

      subroutine co2vmix(updtr, fluxfact )
      include 'newmpar.h'
      include 'arrays.h'
      include 'tracers.h'
      common/work3/vmixarrs(ifull,kl,3),trsrc(ifull,kl),spare(ifull,kl)
      real updtr(ilt*jlt,klt)

      do iq=1,ilt*jlt
        updtr(iq,1)=tr(iq,1,max(1,ico2)) -
     .                fluxfact*trsrc(iq,1)/ps(iq)
      end do

      do k=2,klt     
       do iq=1,ilt*jlt
        updtr(iq,k)=tr(iq,k,max(1,ico2))
       enddo
      enddo
      return
      end

      subroutine radonvmix(updtr, fluxfact )
      include 'const_phys.h'
      include 'newmpar.h'
      include 'arrays.h'
      include 'parm.h'    ! dt
      include 'tracers.h'
      common/work3/vmixarrs(ifull,kl,3),trsrc(ifull,kl),spare(ifull,kl)
      real updtr(ilt*jlt,klt)
      parameter ( fhalflife_ra = 60.*60.*24.*3.8, flog2 = 0.693147 )
      parameter ( flife_ra = fhalflife_ra/flog2 )

      decay=exp(-dt/flife_ra)
      do iq=1,ilt*jlt
        updtr(iq,1)=tr(iq,1,max(1,iradon))*decay -
     .                fluxfact*trsrc(iq,1)/ps(iq)
      end do

      do k=2,klt     
       do iq=1,ilt*jlt
        updtr(iq,k)=tr(iq,k,max(1,iradon))*decay            ! wed  01-04-1995
       enddo
      enddo
      return
      end

      subroutine o2vmix(updtr, fluxfact )
      include 'newmpar.h'
      include 'arrays.h'
      include 'tracers.h'
      common/work3/vmixarrs(ifull,kl,3),trsrc(ifull,kl),spare(ifull,kl)
      real updtr(ilt*jlt,klt)
      do iq=1,ilt*jlt
        updtr(iq,1)=tr(iq,1,max(1,io2)) -
     .                fluxfact*trsrc(iq,2)/ps(iq)
      end do
      do k=2,klt     
       do iq=1,ilt*jlt
        updtr(iq,k)=tr(iq,k,max(1,io2))
       enddo
      enddo
      return
      end

