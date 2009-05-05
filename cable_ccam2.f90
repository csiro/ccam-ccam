module cable_ccam

  ! CABLE interface originally developed by the CABLE group
  ! Subsequently modified by MJT for mosaic

  use define_types, cbm_ms => ms

  private
  public CABLE,sib4,loadcbmparm,savetile

  integer, parameter :: CABLE = 4
  integer, dimension(:), allocatable, save :: cmap
  integer, dimension(5,2), save :: pind  
  real, dimension(:), allocatable, save :: atmco2
  real, dimension(:), allocatable, save :: sv,vl
  integer, parameter :: CO2forcingtype=1   ! 1 constant, 2 prescribed 1900-2004,
                                           ! 3 interactive
  TYPE (air_type)             :: air
  TYPE (bgc_pool_type)        :: bgc
  TYPE (canopy_type)          :: canopy
  TYPE (met_type)             :: met
  TYPE (balances_type)        :: bal
  TYPE (radiation_type)       :: rad
  TYPE (roughness_type)       :: rough
  type (soil_parameter_type)  :: soil       ! soil parameters
  TYPE (soil_snow_type)       :: ssoil
  TYPE (sum_flux_type)        :: sum_flux
  type (veg_parameter_type)   :: veg        ! vegetation parameters
  ! Save these so only have to do the allocation once.
  save air, bgc, canopy, met, bal, rad, rough, soil, ssoil, &
       sum_flux, veg

  contains
  ! ****************************************************************************

  subroutine sib4

      use zenith_m
      USE cbm_module
      use physical_constants, only : umin
  
      implicit none

      include 'newmpar.h'
      include 'arrays.h'
      include 'carbpools.h'
      include 'const_phys.h' ! grav
      include 'dates.h' ! ktime,kdate,timer,timeg,xg,yg
      include 'extraout.h'
      include 'latlong.h'  ! rlatt,rlongg
      include 'morepbl.h'
      include 'nsibd.h'
      include 'parm.h'
      include 'permsurf.h'
      include 'pbl.h'
      include 'prec.h'
      include 'screen.h'   ! tscrn etc
      include 'sigs.h'
      include 'soil.h'     ! ... zmin zolnd zolog sice alb
      include 'soilbal.h'
      include 'soilsnow.h' ! 
      include 'vegpar.h' ! 
!                     met forcing for CBM
      real dirad,dfgdt,degdt,wetfac,degdw,cie
      real factch,qsttg,rho,zo,aft,fh,ri,theta
      real gamm,rg,vmod,dummwk2
      common/work2/dirad(ifull),dfgdt(ifull),degdt(ifull)  &
     & ,wetfac(ifull),degdw(ifull),cie(ifull)              &
     & ,factch(ifull),qsttg(ifull),rho(ifull),zo(ifull)    &
     & ,aft(ifull),fh(ifull),ri(ifull),theta(ifull)        &
     & ,gamm(ifull),rg(ifull),vmod(ifull),dummwk2(ifull)

      ! for calculation of zenith angle
      real fjd, r1, dlt, slag, dhr, coszro2(ifull),taudar2(ifull)
      real bpyear,alp
      real tmps(ifull)

      integer jyear,jmonth,jday,jhour,jmin
      integer mstart,ktauplus,k,mins,kstart
      integer nb

      integer imonth(12)
      data imonth /31,28,31,30,31,30,31,31,30,31,30,31/
      integer ndoy(12)   ! days from beginning of year (1st Jan is 0)
      data ndoy/ 0,31,59,90,120,151,181,212,243,273,304,334/
!
!      set meteorological forcing
!
       jyear=kdate/10000
       jmonth=(kdate-jyear*10000)/100
       jday=kdate-jyear*10000-jmonth*100
       jhour=ktime/100
       jmin=ktime-jhour*100
       ! mins from start of year
       mstart=1440*(ndoy(jmonth)+jday-1) + 60*jhour + jmin
       ktauplus=0
       do k=1,jmonth-1
        ktauplus = ktauplus + imonth(k)*nperday
       enddo

       ! mtimer contains number of minutes since the start of the run.
       mins = mtimer + mstart
       bpyear = 0.
       fjd = float(mod(mins,525600))/1440.  ! 525600 = 1440*365
       call solargh(fjd,bpyear,r1,dlt,alp,slag)
       dhr = dt/3600.
       call zenith(fjd,r1,dlt,slag,rlatt,rlongg,dhr,ifull,coszro2,taudar2)

       call setco2for(jyear)

       kstart = 1

       met%doy=fjd
       met%tk=theta(cmap)
       met%ua=vmod(cmap)
       met%ca=1.e-6*atmco2(cmap)
       met%coszen=max(1.e-8,coszro2(cmap)) ! use instantaneous value
         
       met%fld=-rgsave(cmap)        ! long wave down  
       met%qv=qg(cmap,1)        ! specific humidity in kg/kg
       met%pmb=.01*ps(cmap)     ! pressure in mb at ref height
       met%precip=condx(cmap)
       met%hod=(met%doy-int(met%doy))*24.0 + rlongg(cmap)*180./(15.*pi)
       where (met%hod.gt.24.0) met%hod=met%hod-24.0
       rough%za=-287.*t(cmap,1)*log(sig(1))/grav   ! reference height
       met%fsd=sgsave(cmap)/(1.-0.5*sum(albvisnir(cmap,:),2))! short wave down (positive) W/m^2
       
       veg%vlai(:)=max(vl(:),0.1) ! for updating LAI

       met%tc=met%tk - 273.16
       met%ua=max(met%ua,umin)
       met%precip_s=0. ! in mm not mm/sec
       where ( met%tc < 0.0 ) met%precip_s = met%precip
       
       if (ktau.eq.1.and.all(ssoil%rtsoil.eq.0.)) ssoil%rtsoil=rtsoil(cmap)

!      rml 21/09/07 remove ktauplus+ktau due to change in cable_offline
       CALL cbm(ktau, kstart, ntau, dt, air, bgc, canopy, met, &
     &      bal, rad, rough, soil, ssoil, sum_flux, veg, mxvt, mxst)

      ! average diagnostic fields
      albsoilsn(iperm(1:ipland),:)=0.
      albsav(iperm(1:ipland))=0.
      albnirsav(iperm(1:ipland))=0.
      runoff(iperm(1:ipland))=0.
      rnof1(iperm(1:ipland))=0.
      rnof2(iperm(1:ipland))=0.
      wbtot(iperm(1:ipland))=0.
      tevap(iperm(1:ipland))=0.
      tprecip(iperm(1:ipland))=0.
      totenbal(iperm(1:ipland))=0.
      trnoff(iperm(1:ipland))=0.
      tgg(iperm(1:ipland),:)=0.
      wb(iperm(1:ipland),:)=0.
      wbice(iperm(1:ipland),:)=0.
      tggsn(iperm(1:ipland),:)=0.
      smass(iperm(1:ipland),:)=0.
      ssdn(iperm(1:ipland),:)=0.
      ssdnn(iperm(1:ipland))=0.
      snowd(iperm(1:ipland))=0.
      !osnowd(iperm(1:ipland))=0.
      !osnowd0(iperm(1:ipland))=0.
      snage(iperm(1:ipland))=0.
      rtsoil(iperm(1:ipland))=0.
      rnet(iperm(1:ipland))=0.
      fg(iperm(1:ipland))=0.
      eg(iperm(1:ipland))=0.
      epot(iperm(1:ipland))=0.
      tss(iperm(1:ipland))=0.
      tgf(iperm(1:ipland))=0.
      cansto(iperm(1:ipland))=0.
      gflux(iperm(1:ipland))=0.
      sgflux(iperm(1:ipland))=0.
      fnee(iperm(1:ipland))=0.
      fpn(iperm(1:ipland))=0.
      frd(iperm(1:ipland))=0.
      frp(iperm(1:ipland))=0.
      frpw(iperm(1:ipland))=0.
      frpr(iperm(1:ipland))=0.
      frs(iperm(1:ipland))=0.
      sumpn(iperm(1:ipland))=0.
      sumrp(iperm(1:ipland))=0.
      sumrpw(iperm(1:ipland))=0.
      sumrpr(iperm(1:ipland))=0.
      sumrs(iperm(1:ipland))=0.
      sumrd(iperm(1:ipland))=0.
      cplant(iperm(1:ipland),:)=0.
      csoil(iperm(1:ipland),:)=0.
      zo(iperm(1:ipland))=0.
      tmps=0.
      do nb=1,5
        do k=1,ms
          tgg(cmap(pind(nb,1):pind(nb,2)),k)=tgg(cmap(pind(nb,1):pind(nb,2)),k) &
                                            +sv(pind(nb,1):pind(nb,2))*ssoil%tgg(pind(nb,1):pind(nb,2),k)
          wb(cmap(pind(nb,1):pind(nb,2)),k)=wb(cmap(pind(nb,1):pind(nb,2)),k) &
                                            +sv(pind(nb,1):pind(nb,2))*ssoil%wb(pind(nb,1):pind(nb,2),k)
          wbice(cmap(pind(nb,1):pind(nb,2)),k)=wbice(cmap(pind(nb,1):pind(nb,2)),k) &
                                               +sv(pind(nb,1):pind(nb,2))*ssoil%wbice(pind(nb,1):pind(nb,2),k)
        end do
        do k=1,ncp
          cplant(cmap(pind(nb,1):pind(nb,2)),k)=cplant(cmap(pind(nb,1):pind(nb,2)),k) &
                                                +sv(pind(nb,1):pind(nb,2))*bgc%cplant(pind(nb,1):pind(nb,2),k)
        enddo
        do k=1,ncs
          csoil(cmap(pind(nb,1):pind(nb,2)),k)=csoil(cmap(pind(nb,1):pind(nb,2)),k) &
                                               +sv(pind(nb,1):pind(nb,2))*bgc%csoil(pind(nb,1):pind(nb,2),k)
        enddo
        albsoilsn(cmap(pind(nb,1):pind(nb,2)),1)=albsoilsn(cmap(pind(nb,1):pind(nb,2)),1) &
                                                +sv(pind(nb,1):pind(nb,2))*ssoil%albsoilsn(pind(nb,1):pind(nb,2),1)
        albsoilsn(cmap(pind(nb,1):pind(nb,2)),2)=albsoilsn(cmap(pind(nb,1):pind(nb,2)),2) &
                                                 +sv(pind(nb,1):pind(nb,2))*ssoil%albsoilsn(pind(nb,1):pind(nb,2),2)
        albsav(cmap(pind(nb,1):pind(nb,2)))=albsav(cmap(pind(nb,1):pind(nb,2))) &
                                            +sv(pind(nb,1):pind(nb,2))*rad%albedo(pind(nb,1):pind(nb,2),1)
        albnirsav(cmap(pind(nb,1):pind(nb,2)))=albnirsav(cmap(pind(nb,1):pind(nb,2))) &
                                               +sv(pind(nb,1):pind(nb,2))*rad%albedo(pind(nb,1):pind(nb,2),2)
        runoff(cmap(pind(nb,1):pind(nb,2)))=runoff(cmap(pind(nb,1):pind(nb,2))) &
                                            +sv(pind(nb,1):pind(nb,2))*ssoil%runoff(pind(nb,1):pind(nb,2))
        rnof1(cmap(pind(nb,1):pind(nb,2)))=rnof1(cmap(pind(nb,1):pind(nb,2))) &
                                           +sv(pind(nb,1):pind(nb,2))*ssoil%rnof1(pind(nb,1):pind(nb,2))
        rnof2(cmap(pind(nb,1):pind(nb,2)))=rnof2(cmap(pind(nb,1):pind(nb,2))) &
                                           +sv(pind(nb,1):pind(nb,2))*ssoil%rnof2(pind(nb,1):pind(nb,2))
        wbtot(cmap(pind(nb,1):pind(nb,2)))=wbtot(cmap(pind(nb,1):pind(nb,2))) &
                                           +sv(pind(nb,1):pind(nb,2))*ssoil%wbtot(pind(nb,1):pind(nb,2))
        tevap(cmap(pind(nb,1):pind(nb,2)))=tevap(cmap(pind(nb,1):pind(nb,2))) &
                                           +sv(pind(nb,1):pind(nb,2))*bal%evap_tot(pind(nb,1):pind(nb,2))
        tprecip(cmap(pind(nb,1):pind(nb,2)))=tprecip(cmap(pind(nb,1):pind(nb,2))) &
                                             +sv(pind(nb,1):pind(nb,2))*bal%precip_tot(pind(nb,1):pind(nb,2))
        totenbal(cmap(pind(nb,1):pind(nb,2)))=totenbal(cmap(pind(nb,1):pind(nb,2))) &
                                              +sv(pind(nb,1):pind(nb,2))*bal%ebal_tot(pind(nb,1):pind(nb,2))
        trnoff(cmap(pind(nb,1):pind(nb,2)))=trnoff(cmap(pind(nb,1):pind(nb,2))) &
                                            +sv(pind(nb,1):pind(nb,2))*bal%rnoff_tot(pind(nb,1):pind(nb,2))
        tmps(cmap(pind(nb,1):pind(nb,2)))=tmps(cmap(pind(nb,1):pind(nb,2))) &
                                           +sv(pind(nb,1):pind(nb,2))*real(ssoil%isflag(pind(nb,1):pind(nb,2)))
        rtsoil(cmap(pind(nb,1):pind(nb,2)))=rtsoil(cmap(pind(nb,1):pind(nb,2))) &
                                            +sv(pind(nb,1):pind(nb,2))/ssoil%rtsoil(pind(nb,1):pind(nb,2))
        rnet(cmap(pind(nb,1):pind(nb,2)))=rnet(cmap(pind(nb,1):pind(nb,2))) &
                                          +sv(pind(nb,1):pind(nb,2))*canopy%rnet(pind(nb,1):pind(nb,2))
        fg(cmap(pind(nb,1):pind(nb,2)))=fg(cmap(pind(nb,1):pind(nb,2))) &
                                        +sv(pind(nb,1):pind(nb,2))*canopy%fh(pind(nb,1):pind(nb,2))
        eg(cmap(pind(nb,1):pind(nb,2)))=eg(cmap(pind(nb,1):pind(nb,2))) &
                                        +sv(pind(nb,1):pind(nb,2))*canopy%fe(pind(nb,1):pind(nb,2))
        epot(cmap(pind(nb,1):pind(nb,2)))=epot(cmap(pind(nb,1):pind(nb,2))) &
                                          +sv(pind(nb,1):pind(nb,2))*ssoil%potev(pind(nb,1):pind(nb,2))
        tss(cmap(pind(nb,1):pind(nb,2)))=tss(cmap(pind(nb,1):pind(nb,2))) &
                                         +sv(pind(nb,1):pind(nb,2))*rad%trad(pind(nb,1):pind(nb,2))
        tgf(cmap(pind(nb,1):pind(nb,2)))=tgf(cmap(pind(nb,1):pind(nb,2))) &
                                         +sv(pind(nb,1):pind(nb,2))*canopy%tv(pind(nb,1):pind(nb,2))
        cansto(cmap(pind(nb,1):pind(nb,2)))=cansto(cmap(pind(nb,1):pind(nb,2))) &
                                            +sv(pind(nb,1):pind(nb,2))*canopy%cansto(pind(nb,1):pind(nb,2))
        gflux(cmap(pind(nb,1):pind(nb,2)))=gflux(cmap(pind(nb,1):pind(nb,2))) &
                                           +sv(pind(nb,1):pind(nb,2))*canopy%ghflux(pind(nb,1):pind(nb,2))
        sgflux(cmap(pind(nb,1):pind(nb,2)))=sgflux(cmap(pind(nb,1):pind(nb,2))) &
                                            +sv(pind(nb,1):pind(nb,2))*canopy%sghflux(pind(nb,1):pind(nb,2))
        fnee(cmap(pind(nb,1):pind(nb,2)))=fnee(cmap(pind(nb,1):pind(nb,2))) &
                                          +sv(pind(nb,1):pind(nb,2))*canopy%fnee(pind(nb,1):pind(nb,2))
        fpn(cmap(pind(nb,1):pind(nb,2)))=fpn(cmap(pind(nb,1):pind(nb,2))) &
                                         +sv(pind(nb,1):pind(nb,2))*canopy%fpn(pind(nb,1):pind(nb,2))
        frd(cmap(pind(nb,1):pind(nb,2)))=frd(cmap(pind(nb,1):pind(nb,2))) &
                                         +sv(pind(nb,1):pind(nb,2))*canopy%frday(pind(nb,1):pind(nb,2))
        frp(cmap(pind(nb,1):pind(nb,2)))=frp(cmap(pind(nb,1):pind(nb,2))) &
                                         +sv(pind(nb,1):pind(nb,2))*canopy%frp(pind(nb,1):pind(nb,2))
        frpw(cmap(pind(nb,1):pind(nb,2)))=frpw(cmap(pind(nb,1):pind(nb,2))) &
                                          +sv(pind(nb,1):pind(nb,2))*canopy%frpw(pind(nb,1):pind(nb,2))
        frs(cmap(pind(nb,1):pind(nb,2)))=frs(cmap(pind(nb,1):pind(nb,2))) &
                                         +sv(pind(nb,1):pind(nb,2))*canopy%frs(pind(nb,1):pind(nb,2))
        sumpn(cmap(pind(nb,1):pind(nb,2)))=sumpn(cmap(pind(nb,1):pind(nb,2))) &
                                           +sv(pind(nb,1):pind(nb,2))*sum_flux%sumpn(pind(nb,1):pind(nb,2))
        sumrp(cmap(pind(nb,1):pind(nb,2)))=sumrp(cmap(pind(nb,1):pind(nb,2))) &
                                           +sv(pind(nb,1):pind(nb,2))*sum_flux%sumrp(pind(nb,1):pind(nb,2))
        sumrpw(cmap(pind(nb,1):pind(nb,2)))=sumrpw(cmap(pind(nb,1):pind(nb,2))) &
                                            +sv(pind(nb,1):pind(nb,2))*sum_flux%sumrpw(pind(nb,1):pind(nb,2))
        sumrpr(cmap(pind(nb,1):pind(nb,2)))=sumrpr(cmap(pind(nb,1):pind(nb,2))) &
                                            +sv(pind(nb,1):pind(nb,2))*sum_flux%sumrpr(pind(nb,1):pind(nb,2))
        sumrs(cmap(pind(nb,1):pind(nb,2)))=sumrs(cmap(pind(nb,1):pind(nb,2))) &
                                           +sv(pind(nb,1):pind(nb,2))*sum_flux%sumrs(pind(nb,1):pind(nb,2))
        sumrd(cmap(pind(nb,1):pind(nb,2)))=sumrd(cmap(pind(nb,1):pind(nb,2))) &
                                           +sv(pind(nb,1):pind(nb,2))*sum_flux%sumrd(pind(nb,1):pind(nb,2))
        zo(cmap(pind(nb,1):pind(nb,2)))=zo(cmap(pind(nb,1):pind(nb,2))) &
                                        +sv(pind(nb,1):pind(nb,2))/log(zmin/rough%z0m(pind(nb,1):pind(nb,2)))**2
      end do
      where (land)
        zo=max(zmin*exp(-sqrt(1./zo)),zobgin)
        rtsoil=1./rtsoil
      end where
      
      where (land.and.tmps.ge.0.5)
        isflag=1
      elsewhere
        isflag=0
      endwhere
      do nb=1,5 ! update snow
        do k=1,3
          where (ssoil%isflag(pind(nb,1):pind(nb,2)).lt.isflag(cmap(pind(nb,1):pind(nb,2))).and.k.eq.1) ! pack 1-layer into 3-layer
            tggsn(cmap(pind(nb,1):pind(nb,2)),k)=tggsn(cmap(pind(nb,1):pind(nb,2)),k) &
                                               +sv(pind(nb,1):pind(nb,2))*ssoil%tgg(pind(nb,1):pind(nb,2),1)
            smass(cmap(pind(nb,1):pind(nb,2)),k)=smass(cmap(pind(nb,1):pind(nb,2)),k) &
                                               +sv(pind(nb,1):pind(nb,2))*0.05*ssoil%ssdn(pind(nb,1):pind(nb,2),1)
            ssdn(cmap(pind(nb,1):pind(nb,2)),k)=ssdn(cmap(pind(nb,1):pind(nb,2)),k) &
                                              +sv(pind(nb,1):pind(nb,2))*ssoil%ssdn(pind(nb,1):pind(nb,2),1)
          elsewhere (ssoil%isflag(pind(nb,1):pind(nb,2)).lt.isflag(cmap(pind(nb,1):pind(nb,2))).and.k.eq.2)
            tggsn(cmap(pind(nb,1):pind(nb,2)),k)=tggsn(cmap(pind(nb,1):pind(nb,2)),k) &
                                               +sv(pind(nb,1):pind(nb,2))*ssoil%tgg(pind(nb,1):pind(nb,2),1)
            smass(cmap(pind(nb,1):pind(nb,2)),k)=smass(cmap(pind(nb,1):pind(nb,2)),k) &
                                               +sv(pind(nb,1):pind(nb,2))*(ssoil%snowd(pind(nb,1):pind(nb,2)) &
					       -0.05*ssoil%ssdn(pind(nb,1):pind(nb,2),1))*0.4
            ssdn(cmap(pind(nb,1):pind(nb,2)),k)=ssdn(cmap(pind(nb,1):pind(nb,2)),k) &
                                              +sv(pind(nb,1):pind(nb,2))*ssoil%ssdn(pind(nb,1):pind(nb,2),1)          
          elsewhere (ssoil%isflag(pind(nb,1):pind(nb,2)).lt.isflag(cmap(pind(nb,1):pind(nb,2))).and.k.eq.3)
            tggsn(cmap(pind(nb,1):pind(nb,2)),k)=tggsn(cmap(pind(nb,1):pind(nb,2)),k) &
                                               +sv(pind(nb,1):pind(nb,2))*ssoil%tgg(pind(nb,1):pind(nb,2),1)
            smass(cmap(pind(nb,1):pind(nb,2)),k)=smass(cmap(pind(nb,1):pind(nb,2)),k) &
                                               +sv(pind(nb,1):pind(nb,2))*(ssoil%snowd(pind(nb,1):pind(nb,2)) &
					       -0.05*ssoil%ssdn(pind(nb,1):pind(nb,2),1))*0.6
            ssdn(cmap(pind(nb,1):pind(nb,2)),k)=ssdn(cmap(pind(nb,1):pind(nb,2)),k) &
                                              +sv(pind(nb,1):pind(nb,2))*ssoil%ssdn(pind(nb,1):pind(nb,2),1)
          elsewhere (ssoil%isflag(pind(nb,1):pind(nb,2)).gt.isflag(cmap(pind(nb,1):pind(nb,2))).and.k.eq.1) ! pack 3-layer into 1-layer
            tggsn(cmap(pind(nb,1):pind(nb,2)),k)=tggsn(cmap(pind(nb,1):pind(nb,2)),k) &
                                               +sv(pind(nb,1):pind(nb,2))*273.16
            smass(cmap(pind(nb,1):pind(nb,2)),k)=smass(cmap(pind(nb,1):pind(nb,2)),k) &
                                               +sv(pind(nb,1):pind(nb,2))*ssoil%snowd(pind(nb,1):pind(nb,2))
            ssdn(cmap(pind(nb,1):pind(nb,2)),k)=ssdn(cmap(pind(nb,1):pind(nb,2)),k) &
                                              +sv(pind(nb,1):pind(nb,2))*ssoil%ssdnn(pind(nb,1):pind(nb,2))
          elsewhere (ssoil%isflag(pind(nb,1):pind(nb,2)).gt.isflag(cmap(pind(nb,1):pind(nb,2))).and.k.ge.2)
            tggsn(cmap(pind(nb,1):pind(nb,2)),k)=tggsn(cmap(pind(nb,1):pind(nb,2)),k) &
                                               +sv(pind(nb,1):pind(nb,2))*273.16
            ssdn(cmap(pind(nb,1):pind(nb,2)),k)=ssdn(cmap(pind(nb,1):pind(nb,2)),k) &
                                              +sv(pind(nb,1):pind(nb,2))*ssoil%ssdn(pind(nb,1):pind(nb,2),k)          
          elsewhere
            tggsn(cmap(pind(nb,1):pind(nb,2)),k)=tggsn(cmap(pind(nb,1):pind(nb,2)),k) &
                                               +sv(pind(nb,1):pind(nb,2))*ssoil%tggsn(pind(nb,1):pind(nb,2),k)
            smass(cmap(pind(nb,1):pind(nb,2)),k)=smass(cmap(pind(nb,1):pind(nb,2)),k) &
                                               +sv(pind(nb,1):pind(nb,2))*ssoil%smass(pind(nb,1):pind(nb,2),k)
            ssdn(cmap(pind(nb,1):pind(nb,2)),k)=ssdn(cmap(pind(nb,1):pind(nb,2)),k) &
                                              +sv(pind(nb,1):pind(nb,2))*ssoil%ssdn(pind(nb,1):pind(nb,2),k)
          end where
        end do
        ssdnn(cmap(pind(nb,1):pind(nb,2)))=ssdnn(cmap(pind(nb,1):pind(nb,2))) &
                                         +sv(pind(nb,1):pind(nb,2))*ssoil%ssdnn(pind(nb,1):pind(nb,2))
        snage(cmap(pind(nb,1):pind(nb,2)))=snage(cmap(pind(nb,1):pind(nb,2))) &
                                         +sv(pind(nb,1):pind(nb,2))*ssoil%snage(pind(nb,1):pind(nb,2))
        snowd(cmap(pind(nb,1):pind(nb,2)))=snowd(cmap(pind(nb,1):pind(nb,2))) &
                                         +sv(pind(nb,1):pind(nb,2))*ssoil%snowd(pind(nb,1):pind(nb,2))
        !osnowd(cmap(pind(nb,1):pind(nb,2)))=osnowd(cmap(pind(nb,1):pind(nb,2))) &
        !                                  +sv(pind(nb,1):pind(nb,2))*ssoil%osnowd(pind(nb,1):pind(nb,2))
        !osnowd0(cmap(pind(nb,1):pind(nb,2)))=osnowd0(cmap(pind(nb,1):pind(nb,2))) &
        !                                   +sv(pind(nb,1):pind(nb,2))*bal%osnowd0(pind(nb,1):pind(nb,2))
      end do

      return
      end subroutine sib4

! *************************************************************************************
      subroutine setco2for(jyear)
!     set co2 forcing for cable
!     constant: atmospheric co2 = 360 ppm 
!     prescribed: atmospheric co2 follows prescribed trend from 1900-2004
!                 based on ice core and Mauna Loa/South Pole data
!     interactive: atmospheric co2 taken from tracer (usually cable+fos+ocean)

      implicit none

      include 'newmpar.h'
      include 'tracers.h'

      integer, intent(in) :: jyear
      integer, parameter :: constantCO2 = 1
      integer, parameter :: prescribedCO2 = 2
      integer, parameter :: interactiveCO2 = 3
      integer ipco2

!     rml added values for 2000-2004 using scripps records from cdiac
!     0.75*mauna loa + 0.25*south pole
!     2005 from gv07, noaa flask record
      real co2for(0:105)
      data co2for/296.0049,296.3785,296.7731,297.1795,297.5887,297.9919 &
     &           ,298.3842,298.7654,299.1349,299.4925,299.838 ,300.1709 &
     &           ,300.491 ,300.801 ,301.106 ,301.4113,301.7205,302.0357 &
     &           ,302.3587,302.6915,303.036,303.3941,303.7677,304.1587 &
     &           ,304.569,304.9971,305.4388,305.8894,306.3444,306.7992 &
     &           ,307.2494,307.6902,308.117,308.521,308.8895,309.2135 &
     &           ,309.4877,309.7068,309.8658,309.9668,310.019,310.0358 &
     &           ,310.035,310.0345,310.0522,310.0989,310.1794,310.2977 &
     &           ,310.4581,310.6661,310.928,311.2503,311.6395,312.1015 &
     &           ,312.6341,313.227,313.8694,314.5506,315.2599,315.9866 &
     &           ,316.7167,317.4268,318.106,318.7638,319.4402,320.1901 &
     &           ,321.0398,321.9713,322.9779,324.045,325.1356,326.2445 &
     &           ,327.3954,328.5711,329.723,330.8865,332.1331,333.5012 &
     &           ,334.9617,336.4614,337.9588,339.4225,340.8724,342.3488 &
     &           ,343.8625,345.421,347.0481,348.7452,350.4397,352.0169 &
     &           ,353.4269,354.6917,355.8849,357.1348,358.5514,360.1441 &
     &           ,361.8578,363.6479,365.4682,367.2137 &
     &           ,368.87,370.35,372.49,374.93,376.69,379.09/

      select case (CO2forcingtype)
        case (constantCO2); atmco2 = 360.
        case (prescribedCO2) 
          ipco2 = jyear-1900
          if (ipco2.lt.0.or.ipco2.gt.105) stop 'year out of range for&
     & co2 forcing for cable'
          atmco2 = co2for(ipco2)
        case (interactiveCO2)
          write(6,*) 'need to replace with tracer sum'
          atmco2 = tr(1:ifull,1,1)
      end select

      return
      end subroutine setco2for

! *************************************************************************************
  subroutine loadcbmparm(fveg)

  use cc_mpi
  use define_dimensions, only : mp
  
  implicit none
  
  include 'newmpar.h'
  include 'const_phys.h'
  include 'carbpools.h'
  include 'latlong.h'  
  include 'nsibd.h'
  include 'parm.h'
  include 'parmgeom.h'  ! rlong0,rlat0,schmidt
  include 'sigs.h'
  include 'soil.h'
  include 'soilsnow.h'
  include 'soilv.h'
  include 'vegpar.h' !   

  integer(i_d), dimension(ifull_g,5) :: ivsg
  integer(i_d), dimension(ifull,5) :: ivs
  integer(i_d) iq,n,nb,k,iad,ipos,ip
  real(r_1) :: totdepth,ra,rb
  real(r_1), dimension(mxvt,ms) :: froot2
  real(r_1), dimension(ifull_g,5) :: svsg,vling
  real(r_1), dimension(ifull,5) :: svs,vlin
  character(len=*), intent(in) :: fveg
  character(len=10), dimension(mxvt) :: vegtype
  integer ilx,jlx
  real rlong0x,rlat0x,schmidtx,dsx
  character*47 header

  if (myid == 0) print *,"Setting CABLE defaults (igbp)"
  
  zmin=-rdry*280.*log(sig(1))/grav ! not yet defined in indata.f

  hc=(/ 17.,35.,15.5,20.,19.25,0.6,0.6,7.0426,14.3379,0.567,0.5,0.55,6.017,0.55,0.2,0.2,0.2 /)
  xfang=(/ 0.01,0.1,0.01,0.25,0.125,-0.3,0.01,-0.3,-0.3,-0.3,0.,-0.3,0.,-0.3,0.,0.1,0. /)
  dleaf=(/ 0.055,0.1,0.04,0.15,0.1,0.1,0.1,0.233,0.129,0.3,0.3,0.3,0.242,0.3,0.03,0.03,0.03 /)
  xalbnir=(/ 0.79,0.96,0.81,1.,1.08,1.14,1.2,1.02,1.23,1.16,0.89,0.98,1.1,1.13,1.,1.15,1. /)
  wai=(/ 1.,1.,1.,1.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0. /)
  canst1=0.1
  shelrb=2.
  vegcf=(/ 0.91,1.95,0.73,1.5,1.55,0.6,2.05,2.8,2.75,2.75,0.,2.8,0.,2.8,0.,0.4,0. /)
  vegtype(1)='forest'
  vegtype(2)='forest'
  vegtype(3)='deciduous'
  vegtype(4)='deciduous'
  vegtype(5)='forest'
  vegtype(6)='not used'
  vegtype(7)='shrub'
  vegtype(8)='shrub'
  vegtype(9)='shrub'
  vegtype(10)='grass'
  vegtype(11)='not used'
  vegtype(12)='crop'
  vegtype(13)='not used'
  vegtype(14)='not used'
  vegtype(15)='no veg'
  vegtype(16)='grass'
  vegtype(17)='not used'
  extkn=0.7
  vcmax=(/ 65.2E-6,65.E-6,70.E-6,85.0E-6,80.E-6,20.E-6,20.E-6,10.E-6,20.E-6,10.E-6,50.E-6,80.E-6,1.E-6,80.E-6, &
           17.E-6,17.E-6,17.E-6 /)
  ejmax=2.*vcmax
  rp20=(/ 3.3039,1.1342,3.2879,1.4733,2.3704,5.,5.,1.05,2.,0.8037,1.,3.,1.,3.,0.1,0.1,0.1 /)
  rpcoef=0.0832
  rs20=(/ 1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,0.,0. /)
  tminvj=(/ -15.,-15.,5.,5.,5.,-15.,-15.,-15.,-15.,-15.,-15.,-15.,-15.,-15.,-15.,-15.,-15. /)
  tmaxvj=(/ -10.,-10.,10.,15.,10.,-10.,-10.,-10.,-10.,-10.,-10.,-10.,-10.,-10.,-10.,-10.,-10. /)
  vbeta=(/ 2.,2.,2.,2.,2.,4.,4.,4.,4.,4.,2.,2.,2.,2.,4.,4.,4. /)
  rootbeta=(/ 0.943,0.962,0.966,0.961,0.966,0.914,0.964,0.972,0.943,0.943,0.961,0.961,0.961,0.961,0.961,0.975,0.975 /)
  tcplant(:,1)=(/ 200.  ,300.  ,200.  ,300.  ,200.  ,150. ,150. ,250. ,250. ,250. ,250.,150.,0.1,150.,0.,1.,0. /)
  tcplant(:,2)=(/ 10217.,16833.,5967. ,12000.,10217.,5000.,5000.,5247.,5247.,0.   ,0.  ,0.  ,0. ,0.  ,0.,0.,0. /)
  tcplant(:,3)=(/ 876.  ,1443. ,511.  ,1029. ,876.  ,500. ,500. ,1124.,1124.,500. ,500.,607.,0.1,607.,0.,1.,0. /)
  tcsoil(:,1)=(/ 184.   ,303.  ,107.  ,216.  ,184.  ,100. ,100. ,275. ,275. ,275. ,275.,149.,0.1,149.,0.,1.,0. /)
  tcsoil(:,2)=(/ 367.   ,606.  ,214.  ,432.  ,367.  ,250. ,250. ,314. ,314. ,314. ,314.,300.,0.1,300.,0.,1.,0. /)
  ratecp(1)=1.
  ratecp(2)=0.03
  ratecp(3)=0.14
  ratecs(1)=2.
  ratecs(2)=0.5

  if (cbm_ms.ne.ms) then
    print *,"ERROR: CABLE and CCAM soil levels do not match"
    stop
  end if

  if (myid==0) then
    print *,"Reading land-use data for CABLE"
    open(87,file=fveg,status='old')
    read(87,*) ilx,jlx,rlong0x,rlat0x,schmidtx,dsx,header
    if(ilx.ne.il_g.or.jlx.ne.jl_g.or.rlong0x.ne.rlong0.or.rlat0x.ne.rlat0.or.schmidtx.ne.schmidt) stop 'wrong data file supplied'
    do iq=1,ifull_g
      read(87,*) iad,ra,rb,ivsg(iq,1),svsg(iq,1),vling(iq,1),ivsg(iq,2),svsg(iq,2),vling(iq,2),ivsg(iq,3),svsg(iq,3),vling(iq,3), &
                 ivsg(iq,4),svsg(iq,4),vling(iq,4),ivsg(iq,5),svsg(iq,5),vling(iq,5)
    end do
    close(87)
    do n=1,5
      call ccmpi_distribute(ivs(:,n),ivsg(:,n))
      call ccmpi_distribute(svs(:,n),svsg(:,n))
      call ccmpi_distribute(vlin(:,n),vling(:,n))
    end do
  else
    do n=1,5
      call ccmpi_distribute(ivs(:,n))
      call ccmpi_distribute(svs(:,n))
      call ccmpi_distribute(vlin(:,n))
    end do
  end if
  do n=1,5
    svs(:,n)=svs(:,n)/sum(svs,2)
  end do

   mp=0
   do iq=1,ifull
     if (land(iq)) then
       mp=mp+count(svs(iq,:).gt.0.)
     end if
   end do
   nb=count(land)
  
  allocate(sv(mp))
  allocate(vl(mp))
  allocate(cmap(mp))
  allocate(atmco2(ifull))
  call alloc_cbm_var(air, mp)
  call alloc_cbm_var(bgc, mp)
  call alloc_cbm_var(canopy, mp)
  call alloc_cbm_var(met, mp)
  call alloc_cbm_var(bal, mp)
  call alloc_cbm_var(rad, mp)
  call alloc_cbm_var(rough, mp)
  call alloc_cbm_var(soil, mp)
  call alloc_cbm_var(ssoil, mp)
  call alloc_cbm_var(sum_flux, mp)
  call alloc_cbm_var(veg, mp)
  
  ! soil parameters
  soil%zse = (/.022, .058, .154, .409, 1.085, 2.872/) ! soil layer thickness
  soil%zshh(1) = 0.5 * soil%zse(1)
  soil%zshh(ms+1) = 0.5 * soil%zse(ms)
  soil%zshh(2:ms) = 0.5 * (soil%zse(1:ms-1) + soil%zse(2:ms))
  
! preferred option
! froot is now calculated from soil depth and the new parameter rootbeta 
! according to Jackson et al. 1996, Oceologica, 108:389-411
  totdepth = 0.0
  do k=1,ms
    totdepth = totdepth + soil%zse(k)*100.0
    froot2(:,k) = min(1.0,1.0-rootbeta(:)**totdepth)
  enddo
  do k = ms, 2, -1
    froot2(:,k) = froot2(:,k) - froot2(:,k-1)
  enddo
  
  ipos=0
  do n=1,5
    ip=0
    call getc4(ifull,ivs(:,n),rlatt(:)*180./pi,c4frac(:))
    pind(n,1)=ipos+1
    do iq=1,ifull
      if (land(iq).and.(svs(iq,n).gt.0.)) then
        ipos=ipos+1
        if (ivs(iq,n).lt.1) then
          print *,"ERROR: Land-type/lsmask mismatch at iq=",iq
          stop
        end if
        cmap(ipos)=iq
        sv(ipos)=svs(iq,n)
        vl(ipos)=vlin(iq,n)
        veg%iveg(ipos)=ivs(iq,n)
        soil%isoilm(ipos)=isoilm(iq)
        veg%frac4(ipos)=c4frac(iq)
      end if
    end do
    pind(n,2)=ipos
  end do
  
  if (ipos.ne.mp) then
    print *,"ERROR: Internal memory allocation error for CABLE set-up"
    stop
  end if
  
  do iq=1,ifull
    vlai(iq)=dot_product(vlin(iq,:),svs(iq,:))
  end do
  
  ! aggregate zom
  zolnd=0.
  do n=1,5
    where (land)
      zolnd=zolnd+svs(:,n)/log(zmin/(0.1*hc(ivs(:,n))))**2
    end where
  end do
  where (land)
    zolnd=max(zmin*exp(-sqrt(1./zolnd)),zobgin)
  elsewhere
    zolnd=zobgin
  end where
  
  ! use dominant veg type
  ivegt=ivs(:,1)
  veg%hc     = hc(veg%iveg)
  veg%canst1 = canst1(veg%iveg)
  veg%ejmax  = ejmax(veg%iveg)
  veg%tminvj = tminvj(veg%iveg)
  veg%tmaxvj = tmaxvj(veg%iveg)
  veg%vbeta  = vbeta(veg%iveg)
  veg%rp20   = rp20(veg%iveg)
  veg%rpcoef = rpcoef(veg%iveg)
  veg%shelrb = shelrb(veg%iveg)
  veg%vcmax  = vcmax(veg%iveg)
  veg%xfang  = xfang(veg%iveg)
  veg%dleaf  = dleaf(veg%iveg)
  veg%wai    = wai(veg%iveg)
  veg%vegcf  = vegcf(veg%iveg)
  veg%extkn  = extkn(veg%iveg)
  veg%xalbnir= xalbnir(veg%iveg)
  soil%rs20  = rs20(veg%iveg)
  do k=1,ms
    veg%froot(:,k)=froot2(veg%iveg,k)
  end do
  veg%deciduous =(vegtype(veg%iveg).eq.'deciduous')
  
  if (all(cplant.eq.0.)) then
    if (myid == 0) print *,"Using default carbpools"
    do k=1,ncp
      where (land)
        cplant(:,ncp)=tcplant(ivegt,ncp)
      end where
    end do
    do k=1,ncs
      where (land)        
        csoil(:,ncs)=tcsoil(ivegt,ncs)
      end where
    end do
  else
    if (myid == 0) print *,"Loading carbpools from ifile"
  end if

  soil%bch     = bch(soil%isoilm)
  soil%css     = css(soil%isoilm)
  soil%rhosoil = rhos(soil%isoilm)
  soil%cnsd    = cnsd(soil%isoilm)
  soil%hyds    = hyds(soil%isoilm)
  soil%sucs    = sucs(soil%isoilm)
  soil%hsbh    = hsbh(soil%isoilm)
  soil%sfc     = sfc(soil%isoilm)
  soil%ssat    = ssat(soil%isoilm)
  soil%swilt   = swilt(soil%isoilm)
  soil%ibp2    = ibp2(soil%isoilm)
  soil%i2bp3   = i2bp3(soil%isoilm)
  
  bgc%ratecp(:) = ratecp(:)
  bgc%ratecs(:) = ratecs(:)
  
  sigmf=(1.-exp(-extkn(ivegt(:))*vlai(:)))

  ! store bare soil albedo and define snow free albedo
  albsoilsn(:,:)=0.08
  where (land)
    albsoil(:)=0.5*sum(albvisnir,2)
  elsewhere   
    albsoil(:)=0.08
  end where
  where ((albsoil.le.0.14).and.land)
    !sfact=0.5 for alb <= 0.14 (see cable_soilsnow.f90)
    albsoilsn(:,1)=(1.00/1.50)*albsoil(:)
    albsoilsn(:,2)=(2.00/1.50)*albsoil(:)
  elsewhere ((albsoil(:).le.0.2).and.land)
    !sfact=0.62 for 0.14 < alb <= 0.20 (see cable_soilsnow.f90)
    albsoilsn(:,1)=(1.24/1.62)*albsoil(:)
    albsoilsn(:,2)=(2.00/1.62)*albsoil(:)
  elsewhere (land)
    !sfact=0.68 for 0.2 < alb (see cable_soilsnow.f90)
    albsoilsn(:,1)=(1.36/1.68)*albsoil(:)
    albsoilsn(:,2)=(2.00/1.68)*albsoil(:)
  end where
 ! MJT suggestion to get an approx inital albedo (before cable is called)
  where (land)
    albvisnir(:,1)=albsoilsn(:,1)*(1.-sigmf)+0.1*sigmf
    albvisnir(:,2)=albsoilsn(:,2)*(1.-sigmf)+0.3*sigmf
  elsewhere
    albvisnir(:,1)=albsoilsn(:,1)
    albvisnir(:,2)=albsoilsn(:,2)
  end where

  soil%albsoil=albsoil(cmap)
  ssoil%albsoilsn(:,1)=albsoilsn(cmap,1) ! overwritten by CABLE
  ssoil%albsoilsn(:,2)=albsoilsn(cmap,2) ! overwritten by CABLE
  ssoil%albsoilsn(:,3)=0.05
  
  ssoil%rtsoil=0. ! either load from tile or define in sib4
    
  rad%latitude=rlatt(cmap)*180./pi
  
  gflux=0. ! MJT suggestion
  sgflux=0. ! MJT suggestion
  canopy%ghflux=0.
  canopy%sghflux=0.  
  
 ! Initialise sum flux variables
  sumpn = 0.
  sumrp = 0.
  sumrs = 0.
  sumrd = 0.
  sum_flux%sumpn=0.
  sum_flux%sumrp=0.
  sum_flux%sumrs=0.
  sum_flux%sumrd=0.
  sum_flux%sumrpw=0.
  sum_flux%sumrpr=0.
  sum_flux%dsumpn=0.
  sum_flux%dsumrp=0.
  sum_flux%dsumrs=0.
  sum_flux%dsumrd=0.
  
  bal%evap_tot=0.
  bal%precip_tot=0.
  bal%ebal_tot=0.
  bal%rnoff_tot=0.

  call loadtile ! load tgg,wb,wbice,snowd,snage,tggsn,smass,ssdn,isflag,rtsoil,cansto,cplant and csoil

  ssoil%osnowd=0.                                         ! overwritten by CABLE
  bal%osnowd0=0.                                          ! overwritten by CABLE
  ssoil%ssdnn=120.                                        ! overwritten by CABLE
  where (ssoil%isflag.gt.0)
    ssoil%sdepth(:,1)=ssoil%smass(:,1)/ssoil%ssdn(:,1)    ! overwritten by CABLE
  elsewhere
    ssoil%sdepth(:,1)=ssoil%snowd/ssoil%ssdn(:,1)         ! overwritten by CABLE
  end where
  do k=2,3
    where (ssoil%isflag.gt.0)
      ssoil%sdepth(:,k)=ssoil%smass(:,k)/ssoil%ssdn(:,k)  ! overwritten by CABLE
    elsewhere
      ssoil%sdepth(:,k)=0.                                ! overwritten by CABLE
    end where
  end do  
  
  return
  end subroutine loadcbmparm

! *************************************************************************************
  subroutine getc4(ifull,ivegt,rlatt,c4frac)
  
  implicit none
  
  integer, intent(in) :: ifull
  integer, dimension(ifull), intent(in) :: ivegt
  real, dimension(ifull), intent(in) :: rlatt
  real, dimension(ifull), intent(out) :: c4frac

  where ((ivegt.eq.7).and.(rlatt.ge.-30.).and.(rlatt.le.30.))
    c4frac=0.95
  elsewhere ((ivegt.eq.8).and.(rlatt.ge.-30.).and.(rlatt.le.0.))
    c4frac=0.5
  elsewhere ((ivegt.eq.8).and.(rlatt.ge.0.).and.(rlatt.le.20.))
    c4frac=0.8
  elsewhere ((ivegt.eq.8).and.(rlatt.ge.20.).and.(rlatt.le.30.))
    c4frac=0.5
  elsewhere (ivegt.eq.9)
    c4frac=0.75
  elsewhere ((ivegt.eq.10).and.(rlatt.ge.-30.).and.(rlatt.le.-20.))
    c4frac=0.5
  elsewhere ((ivegt.eq.10).and.(rlatt.ge.-20.).and.(rlatt.le.20.))
    c4frac=0.95
  elsewhere ((ivegt.eq.10).and.(rlatt.ge.20.).and.(rlatt.le.35.))
    c4frac=0.5
  elsewhere ((ivegt.eq.12).and.(rlatt.ge.0.).and.(rlatt.le.40.))
    c4frac=0.3
  elsewhere
    c4frac=0.
  end where
  
  return
  end subroutine getc4

! *************************************************************************************  
  subroutine loadtile

  use cc_mpi, only : myid
  
  implicit none

  include 'newmpar.h'
  include 'darcdf.h'
  include 'parm.h'
  include 'netcdf.inc'
  include 'mpif.h'
  include 'carbpools.h'
  include 'soilsnow.h'
  include 'vegpar.h'    
  
  integer k,n,iarchi,ierr,ierr2,idv
  real, dimension(ifull) :: dat
  character*9 vname
  real sigin  
  integer ik,jk,kk
  common/sigin/ik,jk,kk,sigin(40)  ! for vertint, infile ! MJT bug  

  iarchi=1 ! assume restart file

  if (io_in.eq.1) then
    if (myid.eq.0) idv = ncvid(ncid,"tgg1_1",ierr)
    call MPI_Bcast(ierr,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr2)    
  else
    ierr=1
  end if
  
  if (ierr.ne.0) then
    if (myid==0) write(6,*) "Use averaged data to initialise CABLE"
    do k = 1,ms
      ssoil%tgg(:,k) = tgg(cmap,k)
      ssoil%wb(:,k) = wb(cmap,k)
      ssoil%wbice(:,k) = wbice(cmap,k)
    enddo
    do k = 1,3
      ssoil%tggsn(:,k) = tggsn(cmap,k)
      ssoil%smass(:,k) = smass(cmap,k)
      ssoil%ssdn(:,k) = ssdn(cmap,k)
    enddo      
    ssoil%isflag=isflag(cmap)
    ssoil%snowd=snowd(cmap)
    ssoil%snage=snage(cmap)
    ssoil%rtsoil=0.
    canopy%cansto=0.
    do k=1,ncp
      bgc%cplant(:,k) = cplant(cmap,k)
    enddo
    do k=1,ncs
      bgc%csoil(:,k) = csoil(cmap,k)
    enddo  
  else
    if (myid==0) write(6,*) "Use tiled data to initialise CABLE"
    do n=1,5
      do k=1,ms
        write(vname,'("tgg",I1.1,"_",I1.1)') k,n
        call histrd1(ncid,iarchi,ierr,vname,ik,jk,dat,ifull)
        ssoil%tgg(pind(n,1):pind(n,2),k)=dat(cmap(pind(n,1):pind(n,2)))
        write(vname,'("wb",I1.1,"_",I1.1)') k,n
        call histrd1(ncid,iarchi,ierr,vname,ik,jk,dat,ifull)
        ssoil%wb(pind(n,1):pind(n,2),k)=dat(cmap(pind(n,1):pind(n,2)))
        write(vname,'("wbice",I1.1,"_",I1.1)') k,n
        call histrd1(ncid,iarchi,ierr,vname,ik,jk,dat,ifull)
        ssoil%wbice(pind(n,1):pind(n,2),k)=dat(cmap(pind(n,1):pind(n,2)))
      end do
      do k=1,3
        write(vname,'("tggsn",I1.1,"_",I1.1)') k,n
        call histrd1(ncid,iarchi,ierr,vname,ik,jk,dat,ifull)
        ssoil%tggsn(pind(n,1):pind(n,2),k)=dat(cmap(pind(n,1):pind(n,2)))
        write(vname,'("smass",I1.1,"_",I1.1)') k,n
        call histrd1(ncid,iarchi,ierr,vname,ik,jk,dat,ifull)
        ssoil%smass(pind(n,1):pind(n,2),k)=dat(cmap(pind(n,1):pind(n,2)))
        write(vname,'("ssdn",I1.1,"_",I1.1)') k,n
        call histrd1(ncid,iarchi,ierr,vname,ik,jk,dat,ifull)
        ssoil%ssdn(pind(n,1):pind(n,2),k)=dat(cmap(pind(n,1):pind(n,2)))
      end do
      write(vname,'("sflag_",I1.1)') n
      call histrd1(ncid,iarchi,ierr,vname,ik,jk,dat,ifull)
      ssoil%isflag(pind(n,1):pind(n,2))=int(dat(cmap(pind(n,1):pind(n,2))))
      write(vname,'("snd_",I1.1)') n
      call histrd1(ncid,iarchi,ierr,vname,ik,jk,dat,ifull)
      ssoil%snowd(pind(n,1):pind(n,2))=dat(cmap(pind(n,1):pind(n,2)))
      write(vname,'("snage_",I1.1)') n
      call histrd1(ncid,iarchi,ierr,vname,ik,jk,dat,ifull)
      ssoil%snage(pind(n,1):pind(n,2))=dat(cmap(pind(n,1):pind(n,2)))
      write(vname,'("rtsoil_",I1.1)') n
      call histrd1(ncid,iarchi,ierr,vname,ik,jk,dat,ifull)
      ssoil%rtsoil(pind(n,1):pind(n,2))=dat(cmap(pind(n,1):pind(n,2)))
      write(vname,'("cansto_",I1.1)') n
      call histrd1(ncid,iarchi,ierr,vname,ik,jk,dat,ifull)
      canopy%cansto(pind(n,1):pind(n,2))=dat(cmap(pind(n,1):pind(n,2)))
      do k=1,ncp
        write(vname,'("cplant",I1.1,"_",I1.1)') k,n
        call histrd1(ncid,iarchi,ierr,vname,ik,jk,dat,ifull)
        bgc%cplant(pind(n,1):pind(n,2),k) = dat(cmap(pind(n,1):pind(n,2))) 
      enddo
      do k=1,ncs
        write(vname,'("csoil",I1.1,"_",I1.1)') k,n
        call histrd1(ncid,iarchi,ierr,vname,ik,jk,dat,ifull)
        bgc%csoil(pind(n,1):pind(n,2),k) = dat(cmap(pind(n,1):pind(n,2)))
      enddo
    end do
  end if
  
  return
  end subroutine loadtile

 ! *************************************************************************************
  subroutine savetile(idnc,local,idim)

  use cc_mpi, only : myid
  
  implicit none

  include 'newmpar.h'
  include 'carbpools.h'
  include 'soilsnow.h'
  include 'vegpar.h'
  
  integer k,n,iarch,ierr
  integer, intent(in) :: idnc
  integer, dimension(3), intent(in) :: idim  
  real, dimension(ifull) :: dat
  character*9 vname
  character*40 lname
  logical, intent(in) :: local
  real sigin  
  integer ik,jk,kk
  common/sigin/ik,jk,kk,sigin(40)  ! for vertint, infile ! MJT bug    
  
  iarch=1 ! assume restart file
  
  if (myid.eq.0) then
    print *,"Storing CABLE tile data"
    call ncredf(idnc,ierr)
    do n=1,5
      do k=1,ms
        write(lname,'("Soil temperature lev ",I1.1," tile ",I1.1)') k,n
        write(vname,'("tgg",I1.1,"_",I1.1)') k,n
        call attrib(idnc,idim,3,vname,lname,'K',100.,400.,0)
        write(lname,'("Soil moisture lev ",I1.1," tile ",I1.1)') k,n
        write(vname,'("wb",I1.1,"_",I1.1)') k,n 
        call attrib(idnc,idim,3,vname,lname,'m3/m3',0.,1.,0)
        write(lname,'("Soil ice lev ",I1.1," tile ",I1.1)') k,n
        write(vname,'("wbice",I1.1,"_",I1.1)') k,n 
        call attrib(idnc,idim,3,vname,lname,'m3/m3',0.,1.,0)
      end do
      do k=1,3
        write(lname,'("Snow temperature lev ",I1.1," tile ",I1.1)') k,n
        write(vname,'("tggsn",I1.1,"_",I1.1)') k,n
        call attrib(idnc,idim,3,vname,lname,'K',100.,400.,0)
        write(lname,'("Snow mass lev ",I1.1," tile ",I1.1)') k,n
        write(vname,'("smass",I1.1,"_",I1.1)') k,n 
        call attrib(idnc,idim,3,vname,lname,'K',0.,400.,0)
        write(lname,'("Snow density lev ",I1.1," tile ",I1.1)') k,n
        write(vname,'("ssdn",I1.1,"_",I1.1)') k,n 
        call attrib(idnc,idim,3,vname,lname,'K',0.,400.,0)
      end do
      write(lname,'("Snow flag tile ",I1.1)') n
      write(vname,'("sflag_",I1.1)') n
      call attrib(idnc,idim,3,vname,lname,'mm',0.,4.,0)
      write(lname,'("Snow depth tile ",I1.1)') n
      write(vname,'("snd_",I1.1)') n
      call attrl (idnc,idim,3,vname,lname,'mm',0.,5000.,0)  ! long
      write(lname,'("Snow age tile ",I1.1)') n
      write(vname,'("snage_",I1.1)') n
      call attrib(idnc,idim,3,vname,lname,'none',0.,20.,0)
      write(lname,'("Soil turbulent resistance tile ",I1.1)') n
      write(vname,'("rtsoil_",I1.1)') n
      call attrib(idnc,idim,3,vname,lname,'none',0.,9.e4,0)
      write(lname,'("cansto tile ",I1.1)') n
      write(vname,'("cansto_",I1.1)') n
      call attrib(idnc,idim,3,vname,lname,'none',0.,10.,0)
      write(lname,'("Carbon leaf pool tile ",I1.1)') n
      write(vname,'("cplant1_",I1.1)') n    
      call attrib(idnc,idim,3,vname,lname,'none',0.,50000.,0)
      write(lname,'("Carbon wood pool tile ",I1.1)') n
      write(vname,'("cplant2_",I1.1)') n
      call attrib(idnc,idim,3,vname,lname,'none',0.,50000.,0)
      write(lname,'("Carbon root pool tile ",I1.1)') n
      write(vname,'("cplant3_",I1.1)') n
      call attrib(idnc,idim,3,vname,lname,'none',0.,50000.,0)
      write(lname,'("Carbon soil fast pool tile ",I1.1)') n
      write(vname,'("csoil1_",I1.1)') n
      call attrib(idnc,idim,3,vname,lname,'none',0.,50000.,0)
      write(lname,'("Carbon soil slow pool tile ",I1.1)') n
      write(vname,'("csoil2_",I1.1)') n
      call attrib(idnc,idim,3,vname,lname,'none',0.,50000.,0)
    end do      
    call ncendf(idnc,ierr)
  end if
  do n=1,5
    do k=1,ms
      dat=tgg(:,k)
      dat(cmap(pind(1,n):pind(2,n)))=ssoil%tgg(pind(n,1):pind(n,2),k)
      write(vname,'("tgg",I1.1,"_",I1.1)') k,n
      call histwrt3(dat,vname,idnc,iarch,local)
      dat=wb(:,k)
      dat(cmap(pind(1,n):pind(2,n)))=ssoil%wb(pind(n,1):pind(n,2),k)
      write(vname,'("wb",I1.1,"_",I1.1)') k,n
      call histwrt3(dat,vname,idnc,iarch,local)
      dat=wbice(:,k)
      dat(cmap(pind(1,n):pind(2,n)))=ssoil%wbice(pind(n,1):pind(n,2),k)
      write(vname,'("wbice",I1.1,"_",I1.1)') k,n
      call histwrt3(dat,vname,idnc,iarch,local)
    end do
    do k=1,3
      dat=tggsn(:,k)
      dat(cmap(pind(1,n):pind(2,n)))=ssoil%tggsn(pind(n,1):pind(n,2),k)
      write(vname,'("tggsn",I1.1,"_",I1.1)') k,n
      call histwrt3(dat,vname,idnc,iarch,local)
      dat=smass(:,k)
      dat(cmap(pind(1,n):pind(2,n)))=ssoil%smass(pind(n,1):pind(n,2),k)
      write(vname,'("smass",I1.1,"_",I1.1)') k,n
      call histwrt3(dat,vname,idnc,iarch,local)
      dat=ssdn(:,k)
      dat(cmap(pind(1,n):pind(2,n)))=ssoil%ssdn(pind(n,1):pind(n,2),k)
      write(vname,'("ssdn",I1.1,"_",I1.1)') k,n
      call histwrt3(dat,vname,idnc,iarch,local)
    end do
    dat=real(isflag)
    dat(cmap(pind(1,n):pind(2,n)))=real(ssoil%isflag(pind(n,1):pind(n,2)))
    write(vname,'("sflag_",I1.1)') n
    call histwrt3(dat,vname,idnc,iarch,local)
    dat=snowd
    dat(cmap(pind(1,n):pind(2,n)))=ssoil%snowd(pind(n,1):pind(n,2))
    write(vname,'("snd_",I1.1)') n
    call histwrt3l(dat,vname,idnc,iarch,local)  ! long write    
    dat=snage
    dat(cmap(pind(1,n):pind(2,n)))=ssoil%snage(pind(n,1):pind(n,2))
    write(vname,'("snage_",I1.1)') n
    call histwrt3(dat,vname,idnc,iarch,local)
    dat=rtsoil
    dat(cmap(pind(1,n):pind(2,n)))=ssoil%rtsoil(pind(n,1):pind(n,2))
    write(vname,'("rtsoil_",I1.1)') n
    call histwrt3(dat,vname,idnc,iarch,local)   
    dat=cansto
    dat(cmap(pind(1,n):pind(2,n)))=canopy%cansto(pind(n,1):pind(n,2))
    write(vname,'("cansto_",I1.1)') n
    call histwrt3(dat,vname,idnc,iarch,local)
    do k=1,ncp
      dat=cplant(:,k)
      dat(cmap(pind(1,n):pind(2,n)))=bgc%cplant(pind(n,1):pind(n,2),k)
      write(vname,'("cplant",I1.1,"_",I1.1)') k,n
      call histwrt3(dat,vname,idnc,iarch,local)    
    end do
    do k=1,ncs
      dat=csoil(:,k)
      dat(cmap(pind(1,n):pind(2,n)))=bgc%csoil(pind(n,1):pind(n,2),k)
      write(vname,'("csoil",I1.1,"_",I1.1)') k,n
      call histwrt3(dat,vname,idnc,iarch,local)
    end do
  end do
  
  return
  end subroutine savetile 

end module cable_ccam

