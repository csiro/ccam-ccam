! UCLEM urban canopy model
    
! Copyright 2015-2025 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
! This file is part of the UCLEM urban canopy model
!
! UCLEM is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! UCLEM is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with UCLEM.  If not, see <http://www.gnu.org/licenses/>.

!------------------------------------------------------------------------------

! This code was originally inspired by the TEB scheme of Masson, Boundary-Layer Meteorology, 94, p357 (2000)
! The snow scheme is based on Douville, Royer and Mahfouf, Climate Dynamics, 12, p21 (1995)
! The in-canyon vegetation is based on Kowalczyk et al, DAR Tech Paper 32 (1994), but simplified by assuming sigmaf=1.

! The important differences from previous models include an alternative formulation of the conduction model by
! Lipson et al (2017), an alternative formulation for in-canyon aerodynamical resistances based on Harman, et al (2004)
! and Kanada et al (2007), combined with a second canyon wall for completeness.  The scheme includes nrefl order
! reflections in the canyon for both longwave and shortwave radiation (in TEB infinite reflections are used for
! shortwave and 1st order reflections in longwave). A big-leaf vegetation tile is included in the canyon using the
! Kowalczyk et al (1994) scheme but with a simplified soil moisture budget and no modelling of the soil temperature
! since sigmaf=1.  Snow is also included in the canyon and on roofs using a single-layer scheme based on
! Douville, et al (1995).  Time dependent traffic heat fluxes are based on Coutts, et al (2007).
    
! The Urban CLimate Energy Model (UCLEM) results when the exterior canyon model is combined with an interior model
! of energy use from Lipson, et al (2018).  This model includes a representation of energy storage and conduction for the
! interior slab, prognostic internal air temperature, as well as dynamic anthropogenic energy generation as a function
! of internal air temperature and a representation of human behaviours over the diurnal cycle based on Thatcher (2007).


module uclem

implicit none

private
public uclem_eval, init_lwcoeff, init_internal, getswcoeff, getnetalbedo
public uclem_calc_alb, uclem_prepinterior

contains                    
                    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! urban flux calculations

! Basic loop is:
!  --- Canyon fluxes (coupled with interior)
!  Canyon short wave flux (nrefl reflections, depending on snow)
!  Canyon long wave flux (nrefl reflections precomputed)
!  Calculate canyon facet aerodynamic resistances
!  Solve canyon sensible and latent heat budgets
!    Solve for canyon snow and vegetation energy budgets
!      Calculate canyon vegetation temperature
!      Calculate canyon snow temperature
!    End canyon snow and vegetation energy budgets loop
!    Estimate interior fluxes into the room
!    Solve for interior air energy budget and interior air temperature
!    Calculate air conditioning fluxes
!    Calculate canyon air temperature
!    Calculate canyon water vapor mixing ratio
!  End canyon sensible and latent heat budgets loop
!  Canyon longwave, sensible and latent heat fluxes
!  --- Roof fluxes   
!  Solve roof snow and vegetation energy budget
!    Roof snow temperature
!    Roof vegetation temperature
!  End roof snow and vegetation energy budgets loop
!  Roof longwave, sensible and latent heat fluxes
!  --- Interior fluxes
!  Solve interior longwave fluxes
!  --- Conduction    
!  Update urban roof, road and wall and slab temperatures with conduction model
!  Update sensible heat flux (implicit with conduction)
!  --- Hydrology
!  Update water on road surfaces
!  Update water on roof surfaces
!  --- Output (combine fluxes from facets)
!  Combine drag terms for momentum fluxes
!  Combine runoff terms
!  Combine transpiration fluxes
!  Combine roof fluxes
!  Combine hydrology   
!  Combine canyon and roof fluxes
!  Calculate bulk snow properties
!  Estimate bulk roughness length for heat
!  Calculate screen level diagnostics
!  --- Conservation
!  Check energy conservation

! progcalc = 0   Update prognostic and diagnostic variables
! progcalc = 1   Update diagnostic variables only
    
subroutine uclem_eval(u_fg,u_eg,u_ts,u_wf,u_rn,u_evspsbl,u_sbl,ddt,a_sg,a_rg,a_rho,a_temp,a_temproof, &
                      a_mixr,a_mixrroof,a_ps,a_umag,a_umagroof,a_udir,a_rnd,a_snd,a_zmin,a_zroof,fp,  &
                      fp_intm,fp_road,fp_roof,fp_slab,fp_wall,intm,pd,rdhyd,rfhyd,rfveg,road,roof,    &
                      room,slab,walle,wallw,cnveg,intl,ufull,progcalc,diag)

use uclem_parameters
use uclem_types

implicit none

integer, intent(in) :: ufull, progcalc, diag
integer k
real, intent(in) :: ddt
real, dimension(ufull), intent(in) :: a_sg,a_rg,a_rho,a_temp,a_mixr,a_ps,a_umag,a_udir,a_rnd,a_snd,a_zmin
real, dimension(ufull), intent(in) :: a_temproof,a_mixrroof,a_umagroof,a_zroof
real, dimension(ufull), intent(out) :: u_fg,u_eg,u_ts,u_wf,u_rn,u_evspsbl,u_sbl
real, dimension(ufull) :: ggint_roof,ggint_walle,ggint_wallw,ggint_road,ggint_slab,ggint_intm
real, dimension(ufull) :: rdsntemp,rfsntemp,rdsnmelt,rfsnmelt,garfsn,gardsn
real, dimension(ufull) :: wallpsi,roadpsi,fgtop,egtop,qsatr,qsata
real, dimension(ufull) :: cu,fgrooftop,egrooftop
real, dimension(ufull) :: we,ww,wr,zolog,a,n,zom,zonet,dis
real, dimension(ufull) :: roofvegwetfac,roadvegwetfac,d_irrigwater
real, dimension(ufull) :: z_on_l, pa, dts, dtt
real, dimension(ufull) :: u_alb, u_melt, evspsbltop, evspsblrooftop, sbltop, sblrooftop
real, dimension(ufull) :: sg_roof,sg_vegr,sg_road,sg_walle,sg_wallw,sg_vegc,sg_rfsn,sg_rdsn
real, dimension(ufull) :: rg_roof,rg_road,rg_walle,rg_wallw,rg_vegc,rg_vegr,rg_rfsn,rg_rdsn
real, dimension(ufull) :: rgint_roof,rgint_walle,rgint_wallw,rgint_slab,zero_flux
real, dimension(ufull) :: fg_roof,fg_road,fg_walle,fg_wallw,fg_vegc,fg_vegr,fg_rfsn,fg_rdsn
real, dimension(ufull) :: eg_roof,eg_road,eg_vegc,eg_vegr,eg_rfsn,eg_rdsn
real, dimension(ufull) :: acond_roof,acond_road,acond_walle,acond_wallw
real, dimension(ufull) :: acond_vegc,acond_vegr,acond_rfsn,acond_rdsn
real, dimension(ufull) :: abase_road,abase_walle,abase_wallw,abase_vegc,abase_rdsn
real, dimension(ufull) :: d_roofdelta,d_roaddelta,d_vegdeltac,d_vegdeltar,d_rfsndelta,d_rdsndelta
real, dimension(ufull) :: d_tempc,d_mixrc,d_tempr,d_mixrr,d_sigd,d_sigr,d_rfdzmin
real, dimension(ufull) :: d_ac_canyon,d_canyonrgout,d_roofrgout,d_tranc,d_evapc,d_tranr,d_evapr
real, dimension(ufull) :: d_totdepth,d_netemiss,d_netrad,d_topu
real, dimension(ufull) :: d_cwa,d_cw0,d_cww,d_cwr,d_cra,d_crr,d_crw
real, dimension(ufull) :: d_canyontemp,d_canyonmix,d_traf
real, dimension(ufull) :: ggext_roof, ggext_walle, ggext_wallw, ggext_road, ggext_slab, ggext_intm
real, dimension(ufull) :: ggext_impl, ggint_impl, tmp
real, dimension(ufull) :: ac_coeff, ggint
real, dimension(ufull) :: d_ac_inside, d_intgains_bld, int_infilflux
real, dimension(ufull) :: cyc_traffic, cyc_basedemand, cyc_proportion, cyc_translation
real, dimension(ufull) :: int_infilfg, delintercept, rf_delsnow, rd_delsnow
real, dimension(ufull) :: iroomtemp
real, dimension(ufull) :: cvcoeff_roof, cvcoeff_walle, cvcoeff_wallw, cvcoeff_slab
real, dimension(ufull) :: cvcoeff_intm1, cvcoeff_intm2
real, dimension(ufull) :: d_roomstor
real, dimension(ufull) :: d_surfwater, d_soilwater, d_leafwater, d_snow, d_den, d_snowalpha
real, dimension(ufull) :: fp_coeffbldheight
real, dimension(ufull,1:4)  :: iskintemp          ! interim internal skin temperature
real(kind=8), dimension(ufull,0:nl) :: inodetemp  ! interim tridiagonal solution
type(facetparams), intent(in) :: fp_intm, fp_road, fp_roof, fp_slab, fp_wall
type(hydrodata), intent(inout) :: rdhyd, rfhyd
type(vegdata), intent(inout) :: rfveg
type(facetdata), intent(inout) :: road, roof, room, slab, walle, wallw, intm
type(vegdata), intent(inout) :: cnveg
type(intldata), intent(in) :: intl
type(fparmdata), intent(in) :: fp
type(pdiagdata), intent(inout) :: pd

if ( diag>=1 ) write(6,*) "Evaluating UCLEM"

!------------------------------------------------------------------------------
! Canyon fluxes (coupled with interior)
!------------------------------------------------------------------------------

! calculate water and snow area cover fractions
d_rfsndelta = rfhyd%snow/(rfhyd%snow+maxrfsn)
d_rdsndelta = rdhyd%snow/(rdhyd%snow+maxrdsn)

! variables for canyon geometry with vegetation
fp_coeffbldheight = fp%effhwratio/fp%hwratio

! calculate shortwave reflections
! Here we modify the effective canyon geometry to account for in-canyon vegetation
call getswcoeff(sg_roof,sg_vegr,sg_road,sg_walle,sg_wallw,sg_vegc,sg_rfsn,sg_rdsn,wallpsi,roadpsi,fp%hwratio,   &
                fp%effhwratio,fp%vangle,fp%hangle,fp%fbeam,cnveg%sigma,fp_road%alpha,cnveg%alpha,fp_wall%alpha, &
                rdhyd%snowalpha,d_rdsndelta)
sg_walle = sg_walle*fp_coeffbldheight ! shadow due to in-canyon vegetation
sg_wallw = sg_wallw*fp_coeffbldheight ! shadow due to in-canyon vegetation
call getnetalbedo(u_alb,sg_roof,sg_vegr,sg_road,sg_walle,sg_wallw,sg_vegc,sg_rfsn,sg_rdsn,  &
                  fp%hwratio,fp%sigmabld,rfveg%sigma,fp_roof%alpha,rfveg%alpha,             &
                  cnveg%sigma,fp_road%alpha,fp_wall%alpha,cnveg%alpha,                      &
                  rfhyd%snowalpha,rdhyd%snowalpha,d_rfsndelta,d_rdsndelta)
sg_roof  = (1.-fp_roof%alpha)*sg_roof*a_sg
sg_vegr  = (1.-rfveg%alpha)*sg_vegr*a_sg
sg_walle = (1.-fp_wall%alpha)*sg_walle*a_sg
sg_wallw = (1.-fp_wall%alpha)*sg_wallw*a_sg
sg_road  = (1.-fp_road%alpha)*sg_road*a_sg
sg_vegc  = (1.-cnveg%alpha)*sg_vegc*a_sg
sg_rfsn  = (1.-rfhyd%snowalpha)*sg_rfsn*a_sg
sg_rdsn  = (1.-rdhyd%snowalpha)*sg_rdsn*a_sg

! calculate long wave reflections (pregenerated before canyonflux subroutine)
call getlwcoeff(d_netemiss,d_cwa,d_cra,d_cw0,d_cww,d_crw,d_crr,d_cwr,d_rdsndelta, &
                wallpsi,roadpsi,cnveg%sigma,fp_road%emiss,cnveg%emiss,fp_wall%emiss)
                
pd%emiss = d_rfsndelta*snowemiss+(1.-d_rfsndelta)*((1.-rfveg%sigma)*fp_roof%emiss+rfveg%sigma*rfveg%emiss)
pd%emiss = fp%sigmabld*pd%emiss+(1.-fp%sigmabld)*(2.*fp_wall%emiss*fp%effhwratio*d_cwa+d_netemiss*d_cra) ! diagnostic only

! estimate bulk in-canyon surface roughness length
dis   = max(max(max(0.1*fp_coeffbldheight*fp%bldheight,zocanyon+0.2),cnveg%zo+0.2),zosnow+0.2)
zolog = 1./sqrt(d_rdsndelta/log(dis/zosnow)**2+(1.-d_rdsndelta)*(cnveg%sigma/log(dis/cnveg%zo)**2  &
       +(1.-cnveg%sigma)/log(dis/zocanyon)**2))
zonet = dis*exp(-zolog)

! estimate overall urban roughness length
zom = zomratio*fp%bldheight
where ( zom*fp%sigmabld<zonet*(1.-fp%sigmabld) ) ! MJT suggestion
  zom = zonet
end where
n   = rdhyd%snow/(rdhyd%snow+maxrdsn+0.408*grav*zom)   ! snow cover for urban roughness calc (Douville, et al 1995)
zom = (1.-n)*zom + n*zosnow                            ! blend urban and snow roughness lengths (i.e., snow fills canyon)

! Calculate distance from atmosphere to displacement height
d_rfdzmin = max(a_zroof-fp%bldheight,zoroof+0.2,rfveg%zo+0.2) ! distance to roof displacement height
pd%cndzmin = max(a_zmin-refheight*fp%bldheight,1.5,zom+0.2)   ! distance to canyon displacement height
pd%lzom    = log(pd%cndzmin/zom)

! calculate canyon wind speed and bulk transfer coefficients
! (i.e., acond = 1/(aerodynamic resistance) )
! some terms are updated when calculating canyon air temperature
select case(resmeth)
  case(0) ! Masson (2000)
    cu=exp(-0.25*fp%effhwratio)
    abase_road =cu ! bulk transfer coefficients are updated in canyonflux
    abase_walle=cu
    abase_wallw=cu
    abase_rdsn =cu
    abase_vegc =cu
  case(2) ! Kusaka et al (2001)
    cu=exp(-0.386*fp%effhwratio)
    abase_road =cu ! bulk transfer coefficients are updated in canyonflux
    abase_walle=cu
    abase_wallw=cu
    abase_rdsn =cu
    abase_vegc =cu
  case(3) ! Harman et al (2004)
    we=0. ! for cray compiler
    ww=0. ! for cray compiler
    wr=0. ! for cray compiler
    call getincanwindb(we,ww,wr,a_udir,zonet,fp,ufull)
    dis=max(0.1*fp_coeffbldheight*fp%bldheight,zocanyon+0.2)
    zolog=log(dis/zocanyon)
    a=vkar*vkar/(zolog*(2.3+zolog))  ! Assume zot=zom/10.
    abase_walle=a*we                 ! east wall bulk transfer
    abase_wallw=a*ww                 ! west wall bulk transfer
    dis=max(0.1*fp_coeffbldheight*fp%bldheight,zocanyon+0.2,cnveg%zo+0.2,zosnow+0.2)
    zolog=log(dis/zocanyon)
    a=vkar*vkar/(zolog*(2.3+zolog))  ! Assume zot=zom/10.
    abase_road=a*wr                  ! road bulk transfer
    zolog=log(dis/cnveg%zo)
    a=vkar*vkar/(zolog*(2.3+zolog))  ! Assume zot=zom/10.
    abase_vegc=a*wr
    zolog=log(dis/zosnow)
    a=vkar*vkar/(zolog*(2.3+zolog))  ! Assume zot=zom/10.
    abase_rdsn=a*wr                  ! road snow bulk transfer
  case default
    write(6,*) "ERROR: Invalid choice for resmeth ",resmeth
    stop
end select
  
call getdiurnal(fp%ctime,fp%weekdayload,cyc_traffic,cyc_basedemand,cyc_proportion,cyc_translation)
! traffic sensible heat flux
pd%traf = fp%trafficfg*cyc_traffic
d_traf = pd%traf/(1.-fp%sigmabld)
! internal gains sensible heat flux
d_intgains_bld = (fp%intmassn+1.)*fp%intgains_flr*cyc_basedemand ! building internal gains 
pd%intgains_full= fp%sigmabld*d_intgains_bld                     ! full domain internal gains

! calculate water and snow area cover fractions
d_roaddelta = max(rdhyd%surfwater/maxrdwater,0.)**(2./3.)
d_vegdeltac = max(rdhyd%leafwater/max(maxvwatf*cnveg%lai,1.E-8),0.)**(2./3.)

! canyon level air temp and water vapor (displacement height at refheight*building height)
pa      = a_ps*exp(-grav*a_zmin/(rd*(a_temp+urbtemp)))
d_sigd  = a_ps
a       = (d_sigd/pa)**(rd/aircp)
!d_tempc = (a_temp+urbtemp)*a - urbtemp
d_tempc = a_temp*a + urbtemp*(a-1.)
call getqsat(qsatr,d_tempc,d_sigd)
call getqsat(qsata,a_temp,pa)
d_mixrc = a_mixr*qsatr/qsata

! calculate soil data
d_totdepth = sum( fp_road%depth, dim=2 )

! calculate canyon fluxes
call solvecanyon(sg_road,rg_road,fg_road,eg_road,acond_road,abase_road,                          &
                 sg_walle,rg_walle,fg_walle,acond_walle,abase_walle,                             &
                 sg_wallw,rg_wallw,fg_wallw,acond_wallw,abase_wallw,                             & 
                 sg_vegc,rg_vegc,fg_vegc,eg_vegc,acond_vegc,abase_vegc,                          &
                 sg_rdsn,rg_rdsn,fg_rdsn,eg_rdsn,acond_rdsn,abase_rdsn,rdsntemp,rdsnmelt,gardsn, &
                 a_umag,a_rho,a_rg,a_rnd,a_snd,we,ww,                                            &
                 d_canyontemp,d_canyonmix,d_tempc,d_mixrc,d_sigd,d_topu,d_netrad,                &
                 d_roaddelta,d_vegdeltac,d_rdsndelta,d_ac_canyon,d_traf,d_ac_inside,             &
                 d_canyonrgout,d_tranc,d_evapc,d_cwa,d_cra,d_cw0,d_cww,d_crw,d_crr,              &
                 d_cwr,d_totdepth,d_intgains_bld,d_roomstor,fgtop,egtop,evspsbltop,sbltop,       &
                 int_infilflux,int_infilfg,cyc_proportion,cyc_translation,ddt,                   &
                 cnveg,fp,fp_intm,fp_road,fp_roof,fp_wall,fp_slab,intm,pd,rdhyd,rfveg,road,      &
                 roof,room,slab,walle,wallw,iroomtemp,cvcoeff_roof,cvcoeff_walle,                &
                 cvcoeff_wallw,cvcoeff_slab,cvcoeff_intm1,cvcoeff_intm2,ufull,diag)


!------------------------------------------------------------------------------
! Roof fluxes
!------------------------------------------------------------------------------

! calculate water and snow area cover fractions
d_roofdelta = max(rfhyd%surfwater/maxrfwater,0.)**(2./3.)
d_vegdeltar = max(rfhyd%leafwater/max(maxvwatf*rfveg%lai,1.E-8),0.)**(2./3.)

! roof level air temperature and water vapor (displacement height at building height)
pa      = a_ps*exp(-grav*a_zroof/(rd*(a_temproof+urbtemp)))
d_sigr  = a_ps*exp(-grav*fp%bldheight*(1.-refheight)/(rd*(a_temproof+urbtemp)))
a       = (d_sigr/pa)**(rd/aircp)
!d_tempr = (a_temproof+urbtemp)*a - urbtemp
d_tempr = a_temproof*a + urbtemp*(a-1.)
call getqsat(qsatr,d_tempr,d_sigr)
call getqsat(qsata,a_temproof,pa)
d_mixrr = a_mixrroof*qsatr/qsata

! calculate roof fluxes (fg_roof updated below)
eg_roof = 0. ! For cray compiler
call solveroof(sg_rfsn,rg_rfsn,fg_rfsn,eg_rfsn,garfsn,rfsnmelt,rfsntemp,acond_rfsn,d_rfsndelta, &
               sg_vegr,rg_vegr,fg_vegr,eg_vegr,acond_vegr,d_vegdeltar,                          &
               sg_roof,rg_roof,eg_roof,acond_roof,d_roofdelta,                                  &
               a_rg,a_umagroof,a_rho,a_rnd,a_snd,d_tempr,d_mixrr,d_rfdzmin,d_tranr,d_evapr,     &
               d_sigr,ddt,fp_roof,rfhyd,rfveg,roof,fp,ufull)


!------------------------------------------------------------------------------
! Interior fluxes
!------------------------------------------------------------------------------

! update interior fluxes (first estimate)
if ( intairtmeth==0 .or. lwintmeth==0 ) then
  rgint_slab  = 0.
  rgint_wallw = 0.
  rgint_roof  = 0.
  rgint_walle = 0.    
else if ( lwintmeth==1 .or. progcalc>0 ) then
  iskintemp(:,1) = real( slab%nodetemp(:,nl) )
  iskintemp(:,2) = real( wallw%nodetemp(:,nl) )
  iskintemp(:,3) = real( roof%nodetemp(:,nl) )
  iskintemp(:,4) = real( walle%nodetemp(:,nl) )
  ! calculate interior lwfluxes using last timestep temps lw (rg) fluxes
  call internal_lwflux(rgint_slab,rgint_wallw,rgint_roof,rgint_walle,    &
                       iskintemp,fp,intl,ufull)
end if  


!------------------------------------------------------------------------------
! Conduction
!------------------------------------------------------------------------------

if ( progcalc==0 ) then

  zero_flux = 0. ! just a constant
    
  ggext_roof = (1.-d_rfsndelta)*(sg_roof+rg_roof-eg_roof+aircp*a_rho*d_tempr*acond_roof) &
                +d_rfsndelta*garfsn
  ggext_walle= sg_walle+rg_walle+aircp*a_rho*d_canyontemp*acond_walle*fp_coeffbldheight
  ggext_wallw= sg_wallw+rg_wallw+aircp*a_rho*d_canyontemp*acond_wallw*fp_coeffbldheight
  ggext_road = (1.-d_rdsndelta)*(sg_road+rg_road-eg_road                        &
               +aircp*a_rho*d_canyontemp*acond_road) + d_rdsndelta*gardsn
  ggext_slab = 0.
  if ( intairtmeth/=0 ) then
    if ( intmassmeth/=0 ) then  
      ggext_intm = cvcoeff_intm1*iroomtemp
    else
      ggext_intm = 0.  
    end if
  else
    ggext_intm = 0.
  end if  

  if ( intairtmeth/=0 .and. lwintmeth==2 ) then
    ! explicit longwave using update skin temperature from intermediate facet temp calc
    rgint_slab  = 0.
    rgint_wallw = 0.
    rgint_roof  = 0.
    rgint_walle = 0.

    !!!! slab !!!!
    ggint_impl = cvcoeff_slab
    ggint = -cvcoeff_slab*iroomtemp
    ggext_impl = 0.  ! zero_flux boundary condition
    call solvetridiag(ggext_slab,ggint,rgint_slab,ggext_impl,ggint_impl,slab%nodetemp,  &
                      inodetemp,fp_slab%depth,fp_slab%volcp,fp_slab%lambda,ddt,ufull)
    iskintemp(:,1) = real(inodetemp(:,nl))

    !!!! wallw !!!!
    ggint_impl = cvcoeff_wallw
    ggint = -cvcoeff_wallw*iroomtemp
    ggext_impl = aircp*a_rho*acond_wallw*fp_coeffbldheight ! implicit, later update fg_wallw with final T
    call solvetridiag(ggext_wallw,ggint,rgint_wallw,ggext_impl,ggint_impl,wallw%nodetemp,  &
                      inodetemp,fp_wall%depth,fp_wall%volcp,fp_wall%lambda,ddt,ufull)
    iskintemp(:,2) = real(inodetemp(:,nl))

    !!!! roof !!!!
    ggint_impl = cvcoeff_roof
    ggint = -cvcoeff_roof*iroomtemp
    ggext_impl = (1.-d_rfsndelta)*aircp*a_rho*acond_roof  ! implicit, later update fg_roof with final T
    call solvetridiag(ggext_roof,ggint,rgint_roof,ggext_impl,ggint_impl,roof%nodetemp,  &
                      inodetemp,fp_roof%depth,fp_roof%volcp,fp_roof%lambda,ddt,ufull)
    iskintemp(:,3) = real(inodetemp(:,nl))

    !!!! walle !!!!
    ggint_impl = cvcoeff_walle
    ggint = -cvcoeff_walle*iroomtemp
    ggext_impl = aircp*a_rho*acond_walle*fp_coeffbldheight ! implicit, later update fg_walle with final T
    call solvetridiag(ggext_walle,ggint,rgint_walle,ggext_impl,ggint_impl,walle%nodetemp,  &
                      inodetemp,fp_wall%depth,fp_wall%volcp,fp_wall%lambda,ddt,ufull)         
    iskintemp(:,4) = real(inodetemp(:,nl))

    ! calculate interior lwfluxes using iskintemp
    call internal_lwflux(rgint_slab,rgint_wallw,rgint_roof,rgint_walle,    &
                         iskintemp,fp,intl,ufull)
  end if ! intairtmeth/=0 .and. lwintmeth==2


  ! solve final conduction through facets with explicit interior, implicit exterior boundary
  ggint_impl = cvcoeff_roof
  ggint = -cvcoeff_roof*iroomtemp
  ggext_impl = (1.-d_rfsndelta)*aircp*a_rho*acond_roof  ! later update fg_roof with final roof skin T
  call solvetridiag(ggext_roof,ggint,rgint_roof,ggext_impl,ggint_impl,roof%nodetemp,     &
                    roof%nodetemp,fp_roof%depth,fp_roof%volcp,fp_roof%lambda,ddt,ufull)

  ggint_impl = cvcoeff_walle
  ggint = -cvcoeff_walle*iroomtemp
  ggext_impl = aircp*a_rho*acond_walle*fp_coeffbldheight ! later update fg_walle with final walle skin T
  call solvetridiag(ggext_walle,ggint,rgint_walle,ggext_impl,ggint_impl,walle%nodetemp,  &
                    walle%nodetemp,fp_wall%depth,fp_wall%volcp,fp_wall%lambda,ddt,ufull)

  ggint_impl = cvcoeff_wallw
  ggint = -cvcoeff_wallw*iroomtemp
  ggext_impl = aircp*a_rho*acond_wallw*fp_coeffbldheight ! later update fg_wallw with final wallw skin T
  call solvetridiag(ggext_wallw,ggint,rgint_wallw,ggext_impl,ggint_impl,wallw%nodetemp,  &
                    wallw%nodetemp,fp_wall%depth,fp_wall%volcp,fp_wall%lambda,ddt,ufull)

  ggint_impl = 0.  ! zero_flux boundary condition
  ggint_road = 0.  ! zero_flux boundary condition
  ggext_impl = (1.-d_rdsndelta)*aircp*a_rho*acond_road ! later update fg_road with final road skin T
  call solvetridiag(ggext_road,ggint_road,zero_flux,ggext_impl,ggint_impl,road%nodetemp,   &
                    road%nodetemp,fp_road%depth,fp_road%volcp,fp_road%lambda,ddt,ufull)

  ! calculate internal facet conduction
  if ( intairtmeth==0 ) then
    slab%nodetemp = 1._8 ! + urbtemp 
  else    
    ggint_impl = cvcoeff_slab
    ggint = -cvcoeff_slab*iroomtemp
    ggext_slab = 0.  ! zero_flux boundary condition
    ggext_impl = 0.  ! zero_flux boundary condition
    call solvetridiag(ggext_slab,ggint,rgint_slab,ggext_impl,ggint_impl,slab%nodetemp,    &
                      slab%nodetemp,fp_slab%depth,fp_slab%volcp,fp_slab%lambda,ddt,ufull)
    if ( intmassmeth/=0 ) then
      ggint_impl = cvcoeff_intm2
      ggint = -cvcoeff_intm2*iroomtemp
      ggext_impl = cvcoeff_intm1
      call solvetridiag(ggext_intm,ggint,zero_flux,ggext_impl,ggint_impl,intm%nodetemp,    &
                        intm%nodetemp,fp_intm%depth,fp_intm%volcp,fp_intm%lambda,ddt,ufull)
    end if
  end if
  
end if ! if progcalc==0 ... 


! Update sensible heat fluxes (implicit with conduction)

fg_roof  = aircp*a_rho*(real(roof%nodetemp(:,0))-d_tempr)*acond_roof
fg_walle = aircp*a_rho*(real(walle%nodetemp(:,0))-d_canyontemp)*acond_walle*fp_coeffbldheight
fg_wallw = aircp*a_rho*(real(wallw%nodetemp(:,0))-d_canyontemp)*acond_wallw*fp_coeffbldheight
fg_road  = aircp*a_rho*(real(road%nodetemp(:,0))-d_canyontemp)*acond_road
if ( intairtmeth/=0 ) then
  if ( intmassmeth/=0 ) then  
    ggext_intm = -cvcoeff_intm1*(real(intm%nodetemp(:,0))-iroomtemp)
  else
    ggext_intm = 0.  
  end if
else
  ggext_intm = 0.  
end if  

ggint_roof  = cvcoeff_roof*(real(roof%nodetemp(:,nl))-iroomtemp)
ggint_walle = cvcoeff_walle*(real(walle%nodetemp(:,nl))-iroomtemp)
ggint_wallw = cvcoeff_wallw*(real(wallw%nodetemp(:,nl))-iroomtemp)
if ( intairtmeth/=0 ) then
  ggint_slab = cvcoeff_slab*(real(slab%nodetemp(:,nl))-iroomtemp)  
  if ( intmassmeth/=0 ) then
    ggint_intm = cvcoeff_intm2*(real(intm%nodetemp(:,nl))-iroomtemp) 
  else
    ggint_intm = 0.  
  end if    
else
  ggint_slab = 0.  
  ggint_intm = 0.
end if

d_ac_inside = -ggint_roof - ggint_slab + (ggext_intm - ggint_intm)*real(fp%intmassn) & 
              - (ggint_walle+ggint_wallw)*(fp%bldheight/fp%bldwidth)                 &
              - int_infilflux - d_intgains_bld + d_roomstor
where ( d_ac_inside<0. )
  ! modified UCLEM COP calculation (based on "Coefficient of performance for finite speed heat pump" - Blanchard 1980)
  ac_coeff = max( acfactor*(d_canyontemp+30.-iroomtemp)/(iroomtemp+30.+urbtemp), 0. )
elsewhere 
  ac_coeff = 0.
end where
select case(acmeth) ! AC heat pump into canyon (0=Off, 1=On, 2=Reversible, COP of 1.0)
  case(0) ! unrealistic cooling (buildings act as heat sink)
    d_ac_canyon = 0.
    pd%bldheat = max(0.,d_ac_inside*fp%sigmabld)
    pd%bldcool = max(0.,d_ac_inside*fp%sigmabld)
  case(1) ! d_ac_canyon pumps conducted heat + ac waste heat back into canyon
    d_ac_canyon = max(0.,-d_ac_inside*(1.+ac_coeff)*fp%sigmabld/(1.-fp%sigmabld))  ! canyon domain W/m/m
    pd%bldheat = max(0.,d_ac_inside*fp%sigmabld)            ! entire domain W/m/m
    pd%bldcool = max(0.,-d_ac_inside*ac_coeff*fp%sigmabld)  ! entire domain W/m/m
  case(2) ! reversible heating and cooling (for testing energy conservation)
    d_ac_canyon = -d_ac_inside*fp%sigmabld/(1.-fp%sigmabld)
    pd%bldheat = 0.
    pd%bldcool = 0.
  case DEFAULT
    write(6,*) "ERROR: Unknown acmeth mode ",acmeth
    stop
end select

! update canyon flux
fgtop = fp%hwratio*(fg_walle+fg_wallw) + (1.-d_rdsndelta)*(1.-cnveg%sigma)*fg_road &
      + (1.-d_rdsndelta)*cnveg%sigma*fg_vegc + d_rdsndelta*fg_rdsn                 &
      + d_traf + d_ac_canyon - int_infilfg


!------------------------------------------------------------------------------
! Hydrology
!------------------------------------------------------------------------------

! Road hydrology

! reset inception of leaf water
pd%delintercept = 0.

d_surfwater = rdhyd%surfwater
d_soilwater = rdhyd%soilwater
d_leafwater = rdhyd%leafwater
d_snow = rdhyd%snow
d_den = rdhyd%den
d_snowalpha = rdhyd%snowalpha

! new snowfall
where ( a_snd>1.e-10 )
  ! update snow density
  d_den = (d_snow*d_den+a_snd*ddt*minsnowden)/(d_snow+ddt*a_snd)
  ! reset snow albedo
  d_snowalpha = maxsnowalpha
end where

! calculate water/snow budgets for road surface
call updatewater(ddt,d_surfwater,d_soilwater,d_leafwater,d_snow,              &
                 d_den,d_snowalpha,rdsnmelt,a_rnd,a_snd,eg_road,              &
                 eg_rdsn,d_tranc,d_evapc,d_totdepth, cnveg%lai,wbrelaxc,      &
                 fp%sfc,fp%swilt,d_irrigwater,delintercept,rd_delsnow,ufull)
! record grid irrigation road component [kg m-2 s-1]
pd%irrig = d_irrigwater*waterden*d_totdepth*cnveg%sigma*(1.-fp%sigmabld)/ddt
! record leaf intercept component [kg m-2 s-1]
pd%delintercept = pd%delintercept + delintercept*cnveg%sigma

pd%road_water_runoff = max(d_surfwater-maxrdwater,0.)
pd%road_snow_runoff  = max(d_snow-maxrdsn,0.)
pd%road_soil_runoff  = max(d_soilwater-fp%ssat,0.)

! Update road prognostic variables
if ( progcalc==0 ) then
  rdhyd%surfwater(1:ufull) = max(d_surfwater - pd%road_water_runoff,0.)
  rdhyd%soilwater(1:ufull) = max(d_soilwater - pd%road_soil_runoff,fp%swilt)
  rdhyd%leafwater          = max(d_leafwater,0.)
  rdhyd%snow(1:ufull)      = max(d_snow - pd%road_snow_runoff,0.)
  rdhyd%den                = min(max(d_den,minsnowden),maxsnowden)
  rdhyd%snowalpha          = min(max(d_snowalpha,minsnowalpha),maxsnowalpha)
end if  

! Roof hydrology

d_surfwater = rfhyd%surfwater
d_soilwater = rfhyd%soilwater
d_leafwater = rfhyd%leafwater
d_snow = rfhyd%snow
d_den = rfhyd%den
d_snowalpha = rfhyd%snowalpha

! new snowfall
where ( a_snd>1.e-10 )
  ! update snow density
  d_den = (rfhyd%snow*rfhyd%den+a_snd*ddt*minsnowden)/(rfhyd%snow+ddt*a_snd)
  ! reset snow albedo
  d_snowalpha = maxsnowalpha
end where

! calculate water/snow budgets for roof surface
call updatewater(ddt,d_surfwater,d_soilwater,d_leafwater,d_snow,             &
                 d_den,d_snowalpha,rfsnmelt,a_rnd,a_snd,eg_roof,             &
                 eg_rfsn,d_tranr,d_evapr,fp%rfvegdepth,rfveg%lai,wbrelaxr,   &
                 fp%sfc,fp%swilt,d_irrigwater,delintercept,rf_delsnow,ufull)
! add grid irrigation roof component [kg m-2 s-1]
pd%irrig = pd%irrig + d_irrigwater*waterden*fp%rfvegdepth*rfveg%sigma*fp%sigmabld/ddt
! record leaf intercept component [kg m-2 s-1]
pd%delintercept = pd%delintercept + delintercept*rfveg%sigma

pd%roof_water_runoff = max(d_surfwater-maxrfwater,0.)
pd%roof_snow_runoff  = max(d_snow-maxrfsn,0.)
pd%roof_soil_runoff  = max(d_soilwater-fp%ssat,0.)

! Update roof prognostic variables
if ( progcalc==0 ) then
  rfhyd%surfwater(1:ufull) = max(d_surfwater - pd%roof_water_runoff,0.)
  rfhyd%soilwater(1:ufull) = max(d_soilwater - pd%roof_soil_runoff,fp%swilt)
  rfhyd%leafwater          = max(d_leafwater,0.)
  rfhyd%snow(1:ufull)      = max(d_snow - pd%roof_snow_runoff,0.)
  rfhyd%den                = min(max(d_den,minsnowden),maxsnowden)
  rfhyd%snowalpha          = min(max(d_snowalpha,minsnowalpha),maxsnowalpha)
end if


!------------------------------------------------------------------------------
! Output (combine fluxes from facets)
!------------------------------------------------------------------------------

! area weighted vegetation aerodynamic conductance
pd%acond_vegw = (acond_vegc*cnveg%sigma*(1.-fp%sigmabld) + acond_vegr*rfveg%sigma*fp%sigmabld)/  & 
               (rfveg%sigma*fp%sigmabld + cnveg%sigma*(1.-fp%sigmabld))

! calculate runoff (leafwater runoff already accounted for in precip reaching canyon floor) [kg m-2]
u_rn = pd%roof_water_runoff*fp%sigmabld*(1.-rfveg%sigma)                     & ! roof surfwater runnoff 
      +pd%road_water_runoff*(1.-fp%sigmabld)*(1.-cnveg%sigma)                & ! road surfwater runnoff
      +pd%roof_snow_runoff*fp%sigmabld                                       & ! roof snow removal     
      +pd%road_snow_runoff*(1.-fp%sigmabld)                                  & ! road snow removal     
      +pd%roof_soil_runoff*waterden*fp%rfvegdepth*rfveg%sigma*fp%sigmabld    & ! roofveg saturation runnoff
      +pd%road_soil_runoff*waterden*d_totdepth*cnveg%sigma*(1.-fp%sigmabld)    ! canyveg saturation runnoff

! diagnose surface drainage rate [kg m-2 s-1]
pd%surfrunoff = u_rn/ddt

! diagnose total vegetation transpiration [kg m-2 s-1]
pd%transveg = (d_tranc*cnveg%sigma*(1.-fp%sigmabld) + d_tranr*rfveg%sigma*fp%sigmabld)/lv

! combine roof fluxes from snow and snow-free
d_roofrgout = a_rg-d_rfsndelta*rg_rfsn-(1.-d_rfsndelta)*((1.-rfveg%sigma)*rg_roof+rfveg%sigma*rg_vegr)
fgrooftop   = d_rfsndelta*fg_rfsn+(1.-d_rfsndelta)*((1.-rfveg%sigma)*fg_roof+rfveg%sigma*fg_vegr)
egrooftop   = d_rfsndelta*eg_rfsn+(1.-d_rfsndelta)*((1.-rfveg%sigma)*eg_roof+rfveg%sigma*eg_vegr)
evspsblrooftop = d_rfsndelta*eg_rfsn/ls+(1.-d_rfsndelta)*((1.-rfveg%sigma)*eg_roof+rfveg%sigma*eg_vegr)/lv
sblrooftop  = d_rfsndelta*eg_rfsn/ls
!fgtop       = d_rdsndelta*fg_rdsn+(1.-d_rdsndelta)*((1.-cnveg%sigma)*fg_road+cnveg%sigma*fg_vegc)   &
!             +fp%hwratio*(fg_walle+fg_wallw)+d_traf+d_ac_canyon
!egtop       = d_rdsndelta*eg_rdsn+(1.-d_rdsndelta)*((1.-cnveg%sigma)*eg_road+cnveg%sigma*eg_vegc)

! calculate wetfac for roof and road vegetation (see sflux.f or cable_canopy.f90)
roofvegwetfac = max(min((rfhyd%soilwater-fp%swilt)/(fp%sfc-fp%swilt),1.),0.)
roadvegwetfac = max(min((rdhyd%soilwater-fp%swilt)/(fp%sfc-fp%swilt),1.),0.)

! area weighted soil water (canyon and roof)
pd%soilwater = (rdhyd%soilwater*cnveg%sigma*(1.-fp%sigmabld) + rfhyd%soilwater*rfveg%sigma*fp%sigmabld) / &
               max(cnveg%sigma*(1.-fp%sigmabld) + rfveg%sigma*fp%sigmabld, 1.e-10)

! record previous soilmoist
tmp = pd%soilmoist
! update soilmoist only in areas which support vegetation (i.e. have soil moisture)
pd%soilmoist = ( rdhyd%soilwater*waterden*d_totdepth*cnveg%sigma*(1.-fp%sigmabld)   &  ! canyon integrated soil moisture
               + rfhyd%soilwater*waterden*fp%rfvegdepth*rfveg%sigma*fp%sigmabld ) / &  ! roof integrated soil moisture
               max(cnveg%sigma*(1.-fp%sigmabld) + rfveg%sigma*fp%sigmabld, 1.e-10)     ! limited to areas with soil
! calculate difference in soilmoist, weighted over entire grid (i.e. include dry area under roads + buildings)
if ( pd%first_call ) then
  pd%delsoilmoist = 0.
else
  pd%delsoilmoist = (pd%soilmoist - tmp)*(rfveg%sigma*fp%sigmabld + cnveg%sigma*(1.-fp%sigmabld))
end if

! rootmoist calculated only for canyon veg subgrid
pd%rootmoistc = rdhyd%soilwater*waterden*d_totdepth

! soilwetness: (soil moisture - wilting) / (saturation - wilting).
tmp = (rdhyd%soilwater*cnveg%sigma*(1.-fp%sigmabld) + rfhyd%soilwater*rfveg%sigma*fp%sigmabld) / &
      max(cnveg%sigma*(1.-fp%sigmabld) + rfveg%sigma*fp%sigmabld, 1.e-10)
pd%soilwetness = (tmp-fp%swilt)/(fp%ssat-fp%swilt)

! calculate longwave, sensible heat latent heat outputs
! estimate surface temp from outgoing longwave radiation
u_ts = ((fp%sigmabld*d_roofrgout+(1.-fp%sigmabld)*d_canyonrgout)/sbconst)**0.25 - urbtemp
u_fg = fp%sigmabld*fgrooftop+(1.-fp%sigmabld)*fgtop+fp%industryfg
u_eg = fp%sigmabld*egrooftop+(1.-fp%sigmabld)*egtop
u_evspsbl = fp%sigmabld*evspsblrooftop + (1.-fp%sigmabld)*evspsbltop
u_sbl = fp%sigmabld*sblrooftop + (1.-fp%sigmabld)*sbltop
u_wf = fp%sigmabld*(1.-d_rfsndelta)*((1.-rfveg%sigma)*d_roofdelta       &
      +rfveg%sigma*((1.-d_vegdeltar)*roofvegwetfac+d_vegdeltar))        &
      +(1.-fp%sigmabld)*(1.-d_rdsndelta)*((1.-cnveg%sigma)*d_roaddelta  &
      +cnveg%sigma*((1.-d_vegdeltac)*roadvegwetfac+d_vegdeltac))

pd%snowmelt = fp%sigmabld*rfsnmelt + (1.-fp%sigmabld)*rdsnmelt
u_melt = lf*(fp%sigmabld*d_rfsndelta*rfsnmelt + (1.-fp%sigmabld)*d_rdsndelta*rdsnmelt)

! calculate average snow temperature and snow water equivilent
pd%snowt = ( (rfsntemp+urbtemp)*d_rfsndelta*(1.-fp%sigmabld) + (rdsntemp+urbtemp)*d_rdsndelta*fp%sigmabld ) / &
           max( d_rfsndelta*(1.-fp%sigmabld) + d_rdsndelta*fp%sigmabld, 1.e-10 )
pd%swe = rfhyd%snow*d_rfsndelta*(1.-fp%sigmabld) + rdhyd%snow*d_rdsndelta*fp%sigmabld
pd%snowfrac = d_rfsndelta*(1.-fp%sigmabld) + d_rdsndelta*fp%sigmabld
pd%salbedo = ( rfhyd%snowalpha*d_rfsndelta*(1.-fp%sigmabld) + rdhyd%snowalpha*d_rdsndelta*fp%sigmabld ) / &
           max( d_rfsndelta*(1.-fp%sigmabld) + d_rdsndelta*fp%sigmabld, 1.e-10 ) 
pd%delswe = rf_delsnow*d_rfsndelta*(1.-fp%sigmabld) + rd_delsnow*d_rdsndelta*fp%sigmabld

! calculate total liquid water storage
pd%surfstor = rdhyd%surfwater*(1.-fp%sigmabld)*(1.-d_rfsndelta)*(1.-rfveg%sigma) &
            + rdhyd%surfwater*fp%sigmabld*(1.-d_rdsndelta)*(1.-cnveg%sigma)

! calculate average vegetation temperature
pd%vegt = ( (cnveg%temp+urbtemp)*(1.-fp%sigmabld)*cnveg%sigma + (rfveg%temp+urbtemp)*fp%sigmabld*rfveg%sigma ) / &
          max( (1.-fp%sigmabld)*cnveg%sigma + fp%sigmabld*rfveg%sigma, 1.e-10 )
pd%calbedo = ( cnveg%alpha*(1.-fp%sigmabld)*cnveg%sigma + rfveg%alpha*fp%sigmabld*rfveg%sigma ) / &
          max( (1.-fp%sigmabld)*cnveg%sigma + fp%sigmabld*rfveg%sigma, 1.e-10 )
pd%ulai = ( cnveg%lai*(1.-fp%sigmabld)*cnveg%sigma + rfveg%lai*fp%sigmabld*rfveg%sigma ) / &
          max( (1.-fp%sigmabld)*cnveg%sigma + fp%sigmabld*rfveg%sigma, 1.e-10 )

! Save canyon air temperature
pd%taircanyon = d_canyontemp


! MOST calculations for bulk properties
! (re)calculate heat roughness length for MOST (diagnostic only)
call getqsat(a,u_ts,d_sigd)
a = a*u_wf + min(a_mixr,a)*(1.-u_wf)
dts = u_ts + (u_ts+urbtemp)*0.61*a
dtt = d_tempc + (d_tempc+urbtemp)*0.61*d_mixrc
select case(zohmeth)
  case(0) ! Use veg formulation
    pd%lzoh = 2.3+pd%lzom
    call getinvres(pd%cdtq,pd%cduv,z_on_l,pd%lzoh,pd%lzom,pd%cndzmin,dts,dtt,a_umag,1)
  case(1) ! Use Kanda parameterisation
    pd%lzoh = 2.3+pd%lzom ! replaced in getlna
    call getinvres(pd%cdtq,pd%cduv,z_on_l,pd%lzoh,pd%lzom,pd%cndzmin,dts,dtt,a_umag,2)
  case(2) ! Use Kanda parameterisation
    pd%lzoh = 6.+pd%lzom
    call getinvres(pd%cdtq,pd%cduv,z_on_l,pd%lzoh,pd%lzom,pd%cndzmin,dts,dtt,a_umag,4)
end select

! calculate screen level diagnostics
call scrncalc(a_mixr,a_umag,a_temp,u_ts,d_tempc,d_rdsndelta,d_roaddelta,d_vegdeltac,d_sigd,a,rdsntemp,zonet, &
              cnveg,fp,pd,rdhyd,road,ufull)


!------------------------------------------------------------------------------
! Conservation
!------------------------------------------------------------------------------

if ( progcalc==0 ) then

  call energyclosure(sg_roof,rg_roof,fg_roof,sg_walle,rg_walle,fg_walle,     &
                     sg_road,rg_road,fg_road,sg_wallw,rg_wallw,fg_wallw,     &
                     rgint_roof,rgint_walle,rgint_wallw,rgint_slab,          &
                     eg_roof,eg_road,garfsn,gardsn,d_rfsndelta,d_rdsndelta,  &
                     a_sg,a_rg,u_ts,u_fg,u_eg,u_alb,u_melt,                  &
                     ggint_roof,ggint_road,ggint_walle,ggint_wallw,          &
                     ggext_intm,ggint_slab,ggint_intm,d_intgains_bld,        &
                     int_infilflux,d_ac_inside,fp,ddt,                       &
                     cnveg,fp_intm,fp_road,fp_roof,fp_slab,fp_wall,intm,pd,  &
                     rfveg,road,roof,room,slab,walle,wallw,ufull)
end if

return
end subroutine uclem_eval

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Tridiagonal solver for temperatures

! This version has an implicit estimate for roof sensible heat flux

! [ ggB ggC         ] [ temp ] = [ ggD ]
! [ ggA ggB ggC     ] [ temp ] = [ ggD ]
! [     ggA ggB ggC ] [ temp ] = [ ggD ]
! [         ggA ggB ] [ temp ] = [ ggD ]

subroutine solvetridiag(ggext,ggint,rgint,ggext_impl,ggint_impl,    &
                        nodetemp,inodetemp,depth,volcp,lambda,ddt,ufull)

use uclem_parameters

implicit none

integer, intent(in) :: ufull
real, dimension(ufull),     intent(in)    :: ggext,ggint,rgint  ! surface energy fluxes
real, dimension(ufull),     intent(in)    :: ggext_impl         ! implicit update for stability
real, dimension(ufull),     intent(in)    :: ggint_impl         ! implicit update for stability
real(kind=8), dimension(ufull,0:nl),intent(in)  :: nodetemp     ! temperature of each node
real(kind=8), dimension(ufull,0:nl),intent(out) :: inodetemp    ! updated temperature solution
real, dimension(ufull,nl),  intent(in)    :: depth,volcp,lambda ! facet depth, heat capacity, conductivity
real(kind=8), dimension(ufull,nl)         :: cap,res            ! layer capacitance & resistance
real(kind=8), dimension(ufull,0:nl)       :: ggA,ggB,ggC,ggD    ! tridiagonal matrices
real(kind=8) ggX                                                ! tridiagonal coefficient
real(kind=8), dimension(ufull)            :: ans                ! tridiagonal solution
real, intent(in)                          :: ddt                ! timestep
integer k, iqu

res = real(depth,8)/real(lambda,8)
cap = real(depth,8)*real(volcp,8)
!!!!!!!!! interface conduction !!!!!!!!!!!
ggA(:,1:nl)   = -1._8/res(:,1:nl)
ggB(:,0)      =  1._8/res(:,1) +0.5_8*cap(:,1)/real(ddt,8) + real(ggext_impl,8)
ggB(:,1:nl-1) =  1._8/res(:,1:nl-1) +1._8/res(:,2:nl) +0.5_8*(cap(:,1:nl-1) &
                   +cap(:,2:nl))/real(ddt,8)
ggB(:,nl)     =  1._8/res(:,nl) + 0.5_8*cap(:,nl)/real(ddt,8) + real(ggint_impl,8)
ggC(:,0:nl-1) = -1._8/res(:,1:nl)
ggD(:,0)      = nodetemp(:,0)*0.5_8*cap(:,1)/real(ddt,8) + real(ggext,8)
ggD(:,1:nl-1) = nodetemp(:,1:nl-1)*0.5_8*(cap(:,1:nl-1)+cap(:,2:nl))/real(ddt,8)
ggD(:,nl)     = nodetemp(:,nl)*0.5_8*cap(:,nl)/real(ddt,8) - real(ggint,8) - real(rgint,8)
! tridiagonal solver (Thomas algorithm) to solve node temperatures
do k = 1,nl
  do iqu = 1,ufull
    ggX = ggA(iqu,k)/ggB(iqu,k-1)
    ggB(iqu,k) = ggB(iqu,k)-ggX*ggC(iqu,k-1)
    ggD(iqu,k) = ggD(iqu,k)-ggX*ggD(iqu,k-1)
  end do  
end do
ans = ggD(:,nl)/ggB(:,nl)
inodetemp(:,nl) = ans
do k = nl-1,0,-1
  ans = (ggD(:,k) - ggC(:,k)*ans)/ggB(:,k)
  inodetemp(:,k) = ans
end do

end subroutine solvetridiag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Conservation of energy check

subroutine energyclosure(sg_roof,rg_roof,fg_roof,sg_walle,rg_walle,fg_walle,     &
                         sg_road,rg_road,fg_road,sg_wallw,rg_wallw,fg_wallw,     &
                         rgint_roof,rgint_walle,rgint_wallw,rgint_slab,          &
                         eg_roof,eg_road,garfsn,gardsn,d_rfsndelta,d_rdsndelta,  &
                         a_sg,a_rg,u_ts,u_fg,u_eg,u_alb,u_melt,                  &
                         ggint_roof,ggint_road,ggint_walle,ggint_wallw,          &
                         ggext_intm,ggint_slab,ggint_intm,d_intgains_bld,      &
                         int_infilflux,d_ac_inside,fp,ddt,cnveg,fp_intm,fp_road, &
                         fp_roof,fp_slab,fp_wall,intm,pd,rfveg,road,roof,room,   &
                         slab,walle,wallw,ufull)

use uclem_parameters
use uclem_types

implicit none

integer, intent(in) :: ufull
real, intent(in) :: ddt
real, dimension(ufull), intent(in) :: sg_roof,rg_roof,fg_roof,sg_walle,rg_walle,fg_walle
real, dimension(ufull), intent(in) :: sg_road,rg_road,fg_road,sg_wallw,rg_wallw,fg_wallw
real, dimension(ufull), intent(in) :: rgint_roof,rgint_walle,rgint_wallw,rgint_slab
real, dimension(ufull), intent(in) :: eg_roof,eg_road,garfsn,gardsn,d_rfsndelta,d_rdsndelta
real, dimension(ufull), intent(in) :: a_sg,a_rg,u_ts,u_fg,u_eg,u_alb,u_melt
real, dimension(ufull), intent(in) :: ggint_roof,ggint_road,ggint_walle,ggint_wallw
real, dimension(ufull), intent(in) :: ggext_intm,ggint_slab,ggint_intm,d_intgains_bld
real, dimension(ufull), intent(in) :: int_infilflux,d_ac_inside
real(kind=8), dimension(ufull) :: d_roofflux,d_walleflux,d_wallwflux,d_roadflux,d_slabflux,d_intmflux,d_roomflux 
real(kind=8), dimension(ufull) :: d_roofstor,d_wallestor,d_wallwstor,d_roadstor,d_slabstor,d_intmstor,d_roomstor
real(kind=8), dimension(ufull) :: d_faceterr
real(kind=8), dimension(ufull) :: d_atmosflux
real(kind=8), dimension(ufull,nl) :: roadstorage_prev, roofstorage_prev, wallestorage_prev, wallwstorage_prev
real(kind=8), dimension(ufull,nl) :: slabstorage_prev, intmstorage_prev
real(kind=8), dimension(ufull,1) :: roomstorage_prev
type(facetparams), intent(in) :: fp_intm, fp_road, fp_roof, fp_slab, fp_wall
type(facetdata), intent(inout) :: intm
type(vegdata), intent(in) :: cnveg, rfveg
type(facetdata), intent(inout) :: road, roof, room, slab, walle, wallw
type(fparmdata), intent(in) :: fp
type(pdiagdata), intent(inout) :: pd

! Store previous calculation to determine flux
roofstorage_prev(:,:)  = roof%storage(:,:)
roadstorage_prev(:,:)  = road%storage(:,:)
wallestorage_prev(:,:) = walle%storage(:,:)
wallwstorage_prev(:,:) = wallw%storage(:,:)
slabstorage_prev(:,:)  = slab%storage(:,:)
intmstorage_prev(:,:)  = intm%storage(:,:)
roomstorage_prev(:,:)  = room%storage(:,:)
pd%surferr = 0.


! assume a_rho=1 for interior
room%storage(:,1) = real(fp%bldheight(:),8)*real(aircp,8)*room%nodetemp(:,1)
! Sum heat stored in urban materials from layer 1 to nl
! interface conduction
roof%storage(:,:)  = 0.5_8*real(fp_roof%depth(:,:),8)*real(fp_roof%volcp(:,:),8)                        & 
                        *(roof%nodetemp(:,0:nl-1)+roof%nodetemp(:,1:nl))
road%storage(:,:)  = 0.5_8*real(fp_road%depth(:,:),8)*real(fp_road%volcp(:,:),8)                        & 
                        *(road%nodetemp(:,0:nl-1)+road%nodetemp(:,1:nl))
walle%storage(:,:) = 0.5_8*real(fp_wall%depth(:,:),8)*real(fp_wall%volcp(:,:),8)                        & 
                        *(walle%nodetemp(:,0:nl-1)+walle%nodetemp(:,1:nl))
wallw%storage(:,:) = 0.5_8*real(fp_wall%depth(:,:),8)*real(fp_wall%volcp(:,:),8)                        & 
                        *(wallw%nodetemp(:,0:nl-1)+wallw%nodetemp(:,1:nl))
slab%storage(:,:)  = 0.5_8*real(fp_slab%depth(:,:),8)*real(fp_slab%volcp(:,:),8)                        & 
                        *(slab%nodetemp(:,0:nl-1)+slab%nodetemp(:,1:nl))
intm%storage(:,:)  = 0.5_8*real(fp_intm%depth(:,:),8)*real(fp_intm%volcp(:,:),8)                        & 
                        *(intm%nodetemp(:,0:nl-1)+intm%nodetemp(:,1:nl))

! first call do not calculate comparison with last timestep (return) 
if ( pd%first_call ) then
  pd%first_call = .false.
  return
end if  
  
d_roofstor = sum(roof%storage-roofstorage_prev,dim=2)/real(ddt,8)
d_roofflux = (1._8-real(d_rfsndelta,8))*(real(sg_roof,8)+real(rg_roof,8)-real(fg_roof,8)-real(eg_roof,8))  &
           + real(d_rfsndelta,8)*real(garfsn,8) - real(ggint_roof,8) - real(rgint_roof,8)
d_faceterr  = d_roofstor - d_roofflux
pd%surferr = pd%surferr + d_faceterr
if (any(abs(d_faceterr)>=energytol)) write(6,*) "aTEB roof facet closure error:", maxval(abs(d_faceterr))
d_roadstor = sum(road%storage-roadstorage_prev,dim=2)/real(ddt,8)
d_roadflux = (1._8-real(d_rdsndelta,8))*(real(sg_road,8)+real(rg_road,8)-real(fg_road,8)-real(eg_road,8)) &
           + real(d_rdsndelta,8)*real(gardsn,8) - real(ggint_road,8)
d_faceterr  = d_roadstor - d_roadflux
pd%surferr = pd%surferr + d_faceterr
if (any(abs(d_faceterr)>=energytol)) write(6,*) "aTEB road facet closure error:", maxval(abs(d_faceterr))
d_wallestor= sum(walle%storage-wallestorage_prev,dim=2)/real(ddt,8)
d_walleflux= real(sg_walle,8)+real(rg_walle,8)-real(fg_walle,8) - real(ggint_walle,8) - real(rgint_walle,8)
d_faceterr = d_wallestor - d_walleflux
pd%surferr = pd%surferr + d_faceterr
if (any(abs(d_faceterr)>=energytol)) write(6,*) "aTEB walle facet closure error:", maxval(abs(d_faceterr))
d_wallwstor= sum(wallw%storage-wallwstorage_prev,dim=2)/real(ddt,8)
d_wallwflux= real(sg_wallw,8)+real(rg_wallw,8)-real(fg_wallw,8) - real(ggint_wallw,8) - real(rgint_wallw,8)
d_faceterr = d_wallwstor - d_wallwflux
pd%surferr = pd%surferr + d_faceterr
if (any(abs(d_faceterr)>=energytol)) write(6,*) "aTEB wallw facet closure error:", maxval(abs(d_faceterr))
d_slabstor = sum(slab%storage-slabstorage_prev,dim=2)/real(ddt,8)
d_slabflux = -real(ggint_slab,8) - real(rgint_slab,8)
d_faceterr = d_slabstor - d_slabflux
pd%surferr = pd%surferr + d_faceterr
if (any(abs(d_faceterr)>=energytol)) write(6,*) "aTEB slab facet closure error:", maxval(abs(d_faceterr))
if (intairtmeth/=0) then
  d_intmstor = sum(intm%storage-intmstorage_prev,dim=2)/real(ddt,8)
  d_intmflux = real(ggext_intm,8) - real(ggint_intm,8)
  d_faceterr = d_intmstor - d_intmflux
  pd%surferr = pd%surferr + d_faceterr
  if (any(abs(d_faceterr)>=energytol)) write(6,*) "aTEB intm facet closure error:", maxval(abs(d_faceterr))
  d_roomstor = (room%storage(:,1)-roomstorage_prev(:,1))/real(ddt,8)
  d_roomflux = real(ggint_roof,8)+real(ggint_slab,8)-real(fp%intmassn,8)*d_intmflux                    & 
            + (real(fp%bldheight,8)/real(fp%bldwidth,8))*(real(ggint_walle,8) + real(ggint_wallw,8))   &
            + real(int_infilflux,8) + real(d_ac_inside,8) + real(d_intgains_bld,8)
  d_faceterr = d_roomstor - d_roomflux
  pd%surferr = pd%surferr + d_faceterr
  if (any(abs(d_faceterr)>=energytol)) write(6,*) "aTEB room volume closure error:", maxval(abs(d_faceterr))
else
  d_intmstor = 0._8
  d_roomstor = 0._8
end if

pd%storage_flux = d_roofstor*real(fp%sigmabld,8)*(1._8-real(rfveg%sigma,8))           &
                + d_roadstor*(1._8-real(fp%sigmabld,8))*(1._8-real(cnveg%sigma,8))    &
                + d_wallestor*(1._8-real(fp%sigmabld,8))*real(fp%hwratio,8)           &
                + d_wallwstor*(1._8-real(fp%sigmabld,8))*real(fp%hwratio,8)           &
                + d_slabstor*real(fp%sigmabld,8)                                      &
                + d_intmstor*real(fp%sigmabld,8)*real(fp%intmassn,8)                  &
                + d_roomstor*real(fp%sigmabld,8)

! atmosphere energy flux = (SWdown-SWup) + (LWdown-LWup) - Turbulent + Anthropogenic
d_atmosflux = (real(a_sg,8)-real(a_sg,8)*real(u_alb,8)) + (real(a_rg,8)-real(sbconst,8)*(real(u_ts,8)+urbtemp)**4) &
            - (real(u_fg,8)+real(u_eg,8)+real(u_melt,8)) + real(pd%bldheat,8) + real(pd%bldcool,8)                 & 
            + real(pd%traf,8) + real(fp%industryfg,8) + real(pd%intgains_full,8)
pd%atmoserr = pd%storage_flux - d_atmosflux

if ( any(abs(pd%atmoserr)>=energytol) ) then
  write(6,*) "UCLEM energy not conserved! Atmos. error:", maxval(abs(pd%atmoserr))
end if

return
end subroutine energyclosure

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update water prognostic variables for roads and roofs
                            
subroutine updatewater(ddt,surfwater,soilwater,leafwater,snow,den,alpha, &
                       snmelt,a_rnd,a_snd,eg_surf,eg_snow,d_tran,d_evap, &
                       d_totdepth,fp_vegrlai,iwbrelax,                   &
                       fp_sfc,fp_swilt,irrigwater,delintercept,          &
                       delsnow,ufull)

use uclem_parameters

implicit none

integer, intent(in) :: ufull
integer, intent(in) :: iwbrelax
real, intent(in) :: ddt
real, dimension(ufull), intent(inout) :: surfwater,soilwater,leafwater,snow,den,alpha
real, dimension(ufull), intent(in) :: snmelt,a_rnd,a_snd,eg_surf,eg_snow
real, dimension(ufull), intent(in) :: d_tran,d_evap,d_totdepth,fp_vegrlai
real, dimension(ufull) :: modrnd, oldleafwater
real, dimension(ufull), intent(in) :: fp_sfc, fp_swilt
real, dimension(ufull), intent(out):: irrigwater, delintercept, delsnow

modrnd = max(a_rnd-d_evap/lv-max(maxvwatf*fp_vegrlai-leafwater,0.)/ddt,0.) ! rainfall reaching the soil under vegetation

! soil moisture based on Kowalczyk et al 1994 eq 12
! note that since sigmaf=1, then there is no soil evaporation, only transpiration.
! Evaporation only occurs from water on leafs.
surfwater = surfwater+ddt*(a_rnd-eg_surf/lv+snmelt)                                    ! surface
soilwater = soilwater+ddt*(modrnd+snmelt*den/waterden-d_tran/lv)/(waterden*d_totdepth) ! soil
oldleafwater = leafwater - ddt*d_evap/lv
leafwater = leafwater+ddt*(a_rnd-d_evap/lv)                                            ! leaf
leafwater = min(max(leafwater,0.),maxvwatf*fp_vegrlai)
delintercept = leafwater - oldleafwater

if (iwbrelax==1) then
  ! increase soil moisture for irrigation
  irrigwater = max(0.75*fp_swilt+0.25*fp_sfc-soilwater,0.)/(86400./ddt+1.)    ! 24h e-fold time
  soilwater=soilwater+irrigwater
else
  irrigwater = 0.
end if

! snow fields
snow  = snow + ddt*(a_snd-eg_snow/ls-snmelt)
delsnow = ddt*(a_snd-eg_snow/ls-snmelt)
den   = den + (maxsnowden-den)/(86400./(0.24*ddt)+1.)
alpha = alpha + (minsnowalpha-alpha)/(86400./(0.24*ddt)+1.)

return
end subroutine updatewater

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Estimate saturation mixing ratio

subroutine getqsatr8(qsat,temp,ps)

implicit none
                  
real(kind=8), dimension(:), intent(in) :: temp
real, dimension(:), intent(in) :: ps
real, dimension(:), intent(out) :: qsat
real, dimension(size(temp)) :: tempr4

tempr4 = real(temp)
call getqsat(qsat,tempr4,ps)

return
end subroutine getqsatr8

subroutine getqsat(qsat,temp,ps)

use uclem_parameters

implicit none

real, dimension(:), intent(in) :: temp
real, dimension(:), intent(in) :: ps
real, dimension(:), intent(out) :: qsat
real, dimension(size(temp)) :: esatf,tdiff,rx
integer, dimension(size(temp)) :: ix

real, dimension(0:220), parameter :: table =                                          &
  (/ 1.e-9, 1.e-9, 2.e-9, 3.e-9, 4.e-9,                                               & !-146C
     6.e-9, 9.e-9, 13.e-9, 18.e-9, 26.e-9,                                            & !-141C
     36.e-9, 51.e-9, 71.e-9, 99.e-9, 136.e-9,                                         & !-136C
     0.000000188, 0.000000258, 0.000000352, 0.000000479, 0.000000648,                 & !-131C
     0.000000874, 0.000001173, 0.000001569, 0.000002090, 0.000002774,                 & !-126C
     0.000003667, 0.000004831, 0.000006340, 0.000008292, 0.00001081,                  & !-121C
     0.00001404, 0.00001817, 0.00002345, 0.00003016, 0.00003866,                      & !-116C
     0.00004942, 0.00006297, 0.00008001, 0.0001014, 0.0001280,                        & !-111C
     0.0001613, 0.0002026, 0.0002538, 0.0003170, 0.0003951,                           & !-106C
     0.0004910, 0.0006087, 0.0007528, 0.0009287, 0.001143,                            & !-101C 
     0.001403, 0.001719, 0.002101, 0.002561, 0.003117, 0.003784,                      & !-95C
     0.004584, 0.005542, 0.006685, 0.008049, 0.009672, 0.01160, 0.01388, 0.01658,     & !-87C
     0.01977, 0.02353, 0.02796, 0.03316, 0.03925, 0.04638, 0.05472, 0.06444, 0.07577, & !-78C
     0.08894, 0.1042, 0.1220, 0.1425, 0.1662, 0.1936, 0.2252, 0.2615, 0.3032,         & !-69C
     0.3511, 0.4060, 0.4688, 0.5406, 0.6225, 0.7159, 0.8223, 0.9432, 1.080,           & !-60C
     1.236, 1.413, 1.612, 1.838, 2.092, 2.380, 2.703, 3.067, 3.476,                   & !-51C
     3.935, 4.449, 5.026, 5.671, 6.393, 7.198, 8.097, 9.098,                          & !-43C
     10.21, 11.45, 12.83, 14.36, 16.06, 17.94, 20.02, 22.33, 24.88,                   & !-34C
     27.69, 30.79, 34.21, 37.98, 42.13, 46.69,51.70,57.20,63.23,69.85,                & !-24C 
     77.09, 85.02, 93.70, 103.20, 114.66, 127.20, 140.81, 155.67,                     & !-16C
     171.69, 189.03, 207.76, 227.96 , 249.67, 272.98, 298.00, 324.78,                 & !-8C
     353.41, 383.98, 416.48, 451.05, 487.69, 526.51, 567.52, 610.78,                  & !0C
     656.62, 705.47, 757.53, 812.94, 871.92, 934.65, 1001.3, 1072.2,                  & !8C
     1147.4, 1227.2, 1311.9, 1401.7, 1496.9, 1597.7, 1704.4, 1817.3,                  & !16C
     1936.7, 2063.0, 2196.4, 2337.3, 2486.1, 2643.0, 2808.6, 2983.1,                  & !24C
     3167.1, 3360.8, 3564.9, 3779.6, 4005.5, 4243.0, 4492.7, 4755.1,                  & !32C
     5030.7, 5320.0, 5623.6, 5942.2, 6276.2, 6626.4, 6993.4, 7377.7,                  & !40C
     7780.2, 8201.5, 8642.3, 9103.4, 9585.5, 10089.0, 10616.0,                        & !47C
     11166.0, 11740.0, 12340.0, 12965.0, 13617.0, 14298.0, 15007.0,                   & !54C
     15746.0, 16516.0, 17318.0, 18153.0, 19022.0, 19926.0, 20867.0,                   & !61C
     21845.0, 22861.0, 23918.0, 25016.0, 26156.0, 27340.0, 28570.0,                   & !68C
     29845.0, 31169.0 /)


tdiff=min(max( temp+(urbtemp-123.16), 0.), 219.)
rx=tdiff-aint(tdiff)
ix=int(tdiff)
esatf=(1.-rx)*table(ix)+ rx*table(ix+1)
qsat=0.622*esatf/max(ps-esatf,0.1)

return
end subroutine getqsat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! Interface for calculating ustar and thetastar

subroutine getinvres(invres,cd,z_on_l,olzoh,ilzom,zmin,sthetav,thetav,a_umag,mode)

use uclem_parameters

implicit none

integer, intent(in) :: mode
real, dimension(:), intent(in) :: ilzom
real, dimension(:), intent(in) :: zmin,sthetav,thetav
real, dimension(:), intent(in) :: a_umag
real, dimension(:), intent(out) :: invres,cd,z_on_l
real, dimension(:), intent(inout) :: olzoh
real, dimension(size(ilzom)) :: lna,thetavstar,integralh

lna=olzoh-ilzom
call dyerhicks(integralh,z_on_l,cd,thetavstar,thetav,sthetav,a_umag,zmin,ilzom,lna,mode)
invres=vkar*sqrt(cd)*a_umag/integralh
olzoh=lna+ilzom

return
end subroutine getinvres

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate stability functions using Dyer-Hicks

subroutine dyerhicks(integralh,z_on_l,cd,thetavstar,thetav,sthetav,umagin,zmin,ilzom,lna,mode)

use uclem_parameters

implicit none

integer, intent(in) :: mode
integer ic
real, dimension(:), intent(in) :: thetav
real, dimension(:), intent(in) :: sthetav,umagin,zmin,ilzom
real, dimension(:), intent(inout) :: lna
real, dimension(:), intent(out) :: cd,thetavstar
real, dimension(:), intent(out) :: integralh,z_on_l
real, dimension(size(thetav)) :: z0_on_l,zt_on_l,olzoh,umag
real, dimension(size(thetav)) :: pm0,ph0,pm1,ph1,integralm
!real, parameter :: aa1 = 3.8
!real, parameter :: bb1 = 0.5
!real, parameter :: cc1 = 0.3

umag = max(umagin, 0.01)
cd=(vkar/ilzom)**2                         ! first guess
call getlna(lna,cd,umag,zmin,ilzom,mode)
olzoh=ilzom+lna
integralh=sqrt(cd)*ilzom*olzoh/vkar        ! first guess
thetavstar=vkar*(thetav-sthetav)/integralh ! first guess

do ic = 1,icmax
  z_on_l=vkar*zmin*grav*thetavstar/((thetav+urbtemp)*cd*umag**2)
  z_on_l=min(z_on_l,10.)
  z0_on_l  = z_on_l*exp(-ilzom)
  zt_on_l  = z0_on_l*exp(-lna)
  where (z_on_l<0.)
    pm0     = (1.-16.*z0_on_l)**(-0.25)
    ph0     = (1.-16.*zt_on_l)**(-0.5)
    pm1     = (1.-16.*z_on_l)**(-0.25)
    ph1     = (1.-16.*z_on_l)**(-0.5)
    integralm = ilzom-2.*log((1.+1./pm1)/(1.+1./pm0))-log((1.+1./pm1**2)/(1.+1./pm0**2)) &
               +2.*(atan(1./pm1)-atan(1./pm0))
    integralh = olzoh-2.*log((1.+1./ph1)/(1.+1./ph0))
  elsewhere
    !--------------Beljaars and Holtslag (1991) momentum & heat
    pm0 = -(a_1*z0_on_l+b_1*(z0_on_l-(c_1/d_1))*exp(-d_1*z0_on_l)+b_1*c_1/d_1)
    pm1 = -(a_1*z_on_l+b_1*(z_on_l-(c_1/d_1))*exp(-d_1*z_on_l)+b_1*c_1/d_1)
    ph0 = -((1.+(2./3.)*a_1*zt_on_l)**1.5+b_1*(zt_on_l-(c_1/d_1))*exp(-d_1*zt_on_l)+b_1*c_1/d_1-1.)
    ph1 = -((1.+(2./3.)*a_1*z_on_l)**1.5+b_1*(z_on_l-(c_1/d_1))*exp(-d_1*z_on_l)+b_1*c_1/d_1-1.)
    integralm = ilzom-(pm1-pm0)
    integralh = olzoh-(ph1-ph0)
  endwhere
  integralm = max( integralm, 1.e-10 )
  integralh = max( integralh, 1.e-10 )
  !where (z_on_l<=0.4)
    cd = (max(0.01,min(vkar*umag/integralm,2.))/umag)**2
  !elsewhere
  !  cd = (max(0.01,min(vkar*umag/(aa1*( ( z_on_l**bb1)*(1.0+cc1* z_on_l**(1.-bb1)) &
  !      -(z0_on_l**bb1)*(1.+cc1*z0_on_l**(1.-bb1)) )),2.))/umag)**2
  !endwhere
  thetavstar= vkar*(thetav-sthetav)/integralh
  call getlna(lna,cd,umag,zmin,ilzom,mode)
end do

return
end subroutine dyerhicks

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Estimate roughness length for heat
!

pure subroutine getlna(lna,cd,umag,zmin,ilzom,mode)

implicit none

integer, intent(in) :: mode
integer iqu
real, dimension(:), intent(out) :: lna
real, dimension(:), intent(in) :: cd,umag,zmin,ilzom
real re
real, parameter :: nu = 1.461E-5
!real, parameter :: eta0 = 1.827E-5
!real, parameter :: t0 = 291.15
!real, parameter :: c = 120.
!eta=eta0*((t0+c)/(theta+c))*(theta/t0)**(2./3.)
!nu=eta/rho

select case(mode) ! roughness length for heat
  case(1) ! zot=zom/10.
    lna=2.3
  case(2) ! Kanda et al 2007
    do iqu = 1,size(lna)
      re=max(sqrt(cd(iqu))*umag(iqu)*zmin(iqu)*exp(-ilzom(iqu))/nu,10.)
      !lna(iqu)=2.46*re**0.25-2. !(Brutsaet, 1982)
      lna(iqu)=1.29*re**0.25-2.  !(Kanda et al, 2007)
    end do  
  case(3) ! zot=zom (neglect molecular diffusion)
    lna=0.
  case(4) ! user defined
    ! no change
end select

return
end subroutine getlna

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate shortwave radiation coefficients (modified to include 2nd wall)

subroutine getswcoeff(sg_roof,sg_vegr,sg_road,sg_walle,sg_wallw,sg_vegc,sg_rfsn,sg_rdsn,wallpsi,roadpsi,fp_hwratio,   &
                      fp_effhwratio,fp_vangle,fp_hangle,fp_fbeam,fp_vegsigmac,fp_roadalpha,fp_vegalphac,fp_wallalpha, &
                      ird_alpha,rdsndelta)

use uclem_parameters

implicit none

integer k
real, dimension(:), intent(in) :: rdsndelta
real, dimension(:), intent(in) :: ird_alpha
real, dimension(:), intent(out) :: wallpsi,roadpsi
real, dimension(:), intent(in) :: fp_hwratio, fp_effhwratio
real, dimension(:), intent(in) :: fp_vangle,fp_hangle,fp_fbeam,fp_vegsigmac,fp_roadalpha,fp_vegalphac
real, dimension(:), intent(in) :: fp_wallalpha
real, dimension(:), intent(out) :: sg_roof,sg_vegr,sg_road,sg_walle,sg_wallw,sg_vegc,sg_rfsn,sg_rdsn
real, dimension(size(rdsndelta)) :: thetazero,walles,wallws,roads,ta,tc,xa,ya,roadnetalpha
real, dimension(size(rdsndelta)) :: nwalles,nwallws,nroads

wallpsi=0.5*(fp_effhwratio+1.-sqrt(fp_effhwratio*fp_effhwratio+1.))/fp_effhwratio
roadpsi=sqrt(fp_effhwratio*fp_effhwratio+1.)-fp_effhwratio

! integrate through 180 deg instead of 360 deg.  Hence partitioning to east and west facing walls
where (fp_vangle>=0.5*pi)
  walles=(1.-fp_fbeam)*wallpsi
  wallws=fp_fbeam/fp_hwratio + (1.-fp_fbeam)*wallpsi
  roads=(1.-fp_fbeam)*roadpsi
elsewhere
  ta=tan(fp_vangle)
  thetazero=asin(1./max(fp_hwratio*ta,1.))
  tc=2.*(1.-cos(thetazero))
  xa=min(max(fp_hangle-thetazero,0.),pi)-max(fp_hangle-pi+thetazero,0.)-min(fp_hangle+thetazero,0.)
  ya=cos(max(min(0.,fp_hangle),fp_hangle-pi))-cos(max(min(thetazero,fp_hangle),fp_hangle-pi)) &
    +cos(min(0.,-fp_hangle))-cos(min(thetazero,-fp_hangle)) &
    +cos(max(0.,pi-fp_hangle))-cos(max(thetazero,pi-fp_hangle))
  ! note that these terms now include the azimuth angle
  walles=fp_fbeam*(xa/fp_hwratio+ta*ya)/pi+(1.-fp_fbeam)*wallpsi
  wallws=fp_fbeam*((pi-2.*thetazero-xa)/fp_hwratio+ta*(tc-ya))/pi+(1.-fp_fbeam)*wallpsi
  roads=fp_fbeam*(2.*thetazero-fp_hwratio*ta*tc)/pi+(1.-fp_fbeam)*roadpsi
end where

! Calculate short wave reflections
roadnetalpha=rdsndelta*ird_alpha+(1.-rdsndelta)*((1.-fp_vegsigmac)*fp_roadalpha+fp_vegsigmac*fp_vegalphac)
sg_walle = walles
sg_wallw = wallws
sg_road = roads
do k = 1,nrefl
  nwalles = roadnetalpha*wallpsi*roads + fp_wallalpha*(1.-2.*wallpsi)*wallws
  nwallws = roadnetalpha*wallpsi*roads + fp_wallalpha*(1.-2.*wallpsi)*walles
  nroads = fp_wallalpha*(1.-roadpsi)*0.5*(walles+wallws)
  walles = nwalles
  wallws = nwallws
  roads = nroads
  sg_walle = sg_walle + walles
  sg_wallw = sg_wallw + wallws
  sg_road = sg_road + roads
end do
sg_roof=1.
sg_vegr=1.
sg_rfsn=1.
sg_rdsn=sg_road
sg_vegc=sg_road

return
end subroutine getswcoeff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate longwave radiation coefficients (modified to include 2nd wall)

subroutine getlwcoeff(d_netemiss,d_cwa,d_cra,d_cw0,d_cww,d_crw,d_crr,d_cwr,d_rdsndelta,wallpsi,roadpsi,fp_vegsigmac, &
                      fp_roademiss,fp_vegemissc,fp_wallemiss)

use uclem_parameters

implicit none

integer k
real, dimension(:), intent(inout) :: d_netemiss
real, dimension(:), intent(inout) :: d_cwa,d_cra,d_cw0,d_cww,d_crw,d_crr,d_cwr,d_rdsndelta
real, dimension(:), intent(in) :: fp_vegsigmac,fp_roademiss,fp_vegemissc,fp_wallemiss
real, dimension(:), intent(in) :: wallpsi,roadpsi
real, dimension(size(d_netemiss)) :: rcwa,rcra,rcwe,rcww,rcrw,rcrr,rcwr
real, dimension(size(d_netemiss)) :: ncwa,ncra,ncwe,ncww,ncrw,ncrr,ncwr

d_netemiss=d_rdsndelta*snowemiss+(1.-d_rdsndelta)*((1.-fp_vegsigmac)*fp_roademiss+fp_vegsigmac*fp_vegemissc)
d_cwa=wallpsi
d_cra=roadpsi
d_cw0=0.
d_cww=1.-2.*wallpsi
d_crw=0.5*(1.-roadpsi)
d_crr=0.
d_cwr=wallpsi
rcwa=d_cwa
rcra=d_cra
rcwe=d_cw0
rcww=d_cww
rcrw=d_crw
rcrr=d_crr
rcwr=d_cwr
do k = 1,nrefl
  ncwa=(1.-d_netemiss)*wallpsi*rcra+(1.-fp_wallemiss)*(1.-2.*wallpsi)*rcwa
  ncra=(1.-fp_wallemiss)*(1.-roadpsi)*rcwa
  ncwe=(1.-d_netemiss)*wallpsi*rcrw+(1.-fp_wallemiss)*(1.-2.*wallpsi)*rcww
  ncww=(1.-d_netemiss)*wallpsi*rcrw+(1.-fp_wallemiss)*(1.-2.*wallpsi)*rcwe
  ncrw=(1.-fp_wallemiss)*(1.-roadpsi)*0.5*(rcww+rcwe)  
  ncwr=(1.-d_netemiss)*wallpsi*rcrr+(1.-fp_wallemiss)*(1.-2.*wallpsi)*rcwr
  ncrr=(1.-fp_wallemiss)*(1.-roadpsi)*rcwr
  rcwa=ncwa
  rcra=ncra
  rcwe=ncwe
  rcww=ncww
  rcrw=ncrw
  rcrr=ncrr
  rcwr=ncwr
  d_cwa=d_cwa+rcwa
  d_cra=d_cra+rcra
  d_cw0=d_cw0+rcwe
  d_cww=d_cww+rcww
  d_crw=d_crw+rcrw
  d_cwr=d_cwr+rcwr
  d_crr=d_crr+rcrr
end do

end subroutine getlwcoeff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! solve for road snow temperature (includes vegetation canopy temperature and canyon temperature)

subroutine solvecanyon(sg_road,rg_road,fg_road,eg_road,acond_road,abase_road,                          &
                       sg_walle,rg_walle,fg_walle,acond_walle,abase_walle,                             &
                       sg_wallw,rg_wallw,fg_wallw,acond_wallw,abase_wallw,                             &
                       sg_vegc,rg_vegc,fg_vegc,eg_vegc,acond_vegc,abase_vegc,                          &
                       sg_rdsn,rg_rdsn,fg_rdsn,eg_rdsn,acond_rdsn,abase_rdsn,rdsntemp,rdsnmelt,gardsn, &
                       a_umag,a_rho,a_rg,a_rnd,a_snd,d_we,d_ww,                                        &
                       d_canyontemp,d_canyonmix,d_tempc,d_mixrc,d_sigd,d_topu,d_netrad,                &
                       d_roaddelta,d_vegdeltac,d_rdsndelta,d_ac_canyon,d_traf,d_ac_inside,             &
                       d_canyonrgout,d_tranc,d_evapc,d_cwa,d_cra,d_cw0,d_cww,d_crw,d_crr,              &
                       d_cwr,d_totdepth,d_intgains_bld,d_roomstor,fgtop,egtop,evspsbltop,sbltop,       &
                       int_infilflux,int_infilfg,cyc_proportion,cyc_translation,ddt,                   &
                       cnveg,fp,fp_intm,fp_road,fp_roof,fp_wall,fp_slab,intm,pd,rdhyd,rfveg,road,      &
                       roof,room,slab,walle,wallw,iroomtemp,cvcoeff_roof,cvcoeff_walle,                &
                       cvcoeff_wallw,cvcoeff_slab,cvcoeff_intm1,cvcoeff_intm2,ufull,diag)

use uclem_parameters
use uclem_types

implicit none

integer, intent(in) :: ufull,diag
integer k,l
real, intent(in)    :: ddt
real, dimension(ufull), intent(inout) :: rg_road,fg_road,eg_road,abase_road
real, dimension(ufull), intent(inout) :: rg_walle,fg_walle,abase_walle
real, dimension(ufull), intent(inout) :: rg_wallw,fg_wallw,abase_wallw
real, dimension(ufull), intent(inout) :: rg_vegc,fg_vegc,eg_vegc,abase_vegc
real, dimension(ufull), intent(inout) :: rg_rdsn,fg_rdsn,eg_rdsn,abase_rdsn,rdsntemp,rdsnmelt,gardsn
real, dimension(ufull), intent(in) :: sg_road,sg_walle,sg_wallw,sg_vegc,sg_rdsn
real, dimension(ufull), intent(in) :: a_umag,a_rho,a_rg,a_rnd,a_snd,d_sigd,d_we,d_ww
real, dimension(ufull), intent(out) :: d_canyontemp,d_canyonmix,d_topu,d_netrad
real, dimension(ufull), intent(out) :: d_ac_inside
real, dimension(ufull), intent(inout) :: d_tempc,d_mixrc
real, dimension(ufull), intent(inout) :: d_roaddelta,d_vegdeltac,d_rdsndelta,d_ac_canyon,d_traf
real, dimension(ufull), intent(inout) :: d_canyonrgout,d_tranc,d_evapc,d_cwa,d_cra,d_cw0,d_cww,d_crw,d_crr,d_cwr
real, dimension(ufull), intent(inout) :: d_totdepth,d_intgains_bld
real, dimension(ufull), intent(out) :: fgtop,egtop,evspsbltop,sbltop,int_infilflux
real, dimension(ufull), intent(out) :: acond_road,acond_walle,acond_wallw,acond_vegc,acond_rdsn
real, dimension(ufull), intent(out) :: int_infilfg
real, dimension(ufull), intent(out) :: iroomtemp
real, dimension(ufull), intent(out) :: cvcoeff_roof, cvcoeff_walle, cvcoeff_wallw, cvcoeff_slab
real, dimension(ufull), intent(out) :: cvcoeff_intm1, cvcoeff_intm2
real, dimension(ufull), intent(out) :: d_roomstor
real, dimension(ufull), intent(in) :: cyc_proportion,cyc_translation
real, dimension(ufull) :: newval,sndepth,snlambda,ldratio,roadqsat,vegqsat,rdsnqsat
real, dimension(ufull) :: cu,topinvres,dts,dtt,cduv,z_on_l,dumroaddelta,dumvegdelta,res
real, dimension(ufull) :: aa,bb,cc,dd,ee,ff
real, dimension(ufull,2) :: evct,evctx,oldval
real, dimension(ufull) :: d_canyontemp_prev, fp_coeffbldheight
type(facetparams), intent(in) :: fp_intm, fp_road, fp_roof, fp_wall, fp_slab
type(hydrodata), intent(in) :: rdhyd
type(vegdata), intent(inout) :: cnveg, rfveg
type(facetdata), intent(inout) :: road, room, walle, wallw, roof, slab, intm
type(fparmdata), intent(in) :: fp
type(pdiagdata), intent(inout) :: pd


! -------Initialise out variables------------------------
d_ac_inside = 0.

! snow conductance
sndepth  = rdhyd%snow*waterden/rdhyd%den
snlambda = icelambda*(rdhyd%den/waterden)**1.88

! first guess for canyon and room air temperature 
! also guess for canyon veg, snow temperatures and water vapour mixing ratio
d_canyontemp      = d_tempc
d_canyontemp_prev = d_tempc
d_canyonmix       = d_mixrc
cnveg%temp        = d_tempc
rdsntemp          = real(road%nodetemp(:,1))
iroomtemp         = real(room%nodetemp(:,1))
rdsnmelt        = 0.
dumvegdelta     = 0. ! cray compiler bug

!  solve for aerodynamical resistance between canyon and atmosphere  
! assume zoh=zom when coupling to canyon air temperature
pd%lzoh = pd%lzom
dts    = d_canyontemp + (d_canyontemp+urbtemp)*0.61*d_canyonmix
dtt    = d_tempc + (d_tempc+urbtemp)*0.61*d_mixrc
call getinvres(topinvres,cduv,z_on_l,pd%lzoh,pd%lzom,pd%cndzmin,dts,dtt,a_umag,3)
call gettopu(d_topu,a_umag,z_on_l,fp%bldheight,cduv,pd%cndzmin,fp%hwratio,ufull)

select case(resmeth)
  case(0)
    acond_road  = (11.8+4.2*sqrt((d_topu*abase_road)**2+cduv*a_umag**2))/(aircp*a_rho)  ! From Rowley, et al (1930)
    acond_walle = acond_road
    acond_wallw = acond_road
    acond_rdsn  = acond_road
    acond_vegc  = acond_road
  case(2)  
    cu = abase_road*d_topu
    where (cu<=5.)
      acond_road = (6.15+4.18*cu)/(aircp*a_rho)
    elsewhere
      acond_road = (7.51*cu**0.78)/(aircp*a_rho)
    end where
    acond_walle = acond_road
    acond_wallw = acond_road
    acond_rdsn  = acond_road
    acond_vegc  = acond_road
  case default
    acond_road  = d_topu*abase_road  
    acond_walle = d_topu*abase_walle
    acond_wallw = d_topu*abase_wallw
    acond_vegc  = d_topu*abase_vegc
    acond_rdsn  = d_topu*abase_rdsn
end select

  
! solve for road snow and canyon veg temperatures -------------------------------
d_canyontemp_prev = d_canyontemp
ldratio  = 0.5*( sndepth/snlambda + fp_road%depth(:,1)/fp_road%lambda(:,1) )
cnveg%temp = cnveg%temp + 0.1 ! 1st guess
rdsntemp   = rdsntemp + 0.1   ! 1st guess
oldval(:,1) = cnveg%temp
oldval(:,2) = rdsntemp
call canyonflux(evct,eg_road,acond_road,roadqsat,dumroaddelta,                                &
                sg_vegc,rg_vegc,fg_vegc,eg_vegc,acond_vegc,vegqsat,res,dumvegdelta,           &
                sg_rdsn,rg_rdsn,fg_rdsn,eg_rdsn,acond_rdsn,rdsntemp,gardsn,rdsnmelt,rdsnqsat, &
                a_rg,a_rho,a_rnd,a_snd,d_canyontemp,d_canyontemp_prev,                        &
                d_canyonmix,d_sigd,d_netrad,d_tranc,d_evapc,                                  &
                d_cra,d_crr,d_crw,d_totdepth,d_vegdeltac,d_rdsndelta,d_roaddelta,d_mixrc,     &
                d_intgains_bld,d_ac_canyon,d_ac_inside,d_roomstor,d_we,d_ww,d_tempc,d_traf,   &
                topinvres,ldratio,cyc_translation,cyc_proportion,iroomtemp,                   &
                int_infilflux,int_infilfg,                                                    &
                cvcoeff_roof,cvcoeff_walle,cvcoeff_wallw,cvcoeff_slab,cvcoeff_intm1,          &
                cvcoeff_intm2,acond_walle,acond_wallw,                                        &
                ddt,cnveg,rfveg,fp,fp_wall,fp_road,fp_roof,fp_slab,fp_intm,rdhyd,             &
                road,walle,wallw,roof,slab,intm,room,                                         &
                pd,ufull,diag)
cnveg%temp = cnveg%temp - 0.1 ! 2nd guess
rdsntemp   = rdsntemp - 0.1   ! 2nd guess
do k = 1,nfgits ! sectant
  evctx = evct
  call canyonflux(evct,eg_road,acond_road,roadqsat,dumroaddelta,                                &
                  sg_vegc,rg_vegc,fg_vegc,eg_vegc,acond_vegc,vegqsat,res,dumvegdelta,           &
                  sg_rdsn,rg_rdsn,fg_rdsn,eg_rdsn,acond_rdsn,rdsntemp,gardsn,rdsnmelt,rdsnqsat, &
                  a_rg,a_rho,a_rnd,a_snd,d_canyontemp,d_canyontemp_prev,                        &
                  d_canyonmix,d_sigd,d_netrad,d_tranc,d_evapc,                                  &
                  d_cra,d_crr,d_crw,d_totdepth,d_vegdeltac,d_rdsndelta,d_roaddelta,d_mixrc,     &
                  d_intgains_bld,d_ac_canyon,d_ac_inside,d_roomstor,d_we,d_ww,d_tempc,d_traf,   &
                  topinvres,ldratio,cyc_translation,cyc_proportion,iroomtemp,                   &
                  int_infilflux,int_infilfg,acond_walle,acond_wallw,                            &
                  cvcoeff_roof,cvcoeff_walle,cvcoeff_wallw,cvcoeff_slab,cvcoeff_intm1,          &
                  cvcoeff_intm2,                                                                &
                  ddt,cnveg,rfveg,fp,fp_wall,fp_road,fp_roof,fp_slab,fp_intm,rdhyd,             &
                  road,walle,wallw,roof,slab,intm,room,                                         &
                  pd,ufull,diag)  
  evctx = evct - evctx
  where ( abs(evctx(:,1))>tol )
    newval      = max(min(cnveg%temp-evct(:,1)*(cnveg%temp-oldval(:,1))/evctx(:,1),400.-urbtemp),200.-urbtemp)
    oldval(:,1) = cnveg%temp
    cnveg%temp  = newval
  end where
  where ( abs(evctx(:,2))>tol )
    newval      = max(min(rdsntemp-evct(:,2)*(rdsntemp-oldval(:,2))/evctx(:,2), 300.-urbtemp),100.-urbtemp)
    oldval(:,2) = rdsntemp
    rdsntemp    = newval
  end where
end do ! k = 1,nfgits

room%nodetemp(:,1) = real(iroomtemp,8)

fp_coeffbldheight = fp%effhwratio/fp%hwratio

! solve for canyon sensible heat flux
fg_walle = aircp*a_rho*(real(walle%nodetemp(:,0))-d_canyontemp)*acond_walle*fp_coeffbldheight 
fg_wallw = aircp*a_rho*(real(wallw%nodetemp(:,0))-d_canyontemp)*acond_wallw*fp_coeffbldheight
fg_road  = aircp*a_rho*(real(road%nodetemp(:,0))-d_canyontemp)*acond_road
fg_vegc  = sg_vegc+rg_vegc-eg_vegc
fg_rdsn  = sg_rdsn+rg_rdsn-eg_rdsn-lf*rdsnmelt-gardsn*(1.-cnveg%sigma)
fgtop = fp%hwratio*(fg_walle+fg_wallw) + (1.-d_rdsndelta)*(1.-cnveg%sigma)*fg_road &
      + (1.-d_rdsndelta)*cnveg%sigma*fg_vegc + d_rdsndelta*fg_rdsn                 &
      + d_traf + d_ac_canyon - int_infilfg


! solve for canyon latent heat flux
egtop = (1.-d_rdsndelta)*(1.-cnveg%sigma)*eg_road + (1.-d_rdsndelta)*cnveg%sigma*eg_vegc &
      + d_rdsndelta*eg_rdsn
evspsbltop = (1.-d_rdsndelta)*(1.-cnveg%sigma)*eg_road/lv + (1.-d_rdsndelta)*cnveg%sigma*eg_vegc/lv &
      + d_rdsndelta*eg_rdsn/ls
sbltop = d_rdsndelta*eg_rdsn/ls


! calculate longwave radiation
rg_walle = fp_wall%emiss*(a_rg*d_cwa+sbconst*(real(walle%nodetemp(:,0))+urbtemp)**4*(fp_wall%emiss*d_cw0-1.)          & 
         + sbconst*(real(wallw%nodetemp(:,0))+urbtemp)**4*fp_wall%emiss*d_cww+d_netrad*d_cwr)
rg_wallw = fp_wall%emiss*(a_rg*d_cwa+sbconst*(real(wallw%nodetemp(:,0))+urbtemp)**4*(fp_wall%emiss*d_cw0-1.)          &
         + sbconst*(real(walle%nodetemp(:,0))+urbtemp)**4*fp_wall%emiss*d_cww+d_netrad*d_cwr)
rg_road = fp_road%emiss*(a_rg*d_cra+(d_netrad*d_crr-sbconst*(real(road%nodetemp(:,0))+urbtemp)**4)                    &
        + sbconst*fp_wall%emiss*((real(walle%nodetemp(:,0))+urbtemp)**4+(real(wallw%nodetemp(:,0))+urbtemp)**4)*d_crw)
! outgoing longwave radiation
! note that eff terms are used for outgoing longwave radiation, whereas rg terms are used for heat conduction
d_canyonrgout = a_rg-d_rdsndelta*rg_rdsn-(1.-d_rdsndelta)*((1.-cnveg%sigma)*rg_road+cnveg%sigma*rg_vegc)              &
              - fp%hwratio*(rg_walle+rg_wallw)

return
end subroutine solvecanyon

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! solve for canyon veg and snow fluxes
                     
subroutine canyonflux(evct,eg_road,acond_road,roadqsat,dumroaddelta,                                 &
                      sg_vegc,rg_vegc,fg_vegc,eg_vegc,acond_vegc,vegqsat,res,dumvegdelta,            &
                      sg_rdsn,rg_rdsn,fg_rdsn,eg_rdsn,acond_rdsn,rdsntemp,gardsn,rdsnmelt,rdsnqsat,  &
                      a_rg,a_rho,a_rnd,a_snd,d_canyontemp,d_canyontemp_prev,                         &
                      d_canyonmix,d_sigd,d_netrad,d_tranc,d_evapc,                                   &
                      d_cra,d_crr,d_crw,d_totdepth,d_vegdeltac,d_rdsndelta,d_roaddelta,d_mixrc,      &
                      d_intgains_bld,d_ac_canyon,d_ac_inside,d_roomstor,d_we,d_ww,d_tempc,d_traf,    &
                      topinvres,ldratio,cyc_translation,cyc_proportion,iroomtemp,                    &
                      int_infilflux,int_infilfg,                                                     &
                      cvcoeff_roof,cvcoeff_walle,cvcoeff_wallw,cvcoeff_slab,cvcoeff_intm1,           &
                      cvcoeff_intm2,acond_walle,acond_wallw,                                         &
                      ddt,cnveg,rfveg,fp,fp_wall,fp_road,fp_roof,fp_slab,fp_intm,rdhyd,              &
                      road,walle,wallw,roof,slab,intm,room,                                          &
                      pd,ufull,diag)

use uclem_parameters
use uclem_types

implicit none

integer, intent(in) :: ufull, diag
real, intent(in) :: ddt
real, dimension(ufull,2), intent(out) :: evct
real, dimension(ufull), intent(inout) :: eg_road,acond_road,roadqsat,dumroaddelta
real, dimension(ufull), intent(inout) :: rg_vegc,fg_vegc,eg_vegc,acond_vegc,vegqsat,res,dumvegdelta
real, dimension(ufull), intent(inout) :: rg_rdsn,fg_rdsn,eg_rdsn,acond_rdsn,rdsntemp,gardsn,rdsnmelt,rdsnqsat
real, dimension(ufull), intent(in) :: sg_vegc,sg_rdsn
real, dimension(ufull), intent(in) :: a_rg,a_rho,a_rnd,a_snd,d_sigd
real, dimension(ufull), intent(in) :: topinvres,ldratio,cyc_translation,cyc_proportion
real, dimension(ufull), intent(out) :: int_infilflux,int_infilfg
real, dimension(ufull), intent(inout) :: d_canyontemp,d_canyontemp_prev,d_canyonmix
real, dimension(ufull), intent(inout) :: d_netrad,d_tranc,d_evapc
real, dimension(ufull), intent(inout) :: d_cra,d_crr,d_crw,d_totdepth,d_vegdeltac,d_rdsndelta
real, dimension(ufull), intent(inout) :: d_roaddelta,d_mixrc,d_intgains_bld,d_ac_canyon,d_ac_inside
real, dimension(ufull), intent(inout) :: d_roomstor,d_tempc,d_traf
real, dimension(ufull), intent(in) :: d_we,d_ww
real, dimension(ufull), intent(inout) :: iroomtemp
real, dimension(ufull), intent(out) :: cvcoeff_roof, cvcoeff_walle, cvcoeff_wallw, cvcoeff_slab
real, dimension(ufull), intent(out) :: cvcoeff_intm1, cvcoeff_intm2
real, dimension(ufull), intent(in) :: acond_walle, acond_wallw
real, dimension(ufull) :: f1,f2,f3,f4,aa,bb,cc,dd,ee,ff
real, dimension(ufull) :: snevap,canyontemp_ave
type(vegdata), intent(in) :: cnveg, rfveg
type(facetparams), intent(in) :: fp_wall, fp_road, fp_roof, fp_slab, fp_intm
type(hydrodata), intent(in) :: rdhyd
type(facetdata), intent(inout) :: road, walle, wallw, roof, slab, intm, room
type(fparmdata), intent(in) :: fp
type(pdiagdata), intent(inout) :: pd

! calculate average canyon temperature over iterations for fluxes
canyontemp_ave = infilalpha*d_canyontemp + (1.-infilalpha)*d_canyontemp_prev

! saturated mixing ratio for road, vegetation and snow
call getqsatr8(roadqsat,road%nodetemp(:,0),d_sigd)
call getqsat(vegqsat,cnveg%temp,d_sigd)
call getqsat(rdsnqsat,rdsntemp,d_sigd)

! correction for dew
where (roadqsat < d_canyonmix)
  dumroaddelta = 1.
elsewhere
  dumroaddelta = d_roaddelta
end where
where (vegqsat < d_canyonmix)
  dumvegdelta = 1.
elsewhere
  dumvegdelta = d_vegdeltac
end where

! calculate canyon road latent heat flux
aa = rdhyd%surfwater/ddt + a_rnd + rdsnmelt
eg_road = lv*min( a_rho*d_roaddelta*(roadqsat-d_canyonmix)*acond_road, aa )
  
! Calculate longwave radiation emitted from the canyon floor
! MJT notes - This could be included within the iterative solver for snow and vegetation temperatures.
! However, it creates a (weak) coupling between these two variables and therefore could require
! a multivariate root finding method (e.g,. Broyden's method). Instead we explicitly solve for d_netrad, 
! which allows us to decouple the solutions for snow and vegetation temperatures.
d_netrad = sbconst*(d_rdsndelta*snowemiss*(rdsntemp+urbtemp)**4                                   &
         + (1.-d_rdsndelta)*((1.-cnveg%sigma)*fp_road%emiss*(real(road%nodetemp(:,0))+urbtemp)**4 &
         + cnveg%sigma*cnveg%emiss*(cnveg%temp+urbtemp)**4))
  
! vegetation transpiration terms (developed by Eva in CCAM sflux.f and CSIRO9)
where (cnveg%zo < 0.5)
  ff = 1.1*sg_vegc/max(cnveg%lai*150.,1.E-8)
elsewhere
  ff = 1.1*sg_vegc/max(cnveg%lai*30.,1.E-8)
end where
f1 = (1.+ff)/(ff+cnveg%rsmin*cnveg%lai/5000.)
f2 = max(0.5*(fp%sfc-fp%swilt)/max(rdhyd%soilwater-fp%swilt,1.E-9),1.)
f3 = max(1.-0.00025*(vegqsat-d_canyonmix)*d_sigd/0.622,0.5) ! increased limit from 0.05 to 0.5 following Mk3.6    
f4 = max(1.-0.0016*(298.-urbtemp-canyontemp_ave)**2,0.05)     ! 0.2 in Mk3.6
res = max(30.,cnveg%rsmin*f1*f2/(f3*f4))

! solve for vegetation and snow sensible heat fluxes
fg_vegc = aircp*a_rho*(cnveg%temp-canyontemp_ave)*acond_vegc
fg_rdsn = aircp*a_rho*(rdsntemp-canyontemp_ave)*acond_rdsn

! calculate longwave radiation for vegetation and snow
rg_vegc = cnveg%emiss*(a_rg*d_cra+(d_netrad*d_crr-sbconst*(cnveg%temp+urbtemp)**4)    &
                    +sbconst*fp_wall%emiss*((real(walle%nodetemp(:,0))+urbtemp)**4    &
                    +(real(wallw%nodetemp(:,0))+urbtemp)**4)*d_crw)
rg_rdsn = snowemiss*(a_rg*d_cra+(d_netrad*d_crr-sbconst*(rdsntemp+urbtemp)**4)        &
                    +sbconst*fp_wall%emiss*((real(walle%nodetemp(:,0))+urbtemp)**4    &
                    +(real(wallw%nodetemp(:,0))+urbtemp)**4)*d_crw)

! estimate snow melt
rdsnmelt = min(max(0.,rdsntemp+(urbtemp-273.16))*icecp*rdhyd%snow/(ddt*lf),rdhyd%snow/ddt)

! calculate latent heat of transpiration and evaporation of in-canyon vegetation
d_tranc = lv*min(max((1.-dumvegdelta)*a_rho*(vegqsat-d_canyonmix)/(1./max(acond_vegc,1.e-10)+res),0.), &
                 max((rdhyd%soilwater-fp%swilt)*d_totdepth*waterden/ddt,0.))
d_evapc = lv*min(dumvegdelta*a_rho*(vegqsat-d_canyonmix)*acond_vegc,rdhyd%leafwater/ddt+a_rnd)
eg_vegc = d_evapc + d_tranc

! calculate canyon snow latent heat and ground fluxes
snevap = min(a_rho*max(0.,rdsnqsat-d_canyonmix)*acond_rdsn,rdhyd%snow/ddt+a_snd-rdsnmelt)
eg_rdsn = lv*snevap
rdsnmelt = rdsnmelt+snevap
gardsn = (rdsntemp-real(road%nodetemp(:,0)))/ldratio ! use road temperature to represent canyon bottom surface temperature
                                                     ! (i.e., we have omitted soil under vegetation temperature)

! vegetation energy budget error term
evct(:,1) = sg_vegc + rg_vegc - fg_vegc - eg_vegc

! road snow energy balance error term
evct(:,2) = sg_rdsn + rg_rdsn - fg_rdsn - eg_rdsn - lf*rdsnmelt - gardsn*(1.-cnveg%sigma)

!!!!!!!!!!!!!!!!!!!! start interior models !!!!!!!!!!!!!!!!!!!!!!!!!!
select case(intairtmeth)
  case(0) ! simple internal physics (per Thatcher and Hurley, 2012)
    call interiorflux_simple(ufull,diag,ddt,                                    &
                             d_intgains_bld,cyc_translation,                    &
                             canyontemp_ave,d_ac_canyon,d_ac_inside,d_roomstor, &
                             iroomtemp,int_infilflux,int_infilfg,               &
                             cvcoeff_roof,cvcoeff_walle,                        &
                             cvcoeff_wallw,cvcoeff_slab,cvcoeff_intm1,          &
                             cvcoeff_intm2,fp,fp_roof,fp_wall,rfveg,roof,       &
                             walle,wallw,slab,intm)

  case(1,2) ! complex internal physics (per Lipson et al., 2018) or with interior frac      
    call interiorflux_complex(ufull,diag,ddt,a_rho,d_we,d_ww,                  &
                              d_intgains_bld,cyc_proportion,canyontemp_ave,    &
                              d_ac_canyon,d_ac_inside,d_roomstor,              &
                              iroomtemp,int_infilflux,int_infilfg,             &
                              cvcoeff_roof,cvcoeff_walle,cvcoeff_wallw,        &
                              cvcoeff_slab,cvcoeff_intm1,cvcoeff_intm2,fp,     &
                              fp_roof,fp_wall,fp_slab,fp_intm,roof,walle,      &
                              wallw,slab,intm,room,pd)
  case default
    write(6,*) "ERROR: Unknown option intairtmeth = ",intairtmeth
    stop
end select

! balance canyon latent heat budget
aa = d_rdsndelta*acond_rdsn
bb = (1.-d_rdsndelta)*(1.-cnveg%sigma)*dumroaddelta*acond_road
cc = (1.-d_rdsndelta)*cnveg%sigma*(dumvegdelta*acond_vegc+(1.-dumvegdelta)/(1./max(acond_vegc,1.e-10)+res))
dd = topinvres
d_canyonmix = (aa*rdsnqsat+bb*roadqsat+cc*vegqsat+dd*d_mixrc)/(aa+bb+cc+dd)
  
! update canyon temperature estimate
d_canyontemp_prev = d_canyontemp  
aa = aircp*a_rho*topinvres
bb = d_rdsndelta*aircp*a_rho*acond_rdsn
cc = (1.-d_rdsndelta)*(1.-cnveg%sigma)*aircp*a_rho*acond_road
dd = (1.-d_rdsndelta)*cnveg%sigma*aircp*a_rho*acond_vegc
ee = fp%effhwratio*aircp*a_rho*acond_walle
ff = fp%effhwratio*aircp*a_rho*acond_wallw
d_canyontemp = (aa*d_tempc+bb*rdsntemp+cc*real(road%nodetemp(:,0))+dd*cnveg%temp+ee*real(walle%nodetemp(:,0)) & 
              +ff*real(wallw%nodetemp(:,0))+d_traf+d_ac_canyon-int_infilfg)/(aa+bb+cc+dd+ee+ff)

return
end subroutine canyonflux

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! solve for canyon/ internal heat fluxes and air temperatures

subroutine interiorflux_simple(ufull,diag,ddt,                                        &
                               d_intgains_bld,cyc_translation,                        &
                               d_canyontemp,d_ac_canyon,d_ac_inside,d_roomstor,       &
                               iroomtemp,int_infilflux,int_infilfg,                   &
                               cvcoeff_roof,cvcoeff_walle,cvcoeff_wallw,              &
                               cvcoeff_slab,cvcoeff_intm1,cvcoeff_intm2,fp,fp_roof,   &
                               fp_wall,rfveg,roof,walle,wallw,slab,intm)

use uclem_parameters
use uclem_types

implicit none

!START------Initialise variables -------------------------------------------------
!--+-+------Canyon/internal air temperature loop (default ncyits=6) -----------+-+
!  | |                                                                         | |
!  | |        (intairtmeth=0) Evaluating with fixed internal air temperature   | |
!  | +---0.1: Estimate internal surface convection coefficients                | |
!  | +---0.2: Calculate internal conducted flux (ggint)                        | |
!  | +---0.3: Calculate cooling/heating fluxes (fixed air temp)                | |
!  | +---0.4: Calculate acflux into canyon ------------------------------------+ |
!---------------------------------------------------------------------------------

! passed variables
integer,                    intent(in) :: diag
integer,                    intent(in) :: ufull
real,                       intent(in) :: ddt
real, dimension(ufull),     intent(in) :: d_intgains_bld
real, dimension(ufull),     intent(in) :: cyc_translation
real, dimension(ufull),     intent(in) :: d_canyontemp
real, dimension(ufull),     intent(out) :: d_ac_canyon
real, dimension(ufull),     intent(inout) :: d_ac_inside
real, dimension(ufull),     intent(inout) :: iroomtemp
real, dimension(ufull),     intent(out) :: int_infilflux
real, dimension(ufull),     intent(out) :: int_infilfg
real, dimension(ufull),     intent(out) :: cvcoeff_roof, cvcoeff_walle, cvcoeff_wallw, cvcoeff_slab
real, dimension(ufull),     intent(out) :: cvcoeff_intm1, cvcoeff_intm2
real, dimension(ufull),     intent(out) :: d_roomstor

type(fparmdata),   intent(in) :: fp
type(facetparams), intent(in) :: fp_roof,fp_wall
type(facetdata),   intent(in) :: roof,walle,wallw,slab,intm
type(vegdata),     intent(in) :: rfveg

! local variables
real, dimension(ufull) :: ggint_roof,ggint_walle,ggint_wallw,ggint_slab
real, dimension(ufull) :: ac_coeff

if ( diag>=1 ) write(6,*) "Evaluating with fixed internal air temperature"

iroomtemp = fp%bldairtemp + cyc_translation
d_roomstor = 0.

! ---0.1: Estimate internal surface convection coefficients -------------
call calc_convcoeff(cvcoeff_roof,cvcoeff_walle,cvcoeff_wallw,cvcoeff_slab,  & 
                    cvcoeff_intm1,cvcoeff_intm2,iroomtemp,roof,slab,intm,ufull)
! ---0.2: Calculate internal conducted flux (ggint) ----------------------
ggint_slab=0.  ! no slab for intairtmeth=0
ggint_roof  = cvcoeff_roof*(real(roof%nodetemp(:,nl))-iroomtemp)
ggint_walle = cvcoeff_walle*(real(walle%nodetemp(:,nl))-iroomtemp)
ggint_wallw = cvcoeff_wallw*(real(wallw%nodetemp(:,nl))-iroomtemp)

! ---0.3: Calculate cooling/heating fluxes (fixed air temp) --------------
! flux into room potentially pumped out into canyon (depends on AC method)
d_ac_inside = -ggint_roof - ggint_slab                         & 
              - (ggint_walle+ggint_wallw)*(fp%bldheight/fp%bldwidth)

! ---0.4: Calculate acflux into canyon -----------------------------------
! update heat pumped into canyon
where ( d_ac_inside<0. )
  ! original UCLEM COP calculation (Lipson et al., 2018)
  ! ac_coeff = max(acfactor*(d_canyontemp-iroomtemp)/(iroomtemp+urbtemp), 1./ac_copmax) ! T&H2012 Eq. 10
  ! modified UCLEM COP calculation (based on "Coefficient of performance for finite speed heat pump" - Blanchard 1980)
  ac_coeff = max( acfactor*(d_canyontemp+30.-iroomtemp)/(iroomtemp+30.+urbtemp), 0. )
elsewhere 
  ac_coeff = 0.
end where
select case(acmeth) ! AC heat pump into canyon (0=Off, 1=On, 2=Reversible, COP of 1.0)
  case(0) ! unrealistic cooling (buildings act as heat sink)
    d_ac_canyon  = 0.
  case(1) ! d_ac_canyon pumps conducted heat + ac waste heat back into canyon
    d_ac_canyon = max(0.,-d_ac_inside*(1.+ac_coeff)*fp%sigmabld/(1.-fp%sigmabld))  ! canyon domain W/m/m
  case(2) ! reversible heating and cooling (for testing energy conservation)
    d_ac_canyon  = -d_ac_inside*fp%sigmabld/(1.-fp%sigmabld)
  case DEFAULT
    write(6,*) "ERROR: Unknown acmeth mode ",acmeth
    stop
end select

! immediately eject intgains as infiltration flux
int_infilflux = -d_intgains_bld  
int_infilfg = int_infilflux*fp%sigmabld/(1.-fp%sigmabld)

return
end subroutine interiorflux_simple

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! solve for canyon/ internal heat fluxes and air temperatures

subroutine interiorflux_complex(ufull,diag,ddt,a_rho,d_weu,d_wwu,                      &
                                d_intgains_bld,cyc_proportion,canyontemp_ave,          &
                                d_ac_canyon,d_ac_inside,d_roomstor,                    &
                                iroomtemp,int_infilflux,int_infilfg,                   &
                                cvcoeff_roof,cvcoeff_walle,cvcoeff_wallw,cvcoeff_slab, &
                                cvcoeff_intm1,cvcoeff_intm2,fp,fp_roof,fp_wall,        &
                                fp_slab,fp_intm,roof,walle,wallw,slab,intm,room,pd)

use uclem_parameters
use uclem_types

implicit none

!START------Initialise variables -------------------------------------------------
!--+-+------Canyon/internal air temperature loop (default ncyits=6) -----------+-+
!  |          (intairtmeth=1) Evaluating with varying internal air temperature   |
!  +-----1.1: Estimate internal surface convection coefficients                  |
!  +-----1.2: Calculate air exchange from infiltration and open windows          |
!  +-----1.3: Calculate cooling/heating fluxes (varying air temp)                |
!  +-----1.4: Calculate internal air temperature implicitly                      |
!  +--+         Final loop                                                       |
!  |  +-----1.5: Final iteration flux and temperature                            |
!  +-----1.8: Calculate acflux into canyon and update canyon temperature --------+
!---------------------------------------------------------------------------------

! passed variables
integer,                    intent(in) :: diag                    ! diagnostic integer
integer,                    intent(in) :: ufull                   ! tile vector length
real,                       intent(in) :: ddt                     ! timestep
real, dimension(ufull),     intent(in) :: a_rho                   ! atmospheric air density
real, dimension(ufull),     intent(in) :: d_weu,d_wwu             ! wind speed wall east and west
real, dimension(ufull),     intent(in) :: d_intgains_bld          ! internal gain flux for building
real, dimension(ufull),     intent(in) :: cyc_proportion          ! cyclic proportion of conditioning units
real, dimension(ufull),     intent(in) :: canyontemp_ave          ! canyon air temperature (diagnostic)
real, dimension(ufull),     intent(out) :: d_ac_canyon            ! conditioning flux into canyon
real, dimension(ufull),     intent(inout) :: d_ac_inside          ! conditioning flux inside
real, dimension(ufull),     intent(inout) :: iroomtemp            ! update of room air temperature
real, dimension(ufull),     intent(out) :: int_infilflux          ! internal infiltration flux
real, dimension(ufull),     intent(out) :: int_infilfg            ! external infiltration flux
real, dimension(ufull),     intent(out) :: cvcoeff_roof,cvcoeff_walle,cvcoeff_wallw,cvcoeff_slab
real, dimension(ufull),     intent(out) :: cvcoeff_intm1,cvcoeff_intm2
real, dimension(ufull),     intent(out) :: d_roomstor

type(fparmdata),   intent(in) :: fp
type(facetparams), intent(in) :: fp_roof,fp_wall,fp_slab,fp_intm
type(facetdata),   intent(in) :: roof,walle,wallw,slab,intm
type(facetdata),   intent(inout) :: room
type(pdiagdata),   intent(in) :: pd

! local variables
integer iqu
real, dimension(ufull) :: ggint_roof, ggint_walle, ggint_wallw ! facet convection fluxes
real, dimension(ufull) :: ggint_slab, ggext_intm, ggint_intm   ! facet convection fluxes
real, dimension(ufull) :: rm,rf,we,ww,sl,im1,im2,infl
real, dimension(ufull) :: ac_coeff,d_ac_cool,d_ac_heat,ac_load
real, dimension(ufull) :: infl_dynamic,itemp
real, dimension(ufull) :: ac_heat_temp,ac_cool_temp
real, dimension(ufull) :: int_airden
real, dimension(ufull) :: d_openwindows
real, dimension(ufull) :: dummy, tempr4


if ( diag>=1 ) write(6,*) "Evaluating with varying internal air temperature"

! ---1.1: Estimate internal surface convection coefficients -------------
tempr4 = real(room%nodetemp(:,1))
call calc_convcoeff(cvcoeff_roof,cvcoeff_walle,cvcoeff_wallw,cvcoeff_slab,       & 
                    cvcoeff_intm1,cvcoeff_intm2,tempr4,roof,slab,intm,ufull)


ac_heat_temp = fp%bldairtemp - ac_deltat
ac_cool_temp = fp%bldairtemp + ac_deltat


! ---1.2: Calculate air exchange from infiltration and open windows -----
call calc_openwindows(d_openwindows,fp,iroomtemp,canyontemp_ave,roof,walle,wallw,slab,ufull)
select case(infilmeth)
  case(0) ! constant
    infl_dynamic = fp%infilach
  case(1) ! EnergyPlus with BLAST coefficients (Coblenz and Achenbach, 1963)
    infl_dynamic = fp%infilach*(0.606 + 0.03636*abs(iroomtemp-canyontemp_ave) + 0.1177*0.5*(d_wwu+d_weu) + 0.)
  !case(2) ! AccuRate (Chen, 2010)
  !  not yet implemented
  case default
    write(6,*) "ERROR: Unknown option infilmeth = ",infilmeth
    stop
end select


! ---1.3: Calculate cooling/heating fluxes (varying air temp)------------
! int_airden = air_pr/(air_rd*(iroomtemp+urbtemp))
int_airden = 1.
infl = aircp*a_rho*fp%bldheight*(infl_dynamic+d_openwindows*fp%ventilach)/3600.
rm = aircp*int_airden*fp%bldheight/ddt
rf = cvcoeff_roof
we = (fp%bldheight/fp%bldwidth)*cvcoeff_walle
ww = (fp%bldheight/fp%bldwidth)*cvcoeff_wallw
sl = cvcoeff_slab
im1 = cvcoeff_intm1*real(fp%intmassn)
im2 = cvcoeff_intm2*real(fp%intmassn)


! ---1.4: Calculate cooling/heating fluxes (varying air temp)------------

! remove acflux component added to implicit formulation of air temp
itemp = (rm*real(room%nodetemp(:,1))              & ! room temperature
       + rf*real(roof%nodetemp(:,nl))             & ! roof conduction
       + we*real(walle%nodetemp(:,nl))            & ! wall conduction east
       + ww*real(wallw%nodetemp(:,nl))            & ! wall conduction west
       + sl*real(slab%nodetemp(:,nl))             & ! slab conduction
       + im1*real(intm%nodetemp(:,0))             & ! mass conduction side 1
       + im2*real(intm%nodetemp(:,nl))            & ! mass conduction side 2
       + infl*canyontemp_ave                      & ! weighted infiltration
       + d_intgains_bld                           & ! internal gains (explicit)
       )/(rm+rf+we+ww+sl+im1+im2+infl)

d_ac_cool = 0.
d_ac_heat = 0.

! heating load
where ( itemp < ac_heat_temp )
  !ac_load = rm*(ac_heat_temp-itemp)
  ac_load = (rm+rf+we+ww+sl+im1+im2+infl)*(ac_heat_temp-itemp)
  d_ac_heat = min(ac_heatcap*fp%bldheight,ac_load)*pd%ac_heat_on
end where
! cooling load
where ( itemp > ac_cool_temp )
  !ac_load = rm*(itemp-ac_cool_temp)
  ac_load = (rm+rf+we+ww+sl+im1+im2+infl)*(itemp-ac_cool_temp)
  d_ac_cool = min(ac_coolcap*fp%bldheight,ac_load)*pd%ac_cool_on
end where

! modulate acfluxes to account for proportional heating/cooling area and time of day
if (intairtmeth==1) then
  d_ac_heat = d_ac_heat*fp%heatprop*cyc_proportion
  d_ac_cool = d_ac_cool*fp%coolprop*cyc_proportion
end if

d_ac_inside = d_ac_heat - d_ac_cool

! ---1.4: Calculate internal air temperature implicitly -----------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! fully implicit - old temps with acflux included

iroomtemp = (rm*real(room%nodetemp(:,1))              & ! room temperature
           + rf*real(roof%nodetemp(:,nl))             & ! roof conduction
           + we*real(walle%nodetemp(:,nl))            & ! wall conduction east
           + ww*real(wallw%nodetemp(:,nl))            & ! wall conduction west
           + sl*real(slab%nodetemp(:,nl))             & ! slab conduction
           + im1*real(intm%nodetemp(:,0))             & ! mass conduction side 1
           + im2*real(intm%nodetemp(:,nl))            & ! mass conduction side 2
           + infl*canyontemp_ave                      & ! weighted infiltration
           + d_ac_inside                              & ! ac flux (iterative)
           + d_intgains_bld                           & ! internal gains (explicit)
           )/(rm+rf+we+ww+sl+im1+im2+infl)
!iroomtemp = itemp + d_ac_inside/rm ! time split version


! calculate infiltration flux with dampening factor from previous timestep (reduce oscillation)
int_infilflux = infl*(canyontemp_ave - iroomtemp)
int_infilfg = (fp%sigmabld/(1.-fp%sigmabld))*int_infilflux

d_roomstor = (iroomtemp - real(room%nodetemp(:,1)))*(aircp*int_airden*fp%bldheight)/ddt

! ---1.5 (Depreciated by MJT)

! ---1.6: Calculate acflux into canyon and update canyon temperature-----

where ( d_ac_cool>0. )
  ! original UCLEM COP calculation (Lipson et al., 2018)
  ! ac_coeff = max(acfactor*(canyontemp_ave-iroomtemp)/(iroomtemp+urbtemp), 1./ac_copmax) ! T&H2012 Eq. 10
  ! modified UCLEM COP calculation (based on "Coefficient of performance for finite speed heat pump" - Blanchard 1980)
  ac_coeff = max( acfactor*(canyontemp_ave+30.-iroomtemp)/(iroomtemp+30.+urbtemp), 0. )
elsewhere 
  ac_coeff = 0.
end where
select case(acmeth) ! AC heat pump into canyon (0=Off, 1=On, 2=Reversible)
  case(0) ! unrealistic cooling (buildings act as heat sink)
    d_ac_canyon = 0.
  case(1) ! d_ac_canyon pumps conducted heat + ac waste heat back into canyon
    d_ac_canyon = max(0.,d_ac_cool*(1. + ac_coeff)*fp%sigmabld/(1.-fp%sigmabld))  ! canyon domain W/m2
  case(2) ! reversible heating and cooling (for testing energy conservation)
    d_ac_canyon = -d_ac_inside*fp%sigmabld/(1.-fp%sigmabld)
  case DEFAULT
    write(6,*) "ERROR: Unknown acmeth mode ",acmeth
    stop
end select

return
end subroutine interiorflux_complex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Solve for roof fluxes

subroutine solveroof(sg_rfsn,rg_rfsn,fg_rfsn,eg_rfsn,garfsn,rfsnmelt,rfsntemp,acond_rfsn,d_rfsndelta, &
                     sg_vegr,rg_vegr,fg_vegr,eg_vegr,acond_vegr,d_vegdeltar,                          &
                     sg_roof,rg_roof,eg_roof,acond_roof,d_roofdelta,                                  &
                     a_rg,a_umag,a_rho,a_rnd,a_snd,d_tempr,d_mixrr,d_rfdzmin,d_tranr,d_evapr,         &
                     d_sigr,ddt,fp_roof,rfhyd,rfveg,roof,fp,ufull)

use uclem_parameters
use uclem_types

implicit none

integer, intent(in) :: ufull
integer k
real, intent(in) :: ddt
real, dimension(ufull), intent(inout) :: rg_rfsn,fg_rfsn,eg_rfsn,garfsn,rfsnmelt,rfsntemp,acond_rfsn
real, dimension(ufull), intent(inout) :: rg_vegr,fg_vegr,eg_vegr,acond_vegr
real, dimension(ufull), intent(inout) :: rg_roof,eg_roof,acond_roof
real, dimension(ufull), intent(in) :: sg_rfsn,sg_vegr,sg_roof
real, dimension(ufull), intent(in) :: a_rg,a_umag,a_rho,a_rnd,a_snd
real, dimension(ufull), intent(inout) :: d_tempr,d_mixrr,d_rfdzmin,d_tranr,d_evapr,d_sigr
real, dimension(ufull), intent(inout) :: d_rfsndelta,d_vegdeltar,d_roofdelta
real, dimension(ufull) :: lzomroof,lzohroof,qsatr,dts,dtt,cdroof,z_on_l,newval,ldratio
real, dimension(ufull) :: aa,dd,ee
real, dimension(ufull,2) :: oldval,evctx,evctveg
type(facetparams), intent(in) :: fp_roof
type(hydrodata), intent(in) :: rfhyd
type(vegdata), intent(inout) :: rfveg
type(facetdata), intent(inout) :: roof
type(fparmdata), intent(in) :: fp

lzomroof=log(d_rfdzmin/zoroof)
lzohroof=2.3+lzomroof
call getqsatr8(qsatr,roof%nodetemp(:,0),d_sigr)
dts=real(roof%nodetemp(:,0)) + (real(roof%nodetemp(:,0))+urbtemp)*0.61*d_roofdelta*qsatr
dtt=d_tempr + (d_tempr+urbtemp)*0.61*d_mixrr
! Assume zot=0.1*zom (i.e., Kanda et al 2007, small experiment)
call getinvres(acond_roof,cdroof,z_on_l,lzohroof,lzomroof,d_rfdzmin,dts,dtt,a_umag,1)

! update green roof and snow temperature
rfveg%temp=d_tempr
rfsntemp  =max(min(real(roof%nodetemp(:,0)),300.-urbtemp),100.-urbtemp)
rg_vegr = fp_roof%emiss*(a_rg-sbconst*(real(roof%nodetemp(:,0))+urbtemp)**4) ! 1st guess
rg_rfsn = fp_roof%emiss*(a_rg-sbconst*(real(roof%nodetemp(:,0))+urbtemp)**4) ! 1st guess
eg_vegr = 0.
eg_rfsn = 0.
rfsnmelt = 0.
garfsn = 0.
d_tranr = 0.
d_evapr = 0.
acond_vegr = acond_roof
acond_rfsn = acond_roof
if ( any( d_rfsndelta>0. .or. rfveg%sigma>0. ) ) then
  evctveg = 0.
  oldval(:,1)=rfveg%temp+0.5
  oldval(:,2)=rfsntemp+0.5
  call roofflux(evctveg,rfsntemp,rfsnmelt,garfsn,sg_vegr,rg_vegr,fg_vegr,eg_vegr,acond_vegr, &
                sg_rfsn,rg_rfsn,fg_rfsn,eg_rfsn,acond_rfsn,a_rg,a_umag,a_rho,a_rnd,a_snd,    &
                d_tempr,d_mixrr,d_rfdzmin,d_tranr,d_evapr,d_sigr,d_vegdeltar,                &
                d_rfsndelta,ddt,fp,fp_roof,rfhyd,rfveg,roof,ufull)
  ! turn off roof snow and roof vegetation if they are not needed
  where ( rfveg%sigma>0. )
    rfveg%temp=rfveg%temp-0.5
  end where
  where ( d_rfsndelta>0. )
    rfsntemp  =rfsntemp-0.5
  end where
  do k=1,nfgits
    evctx=evctveg
    call roofflux(evctveg,rfsntemp,rfsnmelt,garfsn,sg_vegr,rg_vegr,fg_vegr,eg_vegr,acond_vegr, &
                  sg_rfsn,rg_rfsn,fg_rfsn,eg_rfsn,acond_rfsn,a_rg,a_umag,a_rho,a_rnd,a_snd,    &
                  d_tempr,d_mixrr,d_rfdzmin,d_tranr,d_evapr,d_sigr,d_vegdeltar,                &
                  d_rfsndelta,ddt,fp,fp_roof,rfhyd,rfveg,roof,ufull)
    evctx=evctveg-evctx
    where ( abs(evctx(:,1))>tol .and. rfveg%sigma>0. )
      newval=rfveg%temp-evctveg(:,1)*(rfveg%temp-oldval(:,1))/evctx(:,1)
      oldval(:,1)=rfveg%temp
      rfveg%temp=newval
    end where
    where ( abs(evctx(:,2))>tol .and. d_rfsndelta>0. )
      newval=max(min(rfsntemp-evctveg(:,2)*(rfsntemp-oldval(:,2))/evctx(:,2), 300.-urbtemp),100.-urbtemp)
      oldval(:,2)=rfsntemp
      rfsntemp=newval
    end where
  end do
end if
fg_vegr=sg_vegr+rg_vegr-eg_vegr
fg_rfsn=sg_rfsn+rg_rfsn-eg_rfsn-lf*rfsnmelt-garfsn*(1.-rfveg%sigma)

! estimate roof latent heat flux
where ( qsatr<d_mixrr )
  ! dew
  eg_roof=lv*a_rho*(qsatr-d_mixrr)*acond_roof
elsewhere
  ! evaporation
  aa=rfhyd%surfwater/ddt+a_rnd+rfsnmelt
  eg_roof=lv*min(a_rho*d_roofdelta*(qsatr-d_mixrr)*acond_roof,aa)
end where

! calculate net roof longwave radiation
! (sensible heat flux will be updated in solvetridiag)
rg_roof=fp_roof%emiss*(a_rg-sbconst*(real(roof%nodetemp(:,0))+urbtemp)**4)

return
end subroutine solveroof

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Solve for green roof and snow fluxes

subroutine roofflux(evct,rfsntemp,rfsnmelt,garfsn,sg_vegr,rg_vegr,fg_vegr,eg_vegr,acond_vegr,   &
                    sg_rfsn,rg_rfsn,fg_rfsn,eg_rfsn,acond_rfsn,a_rg,a_umag,a_rho,a_rnd,a_snd,   &
                    d_tempr,d_mixrr,d_rfdzmin,d_tranr,d_evapr,d_sigr,d_vegdeltar,               &
                    d_rfsndelta,ddt,fp,fp_roof,rfhyd,rfveg,roof,ufull)

use uclem_parameters
use uclem_types

implicit none

integer, intent(in) :: ufull
real, intent(in) :: ddt
real, dimension(ufull,2), intent(inout) :: evct
real, dimension(ufull), intent(in) :: rfsntemp,sg_vegr,sg_rfsn
real, dimension(ufull), intent(inout) :: rfsnmelt,garfsn
real, dimension(ufull), intent(inout) :: rg_vegr,fg_vegr,eg_vegr,acond_vegr
real, dimension(ufull), intent(inout) :: rg_rfsn,fg_rfsn,eg_rfsn,acond_rfsn
real, dimension(ufull), intent(in) :: a_rg,a_umag,a_rho,a_rnd,a_snd
real, dimension(ufull), intent(inout) :: d_tempr,d_mixrr,d_rfdzmin,d_tranr,d_evapr,d_sigr
real, dimension(ufull), intent(inout) :: d_vegdeltar,d_rfsndelta
real, dimension(ufull) :: lzomvegr,lzohvegr,vwetfac,dts,dtt,z_on_l,ff,f1,f2,f3,f4,cdvegr
real, dimension(ufull) :: vegqsat,dumvegdelta,res,sndepth,snlambda,ldratio,lzosnow,rfsnqsat,cdrfsn
real, dimension(ufull) :: lzotdum, snevap
type(facetparams), intent(in) :: fp_roof
type(hydrodata), intent(in) :: rfhyd
type(vegdata), intent(in) :: rfveg
type(facetdata), intent(in) :: roof
type(fparmdata), intent(in) :: fp

call getqsat(vegqsat,rfveg%temp,d_sigr)
where ( vegqsat<d_mixrr )
  dumvegdelta = 1.
elsewhere
  dumvegdelta = d_vegdeltar
end where

! transpiration terms (developed by Eva in CCAM sflux.f and CSIRO9)
where ( rfveg%zo<0.5 )
  ff = 1.1*sg_vegr/max(rfveg%lai*150.,1.E-8)
elsewhere
  ff = 1.1*sg_vegr/max(rfveg%lai*30.,1.E-8)
end where
f1 = (1.+ff)/(ff+rfveg%rsmin*rfveg%lai/5000.)
f2 = max(0.5*(fp%sfc-fp%swilt)/max(rfhyd%soilwater-fp%swilt,1.E-9),1.)
f3 = max(1.-.00025*(vegqsat-d_mixrr)*d_sigr/0.622,0.5)
f4 = max(1.-0.0016*((298.-urbtemp)-d_tempr)**2,0.05)
res = max(30.,rfveg%rsmin*f1*f2/(f3*f4))

vwetfac = max(min((rfhyd%soilwater-fp%swilt)/(fp%sfc-fp%swilt),1.),0.) ! veg wetfac (see sflux.f or cable_canopy.f90)
vwetfac = (1.-dumvegdelta)*vwetfac+dumvegdelta
lzomvegr = log(d_rfdzmin/rfveg%zo)
! xe is a dummy variable for lzohvegr
lzohvegr = 2.3+lzomvegr
dts = rfveg%temp + (rfveg%temp+urbtemp)*0.61*vegqsat*vwetfac
dtt = d_tempr + (d_tempr+urbtemp)*0.61*d_mixrr
! Assume zot=0.1*zom (i.e., Kanda et al 2007, small experiment)
call getinvres(acond_vegr,cdvegr,z_on_l,lzohvegr,lzomvegr,d_rfdzmin,dts,dtt,a_umag,1)
! acond_vegr is multiplied by a_umag

where ( rfveg%sigma>0. )
  ! longwave radiation    
  rg_vegr=rfveg%emiss*(a_rg-sbconst*(rfveg%temp+urbtemp)**4)
  
  ! sensible heat flux
  fg_vegr=aircp*a_rho*(rfveg%temp-d_tempr)*acond_vegr

  ! calculate transpiration and evaporation of roof vegetation
  d_tranr=lv*min(max((1.-dumvegdelta)*a_rho*(vegqsat-d_mixrr)/(1./acond_vegr+res),0.), &
                 max((rfhyd%soilwater-fp%swilt)*fp%rfvegdepth*waterden/ddt,0.))
  d_evapr=lv*min(dumvegdelta*a_rho*(vegqsat-d_mixrr)*acond_vegr,rfhyd%leafwater/ddt+a_rnd)
  eg_vegr=d_evapr+d_tranr
  
  ! balance green roof energy budget
  evct(:,1)=sg_vegr+rg_vegr-fg_vegr-eg_vegr
end where


! snow conductance
sndepth=rfhyd%snow*waterden/rfhyd%den
snlambda=icelambda*(rfhyd%den/waterden)**1.88
ldratio=0.5*(sndepth/snlambda+fp_roof%depth(:,1)/fp_roof%lambda(:,1))

! Update roof snow energy budget
lzosnow=log(d_rfdzmin/zosnow)
call getqsat(rfsnqsat,rfsntemp,d_sigr)
lzotdum=2.3+lzosnow
dts=rfsntemp + (rfsntemp+urbtemp)*0.61*rfsnqsat
call getinvres(acond_rfsn,cdrfsn,z_on_l,lzotdum,lzosnow,d_rfdzmin,dts,dtt,a_umag,1)
! acond_rfsn is multiplied by a_umag

where ( d_rfsndelta>0. )
  rfsnmelt=min(max(0.,rfsntemp+(urbtemp-273.16))*icecp*rfhyd%snow/(ddt*lf),rfhyd%snow/ddt)
  rg_rfsn=snowemiss*(a_rg-sbconst*(rfsntemp+urbtemp)**4)
  fg_rfsn=aircp*a_rho*(rfsntemp-d_tempr)*acond_rfsn
  snevap=min(a_rho*max(0.,rfsnqsat-d_mixrr)*acond_rfsn,rfhyd%snow/ddt+a_snd-rfsnmelt)
  eg_rfsn=lv*snevap
  rfsnmelt=rfsnmelt+snevap
  garfsn=(rfsntemp-real(roof%nodetemp(:,0)))/ldratio
  
  ! balance snow energy budget
  evct(:,2)=sg_rfsn+rg_rfsn-fg_rfsn-eg_rfsn-lf*rfsnmelt-garfsn*(1.-rfveg%sigma)
end where

return
end subroutine roofflux

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Define weights during the diurnal cycle. Array starts at 1am.

subroutine getdiurnal(fp_ctime,fp_weekdayload,icyc_traffic,icyc_basedemand,icyc_proportion,icyc_translation)                  

use uclem_parameters

implicit none

real, dimension(:), intent(in) :: fp_ctime,fp_weekdayload
real, dimension(:), intent(out) :: icyc_traffic,icyc_basedemand,icyc_proportion,icyc_translation
real, dimension(size(fp_ctime)) :: real_p
integer, dimension(size(fp_ctime)) :: int_p

! ! base temperature translation cycle optimised for Melbourne
! real, dimension(25), save :: trancycle =(/ -1.2, -1.7, -2.4, -2.9, -3.0, -2.1,  &
!                                            -1.1, -0.6, -0.7, -0.8, -1.2, -1.4,  &
!                                            -1.6, -1.8, -1.6, -1.0, -0.2,  1.1,  &
!                                             2.1,  2.7,  3.0,  2.6,  1.6,  0.2, -1.2 /)

int_p=int(24.*fp_ctime)
real_p=24.*fp_ctime-real(int_p)
where (int_p<1) int_p=int_p+24

icyc_traffic     = ((1.-real_p)*cyc_traf(int_p)+real_p*cyc_traf(int_p+1))*fp_weekdayload
icyc_basedemand  = ((1.-real_p)*cyc_base(int_p)+real_p*cyc_base(int_p+1))*fp_weekdayload
icyc_proportion  = ((1.-real_p)*cyc_prop(int_p)+real_p*cyc_prop(int_p+1))*fp_weekdayload
icyc_translation = ((1.-real_p)*cyc_tran(int_p)+real_p*cyc_tran(int_p+1))

! remove statistical energy use diurnal adjustments
if (statsmeth==0) then
  icyc_basedemand=1.
  icyc_proportion=1.
  icyc_translation=0.
end if

return
end subroutine getdiurnal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Define heating and cooling for different interior model tiles

subroutine uclem_prepinterior(ifrac,fp,pd,ufull)

use uclem_parameters
use uclem_types

implicit none

integer, intent(in) :: ifrac, ufull
type(fparmdata),   intent(in) :: fp
type(pdiagdata),   intent(inout) :: pd
real, dimension(ufull) :: cyc_traffic, cyc_basedemand, cyc_proportion, cyc_translation

! fractional interior spaces
call getdiurnal(fp%ctime,fp%weekdayload,cyc_traffic,cyc_basedemand,cyc_proportion,cyc_translation)
select case( intairtmeth + (ifrac-1)*10 )
  case(0,1)
    pd%frac_sigma = 1.
    pd%ac_heat_on = 1.
    pd%ac_cool_on = 1.
  case(2)  ! heating and cooling
    pd%frac_sigma = min( fp%heatprop, fp%coolprop )*cyc_proportion  
    pd%ac_heat_on = 1.
    pd%ac_cool_on = 1.
  case(12) ! heating or cooling
    pd%frac_sigma = (max( fp%heatprop, fp%coolprop ) - min( fp%heatprop, fp%coolprop ))*cyc_proportion
    where ( fp%heatprop>fp%coolprop )
      pd%ac_heat_on = 1.
      pd%ac_cool_on = 0.
    elsewhere
      pd%ac_heat_on = 0.
      pd%ac_cool_on = 1.
    end where      
  case(22) ! no conditioning
    pd%frac_sigma = 1. - max( fp%heatprop, fp%coolprop )*cyc_proportion
    pd%ac_heat_on = 0.
    pd%ac_cool_on = 0.
  case default
    write(6,*) "ERROR: Invalid choice for intairtemth+(ifrac-1)*10 ",intairtmeth + (ifrac-1)*10
    stop
end select

return 
end subroutine uclem_prepinterior   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate in-canyon wind speed for walls and road
! This version holds the eddy size in the canyon to a constant width
! that allows for an analytic solution to the integral

subroutine getincanwindb(ueast,uwest,ufloor,a_udir,z0,fp,ufull)

use uclem_parameters
use uclem_types

implicit none

integer, intent(in) :: ufull
real, dimension(ufull), intent(out) :: ueast,uwest,ufloor
real, dimension(ufull), intent(in) :: z0
real, dimension(ufull) :: wsuma,wsumb,fsum
real, dimension(ufull) :: theta1,wdir,h,w
real, dimension(ufull) :: dufa,dura,duva,ntheta
real, dimension(ufull) :: dufb,durb,duvb
real, dimension(ufull) :: fp_coeffbldheight
real, dimension(ufull), intent(in) :: a_udir
type(fparmdata), intent(in) :: fp

! rotate wind direction so that all cases are between 0 and pi
! walls are fliped at the end of the subroutine to account for additional pi rotation
where (a_udir>=0.)
  wdir=a_udir
elsewhere
  wdir=a_udir+pi
endwhere

fp_coeffbldheight = fp%effhwratio/fp%hwratio

h=fp%bldheight*fp_coeffbldheight
w=fp%bldheight/fp%hwratio
theta1=acos(min(w/(3.*h),1.))

call winda(dufa,dura,duva,h,w,z0,ufull) ! jet on road
call windb(dufb,durb,duvb,h,w,z0,ufull) ! jet on wall
ntheta=2. ! i.e., int_0^pi sin(theta) dtheta = 2.)
where (wdir<theta1.or.wdir>pi-theta1) ! jet on wall
  wsuma=duvb*ntheta
  wsumb=durb*ntheta
  fsum=dufb*ntheta
elsewhere                             ! jet on road
  wsuma=dura*ntheta
  wsumb=duva*ntheta
  fsum=dufa*ntheta
end where

! Correct for rotation of winds at start of subroutine
! 0.5 to adjust for factor of 2 in gettopu
where (a_udir>=0.)
  ueast=0.5*wsuma
  uwest=0.5*wsumb
elsewhere
  ueast=0.5*wsumb
  uwest=0.5*wsuma
end where
ufloor=0.5*fsum      ! floor

return
end subroutine getincanwindb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate canyon wind speeds, jet on road

subroutine winda(uf,ur,uv,h,w,z0,ufull)

use uclem_parameters

implicit none

integer, intent(in) :: ufull
integer iqu
real, dimension(ufull), intent(out) :: uf,ur,uv
real, dimension(ufull), intent(in) :: h,w,z0
real a,u0,cuven,zolog

u0=exp(-0.9*sqrt(13./4.))
do iqu = 1,ufull
  a=0.15*max(1.,3.*h(iqu)/(2.*w(iqu)))
  zolog=log(max(h(iqu),z0(iqu)+0.2)/z0(iqu))
  cuven=log(max(refheight*h(iqu),z0(iqu)+0.2)/z0(iqu))/log(max(h(iqu),z0(iqu)+0.2)/z0(iqu))
  cuven=max(cuven*max(1.-3.*h(iqu)/w(iqu),0.),(u0/a)*(h(iqu)/w(iqu))*(1.-exp(max(-a*max(w(iqu)/h(iqu)-3.,0.),-40.))))
  uf(iqu)=(u0/a)*(h(iqu)/w(iqu))*(1.-exp(-3.*a))+cuven
  !uf(iqu)=(u0/a)*(h(iqu)/w(iqu))*(2.-exp(-a*3.)-exp(-a*(w(iqu)/h(iqu)-3.)))
  ur(iqu)=(u0/a)*exp(-a*3.)*(1.-exp(-a))
  ! MJT suggestion
  cuven=1.-1./zolog
  uv(iqu)=(u0/a)*exp(max(-a*max(w(iqu)/h(iqu)-3.,0.),-40.))*(1.-exp(-a))
  uv(iqu)=max(cuven,uv(iqu))
  !uv(iqu)=(u0/a)*exp(-a*(w(iqu)/h(iqu)-3.))*(1.-exp(-a))
end do

return
end subroutine winda

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate canyon wind speeds, jet on wall

subroutine windb(uf,ur,uv,h,win,z0,ufull)

implicit none

integer, intent(in) :: ufull
integer iqu
real, dimension(ufull), intent(out) :: uf,ur,uv
real, dimension(ufull), intent(in) :: h,win,z0
real a,dh,u0,w
real zolog,cuven

do iqu = 1,ufull
  w=min(win(iqu),1.5*h(iqu))

  a=0.15*max(1.,3.*h(iqu)/(2.*w))
  dh=max(2.*w/3.-h(iqu),0.)
  u0=exp(-0.9*sqrt(13./4.)*dh/h(iqu))

  zolog=log(max(h(iqu),z0(iqu)+0.2)/z0(iqu))
  ! MJT suggestion (cuven is multiplied by dh to avoid divide by zero)
  cuven=h(iqu)-(h(iqu)-dh)*log(max(h(iqu)-dh,z0(iqu)+0.2)/z0(iqu))/zolog-dh/zolog
  ! MJT cuven is back to the correct units of m/s
  cuven=max(cuven/h(iqu),(u0/a)*(1.-exp(-a*dh/h(iqu))))

  uf(iqu)=(u0/a)*(h(iqu)/w)*exp(-a*(1.-dh/h(iqu)))*(1.-exp(-a*w/h(iqu)))
  ur(iqu)=(u0/a)*exp(-a*(1.-dh/h(iqu)+w/h(iqu)))*(1.-exp(-a))
  uv(iqu)=(u0/a)*(1.-exp(-a*(1.-dh/h(iqu))))+cuven
  !uv(iqu)=(u0/a)*(2.-exp(-a*(1.-dh/h(iqu)))-exp(-a*dh/h(iqu)))
end do

return
end subroutine windb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate wind speed at canyon top
subroutine gettopu(d_topu,a_umag,z_on_l,fp_bldheight,pd_cduv,pd_cndzmin,fp_hwratio,ufull)
      
use uclem_parameters

implicit none

integer, intent(in) :: ufull
real, dimension(ufull), intent(in) :: z_on_l
real, dimension(ufull) :: z0_on_l,bldheight
real, dimension(ufull) :: pm0,pm1,integralm
real, dimension(ufull) :: ustar,neutral
real, dimension(ufull), intent(out) :: d_topu
real, dimension(ufull), intent(in) :: a_umag
real, dimension(ufull), intent(in) :: fp_bldheight, fp_hwratio
real, dimension(ufull), intent(inout) :: pd_cduv, pd_cndzmin

bldheight=fp_bldheight*(1.-refheight)
ustar=sqrt(pd_cduv)*a_umag

z0_on_l=min(bldheight,pd_cndzmin)*z_on_l/pd_cndzmin ! calculate at canyon top
z0_on_l=min(z0_on_l,10.)
neutral = log(pd_cndzmin/min(bldheight,pd_cndzmin))
where (z_on_l<0.)
  pm0     = (1.-16.*z0_on_l)**(-0.25)
  pm1     = (1.-16.*z_on_l)**(-0.25)
  integralm = neutral-2.*log((1.+1./pm1)/(1.+1./pm0)) &
                -log((1.+1./pm1**2)/(1.+1./pm0**2)) &
                +2.*(atan(1./pm1)-atan(1./pm0))
elsewhere
  !-------Beljaars and Holtslag (1991) heat function
  pm0 = -(a_1*z0_on_l+b_1*(z0_on_l-(c_1/d_1))*exp(-d_1*z0_on_l)+b_1*c_1/d_1)
  pm1  = -(a_1*z_on_l+b_1*(z_on_l-(c_1/d_1))*exp(-d_1*z_on_l)+b_1*c_1/d_1)
  integralm = neutral-(pm1-pm0)
end where
where (bldheight<pd_cndzmin)
  d_topu=(2./pi)*(a_umag-ustar*integralm/vkar)
elsewhere ! within canyon
  d_topu=(2./pi)*a_umag*exp(0.5*fp_hwratio*(1.-pd_cndzmin/bldheight))
end where
d_topu=max(d_topu,0.1)

return
end subroutine gettopu

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate screen diagnostics

subroutine scrncalc(a_mixr,a_umag,a_temp,u_ts,d_tempc,d_rdsndelta,d_roaddelta,d_vegdeltac,d_sigd,smixr,rdsntemp,zonet, &
                    cnveg,fp,pd,rdhyd,road,ufull)

use uclem_parameters
use uclem_types

implicit none

integer, intent(in) :: ufull
real, dimension(ufull), intent(in) :: smixr,rdsntemp,zonet
real, dimension(ufull) :: cd,thetav,sthetav
real, dimension(ufull) :: thetavstar,z_on_l,z0_on_l
real, dimension(ufull) :: pm0,ph0,pm1,ph1,integralm,integralh
real, dimension(ufull) :: ustar,qstar,z10_on_l
real, dimension(ufull) :: neutral,neutral10,pm10
real, dimension(ufull) :: integralm10,tts,tetp
real, dimension(ufull) :: tstar,lna
real, dimension(ufull) :: utop,ttop,qtop,wf,tsurf,qsurf,n
real, dimension(ufull), intent(in) :: a_mixr,a_umag,a_temp
real, dimension(ufull), intent(in) :: u_ts
real, dimension(ufull), intent(in) :: d_tempc,d_rdsndelta,d_roaddelta,d_vegdeltac,d_sigd
type(vegdata), intent(in) :: cnveg
type(hydrodata), intent(in) :: rdhyd
type(facetdata), intent(in) :: road
type(fparmdata), intent(in) :: fp
type(pdiagdata), intent(inout) :: pd

real, parameter :: z0  = 1.5
real, parameter :: z10 = 10.

select case(scrnmeth)
  case(0) ! estimate screen diagnostics (slab at displacement height approach)
    thetav=d_tempc + (d_tempc+urbtemp)*0.61*a_mixr
    sthetav=u_ts + (u_ts+urbtemp)*0.61*smixr
    lna=pd%lzoh-pd%lzom
    call dyerhicks(integralh,z_on_l,cd,thetavstar,thetav,sthetav,a_umag,pd%cndzmin,pd%lzom,lna,4)
    ustar=sqrt(cd)*a_umag
    qstar=vkar*(a_mixr-smixr)/integralh  
    tstar=vkar*(a_temp-u_ts)/integralh
    
    z0_on_l  = z0*z_on_l/pd%cndzmin
    z10_on_l = z10*z_on_l/pd%cndzmin
    z0_on_l  = min(z0_on_l,10.)
    z10_on_l = min(z10_on_l,10.)
    neutral   = log(pd%cndzmin/z0)
    neutral10 = log(pd%cndzmin/z10)
    where (z_on_l<0.)
      ph0     = (1.-16.*z0_on_l)**(-0.50)
      ph1     = (1.-16.*z_on_l)**(-0.50)
      pm0     = (1.-16.*z0_on_l)**(-0.25)
      pm10    = (1.-16.*z10_on_l)**(-0.25)
      pm1     = (1.-16.*z_on_l)**(-0.25)
      integralh   = neutral-2.*log((1.+1./ph1)/(1.+1./ph0))
      integralm   = neutral-2.*log((1.+1./pm1)/(1.+1./pm0)) &
                    -log((1.+1./pm1**2)/(1.+1./pm0**2))     &
                    +2.*(atan(1./pm1)-atan(1./pm0))
      integralm10 = neutral10-2.*log((1.+1./pm1)/(1.+1./pm10)) &
                    -log((1.+1./pm1**2)/(1.+1./pm10**2))       &
                    +2.*(atan(1./pm1)-atan(1./pm10))     
    elsewhere
      !-------Beljaars and Holtslag (1991) heat function
      ph0  = -((1.+(2./3.)*a_1*z0_on_l)**1.5 &
             +b_1*(z0_on_l-(c_1/d_1))        &
             *exp(-d_1*z0_on_l)+b_1*c_1/d_1-1.)
      ph1  = -((1.+(2./3.)*a_1*z_on_l)**1.5 &
             +b_1*(z_on_l-(c_1/d_1))        &
             *exp(-d_1*z_on_l)+b_1*c_1/d_1-1.)
      pm0 = -(a_1*z0_on_l+b_1*(z0_on_l-(c_1/d_1))*exp(-d_1*z0_on_l)+b_1*c_1/d_1)
      pm10 = -(a_1*z10_on_l+b_1*(z10_on_l-(c_1/d_1))*exp(-d_1*z10_on_l)+b_1*c_1/d_1)
      pm1  = -(a_1*z_on_l+b_1*(z_on_l-(c_1/d_1))*exp(-d_1*z_on_l)+b_1*c_1/d_1)
      integralh   = neutral-(ph1-ph0)
      integralm   = neutral-(pm1-pm0)
      integralm10 = neutral10-(pm1-pm10)
    endwhere
    pd%tscrn = a_temp - tstar*integralh/vkar
    pd%qscrn = a_mixr - qstar*integralh/vkar
    pd%uscrn = max(a_umag-ustar*integralm/vkar,0.)
    pd%u10   = max(a_umag-ustar*integralm10/vkar,0.)
    
  case(1) ! estimate screen diagnostics (two step canopy approach)
    thetav=d_tempc + (d_tempc+urbtemp)*0.61*a_mixr
    sthetav=u_ts + (u_ts+urbtemp)*0.61*smixr
    lna=pd%lzoh-pd%lzom
    call dyerhicks(integralh,z_on_l,cd,thetavstar,thetav,sthetav,a_umag,pd%cndzmin,pd%lzom,lna,4)
    ustar=sqrt(cd)*a_umag
    qstar=vkar*(a_mixr-smixr)/integralh
    tts=vkar*(thetav-sthetav)/integralh
    tstar=vkar*(a_temp-u_ts)/integralh
    
    z0_on_l  = fp%bldheight*(1.-refheight)*z_on_l/pd%cndzmin ! calculate at canyon top
    z10_on_l = max(z10-fp%bldheight*refheight,1.)*z_on_l/pd%cndzmin
    z0_on_l  = min(z0_on_l,10.)
    z10_on_l = min(z10_on_l,10.)
    neutral   = log(pd%cndzmin/(fp%bldheight*(1.-refheight)))
    neutral10 = log(pd%cndzmin/max(z10-fp%bldheight*refheight,1.))
    where (z_on_l<0.)
      ph0     = (1.-16.*z0_on_l)**(-0.50)
      ph1     = (1.-16.*z_on_l)**(-0.50)
      pm0     = (1.-16.*z0_on_l)**(-0.25)
      pm10    = (1.-16.*z10_on_l)**(-0.25)
      pm1     = (1.-16.*z_on_l)**(-0.25)
      integralh = neutral-2.*log((1.+1./ph1)/(1.+1./ph0))
      integralm = neutral-2.*log((1.+1./pm1)/(1.+1./pm0)) &
                    -log((1.+1./pm1**2)/(1.+1./pm0**2))   &
                    +2.*(atan(1./pm1)-atan(1./pm0))
      integralm10 = neutral10-2.*log((1.+1./pm1)/(1.+1./pm10)) &
                    -log((1.+1./pm1**2)/(1.+1./pm10**2))       &
                    +2.*(atan(1./pm1)-atan(1./pm10))     
    elsewhere
      !-------Beljaars and Holtslag (1991) heat function
      ph0  = -((1.+(2./3.)*a_1*z0_on_l)**1.5 &
             +b_1*(z0_on_l-(c_1/d_1))        &
             *exp(-d_1*z0_on_l)+b_1*c_1/d_1-1.)
      ph1  = -((1.+(2./3.)*a_1*z_on_l)**1.5 &
             +b_1*(z_on_l-(c_1/d_1))        &
             *exp(-d_1*z_on_l)+b_1*c_1/d_1-1.)
      pm0 = -(a_1*z0_on_l+b_1*(z0_on_l-(c_1/d_1))*exp(-d_1*z0_on_l)+b_1*c_1/d_1)
      pm10 = -(a_1*z10_on_l+b_1*(z10_on_l-(c_1/d_1))*exp(-d_1*z10_on_l)+b_1*c_1/d_1)
      pm1  = -(a_1*z_on_l+b_1*(z_on_l-(c_1/d_1))*exp(-d_1*z_on_l)+b_1*c_1/d_1)
      integralh   = neutral-(ph1-ph0)
      integralm   = neutral-(pm1-pm0)
      integralm10 = neutral10-(pm1-pm10)
    endwhere
    ttop = thetav - tts*integralh/vkar
    tetp = a_temp - tstar*integralh/vkar
    qtop = a_mixr - qstar*integralh/vkar
    utop = a_umag - ustar*integralm/vkar

    where (fp%bldheight<=z10) ! above canyon
      pd%u10=max(a_umag-ustar*integralm10/vkar,0.)
    end where

    ! assume standard stability functions hold for urban canyon (needs more work)
    tsurf = d_rdsndelta*rdsntemp+(1.-d_rdsndelta)*((1.-cnveg%sigma)*real(road%nodetemp(:,0))+cnveg%sigma*cnveg%temp)
    n=max(min((rdhyd%soilwater-fp%swilt)/(fp%sfc-fp%swilt),1.),0.)
    wf = (1.-d_rdsndelta)*((1.-cnveg%sigma)*d_roaddelta+cnveg%sigma*((1.-d_vegdeltac)*n+d_vegdeltac))
    call getqsat(qsurf,tsurf,d_sigd)
    qsurf=qsurf*wf
    n=log(fp%bldheight/zonet)
    
    thetav=ttop + (ttop+urbtemp)*0.61*qtop
    sthetav=tsurf + (tsurf+urbtemp)*0.61*qsurf
    lna=2.3
    call dyerhicks(integralh,z_on_l,cd,thetavstar,thetav,sthetav,utop,fp%bldheight,n,lna,1)
    ustar=sqrt(cd)*utop
    tstar=vkar*(tetp-tsurf)/integralh
    qstar=vkar*(qtop-qsurf)/integralh
    
    z0_on_l   = z0*z_on_l/fp%bldheight
    z10_on_l  = max(z10,fp%bldheight)*z_on_l/fp%bldheight
    z0_on_l   = min(z0_on_l,10.)
    z10_on_l  = min(z10_on_l,10.)
    neutral   = log(fp%bldheight/z0)
    neutral10 = log(fp%bldheight/max(z10,fp%bldheight))
    where (z_on_l<0.)
      ph0     = (1.-16.*z0_on_l)**(-0.50)
      ph1     = (1.-16.*z_on_l)**(-0.50)
      pm0     = (1.-16.*z0_on_l)**(-0.25)
      pm10    = (1.-16.*z10_on_l)**(-0.25)
      pm1     = (1.-16.*z_on_l)**(-0.25)
      integralh = neutral-2.*log((1.+1./ph1)/(1.+1./ph0))
      integralm = neutral-2.*log((1.+1./pm1)/(1.+1./pm0)) &
                    -log((1.+1./pm1**2)/(1.+1./pm0**2))   &
                    +2.*(atan(1./pm1)-atan(1./pm0))
      integralm10 = neutral10-2.*log((1.+1./pm1)/(1.+1./pm10)) &
                    -log((1.+1./pm1**2)/(1.+1./pm10**2))       &
                    +2.*(atan(1./pm1)-atan(1./pm10))     
    elsewhere
      !-------Beljaars and Holtslag (1991) heat function
      ph0  = -((1.+(2./3.)*a_1*z0_on_l)**1.5 &
             +b_1*(z0_on_l-(c_1/d_1))        &
             *exp(-d_1*z0_on_l)+b_1*c_1/d_1-1.)
      ph1  = -((1.+(2./3.)*a_1*z_on_l)**1.5 &
             +b_1*(z_on_l-(c_1/d_1))        &
             *exp(-d_1*z_on_l)+b_1*c_1/d_1-1.)
      pm0 = -(a_1*z0_on_l+b_1*(z0_on_l-(c_1/d_1))*exp(-d_1*z0_on_l)+b_1*c_1/d_1)
      pm10 = -(a_1*z10_on_l+b_1*(z10_on_l-(c_1/d_1))*exp(-d_1*z10_on_l)+b_1*c_1/d_1)
      pm1  = -(a_1*z_on_l+b_1*(z_on_l-(c_1/d_1))*exp(-d_1*z_on_l)+b_1*c_1/d_1)
      integralh   = neutral-(ph1-ph0)
      integralm   = neutral-(pm1-pm0)
      integralm10 = neutral10-(pm1-pm10)
    endwhere

    pd%tscrn = tetp-tstar*integralh/vkar
    pd%qscrn = qtop-qstar*integralh/vkar
    pd%uscrn = max(utop-ustar*integralm/vkar,0.)
    where (fp%bldheight>z10) ! within canyon
      pd%u10 = max(utop-ustar*integralm10/vkar,0.)
    end where

  case(2) ! calculate screen diagnostics from canyon only
    tsurf=d_rdsndelta*rdsntemp+(1.-d_rdsndelta)*((1.-cnveg%sigma)*real(road%nodetemp(:,0))+cnveg%sigma*cnveg%temp)
    n=max(min((rdhyd%soilwater-fp%swilt)/(fp%sfc-fp%swilt),1.),0.)
    wf=(1.-d_rdsndelta)*((1.-cnveg%sigma)*d_roaddelta+cnveg%sigma*((1.-d_vegdeltac)*n+d_vegdeltac))
    call getqsat(qsurf,tsurf,d_sigd)
    qsurf=qsurf*wf
    n=log(fp%bldheight/zonet)

    thetav=d_tempc + (d_tempc+urbtemp)*0.61*a_mixr
    sthetav=tsurf + (tsurf+urbtemp)*0.61*qsurf
    lna=2.3
    call dyerhicks(integralh,z_on_l,cd,thetavstar,thetav,sthetav,a_umag,pd%cndzmin,n,lna,1)
    ustar=sqrt(cd)*a_umag
    qstar=vkar*(a_mixr-smixr)/integralh
    tstar=vkar*(a_temp-tsurf)/integralh
    
    z0_on_l  = z0*z_on_l/pd%cndzmin
    z10_on_l = z10*z_on_l/pd%cndzmin
    z0_on_l  = min(z0_on_l,10.)
    z10_on_l = min(z10_on_l,10.)
    neutral   = log(pd%cndzmin/z0)
    neutral10 = log(pd%cndzmin/z10)
    where (z_on_l<0.)
      ph0     = (1.-16.*z0_on_l)**(-0.50)
      ph1     = (1.-16.*z_on_l)**(-0.50)
      pm0     = (1.-16.*z0_on_l)**(-0.25)
      pm10    = (1.-16.*z10_on_l)**(-0.25)
      pm1     = (1.-16.*z_on_l)**(-0.25)
      integralh   = neutral-2.*log((1.+1./ph1)/(1.+1./ph0))
      integralm   = neutral-2.*log((1.+1./pm1)/(1.+1./pm0)) &
                    -log((1.+1./pm1**2)/(1.+1./pm0**2))     &
                    +2.*(atan(1./pm1)-atan(1./pm0))
      integralm10 = neutral10-2.*log((1.+1./pm1)/(1.+1./pm10)) &
                    -log((1.+1./pm1**2)/(1.+1./pm10**2))       &
                    +2.*(atan(1./pm1)-atan(1./pm10))     
    elsewhere
      !-------Beljaars and Holtslag (1991) heat function
      ph0  = -((1.+(2./3.)*a_1*z0_on_l)**1.5 &
             +b_1*(z0_on_l-(c_1/d_1))        &
             *exp(-d_1*z0_on_l)+b_1*c_1/d_1-1.)
      ph1  = -((1.+(2./3.)*a_1*z_on_l)**1.5 &
             +b_1*(z_on_l-(c_1/d_1))        &
             *exp(-d_1*z_on_l)+b_1*c_1/d_1-1.)
      pm0 = -(a_1*z0_on_l+b_1*(z0_on_l-(c_1/d_1))*exp(-d_1*z0_on_l)+b_1*c_1/d_1)
      pm10 = -(a_1*z10_on_l+b_1*(z10_on_l-(c_1/d_1))*exp(-d_1*z10_on_l)+b_1*c_1/d_1)
      pm1  = -(a_1*z_on_l+b_1*(z_on_l-(c_1/d_1))*exp(-d_1*z_on_l)+b_1*c_1/d_1)
      integralh   = neutral-(ph1-ph0)
      integralm   = neutral-(pm1-pm0)
      integralm10 = neutral10-(pm1-pm10)
    endwhere
    pd%tscrn = a_temp-tstar*integralh/vkar
    pd%qscrn = a_mixr-qstar*integralh/vkar
    pd%uscrn = max(a_umag-ustar*integralm/vkar,0.)
    pd%u10   = max(a_umag-ustar*integralm10/vkar,0.)
    
end select
pd%qscrn       = max(pd%qscrn,1.E-4)
      
return
end subroutine scrncalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The following subroutines are used for internal varying temperature
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculates net longwave radiation flux (flux_rg) at each surface
! longwave flux is temperature dependent, so this subroutine should be run at each timestep
subroutine internal_lwflux(rgint_slab,rgint_wallw,rgint_roof,rgint_walle, &
                           skintemps,fp,intl,ufull)

use uclem_parameters
use uclem_types

implicit none

integer, intent(in) :: ufull
integer :: j
real, dimension(ufull,4), intent(in) :: skintemps  ! floor, wall, ceiling, wall temperature array
real, dimension(ufull,4) :: epsil     ! floor, wall, ceiling, wall emissivity array
real, dimension(ufull,4) :: radnet    ! net flux density on ith surface (+ve leaving)
real, dimension(ufull)   :: radtot    ! net leaving flux density (B) on ith surface
real, dimension(ufull), intent(out) :: rgint_slab,rgint_wallw,rgint_roof,rgint_walle
type(intldata), intent(in) :: intl
type(fparmdata), intent(in) :: fp

radnet = 0.
epsil  = 0.9

!epsil = reshape((/(fp_slab%emiss,fp_wall%emiss,fp_roof%emiss,fp_wall%emiss, & 
!                    i=1,ufull)/), (/ufull,4/))

do j = 2,4
  radnet(:,j) = epsil(:,j)/(1.-epsil(:,j))*((sbconst*(skintemps(:,j)+urbtemp)**4)  & 
               - sum(intl%psi(:,j,:)*(sbconst*(skintemps(:,:)+urbtemp)**4),dim=2))
end do

radnet(:,1) = -fp%bldheight*(radnet(:,2)+radnet(:,4))/fp%bldwidth - radnet(:,3)

! energy conservation check
radtot(:) = abs(fp%bldwidth(:)*(radnet(:,1)+radnet(:,3)) + fp%bldheight*(radnet(:,2)+radnet(:,4)))

do j = 1,ufull
  if ( radtot(j)>energytol ) write(6,*) "error: radiation energy non-closure: ", radtot(j)
end do

rgint_slab  = radnet(:,1)
rgint_wallw = radnet(:,2)
rgint_roof  = radnet(:,3)
rgint_walle = radnet(:,4)

return
end subroutine internal_lwflux

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine sets internal surface convective heat transfer coefficients
! Compares temperature of innermost layer temperature with air temperature
! Considers horizontal and vertical orientation
! Based on EnergyPlus: Simple Natural Convection Algorithm [W m^-2 K^-1]

subroutine calc_convcoeff(cvcoeff_roof,cvcoeff_walle,cvcoeff_wallw,cvcoeff_slab, & 
                          cvcoeff_intm1,cvcoeff_intm2,roomtemp,roof,slab,intm,ufull)

use uclem_parameters
use uclem_types

implicit none

integer, intent(in) :: ufull
real, dimension(ufull), intent(out) :: cvcoeff_roof,cvcoeff_walle,cvcoeff_wallw
real, dimension(ufull), intent(out) :: cvcoeff_intm1,cvcoeff_intm2,cvcoeff_slab
real, dimension(ufull), intent(in)  :: roomtemp
type(facetdata), intent(in) :: roof
type(facetdata), intent(in) :: slab
type(facetdata), intent(in) :: intm
real, parameter :: tdel = 0.1
real, dimension(ufull) :: x

select case(cvcoeffmeth)
  case(0) ! DOE Simple
    cvcoeff_walle = 3.067    ! vertical surface coefficient constant
    cvcoeff_wallw = 3.067    ! vertical surface coefficient constant
    where ( real(roof%nodetemp(:,nl))>=roomtemp+tdel )
      cvcoeff_roof(:)=0.948  ! reduced convection 
    elsewhere ( real(roof%nodetemp(:,nl))>roomtemp-tdel )
      x = (real(roof%nodetemp(:,nl))-roomtemp+tdel)/(2.*tdel)
      cvcoeff_roof(:)=(1.-x)*4.04 + x*0.948 ! blend
    elsewhere
      cvcoeff_roof(:)=4.040  ! enhanced convection  
    end where    
    cvcoeff_intm1 = 3.067   ! vertical surface coefficient constant
    cvcoeff_intm2 = 3.067   ! vertical surface coefficient constant
    where ( real(slab%nodetemp(:,nl))<=roomtemp-tdel )
      cvcoeff_slab(:)=0.7   ! reduced convection
    elsewhere ( real(slab%nodetemp(:,nl))<roomtemp+tdel )
      x = (real(slab%nodetemp(:,nl))-roomtemp+tdel)/(2.*tdel)
      cvcoeff_slab(:) = (1.-x)*0.7 + x*4.04 ! blend
    elsewhere
      cvcoeff_slab(:)=4.040  ! enhanced convection
    end where
  case(1) ! dynamic, from international standard ISO6946:2007, Annex A
    cvcoeff_walle = 2.5    ! vertical surface coefficient constant
    cvcoeff_wallw = 2.5    ! vertical surface coefficient constant
    where ( real(roof%nodetemp(:,nl))>=roomtemp+tdel )
      cvcoeff_roof(:)=0.7  ! reduced convection (upper surface)
    elsewhere ( real(roof%nodetemp(:,nl))>roomtemp-tdel )
      x = (real(roof%nodetemp(:,nl))-roomtemp+tdel)/(2.*tdel)
      cvcoeff_roof(:) = (1.-x)*5.0 + x*0.7 ! blend
    elsewhere
      cvcoeff_roof(:)=5.0  ! enhanced convection (upper surface)  
    end where
    where ( real(intm%nodetemp(:,nl))>=roomtemp+tdel )
      cvcoeff_intm2 = 0.7+5.7   ! reduced convection (upper surface) + radiation @20 deg C
    elsewhere ( real(intm%nodetemp(:,nl))>roomtemp-tdel )
      x = (real(intm%nodetemp(:,nl))-roomtemp+tdel)/(2.*tdel)
      cvcoeff_intm2 = (1.-x)*(5.+5.7) + x*(0.7+5.7) ! blend
    elsewhere
      cvcoeff_intm2 = 5.0+5.7   ! reduced convection (upper surface) + radiation @20 deg C
    end where
    where ( real(intm%nodetemp(:,0))<=roomtemp-tdel )    
      cvcoeff_intm1 = 0.7+5.7   ! reduced convection (lower surface) + radiation @20 deg C
    elsewhere ( real(intm%nodetemp(:,0))<roomtemp+tdel )
      x = (real(intm%nodetemp(:,nl))-roomtemp+tdel)/(2.*tdel)
      cvcoeff_intm1 = (1.-x)*(0.7+5.7) + x*(5.+5.7) ! blend
    elsewhere
      cvcoeff_intm1 = 5.0+5.7   ! reduced convection (lower surface) + radiation @20 deg C
    end where
    where ( real(slab%nodetemp(:,nl))<=roomtemp-tdel )
      cvcoeff_slab(:)=0.7   ! reduced convection (lower surface)
    elsewhere ( real(slab%nodetemp(:,nl))<roomtemp+tdel )
      x = (real(slab%nodetemp(:,nl))-roomtemp+tdel)/(2.*tdel)
      cvcoeff_slab = (1.-x)*0.7 + x*5.
    elsewhere
      cvcoeff_slab(:)=5.0   ! enhanced convection (lower surface)
    end where
  case(2) ! fixed, from international standard IS6946:2007, 5.2 
    cvcoeff_roof  = 1./0.10 ! fixed coefficient up
    cvcoeff_walle = 1./0.13 ! fixed coefficient horizontal
    cvcoeff_wallw = 1./0.13 ! fixed coefficient horizontal
    cvcoeff_intm1 = 1./0.10 ! fixed coefficient up
    cvcoeff_intm2 = 1./0.17 ! fixed coefficient down
    cvcoeff_slab  = 1./0.17 ! fixed coefficient down
end select

end subroutine calc_convcoeff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculates the proportion of open windows for ventilation
subroutine calc_openwindows(d_openwindows,fp,iroomtemp,d_canyontemp, &
                            roof,walle,wallw,slab,ufull)

use uclem_parameters
use uclem_types

implicit none

integer, intent(in)                 :: ufull
real, dimension(ufull), intent(in)  :: d_canyontemp,iroomtemp
real, dimension(ufull), intent(out) :: d_openwindows
real, dimension(ufull)              :: xtemp, mrt
type(facetdata), intent(in) :: roof, walle, wallw, slab
type(fparmdata), intent(in) :: fp

! mean radiant temperature estimation
mrt = 0.5*(fp%bldheight/(fp%bldwidth+fp%bldheight)*(real(walle%nodetemp(:,nl) + wallw%nodetemp(:,nl)))) & 
    + 0.5*(fp%bldwidth/ (fp%bldwidth+fp%bldheight)*(real(roof%nodetemp(:,nl) + slab%nodetemp(:,nl))))
! globe temperature approximation (average of mrt and air temperature) [Celsius]
xtemp = 0.5*(iroomtemp + mrt)

d_openwindows = 1./(1. + exp( 0.5*(fp%bldairtemp+ac_deltat+5.-xtemp) )) &
    *1./(1. + exp( max(min(d_canyontemp-iroomtemp,40.),-40.) ))

return
end subroutine calc_openwindows

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine getnetalbedo(alb,sg_roof,sg_vegr,sg_road,sg_walle,sg_wallw,sg_vegc,sg_rfsn,sg_rdsn,  &
                        fp_hwratio,fp_sigmabld,fp_vegsigmar,fp_roofalpha,fp_vegalphar,          &
                        fp_vegsigmac,fp_roadalpha,fp_wallalpha,fp_vegalphac,                    &
                        roofalpha,roadalpha,snowdeltar,snowdeltac)

implicit none

integer iqu
real, dimension(:), intent(out) :: alb
real, dimension(:), intent(in) :: sg_roof, sg_vegr, sg_road, sg_walle, sg_wallw, sg_vegc
real, dimension(:), intent(in) :: sg_rfsn, sg_rdsn
real, dimension(:), intent(in) :: fp_hwratio, fp_sigmabld
real, dimension(:), intent(in) :: fp_vegsigmar, fp_roofalpha, fp_vegalphar
real, dimension(:), intent(in) :: fp_vegsigmac, fp_roadalpha, fp_vegalphac, fp_wallalpha
real, dimension(:), intent(in) :: roofalpha, roadalpha, snowdeltar, snowdeltac
real albu, albr

do iqu = 1,size(alb)
  ! canyon
  albu=1.-(fp_hwratio(iqu)*(sg_walle(iqu)+sg_wallw(iqu))*(1.-fp_wallalpha(iqu))         &
      +snowdeltac(iqu)*sg_rdsn(iqu)*(1.-roadalpha(iqu))                                 &
      +(1.-snowdeltac(iqu))*((1.-fp_vegsigmac(iqu))*sg_road(iqu)*(1.-fp_roadalpha(iqu)) &
      +fp_vegsigmac(iqu)*sg_vegc(iqu)*(1.-fp_vegalphac(iqu))))
  ! roof
  albr=(1.-snowdeltar(iqu))*((1.-fp_vegsigmar(iqu))*sg_roof(iqu)*fp_roofalpha(iqu) &
      +fp_vegsigmar(iqu)*sg_vegr(iqu)*fp_vegalphar(iqu))                           &
      +snowdeltar(iqu)*sg_rfsn(iqu)*roofalpha(iqu)
  ! net
  alb(iqu)=fp_sigmabld(iqu)*albr+(1.-fp_sigmabld(iqu))*albu
end do

return
end subroutine getnetalbedo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine init_lwcoeff(fp,intl,ufull)
! This subroutine calculates longwave reflection coefficients (int_psi) at each surface
! longwave coefficients do not change, so this subroutine should only be run once
! Infinite reflections per Harman et al., (2004) "Radiative Exchange in Urban Street Canyons"
! Per method in "Radiation Heat Transfer, Sparrow & Cess 1978, Ch 3-3"
! array surface order is: (1) floor; (2) wallw; (3) ceiling; (4) walle

use uclem_types

implicit none

! local variables
integer, intent(in) :: ufull
real, dimension(ufull,4,4) :: chi
real, dimension(4,4)       :: krondelta
real, dimension(ufull)     :: h, w
real, dimension(ufull,4)   :: epsil   ! floor, wall, ceiling, wall emissivity array
integer :: i, j
integer :: ierr       ! inverse matrix error flag
type(intldata), intent(inout) :: intl
type(fparmdata), intent(in) :: fp


krondelta = 0.
chi = 0.
intl%psi = 0.
h = fp%bldheight
w = fp%sigmabld*(fp%bldheight/fp%hwratio)/(1.-fp%sigmabld)

! set int_vfactors
intl%viewf(:,1,1) = 0.                                    ! floor to self
intl%viewf(:,1,2) = 0.5*(1.+(h/w)-sqrt(1.+(h/w)**2))      ! floor to wallw
intl%viewf(:,1,3) = sqrt(1.+(h/w)**2)-(h/w)               ! floor to ceiling
intl%viewf(:,1,4) = intl%viewf(:,1,2)                     ! floor to walle
intl%viewf(:,2,1) = 0.5*(1.+(w/h)-sqrt(1.+(w/h)**2))      ! wallw to floor
intl%viewf(:,2,2) = 0.                                    ! wallw to self
intl%viewf(:,2,3) = intl%viewf(:,2,1)                     ! wallw to ceiling
intl%viewf(:,2,4) = sqrt(1.+(w/h)**2)-(w/h)               ! wallw to walle
intl%viewf(:,3,1) = intl%viewf(:,1,3)                     ! ceiling to floor
intl%viewf(:,3,2) = intl%viewf(:,1,2)                     ! ceiling to wallw
intl%viewf(:,3,3) = 0.                                    ! ceiling to self
intl%viewf(:,3,4) = intl%viewf(:,1,2)                     ! ceiling walle
intl%viewf(:,4,1) = intl%viewf(:,2,1)                     ! walle to floor
intl%viewf(:,4,2) = intl%viewf(:,2,4)                     ! walle to wallw
intl%viewf(:,4,3) = intl%viewf(:,2,1)                     ! walle to ceiling
intl%viewf(:,4,4) = 0.                                    ! walle to self

!epsil = reshape((/(f_slab%emiss,f_wall%emiss,f_roof%emiss,f_wall%emiss, & 
!                    i=1,ufull_g)/), (/ufull_g,4/))
epsil = 0.9

do i = 1,4
  krondelta(i,i) = 1.
end do
do j = 1,4
  do i = 1,4
    chi(:,i,j) = (krondelta(i,j) - (1.-epsil(:,i))*intl%viewf(:,i,j))/(epsil(:,i))
  end do
end do

! invert matrix
intl%psi = chi
call minverse(intl%psi,ierr)

end subroutine init_lwcoeff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine initialises internal variables; 
! building width and number of internal mass floors
subroutine init_internal(fp)

use uclem_parameters
use uclem_types

implicit none

type(fparmdata), intent(inout) :: fp
real :: floor_height = 3. ! assumed average floor to floor height (m)

fp%bldwidth = fp%sigmabld*(fp%bldheight/fp%hwratio)/max(1.-fp%sigmabld,0.001)
! define number of internal mass floors (based on building height)
select case(intmassmeth)
  case(0) ! no internal mass
    fp%intmassn = 0
  case(1) ! one floor of internal mass
    fp%intmassn = 1
  case(2) ! dynamic floors of internal mass
    fp%intmassn = max((nint(fp%bldheight/floor_height)-1),0)
end select

end subroutine init_internal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculates the inverse of a NxN matrix
! input/output = a, size=s, error flag=ier
 
subroutine minverse(a_inout,ierr)
 
implicit none
 
real, dimension(:,:,:), intent(inout) :: a_inout
real(kind=8), dimension(size(a_inout,1),size(a_inout,2),size(a_inout,3)) :: a
real(kind=8), dimension(size(a_inout,1)) :: det, d, amax
real(kind=8), dimension(size(a_inout,2)) :: x
real(kind=8) :: y
integer, intent(out)  :: ierr
integer s, ns, iq, i, j, nu
integer, dimension(size(a_inout,1),size(a_inout,2)) :: row, col
integer, dimension(size(a_inout,1)) :: prow, pcol
logical, dimension(size(a_inout,1),size(a_inout,2)) :: notpiv

a = real(a_inout,8)

nu = size(a_inout,1)
s = size(a_inout,2)

det = 0.
d = 1.
notpiv = .TRUE.
 
do ns = 1,s
 
  amax(:) = 0.
  do j = 1,s
    do i = 1,s
      where ( notpiv(:,j) .and. notpiv(:,i) .and. amax(:)<abs(a(:,i,j)) )
        amax(:) = abs(a(:,i,j))
        prow(:) = i
        pcol(:) = j
      end where  
    end do
  end do
 
  if ( any(amax<0.) ) then
    ierr = 1 ! cannot invert matrix
    return
  end if

  do iq = 1,nu
    notpiv(iq,pcol(iq)) = .FALSE.
    if ( prow(iq)/=pcol(iq) ) then
      d(iq) = -d(iq)
      x(1:s) = a(iq,prow(iq),1:s)
      a(iq,prow(iq),1:s) = a(iq,pcol(iq),1:s)
      a(iq,pcol(iq),1:s) = x(1:s)
    end if
  end do  
 
  row(:,ns) = prow(:)
  col(:,ns) = pcol(:)
  do iq = 1,nu
    amax(iq) = a(iq,pcol(iq),pcol(iq))
  end do  
  d(:) = d(:)*amax(:)
 
  if ( any(abs(d)<=0.) ) then
    ierr = 1 ! cannot invert matrix
    return
  end if
 
  amax(:) = 1./amax(:)
  do iq = 1,nu
    a(iq,pcol(iq),pcol(iq)) = 1.
    a(iq,pcol(iq),1:s) = a(iq,pcol(iq),1:s)*amax(iq)
  end do  
 
  do i=1,s
    do iq = 1,nu
      if ( i/=pcol(iq) ) then
        y = a(iq,i,pcol(iq))
        a(iq,i,pcol(iq)) = 0.
        do j = 1,s
          a(iq,i,j) = a(iq,i,j) - y*a(iq,pcol(iq),j)
        end do
      end if
    end do  
  end do
 
end do
 
det(:) = d(:)
 
do ns = s,1,-1
  prow(:) = row(:,ns)
  pcol(:) = col(:,ns)
  do iq = 1,nu
    if ( prow(iq)/=pcol(iq) ) then
      do i = 1,s
        y = a(iq,i,prow(iq))
        a(iq,i,prow(iq)) = a(iq,i,pcol(iq))
        a(iq,i,pcol(iq)) = y
      end do
    end if
  end do  
end do

a_inout = real(a)

ierr = 0
 
return
end subroutine minverse

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate albedo
!

subroutine uclem_calc_alb(alb,fp_hwratio,fp_effhwratio,fp_vangle,      &
                          fp_hangle,fp_fbeam,                          &
                          fp_vegsigmac,fp_vegsigmar,fp_sigmabld,       &
                          fp_roadalpha,fp_wallalpha,fp_roofalpha,      &
                          fp_vegalphac,fp_vegalphar,                   &
                          road_snow_alpha,roof_snow_alpha,             &
                          can_snow,roof_snow)

use uclem_parameters

implicit none

real, dimension(:), intent(in) :: fp_hwratio, fp_effhwratio, fp_vangle, fp_hangle, fp_fbeam
real, dimension(:), intent(in) :: fp_roadalpha, fp_wallalpha, fp_roofalpha
real, dimension(:), intent(in) :: fp_vegalphar, fp_vegalphac
real, dimension(:), intent(in) :: road_snow_alpha, roof_snow_alpha
real, dimension(:), intent(in) :: can_snow, roof_snow
real, dimension(:), intent(in) :: fp_sigmabld, fp_vegsigmar, fp_vegsigmac
real, dimension(:), intent(out) :: alb
real, dimension(size(alb)) :: snowdeltar, snowdeltac, wallpsi, roadpsi
real, dimension(size(alb)) :: sg_roof, sg_vegr, sg_road, sg_walle, sg_wallw, sg_vegc, sg_rfsn, sg_rdsn
real, dimension(size(alb)) :: fp_coeffbldheight

fp_coeffbldheight = fp_effhwratio/fp_hwratio

! roof
snowdeltar = roof_snow/(roof_snow+maxrfsn)
! canyon
snowdeltac = can_snow/(can_snow+maxrdsn)
          
! call the subroutine
call getswcoeff(sg_roof,sg_vegr,sg_road,sg_walle,sg_wallw,       &
                sg_vegc,sg_rfsn,sg_rdsn,wallpsi,roadpsi,         &
                fp_hwratio,fp_effhwratio,fp_vangle,fp_hangle,    &
                fp_fbeam,fp_vegsigmac,fp_roadalpha,fp_vegalphac, &
                fp_wallalpha,road_snow_alpha,snowdeltac)
sg_walle = sg_walle*fp_coeffbldheight
sg_wallw = sg_wallw*fp_coeffbldheight

call getnetalbedo(alb,sg_roof,sg_vegr,sg_road,sg_walle,          &
                  sg_wallw,sg_vegc,sg_rfsn,sg_rdsn,              &
                  fp_hwratio,fp_sigmabld,fp_vegsigmar,           &
                  fp_roofalpha,fp_vegalphar,fp_vegsigmac,        &
                  fp_roadalpha,fp_wallalpha,fp_vegalphac,        &
                  roof_snow_alpha,road_snow_alpha,               &
                  snowdeltar,snowdeltac)

return
end subroutine uclem_calc_alb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module uclem
