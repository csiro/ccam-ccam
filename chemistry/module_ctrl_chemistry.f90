! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2026 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
module module_ctrl_chemistry

use aerointerface                          ! Aerosol interface
use timeseries                             ! Tracer time series
use tracermodule                           ! Tracer routines
use tracers_m                              ! Tracer data
use trvmix                                 ! Tracer mixing routines

implicit none

private

! aerointerface
public load_aerosolldr, aerocalc, aerocalc_init, aerodrop, convscav
public ppfprec, ppfmelt, ppfsnow, ppfevap, ppfsubl, pplambs, ppmrate
public ppmaccr, ppqfsedice, pprscav, pprfreeze
public opticaldepth, updateoxidant, oxidant_timer
public aerosol_u10, naerofamilies, aero_split, aeroindir
public ch_dust, zvolcemi, so4mtn, carbmtn, saltsmallmtn, saltlargemtn, enhanceu10
public scav_effl
public dustreff
public naero
public aldrinit,aldrend,aldrloademiss,aldrloaderod
public xtg,xtgsav,xtosav
public itracdu,ndust
public dustdd,dustwd,duste,dust_burden
public itracbc,bce,bcdd,bcwd,bc_burden
public itracoc,oce,ocdd,ocwd,oc_burden
public itracdms,itracso2,itracso4
public dmse,dmsso2o,so2e,so2so4o,so2dd,so2wd,so4e,so4dd,so4wd
public dms_burden,so2_burden,so4_burden
public itracsa,nsalt,salte,saltdd,saltwd,salt_burden
public xtg_solub,zoxidant_g,erod,ndcls,emissfield,vso2

! timeseries
public init_ts, write_ts

! tracermodule
public sitefile,shipfile
public co2em,init_tracer
public tracer_mass,interp_tracerflux,tracerlist
public writetrpm
public readtracerflux,tracunit,tracdaytime
public oh,strloss,mcfdep,jmcf
public trden, trreff, trdep, traclevel
public numtracer

! tracers_m
public ngas,ntrac,nllp
public tr,traver,acloss_g
public trpm,npm
public tractype,tracname
public tracers_init,tracers_end

! trvmix
public tracervmix
    
end module module_ctrl_chemistry