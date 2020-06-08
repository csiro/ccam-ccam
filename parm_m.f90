! Conformal Cubic Atmospheric Model
    
! Copyright 2016-2019 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
module parm_m

implicit none

private
public ngwd,nrungcm,newtop
public kountr,nrad,nvmix,nlocal
public nhstest,namip,nspecial,newrough,nsib
public ntaft,ntsea,ntsur,ntsur2,lgwd,newztsea,nglacier,mbd
public mbd_mlo,nbd,kbotdav,kbotu,nud_p,nud_q,nud_t,nud_uv
public nud_hrs,nudu_hrs,ktau,ndi,ndi2,ntau,nperavg,nperday,nperhr
public nmaxpr,nlv,ia,ib,ja,jb,id,jd,idjd
public io_in,io_out,io_rest
public nwt,nrun,nextout,m_fly,nsemble,tbave
public nurban,nmr,nmlo,ktopdav,nud_sst,nud_sss,kbotmlo,ktopmlo
public mloalpha,nud_ouv,nud_sfh,kblock,rescrn,knh,iaero
public nud_aero,mbd_maxscale,mbd_maxgrid,mbd_maxscale_mlo,nriver
public leap,nbarewet,nsigmf,qg_fix
public qgmin, zo_clearing
public av_vmod, vmodmin, snmin, tss_sh, charnock, chn10, zobgin
public rlongdn, rlongdx, rlatdn, rlatdx, ds, dt, dtin, panfg, panzo
public bpyear, helim, fc2, sigbot_gwd, alphaj, divdamp
public sigramplow, sigramphigh, amxlsq, dvmodmin, siburbanfrac, cqmix
public intsch_mode
public diag, localhist, unlimitedhist, synchist, amipo3
public save_aerosols, save_pbl, save_cloud, save_land, save_maxmin
public save_ocean, save_radiation, save_urban, save_carbon, save_river
public diaglevel_aerosols, diaglevel_pbl, diaglevel_cloud, diaglevel_land, diaglevel_maxmin
public diaglevel_ocean, diaglevel_radiation, diaglevel_urban, diaglevel_carbon, diaglevel_river
public diaglevel_pop
public procmode, compression
public nud_period, mins_rad, nalpha, jalbfix, irest, nwrite
public nstagin, nstaguin
public hp_output, surf_cordex, surf_windfarm
public ensemble_mode, ensemble_period, ensemble_rsfactor

integer, save :: ngwd=-5, nrungcm=-1, newtop=1
integer, save :: kountr=0, nrad=4, nvmix=3, nlocal=6
integer, save :: nhstest=0, namip=0, nspecial=0, newrough=0, nsib=3
integer, save :: ntaft=2, ntsea=6, ntsur=6, ntsur2, lgwd=0, newztsea=1, nglacier=1, mbd=0
integer, save :: mbd_mlo=0, nbd=0, kbotdav=4, kbotu=0, nud_p=0, nud_q=0, nud_t=0, nud_uv=1
integer, save :: nud_hrs=24, nudu_hrs=0, ktau, ndi=1, ndi2=0, ntau, nperavg=-99, nperday, nperhr
integer, save :: nmaxpr=99, nlv=1, ia=1, ib=1, ja=1, jb=1, id=1, jd=1, idjd
integer, save :: io_in=1, io_out=1, io_rest=1
integer, save :: nwt=-99, nrun=0, nextout=3, m_fly=4, nsemble=0, tbave=1
integer, save :: nurban=0, nmr=0, nmlo=0, ktopdav=0, nud_sst=0, nud_sss=0, kbotmlo=-1000, ktopmlo=1
integer, save :: mloalpha=0, nud_ouv=0, nud_sfh=0, kblock=-1, rescrn=0, knh=-1, iaero=0
integer, save :: nud_aero=0, mbd_maxscale=3000, mbd_maxgrid=999999, mbd_maxscale_mlo=3000, nriver=0
integer, save :: leap=0, nbarewet=0, nsigmf=1, qg_fix=2
integer, save :: procmode=0, compression=1
integer, save :: nud_period=-1, mins_rad=-1, nalpha=1, jalbfix=1, irest=1, nwrite=0
integer, save :: nstagin=0, nstaguin=0, intsch_mode=-1
integer, save :: hp_output=0, surf_cordex=0, surf_windfarm=0
integer, save :: ensemble_mode=0, ensemble_period=720
integer, save :: diaglevel_aerosols=0, diaglevel_pbl=0, diaglevel_cloud=0, diaglevel_land=0, diaglevel_maxmin=0
integer, save :: diaglevel_ocean=0, diaglevel_radiation=0, diaglevel_urban=0, diaglevel_carbon=0, diaglevel_river=0
integer, save :: diaglevel_pop=0
!integer, save :: filemode=0, ioreaders=-1
real, save :: qgmin=1.e-6, zo_clearing=0.
real, save :: av_vmod=0.7, vmodmin=0.2, snmin=0.11, tss_sh=1., charnock=0.018, chn10=0.00125, zobgin=0.02
real, save :: rlongdn=0., rlongdx=0., rlatdn=0., rlatdx=0., ds=0., dt=0., dtin=0., panfg=4., panzo=0.001
real, save :: bpyear=0., helim=800., fc2=1., sigbot_gwd=0., alphaj=1.e-6, divdamp=450.
real, save :: sigramplow=0., sigramphigh=0., amxlsq=100., dvmodmin=1., siburbanfrac=1., cqmix=2.5
real, save :: ensemble_rsfactor=0.1
logical, save :: diag=.false., localhist=.false., unlimitedhist=.true., synchist=.false., amipo3=.false.
logical, save :: save_aerosols=.true., save_pbl=.true., save_cloud=.true., save_land=.true., save_maxmin=.true.
logical, save :: save_ocean=.true., save_radiation=.true., save_urban=.true., save_carbon=.true., save_river=.true.
!logical, save :: pio=.false., mpiio=.true., npio=.false., useiobuffer=.false.

!$acc declare create(vmodmin,sigbot_gwd,fc2,dt,alphaj,ngwd,iaero,ds,nmr,diag)
!$acc declare create(qgmin,nlocal,cqmix)

end module parm_m
