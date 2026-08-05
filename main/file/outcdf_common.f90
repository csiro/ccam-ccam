! Conformal Cubic Atmospheric Model
    
! Copyright 2026 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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

! Common data and tools for outcdf
    
module outcdf_common
    
private
public month, cordex_levels, height_levels
public cordex_level_data, height_level_data
public mslp, cordex_name, bisect, outparam
public version

character(len=3), dimension(12), parameter :: month = (/'jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec'/)
integer, parameter :: cordex_levels = 17
integer, dimension(cordex_levels) :: cordex_level_data = &
    (/ 1000, 925, 850, 700, 600, 500, 400, 300, 250, 200, 150, 100, 70, 50, 30, 20, 10 /)
integer, parameter :: height_levels = 6
integer, dimension(height_levels) :: height_level_data = &
    (/ 50, 100, 150, 200, 250, 300 /)

include 'version.h'                              ! Model version data

contains

subroutine cordex_name(lname,stringa,press_level,stringb)

use cc_mpi, only : ccmpi_abort

implicit none

integer, intent(in) :: press_level
character(len=*), intent(out) :: lname
character(len=*), intent(in) :: stringa
character(len=*), intent(in), optional :: stringb

if ( present(stringb) ) then
  if ( press_level>=1000 ) then
    write(lname,'(A,I4.4,A)') stringa,press_level,stringb
  else if (press_level>=100 ) then
    write(lname,'(A,I3.3,A)') stringa,press_level,stringb
  else if ( press_level>=10 ) then
    write(lname,'(A,I2.2,A)') stringa,press_level,stringb
  else if ( press_level>=1 ) then
    write(lname,'(A,I1.1,A)') stringa,press_level,stringb
  else
    write(6,*) "ERROR: Unexpected output pressure level in cordex_name"
    call ccmpi_abort(-1)
  end if
else
  if ( press_level>=1000 ) then
    write(lname,'(A,I4.4)') stringa,press_level
  else if (press_level>=100 ) then
    write(lname,'(A,I3.3)') stringa,press_level
  else if ( press_level>=10 ) then
    write(lname,'(A,I2.2)') stringa,press_level
  else if ( press_level>=1 ) then
    write(lname,'(A,I1.1)') stringa,press_level
  else
    write(6,*) "ERROR: Unexpected output pressure level in cordex_name"
    call ccmpi_abort(-1)
  end if
end if

return
end subroutine cordex_name

subroutine mslp(pmsl,psl,zs,t)

use cc_mpi, only : mydiag, myid
use const_phys
use newmpar_m
use parm_m
use sigs_m

implicit none
! this one will ignore negative zs (i.e. over the ocean)

integer, save :: lev = -1
integer, dimension(1) :: pos
!real c, conr, con
real, dimension(ifull), intent(out) :: pmsl
real, dimension(ifull), intent(in) :: psl, zs
real, dimension(ifull) :: phi1, tsurf, tav,  dlnps
real, dimension(:,:), intent(in) :: t
      
!c = grav/stdlapse
!conr = c/rdry
if ( lev<0 ) then
  pos = minloc(abs(sig-0.9),sig>=0.9)
  lev = pos(1)
  if ( myid==0 .and. nmaxpr==1 ) then
    write(6,*) "Reducing ps to MSLP with lev,sig ",lev,sig(lev) 
  end if
end if
!con = sig(lev)**(rdry/c)/c
      
phi1(:) = t(1:ifull,lev)*rdry*(1.-sig(lev))/sig(lev) ! phi of sig(lev) above sfce
tsurf(:) = t(1:ifull,lev)+phi1(:)*stdlapse/grav
tav(:) = tsurf(:)+zs(1:ifull)*.5*stdlapse/grav
dlnps(:) = zs(1:ifull)/(rdry*tav(:))
pmsl(:) = 1.e5*exp(psl(:)+dlnps(:))
      
if ( nmaxpr==1 .and. mydiag ) then
  write(6,*) 'meth,lev,sig(lev) ',1,lev,sig(lev)
  write(6,*) 'zs,t_lev,psl,pmsl ',zs(idjd),t(idjd,lev),psl(idjd),pmsl(idjd)
end if
      
return
end subroutine mslp

subroutine outparam(idnc)

use aerointerface                          ! Aerosol interface 
use cc_mpi                                 ! CC MPI routines
use estab                                  ! Liquid saturation function
use infile                                 ! Input file routines
use kuocom_m                               ! JLM convection
use mlodynamics                            ! Ocean dynamics
use module_aux_rad                         ! Additional cloud and radiation routines
use module_ctrl_microphysics               ! Interface for cloud microphysics
use nharrs_m                               ! Non-hydrostatic atmosphere arrays
use ozoneread                              ! Ozone input routines
use parm_m                                 ! Model configuration
use parmdyn_m                              ! Dynamics parameters
use parmhdff_m                             ! Horizontal diffusion parameters
use parmhor_m                              ! Horizontal advection parameters
use parmvert_m                             ! Vertical advection parameters
use river                                  ! River routing
use seaesfrad_m                            ! SEA-ESF radiation
use sflux_m                                ! Surface flux routines
use tkeeps                                 ! TKE-EPS boundary layer

implicit none

integer, intent(in) :: idnc
integer namipo3, nalways_mspeca, ndo_co2_10um, ndo_quench
integer nremain_rayleigh_bug, nuse_rad_year

! main
call ccnf_put_attg(idnc,'aeroindir',aeroindir)
if ( always_mspeca ) then
  nalways_mspeca = 1
else
  nalways_mspeca = 0
end if
call ccnf_put_attg(idnc,'always_mspeca',nalways_mspeca)
if ( amipo3 ) then
  namipo3 = 1
else
  namipo3 = 0
end if
call ccnf_put_attg(idnc,'amipo3',namipo3)
call ccnf_put_attg(idnc,'av_vmod',av_vmod)
call ccnf_put_attg(idnc,'ch_dust',ch_dust)
call ccnf_put_attg(idnc,'charnock',charnock)
call ccnf_put_attg(idnc,'chn10',chn10)
call ccnf_put_attg(idnc,'cordex_fix',cordex_fix)
call ccnf_put_attg(idnc,'divdamp',divdamp)
call ccnf_put_attg(idnc,'ensemble_mode',ensemble_mode)
call ccnf_put_attg(idnc,'ensemble_period',ensemble_period)
call ccnf_put_attg(idnc,'ensemble_rsfactor',ensemble_rsfactor)
call ccnf_put_attg(idnc,'epsf',epsf)
call ccnf_put_attg(idnc,'epsh',epsh)
call ccnf_put_attg(idnc,'epsp',epsp)
call ccnf_put_attg(idnc,'epsu',epsu)
call ccnf_put_attg(idnc,'estab_bug_fix',estab_bug_fix)
call ccnf_put_attg(idnc,'helmmeth',helmmeth)
call ccnf_put_attg(idnc,'hp_output',hp_output)
call ccnf_put_attg(idnc,'iaero',iaero)   
call ccnf_put_attg(idnc,'iceradmethod',iceradmethod)
call ccnf_put_attg(idnc,'intsch_mode',intsch_mode)
call ccnf_put_attg(idnc,'jalbfix',jalbfix)
call ccnf_put_attg(idnc,'kblock',kblock)
call ccnf_put_attg(idnc,'kbotdav',kbotdav)
call ccnf_put_attg(idnc,'kbotmlo',kbotmlo)
call ccnf_put_attg(idnc,'khdif',khdif)
call ccnf_put_attg(idnc,'khor',khor)
call ccnf_put_attg(idnc,'knh',knh)
call ccnf_put_attg(idnc,'ktopdav',ktopdav)
call ccnf_put_attg(idnc,'ktopmlo',ktopmlo)
call ccnf_put_attg(idnc,'leap',leap)
call ccnf_put_attg(idnc,'liqradmethod',liqradmethod)    
call ccnf_put_attg(idnc,'lgwd',lgwd)
call ccnf_put_attg(idnc,'m_fly',m_fly)
call ccnf_put_attg(idnc,'maxcolour',maxcolour)
call ccnf_put_attg(idnc,'maxuv',maxuv)
call ccnf_put_attg(idnc,'mbd',mbd)
call ccnf_put_attg(idnc,'mbd_maxgrid',mbd_maxgrid)
call ccnf_put_attg(idnc,'mbd_maxscale',mbd_maxscale)
call ccnf_put_attg(idnc,'mbd_maxscale_mlo',mbd_maxscale_mlo)
call ccnf_put_attg(idnc,'mbd_mlo',mbd_mlo)
call ccnf_put_attg(idnc,'mex',mex)
call ccnf_put_attg(idnc,'mfix',mfix)
call ccnf_put_attg(idnc,'mfix_aero',mfix_aero)
call ccnf_put_attg(idnc,'mfix_qg',mfix_qg)
call ccnf_put_attg(idnc,'mfix_t',mfix_t)
call ccnf_put_attg(idnc,'mfix_tr',mfix_tr)
call ccnf_put_attg(idnc,'mh_bs',mh_bs)
call ccnf_put_attg(idnc,'mloalpha',mloalpha)
call ccnf_put_attg(idnc,'mup',mup)
call ccnf_put_attg(idnc,'nalpha',nalpha)
call ccnf_put_attg(idnc,'namip',namip)
call ccnf_put_attg(idnc,'nbarewet',nbarewet)
call ccnf_put_attg(idnc,'nbd',nbd)
call ccnf_put_attg(idnc,'newrough',newrough)
call ccnf_put_attg(idnc,'newtop',newtop)
call ccnf_put_attg(idnc,'newztsea',newztsea)
call ccnf_put_attg(idnc,'nglacier',nglacier)
call ccnf_put_attg(idnc,'nh',nh)
call ccnf_put_attg(idnc,'nhor',nhor)
call ccnf_put_attg(idnc,'nhorjlm',nhorjlm)
call ccnf_put_attg(idnc,'nhorps',nhorps)
call ccnf_put_attg(idnc,'nhstest',nhstest)
call ccnf_put_attg(idnc,'nlocal',nlocal)
call ccnf_put_attg(idnc,'nmlo',nmlo)
call ccnf_put_attg(idnc,'nrad',nrad)
call ccnf_put_attg(idnc,'nritch_t',nritch_t)
call ccnf_put_attg(idnc,'nrungcm',nrungcm)
call ccnf_put_attg(idnc,'nsemble',nsemble)
call ccnf_put_attg(idnc,'nsib',nsib)
call ccnf_put_attg(idnc,'nsigmf',nsigmf)
call ccnf_put_attg(idnc,'nspecial',nspecial)
call ccnf_put_attg(idnc,'nstagu',nstagu)
call ccnf_put_attg(idnc,'nt_adv',nt_adv)
call ccnf_put_attg(idnc,'ntaft',ntaft)
call ccnf_put_attg(idnc,'ntbar',ntbar)
call ccnf_put_attg(idnc,'ntsea',ntsea)
call ccnf_put_attg(idnc,'ntsur',ntsur)
call ccnf_put_attg(idnc,'ntvd',ntvd)
call ccnf_put_attg(idnc,'nud_aero',nud_aero)
call ccnf_put_attg(idnc,'nud_hrs',nud_hrs)
call ccnf_put_attg(idnc,'nud_ouv',nud_ouv)
call ccnf_put_attg(idnc,'nud_p',nud_p)
call ccnf_put_attg(idnc,'nud_period',nud_period)
call ccnf_put_attg(idnc,'nud_q',nud_q)
call ccnf_put_attg(idnc,'nud_sfh',nud_sfh)
call ccnf_put_attg(idnc,'nud_sss',nud_sss)    
call ccnf_put_attg(idnc,'nud_sst',nud_sst)
call ccnf_put_attg(idnc,'nud_t',nud_t)
call ccnf_put_attg(idnc,'nud_uv',nud_uv)
call ccnf_put_attg(idnc,'nudu_hrs',nudu_hrs)
call ccnf_put_attg(idnc,'nurban',nurban)
call ccnf_put_attg(idnc,'nvmix',nvmix)
call ccnf_put_attg(idnc,'nxtrrho',nxtrrho)
!call ccnf_put_attg(idnc,'ol',ol)
call ccnf_put_attg(idnc,'panfg',panfg)
call ccnf_put_attg(idnc,'panzo',panzo)
call ccnf_put_attg(idnc,'pil_single',pil_single)
call ccnf_put_attg(idnc,'process_rate_mode',process_rate_mode)
call ccnf_put_attg(idnc,'precon',precon)
call ccnf_put_attg(idnc,'qg_fix',qg_fix)
call ccnf_put_attg(idnc,'rescrn',rescrn)
call ccnf_put_attg(idnc,'restol',restol)
call ccnf_put_attg(idnc,'rhsat',rhsat)
call ccnf_put_attg(idnc,'sigramphigh',sigramphigh)
call ccnf_put_attg(idnc,'sigramplow',sigramplow)
call ccnf_put_attg(idnc,'snmin',snmin)
call ccnf_put_attg(idnc,'tbave',tbave)
call ccnf_put_attg(idnc,'tbave10',tbave10)
call ccnf_put_attg(idnc,'tss_sh',tss_sh)
call ccnf_put_attg(idnc,'vmodmin',vmodmin)
call ccnf_put_attg(idnc,'zobgin',zobgin)
call ccnf_put_attg(idnc,'zo_clearing',zo_clearing)

! radiation and aerosol
call ccnf_put_attg(idnc,'aero_split',aero_split)
call ccnf_put_attg(idnc,'aeroindir',aeroindir)
call ccnf_put_attg(idnc,'aerosol_u10',aerosol_u10)    
call ccnf_put_attg(idnc,'bpyear',bpyear)
call ccnf_put_attg(idnc,'carbmtn',carbmtn)
call ccnf_put_attg(idnc,'carbonradmethod',carbonradmethod)
call ccnf_put_attg(idnc,'ch_dust',ch_dust)
call ccnf_put_attg(idnc,'continuum_form',continuum_form)
call ccnf_put_attg(idnc,'csolar',csolar)
if ( do_co2_10um ) then
  ndo_co2_10um = 1
else
  ndo_co2_10um = 0
end if      
call ccnf_put_attg(idnc,'do_co2_10um',ndo_co2_10um)
if ( do_quench ) then
  ndo_quench = 1
else
  ndo_quench = 0
end if
call ccnf_put_attg(idnc,'do_quench',ndo_quench)
call ccnf_put_attg(idnc,'dustradmethod',dustradmethod)
call ccnf_put_attg(idnc,'enhanceu10',enhanceu10)
call ccnf_put_attg(idnc,'iceradmethod',iceradmethod)
call ccnf_put_attg(idnc,'linecatalog_form',linecatalog_form)
call ccnf_put_attg(idnc,'liqradmethod',liqradmethod)
call ccnf_put_attg(idnc,'lwem_form',lwem_form)
call ccnf_put_attg(idnc,'mins_rad',mins_rad)
call ccnf_put_attg(idnc,'o3_vert_interpolate',o3_vert_interpolate)
call ccnf_put_attg(idnc,'qgmin',qgmin)
call ccnf_put_attg(idnc,'rad_year',rad_year)    
if ( remain_rayleigh_bug ) then
  nremain_rayleigh_bug = 1
else
  nremain_rayleigh_bug = 0
end if  
call ccnf_put_attg(idnc,'remain_rayleigh_bug',nremain_rayleigh_bug)
call ccnf_put_attg(idnc,'saltlargemtn',saltlargemtn)
call ccnf_put_attg(idnc,'saltsmallmtn',saltsmallmtn)
call ccnf_put_attg(idnc,'seasaltradmethod',seasaltradmethod)
call ccnf_put_attg(idnc,'siglow',siglow)
call ccnf_put_attg(idnc,'sigmid',sigmid)
call ccnf_put_attg(idnc,'so4mtn',so4mtn)
call ccnf_put_attg(idnc,'so4radmethod',so4radmethod)
call ccnf_put_attg(idnc,'sw_diff_streams',sw_diff_streams)
call ccnf_put_attg(idnc,'sw_resolution',sw_resolution)
if ( use_rad_year ) then
  nuse_rad_year = 1
else
  nuse_rad_year = 0
end if
call ccnf_put_attg(idnc,'use_rad_year',nuse_rad_year)
call ccnf_put_attg(idnc,'zvolcemi',zvolcemi)
    
! convection and cloud microphysics
call ccnf_put_attg(idnc,'acon',acon)
call ccnf_put_attg(idnc,'alflnd',alflnd)
call ccnf_put_attg(idnc,'alfsea',alfsea)
call ccnf_put_attg(idnc,'bcon',bcon)
call ccnf_put_attg(idnc,'cld_decay',cld_decay)
call ccnf_put_attg(idnc,'cldh_lnd',cldh_lnd)
call ccnf_put_attg(idnc,'cldl_lnd',cldl_lnd)
call ccnf_put_attg(idnc,'cldm_lnd',cldm_lnd)
call ccnf_put_attg(idnc,'cldh_sea',cldh_sea)
call ccnf_put_attg(idnc,'cldl_sea',cldl_sea)
call ccnf_put_attg(idnc,'cldm_sea',cldm_sea)
call ccnf_put_attg(idnc,'cloud_aerosol_mode',cloud_aerosol_mode)
call ccnf_put_attg(idnc,'cloud_ice_method',cloud_ice_method)
call ccnf_put_attg(idnc,'convfact',convfact)
call ccnf_put_attg(idnc,'convtime',convtime)
call ccnf_put_attg(idnc,'detrain',detrain)
call ccnf_put_attg(idnc,'detrainx',detrainx)
call ccnf_put_attg(idnc,'dsig2',dsig2)
call ccnf_put_attg(idnc,'dsig4',dsig4)
call ccnf_put_attg(idnc,'entrain',entrain)
call ccnf_put_attg(idnc,'fldown',fldown)
call ccnf_put_attg(idnc,'gf_imid',gf_imid)
call ccnf_put_attg(idnc,'iterconv',iterconv)
call ccnf_put_attg(idnc,'ksc',ksc)
call ccnf_put_attg(idnc,'kscmom',kscmom)
call ccnf_put_attg(idnc,'kscsea',kscsea)
call ccnf_put_attg(idnc,'ldr',ldr)
call ccnf_put_attg(idnc,'leon_snowmeth',leon_snowmeth)
call ccnf_put_attg(idnc,'lin_adv',lin_adv)
call ccnf_put_attg(idnc,'maxlintime',maxlintime)
call ccnf_put_attg(idnc,'mbase',mbase)
call ccnf_put_attg(idnc,'mdelay',mdelay)
call ccnf_put_attg(idnc,'methdetr',methdetr)
call ccnf_put_attg(idnc,'methprec',methprec)
call ccnf_put_attg(idnc,'nbase',nbase)
call ccnf_put_attg(idnc,'ncldia',nclddia)
call ccnf_put_attg(idnc,'ncloud',ncloud)
call ccnf_put_attg(idnc,'ncvcloud',ncvcloud)
call ccnf_put_attg(idnc,'ncvmix',ncvmix)
call ccnf_put_attg(idnc,'nevapcc',nevapcc)
call ccnf_put_attg(idnc,'nevapls',nevapls)
call ccnf_put_attg(idnc,'nkuo',nkuo)
call ccnf_put_attg(idnc,'nmr',nmr)
call ccnf_put_attg(idnc,'nrhcrit',nrhcrit)
call ccnf_put_attg(idnc,'nscheme',nscheme)
call ccnf_put_attg(idnc,'nstab_cld',nstab_cld)
call ccnf_put_attg(idnc,'nuvconv',nuvconv)
call ccnf_put_attg(idnc,'qfg_max',qfg_max)
call ccnf_put_attg(idnc,'qlg_max',qlg_max)
call ccnf_put_attg(idnc,'rcm',rcm)
call ccnf_put_attg(idnc,'rcrit_l',rcrit_l)
call ccnf_put_attg(idnc,'rcrit_s',rcrit_s)
call ccnf_put_attg(idnc,'rhcv',rhcv)
call ccnf_put_attg(idnc,'rhmois',rhmois)
call ccnf_put_attg(idnc,'rhsat',rhsat)
call ccnf_put_attg(idnc,'shaltime',shaltime)
call ccnf_put_attg(idnc,'sig_ct',sig_ct)        
call ccnf_put_attg(idnc,'sigcb',sigcb)
call ccnf_put_attg(idnc,'sigksct',sigcll)
call ccnf_put_attg(idnc,'sigkscb',sigkscb)
call ccnf_put_attg(idnc,'sigksct',sigksct)
call ccnf_put_attg(idnc,'tied_con',tied_con)
call ccnf_put_attg(idnc,'tied_over',tied_over)
call ccnf_put_attg(idnc,'tied_rh',tied_rh)
call ccnf_put_attg(idnc,'tiedtke_form',tiedtke_form)
call ccnf_put_attg(idnc,'vdeposition_mode',vdeposition_mode)

! boundary layer turbulence and gravity wave
call ccnf_put_attg(idnc,'alphaj',alphaj)    
call ccnf_put_attg(idnc,'amxlsq',amxlsq)
call ccnf_put_attg(idnc,'b1',b1)
call ccnf_put_attg(idnc,'b2',b2)
call ccnf_put_attg(idnc,'be',be)
call ccnf_put_attg(idnc,'buoymeth',buoymeth)
call ccnf_put_attg(idnc,'ce0',ce0)
call ccnf_put_attg(idnc,'ce1',ce1)
call ccnf_put_attg(idnc,'ce2',ce2)
call ccnf_put_attg(idnc,'ce3',ce3)
call ccnf_put_attg(idnc,'cm0',cm0)
call ccnf_put_attg(idnc,'cqmix',cqmix)
call ccnf_put_attg(idnc,'dtrc0',dtrc0)
call ccnf_put_attg(idnc,'dvmodmin',dvmodmin)
call ccnf_put_attg(idnc,'ent_min',ent_min)
call ccnf_put_attg(idnc,'ent0',ent0)
call ccnf_put_attg(idnc,'ent1',ent1)
call ccnf_put_attg(idnc,'entc0',entc0)
call ccnf_put_attg(idnc,'ezmin',ezmin)
call ccnf_put_attg(idnc,'fc2',fc2)
call ccnf_put_attg(idnc,'helim',helim)
call ccnf_put_attg(idnc,'m0',m0)
call ccnf_put_attg(idnc,'maxdts',maxdts)
call ccnf_put_attg(idnc,'maxl',maxl)
call ccnf_put_attg(idnc,'mfbeta',mfbeta)
call ccnf_put_attg(idnc,'mineps',mineps)
call ccnf_put_attg(idnc,'minl',minl)
call ccnf_put_attg(idnc,'mintke',mintke)
call ccnf_put_attg(idnc,'ngwd',ngwd)
call ccnf_put_attg(idnc,'qcmf',qcmf)
call ccnf_put_attg(idnc,'sigbot_gwd',sigbot_gwd)    
call ccnf_put_attg(idnc,'stabmeth',stabmeth)
call ccnf_put_attg(idnc,'tcalmeth',tcalmeth)
call ccnf_put_attg(idnc,'tkemeth',tkemeth)
call ccnf_put_attg(idnc,'ugs_meth',ugs_meth)
call ccnf_put_attg(idnc,'wg_prob',wg_prob)
call ccnf_put_attg(idnc,'wg_tau',wg_tau)
    
! land, urban and carbon
call ccnf_put_attg(idnc,'ateb_ac_coolcap',ateb_ac_coolcap)
call ccnf_put_attg(idnc,'ateb_ac_deltat',ateb_ac_deltat)
call ccnf_put_attg(idnc,'ateb_ac_heatcap',ateb_ac_heatcap)
call ccnf_put_attg(idnc,'ateb_acfactor',ateb_acfactor)
call ccnf_put_attg(idnc,'ateb_lwintmeth',ateb_lwintmeth)   
call ccnf_put_attg(idnc,'ateb_cvcoeffmeth',ateb_cvcoeffmeth)
call ccnf_put_attg(idnc,'ateb_infilmeth',ateb_infilmeth)
call ccnf_put_attg(idnc,'ateb_intairtmeth',intairtmeth)
call ccnf_put_attg(idnc,'ateb_intmassmeth',intmassmeth)
call ccnf_put_attg(idnc,'ateb_maxrdsn',ateb_maxrdsn)
call ccnf_put_attg(idnc,'ateb_maxrdwater',ateb_maxrdwater)
call ccnf_put_attg(idnc,'ateb_maxrfsn',ateb_maxrfsn)
call ccnf_put_attg(idnc,'ateb_maxrfwater',ateb_maxrfwater)
call ccnf_put_attg(idnc,'ateb_maxsnowalpha',ateb_maxsnowalpha)
call ccnf_put_attg(idnc,'ateb_maxsnowden',ateb_maxsnowden)
call ccnf_put_attg(idnc,'ateb_maxvwatf',ateb_maxvwatf)
call ccnf_put_attg(idnc,'ateb_minsnowalpha',ateb_minsnowalpha)
call ccnf_put_attg(idnc,'ateb_minsnowden',ateb_minsnowden)
call ccnf_put_attg(idnc,'ateb_ncyits',ateb_ncyits)
call ccnf_put_attg(idnc,'ateb_nfgits',ateb_nfgits)
call ccnf_put_attg(idnc,'ateb_nrefl',ateb_nrefl)
call ccnf_put_attg(idnc,'ateb_refheight',ateb_refheight)
call ccnf_put_attg(idnc,'ateb_resmeth',ateb_resmeth)
call ccnf_put_attg(idnc,'ateb_snowemiss',ateb_snowemiss)
call ccnf_put_attg(idnc,'ateb_statsmeth',ateb_statsmeth)
call ccnf_put_attg(idnc,'ateb_tol',ateb_tol)
call ccnf_put_attg(idnc,'ateb_wbrelaxc',ateb_wbrelaxc)
call ccnf_put_attg(idnc,'ateb_wbrelaxr',ateb_wbrelaxr)
call ccnf_put_attg(idnc,'ateb_zocanyon',zocanyon)
call ccnf_put_attg(idnc,'ateb_zohmeth',ateb_zohmeth)
call ccnf_put_attg(idnc,'ateb_zomratio',ateb_zomratio)
call ccnf_put_attg(idnc,'ateb_zoroof',zoroof)
call ccnf_put_attg(idnc,'ateb_zosnow',ateb_zosnow)
call ccnf_put_attg(idnc,'cable_enablefao',cable_enablefao)
call ccnf_put_attg(idnc,'cable_gw_model',cable_gw_model)
call ccnf_put_attg(idnc,'cable_litter',cable_litter)
call ccnf_put_attg(idnc,'cable_roughness',cable_roughness)
!call ccnf_put_attg(idnc,'cable_runoff',cable_runoff)
!call ccnf_put_attg(idnc,'cable_soilevap',cable_soilevap)
!call ccnf_put_attg(idnc,'cable_thermfix',cable_thermfix)
call ccnf_put_attg(idnc,'cable_potev',cable_potev)
call ccnf_put_attg(idnc,'cable_pop',cable_pop)    
!call ccnf_put_attg(idnc,'cable_version',cable_version)
call ccnf_put_attg(idnc,'ccycle',ccycle)
call ccnf_put_attg(idnc,'freshwaterlake_fix',freshwaterlake_fix)
call ccnf_put_attg(idnc,'fwsoil_switch',fwsoil_switch)
call ccnf_put_attg(idnc,'gs_switch',gs_switch)
call ccnf_put_attg(idnc,'proglai',proglai)
call ccnf_put_attg(idnc,'progvcmax',progvcmax)
call ccnf_put_attg(idnc,'smrf_switch',smrf_switch)
call ccnf_put_attg(idnc,'strf_switch',strf_switch)
call ccnf_put_attg(idnc,'siburbanfrac',siburbanfrac)
call ccnf_put_attg(idnc,'soil_struc',soil_struc)
call ccnf_put_attg(idnc,'wt_transport',wt_transport)
    
! ocean
call ccnf_put_attg(idnc,'alphanir_seaice',alphanir_seaice)
call ccnf_put_attg(idnc,'alphanir_seasnw',alphanir_seasnw)        
call ccnf_put_attg(idnc,'alphavis_seaice',alphavis_seaice)    
call ccnf_put_attg(idnc,'alphavis_seasnw',alphavis_seasnw)    
call ccnf_put_attg(idnc,'basinmd',basinmd)
call ccnf_put_attg(idnc,'delwater',delwater)
call ccnf_put_attg(idnc,'factchseaice',factchseaice)
call ccnf_put_attg(idnc,'fluxwgt',fluxwgt)
call ccnf_put_attg(idnc,'kemaxdt',kemaxdt)
call ccnf_put_attg(idnc,'maxsal',maxsal)
call ccnf_put_attg(idnc,'mindep',mindep)
call ccnf_put_attg(idnc,'minsal',minsal)
call ccnf_put_attg(idnc,'minwater',minwater)
call ccnf_put_attg(idnc,'mlo_adjeta',mlo_adjeta)
call ccnf_put_attg(idnc,'mlo_bs',mlo_bs)
call ccnf_put_attg(idnc,'mlo_step',mlo_step)
call ccnf_put_attg(idnc,'mlo_uvcoupl',mlo_uvcoupl)
call ccnf_put_attg(idnc,'mlodiff',mlodiff)
call ccnf_put_attg(idnc,'mlodiff_numits',mlodiff_numits)
call ccnf_put_attg(idnc,'mlodps',mlodps)
call ccnf_put_attg(idnc,'mloiceadv',mloiceadv)
call ccnf_put_attg(idnc,'mlointschf',mlointschf)
call ccnf_put_attg(idnc,'mlojacobi',mlojacobi)
call ccnf_put_attg(idnc,'mlomfix',mlomfix)
call ccnf_put_attg(idnc,'mlosigma',mlosigma)
call ccnf_put_attg(idnc,'mlontvd',mlontvd)
call ccnf_put_attg(idnc,'mprecond',mprecond)
call ccnf_put_attg(idnc,'mstagf',mstagf)
call ccnf_put_attg(idnc,'mxd',mxd)
call ccnf_put_attg(idnc,'nodrift',nodrift)
call ccnf_put_attg(idnc,'oclosure',oclosure)
call ccnf_put_attg(idnc,'ocnepr',ocnepr)
call ccnf_put_attg(idnc,'ocneps',ocneps)
call ccnf_put_attg(idnc,'ocnsmag',ocnsmag)
call ccnf_put_attg(idnc,'omaxl',omaxl)
call ccnf_put_attg(idnc,'omineps',omineps)
call ccnf_put_attg(idnc,'omink',omink)
call ccnf_put_attg(idnc,'ominl',ominl)
call ccnf_put_attg(idnc,'otaumode',otaumode)
call ccnf_put_attg(idnc,'rivercoeff',rivercoeff)
call ccnf_put_attg(idnc,'rivermd',rivermd)
call ccnf_put_attg(idnc,'usepice',usepice)
call ccnf_put_attg(idnc,'usetide',usetide)
call ccnf_put_attg(idnc,'zomode',zomode)
call ccnf_put_attg(idnc,'zoseaice',zoseaice)
    
! ensemble data
if ( driving_model_id /= ' ' ) then
  call ccnf_put_attg(idnc,'driving_model_id',trim(driving_model_id))
end if
if ( driving_institution_id /= ' ' ) then
  call ccnf_put_attg(idnc,'driving_institution_id',trim(driving_institution_id))
end if
if ( driving_model_ensemble_number /= ' ' ) then
  call ccnf_put_attg(idnc,'driving_model_ensemble_number',trim(driving_model_ensemble_number))
end if
if ( driving_experiment_name /= ' ' ) then
  call ccnf_put_attg(idnc,'driving_experiment_name',trim(driving_experiment_name))
end if 

end subroutine outparam

! Find pressure level
pure function bisect(press_target, ps, sig) result(ans)
real, intent(in) :: press_target, ps
real, dimension(:), intent(in) :: sig
integer :: ans
integer a, b, i, kx

kx = size(sig)
a = 1
b = kx
do while ( b-a > 1 )
  i = (a+b)/2
  if ( press_target > ps*sig(i) ) then
    b = i
  else
    a = i
  end if
end do
if ( ps*sig(a)>=press_target .and. ps*sig(b)>=press_target ) then
  ans = b
else
  ans = a
end if
ans = min( ans, kx-1 )

end function bisect

end module outcdf_common
