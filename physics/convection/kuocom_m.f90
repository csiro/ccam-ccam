! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2025 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
module kuocom_m

implicit none

private
public iterconv,ksc,kscmom,kscsea,kuocb,ldr,mbase
public mdelay,methdetr,methprec,nbase,nclddia,ncvcloud,ncvmix
public nevapcc,nevapls,nkuo,nrhcrit,nstab_cld,nuvconv,ncloud
public alflnd,alfsea ,cldh_lnd,cldm_lnd,cldl_lnd,cldh_sea,cldm_sea
public cldl_sea,convfact,convtime
public detrain,detrainx,dsig2,dsig4,entrain
public fldown,rhcv,rhmois,rhsat,shaltime
public sigcb,sigcll,sigkscb,sig_ct,sigksct
public tied_con,tied_over,tied_rh
public nscheme
public acon,bcon,rcm,rcrit_l,rcrit_s,cld_decay
public vdeposition_mode,tiedtke_form
public kbsav,ktsav
public convpsav,fluxtot
public mconv_save
public kuocom_init,kuocom_end

public mse_t1, mse_t2
public dmsedt_adv, dmsedt_rad, dmsedt_pbl

integer, save :: iterconv=3,ksc=-95,kscmom=1,kscsea=0,kuocb,ldr=1,mbase=101
integer, save :: mdelay=-1,methdetr=2,methprec=8,nbase=-4,nclddia=1,ncvcloud=0,ncvmix=0
integer, save :: nevapcc=0,nevapls=-4,nkuo=23,nrhcrit=10,nstab_cld=0,nuvconv=0,ncloud=0
integer, save :: vdeposition_mode=0,tiedtke_form=0, nscheme=1
integer, dimension(:), allocatable, target, save :: kbsav, ktsav
real, save :: alflnd=1.1,alfsea=1.1,cldh_lnd=95.,cldm_lnd=85.,cldl_lnd=75.,cldh_sea=95.,cldm_sea=90.
real, save :: cldl_sea=80.,convfact=1.02,convtime=0.33
real, save :: detrain=0.15,detrainx=0.,dsig2=0.15,dsig4=0.4,entrain=0.05
real, save :: fldown=0.6,rhcv=0.,rhmois=0.1,rhsat=1.,shaltime=0.
real, save :: sigcb=1.,sigcll=0.95,sigkscb=0.95,sig_ct=1.,sigksct=0.8
real, save :: tied_con=2.,tied_over=0.,tied_rh=0.75
real, save :: acon=0.2,bcon=0.07,rcm=0.92e-5,rcrit_l=0.75,rcrit_s=0.85,cld_decay=7200.
real, dimension(:,:), allocatable, save :: fluxtot
real, dimension(:,:), allocatable, save :: mconv_save
real, dimension(:), allocatable, save :: convpsav

real, dimension(:,:), allocatable, save :: mse_t1, mse_t2
real, dimension(:,:), allocatable, save :: dmsedt_adv, dmsedt_rad, dmsedt_pbl

contains

subroutine kuocom_init(ifull,kl)

implicit none

integer, intent(in) :: ifull,kl

allocate(kbsav(ifull),ktsav(ifull))
allocate(convpsav(ifull),fluxtot(ifull,kl))
allocate(mconv_save(ifull,kl))
allocate(mse_t1(ifull,kl))
allocate(mse_t2(ifull,kl))
allocate(dmsedt_adv(ifull,kl), dmsedt_rad(ifull,kl), dmsedt_pbl(ifull,kl))
kbsav=kl-1
ktsav=kl-1
convpsav=0.
fluxtot=0.
mconv_save=0.
mse_t1=0.
mse_t2=0.
dmsedt_adv=0.
dmsedt_rad=0.
dmsedt_pbl=0.

return
end subroutine kuocom_init

subroutine kuocom_end

implicit none

deallocate(kbsav,ktsav)
deallocate(convpsav,fluxtot)
deallocate(mconv_save)
deallocate(mse_t1,mse_t2)
deallocate(dmsedt_adv, dmsedt_rad, dmsedt_pbl)
return
end subroutine kuocom_end

end module kuocom_m
