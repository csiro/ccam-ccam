
! Conformal Cubic Atmospheric Model
    
! Copyright 2016 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
public vdeposition_mode,tiedtke_form

public alflnd,alfsea ,cldh_lnd,cldm_lnd,cldl_lnd,cldh_sea,cldm_sea
public cldl_sea,convfact,convtime
public detrain,detrainx,dsig2,dsig4,entrain
public fldown,rhcv,rhmois,rhsat,shaltime
public sigcb,sigcll,sigkscb,sig_ct,sigksct
public tied_con,tied_over,tied_rh
public acon,bcon,rcm,rcrit_l,rcrit_s,cld_decay

real, save :: alflnd=1.1                               ! land-weighting ratio for cloud bases         [1.15]
real, save :: alfsea=1.1                               ! sea-weighting ratio for cloud bases          [1.05]
real, save :: cldh_lnd=95.,cldm_lnd=85.,cldl_lnd=75.   !      for old nrhcrit=5     [95.,85.,75.]
real, save :: cldh_sea=95.,cldm_sea=90.,cldl_sea=80.   !      for old nrhcrit=5     [95.,90.,80.]
real, save :: convfact=1.02                            ! overshooting factor for mass flux            [1.]
!real, save :: convrh                                  ! convective rh    inhibitor                   [0.]
real, save :: convtime=.33                             ! adjustment time (h) of cu scheme             [.3]
real, save :: detrain=.15                              ! fraction of precip into detrainment          [.1]
real, save :: detrainx=0.                              ! fraction into detrainment for shallow clouds [1.]
real, save :: dsig2=.15                                ! delta-sigma2 for end of shallow clouds       [.1]
real, save :: dsig4=.4                                 ! delta-sigma4 for start of deep clouds        [.7]
real, save :: entrain=.05                              ! entrainment factor                           [0.]
real, save :: fldown=.6                                ! fraction of convective flux into downdraft   [.6]
integer, save :: iterconv=3                            ! number of iterations in convjlm              [1]
integer, save :: ksc=-95                               ! shallow convection switch (99 for Tiedtke on)[0]
integer, save :: kscmom=1                              ! shallow convection momentum switch (1 for on)[0]
integer, save :: kscsea=0                              ! 1 for doing Tiedtke only over sea            [0]
integer, save :: kuocb                                 ! level of min. cloud base, calc. from sigcb   [1]
integer, save :: ldr=1                                 ! ldr scheme options; 0 for off                [?]
integer, save :: mbase=101                             ! base test: 1 cfrac; 2 cfrac and omega        [0]
integer, save :: mdelay=-1                             ! convective delay time in secs                [0]
integer, save :: methdetr=2                            ! meth_shallow_detrainment for convjlm, 2 off  [2]
integer, save :: methprec=8                            ! meth_precip (deep_detrainment) for convjlm   [8]
integer, save :: nbase=-4                              ! type of base: 1 simple; 2 linear to sfce     [1]
integer, save :: nclddia=1                             ! conversion of RH to cloudiness, 0, 3, or     [5]
integer, save :: ncvcloud=0                            ! convective cloud enhancement in radrive      [0]
integer, save :: ncvmix=0                              ! cumulus mixing in vertmix                    [0]
integer, save :: nevapcc=0                             ! evap scheme of convective rain               [0]
integer, save :: nevapls=-4                            ! evap scheme of large-scale rain              [5]
integer, save :: nkuo=23                               ! convective scheme                            [23]
integer, save :: nrhcrit=10                            ! Hal's original 0; for jlm 7, 8               [8 ]
integer, save :: nstab_cld=0                           ! 0 off, 3 for stability-enhanced cll          [0]
integer, save :: nuvconv=0                             ! 0 off, 1 to turn on momentum mixing          [0]
real, save :: rhcv=0.                                  ! RH trigger for convective scheme             [.75]
real, save :: rhmois=.1                                ! used by conjob, convjlm for nevapcc=5        [.1]
real, save :: rhsat=1.                                 ! saturation trigger for large-scale rain      [1.]
real, save :: shaltime=0.                              ! shallow convection time scale (h)            [0.]
real, save :: sigcb=1.                                 ! sig value for base of convection scheme      [1.]
real, save :: sigcll=.95                               ! sig value of low cloud base (for cll)        [.95]
real, save :: sig_ct=1.                                ! min sig value of cloud tops (for convjlm)    [.8]
real, save :: sigkscb=.95                              ! for tiedtke shallow convection               [.98]
real, save :: sigksct=.8                               ! for tiedtke shallow convection               [.75]
real, save :: tied_con=2.                              ! tiedtke diffsn const. e.g. 25., 20. or       [6.]
real, save :: tied_over=0.                             ! tiedtke overshooting constant                [2.]
real, save :: tied_rh=.75                              ! tiedtke RH trigger:                          [.75]
integer, save :: ncloud=0                              ! New options for microphysics 0 = original    [0]
integer, save :: vdeposition_mode=0                    ! New options for ice deposition 0 = original  [0]
integer, save :: tiedtke_form=0                        ! Options for Tiedtke cloud formation          [0]

real, save :: acon=.2                                  ! Cloud fraction for non-precipitating convection  [.2]
real, save :: bcon=.07                                 ! Rate at which conv cloud frac increases with R   [.07]
real, save :: rcm=.92e-5                               ! Threshold cloud dropl R for coalescence to begin [1e-5]
real, save :: rcrit_l=.75                              ! rcrit_land for ldr newcloud                      [.75]
real, save :: rcrit_s=.85                              ! rcrit_sea  for ldr newcloud                      [.85]
real, save :: cld_decay=7200.                          ! time-decay factor for cirrus                     [7200.]

!$acc declare create(acon,bcon,ldr)

end module kuocom_m
