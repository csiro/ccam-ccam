c $Log$
c Revision 1.6  2006/02/23 22:49:54  mcg130
c CCAM-5-11
c
c Revision 1.2  2004/05/05 07:09:51  dix043
c JMcG change - May 2004
c
c Revision 1.1.1.1  2003/08/13 01:24:20  dix043
c Imported sources
c
c Revision 1.34  2003/02/26 04:01:15  rot032
c Bring Rcm value back to same as in HBG's version (7.5 um).
c
c Revision 1.33  2003/02/25 05:56:14  rot032
c Changes from LDR for aerosols and radiation
c
c Revision 1.32  2002/11/21 03:17:51  rot032
c Ice-cloud and aerosol radiative changes and new BC/OC emissions from LDR
c
c Revision 1.31  2001/10/12 02:13:45  rot032
c LDR changes, to bring sulfur cycle to run P76, and Mk2 OGCM fixes
c
c Revision 1.30  2000/11/14 03:11:37  rot032
c Energy/water balance for land-surface and sea-ice, plus tidy ups from HBG
c
c Revision 1.29  1999/06/30 05:29:37  rot032
c Add Ecol to this file.
c
c Revision 1.28  1998/12/10 00:55:39  ldr
c HBG changes to V5-1-21
c
c Revision 1.27  1998/05/26  06:08:06  ldr
c New mixed-phase cloud scheme from LDR (use plates option).
c
c Revision 1.26  1998/05/26  05:42:26  ldr
c Qcloud changes from LDR: new icefall Kessler option, Platt's (1996)
c optical properties and emissivity "fix" for low clouds.
c
c Revision 1.25  1997/11/24  23:25:27  ldr
c Indirect sulfates changes to give version used for GRL paper
c
c Revision 1.24  1997/07/24  05:59:17  ldr
c Mods to re-tune qcloud for 9L model from LDR.
c
c Revision 1.23  1996/10/24  01:02:46  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.22  1996/06/13  01:50:56  ldr
c Update qcloud to run 24E (works for 24L and 18L)
c
c Revision 1.21  1996/03/21  03:41:39  ldr
c Merge of TIE and LDR changes.
c
c Revision 1.20  1996/03/21  03:18:43  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.19.1.1  1996/03/21  03:30:58  ldr
c Reduce conv cloud and increase rcritl to 0.85 (for 24 levels).
c
c Revision 1.19  1996/01/09  06:18:16  ldr
c This is the version used for the long run H06 and write-up.
c
c Revision 1.18  1995/11/30  04:40:16  ldr
c Final updates to qcloud scheme to bring it to long run h63 used for
c BMRC workshop and write-up.
c
c Revision 1.17  1995/11/23  06:03:29  ldr
c Update qcloud to run h32 (used for BMRC & long run)
c
c Revision 1.16  1995/08/29  01:46:11  ldr
c Update qcloud to run g61.
c
c Revision 1.15  1995/08/14  07:32:40  ldr
c Updated qcloud to run g39, run for several years on kaos.
c
c Revision 1.14  1995/08/08  02:02:15  ldr
c Changes to qcloud (run g28) and tidied some diagnostics
c
c Revision 1.13  1995/06/30  02:44:40  ldr
c Changes to qcloud (run F88) and put qsat formula into ESTABL.f
c
c Revision 1.12  1995/05/02  06:10:36  ldr
c Changes to V4-6-16l from LDR - mainly affecting qcloud scheme.
c This version used for run F35.
c
c Revision 1.11  1995/03/14  01:21:44  ldr
c Changes to qcloud: incorporate hlfusion in snow that hits ground; use
c preset vice in slow icefall calc. - this is run E85, still with small
c fluxi array in icefall.f.
c
c Revision 1.10  1994/12/15  07:04:23  ldr
c Corrected Ecr to be Ec, consostent with V4-6-5l.
c
c Revision 1.9  1994/12/15  06:27:05  ldr
c Further mods to LDR's cloud scheme - this version is as used in E58 run
c and as described in seminar of 27/10/94.
c
c Revision 1.8  1994/12/07  23:53:36  mrd
c Explicitly declare all variables in include files.
c
c Revision 1.7  1994/09/12  12:50:06  ldr
c Changes to qcloud scheme of V4-6 to create run E32, except that E32 has
c enhanced collection of cloud water by falling ice.
c
c Revision 1.6  94/08/04  11:07:41  ldr
c A host of minor changes - this version used for run E03.
c 
c Revision 1.5  94/07/01  14:43:34  ldr
c New stuff for Manton-Cotton precipitation parameterization - these values
c used for run e58.
c 
c Revision 1.4  94/06/17  15:49:24  ldr
c Version used for run e37.
c 
c Revision 1.3  94/05/13  14:55:43  ldr
c Changes to introduce separate prognostic variable for cloud ice.
c 
c Revision 1.2  94/05/03  10:00:35  ldr
c Increase qcrl/qcrs back to 8/4e-4.
c 
c Revision 1.1  94/03/30  09:50:41  ldr
c Initial revision
c 

c Physical constants

      real rhow, rhoice, um, Dva, rKa
      parameter (rhow=1000.) !Density of water
      parameter (rhoice=917.)!Density of ice
      parameter (um=1.717e-5) !Dynamic viscosity of air (at 0 deg. C)
      parameter (Dva=2.21) !Diffusivity of qv in air (0 deg. and 1 Pa)
      parameter (rKa=2.4e-2) !Thermal conductivity of air (0 deg)


c Tunable parameters for qcloud scheme

      real ti, tice, aa, bb
      parameter (ti = -40.)  ! Min T for liquid water clouds in Celsius
      parameter (tice=273.15+ti)  !Convert ti to Kelvin
      parameter (aa=-2/ti**3, bb=3/ti**2) ! Coeffs for cubic interp of fracice

c The following are used in the Manton-Cotton rain parameterization

      real cdrops_nh, cdropl_nh, cdrops_sh, cdropl_sh, Ec, Aurate
      parameter (cdrops_nh=1.e8, cdropl_nh=3.e8) !Cloud droplet conc sea/land nh
!     parameter (cdrops_sh=1.e8, cdropl_sh=3.e8) !Cloud droplet conc sea/land sh
      parameter (cdrops_sh=.5e8, cdropl_sh=1.e8) !Cloud droplet conc sea/land sh
!     parameter (Rcm=7.5e-6) !Threshold cloud dropl R for coalescence to begin
!!    parameter (Rcm=10.e-6) !Threshold cloud dropl R for coalescence to begin
      parameter (Ec=0.55) !Mean collection efficiency for cloud drops
      parameter (Aurate=0.104*grav*Ec/um) !Part of rate constant

c Parameters related to snow

      real rhosno, Eac
      parameter (rhosno=100.) !Assumed density of snow in kg/m^3
      parameter (Eac=0.7)  !Mean collection efficiency of ql by snow

c Parameters related to rain

      real Nor, rk2, rho0, Ecol
      parameter (Nor=8.0e6) !Marshall-Palmer intercept parameter
      parameter (rk2=142.) !Fall speed of rain: V(D)=rk2*sqrt(D)*sqrt(rho0/rhoa)
      parameter (rho0=1.2)  !Standard air density
      parameter (Ecol=0.7)  !Mean collection efficiency of ql by rain

c Parameters related to diagnosed convective cloud

      real wlc,wls,ticon
      parameter (wlc=0.2e-3)  !LWC of deep conv cloud (kg/m**3)
      parameter (wls=0.35e-3) !LWC of shallow (non-preciptating) conv cloud
      parameter (ticon=238.15)!Temp at which conv cloud becomes ice

c Parameters related to cloud radiative properties

      real aice, bice
c     parameter (aice=0.032)!Constant in Platt optical depth for ice (SI units)
c     parameter (bice=0.333)!Constant in Platt optical depth for ice (SI units)
      parameter (aice=1.016)!Constant in Platt optical depth for ice (SI units)
      parameter (bice=0.68)!Constant in Platt optical depth for ice (SI units)
cldr  parameter (aice=0.23)!Constant in H-P optical depth for ice (SI units)
cldr  parameter (bice=0.54)!Constant in H-P optical depth for ice (SI units)
