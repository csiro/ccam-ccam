c Imported from CSIRO GCM
      real stefbo, erad, eradsq, cp, rdry, epsil, rvap, hlf, hl
     &, hls, hlcp, grav, hlfcp, hlscp, sq2, cappa, tomg, pi

      parameter (stefbo=5.67e-8) !Stefan-Boltzmann constant
      parameter (erad=6.37122e6, eradsq=erad*erad) !Radius of earth
      parameter (cp=1004.64) ! Specific heat of dry air at const P
      parameter (rdry=287.04) ! Specific gas const for dry air
      parameter (epsil=0.622) ! Ratio molec wt of H2O vapour to dry air
      parameter (rvap=461.) ! Gas constant for water vapour
      parameter (hlf=3.35e5) ! Latent heat of fusion (at 0 deg C)
      parameter (hl=2.5e6) !Latent heat of vaporization (at 0 deg. C)
      parameter (hls=hl+hlf) !  "      "   " sublimation
      parameter (hlcp=hl/cp, hlfcp=hlf/cp, hlscp=hlcp+hlfcp)
      parameter (grav=9.80616) ! Acceleration of gravity

      parameter (sq2=1.414213562373092) ! Square root of 2
      parameter (cappa=rdry/cp) 
      parameter (tomg=2*7.2921233e-5) ! 2*omega
      parameter (pi=3.14159265) !Good ol' pi

      real hcap50 ! Heat capacity of sea water * 50m (J/m**2/K)
      real dzmlo  ! Depth of mixed layer ocean points
      parameter (hcap50=2.095e8,dzmlo=100.0)
      real hcap   ! Heat capacity of 1m of sea water (J/m**2/K)
      parameter (hcap=hcap50/50.0)

      real tfrz
      parameter (tfrz=273.15)

