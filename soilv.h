      !integer, parameter :: mxst=13       ! max_no_of_soil_types ! MJT cable (see newmpar.h)
      real, dimension(0:mxst) :: swilt, ssat, sfc
      real, dimension(mxst)   :: bch, cnsd, css, hsbh, hyds, rhos, sucs,  &
     &                           clay, sand, silt
      real, dimension(mxvt)   :: rs20 ! MJT cable
      integer, dimension(mxst) :: i2bp3, ibp2
      real, dimension(44) :: rlaim44,rlais44,scveg44,rsunc44,slveg44
      real :: froot(5), zse(ms)
      common/soilpr/swilt,ssat,sfc,bch,cnsd,css,hsbh,hyds,i2bp3,ibp2,     &
     &              rhos,sucs,clay,sand,silt,rlaim44,rlais44,scveg44,     &
     &              rsunc44,slveg44,froot,zse,rs20 ! MJT cable

