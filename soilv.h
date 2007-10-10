! eak version 16/03/06
!      integer, parameter :: mxst=16       ! max_no_of_soil_types  ! BP changed 13 to 16 due to IGBP types (OCT 2007)
      real, dimension(0:mxst) :: swilt, ssat, sfc
      real, dimension(mxst)   :: bch, cnsd, css, hsbh, hyds, rhos, sucs
      real, dimension(mxst)   :: clay, sand, silt, c3
      integer, dimension(mxst) :: i2bp3, ibp2
      real, dimension(mxst)   :: rs20
      real, dimension(ms) :: zse
      real, dimension(ms+1) :: zshh
      real, dimension(mxst,ms) :: froot
      real, dimension(44) :: rlaim44,rlais44,scveg44,rsunc44,slveg44

      common/soilpr/swilt,ssat,sfc,bch,cnsd,css,hsbh,hyds
      common/soilpr/i2bp3,ibp2,rhos,sucs,clay,sand,silt,c3
      common/soilpr/rs20
      common/soilpr/froot,zse,zshh
      common/soilpr/rlaim44,rlais44,scveg44,rsunc44,slveg44

