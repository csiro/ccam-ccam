      integer, parameter :: nstnmax=47
      integer :: nstn
      integer, dimension(nstnmax) :: istn,jstn,iunp
      real, dimension(nstnmax)    :: slat,slon,zstn
      logical, dimension(nstnmax) :: mystn
      character(len=3), dimension(nstnmax) :: name_stn
      common /trcom2/nstn,istn,jstn,iunp,slat,slon,zstn,                  &
     &               mystn,name_stn
