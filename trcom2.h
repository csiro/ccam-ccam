      integer, parameter :: nstnmax=20, nstn2=3
      integer :: mstn,nstn,iaustw,iauste,iaustn,iausts
      integer, dimension(nstnmax) :: istn,jstn,iunp
      integer, dimension(nstn2)   :: istn2,jstn2,iunp2
      real, dimension(nstnmax)    :: slat,slon,zstn
      real, dimension(nstn2)      :: slat2,slon2
      character(len=3) name_stn
      common /trcom2/ mstn,nstn,slat,slon,istn,jstn,iunp,zstn,
     &                slat2,slon2,istn2,jstn2,iunp2,
     &                iaustw,iauste,iaustn,iausts,name_stn
