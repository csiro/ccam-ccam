module utilities
   implicit none
contains
   function calc_rotpole(rlong0,rlat0) result (rotpole)
!     Coefficients for rotation of pole
!     rotpole(1,) is x-axis of rotated coords in terms of orig Cartesian
!     rotpole(2,) is y-axis of rotated coords in terms of orig Cartesian
!     rotpole(3,) is z-axis of rotated coords in terms of orig Cartesian
      include 'const_phys.h'
      real, intent(in) :: rlong0, rlat0
      real, dimension(3,3) :: rotpole
      real ::  coslong, sinlong, coslat, sinlat
      coslong=cos(rlong0*pi/180.)
      sinlong=sin(rlong0*pi/180.)
      coslat=cos(rlat0*pi/180.)
      sinlat=sin(rlat0*pi/180.)
      rotpole(1,1)=coslong*sinlat
      rotpole(1,2)=-sinlong
      rotpole(1,3)=coslong*coslat
      rotpole(2,1)=sinlong*sinlat
      rotpole(2,2)=coslong
      rotpole(2,3)=sinlong*coslat
      rotpole(3,1)=-coslat
      rotpole(3,2)=0.
      rotpole(3,3)=sinlat
   end function calc_rotpole
end module utilities

