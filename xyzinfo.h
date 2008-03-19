!     integer, parameter :: mydoubl1 = selected_real_kind (10,100)
!     real (kind=mydoubl1) , dimension(ifull) :: x, y, z
      real*8 x(ifull),y(ifull),z(ifull)
      real wts(ifull)
      common/xyz/x,y,z,wts ! rlat and rlong gone
