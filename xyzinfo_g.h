!     integer, parameter :: mydoubl3 = selected_real_kind (10,100)
!     real (kind=mydoubl3) , dimension(ifull_g) :: x_g, y_g, z_g
      real*8, dimension(ifull_g) :: x_g, y_g, z_g
      real, dimension(ifull_g) :: wts_g
      common /xyz_g/ x_g,y_g,z_g,wts_g ! rlat and rlong gone
