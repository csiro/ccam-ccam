MODULE math_constants
  USE define_dimensions, ONLY : i_d, r_1
  IMPLICIT NONE
  PRIVATE i_d, r_1
  REAL(r_1), PARAMETER :: pi_c = 3.1415927
  REAL(r_1), PARAMETER :: pi180 = pi_c / 180.0 ! radians / degree
  REAL(r_1), PARAMETER :: two_pi = 2.0 * pi_c
END MODULE math_constants
