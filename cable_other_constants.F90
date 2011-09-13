MODULE other_constants
  USE define_dimensions, ONLY : i_d, r_1, nrb
  REAL(r_1), PARAMETER, DIMENSION(nrb) :: gauss_w=(/0.308,0.514,0.178/) ! Gaussian integ. weights
!  REAL(r_1), PARAMETER, DIMENSION(nrb) :: refl = (/ 0.1, 0.425, 0.05 /) ! leaf reflectance
!  REAL(r_1), PARAMETER, DIMENSION(nrb) :: taul = (/ 0.1, 0.425, 0.05/)  ! leaf transmittance
END MODULE other_constants
