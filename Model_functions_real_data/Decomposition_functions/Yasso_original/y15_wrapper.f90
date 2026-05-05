! y15_wrapper.f90 -- plain subroutine wrapper for R
SUBROUTINE mod5c_r(theta, time, climate, init, b, d, leac, xt, steady)
  USE yasso
  IMPLICIT NONE
  REAL, DIMENSION(35), INTENT(IN)  :: theta
  REAL, INTENT(IN)                 :: time, d, leac
  REAL, DIMENSION(3),  INTENT(IN)  :: climate
  REAL, DIMENSION(5),  INTENT(IN)  :: init, b
  REAL, DIMENSION(5),  INTENT(OUT) :: xt
  LOGICAL,             INTENT(IN)  :: steady
  CALL mod5c(theta, time, climate, init, b, d, leac, xt, steady)
END SUBROUTINE mod5c_r