! rothc_step_f.f90
!
! Single-month RothC decomposition step for use with R's .Fortran() interface.
!
! Accepts a pre-computed combined rate modifier xi (= RM_Tmp * RM_Moist * RM_PC)
! rather than raw climate data. Climate modifier computation stays in R
! (compute_xi_rothc); Fortran handles only the decomposition arithmetic.
!
! Rate constants k_DPM, k_RPM, k_Bio, k_Hum and clay-scaled partitioning
! fractions f_bio, f_hum are now passed as arguments rather than hardcoded,
! matching the parameterisation structure of the Yasso wrappers.
!
! The clay-scaled effective fractions f_bio and f_hum must be pre-computed
! in R before calling this subroutine:
!   x     = 1.67 * (1.85 + 1.60 * exp(-0.0786 * clay))
!   f_bio = f_bio_param / (x + 1)
!   f_hum = f_hum_param / (x + 1)
!
! Pool update (RothC manual Section 1.4):
!   pool1  = pool0 * exp(-xi * k * dt)
!   delta  = pool0 - pool1
!   Bio_in = f_bio * sum(delta)
!   Hum_in = f_hum * sum(delta)
!
! Timestep: dt = 1/12 yr (monthly)
! IOM is passed through unchanged.

subroutine rothc_step_f(DPM, RPM, Bio, Hum, IOM,    &
                        xi,                           &
                        k_DPM, k_RPM, k_Bio, k_Hum,  &
                        f_bio, f_hum,                 &
                        C_in_DPM, C_in_RPM)

  implicit none

  ! Pool states (in/out)
  real*8, intent(inout) :: DPM, RPM, Bio, Hum, IOM

  ! Combined rate modifier (precomputed in R: RM_Tmp * RM_Moist * RM_PC)
  real*8, intent(in) :: xi

  ! Decomposition rate constants [yr^-1] -- calibration parameters
  real*8, intent(in) :: k_DPM, k_RPM, k_Bio, k_Hum

  ! Clay-scaled partitioning fractions -- pre-computed in R from calibration
  ! parameters f_bio, f_hum and the clay texture factor x:
  !   f_bio_eff = f_bio_param / (x + 1)
  !   f_hum_eff = f_hum_param / (x + 1)
  real*8, intent(in) :: f_bio, f_hum

  ! Monthly litter inputs [t C ha^-1 month^-1]
  real*8, intent(in) :: C_in_DPM, C_in_RPM

  ! Fixed timestep
  real*8, parameter :: dt = 1.0d0 / 12.0d0

  ! Local working variables
  real*8 :: DPM1, RPM1, Bio1, Hum1
  real*8 :: dDPM, dRPM, dBio, dHum
  real*8 :: Bio_in, Hum_in

  ! Decay remaining in each pool after one timestep (Section 1.4)
  DPM1 = DPM * exp(-xi * k_DPM * dt)
  RPM1 = RPM * exp(-xi * k_RPM * dt)
  Bio1 = Bio * exp(-xi * k_Bio * dt)
  Hum1 = Hum * exp(-xi * k_Hum * dt)

  ! Amount decomposed from each pool
  dDPM = DPM - DPM1
  dRPM = RPM - RPM1
  dBio = Bio - Bio1
  dHum = Hum - Hum1

  ! Recycled fractions routed to Bio and Hum from all decomposed material
  Bio_in = f_bio * (dDPM + dRPM + dBio + dHum)
  Hum_in = f_hum * (dDPM + dRPM + dBio + dHum)

  ! Update pools: remaining + recycled + fresh litter input
  DPM = DPM1 + C_in_DPM
  RPM = RPM1 + C_in_RPM
  Bio = Bio1 + Bio_in
  Hum = Hum1 + Hum_in
  ! IOM unchanged

end subroutine rothc_step_f
