! rothc_step_f.f90
!
! Single-month RothC decomposition step for use with R's .Fortran() interface.
!
! Accepts a pre-computed combined rate modifier xi (= RM_Tmp * RM_Moist * RM_PC)
! rather than raw climate data. Climate modifier computation stays in R
! (compute_xi_rothc); Fortran handles only the decomposition arithmetic.
!
! Pool update (RothC manual Section 1.4):
!   pool1  = pool0 * exp(-xi * k * dt)
!   delta  = pool0 - pool1
!   Bio_in = 0.46/(x+1) * sum(delta)
!   Hum_in = 0.54/(x+1) * sum(delta)
!
! Clay texture factor (Section 1.7):
!   x = 1.67 * (1.85 + 1.60 * exp(-0.0786 * clay))
!
! Rate constants [yr^-1]: DPM=10.0, RPM=0.3, Bio=0.66, Hum=0.02
! Timestep: dt = 1/12 yr (monthly)
! IOM is passed through unchanged.

subroutine rothc_step_f(DPM, RPM, Bio, Hum, IOM, &
                        xi, clay, C_in_DPM, C_in_RPM)

  implicit none

  ! Pool states (in/out)
  real*8, intent(inout) :: DPM, RPM, Bio, Hum, IOM

  ! Combined rate modifier (precomputed in R)
  real*8, intent(in) :: xi

  ! Site property
  real*8, intent(in) :: clay

  ! Monthly litter inputs (already split by PL fractions)
  real*8, intent(in) :: C_in_DPM, C_in_RPM

  ! Fixed rate constants [yr^-1]
  real*8, parameter :: DPM_k = 10.0d0
  real*8, parameter :: RPM_k =  0.3d0
  real*8, parameter :: Bio_k =  0.66d0
  real*8, parameter :: Hum_k =  0.02d0
  real*8, parameter :: dt    =  1.0d0 / 12.0d0

  ! Local working variables
  real*8 :: DPM1, RPM1, Bio1, Hum1
  real*8 :: dDPM, dRPM, dBio, dHum
  real*8 :: x, f_bio, f_hum
  real*8 :: Bio_in, Hum_in

  ! Clay texture factor (Section 1.7, eq. 13)
  x     = 1.67d0 * (1.85d0 + 1.60d0 * exp(-0.0786d0 * clay))
  f_bio = 0.46d0 / (x + 1.0d0)
  f_hum = 0.54d0 / (x + 1.0d0)

  ! Decay remaining in each pool (Section 1.4)
  DPM1 = DPM * exp(-xi * DPM_k * dt)
  RPM1 = RPM * exp(-xi * RPM_k * dt)
  Bio1 = Bio * exp(-xi * Bio_k * dt)
  Hum1 = Hum * exp(-xi * Hum_k * dt)

  ! Amount decomposed from each pool
  dDPM = DPM - DPM1
  dRPM = RPM - RPM1
  dBio = Bio - Bio1
  dHum = Hum - Hum1

  ! Recycled fractions to Bio and Hum from all decomposed material
  Bio_in = f_bio * (dDPM + dRPM + dBio + dHum)
  Hum_in = f_hum * (dDPM + dRPM + dBio + dHum)

  ! Update pools: remaining + recycled + fresh litter input
  DPM = DPM1 + C_in_DPM
  RPM = RPM1 + C_in_RPM
  Bio = Bio1 + Bio_in
  Hum = Hum1 + Hum_in
  ! IOM unchanged

end subroutine rothc_step_f
