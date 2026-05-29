module mod_params
  implicit none

  integer,  parameter :: dp     = kind(1.0d0)
  integer,  parameter :: sp     = kind(1.0e0)
  integer,  parameter :: MAXPC  = 2

  real(dp), parameter :: PI     = 3.1415926535897932384626433832795_dp
  real(dp), parameter :: TWOPI  = 2.0_dp * PI
  real(dp), parameter :: FPIOT  = 4.1887902047863909846168578443727_dp  ! 4π/3

end module mod_params
