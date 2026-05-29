module mod_rng
  ! Marsaglia-Zaman-Tsang pseudo-random number generator.
  ! Ported from marsaglia.f: COMMON replaced by module-level save variables.
  ! Period ~ 2^144. Uniform output in (0,1).
  use mod_params, only: dp
  implicit none
  private

  public :: rng_init, rng_save, rng_load, rng_uniform, rng_unit_vector

  real(dp), save :: u(97), c, cd, cm
  integer,  save :: i_idx, j_idx

contains

  subroutine rng_init(ij, kl)
    ! Initialise from integer seeds. Defaults match marsaglia.f: ij=1802, kl=9373.
    integer, intent(in), optional :: ij, kl
    integer :: ij_, kl_, ii, jj, k, m, n, i_, j_
    real(dp) :: s, t

    ij_ = 1802; if (present(ij)) ij_ = ij
    kl_ = 9373; if (present(kl)) kl_ = kl

    i_ = mod(ij_/177, 177) + 2
    j_ = mod(ij_,    177) + 2
    k  = mod(kl_/169, 178) + 1
    m  = mod(kl_,    169)

    do ii = 1, 97
      s = 0.0_dp
      t = 0.5_dp
      do jj = 1, 24
        n = mod(mod(i_*j_, 179)*k, 179)
        i_ = j_; j_ = k; k = n
        m = mod(53*m + 1, 169)
        if (mod(m*n, 64) >= 32) s = s + t
        t = 0.5_dp * t
      end do
      u(ii) = s
    end do

    c      = 362436.0_dp  / 16777216.0_dp
    cd     = 7654321.0_dp / 16777216.0_dp
    cm     = 16777213.0_dp/ 16777216.0_dp
    i_idx  = 97
    j_idx  = 33

    write(*,'(a,2i8)') ' RNG initialised (seeds): ', ij_, kl_
  end subroutine

  subroutine rng_save(filename)
    character(len=*), intent(in) :: filename
    integer :: u_unit
    open(newunit=u_unit, file=filename, status='unknown', form='unformatted')
    write(u_unit) u, c, cd, cm, i_idx, j_idx
    close(u_unit)
  end subroutine

  subroutine rng_load(filename)
    character(len=*), intent(in) :: filename
    integer :: u_unit
    open(newunit=u_unit, file=filename, status='old', form='unformatted')
    read(u_unit) u, c, cd, cm, i_idx, j_idx
    close(u_unit)
    write(*,'(a,a)') ' RNG state loaded from ', trim(filename)
  end subroutine

  subroutine rng_uniform(x)
    real(dp), intent(out) :: x
    real(dp) :: uni
    uni = u(i_idx) - u(j_idx)
    if (uni < 0.0_dp) uni = uni + 1.0_dp
    u(i_idx) = uni
    i_idx = i_idx - 1; if (i_idx == 0) i_idx = 97
    j_idx = j_idx - 1; if (j_idx == 0) j_idx = 97
    c = c - cd;        if (c < 0.0_dp)  c = c + cm
    uni = uni - c
    if (uni < 0.0_dp) uni = uni + 1.0_dp
    x = uni
  end subroutine

  ! Uniformly distributed random unit vector (Marsaglia rejection method).
  subroutine rng_unit_vector(a)
    real(dp), intent(out) :: a(3)
    real(dp) :: rt, zz
10  call rng_uniform(zz); a(1) = zz - 0.5_dp
    call rng_uniform(zz); a(2) = zz - 0.5_dp
    rt = a(1)*a(1) + a(2)*a(2)
    if (rt > 0.125_dp) goto 10
    call rng_uniform(zz); a(3) = zz - 0.5_dp
    rt = rt + a(3)*a(3)
    if (rt > 0.125_dp) goto 10
    rt = 1.0_dp / sqrt(rt)
    a  = a * rt
  end subroutine

end module mod_rng
