module mod_pivot
  ! Pivot MC move for linear polymers (THS model).
  !
  ! Algorithm:
  !   1. Pick a random pivot IPIV (not endpoint: 2..NP-1).
  !   2. 50/50: rotate LEFT arm (1..IPIV-1) or RIGHT arm (IPIV+1..NP).
  !   3. Apply a uniformly random rotation around IPIV.
  !   4. Check moved arm against fixed arm incrementally (smart skip).
  !   5. On overlap at aux_ip: restore only istart..aux_ip.
  !
  ! Bonded pair (moved monomer adjacent to IPIV) is skipped in check.
  use mod_params,  only: dp, TWOPI
  use mod_system,  only: sys, chain_t
  use mod_rng,     only: rng_uniform, rng_unit_vector
  implicit none
  private

  public :: pivot_move, pivot_init, pivot_free

  real(dp), allocatable, save :: rold_piv(:,:)   ! (3, NP)

contains

  subroutine pivot_init(np)
    integer, intent(in) :: np
    if (.not. allocated(rold_piv)) allocate(rold_piv(3, np))
  end subroutine

  subroutine pivot_free()
    if (allocated(rold_piv)) deallocate(rold_piv)
  end subroutine

  subroutine pivot_move(ch, acc)
    type(chain_t), intent(inout) :: ch
    real(dp),      intent(out)   :: acc

    integer  :: ipiv, istart, iend, istep, ip, jp, ii, jj, kk, np_, aux_ip
    real(dp) :: zz, rr
    real(dp) :: rm(3,3), bvec(3), rpiv(3), ri(3), rj(3)
    logical  :: ov

    acc    = 0.0_dp
    ov     = .false.
    np_    = ch%np
    aux_ip = 0

    ! --- 1. Pick pivot (not endpoint: ipiv in 2..NP-1) ---
100 call rng_uniform(zz)
    ipiv = int((np_ - 2) * zz) + 2
    if (ipiv < 2 .or. ipiv > np_ - 1) goto 100

    rpiv = ch%r(:, ipiv)

    ! --- 2. Choose arm: LEFT (1..ipiv-1) or RIGHT (ipiv+1..NP) ---
    call rng_uniform(zz)
    if (zz <= 0.5_dp) then
      istart = 1;      iend = ipiv - 1; istep =  1
    else
      istart = np_;    iend = ipiv + 1; istep = -1
    end if

    if (abs(iend - istart) < 0) goto 100   ! arm too short (NP=2 edge case)

    ! --- 3. Random rotation matrix (uniformly random axis + angle) ---
    call random_rot_matrix(rm)

    ! --- 4. Rotate arm one monomer at a time, check overlap ---
    aux_ip = istart - istep   ! nothing moved yet

    do ip = istart, iend, istep
      ii = ip

      rold_piv(:, ii) = ch%r(:, ii)
      bvec = ch%r(:, ii) - rpiv
      ch%r(:, ii) = rpiv + matmul(rm, bvec)
      aux_ip = ip   ! ii is now in rotated state

      ri = ch%r(:, ii)

      ! Check ii against the fixed arm (on the opposite side of ipiv)
      if (istep == 1) then
        ! moved arm: 1..ii (growing rightward)
        ! fixed arm: ipiv+1..NP (traverse forward with smart skip)
        jp = ipiv + 1
        do while (jp <= np_)
          jj = jp
          if (abs(ii - jj) <= 1) then; jp = jp + 1; cycle; end if
          rj = ch%r(:, jj)
          rr = sum((ri - rj) * (ri - rj))
          if (rr < sys%diam2) then; ov = .true.; goto 555; end if
          kk = max(1, int((sqrt(rr) - sys%diam) / sys%bondl))
          jp = jp + kk
        end do
      else
        ! moved arm: NP..ii (growing leftward)
        ! fixed arm: 1..ipiv-1 (traverse backward with smart skip)
        jp = ipiv - 1
        do while (jp >= 1)
          jj = jp
          if (abs(ii - jj) <= 1) then; jp = jp - 1; cycle; end if
          rj = ch%r(:, jj)
          rr = sum((ri - rj) * (ri - rj))
          if (rr < sys%diam2) then; ov = .true.; goto 555; end if
          kk = max(1, int((sqrt(rr) - sys%diam) / sys%bondl))
          jp = jp - kk
        end do
      end if
    end do

555 continue

    if (ov) then
      ! Restore only monomers that were actually saved + rotated: istart..aux_ip
      do ip = istart, aux_ip, istep
        ch%r(:, ip) = rold_piv(:, ip)
      end do
    else
      acc = 1.0_dp
    end if

  end subroutine

  ! Uniformly random rotation matrix: random unit axis + random angle in [0,2π].
  subroutine random_rot_matrix(r)
    real(dp), intent(out) :: r(3,3)
    real(dp) :: a(3), theta, zz, c, s, t
    call rng_unit_vector(a)
    call rng_uniform(zz)
    theta = TWOPI * zz
    c = cos(theta); s = sin(theta); t = 1.0_dp - c
    r(1,1) = c + a(1)*a(1)*t;       r(1,2) = a(1)*a(2)*t - a(3)*s; r(1,3) = a(1)*a(3)*t + a(2)*s
    r(2,1) = a(2)*a(1)*t + a(3)*s;  r(2,2) = c + a(2)*a(2)*t;      r(2,3) = a(2)*a(3)*t - a(1)*s
    r(3,1) = a(3)*a(1)*t - a(2)*s;  r(3,2) = a(3)*a(2)*t + a(1)*s; r(3,3) = c + a(3)*a(3)*t
  end subroutine

end module mod_pivot
