module mod_crank
  ! Crank-shaft MC move for ring polymers (THS model).
  !
  ! Algorithm:
  !   1. Pick FIRST, LAST in 1..NP with topological distance >= 2.
  !   2. 50/50 coin flip: swap(FIRST,LAST) to choose which arc to rotate.
  !   3. Normalise so LAST > FIRST (add NP if needed).
  !   4. Rotate interior monomers FIRST+1..LAST-1 around axis r(IFI)→r(ILA),
  !      checking overlap incrementally against complementary arc LAST..FIRST+NP.
  !   5. On overlap at monomer AUX_IP: restore only first+1..aux_ip.
  !
  ! Bonded pairs (topo dist = 1) are skipped in the check loop.
  ! Smart skip: KK = max(1, floor((dist - diam) / bondl))
  use mod_params,  only: dp, TWOPI
  use mod_system,  only: sys, chain_t, update_ghosts
  use mod_overlap, only: overlap_pair, check_moved_vs_fixed
  use mod_rng,     only: rng_uniform
  implicit none
  private

  public :: crank_move, crank_init, crank_free

  ! Set DEBUG_CRANK=.true. to validate every accepted move (slow — debug only).
  logical, parameter :: DEBUG_CRANK = .false.

  ! Module-level workspace: allocated once, avoids per-call allocation in hot loop.
  real(dp), allocatable, save :: rold(:,:)   ! (3, NP)

contains

  subroutine crank_init(np)
    integer, intent(in) :: np
    if (.not. allocated(rold)) allocate(rold(3, np))
  end subroutine

  subroutine crank_free()
    if (allocated(rold)) deallocate(rold)
  end subroutine

  subroutine crank_move(ch, acc)
    type(chain_t), intent(inout) :: ch
    real(dp),      intent(out)   :: acc

    integer  :: first, last, ifi, ila, ip, jp, ii, jj, kk, np_
    integer  :: n_move, topo_comp, tmp_i, aux_ip
    real(dp) :: zz, rr, theta
    real(dp) :: axis(3), rm(3,3), bvec(3), ri(3), rj(3)
    logical  :: ov

    acc    = 0.0_dp
    ov     = .false.
    np_    = ch%np
    aux_ip = 0   ! last ip that has been saved + rotated (0 = nothing moved)

    ! --- 1. Pick FIRST and LAST with topological distance >= 2 ---
100 call rng_uniform(zz); first = int(np_ * zz) + 1
    call rng_uniform(zz); last  = int(np_ * zz) + 1
    if (min(abs(last - first), np_ - abs(last - first)) < 2) goto 100

    ! --- 2. 50/50: swap to choose which arc we rotate ---
    call rng_uniform(zz)
    if (zz > 0.5_dp) then; tmp_i = first; first = last; last = tmp_i; end if

    ! --- 3. Normalise: ensure LAST > FIRST ---
    if (last <= first) last = last + np_

    n_move    = last - first - 1
    topo_comp = np_ - (last - first) - 1

    if (n_move < 1) goto 100

    ! Physical indices of pivot endpoints (1..NP)
    ifi = mod(first - 1, np_) + 1
    ila = mod(last  - 1, np_) + 1

    ! --- 4. Rotation matrix around axis r(ifi)→r(ila) ---
    axis = ch%r(:, ila) - ch%r(:, ifi)
    rr   = sum(axis * axis)
    if (rr < 1.0e-12_dp) goto 100
    axis = axis / sqrt(rr)

    call rng_uniform(zz)
    theta = TWOPI * zz
    call rot_matrix(axis, theta, rm)

    ! --- 5. Rotate monomers one by one, check each vs complementary arc ---
    do ip = first + 1, last - 1
      ii = mod(ip - 1, np_) + 1

      ! Save BEFORE rotating (needed for rollback)
      rold(:, ii) = ch%r(:, ii)

      bvec = ch%r(:, ii) - ch%r(:, ifi)
      ch%r(:, ii) = ch%r(:, ifi) + matmul(rm, bvec)
      aux_ip = ip   ! this monomer is now in the rotated state

      ri = ch%r(:, ii)

      ! Check against complementary arc: LAST .. FIRST+NP (inclusive)
      ! Includes both shared endpoints IFI (=FIRST) and ILA (=LAST).
      jp = last
      do while (jp <= first + np_)
        jj = mod(jp - 1, np_) + 1

        ! Skip bonded pair (topo dist 1): only boundary bonds trigger this
        if (abs(ii - jj) == 1 .or. abs(ii - jj) == np_ - 1) then
          jp = jp + 1
          cycle
        end if

        rj = ch%r(:, jj)
        rr = sum((ri - rj) * (ri - rj))

        if (rr < sys%diam2) then
          ov = .true.
          goto 555
        end if

        ! Smart skip along complementary arc
        kk = max(1, int((sqrt(rr) - sys%diam) / sys%bondl))
        jp = jp + kk
      end do
    end do

555 continue

    if (ov) then
      ! Restore only the monomers that were actually rotated (first+1 .. aux_ip)
      do ip = first + 1, aux_ip
        ii = mod(ip - 1, np_) + 1
        ch%r(:, ii) = rold(:, ii)
      end do
    else
      acc = 1.0_dp
      call update_ghosts(ch)

      ! Debug: brute-force check that the accepted move introduced no overlaps.
      ! Catches any pair the crank inner loop may have missed.
      if (DEBUG_CRANK) then
        block
          integer :: moved_list(n_move), k
          do k = 1, n_move
            moved_list(k) = mod(first + k, np_) + 1
          end do
          if (.not. check_moved_vs_fixed(ch, moved_list, n_move)) then
            write(*,'(a,2i6)') ' BUG in crank_move! first=', first, last
          end if
        end block
      end if
    end if

  end subroutine

  ! Rodrigues rotation matrix around unit axis a by angle theta.
  pure subroutine rot_matrix(a, theta, r)
    real(dp), intent(in)  :: a(3), theta
    real(dp), intent(out) :: r(3,3)
    real(dp) :: c, s, t
    c = cos(theta); s = sin(theta); t = 1.0_dp - c
    r(1,1) = c + a(1)*a(1)*t;       r(1,2) = a(1)*a(2)*t - a(3)*s; r(1,3) = a(1)*a(3)*t + a(2)*s
    r(2,1) = a(2)*a(1)*t + a(3)*s;  r(2,2) = c + a(2)*a(2)*t;      r(2,3) = a(2)*a(3)*t - a(1)*s
    r(3,1) = a(3)*a(1)*t - a(2)*s;  r(3,2) = a(3)*a(2)*t + a(1)*s; r(3,3) = c + a(3)*a(3)*t
  end subroutine

end module mod_crank
