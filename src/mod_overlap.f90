module mod_overlap
  ! THS (Tangent Hard Sphere) overlap detection.
  ! Overlap criterion: |r_i - r_j|² < diam²
  ! Bonded pairs (topological distance = 1) are NOT checked:
  !   the user guarantees b >= d, so bonded monomers never overlap.
  use mod_params, only: dp
  use mod_system, only: sys, chain_t
  implicit none
  private

  public :: overlap_pair, check_chain_valid, check_moved_vs_fixed

contains

  ! Pure sphere-sphere overlap test. No sqrt needed.
  pure logical function overlap_pair(ri, rj, diam2)
    real(dp), intent(in) :: ri(3), rj(3), diam2
    real(dp) :: dx, dy, dz
    dx = ri(1) - rj(1)
    dy = ri(2) - rj(2)
    dz = ri(3) - rj(3)
    overlap_pair = (dx*dx + dy*dy + dz*dz) < diam2
  end function

  ! Full O(N²) check — used at startup/restart and in debug mode.
  ! Prints first 5 overlapping pairs with their distances.
  logical function check_chain_valid(ch)
    type(chain_t), intent(in) :: ch
    integer  :: i, j, topo, nfound
    real(dp) :: dist

    nfound = 0
    do i = 1, ch%np - 1
      do j = i + 1, ch%np
        topo = j - i
        if (ch%is_ring) topo = min(topo, ch%np - topo)
        if (topo <= 1) cycle

        if (overlap_pair(ch%r(:,i), ch%r(:,j), sys%diam2)) then
          nfound = nfound + 1
          dist = sqrt(sum((ch%r(:,i)-ch%r(:,j))**2))
          write(*,'(a,2i6,a,f10.6,a,i4)') &
            '  OVERLAP pair: (', i, j, ')  dist=', dist, '  topo=', topo
          if (nfound >= 5) then
            write(*,*) '  ... (more overlaps not shown)'
            check_chain_valid = .false.
            return
          end if
        end if
      end do
    end do

    check_chain_valid = (nfound == 0)
  end function

  ! Debug: brute-force check of moved arc vs complementary arc.
  ! Call after an accepted crank move to catch bugs the crank check missed.
  ! moved(1:n_move) = physical indices of rotated monomers.
  ! Returns .false. and prints details if any new overlap found.
  logical function check_moved_vs_fixed(ch, moved, n_move)
    type(chain_t), intent(in) :: ch
    integer,       intent(in) :: moved(:), n_move
    integer :: k, j, topo
    real(dp) :: dist

    check_moved_vs_fixed = .true.
    do k = 1, n_move
      do j = 1, ch%np
        ! skip if j is also in the moved arc
        if (any(moved(1:n_move) == j)) cycle
        topo = abs(moved(k) - j)
        if (ch%is_ring) topo = min(topo, ch%np - topo)
        if (topo <= 1) cycle

        if (overlap_pair(ch%r(:, moved(k)), ch%r(:, j), sys%diam2)) then
          dist = sqrt(sum((ch%r(:,moved(k)) - ch%r(:,j))**2))
          write(*,'(a,2i6,a,f10.6)') &
            '  BUG: missed overlap (', moved(k), j, ')  dist=', dist
          check_moved_vs_fixed = .false.
        end if
      end do
    end do
  end function

end module mod_overlap
