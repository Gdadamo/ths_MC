module mod_system
  use mod_params, only: dp, MAXPC
  implicit none

  ! Single polymer chain.
  ! Positions: r(3, 0:np+1)
  !   r(:,1..np)   physical monomers
  !   r(:,0)       ghost = r(:,np)   \  updated after every accepted move
  !   r(:,np+1)    ghost = r(:,1)   /  (ring only — eliminates MODULO in hot loop)
  type :: chain_t
    integer               :: np
    logical               :: is_ring
    real(dp), allocatable :: r(:,:)   ! (3, 0:np+1)
  end type

  ! Full system state, shared by all routines via use association.
  type :: system_t
    integer        :: nc              ! 1 or 2
    integer        :: np              ! monomers per chain (global)
    type(chain_t)  :: chain(MAXPC)
    real(dp)       :: diam            ! sphere diameter d
    real(dp)       :: diam2           ! d²
    real(dp)       :: bondl           ! bond length b
  end type

  ! Global singleton — every module uses this via:
  !   use mod_system, only: sys
  type(system_t), target :: sys

contains

  subroutine alloc_chain(ch, np, is_ring)
    type(chain_t), intent(inout) :: ch
    integer,       intent(in)    :: np
    logical,       intent(in)    :: is_ring
    ch%np      = np
    ch%is_ring = is_ring
    allocate(ch%r(3, 0:np+1))
    ch%r = 0.0_dp
  end subroutine

  subroutine dealloc_chain(ch)
    type(chain_t), intent(inout) :: ch
    if (allocated(ch%r)) deallocate(ch%r)
  end subroutine

  ! Refresh ghost nodes of a ring chain after a move.
  subroutine update_ghosts(ch)
    type(chain_t), intent(inout) :: ch
    if (.not. ch%is_ring) return
    ch%r(:, 0)       = ch%r(:, ch%np)
    ch%r(:, ch%np+1) = ch%r(:, 1)
  end subroutine

end module mod_system
