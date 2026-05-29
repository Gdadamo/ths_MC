module mod_setup
  use mod_params, only: dp, PI, TWOPI
  use mod_system, only: sys, chain_t, update_ghosts
  use mod_rng,    only: rng_uniform
  implicit none
  private

  public :: setup_system

contains

  subroutine setup_system(basename)
    character(len=*), intent(in) :: basename
    logical  :: exists
    integer  :: ic

    inquire(file=trim(basename)//'.RS', exist=exists)

    if (exists) then
      call read_restart(basename)
    else
      do ic = 1, sys%nc
        if (sys%chain(ic)%is_ring) then
          call init_ring(sys%chain(ic))
        else
          call init_linear(sys%chain(ic))
        end if
      end do
    end if

    ! Refresh ghost nodes for all chains
    do ic = 1, sys%nc
      call update_ghosts(sys%chain(ic))
    end do

    write(*,'(a)') ' Configuration ready.'
  end subroutine

  subroutine init_ring(ch)
    type(chain_t), intent(inout) :: ch
    integer  :: i
    real(dp) :: radius, theta, dtheta, zz

    dtheta = TWOPI / real(ch%np, dp)
    radius = sys%bondl / (2.0_dp * sin(PI / real(ch%np, dp)))

    call rng_uniform(zz)
    theta = TWOPI * zz

    do i = 1, ch%np
      ch%r(1, i) = radius * cos(theta + dtheta * real(i-1, dp))
      ch%r(2, i) = radius * sin(theta + dtheta * real(i-1, dp))
      ch%r(3, i) = 0.0_dp
    end do
  end subroutine

  subroutine init_linear(ch)
    type(chain_t), intent(inout) :: ch
    integer  :: i
    real(dp) :: dx, dz, sign_z

    dx     = sys%bondl * 0.5_dp
    dz     = sys%bondl * sqrt(3.0_dp)*0.5_dp
    sign_z = 1.0_dp

    ch%r(:, 1) = 0.0_dp

    do i = 2, ch%np
      ch%r(1, i) = ch%r(1, i-1) + dx
      ch%r(2, i) = 0.0_dp
      ch%r(3, i) = ch%r(3, i-1) + sign_z * dz
      sign_z = -sign_z
    end do

    ch%r(1, :) = ch%r(1, :) - ch%r(1, (ch%np+1)/2)
    ch%r(3, :) = ch%r(3, :) - ch%r(3, (ch%np+1)/2)
  end subroutine

  subroutine read_restart(basename)
    character(len=*), intent(in) :: basename
    integer  :: u_rs, ic, i, ierr
    real(dp) :: rx, ry, rz

    open(newunit=u_rs, file=trim(basename)//'.RS', status='old', &
         action='read', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) 'ERROR: cannot open ', trim(basename)//'.RS'; stop 1
    end if

    do ic = 1, sys%nc
      do i = 1, sys%np
        read(u_rs, *) rx, ry, rz
        sys%chain(ic)%r(1, i) = rx
        sys%chain(ic)%r(2, i) = ry
        sys%chain(ic)%r(3, i) = rz
      end do
    end do

    close(u_rs)
    write(*,'(a)') ' Restart loaded from '//trim(basename)//'.RS'
  end subroutine

end module mod_setup
