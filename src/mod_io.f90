module mod_io
  ! I/O management for THS MC simulation.
  !
  ! Files managed:
  !   .OUT  — human-readable log (header written by mod_input, appended here)
  !   .DAT  — per-block observables (text, one line per block)
  !   .traj — binary trajectory (float32 coordinates, sampled every KSAMPLE blocks)
  !   .RS   — restart file (text, last configuration)
  !   .SEED — RNG state (binary, written at end)
  !
  ! All units are opened once in io_open, closed in io_close.
  use mod_params,      only: dp, sp
  use mod_system,      only: sys
  use mod_observables, only: obs_t
  use mod_rng,         only: rng_save
  implicit none
  private

  public :: io_open, io_close, io_write_block, io_write_traj, io_write_restart

  ! File unit numbers (module-level, persistent)
  integer, save :: u_log  = -1
  integer, save :: u_dat  = -1
  integer, save :: u_traj = -1

  ! Trajectory magic header string
  character(len=6), parameter :: TRAJ_MAGIC = 'THSMC1'

contains

  ! ----------------------------------------------------------------
  ! Open all output files, write headers.
  ! Call after read_input and setup_system.
  ! ----------------------------------------------------------------
  subroutine io_open(basename, ntrials, ksample)
    character(len=*), intent(in) :: basename
    integer,          intent(in) :: ntrials, ksample

    integer  :: ic
    integer(kind=4) :: i4
    real(kind=8)    :: r8

    ! --- .OUT log (append to header written by mod_input) ---
    open(newunit=u_log, file=trim(basename)//'.OUT', status='old', &
         position='append', action='write')

    ! --- .DAT observables ---
    open(newunit=u_dat, file=trim(basename)//'.DAT', status='replace', action='write')
    write(u_dat, '(a)') 'IB Rg2_1 Rg2_2 Ree2_1 Ree2_2 B2 acc_1 acc_2'

    ! --- .traj binary ---
    open(newunit=u_traj, file=trim(basename)//'.traj', status='replace', &
         action='write', form='unformatted', access='stream')

    ! Write header
    write(u_traj) TRAJ_MAGIC
    i4 = int(sys%nc, 4);   write(u_traj) i4
    i4 = int(sys%np, 4);   write(u_traj) i4
    r8 = sys%diam;          write(u_traj) r8
    r8 = sys%bondl;         write(u_traj) r8
    do ic = 1, sys%nc
      if (sys%chain(ic)%is_ring) then; i4 = 1; else; i4 = 0; end if
      write(u_traj) i4
    end do
    i4 = int(ksample, 4);  write(u_traj) i4

  end subroutine

  ! ----------------------------------------------------------------
  ! Write one line to .DAT and one line to .OUT (acceptance ratios).
  ! ----------------------------------------------------------------
  subroutine io_write_block(ib, obs, acc)
    integer,     intent(in) :: ib
    type(obs_t), intent(in) :: obs
    real(dp),    intent(in) :: acc(2)   ! acceptance ratios per chain

    ! .DAT: one row per block
    write(u_dat, '(i8, 6es20.10, 2f10.6)') &
      ib, obs%rg2(1), obs%rg2(2), obs%ree2(1), obs%ree2(2), obs%b2, acc(1), acc(2)
    flush(u_dat)

    ! .OUT: block summary
    write(u_log, '(a,i8,a,2f10.6)') &
      'Block ', ib, '  acc= ', acc(1), acc(2)
    flush(u_log)

  end subroutine

  ! ----------------------------------------------------------------
  ! Write a trajectory frame to .traj (float32 coordinates, CM-centred).
  ! Called only every KSAMPLE blocks.
  ! ----------------------------------------------------------------
  subroutine io_write_traj(ib)
    integer, intent(in) :: ib

    integer  :: ic, i
    integer(kind=4)  :: i4
    real(kind=4) :: x4

    i4 = int(ib, 4)
    write(u_traj) i4

    do ic = 1, sys%nc
      do i = 1, sys%np
        x4 = real(sys%chain(ic)%r(1, i), 4); write(u_traj) x4
        x4 = real(sys%chain(ic)%r(2, i), 4); write(u_traj) x4
        x4 = real(sys%chain(ic)%r(3, i), 4); write(u_traj) x4
      end do
    end do

  end subroutine

  ! ----------------------------------------------------------------
  ! Write configuration to .RS and RNG state to .SEED (both every block).
  ! Keeping them in sync guarantees a crash-safe restart: the RNG stream
  ! and the configuration always correspond to the same block.
  ! ----------------------------------------------------------------
  subroutine io_write_restart(basename)
    character(len=*), intent(in) :: basename

    integer :: u_rs, ic, i

    open(newunit=u_rs, file=trim(basename)//'.RS', status='replace', action='write')
    do ic = 1, sys%nc
      do i = 1, sys%np
        write(u_rs, '(3es22.14)') sys%chain(ic)%r(1,i), &
                                    sys%chain(ic)%r(2,i), &
                                    sys%chain(ic)%r(3,i)
      end do
    end do
    close(u_rs)

    call rng_save(trim(basename)//'.SEED')

  end subroutine

  ! ----------------------------------------------------------------
  ! Close all output files. RNG state is already in .SEED (written
  ! by io_write_restart every block), no need to save again here.
  ! ----------------------------------------------------------------
  subroutine io_close(basename)
    character(len=*), intent(in) :: basename

    if (u_log  > 0) close(u_log)
    if (u_dat  > 0) close(u_dat)
    if (u_traj > 0) close(u_traj)

  end subroutine

end module mod_io
