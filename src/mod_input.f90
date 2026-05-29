module mod_input
  use mod_params, only: dp
  use mod_system, only: sys, alloc_chain
  use mod_rng,    only: rng_init, rng_load
  implicit none
  private

  public :: mc_params_t, mc_params, read_input, rng_init_default, rng_load_seed

  type :: mc_params_t
    integer  :: nblock
    integer  :: nstep
    integer  :: ksample
    integer  :: ntrials
    integer  :: iseed
  end type

  type(mc_params_t) :: mc_params

contains

  subroutine read_input(basename)
    character(len=*), intent(in) :: basename
    integer  :: u_in, u_log, ic
    real(dp) :: diam
    character(len=16) :: ctype(2)

    open(newunit=u_in, file=trim(basename)//'.IN', status='old', &
         action='read', iostat=ic)
    if (ic /= 0) then
      write(*,*) 'ERROR: cannot open ', trim(basename)//'.IN'
      stop 1
    end if

    read(u_in, *) sys%nc
    read(u_in, *) sys%np

    if (sys%nc < 1 .or. sys%nc > 2) then
      write(*,*) 'ERROR: NC must be 1 or 2, got ', sys%nc; stop 1
    end if

    do ic = 1, sys%nc
      read(u_in, *) ctype(ic)
      ctype(ic) = trim(adjustl(ctype(ic)))
      if (ctype(ic) /= 'ring' .and. ctype(ic) /= 'linear') then
        write(*,*) 'ERROR: chain type must be ring or linear, got: ', trim(ctype(ic))
        stop 1
      end if
    end do
    if (sys%nc == 1) ctype(2) = 'linear'

    read(u_in, *) sys%diam
    read(u_in, *) sys%bondl

    if (sys%diam <= 0.0_dp .or. sys%bondl <= 0.0_dp) then
      write(*,*) 'ERROR: DIAM and BONDL must be positive'; stop 1
    end if
    if (sys%diam > sys%bondl) then
      write(*,*) 'WARNING: DIAM > BONDL — bonded monomers will overlap!'
    end if

    sys%diam2 = sys%diam * sys%diam

    read(u_in, *) mc_params%nblock
    read(u_in, *) mc_params%nstep
    read(u_in, *) mc_params%ksample
    read(u_in, *) mc_params%ntrials
    read(u_in, *) mc_params%iseed

    close(u_in)

    do ic = 1, sys%nc
      call alloc_chain(sys%chain(ic), sys%np, ctype(ic) == 'ring')
    end do

    open(newunit=u_log, file=trim(basename)//'.OUT', status='replace', action='write')
    write(u_log, '(a,i4)')      'NC      = ', sys%nc
    write(u_log, '(a,i8)')      'NP      = ', sys%np
    do ic = 1, sys%nc
      if (sys%chain(ic)%is_ring) then
        write(u_log, '(a,i2,a)') 'chain ', ic, ' = ring'
      else
        write(u_log, '(a,i2,a)') 'chain ', ic, ' = linear'
      end if
    end do
    write(u_log, '(a,f12.6)')   'DIAM    = ', sys%diam
    write(u_log, '(a,f12.6)')   'BONDL   = ', sys%bondl
    write(u_log, '(a,f12.6)')   'd/b     = ', sys%diam / sys%bondl
    write(u_log, '(a,i8)')      'NBLOCK  = ', mc_params%nblock
    write(u_log, '(a,i8)')      'NSTEP   = ', mc_params%nstep
    write(u_log, '(a,i8)')      'KSAMPLE = ', mc_params%ksample
    write(u_log, '(a,i8)')      'NTRIALS = ', mc_params%ntrials
    write(u_log, '(a,i8)')      'ISEED   = ', mc_params%iseed
    flush(u_log)
    close(u_log)

    write(*,'(a,i4,a,i8,a,f8.4,a,f8.4)') &
      'NC=', sys%nc, '  NP=', sys%np, '  d=', sys%diam, '  b=', sys%bondl

  end subroutine

  ! Step 1 (before setup_system): always initialise RNG from scratch.
  ! This ensures rng_uniform works when init_ring is called.
  subroutine rng_init_default(basename)
    character(len=*), intent(in) :: basename
    call rng_init(kl = 9373 + mc_params%iseed)
  end subroutine

  ! Step 2 (after setup_system): if a SEED file exists, override the RNG
  ! state with it.  The unformatted read in rng_load corrupts gfortran 15's
  ! list-directed float parser, so it must happen AFTER all text reads.
  subroutine rng_load_seed(basename)
    character(len=*), intent(in) :: basename
    logical :: exists
    inquire(file=trim(basename)//'.SEED', exist=exists)
    if (exists) call rng_load(trim(basename)//'.SEED')
  end subroutine

end module mod_input
