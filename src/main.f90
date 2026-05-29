program ths_mc
  use mod_params,      only: dp
  use mod_system,      only: sys, dealloc_chain
  use mod_input,       only: mc_params, read_input, rng_init_default, rng_load_seed
  use mod_setup,       only: setup_system
  use mod_overlap,     only: check_chain_valid
  use mod_crank,       only: crank_move, crank_init, crank_free
  use mod_pivot,       only: pivot_move, pivot_init, pivot_free
  use mod_observables, only: obs_t, compute_observables
  use mod_io,          only: io_open, io_close, io_write_block, &
                              io_write_traj, io_write_restart
  implicit none

  character(len=300) :: basename
  integer  :: ib, is, ic
  real(dp) :: acc_move, acc_sum(2)
  type(obs_t) :: obs
  real :: t0, t1

  call get_command_argument(1, basename)
  if (len_trim(basename) == 0) then
    write(*,*) 'Usage: ths_mc <basename>  (reads <basename>.IN)'
    stop 1
  end if

  call cpu_time(t0)

  ! ---- Initialisation ----
  call read_input(trim(basename))
  call rng_init_default(trim(basename))  ! must be before setup_system (init_ring needs RNG)
  call setup_system(trim(basename))
  call rng_load_seed(trim(basename))     ! override with saved state AFTER text reads (gfortran 15 workaround)

  do ic = 1, sys%nc
    if (.not. check_chain_valid(sys%chain(ic))) then
      write(*,'(a,i2,a)') ' WARNING: chain ', ic, &
        ' has overlaps in initial config — run will self-heal via MC rejection'
    end if
  end do

  call crank_init(sys%np)
  call pivot_init(sys%np)

  call io_open(trim(basename), mc_params%ntrials, mc_params%ksample)

  write(*,'(a)') ' Starting MC run...'

  do ib = 1, mc_params%nblock

    acc_sum = 0.0_dp

    do ic = 1, sys%nc
      do is = 1, mc_params%nstep

        if (sys%chain(ic)%is_ring) then
          call crank_move(sys%chain(ic), acc_move)
        else
          call pivot_move(sys%chain(ic), acc_move)
        end if

        acc_sum(ic) = acc_sum(ic) + acc_move

      end do
    end do

    acc_sum = acc_sum / real(mc_params%nstep, dp)

    call compute_observables(mc_params%ntrials, obs)

    call io_write_block(ib, obs, acc_sum)

    if (mod(ib, mc_params%ksample) == 0) then
      call io_write_traj(ib)
    end if

    call io_write_restart(trim(basename))

    if (mod(ib, max(1, mc_params%nblock/10)) == 0) then
      write(*,'(a,i8,a,i8,a,2f8.4)') &
        ' Block ', ib, '/', mc_params%nblock, '  acc= ', acc_sum(1), acc_sum(2)
    end if

  end do

  call crank_free()
  call pivot_free()
  do ic = 1, sys%nc
    call dealloc_chain(sys%chain(ic))
  end do
  call io_close(trim(basename))

  call cpu_time(t1)
  write(*,'(a,f10.2,a)') ' Done. CPU time: ', t1 - t0, ' s'

end program ths_mc
