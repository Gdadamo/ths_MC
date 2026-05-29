module mod_observables
  ! Per-block observables for THS ring/linear polymers.
  !
  !  Rg2  = (1/NP) Σ |r_i - r_cm|²        (both topologies)
  !  Ree2 = |r_1 - r_NP|²                  (linear only; 0 for ring)
  !  Rmax = max_i |r_i - r_cm|              (used as sampling radius for B2)
  !  B2   = -0.5 * <f_ij> * V_sphere        (NC=2 only, MC sampling)
  !
  ! B2 is computed by placing chain 2 at random positions B within a sphere
  ! of radius R12 = Rmax(1)+Rmax(2)+diam centred on the CM of chain 1,
  ! then counting overlapping configurations.
  !   B2 = -0.5 * (n_overlap / n_trials) * V_sphere
  use mod_params, only: dp, FPIOT
  use mod_system, only: sys, chain_t
  use mod_overlap, only: overlap_pair
  use mod_rng,    only: rng_uniform, rng_unit_vector
  implicit none
  private

  public :: compute_observables, obs_t

  type, public :: obs_t
    real(dp) :: rg2(2)    ! radius of gyration² per chain
    real(dp) :: ree2(2)   ! end-to-end² per chain (0 for rings)
    real(dp) :: rmax(2)   ! max distance from CM per chain
    real(dp) :: b2        ! second virial coefficient (0 if NC=1)
  end type

contains

  subroutine compute_observables(ntrials, obs)
    integer,    intent(in)  :: ntrials
    type(obs_t), intent(out) :: obs

    integer  :: ic
    real(dp) :: rcm(3)

    obs%rg2  = 0.0_dp
    obs%ree2 = 0.0_dp
    obs%rmax = 0.0_dp
    obs%b2   = 0.0_dp

    ! Centre and compute Rg², Ree², Rmax per chain
    do ic = 1, sys%nc
      call chain_observables(sys%chain(ic), obs%rg2(ic), obs%ree2(ic), obs%rmax(ic))
    end do

    ! Second virial coefficient (NC=2 only)
    if (sys%nc == 2) then
      call compute_b2(ntrials, obs%rmax, obs%b2)
    end if

  end subroutine

  ! ----------------------------------------------------------------
  ! Rg², Ree², Rmax for a single chain (recentered on CM in place).
  ! ----------------------------------------------------------------
  subroutine chain_observables(ch, rg2, ree2, rmax)
    type(chain_t), intent(inout) :: ch   ! positions shifted to CM frame
    real(dp),      intent(out)   :: rg2, ree2, rmax

    integer  :: i
    real(dp) :: rcm(3), rr, dr(3)

    ! Centre of mass
    rcm = 0.0_dp
    do i = 1, ch%np
      rcm = rcm + ch%r(:, i)
    end do
    rcm = rcm / real(ch%np, dp)

    ! Shift positions to CM frame, compute Rg² and Rmax
    rg2  = 0.0_dp
    rmax = 0.0_dp
    do i = 1, ch%np
      ch%r(:, i) = ch%r(:, i) - rcm
      rr = sum(ch%r(:, i)**2)
      rg2 = rg2 + rr
      if (rr > rmax) rmax = rr
    end do
    rg2  = rg2 / real(ch%np, dp)
    rmax = sqrt(rmax)

    ! Update ghost nodes after CM shift
    if (ch%is_ring) then
      ch%r(:, 0)       = ch%r(:, ch%np)
      ch%r(:, ch%np+1) = ch%r(:, 1)
    end if

    ! End-to-end distance (0 for rings)
    if (.not. ch%is_ring) then
      dr   = ch%r(:, 1) - ch%r(:, ch%np)
      ree2 = sum(dr * dr)
    else
      ree2 = 0.0_dp
    end if

  end subroutine

  ! ----------------------------------------------------------------
  ! B2 via MC sampling (Mayer f-function).
  !   Place chain 2 at random displacement B (uniform in sphere of radius R12)
  !   relative to chain 1. B2 = -0.5 * <overlap> * V_sphere.
  ! ----------------------------------------------------------------
  subroutine compute_b2(ntrials, rmax, b2)
    integer,  intent(in)  :: ntrials
    real(dp), intent(in)  :: rmax(2)
    real(dp), intent(out) :: b2

    integer  :: it, ip, jp
    real(dp) :: r12, vol, zz, rr, b_vec(3), rip(3), rjp(3), n_ov
    logical  :: ov

    ! Sampling sphere radius: maximum possible contact distance
    r12 = rmax(1) + rmax(2) + sys%diam
    vol = FPIOT * r12 * r12 * r12   ! (4/3)π r12³

    n_ov = 0.0_dp

    do it = 1, ntrials
      ! Random displacement B uniformly in sphere of radius r12
      call rng_uniform(zz)
      rr = r12 * zz**(1.0_dp/3.0_dp)   ! radial sampling with r² weight
      call rng_unit_vector(b_vec)
      b_vec = b_vec * rr

      ! Check if chain 2 displaced by B overlaps with chain 1
      ov = .false.
      outer: do ip = 1, sys%chain(1)%np
        rip = sys%chain(1)%r(:, ip)
        do jp = 1, sys%chain(2)%np
          rjp = sys%chain(2)%r(:, jp) + b_vec
          if (overlap_pair(rip, rjp, sys%diam2)) then
            ov = .true.
            exit outer
          end if
        end do
      end do outer

      if (ov) n_ov = n_ov + 1.0_dp
    end do

    ! B2 = -½ ∫ f dr  with f=-1 on overlap, 0 elsewhere
    ! <f> = -n_ov/N  →  B2 = -½ × (-n_ov/N) × V = +½ × (n_ov/N) × V > 0
    b2 = 0.5_dp * (n_ov / real(ntrials, dp)) * vol

  end subroutine

end module mod_observables
