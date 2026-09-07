program efold_cmap_test
! Solves the set of coupled equations for an inflation model with two scalar fields with a non-trivial field space metric.
! This script evolves 6 functions of time stored in the entries of the array y(1:6).
! y(1:2) stores $\phi^{A}$
! y(3:4) stores $\dot{\phi}^{1}$
! y(5)   stores $H$
! y(6)   stores $N$
!
! For each initial condition the background trajectory is evolved until one
! of four outcomes is reached:
!   (1) $\epsilon = 1$                                  -> efold = N
!   (2) $V(phi) \le 0$ is encountered                   -> efold = negV_flag
!   (3) the trajectory is certified to be captured at
!       an attracting local minimum of $V$ (genuine
!       eternal-inflation candidate); see the
!       attractivity + persistence test below           -> efold = eternal_flag
!   (4) the hard cap $N_{\rm max}$ is reached without
!		(1)-(3) triggering (safety net for
!		pathological/very long saddle transients)       -> efold = eternal_flag
!
! Attractivity + persistence test (outcome 3):
!   At an attracting critical point of $V$, linearizing the slow-roll flow
!   $phi^A(t)$ around it gives $\delta\phi \sim \exp(-\lambda N)$ along each
!   eigen-direction of the Hessian $\partial_{B}\partial_{A} V$. If all
!   eigenvalues are positive (a true local minimum, not a saddle) and V > 0
!   there, epsilon decays to zero forever and $\epsilon = 1$ is never reached. A
!   saddle point instead has a mixed-sign spectrum: trajectories stall near it
!   but eventually escape along the unstable direction, so requiring
!   positive-definiteness (not just small $\epsilon$) is what distinguishes
!   genuine eternal inflation from a long transient. The persistence
!   window (in e-folds, scaled by the local relaxation rate $1/\lambda_{\rm min}$)
!   then guards against a trajectory merely grazing close to a saddle
!   before being flung away.

use multifield_globals
use multifield_grf
use multifield_utils

implicit none

real, parameter      :: N_max = 100.0                       ! Hard cap on e-folds (safety net)
real, parameter      :: eps_th = 1.0d-3                     ! "epsilon is small" threshold
real, parameter      :: persist_factor = 8.0                ! Persistence window, in units of 1/lambda_min
real, parameter      :: deltaN_min = 20.0                   ! Floor on the persistence window (e-folds)
real, parameter      :: deltaN_max = 100.0                  ! Ceiling on the persistence window (e-folds)
real, parameter      :: eternal_flag = 1.0d3                ! Output value for eternal-inflation candidates
real, parameter      :: negV_flag = -1.0d3                  ! Output value for V <= 0 excursions
integer, parameter   :: check_stride = 1                    ! Evaluate the stopping criterion every this many GL8 steps
integer, parameter   :: ngrid = 1024                        ! Number of grid points per axis
real, allocatable    :: efold(:,:)                          ! E-fold number grid
real, dimension(6)   :: y                                   ! State array
real, dimension(2)   :: phi, phidot, phi_min, phi_max, dphi ! Colormap grid parameters
real                 :: H, Hdot, N, epsilon                 ! Variables
real                 :: Vpot, disc, trM                     ! Potential and Hessian discriminant and trace
real, dimension(2,2) :: DDV                                 ! h^{AB}, covariant Hessian M_{AB}, and h^{AC}M_{CB}
real                 :: lambda1, lambda2, lambda_min        ! Minimum eigenvalue of the Hessian of V
real                 :: N_window_start, deltaN_required     ! Persistence-window bookkeeping
logical              :: attracting, window_active           ! Persistence-window bookkeeping
integer              :: i, j, step                          ! Loop indices
integer              :: status                              ! Status parameter
character(len=32)    :: arg                                 ! Command-line argument

! Load FITS file containing random potential
call get_command_argument(1, potential_filename)
call load_grf(potential_filename)

! Create destination directory
!call execute_command_line('rm -rf trajs; mkdir trajs')

! Define the range for initial conditions
phi_min = [-lx, -ly]
phi_max = [ lx,  ly]

! Define the step size (NGRID-1 intervals between NGRID points)
dphi = [(phi_max(1) - phi_min(1)) / real(ngrid),  (phi_max(2) - phi_min(2)) / real(ngrid)]

allocate(efold(ngrid,ngrid))

!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(dynamic) &
!$OMP PRIVATE(i, j, step, status, phi, phidot, H, N, epsilon, y, Vpot, DDV, &
!$OMP         trM, disc, lambda1, lambda2, lambda_min, N_window_start, &
!$OMP         deltaN_required, attracting, window_active) &
!$OMP SHARED(phi_min, phi_max, dphi, efold)
do j = 1, ngrid
	do i = 1, ngrid
		! Background initial conditions
		phi    = [phi_min(1) + (real(i)-0.5) * dphi(1), phi_min(2) + (real(j)-0.5) * dphi(2)] ! $\phi^{A}(t_{0})$
		phidot = terminal_phidot(phi)                                                         ! $\dot{\phi}^{A}(t_{0})$
		H      = Hubble(phi, phidot)                                                          ! $H(t_{0})$
		N      = 0.0                                                                          ! $N(t_{0})$

		! Initialize background arrays
		call pack_state_background(y, phi, phidot, H, N)

		! Initialize slow-roll parameter $\epsilon(t_{0})$
		epsilon = 0.0

		! Persistence-window bookkeeping for the eternal-inflation test
		window_active  = .false.
		N_window_start = 0.0
		step           = 0
		
		! Status parameter
		! 0 = undetermined, 1 = genuine end of inflation, 2 = eternal inflation, 3 = negative potential local minimum
		status = 0

		do while (epsilon <= 1.0 .and. N <= N_max)
			! Update functions of time
			call unpack_state_background(y, phi, phidot, H, N)

			! Update slow-roll parameter (exact $epsilon = -Hdot/H^2$)
			epsilon = epsilon_sr(phi, phidot, H)

			step = step + 1
			if (mod(step, check_stride) == 0) then

				! A negative potential excursion makes epsilon/H ill-defined;
				! flag it separately from both endpoints and eternal candidates
				Vpot = potential(phi)
				if (Vpot <= 0.0) then
					status = 3
					exit
				end if

				! Eigenvalues of the mixed tensor $\partial_{B}\partial_{A} V$. These are
				DDV        = HessianV(phi)
				trM        = DDV(1,1) + DDV(2,2)             ! Trace
				disc       = max(trM**2 - 4.0 * det(DDV), 0.0) ! Clamp roundoff; analytically disc >= 0
				lambda1    = 0.5 * (trM + sqrt(disc))
				lambda2    = 0.5 * (trM - sqrt(disc))
				lambda_min = min(lambda1, lambda2)

				! Attractivity test
				! $\epsilon$ small and both Hessian eigenvalues positive (true local minimum, not a saddle)
				! $V > 0$ already guaranteed by the check above
				attracting = (epsilon < eps_th) .and. (lambda_min > 0.0)

				if (attracting) then
					if (.not. window_active) then
						window_active  = .true.
						N_window_start = N
					end if
					! Window length scales with the local relaxation time, clamped to a sane range
					deltaN_required = min(max(persist_factor / lambda_min, deltaN_min), deltaN_max)
					if (N - N_window_start >= deltaN_required) then
						status = 2
						exit
					end if
				else
					! Criterion (a) failed: reset the persistence window
					window_active = .false.
				end if

			end if

			call gl8_background(y, dt_back)
		end do

		! Loop ended without an early exit: classify by how it ended
		if (status == 0) then
			if (epsilon > 1.0) then
				status = 1 ! Genuine end of inflation
			else
				status = 2 ! Hit $N_{\rm max}$ without confirmed attractivity: flag for review
			end if
		end if

		select case (status)
		case (1)
			efold(i,j) = N
		case (2)
			efold(i,j) = eternal_flag
		case (3)
			efold(i,j) = negV_flag
		end select

	end do
end do
!$OMP END PARALLEL DO

do j = 1, ngrid
	write (*, '(*(f14.4))') (efold(i,j), i = 1, ngrid)
end do

deallocate(efold)

end program efold_cmap_test
