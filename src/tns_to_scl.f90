program tns_to_scl
! Evolves background, scalar amplitudes and a tensor amplitude for each injected mode.
! Calculates the tensor-to-scalar ratio as a function of log(k/(aH)) at the end of inflation.
! Reduced Planck mass = 1. Both tensor polarizations are included in the tensor spectrum.
! y(1:2) stores the background fields
! y(3:4) stores their cosmic-time derivatives
! y(5)   stores $H$
! y(6)   stores $N$
! y(7)   stores $r_{T}$
! y(8)   stores $r_{T}^{\prime}$ (conformal-time derivative)

use multifield_globals
use multifield_grf
use multifield_utils

implicit none

real, parameter    :: N_start = 2.5, N_stop = 100.0
real, parameter    :: N_step = 0.25
integer, parameter :: ngrid_max = int((N_stop - N_start) / N_step) + 1
real               :: y_back(6), y_pert(8)
real               :: phi(2), phidot(2)
real               :: phi_grid(2,ngrid_max), phidot_grid(2,ngrid_max)
real               :: H_grid(ngrid_max), N_grid(ngrid_max)
real               :: x_grid(ngrid_max), PT(ngrid_max)
real               :: H, H_end, N, slowroll
real               :: k_mode, k_phys, N_trigger
real               :: r_T, r_T_p, W2_T, omega_T
integer            :: i, ngrid

! Load FITS file containing random potential
call get_command_argument(1, potential_filename)
call load_grf(potential_filename)

!========================================================================================================
! Background
!========================================================================================================
! Background initial conditions
phi    = initial_phi
phidot = terminal_phidot(phi)
H      = Hubble(phi, phidot)
N      = 0.0

! Initialize background arrays
call pack_state_background(y_back, phi, phidot, H, N)

i         = 0
N_trigger = N_start
condition = .true.

do while (condition)
    call unpack_state_background(y_back, phi, phidot, H, N)

    slowroll = -Hubbledot(phi, phidot) / H**2

    ! Same injection ordering as ps_full.
    if (N >= N_trigger .and. N <= N_stop .and. i < ngrid_max) then
        i                = i + 1
        phi_grid(:,i)    = phi
        phidot_grid(:,i) = phidot
        H_grid(i)        = H
        N_grid(i)        = N
        N_trigger        = N_trigger + N_step
    end if

    call gl8_background(y_back, dt_back)

    if (slowroll >= 1.0) condition = .false.
end do

N_end = N
H_end = H
ngrid = i

!========================================================================================================
! Perturbations
!========================================================================================================
! Same injection wavenumber prescription as ps_full.
k_phys = 1.0d5 * H_end

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i, phi, phidot, H, N, k_mode, y_pert, r_T, r_T_p, W2_T, omega_T) SHARED(ngrid, phi_grid, phidot_grid, H_grid, N_grid, k_phys, N_end, x_grid, PT) SCHEDULE(dynamic)
do i = 1, ngrid
    phi    = phi_grid(:,i)
    phidot = phidot_grid(:,i)
    H      = H_grid(i)
    N      = N_grid(i)
    k_mode = exp(N) * k_phys

    ! Initial tensor vacuum, identical to the original program.
    W2_T = k_mode**2 - exp(2.0*N) * (2.0*H**2 + Hubbledot(phi, phidot))

    if (W2_T <= 0.0) then
        error stop 'Initial tensor frequency squared must be positive'
    end if

    omega_T = sqrt(W2_T)
    r_T     = 1.0 / sqrt(2.0 * omega_T)
    r_T_p   = 0.0

    call pack_state_tns_to_scl(y_pert, phi, phidot, H, N, r_T, r_T_p)

    ! Preserve ps_full's stepping and stopping convention.
    do
        call gl8_tns_to_scl(y_pert, dt_pert, k_mode)
        if (y_pert(6) >= N_end) exit
    end do

    call unpack_state_tns_to_scl(y_pert, phi, phidot, H, N, r_T, r_T_p)

    x_grid(i) = log(k_mode/(exp(N)*H))
    PT(i)     = 4.0 * k_mode**3 * exp(-2.0*N) * r_T**2 / pi**2
end do
!$OMP END PARALLEL DO

! Columns: $\log(k/(aH))$, $\mathcal{P}_{\rm T}$ including both tensor polarizations.
do i = 1, ngrid
    write (*,'(2e25.10e3)') x_grid(i), PT(i)
end do

end program tns_to_scl
