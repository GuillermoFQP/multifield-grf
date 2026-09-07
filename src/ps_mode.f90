program ps_mode
! Solves the set of coupled background and perturbation equations for an inflation model with two scalar fields with a non-trivial field space metric and calculates the curvature power spectrum of the perturbation fields $\mathcal{P}_{\mathcal{R}}$ as a function of $N$. This script evolves 28 functions of time stored in the entries of the array y(1:28).
! y(1:2)   stores $\phi^{A}$
! y(3:4)   stores $\dot{\phi}^{1}$
! y(5)     stores $H$
! y(6)     stores $N$
! y(7:10)  stores $\Lambda^{A}_{j}$
! y(11:12) stores $\mathrm{Re} \, \mathbf{r}_{1}$
! y(13:14) stores $\mathrm{Im} \, \mathbf{r}_{1}$
! y(15:16) stores $\mathrm{Re} \, \mathbf{r}_{2}$
! y(17:18) stores $\mathrm{Im} \, \mathbf{r}_{2}$
! y(19:20) stores $\mathrm{Re} \, \mathbf{r}_{1}^{\prime}$
! y(21:22) stores $\mathrm{Im} \, \mathbf{r}_{1}^{\prime}$
! y(23:24) stores $\mathrm{Re} \, \mathbf{r}_{2}^{\prime}$
! y(25:26) stores $\mathrm{Im} \, \mathbf{r}_{2}^{\prime}$

use multifield_globals
use multifield_grf
use multifield_utils

implicit none

real, parameter      :: N_trigger = 2.5                             ! Mode trigger
real, dimension(6)   :: y_back                                      ! Background state array
real, dimension(28)  :: y_pert                                      ! Background + perturbation state array
real, dimension(2)   :: phi, phidot                                 ! Field multiplet and its time derivative
real, dimension(2)   :: Re_r1, Im_r1, Re_r2, Im_r2                  ! Complex amplitude vector
real, dimension(2)   :: Re_r1_p, Im_r1_p, Re_r2_p, Im_r2_p          ! Time derivative of complex amplitude vector
real, dimension(2)   :: phi_0, phidot_0                             ! Initial conditions for mode-injection
real                 :: H_0, N_0                                    ! Initial conditions for mode-injection
real                 :: H, H_end, Hdot, N, slowroll, k_mode, k_phys ! Variables
real, dimension(2,2) :: phidotphidot, vbein_PT, W2_ij               ! Variables
real                 :: W_11, W_22, N_flush                         ! Variables
logical              :: trigger                                     ! Mode trigger and loop condition
character(len=32)    :: arg                                         ! Command-line argument

! Load FITS file containing random potential
call get_command_argument(1, potential_filename)
call load_grf(potential_filename)

!========================================================================================================
! Background
!========================================================================================================
! Background initial conditions
phi    = initial_phi          ! $\phi^{A}(t_{0})$
phidot = terminal_phidot(phi) ! $\dot{\phi}^{A}(t_{0})$
H      = Hubble(phi, phidot)  ! $H(t_{0})$
N      = 0.0                  ! $N(t_{0})$

! Initialize background arrays
call pack_state_background(y_back, phi, phidot, H, N)

trigger   = .true. ! Counter
slowroll  = 0.0    ! Initialize slow-roll parameter $\epsilon(t_{0})$
condition = .true. ! Initialize condition

do while (condition)
	! Update functions of time
	call unpack_state_background(y_back, phi, phidot, H, N)
	
	! Update slow-roll parameter
	slowroll = - Hubbledot(phi, phidot) / H**2
	
	if (N >= N_trigger .and. trigger) then
		trigger  = .false.
		phi_0    = phi
		phidot_0 = phidot
		H_0      = H
		N_0      = N
	end if
	
	call gl8_background(y_back, dt_back)
	
	! Update condition
!	if (slowroll >= 1.0 .and. convergence(y_back(6), N)) condition = .false.
	if (slowroll >= 1.0) condition = .false.
end do

N_end = N ! E-fold number at the end of inflation
H_end = H ! Hubble parameter at the end of inflation

!========================================================================================================
! Perturbations
!========================================================================================================
k_phys = 1.0d5 * H_end ! Physical wave vector magnitude

! Background initial conditions
phi      = phi_0                           ! $\phi^{A}(t_{0})$
phidot   = phidot_0                        ! $\dot{\phi}^{A}(t_{0})$
H        = H_0                             ! $H(t_{0})$
N        = N_0                             ! $N(t_{0})$
slowroll = - Hubbledot(phi, phidot) / H**2 ! $\epsilon(t_{0})$
vbein_PT = vielbein_ad_is(phi, phidot)     ! $\Lambda^{A}_{i}(t_{0})$

! Comoving wave vector magnitude
k_mode = exp(N) * k_phys

! Initial conditions for perturbations
W2_ij   = freq_matrix(phi, phidot, H, N, vbein_PT, k_mode) ! $\Omega^{2}_{ij}(t_{0})$
W_11    = sqrt(W2_ij(1,1))                                 ! $\sqrt{\Omega^{2}_{11}}$
W_22    = sqrt(W2_ij(2,2))                                 ! $\sqrt{\Omega^{2}_{22}}$
Re_r1   = [sqrt(2.0 * W_11)**(-1) , 0.0]
Re_r2   = [0.0, sqrt(2.0 * W_22)**(-1)]
Im_r1   = 0.0
Im_r2   = 0.0
Re_r1_p = 0.0
Re_r2_p = 0.0
Im_r1_p = 0.0
Im_r2_p = 0.0

! Initialize perturbation functions
call pack_state_perturbations(y_pert, phi, phidot, H, N, vbein_PT, Re_r1, Im_r1, Re_r2, Im_r2, Re_r1_p, Im_r1_p, Re_r2_p, Im_r2_p)

N_flush   = N      ! Initialize data writing trigger
condition = .true. ! Initialize condition

do while (condition)
	! Update functions of time
	call unpack_state_perturbations(y_pert, phi, phidot, H, N, vbein_PT, Re_r1, Im_r1, Re_r2, Im_r2, Re_r1_p, Im_r1_p, Re_r2_p, Im_r2_p)
	
	! Update slow-roll parameter
	slowroll = - Hubbledot(phi, phidot) / H**2
	
	if (N >= N_flush) then
		write (*,'(7(6e25.10e3))') N, powerspectrum(phi, phidot, H, N, vbein_PT, Re_r1, Im_r1, Re_r2, Im_r2, k_mode), N_end
		N_flush = N_flush + 0.01
	end if
	
	call gl8_perturbations(y_pert, dt_pert, k_mode)
	
	! Update condition
	if (y_pert(6) >= N_end) condition = .false.
end do

end program ps_mode
