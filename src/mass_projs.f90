program mass_projs
! Solves the set of coupled equations for an inflation model with two scalar fields with a non-trivial field space metric. This script evolves 6 functions of time stored in the entries of the array y(1:6).
! y(1:2)   stores $\phi^{A}$
! y(3:4)   stores $\dot{\phi}^{1}$
! y(5)     stores $H$
! y(6)     stores $N$

use multifield_globals
use multifield_grf
use multifield_utils

implicit none

real, parameter      :: N_bound = 200.0               ! Upper bound in N
real, parameter      :: dN = 0.01                     ! Data flushing period
real, dimension(6)   :: y_back                        ! State array
real, dimension(2)   :: phi, phidot                   ! Field multiplet
real                 :: H, Hdot, N, slowroll, N_flush ! Variables
real, dimension(2,2) :: M_XY                          ! Variables
real, dimension(2,2) :: M_AB, e                       ! Mass matrix $\mathcal{M}^{2}_{AB}$
integer              :: i, j, k                       ! Indices
character(len=32)    :: arg                           ! Command-line argument
character(len=100)   :: filename                      ! Output file name

! Load FITS file containing random potential
call get_command_argument(1, potential_filename)
call load_grf(potential_filename)

! Background initial conditions
phi    = initial_phi          ! $\phi^{A}(t_{0})$
phidot = terminal_phidot(phi) ! $\dot{\phi}^{A}(t_{0})$
H      = Hubble(phi, phidot)  ! $H(t_{0})$
N      = 0.0                  ! $N(t_{0})$

! Initialize background arrays
call pack_state_background(y_back, phi, phidot, H, N)

slowroll = 0.0     ! Initialize slow-roll parameter $\epsilon(t_{0})$
N_flush  = 0.0     ! Data writing trigger
condition = .true. ! Loop condition

do while (condition)
	! Update functions of time
	call unpack_state_background(y_back, phi, phidot, H, N)
	
	! Update slow-roll parameter
	slowroll = - Hubbledot(phi, phidot) / H**2
	
	! Update mass matrix
	M_AB = mass_matrix(phi, phidot, H)
	
	! Update vielbein adiabatic-isocurvature
	e = vielbein_ad_is(phi, phidot)
	
	! Squared mass matrix projections
	M_XY = matmul(transpose(e) , matmul(M_AB, e))
	
	if (N >= N_flush) then
		write (*, '(6(6e25.10e3))') N, M_XY, slowroll
		N_flush = N_flush + dN
	end if
	
	call gl8_background(y_back, dt_back)
	
	! Update condition
	if (slowroll >= 1.0) condition = .false.
end do

end program mass_projs
