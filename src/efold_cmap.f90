program efold_cmap
! Solves the set of coupled equations for an inflation model with two scalar fields with a non-trivial field space metric. This script evolves 6 functions of time stored in the entries of the array y(1:6).
! y(1:2)   stores $\phi^{A}$
! y(3:4)   stores $\dot{\phi}^{1}$
! y(5)     stores $H$
! y(6)     stores $N$

use multifield_globals
use multifield_grf
use multifield_utils

implicit none

real, parameter    :: N_bound = 100.0                     ! Upper bound in N
integer, parameter :: ngrid = 200                         ! Number of grid points
real               :: efold(ngrid, ngrid)                 ! E-fold number grid
real, dimension(6) :: y                                   ! State array
real, dimension(2) :: phi, phidot, phi_min, phi_max, dphi ! Colormap grid parameters
real               :: H, Hdot, N, slowroll                ! Variables
integer            :: i, j                                ! Indices
character(len=32)  :: arg                                 ! Command-line argument

! Parse argument
if (field_space_geometry == "RPM") then
	call get_command_argument(1, arg)
	read (arg, *) M_RPM
end if

! Define the range for initial conditions
phi_min = [-lx, -ly]
phi_max = [ lx,  ly]

! Define the step size (NGRID-1 intervals between NGRID points)
dphi = [(phi_max(1) - phi_min(1)) / real(ngrid),  (phi_max(2) - phi_min(2)) / real(ngrid)]

!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(dynamic) PRIVATE(i, j, phi, phidot, H, N, slowroll, y) SHARED(phi_min, phi_max, dphi, efold)
do j = 1, ngrid
	do i = 1, ngrid
		! Background initial conditions
		phi    = [phi_min(1) + (real(i)-0.5) * dphi(1), phi_min(2) + (real(j)-0.5) * dphi(2)] ! $\phi^{A}(t_{0})$
		phidot = [0.0, 0.0]                                                                   ! $\dot{\phi}^{A}(t_{0})$
		H      = Hubble(phi, phidot)                                                          ! $H(t_{0})$
		N      = 0.0                                                                          ! $N(t_{0})$

		! Initialize background arrays
		call pack_state_background(y, phi, phidot, H, N)

		! Initialize slow-roll parameter $\epsilon(t_{0})$
		slowroll = 0.0

!		do while (slowroll <= 1.0 .and. N <= N_bound)
		do while (slowroll <= 1.0)
			! Update functions of time
			call unpack_state_background(y, phi, phidot, H, N)
			
			! Update slow-roll parameter
			slowroll = - Hubbledot(phi, phidot) / H**2
			
			call gl8_background(y, dt_back)
		end do
		
		efold(i,j) = N
	end do
end do
!$OMP END PARALLEL DO

do j = 1, ngrid
	write (*, '(800f10.4)') (efold(i,j), i = 1, ngrid)
end do

end program efold_cmap
