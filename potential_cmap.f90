program potential_cmap
! Solves the set of coupled equations for an inflation model with two scalar fields with a non-trivial field space metric. This script evolves 6 functions of time stored in the entries of the array y(1:6).
! y(1:2)   stores $\phi^{A}$
! y(3:4)   stores $\dot{\phi}^{1}$
! y(5)     stores $H$
! y(6)     stores $N$

use multifield_globals
use multifield_utils

implicit none

real, parameter    :: dt = 5000.0                 ! Time-step
integer, parameter :: ngrid = 500                 ! Number of grid points
real, parameter    :: half_square_side = 15.0     ! Grid size
real               :: Vgrid(ngrid, ngrid)         ! E-fold number grid
real, dimension(2) :: phi, phi_min, phi_max, dphi ! Colormap grid parameters
integer            :: i, j                        ! Indices
character(len=32)  :: arg                         ! Command-line argument

! Load FITS file containing random potential
call get_command_argument(1, potential_filename)
call load_grf(potential_filename)

! Define the grid domain
!phi_min = [-half_square_side, -half_square_side]
!phi_max = [ half_square_side,  half_square_side]
phi_min = [-lx, -ly]
phi_max = [ lx,  ly]

! Define the step size (NGRID-1 intervals between NGRID points)
dphi = [(phi_max(1) - phi_min(1)) / real(ngrid),  (phi_max(2) - phi_min(2)) / real(ngrid)]

!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(dynamic) PRIVATE(i, j, phi) SHARED(phi_min, phi_max, dphi)
do j = 1, ngrid
	do i = 1, ngrid
		phi        = [phi_min(1) + (real(i)-0.5) * dphi(1), phi_min(2) + (real(j)-0.5) * dphi(2)] ! $\phi^{A}(t_{0})$
		Vgrid(i,j) = potential(phi)
	end do
end do
!$OMP END PARALLEL DO

do j = 1, ngrid
	write (*,'(500(6e25.10e3))') (Vgrid(i,j), i = 1, ngrid)
end do

end program potential_cmap
