module multifield_grf

use multifield_globals

implicit none

integer :: nx, ny
real    :: lx, ly
real    :: dx, dy

real :: initial_phi(2)

real, allocatable :: Vmap(:,:)
real, allocatable :: Vxmap(:,:)
real, allocatable :: Vymap(:,:)
real, allocatable :: Vxxmap(:,:)
real, allocatable :: Vxymap(:,:)
real, allocatable :: Vyymap(:,:)

contains

! Downsample a (nx_full x ny_full) map to (nx_out x ny_out) by box-averaging.
! nx_full must be an integer multiple of nx_out, same for ny.
subroutine downsample(full, nx_full, ny_full, out, nx_out, ny_out)
	integer, intent(in) :: nx_full, ny_full, nx_out, ny_out
	real,    intent(in) :: full(nx_full, ny_full)
	real,    intent(out):: out(nx_out, ny_out)
	integer :: i, j, ii, jj, bx, by
	real    :: s
	
	bx = nx_full / nx_out   ! block size in x  (e.g. 2)
	by = ny_full / ny_out   ! block size in y  (e.g. 2)
	
	do j = 1, ny_out
		do i = 1, nx_out
			s = 0.0
			do jj = 1, by
				do ii = 1, bx
					s = s + full((i-1)*bx + ii, (j-1)*by + jj)
				end do
			end do
			out(i,j) = s / real(bx * by)
		end do
	end do
end subroutine downsample

! Load a FITS file containing a random potential and its gradient and Hessian
subroutine load_grf(filename)
	character(*), intent(in) :: filename
	real                     :: nullval
	integer                  :: status, blocksize, unit, anynull, hdutype, naxes(2)
	integer                  :: nx_full, ny_full
    integer                  :: scaled, scaled_status
    real                     :: spline_shift_x, spline_shift_y
	integer, parameter       :: downgrade_factor = 1
	real, allocatable        :: tmp(:,:)
	! For reading the spline curve from HDU 8 (index 9 in CFITSIO 1-based HDU count)
	integer                  :: nrows, felem
	real, allocatable        :: col_x(:), col_y(:)
	integer                  :: colnum_x, colnum_y
    real, parameter          :: potential_factor = 1.d-9
	
	status    = 0
	blocksize = 0
	hdutype   = 0
	nullval   = 0.0
	
	call ftgiou(unit, status)
	call ftopen(unit, filename, 0, blocksize, status)
	
	call ftgkyd(unit, 'LX', lx, '', status)
	call ftgkyd(unit, 'LY', ly, '', status)
	call ftgkyj(unit, 'NX', nx_full, '', status)
	call ftgkyj(unit, 'NY', ny_full, '', status)
    
    ! Fourier refinement factor; an absent SCALED means no refinement.
    scaled = 1
    scaled_status = 0
    call ftgkyj(unit, 'SCALED', scaled, '', scaled_status)

    if (scaled_status == 202) then
        ! CFITSIO: keyword not found.
        scaled = 1
    else if (scaled_status /= 0) then
        error stop 'Could not read SCALED from the FITS header.'
    end if

    if (scaled < 1) then
        error stop 'SCALED must be a positive integer.'
    end if

    ! Express the original spline coordinates in the pixel-center frame.
    ! Use the FITS resolution BEFORE any optional downsampling.
    spline_shift_x = real(scaled - 1) * lx / real(nx_full)
    spline_shift_y = real(scaled - 1) * ly / real(ny_full)
	
	! Set module resolution to the downsampled grid
	nx = nx_full / downgrade_factor
	ny = ny_full / downgrade_factor
	dx = 2.0 * lx / nx
	dy = 2.0 * ly / ny
	
	allocate(tmp(nx_full, ny_full))
	allocate(Vmap(nx,ny))
	allocate(Vxmap(nx,ny))
	allocate(Vymap(nx,ny))
	allocate(Vxxmap(nx,ny))
	allocate(Vxymap(nx,ny))
	allocate(Vyymap(nx,ny))
	
	! HDU 1: V
	call ftmahd(unit, 2, hdutype, status)
	call ftg2dd(unit, 1, nullval, nx_full, nx_full, ny_full, tmp, anynull, status)
	call downsample(tmp, nx_full, ny_full, Vmap, nx, ny)
	
	! HDU 2: Vx
	call ftmahd(unit, 3, hdutype, status)
	call ftg2dd(unit, 1, nullval, nx_full, nx_full, ny_full, tmp, anynull, status)
	call downsample(tmp, nx_full, ny_full, Vxmap, nx, ny)
	
	! HDU 3: Vy
	call ftmahd(unit, 4, hdutype, status)
	call ftg2dd(unit, 1, nullval, nx_full, nx_full, ny_full, tmp, anynull, status)
	call downsample(tmp, nx_full, ny_full, Vymap, nx, ny)
	
	! HDU 4: Vxx
	call ftmahd(unit, 5, hdutype, status)
	call ftg2dd(unit, 1, nullval, nx_full, nx_full, ny_full, tmp, anynull, status)
	call downsample(tmp, nx_full, ny_full, Vxxmap, nx, ny)
	
	! HDU 5: Vxy
	call ftmahd(unit, 6, hdutype, status)
	call ftg2dd(unit, 1, nullval, nx_full, nx_full, ny_full, tmp, anynull, status)
	call downsample(tmp, nx_full, ny_full, Vxymap, nx, ny)
	
	! HDU 6: Vyy
	call ftmahd(unit, 7, hdutype, status)
	call ftg2dd(unit, 1, nullval, nx_full, nx_full, ny_full, tmp, anynull, status)
	call downsample(tmp, nx_full, ny_full, Vyymap, nx, ny)
	
	deallocate(tmp)
	
	! HDU 8 (CFITSIO HDU 9): SPLINE binary table
	! Contains the spline curve sampled at 1024 points; first row is s=0 (trajectory start)
	call ftmahd(unit, 9, hdutype, status)
	call ftgkyj(unit, 'NAXIS2', nrows, '', status)
	
	allocate(col_x(nrows))
	allocate(col_y(nrows))
	
	call ftgcno(unit, .false., 'X', colnum_x, status)
	call ftgcno(unit, .false., 'Y', colnum_y, status)
	
	felem = 1
	call ftgcvd(unit, colnum_x, 1, felem, nrows, nullval, col_x, anynull, status)
	call ftgcvd(unit, colnum_y, 1, felem, nrows, nullval, col_y, anynull, status)
	
    ! Start at the spline origin, expressed in the pixel-center frame.
    initial_phi(1) = col_x(1) - spline_shift_x
    initial_phi(2) = col_y(1) - spline_shift_y
	
	deallocate(col_x, col_y)
	
	call ftclos(unit, status)
	call ftfiou(unit, status)
	
	Vmap   = potential_factor * Vmap
	Vxmap  = potential_factor * Vxmap
	Vymap  = potential_factor * Vymap
	Vxxmap = potential_factor * Vxxmap
	Vxymap = potential_factor * Vxymap
	Vyymap = potential_factor * Vyymap
	
end subroutine load_grf

! Bilinear interpolation
pure function interp2(map, x, y) result(f)
	real, intent(in) :: map(:,:)
	real, intent(in) :: x, y
	real             :: f

	integer :: i, j
	real    :: tx, ty

	! --- locate cell ---
	tx = (x + lx - 0.5*dx) / dx
	ty = (y + ly - 0.5*dy) / dy

	i = floor(tx) + 1
	j = floor(ty) + 1

	! keep interpolation cell inside array
	i = max(1, min(i, nx-1))
	j = max(1, min(j, ny-1))

	tx = tx - floor(tx)   ! fractional part in [0,1)
	ty = ty - floor(ty)

	! --- bilinear interpolation ---
	f = (1.0-tx)*(1.0-ty)*map(i  ,j  ) &
	  + tx      *(1.0-ty)*map(i+1,j  ) &
	  + (1.0-tx)*ty      *map(i  ,j+1) &
	  + tx      *ty      *map(i+1,j+1)

end function interp2

! Bicubic interpolation with periodic boundaries (Fourier-generated maps).
pure function interp3(map, x, y) result(f)
	real, intent(in) :: map(:,:)
	real, intent(in) :: x, y
	real             :: f

	integer :: i, j, a, b, ix(4), iy(4), mx, my
	real    :: tx, ty, wx(4), wy(4)

	mx = size(map, 1)
	my = size(map, 2)

	! Locate cell in the pixel-center frame and wrap to the map.
	tx = modulo((x + lx - 0.5*dx) / dx, real(mx))
	ty = modulo((y + ly - 0.5*dy) / dy, real(my))
	i = floor(tx)
	j = floor(ty)

	wx = catmull_rom_weights(tx - real(i))
	wy = catmull_rom_weights(ty - real(j))

	do a = 1, 4
		ix(a) = modulo(i + a - 2, mx) + 1
		iy(a) = modulo(j + a - 2, my) + 1
	end do

	f = 0.0
	do b = 1, 4
		do a = 1, 4
			f = f + wx(a)*wy(b)*map(ix(a),iy(b))
		end do
	end do

end function interp3

! Catmull-Rom cardinal spline weights for the four neighbours.
! Returns w(1..4) for offsets -1, 0, +1, +2 given t in [0,1).
pure function catmull_rom_weights(t) result(w)
	real, intent(in) :: t
	real             :: w(4)
	real             :: t2, t3
	
	t2 = t*t
	t3 = t2*t
	
	! alpha = 0.5  (standard Catmull-Rom)
	w(1) =  0.5*(-t3 + 2.0*t2 - t)
	w(2) =  0.5*( 3.0*t3 - 5.0*t2 + 2.0)
	w(3) =  0.5*(-3.0*t3 + 4.0*t2 + t)
	w(4) =  0.5*( t3 - t2)
	
end function catmull_rom_weights

end module multifield_grf
