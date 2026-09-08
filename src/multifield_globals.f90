module multifield_globals

implicit none

!===========================================================================================================================
! Constants
real, parameter :: pi = 4.0 * atan(1.0)
!===========================================================================================================================
! Time steps
real, parameter :: dt_back = 1.0d2 ! Background time step
real, parameter :: dt_pert = 1.0d1 ! Perturbation time step
!===========================================================================================================================
! Field-space metric
! "EUM" for Euclidean metric
! "RPM" for Renaux-Petel metric
! "AAM" for alpha-attractor metric
! "OWM" for oscillating warp metric
! "GBM" for Gaussian bumps metric
character(len=*), parameter :: field_space_geometry = "EUM"
!===========================================================================================================================
! Potential
! "ELP" for elliptic potential
! "NLP" for non-linear potential
! "HYP" for hybrid potential
! "MVP" for modulated valley potential
! "DEP" for deformed elliptic potential
character(len=*), parameter :: potential_shape = "ELP"
!===========================================================================================================================
! Elliptic potential parameters
real, parameter :: m1_ELP = 1.4d-6 ! 1.4d-6
real, parameter :: m2_ELP = 15.0 * m1_ELP
!===========================================================================================================================
! Non-linear potential parameters
real, parameter :: lambda_NLP = 1.d-14
real, parameter :: g_NLP = lambda_NLP
real, parameter :: mu_NLP = 0.0
!===========================================================================================================================
! Hybrid inflation potential parameters
real, parameter :: chi0_HYP = 2.5
real, parameter :: m2_HYP = 1.47d-5
real, parameter :: m1_HYP = 0.3 * m2_HYP
real, parameter :: g_HYP = 0.8 * m2_HYP
real, parameter :: mu_HYP = 0.0
!===========================================================================================================================
! Modulated valley potential parameters
real, parameter :: m_MVP = 1.4d-6
real, parameter :: f1_MVP = 5.00 ! 1.00
real, parameter :: f2_MVP = 12.0 ! 4.00
real, parameter :: f3_MVP = 0.50 ! 2.00
real            :: lambda1_MVP ! = 0.05 * sqrt(m_mvp) ! 2.0d-3
!===========================================================================================================================
! Deformed elliptic potential parameters
real, parameter :: m1_DEP = 1.4d-6
real, parameter :: m2_DEP = 1.00 * m1_DEP
real, parameter :: amp_DEP = 5.00
real, parameter :: frq_DEP = 1.00
!===========================================================================================================================
! $\alpha$-attractor metric parameters
real, parameter :: alpha_aa = 1.0
real, parameter :: beta_aa = 1.0d10
!===========================================================================================================================
! Oscillatory warp metric parameters
real, parameter :: amp_OWM = 0.90 ! 0.92
real, parameter :: frq_OWM = 5.00 ! 5.00
!===========================================================================================================================
! Localized gaussian bumps metric number of bumps (assumes a domain in $\varphi^{1}$ from -20 to 20)
integer            :: c_GBM
integer, parameter :: n_GBM = 2 * 30
real, parameter    :: spp_GBM = 5.00
real               :: left_GBM(n_GBM/2), amp_GBM(n_GBM)
real, parameter    :: std_GBM(n_GBM) = 0.06
real, parameter    :: x_GBM(n_GBM) = [ ( -19.0 + (c_GBM-1) * (38.0 / real(n_GBM-1)), c_GBM = 1, n_GBM ) ]
!===========================================================================================================================
! Weierstrass metric parameters
real, parameter :: dsp_WNM = 1.9
real, parameter :: amp_WNM = 0.2
real, parameter :: a_WNM = 1.3
real, parameter :: b_WNM = 5.0
real, parameter :: c_WNM = 0.25
real, parameter :: d_WNM = 0.0
!===========================================================================================================================
! Global variables set at runtime
real    :: M_RPM       ! Renaux-Petel metric parameter
real    :: N_end       ! E-fold number at the end of inflation
logical :: condition   ! Loop condition
real    :: amp_factor  ! Variable amplitude factor for metrics and potentials
!===========================================================================================================================
! Random field potential file name
character(len=1024) :: potential_filename

end module multifield_globals
