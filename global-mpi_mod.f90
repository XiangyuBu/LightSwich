MODULE global_parameters
IMPLICIT NONE

DOUBLE PRECISION, PARAMETER :: pi=3.1415926535897932384D0
DOUBLE PRECISION, PARAMETER :: eps = 1d-10

type node
	DOUBLE PRECISION :: x
	DOUBLE PRECISION :: y
	DOUBLE PRECISION :: z
end type
type(node),allocatable :: polymer(:,:)     !polymer in the ball 
type(node),allocatable :: bond_vector(:)
type(node),allocatable :: azo(:,:)         !azo-polymer in the substrate
type(node),allocatable :: sub(:,:)

Integer :: Nm, N_chain, N_azo, jj, i_azo_temp_1, i_azo_temp_2 ! jj is the chain will be rotated
Integer :: Nm_sub, N_sub, Nm_pol
Integer :: Nz, Nr, i_interface
Integer :: Movestep    !each MC attempt move
Integer :: num
Double Precision :: dr, Lr, Lz, Lz_2, dz, L_interface, np_r
Double Precision :: r_sphere_2,r_sphere, roL, rho_0
Double Precision :: move_max, rotate, rotate_s,hahah,csoL

Integer, DIMENSION(:,:), ALLOCATABLE :: ir, iz, ir_azo, iz_azo, ir_sub, iz_sub

DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: w, w_new
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: eta, eta_new,eta_azo,eta_azo_new,eta_sub,eta_sub_new

Integer, DIMENSION(:), ALLOCATABLE :: i_azo
Double PRECISION :: r_dr,r_dz, p_sphere_2
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: r_a, z_i, phi_rtot, phi_r, rr_r, zz_r
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: trial_move_rate

DOUBLE PRECISION :: nu, tau, epsilon, epsilon_azo, Lbox,azo_position_1,azo_position_2
DOUBLE PRECISION ::  Loa, deltaS,lambda
Integer :: n_iter, Max_iter, N_pre, Npre, Nmove, moves, NMCs,NMCstot, MCS, ncount
Integer :: seed
Integer :: comm, ierr, myid, numprocs  !MPI parameter
logical ::keepruning   !if keep MC simulating in the w field ,keepruning=.true.
!while the Npre will times 10 to keep a longer MC simulation time ,you can times
!a much longer time if you want


END MODULE global_parameters