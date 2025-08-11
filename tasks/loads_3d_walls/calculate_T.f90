module heat_diffusion
  implicit none

#ifdef BE_WALL
    real*8,  parameter :: mat_density = 1.750d3 ! kg/m3
#else
    real*8,  parameter :: mat_density = 19.254d3 ! kg/m3
#endif

  contains

  function heat_diffusion_step(T_curr, nx, dx, dt, q) result(T_next)
    ! Input parameters
    integer, intent(in) :: nx        ! Number of spatial grid points
    real(8), intent(in) :: dx, dt    ! Spatial and temporal step sizes
    real(8), intent(in) :: q         ! Heat flux at x=0
    real(8), intent(in), dimension(nx) :: T_curr  ! Current temperature distribution

    ! Output
    real(8), dimension(nx) :: T_next ! Temperature distribution after one time step

    ! Local variables
    real(8) :: r, dTdx, alpha, T0, TL, dalpha_dT, Ti
    integer :: i

    ! Boundary conditions (x=0)
    T0         = T_curr(1)
    dTdx       = -q / heat_conductivity(T0)
    alpha      = heat_conductivity(T0) / (mat_density*spec_heat_capacity(T0))
    dalpha_dT  = alpha * ( dheat_conductivity_dT(T0)/heat_conductivity(T0)) ! - dspec_heat_capacity_dT(T0)/spec_heat_capacity(T0) )
    r          = alpha * dt / (dx * dx)
    T_next(1)  = T_curr(1) + 2.d0* r * (T_curr(2) - T_curr(1) - dTdx*dx)  & ! Neumann BC: heat flux at x=0
                 + dalpha_dT * dTdx**2 * dt

    ! Boundary conditions (x=L)
    TL         = T_curr(nx)
    alpha      = heat_conductivity(TL) / (mat_density*spec_heat_capacity(TL))
    r          = alpha * dt / (dx * dx)
    T_next(nx) = T_curr(nx) + 2.d0* r * (T_curr(nx-1) - T_curr(nx) )       ! Neumann BC: heat flux=0 at x=L

    ! Update temperature for internal points
    do i = 2, nx - 1
      Ti        = T_curr(i)
      alpha     = heat_conductivity( Ti ) / (mat_density*spec_heat_capacity( Ti ))
      dalpha_dT = alpha * ( dheat_conductivity_dT(Ti)/heat_conductivity(Ti) )!- dspec_heat_capacity_dT(Ti)/spec_heat_capacity(Ti) )
      r         = alpha * dt / (dx * dx)
      T_next(i) = T_curr(i) + r * (T_curr(i+1) - 2.0 * T_curr(i) + T_curr(i-1))   &
                  + dt* dalpha_dT / (4.d0*dx*dx) * (T_curr(i+1)-T_curr(i-1))**2
    end do

  end function heat_diffusion_step


  pure function spec_heat_capacity(T) result(cp)

    implicit none
    
    ! Input parameters
    real(8), intent(in) :: T        ! Temperature in K
    real(8) :: cp, TC   ! specific capacity [J / (kg K)]

#ifdef BE_WALL
    TC = T - 273.15
    cp = -1.338948d-11*TC**4 + 1.310600d-6*TC**3 - 3.143247d-3*TC**2 + 3.344647d0*TC + 1.741434d3
    if (TC > 1283) cp = 3590
    if (TC <   20) cp = 1807
#else
    ! DOI: 10.1007/978-94-007-7587-9_3
    cp = 135.76 + 9.1159d-3*T + 2.3134d-9*T**3 - 6.5233d5 /T**2
    cp = max(cp,129.d0) ! Min for 273
    cp = min(cp,248.d0) ! Max for 3300
#endif

    !cp = 133.d0
  
  end function spec_heat_capacity


  pure function dspec_heat_capacity_dT(T) result(cp)

    implicit none
    
    ! Input parameters
    real(8), intent(in) :: T        ! Temperature in K
    real(8) :: cp, TC   ! specific capacity [J / (kg K)]

#ifdef BE_WALL
    TC = T - 273.15
    cp = -5.355792d-11*TC**3 + 3.9318d-6*TC**2 - 6.286494d-3*TC + 3.344647d0
    if (TC > 1283) cp = 0
    if (TC <   20) cp = 0
#else
    ! DOI: 10.1007/978-94-007-7587-9_3
    cp = 9.1159d-3 + 6.9402d-9*T**2 + 13.0466d5 /T**3
    if (T < 273)   cp = 0.d0
    if (T > 3300)  cp = 0.d0
#endif
  
  end function dspec_heat_capacity_dT


  ! --- For W
  pure function heat_conductivity(T) result(K)

    implicit none
    
    ! Input parameters
    real(8), intent(in) :: T        ! Temperature in K
    real(8) :: K, TC   !  [W / (m K)]

#ifdef BE_WALL
    TC = T - 273.15
    K = 2.472945d-10*TC**4 - 7.354535d-7*TC**3 + 7.959136d-4*TC**2 - 4.470410d-1*TC + 2.073906d2 
    K = min(K,200.d0) !20 C
    K = max(K,60.d0)  !1283 C
    if (TC > 1283) K = 84.7 
#else
    ! DOI: 10.1007/978-94-00 7-7587-9_3
    K = 108.34 - 1.052d-2*T + 23419.9/T
    K = min(K,173.d0) !340
    K = max(K,77.d0)  !3600
#endif

    !K = 173.d0

  end function heat_conductivity


  pure function dheat_conductivity_dT(T) result(K)

    implicit none
    
    ! Input parameters
    real(8), intent(in) :: T        ! Temperature in K
    real(8) :: K, TC   ! [W / (m K)]

#ifdef BE_WALL
    TC = T - 273.15
    K = 9.89178d-10*TC**3 - 22.063605d-7*TC**2 + 15.918272d-4*TC - 4.470410d-1
    if (TC > 1283) K = 0
    if (TC <   20) K = 0
#else
    ! DOI: 10.1007/978-94-00 7-7587-9_3
    K = -1.052d-2 - 23419.9/T**2
    if (T < 340)   K = -0.2131142d0
    if (T > 3600)  K =-0.01232709d0
#endif

    !K = 0.d0

  end function dheat_conductivity_dT


end module heat_diffusion



program calculate_T

  use mod_utils
  use heat_diffusion

  implicit none

  ! Declarations
  integer :: n_nodes, n_tri
  integer :: i_begin, i_end, i_step, i_jump_steps, i_tri, n_qperp
  character(len=64) :: file_name, wall_f_name
  real(8), allocatable :: nodes_xyz(:,:)
  integer, allocatable :: indices(:,:), ind_qperp(:)
  real(8), allocatable :: Jperp(:), l_part(:), q_heat_perp_3d(:), field_wall_angle(:), &
                          T_tri(:), T_tri_hist(:), T_max(:), dT_tri(:), dT_tri_max(:), T_tri_tmp(:)

  ! --- Parameters for heat diffusion
  integer, parameter :: nx=120
  integer            :: nt, i_time, ierr, i
  real*8,  parameter :: L=0.012d0
  real*8,  parameter :: stab_param = 0.01d0, T_init=473.d0, alpha_max = 6.75d-5
  real(8) ::  dx, dt, t_interval, time_now, time_before
  real(8), allocatable :: T_curr(:,:), q_perp(:)

  call read_namelist(i_begin, i_end, i_jump_steps, wall_f_name)

  dx    = L / real(nx - 1)
  write(*,*) ' heat diffusivity = ', alpha_max, ' m2/s'
  write(*,*) ' Radial grid width for T = ', dx

  ! -- Initialize parameters for temperature time evolution
  ! dTdx = 1.d9 / heat_conduct
  ! dt = stab_param * dx**2 / alpha
  ! nt = int(1d-3/dt) 
  ! allocate(T_tri_hist(nx))
  ! write(*,*) 'nt = ', nt
  ! write(*,*) 'dt = ', dt
  ! T_tri_hist = 0.d0
  
  ! do i_time=1, nt
  !   T_tri_hist(1) = T_tri_hist(1) + 2.d0* stab_param * (T_tri_hist(2) - T_tri_hist(1) + dTdx*dx) 
  !   do i = 2, nx - 1
  !     T_tri_hist(i) = T_tri_hist(i) + stab_param * (T_tri_hist(i+1) - 2.0d0 * T_tri_hist(i) + T_tri_hist(i-1))
  !   end do
  ! enddo
  ! write(*,*) T_tri_hist(1) / (2.d0/sqrt(3.1415)*1d9*sqrt(1d-3)/sqrt(heat_conduct*spec_heat_cap*mat_density))
  

  ! Read all vtk files available and check wetted triangles
  write(*,*) 'Read all vtks to calculate number of wetted triangles'
  do i_step = i_begin, i_end, i_jump_steps
    
    write(file_name,'(a,i6.6,a)')   trim(wall_f_name), i_step, '.dat'
    write(*,*) 'Reading ', file_name
    call read_data(file_name, n_nodes, n_tri, nodes_xyz, indices, Jperp, l_part, q_heat_perp_3d, field_wall_angle, time_now, T_tri_tmp, ierr)
    if (ierr .ne. 0) cycle

    if(i_step==i_begin) then
      allocate(ind_qperp(n_tri))
      ind_qperp = 0
    endif

    do i_tri=1, n_tri
      if (q_heat_perp_3d(i_tri)>1.d0) then
        ind_qperp(i_tri) = 1
      endif
    enddo
  end do

  ! -- Count number of wetted triangles
  n_qperp = 0
  do i_tri=1, n_tri
    if (ind_qperp(i_tri)>0) then
      n_qperp = n_qperp + 1
      ind_qperp(i_tri) = n_qperp
    endif
  enddo
  write(*,*) 'Number of total wetted triangles over time = ', n_qperp

  write(*,*) ' '
  write(*,*) ' '

  ! Allocate and define initial temperature
  allocate(T_curr(n_qperp,nx), q_perp(n_qperp))
  allocate(T_tri_hist(nx))
  allocate(dT_tri(n_tri))
  T_curr = T_init
  T_tri  = T_init

  do i_step = i_begin, i_end, i_jump_steps
    
    write(file_name,'(a,i6.6,a)')   trim(wall_f_name), i_step, '.dat'
    write(*,*) 'Reading ', file_name
    call read_data(file_name, n_nodes, n_tri, nodes_xyz, indices, Jperp, l_part, q_heat_perp_3d, field_wall_angle, time_now, T_tri, ierr)
    if (ierr .ne. 0) cycle

    write(*,*) 'At time = ', time_now, ' s'

    if (i_step == i_begin) time_before = time_now

    ! --- Get q_perp for wetted triangles
    do i_tri=1, n_tri
      if (ind_qperp(i_tri)>0) then
        q_perp(ind_qperp(i_tri)) = q_heat_perp_3d(i_tri)*(1.d0-0.22d0)
      endif
    enddo

    t_interval = (time_now-time_before)
    dt = stab_param * dx**2 / alpha_max
    nt = int(t_interval/dt) ! --- number of steps to take

    write(*,*) 't_interval = ', t_interval
    write(*,*) 'ntl = ', nt
    write(*,*) 'dt = ', dt


    if (t_interval < 0) then
      write(*,*) 'Issue with data files, time not monotonic'
      cycle
    endif

    ! --- Advance the temperaure for each triangle
    !$omp parallel shared(T_curr, dx, dt, q_perp, nt, n_qperp) &
    !$omp private(i_tri, i_time, i) firstprivate(T_tri_hist)
    !$omp do
    do i_tri=1, n_qperp
      T_tri_hist = T_curr(i_tri,:)
      do i_time=1, nt
        T_tri_hist = heat_diffusion_step(T_tri_hist, nx, dx, dt, q_perp(i_tri)) 
      enddo
      T_curr(i_tri,:) = T_tri_hist
    enddo
    !$omp end do
    !$omp end parallel

    do i_tri=1, n_tri
      if (ind_qperp(i_tri)>0) then
        T_tri(i_tri)  = T_curr(ind_qperp(i_tri),1)
        !dT_tri(i_tri) = T_tri(i_tri) - T_init
      endif
    enddo

    ! if (i_step==i_begin) then
    !   allocate(T_max(n_tri))
    !   allocate(dT_tri_max(n_tri))
    !   T_max = 0.d0
    ! endif

    ! do i=1, n_tri
    !   if (T_tri(i)>T_max(i)) then
    !     T_max(i)      = T_tri(i)
    !     dT_tri_max(i) = dT_tri(i)
    !   endif
    ! enddo

    ! write(file_name,'(a,i5.5,a)')   'full_wall_with_Tmax.', i_step, '.vtk'
    ! call write_vtk(file_name, n_nodes, n_tri, nodes_xyz, indices, Jperp, l_part, &
    !                q_heat_perp_3d, field_wall_angle, time_now, T_tri, dT_tri_max)  
    
    write(file_name,'(a,i6.6,a)')   'files_with_T.', i_step, '.dat'

    call write_formatted_data(file_name, n_nodes, n_tri, nodes_xyz, indices, Jperp, l_part, &
                  q_heat_perp_3d, field_wall_angle, time_now, T_tri) 

    time_before = time_now

  end do

  ! Deallocate arrays after usage
  deallocate(nodes_xyz, indices, Jperp, l_part, q_heat_perp_3d, field_wall_angle, T_tri, dT_tri)

end program calculate_T

