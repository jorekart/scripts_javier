module mod_vtk
  implicit none

  contains

  subroutine read_vtk(filename, n_nodes, n_tri, nodes_xyz, indices, Jperp, l_part, & 
                      q_heat_perp_3d, field_wall_angle, time, ierr)
    implicit none
  

    ! Subroutine arguments
    character(len=*), intent(in) :: filename
    integer, intent(out) :: n_nodes, n_tri
    real(8), allocatable, intent(out) :: nodes_xyz(:,:)
    integer, allocatable, intent(out) :: indices(:,:)
    real(8), allocatable, intent(out) :: Jperp(:), l_part(:), q_heat_perp_3d(:), field_wall_angle(:)
    real(8), intent(inout) :: time
    integer, intent(out)   :: ierr

    ! Local variables
    integer :: i, n_int, j
    character(len=256) :: line
    character(len=256) :: keyword

    ! Open the VTK file for reading
    open(unit=10, file=filename, status='old', action='read', IOSTAT=ierr)
    if (ierr .ne. 0) return

    ! Read and ignore the header lines
    do i = 1, 4
      read(10, '(A)') line
    end do

    read(10, '(A)') line
    read(10, '(A)') line
    read(10, '(3ES16.8)') time
    
    ! Read the POINTS line
    read(10, '(A)') line

    ! Manually parse the number of nodes
    read(line, '(A6, I8)') keyword, n_nodes
   
    ! Allocate memory for the points
    allocate(nodes_xyz(3, n_nodes))

    ! Read the points
    do i = 1, n_nodes
      read(10, '(3ES16.8)') nodes_xyz(:, i)
    end do

    ! Read the POLYGONS line
    read(10, '(A)') line

    ! Manually parse the number of triangles and total integers
    read(line, '(A8, I8, I8)') keyword, n_tri, n_int
 
    ! Allocate memory for the triangle indices
    allocate(indices(3, n_tri))

    ! Read the triangle indices
    do i = 1, n_tri
      read(10, '(I8, 3I8)') j, indices(:, i)
    end do
  
    ! Read CELL_DATA header
    read(10, '(A, I8)') line, j
    
    ! Allocate memory for the cell data
    allocate(Jperp(n_tri), l_part(n_tri), q_heat_perp_3d(n_tri), field_wall_angle(n_tri))
    
    ! Skip SCALARS Jperp[A/m2] header
    read(10, '(A)') line
    read(10, '(A)') line
    
    ! Read Jperp values
    do i = 1, n_tri
      read(10, '(3ES16.8)') Jperp(i)
    end do
    
    ! Skip SCALARS L_pre_collision[m] header
    read(10, '(A)') line
    read(10, '(A)') line
    
    ! Read L_pre_collision values
    do i = 1, n_tri
      read(10, '(3ES16.8)') l_part(i)
    end do
    
    ! Skip SCALARS q_perp[W/m2] header
    read(10, '(A)') line
    read(10, '(A)') line
    
    ! Read q_perp values
    do i = 1, n_tri
      read(10, '(3ES16.8)') q_heat_perp_3d(i)
    end do
    
    ! Skip SCALARS B_wall_angle[deg] header
    read(10, '(A)') line
    read(10, '(A)') line
    
    ! Read B_wall_angle values
    do i = 1, n_tri
      read(10, '(3ES16.8)') field_wall_angle(i)
    end do
    
    ! Close the file
    close(10)
    
  end subroutine read_vtk



  subroutine write_vtk(filename, n_nodes, n_tri, nodes_xyz, indices, Jperp, l_part, &
                q_heat_perp_3d, field_wall_angle, time, T_tri)

  implicit none

  ! Subroutine arguments
  character(len=*), intent(in) :: filename
  integer, intent(in) :: n_nodes, n_tri
  real(8), intent(in) :: nodes_xyz(:,:)
  integer, intent(in) :: indices(:,:)
  real(8), intent(in) :: Jperp(:), l_part(:), q_heat_perp_3d(:), field_wall_angle(:), T_tri(:)
  real(8), intent(in) :: time

  ! Local variables
  integer :: i, j
  character(80) :: buffer
  character(12) :: str1, str2
  character(1)  :: lf
  integer, parameter :: ivtk = 22 ! Unit number for the VTK output file

  lf = char(10) ! Line feed character

  open(unit=ivtk,file=filename,form='unformatted',status="unknown",access='stream',convert='BIG_ENDIAN')

  ! Write the header
  buffer = '# vtk DataFile Version 3.0'//lf     ; write(ivtk) trim(buffer)
  buffer = 'vtk output'//lf                     ; write(ivtk) trim(buffer)
  buffer = 'BINARY'//lf                         ; write(ivtk) trim(buffer)
  buffer = 'DATASET POLYDATA'//lf               ; write(ivtk) trim(buffer)

  buffer = 'FIELD FieldData 1'//lf              ; write(ivtk) trim(buffer)
  buffer = 'TIME 1 1 float'//lf                 ; write(ivtk) trim(buffer)
  write(ivtk) real(time,4)

  ! Write the points
  write(str1, '(I8)') n_nodes
  buffer = 'POINTS '//str1//' float'//lf        ; write(ivtk) trim(buffer)
  write(ivtk) ((real(nodes_xyz(j,i),4), j=1,3), i=1,n_nodes)

  ! Write the polygons
  write(str1, '(I8)') n_tri
  write(str2, '(I8)') n_tri * 4
  buffer = 'POLYGONS '//str1//' '//str2//lf     ; write(ivtk) trim(buffer)
  write(ivtk) (int(3,4), (int(indices(j,i),4), j=1,3), i=1,n_tri)

  ! Write the cell data
  buffer = 'CELL_DATA '//str1//lf              ; write(ivtk) trim(buffer)
 
  buffer = 'SCALARS Jperp[A/m2] float'//lf     ; write(ivtk) trim(buffer)
  buffer = 'LOOKUP_TABLE default'//lf          ; write(ivtk) trim(buffer)
  write(ivtk) (real(Jperp(i),4), i=1,n_tri)

  buffer = 'SCALARS L_pre_collision[m] float'//lf ; write(ivtk) trim(buffer)
  buffer = 'LOOKUP_TABLE default'//lf          ; write(ivtk) trim(buffer)
  write(ivtk) (real(l_part(i),4), i=1,n_tri)

  buffer = 'SCALARS q_perp[W/m2] float'//lf    ; write(ivtk) trim(buffer)
  buffer = 'LOOKUP_TABLE default'//lf          ; write(ivtk) trim(buffer)
  write(ivtk) (real(q_heat_perp_3d(i),4), i=1,n_tri)

  buffer = 'SCALARS B_wall_angle[deg] float'//lf ; write(ivtk) trim(buffer)
  buffer = 'LOOKUP_TABLE default'//lf          ; write(ivtk) trim(buffer)
  write(ivtk) (real(field_wall_angle(i),4), i=1,n_tri)

  buffer = 'SCALARS DeltaT_tot[K] float'//lf   ; write(ivtk) trim(buffer)
  buffer = 'LOOKUP_TABLE default'//lf          ; write(ivtk) trim(buffer)
  write(ivtk) (real(T_tri(i),4), i=1,n_tri)

  ! Close the file
  close(ivtk)

end subroutine write_vtk




subroutine write_formatted_data(filename, n_nodes, n_tri, nodes_xyz, indices, Jperp, l_part, &
  q_heat_perp_3d, field_wall_angle, time, T_tri)

implicit none

! Subroutine arguments
character(len=*), intent(in) :: filename
integer, intent(in) :: n_nodes, n_tri
real(8), intent(in) :: nodes_xyz(:,:)
integer, intent(in) :: indices(:,:)
real(8), intent(in) :: Jperp(:), l_part(:), q_heat_perp_3d(:), field_wall_angle(:), T_tri(:)
real(8), intent(in) :: time

! Local variables
integer :: i, j
integer, parameter :: ivtk = 22 ! Unit number for the VTK output file

open(unit=ivtk,file=filename,form='unformatted',status='replace', action='write' )

write(ivtk) time
write(ivtk) n_nodes
write(ivtk) ((nodes_xyz(j,i), j=1,3), i=1,n_nodes)
write(ivtk) n_tri
write(ivtk) ((indices(j,i), j=1,3), i=1,n_tri)
write(ivtk) (Jperp(i), i=1,n_tri)
write(ivtk) (l_part(i), i=1,n_tri)
write(ivtk) (q_heat_perp_3d(i), i=1,n_tri)
write(ivtk) (field_wall_angle(i), i=1,n_tri)
write(ivtk) (T_tri(i), i=1,n_tri)

! Close the file
close(ivtk)

end subroutine write_formatted_data



end module mod_vtk


module heat_diffusion
  implicit none
  contains

  pure function heat_diffusion_step(T_curr, nx, dx, dt, alpha, dTdx, T_left) result(T_next)
    ! Input parameters
    integer, intent(in) :: nx        ! Number of spatial grid points
    real(8), intent(in) :: dx, dt    ! Spatial and temporal step sizes
    real(8), intent(in) :: alpha     ! Thermal diffusivity
    real(8), intent(in) :: dTdx      ! Temperature gradtient at x=0
    real(8), intent(in) :: T_left    ! Dirichlet boundary condition at x=L
    real(8), intent(in), dimension(nx) :: T_curr  ! Current temperature distribution

    ! Output
    real(8), dimension(nx) :: T_next ! Temperature distribution after one time step

    ! Local variables
    real(8) :: r
    integer :: i

    ! Calculate stability parameter
    r = alpha * dt / (dx * dx)

    ! Boundary conditions
    !T_next(1) = T_curr(1) + r * (T_curr(2) - T_curr(1)) + dTdx * dt / dx  ! Neumann BC: heat flux at x=0
    T_next(1) = T_curr(1) + 2.d0* r * (T_curr(2) - T_curr(1) + dTdx*dx)  ! Neumann BC: heat flux at x=0
    T_next(nx) = T_left                                               ! Dirichlet BC: fixed temperature at x=L

    ! Update temperature for internal points
    do i = 2, nx - 1
        T_next(i) = T_curr(i) + r * (T_curr(i+1) - 2.0 * T_curr(i) + T_curr(i-1))
    end do

  end function heat_diffusion_step
end module heat_diffusion



program calculate_T_rise_in_3D_wall

  use mod_vtk
  use heat_diffusion

  implicit none

  ! Declarations
  integer :: n_nodes, n_tri
  integer :: i_begin=1180, i_end=2140, i_step, i_jump_steps=20, i_tri, n_qperp
  character(len=64) :: file_name
  real(8), allocatable :: nodes_xyz(:,:)
  integer, allocatable :: indices(:,:), ind_qperp(:)
  real(8), allocatable :: Jperp(:), l_part(:), q_heat_perp_3d(:), field_wall_angle(:), T_tri(:), T_tri_hist(:)

  ! --- Parameters for heat diffusion
  integer, parameter :: nx=100
  integer            :: nt, i_time, i, ierr
  real*8,  parameter :: L=0.01d0, spec_heat_cap =133.d0 , mat_density = 19.254d3, heat_conduct=173.d0, fscale=60.d0
  real*8,  parameter :: stab_param = 0.1d0, T_left=0.d0!373.15
  real(8) ::  dx, dt, alpha, t_interval, time_now, time_before, dTdx
  real(8), allocatable :: T_curr(:,:), q_perp(:)

  alpha = heat_conduct / (spec_heat_cap*mat_density)
  dx    = L / real(nx - 1)
  write(*,*) ' heat diffusivity = ', alpha, ' m2/s'
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
    
    write(file_name,'(a,i5.5,a)')   'full_wall_corrected', i_step, '.vtk'
    write(*,*) 'Reading ', file_name
    call read_vtk(file_name, n_nodes, n_tri, nodes_xyz, indices, Jperp, l_part, q_heat_perp_3d, field_wall_angle, time_now, ierr)
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
  allocate(T_curr(n_qperp,nx), q_perp(n_qperp), T_tri(n_tri))
  allocate(T_tri_hist(nx))

  do i_step = i_begin, i_end, i_jump_steps
    
    write(file_name,'(a,i5.5,a)')   'full_wall_corrected', i_step, '.vtk'
    write(*,*) 'Calculate T rise for ', file_name
  
    call read_vtk(file_name, n_nodes, n_tri, nodes_xyz, indices, Jperp, l_part, q_heat_perp_3d, field_wall_angle, time_now, ierr)
    if (ierr .ne. 0) cycle
    time_now = time_now*fscale
    q_heat_perp_3d = q_heat_perp_3d/fscale

    write(*,*) 'At time = ', time_now, ' s'

    if (i_step == i_begin) time_before = time_now

    ! --- Get q_perp for wetted triangles
    do i_tri=1, n_tri
      if (ind_qperp(i_tri)>0) then
        q_perp(ind_qperp(i_tri)) = q_heat_perp_3d(i_tri) 
      endif
    enddo

    t_interval = (time_now-time_before)
    dt = stab_param * dx**2 / alpha
    nt = int(t_interval/dt) ! --- number of steps to take

    write(*,*) 't_interval = ', t_interval
    write(*,*) 'ntl = ', nt
    write(*,*) 'dt = ', dt


    if (t_interval < 0) then
      write(*,*) 'Issue with restart files, time not monotonic'
      cycle
    endif


    T_tri_hist = 0.d0

    ! --- Advance the temperaure for each triangle
    !$omp parallel shared(T_curr, dx, dt, alpha, q_perp, nt, n_qperp) &
    !$omp private(i_tri, i_time, dTdx, i) firstprivate(T_tri_hist)
    !$omp do
    do i_tri=1, n_qperp
      T_tri_hist = T_curr(i_tri,:)
      dTdx       = q_perp(i_tri) / heat_conduct
      do i_time=1, nt
        !T_tri_hist = heat_diffusion_step(T_tri_hist, nx, dx, dt, alpha, dTdx, T_left) 
        T_tri_hist(1) = T_tri_hist(1) + 2.d0* stab_param * (T_tri_hist(2) - T_tri_hist(1) + dTdx*dx) 
        do i = 2, nx - 1
          T_tri_hist(i) = T_tri_hist(i) + stab_param * (T_tri_hist(i+1) - 2.0d0 * T_tri_hist(i) + T_tri_hist(i-1))
        end do
      enddo
      T_curr(i_tri,:) = T_tri_hist
    enddo
    !$omp end do
    !$omp end parallel

    do i_tri=1, n_tri
      if (ind_qperp(i_tri)>0) then
        T_tri(i_tri) = T_curr(ind_qperp(i_tri),1)
      endif
    enddo

    ! write(*,*) 'Write vtk'
    ! write(file_name,'(a,i5.5,a)')   'full_wall_with_Trise', i_step, '.vtk'

    ! call write_vtk(file_name, n_nodes, n_tri, nodes_xyz, indices, Jperp, l_part, &
    !                q_heat_perp_3d, field_wall_angle, time_now, T_tri)  

    write(file_name,'(a,i5.5,a)')   'full_wall_with_Trise', i_step, '.dat'

    call write_formatted_data(file_name, n_nodes, n_tri, nodes_xyz, indices, Jperp, l_part, &
                   q_heat_perp_3d, field_wall_angle, time_now, T_tri) 

    time_before = time_now

  end do

  ! Deallocate arrays after usage
  deallocate(nodes_xyz, indices, Jperp, l_part, q_heat_perp_3d, field_wall_angle)

end program calculate_T_rise_in_3D_wall

