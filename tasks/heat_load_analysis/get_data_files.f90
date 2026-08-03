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




program get_data_files

  use mod_vtk

  implicit none

  ! Declarations
  integer :: n_nodes, n_tri
  integer :: i_begin=5000, i_end=24000, i_step, i_jump_steps=100
  character(len=64) :: file_name
  real(8), allocatable :: nodes_xyz(:,:)
  integer, allocatable :: indices(:,:), ind_qperp(:)
  real(8), allocatable :: Jperp(:), l_part(:), q_heat_perp_3d(:), field_wall_angle(:)

  integer            :: ierr
  real*8,  parameter :: fscale=1.d0
  real(8) ::  time_now

  do i_step = i_begin, i_end, i_jump_steps
    
    write(file_name,'(a,i5.5,a)')   'wall_loads', i_step, '.vtk'
    write(*,*) 'Create data file for ', file_name
  
    call read_vtk(file_name, n_nodes, n_tri, nodes_xyz, indices, Jperp, l_part, q_heat_perp_3d, field_wall_angle, time_now, ierr)
    
    if (ierr .ne. 0) cycle
    time_now = time_now*fscale
    q_heat_perp_3d = q_heat_perp_3d/fscale

    write(file_name,'(a,i5.5,a)')   'wall_loads', i_step, '.dat'

    call write_formatted_data(file_name, n_nodes, n_tri, nodes_xyz, indices, Jperp, l_part, &
                   q_heat_perp_3d, field_wall_angle, time_now, field_wall_angle*0.d0) 

  end do

  ! Deallocate arrays after usage
  deallocate(nodes_xyz, indices, Jperp, l_part, q_heat_perp_3d, field_wall_angle)

end program get_data_files

