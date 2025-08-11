! Just transform the vtk files in faster files to work with
program get_data_files

  use mod_utils

  implicit none

  ! Declarations
  integer :: n_nodes, n_tri
  integer :: i_begin, i_end, i_step, i_jump_steps
  character(len=64) :: file_name, wall_f_name
  real(8), allocatable :: nodes_xyz(:,:)
  integer, allocatable :: indices(:,:), ind_qperp(:)
  real(8), allocatable :: Jperp(:), l_part(:), q_heat_perp_3d(:), field_wall_angle(:)

  integer            :: ierr
  real*8,  parameter :: fscale=1.d0
  real(8) ::  time_now

  call read_namelist(i_begin, i_end, i_jump_steps, wall_f_name)

  do i_step = i_begin, i_end, i_jump_steps
    
    write(file_name,'(a,i6.6,a)')   trim(wall_f_name), i_step, '.vtk'
    write(*,*) 'Create data file for ', file_name
  
    call read_vtk(file_name, n_nodes, n_tri, nodes_xyz, indices, Jperp, l_part, q_heat_perp_3d, field_wall_angle, time_now, ierr)
    
    if (ierr .ne. 0) cycle

    time_now = time_now*fscale
    q_heat_perp_3d = q_heat_perp_3d/fscale

    write(file_name,'(a,i6.6,a)')   trim(wall_f_name), i_step, '.dat'

    call write_formatted_data(file_name, n_nodes, n_tri, nodes_xyz, indices, Jperp, l_part, &
                   q_heat_perp_3d, field_wall_angle, time_now, field_wall_angle*0.d0) 

  end do

  ! Deallocate arrays after usage
  deallocate(nodes_xyz, indices, Jperp, l_part, q_heat_perp_3d, field_wall_angle)

end program get_data_files

