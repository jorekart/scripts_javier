module mod_vtk
  implicit none

  contains


  subroutine read_data(filename, n_nodes, n_tri, nodes_xyz, indices, Jperp, l_part, &
                       q_heat_perp_3d, field_wall_angle, time, T_tri)

    implicit none

    ! Subroutine arguments
    character(len=*), intent(in) :: filename
    integer, intent(out) :: n_nodes, n_tri
    real(8), allocatable, intent(out) :: nodes_xyz(:,:)
    integer, allocatable, intent(out) :: indices(:,:)
    real(8), allocatable, intent(out) :: Jperp(:), l_part(:), q_heat_perp_3d(:), field_wall_angle(:), T_tri(:)
    real(8), intent(out) :: time

    ! Local variables
    integer :: i, j, ierr
    integer :: ivtk=22 ! Unit number for the VTK input file


    open(unit=ivtk,file=filename,form='unformatted',status='old', action='read',iostat=ierr )
    if (ierr /= 0) then
      print *, "Error opening file: ", filename
      stop
    end if

    read(ivtk) time
    read(ivtk) n_nodes
    if (allocated(nodes_xyz)) deallocate(nodes_xyz)
    allocate(nodes_xyz(3,n_nodes))
    read(ivtk) ((nodes_xyz(j,i), j=1,3), i=1,n_nodes)

    read(ivtk) n_tri

    if (allocated(indices)) deallocate(indices)
    allocate(indices(3,n_tri))
    read(ivtk) ((indices(j,i), j=1,3), i=1,n_tri)

    if (allocated(Jperp)) deallocate(Jperp)
    allocate(Jperp(n_tri))
    read(ivtk) (Jperp(i), i=1,n_tri)
    
    if (allocated(l_part)) deallocate(l_part)
    allocate(l_part(n_tri))
    read(ivtk) (l_part(i), i=1,n_tri)

    if (allocated(q_heat_perp_3d)) deallocate(q_heat_perp_3d)
    allocate(q_heat_perp_3d(n_tri))
    read(ivtk) (q_heat_perp_3d(i), i=1,n_tri)

    if (allocated(field_wall_angle)) deallocate(field_wall_angle)
    allocate(field_wall_angle(n_tri))
    read(ivtk) (field_wall_angle(i), i=1,n_tri)

    if (allocated(T_tri)) deallocate(T_tri)
    allocate(T_tri(n_tri))
    read(ivtk) (T_tri(i), i=1,n_tri)

    close(ivtk)

  end subroutine read_data




  
  subroutine write_vtk(filename, n_nodes, n_tri, nodes_xyz, indices, Jperp, l_part, &
                q_heat_perp_3d, field_wall_angle, time, T_tri, T_max)

  implicit none

  ! Subroutine arguments
  character(len=*), intent(in) :: filename
  integer, intent(in) :: n_nodes, n_tri
  real(8), intent(in) :: nodes_xyz(:,:)
  integer, intent(in) :: indices(:,:)
  real(8), intent(in) :: Jperp(:), l_part(:), q_heat_perp_3d(:), field_wall_angle(:), T_tri(:), T_max(:)
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

  buffer = 'SCALARS DeltaT_max[K] float'//lf   ; write(ivtk) trim(buffer)
  buffer = 'LOOKUP_TABLE default'//lf          ; write(ivtk) trim(buffer)
  write(ivtk) (real(T_max(i),4), i=1,n_tri)

  ! Close the file
  close(ivtk)

end subroutine write_vtk



end module mod_vtk




program analyze_T_rise

  use mod_vtk

  implicit none

  ! Declarations
  integer :: n_nodes, n_tri
  integer :: i_begin=100, i_end=42000, i_step, i_jump_steps=100
  character(len=64) :: file_name
  real(8), allocatable :: nodes_xyz(:,:)
  integer, allocatable :: indices(:,:)
  real(8), allocatable :: Jperp(:), l_part(:), q_heat_perp_3d(:), field_wall_angle(:), T_tri(:)

  ! --- Parameters for heat diffusion
  integer, parameter :: nx=100
  integer            ::  i
  real*8,  parameter :: L=0.01d0, spec_heat_cap =133.d0 , mat_density = 19.254d3, heat_conduct=173.d0, fscale=60.d0
  real*8,  parameter :: stab_param = 0.1d0, T_left=0.d0!373.15
  real(8) :: time_now, tri_x, tri_y, tri_z, tri_R, tri_phi, tri_theta
  real(8), allocatable :: T_max(:)

  do i_step = i_begin, i_end, i_jump_steps
    
    write(file_name,'(a,i5.5,a)')   'full_wall_with_Trise', i_step, '.dat'
    write(*,*) 'Reading ', file_name
    call read_data(file_name, n_nodes, n_tri, nodes_xyz, indices, Jperp, l_part, q_heat_perp_3d, field_wall_angle, time_now, T_tri)
    
    if (i_step==i_begin) then
      allocate(T_max(n_tri))
      T_max = 0.d0
    endif

    do i=1, n_tri
      if (T_tri(i)>T_max(i)) T_max(i) = T_tri(i)
    enddo


  end do

  ! --- Write triangles to theta and phi coordinates
  do i=1, n_tri
    tri_x = (nodes_xyz(1,indices(1,i)) + nodes_xyz(1,indices(2,i)) + nodes_xyz(1,indices(3,i)) ) / 3.d0
    tri_y = (nodes_xyz(2,indices(1,i)) + nodes_xyz(2,indices(2,i)) + nodes_xyz(2,indices(3,i)) ) / 3.d0
    tri_z = (nodes_xyz(3,indices(1,i)) + nodes_xyz(3,indices(2,i)) + nodes_xyz(3,indices(3,i)) ) / 3.d0
    tri_R = sqrt(tri_x**2 +tri_y**2)

    tri_phi   = atan2(-tri_y,tri_x)
    tri_theta = atan2(-(tri_Z-1.0d0),(tri_R-6.2d3))

    write(78, '(3ES14.6)') tri_phi, tri_theta, T_max(i)

  enddo

    ! write(*,*) 'Write vtk'
  write(file_name,'(a,i5.5,a)')   'files_with_Tmax', i_step, '.vtk'
  call write_vtk(file_name, n_nodes, n_tri, nodes_xyz, indices, Jperp, l_part, &
                 q_heat_perp_3d, field_wall_angle, time_now, T_tri, T_max)  



  ! Deallocate arrays after usage
  deallocate(nodes_xyz, indices, Jperp, l_part, q_heat_perp_3d, field_wall_angle)

end program analyze_T_rise

