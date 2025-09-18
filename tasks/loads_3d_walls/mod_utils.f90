module mod_utils
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


  subroutine read_data(filename, n_nodes, n_tri, nodes_xyz, indices, Jperp, l_part, &
    q_heat_perp_3d, field_wall_angle, time, T_tri, ierr)

    implicit none

    ! Subroutine arguments
    character(len=*), intent(in) :: filename
    integer, intent(out) :: n_nodes, n_tri, ierr
    real(8), allocatable, intent(out) :: nodes_xyz(:,:)
    integer, allocatable, intent(out) :: indices(:,:)
    real(8), allocatable, intent(out) :: Jperp(:), l_part(:), q_heat_perp_3d(:), field_wall_angle(:), T_tri(:)
    real(8), intent(out) :: time

    ! Local variables
    integer :: i, j
    integer :: ivtk=22 ! Unit number for the VTK input file


    open(unit=ivtk,file=filename,form='unformatted',status='old', action='read',iostat=ierr )
    if (ierr /= 0) then
      print *, "Error opening file: ", filename
    return
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




  subroutine read_namelist(i_begin, i_end, i_jump_steps, file_name, frad, R_geo, Z_geo)
    !! Reads Namelist from given file.
    integer,           intent(inout) :: i_begin, i_end, i_jump_steps
    integer                          :: fu, rc
    logical                          :: file_exists
    real*8, optional,intent(inout)   :: frad             !> radiation fraction (if simulation has no radiation,
                                                         !> but a radiation fraction is assumed in reality)
    real*8, optional,intent(inout)   :: R_geo, Z_geo     !> for theta calculation
    character(len=64),intent(inout)  :: file_name
    character(len=64)                :: wall_f_name

    ! Namelist definition.
    namelist /restart_inputs/ i_begin, i_end, i_jump_steps, wall_f_name, frad, R_geo, Z_geo

    ! Check whether file exists.
    inquire (file='heat_load.nml', exist=file_exists)

    if (.not. file_exists) then
        write(*,*) 'heat_load.nml mamelist file does not exist'
        stop
    end if

    ! Open and read Namelist file.
    open (action='read', file='heat_load.nml', iostat=rc, newunit=fu)
    read (nml=restart_inputs, iostat=rc, unit=fu)
    close(fu)

    write(*,*) wall_f_name

    file_name = wall_f_name

  end subroutine read_namelist



end module mod_utils
