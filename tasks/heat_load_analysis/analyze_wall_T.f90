module mod_vtk
  implicit none

  contains


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
    !print *, "Error opening file: ", filename
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

  buffer = 'SCALARS T[K] float'//lf   ; write(ivtk) trim(buffer)
  buffer = 'LOOKUP_TABLE default'//lf          ; write(ivtk) trim(buffer)
  write(ivtk) (real(T_tri(i),4), i=1,n_tri)

  buffer = 'SCALARS T_max[K] float'//lf   ; write(ivtk) trim(buffer)
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
  integer :: n_nodes, n_tri, i, i1,i2, i3, ierr
  integer :: i_begin=5000, i_end=52700, i_step, i_jump_steps=100
  character(len=64)  :: file_name
  real(8), parameter :: R_geo=2.8, T_melt=3700
  real(8), allocatable :: nodes_xyz(:,:)
  integer, allocatable :: indices(:,:)
  real(8), allocatable :: Jperp(:), l_part(:), q_heat_perp_3d(:), field_wall_angle(:), T_tri(:)

  real(8) :: time_now, tri_x, tri_y, tri_z, tri_R, tri_phi, tri_theta, dt, time_old, E_dep_tot
  real(8), dimension(3):: v21, v31, norm_tria
  real(8), allocatable :: T_max(:), melt_duration(:), E_dep(:), tri_areas(:)

  do i_step = i_begin, i_end, i_jump_steps
    
    write(file_name,'(a,i5.5,a)')   'fluxes_with_T_rec.', i_step, '.dat'
    write(*,*) 'Reading ', file_name
    call read_data(file_name, n_nodes, n_tri, nodes_xyz, indices, Jperp, l_part, q_heat_perp_3d, field_wall_angle, time_now, T_tri, ierr)

    
    nodes_xyz = nodes_xyz * 1.d-3
    if (i_step==i_begin) then
      allocate(T_max(n_tri), melt_duration(n_tri), E_dep(n_tri))
      allocate(tri_areas(n_tri))
      T_max = 0.d0; melt_duration = 0.d0; E_dep = 0.d0
      
      time_old = time_now

      ! get triangle areas
      do i=1, n_tri
        i1 = indices(1,i) + 1
        i2 = indices(2,i) + 1
        i3 = indices(3,i) + 1
        v21 = nodes_xyz(:,i2) - nodes_xyz(:,i1)
        v31 = nodes_xyz(:,i3) - nodes_xyz(:,i1)

        norm_tria = [v21(2)*v31(3)- v21(3)*v31(2), v21(3)*v31(1)- v21(1)*v31(3), v21(1)*v31(2)- v21(2)*v31(1)]

        tri_areas(i) = norm2(norm_tria)*0.5d0
      enddo
      write(*,*) sum(tri_areas)
    endif
    dt = time_now - time_old

    do i=1, n_tri
      if (T_tri(i)>T_max(i)) T_max(i) = T_tri(i)
      if (T_tri(i)>T_melt) then
        melt_duration(i) = melt_duration(i) + dt
      endif
      E_dep(i) = E_dep(i) + dt*q_heat_perp_3d(i)
    enddo

    ! write(*,*) 'Write vtk'
    write(file_name,'(a,i5.5,a)')   'final_files', i_step, '.vtk'
    call write_vtk(file_name, n_nodes, n_tri, nodes_xyz, indices, Jperp, l_part, &
                   q_heat_perp_3d, field_wall_angle, time_now, T_tri, T_max)  

    time_old = time_now

  end do

  ! --- Write triangles to theta and phi coordinates
  do i=1, n_tri
    i1 = indices(1,i) + 1
    i2 = indices(2,i) + 1
    i3 = indices(3,i) + 1
    tri_x = (nodes_xyz(1,i1) + nodes_xyz(1,i2) + nodes_xyz(1,i3) ) / 3.d0
    tri_y = (nodes_xyz(2,i1) + nodes_xyz(2,i2) + nodes_xyz(2,i3) ) / 3.d0
    tri_z = (nodes_xyz(3,i1) + nodes_xyz(3,i2) + nodes_xyz(3,i3) ) / 3.d0
    tri_R = sqrt(tri_x**2 +tri_y**2)

    tri_phi   = atan2(-tri_y,tri_x)
    tri_theta = atan2(-(tri_Z-0.0d0),(tri_R-R_geo))

    write(78, '(6ES14.6)') tri_phi, tri_theta, T_max(i), melt_duration(i), E_dep(i), tri_areas(i)

  enddo

  E_dep_tot = sum(E_dep*tri_areas)
  write(*,*) ' E_dep_tot[MJ] = ', E_dep_tot*1e-6

  ! Deallocate arrays after usage
  deallocate(nodes_xyz, indices, Jperp, l_part, q_heat_perp_3d, field_wall_angle)

end program analyze_T_rise

