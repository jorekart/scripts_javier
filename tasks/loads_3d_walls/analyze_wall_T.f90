! get deposited energy and current
program analyze_wall_T

  use mod_utils

  implicit none

  ! Declarations
  integer :: n_nodes, n_tri, i, i1,i2, i3, ierr
  integer :: i_begin, i_end, i_step, i_jump_steps
  real(8), parameter :: R_geo=2.8, T_melt=1556
  character(len=64)    :: file_name, wall_f_name
  real(8), allocatable :: nodes_xyz(:,:)
  integer, allocatable :: indices(:,:)
  real(8), allocatable :: Jperp(:), l_part(:), q_heat_perp_3d(:), field_wall_angle(:), T_tri(:)

  real(8) :: time_now, tri_x, tri_y, tri_z, tri_R, tri_phi, tri_theta, dt, time_old, E_dep_tot
  real(8), dimension(3):: v21, v31, norm_tria
  real(8), allocatable :: T_max(:), melt_duration(:), E_dep(:), tri_areas(:)

  call read_namelist(i_begin, i_end, i_jump_steps, wall_f_name)

  do i_step = i_begin, i_end, i_jump_steps
    
    write(file_name,'(a,i6.6,a)')   trim(wall_f_name), i_step, '.dat'
    write(*,*) 'Reading ', file_name
    call read_data(file_name, n_nodes, n_tri, nodes_xyz, indices, Jperp, l_part, q_heat_perp_3d, field_wall_angle, time_now, T_tri, ierr)
    if (ierr .ne. 0) cycle
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
    ! write(file_name,'(a,i6.6,a)')   trim(wall_f_name), i_step, '.vtk'
    ! call write_vtk(file_name, n_nodes, n_tri, nodes_xyz, indices, Jperp, l_part, &
    !                q_heat_perp_3d, field_wall_angle, time_now, T_tri, T_max)  

    time_old = time_now

  end do

  open(unit=20, file="T_tria.dat", status="unknown", action="write", iostat=ierr)

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

    tri_theta = atan2(-(tri_Z-1.49),(tri_R-2.6)) ! UDP arc center

    write(20, '(6ES14.6)') tri_phi, tri_theta, T_max(i), melt_duration(i), E_dep(i), tri_areas(i)

  enddo

  E_dep_tot = sum(E_dep*tri_areas)
  write(*,*) ' E_dep_tot[MJ] = ', E_dep_tot*1e-6

  ! Deallocate arrays after usage
  deallocate(nodes_xyz, indices, Jperp, l_part, q_heat_perp_3d, field_wall_angle)

  close(20)

end program analyze_wall_T

