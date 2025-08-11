! get quantities in JET mushrooms
program get_data_files

  use mod_utils

  implicit none

  ! Declarations
  integer :: n_nodes, n_tri, phi_idx, theta_idx
  integer :: i_begin, i_end, i_step, i_jump_steps, i1, i2, i3, i, j
  integer, parameter :: n_tor = 8, n_pol = 2
  character(len=64) :: file_name, wall_f_name
  real(8), allocatable :: nodes_xyz(:,:)
  integer, allocatable :: indices(:,:), ind_qperp(:)
  real(8), allocatable :: Jperp(:), l_part(:), q_heat_perp_3d(:), field_wall_angle(:)

  real(8) :: time_now, tri_x, tri_y, tri_z, tri_R, tri_phi, tri_theta, dt, time_old, E_dep_tot
  real(8) :: R_geo=2.8d0, Z_geo = 0.7d0, phi0=0, phif = 2.d0 * 3.14159265d0
  real(8), dimension(n_pol+1) :: theta_ranges
  real(8), dimension(n_tor+1) :: phi_ranges
  real(8), dimension(n_pol,n_tor) :: I_curr
  real(8), dimension(3):: v21, v31, norm_tria
  real(8), allocatable :: T_max(:), melt_duration(:), E_dep(:), tri_areas(:), T_tri(:)
  integer            :: ierr
  logical :: valid


  ! --- Create phi and theta ranges to separate the different mushrooms (it assumes all of them are together in the same file)
  do i = 1, n_tor+1
    phi_ranges(i) = phi0 + (phif - phi0) * float(i - 1) / float(n_tor)
  end do
  theta_ranges(:) = (/ -2.0, -1.2, 0.0 /)


  call read_namelist(i_begin, i_end, i_jump_steps, wall_f_name)

  open(unit=20, file="I_curr_output.dat", status="unknown", action="write", iostat=ierr)

  do i_step = i_begin, i_end, i_jump_steps
    
    write(file_name,'(a,i6.6,a)')   trim(wall_f_name), i_step, '.dat'
    write(*,*) 'Reading ', file_name
    call read_data(file_name, n_nodes, n_tri, nodes_xyz, indices, Jperp, l_part, q_heat_perp_3d, field_wall_angle, time_now, T_tri, ierr)

    
    nodes_xyz = nodes_xyz * 1.d-3

    if (i_step==i_begin) then
      allocate(T_max(n_tri), melt_duration(n_tri), E_dep(n_tri))
      allocate(tri_areas(n_tri))
      T_max = 0.d0; melt_duration = 0.d0; E_dep = 0.d0
    endif
    
    time_old = time_now
    I_curr = 0.d0

    ! get triangle areas
    do i=1, n_tri
      i1 = indices(1,i) + 1
      i2 = indices(2,i) + 1
      i3 = indices(3,i) + 1
      v21 = nodes_xyz(:,i2) - nodes_xyz(:,i1)
      v31 = nodes_xyz(:,i3) - nodes_xyz(:,i1)

      norm_tria = [v21(2)*v31(3)- v21(3)*v31(2), v21(3)*v31(1)- v21(1)*v31(3), v21(1)*v31(2)- v21(2)*v31(1)]

      tri_areas(i) = norm2(norm_tria)*0.5d0

      tri_x = (nodes_xyz(1,i1) + nodes_xyz(1,i2) + nodes_xyz(1,i3) ) / 3.d0
      tri_y = (nodes_xyz(2,i1) + nodes_xyz(2,i2) + nodes_xyz(2,i3) ) / 3.d0
      tri_z = (nodes_xyz(3,i1) + nodes_xyz(3,i2) + nodes_xyz(3,i3) ) / 3.d0
      tri_R = sqrt(tri_x**2 +tri_y**2)

      tri_phi   = atan2(-tri_y,tri_x)
      if (tri_phi < 0.d0) tri_phi = tri_phi + 2.d0 * 3.14159265d0
      tri_theta = atan2(-(tri_Z-Z_geo),(tri_R-R_geo))   
      

      ! Find to which mush indices this triangle belongs to
      phi_idx = find_bin(tri_phi, phi_ranges, n_tor)
      theta_idx = find_bin(tri_theta, theta_ranges, n_pol)
      valid = phi_idx >= 1 .and. phi_idx <= n_tor .and. &
              theta_idx >= 1 .and. theta_idx <= n_pol

      if (.not. valid) then
        write(*,*) 'Not valid!'
      endif

      if (valid) then
        I_curr(theta_idx, phi_idx) = I_curr(theta_idx, phi_idx) + Jperp(i)*tri_areas(i)
      end if

    enddo
    write(*,*) sum(tri_areas)

    write(20, '(E16.5)', advance='no') time_now

    do i = 1, n_pol
      do j = 1, n_tor
        write(20, '(E16.5)', advance='no') I_curr(i,j)
      !  write(20, '(" ")', advance='no')  ! add space between values
      end do
    end do
    write(20,*)




    dt = time_now - time_old

    do i=1, n_tri
      E_dep(i) = E_dep(i) + dt*q_heat_perp_3d(i)
    enddo

    ! write(*,*) 'Write vtk'
    ! write(file_name,'(a,i5.5,a)')   'final_files', i_step, '.vtk'
    ! call write_vtk(file_name, n_nodes, n_tri, nodes_xyz, indices, Jperp, l_part, &
    !                q_heat_perp_3d, field_wall_angle, time_now, T_tri, T_max)  

    time_old = time_now

  end do

  ! ! --- Write triangles to theta and phi coordinates
  ! do i=1, n_tri

  !   write(78, '(7ES14.6)') tri_phi, tri_theta, T_max(i), melt_duration(i), E_dep(i), tri_areas(i), tri_z

  ! enddo

  ! E_dep_tot = sum(E_dep*tri_areas)
  ! write(*,*) ' E_dep_tot[MJ] = ', E_dep_tot*1e-6

  ! Deallocate arrays after usage
  deallocate(nodes_xyz, indices, Jperp, l_part, q_heat_perp_3d, field_wall_angle)



  contains

  integer function find_bin(x, edges, nbins)
    real(8), intent(in) :: x
    real(8), intent(in) :: edges(:)
    integer, intent(in) :: nbins
    integer :: j
    find_bin = -1
    do j = 1, nbins
       if (x >= edges(j) .and. x < edges(j+1)) then
          find_bin = j
          return
       end if
    end do
    ! Special case if x == last edge
    if (x == edges(nbins+1)) find_bin = nbins
  end function find_bin

end program get_data_files

