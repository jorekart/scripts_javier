

1) cp ~/scripts_hub/scripts_javier/tasks/loads_3d_walls/job_fluid_loads_3D .
   modify input parameters
   sbatch job_fluid_loads_3D

2) Cp vtk files and namelist to analysis folder

3) Transfor them into .dat files for faster analysis

  cd ~/scripts_hub/scripts_javier/tasks/loads_3d_walls/
  ifort -fpp -DBE_WALL -qopenmp -o get_data_files mod_utils.f90 get_data_files.f90
  cd -
  ~/scripts_hub/scripts_javier/tasks/loads_3d_walls/get_data_files

4) Analyze the data, for example
  ifort -fpp -DBE_WALL -qopenmp -o currents_in_mush mod_utils.f90 currents_in_mush.f90
  ~/scripts_hub/scripts_javier/tasks/loads_3d_walls/currents_in_mush

