This is a set of instructions, to put JOREK simulations into IMAS DBs.

1) Go to the parent folder, where the JOREK simulation folders are

2) Prepare each simulation inside each folder.
   - The phase_file.yaml, that must be included in each sim.
     folder, in case multiple input files have been used. There we set the times
     and restart files that you need for each sim. phase (e.g TQ, CQ)

3) Go back to parent folder, and create "folder_list.txt", which is simply
   a list of sim. folder names. These are the simulations to be exported to IMAS.

4) Create your master "imas.nml" file, to indicate IDSs to export, the DB and the run number.
   Also add general info about the cases as a comment

4) Set input parameters in "sim_to_IMASBD.py". 
   IMPORTANT: CHOOSE THE CORRECT SHOT #!!! otherwise you may overwrite old entries
   other input params, are just job_IDS and the executable jorek2_DS

5) Execute sim_to_IMASBD.py in the parent folder

POTENTIAL ISSUES: jorek2_postproc fails sometimes with averaging profiles.
   Using a smaller radial coordinate helps sometimes. To be changed in
   initialise_postproc_settings in communication/IMAS/mod_jorek2IDS.f90
      call set_setting('rad_range_min',   '0.05',  ierr, 'numerical parameter for field line tracing'         )
      call set_setting('rad_range_max',   '0.95', ierr, 'numerical parameter for field line tracing'         )
