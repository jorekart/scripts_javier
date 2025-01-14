import os
from dina_old_2_IMAS import dina_ds_2_imas
database,user_or_path = 'DINA_2010',os.getenv('USER')

def find_all(name, path):
    result = []
    for root, dirs, files in os.walk(path):
        if name in files:
            #result.append(os.path.join(root, name))
            result.append(root)
    return result

# EACH FOLDER WHICH CONTAINS PLASMA.DAT IS ASSUMED TO BE A DINA SIMULATION

folders = find_all('plasma.dat',os.getcwd())

shot = 100499
initial_folder = os.getcwd()
for ifolder in range(len(folders)):
    shot = shot + 1
    run  = 1
    ref_name = folders[ifolder].replace(os.getcwd()+'/','')
    #print('Shot = '+str(shot)+', Run = '+str(run)+', Folder = '+ref_name)
    print(str(shot)+'  '+ref_name)
    #os.chdir(folders[ifolder])
    #dina_ds_2_imas(database,user_or_path,shot,run,ref_name)
    #os.chdir(initial_folder)



    
