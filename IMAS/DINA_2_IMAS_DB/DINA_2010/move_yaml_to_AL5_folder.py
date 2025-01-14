import os
import shutil

run = 1

def move_files_to_folders(src_folder, dest_folder):
    try:

        # List all files in the source folder
        for file_name in os.listdir(src_folder):
            if file_name.endswith(".watcher") or file_name.endswith(".yaml"):
                # Extract the 6-digit folder name from the file name
                folder_number = file_name.split('_')[1][:6]
                
                folder_path = os.path.join(dest_folder, folder_number+"/"+str(run))

                # Move the file into the corresponding folder
                src_file_path = os.path.join(src_folder, file_name)
                dest_file_path = os.path.join(folder_path, file_name)
                shutil.move(src_file_path, dest_file_path)

                print(f"Moved {file_name} to {folder_path}")

        print("Files have been successfully organized.")

    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    src_folder  = "/home/ITER/artolaj/imas_work/DINA_old_2_IMAS/DINA2010"
    dest_folder = "/home/ITER/artolaj/public/imasdb/DINA_2010/4"

    if os.path.isdir(src_folder):
        move_files_to_folders(src_folder, dest_folder)
    else:
        print("The source folder path is not valid.")