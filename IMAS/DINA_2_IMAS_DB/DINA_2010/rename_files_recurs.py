import os

def rename_files_in_folder(folder_path):
    # Define mapping of file name prefixes to their new names
    rename_mapping = {
        "plasma": "plasma.dat",
        "Plasma": "plasma.dat",
        "halo": "halo.dat",
        "Halo": "halo.dat",
        "psi": "psi_data",
        "tor": "toroidal_currents.dat"
    }

    for root, dirs, files in os.walk(folder_path):
        for file_name in files:
            for prefix, new_name in rename_mapping.items():
                if file_name.startswith(prefix):
                    old_path = os.path.join(root, file_name)
                    new_path = os.path.join(root, new_name)

                    # Rename the file
                    try:
                        os.rename(old_path, new_path)
                        print(f"Renamed: {old_path} -> {new_path}")
                    except FileExistsError:
                        print(f"Skipped: {new_path} already exists")
                    except Exception as e:
                        print(f"Error renaming {old_path}: {e}")

if __name__ == "__main__":
    folder_path = input("Enter the path to the folder: ").strip()
    if os.path.isdir(folder_path):
        rename_files_in_folder(folder_path)
    else:
        print("The provided path is not a valid directory.")
