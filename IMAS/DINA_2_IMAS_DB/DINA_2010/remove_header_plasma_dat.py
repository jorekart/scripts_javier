import os

def process_plasma_file(file_path):
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()

        # Find the line containing the specific string
        target_string = "Res[mkmOm*m]"
        target_index = next((i for i, line in enumerate(lines) if target_string in line), None)

        if target_index is not None:
            # Keep the target line and all lines after it
            updated_lines = lines[target_index:]

            # Write back the updated content to the file
            with open(file_path, 'w') as file:
                file.writelines(updated_lines)

            print(f"Processed file: {file_path}")
        else:
            print(f"Target string not found in file: {file_path}")
    except Exception as e:
        print(f"Error processing file {file_path}: {e}")

def process_all_plasma_files(folder_path):
    for root, dirs, files in os.walk(folder_path):
        for file_name in files:
            if file_name == "plasma.dat":
                file_path = os.path.join(root, file_name)
                process_plasma_file(file_path)

if __name__ == "__main__":
    folder_path = input("Enter the path to the folder: ").strip()
    if os.path.isdir(folder_path):
        process_all_plasma_files(folder_path)
    else:
        print("The provided path is not a valid directory.")