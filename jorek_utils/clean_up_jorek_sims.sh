#!/bin/bash

# Define your filename patterns here
#patterns=("*.out" "*.vtk" "grid_*" "STARWALL" "jorek_model*" "jorek2vtk" "jorek2_postproc" "*.ps" "*.h5" "starwall-response.dat")
#patterns=("*.out" "*.vtk" "grid_*" "STARWALL" "jorek_model*" "jorek2vtk" "jorek2_postproc")
patterns=("*.out" "*.vtk" "grid_*" "STARWALL" "jorek_model*" "jorek2vtk" "jorek2_postproc")

echo "You are in: $(pwd)"
echo "Files matching these patterns may be deleted:"
for pattern in "${patterns[@]}"; do
    echo "  - $pattern"
done

matches=()
for pattern in "${patterns[@]}"; do
    while IFS= read -r file; do
        matches+=("$file")
    done < <(find . -type f -name "$pattern")
done

# Show summary and confirm
count=${#matches[@]}
if (( count == 0 )); then
    echo "No matching files found in $(pwd)."
    exit 0
fi

echo "Found $count matching file(s) in directory: $(pwd)"

read -p "Proceed to delete matching files? [y/N]: " confirm
if [[ "$confirm" =~ ^[Yy]$ ]]; then
    for pattern in "${patterns[@]}"; do
        find . -type f -name "$pattern" -delete
    done
    echo "Matching files deleted."
else
    echo "Aborted. No files deleted."
fi

