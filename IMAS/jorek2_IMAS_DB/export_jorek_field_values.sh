#!/bin/bash

# Check if the correct number of arguments are provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 input_file"
    exit 1
fi

# Input and output files
input_file="$1"
output_file="filtered_file"
output_file2="filtered_2"
output_file3="filtered_3"
output_file4="filtered_4"

# Extract lines containing exactly one "=" and at least one "%"
awk '{
    if (gsub("=", "=") == 1 && index($0, "%") > 0 ) {
        line = $0
        gsub("=>", "-->", line)       # Replace
        gsub("=.*", "", line)       # Replace
        gsub(" ", "", line)        # Remove all spaces
        gsub("ggd_scalar-->", "", line)        # Remove all spaces
        gsub("ggd_vector-->", "", line)        # Remove all spaces
        gsub("grid_sub_ind", "", line)        # Remove all spaces
        gsub("\\([^)]*\\)", "", line)
        gsub("(*)", "[]", line)        # Remove all spaces
        gsub("%", ".", line)       # Replace "%" with "."
       # gsub("\\(.*\\)", "[:]", line)
        print line
    }
}' "$input_file" > "$output_file"

awk '{
    if (index($0, ".") > 0 && index($0,"command_tmp") == 0) {
        line = $0
        print line
    }
}' "$output_file" > "$output_file2"

awk '!seen[$0]++' "$output_file2" > "$output_file3"

echo "**summary**" > "$output_file4"
echo "" >> "$output_file4"
awk '{
    if (index($0, "summary") > 0 ) {
        print $0
    }
}' "$output_file3" >> "$output_file4"

echo "" >> "$output_file4"
echo "**equilibrium**" >> "$output_file4"
echo "" >> "$output_file4"
awk '{
    if (index($0, "equilibrium") > 0 ) {
        print $0
    }
}' "$output_file3" >> "$output_file4"

echo "" >> "$output_file4"
echo "**disruption**" >> "$output_file4"
echo "" >> "$output_file4"
awk '{
    if (index($0, "disruption") > 0 ) {
        print $0
    }
}' "$output_file3" >> "$output_file4"

echo "" >> "$output_file4"
echo "**spi**" >> "$output_file4"
echo "" >> "$output_file4"
awk '{
    if (index($0, "spi") > 0 ) {
        print $0
    }
}' "$output_file3" >> "$output_file4"

echo "" >> "$output_file4"
echo "**pf_passive**" >> "$output_file4"
echo "" >> "$output_file4"
awk '{
    if (index($0, "pf_passive") > 0 ) {
        print $0
    }
}' "$output_file3" >> "$output_file4"

echo "" >> "$output_file4"
echo "**pf_active**" >> "$output_file4"
echo "" >> "$output_file4"
awk '{
    if (index($0, "pf_active") > 0 ) {
        print $0
    }
}' "$output_file3" >> "$output_file4"

echo "" >> "$output_file4"
echo "**wall**" >> "$output_file4"
echo "" >> "$output_file4"
awk '{
    if (index($0, "wall_ids") > 0 ) {
        print $0
    }
}' "$output_file3" >> "$output_file4"


echo "" >> "$output_file4"
echo "**plasma_profiles**" >> "$output_file4"
echo "" >> "$output_file4"
awk '{
    if (index($0, "plasma_profiles") > 0 ) {
        print $0
    }
}' "$output_file3" >> "$output_file4"



# # Process the file and write the cleaned-up content to the output file
# awk '{
#     line = $0                  # Get the current line
#     gsub("=>", "-->", line)       # Replace "%" with "."
#     sub("=.*", "="); # Remove everything after and including the "="
#     gsub("%", ".", line)       # Replace "%" with "."
#     gsub("\\(.*\\)", "[:]", line) # Replace anything inside parentheses with ":"
#     gsub(" ", "", line)        # Remove all spaces
#     print line                 # Print the cleaned line
# }' "$output_file" > "$output_file2"

# # Inform the user
# echo "Cleanup complete. Cleaned output written to $output_file2"
