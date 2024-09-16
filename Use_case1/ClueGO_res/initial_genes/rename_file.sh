#!/bin/bash

# Loop through each folder in the current directory
for dir in */; do
    # Get the folder name without the trailing slash
    folder_name=$(basename "$dir")

    # Loop through all files in the current folder
    for file in "$dir"*; do
        # Get the file name without the path
        file_name=$(basename "$file")

        # Check if the file name starts with the folder name
        if [[ "$file_name" != "$folder_name"* ]]; then
            # Rename the file to prepend the folder name
            extension="${file_name##*.}"
            new_name="$folder_name ${file_name#* }"
            mv "$file" "$dir$new_name"
            echo "Renamed '$file_name' to '$new_name' in folder '$folder_name'"
        fi
    done
done
