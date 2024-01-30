#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <input_file>"
    exit 1
fi

# Set the input directory and output directory
input_dir="/Users/bordonis/ResearchActivities/Perovskite/data"
output_dir="/Users/bordonis/ResearchActivities/PeroAna/rootinputfiles"

# Get the input file name from the command line argument
input_file="$1"

# Extract the base name (without extension) from the input file
base_name=$(basename -- "$input_file")
base_name_no_ext="${base_name%.*}"

# Set the output file name with the output directory
output_file="$output_dir/${base_name_no_ext}_output.root"


# Create build directory if not exists
mkdir -p build
cd build

# Run CMake from the project root
cmake ..

# Build the project
make

# Check if compilation was successful
if [ $? -eq 0 ]; then
    # Run the compiled program with the input and output file names
    ./readData "$input_dir/$input_file" "$output_file"

    # Check if execution was successful
    if [ $? -eq 0 ]; then
        echo "Program executed successfully. Output saved to $output_file"
    else
        echo "Error: Program execution failed."
    fi
else
    echo "Error: Compilation failed."
fi
