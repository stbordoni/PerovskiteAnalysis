#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <input_file>"
    exit 1
fi

# Set the input directory and output directory
#input_dir is the path where the raw data file (from the DAQ) are stored
input_dir="/Users/bordonis/ResearchActivities/Perovskite/data/NewSiPM_KEPboard_DPNCboard/" 
#input_dir="/Users/bordonis/ResearchActivities/Perovskite/data/"

#rootfile_dir is the location where the rootfile built by the readData code is saved 
rootfile_dir="/Users/bordonis/ResearchActivities/PerovskiteAnalysis/PeroAna/rootinputfiles"

# Get the input file name from the command line argument
input_file="$1"

# Extract the base name (without extension) from the input file
base_name=$(basename -- "$input_file")
base_name_no_ext="${base_name%.*}"

# Set the output file name with the output directory
rootfile_path="$rootfile_dir/${base_name_no_ext}_output.root"


# Create build directory if not exists
mkdir -p build
cd build

# Run CMake from the project root
cmake ..

# Build the project
make

# Check if compilation was successful
if [ $? -eq 0 ]; then
    ################################################
    # Execute the code
    cd src/ || exit; # executable is here

    # Run the compiled program with the input and output file names
    ./readData "$input_dir/$input_file" "$rootfile_path"

    # Check if execution was successful
    if [ $? -eq 0 ]; then
        echo "Program executed successfully. Output saved to $rootfile_path"
    else
        echo "Error: Program execution failed."
    fi
else
    echo "Error: Compilation failed."
fi

echo "Currently we are in"
pwd
echo ""

echo "Now executing the analysis code"



if [ $? -eq 0 ]; then
    ################################################
    # Execute the code
    #cd src/ || exit; # executable is here
    
    #here are analysis options
    display_waveforms=true


    # Run the analysis
    #echo $display_waveforms

    ./PerovAna $rootfile_path $display_waveforms
    #./PerovAna
    echo "Analysis code successfully executed. All done! "
else
    echo "Error: Execution failed."
fi