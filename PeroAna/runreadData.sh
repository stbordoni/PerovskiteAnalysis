#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <input_file>"
    exit 1
fi

# Set the input directory and output directory
#input_dir is the path where the raw data file (from the DAQ) are stored
#input_dir="/Users/bordonis/ResearchActivities/Perovskite/data/NewSiPM_KEPboard_DPNCboard"
input_dir="/Users/bordonis/ResearchActivities/Perovskite/data/NewSiPM_KEPboard_DPNCboard/CsPbBr3/Cs104-A" 
#input_dir="/Users/bordonis/ResearchActivities/Perovskite/data/NewSiPM_KEPboard_DPNCboard/MAPbBr3/SC67-B" 
#input_dir="/Users/bordonis/ResearchActivities/Perovskite/data/NewSiPM_KEPboard_DPNCboard/CERN/MAPbBr3-SC67_B" 
#input_dir="/eos/home-s/sbordoni/Perovskite/noiseTests/CERN/Cs104-A/"

#input_dir="/Users/bordonis/ResearchActivities/Perovskite/data/PCB_box/UniGe/BulkPerovskite"

#rootfile_dir is the location where the rootfile built by the readData code is saved 
#rootfile_dir="/afs/cern.ch/work/s/sbordoni/PerovskiteReD/PerovskiteAnalysis/PeroAna/rootinputfiles"
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

 # check if cmake was successful
if [ $? -ne 0 ]
then
    echo "  CMake failed. Stopping execution."
    echo ""
    exit 0
fi


# Build the project
make

 # check if make was successful
if [ $? -ne 0 ]
then
    echo "  Make failed. Stopping execution."
    echo ""
    exit 0
fi

echo "  Compilation successful!"
echo ""





# Check if compilation was successful



################################################
# Execute the code
cd src/ || exit; # executable is here


    
# Run the compiled program with the input and output file names
echo "Converting the raw data in to a rootfile. The following command is exected: "
command="./readData $input_dir/$input_file $rootfile_path"
echo $command
echo ""

./readData "$input_dir/$input_file" "$rootfile_path"

echo ""


# Check if execution was successful
if [ $? -eq 0 ]; then
    echo "Program executed successfully. Output saved to $rootfile_path"
else
    echo "Error: Program execution failed."
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
    display_waveforms=false


    # Run the analysis
    #echo $display_waveforms

    ./PerovAna $rootfile_path $display_waveforms
    #./PerovAna
    echo "Analysis code successfully executed. All done! "
else
    echo "Error: Execution failed."
fi
