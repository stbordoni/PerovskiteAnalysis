# PerovskiteAnalysis

The script runReadData handle several things:
- read of the raw data file and storage of the information on a rootfile (readData.sh)
- compilation of the analysis code
- execution of the code  PerovAna 

How to run: 
From a terminal you run the line: 
./runreadData.sh filename 

where filename is just the name (no path since this is already given inside the runreadData script)

the script first compile the code and move to the execution of the code. 
In the script there is a variable called "display_waveforms" which is a boolean. This value is taken in consideration when the code executes. 
 - If set to true, the code run event by event and show in a canvas the waveform of that event. This is useful to have an eye-scan of how the waveforms look like (to check noise, peaks, .. )
 - if set to false, the code run over all events and plots some overall distributions (baseline, integral, noise .. )


 NOTE: the script so far does not stop if the compilation fails. It try to execute the code and it fails since in case of failed compilation the executable is not available. TO DO: add a check on the compilation to stop the script.   