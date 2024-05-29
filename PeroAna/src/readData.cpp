#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <TFile.h>
#include <TTree.h>

// Define the structure for header information
struct HeaderInfo {
    std::string softwareVersion;
    int tdcValue;
    int NChannels;
    double samplingPeriod;
    // Add more header information if needed
};

// Define the structure for event data
struct EventData {
    int event;
    //int eventId;
    int channelId;
    int fcr;
    double baseline;
    double amplitude;
    double charge;
    double leadingEdgeTime;
    double trailingEdgeTime;
    double rateCounter;
    std::vector<double> dataSamples;
};


bool readHeaderInfo(std::ifstream& fileStream, HeaderInfo& headerInfo) {
    std::string line;
    while (std::getline(fileStream, line)) {
        //std::cout << line << std::endl;
        if (line.find("=== DATA FILE SAVED WITH SOFTWARE VERSION: V") != std::string::npos) {
            sscanf(line.c_str(), "=== DATA FILE SAVED WITH SOFTWARE VERSION: V%s ===", &headerInfo.softwareVersion);
        }
        //if (line.find("=== UnixTime") != std::string::npos) {
        //    sscanf(line.c_str(), "=== UnixTime = %*f date = %*f time = %*s == TDC = %d =", &headerInfo.tdcValue);
        //    return true; // Header information is successfully read
        //}
        if (line.find("=== DATA SAMPLES") != std::string::npos) {
            sscanf(line.c_str(), "=== DATA SAMPLES [%*d] in Volts == NB OF CHANNELS ACQUIRED: %d == Sampling Period: %lf ===",
                   &headerInfo.NChannels, &headerInfo.samplingPeriod);
            return true;
        }
        
    }
    return false; // Header information not found
}




bool readEventData(std::ifstream& fileStream, EventData& eventData) {
    std::string line;
    
    while (std::getline(fileStream, line)) {
        
        if (line.find("=== EVENT") != std::string::npos) {
            sscanf(line.c_str(), "=== EVENT %d ===", &eventData.event);
            

            std::getline(fileStream, line);
            //std::cout << " read this line " << std::endl;
            //std::cout<< line << std::endl;
            // Add code to parse UnixTime if needed

            // Read the line with channelId separately
            if (std::getline(fileStream, line)) {
                sscanf(line.c_str(), "=== CH: %d EVENTID: %*d FCR: %d Baseline: %lf Amplitude: %lf Charge: %lf LeadingEdgeTime: %lf TrailingEdgeTime: %lf RateCounter %lf ===",
                    &eventData.channelId, &eventData.fcr, &eventData.baseline, &eventData.amplitude,
                    &eventData.charge, &eventData.leadingEdgeTime, &eventData.trailingEdgeTime, &eventData.rateCounter);

                //std::cout << "  event " << eventData.event << "  ch " << eventData.channelId  <<   std::endl;

                // Assuming your data samples are on the next line, modify as needed
                if (std::getline(fileStream, line)) {
                    std::istringstream iss(line);
                    double sample;
                    while (iss >> sample) {
                        eventData.dataSamples.push_back(sample);
                    }
                }

                //std::cout <<  eventData.dataSamples.size() << std::endl;
                //if (eventData.event < 4 ){
                //    for (int isample =0; isample < eventData.dataSamples.size(); isample++)
                //        std::cout << eventData.dataSamples.at(isample) << " ";
                
                //std::cout << std::endl;
                //}

                // Successfully read an event, return true
                return true;
            } else {
                // End of file reached
                return false;
            }
        }
        

    }

   
    // No event found
    return false;
}




int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <inputfile> <outputfile>" << std::endl;
        return 1;
    }

    const char* inputFileName = argv[1];
    const char* outputFileName = argv[2];

    std::ifstream inputFile(inputFileName);
    if (!inputFile.is_open()) {
        std::cerr << "Error: Could not open input file " << inputFileName << std::endl;
        return 1;
    }

    HeaderInfo headerInfo;
    if (!readHeaderInfo(inputFile, headerInfo)) {
        std::cerr << "Error: Could not read header information" << std::endl;
        return 1;
    }

    TFile outputFile(outputFileName, "RECREATE");

    // Create header tree and fill with header information
    TTree headerTree("headerTree", "Header Information");
    headerTree.Branch("softwareVersion", &headerInfo.softwareVersion);
    headerTree.Branch("NChannels", &headerInfo.NChannels);
    headerTree.Branch("tdcValue", &headerInfo.tdcValue);
    headerTree.Branch("SamplingPeriod", &headerInfo.samplingPeriod);
    headerTree.Fill();

    // Create event tree
    EventData eventData;
    TTree eventTree("eventTree", "Event Information");
    eventTree.Branch("event", &eventData.event);
    eventTree.Branch("channelId", &eventData.channelId);
    eventTree.Branch("fcr", &eventData.fcr);
    eventTree.Branch("baseline", &eventData.baseline);
    eventTree.Branch("amplitude", &eventData.amplitude);
    eventTree.Branch("charge", &eventData.charge);
    eventTree.Branch("leadingEdgeTime", &eventData.leadingEdgeTime);
    eventTree.Branch("trailingEdgeTime", &eventData.trailingEdgeTime);
    eventTree.Branch("rateCounter", &eventData.rateCounter);
    eventTree.Branch("dataSamples", &eventData.dataSamples);

    std::cout << " Reading the input file   "  << std::endl;
    // Read and fill event data
    int eventcounter = 0;
    while (readEventData(inputFile, eventData)) {
        eventcounter++;

        eventTree.Fill();
        // Clear data samples for the next event
        eventData.dataSamples.clear();
    }
    std::cout << " Found  " << eventcounter << " events in the input file " << std::endl;

    // Write trees to the output file
    headerTree.Write();
    eventTree.Write();

    // Close the output file
    outputFile.Close();

    return 0;
}
