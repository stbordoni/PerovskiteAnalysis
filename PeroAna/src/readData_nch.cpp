#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <string>

// Define the structure for header information
struct HeaderInfo {
    std::string softwareVersion;
    int tdcValue;
    int NChannels;
    double samplingPeriod;
    // Add more header information if needed
};

// Define the structure for event data (now a vector of channels) to write the ROOT tree
struct EventTreeData {
    int event;
    std::vector<int> channelId;
    std::vector<int> fcr;
    std::vector<double> baseline;
    std::vector<double> amplitude;
    std::vector<double> charge;
    std::vector<double> leadingEdgeTime;
    std::vector<double> trailingEdgeTime;
    std::vector<double> rateCounter;
    std::vector<std::vector<double>> dataSamples;
};

// Define the structure for channel data (one waveform & measurements per channel)
struct ChannelData {
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

struct Event {
    int eventId;
    std::vector<ChannelData> channels;
}; 


bool readHeaderInfo(std::ifstream& fileStream, HeaderInfo& headerInfo) {
    std::string line;
    while (std::getline(fileStream, line)) {
        //std::cout << line << std::endl;
        if (line.find("=== DATA FILE SAVED WITH SOFTWARE VERSION: V") != std::string::npos) {
	  //  sscanf(line.c_str(), "=== DATA FILE SAVED WITH SOFTWARE VERSION: V%s ===", &headerInfo.softwareVersion);
	  size_t prefixLen = std::string("=== DATA FILE SAVED WITH SOFTWARE VERSION: V").length();
	  size_t suffixPos = line.find(" ===", prefixLen);
	  if (suffixPos != std::string::npos) {
	    headerInfo.softwareVersion = line.substr(prefixLen, suffixPos - prefixLen);
	  }
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



bool readEventData(std::ifstream& fileStream, Event& event) {
    std::string line;

    // Find an event line
    while (std::getline(fileStream, line)) {
        if (line.find("=== EVENT") != std::string::npos) {
            sscanf(line.c_str(), "=== EVENT %d ===", &event.eventId);

            // We now start reading channels for this event
            while (true) {
                std::streampos pos = fileStream.tellg();

                if (!std::getline(fileStream, line)) {
                    // EOF
                    return true;
                }

                // Check if the line is a new event marker
                if (line.find("=== EVENT") != std::string::npos) {
                    // We found next event, rewind and exit
                    fileStream.seekg(pos);
                    return true;
                }

                if (line.find("=== CH:") != std::string::npos) {
                    ChannelData channelData;
                    //channelData.event = event.eventId;

                    sscanf(line.c_str(),
                        "=== CH: %d EVENTID: %*d FCR: %d Baseline: %lf Amplitude: %lf Charge: %lf LeadingEdgeTime: %lf TrailingEdgeTime: %lf RateCounter %lf ===",
                        &channelData.channelId,
                        &channelData.fcr,
                        &channelData.baseline,
                        &channelData.amplitude,
                        &channelData.charge,
                        &channelData.leadingEdgeTime,
                        &channelData.trailingEdgeTime,
                        &channelData.rateCounter);

                    // Next line = data samples
                    if (std::getline(fileStream, line)) {
                        std::istringstream iss(line);
                        double sample;
                        while (iss >> sample) {
                            channelData.dataSamples.push_back(sample);
                        }
                    }

                    event.channels.push_back(channelData);
                }
            }
        }
    }

    return false;
}






int main(int argc, char *argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <inputfile> <outputfile>" << std::endl;
        return 1;
    }

    //const char* inputFileName = argv[1];
    //const char* outputFileName = argv[2];

    std::string inputFileName = argv[1];
    std::string outputFileName = argv[2];

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

    TFile outputFile(outputFileName.c_str(), "RECREATE");

    // Create header tree and fill with header information
    TTree headerTree("headerTree", "Header Information");
    headerTree.Branch("softwareVersion", &headerInfo.softwareVersion);
    headerTree.Branch("NChannels", &headerInfo.NChannels);
    headerTree.Branch("tdcValue", &headerInfo.tdcValue);
    headerTree.Branch("SamplingPeriod", &headerInfo.samplingPeriod);
    headerTree.Fill();

    // Create event tree
    EventTreeData eventTreeData;
    TTree eventTree("eventTree", "Event Information");
    eventTree.Branch("event", &eventTreeData.event);
    eventTree.Branch("channelId", &eventTreeData.channelId);
    eventTree.Branch("fcr", &eventTreeData.fcr);
    eventTree.Branch("baseline", &eventTreeData.baseline);
    eventTree.Branch("amplitude", &eventTreeData.amplitude);
    eventTree.Branch("charge", &eventTreeData.charge);
    eventTree.Branch("leadingEdgeTime", &eventTreeData.leadingEdgeTime);
    eventTree.Branch("trailingEdgeTime", &eventTreeData.trailingEdgeTime);
    eventTree.Branch("rateCounter", &eventTreeData.rateCounter);
    eventTree.Branch("dataSamples", &eventTreeData.dataSamples);


    std::cout << " Reading the input file   "  << std::endl;


    // Read and fill event data

    Event event;
    int eventcounter = 0;

    while (readEventData(inputFile, event)) {
        eventcounter++;

        eventTreeData.event = event.eventId;
       
        eventTreeData.channelId.clear();
        eventTreeData.fcr.clear();
        eventTreeData.baseline.clear();
        eventTreeData.amplitude.clear();
        eventTreeData.charge.clear();
        eventTreeData.leadingEdgeTime.clear();
        eventTreeData.trailingEdgeTime.clear();
        eventTreeData.rateCounter.clear();
        eventTreeData.dataSamples.clear();

        for (const auto& channel : event.channels) {
            eventTreeData.channelId.push_back(channel.channelId);
            eventTreeData.fcr.push_back(channel.fcr);
            eventTreeData.baseline.push_back(channel.baseline);
            eventTreeData.amplitude.push_back(channel.amplitude);
            eventTreeData.charge.push_back(channel.charge);
            eventTreeData.leadingEdgeTime.push_back(channel.leadingEdgeTime);
            eventTreeData.trailingEdgeTime.push_back(channel.trailingEdgeTime);
            eventTreeData.rateCounter.push_back(channel.rateCounter);
            eventTreeData.dataSamples.push_back(channel.dataSamples);
        }

        eventTree.Fill();
        // Clear data samples for the next event
        event.channels.clear();
    }
    std::cout << " Found  " << eventcounter << " events in the input file " << std::endl;

    // Write trees to the output file
    headerTree.Write();
    eventTree.Write();

    // Close the output file
    outputFile.Close();

    return 0;
    
}
