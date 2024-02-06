#include <TFile.h>
#include <TTree.h>
#include <iostream>

#include "Event.h"

struct HeaderInfo {
    // Define the structure for header information
    // Add the necessary members based on your actual header structure

    std::string softwareVersion;
    int tdcValue;

};

struct EventData {
    // Define the structure for event data
    // Add the necessary members based on your actual event data structure

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
    std::vector<double> *dataSamples =0;
    
};

void analyzeData(const char* filename) {
    // Open the ROOT file
    TFile *file = new TFile(filename, "READ");
    
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        return;
    }

    
    // Read the header tree
    TTree *headerTree = dynamic_cast<TTree*>(file->Get("headerTree"));
    if (!headerTree) {
        std::cerr << "Error: Header tree not found in file " << filename << std::endl;
        file->Close();
        return;
    }

    HeaderInfo headerInfo;
    // Attach branches to the header tree, assuming you've set up branches in a similar way
    headerTree->SetBranchAddress("tdcValue", &headerInfo.tdcValue);
    // headerTree->SetBranchAddress("headerInfoMember2", &headerInfo.member2);
    // ...

    // Read the header information from the first entry in the tree
    headerTree->GetEntry(0);

    // Print or use headerInfo as needed

    // Read the event tree
    TTree *eventTree = dynamic_cast<TTree*>(file->Get("eventTree"));
    if (!eventTree) {
        std::cerr << "Error: Event tree not found in file " << filename << std::endl;
        file->Close();
        return;
    }

    EventData eventData;

    // Attach branches to the event tree, assuming you've set up branches in a similar way  
    eventTree->SetBranchAddress("event", &eventData.event);
    eventTree->SetBranchAddress("channelId", &eventData.channelId);
    eventTree->SetBranchAddress("dataSamples", &eventData.dataSamples );

    // ...

    // Loop over entries in the event tree
    Long64_t nEntries = eventTree->GetEntries();
    std::cout << "entries " << nEntries << std::endl;

    for (Long64_t ievt = 0; ievt < 3; ievt++) {
         std::cout << "====== EVENT " << ievt << "======" << std::endl;
    //for (Long64_t ievt = 0; ievt < nEntries; ievt++) {
        if ((ievt+1)%100 == 0) 
            std::cout << "====== EVENT " << ievt+1 << "======" << std::endl;
        
        eventTree->GetEntry(ievt);
        
        //if (ievt ==0)
        Event myevent(eventData.event, eventData.channelId, *eventData.dataSamples);
        
        std::cout << " evt  " << myevent.GetEventId() << "  ch " <<  myevent.GetChannelId() << " size  "<< myevent.GetRawWaveform()->size()<< std::endl;
        
        myevent.ComputeBaseline();
        myevent.SubtractBaseline();


        // to printout the waveform
        for (int isample =0; isample<myevent.GetRawWaveform()->size(); isample++ )
            std::cout << myevent.GetRawWaveform()->at(isample) << " " ;
        std::cout << "\n " << std::endl;


        for (int isample =0; isample<myevent.GetWaveform()->size(); isample++ )
            std::cout << myevent.GetWaveform()->at(isample) << " " ;
        std::cout << "\n " << std::endl;


    }

    // Close the file
    file->Close();
}

int main() {
    // Replace "output.root" with the actual name of your output file
    analyzeData("/Users/bordonis/ResearchActivities/PeroAna/rootinputfiles/Run_BR300_Data_1_24_2024_Ascii_Am003_output.root");
    
    return 0;
}
