#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1.h>
#include <iostream>
#include "TApplication.h"

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




//void analyzeData(const char* filename, bool displaywave) {
int main(int argc, char *argv[]){
    
    std::cout << "found " << argc << " arguments " << std::endl;
    for (int iarg=0; iarg <argc; iarg++)
        std::cout << "argv["  << iarg << "]" <<argv[iarg] << std::endl; 

    char* filename=argv[1];
    bool displaywave;
    if (argv[2]=="true")
        displaywave=true;
    else
        displaywave=false;
    //char* filename="/Users/bordonis/ResearchActivities/PerovskiteAnalysis/PeroAna/rootinputfiles/Run_BR300_Data_1_24_2024_Ascii_Am003_output.root";
    //char* filename="/Users/bordonis/ResearchActivities/PerovskiteAnalysis/PeroAna/rootinputfiles/Run_000_bulkMAPbBr3_Data_3_15_2024_Ascii_output.root";
    //bool displaywave=false;
    
    std::cout << "Running analysis on file " << filename <<std::endl;
    std::cout << "Drawing waveforms " << displaywave << std::endl;

    //////////////////////////////////////////////////////////////
    // ROOT app and objects to read the data in the Run
    auto *app = new TApplication("myapp", &argc, argv);
    //auto *app = new TApplication();


    std::cout << "display wave " << displaywave << std::endl;

    // Open the ROOT file
    TFile *file = new TFile(filename, "READ");
    
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        return 1;
    }

    
    // Read the header tree
    TTree *headerTree = dynamic_cast<TTree*>(file->Get("headerTree"));
    if (!headerTree) {
        std::cerr << "Error: Header tree not found in file " << filename << std::endl;
        file->Close();
        return 1;
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
        return 1;
    }

    EventData eventData;

    // Attach branches to the event tree, assuming you've set up branches in a similar way  
    eventTree->SetBranchAddress("event", &eventData.event);
    eventTree->SetBranchAddress("channelId", &eventData.channelId);
    eventTree->SetBranchAddress("dataSamples", &eventData.dataSamples );

    // ...

    TH1F* h_baseline = new TH1F("h_baseline","", 250, 10, 10);
    TH1F* h_maxAmp   = new TH1F("h_maxAmp","", 250, 10, 10);
    TH1F* h_integral = new TH1F("h_integral","", 100, 10,10);

    // Loop over entries in the event tree
    Long64_t nEntries = eventTree->GetEntries();
    std::cout << "entries " << nEntries << std::endl;

    for (Long64_t ievt = 0; ievt < nEntries; ievt++) {
         std::cout << "====== EVENT " << ievt << "======" << std::endl;
    //for (Long64_t ievt = 0; ievt < nEntries; ievt++) {
        if ((ievt+1)%100 == 0) 
            std::cout << "====== EVENT " << ievt+1 << "======" << std::endl;
        
        eventTree->GetEntry(ievt);
        
        //if (ievt ==0)
        Event myevent(eventData.event, eventData.channelId, *eventData.dataSamples);
        
        //std::cout << " evt  " << myevent.GetEventId() << "  ch " <<  myevent.GetChannelId() << " size  "<< myevent.GetRawWaveform()->size()<< std::endl;
        
        myevent.ComputeBaseline();
        myevent.SubtractBaseline();

        myevent.ComputeIntegral();
        myevent.FindMaxAmp();

        // to printout the waveform
        //for (int isample =0; isample<myevent.GetRawWaveform()->size(); isample++ )
        //    std::cout << myevent.GetRawWaveform()->at(isample) << " " ;
        //std::cout << "\n " << std::endl;


        //for (int isample =0; isample<myevent.GetWaveform()->size(); isample++ )
        //    std::cout << myevent.GetWaveform()->at(isample) << " " ;
        //std::cout << "\n " << std::endl;

        std::cout << "integral " << myevent.integral << " maxAmp  " << myevent.maxAmp << std::endl; 
        h_baseline->Fill(myevent.baseline);
        h_maxAmp->Fill(myevent.maxAmp);
        h_integral->Fill(myevent.integral);

        if (displaywave){
                std::cout << "--- CH " << myevent.GetChannelId() << " ---" << std::endl;
                TH1F* h_waveform = new TH1F("h_waveform","", 1024,0,1023);
                for(int isample=0; isample<1024; isample++){
                    if (isample>10) h_waveform->SetBinContent(isample,myevent.GetRawWaveform()->at(isample));
                    else      h_waveform->SetBinContent(isample,0); // prevent first bins with strange values to be shown.
                    //if(w<10) cout << w << ", " << kv.second->Waveform[w] << endl;
                }
            
                TCanvas* c = new TCanvas("c", "c", 600, 500);
                c->cd();
                h_waveform->Draw("HIST");
                //h_waveform->GetXaxis()->SetRangeUser(0,200);
                c->Update();
                c->WaitPrimitive();

                delete h_waveform; delete c;
            
            }   

    }


    
    TCanvas* c1 = new TCanvas("c1", "c1", 1200, 1000);
    c1->Divide(3,2);   
    c1->cd(1);
    h_baseline->Draw();
    c1->cd(2);
    h_maxAmp->Draw();
    c1->cd(3);
    h_integral->Draw();

    // Close the file
    //file->Close();
    

    app->Run();    
    return 0;
}

//int main() {

//    bool displaywave = true;
    // Replace "output.root" with the actual name of your output file
//    analyzeData("/Users/bordonis/ResearchActivities/PerovskiteAnalysis/PeroAna/rootinputfiles/Run_BR300_Data_1_24_2024_Ascii_Am003_output.root", displaywave);
    
//    return 0;
//}
