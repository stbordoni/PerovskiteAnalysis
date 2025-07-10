#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TLine.h>
#include <TLegend.h>


#include <iostream>
#include "TApplication.h"
#include <TSpectrum.h>
#include <TPolyMarker.h>
#include "TVirtualFFT.h"

#include <filesystem>
#include <string>

#include "Event.h"

struct HeaderInfo {
    // Define the structure for header information
    // Add the necessary members based on your actual header structure

    std::string softwareVersion;
    int tdcValue;

};

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


struct EventData {
    int event;
    std::vector<ChannelData> channels;
};

struct MultiChannelEvent {
    int eventId;
    std::vector<Event> channels;
};


int main(int argc, char *argv[]){
         
    std::cout << "found " << argc << " arguments " << std::endl;
    for (int iarg=0; iarg <argc; iarg++)
        std::cout << "argv["  << iarg << "]" <<argv[iarg] << std::endl; 

    //char* filename=argv[1];
    std::string filename=argv[1];
    bool displaywave;
    if (std::string(argv[2])=="true")
        displaywave=true;
    else
        displaywave=false;
    
    std::cout << "Running analysis on file " << filename <<std::endl;
    std::cout << "Drawing waveforms " << displaywave << std::endl;

    bool mitigate_noise = true;
    bool verbose = false; // to have more printout in the various steps of the analysis

    if (mitigate_noise)
       std::cout << " Remove noise routine  applied : moving average " << std::endl;
    else
        std::cout << " Using raw waveforms" << std::endl;
    
    if (displaywave)
        std::cout << "display waveform is active. Waveform will be displayed one by one  " << displaywave << std::endl;
    

    //////////////////////////////////////////////////////////////
    // ROOT app and objects to read the data in the Run
    auto *app = new TApplication("myapp", &argc, argv);
    
    
    Int_t nSamples = 1024;    


    TFile* FileOutput   ; 
    std::string outputana_filename;
    //create the filename of the output file from the current input filename. Just replacing _output.root to _ana.root 
    std::size_t found = filename.find("_output.root");
    if (found != std::string::npos) {
        // Create the new filename by replacing "_output.root" with "_ana.root"
        outputana_filename = filename.substr(0, found) + "_ana.root";

        // Rename the file
        //fs::rename(entry.path(), entry.path().parent_path() / outputana_filename);

        std::cout << "Renamed: " << filename << " to " << outputana_filename << std::endl;
    }

    // Extract the details of the run to be used to save plots..
    std::string Runname; 
    // Find the last slash
    size_t lastSlash = outputana_filename.find_last_of('/');
    
    // Find the first underscore after "Run_"
    size_t runStart = outputana_filename.find("Run_", lastSlash); // Start of "Run_"
    size_t runEnd = outputana_filename.find("_Ascii_", runStart);  // End of "Run_*_Cs104-A"

    // Extract the substring
    if (runStart != std::string::npos && runEnd != std::string::npos) {
        Runname = outputana_filename.substr(runStart, runEnd - runStart);
        std::cout << "Extracted part: " << Runname << std::endl;
    } else {
        std::cout << "Unable to extract the desired part of the string." << std::endl;
    }

    FileOutput = new TFile(outputana_filename.c_str(),"RECREATE");
    if(!FileOutput) {std::cout << "no output file! exit." << std::endl; exit(1);}

    TTree *output_tree; 
    //storing variables in another tree
    output_tree = new TTree("output_tree", "Tree contanining all pre-computed quatitities");


    // Multi-channel branches (one value per channel)
    std::vector<double> out_baseline;
    std::vector<double> out_pulseInt;
    std::vector<double> out_tailInt;
    std::vector<int> out_npeaks;
    std::vector<int> out_npeaks_specut;
    std::vector<int> out_distmaxAmppeaks_right;
    std::vector<int> out_distmaxAmppeaks_left;
    std::vector<int> out_nphotons_tailInt;

    // Multi-channel branches (vector-of-vectors: one vector per channel)
    std::vector<std::vector<double>> out_peakInt;
    std::vector<std::vector<double>> out_asympeakInt;
    std::vector<std::vector<double>> out_peakAmp;
    //std::vector<std::vector<std::vector<double>>> out_peak_interdistance;


    // Output tree branches
    output_tree->Branch("baseline", &out_baseline);
    output_tree->Branch("pulseInt", &out_pulseInt);
    output_tree->Branch("tailInt", &out_tailInt);
    output_tree->Branch("npeaks", &out_npeaks);
    output_tree->Branch("npeaks_specut", &out_npeaks_specut);
    output_tree->Branch("distmaxAmppeaks_right", &out_distmaxAmppeaks_right);
    output_tree->Branch("distmaxAmppeaks_left", &out_distmaxAmppeaks_left);
    output_tree->Branch("nphotons_tailInt", &out_nphotons_tailInt);
    output_tree->Branch("peakInt", &out_peakInt);
    output_tree->Branch("asympeakInt", &out_asympeakInt);
    output_tree->Branch("peakAmp", &out_peakAmp);
    //output_tree->Branch("peak_interdistance", &out_peak_interdistance);




    
    // Open the ROOT file
    TFile *file = new TFile(filename.c_str(), "READ");
    
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        return 1;
    }

    // Read the event tree
    TTree *eventTree = dynamic_cast<TTree*>(file->Get("eventTree"));
    if (!eventTree) {
        std::cerr << "Error: Event tree not found in file " << filename << std::endl;
        file->Close();
        return 1;
    }


    int event;
    std::vector<int>* channelId = nullptr;
    std::vector<int>* fcr = nullptr;
    std::vector<double>* baseline = nullptr;
    std::vector<double>* amplitude = nullptr;
    std::vector<double>* charge = nullptr;
    std::vector<double>* leadingEdgeTime = nullptr;
    std::vector<double>* trailingEdgeTime = nullptr;
    std::vector<double>* rateCounter = nullptr;
    std::vector<std::vector<double>>* dataSamples = nullptr;


    eventTree->SetBranchAddress("event", &event);
    eventTree->SetBranchAddress("channelId", &channelId);
    eventTree->SetBranchAddress("fcr", &fcr);
    eventTree->SetBranchAddress("baseline", &baseline);
    eventTree->SetBranchAddress("amplitude", &amplitude);
    eventTree->SetBranchAddress("charge", &charge);
    eventTree->SetBranchAddress("leadingEdgeTime", &leadingEdgeTime);
    eventTree->SetBranchAddress("trailingEdgeTime", &trailingEdgeTime);
    eventTree->SetBranchAddress("rateCounter", &rateCounter);
    eventTree->SetBranchAddress("dataSamples", &dataSamples);

    //declare the vector of MultiChannelEvent to store the events
    //std::vector<MultiChannelEvent> eventList;


    Long64_t nEntries = eventTree->GetEntries();

    //for (Long64_t ievt = 0; ievt < nEntries-1; ievt++) {
    for (Long64_t ievt = 0; ievt < 5; ievt++) {
        eventTree->GetEntry(ievt);

        std::cout << "====== EVENT " << ievt+1 << "======" << std::endl;

        MultiChannelEvent multiCh;
        multiCh.eventId = event;

        if ((ievt+1)%100 == 0) 
            std::cout << "====== EVENT " << ievt+1 << "======" << std::endl;

        size_t nChannels = channelId->size();    
       

        std::vector<TH1F*> h_waveform(nChannels);
        std::vector<TH1F*> h_waveforminTime(nChannels);
        std::vector<TH1F*> h_AvgMeanwaveform(nChannels);
    
        for (size_t ich = 0; ich < nChannels; ich++) {
            int chId = channelId->at(ich);
            std::vector<double> waveform = dataSamples->at(ich);

            // Create your Event object
            Event myevent(event, chId, waveform);
            
        
            if (myevent.GetRawWaveform().size()<nSamples){
                std::cout << "Error on Event " <<  myevent.GetEventId() << ": waveform corrupted (<1024 samples) --> "<< myevent.GetRawWaveform().size()  << std::endl;
                continue;
            } 
            else{

                myevent.ComputeMovingAverage(20, verbose);
                myevent.ComputeBaseline(mitigate_noise);
                myevent.SubtractBaseline(mitigate_noise);



            } // end if raw waveform size


           
            multiCh.channels.push_back(myevent);
        } // end loop over channels

        //fill the output tree with the computed quantities
        std::cout << "Filling the variables to be recorded in the tree for all channels" << multiCh.eventId << std::endl;
        for (const auto& chEvt : multiCh.channels) {
            out_baseline.push_back(chEvt.GetBaseline());
            //out_pulseInt.push_back(chEvt.GetPulseIntegral());
            // ... other branches
            //pulseInt.push_back(chEvt.GetPulseIntegral());
            //tailInt.push_back(chEvt.GetTailIntegral());
            //npeaks.push_back(chEvt.GetNumberOfPeaks());
            //npeaks_specut.push_back(chEvt.GetNumberOfPeaksWithSpecut());
            //distmaxAmppeaks_right.push_back(chEvt.GetRightDistanceFromMax());
            //distmaxAmppeaks_left.push_back(chEvt.GetLeftDistanceFromMax());
            //nphotons_tailInt.push_back(chEvt.GetPhotonCountInTail());

            //peakInt.push_back(chEvt.GetPeakIntegrals()); // returns std::vector<double>
            //asympeakInt.push_back(chEvt.GetAsymPeakIntegrals());
            //peakAmp.push_back(chEvt.GetPeakAmplitudes());
            //peak_interdistance.push_back(chEvt.GetPeakInterDistances());
        }
    
        std::cout<< " fill the tree " << std::endl;
        output_tree->Fill();

        std::cout<< " clear the vectors " << std::endl;
        out_baseline.clear();
        out_pulseInt.clear();
        out_tailInt.clear();
        out_npeaks.clear();
        out_npeaks_specut.clear();
        out_distmaxAmppeaks_right.clear();
        out_distmaxAmppeaks_left.clear();
        out_nphotons_tailInt.clear();
        out_peakInt.clear();
        //out_peak_interdistance.clear();
     //eventList.push_back(multiCh);
    } // end loop over events
   

    std::cout << "All events processed. Now writing the output tree to the file." << std::endl;
    // Write the output tree to the file
    if (FileOutput) {
        FileOutput->cd();
        output_tree->Write();
        FileOutput->Close();
    } else {
        std::cerr << "Error: Unable to create output file " << outputana_filename << std::endl;
    }

    std::cout << "Output tree written successfully." << std::endl;

    app->Run(); // Start the ROOT application event loop
    // Clean up
    //delete app;
    
    std::cout << "Analysis completed successfully." << std::endl;

    return 0;

}// end main function