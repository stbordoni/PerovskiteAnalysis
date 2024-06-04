#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1.h>
#include <iostream>
#include "TApplication.h"
#include <TSpectrum.h>

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
    if (std::string(argv[2])=="true")
        displaywave=true;
    else
        displaywave=false;
    
    std::cout << "Running analysis on file " << filename <<std::endl;
    std::cout << "Drawing waveforms " << displaywave << std::endl;

    bool remove_noise = true;
    bool verbose = false; // to have more printout in the various steps of the analysis

    if (remove_noise)
       std::cout << " Remove noise routine  applied : moving average " << std::endl;
    else
        std::cout << " Using raw waveforms" << std::endl;
    
    if (displaywave)
        std::cout << "display waveform is active. Waveform will be displayed one by one  " << displaywave << std::endl;
    

    
    
    //////////////////////////////////////////////////////////////
    // ROOT app and objects to read the data in the Run
    auto *app = new TApplication("myapp", &argc, argv);
    
    
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

    for (Long64_t ievt = 0; ievt < 500; ievt++) {
         //std::cout << "====== EVENT " << ievt << "======" << std::endl;
        //for (Long64_t ievt = 0; ievt < nEntries; ievt++) {
        if ((ievt+1)%100 == 0) 
            std::cout << "====== EVENT " << ievt+1 << "======" << std::endl;
        
        eventTree->GetEntry(ievt);
        
        //if (ievt ==0)
        Event myevent(eventData.event, eventData.channelId, *eventData.dataSamples);
        
        if (myevent.GetRawWaveform()->size()<1024){
            std::cout << "Error on Event " <<  myevent.GetEventId() << ": waveform corrupted (<1024 samples) --> "<< myevent.GetRawWaveform()->size()  << std::endl;
            continue;

        } 
        else{
        
            
            myevent.ComputeMovingAverage(10, verbose);
            myevent.ComputeBaseline(remove_noise);
            myevent.SubtractBaseline(remove_noise);

            myevent.ComputeIntegral();
            myevent.FindMaxAmp();

            // to printout the waveform
            //for (int isample =0; isample<myevent.GetRawWaveform()->size(); isample++ )
            //    std::cout << myevent.GetRawWaveform()->at(isample) << " " ;
            //std::cout << "\n " << std::endl;

        

            //std::cout << "integral " << myevent.integral << " maxAmp  " << myevent.maxAmp << std::endl; 
            h_baseline->Fill(myevent.baseline);
            h_maxAmp->Fill(myevent.maxAmp);
            h_integral->Fill(myevent.integral);

            if (displaywave){
                std::cout << "--- CH " << myevent.GetChannelId() << " ---" << std::endl;
                TH1F* h_waveform = new TH1F("h_waveform","", 1024,0,1023);
                TH1F* h_AvgMeanwaveform = new TH1F("h_AvgMeanwaveform","", 1024,0,1023);
                h_AvgMeanwaveform->SetLineColor(1);
                h_AvgMeanwaveform->SetLineWidth(3);
                

                for(int isample=0; isample<1024; isample++){
                    if (isample>10) {
                    h_waveform->SetBinContent(isample,myevent.GetRawWaveform()->at(isample)); //raw waveform
                    h_AvgMeanwaveform->SetBinContent(isample,myevent.GetAvgMeanWaveform()->at(isample)); //averaged mean waveform
                    }
                    else  {
                    h_waveform->SetBinContent(isample,0); // prevent first bins with strange values to be shown.
                    h_AvgMeanwaveform->SetBinContent(isample,0);
                    //if(w<10) cout << w << ", " << kv.second->Waveform[w] << endl;
                    }
                }
                
                TCanvas* c = new TCanvas("c", "c", 1000, 700);
                c->cd();
                h_waveform->Draw("HIST");
                h_waveform->GetYaxis()->SetRangeUser(0,0.1);
                h_AvgMeanwaveform->Draw("HISTsame");

                //search for peaks in the waveform 
                Int_t np=20;
                Int_t npeaks = TMath::Abs(np);

                // Use TSpectrum to find the peak candidates
                TSpectrum *s = new TSpectrum(2*npeaks);
                
                Int_t nfound = s->Search(h_AvgMeanwaveform, 20, " ", 0.1);
                printf("Found %d candidate peaks to fit\n",nfound);

                //TH1 *hb = s->Background(h_AvgMeanwaveform,20,"same");
                h_waveform->Draw("HISTsame");
                
                // this is to estimate the background but I don't think it's needed here. 
                //if (hb) c->Update();    
                //TCanvas* c = new TCanvas("c", "c", 1000, 700);
                //c->cd();
                
                c->Update();
                c->WaitPrimitive();

                delete h_waveform; 
                delete h_AvgMeanwaveform; 
                delete c;
                
            }
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
