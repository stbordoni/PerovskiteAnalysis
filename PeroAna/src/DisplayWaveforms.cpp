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



    Long64_t nEntries = eventTree->GetEntries();

    for (Long64_t ievt = 0; ievt < nEntries-1; ievt++) {
        eventTree->GetEntry(ievt);

        if ((ievt+1)%100 == 0) 
            std::cout << "====== EVENT " << ievt+1 << "======" << std::endl;

        size_t nChannels = channelId->size();    


        std::vector<TH1F*> h_waveform(nChannels);
        std::vector<TH1F*> h_waveforminTime(nChannels);
        std::vector<TH1F*> h_AvgMeanwaveform(nChannels);
    
        // loop to declare histos
        for (size_t ich = 0; ich < nChannels; ich++) {
            TString name = Form("h_waveform_evt%d_ch%d", event-1, channelId->at(ich));
            h_waveform[ich] = new TH1F(name, name, nSamples, 0, 1023);

            //name = Form("h_waveforminTime_evt%d_ch%d", event, channelId->at(ich));
            //h_waveforminTime[ich] = new TH1F(name, name, nSamples, 0*0.3125, 1023*0.3125);

            name = Form("h_AvgMeanwaveform_evt%d_ch%d", event-1, channelId->at(ich));
            h_AvgMeanwaveform[ich] = new TH1F(name, name, nSamples, 0, 1023);
            //h_AvgMeanwaveform[ich]->SetLineColor(1);
            //h_AvgMeanwaveform[ich]->SetLineWidth(3);
        }



        for (size_t ich = 0; ich < nChannels; ich++) {
            int chId = channelId->at(ich);
            std::vector<double> waveform = dataSamples->at(ich);

            // Create your Event object
            Event myevent(event, chId, waveform);
            //Event myevent(eventData.event, eventData.channelId, *eventData.dataSamples);
        
            if (myevent.GetRawWaveform()->size()<nSamples){
                std::cout << "Error on Event " <<  myevent.GetEventId() << ": waveform corrupted (<1024 samples) --> "<< myevent.GetRawWaveform()->size()  << std::endl;
                continue;
            } 
            else{

                myevent.ComputeMovingAverage(20, verbose);
                myevent.ComputeBaseline(mitigate_noise);
                myevent.SubtractBaseline(mitigate_noise);
                //myevent.FindMaxAmp();

                for (int isample = 0; isample < nSamples; isample++) {
                    double val = waveform.at(isample);
                    double avg_val = myevent.GetAvgMeanWaveform()->at(isample);

                    if (isample > 10) {
                        h_waveform[ich]->SetBinContent(isample+1, val);
                        //h_waveforminTime[ich]->SetBinContent(isample+1, val);
                        h_AvgMeanwaveform[ich]->SetBinContent(isample+1, avg_val);
                    } else {
                        h_waveform[ich]->SetBinContent(isample+1, 0);
                        //h_waveforminTime[ich]->SetBinContent(isample+1, 0);
                        h_AvgMeanwaveform[ich]->SetBinContent(isample+1, 0);
                    }
                }
                
            }
        }//end loop over channels


        // now drawing for this event 
        // routine to display the waveforms one by one and monitor the peak searches
        if (displaywave){
            
            std::cout << " Event "<< ievt << std::endl; //counting of events in Wavecatcher starts from 1, not 0. 
            
            TCanvas* c = new TCanvas("c", "c", 1000, 700);
            c->cd(1);
            c->SetTitle(Form("Event %d", (int)ievt+1));

            std::vector <Color_t> color = {kBlue, kOrange+5, kOrange-3, kOrange-9};


            for (size_t ich = 0; ich < nChannels; ich++) {
                h_AvgMeanwaveform[ich]->SetLineColor(color.at(ich)); 
                h_AvgMeanwaveform[ich]->SetLineStyle(1);  
                h_AvgMeanwaveform[ich]->SetLineWidth(3);    
                h_waveform[ich]->SetLineColor(color.at(ich));    



                if (ich==0){    
                    h_AvgMeanwaveform[ich]->SetLineStyle(1);        
                    h_AvgMeanwaveform[ich]->SetLineWidth(3);
                    h_AvgMeanwaveform[ich]->SetTitle(Form("Event %d ", event));
                    h_AvgMeanwaveform[ich]->GetXaxis()->SetTitle("Sample number");
                    h_AvgMeanwaveform[ich]->GetYaxis()->SetTitle("Amplitude (V)");
                    h_AvgMeanwaveform[ich]->SetMaximum(h_waveform[ich]->GetMaximum()*1.2);

                    h_waveform[ich]->SetLineWidth(1);
                    
                    
                    h_AvgMeanwaveform[ich]->Draw("HIST");
                    h_waveform[ich]->Draw("HISTsame");
                }
                else{
                    h_AvgMeanwaveform[ich]->SetLineStyle(1);        
                    h_AvgMeanwaveform[ich]->SetLineWidth(3);
                    h_waveform[ich]->SetLineWidth(1);
                    
                    h_AvgMeanwaveform[ich]->Draw("HISTsame");
                    h_waveform[ich]->Draw("HISTsame");
                }
            }//loop over channels

            c->Update();
            c->WaitPrimitive();

            delete c;
                    
        } //displaywave

    }//loop over events

}