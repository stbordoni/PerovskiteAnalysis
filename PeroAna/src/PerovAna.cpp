#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1.h>
#include <iostream>
#include "TApplication.h"
#include <TSpectrum.h>
#include <TPolyMarker.h>
#include "TVirtualFFT.h"

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
    TH1F* h_peakInt  = new TH1F("h_peakInt", "", 100, 10,10);
    TH1F *h_goodpeaksperevt = new TH1F("h_goodpeaksperevt", "", 15, -0.5, 14.5);    


    TH1F* h_waveform = new TH1F("h_waveform","", nSamples,0,1023);
    TH1F* h_waveforminTime = new TH1F("h_waveforminTime","", nSamples,0*0.3125,1023*0.3125);
    TH1F* h_AvgMeanwaveform = new TH1F("h_AvgMeanwaveform","", nSamples,0,1023);

    std::cout << "  here " << std::endl;
     // initiate the FFT
    TVirtualFFT::SetTransform(0);

    // get the magnitude (i.e., power) of the transform of f(x)
    TH1* h_1dMagRaw = NULL;
    TH1* h_noisefreq =NULL;
    //hnoisefreq->SetName("hnoisefreq");

    

    // Loop over entries in the event tree
    Long64_t nEntries = eventTree->GetEntries();
    std::cout << "entries " << nEntries << std::endl;

    
    //for (Long64_t ievt = 0; ievt < 50; ievt++) {
    for (Long64_t ievt = 0; ievt < nEntries-1; ievt++) {
        //std::cout << "====== EVENT " << ievt+1 << "======" << std::endl;

        if ((ievt+1)%100 == 0) 
            std::cout << "====== EVENT " << ievt+1 << "======" << std::endl;
        
        eventTree->GetEntry(ievt);
        if (ievt>0)
            h_1dMagRaw->Clear();
        

        Event myevent(eventData.event, eventData.channelId, *eventData.dataSamples);
        
        if (myevent.GetRawWaveform()->size()<nSamples){
            std::cout << "Error on Event " <<  myevent.GetEventId() << ": waveform corrupted (<1024 samples) --> "<< myevent.GetRawWaveform()->size()  << std::endl;
            continue;
        } 
        else{
        
            
            myevent.ComputeMovingAverage(50, verbose);
            myevent.ComputeBaseline(mitigate_noise);
            myevent.SubtractBaseline(mitigate_noise);

            myevent.ComputeIntegral();
            myevent.FindMaxAmp();

            // to printout the waveform
            if (verbose){
                for (int isample =0; isample<myevent.GetRawWaveform()->size(); isample++ )
                    std::cout << myevent.GetRawWaveform()->at(isample) << " " ;
                std::cout << "\n " << std::endl;
            }
        
            
            h_baseline->Fill(myevent.baseline);
            h_maxAmp->Fill(myevent.maxAmp);
            h_integral->Fill(myevent.integral);

            
	  


            //ToDO: if it is possible to use TSpectrum over the array of values and not the histo 
            //then remove all this histo part to the displaywave block
            h_AvgMeanwaveform->SetLineColor(1);
            h_AvgMeanwaveform->SetLineWidth(3);
            

            for(int isample=0; isample<nSamples; isample++){
                if (isample>10) {
                    h_waveform->SetBinContent(isample,myevent.GetRawWaveform()->at(isample)); //raw waveform
                    h_waveforminTime->SetBinContent(isample,myevent.GetRawWaveform()->at(isample));
                    h_AvgMeanwaveform->SetBinContent(isample,myevent.GetAvgMeanWaveform()->at(isample)); //averaged mean waveform
                }
                else  {
                    h_waveform->SetBinContent(isample,0); // prevent first bins with strange values to be shown.
                    h_waveforminTime->SetBinContent(isample*0.3125,0);
                    h_AvgMeanwaveform->SetBinContent(isample,0);
              
                }
            }


            
	


            //search for peaks in the waveform 
            Int_t np=20;
            Int_t npeaks = TMath::Abs(np);

            Int_t ngoodpeaks =0;

            // Use TSpectrum to find the peak candidates
            TSpectrum *s = new TSpectrum(2*npeaks);
            
            Int_t nfound = s->Search(h_AvgMeanwaveform, 20, "goff", 0.1);
            if (verbose) printf("Found %d candidate peaks to fit\n",nfound);

            Double_t *xpeaks;
            xpeaks = s->GetPositionX();

            for (Int_t p=0;p<nfound;p++) {
                Double_t xp = xpeaks[p];
                if (verbose) std::cout << "p  " << p << " : xp " << xp << std::endl;     

                Int_t bin = h_AvgMeanwaveform->GetXaxis()->FindBin(xp);
                Double_t yp = h_AvgMeanwaveform->GetBinContent(bin);

                if (yp > 5*myevent.baseline) {
                    //couont and record the information about the good peaks    
                    ngoodpeaks++;
                    myevent.x_peak.push_back(xp);
                    myevent.y_peak.push_back(yp);
                    if (verbose) std::cout << " found a good peak  x: " << xp << " ; y: " << yp << std::endl;

                    // compute integral around the peak
                    double localI=0;
                    
                    localI = myevent.ComputeLocalIntegral(bin, 10 );
                    h_peakInt->Fill(localI);
                }
            }

            if (verbose) 
            std::cout << "found " << ngoodpeaks << " good peaks (i.e. x5*baseline)" << std::endl;

            h_goodpeaksperevt->Fill(ngoodpeaks);

            
            //noise study
           
            h_1dMagRaw = h_waveforminTime->FFT(h_1dMagRaw, "MAG"); // this has units of 1/f_max
            Int_t nbinsx = h_1dMagRaw->GetNbinsX();
            Double_t x_low =0;
            Double_t x_up = h_1dMagRaw->GetXaxis()->GetXmax()/h_waveforminTime->GetXaxis()->GetXmax();
            
            
            TH1D * h_1dMag= new TH1D( "h_1dMag", "Magnitude", nbinsx, x_low, x_up);
            
            // rescale axis to get real units
            for (Int_t bin = 1; bin <= nSamples; bin++){
                h_1dMag->SetBinContent(bin, h_1dMagRaw->GetBinContent(bin)/sqrt(h_1dMagRaw->GetBinContent(bin)));
            }
            
            if (ievt==0){
                h_noisefreq = (TH1D*)h_1dMag->Clone("h_noisefreq");
            }
            else
                h_noisefreq->Add(h_1dMag);

            delete h_1dMag; 
            

            // routine to display the waveforms one by one and monitor the peak searches
            if (displaywave){
                
                
                TCanvas* c = new TCanvas("c", "c", 1000, 700);
                //c->Divide(2,1);
                c->cd(1);
                
                
                h_AvgMeanwaveform->Draw("HIST");
                //h_AvgMeanwaveform->GetYaxis()->SetRangeUser(-0.002,0.03);

                

                
                
                // to draw the markers for the found peaks 
                Double_t xp[100]; //= new Double_t[100];
                Double_t yp[100];// = new Double_t[100]; 

                for (int ip=0; ip<nfound;ip++) {
                    xp[ip] = xpeaks[ip];
                    Int_t bin = h_AvgMeanwaveform->GetXaxis()->FindBin(xp[ip]);
                    yp[ip] = h_AvgMeanwaveform->GetBinContent(bin);
                }

                TPolyMarker *pm = (TPolyMarker*) h_AvgMeanwaveform->GetListOfFunctions()->FindObject("TPolyMarker");
                if (pm){
                    h_AvgMeanwaveform->GetListOfFunctions()->Remove(pm);
                    delete pm;
                    }

                pm = new TPolyMarker(npeaks, xp, yp);
                
                h_AvgMeanwaveform->GetListOfFunctions()->Add(pm);
                pm->SetMarkerStyle(23);
                pm->SetMarkerColor(kRed);
                pm->SetMarkerSize(1.3);
                h_AvgMeanwaveform->Draw("");

                
                
                
                h_waveform->Draw("HISTsame");
                
                // this is to estimate the background but I don't think it's needed here. 
                //TH1 *hb = s->Background(h_AvgMeanwaveform,20,"same");
                //if (hb) c->Update();    
                //TCanvas* c = new TCanvas("c", "c", 1000, 700);
                //c->cd();
                
                

                //delete h_waveform; 
                //delete h_AvgMeanwaveform; 

                

                c->Update();
                c->WaitPrimitive();

                delete c;
                //delete h_1dMagRaw; 
                 
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
    c1->cd(4);
    h_peakInt->Draw();
    c1->cd(5);
    h_goodpeaksperevt->Draw();

    TCanvas* c2 = new TCanvas("c2", "c2", 1000, 700);
    c2->cd();
    h_noisefreq->Draw();
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
