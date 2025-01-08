#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
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
    size_t runEnd = outputana_filename.find("_Data_", runStart);  // End of "Run_*_Cs104-A"

    // Extract the substring
    if (runStart != std::string::npos && runEnd != std::string::npos) {
        Runname = outputana_filename.substr(runStart, runEnd - runStart);
        std::cout << "Extracted part: " << Runname << std::endl;
    } else {
        std::cout << "Unable to extract the desired part of the string." << std::endl;
    }


    //FileOutput = new TFile(fileOut.Data(),"RECREATE");
    FileOutput = new TFile(outputana_filename.c_str(),"RECREATE");
    if(!FileOutput) {std::cout << "no output file! exit." << std::endl; exit(1);}
    //dataOut    = new TTree("Events","Events");
    //storing variables in another tree
    
    TTree *output_tree; 
    //storing variables in another tree
    output_tree = new TTree("output_tree", "Tree contanining all pre-computed quatitities");

    

    double baseline;
    double pulseInt;
    double tailInt;
    int npeaks;
    int distmaxAmppeaks_right;
    int distmaxAmppeaks_left;
    std::vector<double> peakInt;
    std::vector<double> peakAmp;

    output_tree->Branch("baseline", &baseline);
    output_tree->Branch("pulseInt", &pulseInt);   
    output_tree->Branch("tailInt", &tailInt);   
    output_tree->Branch("npeaks", &npeaks);   
    output_tree->Branch("distmaxAmppeaks_right", &distmaxAmppeaks_right);  
    output_tree->Branch("distmaxAmppeaks_left", &distmaxAmppeaks_left);   
     
    output_tree->Branch("peakInt", &peakInt);
    output_tree->Branch("peakAmp", &peakAmp);
    
    //output_tree->Branch("ToT", &ToT);
   

    

    // Open the ROOT file
    TFile *file = new TFile(filename.c_str(), "READ");
    
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
    TH1F* h_maxAmp   = new TH1F("h_maxAmp","", 100, 10, 10);
    TH1F* h_integral = new TH1F("h_integral","", 100, 10,10);
    TH1F* h_peakInt  = new TH1F("h_peakInt", "", 100, 10,10);
    TH1F* h_tailInt  = new TH1F("h_tailInt", "", 100, 10,10);
    TH1F* h_goodpeaksperevt = new TH1F("h_goodpeaksperevt", "", 15, -0.5, 14.5);    
    //TH1F* h_distancetomaxAmppeak_left = new TH1F("h_distancetomaxAmppeak_left", "", nSamples, -500, 524 );
    //TH1F* h_distancetomaxAmppeak_right = new TH1F("h_distancetomaxAmppeak_right","", nSamples, -500, 524 );
    TH1F* h_distancetomaxAmppeak_left = new TH1F("h_distancetomaxAmppeak_left", "", nSamples, -200, 824 );
    TH1F* h_distancetomaxAmppeak_right = new TH1F("h_distancetomaxAmppeak_right","", nSamples, -200, 824 );


    std::vector <TH1F*> h_position_peak;
    std::vector <TH1F*> h_distance_peak;

    h_position_peak.reserve(10);
    h_distance_peak.reserve(10);
    for (int i = 0; i<10; i++){
        h_position_peak.emplace_back(new TH1F(Form("h_position_peak_%d", i), Form("position peak %d",i), nSamples/4, 0, 1023 ));
        h_distance_peak.emplace_back(new TH1F(Form("h_distance_peak_%d", i), Form("distance peak %d",i), nSamples/4, -500, 524 ));
    }

    TH1F* h_waveform = new TH1F("h_waveform","raw waveform", nSamples,0,1023);
    TH1F* h_waveforminTime = new TH1F("h_waveforminTime","", nSamples,0*0.3125,1023*0.3125);
    TH1F* h_AvgMeanwaveform = new TH1F("h_AvgMeanwaveform","", nSamples,0,1023);
    h_AvgMeanwaveform->SetLineColor(1);
    h_AvgMeanwaveform->SetLineWidth(3);
    TH1F* h_tailIntegral = new TH1F("h_tailIntegral","", nSamples,0,1023);

    TH1F* h_peakAmp   = new TH1F("h_peakAmp","", 100, 10, 10);
    TH2F* h_peakAmp_vs_peakI = new TH2F("h_peakAmp_vs_peakI","" ,100, 10, 10, 100, 10,10);

    TH1F *h_count_peaks_entries = new TH1F("h_count_peaks_entries", "counting peaks entries", 5, -0.5, 4.5); 
    
     // initiate the FFT
    TVirtualFFT::SetTransform(0);

    // get the magnitude (i.e., power) of the transform of f(x)
    TH1* h_1dMagRaw = NULL;
    TH1* h_noisefreq =NULL;
    //hnoisefreq->SetName("hnoisefreq");

    TH1* h_noisedB =NULL;

    // Loop over entries in the event tree
    Long64_t nEntries = eventTree->GetEntries();
    std::cout << "entries " << nEntries << std::endl;

    
    //for (Long64_t ievt = 2350; ievt < nEntries-1; ievt++) {
    for (Long64_t ievt = 0; ievt < nEntries-1; ievt++) {
        //std::cout << "====== EVENT " << ievt+1 << "======" << std::endl;

        if ((ievt+1)%100 == 0) 
            std::cout << "====== EVENT " << ievt+1 << "======" << std::endl;

        h_waveform->Reset();
        h_AvgMeanwaveform->Reset();
        h_waveforminTime->Reset();
        h_tailIntegral->Reset();

        eventTree->GetEntry(ievt);
        if (ievt>0)
            h_1dMagRaw->Clear();
        

        Event myevent(eventData.event, eventData.channelId, *eventData.dataSamples);
    
       


        if (myevent.GetRawWaveform()->size()<nSamples){
            std::cout << "Error on Event " <<  myevent.GetEventId() << ": waveform corrupted (<1024 samples) --> "<< myevent.GetRawWaveform()->size()  << std::endl;
            continue;
        } 
        else{
        
                   
            myevent.ComputeMovingAverage(20, verbose);
            myevent.ComputeBaseline(mitigate_noise);
            myevent.SubtractBaseline(mitigate_noise);
            
            myevent.ComputeIntegral();
            myevent.FindMaxAmp();



            // to printout the waveform
            //if (verbose){
            //    for (int isample =0; isample<myevent.GetRawWaveform()->size(); isample++ )
            //        std::cout << myevent.GetRawWaveform()->at(isample) << " " ;
            //    std::cout << "\n " << std::endl;
            //}
        
            
            h_baseline->Fill(myevent.baseline);
            h_maxAmp->Fill(myevent.maxAmp);
            h_integral->Fill(myevent.integral);

            
            


            //ToDO: if it is possible to use TSpectrum over the array of values and not the histo 
            //then remove all this histo part to the displaywave block
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

            myevent.ngoodpeaks =0;

            // set value to define a good peak

            double peak_thsld = 0.002; // hard cut on the assumed value of a single p.e of ~3mV
            //  double peak_thsld = 5*fabs(myevent.baseline);


            // Use TSpectrum to find the peak candidates
            TSpectrum *s = new TSpectrum(2*npeaks);
            
            Int_t nfound = s->Search(h_AvgMeanwaveform, 10, "goff", 0.05);
            if (verbose) printf("Found %d candidate peaks to fit\n",nfound);

            Double_t *xpeaks;
            xpeaks = s->GetPositionX();

            for (Int_t p=0;p<nfound;p++) {
                Double_t xp = xpeaks[p];
                if (verbose) std::cout << "p  " << p << " : xp " << xp << std::endl;     

                Int_t bin = h_AvgMeanwaveform->GetXaxis()->FindBin(xp);
                Double_t yp = h_AvgMeanwaveform->GetBinContent(bin);

                //if (yp > 5*fabs(myevent.baseline)) {
                //if (yp > 0.002) { //hard cut on 2mV threshold (assuming 1 p.e. to be of 3mV)
                if (yp > peak_thsld){ 
                    //couont and record the information about the good peaks    
                    myevent.ngoodpeaks++;         
                    myevent.x_peak.push_back(xp);
                    myevent.y_peak.push_back(yp);
                    if (verbose) std::cout << " found a good peak  x: " << xp << " ; y: " << yp << std::endl;

                    // compute integral around the peak
                    double localI=0;
                    
                    localI = myevent.ComputeLocalIntegral(bin, 10 );
                    h_peakInt->Fill(localI);
                    h_peakAmp->Fill(yp);
                    h_peakAmp_vs_peakI->Fill(yp, localI); 
                  
                }
            }

            
            if (verbose) 
            std::cout << "found " << myevent.ngoodpeaks << " good peaks (i.e. x3*baseline) with baseline =" << myevent.baseline << std::endl;
            
            if (! myevent.y_peak.empty()){
        
                int maxElementIndex = std::max_element(myevent.y_peak.begin(),myevent.y_peak.end()) - myevent.y_peak.begin();
                int maxElement = *std::max_element(myevent.y_peak.begin(), myevent.y_peak.end());
                //double tailIntegral = 0;
                if (verbose) 
                    std::cout << "maxElementIndex:" << maxElementIndex << ", maxElement:" << maxElement << '\n';

                std::vector <double> x_peak_right;                    
                int countpeaks_right=0;
                int countpeaks_left=0;
                for (auto ix: myevent.x_peak ){
                    if (ix < myevent.x_peak.at(maxElementIndex) ){
                        h_count_peaks_entries->Fill(1); //if peak is on the left of the max peak
                        countpeaks_left++;
                        h_distancetomaxAmppeak_left->Fill(ix-myevent.x_peak.at(maxElementIndex));  //the histo should give only negative numbers
                        distmaxAmppeaks_left = ix-myevent.x_peak.at(maxElementIndex);
                        //std::cout<< " left peak - distance: " << ix-myevent.x_peak.at(maxElementIndex)<< std::endl;
                        if (verbose) std::cout <<  ix  << " : filling 1 " << std::endl;
                        }
                    else if (ix > myevent.x_peak.at(maxElementIndex) ){
                        h_count_peaks_entries->Fill(2);  //if peak is on the right of the max peak
                        countpeaks_right++;
                        h_distancetomaxAmppeak_right->Fill(ix-myevent.x_peak.at(maxElementIndex));//the histo should give only positive numbers
                        distmaxAmppeaks_right = ix-myevent.x_peak.at(maxElementIndex); 
                        if (countpeaks_right < 10){
                            h_position_peak.at(countpeaks_right)->Fill(ix);
                            h_distance_peak.at(countpeaks_right)->Fill(ix-myevent.x_peak.at(maxElementIndex));
                        }
                        
                        //std::cout<< " right peak - distance: " << ix-myevent.x_peak.at(maxElementIndex) << std::endl;
                        if (verbose) std::cout <<  ix << " : filling 2 "<< std::endl;
                        }
                    else {
                        h_count_peaks_entries->Fill(3);  //if peak is the same as the max peak (should be one per event)
                        if (verbose) std::cout <<  ix  << " : filling 3 "<< std::endl;
                        h_position_peak.at(0)->Fill(ix); // filling the position of the main peak
                        h_distance_peak.at(0)->Fill(ix-myevent.x_peak.at(maxElementIndex));

                        // compute waveform integral starting from the max peak to the end of the waveform
                        for (int isample = ix -25; isample < myevent.GetAvgMeanWaveform()->size(); isample++) {
                            //tailIntegral += myevent.GetAvgMeanWaveform()->at(isample);
                            myevent.tailIntegral += myevent.GetAvgMeanWaveform()->at(isample);
                            h_tailIntegral->SetBinContent(isample,myevent.GetAvgMeanWaveform()->at(isample));
                        }

                        for (int isample = 0; isample < ix -25; isample++)
                            h_tailIntegral->SetBinContent(isample,0);
                    }
                    //countpeaks_right++;
                }

                h_tailInt->Fill(myevent.tailIntegral);
                //tailInt=tailIntegral; 
            }                

            h_goodpeaksperevt->Fill(myevent.ngoodpeaks);
            
            
            //noise study
           
            h_1dMagRaw = h_waveforminTime->FFT(h_1dMagRaw, "MAG R2C M"); // this has units of 1/f_max
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
                
                std::cout << " Event "<< ievt << std::endl;
                
                TCanvas* c = new TCanvas("c", "c", 1000, 700);
                c->cd(1);
                c->SetTitle(Form("Event %d", (int)ievt));

                h_AvgMeanwaveform->Draw("HIST");
                

                // to draw the markers for the found peaks 
                Double_t *xp = new Double_t[50];
                Double_t *yp = new Double_t[50]; 

            
                //for (int ip=0; ip<nfound;ip++) {
                for (int ip=0; ip<myevent.ngoodpeaks;ip++) {
                    xp[ip] = myevent.x_peak[ip];
                    yp[ip] = myevent.y_peak[ip];
                    if (verbose) std::cout << " found a good peak  x: " << xp[ip] << " ; y: " << yp[ip] << std::endl;
                }

                TPolyMarker *pm = (TPolyMarker*) h_AvgMeanwaveform->GetListOfFunctions()->FindObject("TPolyMarker");
                if (pm){
                    h_AvgMeanwaveform->GetListOfFunctions()->Remove(pm);
                    delete pm;
                }

                pm = new TPolyMarker(myevent.ngoodpeaks, xp, yp);
                
                h_AvgMeanwaveform->GetListOfFunctions()->Add(pm);
                pm->SetMarkerStyle(23);
                pm->SetMarkerColor(kRed);
                pm->SetMarkerSize(1.3);
                h_AvgMeanwaveform->Draw("");

                h_waveform->Draw("HISTsame");

                h_tailIntegral->SetLineColor(kRed);
                h_tailIntegral->SetFillStyle(3352);
                h_tailIntegral->SetFillColor(kRed);
                h_tailIntegral->Draw("same");

                TLine *l_baseline = new TLine(0, myevent.baseline, 1023, myevent.baseline);
                l_baseline->SetLineColor(kBlue);
                l_baseline->SetLineStyle(9);
                l_baseline->SetLineWidth(2);

                //TLine *l_peakthrsld = new TLine(0, 5*fabs(myevent.baseline), 1023, 5*fabs(myevent.baseline));
                //TLine *l_peakthrsld = new TLine(0, 0.002, 1023, 0.002);
                TLine *l_peakthrsld = new TLine(0, peak_thsld, 1023, peak_thsld);
                l_peakthrsld->SetLineColor(kBlack);
                l_peakthrsld->SetLineStyle(4);
                l_peakthrsld->SetLineWidth(2);

                l_baseline->Draw("same");
                l_peakthrsld->Draw("same");

                
                
                //TLegend *leg = new TLegend(0.1,0.7,0.48,0.9);
                TLegend *leg = new TLegend();
                leg->AddEntry("h_waveform", "raw waveform ", "l");
                leg->AddEntry("h_AvgMeanwaveform", "waveform (after running average)", "l");
                leg->AddEntry("h_tailIntegral", "tail Integral", "f");
                leg->AddEntry(l_baseline, "baseline", "l");
                leg->AddEntry(l_peakthrsld, "peak threshold", "l");
                leg->Draw("same");

                delete [] xp;
                delete [] yp;

                c->Update();
                c->WaitPrimitive();

                delete c;
                     
            }//end displaywave
    
            /*
            baseline=myevent.baseline;
            pulseInt=myevent.integral;
            npeaks=myevent.ngoodpeaks;
            peakAmp=myevent.y_peak;
            peakInt=myevent.peak_integral;
            tailInt=myevent.tailIntegral;
            
            output_tree->Fill();
            */
        }   

    
    baseline=myevent.baseline;
    pulseInt=myevent.integral;
    npeaks=myevent.ngoodpeaks;
    peakAmp=myevent.y_peak;
    peakInt=myevent.peak_integral;
    tailInt=myevent.tailIntegral;
    
    output_tree->Fill();
    }


   //output path for plots
    TString plotspath="/Users/bordonis/ResearchActivities/PerovskiteAnalysis/PeroAna/Plots/";

    TCanvas* c1 = new TCanvas("c1", "c1", 1200, 1000);
    c1->Divide(3,2);   
    c1->cd(1);
    h_baseline->Draw();
    c1->cd(2);
    h_maxAmp->Draw();
    h_peakAmp->SetLineColor(2);
    h_peakAmp->Draw("");   
    h_maxAmp->Draw("SAME"); 
    c1->cd(3);
    h_integral->Draw();
    c1->cd(4);
    h_peakInt->Draw();
    c1->cd(5);
    h_goodpeaksperevt->Draw();
    c1->cd(6);
    h_peakAmp_vs_peakI->GetXaxis()->SetTitle("peak Amplitude [V]");
    h_peakAmp_vs_peakI->GetYaxis()->SetTitle("peak Integral [V]");
    h_peakAmp_vs_peakI->Draw("COLZ");

    c1->SaveAs(plotspath+Runname+"_spectra_2mVthr.pdf");

    TCanvas* c2 = new TCanvas("c2", "c2", 1000, 700);
    c2->cd();
    //h_noisefreq->GetYaxis()->SetRangeUser(1000,12000);
    h_noisefreq->Draw();
    // Close the file
    //file->Close();

    c2->SaveAs(plotspath+Runname+"_noisefreq.pdf");
    
    TCanvas* c3 = new TCanvas("c3", "c3", 1500, 500);
    c3->Divide(3,1);
    c3->cd(1);
    h_count_peaks_entries->Draw();
    c3->cd(2);
    gPad->SetLogy();
    h_tailInt->Draw();

    c3->cd(3);
    h_distancetomaxAmppeak_right->SetLineColor(kRed);
    h_distancetomaxAmppeak_right->Draw();
    h_distancetomaxAmppeak_left->SetLineColor(kBlue);
    h_distancetomaxAmppeak_left->Draw("same");
    c3->SaveAs(plotspath+Runname+"_secondarypeaks_2mVthr.pdf");
    

    TCanvas* c4 = new TCanvas("c4", "c4", 1500, 500);
    c4->Divide(2,1);
    
    std::vector <Color_t> color = {kBlack, kRed, kOrange, kOrange-3, kOrange-9};

    TLegend *leg_peaks = new TLegend(0.1,0.7,0.48,0.9);
    
    //for (int i=0; i<h_position_peak.size(); i++){
    for (int i=0; i<5; i++){
        h_position_peak.at(i)->SetLineColor(color.at(i));
        h_distance_peak.at(i)->SetLineColor(color.at(i));
        if (i>0) {
            c4->cd(1);
            h_position_peak.at(i)->Draw("same");
            leg_peaks->AddEntry(h_position_peak.at(i), Form("position of %d secondary peak",i), "l");
            c4->cd(2);
            h_distance_peak.at(i)->Draw("same");
        }
        else {
            c4->cd(1);
            gPad->SetLogy();
            h_position_peak.at(i)->Draw();
            leg_peaks->AddEntry("h_position_peak.at(i)", "Position of main peak", "l");
            c4->cd(2);
            gPad->SetLogy();
            h_distance_peak.at(i)->Draw();
        }
    }
    c4->cd(1);
    leg_peaks->Draw("same");

    c4->SaveAs(plotspath+Runname+"_secondarypeaks_distance_2mVthr.pdf");

    //before closing produce an output tree
    FileOutput->cd();
    output_tree->Write();
    FileOutput->Close();

    app->Run();    
    return 0;
}


