#include "Event.h"
#include <iostream>
#include <algorithm>
#include <TSpectrum.h>
#include "Config.h"




//Event::Event(int evt, int chan, const std::vector<double>& wave)
//    : eventId(evt), channelId(chan), rawWaveform(new std::vector<double>(wave)) 
Event::Event(int evt, int chan, const std::vector<double>& wave)
    : eventId(evt), channelId(chan), rawWaveform(wave)
{
    Init();
    // Additional initialization if needed
}


Event::~Event() {
    //delete rawWaveform; 
    //delete Waveform; 
}

void Event::Init(){

    //std::cout<< " Initialize event " << std::endl;
    baseline = 0;
    integral = 0;
    tailIntegral = 0;
    maxAmp = 0;
    mainPeakIdx = -1; // Initialize to -1 to indicate no peak found yet
    ToT = 0;
    //distmaxAmppeaks_right = 0;
    //distmaxAmppeaks_left = 0;
    //pulseInt = 0;
    //Waveform.clear();
    
    ngoodpeaks=0;
    ngoodpeaks_specut=0;
    peak_integral.clear();
    asympeak_integral.clear();
    peak_interdistance.clear();

    //setEventId();
    //channelId();
}



void Event::SetEventId(int id) {
    eventId = id;
}

void Event::SetChannelId(int ch) {
    channelId = ch;
}



//
int Event::GetEventId() const {
    return eventId;
}

int Event::GetChannelId() const {
    return channelId;
}

//const std::vector<double>* Event::GetRawWaveform() const {
const std::vector<double> Event::GetRawWaveform() const {
    return rawWaveform;
}

//const std::vector<double>* Event::GetAvgMeanWaveform() const {
const std::vector<double> Event::GetAvgMeanWaveform() const {
    return avgWaveform;
}


const std::vector<double> Event::GetWaveform() const {
     return Waveform;
}






void Event::ComputeMovingAverage(int step, bool debug){
   //if (debug) std::cout << "entering the Compute Moving Average" <<std::endl;

    //avgWaveform = new std::vector<double>(*rawWaveform); // create a new WF to recorde the one after moving average
    //avgWaveform = new std::vector<double>(*rawWaveform);
    //avgWaveform = rawWaveform; 
    avgWaveform.resize(rawWaveform.size()); // resize the vector to the size of the raw waveform

    double sum = 0;
    int start = step;
    int end = rawWaveform.size();

    for (int isample=0; isample < end; isample++){
    
        //if (debug) std::cout << "--- isample " << isample << std::endl;
        if ( ((isample - step/2) <0) || ( (isample + step/2) >end )) {
            continue;
        }
        else{
            //if (debug) std::cout << "start:  isample-step/2 : " <<  isample-step/2 <<  " stop : isample+step/2 : " << isample+step/2 << std::endl;
            sum = 0;
            
            for(int isum =isample-step/2; isum<isample+step/2; isum++){ 
                //if (debug) std::cout << "isum " << isum << std::endl;
                //sum += rawWaveform->at(isum);
                sum += rawWaveform.at(isum);
                //if (debug) std::cout << "sum " << sum << std::endl;
            }
            //avgWaveform->at(isample) = (float)sum/step;
            avgWaveform.at(isample) = (float)sum/step;
            //if (debug) std::cout << " isample  " << isample << " rawWaveform->at(isample)  " << rawWaveform->at(isample)  <<  ";  avgWaveform  " << avgWaveform->at(isample)<<std::endl;
        
        }

    }
   

}


void Event::ComputeBaseline(bool _removenoise){

    double sum = 0;
    int start = 5;
    int end = start+50;
    for(int i=start;i<end;i++) {
        if (! _removenoise)
            //sum += rawWaveform->at(i);
            sum += rawWaveform.at(i);
        else 
            //sum += avgWaveform->at(i);
            sum += avgWaveform.at(i);
    }
    baseline = sum/(end-start);
    
}


void Event::SubtractBaseline(bool _removenoise){

    std::vector<double> Waveform;

    if (! _removenoise)
            Waveform = rawWaveform;
        else 
            Waveform = avgWaveform;

        //for(int i=0;i<Waveform->size();i++) {
        for(int i=0;i<Waveform.size();i++) {
            (Waveform)[i] -= baseline;
        }


}



void Event::ComputeIntegral(){
    //integral = 0;
    //for(int i=0;i<Waveform->size();i++) integral += Waveform->at(i);  // for the moment without conversion
    for(int i=0;i<rawWaveform.size();i++) integral += rawWaveform.at(i);  // for the moment without conversion
}




double Event::ComputeLocalIntegral(int xbin, int Irange ){
    
    double localI = 0;
   
    // compute integral around the peak
    int startI;
    int stopI;
    //double localI=0;
    //Int_t Irange =10; 

    if ((xbin- Irange/2) > 0)
        startI = xbin-Irange/2;
    else
        startI = 0;

    //if ((xbin+Irange/2) < avgWaveform->size())
    if ((xbin+Irange/2) < avgWaveform.size())
        stopI = xbin+Irange/2;
    else
        stopI = avgWaveform.size();
        //stopI = avgWaveform->size();

    //if (verbose) std::cout  << " bin " << xbin << " startI " << startI << "  stopI " << stopI <<std::endl;
    for (int isample = startI; isample<stopI; isample++){
        //localI += avgWaveform->at(isample);
        localI += avgWaveform.at(isample);
    } 

    peak_integral.push_back(localI);
    return localI; 
}



void Event::ComputeTailIntegral(int preSamples)
{
    tailIntegral = 0.0;

    if (mainPeakIdx < 0 || mainPeakIdx >= (int)peakPositions.size())
        return;

    int peakSample = static_cast<int>(std::round(peakPositions[mainPeakIdx]));
    int startSample = (peakSample - preSamples >= 0) ? (peakSample - preSamples) : 0;

    for (size_t i = startSample; i < avgWaveform.size(); ++i) {
        tailIntegral += avgWaveform.at(i);
    }
}




void Event::FindMainPeak()
{
    /*
    if (peakAmplitudes.empty()) {
        mainPeakIdx = -1;
        return;
    }

    size_t imax = 0;
    maxAmp = peakAmplitudes[0];

    for (size_t i = 1; i < peakAmplitudes.size(); ++i) {
        if (peakAmplitudes[i] > maxAmp) {
            maxAmp = peakAmplitudes[i];
            imax = i;
        }
    }

    mainPeakIdx = static_cast<int>(imax);*/


    mainPeakIdx = -1;
    maxAmp = -1.;

    if (goodPeakAmplitudes.empty()) {
        std::cerr << "[FindMainPeak]:  Channel " << channelId << ":  No good peaks found, skipping.\n";
        return;
    }

    for (size_t i = 0; i < goodPeakAmplitudes.size(); ++i) {
        if (goodPeakAmplitudes[i] > maxAmp) {
            maxAmp = goodPeakAmplitudes[i];
            mainPeakIdx = static_cast<int>(i);
        }
    }

}


void Event::FindMaxAmp(){
    //maxAmp = *max_element(Waveform->begin(), Waveform->end());
    maxAmp = *max_element(Waveform.begin(), Waveform.end());

}


void Event::FindPeaks(int maxPeaks, float sigma, float threshold, bool verbose) {
    peakPositions.clear();
    peakAmplitudes.clear();

    // Sanity check
    if (avgWaveform.empty()) {
        std::cerr << "[FindPeaks] Empty waveform, skipping.\n";
        return;
    }

    int NSamples = avgWaveform.size();

    // Create histogram from waveform
    TH1F* h_AvgMeanwaveform = new TH1F("h_AvgMeanwaveform","", NSamples,0,1023);

    for(int isample=0; isample<NSamples; isample++){
        if (isample>10) {
            h_AvgMeanwaveform->SetBinContent(isample+1, avgWaveform.at(isample));
        } else {
            h_AvgMeanwaveform->SetBinContent(isample+1, 0);
        }
    }   
    

    // Run TSpectrum
    TSpectrum *spectrum = new TSpectrum(2 * maxPeaks);
    int nFound = spectrum->Search(h_AvgMeanwaveform, sigma, "goff", threshold);

    if (verbose)
        std::cout << "[FindPeaks] Found " << nFound << " candidate peaks\n";

    Double_t *xPeaks = spectrum->GetPositionX();

    for (int p = 0; p < nFound; p++) {
        Double_t xp = xPeaks[p];
        Int_t bin = h_AvgMeanwaveform->GetXaxis()->FindBin(xp);
        Double_t yp = h_AvgMeanwaveform->GetBinContent(bin);

        if (verbose)
            std::cout << "Peak " << p << ": x = " << xp << ", y = " << yp << std::endl;

        peakPositions.push_back(xp);
        peakAmplitudes.push_back(yp);
    }

    delete spectrum;
    delete h_AvgMeanwaveform;
}



void Event::AnalyzePeaks(double SPE_INTEGRAL, TH1* h_AsymLocalIntegral, TH1* h_AsympeakInt, bool verbose)
{
    int Nsamples = avgWaveform.size();

    //std::cout << "Found NPeaks: " << peakPositions.size() << std::endl;

    for (size_t ip = 0; ip < peakPositions.size(); ++ip) {
        double xp = peakPositions[ip];
        int bin = static_cast<int>(xp + 0.5);

        int range_low = 15;
        int range_up = 40;

        int lowedge_asymInt = std::max(0, bin - range_low);
        int upedge_asymInt = std::min(Nsamples - 1, bin + range_up);

        double asymlocalI = 0;
        for (int isample = lowedge_asymInt; isample < upedge_asymInt; ++isample) {
            double val = avgWaveform.at(isample);
            asymlocalI += val;

            if (h_AsymLocalIntegral)
                h_AsymLocalIntegral->SetBinContent(isample, val);
        }

        if (h_AsympeakInt)
            h_AsympeakInt->Fill(asymlocalI);

        // Save integral for this peak
        asympeak_integral.push_back(asymlocalI);

        if (asymlocalI > SPE_INTEGRAL) {
            ngoodpeaks_specut++;

            // Save only good peaks
            goodPeakPositions.push_back(xp);
            goodPeakAmplitudes.push_back(peakAmplitudes[ip]);
            goodPeakIntegrals.push_back(asymlocalI);
    
            if (verbose) {
                std::cout << "Peak at x=" << xp
                          << " passed SPE integral cut: " << asymlocalI << std::endl;
            }
        } else {
            if (verbose) {
                std::cout << "Peak at x=" << xp
                          << " rejected (integral too small): " << asymlocalI << std::endl;
            }
        }
    }
    //std::cout << "[Event::AnalyzePeaks] Found " 
    //          << goodPeakPositions.size() << " good peaks after SPE integral cut." << std::endl;

}


void Event::SeparateLeftRightPeaks(bool verbose)
{
    leftPeakPositions.clear();
    leftPeakAmplitudes.clear();
    rightPeakPositions.clear();
    rightPeakAmplitudes.clear();

    if (mainPeakIdx < 0 || mainPeakIdx >= static_cast<int>(peakPositions.size())) {
        std::cerr << "[Event::SeparateLeftRightPeaks] Invalid main peak index." << std::endl;
        return;
    }

    double mainPeakPos = goodPeakPositions[mainPeakIdx];

    for (size_t i = 0; i < goodPeakPositions.size(); ++i) {
        if (i == static_cast<size_t>(mainPeakIdx))
            continue;

        if (goodPeakPositions[i] < mainPeakPos) {
            leftPeakPositions.push_back(goodPeakPositions[i]);
            leftPeakAmplitudes.push_back(goodPeakAmplitudes[i]);
        } else {
            rightPeakPositions.push_back(goodPeakPositions[i]);
            rightPeakAmplitudes.push_back(goodPeakAmplitudes[i]);
            if (goodPeakPositions[i] < 0){
                std::cout << "main peak position: " << peakPositions.at(mainPeakIdx) << std::endl;
                std::cout << "right peak position: " << goodPeakPositions[i] << std::endl;
                
            }
            
        }
    }
    if (verbose) {
        
    std::cout << "[Event::SeparateLeftRightPeaks] Found "
              << leftPeakPositions.size() << " peaks left, "
              << rightPeakPositions.size() << " peaks right of main peak." << std::endl;
    }
}


void Event::ComputePeakDistances()
{
    distancesRight.clear();

    double mainPos = peakPositions[mainPeakIdx];

    for (const auto& pos : rightPeakPositions) {
        distancesRight.push_back(pos - mainPos);
        if ((pos - mainPos) <0){
            std::cout << "main peak position: " << peakPositions.at(mainPeakIdx) << std::endl;
            std::cout << "right peak position: " << pos << std::endl;
                
        }
    }
}



void Event::ComputeSecPeakInterdistance(int maxpeak_x, std::vector <double> peak_x, std::vector<TH1F*> &h_tmp)
{

    //peak_interdistance.clear();
    peak_interdistance.reserve(peak_x.size());
    for (int i=0; i<peak_interdistance.size(); i++){
        peak_interdistance.emplace_back(std::vector<double>());
    }
    
    if (peak_x.size()){
    //compute first the distance between the first secondary peak and the maxAmp peak
        double maxAmp_dist = peak_x.at(0) - maxpeak_x;
        h_tmp[0]->Fill(maxAmp_dist);
        peak_interdistance[0].emplace_back(maxAmp_dist);

        double tmp_dist = -999;
        int maxNpeaks = peak_x.size();

        if (maxNpeaks >10)
            maxNpeaks = 10; // limit the number of peaks to be considered to limit the size of the array of histogram.

        for (int ix=1; ix<maxNpeaks; ix++){
            tmp_dist = peak_x.at(ix) - peak_x.at(ix-1);
            h_tmp[ix]->Fill(tmp_dist);
            //std::cout << "before filling the vector , ix " << ix << std::endl;
            peak_interdistance[ix].emplace_back(tmp_dist);
            //std::cout << peak_interdistance[ix].size() << " elements in the vector " << ix << std::endl;
        }
    }

    return;
}