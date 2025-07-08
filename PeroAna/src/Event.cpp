#include "Event.h"
#include <iostream>
#include <algorithm>




Event::Event(int evt, int chan, const std::vector<double>& wave)
    : eventId(evt), channelId(chan), rawWaveform(new std::vector<double>(wave)) 
{
    Init();
    // Additional initialization if needed
}


Event::~Event() {
    delete rawWaveform; 
    delete Waveform; 
}

void Event::Init(){

    //std::cout<< " Initialize event " << std::endl;
    tailIntegral=0; 
    ngoodpeaks=0;
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


int Event::GetEventId() const {
    return eventId;
}

int Event::GetChannelId() const {
    return channelId;
}

const std::vector<double>* Event::GetRawWaveform() const {
    return rawWaveform;
}

const std::vector<double>* Event::GetAvgMeanWaveform() const {
    return avgWaveform;
}


const std::vector<double>* Event::GetWaveform() const {
     return Waveform;
}

void Event::ComputeMovingAverage(int step, bool debug){
   //if (debug) std::cout << "entering the Compute Moving Average" <<std::endl;

    avgWaveform = new std::vector<double>(*rawWaveform); // create a new WF to recorde the one after moving average

    double sum = 0;
    int start = step;
    int end = rawWaveform->size();

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
                sum += rawWaveform->at(isum);
                //if (debug) std::cout << "sum " << sum << std::endl;
            }
            avgWaveform->at(isample) = (float)sum/step;
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
            sum += rawWaveform->at(i);
        else 
            sum += avgWaveform->at(i);
    }
    baseline = sum/(end-start);
    
}


void Event::SubtractBaseline(bool _removenoise){
        if (! _removenoise)
            Waveform = new std::vector<double>(*rawWaveform);
        else 
            Waveform = new std::vector<double>(*avgWaveform);

        for(int i=0;i<Waveform->size();i++) {
            (*Waveform)[i] -= baseline;
        }


}



void Event::ComputeIntegral(){
    integral = 0;
    //for(int i=0;i<1024;i++) integral += Waveform->at(i)*sampling_rate/50.; //this is pC 
    for(int i=0;i<Waveform->size();i++) integral += Waveform->at(i);  // for the moment without conversion
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

    if ((xbin+Irange/2) < avgWaveform->size())
        stopI = xbin+Irange/2;
    else
        stopI = avgWaveform->size();

    //if (verbose) std::cout  << " bin " << xbin << " startI " << startI << "  stopI " << stopI <<std::endl;
    for (int isample = startI; isample<stopI; isample++){
        localI += avgWaveform->at(isample);
    } 

    peak_integral.push_back(localI);
    return localI; 
}


void Event::FindMaxAmp(){
    maxAmp = *max_element(Waveform->begin(), Waveform->end());

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