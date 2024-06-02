#include "Event.h"
#include <iostream>




Event::Event(int evt, int chan, const std::vector<double>& wave)
    : eventId(evt), channelId(chan), rawWaveform(new std::vector<double>(wave)) 
{
    // Additional initialization if needed
}


Event::~Event() {
    delete rawWaveform; 
    delete Waveform; 
}

void Event::Init(){

    std::cout<< " Initialize event " << std::endl;
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

void Event::ComputeMovingAverage(int step){
    std::cout << "entering the Compute Moving Average" <<std::endl;

    avgWaveform = new std::vector<double>(*rawWaveform); // create a new WF to recorde the one after moving average

    double sum = 0;
    int start = step;
    int end = rawWaveform->size();

    for (int isample=0; isample < end; isample++){
    
        std::cout << "--- isample " << isample << std::endl;
        if ( ((isample - step/2) <0) || ( (isample + step/2) >end )) {
            continue;
        }
        else{
            std::cout << "start:  isample-step/2 : " <<  isample-step/2 <<  " stop : isample+step/2 : " << isample+step/2 << std::endl;
            sum = 0;
            
            for(int isum =isample-step/2; isum<isample+step/2; isum++){ 
                std::cout << "isum " << isum << std::endl;
                sum += rawWaveform->at(isum);
                std::cout << "sum " << sum << std::endl;
            }
            avgWaveform->at(isample) = (float)sum/step;
            std::cout << " isample  " << isample << " rawWaveform->at(isample)  " << rawWaveform->at(isample)  <<  ";  avgWaveform  " << avgWaveform->at(isample)<<std::endl;
        
        }

    }
   

}


void Event::ComputeBaseline(){

    double sum = 0;
    int start = 5;
    int end = start+50;
    for(int i=start;i<end;i++) 
        sum += rawWaveform->at(i);
    
    baseline = sum/(end-start);
    
}


void Event::SubtractBaseline(){
        Waveform = new std::vector<double>(*rawWaveform);

        for(int i=0;i<Waveform->size();i++) {
            (*Waveform)[i] -= baseline;
        }


}



void Event::ComputeIntegral(){
    integral = 0;
    //for(int i=0;i<1024;i++) integral += Waveform->at(i)*sampling_rate/50.; //this is pC 
    for(int i=0;i<Waveform->size();i++) integral += Waveform->at(i);  // for the moment without conversion
}

void Event::FindMaxAmp(){
    maxAmp = *max_element(Waveform->begin(), Waveform->end());

}