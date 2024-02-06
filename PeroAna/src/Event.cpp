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

const std::vector<double>* Event::GetWaveform() const {
     return Waveform;
}



void Event::ComputeBaseline(){

    double sum = 0;
    int start = 5;
    int end = start+100;
    for(int i=start;i<end;i++) 
        sum += rawWaveform->at(i);
    
    baseline = sum/(end-start);
    
}


void Event::SubtractBaseline(){
        Waveform = new std::vector<double>(*rawWaveform);

        for(int i=0;i<1024;i++) {
            (*Waveform)[i] -= baseline;
        }


}



void Event::ComputeIntegral(){
    integral = 0;
    for(int i=0;i<1024;i++) integral += Waveform->at(i)*sampling_rate/50.; //this is pC 
}

void Event::FindMaxAmp(){
    maxAmp = *max_element(Waveform->begin(), Waveform->end());

}