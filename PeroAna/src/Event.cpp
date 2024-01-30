#include "Event.h"
#include <iostream>


//Event::Event() : eventId(0), channelId(0) {}

Event::Event() {
eventId =0; 
channelId = 0; 
Init();

}


Event::~Event() {}

void Event::Init(){

    std::cout<< " Initialize event " << std::endl;
    //setEventId();
    //channelId();


}

void Event::setEventId(int id) {
    eventId = id;
}

void Event::setChannelId(int ch) {
    channelId = ch;
}

int Event::getEventId() const {
    return eventId;
}

int Event::getChannelId() const {
    return channelId;
}
