#ifndef EVENT_H
#define EVENT_H

#include <vector>


class Event {
public:
    Event();
    ~Event();

    void Init();
    // Setter methods
    void setEventId(int id);
    void setChannelId(int ch);

    //void setChannels(const std::vector<int>& channels);

    // Getter methods
    int getEventId() const;
    int getChannelId() const;
    //const std::vector<int>& getChannels() const;

private:
    int eventId;
    int channelId;
    //std::vector<double> waveform;
};

#endif // EVENT_H
