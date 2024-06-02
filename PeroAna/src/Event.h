#ifndef EVENT_H
#define EVENT_H

#include <vector>



class Event {
public:
    //Event();
   
    Event(int evt, int chan, const std::vector<double>& wave);
    ~Event();

    
    void Init();
    // Setter methods
    void SetEventId(int id);
    void SetChannelId(int ch);
    

    //void setChannels(const std::vector<int>& channels);

    // Getter methods
    int GetEventId() const;
    int GetChannelId() const;
    const std::vector<double>* GetWaveform() const; 
    const std::vector<double>* GetRawWaveform() const; 
    const std::vector<double>* GetAvgMeanWaveform() const; 
    void ComputeMovingAverage(int step=5);
    void ComputeBaseline();
    void SubtractBaseline();
    void ComputeIntegral();
    void FindMaxAmp();
    //const std::vector<int>& getChannels() const;


    double baseline;
    double integral;
    double maxAmp;
    double ToT; 
    std::vector<double>* Waveform;
    std::vector<double>* avgWaveform;


private:
    
    int eventId;
    int channelId;
    float sampling_rate; 
   
    //int fcr;
    //double baseline;
    //double amplitude;
    //double charge;
    //double leadingEdgeTime;
    //double trailingEdgeTime;
    //double rateCounter;
    std::vector<double>* rawWaveform =0;
    
};

#endif // EVENT_H
