#ifndef EVENT_H
#define EVENT_H

#include <vector>
#include "TH1F.h"
#include "Config.h"


class Event {
public:
    //Event();
   
    //Event(int evt, int chan, const std::vector<double>& wave);
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
    double GetBaseline() const {return baseline;}
    double GetPulseIntegral() const {return integral;}
    double GetTailIntegral() const {return tailIntegral;}
    double GetMaxAmp() const {return maxAmp;}
    int GetMainPeakIdx() const { return mainPeakIdx; }
    int GetNumberOfPeaks() const {return GetPeakPositions().size();}
    int GetNumberOfPeaksWithSpecut() const {return ngoodpeaks_specut;}
    double GetPhotonCountInTail() const {return tailIntegral / SPE_INTEGRAL;}
    int GetNRightPeaks() const { return rightPeakPositions.size(); }
    int GetNLeftPeaks() const { return leftPeakPositions.size(); }
    //double GetLeftDistanceFromMax() const {return distmaxAmppeaks_left;}
    
    
    //int GetRightDistanceFromMax() const {return distmaxAmppeaks_right
    //const std::vector<double>* GetWaveform() const; 
    //const std::vector<double>* GetRawWaveform() const; 
    //const std::vector<double>* GetAvgMeanWaveform() const; 
    const std::vector<double> GetWaveform() const; 
    const std::vector<double> GetRawWaveform() const; 
    const std::vector<double> GetAvgMeanWaveform() const; 
    const std::vector<double>& GetRightDistanceFromMax() const { return distancesRight; }
    

    void ComputeMovingAverage(int step=5, bool debug=false);
    void ComputeBaseline(bool _removenoise);
    void SubtractBaseline(bool _removenoise);
    void ComputeIntegral();
    double ComputeLocalIntegral(int xbin, int Irange=10 );
    void FindMainPeak();
    
    //double ComputeTailIntegral(int mainPeakIdx, int preSamples) const
    void ComputeTailIntegral(int preSamples = 25);
    void SeparateLeftRightPeaks(bool verbose = false);
    void ComputePeakDistances();

    void FindMaxAmp();
    //void ComputeSecPeakInterdistance(int maxpeak_x, std::vector <double> peak_x, TH1F &h_tmp);
    void ComputeSecPeakInterdistance(int maxpeak_x, std::vector <double> peak_x, std::vector<TH1F*> &h_tmp);

    // Finds peaks and stores their positions (x) and amplitudes (y)
    void FindPeaks(int maxPeaks, float sigma, float threshold, bool verbose = false);
    void AnalyzePeaks(double SPE_INTEGRAL, TH1* h_AsymLocalIntegral = nullptr, TH1* h_AsympeakInt = nullptr, bool verbose = false);

    // Accessors to get peak results
    const std::vector<double>& GetPeakPositions() const { return peakPositions; }
    const std::vector<double>& GetPeakAmplitudes() const { return peakAmplitudes; }

    //const std::vector<int>& getChannels() const;


    double baseline;
    double integral;
    double tailIntegral;
    double maxAmp;
    int mainPeakIdx; // Index of the main peak in the waveform
    double ToT; 
    
    int ngoodpeaks;
    int ngoodpeaks_specut;
    //int distmaxAmppeaks_right;
    //int distmaxAmppeaks_left;
    std::vector<double> x_peak;
    std::vector<double> y_peak;
    std::vector<double> peak_integral;
    std::vector<double> asympeak_integral;

    std::vector<double> peakPositions;
    std::vector<double> peakAmplitudes;

    std::vector<double> leftPeakPositions;
    std::vector<double> leftPeakAmplitudes;
    std::vector<double> rightPeakPositions;
    std::vector<double> rightPeakAmplitudes;
    std::vector<std::vector<double>> peak_interdistance;

    

    //std::vector<double>* Waveform;
    //std::vector<double>* avgWaveform;

    std::vector<double> Waveform;
    std::vector<double> avgWaveform;

    std::vector<double> distancesRight;


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
    std::vector<double> rawWaveform;

    

    

    
};

#endif // EVENT_H
