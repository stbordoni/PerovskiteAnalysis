#ifndef CONFIG_H
#define CONFIG_H

// Single photoelectron reference integral value
constexpr double SPE_INTEGRAL = 0.09636; // in V (estimated from run6 SiPM only )

// Number of samples in the waveform
constexpr int nSamples = 1024; // Number of samples in the waveform

//NChannels = 2; // Number of channels in the event
constexpr int NChannels = 2; // Number of channels in the event

//NmaxRightPeaks = 5; // Maximum number of right peaks to consider
constexpr int NmaxRightPeaks = 5; // Maximum number of right peaks to consider
#endif
