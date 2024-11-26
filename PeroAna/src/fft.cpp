
#include "TApplication.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TMath.h"
#include "TVirtualFFT.h"
#include "TRandom3.h"

#include <iostream>
using std::cout;
using std::endl;

const Double_t pi = TMath::Pi();
const Double_t twoPi = TMath::TwoPi();

int main(int argc, char **argv)
{
	TApplication theApp("App", &argc, argv);
	
	TCanvas *c = new TCanvas("c", "FFTW Practice (forward transform)", 800, 800);
	c->Connect("Closed()", "TApplication", &theApp, "Terminate()");
	
	TPad* c1 = new TPad("c1", "f(x)", 0.005, 0.755, 0.745, 0.995); c1->Draw();
	TPad* c2 = new TPad("c2", "f(w)", 0.005, 0.505, 0.745, 0.745); c2->Draw();
	TPad* c3 = new TPad("c3", "f(x)", 0.005, 0.255, 0.745, 0.495); c3->Draw();
	TPad* c4 = new TPad("c4", "f(x)", 0.005, 0.005, 0.745, 0.245); c4->Draw();
	
	TPad* c5 = new TPad("c5", "f(x)", 0.755, 0.755, 0.995, 0.995); c5->Draw();
	TPad* c6 = new TPad("c6", "f(w)", 0.755, 0.505, 0.995, 0.745); c6->Draw();
	TPad* c7 = new TPad("c7", "f(x)", 0.755, 0.255, 0.995, 0.495); c7->Draw();
	TPad* c8 = new TPad("c8", "f(x)", 0.755, 0.005, 0.995, 0.245); c8->Draw();


	
	Int_t nSamples = 16384;
	

	/**************************************************
	 
	 make a times series (i.e., a digitized waveform)
	 
	 **************************************************/
	
	// nSamples samples of f(x) = a sawtooth wave with period 2pi, defined from 0--20pi.
	/*
	TH1D h1dInit("h1dInit", "Sawtooth Wave", nSamples, 0, 10.*twoPi);
	for (Int_t bin = 1; bin <= nSamples; bin++){
		Double_t x = h1dInit.GetBinLowEdge(bin);
		Double_t value = x/twoPi - TMath::Floor(x/twoPi) - 0.5;
		h1dInit.SetBinContent(bin , value);
	}
	*/
	
	// white noise w/ rms = RMS
	const Double_t RMS = 1.;
	TRandom3 rand(0);
	TH1D h1dInit("h1dInit", "White Noise", nSamples, 0, 10.*twoPi);
	for (Int_t bin = 1; bin <= nSamples; bin++){
		h1dInit.SetBinContent(bin, rand.Gaus(0, RMS));
	}
	
	c1->cd();
	h1dInit.Draw();

	TH1D h1dInitProj("h1dInitProj", "", 100, -5, 5);
	for (Int_t i = 1; i <= h1dInit.GetNbinsX(); i++)
		h1dInitProj.Fill(h1dInit.GetBinContent(i));
	c5->cd();
	h1dInitProj.Draw();
	
	// initiate the FFT
	TVirtualFFT::SetTransform(0);
	
	
	/**************************************************
	 
	 discrete fourier transform to frequency space
	 
	 **************************************************/
	
	// get the magnitude (i.e., power) of the transform of f(x)
	TH1* h1dMagRaw = NULL;
	h1dMagRaw = h1dInit.FFT(h1dMagRaw, "MAG"); // this has units of 1/f_max
	TH1D h1dMag("h1dMag", "Magnitude (i.e., Power Spectrum)", h1dMagRaw->GetNbinsX(), 0, h1dMagRaw->GetXaxis()->GetXmax()/h1dInit.GetXaxis()->GetXmax()); 
	// rescale axis to get real units
	for (Int_t bin = 1; bin <= nSamples; bin++){
		h1dMag.SetBinContent(bin, h1dMagRaw->GetBinContent(bin));
	}
	
	
	c2->cd();
	h1dMag.Draw();
	
	TH1D h1dMagProj("h1dMagProj", "", 50, 0, 500);
	for (Int_t i = 1; i <= h1dMag.GetNbinsX(); i++)
		h1dMagProj.Fill(h1dMag.GetBinContent(i));
	c6->cd();
	h1dMagProj.Draw();
	
	
	// get the phase of the transform of f(x)
	TH1* h1dPhaseRaw = NULL;
	h1dPhaseRaw = h1dInit.FFT(h1dPhaseRaw, "PH"); // this has units of 1/f_max
	TH1D h1dPhase("h1dPhase", "Phase", h1dPhaseRaw->GetNbinsX(), 0, h1dPhaseRaw->GetXaxis()->GetXmax()/h1dInit.GetXaxis()->GetXmax()); 
	// rescale axis to get real units
	for (Int_t bin = 1; bin <= nSamples; bin++){
		h1dPhase.SetBinContent(bin, h1dPhaseRaw->GetBinContent(bin));
	}
	
	c3->cd();
	h1dPhase.Draw();
	
	TH1D h1dPhaseProj("h1dPhaseProj", "", 100, -pi, pi);
	for (Int_t i = 1; i <= h1dPhase.GetNbinsX(); i++)
		h1dPhaseProj.Fill(h1dPhase.GetBinContent(i));
	c7->cd();
	h1dPhaseProj.Draw("e");
	
	// get full re + i*im contents of the transform
	TVirtualFFT* fft = TVirtualFFT::GetCurrentTransform();
	Double_t re[nSamples];
	Double_t im[nSamples];
	fft->GetPointsComplex(re, im);
	
	
	/**************************************************
	 
	 now transform it back to check your head (i.e., 
	 see if it matches the original).
	 
	 **************************************************/
	
	//inverse transform
	TVirtualFFT* fft_inv = TVirtualFFT::FFT(1, &nSamples, "C2R M K");
	fft_inv->SetPointsComplex(re, im);
	fft_inv->Transform();
	TH1* h1dReconRaw = 0;
	h1dReconRaw = TH1::TransformHisto(fft_inv, h1dReconRaw, "Re");
	TH1D h1dRecon("h1dRecon", "Reconstructed Signal", h1dReconRaw->GetNbinsX(), 0, h1dInit.GetXaxis()->GetXmax()); 
	// rescale axis to get real units
	for (Int_t bin = 1; bin <= nSamples; bin++){
		h1dRecon.SetBinContent(bin, h1dReconRaw->GetBinContent(bin));
	}
	h1dRecon.Scale(1./((Double_t) nSamples));
	
	c4->cd();
	h1dRecon.Draw();
	
	TH1D h1dReconProj("h1dReconProj", "", 100, -5, 5);
	for (Int_t i = 1; i <= h1dRecon.GetNbinsX(); i++)
		h1dReconProj.Fill(h1dRecon.GetBinContent(i));
	c8->cd();
	h1dReconProj.Draw();
	

	
	c->Update();
	theApp.Run();
	
	delete c;
	delete c1;
	delete c2;
	delete c3;
	delete c4;
	
	return 0;
	
}