#include "TFile.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TF1.h"
#include "TTree.h"
#include "TPave.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TLatex.h"


#include "TMath.h"
#include <numeric>


#include "TF1.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"


void compareRuns(){
    
    //gStyle->SetLegendFont(54);


    TString PATH = "/Users/bordonis/ResearchActivities/PerovskiteAnalysis/PeroAna/rootinputfiles/";

    
    //TString f_Cs ="Run_012_Cs104-A_Data_9_24_2024_Ascii_ana.root";   //10k
    //TString f_Cs_Sr90 ="Run_013_Cs104-A_Data_9_24_2024_Ascii_ana.root";  //10k

    TString f_Bkg ="Run_21_Data_9_2_2024_Ascii_ana.root";   //10k
    TFile *file_Bkg= new TFile(PATH+f_Bkg, "OPEN");
    TTree *t_Bkg = (TTree*)file_Bkg->Get("output_tree");

    TString f_Cs ="Run_018_Cs104-A_Data_9_27_2024_Ascii_ana.root";   //50k
    TString f_Cs_Sr90 ="Run_019_Cs104-A_Data_9_27_2024_Ascii_ana.root";  //50k

    TFile *file_Cs= new TFile(PATH+f_Cs, "OPEN");
    TTree *t_Cs = (TTree*)file_Cs->Get("output_tree");
    t_Cs->SetLineColor(2);
    t_Cs->SetLineWidth(3);
    //t_Cs->Print();
    

    TFile *file_Cs_Sr90= new TFile(PATH+f_Cs_Sr90, "OPEN");
    TTree *t_Cs_Sr90 = (TTree*)file_Cs_Sr90->Get("output_tree");
    t_Cs_Sr90->SetLineColor(2);
    t_Cs_Sr90->SetLineWidth(3);
    t_Cs_Sr90->SetLineStyle(3);
    //t_Cs_Sr90->Print();

    
    //TString f_MA ="Run_001_SC67-B_Data_9_27_2024_Ascii_ana.root";
    //TString f_MA_Sr90 ="Run_002_SC67-B_Data_9_27_2024_Ascii_ana.root";

    //TString f_MA ="Run_004_SC67-B_Data_9_27_2024_Ascii_ana.root";  //50k no source
    //TString f_MA_Sr90 ="Run_003_SC67-B_Data_9_27_2024_Ascii_ana.root"; //50k w/ source
    
    TString f_MA ="Run_007_SC67-B_Data_10_3_2024_Ascii_ana.root";  //10k no source thsdl 10mV
    TString f_MA_Sr90 ="Run_006_SC67-B_Data_10_3_2024_Ascii_ana.root"; //10k w/ source thsdl 10mV
    

    TFile *file_MA = new TFile(PATH+f_MA, "OPEN");
    TTree *t_MA = (TTree*)file_MA->Get("output_tree");
    t_MA->SetLineColor(4);
    t_MA->SetLineWidth(3);
    
    TFile *file_MA_Sr90 = new TFile(PATH+f_MA_Sr90, "OPEN");
    TTree *t_MA_Sr90 = (TTree*)file_MA_Sr90->Get("output_tree");
    t_MA_Sr90->SetLineColor(4);
    t_MA_Sr90->SetLineWidth(3);
    t_MA_Sr90->SetLineStyle(3);
    


    TCanvas *c0 = new TCanvas("c0", "peak Amplitudes", 800, 600);
    c0->cd();

    TH1F *h_peakAmp_Bkg = new TH1F ("h_peakAmp_Bkg", "peak Amp",  250, -0.001, 0.06); 
    //h_peakAmp_Bkg->SetTitle("CsPbBr3");
    h_peakAmp_Bkg->SetLineColor(1);
    h_peakAmp_Bkg->SetLineWidth(3);


    TH1F *h_peakAmp_Cs = new TH1F ("h_peakAmp_Cs", "peak Amp",  250, -0.001, 0.06); 
    h_peakAmp_Cs->SetTitle("CsPbBr3");
    h_peakAmp_Cs->SetLineColor(4);
    h_peakAmp_Cs->SetLineWidth(3);

    TH1F *h_peakAmp_Cs_Sr90 = new TH1F ("h_peakAmp_Cs_Sr90", " peak Amp w/ Sr90",  250, -0.001, 0.06); 
    h_peakAmp_Cs_Sr90->SetTitle("CsPbBr3");
    h_peakAmp_Cs_Sr90->SetLineColor(4);
    h_peakAmp_Cs_Sr90->SetLineWidth(3);
    h_peakAmp_Cs_Sr90->SetLineStyle(3);


    TLegend *leg_Cs = new TLegend();
    
    leg_Cs->AddEntry(h_peakAmp_Cs, "Cs104-A ", "l");
    leg_Cs->AddEntry(h_peakAmp_Cs_Sr90, "Cs104-A w/ Sr90 ", "l");
    
    
    t_Cs->Draw("peakAmp>>h_peakAmp_Cs", "", "");
    t_Cs_Sr90->Draw("peakAmp>>h_peakAmp_Cs_Sr90", "", "sames");
    leg_Cs->Draw("same");





    

    TH1F *h_peakAmp_MA = new TH1F ("h_peakAmp_MA", "peak Amp",  250, -0.001, 0.06); 
    h_peakAmp_MA->SetTitle("MAPbBr3");
    h_peakAmp_MA->SetLineColor(2);
    h_peakAmp_MA->SetLineWidth(3);

    TH1F *h_peakAmp_MA_Sr90 = new TH1F ("h_peakAmp_MA_Sr90", " peak Amp w/ Sr90",  250, -0.001, 0.06); 
    h_peakAmp_MA_Sr90->SetTitle("MAPbBr3");
    h_peakAmp_MA_Sr90->SetLineColor(2);
    h_peakAmp_MA_Sr90->SetLineWidth(3);
    h_peakAmp_MA_Sr90->SetLineStyle(3);

    TLegend *leg_MA = new TLegend();
    leg_MA->AddEntry(h_peakAmp_MA, "SR67-B ", "l");
    leg_MA->AddEntry(h_peakAmp_MA_Sr90, "SR67-B w/ Sr90 ", "l");
    
    TCanvas *c1 = new TCanvas("c1", "peak Amplitudes", 800, 600);
    c1->cd();
    t_MA->Draw("peakAmp>>h_peakAmp_MA", "", "");
    t_MA_Sr90->Draw("peakAmp>>h_peakAmp_MA_Sr90", "", "sames");
    leg_MA->Draw("same");




    TCanvas *c2 = new TCanvas("c2", "peak Amplitudes no sources", 800, 600);
    c2->cd();
    h_peakAmp_MA->Draw("");
    h_peakAmp_Cs->Draw("sames");
    t_Bkg->Draw("peakAmp>>h_peakAmp_Bkg", "", "sames");

    TLegend *leg = new TLegend();
    leg->AddEntry(h_peakAmp_MA, "SR67-B ", "l");
    leg->AddEntry(h_peakAmp_Cs, "Cs104-A ", "l");
    leg->AddEntry(h_peakAmp_Bkg, "Bkg ", "l");
    leg->Draw("same");


    TCanvas *c3 = new TCanvas("c3", "peak Amplitudes W/ sources", 800, 600);
    c3->cd();
    h_peakAmp_MA_Sr90->Draw("");
    h_peakAmp_Cs_Sr90->Draw("sames");

    TLegend *leg_Sr90 = new TLegend();
    leg_Sr90->AddEntry(h_peakAmp_MA_Sr90, "SR67-B w/ Sr90 ", "l");
    leg_Sr90->AddEntry(h_peakAmp_Cs_Sr90, "Cs104-A w/Sr90 ", "l");
    leg_Sr90->Draw("same");

}