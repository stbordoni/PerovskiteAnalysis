//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed May 22 15:45:27 2024 by ROOT version 6.24/06
// from TTree eventTree/Event Information
// found on file: /Users/bordonis/ResearchActivities/PerovskiteAnalysis/PeroAna/rootinputfiles/Run_000_Data_2_21_2024_Ascii_output.root
//////////////////////////////////////////////////////////

#ifndef myAna_h
#define myAna_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

class myAna {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           event;
   Int_t           channelId;
   Int_t           fcr;
   Double_t        baseline;
   Double_t        amplitude;
   Double_t        charge;
   Double_t        leadingEdgeTime;
   Double_t        trailingEdgeTime;
   Double_t        rateCounter;
   vector<double>  *dataSamples;

   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_channelId;   //!
   TBranch        *b_fcr;   //!
   TBranch        *b_baseline;   //!
   TBranch        *b_amplitude;   //!
   TBranch        *b_charge;   //!
   TBranch        *b_leadingEdgeTime;   //!
   TBranch        *b_trailingEdgeTime;   //!
   TBranch        *b_rateCounter;   //!
   TBranch        *b_dataSamples;   //!

   myAna(TTree *tree=0);
   myAna(string filename);
   virtual ~myAna();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef myAna_cxx
myAna::myAna(string filename) : fChain(0) 
{
TTree *tree; 
TFile *f = new TFile(filename.c_str());
  
  f->GetObject("eventTree",tree);
  
  Init(tree);
}


myAna::myAna(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/Users/bordonis/ResearchActivities/PerovskiteAnalysis/PeroAna/rootinputfiles/Run_000_Data_2_21_2024_Ascii_output.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/Users/bordonis/ResearchActivities/PerovskiteAnalysis/PeroAna/rootinputfiles/Run_000_Data_2_21_2024_Ascii_output.root");
      }
      f->GetObject("eventTree",tree);

   }
   Init(tree);
}

myAna::~myAna()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t myAna::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t myAna::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void myAna::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   dataSamples = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("channelId", &channelId, &b_channelId);
   fChain->SetBranchAddress("fcr", &fcr, &b_fcr);
   fChain->SetBranchAddress("baseline", &baseline, &b_baseline);
   fChain->SetBranchAddress("amplitude", &amplitude, &b_amplitude);
   fChain->SetBranchAddress("charge", &charge, &b_charge);
   fChain->SetBranchAddress("leadingEdgeTime", &leadingEdgeTime, &b_leadingEdgeTime);
   fChain->SetBranchAddress("trailingEdgeTime", &trailingEdgeTime, &b_trailingEdgeTime);
   fChain->SetBranchAddress("rateCounter", &rateCounter, &b_rateCounter);
   fChain->SetBranchAddress("dataSamples", &dataSamples, &b_dataSamples);
   Notify();
}

Bool_t myAna::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void myAna::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t myAna::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef myAna_cxx
