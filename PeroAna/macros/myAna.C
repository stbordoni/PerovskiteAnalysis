#define myAna_cxx
#include "myAna.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>



bool displaywave = true;

void myAna::Loop()
{
//   In a ROOT session, you can do:
//      root> .L myAna.C
//      root> myAna t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;


      // printout for the first events to see if the rootfile is read correctly and check the data consistency with the raw file
      if (jentry < 100){
      std::cout <<" entry "  << jentry  <<  " ;  EVENT " << event << std::endl;
     

      // to printout the waveform
        for (int isample =0; isample<dataSamples->size(); isample++ )
            std::cout << dataSamples->at(isample) << " " ;
        std::cout << "\n " << std::endl;

      


      if (displaywave){
         std::cout << "--- CH " << channelId << " ---" << std::endl;
         TH1F* h_waveform = new TH1F("h_waveform","", 1024,0,1023);
         for(int isample=0; isample<1024; isample++){
            if (isample>10) h_waveform->SetBinContent(isample,dataSamples->at(isample));
            else      
            h_waveform->SetBinContent(isample,0); // prevent first bins with strange values to be shown.
            
         }

         TCanvas* c = new TCanvas("c", "c", 600, 500);
         c->cd();
         h_waveform->Draw("HIST");
         //h_waveform->GetXaxis()->SetRangeUser(0,200);
         c->Update();
         c->WaitPrimitive();

         delete h_waveform; delete c;

      }   
   }


   }
}
