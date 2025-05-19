#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TSpectrum.h>
#include <TPolyMarker.h>

TTree* GetTree(TString filename){   
    TString PATH = "/Users/bordonis/ResearchActivities/PerovskiteAnalysis/PeroAna/rootinputfiles/";
    TTree *tree_0;

    if (filename){
        TFile *file_0= new TFile(PATH+filename, "OPEN");
        tree_0 = (TTree*)file_0->Get("output_tree");
    }
    else
        tree_0;

    return tree_0;
}

Int_t npeaks = 20;
Double_t fpeaks(Double_t *x, Double_t *par) {
   Double_t result = par[0] + par[1]*x[0];
   for (Int_t p=0;p<npeaks;p++) {
      Double_t norm  = par[3*p+2]; // "height" or "area"
      Double_t mean  = par[3*p+3];
      Double_t sigma = par[3*p+4];
#if defined(__PEAKS_C_FIT_AREAS__)
      norm /= sigma * (TMath::Sqrt(TMath::TwoPi())); // "area"
#endif /* defined(__PEAKS_C_FIT_AREAS__) */
      result += norm*TMath::Gaus(x[0],mean,sigma);
   }
   return result;
}


void fingerplot(TString f_1 ="file", Bool_t verbose=true){

    //get the tree 
    TTree *tree1;

    if (f_1.IsNull()) return;
       
    tree1 = GetTree(f_1);
    tree1->SetLineColor(1);
    tree1->SetLineWidth(3);


    // get the histogram with the maximum amplitudes and draw it
    TH1F* h_peakAmp   = new TH1F("h_peakAmp","", 120, 10, 10);

    TCanvas *c0 = new TCanvas("c0", "peak Amplitudes", 1600, 600);
    c0->Divide(2,1);
    c0->cd(1);
    gPad->SetLogy();
    tree1->Draw("peakAmp>> h_peakAmp");

    // find the peaks of the spectrum  
    Int_t np=20;
    Int_t npeaks = TMath::Abs(np);
    Double_t par[100];

    // Use TSpectrum to find the peak candidates
    TSpectrum *s = new TSpectrum(2*npeaks);
    
    Int_t nfound = s->Search(h_peakAmp, 3, "", 0.0001);
    if (verbose) printf("Found %d candidate peaks to fit\n",nfound);

    TF1 *fline = new TF1("fline","expo",0.001,0.02);
    h_peakAmp->Fit("fline","");
    par[0] = fline->GetParameter(0);
    par[1] = fline->GetParameter(1);

    //Estimate background using TSpectrum::Background
    TH1 *hb = s->Background(h_peakAmp,20,"same");
    hb->SetLineStyle(5);
    if (hb) c0->Update();

    c0->Update();

    npeaks = 0;
    Double_t *xpeaks;
    xpeaks = s->GetPositionX();

    for (Int_t p=0;p<nfound;p++) {
        Double_t xp = xpeaks[p];
        if (verbose) std::cout << "p  " << p << " : xp " << xp << std::endl;     

            
        Int_t bin = h_peakAmp->GetXaxis()->FindBin(xp);
        Double_t yp = h_peakAmp->GetBinContent(bin);
        //if (yp-TMath::Sqrt(yp) < fline->Eval(xp)) continue;

        par[3*npeaks+2] = yp; // "height"
        par[3*npeaks+3] = xp; // "mean"
        par[3*npeaks+4] = 3.5e-4; // "sigma"

        //if (verbose){
        //    std::cout << "par["<<3*npeaks+2<<"] = " << par[3*npeaks+2] << std::endl;
        //    std::cout << "par["<<3*npeaks+3<<"] = " << par[3*npeaks+3] << std::endl;
        //    std::cout << "par["<<3*npeaks+4<<"] = " << par[3*npeaks+4] << std::endl;
        //}

        #if defined(__PEAKS_C_FIT_AREAS__)
            par[3*npeaks+2] *= par[3*npeaks+4] * (TMath::Sqrt(TMath::TwoPi())); // "area"
        #endif /* defined(__PEAKS_C_FIT_AREAS__) */
            npeaks++;
    }
        
    printf("Found %d useful peaks to fit\n",npeaks);
    printf("Now fitting: Be patient\n");
    c0->cd(2);
    TH1F *h2 = (TH1F*)h_peakAmp->Clone("h2");
    TF1 *fit = new TF1("fit",fpeaks,0,0.02,2+3*npeaks);
    fit->SetLineColor(4);

    TVirtualFitter::Fitter(h2,10+3*npeaks);
    fit->SetParameters(par);
    //fit->SetNpx(1000);
    //gPad->SetLogy();             
    h2->Fit("fit");
    //h2->Draw();


     
        
    

    Double_t FittedPeaks_x[10];

    for (int ipeak =0; ipeak<npeaks; ipeak++){
        std::cout << "peak " << ipeak << ": " << fit->GetParameter(3*ipeak+3) << std::endl;
        FittedPeaks_x[ipeak] = fit->GetParameter(3*ipeak+3);

        std::cout << "distance to the previous peaks "<< std::endl;

        if (ipeak!=0){
            std::cout << "Distance peaks " << ipeak << " -  " << ipeak-1 << " : " << FittedPeaks_x[ipeak] - FittedPeaks_x[ipeak-1] << std::endl;
        }
    }

    
    c0->Update();      
}