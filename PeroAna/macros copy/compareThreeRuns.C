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

void compareThreeRuns(TString f_1 ="run1",
                     TString f_2 ="run2",
                     TString f_3 = "run3 ", 
                     TString label1 = "label1 ",
                     TString label2 = "label2 ",
                     TString label3 = "label3 ")
{
    
    TTree *tree1;
    TTree *tree2;
    TTree *tree3;

    TLegend *leg = new TLegend();

    Bool_t tree1_ok = false;
    Bool_t tree2_ok = false;
    Bool_t tree3_ok = false;

    if (! f_1.IsNull()) {
        //TTree *tree1 = GetTree(f_1);
        tree1 = GetTree(f_1);
        tree1->Print();
        tree1->SetLineColor(1);
        tree1->SetLineWidth(3);
        tree1_ok = true;
        //tree1->SetLineStyle(3);
    }

    if (! f_2.IsNull()) {
        tree2 = GetTree(f_2);
        tree2->Print();
        tree2->SetLineColor(4);
        tree2->SetLineWidth(3);
        tree2_ok = true;
    }

    if (! f_3.IsNull()) {     
        tree3 = GetTree(f_3);
        tree3->Print();
        tree3->SetLineColor(kGray);
        tree3->SetLineWidth(3);
        tree3_ok = true;
    }
   

    TCanvas *c0 = new TCanvas("c0", "peak Amplitudes", 800, 600);
    c0->cd();
    if (tree1_ok) tree1->Draw("peakAmp>>htemp1(100, 0, 0.05)");     
    //htemp1->GetXaxis()->SetTitle("[V]");
    //htemp1->GetYaxis()->SetTitle("counts");
     if (tree2_ok) tree2->Draw("peakAmp>>htemp2(100, 0, 0.05)", "", "sames");
     if (tree3_ok) tree3->Draw("peakAmp>>htemp3(100, 0, 0.05)", "", "sames");

    
    if (tree1_ok) leg->AddEntry("htemp1", label1);
    if (tree2_ok) leg->AddEntry("htemp2", label2);
    if (tree3_ok) leg->AddEntry("htemp3", label3);
    leg->Draw("same");

    

    TCanvas *c1 = new TCanvas("c1", "Local integral", 800, 600);
    c1->cd();
    if (tree1_ok) tree1->Draw("peakInt>>htemp1int(100, 0, 0.5)");     
    //htemp1->GetXaxis()->SetTitle("[V]");
    //htemp1->GetYaxis()->SetTitle("counts");
    if (tree2_ok) tree2->Draw("peakInt>>htemp2int(100, 0, 0.5)", "", "sames");
    if (tree3_ok) tree3->Draw("peakInt>>htemp3int(100, 0, 0.5)", "", "sames");
    leg->Draw("same");

    TCanvas *c2 = new TCanvas("c2", " ", 1600, 600);
    c2->Divide(2,1);
    c2->cd(1);
    if (tree1_ok) tree1->Draw("tailInt>>htemp1tail(100, 0, 10)");     
    //htemp1->GetXaxis()->SetTitle("[V]");
    //htemp1->GetYaxis()->SetTitle("counts");
    if (tree2_ok) tree2->Draw("tailInt>>htemp2tail(100, 0, 10)", "", "sames");
    if (tree3_ok) tree3->Draw("tailInt>>htemp3tail(100, 0, 10)", "", "sames");
    leg->Draw("same");
    
    c2->cd(2);
    if (tree1_ok) tree1->Draw("npeaks>>htemp1npeaks(12, -0.5, 11.5)");     
    //htemp1->GetXaxis()->SetTitle("[V]");
    //htemp1->GetYaxis()->SetTitle("counts");
    if (tree2_ok) tree2->Draw("npeaks>>htemp2npeaks(12, -0.5, 11.5)", "", "sames");
    if (tree3_ok) tree3->Draw("npeaks>>htemp3npeaks(12, -0.5, 11.5)", "", "sames");

}