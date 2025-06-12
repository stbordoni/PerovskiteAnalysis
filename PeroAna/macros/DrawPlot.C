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

void DrawPlot(TString f_1 ="run1",
            TString f_2 ="run2", 
            TString varname = "peakInt", 
            TString label1 = "label1 ",
            TString label2 = "label2 ")
{

    TTree *tree1;
    TTree *tree2;

    TLegend *leg = new TLegend(0.57,0.73,0.87, 0.88);

    tree1 = GetTree(f_1);
    tree1->Print();
    tree1->SetFillColor(1);
    tree1->SetLineWidth(3);


    tree2 = GetTree(f_2);
    tree2->Print();
    //tree2->SetFillColor(4);
    //tree2->SetLineWidth(3);

    TCanvas *c0 = new TCanvas("c0", "peak Amplitudes", 800, 600);
    c0->cd();
    gStyle->SetOptStat(0);
    c0->SetLogy();
    TH1F *htemp1 = new TH1F("htemp1", " ", 100, 10, 10);
    
    tree1->Draw(Form("%s>>htemp1", varname.Data()));
    
    htemp1->SetFillColor(1);
    htemp1->SetLineWidth(3);
    htemp1->SetLineColor(1);

    TH1F *htemp2 = new TH1F("htemp2", " ", 100, 10, 10);
    
    tree2->Draw(Form("%s>>htemp2", varname.Data()));
    
    htemp2->SetFillColor(kOrange-2);
    htemp2->SetLineWidth(3);
    htemp2->SetLineColor(kOrange-2);

    htemp2->GetXaxis()->SetTitle("Deposited Energy [a.u.]");
    htemp2->GetYaxis()->SetTitle("Entries");

    htemp2->Draw("");
    htemp1->Draw("same");

    //leg->AddEntry(htemp1, "background", "f");
    //leg->AddEntry(htemp2, "irradiated crystal", "f");
    leg->AddEntry(htemp1, label1, "f");
    leg->AddEntry(htemp2, label2, "f");

    leg->Draw("same");

    //TH1F *htemp3 = new TH1F("htemp3", " ", 100, 0, 0.05);
    
    //tree2->Draw(Form("%s>>htemp3", varname.Data()));

    TPaveText *pt = new TPaveText(0.49,0.66,0.90,0.70, "NBNDC");
    pt->SetFillColor(0);
    pt->AddText("Threshold at 15mV - CERN data 2025");
    pt->SetBorderSize(0);
    pt->Draw();
}