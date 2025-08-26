#include <TH1.h>
#include <TFitResult.h>

TTree* GetTree(TString filename){

    TString PATH = "/Users/bordonis/ResearchActivities/PerovskiteAnalysis/PeroAna/rootinputfiles/";
    TTree *tree_0;

    if (filename){
        TFile *file_0= new TFile(PATH+filename, "OPEN");
        tree_0 = (TTree*)file_0->Get("flat_tree");
    }
    else
        tree_0;

    return tree_0;
}


Color_t GetChannelColor(int ich) {
    static Color_t colors[] = { kRed, kBlue, kGreen+2, kMagenta, kCyan+2, kOrange };
    int Ncolors = sizeof(colors)/sizeof(colors[0]);
    return colors[ich % Ncolors];
}


int GetNChannels(TTree* tree) {
    if (!tree) return 0;

    tree->Draw("channel", "", "goff");
    int n = tree->GetSelectedRows();
    double* ch_array = tree->GetV1();

    int max_channel = -1;
    for (int i = 0; i < n; ++i) {
        int ch = static_cast<int>(ch_array[i]);
        if (ch > max_channel)
            max_channel = ch;
    }

    return max_channel + 1;
}

void ExtractSignal(TString f_1 ="run1",
                    TString f_2 ="run2",
                    TString label1 = "label1 ",
                    TString label2 = "label2 ")
{

    TTree *tree1;
    TTree *tree2;

    Bool_t tree1_ok = false;
    Bool_t tree2_ok = false;


     if (! f_1.IsNull()) {
        //TTree *tree1 = GetTree(f_1);
        tree1 = GetTree(f_1);
        tree1->Print();
        //tree1->SetLineColor(1);
        //tree1->SetLineWidth(3);
        tree1_ok = true;
        //tree1->SetLineStyle(3);
    }

    if (! f_2.IsNull()) {
        tree2 = GetTree(f_2);
        tree2->Print();
        //tree2->SetLineColor(4);
        //tree2->SetLineWidth(3);
        tree2_ok = true;
    }


    std::vector<TH1F*> h_tree1_pulseInt;

    std::vector<TH1F*> h_tree2_pulseInt;

    std::vector<TH1F*> h_diff_pulseInt;

    if (tree1_ok) {
        int nCh1 = GetNChannels(tree1);
        for (int ich = 0; ich < nCh1; ++ich) {
            TString hname = Form("h_tree1_pulseInt_ch%d", ich);
            TH1F *h = new TH1F(hname, Form("pulseInt - Tree1 - Ch %d", ich),
                            100, 0, 5.);
            tree1->Draw(Form("pulseInt >> %s", hname.Data()), 
                        Form("channel == %d", ich), 
                        "goff");
            h->SetLineColor(GetChannelColor(ich));
            h->SetLineStyle(1); // solid
            h->SetLineWidth(2);
            h_tree1_pulseInt.push_back(h);
        }   
    }


    if (tree2_ok) {
        int nCh2 = GetNChannels(tree2);
        for (int ich = 0; ich < nCh2; ++ich) {
            TString hname = Form("h_tree2_pulseInt_ch%d", ich);
            TH1F *h = new TH1F(hname, Form("pulseInt - Tree2 - Ch %d", ich),
                            100, 0, 5.);
            tree2->Draw(Form("pulseInt >> %s", hname.Data()), 
                        Form("channel == %d", ich), 
                        "goff");
            h->SetLineColor(GetChannelColor(ich));
            h->SetLineStyle(2); // dashed
            h->SetLineWidth(2);
            h_tree2_pulseInt.push_back(h);
        }   
    }

    if (tree1_ok && tree2_ok) {
        int nCh = std::min(h_tree1_pulseInt.size(), h_tree2_pulseInt.size());
        for (int ich = 0; ich < nCh; ++ich) {
            TString hname = Form("h_ratio_pulseInt_ch%d", ich);
            TH1F *h_diff = (TH1F*)h_tree2_pulseInt[ich]->Clone(hname);
            h_diff->SetTitle(Form("Ratio pulseInt - Ch %d", ich));
            h_diff->Add(h_tree1_pulseInt[ich], -1);
            h_diff->SetLineColor(GetChannelColor(ich));
            h_diff->SetLineStyle(1); // solid
            h_diff->SetLineWidth(2);
            h_diff_pulseInt.push_back(h_diff);
        }
    }


    // Now create canvases to display the histograms
    std::vector<TCanvas*> canvases;
    for (int ich = 0; ich < std::max(h_tree1_pulseInt.size(), h_tree2_pulseInt.size()); ++ich) {
        TString h_cname = Form("c_ch%d", ich);
        TCanvas *c = new TCanvas(h_cname, Form("pulseInt comparison - Ch %d", ich), 800, 600);         
        canvases.push_back(c);
    }

    bool firstHist = true;
    TLegend *leg = new TLegend(0.6,0.6,0.88,0.88);

    std::vector<TF1*> fit_funcs_tree1_gaus;
    std::vector<TF1*> fit_funcs_tree2_gaus;
   
    int nCh = std::min(h_tree1_pulseInt.size(), h_tree2_pulseInt.size());
    for (int ich = 0; ich < nCh; ++ich) {
        TString f_name_tree1 = Form("f_tree1_gaus_ch%d", ich);
        TString f_name_tree2 = Form("f_tree2_gaus_ch%d", ich);
        TF1 *f_tree1 = new TF1(f_name_tree1, "gaus", 0, 5);
        TF1 *f_tree2 = new TF1(f_name_tree2, "gaus", 0, 5);
        fit_funcs_tree1_gaus.push_back(f_tree1);
        fit_funcs_tree2_gaus.push_back(f_tree2);
    }
        

    for (int ch = 0; ch < std::max(h_tree1_pulseInt.size(), h_tree2_pulseInt.size()); ++ch) {
        TCanvas *c_ch = canvases[ch];
        c_ch->cd();
        std::cout << "Drawing channel " << ch << std::endl;
        if (tree1_ok && ch < h_tree1_pulseInt.size()) {
            h_tree1_pulseInt[ch]->Draw(firstHist ? "hist" : "hist same");
            leg->AddEntry(h_tree1_pulseInt[ch], Form("%s - Ch%d", label1.Data(), ch), "l");
            firstHist = false;  

            fit_funcs_tree1_gaus[ch]->SetLineColor(h_tree1_pulseInt[ch]->GetLineColor()+2);
            TFitResultPtr r = h_tree1_pulseInt[ch]->Fit(fit_funcs_tree1_gaus[ch],"S");
           
            
            fit_funcs_tree1_gaus[ch]->Draw("same");
            std::cout << endl;
            std::cout << "Tree1 - Ch " << ch << ": mean = " << fit_funcs_tree1_gaus[ch]->GetParameter(1) << " ± " << fit_funcs_tree1_gaus[ch]->GetParError(1) 
                      << ", sigma = " << fit_funcs_tree1_gaus[ch]->GetParameter(2) << " ± " << fit_funcs_tree1_gaus[ch]->GetParError(2) << std::endl;
            std::cout << endl;

        }


        if (tree2_ok && ch < h_tree2_pulseInt.size()) {
            h_tree2_pulseInt[ch]->Draw(firstHist ? "hist" : "hist same");
            leg->AddEntry(h_tree2_pulseInt[ch], Form("%s - Ch%d", label2.Data(), ch), "l");
            firstHist = false;

            fit_funcs_tree2_gaus[ch]->SetLineColor(h_tree2_pulseInt[ch]->GetLineColor()+2);
            TFitResultPtr r = h_tree2_pulseInt[ch]->Fit(fit_funcs_tree2_gaus[ch],"S");            
            fit_funcs_tree2_gaus[ch]->Draw("same");

            std::cout << endl;
            std::cout << "Tree2 - Ch " << ch << ": mean = " << fit_funcs_tree2_gaus[ch]->GetParameter(1) << " ± " << fit_funcs_tree2_gaus[ch]->GetParError(1) 
                      << ", sigma = " << fit_funcs_tree2_gaus[ch]->GetParameter(2) << " ± " << fit_funcs_tree2_gaus[ch]->GetParError(2) << std::endl;
            std::cout << endl;

        }

        leg->Draw();
        c_ch->Update();
   }

    

   /*
    TCanvas *c = new TCanvas("c_ratio", "ratio ", 800, 1200);
    c->Divide(2, 1);  
    c->cd(1);       
    h_diff_pulseInt[0]->Draw("hist");

    c->cd(2);       
    h_diff_pulseInt[1]->Draw("hist");
    c->Update();*/



}