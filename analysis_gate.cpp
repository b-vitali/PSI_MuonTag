/*
--------------------------------------------------------------------------------------
------------------ Make stack plots from TTrees in different files -------------------
--------------------------------------------------------------------------------------
*/

#include <stdlib.h>
#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TChain.h"
#include "TStyle.h"
#include "THStack.h"
#include "TH1.h"
#include "TH2.h"
#include <TMath.h>
#include "TString.h"
#include <TStyle.h>

/*
Choose:
- colors
- main path
- name of the folder
- file names
- what plots you would like to have
- do you want to apply any cuts?
- range for those plots
- fill with one every N TTree
- some printout to see what is happening
*/

// int my_colors[] = {kBlue, kRed, kViolet, kOrange, kYellow, kMagenta, kPink, kSpring, kGreen, kTeal, kCyan, kAzure};
int my_colors[] = {30, 40, 31, 41, 32, 42, 33, 43, 34, 44, 35, 45, 36, 46, 37, 47, 38, 48, 39, 49};

TString path="/home/bastiano/Documents/Geant4/PSI/insertion/data/";

TString folder="28MeV_1k";

TString file_names[] = {
        "28MeV_0.01.root"
    ,   "28MeV_0.02.root"
    ,   "28MeV_0.03.root"
    ,   "28MeV_0.04.root"
    ,   "28MeV_0.05.root"
    ,   "28MeV_0.06.root"
    ,   "28MeV_0.07.root"
    ,   "28MeV_0.08.root"
    ,   "28MeV_0.09.root"
    ,   "28MeV_0.1.root"
    ,   "28MeV_0.2.root"
    ,   "28MeV_0.3.root"
    ,   "28MeV_0.4.root"
    ,   "28MeV_0.5.root"
    ,   "28MeV_0.6.root"
    ,   "28MeV_0.7.root"
    ,   "28MeV_0.8.root"
    ,   "28MeV_0.9.root"
    ,   "28MeV_1.root"
};

double_t thickness[] = {
        0.01
    ,   0.02
    ,   0.03
    ,   0.04
    ,   0.05
    ,   0.06
    ,   0.07
    ,   0.08
    ,   0.09
    ,   0.1
    ,   0.2
    ,   0.3
    ,   0.4
    ,   0.5
    ,   0.6
    ,   0.7
    ,   0.8
    ,   0.9
    ,   1.
};

// these are the commands you would give to TTree->Draw() with the branch names
std::vector <char *> plots {
        (char*)"fEin"
    ,   (char*)"fEout"
    ,   (char*)"fEdep"
    ,   (char*)"fThetaOut"
    ,   (char*)"fTrackLength"
    ,   (char*)"fNgamma"
};

// in principle you can add also some cuts, on the same or on a different branch
std::vector <char *> cuts {
        (char*)""
    ,   (char*)""
    ,   (char*)""
    ,   (char*)""
    ,   (char*)""
    ,   (char*)""
};

std::vector <char *> ranges {
        (char*)"(100,0,40)"
    ,   (char*)"(400,0,40)"
    ,   (char*)"(250,0,5)"
    ,   (char*)"(300,0,30)"
    ,   (char*)"(1000,0,1.3)"
    ,   (char*)"(1000,0,100000)"
}; 

/*
basically the first plot will be:
    the variable "currentleft/currentback"; requiring ""thetapositron>0.01; in range 0-0.5 with 100 bin
*/

int skim = 1;

bool debug = false;

//* Some things are char * instead of strings because I think they work better with TTree->Draw()
//* To improve might be nice to have one obj with both names and ranges for the plots, something like
/*
std::vector <std::pair<char*, char*>> ppp {
    {(char*)"currentleft/currentback",  (char*)"(100,0,0.5)"},
    {(char*)"edep",                     (char*)"(100,0,0.5)"}
};
*/

/*
--------------------------------------------------------------------------------------
------------------ From here on everything should be automatic -----------------------
--------------------------------------------------------------------------------------
*/

// creates legend given a vector of histograms
void fai_legenda(TLegend *legendina, std::vector<TH1F *> h){
    for(int i=0; i<h.size();i++){
        if(i%skim==0) legendina->AddEntry(h[i],file_names[i],"lf");    
    }
}

// set aspect of histogram in a vector
void color_histos(std::vector<TH1F *> h_v){
    for(int i=0; i<h_v.size();i++){
        h_v[i]->SetLineColor(my_colors[i]);
        h_v[i]->SetLineWidth(3);
        h_v[i]->SetFillColor(my_colors[i]);
        h_v[i]->SetFillStyle(3002);
    }
}

void test(){
    cout<<"Data from the files in path : "<<path<<endl<<endl;

    std::vector<TTree *> tree_v;

    std::vector<std::vector<TH1F *>> h1;

    std::vector<THStack *> sh_v; 

    TH1F * h1_tmp;
    TFile *file_tmp;
    std::vector<TGraphErrors *> gr_v;

    // open the files and collect the TTrees in a vector
    for(const TString &file_name : file_names){
        file_tmp =  TFile::Open(path+folder+"/"+file_name);
        tree_v.push_back(file_tmp->Get<TTree>("Scint_gate"));
    }
    
    if(debug) cout<<"file importend and vector<TTree *> filled"<<endl<<endl;

    // loop on every plot you want
    for(int j = 0; j < plots.size(); j++){

        std::vector<TH1F *> v_tmp;
        sh_v.push_back(new THStack);

        gr_v.push_back(new TGraphErrors);
        gr_v[j]->SetTitle(folder+"/"+TString::Format("%s", plots[j]));
        gr_v[j]->SetName("gr_"+TString::Format("%s", plots[j]));

        // for every plot loop on all the TTrees
        for (int i=0; i<tree_v.size();i++)
        {
            // just a check on which file we are right now
            if(debug) {
                std::cout<<file_names[i]<<endl;
                std::cout<<TString::Format("%s>>h1_%d%d%s", plots[j], j,i,ranges[j])<<endl;
            }

            /*
                if you want to see all the plots un-comment new TCanvas and "goff"->"" in the Draw()
            */

            //new TCanvas("",file_names[i]);
            tree_v[i]->Draw(TString::Format("%s>>h1_%d%d%s", plots[j], j,i,ranges[j]),cuts[j],"goff");
            h1_tmp = (TH1F*)gDirectory->Get(TString::Format("h1_%d%d",j,i));

            v_tmp.push_back(h1_tmp);
            
            // make a tgrapherror
            gr_v[j]->SetPoint(gr_v[j]->GetN(),thickness[i],h1_tmp->GetMean());
            gr_v[j]->SetPointError(gr_v[j]->GetN()-1,0,h1_tmp->GetStdDev()/sqrt(h1_tmp->GetEntries())); //

            // just a check on the lenght of the two arrays
            if(debug) {
                std::cout<<"Histogram number: "<<v_tmp.size()<<std::endl;
                std::cout<<"Stack number: "<<sh_v.size()<<std::endl<<std::endl;
            }

            if(i%skim==0) sh_v[j]->Add(v_tmp[i]);
        }
        
        h1.push_back(v_tmp);

    }
    
    std::vector<TLegend*> legende;    

    // make a canvas for every stack and put proper legend
    for(int j=0; j<h1.size();j++){
        color_histos(h1[j]);
        new TCanvas("",folder);
        sh_v[j]->SetTitle(folder+"/"+TString::Format("%s", plots[j]));
        sh_v[j]->SetName("sh_"+TString::Format("%s", plots[j]));
        sh_v[j]->Draw("ehist nostack");
        legende.push_back(new TLegend(0.83,0.3,0.98,0.7));
        fai_legenda(legende[j],h1[j]);
        legende[j]->Draw();

        new TCanvas("",folder);
        gr_v[j]->Draw("AP0");
    }

    TFile * output_file = new TFile(path+folder+"_graphs.root", "recreate");
    for(int j=0; j<h1.size();j++){
        gr_v[j]->Write();
    }
    for(int j=0; j<h1.size();j++){
        sh_v[j]->Write();
    }
    output_file->Close();

}

void make_graphs(std::vector<TTree *> tree_v, TString nani){
    TGraphErrors * gr = new TGraphErrors();
    for (int i=0; i<tree_v.size();i++){
        gr->SetPoint(gr->GetN(), 1,1);
    }
   gr->Draw("ALP");
}