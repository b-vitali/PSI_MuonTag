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

// https://root.cern.ch/doc/master/classTColor.html#C06
// int my_colors[] = {kBlue, kRed, kViolet, kOrange, kYellow, kMagenta, kPink, kSpring, kGreen, kTeal, kCyan, kAzure};
int my_colors[] = {30, 40, 31, 41, 32, 42, 33, 43, 34, 44, 35, 45, 36, 46, 37, 47, 38, 48, 39, 49};

TString path="/home/bastiano/Documents/Geant4/PSI/insertion/data/";

TString folder="new1k"; //28MeV_10k2

TString energy = "28MeV"; //125MeV

std::vector<TString> file_names;

std::vector<double_t> thickness{
        0.01
    ,   0.02
    ,   0.03
    ,   0.04
    ,   0.05
    ,   0.06
    ,   0.07
    ,   0.08
    ,   0.09
    ,   0.10
    ,   0.11
    ,   0.12
    ,   0.13
    ,   0.14
    ,   0.15
    ,   0.16
    ,   0.17
    ,   0.18
    ,   0.19
    ,   0.20
    ,   0.30
    ,   0.40
    ,   0.50
    ,   0.60
    ,   0.70
    ,   0.80
    ,   0.90
    ,   1.00
};

// these are the commands you would give to TTree->Draw() with the branch names
// Choose *variable* *some cut* *histogram range*
std::vector< std::tuple<char*, char*, char*> > plots {
        {   (char*)"fEin",          (char*)"",  (char*)"(300,0,61)"     }
    ,   {   (char*)"fEout",         (char*)"",  (char*)"(300,0,5)"      }
    ,   {   (char*)"fEdep",         (char*)"",  (char*)"(500,0,5)"      }
    ,   {   (char*)"fThetaOut",     (char*)"",  (char*)"(250,0,100)"    }
    ,   {   (char*)"fTrackLength",  (char*)"",  (char*)"(600,0,1.2)"    }
    ,   {   (char*)"fNgamma",       (char*)"",  (char*)"(100,0,1000000)"}
    ,   {   (char*)"fFront",        (char*)"",  (char*)"(1000,0,100000)"}

};

/*
basically the first plot will be:
    the variable "currentleft/currentback"; requiring ""thetapositron>0.01; in range 0-0.5 with 100 bin
*/

int skim = 1;

bool debug = false;

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
        // h_v[i]->SetLineColor(my_colors[i]);
        h_v[i]->SetLineWidth(3);
        // h_v[i]->SetFillColor(my_colors[i]);
        h_v[i]->SetFillStyle(3002);
    }
}

void test(){
    cout<<"Data from the files in path : "<<path<<endl<<endl;
    
    std::stringstream ss;
    for(int i = 0; i< thickness.size(); i++){
        ss << std::fixed<<std::setprecision(2)<<thickness[i];
        TString str = ss.str();

        file_names.push_back(energy+"_"+ str+".root");
        std::cout<<file_names[i]<<std::endl;
        
        ss.str(std::string());
    }

    std::vector<TTree *> tree_v;

    std::vector<std::vector<TH1F *>> h1;

    std::vector<THStack *> sh_v; 

    TH1F * h1_tmp;
    TFile *file_tmp;
    std::vector<TGraphErrors *> gr_v;

    // open the files and collect the TTrees in a vector
    for(const TString &file_name : file_names){
        file_tmp =  TFile::Open(path+folder+"/"+file_name);
        tree_v.push_back(file_tmp->Get<TTree>("Scint_telescope"));
    }
    
    if(debug) cout<<"file importend and vector<TTree *> filled"<<endl<<endl;
    if(debug) cout<<"plot size() = "<<plots.size()<<endl;
    if(debug) cout<<"tree_v size() = "<<tree_v.size()<<endl;

    // loop on every plot you want
    for(int j = 0; j < plots.size(); j++){

        std::vector<TH1F *> v_tmp;
        sh_v.push_back(new THStack);

        gr_v.push_back(new TGraphErrors);
        gr_v[j]->SetTitle(folder+"/"+TString::Format("%s", std::get<0>(plots[j]) ));
        gr_v[j]->SetName("gr_"+TString::Format("%s", std::get<0>(plots[j]) ));

        // for every plot loop on all the TTrees
        for (int i=0; i<tree_v.size();i++)
        {
            // just a check on which file we are right now
            if(debug) {
                std::cout<<file_names[i]<<endl;
                std::cout<<TString::Format("%s>>h1_%d%d%s", std::get<0>(plots[j]), j,i,std::get<2>(plots[j]))<<endl;
            }

            /*
                if you want to see all the plots un-comment new TCanvas and "goff"->"" in the Draw()
            */

            if(tree_v[i]->GetEntries()==0) return;

            //new TCanvas("",file_names[i]);
            tree_v[i]->Draw(TString::Format("%s>>h1_%d%d%s", std::get<0>(plots[j]), j,i,std::get<2>(plots[j])),std::get<1>(plots[j]),"goff");
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

    gStyle->SetPalette(kRainBow);

    // make a canvas for every stack and put proper legend
    for(int j=0; j<h1.size();j++){
        color_histos(h1[j]);
        new TCanvas("",folder);
        sh_v[j]->SetTitle(folder+"/"+TString::Format("%s", std::get<0>(plots[j]) ));
        sh_v[j]->SetName("sh_"+TString::Format("%s", std::get<0>(plots[j]) ));
        sh_v[j]->Draw("ehist nostack pfc plc");
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