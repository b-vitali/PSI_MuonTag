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

int my_colors[] = {kBlue, kRed, kViolet, kOrange, kYellow, kMagenta, kPink, kSpring, kGreen, kTeal, kCyan, kAzure};

THStack * sh = new THStack;
THStack * sh2 = new THStack;


TString path="/home/bastiano/Documents/Geant4/PSI/A_NewMuonTag/build/";

TString folder="mu+28"; //mu+28

TString file_names[] = {
        "0.05mm.root"
    ,   "0.1mm.root"
    ,   "0.25mm.root"
    ,   "0.5mm.root"
};

void fai_legenda(TLegend *legendina, std::vector<TH1F *> h){ //crea la legenda per i vari sample_type

    for(int i=0; i<h.size();i++){
        legendina->AddEntry(h[i],file_names[i],"l");    
    }

}

void test(){
    cout<<"Data from the files in path : "<<path<<endl;
    
    std::vector<TTree *> tree_v;
    std::vector<TH1F *> h1_v1;
    std::vector<TH1F *> h1_v2;   //TString::Format("h_%s", file_names[i])

    TH1F * h_tmp;
    
    TFile *file_tmp;
    for(const TString &file_name : file_names){
        file_tmp =  TFile::Open(path+folder+"/"+file_name);
        tree_v.push_back(file_tmp->Get<TTree>("T"));
        //cout<<file_name.Data()<<endl;
    }
    cout<<"ok"<<endl;

    for(int i=0; i<tree_v.size();i++){
        
        //cout<<tree_v[i]<<endl; 
        //tree_v[i]->Print();
        
        new TCanvas("",folder+file_names[i]);
        tree_v[i]->Draw("currentleft:currentback","","colz");

        new TCanvas("",file_names[i]);
        tree_v[i]->Draw(TString::Format("thetapositron>>hL_%d",i));
        h_tmp = (TH1F*)gDirectory->Get(TString::Format("hL_%d",i)); 
        h1_v2.push_back(h_tmp);
        sh2->Add(h1_v2[i]);


        tree_v[i]->Draw(TString::Format("currentleft/currentback>>h_%d",i),"","goff");
        h_tmp = (TH1F*)gDirectory->Get(TString::Format("h_%d",i)); 
        h1_v1.push_back(h_tmp);
        sh->Add(h1_v1[i]);
        
        new TCanvas("",folder+file_names[i]);
        tree_v[i]->Draw("currentleft:currentright:currentup:currentdown","","candle");
    }

    for(int i=0; i<tree_v.size();i++){
        h1_v1[i]->SetLineColor(my_colors[i]);
        h1_v2[i]->SetLineColor(my_colors[i]);
    }
    
    new TCanvas("",folder);
    sh->Draw("nostack");
    TLegend *legendina = new TLegend(0.83,0.3,0.98,0.7);
    fai_legenda(legendina,h1_v1);
    legendina->Draw();

    new TCanvas("",folder);
    sh2->Draw("nostack");
    TLegend *legendina2 = new TLegend(0.83,0.3,0.98,0.7);
    fai_legenda(legendina2,h1_v2);
    legendina2->Draw();
}