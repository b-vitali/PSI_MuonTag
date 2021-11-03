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

//int my_colors[] = {kBlue, kRed, kViolet, kOrange, kYellow, kMagenta, kPink, kSpring, kGreen, kTeal, kCyan, kAzure};

int my_colors[] = {30, 40, 31, 41, 32, 42, 33, 43, 34, 44, 35, 45, 36, 46, 37, 47, 38, 48, 39, 49};

TString path="/home/bastiano/Documents/Geant4/PSI/A_NewMuonTag/build/data/";

/*
TString folder="e+28_1k_300x50xZ"; //mu+28 mu+28_1k e+28_5k_50x50xZ

TString file_names[] = {
        "1mm.root"
    ,   "2mm.root"
    ,   "5mm.root"
    ,   "10mm.root"
};
*/

TString folder="e+28_5k_50x50xZ"; //mu+28 mu+28_1k e+28_5k_50x50xZ

TString file_names[] = {
        "0.025mm.root"
    ,   "0.05mm.root"
    ,   "0.075mm.root"
    ,   "0.1mm.root"
    ,   "0.15mm.root"
    ,   "0.2mm.root"
    ,   "0.25mm.root"
    ,   "0.5mm.root"
};

float thickness[] = {
        0.025
    ,   0.05
    ,   0.075
    ,   0.1
    ,   0.15
    ,   0.2
    ,   0.25
    ,   0.5
};

/*
TString folder="e+E_5k_50x50x2"; //mu+28 mu+28_1k

TString file_names[] = {
        "2MeV.root"
    ,   "28MeV.root"
    ,   "52.8MeV.root"
    ,   "128MeV.root"
};
*/
std::vector<THStack *> sh_v; 

THStack * sh = new THStack;
THStack * sh2 = new THStack;
THStack * sh3 = new THStack;
THStack * sh4 = new THStack;


char * plot_1 = (char*)"currentleft/currentback"; // 
char * plot_2 = (char*)"thetapositron";
char * plot_3 = (char*)"eout";
char * plot_4 = (char*)"edep";


void fai_legenda(TLegend *legendina, std::vector<TH1F *> h){ //crea la legenda per i vari sample_type

    for(int i=0; i<h.size();i++){
        legendina->AddEntry(h[i],file_names[i],"lf");    
    }

}

void color_histos(std::vector<TH1F *> h_v){
    for(int i=0; i<h_v.size();i++){
        h_v[i]->SetLineColor(my_colors[i]);
        h_v[i]->SetLineWidth(3);
        h_v[i]->SetFillColor(my_colors[i]);
        h_v[i]->SetFillStyle(3002);
    }
}

void test(){
    cout<<"Data from the files in path : "<<path<<endl;
    
    std::vector<TTree *> tree_v;
    std::vector<TH1F *> h1_v1;
    std::vector<TH1F *> h1_v2;   
    std::vector<TH1F *> h1_v3;   
    std::vector<TH1F *> h1_v4;   

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
        
        new TCanvas("",folder+"-"+file_names[i]);
        tree_v[i]->Draw("currentleft:edep","","colz");

        new TCanvas("",file_names[i]);
        tree_v[i]->Draw(TString::Format("%s>>h1_%d(100,0,1)", plot_1, i),"",""); //(500,0,2000)
        h_tmp = (TH1F*)gDirectory->Get(TString::Format("h1_%d",i)); 
        h1_v1.push_back(h_tmp);
        sh->Add(h1_v1[i]);

        new TCanvas("",file_names[i]);
        tree_v[i]->Draw(TString::Format("%s>>h2_%d(160*10,0,1.6)", plot_2, i),"","gott");
        h_tmp = (TH1F*)gDirectory->Get(TString::Format("h2_%d",i)); 
        h1_v2.push_back(h_tmp);
        sh2->Add(h1_v2[i]);

        //new TCanvas("",file_names[i]);
        tree_v[i]->Draw(TString::Format("%s>>h3_%d(200,0,4)", plot_3, i),"","goff"); //(5000,0,0.5)
        h_tmp = (TH1F*)gDirectory->Get(TString::Format("h3_%d",i)); 
        h1_v3.push_back(h_tmp);
        sh3->Add(h1_v3[i]);

        tree_v[i]->Draw(TString::Format("%s>>h4_%d(500,0,4)", plot_4, i),"","goff"); //(5000,0,0.5)
        h_tmp = (TH1F*)gDirectory->Get(TString::Format("h4_%d",i)); 
        h1_v4.push_back(h_tmp);
        sh4->Add(h1_v4[i]);

//        new TCanvas("",folder+"-"+file_names[i]);
//        tree_v[i]->Draw("currentleft:currentright:currentup:currentdown","","candle");
//        tree_v[i]->Draw("log10(currentleft):log10(currentright):log10(currentup):log10(currentdown)","","candle");
    }

    color_histos(h1_v1);
    new TCanvas("",folder);
    sh->SetTitle("Number of #gamma left side"); //plot_1
    sh->Draw("nostack");
    TLegend *legendina = new TLegend(0.83,0.3,0.98,0.7);
    fai_legenda(legendina,h1_v1);
    legendina->Draw();

    color_histos(h1_v2);
    new TCanvas("",folder);
    sh2->SetTitle(plot_2);
    sh2->Draw("nostack");
    TLegend *legendina2 = new TLegend(0.83,0.3,0.98,0.7);
    fai_legenda(legendina2,h1_v2);
    legendina2->Draw();

    color_histos(h1_v3);
    new TCanvas("",folder);
    sh3->SetTitle(plot_3);
    sh3->Draw("nostack");
    TLegend *legendina3 = new TLegend(0.83,0.3,0.98,0.7);
    fai_legenda(legendina3,h1_v3);
    legendina3->Draw();

    color_histos(h1_v4);
    new TCanvas("",folder);
    sh4->SetTitle(plot_4);
    sh4->Draw("nostack");
    TLegend *legendina4 = new TLegend(0.83,0.3,0.98,0.7);
    fai_legenda(legendina4,h1_v4);
    legendina4->Draw();

}

Double_t angle(Double_t *x, Double_t *par)
{
    Double_t m  = par[0];
    Double_t p  = par[1];
    Double_t X0 = par[2];

    Double_t xx     = x[0];
    Double_t E      = sqrt(p*p*+m*m);
    Double_t gamma  = E/m;
    Double_t beta   = p/E;

    Double_t f = ( 13.6/(beta*p) ) * sqrt(xx/X0) * ( 1+0.038*log( xx / (X0*beta*beta) ) );
    return f;
}

void multiplescattering(){
        cout<<"Data from the files in path : "<<path<<endl;
    
    std::vector<TTree *> tree_v;
    std::vector<TH1F *> h1_v1;   

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

        //new TCanvas("",file_names[i]);
        tree_v[i]->Draw(TString::Format("thetapositron>>h1_%d(1600,0,1.6)", i),"","goff"); //(500,0,2000)
        h_tmp = (TH1F*)gDirectory->Get(TString::Format("h1_%d",i)); 
        h1_v1.push_back(h_tmp);
        sh->Add(h1_v1[i]);
    }

    color_histos(h1_v1);
    new TCanvas("",folder);
    sh->SetTitle("Angle exit"); //plot_1
    sh->Draw("nostack");
    TLegend *legendina = new TLegend(0.83,0.3,0.98,0.7);
    fai_legenda(legendina,h1_v1);
    legendina->Draw();

    TGraphErrors *g = new TGraphErrors();
    g->SetMarkerStyle(21);
    for(int i=0; i<tree_v.size();i++){
        g->SetPoint(g->GetN(),thickness[i]*0.001,h1_v1[i]->GetMean()); // h1_v1[i]->GetStdDev()/h1_v1[i]->GetEntries(),thickness[i]/1000.
        g->SetPointError(i,0.000001,h1_v1[i]->GetStdDev()/h1_v1[i]->GetEntries());
    }

    TF1 * f = new TF1("Scattering","angle",0,10,3);
    f->SetParameters(1,1,1);
    f->FixParameter(0,0.511); //0.511 105.66
    f->FixParameter(1,28);
    new TCanvas;
    f->Draw();
    
    new TCanvas;
    g->Fit(f,"EM","",0,0.5);
    g->Draw("AP");
}