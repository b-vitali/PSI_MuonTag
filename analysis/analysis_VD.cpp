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

TString folder="28MeV_10k2"; //28MeV_10k2 new1k

TString energy = "28MeV"; //125MeV

int skim = 1;

bool debug = false;

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

//############################################################
//####################### 1D Plots ###########################
//############################################################

// these are the commands you would give to TTree->Draw() with the branch names
// Choose *variable* *some cut* *histogram range*
std::vector< std::tuple<char*, char*, char*> > plots {
        {   (char*)"fParticleID",   (char*)"fVDNo == 1 && fParticleID == -13",  (char*)"(200,-100,+100)"}
    // ,   {   (char*)"fMom",          (char*)"fVDNo == 1",  (char*)"(1300,0, 130 )"   }
    // ,   {   (char*)"fPosX",         (char*)"fVDNo == 1 && fParticleID == -13",  (char*)"(880,-55,55)" }
    // ,   {   (char*)"fPosY",         (char*)"fVDNo == 1 && fParticleID == -13",  (char*)"(880,-55,55)" }
    // ,   {   (char*)"fMomX",         (char*)"fVDNo == 1",  (char*)"(200,-5,5)"     }
    // ,   {   (char*)"fMomY",         (char*)"fVDNo == 1",  (char*)"(200,-5,5)"     }
    ,   {   (char*)"fMomZ",         (char*)"fVDNo == 1 && fParticleID == -13",  (char*)"(500,10,30)"     }

        // {   (char*)"fParticleID",   (char*)"fVDNo == 1 && fParticleID == -13",  (char*)"(200,-100,+100)"    }
    // ,   {   (char*)"fMom",          (char*)"fVDNo == 1 && fParticleID == -13",  (char*)"(300,0, 30 )"    }
    // ,   {   (char*)"fPosX",         (char*)"fVDNo == 1 && fParticleID == -13",  (char*)"(880,-55,55)"       }
    // ,   {   (char*)"fPosY",         (char*)"fVDNo == 1 && fParticleID == -13",  (char*)"(880,-55,55)"       }
    // ,   {   (char*)"fMomX",         (char*)"fVDNo == 1 && fParticleID == -13",  (char*)"(200,-5,5)"   }
    // ,   {   (char*)"fMomY",         (char*)"fVDNo == 1 && fParticleID == -13",  (char*)"(200,-5,5)" }
    // ,   {   (char*)"fMomZ",         (char*)"fVDNo == 1 && fParticleID == -13",  (char*)"(300,126,129)"     }

};

//############################################################
//####################### 2D Plots ###########################
//############################################################

// these are the commands you would give to TTree->Draw() with the branch names
std::vector <char *> plots2D {
    //     (char*)"fMomX:fPosX"
    // ,   (char*)"fMomY:fPosY"
};

// in principle you can add also some cuts, on the same or on a different branch
std::vector <char *> cuts2D {
        // (char*)""
    // ,   (char*)""


        (char*)"fVDNo == 1 && fParticleID == -13"
    ,   (char*)"fVDNo == 1 && fParticleID == -13"

};

std::vector <char *> ranges2D {
        // (char*)"(300,58,61)"
    // ,   (char*)"(300,58,61)"

        (char*)"(100,-20,20, 100,-5,5)"
    ,   (char*)"(100,-20,20, 100,-5,5)"
}; 

double_t covariance(std::vector<double_t> v1, std::vector<double_t> v2){
if(debug) std::cout<< "covariance function with v1.size = "<<v1.size()<<std::endl;

    if(v1.size() != v2.size()) std::cout<< "must have same size"<<std::endl;

    double mean1 = std::accumulate(v1.begin(), v1.end(), 0.0);
if(debug) std::cout<<mean1<<std::endl;
    mean1 = mean1 / v1.size();
if(debug) std::cout<<mean1<<std::endl;
    double mean2 = std::accumulate(v2.begin(), v2.end(), 0.0);
//  std::cout<<mean2<<std::endl;
   mean2 = mean2 / v2.size();
if(debug) std::cout<<mean2<<std::endl;

    double_t var = 0;
    for (int i=0; i<v1.size();i++){
        var = var + (v1[i]-mean1)*(v2[i]-mean2);
    }
if(debug) std::cout<<var<<std::endl;    
    var = var / (v1.size()-1);
if(debug) std::cout<<var<<std::endl;
    return var;
}

void emittance(){
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
    std::vector<TTree *> tree_telescope_v;
    TFile *file_tmp;

    TGraphErrors * grx = new TGraphErrors;
    grx->SetMarkerStyle(20);
    grx->SetMarkerColor(kBlue);
    grx->SetName("#varepsilon_{x}");
    grx->SetTitle("#varepsilon_{x}");
    TGraphErrors * gry = new TGraphErrors;
    gry->SetMarkerStyle(20);
    gry->SetMarkerColor(kRed);
    gry->SetName("#varepsilon_{y}");
    gry->SetTitle("#varepsilon_{y}");

    double_t x_tmp, px_tmp, y_tmp, py_tmp, pz_tmp;
    int fVDNo, fParticleID, fEvent;
    int fEvent_telescope, fNgamma;
    std::vector <double_t> x;
    std::vector <double_t> px;
    std::vector <double_t> y;
    std::vector <double_t> py;

    std::vector <double_t> pz;


    double_t emittance, q_var, p_var, co_var;

    // open the files and collect the TTrees in a vector
    for(const TString &file_name : file_names){
        file_tmp =  TFile::Open(path+folder+"/"+file_name);
        tree_v.push_back(file_tmp->Get<TTree>("VD"));
        tree_telescope_v.push_back(file_tmp->Get<TTree>("Scint_telescope"));
    }
if(debug) std::cout<<"filled vector<ttree> : "<<tree_v.size()<<std::endl;
    
    bool telescope_hit = false;
    for (int i=0; i<tree_v.size();i++){
        tree_v[i]->SetBranchAddress("fPosX",&x_tmp);
        tree_v[i]->SetBranchAddress("fMomX",&px_tmp);
        tree_v[i]->SetBranchAddress("fPosY",&y_tmp);
        tree_v[i]->SetBranchAddress("fMomY",&py_tmp);

        tree_v[i]->SetBranchAddress("fMomZ",&pz_tmp);

        tree_v[i]->SetBranchAddress("fVDNo",&fVDNo);
        tree_v[i]->SetBranchAddress("fParticleID",&fParticleID);
        tree_v[i]->SetBranchAddress("fEvent",&fEvent);

        tree_telescope_v[i]->SetBranchAddress("fEvent",&fEvent_telescope);
        tree_telescope_v[i]->SetBranchAddress("fNgamma",&fNgamma);

        for (int j=0; j<tree_v[i]->GetEntries(); j++){
            tree_v[i]->GetEntry(j);
            
            if(fVDNo == 1 && fParticleID == -13){
            
                for(int k = 0; k<tree_telescope_v[i]->GetEntries(); k++){
                    tree_telescope_v[i]->GetEntry(k);
                    if(fEvent == fEvent_telescope && fNgamma>0){
                        telescope_hit = true;
                    }
                }
                if(!telescope_hit){
                    x.push_back(x_tmp);
                    px.push_back(px_tmp);
                    y.push_back(y_tmp);
                    py.push_back(py_tmp);
                }
            }
            telescope_hit = false;
        }
        cout<<"i "<<i<<endl;
if(debug) std::cout<<"filled vectors <double_t> : "<<x.size()<<std::endl;
if(debug) std::cout<<"filled vectors <double_t> : "<<x[0]<<" "<<px[0]<<" "<<y[0]<<" "<<py[0]<<std::endl;

        if(x.size() != 0) {
if(debug) std::cout<<"!empty"<<std::endl;

        q_var = covariance(x,x);
        p_var = covariance(px,px);
        co_var = covariance(x,px);
if(debug) std::cout<<"var : "<<q_var<<" "<<p_var<<" "<<co_var<<std::endl;

        emittance =  q_var*p_var-co_var*co_var; //
if(debug) std::cout<<"emittance : "<<emittance<<std::endl;
        emittance = sqrt(emittance);
if(debug) std::cout<<"emittance : "<<emittance<<std::endl;

        grx->SetPoint(grx->GetN(),thickness[i],emittance);
        grx->SetPointError(grx->GetN()-1,0,0);

        q_var = covariance(y,y);
        p_var = covariance(py,py);
        co_var = covariance(y,py);
if(debug) std::cout<<"var : "<<q_var<<" "<<p_var<<" "<<co_var<<std::endl;


        emittance =  q_var*p_var-co_var*co_var ;
if(debug) std::cout<<"emittance : "<<emittance<<std::endl;
        emittance = sqrt( emittance );
if(debug) std::cout<<"emittance : "<<emittance<<std::endl;

        gry->SetPoint(gry->GetN(),thickness[i],emittance);
        gry->SetPointError(gry->GetN()-1,0,0);
        }
        x.clear(); px.clear(); y.clear(); py.clear();
    }

    new TCanvas;
    grx->Draw("AP0");
    new TCanvas;
    gry->Draw("AP0");

    TMultiGraph *mgr = new TMultiGraph;
    mgr->SetTitle("28 MeV : Emittance (x,px) and (y,py); Thickness [mm]; #varepsilon [mm*MeV/c]");
    mgr->Add(grx);
    mgr->Add(gry);
    TCanvas * c_mgr = new TCanvas;
    mgr->Draw("AP0");
    c_mgr->BuildLegend();
}

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
    std::vector<std::vector<TH2F *>> h2;

    std::vector<THStack *> sh_v; 

    TH1F * h1_tmp;
    TH2F * h2_tmp;
    TFile *file_tmp;
    std::vector<TGraphErrors *> gr_v;

    // open the files and collect the TTrees in a vector
    for(const TString &file_name : file_names){
        file_tmp =  TFile::Open(path+folder+"/"+file_name);
        tree_v.push_back(file_tmp->Get<TTree>("VD"));
    }
    
    if(debug) cout<<"file importend and vector<TTree *> filled"<<endl<<endl;

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
                std::cout<<TString::Format("%s>>h1_%d%d%s", std::get<0>(plots[j]) , j,i, std::get<2>(plots[j]) )<<endl;
            }

            /*
                if you want to see all the plots un-comment new TCanvas and "goff"->"" in the Draw()
            */

            //new TCanvas("",file_names[i]);
            tree_v[i]->Draw(TString::Format("%s>>h1_%d%d%s", std::get<0>(plots[j]) , j,i,std::get<2>(plots[j])),std::get<1>(plots[j]),""); //goff
            h1_tmp = (TH1F*)gDirectory->Get(TString::Format("h1_%d%d",j,i));
            
            v_tmp.push_back(h1_tmp);
            
            // make a tgrapherror
            if(thickness[i]<0.5) {
            gr_v[j]->SetPoint(gr_v[j]->GetN(),thickness[i],h1_tmp->GetStdDev()/h1_tmp->GetMean()); //h1_tmp->GetStdDev() h1_tmp->Integral(0,350)/h1_tmp->Integral()
            gr_v[j]->SetPointError(gr_v[j]->GetN()-1,0,0); ///sqrt(h1_tmp->GetEntries())
            }
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
        sh_v[j]->SetName("sh_"+TString::Format("%s",std::get<0>(plots[j]) ));
        sh_v[j]->Draw("ehist nostack pfc plc");
        legende.push_back(new TLegend(0.83,0.3,0.98,0.7));
        fai_legenda(legende[j],h1[j]);
        legende[j]->Draw();

        new TCanvas("",folder);
        gr_v[j]->Draw("AP0");
    }


// loop on every 2D plot you want

    double_t emittance_tmp;

    for(int j = 0; j < plots2D.size(); j++){

        std::vector<TH2F *> v_tmp2D;

        // for every plot loop on all the TTrees
        for (int i=0; i<tree_v.size();i++)
        {
            // just a check on which file we are right now
            if(debug) {
                std::cout<<file_names[i]<<endl;
                std::cout<<TString::Format("%s>>h1_%d%d%s", plots2D[j], j,i,ranges2D[j])<<endl;
            }

            /*
                if you want to see all the plots un-comment new TCanvas and "goff"->"" in the Draw()
            */

            new TCanvas("",file_names[i]);
            tree_v[i]->Draw(TString::Format("%s>>h2_%d%d%s", plots2D[j], j,i,ranges2D[j]),cuts2D[j],"COLZ");
            h2_tmp = (TH2F*)gDirectory->Get(TString::Format("h2_%d%d",j,i));

            v_tmp2D.push_back(h2_tmp);
                
            // make a tgrapherror

            emittance_tmp = sqrt (
                pow (h2_tmp->ProjectionX()->GetStdDev(),2) 
                *pow (h2_tmp->ProjectionY()->GetStdDev(),2) 
                - h2_tmp->GetCovariance()
                );
  
            // just a check on the lenght of the two arrays
            if(debug) {
                std::cout<<"Histogram number: "<<v_tmp2D.size()<<std::endl;
            }
        }
        h2.push_back(v_tmp2D);
    }

    TFile * output_file = new TFile(path+folder+"_VDgraphs.root", "recreate");
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