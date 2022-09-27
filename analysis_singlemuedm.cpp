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

TString path="/home/bastiano/Documents/Geant4/PSI/FromScratch/";

TString folder="100_200x10x3"; //28MeV_10k2

std::vector<TString> file_names;

std::vector<double_t> thickness{
        0.0
        ,5.0
        ,10.0
        ,15.0
        ,20.0
        ,25.0
        ,30.0
        ,35.0
        ,40.0
        ,45.0
        ,50.0
        ,55.0
        ,60.0
        ,65.0
        ,70.0
        ,75.0
        ,80.0
        ,85.0
        ,90.0
        ,95.0
};

// these are the commands you would give to TTree->Draw() with the branch names
// Choose *variable* *some cut* *histogram range*
std::vector< std::tuple<char*, char*, char*> > plots {
        {   (char*)"fSiPMNo",          (char*)"",  (char*)"(10,0,10)"     }
       ,{   (char*)"fPosSiPMInX",          (char*)"",  (char*)"(300,-150,150)"     }
       ,{   (char*)"fTimeIn",          (char*)"fSiPMNo>3",  (char*)"(300,-150,150)"     }
       ,{   (char*)"fTimeIn",          (char*)"fSiPMNo<=3",  (char*)"(300,-150,150)"     }

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

void epos(){
    //? Compile a list of the files you wnat to use
    cout<<"Data from the files in path : "<<path<<endl<<endl;
    
    std::stringstream ss;
    for(int i = 0; i< thickness.size(); i++){
        ss << std::fixed<<std::setprecision(1)<<thickness[i];
        TString str = ss.str();

        file_names.push_back(str+".root");
        std::cout<<file_names[i]<<std::endl;
        
        ss.str(std::string());
    }
    
    //? Create a vector to keep all the ttrees you will extract
    std::vector<TTree *> tree_v;

    //? Create the TGraph to be filled and the variables to calculate what you need
    TGraphErrors * grx = new TGraphErrors;
    grx->SetMarkerStyle(20);
    grx->SetMarkerColor(kBlue);
    grx->SetName("epos_{x}");
    grx->SetTitle("Positron x position");

    TGraphErrors * gry = new TGraphErrors;
    gry->SetMarkerStyle(20);
    gry->SetMarkerColor(kRed);
    gry->SetName("epos_{y}");
    gry->SetTitle("Positron y position");

    double_t mean_x=0;
    double_t mean_x_ev=0;
    double_t mean_y=0;
    double_t mean_sipmno=0;
    double_t mean_sipmno_ev=0;

    double_t nhits=0;
    double_t nhits_ev=0;

    TH2F * h2 = new TH2F("h2", "h2", 20,0,100,200,-100,100);

    //? Reading vector from ttree is not straightforward
    //? std::vector <double_t> * x_tmp; doesnt work. you need to put = 0; 
    std::vector <double_t> * x_tmp=0;
    std::vector <double_t> * y_tmp=0;
    std::vector <double_t> * z_tmp=0;
    std::vector <int> * fSiPMNo_tmp=0;

    std::vector < std::vector <double_t> > x;
    std::vector < std::vector <double_t> > y;
    std::vector < std::vector <double_t> > z;
    std::vector < std::vector <int> > fSiPMNo;

    //? open the files and collect the TTrees in a vector
    TFile *file_tmp;
    for(const TString &file_name : file_names){
        file_tmp =  TFile::Open(path+folder+"/"+file_name);
        tree_v.push_back(file_tmp->Get<TTree>("SiPM"));
    }

    //? Loop over the vector of ttrees:
    /*
    What eactually happens?
        - For each file read the vectors of each event
        - Compile them in vectors of vectors
        - Evaluate whatever you need
        - Feed this variable to the TGraphs
    */
    for (int i=0; i<tree_v.size();i++){
        cout<< "tree "<<i<<endl;
        tree_v[i]->SetBranchAddress("fPosInX",&x_tmp);
        tree_v[i]->SetBranchAddress("fPosInY",&y_tmp);
        //tree_v[i]->SetBranchAddress("fPosInZ",&z_tmp);
        tree_v[i]->SetBranchAddress("fSiPMNo",&fSiPMNo_tmp);
        
        for (int j=0; j<tree_v[i]->GetEntries(); j++){
            tree_v[i]->GetEntry(j);

            //cout<<"xtmp size "<<x_tmp->size()<<endl;

            x.push_back(*x_tmp);
            y.push_back(*y_tmp);
            // z.push_back(z_tmp);
            fSiPMNo.push_back(*fSiPMNo_tmp);
        }

        //cout<<"x size "<<x.size()<<endl;

        //? For every file, if the vector is not empty
        if(x.size() != 0) {

            //? Loop on the events
            for(int j=0; j < x.size(); j++){
                //cout<<"x ["<<i<<"] size"<<x[j].size()<<endl;
                //cout<<"x size "<<x.size()<<endl;
                nhits +=x[j].size();

                //? Loop on the arrays in every event
                for(int k=0; k < x[j].size(); k++){
                    nhits_ev+=1;

                    mean_x+=x[j][k];
                    mean_x_ev+=x[j][k];
                    mean_y+=y[j][k];

                    mean_sipmno+=fSiPMNo[j][k];
                    mean_sipmno_ev+=fSiPMNo[j][k];
                    //cout<<fSiPMNo[j][k]<<endl;
                    //cout<<mean_sipmno<<endl;

                //? Here you will fill graphs "PER EVENT"
                // gry->SetPoint(gry->GetN(),thickness[i],fSiPMNo[i]);
                // gry->SetPointError(gry->GetN()-1,x[i][j],y[i][j]);
                // cout<<"k "<<k<<endl;
                }
                h2->Fill(thickness[i], mean_x_ev/nhits_ev);
                mean_sipmno_ev = 0; mean_x_ev=0; nhits_ev = 0;
            }
            mean_x=mean_x/nhits;
            mean_y=mean_y/nhits;
            mean_sipmno=mean_sipmno/nhits;

            //? Here you will fill graphs "PER FILE"
            grx->SetPoint(grx->GetN(),thickness[i],mean_x);
            grx->SetPointError(grx->GetN()-1,0,0);

            gry->SetPoint(gry->GetN(),thickness[i],mean_y);
            gry->SetPointError(grx->GetN()-1,0,0);
        }

        //? Before looping to the next file clear the variables
        x.clear(); y.clear(); z.clear(); fSiPMNo.clear();
        nhits = 0; mean_x=0; mean_y=0; mean_sipmno=0;
    }
    new TCanvas; h2->Draw();
    new TCanvas; grx->Draw();
    new TCanvas; gry->Draw();


    TFile * output_file = new TFile(folder+"/"+"epos-graphs.root", "recreate");
    grx->Write();
    gry->Write();
    h2->Write();
    output_file->Close();
}


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
        ss << std::fixed<<std::setprecision(1)<<thickness[i];
        TString str = ss.str();

        file_names.push_back(str+".root");
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
        tree_v.push_back(file_tmp->Get<TTree>("SiPM"));
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

    TFile * output_file = new TFile(folder+"/"+"test-graphs.root", "recreate");
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