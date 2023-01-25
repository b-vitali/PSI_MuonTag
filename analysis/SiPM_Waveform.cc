#include <iostream>
#include <cmath>
#include "TCanvas.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TMath.h"
#include <stdlib.h>
#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TChain.h"
#include "TStyle.h"
#include "THStack.h"
#include "TH1.h"
#include "TH2.h"
#include "TString.h"


//? Name of the file
TString filename="test";
bool print = false;
double PDE = 0.4;
double thr = -3;
double tmax = 30;
double tmin = -10;
double DarkRate = 100e3;


//? Ask Giovanni
double Gp[5] = {8.88913032e-01,  6.46867235e+01,  4.16687779e-01, 212.49027149107357, 1.5};

double C(double x,double a,double b){
	return b * sqrt(TMath::Pi()/2) * exp(b*b/2/a/a - x/a)*(TMath::Erf((a*x - b*b)/(a*b*sqrt(2))) + 1);
}

double OneWave(double x){
	return Gp[3]*(C(x - Gp[4], Gp[0], Gp[2]) - C(x - Gp[4], Gp[0]*Gp[1]/(Gp[0]+Gp[1]), Gp[2]));
}

//? Sum a vector of functions std::vector<TF1*> fFuncList; 
struct SumTF1 { 
    SumTF1(const std::vector<TF1 *> & flist) : fFuncList(flist) {}
   
    double operator() (const double * x, const double *p) {
        double result = 0;
        for (unsigned int i = 0; i < fFuncList.size(); ++i) 
           result += -fFuncList[i]->EvalPar(x,p); 
        return result; 
    } 
   
    std::vector<TF1*> fFuncList; 
};

std::vector<TF1 *> AddDark(std::vector<TF1 *> v){
    TRandom3 r;
    r.SetSeed();
    double t;
    double prob = DarkRate*1e-6/tmax;
    // cout<<prob<<endl;

    float whole, fractional;
    fractional = std::modf(prob, &whole);

    // cout<<whole<<" "<<fractional<<endl;

    r.Uniform(tmin,prob);
    for(int i = 0; i < whole; i++){
        t = r.Uniform(tmin,tmax);
        v.push_back(new TF1("f",TString::Format("OneWave(x-%f)", t), tmin, tmax));
    }
    if(r.Rndm() < fractional) v.push_back(new TF1("f",TString::Format("OneWave(x-%f)", t), tmin, tmax));
    return v;
}

vector<TF1* > Waves(std::vector<std::pair<int, double>> data, double thr){
    //? Sort the pair according to increasing SiPMNo
    sort(data.begin(), data.end());

    //? Vector for the resulting functions
    vector<TF1 *> v_fsum;
    std::vector<TF1 *> v; 

    TF1* f_tmp;

    //? to implement the PDE
    TRandom3 r;
    r.SetSeed();
    double p;

    int channel = -1;
    int index = -1;
    for(int i = 0; i < data.size(); i++){
        if(print) cout<<std::get<0>(data[i])<< endl;

        //? If it is a new channel conclude the previous, clear v and increase index
        if(channel != std::get<0>(data[i])) {
            if(print) cout<< "new channel" << endl;
            channel = std::get<0>(data[i]);
            if(index > -1) {
                v = AddDark(v);
                f_tmp = new TF1(TString::Format("f_ch%d", std::get<0>(data[i])),SumTF1(v),tmin,tmax,0);
                if(f_tmp->GetMinimum() < thr){
                    v_fsum.push_back(f_tmp);
                    // new TCanvas;
                    // v_fsum[index]->Draw();
                }
            }
            v.clear();
            index += 1;
        }

        //? PDE SiPM
        p = r.Rndm();
        if(p < PDE){
            if(print) cout<<setprecision(3)<<p<<"<"<<PDE<<" => New photoelectron"<<endl;
            v.push_back(new TF1("f",TString::Format("OneWave(x-%f)", std::get<1>(data[i])), tmin, tmax));
        }
        else{if(print) cout<<setprecision(3)<<p<<">"<<PDE<<" => No photoelectron"<<endl;}

    }

    return v_fsum;
}

TBrowser *OpenBrowser() { return new TBrowser; }

int CreateWaveforms(){
    //? Open the file and import the TTree
    TFile *_file0 = TFile::Open("../build/"+filename+".root");
    TTree *T = _file0->Get<TTree>("SiPM_out");
    
    vector<TF1*> waves ;
    vector<TH1F *> h_charges;
    double charge;

    //? Variables to read the necessary branches
    std::vector<int> *Event=0;
    std::vector<int> *SiPMNo=0;
    std::vector<double> *Time=0;

    T->SetBranchAddress("fEvent",&(Event));
    T->SetBranchAddress("fSiPMNo",&(SiPMNo));
    T->SetBranchAddress("fTimeIn",&(Time));

    //? Create the file in which we will store the waveforms
	TFile *_file1 = TFile::Open("wavetest.root","recreate");
    int ev;

    //? Pairs with time and SiPMNo to generate the signals
    std::vector<std::pair<int, double>> data;

    //? Loop on the entries (each contains vectors of time and sipmNo)
    for(int i = 0; i< T->GetEntries(); i++){

        cout<<"*****************************************"<<endl;
        cout<<"*************** Event N : " <<i<<"**************"<<endl;
        cout<<"*****************************************"<<endl;
        //? Get the i-entry of the TTree
        T->GetEntry(i);

        cout<<Event[0][0]<<endl;
        cout<<SiPMNo[0][0]<<endl;
        cout<<SiPMNo[0].size()<<endl;
        cout<<"*****************************************"<<endl<<endl;

        //? Put SiPMNo and Time in a pair
        data.reserve(SiPMNo[0].size());
        std::transform(SiPMNo[0].begin(), SiPMNo[0].end(), Time[0].begin(), std::back_inserter(data),
               [](int a, double b) { return std::make_pair(a, b); });

        //? Call the Waves routine which gives you an array of TF1
        waves = Waves(data, thr);

        //? Insert the TF1 in folders accordin to the Ev number
        ev = Event[0][i];
        _file1->mkdir(TString::Format("Ev%d", ev));
        _file1->cd(TString::Format("Ev%d", ev));
        h_charges.push_back(new TH1F(TString::Format("h%d", ev),TString::Format("h%d", ev),100,-300,0));
        for(int j = 0; j < waves.size(); j++){
            charge = waves[j]->Integral(tmin,tmax);
            h_charges[i]->Fill(charge);
            waves[j]->Write("");
        }
        h_charges[i]->Write("");
    }
    OpenBrowser();
    return 1;
}

int SiPM_Waveform(){
    cout<<endl;
    cout<<"*****************************************"<<endl;
    cout<<"*********** Waveform analysis ***********"<<endl;
    cout<<"*****************************************"<<endl;
    cout<<"* File name\t: "<<   filename<<"\t\t\t*"<<endl;
    cout<<"* Debug\t\t: "<<     print<<"\t\t\t*"<<endl;
    cout<<"* PDE\t\t: "<<       PDE<<"\t\t\t*"<<endl;
    cout<<"* Threshold\t: "<<   thr<<"\t\t\t*"<<endl;
    cout<<"* Time window\t: ("<<tmin<<", "<<tmax<<")\t\t*"<<endl;
    std::cout.precision(2);
    cout<<"* Dark rate\t: "<<std::scientific<<   DarkRate<<"\t\t*"<<endl;
    cout<<"*****************************************"<<endl;
    cout<<"Proceede (1/0)?  ";
    
    bool choice;
    cin>>choice;

    if(choice) CreateWaveforms();
    return 1;
}

void binary(){
    TRandom3 r;
    r.SetSeed();
    short int a = r.Uniform(0,8);   // N bit which SciFi
    short int b = r.Uniform(0,2);   // 1 bit Up/Down Left/Right
    short int c = r.Uniform(0,128); // 8 bit which SiPM

    cout<<a<<" "<<b<<" "<<c<<" "<<endl<<endl;

    int bUD=1;
    int bSiPM=8;

    short int q = 0;
    cout<<"q = "<<q<<endl;
    q = a;
    cout<<"q = "<<q<<endl;
    q = q<<bUD;
    cout<<"q = "<<q<<endl;
    q = q+b;
    cout<<"q = "<<q<<endl;
    q = q<<bSiPM;
    cout<<"q = "<<q<<endl;
    q = q+c;
    cout<<"q = "<<q<<endl;


    short int pUD = 0b1;            // 1 bit Up/Down Left/Right
    short int pSiPM = 0b11111111;   // 8 bit which SiPM

    c = q&pSiPM;
    b = q>>bSiPM&pUD;
    a = (q>>bSiPM)>>bUD;

    cout<<endl<<a<<" "<<b<<" "<<c<<" "<<endl;

}