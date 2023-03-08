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

#include "TNtuple.h"

//? Name of the file
TString filename="data";
int NEvents=1e5;
bool print = false;
double PDE = 0.4;
double thr = 0;
double thr_tot = -3;
double tmax = 30;
double tmin = -10;
double DarkRate = 100e3;

/*
    This class is to store per each event all the gammas which arrived to SiPMs
    It contains ev, sipmno, position (x,y,z) and time
*/
class Data {
    public:           
        Data() {}     // Empty Constructor
        Data(int fev, int fsipm, double fx, double fy, double fz, double ft) { // Constructor
            ev = fev;
            sipm = fsipm;
            x = fx;
            y = fy;
            z = fz;
            t = ft;
        }     
        int ev, sipm;
        double x,y,z,t;
};

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
    if(print) cout<<prob<<endl;

    float whole, fractional;
    fractional = std::modf(prob, &whole);

    if(print) cout<<whole<<" "<<fractional<<endl;

    r.Uniform(tmin,prob);
    for(int i = 0; i < whole; i++){
        t = r.Uniform(tmin,tmax);
        v.push_back(new TF1("f",TString::Format("OneWave(x-%f)", t), tmin, tmax));
    }
    if(r.Rndm() < fractional) v.push_back(new TF1("f",TString::Format("OneWave(x-%f)", t), tmin, tmax));
    return v;
}


bool compareData(Data* d1, Data* d2)
{
    return (d1->sipm < d2->sipm);
}

void Waves(vector<TF1* > * waves, TNtuple* T_new, std::vector<Data *> data, double thr){
    
    double Tx,Ty,Tz,Tt,Tcharge,Tamplitude;
    int Tev, Tsipm;

    //? Sort the pair according to increasing SiPMNo
    sort(data.begin(), data.end(), compareData);

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
        if(print) cout<<data[i]->sipm<< endl;

        //? If it is a new channel conclude the previous, clear v and increase index
        if(channel != data[i]->sipm) {
            if(print) cout<< "new channel" << endl;
            channel = data[i]->sipm;
            if(index > -1) {
                v = AddDark(v);
                f_tmp = new TF1(TString::Format("f_ch%d", data[i]->sipm),SumTF1(v),tmin,tmax,0);
                Tamplitude = f_tmp->GetMinimum();
                if(Tamplitude < thr){
                    v_fsum.push_back(f_tmp);
                    // new TCanvas;
                    // v_fsum[index]->Draw();
                    Tev = data[i]->ev;
                    Tsipm = data[i]->sipm;
                    Tx = data[i]->x;
                    Ty = data[i]->y;
                    Tz = data[i]->z;
                    Tt = data[i]->t;
                    Tcharge = f_tmp->Integral(tmin,tmax);
                    T_new->Fill(Tev, Tsipm, Tx, Ty, Tz, Tt, Tcharge, Tamplitude);
                }
            }
            v.clear();
            index += 1;
        }

        //? PDE SiPM
        p = r.Rndm();
        if(p < PDE){
            if(print) cout<<setprecision(3)<<p<<"<"<<PDE<<" => New photoelectron"<<endl;
            v.push_back(new TF1("f",TString::Format("OneWave(x-%f)", data[i]->t), tmin, tmax));
        }
        else{if(print) cout<<setprecision(3)<<p<<">"<<PDE<<" => No photoelectron"<<endl;}
    }
    
    *waves = v_fsum;
}

TBrowser *OpenBrowser() { return new TBrowser; }

int CreateWaveforms(TString Tree, TFile * _file1){
    //? Open the file and import the TTree
    TFile *_file0 = TFile::Open("../build/"+filename+".root");
    TTree *T = _file0->Get<TTree>(Tree);
    cout<<Tree<<endl;

    //? Variables needed for the output
    vector<TF1*> * waves = new vector<TF1*>;
    double charge;

    TString treename= Tree + "_hits";
    TNtuple *T_new = new TNtuple(treename,"","ev:sipm:x:y:z:t:charge:amplitude");
    
    //? Variables to read the necessary branches
    std::vector<int> *Event   =0;
    std::vector<int> *SiPMNo  =0;
    std::vector<double> *X =0;
    std::vector<double> *Y =0;
    std::vector<double> *Z =0;
    std::vector<double> *Time =0;

    T->SetBranchAddress("fEvent",&(Event));
    T->SetBranchAddress("fSiPMNo",&(SiPMNo));
    T->SetBranchAddress("fPosSiPMInX",&(X));
    T->SetBranchAddress("fPosSiPMInY",&(Y));
    T->SetBranchAddress("fPosSiPMInZ",&(Z));
    T->SetBranchAddress("fTimeIn",&(Time));

    std::vector<Data *> data; 
    Data* data_tmp; 

    int ev;

    //? Loop on the entries (each contains vectors of time and sipmNo)
    if(NEvents > T->GetEntries()) NEvents = T->GetEntries();
    cout<< NEvents <<endl;
    for(int i = 0; i < NEvents; i++){ //i< T->GetEntries()

        cout<<"*****************************************"<<endl;
        cout<<"*************** Event N : " <<i<<"**************"<<endl;
        cout<<"*****************************************"<<endl;
        //? Get the i-entry of the TTree
        T->GetEntry(i);

        cout<<Event[0][0]<<endl;
        cout<<(*Event)[0]<<endl;
        cout<<SiPMNo[0][0]<<endl;
        cout<<SiPMNo[0].size()<<endl;
        cout<<"*****************************************"<<endl<<endl;

        data_tmp = new Data();
        
        //? Put info (for the i event) in array of Data
        data.reserve(SiPMNo[0].size());
        for(int j=0; j<SiPMNo[0].size(); j++){
            data.push_back( new Data(Event[0][j], SiPMNo[0][j], X[0][j], Y[0][j], Z[0][j], Time[0][j]));
        }

        //? Call the Waves routine which gives you an array of TF1
        Waves(waves, T_new, data, thr);
        if(print) cout<<"waves.size() ="<< waves->size()<<endl;

        //? Insert the TF1 in folders accordin to the Ev number
        ev = Event[0][i];
        TString dirname = TString::Format("Ev_%d_", ev)+Tree;
        _file1->mkdir(dirname);
        _file1->cd(dirname);
        for(int j = 0; j < waves->size(); j++){
            waves->at(j)->Write("");
        }
    }
    _file1->cd("");
    T_new->Write("");
    return 1;
}

int SiPM_Waveform(){
    cout<<endl;
    cout<<"*****************************************"<<endl;
    cout<<"*********** Waveform analysis ***********"<<endl;
    cout<<"*****************************************"<<endl;
    cout<<"* File name\t: "<<   filename<<"\t\t\t*"<<endl;
    cout<<"* N Events\t: "<<   NEvents<<"\t\t\t*"<<endl;
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

    //? Create the file in which we will store the waveforms
    TString  newfile = "wavetest_"+filename+".root";
	TFile *_file1 = TFile::Open(newfile,"recreate");
    if(choice) CreateWaveforms("SiPM_gate", _file1);
    // if(choice) CreateWaveforms("SiPM_in", _file1);

    OpenBrowser();
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