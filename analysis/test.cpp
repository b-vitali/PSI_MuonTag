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

#include "TRandom2.h"

struct Hit{
    int fEvent;
    int fParticleID;
    double x,dx,y,dy,z,t,phi;
};

double x_sigma=2;
double z_sigma=2;

double x_sigma_proj;
double z_sigma_proj;

TRandom2 * r = new TRandom2();
void smear (vector<double>  * v_fx, vector<double>  *v_fz, vector<double>  *v_fphi){
    for(int i=0; i< v_fx->size(); i++){
        if(i == 0)cout<<v_fx->at(i)<<"goes to ";
        // cout<<v_fphi->at(i)<<"goes to "<<cos(v_fphi->at(i)*TMath::Pi()/180)<<endl;

        // x_sigma_proj = abs(cos(v_fphi->at(i))*x_sigma);
        // z_sigma_proj = abs(sin(v_fphi->at(i))*z_sigma);

        v_fx->at(i) = r->Gaus(v_fx->at(i),x_sigma);
        v_fz->at(i) = r->Gaus(v_fz->at(i),z_sigma);
        if(i == 0)cout<<v_fx->at(i)<<endl;
    }
}

void asd (bool draw)
{
    TFile *_file0 = TFile::Open("build/data.root");
    TTree *T = _file0->Get<TTree>("VD");

    int fEvent;
    int fParticleID;
    double x,y,z,t,phi;
    
    vector<int> v_fEvent;
    vector<int> v_fParticleID;
    vector<double>  v_fx,v_fy,v_fz,v_ft, v_fphi;

    vector<Hit> Hits;

    T->SetBranchAddress("fEvent",&fEvent);
    T->SetBranchAddress("fParticleID",&fParticleID);
    T->SetBranchAddress("fPosX",&x);
    T->SetBranchAddress("fPosY",&y);
    T->SetBranchAddress("fPosZ",&z);
    T->SetBranchAddress("fVDTime",&t);
    T->SetBranchAddress("fRotY",&phi);


    TGraphErrors * grxy = new TGraphErrors;
    grxy->SetMarkerStyle(20);
    grxy->SetMarkerColor(kBlue);
    grxy->SetName("x");
    grxy->SetTitle("z");

    TH1F *h_d = new TH1F("h_d","h_d",100,0,30);   
    
    auto chi2Function = [&](const double *par) {
        //minimisation function computing the sum of squares of residuals
        // looping at the graph points
        int np = grxy->GetN();
        double f = 0;
        double *x = grxy->GetX();
        double *y = grxy->GetY();
        for (int i=0;i<np;i++) {
            double u = x[i] - par[0];
            double v = y[i] - par[1];
            double dr = par[2] - std::sqrt(u*u+v*v);
            f += dr*dr;
        } ///(x_sigma*x_sigma) /(z_sigma*z_sigma)
        return f;
    };

    // wrap chi2 function in a function object for the fit
    // 3 is the number of fit parameters (size of array par)
    ROOT::Math::Functor fcn(chi2Function,3);
    ROOT::Fit::Fitter  fitter;

    Hit tmp_hit;
    for(int i = 0; i< T->GetEntries(); i++){
        T->GetEntry(i);
        
        tmp_hit.fEvent=fEvent;
        tmp_hit.fParticleID=fParticleID;
        tmp_hit.x=x;
        tmp_hit.y=y;
        tmp_hit.z=z;
        tmp_hit.t=t;
        tmp_hit.phi=phi;
        Hits.push_back(tmp_hit);

        v_fEvent.push_back(fEvent);
        v_fParticleID.push_back(fParticleID);
        v_fx.push_back(x);
        v_fy.push_back(y);
        v_fz.push_back(z);
        v_ft.push_back(t);
        v_fphi.push_back(phi);
    }

    smear(&v_fx,&v_fz, &v_fphi);
cout << v_fx[0];
    int prev = 0;
        
    for(int i = 0; i< T->GetEntries(); i++){
        cout<<Hits.at(i).fEvent<<endl;
        // cout<< i <<" "<<v_fEvent[i]<<endl;
        if(v_fEvent[i]==prev){
            if(v_ft[i]<1 && v_fParticleID[i] == -11){
                cout<< i <<" "<<v_fEvent[i]<<endl;
       
                grxy->SetPoint(grxy->GetN(),v_fx[i],v_fz[i]);

                // x_sigma_proj = abs(cos(v_fphi[i]*x_sigma));
                // z_sigma_proj = abs(sin(v_fphi[i]*z_sigma));
                grxy->SetPointError(grxy->GetN()-1,x_sigma,z_sigma);
            }
        }
        else{ 
            prev=v_fEvent[i];           
            double pStart[3] = {0,0,1};
            fitter.SetFCN(fcn, pStart);
            fitter.Config().ParSettings(0).SetName("x0");
            fitter.Config().ParSettings(1).SetName("y0");
            fitter.Config().ParSettings(2).SetName("R");

            // do the fit
            bool ok = fitter.FitFCN();
            if (!ok) {
               Error("line3Dfit","Line3D Fit failed");
            }


            const ROOT::Fit::FitResult & result = fitter.Result();
            result.Print(std::cout);

            if(draw) grxy->Draw("AP");

            //Draw the muon circle on top
            TArc *mu_arc = new TArc(0,0,30);
            mu_arc->SetLineColor(kBlue);
            mu_arc->SetLineWidth(3);
            mu_arc->SetFillStyle(0);
            if(draw) mu_arc->Draw();

            //Draw the circle on top of the points
            TArc *arc = new TArc(result.Parameter(0),result.Parameter(1),result.Parameter(2));
            arc->SetLineColor(kRed);
            arc->SetLineWidth(4);
            arc->SetFillStyle(0);
            if(draw) arc->Draw();

            if(draw) {gPad->Modified(); gPad->Update();}
            if(draw) cin.get();

            double d = sqrt(result.Parameter(0)*result.Parameter(0)+result.Parameter(1)*result.Parameter(1))-result.Parameter(2)-30;
            cout<<"Distance "<<d<<endl;
            h_d->Fill(abs(d));

            grxy->Set(0);

            if(v_ft[i]<1 && v_fParticleID[i] == -11){
            cout<< i <<" "<<v_fEvent[i]<<endl;
            //fill first of new one
            grxy->SetPoint(grxy->GetN(),v_fx[i],v_fz[i]);
            grxy->SetPointError(grxy->GetN()-1,x_sigma,z_sigma);
            }
        }
    }
    h_d->Draw();
}