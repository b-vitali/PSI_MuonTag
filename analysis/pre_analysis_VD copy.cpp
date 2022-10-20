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
    int Event;
    int ParticleID;
    double x,dx,y,dy,z,dz,t,dt,phi;
    double px,py,pz;
};

// Name of the file
TString filename="type2_5deg";

// Time cut on the hit to fit
double t_cut = 2;

double x_sigma=1;   //5; gaus 7 bin da 250 micron => 1 mm
double y_sigma=10;  //10; dipende da se fai una X o no
double z_sigma=1;   //1; spessore della fibra
double t_sigma=0.3;   // 500 ps come risoluzione in tempo?

double x_sigma_proj;
double y_sigma_proj;
double z_sigma_proj;
double t_sigma_proj;

// Function to smear my points with the experimental resolutions (projected on the axis)
TRandom2 * r = new TRandom2();
void smear (vector<Hit>  * Hits){
    Hit * hit;
    for(int i=0; i< Hits->size(); i++){
        hit = &Hits->at(i);
        x_sigma_proj = sqrt(pow(sin(hit->phi)*z_sigma,2) + pow(cos(hit->phi)*x_sigma,2));
        y_sigma_proj = y_sigma;
        z_sigma_proj = sqrt(pow(cos(hit->phi)*z_sigma,2) + pow(sin(hit->phi)*x_sigma,2));
        t_sigma_proj = t_sigma;
        hit->x = r->Gaus(hit->x,x_sigma_proj);
        hit->y = r->Gaus(hit->y,y_sigma_proj);
        hit->z = r->Gaus(hit->z,z_sigma_proj);
        hit->t = r->Gaus(hit->t,t_sigma_proj);
        hit->dx = x_sigma_proj;
        hit->dy = y_sigma_proj;
        hit->dz = z_sigma_proj;
        hit->dt = t_sigma_proj;
    }
}

void pre_analysis_VD(){
    cout<<"____________________________________________________________"<<endl;
    cout<<"Hi, this is a simple analysis of the MuEDM e+ reconstruction"<<endl;
    cout<<"The main function is asd(bool display):"<<endl;
    cout<<"\t display:"<<endl;
    cout<<"\t\t 0: shows only the overall plots"<<endl;
    cout<<"\t\t 1: shows the reconstruction event by event"<<endl;
}


void ytlinear(TGraphErrors * gryt){
    int to = 0;
    double chi2=1e6;
    double chi2_tmp = 0;
    TH1F * chi = new TH1F("chi","chi",1000,0,100);
    TF1 * fit = new TF1("fit","pol1", -10, 10);
    TFitResultPtr result;
    double * t = gryt->GetX();

    // double x_min = TMath::MinElement(gryt->GetN(),gryt->GetX());
    double x_min = TMath::MinElement(gryt->GetN(),gryt->GetX());
    cout<<x_min<<endl;


    for(int i = 0; i<gryt->GetN(); i++){
        fit->SetParameters(0,0);
        cout<<x_min<<" "<<t[i]<<endl;
        result = gryt->Fit(fit,"S","", x_min,t[i]);
        // gryt->Draw();
        // gPad->Modified(); gPad->Update();
        // cin.get();
        chi2_tmp = result->MinFcnValue()/result->Ndf();
        chi->Fill(chi2_tmp);
        cout<<chi2_tmp<<endl;
        if(result->Ndf()>3 && chi2_tmp<chi2) {cout<<"less"<<endl; chi2 = chi2_tmp; to = i;}
    }

    cout<<x_min<<" "<<t[to]<<endl;
    result = gryt->Fit(fit,"S","",x_min,t[to]);
    gPad->Modified(); gPad->Update();
    // new TCanvas();
    // gryt->Draw();
    // 
    // new TCanvas();chi->Draw();
    // gPad->Modified(); gPad->Update();
    // cin.get();


}

void asd (bool display)
{
    TFile *_file0 = TFile::Open("build/"+filename+".root");
    TTree *T = _file0->Get<TTree>("VD");

    int fEvent;
    int fParticleID;
    double x,y,z,t,phi;
    double px,py,pz;
    
    vector<Hit> Hits;
    vector<double> angles;

    T->SetBranchAddress("fEvent",&fEvent);
    T->SetBranchAddress("fParticleID",&fParticleID);
    T->SetBranchAddress("fPosX",&x);
    T->SetBranchAddress("fPosY",&y);
    T->SetBranchAddress("fPosZ",&z);
    T->SetBranchAddress("fMomX",&px);
    T->SetBranchAddress("fMomY",&py);
    T->SetBranchAddress("fMomZ",&pz);
    T->SetBranchAddress("fVDTime",&t);
    T->SetBranchAddress("fRotY",&phi);

    TCanvas *c_xz=new TCanvas("xz");
    c_xz->SetWindowSize(500,500);
    c_xz->SetWindowPosition(0,0);
    TGraphErrors * grxz = new TGraphErrors;
    grxz->SetMarkerStyle(20);
    grxz->SetMarkerColor(kBlue);
    grxz->SetName("XZ circle");
    grxz->SetTitle("XZ circle");


    TCanvas *c_yt=new TCanvas("yt");
    c_yt->SetWindowSize(500,500);
    c_yt->SetWindowPosition(0,570);
    TGraphErrors * gryt = new TGraphErrors;
    gryt->SetMarkerStyle(20);
    gryt->SetMarkerColor(kBlue);
    gryt->SetName("YT Line");
    gryt->SetTitle("YT Line");

    TH1F *h_d = new TH1F("h_d","h_d",100,0,30);   
    TH1F *h_Dp = new TH1F("h_Dp","h_Dp",40,-10,10);   
    TH2F *h_yt = new TH2F("h_yt","h_yt",100,0,10, 300,-150,150);   
    h_yt->SetOption("colz");

    auto chi2Function = [&](const double *par) {
        //minimisation function computing the sum of squares of residuals
        // looping at the graph points
        int np = grxz->GetN();
        double f = 0;
        double *x = grxz->GetX();
        double *y = grxz->GetY();
        double *dx = grxz->GetEX();
        double *dy = grxz->GetEY();
        for (int i=0;i<np;i++) {
            double u = x[i] - par[0];
            double v = y[i] - par[1];
            double sigma_dr = 2*(dx[i]+dy[i])/std::sqrt(u*u+v*v);
            double dr = par[2] - std::sqrt(u*u+v*v);
            f += dr*dr;///(sigma_dr*sigma_dr) ;
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
        
        tmp_hit.Event=fEvent;
        tmp_hit.ParticleID=fParticleID;
        tmp_hit.x=x;
        tmp_hit.y=y;
        tmp_hit.z=z;
        tmp_hit.t=t;
        tmp_hit.px=px;
        tmp_hit.py=py;
        tmp_hit.pz=pz;
        tmp_hit.phi=phi;
        Hits.push_back(tmp_hit);
        angles.push_back(acos(phi));
    }

    smear(&Hits);
    
    int j=0; //first hit of the stack
    int prev = 0;
    Hit hit;   
    for(int i = 0; i< T->GetEntries()+1; i++){ //300
        if( i == T->GetEntries() ) goto last;
        hit = Hits.at(i);
        if(hit.Event==prev){
            // cout<<"prev "<<hit.Event<<endl;
            if(hit.ParticleID == -11){
                gryt->SetPoint(gryt->GetN(),hit.t,hit.y);
                gryt->SetPointError(gryt->GetN()-1,hit.dt,hit.dy);
                h_yt->Fill(hit.t,hit.y);
            }
            if(hit.t<t_cut && hit.ParticleID == -11){       
                grxz->SetPoint(grxz->GetN(),hit.x,hit.z);
                grxz->SetPointError(grxz->GetN()-1,hit.dx,hit.dz);
            }
        }
        else if(hit.Event!=prev){ 
            last:
            // cout<<"post "<<hit.Event<<endl;
            prev=hit.Event;           
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

            if(display) {
                c_yt->cd();
                gryt->Draw("AP");
                ytlinear(gryt);
                gryt->GetYaxis()->SetRangeUser( -150, 150);

                TLine *line_tcut = new TLine(t_cut,-150,t_cut,150);
                line_tcut->Draw();     
                c_yt->Modified(); c_yt->Update();

                c_xz->cd();
                grxz->Draw("AP");            
                grxz->GetXaxis()->SetLimits( -100, 100);
                grxz->GetYaxis()->SetRangeUser( -100, 100);
                c_xz->Modified(); c_xz->Update();

                //Draw the muon circle on top
                TArc *mu_arc = new TArc(0,0,30);
                mu_arc->SetLineColor(kBlue);
                mu_arc->SetLineWidth(3);
                mu_arc->SetFillStyle(0);
                mu_arc->Draw();

                //Draw the circle on top of the points
                TArc *arc = new TArc(result.Parameter(0),result.Parameter(1),result.Parameter(2));
                arc->SetLineColor(kRed);
                arc->SetLineWidth(4);
                arc->SetFillStyle(0);
                arc->Draw();

                //Draw the circle on top of the points
                TLine *line = new TLine(0,0,result.Parameter(0),result.Parameter(1));
                line->SetLineColor(kRed);
                line->SetLineWidth(1);
                line->Draw();

                Double_t pl_x[4];
                Double_t pl_y[4];
                double Z = 5;
                double X = 30;
                for(int j=0; j<angles.size(); j++){
                    // Double_t pl_x[4]={30*cos(angles.at(j))-0.5*Z*sin(angles.at(j)),30*cos(angles.at(j))+0.5*Z*sin(angles.at(j)),90*cos(angles.at(j))+0.5*Z*sin(angles.at(j)),90*cos(angles.at(jBlue))-0.5*Z*sin(angles.at(j))};
                    // Double_t pl_y[4]={30*sin(angles.at(j))+0.5*Z*cos(angles.at(j)),30*sin(angles.at(j))-0.5*Z*cos(angles.at(j)),90*sin(angles.at(j))-0.5*Z*cos(angles.at(j)),90*sin(angles.at(j))+0.5*Z*cos(angles.at(j))};
                    // TPolyLine *pline = new TPolyLine(4,pl_x, pl_y);
                    // pline->SetFillColor(38);
                    // pline->SetLineColor(2);
                    //pline->Draw("l");
                }

                double p = sqrt(Hits.at(j).px*Hits.at(j).px+Hits.at(j).py*Hits.at(j).py+Hits.at(j).pz*Hits.at(j).pz);
                auto legend = new TLegend(0.1,0.7,0.4,0.9);
                //cout<<result.MinFcnValue()<<" "<<grxz->GetN()-result.NFreeParameters()<<endl;
                //legend->AddEntry("",TString::Format("Fit: %0.2f/%d",result.MinFcnValue(),grxz->GetN()-result.NFreeParameters()),"");
                legend->AddEntry("",TString::Format("Ev: %d",Hits.at(j).Event),"");
                legend->AddEntry(grxz,"Smeared Points","LP");
                double p_rec =0.3*3*result.Parameter(2);
                double dp_rec=0.3*3*result.Error(2);
                
                legend->AddEntry(arc,TString::Format("%0.1f#pm%0.1f/%0.1f MeV", p_rec, dp_rec, p),"l");
                cout<<result.Error(2)<<endl;
                legend->AddEntry(mu_arc,"muon orbit","l");
                legend->SetTextSize(0.03);
                legend->Draw();

                c_xz->Modified(); 
                c_xz->Update();

                j=i;
                cin.get();
            }

            double d = sqrt(result.Parameter(0)*result.Parameter(0)+result.Parameter(1)*result.Parameter(1))-result.Parameter(2)-30;
            cout<<"Distance "<<d<<endl;
            h_d->Fill(abs(d));

            // Clear the TGraphErrors
            grxz->Set(0);
            gryt->Set(0);

            double p = sqrt(Hits.at(j).px*Hits.at(j).px+Hits.at(j).py*Hits.at(j).py+Hits.at(j).pz*Hits.at(j).pz);
            h_Dp->Fill(0.3*3*result.Parameter(2)-p);

        //! is this ok?
            if( i != T->GetEntries() )i = i-1;
            //fill first of new one
            // if(hit.t<1 && hit.ParticleID == -11){
            //     grxz->SetPoint(grxz->GetN(),hit.x,hit.z);
            //     grxz->SetPointError(grxz->GetN()-1,x_sigma,z_sigma);
            // }
        }
    }
            // double d = sqrt(result.Parameter(0)*result.Parameter(0)+result.Parameter(1)*result.Parameter(1))-result.Parameter(2)-30;
            // cout<<"Distance "<<d<<endl;
            // h_d->Fill(abs(d));
            // double p = sqrt(Hits.at(j).px*Hits.at(j).px+Hits.at(j).py*Hits.at(j).py+Hits.at(j).pz*Hits.at(j).pz);
            // h_Dp->Fill(0.3*3*result.Parameter(2)-p);
    
        if(!display){
            new TCanvas;
            h_d->Draw();
            gPad->Modified(); 
            gPad->Update();

            new TCanvas;
            h_Dp->Draw();
            gPad->Modified(); 
            gPad->Update();

            new TCanvas;
            h_yt->Draw();
            gPad->Modified(); 
            gPad->Update();

            new TCanvas;
            h_yt->ProfileX()->Draw();
            gPad->Modified(); 
            gPad->Update();

            cin.get();

            TFile *out_file = TFile::Open("data/"+filename+".root","recreate");
            out_file->Add(h_d);
            out_file->Add(h_Dp);
            out_file->Add(h_yt);
            out_file->Add(h_yt->ProfileX());
            out_file->Write();
            out_file->Close();        
        }
}


double f_dx, f_dy;
double dx = 1, dy = 10;
bool both=false;
bool squared=false;
Double_t Xerrors(double phi) { 
        if(!both){
            double Dx=sqrt(pow(dx*cos(phi),2) + pow(dy*sin(phi),2));
            double Dy=sqrt(pow(dx*sin(phi),2) + pow(dy*cos(phi),2));
            f_dx = sqrt(Dx*Dx + dx*dx);
            if(!squared) f_dx = min(Dx, dx);
        }
        else{
            double Dx=sqrt(pow(dx*cos(phi),2) + pow(dy*sin(phi),2));
            double Dy=sqrt(pow(dx*sin(phi),2) + pow(dy*cos(phi),2));
            f_dx = sqrt(Dx*Dx + Dx*Dx); 
            if(!squared) f_dx = Dx;
        }       
        return f_dx; 
    }
Double_t Yerrors(double phi) { 
        if(!both){
            double Dx=sqrt(pow(dx*cos(phi),2) + pow(dy*sin(phi),2));
            double Dy=sqrt(pow(dx*sin(phi),2) + pow(dy*cos(phi),2));
            f_dy = sqrt(Dy*Dy + dy*dy);
            if(!squared) f_dy = min(Dy, dy);
        }
        else{
            double Dx=sqrt(pow(dx*cos(phi),2) + pow(dy*sin(phi),2));
            double Dy=sqrt(pow(dx*sin(phi),2) + pow(dy*cos(phi),2));
            f_dy = sqrt(Dy*Dy + Dy*Dy); 
            if(!squared) f_dy = Dy;
        }   
        return f_dy;  
    }

void Errors(bool a, bool b){
    both = a;
    squared = b;
    TF1 *f_xerrors = new TF1("f_xerrors","Xerrors(x)",0,3.14);
    TF1 *f_yerrors = new TF1("f_yerrors","Yerrors(x)",0,3.14);
    new TCanvas;
    f_xerrors->Draw();
    new TCanvas;
    f_yerrors->Draw();
}