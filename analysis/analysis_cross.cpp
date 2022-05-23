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

void test(){    
    TFile * f = TFile::Open("data/28MeV_10k2/28MeV_0.10.root");

    TTree * T_VD = f->Get<TTree>("VD");
    TTree * T_telescope = f->Get<TTree>("Scint_telescope");

    int fEvent_telescope, fEvent_VD, fNgamma_telescope, fVDNo, fParticleID;
    double fMomZ_VD;
    T_telescope->SetBranchAddress("fEvent",&fEvent_telescope);
    T_telescope->SetBranchAddress("fNgamma",&fNgamma_telescope);
    T_VD->SetBranchAddress("fEvent",&fEvent_VD);
    T_VD->SetBranchAddress("fMomZ",&fMomZ_VD);
    T_VD->SetBranchAddress("fVDNo",&fVDNo);
    T_VD->SetBranchAddress("fParticleID",&fParticleID);

    TH2F * h = new TH2F("cross","cross",10000,0,30000,300,0,30);

    //to keep track of where we are in the ttree;
    int j=0;
    int k=0; 
    int N_telescope = T_telescope->GetEntries();
    int N_VD = T_VD->GetEntries();
            
    cout<<N_telescope<<" "<<N_VD<<endl;

    for(j=0; j<N_telescope;j++){
        T_telescope->GetEntry(j);
        for(int k =0; k<N_VD; k++){
            T_VD->GetEntry(k);
            if(fEvent_telescope==fEvent_VD){
                if(fParticleID == -13 && fVDNo == 1) h->Fill(fNgamma_telescope,fMomZ_VD);
            } 
            // if(k%1000==0) cout<<k<<" "<<j<<endl;
        }
        if(j%1000==0)cout<<j<<endl;
    }

    h->Draw("colz");
}