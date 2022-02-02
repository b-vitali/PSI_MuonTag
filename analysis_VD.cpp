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

void asd ()
{
    TFile *_file0 = TFile::Open("build/output.root");
    TTree *T = _file0->Get<TTree>("VD");

    int Row;
    int fEvent;
    int fParticleID;
    double fEout;

    vector<int> v_fEvent;
    vector<int> v_fParticleID;
    vector<double> v_fEout;

    T->SetBranchAddress("Row",&Row);
    T->SetBranchAddress("fEvent",&fEvent);
    T->SetBranchAddress("fParticleID",&fParticleID);
    T->SetBranchAddress("fEout",&fEout);

    for(int i = 0; i< T->GetEntries(); i++){
        T->GetEntry(i);
        v_fEvent.push_back(fEvent);
        v_fParticleID.push_back(fParticleID);
        v_fEout.push_back(fEout);
    }

   TH2F *h2 = new TH2F("h2","",40,0,20,40,0,20);

    for(int i = 0; i < T->GetEntries(); i++){
        //std::cout<<v_fEvent[i]<<" "<<v_fParticleID[i]<<" "<<v_fEout[i]<<" "<<endl;
        if(v_fEvent[i]==v_fEvent[i+1] && v_fParticleID[i] == -v_fParticleID [i+1] && std::abs (v_fParticleID[i]) == 1)
        {
            std::cout<<v_fEvent[i]<<" "<<v_fParticleID[i]<<" "<<v_fEout[i]<<" "<<endl;
            std::cout<<v_fEvent[i+1]<<" "<<v_fParticleID[i+1]<<" "<<v_fEout[i+1]<<" "<<endl;
            if(v_fParticleID[i]==-1) h2->Fill(v_fEout[i],v_fEout[i+1]);
            else h2->Fill(v_fEout[i+1],v_fEout[i]);
        }
    }
    h2->Draw("colz");
    new TCanvas;
    h2->ProjectionX()->Draw();
    new TCanvas;
    h2->ProjectionY()->Draw();

}