
#include <string.h>
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"

int runs = 28;

double_t thickness[] = {
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

void sum_files(){
    int threads = 8;
    
    for(int i = 0; i< runs; i++){
        ; ;
        std::cout<<std::fixed<<std::setprecision(2)<<"hadd 28MeV_"<<thickness[i]<<".root "<<TString::Format("output_%i_t{0..%i}.root; ",i , threads-1);
    }
}
