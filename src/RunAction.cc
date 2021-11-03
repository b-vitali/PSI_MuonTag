/// \file  RunAction.cc
/// \brief Implementation of the RunAction class

#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunActionMessenger.hh"

#include "TFile.h"
#include "TTree.h"
#include "globals.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "G4GenericMessenger.hh"

RunAction::RunAction() : 
	G4UserRunAction(), fData(nullptr), fTree(nullptr), fCmdOCT(false), 
	fCmdDN(false), fCmdPhotons(1), fCmdTracks(1), fRight(0), fLeft(0), 
	fDown(0), fUp(0), fBack(0), fFront(0), fSiPM(0), fGunTime(0), fDNTime(0), 
	fGunTimeMean(1/(1.9e9*CLHEP::hertz)), fDNTimeMean(1/(90*CLHEP::kilohertz)), 
	fName("./data.root"){
	//DefineCommands();
	fMessenger = new RunActionMessenger(this);
}

RunAction::~RunAction(){
	delete fMessenger;
}

void RunAction::BeginOfRunAction(const G4Run*){
	fGunTime = 0;
	fDNTime = 0;
	this->AdvanceDNTime();

	fData = TFile::Open(fName, "RECREATE");
	fTree = new TTree("T","A tree containing simulation values");
	
	// Defining tree branches
	fTree->Branch( "ein",  &fEin);
	fTree->Branch("edep", &fEdep);
	fTree->Branch("eout", &fEout);
	fTree->Branch("delta", &fDelta);
	fTree->Branch("ThetaIn", &fThetaIn);
	fTree->Branch("TrackLength", &fTrackLength);
	fTree->Branch("thetapositron", &fThetaPositron);
	fTree->Branch("eventID", &fID);
	fTree->Branch("bounce",  &fBounce);
	fTree->Branch("trackX", &fPosX);
	fTree->Branch("trackY", &fPosY);
	fTree->Branch("trackZ", &fPosZ);
	fTree->Branch("trackT", &fTime);
	fTree->Branch("Ngamma", &fNgamma);
	fTree->Branch("NgammaSec", &fNgammaSec);
	if(fCmdPhotons == 0){
		fTree->Branch("CerTag", &fCer);
		fTree->Branch("costhetagamma", &fThetaGamma);
		fTree->Branch("timegamma", &fTimeGamma);
		fTree->Branch("egamma", &fEGamma);
	}
	
	fTree->Branch("CerNumber", &fNCer);
	fTree->Branch("currentright", &fRight);
	fTree->Branch("currentleft", &fLeft);
	fTree->Branch("currentdown", &fDown);
	fTree->Branch("currentup", &fUp);
	fTree->Branch("currentback", &fBack);
	fTree->Branch("currentfront", &fFront);
	fTree->Branch("SiPM", &fSiPM);

	fTree->Branch("NCells", &fNCells);
	fTree->Branch("NPhotoElectrons", &fNPhotoElectrons);
	fTree->Branch("Cells", &fCells);
	fTree->Branch("CellTime", &fCellTime);
	fTree->Branch("OCTflag", &fOCTflag);
	fTree->Branch("DNflag", &fDNflag);
	fTree->Branch("GunTime", &fGunTime);
	fTree->Branch("DecayTime", &fDecayTime);
}


void RunAction::EndOfRunAction(const G4Run*){
	fData->cd();
	//fTree->Print();
	fTree->Write();
	fData->Close();
}


