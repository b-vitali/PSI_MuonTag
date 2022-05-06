/// \file  RunAction.cc
/// \brief Implementation of the RunAction class, takes care of ROOT output

#include "RunAction.hh"

RunAction::RunAction()
{
	//DefineCommands();
	
	G4AnalysisManager *man = G4AnalysisManager::Instance();
	
	// Ntuple for the VirtualDetector
	man->CreateNtuple("VD", "VD");
	man->CreateNtupleIColumn("fEvent");
	man->CreateNtupleIColumn("fVDNo");
	man->CreateNtupleIColumn("fParticleID");
	man->CreateNtupleDColumn("fVDTime");
	man->CreateNtupleDColumn("fMom");
	man->FinishNtuple(0);

}

RunAction::~RunAction()
{}

void RunAction::BeginOfRunAction(const G4Run* run)
{
	G4AnalysisManager *man = G4AnalysisManager::Instance();

	G4int runID = run->GetRunID();
	std::stringstream strRunID;
	strRunID << runID;
	man->OpenFile("output_"+strRunID.str()+".root");

	// man->CreateH2("Scint_PhotonTime","PhotonTime", 200, 0., 200, 5, 0,5);
	// man->CreateH3("Scint_PhotonPos_Front","PhotonPositionFront",200, -50., 50, 200, -50., 50, 5, 0,5);
	// man->CreateH3("Scint_PhotonPos_Back","PhotonPositionBack",200, -50., 50, 200, -50., 50, 5, 0,5);
	// man->CreateH3("Scint_PhotonPos_Left","PhotonPositionLeft",200, -5., 5, 200, -50., 50, 5, 0,5);
	// man->CreateH3("Scint_PhotonPos_Right","PhotonPositionRight",200, -5., 5, 200, -50., 50, 5, 0,5);
	// man->CreateH3("Scint_PhotonPos_Up","PhotonPositionUp",200, -50., 50, 200, -50., 50, 5, 0,5);
	// man->CreateH3("Scint_PhotonPos_Down","PhotonPositionDown",200, -50., 50, 200, -50., 50, 5, 0,5);


	/*
	fData = TFile::Open(fName, "RECREATE");
	fTree = new TTree("T","A tree containing simulation values");
	
	// Defining tree branches
	fTree->Branch( "VDTime",  &fVDTime);
	*/
}


void RunAction::EndOfRunAction(const G4Run*)
{
	// Write and Close the Root file
	G4AnalysisManager *man = G4AnalysisManager::Instance();
	man->Write();
	man->CloseFile();

	/*
	fData->cd();
	//fTree->Print();
	fTree->Write();
	fData->Close();
	*/
}