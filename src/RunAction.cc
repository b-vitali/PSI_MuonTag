/// \file  RunAction.cc
/// \brief Implementation of the RunAction class, takes care of ROOT output

#include "RunAction.hh"

RunAction::RunAction():fName("data.root")
{
	//DefineCommands();
	
	G4AnalysisManager *man = G4AnalysisManager::Instance();
	
	// Ntuple for the VirtualDetector
	man->CreateNtuple("VD", "VD");
	man->CreateNtupleIColumn("fEvent");
	man->CreateNtupleIColumn("fVDNo");
	man->CreateNtupleIColumn("fInOut");
	man->CreateNtupleIColumn("fParticleID");
	man->CreateNtupleDColumn("fVDTime");
	man->CreateNtupleDColumn("fMom");
	man->CreateNtupleDColumn("fPosX");
	man->CreateNtupleDColumn("fPosY");
	man->CreateNtupleDColumn("fPosZ");
	man->CreateNtupleDColumn("fMomX");
	man->CreateNtupleDColumn("fMomY");
	man->CreateNtupleDColumn("fMomZ");
	man->FinishNtuple(0);

	// Ntuple for the PrimaryGenerator
	man->CreateNtuple("gun", "gun");
	man->CreateNtupleIColumn("fEvent");
	man->CreateNtupleDColumn("fMom");
	man->CreateNtupleDColumn("fPosX");
	man->CreateNtupleDColumn("fPosY");
	man->CreateNtupleDColumn("fPosZ");
	man->CreateNtupleDColumn("fAngX");
	man->CreateNtupleDColumn("fAngY");
	man->CreateNtupleDColumn("fMomX");
	man->CreateNtupleDColumn("fMomY");
	man->CreateNtupleDColumn("fMomZ");
	man->FinishNtuple(1);

	fMessenger = new RunActionMessenger(this);
}

RunAction::~RunAction()
{
	delete fMessenger;	
}

void RunAction::BeginOfRunAction(const G4Run* run)
{
	G4AnalysisManager *man = G4AnalysisManager::Instance();

	G4int runID = run->GetRunID();
	std::stringstream strRunID;
	strRunID << runID;
	man->OpenFile(fName);

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