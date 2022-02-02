/// \file  RunAction.cc
/// \brief Implementation of the RunAction class, takes care of ROOT output

#include "RunAction.hh"

RunAction::RunAction()
{
	//DefineCommands();
}

RunAction::~RunAction()
{}

void RunAction::BeginOfRunAction(const G4Run*)
{
	G4AnalysisManager *man = G4AnalysisManager::Instance();
	man->OpenFile("output.root");

	// Ntuple for the VirtualDetector
	man->CreateNtuple("VD", "VD");
	man->CreateNtupleIColumn("fEvent");
	man->CreateNtupleDColumn("fVDTime");
	man->CreateNtupleIColumn("fParticleID");
	man->CreateNtupleDColumn("fMomOut");
	man->FinishNtuple(0);
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