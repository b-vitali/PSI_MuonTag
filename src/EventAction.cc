/// \file  EventAction.cc
/// \brief Implementation of the EventAction class

#include "TTree.h"

#include "RunAction.hh"
#include "EventAction.hh"
#include "ScintHit.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"

#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"
#include "G4THitsMap.hh"
#include "G4SystemOfUnits.hh"

EventAction::EventAction(RunAction* runAction) : 
	G4UserEventAction(), fRunAction(runAction), fCollIDScint(-1)
	, fEvID(-1){}

EventAction::~EventAction(){}

void EventAction::BeginOfEventAction(const G4Event*){}

void EventAction::EndOfEventAction(const G4Event* event){
	// Hits collections
	G4HCofThisEvent*HCE = event->GetHCofThisEvent();
	if(!HCE) return;
	// Get hits collections IDs
	if(fCollIDScint < 0){
		G4SDManager* SDMan = G4SDManager::GetSDMpointer();
		fCollIDScint = SDMan->GetCollectionID("scintCollection");
	}
	
	ScintHitsCollection* ScintHitCollection = (ScintHitsCollection*) (HCE->GetHC(fCollIDScint));
	
	G4AnalysisManager *man = G4AnalysisManager::Instance();

	ScintHit* scintHit;
	G4int N = ScintHitCollection->entries();
	for(int i = 0; i < N; i++){
		scintHit = (*ScintHitCollection)[i];

		fEvID = event->GetEventID();
		man->FillNtupleIColumn(1, 0, scintHit->GetEvent());
		man->FillNtupleIColumn(1, 1, scintHit->GetScintNo());
		man->FillNtupleIColumn(1, 2, scintHit->GetParticleID());
       	man->FillNtupleIColumn(1, 3, scintHit->GetCurrentFront());
       	man->FillNtupleDColumn(1, 4, scintHit->GetThetaIn());
       	man->FillNtupleDColumn(1, 5, scintHit->GetThetaOut());
       	man->FillNtupleDColumn(1, 6, scintHit->GetEin());
       	man->FillNtupleDColumn(1, 7, scintHit->GetEout());
       	man->FillNtupleDColumn(1, 8, scintHit->GetEdep());

       	man->FillNtupleDColumn(1, 9,  scintHit->GetPosIn().x());
       	man->FillNtupleDColumn(1, 10, scintHit->GetPosIn().y());
       	man->FillNtupleDColumn(1, 11, scintHit->GetPosIn().z());
       	man->FillNtupleDColumn(1, 12, scintHit->GetTimeIn());
       	man->FillNtupleDColumn(1, 13, scintHit->GetPosOut().x());
       	man->FillNtupleDColumn(1, 14, scintHit->GetPosOut().y());
       	man->FillNtupleDColumn(1, 15, scintHit->GetPosOut().z());
       	man->FillNtupleDColumn(1, 16, scintHit->GetTimeOut());


		man->AddNtupleRow(1);

		scintHit->Clear();
	}
	if(fEvID % 100 == 0 || (fEvID & (fEvID - 1)) == 0 ) 
	std::cout << "Event n. " << fEvID << std::endl;
}
