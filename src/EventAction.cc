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
	, fEvID(-1){
		tmp_scint = new ScintSD("Scint",1);
		tmp_scint2 = new ScintSD("Scint2",2);

	}

EventAction::~EventAction(){}

void EventAction::BeginOfEventAction(const G4Event*){}

void EventAction::EndOfEventAction(const G4Event* event){
	G4bool fmuEDM = false;
	// Hits collections
	G4HCofThisEvent*HCE = event->GetHCofThisEvent();
	if(!HCE) return;
	// Get hits collections IDs
	if(fCollIDScint < 0){
		G4SDManager* SDMan = G4SDManager::GetSDMpointer();
		fCollIDScint = SDMan->GetCollectionID("Scint/scintCollection");
	}
	
	ScintHitsCollection* ScintHitCollection = (ScintHitsCollection*) (HCE->GetHC(fCollIDScint));
	
	G4AnalysisManager *man = G4AnalysisManager::Instance();

	ScintHit* scintHit;
	G4int N = ScintHitCollection->entries();
	for(int i = 0; i < N; i++){
		scintHit = (*ScintHitCollection)[i];

		fEvID = event->GetEventID();

		tmp_scint->FillNtupla(man, scintHit,1);

		scintHit->Clear();
	}

	// if(!fmuEDM){ 
		// G4SDManager* SDMan = G4SDManager::GetSDMpointer();
		// fCollIDScint2 = SDMan->GetCollectionID("Scint2/scintCollection");
		// ScintHitsCollection* ScintHitCollection2 = (ScintHitsCollection*) (HCE->GetHC(fCollIDScint2));
		// G4int M = ScintHitCollection2->entries();
		// G4cout<<"N "<<N<<" M "<<M<<G4endl;

		// for(int i = 0; i < M; i++){
			// scintHit = (*ScintHitCollection2)[i];

			// fEvID = event->GetEventID();

			// tmp_scint2->FillNtupla(man, scintHit,2);

			// scintHit->Clear();
		// }
	// }
	if(fEvID % 100 == 0 || (fEvID & (fEvID - 1)) == 0 ) 
	std::cout << "Event n. " << fEvID << std::endl;
}
