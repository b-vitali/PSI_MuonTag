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
	G4UserEventAction(), fRunAction(runAction), fCollIDScint_out(-1), fCollIDScint_in(-1)
	, fEvID(-1){
		tmp_scint_out = new ScintSD("Scint_out",1);
		tmp_scint_in = new ScintSD("Scint_in",2);

	}

EventAction::~EventAction(){}

void EventAction::BeginOfEventAction(const G4Event*){}

void EventAction::EndOfEventAction(const G4Event* event){
	// Hits collections
	G4HCofThisEvent*HCE = event->GetHCofThisEvent();
	if(!HCE) return;

	// Get hits collections IDs
	if(fCollIDScint_out < 0 && fCollIDScint_in < 0 ){
		G4SDManager* SDMan = G4SDManager::GetSDMpointer();
		fCollIDScint_out = SDMan->GetCollectionID("Scint_out/scintCollection");
		fCollIDScint_in = SDMan->GetCollectionID("Scint_in/scintCollection");
	}
	else{
		G4cout<<"Something is wrong with the hit collections"<<G4endl;
	}
	
	ScintHitsCollection* ScintHitCollection_out = (ScintHitsCollection*) (HCE->GetHC(fCollIDScint_out));
	ScintHitsCollection* ScintHitCollection_in = (ScintHitsCollection*) (HCE->GetHC(fCollIDScint_in));
	
	G4AnalysisManager *man = G4AnalysisManager::Instance();

	ScintHit* scintHit_out;
	G4int N_out = ScintHitCollection_out->entries();
	for(int i = 0; i < N_out; i++){
		scintHit_out = (*ScintHitCollection_out)[i];

		fEvID = event->GetEventID();

		tmp_scint_out->FillNtupla(man, scintHit_out,1);

		scintHit_out->Clear();
	}

	ScintHit* scintHit_in;
	G4int N_in = ScintHitCollection_in->entries();
	for(int i = 0; i < N_in; i++){
		scintHit_in = (*ScintHitCollection_in)[i];

		fEvID = event->GetEventID();

		tmp_scint_in->FillNtupla(man, scintHit_in,2);

		scintHit_in->Clear();
	}

	if(fEvID % 100 == 0 || (fEvID & (fEvID - 1)) == 0 ) 
	std::cout << "Event n. " << fEvID << std::endl;
}
