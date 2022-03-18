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

#include "ScintSD.hh"

EventAction::EventAction(RunAction* runAction) : 
	G4UserEventAction(), fRunAction(runAction), fCollIDScint_gate(-1), fCollIDScint_telescope(-1)
	, fEvID(-1){}

EventAction::~EventAction(){}

void EventAction::BeginOfEventAction(const G4Event*){}

void EventAction::EndOfEventAction(const G4Event* event){
	G4AnalysisManager *man = G4AnalysisManager::Instance();

	// Hits collections
	G4HCofThisEvent*HCE = event->GetHCofThisEvent();
	if(!HCE) return;

	// Get hits collections IDs
	if(fCollIDScint_telescope < 0){
		G4SDManager* SDMan = G4SDManager::GetSDMpointer();
		fCollIDScint_telescope = SDMan->GetCollectionID("Scint_gate/scintCollection");
	}
	
	ScintHitsCollection* ScintHitCollection_telescope = (ScintHitsCollection*) (HCE->GetHC(fCollIDScint_telescope));

	ScintHit* scintHit;
	G4int N = ScintHitCollection_telescope->entries();
	for(int i = 0; i < N; i++){
		scintHit = (*ScintHitCollection_telescope)[i];

		fEvID = event->GetEventID();

		FillScintNtupla(man, scintHit, 1);

		scintHit->Clear();
	}

	if(fEvID % 100 == 0 || (fEvID & (fEvID - 1)) == 0 ) 
	std::cout << "Filled first ntupla" << std::endl;

	//###################################à
	if(fCollIDScint_gate < 0){
		G4SDManager* SDMan = G4SDManager::GetSDMpointer();
		fCollIDScint_gate = SDMan->GetCollectionID("Scint_telescope/scintCollection");
	}
	
	ScintHitsCollection* ScintHitCollection_gate = (ScintHitsCollection*) (HCE->GetHC(fCollIDScint_gate));
	
	N = ScintHitCollection_gate->entries();
	for(int i = 0; i < N; i++){
		scintHit = (*ScintHitCollection_gate)[i];

		fEvID = event->GetEventID();

		FillScintNtupla(man, scintHit, 2);

		scintHit->Clear();
	}


	if(fEvID % 100 == 0 || (fEvID & (fEvID - 1)) == 0 ) 
	{std::cout << "Filled second ntupla" << std::endl;
	std::cout << "Event n. " << fEvID << std::endl;}
}
