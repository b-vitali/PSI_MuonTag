/// \file  EventAction.cc
/// \brief Implementation of the EventAction class

#include "TTree.h"

#include "RunAction.hh"
#include "EventAction.hh"
#include "ScintHit.hh"
#include "SiPMHit.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"

#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"
#include "G4THitsMap.hh"
#include "G4SystemOfUnits.hh"

EventAction::EventAction(RunAction* runAction) : 
	G4UserEventAction(), fRunAction(runAction), fCollIDScint(-1), fCollIDSiPM(-1)
	, fEvID(-1){
		tmp_scint = new ScintSD("Scint",1);
		//tmp_scint2 = new ScintSD("Scint2",2);
		tmp_sipm = new SiPMSD("SiPM",2);
	}

EventAction::~EventAction(){}

void EventAction::BeginOfEventAction(const G4Event*){}

void EventAction::EndOfEventAction(const G4Event* event){
	G4bool fmuEDM = false;
	// Hits collections
	G4HCofThisEvent*HCE = event->GetHCofThisEvent();
	if(!HCE) return;
	// Get hits collections IDs
	if(fCollIDScint < 0 || fCollIDSiPM < 0){
		G4SDManager* SDMan = G4SDManager::GetSDMpointer();
		fCollIDScint = SDMan->GetCollectionID("Scint/scintCollection");
		fCollIDSiPM = SDMan->GetCollectionID("SiPM/sipmCollection");
		G4cout<<"collections: "<<fCollIDScint<<" and "<<fCollIDSiPM<<G4endl;
	}
	
	ScintHitsCollection* ScintHitCollection = (ScintHitsCollection*) (HCE->GetHC(fCollIDScint));
	
	G4AnalysisManager *man = G4AnalysisManager::Instance();

	ScintHit* scintHit;
	G4int N = ScintHitCollection->entries();
	for(int i = 0; i < N; i++){
		scintHit = (*ScintHitCollection)[i];

		fEvID = event->GetEventID();
		G4cout<<fEvID<<G4endl;
		tmp_scint->FillNtupla(man, scintHit,1);

		scintHit->Clear();
	}

	//! //////////////////////////////////////////////////////////////////////////////////////

	SiPMHitsCollection* SiPMHitCollection = (SiPMHitsCollection*) (HCE->GetHC(fCollIDSiPM));

	SiPMHit* sipmHit;
	G4int M = SiPMHitCollection->entries();
	G4cout<<"entries "<<M<<G4endl;
	for(int i = 0; i < M; i++){
		G4cout<<"LOL"<<G4endl;
		sipmHit = (*SiPMHitCollection)[i];

		fEvID = event->GetEventID();

		tmp_sipm->FillNtupla(man, sipmHit,2);

		sipmHit->Clear();
	}
	
	if(fEvID % 100 == 0 || (fEvID & (fEvID - 1)) == 0 ) 
	std::cout << "Event n. " << fEvID << std::endl;
}
