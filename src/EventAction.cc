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
	G4UserEventAction(), fRunAction(runAction), fCollIDScint_out(-1), fCollIDScint_in(-1), fCollIDSiPM_in(-1),fCollIDSiPM_out(-1)
	, fEvID(-1){
		// tmp_scint_out 	= new ScintSD("Scint_out",1);
		// tmp_scint_in 	= new ScintSD("Scint_in",2);
		// tmp_sipm_out 	= new SiPMSD("SiPM_out",2);
		// tmp_sipm_in 	= new SiPMSD("SiPM_in",3);
	}

EventAction::~EventAction(){}

void EventAction::BeginOfEventAction(const G4Event*){}

void EventAction::EndOfEventAction(const G4Event* event){
	
	// Hits collections
	G4HCofThisEvent*HCE = event->GetHCofThisEvent();
	if(!HCE) return;
	
	fRunAction->GetTmpCyFi()->FillNTuples(HCE, event);
/*
	// // Get hits collections IDs
	// if(fCollIDScint_out < 0 && fCollIDScint_in < 0 && fCollIDSiPM_out < 0 && fCollIDSiPM_in < 0 ){
		G4SDManager* SDMan = G4SDManager::GetSDMpointer();
	// 	fCollIDScint_out = SDMan->GetCollectionID("Scint_out/scintCollection");
	// 	fCollIDScint_in = SDMan->GetCollectionID("Scint_in/scintCollection");
		fCollIDSiPM_out = SDMan->GetCollectionID("SiPM_out/sipmCollection");
		fCollIDSiPM_in = SDMan->GetCollectionID("SiPM_in/sipmCollection");
	// }
	// else{
	// 	G4cout<<"Something is wrong with the hit collections"<<G4endl;
	// }
	
	// ScintHitsCollection* ScintHitCollection_out = (ScintHitsCollection*) (HCE->GetHC(fCollIDScint_out));
	// ScintHitsCollection* ScintHitCollection_in = (ScintHitsCollection*) (HCE->GetHC(fCollIDScint_in));
	
	G4AnalysisManager *man = G4AnalysisManager::Instance();

	// ScintHit* scintHit_out;
	// G4int N_out = ScintHitCollection_out->entries();
	// for(int i = 0; i < N_out; i++){
	// 	scintHit_out = (*ScintHitCollection_out)[i];

	// 	fEvID = event->GetEventID();

	// 	tmp_scint_out->FillNtupla(man, scintHit_out,1);

	// 	scintHit_out->Clear();
	// }

	// ScintHit* scintHit_in;
	// G4int N_in = ScintHitCollection_in->entries();
	// for(int i = 0; i < N_in; i++){
	// 	scintHit_in = (*ScintHitCollection_in)[i];

	// 	fEvID = event->GetEventID();

	// 	tmp_scint_in->FillNtupla(man, scintHit_in,2);

	// 	scintHit_in->Clear();
	// }

	SiPMHitsCollection* SiPMHitCollection_out = (SiPMHitsCollection*) (HCE->GetHC(fCollIDSiPM_out));

	SiPMHit* sipmHit;
	G4int M_out = SiPMHitCollection_out->entries();
	G4cout<<"entries "<<M_out<<G4endl;
	for(int i = 0; i < M_out; i++){
		sipmHit = (*SiPMHitCollection_out)[i];
		G4cout<<"subhits : "<<sipmHit->GetEvent().size()<<G4endl;
		fEvID = event->GetEventID();

		tmp_sipm_out->FillNtupla(man, sipmHit,2);

		sipmHit->Clear();
	}

	SiPMHitsCollection* SiPMHitCollection_in = (SiPMHitsCollection*) (HCE->GetHC(fCollIDSiPM_in));

	G4int M_in = SiPMHitCollection_in->entries();
	G4cout<<"entries "<<M_in<<G4endl;
	for(int i = 0; i < M_in; i++){
		sipmHit = (*SiPMHitCollection_in)[i];

		fEvID = event->GetEventID();

		tmp_sipm_in->FillNtupla(man, sipmHit,3);

		sipmHit->Clear();
	}

	if(fEvID % 100 == 0 || (fEvID & (fEvID - 1)) == 0 ) 
	std::cout << "Event n. " << fEvID << std::endl;
*/
}
