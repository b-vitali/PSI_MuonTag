/// \file  EventAction.cc
/// \brief Implementation of the EventAction class

#include "TTree.h"

#include "RunAction.hh"
#include "EventAction.hh"
#include "ScintHit.hh"
#include "PixelHit.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"

#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"
#include "G4THitsMap.hh"
#include "G4SystemOfUnits.hh"

EventAction::EventAction(RunAction* runAction) : 
	G4UserEventAction(), fRunAction(runAction), fCollIDScint(-1), 
	fCollIDSiPM(-1), fCollIDSiPMDraw(-1), fEvID(-1){}

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

	if(fCollIDSiPM < 0){
		G4SDManager* SDMan = G4SDManager::GetSDMpointer();
		fCollIDSiPM = SDMan->GetCollectionID("pixelCollection");
	}
	
	PixelHitsCollection* PixelHitCollection = (PixelHitsCollection*) (HCE->GetHC(fCollIDSiPM));

	if(fCollIDSiPMDraw < 0){
		G4SDManager* SDMan = G4SDManager::GetSDMpointer();
		fCollIDSiPMDraw = SDMan->GetCollectionID("pixelCollectionDraw");
	}
	
	PixelHitsCollection* PixelHitCollectionDraw = (PixelHitsCollection*) (HCE->GetHC(fCollIDSiPMDraw));
	G4int N1 = PixelHitCollectionDraw->entries();
	if(PixelHitCollectionDraw){
		for(int i = 0; i < N1; i++){
			(*PixelHitCollectionDraw)[i]->Draw();
		}
	}
	
	
	ScintHit* scintHit;
	G4int N = ScintHitCollection->entries();
	PixelHit* pixelHit;
	assert(N == (int) PixelHitCollection->entries());
	for(int i = 0; i < N; i++){
		scintHit = (*ScintHitCollection)[i];
		pixelHit = (*PixelHitCollection)[i];
		if(fEvID < 0){
			fEvID = event->GetEventID();
			fRunAction->SetEin(scintHit->GetEin());
			fRunAction->SetEdep(scintHit->GetEdep());
			fRunAction->SetEout(scintHit->GetEout());
			fRunAction->SetEdelta(scintHit->GetEdelta());
			fRunAction->SetThetaIn(scintHit->GetThetaIn());
			fRunAction->SetTrackLength(scintHit->GetTrackLength());
			fRunAction->SetThetaPositron(scintHit->GetThetaPositron());
			fRunAction->SetID(fEvID);
			fRunAction->SetBounce(scintHit->GetBounce());
			std::vector<G4double> X, Y, Z, T;
			
			if(fRunAction->GetCmdTracks() > 1 && fRunAction->GetCmdTracks() < int(scintHit->GetPosX().size())){
				G4int npoints = fRunAction->GetCmdTracks();
				G4int ntot = scintHit->GetPosX().size();
				G4int temp = ntot/(npoints - 1);
				for(int j = 0; j < npoints; j++){
					if(j == npoints - 1){
						X.push_back(scintHit->GetPosX().at(ntot - 1));
						Y.push_back(scintHit->GetPosY().at(ntot - 1));
						Z.push_back(scintHit->GetPosZ().at(ntot - 1));
						T.push_back(scintHit->GetTime().at(ntot - 1));
					}
					else{
						X.push_back(scintHit->GetPosX().at(j * temp));
						Y.push_back(scintHit->GetPosY().at(j * temp));
						Z.push_back(scintHit->GetPosZ().at(j * temp));
						T.push_back(scintHit->GetTime().at(j * temp));
					}
				}
			}
			else{
				X = scintHit->GetPosX();
				Y = scintHit->GetPosY();
				Z = scintHit->GetPosZ();
				T = scintHit->GetTime();
			}
			fRunAction->SetPosX(X);
			fRunAction->SetPosY(Y);
			fRunAction->SetPosZ(Z);
			fRunAction->SetTime(T);
			fRunAction->SetNgamma(scintHit->GetNgamma());
			fRunAction->SetNgammaSec(scintHit->GetNgammaSec());
			if(fRunAction->GetCmdPhotons() == 0){
				fRunAction->SetCer(scintHit->GetCer());
				fRunAction->SetThetaGamma(scintHit->GetThetaGamma());
				fRunAction->SetTimeGamma(scintHit->GetTimeGamma());
				fRunAction->SetEGamma(scintHit->GetEGamma());
			}
			fRunAction->SetNCer(scintHit->GetNCer());
			fRunAction->SetCurrentRight(scintHit->GetCurrentRight());
			fRunAction->SetCurrentLeft(scintHit->GetCurrentLeft());
			fRunAction->SetCurrentDown(scintHit->GetCurrentDown());
			fRunAction->SetCurrentUp(scintHit->GetCurrentUp());
			fRunAction->SetCurrentBack(scintHit->GetCurrentBack());
			fRunAction->SetCurrentFront(scintHit->GetCurrentFront());
			fRunAction->SetSiPM(scintHit->GetSiPM());
			fRunAction->SetDecayTime(scintHit->GetDecayTime());

			fRunAction->SetNCells(pixelHit->GetNCells());
			fRunAction->SetNPhotoElectrons(pixelHit->GetNPhotoElectrons());
			fRunAction->SetCells(pixelHit->GetCells());
			fRunAction->SetCellTime(pixelHit->GetCellTime());
			fRunAction->SetOCTFlag(pixelHit->GetOCTFlag());
			fRunAction->SetDNFlag(pixelHit->GetDNFlag());
			(fRunAction->GetTreePtr())->Fill();
		}
		scintHit->Clear();
		pixelHit->Clear();
	}
	if(fEvID % 100 == 0 || (fEvID & (fEvID - 1)) == 0 ) std::cout << "Event n. " << fEvID << std::endl;
	fRunAction->AdvanceGunTime();
	fEvID = -1;
}
