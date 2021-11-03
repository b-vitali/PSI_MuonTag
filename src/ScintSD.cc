/// \file  ScinSD.cc
/// \brief Implementation of the ScintSD class

#include "ScintSD.hh"
#include "ScintHit.hh"
#include "RunAction.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4VTouchable.hh"
#include "G4ios.hh"
#include "G4VProcess.hh"
#include "G4RunManager.hh"

#include "G4Box.hh"
#include "G4GeometryTolerance.hh"
#include "G4VSolid.hh"


#include "TVector3.h"

ScintSD::ScintSD(G4String name) : 
G4VSensitiveDetector(name), fEin(0), fEdep(0), fEout(0), fDelta(0), fThetaIn(0), 
fTrackLength(0), fThetaPositron(0), fBounce(0), fDirIN(G4ThreeVector()), 
fDirOUT(G4ThreeVector()), fNgamma(0), fNgammaSec(0), fNCer(0), fRight(0), fLeft(0), 
fDown(0), fUp(0), fBack(0), fFront(0), fSiPM(0), fDecayTime(-1){
	fScintCollection = nullptr;
	collectionName.insert("scintCollection");
}

ScintSD::~ScintSD(){}

void ScintSD::Initialize(G4HCofThisEvent* hitsCE){
	fScintCollection = new ScintHitsCollection(SensitiveDetectorName, collectionName[0]);

	// Putting all the hits in the same place
	static G4int hitsCID = -1;
	if(hitsCID<0){
		hitsCID = GetCollectionID(0);
	}
	hitsCE->AddHitsCollection(hitsCID, fScintCollection);
	RunAction* runAction = (RunAction*) G4RunManager::GetRunManager()->GetUserRunAction();
	fPhotonsCmd = runAction->GetCmdPhotons();
	fTracksCmd = runAction->GetCmdTracks();
}

///aStep->GetTrack()->GetParticleDefinition()->GetParticleName() == "e+" && 
G4bool ScintSD::ProcessHits(G4Step *aStep, G4TouchableHistory*){
	if(aStep->GetTrack()->GetTrackID() == 1){
		G4double edep = aStep->GetTotalEnergyDeposit();
		G4double delta = aStep->GetPostStepPoint()->GetKineticEnergy() - aStep->GetPreStepPoint()->GetKineticEnergy() + edep;
		
		G4TouchableHistory* thePreTouchable = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
		G4VPhysicalVolume* thePrePV = thePreTouchable->GetVolume();
		G4TouchableHistory* thePostTouchable = (G4TouchableHistory*)(aStep->GetPostStepPoint()->GetTouchable());
		G4VPhysicalVolume* thePostPV = thePostTouchable->GetVolume();
		if(thePrePV->GetName() != "Crystal" && thePostPV->GetName() != "Crystal"){
			return false;
		} 

		fEdep += edep;
		fDelta -= delta;
		fTrackLength += aStep->GetStepLength();
		
		G4ThreeVector temppos = aStep->GetPreStepPoint()->GetPosition();
		
		if(fTracksCmd != 1 || fPosX.size() == 0){
			fPosX.push_back(temppos.getX());
			fPosY.push_back(temppos.getY());
			fPosZ.push_back(temppos.getZ());
			fTime.push_back(aStep->GetPreStepPoint()->GetGlobalTime());
		}
		G4double ein = 0, eout = 0;
		G4StepPoint* preStep = aStep->GetPreStepPoint();
		G4StepPoint* postStep = aStep->GetPostStepPoint();

		// Counting and classifying photons
		const std::vector<const G4Track*>* secondaries = aStep->GetSecondaryInCurrentStep();
		if(secondaries->size() > 0){
			for(unsigned int i = 0; i < secondaries->size(); i++){
				if(secondaries->at(i)->GetParentID() > 0){
					if(secondaries->at(i)->GetDynamicParticle()->GetParticleDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()){
						fNgamma += 1;
						fNgammaSec += 1;
						if(secondaries->at(i)->GetCreatorProcess()->GetProcessName() == "Cerenkov"){
							if(fPhotonsCmd == 0) fCer.push_back(1);
							fNCer += 1;
						}
						else{
							if(fPhotonsCmd == 0) fCer.push_back(0);
						}
						
						if(fPhotonsCmd == 0){
							fThetaGamma.push_back(secondaries->at(i)->GetDynamicParticle()->GetMomentumDirection().dot(aStep->GetPreStepPoint()->GetMomentumDirection()));
							fTimeGamma.push_back(secondaries->at(i)->GetGlobalTime());
							fEGamma.push_back(secondaries->at(i)->GetKineticEnergy());
						}
					}
					else if(secondaries->at(i)->GetParticleDefinition()->GetParticleName() == "e+"){
						if(secondaries->at(i)->GetCreatorProcess()->GetProcessName() == "Decay"){
							fDecayTime = aStep->GetTrack()->GetGlobalTime();
						}
					}
				}
			}
		}
		// Saving e+ characteristics				
		if(aStep->IsFirstStepInVolume() && fEin == 0){
			ein = preStep->GetKineticEnergy();
			fEin = ein;
			G4double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
			G4double dimensionx = ((G4Box*) thePreTouchable->GetVolume(0)->GetLogicalVolume()->GetSolid())->GetXHalfLength();
			G4double dimensiony = ((G4Box*) thePreTouchable->GetVolume(0)->GetLogicalVolume()->GetSolid())->GetYHalfLength();
			G4double dimensionz = ((G4Box*) thePreTouchable->GetVolume(0)->GetLogicalVolume()->GetSolid())->GetZHalfLength();
			G4ThreeVector worldPos = aStep->GetPreStepPoint()->GetPosition();
			G4ThreeVector localPos = thePreTouchable->GetHistory()->GetTopTransform().TransformPoint(worldPos);
			G4AffineTransform momentumTransform = thePreTouchable->GetHistory()->GetTopTransform();
			momentumTransform.SetNetTranslation(G4ThreeVector(0,0,0));
			G4ThreeVector momentumDir = momentumTransform.TransformPoint(aStep->GetPreStepPoint()->GetMomentumDirection());
			fDirIN = momentumDir;
			

			/// Surfaces:
			///          - x>0: DS right. Surf = 0
			///          - x<0: DS left.  Surf = 1
			///          - y>0: DS up.    Surf = 2
			///          - y<0: DS down.  Surf = 3
			///          - z>0: DS front. Surf = 4
			///          - z<0: DS back.  Surf = 5
			if(std::fabs(localPos.x() - dimensionx) < kCarTolerance &&
			   momentumDir.x() < 0){
				fThetaIn = acos(-momentumDir.x());
			}
			else if(std::fabs(localPos.x() + dimensionx) < kCarTolerance &&
			   momentumDir.x() > 0){
				fThetaIn = acos(momentumDir.x());
			}
			else if(std::fabs(localPos.y() - dimensiony) < kCarTolerance &&
			   momentumDir.y() < 0){
				fThetaIn = acos(-momentumDir.y());
			}
			else if(std::fabs(localPos.y() + dimensiony) < kCarTolerance &&
			   momentumDir.y() > 0){
				fThetaIn = acos(momentumDir.y());
			}
			else if(std::fabs(localPos.z() - dimensionz) < kCarTolerance &&
			   momentumDir.z() < 0){
				fThetaIn = acos(-momentumDir.z());
			}
			else if(std::fabs(localPos.z() + dimensionz) < kCarTolerance &&
			   momentumDir.z() > 0){
				fThetaIn = acos(momentumDir.z());
			}

			return false;
		}
		
		eout = postStep->GetKineticEnergy();
		if(postStep->GetStepStatus() == fGeomBoundary){
			fEout = eout;
			fDirOUT = postStep->GetMomentumDirection();
			fThetaPositron = fDirIN.dot(fDirOUT);
			fBounce += 1;
			temppos = aStep->GetPostStepPoint()->GetPosition();
			fPosX.push_back(temppos.getX());
			fPosY.push_back(temppos.getY());
			fPosZ.push_back(temppos.getZ());
			fTime.push_back(aStep->GetPostStepPoint()->GetGlobalTime());
			//aStep->GetTrack()->SetTrackStatus(fStopAndKill);
			return true;
		}
		else if (eout == 0.){
			fEout = eout;
			fDirOUT = preStep->GetMomentumDirection();
			fThetaPositron = fDirIN.dot(fDirOUT);
			if(fBounce > 0) fBounce += 1;
			temppos = aStep->GetPreStepPoint()->GetPosition();
			fPosX.push_back(temppos.getX());
			fPosY.push_back(temppos.getY());
			fPosZ.push_back(temppos.getZ());
			fTime.push_back(aStep->GetPostStepPoint()->GetGlobalTime());
			//aStep->Getrack()->SetTrackStatus(fStopAndKill);
			return true;
		}
		
		return false;
	}
	else if(aStep->GetTrack()->GetDynamicParticle()->GetParticleDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()){
		G4StepPoint* preStep = aStep->GetPreStepPoint();
		G4StepPoint* postStep = aStep->GetPostStepPoint();
		G4VPhysicalVolume* physVol = preStep->GetPhysicalVolume();
		G4VSolid* solid = physVol->GetLogicalVolume()->GetSolid();
		if(physVol->GetName() != "Crystal" || postStep->GetPhysicalVolume()->GetName() == "Crystal") return false;

		G4Box* boxSolid = (G4Box*)(solid);
		if(postStep->GetStepStatus() == fGeomBoundary){
			G4TouchableHandle theTouchable = aStep->GetPreStepPoint()->GetTouchableHandle();
			G4double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
			G4ThreeVector stppos = postStep->GetPosition();
			G4ThreeVector localpos = theTouchable->GetHistory()->GetTopTransform().TransformPoint(stppos);
			G4double dimensionX = boxSolid->GetXHalfLength();
			G4double dimensionY = boxSolid->GetYHalfLength();
			G4double dimensionZ = boxSolid->GetZHalfLength();
			
			if(std::fabs(localpos.x() + dimensionX) < kCarTolerance && postStep->GetMomentumDirection().getX() < 0){
				fLeft += 1;
				aStep->GetTrack()->SetTrackStatus(fStopAndKill);
			}
			else if(std::fabs(localpos.x() - dimensionX) < kCarTolerance && postStep->GetMomentumDirection().getX() > 0){
				fRight += 1;
				aStep->GetTrack()->SetTrackStatus(fStopAndKill);
			}
			else if(std::fabs(localpos.y() + dimensionY) < kCarTolerance && postStep->GetMomentumDirection().getY() < 0){
				fDown += 1;
				aStep->GetTrack()->SetTrackStatus(fStopAndKill);
			}
			else if(std::fabs(localpos.y() - dimensionY) < kCarTolerance && postStep->GetMomentumDirection().getY() > 0){
				fUp += 1;
				aStep->GetTrack()->SetTrackStatus(fStopAndKill);
			}
			else if(std::fabs(localpos.z() + dimensionZ) < kCarTolerance && postStep->GetMomentumDirection().getZ() < 0) {
				fBack += 1; 
				if (fabs(localpos.x()) < 0.65*CLHEP::mm && fabs(localpos.y()) < 0.65*CLHEP::mm) fSiPM += 1;
				else{
				aStep->GetTrack()->SetTrackStatus(fStopAndKill);
				}
			}
			else if(std::fabs(localpos.z() - dimensionZ) < kCarTolerance && postStep->GetMomentumDirection().getZ() > 0){
				fFront += 1;
				aStep->GetTrack()->SetTrackStatus(fStopAndKill);
			}
		}
		return false;
	}

	else{
		const std::vector<const G4Track*>* secondaries = aStep->GetSecondaryInCurrentStep();
		if(secondaries->size() > 0){
			for(unsigned int i = 0; i < secondaries->size(); i++){
				if(secondaries->at(i)->GetParentID() > 0){
					if(secondaries->at(i)->GetDynamicParticle()->GetParticleDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()){
						fNgammaSec += 1;
						if(secondaries->at(i)->GetCreatorProcess()->GetProcessName() == "Cerenkov"){
							if(fPhotonsCmd == 0) fCer.push_back(1);
							fNCer += 1;
						}
						else{
							if(fPhotonsCmd == 0) fCer.push_back(0);
						}
						
						if(fPhotonsCmd == 0) {
							fThetaGamma.push_back(secondaries->at(i)->GetDynamicParticle()->GetMomentumDirection().dot(aStep->GetPreStepPoint()->GetMomentumDirection()));
							fTimeGamma.push_back(secondaries->at(i)->GetGlobalTime());
							fEGamma.push_back(secondaries->at(i)->GetKineticEnergy());
						}
					}
				}
			}
		}
		return false;
	}

}

void ScintSD::EndOfEvent(G4HCofThisEvent*){
	ScintHit* Hit = new ScintHit();
	Hit->SetEin(fEin);
	Hit->SetEdep(fEdep);
	Hit->SetEout(fEout);
	Hit->SetEdelta(fDelta);
	Hit->SetThetaIn(fThetaIn);
	Hit->SetTrackLength(fTrackLength);
	Hit->SetThetaPositron(std::acos(fThetaPositron));
	Hit->SetBounce(fBounce);
	Hit->SetPosX(fPosX);
	Hit->SetPosY(fPosY);
	Hit->SetPosZ(fPosZ);
	Hit->SetTime(fTime);
	Hit->SetNgamma(fNgamma);
	Hit->SetNgammaSec(fNgammaSec);
	Hit->SetCer(fCer);
	Hit->SetThetaGamma(fThetaGamma);
	Hit->SetTimeGamma(fTimeGamma);
	Hit->SetEGamma(fEGamma);
	Hit->SetNCer(fNCer);
	Hit->SetCurrentRight(fRight);
	Hit->SetCurrentLeft(fLeft);
	Hit->SetCurrentDown(fDown);
	Hit->SetCurrentUp(fUp);
	Hit->SetCurrentBack(fBack);
	Hit->SetCurrentFront(fFront);
	Hit->SetSiPM(fSiPM);
	Hit->SetDecayTime(fDecayTime);
	fScintCollection->insert(Hit);
	fEdep = 0;
	fEin = 0;
	fDelta = 0;
	fTrackLength = 0;
	fBounce = 0;
	fDirIN = fDirOUT = G4ThreeVector();
	fPosX.clear();
	fPosY.clear();
	fPosZ.clear();
	fTime.clear();
	fNgamma = 0;
	fNgammaSec = 0;
	fThetaGamma.clear();
	fTimeGamma.clear();
	fEGamma.clear();
	fNCer = 0;
	fNCer = 0;
	fRight = 0;
	fLeft = 0;
	fDown = 0;
	fUp = 0;
	fBack = 0;
	fFront = 0;
	fSiPM = 0;
	fDecayTime = -1;
	fCer.clear();
}

void ScintSD::clear(){}

void ScintSD::DrawAll(){}

void ScintSD::PrintAll(){}