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
ScintSD::ScintSD(G4String name, G4int ntuple) : 
G4VSensitiveDetector(name), fEvent(-1), fScintNo(-1), fParticleID(0),fEin(0), fEdep(0), fEout(0), fDelta(0), fThetaIn(0), 
fTrackLength(0), fThetaOut((Double_t)(std::numeric_limits<double>::infinity())), fBounce(0), fDirIN(G4ThreeVector()), 
fDirOUT(G4ThreeVector()), fNgamma(0), fNgammaSec(0), fNCer(0), fRight(0), fLeft(0), 
fDown(0), fUp(0), fBack(0), fFront(0), fSiPM(0), fDecayTime(-1){
	fScintCollection = nullptr;
	collectionName.insert("scintCollection");
hitsCID = -1;

	G4AnalysisManager *man = G4AnalysisManager::Instance();

	// Ntuple for the Scintillator
	man->CreateNtuple(name, name);
	man->CreateNtupleIColumn("fEvent");
	man->CreateNtupleIColumn("fScintNo");
	man->CreateNtupleIColumn("fParticleID");
	man->CreateNtupleIColumn("fFront");
	man->CreateNtupleDColumn("fThetaIn");
	man->CreateNtupleDColumn("fThetaOut");
	man->CreateNtupleDColumn("fEin");
	man->CreateNtupleDColumn("fEout");
	man->CreateNtupleDColumn("fEdep");

	man->CreateNtupleDColumn("fPosInX");
	man->CreateNtupleDColumn("fPosInY");
	man->CreateNtupleDColumn("fPosInZ");
	man->CreateNtupleDColumn("fTimeIn");
	man->CreateNtupleDColumn("fPosOutX");
	man->CreateNtupleDColumn("fPosOutY");
	man->CreateNtupleDColumn("fPosOutZ");
	man->CreateNtupleDColumn("fTimeOut");

	G4cout<<"Createntupla "<<ntuple<<" for scint "<<name<<G4endl;

	man->FinishNtuple(ntuple);

}

ScintSD::~ScintSD(){}

void ScintSD::Initialize(G4HCofThisEvent* hitsCE){
	fScintCollection = new ScintHitsCollection(SensitiveDetectorName, collectionName[0]);

	// Putting all the hits in the same place
	
	if(hitsCID<0){
		hitsCID = G4SDManager::GetSDMpointer()->GetCollectionID(fScintCollection); 
	}
	hitsCE->AddHitsCollection(hitsCID, fScintCollection);
}
int t = 0;

G4bool ScintSD::ProcessHits(G4Step *aStep, G4TouchableHistory* ROhist){	
// no printout == 0 ; 
int debug	=	0;

	// Take start and end of the G4Step
	G4StepPoint* preStep = aStep->GetPreStepPoint();
	G4StepPoint* postStep = aStep->GetPostStepPoint();

	// Take the G4VPhysicalVolume for both start and end
    G4TouchableHistory* thePreTouchable = (G4TouchableHistory*)(preStep->GetTouchable());
    G4VPhysicalVolume* thePrePV = thePreTouchable->GetVolume();
    G4TouchableHistory* thePostTouchable = (G4TouchableHistory*)(postStep->GetTouchable());
    G4VPhysicalVolume* thePostPV = thePostTouchable->GetVolume();

if(debug>0)std::cout<<thePrePV->GetName()<<" "<<thePostPV->GetName()<< std::endl;

	// Add a check that the step is in the volume we are interested in

	if(aStep->GetTrack()->GetTrackID() == 1){
if(debug>0)std::cout<<"track id "<< aStep->GetTrack()->GetTrackID()<< std::endl;

		G4double edep = aStep->GetTotalEnergyDeposit();
		G4double delta = aStep->GetPostStepPoint()->GetKineticEnergy() - aStep->GetPreStepPoint()->GetKineticEnergy() + edep;

		fEdep += edep;
		fDelta -= delta;
		fTrackLength += aStep->GetStepLength();
				
		G4double ein = 0, eout = 0;

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

		// Saving incoming primary characteristics	
		if(aStep->IsFirstStepInVolume() && fEin == 0){	

			// Track, pdgID, event number, which scintillator
			G4Track * track = aStep->GetTrack();
    		fParticleID = track->GetParticleDefinition()->GetPDGEncoding();
    		fEvent = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
    		fScintNo = thePrePV->GetCopyNo();

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

			fPosIn = aStep->GetPreStepPoint()->GetPosition();
			fTimeIn = aStep->GetPreStepPoint()->GetGlobalTime();

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
			if (!aStep->IsLastStepInVolume()) return false;
		}
		
		eout = postStep->GetKineticEnergy();
if(debug>0)std::cout<<"eout "<< eout <<std::endl;

		if (eout == 0.){
			fEout = eout;
			fDirOUT = preStep->GetMomentumDirection();
			fThetaOut = (Double_t)(std::numeric_limits<double>::infinity());
			if(fBounce > 0) fBounce += 1;
			fPosOut = aStep->GetPostStepPoint()->GetPosition();
			fTimeOut = aStep->GetPostStepPoint()->GetGlobalTime();

			//aStep->Getrack()->SetTrackStatus(fStopAndKill);

			// std::cout<<"STOPPED "<< fThetaOut<<std::endl;
			// std::cout<<"fDirOUT "<< fDirOUT.getX()<<" "<<fDirOUT.getY()<<std::endl;
			// std::cout<<"fDirIN "<< fDirIN.getX()<<" "<<fDirIN.getY()<<std::endl;
			// std::cout<<"fThetaOut "<< fThetaOut<<std::endl;

if(debug>0)std::cout<<"c"<<std::endl;
			FillHit();
			return true;
		}

		else if(postStep->GetStepStatus() == fGeomBoundary){
			fEout = eout;
			fDirOUT = postStep->GetMomentumDirection();
if(debug>0)std::cout<<fDirOUT.x()<<" "<<fDirOUT.y()<<" "<<fDirOUT.z()<<std::endl;
if(debug>0)std::cout<<fDirIN.x()<<" "<<fDirIN.y()<<" "<<fDirIN.z()<<std::endl;

			fThetaOut = fDirIN.dot(fDirOUT);
if(debug>0)std::cout<<"boh "<<fThetaOut<<std::endl;
						//fThetaOut = (Double_t)(std::numeric_limits<double>::infinity());

			fBounce += 1;
			fPosOut = aStep->GetPostStepPoint()->GetPosition();
			fTimeOut = aStep->GetPostStepPoint()->GetGlobalTime();
			//aStep->GetTrack()->SetTrackStatus(fStopAndKill);

if(debug>0)std::cout<<"NOT STOPPED "<< fThetaOut<<std::endl;
if(debug>0)std::cout<<"fDirOUT "<< fDirOUT.getX()<<" "<<fDirOUT.getY()<<std::endl;
if(debug>0)std::cout<<"fDirIN "<< fDirIN.getX()<<" "<<fDirIN.getY()<<std::endl;
if(debug>0)std::cout<<"fThetaOut "<< fThetaOut<<std::endl;

if(debug>0)std::cout<<"d"<<std::endl;
			FillHit();
			return true;
		}
if(debug>0)std::cout<<"e"<<std::endl;
		return false;
	}
	else if(aStep->GetTrack()->GetDynamicParticle()->GetParticleDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()){
		G4StepPoint* preStep = aStep->GetPreStepPoint();
		G4StepPoint* postStep = aStep->GetPostStepPoint();
		G4VPhysicalVolume* physVol = preStep->GetPhysicalVolume();
		G4VSolid* solid = physVol->GetLogicalVolume()->GetSolid();

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

				//std::cout<<secondaries->at(i)->GetDynamicParticle()->GetParticleDefinition()<<std::endl;

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

void ScintSD::FillHit(){
	ScintHit* Hit = new ScintHit();
	
	Hit->SetEvent(fEvent);
	Hit->SetScintNo(fScintNo);
	Hit->SetParticleID(fParticleID);
	
	Hit->SetEin(fEin);
	Hit->SetEdep(fEdep);
	Hit->SetEout(fEout);
	Hit->SetEdelta(fDelta);
	Hit->SetThetaIn(fThetaIn * 180/CLHEP::pi);
	Hit->SetTrackLength(fTrackLength);
	Hit->SetThetaOut(std::acos(fThetaOut) * 180/CLHEP::pi);
	Hit->SetBounce(fBounce);
	Hit->SetPosIn(fPosIn);
	Hit->SetTimeIn(fTimeIn);
	Hit->SetPosOut(fPosOut);
	Hit->SetTimeOut(fTimeOut);
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

	fEvent = 0;
	fScintNo = 0;
	fParticleID = 0;

	fEdep = 0;
	fEin = 0;
	fDelta = 0;
	fTrackLength = 0;
	fBounce = 0;
	fDirIN = fDirOUT = G4ThreeVector();
fPosIn = CLHEP::Hep3Vector();
fTimeIn = 0;
fPosOut = CLHEP::Hep3Vector();
fTimeOut = 0;
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

void ScintSD::EndOfEvent(G4HCofThisEvent* hitsCE){

}

void ScintSD::clear(){}

void ScintSD::DrawAll(){}

void ScintSD::PrintAll(){}

void FillScintNtupla(G4AnalysisManager *man, ScintHit* scintHit, G4int ntupla){
//G4cout<<"FillScintNtupla "<<ntupla<<G4endl;

	man->FillNtupleIColumn(ntupla, 0, scintHit->GetEvent());
	man->FillNtupleIColumn(ntupla, 1, scintHit->GetScintNo());
	man->FillNtupleIColumn(ntupla, 2, scintHit->GetParticleID());
    man->FillNtupleIColumn(ntupla, 3, scintHit->GetCurrentFront());
    man->FillNtupleDColumn(ntupla, 4, scintHit->GetThetaIn());
    man->FillNtupleDColumn(ntupla, 5, scintHit->GetThetaOut());
    man->FillNtupleDColumn(ntupla, 6, scintHit->GetEin());
    man->FillNtupleDColumn(ntupla, 7, scintHit->GetEout());
    man->FillNtupleDColumn(ntupla, 8, scintHit->GetEdep());
    man->FillNtupleDColumn(ntupla, 9,  scintHit->GetPosIn().x());
    man->FillNtupleDColumn(ntupla, 10, scintHit->GetPosIn().y());
    man->FillNtupleDColumn(ntupla, 11, scintHit->GetPosIn().z());
    man->FillNtupleDColumn(ntupla, 12, scintHit->GetTimeIn());
    man->FillNtupleDColumn(ntupla, 13, scintHit->GetPosOut().x());
    man->FillNtupleDColumn(ntupla, 14, scintHit->GetPosOut().y());
    man->FillNtupleDColumn(ntupla, 15, scintHit->GetPosOut().z());
    man->FillNtupleDColumn(ntupla, 16, scintHit->GetTimeOut());

	man->AddNtupleRow(ntupla);
}