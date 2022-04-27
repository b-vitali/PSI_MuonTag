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
G4VSensitiveDetector(name){ //fThetaOut((Double_t)(std::numeric_limits<double>::infinity()))
	fScintCollection = nullptr;
	collectionName.insert("scintCollection");
	hitsCID = -1;
	G4cout<<"create ScintSD: "<<name<<G4endl;
}

ScintSD::ScintSD(G4String name, G4int ntuple):
G4VSensitiveDetector("tmp"){
	G4cout<<"Create a tmp ScintSD to fill the Ntupla"<<G4endl;

	G4AnalysisManager *man = G4AnalysisManager::Instance();

	// Ntuple for the Scintillator
	man->CreateNtuple(name, name);
	// fEvent_ID 		= man->CreateNtupleIColumn("fEvent");
	// fScintNo_ID 	= man->CreateNtupleIColumn("fScintNo");
	// fParticleID_ID	= man->CreateNtupleIColumn("fParticleID");
	// fNgamma_ID 		= man->CreateNtupleIColumn("fNgamma");
	// fNgammaSec_ID 	= man->CreateNtupleIColumn("fNgammaSec");
	// fNCer_ID 		= man->CreateNtupleIColumn("fNCer");

	// fRight_ID 		= man->CreateNtupleIColumn("fRight");
	// fLeft_ID 		= man->CreateNtupleIColumn("fLeft");
	// fDown_ID 		= man->CreateNtupleIColumn("fDown");
	// fUp_ID 			= man->CreateNtupleIColumn("fUp");
	// fBack_ID 		= man->CreateNtupleIColumn("fBack");
	// fFront_ID 		= man->CreateNtupleIColumn("fFront");

	// fEin_ID 		= man->CreateNtupleDColumn("fEin");
	// fEdep_ID 		= man->CreateNtupleDColumn("fEdep");
	// fEout_ID 		= man->CreateNtupleDColumn("fEout");
	// fThetaIn_ID 	= man->CreateNtupleDColumn("fThetaIn");
	// fTrackLength_ID	= man->CreateNtupleDColumn("fTrackLength");
	// fThetaOut_ID 	= man->CreateNtupleDColumn("fThetaOut");
	// fDecayTime_ID 	= man->CreateNtupleDColumn("fDecayTime");

	// fPosInX_ID	= man->CreateNtupleDColumn("fPosInX");
	// fPosInY_ID	= man->CreateNtupleDColumn("fPosInY");
	// fPosInZ_ID	= man->CreateNtupleDColumn("fPosInZ");
	// fMomInX_ID	= man->CreateNtupleDColumn("fMomInX");
	// fMomInY_ID	= man->CreateNtupleDColumn("fMomInY");
	// fMomInZ_ID	= man->CreateNtupleDColumn("fMomInZ");	
	// fTimeIn_ID	= man->CreateNtupleDColumn("fTimeIn");
	// fPosOutX_ID	= man->CreateNtupleDColumn("fPosOutX");
	// fPosOutY_ID	= man->CreateNtupleDColumn("fPosOutY");
	// fPosOutZ_ID	= man->CreateNtupleDColumn("fPosOutZ");
	// fMomOutX_ID	= man->CreateNtupleDColumn("fMomOutX");
	// fMomOutY_ID	= man->CreateNtupleDColumn("fMomOutY");
	// fMomOutZ_ID	= man->CreateNtupleDColumn("fMomOutZ");
	// fTimeOut_ID	= man->CreateNtupleDColumn("fTimeOut");

	man->CreateNtupleIColumn("fEvent",fEvent);
	man->CreateNtupleIColumn("fScintNo",fScintNo);
	man->CreateNtupleIColumn("fParticleID",fParticleID);
	man->CreateNtupleIColumn("fNgamma",fNgamma);
	man->CreateNtupleIColumn("fNgammaSec",fNgammaSec);
	man->CreateNtupleIColumn("fNCer",fNCer);

	man->CreateNtupleIColumn("fRight",fRight);
	man->CreateNtupleIColumn("fLeft",fLeft);
	man->CreateNtupleIColumn("fDown",fDown);
	man->CreateNtupleIColumn("fUp",fUp);
	man->CreateNtupleIColumn("fBack",fBack);
	man->CreateNtupleIColumn("fFront",fFront);

	man->CreateNtupleDColumn("fEin",fEin);
	man->CreateNtupleDColumn("fEdep",fEdep);
	man->CreateNtupleDColumn("fEout",fEout);
	man->CreateNtupleDColumn("fThetaIn",fThetaIn);
	man->CreateNtupleDColumn("fTrackLength",fTrackLength);
	man->CreateNtupleDColumn("fThetaOut",fThetaOut);
	man->CreateNtupleDColumn("fDecayTime",fDecayTime);

	man->CreateNtupleDColumn("fPosInX",fPosInX);
	man->CreateNtupleDColumn("fPosInY",fPosInY);
	man->CreateNtupleDColumn("fPosInZ",fPosInZ);
	man->CreateNtupleDColumn("fMomInX",fMomInX);
	man->CreateNtupleDColumn("fMomInY",fMomInY);
	man->CreateNtupleDColumn("fMomInZ",fMomInZ);	
	man->CreateNtupleDColumn("fTimeIn",fTimeIn);
	man->CreateNtupleDColumn("fPosOutX",fPosOutX);
	man->CreateNtupleDColumn("fPosOutY",fPosOutY);
	man->CreateNtupleDColumn("fPosOutZ",fPosOutZ);
	man->CreateNtupleDColumn("fMomOutX",fMomOutX);
	man->CreateNtupleDColumn("fMomOutY",fMomOutY);
	man->CreateNtupleDColumn("fMomOutZ",fMomOutZ);
	man->CreateNtupleDColumn("fTimeOut",fTimeOut);
	
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
G4int nReflection;
G4bool ScintSD::ProcessHits(G4Step *aStep, G4TouchableHistory* ROhist){	
// no printout == ""; p = pre-info, i = in, o = out, g = gamma, e = else; 
G4String debug	= "p1 i o n e";

if(aStep->GetStepLength() == 0 && aStep->GetTotalEnergyDeposit() == 0) {G4cout<<"step of lenght 0 and edep 0"<<G4endl; return false;}
// if(nReflection%100000==0)G4cout<<"nReflection "<<nReflection<<G4endl;
G4bool TrackOneIn = false; 
	//? Take start and end of the G4Step
	G4StepPoint* preStep = aStep->GetPreStepPoint();
	G4StepPoint* postStep = aStep->GetPostStepPoint();

	//? Take the G4VPhysicalVolume for both start and end
    G4TouchableHistory* thePreTouchable = (G4TouchableHistory*)(preStep->GetTouchable());
    G4VPhysicalVolume* thePrePV = thePreTouchable->GetVolume();
    G4TouchableHistory* thePostTouchable = (G4TouchableHistory*)(postStep->GetTouchable());
    G4VPhysicalVolume* thePostPV = thePostTouchable->GetVolume();

if(debug.contains("p ")) G4cout<<thePrePV->GetName()<<" "<<thePostPV->GetName()<< G4endl;
if(debug.contains("p ")) G4cout<<"track id "<< aStep->GetTrack()->GetTrackID()<< G4endl;
	
	//? If it is the particle I generated
	if(aStep->GetTrack()->GetTrackID() == 1){
nReflection = 0;
TrackOneIn = true;
		//? Debug on StepProcess and TrackStatus
		if (preStep->GetProcessDefinedStep()){
			G4String StepProcessName = preStep->GetProcessDefinedStep()->GetProcessName();
if(debug.contains("p1"))G4cout<<"StepProcessName " <<StepProcessName<<G4endl;	
		} 
		if(aStep->GetTrack()->GetTrackStatus()==fStopAndKill){
if(debug.contains("p1"))G4cout<< "fStopAndKill"<<G4endl;
		}

		//? Saving incoming primary characteristics	
		if(aStep->IsFirstStepInVolume()){	

			//? Track, pdgID, event number, which scintillator, energy, momentum
			G4Track * track = aStep->GetTrack();
    		fEvent.push_back(G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID());
			fParticleID.push_back(track->GetParticleDefinition()->GetPDGEncoding());
			fScintNo.push_back(preStep->GetTouchable()->GetCopyNumber(1)+1);
			fEin.push_back(preStep->GetKineticEnergy());
			fMomInX.push_back(preStep->GetMomentum().getX());
			fMomInY.push_back(preStep->GetMomentum().getY());
			fMomInZ.push_back(preStep->GetMomentum().getZ());

			//? Position and time
			fPosInX.push_back(aStep->GetPreStepPoint()->GetPosition().getX());
			fPosInY.push_back(aStep->GetPreStepPoint()->GetPosition().getY());
			fPosInZ.push_back(aStep->GetPreStepPoint()->GetPosition().getZ());
			fTimeIn.push_back(aStep->GetPreStepPoint()->GetGlobalTime());

			//? Angle (transform the momentum direction to the volume's reference system)
			// G4double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
			// G4double dimensionx = ((G4Box*) thePreTouchable->GetVolume(0)->GetLogicalVolume()->GetSolid())->GetXHalfLength();
			// G4double dimensiony = ((G4Box*) thePreTouchable->GetVolume(0)->GetLogicalVolume()->GetSolid())->GetYHalfLength();
			// G4double dimensionz = ((G4Box*) thePreTouchable->GetVolume(0)->GetLogicalVolume()->GetSolid())->GetZHalfLength();
			G4ThreeVector worldPos = aStep->GetPreStepPoint()->GetPosition();
			G4ThreeVector localPos = thePreTouchable->GetHistory()->GetTopTransform().TransformPoint(worldPos);
			G4AffineTransform momentumTransform = thePreTouchable->GetHistory()->GetTopTransform();
			momentumTransform.SetNetTranslation(G4ThreeVector(0,0,0));
			fDirIn = aStep->GetPreStepPoint()->GetMomentumDirection();
			fDirIn_trans = momentumTransform.TransformPoint(aStep->GetPreStepPoint()->GetMomentumDirection());
if(debug.contains("i+")) G4cout<<"fDirIn [pre trasform] : "<<fDirIn.x()<<" "<<fDirIn.y()<<" "<<fDirIn.z()<<G4endl;
if(debug.contains("i+")) G4cout<<"fDirIn [Volume's reference] :"<<fDirIn_trans.x()<<" "<<fDirIn_trans.y()<<" "<<fDirIn_trans.z()<<G4endl;

			G4ThreeVector norm = -thePreTouchable->GetVolume(0)->GetLogicalVolume()->GetSolid()->SurfaceNormal(localPos);
			fThetaIn.push_back(norm.dot(fDirIn_trans));
if(debug.contains("i+")) G4cout<<"norm: "<<norm.x()<<" "<<norm.y()<<" "<<norm.z()<<G4endl;
if(debug.contains("i")) G4cout<<"cos(fThetaIn) = "<<norm.dot(fDirIn_trans)<<" and fThetaIn [deg] = "<<std::acos(norm.dot(fDirIn_trans)) * 180/CLHEP::pi<<G4endl;

			//? create the entry in the arrays to be filled;
			fEdep.push_back(0);
			fDelta.push_back(0);
			fTrackLength.push_back(0);
			fNgamma.push_back(0);
			fNgammaSec.push_back(0);
			fNCer.push_back(0);
			fRight.push_back(0);
			fLeft.push_back(0);
			fDown.push_back(0);
			fUp.push_back(0);
			fBack.push_back(0);
			fFront.push_back(0);

			//? return only if it is NOT also the last step
			// if (!aStep->IsLastStepInVolume()) return false;
		}
		
		//? Deposited energy, delta and lenght of the track
		G4double edep = aStep->GetTotalEnergyDeposit();
		G4double delta = aStep->GetPostStepPoint()->GetKineticEnergy() - aStep->GetPreStepPoint()->GetKineticEnergy() + edep;
		
		fEdep[fEdep.size()-1] += edep;
		fDelta[fDelta.size()-1] -= delta;
		fTrackLength[fTrackLength.size()-1] += aStep->GetStepLength();
				
		G4double eout = 0;

		//? Eout to see if the particle was stopped
		eout = postStep->GetKineticEnergy();
if(debug.contains("p")) G4cout<<"eout "<< eout <<G4endl;

		//? Stopped particle (set theta to infinity)
		//! DirOut?
		if (eout == 0. && TrackOneIn){
TrackOneIn = false;
if(debug.contains("o")) G4cout<<"Particle stopped!"<<G4endl;
			fEout.push_back(eout);
			fThetaOut.push_back((Double_t)(std::numeric_limits<double>::infinity()));
			// if(fBounce > 0) fBounce += 1;
		
			fMomOutX.push_back(0);
			fMomOutY.push_back(0);
			fMomOutZ.push_back(0);

			fPosOutX.push_back(aStep->GetPostStepPoint()->GetPosition().getX());
			fPosOutY.push_back(aStep->GetPostStepPoint()->GetPosition().getY());
			fPosOutZ.push_back(aStep->GetPostStepPoint()->GetPosition().getZ());
			fTimeOut.push_back(aStep->GetPostStepPoint()->GetGlobalTime());
			
			//? Should I kill the track?
			//aStep->Getrack()->SetTrackStatus(fStopAndKill);
			// return true;
		}

		//? Exiting particle
		else if(postStep->GetStepStatus() == fGeomBoundary  && eout != 0 && TrackOneIn){
TrackOneIn = false;
if(debug.contains("o")) G4cout<<"Particle NOT stopped!"<<G4endl;
			fEout.push_back(eout);
			
			fMomOutX.push_back(postStep->GetMomentum().getX());
			fMomOutY.push_back(postStep->GetMomentum().getY());
			fMomOutZ.push_back(postStep->GetMomentum().getZ());

			fPosOutX.push_back(postStep->GetPosition().getX());
			fPosOutY.push_back(postStep->GetPosition().getY());
			fPosOutZ.push_back(postStep->GetPosition().getZ());
			fTimeOut.push_back(postStep->GetGlobalTime());
			fDirOut = postStep->GetMomentumDirection();
if(debug.contains("o")) G4cout<<"p out x"<<postStep->GetPosition().getZ()<<G4endl;

			//? Angle (transform the momentum direction to the volume's reference system)
			G4ThreeVector worldPos = postStep->GetPosition();
			G4ThreeVector localPos = thePreTouchable->GetHistory()->GetTopTransform().TransformPoint(worldPos);
			G4AffineTransform momentumTransform = thePreTouchable->GetHistory()->GetTopTransform();
			momentumTransform.SetNetTranslation(G4ThreeVector(0,0,0));
			fDirOut_trans = momentumTransform.TransformPoint(postStep->GetMomentumDirection());
if(debug.contains("o+")) G4cout<<"fDirOut [pre trasform] : "<<fDirOut.x()<<" "<<fDirOut.y()<<" "<<fDirOut.z()<<G4endl;
if(debug.contains("o+")) G4cout<<"fDirOut [Volume's reference] :"<<fDirOut_trans.x()<<" "<<fDirOut_trans.y()<<" "<<fDirOut_trans.z()<<G4endl;

			//? If it is the first step filp the norm got from thePreTouchable 
			G4ThreeVector norm = thePreTouchable->GetVolume(0)->GetLogicalVolume()->GetSolid()->SurfaceNormal(localPos);
			fThetaOut.push_back(norm.dot(fDirOut_trans));
			if(fThetaOut[fThetaOut.size()-1]<0){
				if(aStep->IsFirstStepInVolume()) norm = - norm;
				fThetaOut[fThetaOut.size()-1] = norm.dot(fDirOut_trans);
			}
if(debug.contains("o+")) G4cout<<"norm: "<<norm.x()<<" "<<norm.y()<<" "<<norm.z()<<G4endl;
if(debug.contains("o")) G4cout<<"cos(fThetaOut) = "<<norm.dot(fDirOut_trans)<<" and fThetaOut [deg] = "<<std::acos(norm.dot(fDirOut_trans)) * 180/CLHEP::pi<<G4endl;

			
			//? Should I kill the track?
			// aStep->GetTrack()->SetTrackStatus(fStopAndKill);
			// FillHit();
			// return true;
		}

		//? Counting and classifying photons
		const std::vector<const G4Track*>* secondaries = aStep->GetSecondaryInCurrentStep();
		if(secondaries->size() > 0){
			for(unsigned int i = 0; i < secondaries->size(); i++){
				if(secondaries->at(i)->GetParentID() > 0){
					if(secondaries->at(i)->GetDynamicParticle()->GetParticleDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()){
						fNgamma[fNgamma.size()-1]  += 1;
						fNgammaSec[fNgammaSec.size()-1]  += 1;
					}

					//! Better definition for the decay to keep it general purpose
					else if(secondaries->at(i)->GetParticleDefinition()->GetParticleName() == "e+"){
						if(secondaries->at(i)->GetCreatorProcess()->GetProcessName() == "Decay"){
							fDecayTime.push_back(aStep->GetTrack()->GetGlobalTime());
						}
					}
				}
			}
		}
if(debug.contains("p"))G4cout<<"fNgamma " <<fNgamma[fNgamma.size()-1]<<G4endl;	

if(debug.contains("p")) G4cout<<"End of TrackID = 1"<<G4endl;
		return false;
	}

	//? If it is an G4OpticalPhoton
	//! fDirOut_trans = momentumTransform.TransformPoint(postStep->GetMomentumDirection());
	else if(aStep->GetTrack()->GetDynamicParticle()->GetParticleDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()){
		G4VSolid* solid = thePrePV->GetLogicalVolume()->GetSolid();

		G4Box* boxSolid = (G4Box*)(solid);
		//if it is at the border
		if(postStep->GetStepStatus() == fGeomBoundary){
			G4double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
			G4ThreeVector worldPos = postStep->GetPosition();
			G4ThreeVector localpos = thePreTouchable->GetHistory()->GetTopTransform().TransformPoint(worldPos);
			G4double dimensionX = boxSolid->GetXHalfLength();
			G4double dimensionY = boxSolid->GetYHalfLength();
			G4double dimensionZ = boxSolid->GetZHalfLength();
			
			//transform the direction to match the rotation of the faces
			G4AffineTransform momentumTransform = thePreTouchable->GetHistory()->GetTopTransform();
			momentumTransform.SetNetTranslation(G4ThreeVector(0,0,0));
			fDirOut_trans = momentumTransform.TransformPoint(postStep->GetMomentumDirection());
if(debug.contains("g"))G4cout<<"track id "<< aStep->GetTrack()->GetTrackID()<< G4endl;		
if(debug.contains("g"))G4cout<<"pos: "<<localpos.x()<<" "<<localpos.y()<<" "<<localpos.z()<<G4endl;
if(debug.contains("g"))G4cout<<"dir: "<<fDirOut_trans.x()<<" "<<fDirOut_trans.y()<<" "<<fDirOut_trans.z()<<G4endl;


			G4AnalysisManager *man = G4AnalysisManager::Instance();
			// this checks if the photon GOES OUT from the scint
			if(std::fabs(localpos.x() + dimensionX) < kCarTolerance && fDirOut_trans.getX() < 0){
if(debug.contains("g"))G4cout<<"Left"<<G4endl;
				fLeft[fLeft.size()-1]  += 1;
				man->FillH2(0, postStep->GetGlobalTime(), fScintNo[fScintNo.size()-1]);
				man->FillH3(2, localpos.z(), localpos.y(), fScintNo[fScintNo.size()-1]);
				aStep->GetTrack()->SetTrackStatus(fStopAndKill);
			}
			else if(std::fabs(localpos.x() - dimensionX) < kCarTolerance && fDirOut_trans.getX() > 0){
if(debug.contains("g"))G4cout<<"Right"<<G4endl;
				fRight[fRight.size()-1]  += 1;
				man->FillH2(0, postStep->GetGlobalTime(), fScintNo[fScintNo.size()-1]);
				man->FillH3(3, localpos.z(), localpos.y(), fScintNo[fScintNo.size()-1]);
				aStep->GetTrack()->SetTrackStatus(fStopAndKill);
			}
			else if(std::fabs(localpos.y() + dimensionY) < kCarTolerance && fDirOut_trans.getY() < 0){
if(debug.contains("g"))G4cout<<"Down"<<G4endl;
				fDown[fDown.size()-1]  += 1;
				man->FillH2(0, postStep->GetGlobalTime(), fScintNo[fScintNo.size()-1]);
				man->FillH3(5, localpos.z(), localpos.x(), fScintNo[fScintNo.size()-1]);
				aStep->GetTrack()->SetTrackStatus(fStopAndKill);
			}
			else if(std::fabs(localpos.y() - dimensionY) < kCarTolerance && fDirOut_trans.getY() > 0){
if(debug.contains("g"))G4cout<<"Up"<<G4endl;
				fUp[fUp.size()-1]  += 1;
				man->FillH2(0, postStep->GetGlobalTime(), fScintNo[fScintNo.size()-1]);
				man->FillH3(4, localpos.z(), localpos.x(), fScintNo[fScintNo.size()-1]);
				aStep->GetTrack()->SetTrackStatus(fStopAndKill);
			}
			else if(std::fabs(localpos.z() + dimensionZ) < kCarTolerance && fDirOut_trans.getZ() < 0) {
if(debug.contains("g"))G4cout<<"Back"<<G4endl;
				fBack[fBack.size()-1]  += 1;
				man->FillH2(0, postStep->GetGlobalTime(), fScintNo[fScintNo.size()-1]);
				man->FillH3(1, localpos.x(),localpos.y(), fScintNo[fScintNo.size()-1]);
				aStep->GetTrack()->SetTrackStatus(fStopAndKill);
			}
			else if(std::fabs(localpos.z() - dimensionZ) < kCarTolerance && fDirOut_trans.getZ() > 0){
if(debug.contains("g"))G4cout<<"Front"<<G4endl;
				fFront[fFront.size()-1]  += 1;
				man->FillH2(0, postStep->GetGlobalTime(), fScintNo[fScintNo.size()-1]);
				man->FillH3(0, localpos.x(),localpos.y(), fScintNo[fScintNo.size()-1]);
				aStep->GetTrack()->SetTrackStatus(fStopAndKill);
			}
			else {
				nReflection += 1; 
if(debug.contains("g"))G4cout<<"Reflect"<<G4endl;
			}
		// return false;
		}
	}

	//? Everything else: ionization, decay etc
	else{
if(debug.contains("e"))G4cout<<"else particle id "<<aStep->GetTrack()->GetDynamicParticle()->GetParticleDefinition()->GetPDGEncoding()<<G4endl;
		if (preStep->GetProcessDefinedStep()){
			G4String StepProcessName = preStep->GetProcessDefinedStep()->GetProcessName();
if(debug.contains("e"))G4cout<<"StepProcessName " <<StepProcessName<<G4endl;	
		} 
		//? save if it generated photons
		const std::vector<const G4Track*>* secondaries = aStep->GetSecondaryInCurrentStep();
		if(secondaries->size() > 0){
			for(unsigned int i = 0; i < secondaries->size(); i++){
if(debug.contains("e")) G4cout<<secondaries->at(i)->GetDynamicParticle()->GetParticleDefinition()->GetPDGEncoding()<<G4endl;
				if(secondaries->at(i)->GetParentID() > 0){
					if(secondaries->at(i)->GetDynamicParticle()->GetParticleDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()){
						fNgammaSec[fNgammaSec.size()-1]  += 1;
						if(secondaries->at(i)->GetCreatorProcess()->GetProcessName() == "Cerenkov"){
							fNCer[fNCer.size()-1]  += 1;
						}
					}
				}
			}
		}
		// return false;
	}
}

void ScintSD::FillHit(){
	ScintHit* Hit = new ScintHit();

	Hit->SetEvent(fEvent);
	// Hit->SetBounce(fBounce);
	
	Hit->SetScintNo(fScintNo);
	Hit->SetParticleID(fParticleID);
	Hit->SetNgamma(fNgamma);
	Hit->SetNgammaSec(fNgammaSec);
	Hit->SetNCer(fNCer);
	
	Hit->SetCurrentRight(fRight);
	Hit->SetCurrentLeft(fLeft);
	Hit->SetCurrentDown(fDown);
	Hit->SetCurrentUp(fUp);
	Hit->SetCurrentBack(fBack);
	Hit->SetCurrentFront(fFront);

	Hit->SetEin(fEin);
	Hit->SetEdep(fEdep);
	Hit->SetEout(fEout);
	Hit->SetEdelta(fDelta);
	Hit->SetThetaIn(fThetaIn); //std::acos(fThetaIn) * 180/CLHEP::pi
	Hit->SetTrackLength(fTrackLength);
	Hit->SetThetaOut(fThetaOut); //std::acos(fThetaOut) * 180/CLHEP::pi
	Hit->SetDecayTime(fDecayTime);
	
	Hit->SetPosInX(fPosInX);
	Hit->SetPosInY(fPosInY);
	Hit->SetPosInZ(fPosInZ);
	Hit->SetPosOutX(fPosOutX);
	Hit->SetPosOutY(fPosOutY);
	Hit->SetPosOutZ(fPosOutZ);
	Hit->SetMomInX(fMomInX);
	Hit->SetMomInY(fMomInY);
	Hit->SetMomInZ(fMomInZ);
	Hit->SetMomOutX(fMomOutX);
	Hit->SetMomOutY(fMomOutY);
	Hit->SetMomOutZ(fMomOutZ);
	Hit->SetTimeIn(fTimeIn);
	Hit->SetTimeOut(fTimeOut);

	fScintCollection->insert(Hit);

	// fBounce	= 0;

	fEvent.clear();
	fScintNo.clear();
	fParticleID.clear();
	fNgamma.clear();
	fNgammaSec.clear();
	fNCer.clear();

	fRight.clear();
	fLeft.clear();
	fDown.clear();
	fUp.clear();
	fBack.clear();
	fFront.clear();

	fEin.clear();
	fEdep.clear();
	fEout.clear();
	fDelta.clear();
	fThetaIn.clear();
	fTrackLength.clear();
	fThetaOut.clear();
	fDecayTime.clear();

	
	fPosInX.clear();
	fPosInY.clear();
	fPosInZ.clear();
	fPosOutX.clear();
	fPosOutY.clear();
	fPosOutZ.clear();
	fMomInX.clear();
	fMomInY.clear();
	fMomInZ.clear();
	fMomOutX.clear();
	fMomOutY.clear();
	fMomOutZ.clear();
	fTimeIn.clear();
	fTimeOut.clear();
	
	fDirIn_trans = fDirOut_trans = G4ThreeVector();
	fDirIn = fDirOut = G4ThreeVector();
}

void ScintSD::EndOfEvent(G4HCofThisEvent* hitsCE){
	G4cout<<"number of sub hits "<<fEvent.size()<<G4endl;
	FillHit();
}

void ScintSD::clear(){}

void ScintSD::DrawAll(){}

void ScintSD::PrintAll(){}

void ScintSD::FillNtupla(G4AnalysisManager *man, ScintHit* scintHit, G4int ntupla){
//G4cout<<"FillScintNtupla "<<ntupla<<G4endl;
	
	fEvent 	=  scintHit->GetEvent();
	fScintNo 	=  scintHit->GetScintNo();
	fParticleID 	=  scintHit->GetParticleID();
	fNgamma 	=  scintHit->GetNgamma();
	fNgammaSec 	=  scintHit->GetNgammaSec();
	fNCer 	=  scintHit->GetNCer();
	fRight 	=  scintHit->GetCurrentRight();
	fLeft 	=  scintHit->GetCurrentLeft();
	fDown 	=  scintHit->GetCurrentDown();
	fUp 	=  scintHit->GetCurrentUp();
	fBack 	=  scintHit->GetCurrentBack();
	fFront 	=  scintHit->GetCurrentFront();
	fEin 	=  scintHit->GetEin();
	fEdep 	=  scintHit->GetEdep();
	fEout 	=  scintHit->GetEout();
	fThetaIn 	=  scintHit->GetThetaIn();
	fTrackLength 	=  scintHit->GetTrackLength();
	fThetaOut 	=  scintHit->GetThetaOut();
	fDecayTime 	=  scintHit->GetDecayTime();
	fPosInX 	=  scintHit->GetPosInX();
	fPosInY 	=  scintHit->GetPosInY();
	fPosInZ 	=  scintHit->GetPosInZ();
	fMomInX 	=  scintHit->GetMomInX();
	fMomInY 	=  scintHit->GetMomInY();
	fMomInZ 	=  scintHit->GetMomInZ();
	fTimeIn 	=  scintHit->GetTimeIn();
	fPosOutX 	=  scintHit->GetPosOutX();
	fPosOutY 	=  scintHit->GetPosOutY();
	fPosOutZ 	=  scintHit->GetPosOutZ();
	fMomOutX 	=  scintHit->GetMomOutX();
	fMomOutY 	=  scintHit->GetMomOutY();
	fMomOutZ 	=  scintHit->GetMomOutZ();
	fTimeOut 	=  scintHit->GetTimeOut();
	
	// //! Smarter way to fill it?
	// man->FillNtupleIColumn(ntupla, fEvent_ID		, scintHit->GetEvent());
	// man->FillNtupleIColumn(ntupla, fScintNo_ID		, scintHit->GetScintNo());
	// man->FillNtupleIColumn(ntupla, fParticleID_ID	, scintHit->GetParticleID());
    // man->FillNtupleIColumn(ntupla, fNgamma_ID		, scintHit->GetNgamma());
    // man->FillNtupleIColumn(ntupla, fNgammaSec_ID	, scintHit->GetNgammaSec());
    // man->FillNtupleIColumn(ntupla, fNCer_ID			, scintHit->GetNCer());
    // man->FillNtupleIColumn(ntupla, fRight_ID		, scintHit->GetCurrentRight());
    // man->FillNtupleIColumn(ntupla, fLeft_ID		, scintHit->GetCurrentLeft());
    // man->FillNtupleIColumn(ntupla, fDown_ID		, scintHit->GetCurrentDown());
    // man->FillNtupleIColumn(ntupla, fUp_ID		, scintHit->GetCurrentUp());
    // man->FillNtupleIColumn(ntupla, fBack_ID		, scintHit->GetCurrentBack());
    // man->FillNtupleIColumn(ntupla, fFront_ID		, scintHit->GetCurrentFront());
    // man->FillNtupleDColumn(ntupla, fThetaIn_ID		, scintHit->GetThetaIn());
    // man->FillNtupleDColumn(ntupla, fThetaOut_ID		, scintHit->GetThetaOut());
    // man->FillNtupleDColumn(ntupla, fEin_ID			, scintHit->GetEin());
    // man->FillNtupleDColumn(ntupla, fEout_ID			, scintHit->GetEout());
    // man->FillNtupleDColumn(ntupla, fEdep_ID			, scintHit->GetEdep());
    // man->FillNtupleDColumn(ntupla, fTrackLength_ID	, scintHit->GetTrackLength());
	// man->FillNtupleDColumn(ntupla, fDecayTime_ID	, scintHit->GetDecayTime());
    // man->FillNtupleDColumn(ntupla, fPosInX_ID, scintHit->GetPosIn().x());
    // man->FillNtupleDColumn(ntupla, fPosInY_ID, scintHit->GetPosIn().y());
    // man->FillNtupleDColumn(ntupla, fPosInZ_ID, scintHit->GetPosIn().z());
	// man->FillNtupleDColumn(ntupla, fMomInX_ID, scintHit->GetMomIn().x());
    // man->FillNtupleDColumn(ntupla, fMomInY_ID, scintHit->GetMomIn().y());
    // man->FillNtupleDColumn(ntupla, fMomInZ_ID, scintHit->GetMomIn().z());
    // man->FillNtupleDColumn(ntupla, fTimeIn_ID, scintHit->GetTimeIn());
    // man->FillNtupleDColumn(ntupla, fPosOutX_ID, scintHit->GetPosOut().x());
    // man->FillNtupleDColumn(ntupla, fPosOutY_ID, scintHit->GetPosOut().y());
    // man->FillNtupleDColumn(ntupla, fPosOutZ_ID, scintHit->GetPosOut().z());
    // man->FillNtupleDColumn(ntupla, fMomOutX_ID, scintHit->GetMomOut().x());
    // man->FillNtupleDColumn(ntupla, fMomOutY_ID, scintHit->GetMomOut().y());
    // man->FillNtupleDColumn(ntupla, fMomOutZ_ID, scintHit->GetMomOut().z());
    // man->FillNtupleDColumn(ntupla, fTimeOut_ID, scintHit->GetTimeOut());

	man->AddNtupleRow(ntupla);
}