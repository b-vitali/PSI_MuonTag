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
G4VSensitiveDetector(name), fEvent(-1), fScintNo(-1), fParticleID(0),fEin(0), fEdep(0), fEout(0), fDelta(0), fThetaIn(0), 
fTrackLength(0), fThetaOut((Double_t)(std::numeric_limits<double>::infinity())), fBounce(0), fDirIn(G4ThreeVector()), 
fDirOut(G4ThreeVector()), fNgamma(0), fNgammaSec(0), fNCer(0), fRight(0), fLeft(0), 
fDown(0), fUp(0), fBack(0), fFront(0), fSiPM(0), fDecayTime(-1){
	fScintCollection = nullptr;
	collectionName.insert("scintCollection");
hitsCID = -1;

}

ScintSD::ScintSD(G4String name, G4int ntuple):
G4VSensitiveDetector("tmp"), fEvent(-1), fScintNo(-1), fParticleID(0),fEin(0), fEdep(0), fEout(0), fDelta(0), fThetaIn(0), 
fTrackLength(0), fThetaOut((Double_t)(std::numeric_limits<double>::infinity())), fBounce(0), fDirIn(G4ThreeVector()), 
fDirOut(G4ThreeVector()), fNgamma(0), fNgammaSec(0), fNCer(0), fRight(0), fLeft(0), 
fDown(0), fUp(0), fBack(0), fFront(0), fSiPM(0), fDecayTime(-1){
	G4cout<<"Create a tmp ScintSD to fill the Ntupla"<<G4endl;

	G4AnalysisManager *man = G4AnalysisManager::Instance();



	// Ntuple for the Scintillator
	man->CreateNtuple(name, name);
	fEvent_ID = 	man->CreateNtupleIColumn("fEvent");
	fScintNo_ID = 	man->CreateNtupleIColumn("fScintNo");
	fParticleID_ID= man->CreateNtupleIColumn("fParticleID");
	fNgamma_ID = 	man->CreateNtupleIColumn("fNgamma");
	fNgammaSec_ID = man->CreateNtupleIColumn("fNgammaSec");
	fNCer_ID = 		man->CreateNtupleIColumn("fNCer");
	fRight_ID = 	man->CreateNtupleIColumn("fRight");
	fLeft_ID = 		man->CreateNtupleIColumn("fLeft");
	fDown_ID = 		man->CreateNtupleIColumn("fDown");
	fUp_ID = 		man->CreateNtupleIColumn("fUp");
	fBack_ID = 		man->CreateNtupleIColumn("fBack");
	fFront_ID = 	man->CreateNtupleIColumn("fFront");
	fTrackLength_ID = man->CreateNtupleDColumn("fTrackLength");
	fThetaIn_ID = 	man->CreateNtupleDColumn("fThetaIn");
	fThetaOut_ID = 	man->CreateNtupleDColumn("fThetaOut");
	fEin_ID = 		man->CreateNtupleDColumn("fEin");
	fEout_ID = 		man->CreateNtupleDColumn("fEout");
	fEdep_ID = 		man->CreateNtupleDColumn("fEdep");
	fDecayTime_ID = man->CreateNtupleDColumn("fDecayTime");

	fPosInX_ID	= man->CreateNtupleDColumn("fPosInX");
	fPosInY_ID	= man->CreateNtupleDColumn("fPosInY");
	fPosInZ_ID	= man->CreateNtupleDColumn("fPosInZ");
	fMomInX_ID	= man->CreateNtupleDColumn("fMomInX");
	fMomInY_ID	= man->CreateNtupleDColumn("fMomInY");
	fMomInZ_ID	= man->CreateNtupleDColumn("fMomInZ");	
	fTimeIn_ID	= man->CreateNtupleDColumn("fTimeIn");
	fPosOutX_ID	= man->CreateNtupleDColumn("fPosOutX");
	fPosOutY_ID	= man->CreateNtupleDColumn("fPosOutY");
	fPosOutZ_ID	= man->CreateNtupleDColumn("fPosOutZ");
	fMomOutX_ID	= man->CreateNtupleDColumn("fMomOutX");
	fMomOutY_ID	= man->CreateNtupleDColumn("fMomOutY");
	fMomOutZ_ID	= man->CreateNtupleDColumn("fMomOutZ");
	fTimeOut_ID	= man->CreateNtupleDColumn("fTimeOut");
	man->CreateNtupleDColumn("fTest",fTest);

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
int debug	=	3;

	//? Take start and end of the G4Step
	G4StepPoint* preStep = aStep->GetPreStepPoint();
	G4StepPoint* postStep = aStep->GetPostStepPoint();

	//? Take the G4VPhysicalVolume for both start and end
    G4TouchableHistory* thePreTouchable = (G4TouchableHistory*)(preStep->GetTouchable());
    G4VPhysicalVolume* thePrePV = thePreTouchable->GetVolume();
    G4TouchableHistory* thePostTouchable = (G4TouchableHistory*)(postStep->GetTouchable());
    G4VPhysicalVolume* thePostPV = thePostTouchable->GetVolume();

if(debug>4) G4cout<<thePrePV->GetName()<<" "<<thePostPV->GetName()<< G4endl;
if(debug>4) G4cout<<"track id "<< aStep->GetTrack()->GetTrackID()<< G4endl;
	
	//? If it is the particle I generated
	if(aStep->GetTrack()->GetTrackID() == 1){

		//? Debug on StepProcess and TrackStatus
		if (preStep->GetProcessDefinedStep()){
			G4String StepProcessName = preStep->GetProcessDefinedStep()->GetProcessName();
if(debug>5)G4cout<<"StepProcessName " <<StepProcessName<<G4endl;	
		} 
		if(aStep->GetTrack()->GetTrackStatus()==fStopAndKill){
if(debug>5)G4cout<< "fStopAndKill"<<G4endl;
		}

		//? Deposited energy, delta and lenght of the track
		G4double edep = aStep->GetTotalEnergyDeposit();
		G4double delta = aStep->GetPostStepPoint()->GetKineticEnergy() - aStep->GetPreStepPoint()->GetKineticEnergy() + edep;

		fEdep += edep;
		fDelta -= delta;
		fTrackLength += aStep->GetStepLength();
				
		G4double ein = 0, eout = 0;

		//? Counting and classifying photons
		const std::vector<const G4Track*>* secondaries = aStep->GetSecondaryInCurrentStep();
		if(secondaries->size() > 0){
			for(unsigned int i = 0; i < secondaries->size(); i++){
				if(secondaries->at(i)->GetParentID() > 0){
					if(secondaries->at(i)->GetDynamicParticle()->GetParticleDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()){
						fNgamma += 1;
						fNgammaSec += 1;
					}

					//! Better definition for the decay to keep it general purpose
					else if(secondaries->at(i)->GetParticleDefinition()->GetParticleName() == "e+"){
						if(secondaries->at(i)->GetCreatorProcess()->GetProcessName() == "Decay"){
							fDecayTime = aStep->GetTrack()->GetGlobalTime();
						}
					}
				}
			}
		}

		//? Saving incoming primary characteristics	
		if(aStep->IsFirstStepInVolume() && fEin == 0){	

			//? Track, pdgID, event number, which scintillator, energy, momentum
			G4Track * track = aStep->GetTrack();
    		fParticleID = track->GetParticleDefinition()->GetPDGEncoding();
    		fEvent = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
    		fScintNo = thePrePV->GetCopyNo();
			ein = preStep->GetKineticEnergy();
			fEin = ein;
			fMomIn = preStep->GetMomentum();

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
if(debug>0) G4cout<<"fDirIn [pre trasform] : "<<fDirIn.x()<<" "<<fDirIn.y()<<" "<<fDirIn.z()<<G4endl;
if(debug>0) G4cout<<"fDirIn [Volume's reference] :"<<fDirIn_trans.x()<<" "<<fDirIn_trans.y()<<" "<<fDirIn_trans.z()<<G4endl;

			G4ThreeVector norm = -thePreTouchable->GetVolume(0)->GetLogicalVolume()->GetSolid()->SurfaceNormal(localPos);
			fThetaIn = norm.dot(fDirIn_trans);
if(debug>0) G4cout<<"norm: "<<norm.x()<<" "<<norm.y()<<" "<<norm.z()<<G4endl;
if(debug>0) G4cout<<"cos(fThetaIn) = "<<fThetaIn<<" and fThetaIn [deg] = "<<std::acos(fThetaIn) * 180/CLHEP::pi<<G4endl;

			//? Position and time
			fPosIn = aStep->GetPreStepPoint()->GetPosition();
			fTimeIn = aStep->GetPreStepPoint()->GetGlobalTime();

			//? return only if it is NOT also the last step
			if (!aStep->IsLastStepInVolume()) return false;
		}
		
		//? Eout to see if the particle was stopped
		eout = postStep->GetKineticEnergy();
if(debug>3) G4cout<<"eout "<< eout <<G4endl;

		//? Stopped particle (set theta to infinity)
		//! DirOut?
		if (eout == 0.){
if(debug>0) G4cout<<"Particle stopped!"<<G4endl;
			fEout = eout;
			fDirOut = preStep->GetMomentumDirection();
			fThetaOut = (Double_t)(std::numeric_limits<double>::infinity());
			if(fBounce > 0) fBounce += 1;
			fPosOut = aStep->GetPostStepPoint()->GetPosition();
			fTimeOut = aStep->GetPostStepPoint()->GetGlobalTime();
			
			//? Should I kill the track?
			//aStep->Getrack()->SetTrackStatus(fStopAndKill);
			FillHit();
			return true;
		}

		//? Exiting particle
		else if(postStep->GetStepStatus() == fGeomBoundary  && eout != 0){
if(debug>0) G4cout<<"Particle NOT stopped!"<<G4endl;
			fEout = eout;
			fMomOut = postStep->GetMomentum();
			fDirOut = postStep->GetMomentumDirection();

			//? Angle (transform the momentum direction to the volume's reference system)
			G4ThreeVector worldPos = postStep->GetPosition();
			G4ThreeVector localPos = thePreTouchable->GetHistory()->GetTopTransform().TransformPoint(worldPos);
			G4AffineTransform momentumTransform = thePreTouchable->GetHistory()->GetTopTransform();
			momentumTransform.SetNetTranslation(G4ThreeVector(0,0,0));
			fDirOut_trans = momentumTransform.TransformPoint(postStep->GetMomentumDirection());
if(debug>0) G4cout<<"fDirOut [pre trasform] : "<<fDirOut.x()<<" "<<fDirOut.y()<<" "<<fDirOut.z()<<G4endl;
if(debug>0) G4cout<<"fDirOut [Volume's reference] :"<<fDirOut_trans.x()<<" "<<fDirOut_trans.y()<<" "<<fDirOut_trans.z()<<G4endl;

			//? If it is the first step filp the norm got from thePreTouchable 
			G4ThreeVector norm = thePreTouchable->GetVolume(0)->GetLogicalVolume()->GetSolid()->SurfaceNormal(localPos);
			fThetaOut = norm.dot(fDirOut_trans);
			if(fThetaOut<0){
				if(aStep->IsFirstStepInVolume()) norm = - norm;
				fThetaOut = norm.dot(fDirOut_trans);
			}
if(debug>0) G4cout<<"norm: "<<norm.x()<<" "<<norm.y()<<" "<<norm.z()<<G4endl;
if(debug>0) G4cout<<"cos(fThetaOut) = "<<fThetaOut<<" and fThetaOut [deg] = "<<std::acos(fThetaOut) * 180/CLHEP::pi<<G4endl;

			//? Position and time
			fPosOut = aStep->GetPostStepPoint()->GetPosition();
			fTimeOut = aStep->GetPostStepPoint()->GetGlobalTime();
			
			//? Should I kill the track?
			//aStep->GetTrack()->SetTrackStatus(fStopAndKill);
			FillHit();
			return true;
		}
if(debug>0) G4cout<<"End of TrackID = 1"<<G4endl;
		return false;
	}

	//? If it is an G4OpticalPhoton
	//! CHECK REQUIREMENTS ON THE MOMENTUMDIRECTION()
	else if(aStep->GetTrack()->GetDynamicParticle()->GetParticleDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()){
		G4VSolid* solid = thePrePV->GetLogicalVolume()->GetSolid();

		G4Box* boxSolid = (G4Box*)(solid);
		if(postStep->GetStepStatus() == fGeomBoundary){
			G4double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
			G4ThreeVector stppos = postStep->GetPosition();
			G4ThreeVector localpos = thePreTouchable->GetHistory()->GetTopTransform().TransformPoint(stppos);
			G4double dimensionX = boxSolid->GetXHalfLength();
			G4double dimensionY = boxSolid->GetYHalfLength();
			G4double dimensionZ = boxSolid->GetZHalfLength();
			
			if(std::fabs(localpos.x() + dimensionX) < kCarTolerance && postStep->GetMomentumDirection().getX() < 0){
				fLeft += 1;
			}
			else if(std::fabs(localpos.x() - dimensionX) < kCarTolerance && postStep->GetMomentumDirection().getX() > 0){
				fRight += 1;
			}
			else if(std::fabs(localpos.y() + dimensionY) < kCarTolerance && postStep->GetMomentumDirection().getY() < 0){
				fDown += 1;
			}
			else if(std::fabs(localpos.y() - dimensionY) < kCarTolerance && postStep->GetMomentumDirection().getY() > 0){
				fUp += 1;
			}
			else if(std::fabs(localpos.z() + dimensionZ) < kCarTolerance && postStep->GetMomentumDirection().getZ() < 0) {
				fBack += 1; 
			}
			else if(std::fabs(localpos.z() - dimensionZ) < kCarTolerance && postStep->GetMomentumDirection().getZ() > 0){
				fFront += 1;
			}
			aStep->GetTrack()->SetTrackStatus(fStopAndKill);
		}
		return false;
	}

	//? Everything else: ionization, decay etc
	else{
		if(debug>10)G4cout<<"else particle id "<<aStep->GetTrack()->GetDynamicParticle()->GetParticleDefinition()->GetPDGEncoding()<<G4endl;
		if (preStep->GetProcessDefinedStep()){
			G4String StepProcessName = preStep->GetProcessDefinedStep()->GetProcessName();
if(debug>5)G4cout<<"StepProcessName " <<StepProcessName<<G4endl;	
		} 
		//? save if it generated photons
		const std::vector<const G4Track*>* secondaries = aStep->GetSecondaryInCurrentStep();
		if(secondaries->size() > 0){
			for(unsigned int i = 0; i < secondaries->size(); i++){
if(debug>4) G4cout<<secondaries->at(i)->GetDynamicParticle()->GetParticleDefinition()->GetPDGEncoding()<<G4endl;
				if(secondaries->at(i)->GetParentID() > 0){
					if(secondaries->at(i)->GetDynamicParticle()->GetParticleDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()){
						fNgammaSec += 1;
						if(secondaries->at(i)->GetCreatorProcess()->GetProcessName() == "Cerenkov"){
							fNCer += 1;
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
	Hit->SetThetaIn(std::acos(fThetaIn) * 180/CLHEP::pi);
	Hit->SetTrackLength(fTrackLength);
	Hit->SetThetaOut(std::acos(fThetaOut) * 180/CLHEP::pi);
	Hit->SetBounce(fBounce);
	Hit->SetPosIn(fPosIn);
	Hit->SetMomIn(fMomIn);
	Hit->SetTimeIn(fTimeIn);
	Hit->SetPosOut(fPosOut);
	Hit->SetMomOut(fMomOut);
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
	fMomIn = fMomOut = G4ThreeVector();
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

void ScintSD::FillNtupla(G4AnalysisManager *man, ScintHit* scintHit, G4int ntupla){
//G4cout<<"FillScintNtupla "<<ntupla<<G4endl;
	fTest.push_back(scintHit->GetPosIn().getX());
	fTest.push_back(scintHit->GetPosIn().getY());
	fTest.push_back(scintHit->GetPosIn().getZ());
	//! Smarter way to fill it?
	man->FillNtupleIColumn(ntupla, fEvent_ID		, scintHit->GetEvent());
	man->FillNtupleIColumn(ntupla, fScintNo_ID		, scintHit->GetScintNo());
	man->FillNtupleIColumn(ntupla, fParticleID_ID	, scintHit->GetParticleID());
    man->FillNtupleIColumn(ntupla, fNgamma_ID		, scintHit->GetNgamma());
    man->FillNtupleIColumn(ntupla, fNgammaSec_ID	, scintHit->GetNgammaSec());
    man->FillNtupleIColumn(ntupla, fNCer_ID			, scintHit->GetNCer());
    man->FillNtupleIColumn(ntupla, fRight_ID		, scintHit->GetCurrentFront());
    man->FillNtupleIColumn(ntupla, fLeft_ID		, scintHit->GetCurrentLeft());
    man->FillNtupleIColumn(ntupla, fDown_ID		, scintHit->GetCurrentDown());
    man->FillNtupleIColumn(ntupla, fUp_ID		, scintHit->GetCurrentUp());
    man->FillNtupleIColumn(ntupla, fBack_ID		, scintHit->GetCurrentBack());
    man->FillNtupleIColumn(ntupla, fFront_ID		, scintHit->GetCurrentFront());
    man->FillNtupleDColumn(ntupla, fThetaIn_ID		, scintHit->GetThetaIn());
    man->FillNtupleDColumn(ntupla, fThetaOut_ID		, scintHit->GetThetaOut());
    man->FillNtupleDColumn(ntupla, fEin_ID			, scintHit->GetEin());
    man->FillNtupleDColumn(ntupla, fEout_ID			, scintHit->GetEout());
    man->FillNtupleDColumn(ntupla, fEdep_ID			, scintHit->GetEdep());
    man->FillNtupleDColumn(ntupla, fTrackLength_ID	, scintHit->GetTrackLength());
	man->FillNtupleDColumn(ntupla, fDecayTime_ID	, scintHit->GetDecayTime());
    man->FillNtupleDColumn(ntupla, fPosInX_ID, scintHit->GetPosIn().x());
    man->FillNtupleDColumn(ntupla, fPosInY_ID, scintHit->GetPosIn().y());
    man->FillNtupleDColumn(ntupla, fPosInZ_ID, scintHit->GetPosIn().z());
	man->FillNtupleDColumn(ntupla, fMomInX_ID, scintHit->GetMomIn().x());
    man->FillNtupleDColumn(ntupla, fMomInY_ID, scintHit->GetMomIn().y());
    man->FillNtupleDColumn(ntupla, fMomInZ_ID, scintHit->GetMomIn().z());
    man->FillNtupleDColumn(ntupla, fTimeIn_ID, scintHit->GetTimeIn());
    man->FillNtupleDColumn(ntupla, fPosOutX_ID, scintHit->GetPosOut().x());
    man->FillNtupleDColumn(ntupla, fPosOutY_ID, scintHit->GetPosOut().y());
    man->FillNtupleDColumn(ntupla, fPosOutZ_ID, scintHit->GetPosOut().z());
    man->FillNtupleDColumn(ntupla, fMomOutX_ID, scintHit->GetMomOut().x());
    man->FillNtupleDColumn(ntupla, fMomOutY_ID, scintHit->GetMomOut().y());
    man->FillNtupleDColumn(ntupla, fMomOutZ_ID, scintHit->GetMomOut().z());
    man->FillNtupleDColumn(ntupla, fTimeOut_ID, scintHit->GetTimeOut());

	man->AddNtupleRow(ntupla);

	fTest.clear();
}