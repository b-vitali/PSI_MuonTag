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
	ScintName = name;
	fScintCollection = nullptr;
	collectionName.insert("scintCollection");
	hitsCID = -1;
	G4cout<<"create ScintSD: "<<name<<G4endl;

	G4AnalysisManager *man = G4AnalysisManager::Instance();

	double x = 20. + 2; 	x= x/2.;
	double y = 20. + 2;		y= y/2.;
	double z = 2.5 + 0.5;	z= z/2.;
	PhotonTime_ID			= man->CreateH2(name+"_PhotonTime","PhotonTime", 200, 0., 200, 5, 0,5);
	PhotonPositionFront_ID	= man->CreateH3(name+"_PhotonPos_Front","PhotonPositionFront",220, -11., 11, 220, -11., 11, 5, 0,5);
	PhotonPositionBack_ID	= man->CreateH3(name+"_PhotonPos_Back","PhotonPositionBack",220, -11., 11, 220, -11., 11, 5, 0,5);
	PhotonPositionLeft_ID	= man->CreateH3(name+"_PhotonPos_Left","PhotonPositionLeft",120, -3., 3., 220, -11., 11, 5, 0,5);
	PhotonPositionRight_ID	= man->CreateH3(name+"_PhotonPos_Right","PhotonPositionRight",120, -3., 3., 220, -11., 11, 5, 0,5);
	PhotonPositionUp_ID		= man->CreateH3(name+"_PhotonPos_Up","PhotonPositionUp",220, -11., 11, 120, -3., 3., 5, 0,5);
	PhotonPositionDown_ID	= man->CreateH3(name+"_PhotonPos_Down","PhotonPositionDown",220, -11., 11, 120, -3., 3., 5, 0,5);
	
}

ScintSD::ScintSD(G4String name, G4int ntuple):
G4VSensitiveDetector(name){
	G4cout<<"Create a tmp ScintSD to fill the Ntupla"<<G4endl;

	G4AnalysisManager *man = G4AnalysisManager::Instance();

	// Ntuple for the Scintillator
	man->CreateNtuple(name, name);

	man->CreateNtupleIColumn("fEvent",fEvent);
	man->CreateNtupleIColumn("fScintNo",fScintNo);
	man->CreateNtupleIColumn("fParticleID",fParticleID);
	man->CreateNtupleIColumn("fNgamma",fNgamma);
	man->CreateNtupleIColumn("fNgammaSec",fNgammaSec);
	man->CreateNtupleIColumn("fNCer",fNCer);
	man->CreateNtupleIColumn("fAbsorption",fAbsorption);	
	man->CreateNtupleIColumn("fReflection",fReflection);

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

	// aid variables just to check
	Trk=-5;
	TrkParent=-5;
	TrackOneIn = false;
	EntryCreated = false;

	/*
		Debug feature:
		use a G4String and "contains"
		no printout == ""; p = pre-info, i = in, o = out, g = gamma,
		e = else, 1 = trk ==1, + = additionla info; 
	*/
	// G4String debug	= "p1 i o n g e";
	debug	= "";
}


void ScintSD::CreateEntry(G4Step *aStep){
	G4cout<<"CreateEntry trk : "<<aStep->GetTrack()->GetTrackID()<<" -> "<<Trk<<G4endl;
	fEdep.push_back(0);
	fDelta.push_back(0);
	fTrackLength.push_back(0);
	fNgamma.push_back(0);
	fNgammaSec.push_back(0);
	fNCer.push_back(0);
	fAbsorption.push_back(0);
	fReflection.push_back(0);

	fRight.push_back(0);
	fLeft.push_back(0);
	fDown.push_back(0);
	fUp.push_back(0);
	fBack.push_back(0);
	fFront.push_back(0);

	G4Track * track = aStep->GetTrack();
	G4StepPoint* preStep = aStep->GetPreStepPoint();
    G4TouchableHistory* thePreTouchable = (G4TouchableHistory*)(preStep->GetTouchable());

	fEvent.push_back(G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID());
	fParticleID.push_back(track->GetParticleDefinition()->GetPDGEncoding());

	//! GetCopyNumber 0 or 1 if there is an "element volume" with "read" and "scint in"
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
	G4ThreeVector worldPos = preStep->GetPosition();
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

	EntryCreated = true;
	EntryTrk = aStep->GetTrack()->GetTrackID();
}



G4bool ScintSD::ProcessHits(G4Step *aStep, G4TouchableHistory* ROhist){	

	//if(debug.contains("p"))	G4cout<<"Ev : "<<G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID()<<G4endl;

	if(aStep->GetStepLength() == 0 && aStep->GetTotalEnergyDeposit() == 0) {
		// G4cout<<"step of lenght 0 and edep 0"<<G4endl; 
		return false;
	}

	//? Take start and end of the G4Step
	G4StepPoint* preStep = aStep->GetPreStepPoint();
	G4StepPoint* postStep = aStep->GetPostStepPoint();

	//? Take the G4VPhysicalVolume for both start and end
    G4TouchableHistory* thePreTouchable = (G4TouchableHistory*)(preStep->GetTouchable());
    G4VPhysicalVolume* thePrePV = thePreTouchable->GetVolume();
    G4TouchableHistory* thePostTouchable = (G4TouchableHistory*)(postStep->GetTouchable());
    G4VPhysicalVolume* thePostPV = thePostTouchable->GetVolume();

	if(debug.contains("p+")) G4cout<<thePrePV->GetName()<<" "<<thePostPV->GetName()<< G4endl;
	if(debug.contains("p+")) G4cout<<"track id "<< aStep->GetTrack()->GetTrackID()<< G4endl;


	if(aStep->GetTrack()->GetTrackID() != Trk) {
		Trk = aStep->GetTrack()->GetTrackID();
		TrkParent = aStep->GetTrack()->GetParentID();
       	if(debug.contains("p+")) G4cout<<"new Trk : "<< Trk <<" child of :"<<TrkParent<<G4endl;

		//? If the current hit is in another scintillator reset the EntryCreated
		//? (first check that fScintNo is not empty)
		if(EntryCreated){
			//G4cout<<"the Hit stored scintNo : "<< fScintNo.at(fScintNo.size()-1) << " and we are in "<<preStep->GetTouchable()->GetCopyNumber(1)+1<<G4endl;

			if(fScintNo.at(fScintNo.size()-1) != preStep->GetTouchable()->GetCopyNumber(1)+1){
				if(debug.contains("p+")) G4cout<<"new Trk : "<< Trk <<" child of :"<<TrkParent<<G4endl;
				if(debug.contains("p+")) G4cout<< "fScintNo changed, EntryCreated = false"<<G4endl;
				EntryCreated = false;
			}

			if(Trk==1 && !TrackOneIn){
				if(debug.contains("p+")) G4cout<< "Same fScintNo but the TrackOneIn == false"<<G4endl;
				EntryCreated = false;
			}

		}
	}

	//? ------------------------------ ?//
	//? First step in the scintillator ?//
	//? ------------------------------ ?//
	if(!EntryCreated){

		int pdgid = aStep->GetTrack()->GetDynamicParticle()->GetParticleDefinition()->GetPDGEncoding();
		if(debug.contains("p+")) G4cout<<"pdgid : "<<pdgid<<G4endl;
		
		//? A photon is the first entry? 
		if( pdgid == 22 || pdgid == -22 ){ //photon pdgid == 22 or -22
			if(debug.contains("p+")) G4cout<<"Photon with no entry"<<G4endl;
			aStep->GetTrack()->SetTrackStatus(fStopAndKill);
			return false;
		}

		//? A (anti)neutrino is the first entry? 
		else if( pdgid == 12 || pdgid == 14 || pdgid == 16 || pdgid == -12 || pdgid == -14 || pdgid == -16 ){
			if(debug.contains("p+")) G4cout<<"Neutrino with no entry"<<G4endl;
			aStep->GetTrack()->SetTrackStatus(fStopAndKill);
			return false;
		}

		//? A NON photon is the first entry? 
		else if(aStep->GetTrack()->GetDynamicParticle()->GetParticleDefinition() != G4OpticalPhoton::OpticalPhotonDefinition()){ //aStep->GetTrack()->GetParentID() == 1 &&  || (!EntryCreated)
			
			//? double check that it's the first step in volume
			if(aStep->IsFirstStepInVolume()){
				if(debug.contains("p+")) G4cout<<"First step"<<G4endl;
				CreateEntry(aStep);
				if(aStep->GetTrack()->GetTrackID() == 1) TrackOneIn = true;	
			}
			else {G4cout<<"NOT First step, how?"<<G4endl;}
		}
	}

	//? A step from the first particle entering?
	if(aStep->GetTrack()->GetTrackID() == EntryTrk){
		
		//? Counting and classifying photons
		const std::vector<const G4Track*>* secondaries = aStep->GetSecondaryInCurrentStep();
		if(secondaries->size() > 0){
			for(unsigned int i = 0; i < secondaries->size(); i++){
				if(secondaries->at(i)->GetParentID() > 0){
					if(secondaries->at(i)->GetDynamicParticle()->GetParticleDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()){
						fNgamma.at(fNgamma.size()-1)  += 1;
						if(secondaries->at(i)->GetCreatorProcess()->GetProcessName() == "Cerenkov"){
							fNCer.at(fNCer.size()-1)  += 1;
						}	
					}
					
					//! Better definition for the decay to keep it general purpose
					//? If a secondary is a e+ from "Decay" and the hit was created by trk == 1
					else if(secondaries->at(i)->GetParticleDefinition()->GetParticleName() == "e+" && TrackOneIn){
						if(secondaries->at(i)->GetCreatorProcess()->GetProcessName() == "Decay"){
							fDecayTime.push_back(aStep->GetTrack()->GetGlobalTime());
							TrkDecay = secondaries->at(i)->GetTrackID();
						}
					}
				}
			}
		}	
		const G4VProcess* pds = postStep->GetProcessDefinedStep();
    	G4String procname     = pds->GetProcessName();
		if(debug.contains("p "))    if(procname.compare("Transportation") != 0) {	G4cout<<procname<<G4endl;}
		// if(fAbsorption.size() == 0) G4cout<<"STOP fAbsorption"<<G4endl;
		//    if(procname.compare("OpAbsorption") == 0) {	fAbsorption.at(fAbsorption.size()-1) += 1; }
		//? Deposited energy, delta and lenght of the track
		G4double edep = aStep->GetTotalEnergyDeposit();
		G4double delta = aStep->GetPostStepPoint()->GetKineticEnergy() - aStep->GetPreStepPoint()->GetKineticEnergy() + edep;
		fEdep.at(fEdep.size()-1) += edep;
		fDelta.at(fDelta.size()-1) -= delta;
		fTrackLength.at(fTrackLength.size()-1) += aStep->GetStepLength();
		//? Eout to see if the particle was stopped
		G4double eout = 0;
		eout = postStep->GetKineticEnergy();
		if(debug.contains("p")) G4cout<<"eout "<< eout <<G4endl;
		//? Stopped particle (set theta to infinity)
		if (eout == 0.){
		
			if(debug.contains("o")) G4cout<<"Particle stopped!"<<G4endl;
			fEout.push_back(eout);
			fThetaOut.push_back((Double_t)(std::numeric_limits<double>::infinity()));
			fMomOutX.push_back(0);
			fMomOutY.push_back(0);
			fMomOutZ.push_back(0);
			fPosOutX.push_back(aStep->GetPostStepPoint()->GetPosition().getX());
			fPosOutY.push_back(aStep->GetPostStepPoint()->GetPosition().getY());
			fPosOutZ.push_back(aStep->GetPostStepPoint()->GetPosition().getZ());
			fTimeOut.push_back(aStep->GetPostStepPoint()->GetGlobalTime());
			//! WHEN DO I RESET  EntryCreated ???
			TrackOneIn = false;
			//EntryCreated = false;
			return true;
		}
		//? Exiting particle
		else if(postStep->GetStepStatus() == fGeomBoundary  && eout != 0){
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
			if(fThetaOut.at(fThetaOut.size()-1)<0){
				if(aStep->IsFirstStepInVolume()) norm = - norm;
				fThetaOut.at(fThetaOut.size()-1) = norm.dot(fDirOut_trans);
			}
			if(debug.contains("o+")) G4cout<<"norm: "<<norm.x()<<" "<<norm.y()<<" "<<norm.z()<<G4endl;
			if(debug.contains("o")) G4cout<<"cos(fThetaOut) = "<<norm.dot(fDirOut_trans)<<" and fThetaOut [deg] = "<<std::acos(norm.dot(fDirOut_trans)) * 180/CLHEP::pi<<G4endl;
			//! WHEN DO I RESET  EntryCreated ???
			TrackOneIn = false;
			//EntryCreated = false;
			return true;
		}
		if(debug.contains("p"))G4cout<<"fNgamma " <<fNgamma.at(fNgamma.size()-1)<<G4endl;	
		//return false;
	}

	//? If the entry was created and it is an G4OpticalPhoton 
	else if(aStep->GetTrack()->GetDynamicParticle()->GetParticleDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()){
		G4VSolid* solid = thePrePV->GetLogicalVolume()->GetSolid();
		// if(postStep->GetPhysicalVolume()->GetName() != "Element") {return false;}
		G4Box* boxSolid = (G4Box*)(solid);
		
		//?if it is at the border
		if(postStep->GetStepStatus() == fGeomBoundary){
			G4double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
			G4ThreeVector worldPos = postStep->GetPosition();
			G4ThreeVector localpos = thePreTouchable->GetHistory()->GetTopTransform().TransformPoint(worldPos);
			G4double dimensionX = boxSolid->GetXHalfLength();
			G4double dimensionY = boxSolid->GetYHalfLength();
			G4double dimensionZ = boxSolid->GetZHalfLength();
			//?transform the direction to match the rotation of the faces
			G4AffineTransform momentumTransform = thePreTouchable->GetHistory()->GetTopTransform();
			momentumTransform.SetNetTranslation(G4ThreeVector(0,0,0));
			fDirOut_trans = momentumTransform.TransformPoint(postStep->GetMomentumDirection());
			if(debug.contains("g+"))G4cout<<"track id "<< aStep->GetTrack()->GetTrackID()<< G4endl;		
			if(debug.contains("g+"))G4cout<<"pos: "<<localpos.x()<<" "<<localpos.y()<<" "<<localpos.z()<<G4endl;
			if(debug.contains("g+"))G4cout<<"dir: "<<fDirOut_trans.x()<<" "<<fDirOut_trans.y()<<" "<<fDirOut_trans.z()<<G4endl;
			if(debug.contains("g+"))G4cout<<"dir: "<<aStep->GetPostStepPoint()->GetKineticEnergy()<<G4endl;
			//! GetCopyNumber 0 or 1 if there is an "element volume" with "read" and "scint in"
			fScintNo.at(fScintNo.size()-1) = preStep->GetTouchable()->GetCopyNumber(1)+1;
			G4AnalysisManager *man = G4AnalysisManager::Instance();
			// this checks if the photon GOES OUT from the scint
			if(std::fabs(localpos.x() + dimensionX) < kCarTolerance && fDirOut_trans.getX() < 0){
				if(debug.contains("g"))G4cout<<"Left"<<G4endl;
				fLeft.at(fLeft.size()-1)  += 1;
				man->FillH2(PhotonTime_ID, postStep->GetGlobalTime(), fScintNo.at(fScintNo.size()-1));
				man->FillH3(PhotonPositionLeft_ID, localpos.z(), localpos.y(), fScintNo.at(fScintNo.size()-1));
				// aStep->GetTrack()->SetTrackStatus(fStopAndKill);
			}
			else if(std::fabs(localpos.x() - dimensionX) < kCarTolerance && fDirOut_trans.getX() > 0){
				if(debug.contains("g"))G4cout<<"Right"<<G4endl;
				fRight.at(fRight.size()-1)  += 1;
				man->FillH2(PhotonTime_ID, postStep->GetGlobalTime(), fScintNo.at(fScintNo.size()-1));
				man->FillH3(PhotonPositionRight_ID, localpos.z(), localpos.y(), fScintNo.at(fScintNo.size()-1));
				//aStep->GetTrack()->SetTrackStatus(fStopAndKill);
			}
			else if(std::fabs(localpos.y() + dimensionY) < kCarTolerance && fDirOut_trans.getY() < 0){
				if(debug.contains("g"))G4cout<<"Down"<<G4endl;
				fDown.at(fDown.size()-1)  += 1;
				man->FillH2(PhotonTime_ID, postStep->GetGlobalTime(), fScintNo.at(fScintNo.size()-1));
				man->FillH3(PhotonPositionDown_ID, localpos.x(), localpos.z(), fScintNo.at(fScintNo.size()-1));
				// aStep->GetTrack()->SetTrackStatus(fStopAndKill);
			}
			else if(std::fabs(localpos.y() - dimensionY) < kCarTolerance && fDirOut_trans.getY() > 0){
				if(debug.contains("g"))G4cout<<"Up"<<G4endl;
				fUp.at(fUp.size()-1)  += 1;
				man->FillH2(PhotonTime_ID, postStep->GetGlobalTime(), fScintNo.at(fScintNo.size()-1));
				man->FillH3(PhotonPositionUp_ID, localpos.x(), localpos.z(), fScintNo.at(fScintNo.size()-1));
				// aStep->GetTrack()->SetTrackStatus(fStopAndKill);
			}
			else if(std::fabs(localpos.z() + dimensionZ) < kCarTolerance && fDirOut_trans.getZ() < 0) {
				if(debug.contains("g"))G4cout<<"Back"<<G4endl;
				fBack.at(fBack.size()-1)  += 1;
				man->FillH2(PhotonTime_ID, postStep->GetGlobalTime(), fScintNo.at(fScintNo.size()-1));
				man->FillH3(PhotonPositionBack_ID, localpos.x(),localpos.y(), fScintNo.at(fScintNo.size()-1));
				// aStep->GetTrack()->SetTrackStatus(fStopAndKill);
			}
			else if(std::fabs(localpos.z() - dimensionZ) < kCarTolerance && fDirOut_trans.getZ() > 0){
				if(debug.contains("g"))G4cout<<"Front"<<G4endl;
				fFront.at(fFront.size()-1)  += 1;
				man->FillH2(PhotonTime_ID, postStep->GetGlobalTime(), fScintNo.at(fScintNo.size()-1));
				man->FillH3(PhotonPositionFront_ID, localpos.x(),localpos.y(), fScintNo.at(fScintNo.size()-1));
				// aStep->GetTrack()->SetTrackStatus(fStopAndKill);
			}
			else {
				fReflection.at(fReflection.size()-1) += 1; 
				if(debug.contains("g") && fReflection.at(fReflection.size()-1)%1000 == 0)G4cout<<"1k Reflect"<<G4endl;
			}
		return false;
		}
	}
	//! Should I count the energy deposit of these secondaries??
	//?	Everything else: ionization, decay etc
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
						fNgamma.at(fNgamma.size()-1)  += 1;
						fNgammaSec.at(fNgammaSec.size()-1)  += 1;
						if(secondaries->at(i)->GetCreatorProcess()->GetProcessName() == "Cerenkov"){
							fNCer.at(fNCer.size()-1)  += 1;
						}
					}
				}
			}
		}
		return true;
	}
	
}

void ScintSD::FillHit(){
	ScintHit* Hit = new ScintHit();

	Hit->SetEvent(fEvent);
	
	Hit->SetScintNo(fScintNo);
	Hit->SetParticleID(fParticleID);
	Hit->SetNgamma(fNgamma);
	Hit->SetNgammaSec(fNgammaSec);
	Hit->SetNCer(fNCer);
	Hit->SetNAbsorption(fAbsorption);
	Hit->SetNReflection(fReflection);

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

	fEvent.clear();
	fScintNo.clear();
	fParticleID.clear();
	fNgamma.clear();
	fNgammaSec.clear();
	fNCer.clear();
	fAbsorption.clear();
	fReflection.clear();

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
	G4cout<<G4endl<<"End of event n: "<<G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID()<<G4endl;
	if(fAbsorption.size()>0){
		G4cout<<"Number of absorption "<<fAbsorption.at(fAbsorption.size()-1)<<G4endl;
	}
	if(fReflection.size()>0){
		G4cout<<"Number of reflection "<<fReflection.at(fReflection.size()-1)<<G4endl;
	}
	if(fEvent.size()>0){
		G4cout<<"Number of sub-hits "<<fEvent.size()<<G4endl;
		FillHit();
	}
	
	TrackOneIn = false;
	EntryCreated = false;
}

void ScintSD::clear(){}

void ScintSD::DrawAll(){}

void ScintSD::PrintAll(){}

void ScintSD::FillNtupla(G4AnalysisManager *man, ScintHit* scintHit, G4int ntupla){

	G4cout<<"Fill Ntupla: "<<ntupla<<G4endl;
	
	fEvent 	=  scintHit->GetEvent();
	fScintNo 	=  scintHit->GetScintNo();
	fParticleID 	=  scintHit->GetParticleID();
	fNgamma 	=  scintHit->GetNgamma();
	fNgammaSec 	=  scintHit->GetNgammaSec();
	fNCer 	=  scintHit->GetNCer();
	fAbsorption 	=  scintHit->GetNAbsorption();
	fReflection 	=  scintHit->GetNReflection();
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

	man->AddNtupleRow(ntupla);
}