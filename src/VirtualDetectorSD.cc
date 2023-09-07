/// \file  VirtualDetectorSD.hh
/// \brief Implementation of the SensitiveDetector for the VirtualDetector

// Basic class
#include "VirtualDetectorSD.hh"

VirtualDetectorSD::VirtualDetectorSD(G4String name):
G4VSensitiveDetector(name)
{
        G4cout << "VD created"<< G4endl;
}

VirtualDetectorSD::~VirtualDetectorSD()
{}


G4bool VirtualDetectorSD::ProcessHits(G4Step * aStep, G4TouchableHistory * ROHist)
{
    G4Track * track = aStep->GetTrack();
    fParticleID = track->GetParticleDefinition()->GetPDGEncoding();

    // Reject the OpticalPhotons and "not the first step" form the VD
    if(!aStep->IsFirstStepInVolume() || fParticleID == -22) return true;
   	//if(aStep->GetTrack()->GetTrackID() == 1)
    //{ 
        //track->SetTrackStatus(fStopAndKill);
    
        G4StepPoint *preStepPoint  = aStep->GetPreStepPoint();
        G4StepPoint *postStepPoint = aStep->GetPostStepPoint();
    
        G4TouchableHistory* thePreTouchable = (G4TouchableHistory*)(preStepPoint->GetTouchable());
    	G4VPhysicalVolume* thePrePV = thePreTouchable->GetVolume();
    	G4TouchableHistory* thePostTouchable = (G4TouchableHistory*)(postStepPoint->GetTouchable());
    	G4VPhysicalVolume* thePostPV = thePostTouchable->GetVolume();
    
        fVDNo = thePrePV->GetCopyNo();

    //    G4cout << thePrePV->GetName() << " " <<preStepPoint->GetPosition()<< " && " << thePostPV->GetName()<< " " <<postStepPoint->GetPosition() << G4endl;
    //    G4cout << thePrePV->GetName() << " " <<preStepPoint->GetMomentum()<< " && " << thePostPV->GetName()<< " " <<postStepPoint->GetMomentum() << G4endl;
    //    G4cout << thePrePV->GetName() << " " <<preStepPoint->GetLocalTime()<< " && " << thePostPV->GetName()<< " " <<postStepPoint->GetLocalTime() << G4endl;
    
        fVDTime = preStepPoint->GetLocalTime();
        G4int evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

        G4AnalysisManager *man = G4AnalysisManager::Instance();
        man->FillNtupleIColumn(0, 0, evt);
        man->FillNtupleIColumn(0, 1, fVDNo);
        man->FillNtupleIColumn(0, 2, fParticleID);
        man->FillNtupleDColumn(0, 3, fVDTime);
        man->FillNtupleDColumn(0, 4, preStepPoint->GetMomentum().mag());
        man->FillNtupleDColumn(0, 5, preStepPoint->GetPosition().x());
        man->FillNtupleDColumn(0, 6, preStepPoint->GetPosition().y());
        man->FillNtupleDColumn(0, 7, preStepPoint->GetPosition().z());
        man->FillNtupleDColumn(0, 8, preStepPoint->GetMomentum().x());
        man->FillNtupleDColumn(0, 9, preStepPoint->GetMomentum().y());
        man->FillNtupleDColumn(0, 10, preStepPoint->GetMomentum().z());
        man->AddNtupleRow(0);
    //}

    //else return false;  

}