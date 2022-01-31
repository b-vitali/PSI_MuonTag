/// \file  VirtualDetectorSD.hh
/// \brief Definition the SensitiveDetector for the VirtualDetector

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
   	if(aStep->GetTrack()->GetTrackID() == 1){
           
        track->SetTrackStatus(fStopAndKill);
    
        G4StepPoint *preStepPoint  = aStep->GetPreStepPoint();
        G4StepPoint *postStepPoint = aStep->GetPostStepPoint();
    
        G4TouchableHistory* thePreTouchable = (G4TouchableHistory*)(preStepPoint->GetTouchable());
    	G4VPhysicalVolume* thePrePV = thePreTouchable->GetVolume();
    	G4TouchableHistory* thePostTouchable = (G4TouchableHistory*)(postStepPoint->GetTouchable());
    	G4VPhysicalVolume* thePostPV = thePostTouchable->GetVolume();
    
        G4cout << thePrePV->GetName() << " " <<preStepPoint->GetPosition()<< "&&" << thePostPV->GetName()<<postStepPoint->GetPosition() << G4endl;
    
    
        if(thePrePV->GetName() == "World" && thePostPV->GetName() == "VD")
        {
            G4ThreeVector posIn = preStepPoint->GetPosition();
            G4cout << "position in :" << posIn << G4endl;
        }
    
        else return false;  
    }
}