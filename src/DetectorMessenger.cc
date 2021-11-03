/// \file DetectorMessenger.cc
/// \brief Implementation of the DetectorMessenger class

#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithAString.hh"


DetectorMessenger::DetectorMessenger(DetectorConstruction* Det) : G4UImessenger(), fDetectorConstruction(Det){
	fElementDirectory = new G4UIdirectory("/Element/");
	fElementDirectory->SetGuidance("UI commands specific to this simulation.");
	
	fDetDirectory = new G4UIdirectory("Element/det/");
	fDetDirectory->SetGuidance("Detector construction control");
	
	fCrysSizeCmd = new G4UIcmdWithADoubleAndUnit("/Element/det/CrysSize", this);
	fCrysSizeCmd->SetGuidance("Define BC400 crystal size");
	fCrysSizeCmd->SetParameterName("CrysSize", false);
	fCrysSizeCmd->SetUnitCategory("Length");
	fCrysSizeCmd->AvailableForStates(G4State_Idle);
	
	fCrysSizeCmd3 = new G4UIcmdWith3VectorAndUnit("/Element/det/CrysSize3", this);
	fCrysSizeCmd3->SetGuidance("Define BC400 crystal size");
	fCrysSizeCmd3->SetParameterName("CrysSizeX", "CrysSizeY", "CrysSizeZ", false);
	fCrysSizeCmd3->SetUnitCategory("Length");
	fCrysSizeCmd3->AvailableForStates(G4State_Idle);
	
	fGround = new G4UIcmdWithADouble("/Element/det/Ground", this);
	fGround->SetGuidance("Set ground value in glisur value for scintillator");
	fGround->SetParameterName("Ground", false);
	fGround->AvailableForStates(G4State_Idle);
	
	fAngle = new G4UIcmdWithADoubleAndUnit("/Element/det/Angle", this);
	fAngle->SetGuidance("Set angle with SiPm window");
	fAngle->SetParameterName("Angle", false);
	fAngle->SetUnitCategory("Angle");
	fAngle->AvailableForStates(G4State_Idle);
	
	fAngleWithOpticalGrease = new G4UIcmdWithADoubleAndUnit("/Element/det/AngleWOG", this);
	fAngleWithOpticalGrease->SetGuidance("Set angle with SiPm window filling gap with optical grease");
	fAngleWithOpticalGrease->SetParameterName("AngleWOG", false);
	fAngleWithOpticalGrease->SetUnitCategory("Angle");
	fAngleWithOpticalGrease->AvailableForStates(G4State_Idle);
	
	fTilt = new G4UIcmdWithADoubleAndUnit("/Element/det/Tilt", this);
	fTilt->SetGuidance("Set scintillator's tilt");
	fTilt->SetParameterName("Tilt", false);
	fTilt->SetUnitCategory("Length");
	fTilt->AvailableForStates(G4State_Idle);
	
	fSiPMmodelCmd = new G4UIcmdWithAString("/Element/det/SiPMmodel", this);
	fSiPMmodelCmd->SetGuidance("Set SiPM model");
	fSiPMmodelCmd->SetParameterName("model", false);
	fSiPMmodelCmd->SetCandidates("75PE || 50PE || 25PE || 75CS || 50CS || 25CS");
	fSiPMmodelCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

	fCrysMaterialCmd = new G4UIcmdWithAString("/Element/det/CrysMaterial", this);
	fCrysMaterialCmd->SetGuidance("Set scintillating material");
	fCrysMaterialCmd->SetParameterName("material", false);
	fCrysMaterialCmd->SetCandidates("BC400 || LYSO");
	fCrysMaterialCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

}

DetectorMessenger::~DetectorMessenger(){
	delete fCrysSizeCmd;
	delete fCrysSizeCmd3;
	delete fGround;
	delete fAngle;
	delete fAngleWithOpticalGrease;
	delete fTilt;
	delete fSiPMmodelCmd;
	delete fCrysMaterialCmd;
	delete fDetDirectory;
	delete fElementDirectory;
}

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue){
/*
	if(command == fCrysSizeCmd){
		if(fCrysSizeCmd3->GetNew3VectorValue(newValue)){
			fDetectorConstruction->SetCrystalSize(fCrysSizeCmd3->GetNew3VectorValue(newValue));
		}
		else{
			fDetectorConstruction->SetCrystalSize(fCrysSizeCmd->GetNewDoubleValue(newValue));
		}
	}
*/

	if(command == fCrysSizeCmd){
		fDetectorConstruction->SetCrystalSize(fCrysSizeCmd->GetNewDoubleValue(newValue));
	}
	
	else if (command == fCrysSizeCmd3){
		fDetectorConstruction->SetCrystalSize(fCrysSizeCmd3->GetNew3VectorValue(newValue));
	}
	
	else if (command == fGround){
		fDetectorConstruction->SetGround(fGround->GetNewDoubleValue(newValue));
	}
	
	else if (command == fAngle){
		fDetectorConstruction->SetAngle(fAngle->GetNewDoubleValue(newValue));
	}
	
	else if (command == fAngleWithOpticalGrease){
		fDetectorConstruction->SetAngleWithOpticalGrease(fAngleWithOpticalGrease->GetNewDoubleValue(newValue));
	}
	
	else if (command == fTilt){
		fDetectorConstruction->SetTilt(fTilt->GetNewDoubleValue(newValue));
	}

	else if(command == fSiPMmodelCmd){
		fDetectorConstruction->SetSiPMmodel(newValue);
	}

	else if(command == fCrysMaterialCmd){
		fDetectorConstruction->SetCrystalMaterial(newValue);
	}
}



