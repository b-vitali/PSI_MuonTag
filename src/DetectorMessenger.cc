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
	fMuEDMDirectory = new G4UIdirectory("/MuEDM/");
	fMuEDMDirectory->SetGuidance("UI commands specific to this simulation.");
	
	fScintDirectory = new G4UIdirectory("/MuEDM/Scint/");
	fScintDirectory->SetGuidance("Scintillator construction control");
	
	fScintSizeCmd = new G4UIcmdWithADoubleAndUnit("/MuEDM/Scint/ScintSize", this);
	fScintSizeCmd->SetGuidance("Define BC400 Scintillator size");
	fScintSizeCmd->SetParameterName("ScintSize", false);
	fScintSizeCmd->SetUnitCategory("Length");
	fScintSizeCmd->AvailableForStates(G4State_Idle);
	
	fScintSizeCmd3 = new G4UIcmdWith3VectorAndUnit("/MuEDM/Scint/ScintSize3", this);
	fScintSizeCmd3->SetGuidance("Define BC400 Scintillator size");
	fScintSizeCmd3->SetParameterName("ScintSizeX", "ScintSizeY", "ScintSizeZ", false);
	fScintSizeCmd3->SetUnitCategory("Length");
	fScintSizeCmd3->AvailableForStates(G4State_Idle);

	fScintMaterialCmd = new G4UIcmdWithAString("/MuEDM/Scint/ScintMaterial", this);
	fScintMaterialCmd->SetGuidance("Set scintillating material");
	fScintMaterialCmd->SetParameterName("material", false);
	fScintMaterialCmd->SetCandidates("BC400 || LYSO");
	fScintMaterialCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

DetectorMessenger::~DetectorMessenger(){
	delete fScintSizeCmd;
	delete fScintSizeCmd3;
	delete fScintMaterialCmd;
	delete fMuEDMDirectory;
	delete fScintDirectory;
}

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue){
/*
	if(command == fScintSizeCmd){
		if(fScintSizeCmd3->GetNew3VectorValue(newValue)){
			fDetectorConstruction->SetScinttalSize(fScintSizeCmd3->GetNew3VectorValue(newValue));
		}
		else{
			fDetectorConstruction->SetScinttalSize(fScintSizeCmd->GetNewDoubleValue(newValue));
		}
	}
*/

	if(command == fScintSizeCmd){
		fDetectorConstruction->SetScintSize(fScintSizeCmd->GetNewDoubleValue(newValue));
	}
	
	else if (command == fScintSizeCmd3){
		fDetectorConstruction->SetScintSize(fScintSizeCmd3->GetNew3VectorValue(newValue));
	}

	else if(command == fScintMaterialCmd){
		fDetectorConstruction->SetScintMaterial(newValue);
	}
}



