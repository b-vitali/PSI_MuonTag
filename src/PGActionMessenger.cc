/// \file  PGActionMessenger.cc
/// \brief Implementation of the PGActionMessenger class

#include "PGActionMessenger.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

PGActionMessenger::PGActionMessenger(PrimaryGeneratorAction* action) : G4UImessenger(), fPGAction(action){
	fPrimary = new G4UIdirectory("/Primary/");
	fPrimary->SetGuidance("UI commands to specify primary generator.");

	fCmdOCT = new G4UIcmdWithABool("/Primary/OCT", this);
	fCmdOCT->SetGuidance("Activate primary generator for OCT testing.");
	fCmdOCT->SetParameterName("OCT", false);
	fCmdOCT->AvailableForStates(G4State_Idle);

	fCmdDivergence = new G4UIcmdWithADoubleAndUnit("/Primary/Divergence", this);
	fCmdDivergence->SetGuidance("Set angular divergence.");
	fCmdDivergence->SetUnitCategory("angle");
	fCmdDivergence->SetParameterName("divergence", false);
	fCmdDivergence->AvailableForStates(G4State_Idle);
}

PGActionMessenger::~PGActionMessenger(){
	delete fCmdOCT;
	delete fCmdDivergence;
	delete fPrimary;
}

void PGActionMessenger::SetNewValue(G4UIcommand* command, G4String newValue){
	if(command == fCmdOCT){
		fPGAction->SetOCT(fCmdOCT->GetNewBoolValue(newValue));
	}
	else if(command == fCmdDivergence){
		fPGAction->SetDivergence(fCmdDivergence->GetNewDoubleValue(newValue));
	}
}
