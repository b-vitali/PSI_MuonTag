/// \file  RunActionMessenger.cc
/// \brief Implemenatation of the RunActionMessenger class

#include "RunActionMessenger.hh"
#include "RunAction.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

RunActionMessenger::RunActionMessenger(RunAction* action) : G4UImessenger(), fRunAction(action){
	fAnalysisDirectory = new G4UIdirectory("/Analysis/");
	fAnalysisDirectory->SetGuidance("UI commands to specify analysis options.");

	fCmdFileName = new G4UIcmdWithAString("/Analysis/SetFileName", this);
	fCmdFileName->SetGuidance("Choose analysis file name.");
	fCmdFileName->SetParameterName("Name", false);
	fCmdFileName->AvailableForStates(G4State_Idle);

	fCmdPhotons = new G4UIcmdWithAnInteger("/Analysis/Photons", this);
	fCmdPhotons->SetGuidance("Choose to save photons information or not.");
	fCmdPhotons->SetParameterName("photons", false);
	fCmdPhotons->SetRange("photons >= 0 && photons <= 1");
	fCmdPhotons->SetDefaultValue(0);
	fCmdPhotons->AvailableForStates(G4State_Idle);

	fCmdTracks = new G4UIcmdWithAnInteger("/Analysis/TracksPoints", this);
	fCmdTracks->SetGuidance("Choose the number of points saved for e+ tracks. If 0 they will be completly saved. If 1 they won't be saved, just the starting and the ending point. You can choose to save from 2 to 10 points.");
	fCmdTracks->SetParameterName("tracks", false);
	fCmdTracks->SetRange("tracks >= 0 && tracks <= 10");
	fCmdTracks->AvailableForStates(G4State_Idle);

	fCmdOCT = new G4UIcmdWithABool("/Element/det/OCT", this);
	fCmdOCT->SetGuidance("Activate optical cross talk among pixels in SiPMs");
	fCmdOCT->SetParameterName("OCT", false);
	fCmdOCT->AvailableForStates(G4State_Idle);

	fCmdDN = new G4UIcmdWithABool("/Element/det/DN", this);
	fCmdDN->SetGuidance("Activate dark noise in SiPMs");
	fCmdDN->SetParameterName("DN", false);
	fCmdDN->AvailableForStates(G4State_Idle);

	fCmdGunTime = new G4UIcmdWithADoubleAndUnit("/Primary/Rate", this);
	fCmdGunTime->SetGuidance("Set beam rate.");
	fCmdGunTime->SetParameterName("rate", false);
	fCmdGunTime->SetUnitCategory("Rate");
	fCmdGunTime->SetDefaultValue(1.9e9*hertz);
	fCmdGunTime->AvailableForStates(G4State_Idle);
}

RunActionMessenger::~RunActionMessenger(){
	delete fCmdGunTime;
	delete fCmdFileName;
	delete fCmdOCT;
	delete fCmdDN;
	delete fCmdPhotons;
	delete fCmdTracks;
	delete fAnalysisDirectory;
}

void RunActionMessenger::SetNewValue(G4UIcommand* command, G4String newValue){
	if(command == fCmdFileName){
		fRunAction->SetFileName(newValue);
	}
	else if (command == fCmdOCT){
		fRunAction->SetCmdOCT(fCmdOCT->GetNewBoolValue(newValue));
	}
	else if (command == fCmdDN){
		fRunAction->SetCmdDN(fCmdDN->GetNewBoolValue(newValue));
	}
	else if (command == fCmdPhotons){
		fRunAction->SetCmdPhotons(fCmdPhotons->GetNewIntValue(newValue));
	}
	else if (command == fCmdTracks){
		fRunAction->SetCmdTracks(fCmdTracks->GetNewIntValue(newValue));
	}
	else if (command == fCmdGunTime){
		fRunAction->SetGunTimeMean(1 / fCmdGunTime->GetNewDoubleValue(newValue));
	}
}
