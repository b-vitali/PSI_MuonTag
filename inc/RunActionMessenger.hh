/// \file  RunActionMessenger.hh
/// \brief Definition of the RunActionMessenger class

#ifndef RunActionMessenger_h
#define RunActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class RunAction;
class G4UIdirectory;
class G4UIcmdWithABool;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;

/// it implements command:
///  - /analysis/SetFileName name.root

class RunActionMessenger : public G4UImessenger{
	public:
		RunActionMessenger(RunAction*);
		virtual ~RunActionMessenger();

		virtual void SetNewValue(G4UIcommand*, G4String);
	
	private:
		RunAction* fRunAction;
		
		G4UIdirectory* fAnalysisDirectory;
		
		G4UIcmdWithAString*   fCmdFileName;
		G4UIcmdWithAString*   fCmdFolderName;
};

#endif


