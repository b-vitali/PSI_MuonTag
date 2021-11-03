/// \file  PGActionMessenger.hh
/// \brief Definition of the PGActionMessenger class

#ifndef PGActionMessenger_h
#define PGActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class PrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcmdWithABool;
class G4UIcmdWithADoubleAndUnit;

class PGActionMessenger : public G4UImessenger{
	public:
		PGActionMessenger(PrimaryGeneratorAction*);
		virtual ~PGActionMessenger();

		virtual void SetNewValue(G4UIcommand*, G4String);

	private:
		PrimaryGeneratorAction* fPGAction;

		G4UIdirectory* fPrimary;

		G4UIcmdWithABool* fCmdOCT;
		G4UIcmdWithADoubleAndUnit* fCmdDivergence;
};

#endif


