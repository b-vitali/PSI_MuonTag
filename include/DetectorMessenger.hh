/// \file DetectorMessenger.hh
/// \brief Definition of the DetectorMessenger class

#ifndef DetectorMessenger_h
#define DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithADouble;

/// it implements command:
/// - /Element/det/SetCrysSize value unit
/// - /Element/det/SiPMmodel string

class DetectorMessenger : public G4UImessenger{
	public:
		DetectorMessenger(DetectorConstruction*);
		virtual ~DetectorMessenger();
		
		virtual void SetNewValue(G4UIcommand*, G4String);
	
	private:
		DetectorConstruction* fDetectorConstruction;
		
		G4UIdirectory* fElementDirectory;
		G4UIdirectory* fDetDirectory;
		
		G4UIcmdWithADoubleAndUnit* fCrysSizeCmd;
		G4UIcmdWithADouble*        fGround;
		G4UIcmdWithADoubleAndUnit* fAngle;
		G4UIcmdWithADoubleAndUnit* fAngleWithOpticalGrease;
		G4UIcmdWithADoubleAndUnit* fTilt;
		G4UIcmdWith3VectorAndUnit* fCrysSizeCmd3;
		G4UIcmdWithAString* fCrysMaterialCmd;
		G4UIcmdWithAString* fSiPMmodelCmd;
};

#endif


