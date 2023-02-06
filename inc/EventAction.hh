/// \file  EventAction.hh
/// \brief Definition of the EventAction class

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

#include "ScintSD.hh"
#include "SiPMSD.hh"

class RunAction;
class TTree;

class EventAction : public G4UserEventAction{
	public:
		EventAction(RunAction* runAction);
		virtual ~EventAction();

		virtual void BeginOfEventAction(const G4Event*);
		virtual void   EndOfEventAction(const G4Event*);

	private:
		RunAction* fRunAction;
		
		ScintSD * tmp_scint_gate;
		ScintSD * tmp_scint_telescope;
		SiPMSD * tmp_sipm_telescope;
		SiPMSD * tmp_sipm_gate;

		G4int fCollIDScint_gate;
		G4int fCollIDScint_telescope;
		G4int fCollIDSiPM_gate;
		G4int fCollIDSiPM_telescope;
		
		G4int fEvID; // to register each event just once
};

#endif


