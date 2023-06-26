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
		
		ScintSD * tmp_scint_out;
		ScintSD * tmp_scint_in;

		SiPMSD * tmp_sipm_out;
		SiPMSD * tmp_sipm_in;

		G4int fCollIDScint_out;
		G4int fCollIDScint_in;

		G4int fCollIDSiPM_out;
		G4int fCollIDSiPM_in;
		
		G4int fEvID; // to register each event just once
};

#endif


