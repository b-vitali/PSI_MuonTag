/// \file  EventAction.hh
/// \brief Definition of the EventAction class

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

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
		
		G4int fCollIDScint;
		G4int fCollIDSiPM;
		G4int fCollIDSiPMDraw;
		
		G4int fEvID; // to register each event just once

		TTree* fTree;
};

#endif


