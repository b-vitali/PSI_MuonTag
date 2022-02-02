/// \file  RunAction.cc
/// \brief Definition of the RunAction class, takes care of ROOT output

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "G4Run.hh"				// to have the run number in the file name
#include "g4root.hh"			// to access root stuff

/// Run action class
///
/// It prints out the data Tree

class RunAction : public G4UserRunAction 
{
	public:
		RunAction();
		virtual ~RunAction();

        // Needed minimal functions from G4UserRunAction
        //virtual G4Run* GenerateRun();
		virtual void BeginOfRunAction(const G4Run*);
		virtual void   EndOfRunAction(const G4Run*);
		
		//inline TFile* GetFilePtr(){return fData;}
		//inline TTree* GetTreePtr(){return fTree;}

	private:

};

#endif