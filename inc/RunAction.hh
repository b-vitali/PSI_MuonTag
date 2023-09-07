/// \file  RunAction.cc
/// \brief Definition of the RunAction class, takes care of ROOT output

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "G4Run.hh"				// to have the run number in the file name
#include "G4AnalysisManager.hh"			// to access root stuff

#include "RunActionMessenger.hh"

/// Run action class
///
/// It prints out the data Tree

class RunActionMessenger;

class RunAction : public G4UserRunAction 
{
	public:
		RunAction();
		virtual ~RunAction();

		void SetFileName(G4String name){fName = name;}
		void SetFolderName(G4String name){fFolder = name;}

        // Needed minimal functions from G4UserRunAction
        //virtual G4Run* GenerateRun();
		virtual void BeginOfRunAction(const G4Run*);
		virtual void   EndOfRunAction(const G4Run*);
		
		//inline TFile* GetFilePtr(){return fData;}
		//inline TTree* GetTreePtr(){return fTree;}

	private:
		G4String fFolder, fName;
		RunActionMessenger* fMessenger;

};

#endif