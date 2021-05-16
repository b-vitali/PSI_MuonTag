
#include "YourRunAction.hh"

#include "YourDetectorConstruction.hh"
#include "YourPrimaryGeneratorAction.hh"
#include "YourRun.hh"
#include "G4Run.hh"
#include "YourRunActionMessenger.hh"

#include "Randomize.hh"

#include "MyAnalysis.hh"

YourRunAction::YourRunAction(YourDetectorConstruction* det, YourPrimaryGeneratorAction* prim) 
:   G4UserRunAction(),
    fYourDetector(det),
    fYourPrimary(prim),
    fYourRun(nullptr),
    fIsEdepHistogramUICmdInvoked(false) { 
  /* histo manager could be created here */ 
  // Create our own UI messenger object that will interact to this Run-Action
  // to set some properties that will be used to update YourRun object (generated 
  // by calling this Run-Action::GenerateRun method) in the BeginOfRunAction method.
  fMessenger = new YourRunActionMessenger(this);
}

YourRunAction::~YourRunAction() { 
  /* histo manager must be deleted here then*/ 
  // delete all dynamically allocated objects here 
  delete fMessenger;
}


G4Run* YourRunAction::GenerateRun() {
	fYourRun = new YourRun(fYourDetector, fYourPrimary);
	return fYourRun;
}


void YourRunAction::BeginOfRunAction(const G4Run* /*run*/) {
  // Show Rndm status (only for the Master thread)
  //if ( IsMaster() ) G4Random::showEngineStatus();
  //
  // Make sure that the Gun position is correct: the user can change the target
  // thickness between construction of objects and start of the run.  
  // note: primary generator is set in the CTR only for the Worker threads in the 
  //       ActionInitialization (left null for Master in the BuildForMaster())
  if ( fYourPrimary ) { 
      fYourPrimary->UpdatePosition();
  }
  // Update the properties of the Energy-deposit histogram member of YourRun, 
  // that is already available at this point: Only if the user invoked the UI 
  // command /yourApp/runAction/edepHistoto set properties of the Edep-histo.
  if ( fIsEdepHistogramUICmdInvoked ) {
    // user defined the properties of the Edep-histo by invoking the UI command 
    fYourRun->SetEdepHisto(fEdepHistFileName, fEdepHistMinEnergy, fEdepHistMaxEnergy, fEdepHistNumBins);
  }                        
  //
  // G4AnalysisManager* analysisManager OpenFile

  // Create/get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();  
  analysisManager->SetVerboseLevel(1);
 
  analysisManager->SetHistoDirectoryName("histo");
  analysisManager->SetNtupleDirectoryName("ntuple"); 
  
  // Open an output file  
  analysisManager->OpenFile("MyApplication"); 


  // Create histograms  
  analysisManager->CreateH1("Edep","Energy deposit", 1250, 0., 125, "MeV");  
  analysisManager->CreateH1("Tlen","Track length", 10000, 0., 1000,"mm");
  analysisManager->CreateH2("Edep_vs_Tlen","Energy deposit vs Track length", 10000, 0., 1000, 1250, 0., 125, "mm","MeV");

  analysisManager->CreateNtuple("MyNtuple", "Edep and TrackLength");
  analysisManager->CreateNtupleDColumn("Eabs");
  analysisManager->CreateNtupleDColumn("Labs");
  analysisManager->FinishNtuple();

  analysisManager->CreateNtuple("Positions", "In Out x,y,z");
  analysisManager->CreateNtupleDColumn("in_x");
  analysisManager->CreateNtupleDColumn("in_y");
  analysisManager->CreateNtupleDColumn("in_z");
  analysisManager->CreateNtupleDColumn("out_x");
  analysisManager->CreateNtupleDColumn("out_y");
  analysisManager->CreateNtupleDColumn("out_z");
  analysisManager->FinishNtuple();

}


void YourRunAction::EndOfRunAction(const G4Run*) {
  // Print Run summary (only for the Master thread)
  if ( IsMaster() ) { 
  	fYourRun->EndOfRunSummary();    
  }
  //
  // Show Rndm status (only for the Master thread)
  //  if ( IsMaster() ) G4Random::showEngineStatus();
	//
	// G4AnalysisManager* analysisManager Write and CloseFile

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  // Write and close the output file
  analysisManager->Write();  
  analysisManager->CloseFile();
}

