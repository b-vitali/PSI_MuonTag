#include "YourEventAction.hh"

#include "G4RunManager.hh"
#include "YourRun.hh"

#include "MyAnalysis.hh"

YourEventAction::YourEventAction() 
: G4UserEventAction(),
	fEdepPerEvt(0.0),
  fChTrackLengthPerEvt(0.0),
  in_x(0.0), in_y(0.0), in_z(0.0) {}


YourEventAction::~YourEventAction() {}


// Beore each event: reset per-event variables 
void YourEventAction::BeginOfEventAction(const G4Event* /*anEvent*/) {
	fEdepPerEvt           = 0.0;
  fChTrackLengthPerEvt  = 0.0;	
  in_x                  = 0.0;
  in_y                  = 0.0;
  in_z                  = 0.0;
  out_x                 = 0.0;
  out_y                 = 0.0;
  out_z                 = 0.0;
}


// After each event:
// fill the data collected for this event into the Run global (thread local)
// data Run data object (i.e. into YourRun)  
void YourEventAction::EndOfEventAction(const G4Event* /*anEvent*/) {
	// get the current Run object and cast it to YourRun (because for sure this is its type)
	YourRun* currentRun = static_cast< YourRun* > ( G4RunManager::GetRunManager()->GetNonConstCurrentRun() );
  // add the quantities to the (thread local) run global YourRun object 
  currentRun->AddEnergyDepositPerEvent( fEdepPerEvt );
  currentRun->AddChTrackLengthPerEvent( fChTrackLengthPerEvt );
	currentRun->FillEdepHistogram( fEdepPerEvt );

  // Get analysis manager  
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();  
  // Fill histograms  
  analysisManager->FillH1(0, fEdepPerEvt);  
  analysisManager->FillH1(1, fChTrackLengthPerEvt);
  analysisManager->FillH2(0, fChTrackLengthPerEvt, fEdepPerEvt);

  analysisManager->FillNtupleDColumn(0,0, fEdepPerEvt);
  analysisManager->FillNtupleDColumn(0,1, fChTrackLengthPerEvt);

  analysisManager->AddNtupleRow(0);

  analysisManager->FillNtupleDColumn(1,0, in_x);
  analysisManager->FillNtupleDColumn(1,1, in_y);
  analysisManager->FillNtupleDColumn(1,2, in_z);
  analysisManager->FillNtupleDColumn(1,3, out_x);
  analysisManager->FillNtupleDColumn(1,4, out_y);
  analysisManager->FillNtupleDColumn(1,5, out_z);

  analysisManager->AddNtupleRow(1);

}