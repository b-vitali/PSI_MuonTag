#ifndef YOURRUN_HH
#define YOURRUN_HH

#include "G4Run.hh"

#include "Hist.hh"

// forward declarations
class YourDetectorConstruction;
class YourPrimaryGeneratorAction;

class YourRun : public G4Run {

  // Method declaration:
  public:
    
    // CTR: 
    YourRun(YourDetectorConstruction* det, YourPrimaryGeneratorAction* prim);
    virtual ~YourRun();

    // Virtual method to be implemented to define the way of merging the underlying 
    // (thread local) Run global data structures into one global instance.
    virtual void Merge(const G4Run*);

   
    // Method to be called by the Master to compute final quantities using the 
    // its run global Run object into which all thread local run global run 
    // object instance has already been merged
    void  EndOfRunSummary();


    // Public methods (will be used by the EventAction) to add energy deposit and 
    // track length data after each events.
    void AddEnergyDepositPerEvent(const G4double edep)   { 
    	fEdepInTarget  += edep;  
    	fEdepInTarget2 += edep*edep; 
    }

    void AddChTrackLengthPerEvent(const G4double length) {
        fChargedTrackLengthInTarget  += length;
        fChargedTrackLengthInTarget2 += length*length;
    }
    
    // Public method (will be used by the EventAction) to fill the energy deposit
    // histogram after each event (in the EndOfEventAction)
    void FillEdepHistogram(const G4double edep) {
        if ( fIsActiveEdepHistogram ) {
          fEdepHistogram->Fill(edep);
        }
    }

    // Public method to set properties of the energy deposit histogram: this will 
    // be invoked by the RunAction in its BeginOfRunAction method to set all 
    // properties of the histogram before the simulation starts (only if the user 
    // invoked the corresponding UI command i.e. /yourApp/runAction/edepHisto ).
    // A flag (fIsActiveEdepHisto) will also set to indicate that the user invoked 
    // the UI command /yourApp/runAction/edepHisto ... that activates the energy 
    // deposit histogram. The flag stays to be 'false' otherwise (as initialised).
    void SetEdepHisto(const G4String& filename, G4double emin, G4double emax, G4int numbins) {
        fIsActiveEdepHistogram = true; 
        fEdepHistogram->ReSet(filename, emin, emax, numbins);
    }


  // Data member declarations:
  private:
    // data members to obtain some information needed at the end for the summary
    // Note: that the PrimaryGeneratorAction is not set to the master's RunAction
    //       so it will be nullptr for the master's YourRun. We will set it in 
    //       the Merge.
  	YourDetectorConstruction*    fYourDetector;
  	YourPrimaryGeneratorAction*  fYourPrimary;   
  	//
  	// Run global data members: 
  	// - sum of the energy deposit (charged track length) in the target per event 
  	// - and sum of the squared quantity per event for rms computation
  	G4double fEdepInTarget;
  	G4double fEdepInTarget2;

  	G4double fChargedTrackLengthInTarget;
  	G4double fChargedTrackLengthInTarget2;
    
    // Energy deposit per event (in the target) histogram and a flag to indicate 
    // if this histogram was activated (requested) by the user invoking the 
    // corresponding custum /yourApp/runAction/edepHisto UI command
    G4bool   fIsActiveEdepHistogram;
    Hist*    fEdepHistogram;
    
    
};

#endif