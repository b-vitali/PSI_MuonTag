#ifndef YOUREVENTACTION_HH
#define YOUREVENTACTION_HH

#include "G4UserEventAction.hh"

#include "globals.hh"

class YourEventAction : public G4UserEventAction {

  // Method declaration:
  public:
    
    // CTR: 
    YourEventAction();
    virtual ~YourEventAction();

    // Two virtual method called at the beginning and at the end of each event 
    virtual void BeginOfEventAction(const G4Event* anEvent);
    virtual void EndOfEventAction(const G4Event* anEvent);

    // Public method (used by the SteppingAction) to add energy deposit and 
    // charged particle track length information to this object after each step.
    void AddEnergyDepositPerStep(const G4double edepAlongTheStep) { 
        fEdepPerEvt += edepAlongTheStep;
    }
    void AddChargedTrackLengthPerStep(const G4double trackLAlongTheStep) { 
    	fChTrackLengthPerEvt += trackLAlongTheStep;
    }

    void PositionIn(G4double a, G4double b, G4double c) { 
      in_x = a;
      in_y = b;
      in_z = c;
    }

    void PositionOut(G4double a, G4double b, G4double c) { 
      out_x = a;
      out_y = b;
      out_z = c;
    }

  // Data member declarations:
  private:
    G4double   fEdepPerEvt;
    G4double   fChTrackLengthPerEvt;
    G4double   in_x, in_y, in_z;
    G4double   out_x, out_y, out_z;

};

#endif