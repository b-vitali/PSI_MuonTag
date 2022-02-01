/// \file  ActionInitialization.hh
/// \brief Definition of the ActionInitialization class

#ifndef ActionInitialization_h
#define ActionInitialization_h 1

// Basic class
#include "G4VUserActionInitialization.hh"

#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"

/*
// Dependences wit Run Event and PrimaryGenerator
#include "RunAction.hh"
#include "EventAction.hh"
*/

class ActionInitialization : public G4VUserActionInitialization{
    public:
        ActionInitialization();
        virtual ~ActionInitialization();

        // This is basically the main function of our simulation. 
        //Takes care of a lot of stuff
        virtual void BuildForMaster() const;
        virtual void Build() const;
};

#endif

