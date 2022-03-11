/// \file  ActionInitialization.hh
/// \brief Definition of the ActionInitialization class

#ifndef ActionInitialization_h
#define ActionInitialization_h 1

// Basic class
#include "G4VUserActionInitialization.hh"

// Dependences wit PrimaryGenerator Run and Event
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"


class ActionInitialization : public G4VUserActionInitialization{
    public:
        ActionInitialization();
        virtual ~ActionInitialization();

        // This is basically the main function of our simulation. 
        // Takes care of a lot of stuff
        virtual void BuildForMaster() const;
        virtual void Build() const;
};

#endif

