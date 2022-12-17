/// \file  ActionInitialization.cc
/// \brief Implementation of the ActionInitialization class

#include "ActionInitialization.hh"

ActionInitialization::ActionInitialization() : G4VUserActionInitialization(){}

ActionInitialization::~ActionInitialization(){}

void ActionInitialization::BuildForMaster() const{
    RunAction* runAction = new RunAction();
    SetUserAction(runAction);
}

void ActionInitialization::Build() const{
    SetUserAction(new PrimaryGeneratorAction);

    
    RunAction* runAction = new RunAction();
    SetUserAction(runAction);
    
    SetUserAction(new EventAction(runAction));
}


