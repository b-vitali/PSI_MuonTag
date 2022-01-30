/// \file  ActionInitialization.cc
/// \brief Implementation on the ActionInitialization class

#include "ActionInitialization.hh"

ActionInitialization::ActionInitialization() : G4VUserActionInitialization(){}

ActionInitialization::~ActionInitialization(){}

void ActionInitialization::BuildForMaster() const{
    //SetUserAction(new RunAction);
}

void ActionInitialization::Build() const{
    
    /*
    RunAction* runAction = new RunAction();
    SetUserAction(runAction);

    SetUserAction(new EventAction(runAction));
     
    PrimaryGeneratorAction *generator = new PrimaryGeneratorAction();
    SetUserAction(generator);
    */
}


