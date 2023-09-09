/// \file  scintsim.cc
/// \brief Main of the scintillator simulation

// Stuff needed to have some handle
#include <TTree.h>
#include <iostream>
#include <cmath>
//#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
//#else
#include "G4RunManager.hh"
//#endif

#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4VisManager.hh"
#include "G4VisExecutive.hh"

// My files to define the detector, the action and the physics list
#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"
#include "Physics.hh"


int main(int argc, char** argv){
    // Detect interactive mode (if no macro) and define UI session
    G4UIExecutive* ui = 0;

	// Construct the default run manager
    #ifdef G4MULTITHREADED
      G4MTRunManager* runManager = new G4MTRunManager;
    #else
      G4RunManager* runManager = new G4RunManager;
    #endif


    // Set mandatory initialization classes
    MyPhysicsList * myphysicslist = new MyPhysicsList();
    runManager->SetUserInitialization(new DetectorConstruction());
    runManager->SetUserInitialization(myphysicslist->GetPhysicsList());
    runManager->SetUserInitialization(new ActionInitialization());

    if(argc == 1){
        ui = new G4UIExecutive(argc, argv);
    }

    //Initialize visualization
    G4VisManager* visManager = new G4VisExecutive();
    visManager->Initialize();

    // Get the pointer to the User Interface manager
    G4UImanager* UImanager = G4UImanager::GetUIpointer();

    // Process macro or start UI session
    if(!ui){
        // batch mode
        G4String command  = "/control/execute ";
		G4String fileName = argv[1];
        UImanager->ApplyCommand(command + fileName);
    }
    else{
        // interactive mode
        UImanager->ApplyCommand("/control/execute mac/vis.mac");
        /*
        UImanager->ApplyCommand("/run/initialize");
        UImanager->ApplyCommand("/vis/open OGL 600x600-0+0");
        UImanager->ApplyCommand("/vis/drawVolume");
        UImanager->ApplyCommand("/vis/viewer/set/viewpointVector -1 1 1");
        //UImanager->ApplyCommand("/vis/scene/add/axes");
        //UImanager->ApplyCommand("/vis/scene/add/scale");
        UImanager->ApplyCommand("/vis/viewer/set/style s");
        UImanager->ApplyCommand("/vis/viewer/set/autorefresh true");
        UImanager->ApplyCommand("/vis/scene/add/trajectories smooth");
        UImanager->ApplyCommand("/vis/scene/endOfEventAction accumulate");
        */
        ui->SessionStart();
    }

    // Job termination
    delete ui;
    delete visManager;
    delete runManager;
}


