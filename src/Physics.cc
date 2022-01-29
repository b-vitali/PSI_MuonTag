/// \file  Physics.cc
/// \brief Define our PhysicsList : how the particles interact

#include "Physics.hh"

MyPhysicsList::MyPhysicsList()
{
	
	G4VModularPhysicsList* physicsList = new FTFP_BERT;
	
	physicsList->ReplacePhysics(new G4EmStandardPhysics_option4());
	G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();

	//opticalPhysics->SetWLSTimeProfile("delta");

	//opticalPhysics->SetScintillationYieldFactor(1.0);
	//opticalPhysics->SetScintillationExcitationRatio(0.0);

	//opticalPhysics->SetMaxNumPhotonsPerStep(100);
	//opticalPhysics->SetMaxBetaChangePerStep(10.0);

	//opticalPhysics->SetTrackSecondariesFirst(kCerenkov, true);
	opticalPhysics->SetTrackSecondariesFirst(kScintillation, true);

	physicsList->RegisterPhysics(opticalPhysics);
	
}

MyPhysicsList::~MyPhysicsList(){}