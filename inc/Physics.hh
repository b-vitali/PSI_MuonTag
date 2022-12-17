/// \file  Physics.cc
/// \brief Define our PhysicsList : how the particles interact

#ifndef Physics_h
#define Physics_h 1

#include "G4VModularPhysicsList.hh"
#include "FTFP_BERT.hh"
#include "G4OpticalPhysics.hh"
#include "G4EmStandardPhysics_option4.hh"

#include "G4EmStandardPhysics.hh"
#include "G4OpticalPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4DecayPhysics.hh"

class MyPhysicsList : public G4VModularPhysicsList
{
	public:
		MyPhysicsList();
		~MyPhysicsList();
		G4VModularPhysicsList * GetPhysicsList(){return physicsList;};
	private:	
		G4VModularPhysicsList* physicsList;
};

#endif
