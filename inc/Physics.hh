/// \file  Physics.cc
/// \brief Define our PhysicsList : how the particles interact

#ifndef Physics_h
#define Physics_h 1

#include "G4VModularPhysicsList.hh"
#include "FTFP_BERT.hh"
#include "G4OpticalPhysics.hh"
#include "G4EmStandardPhysics_option4.hh"

class MyPhysicsList : public G4VModularPhysicsList
{
	public:
		MyPhysicsList();
		~MyPhysicsList();
};

#endif
