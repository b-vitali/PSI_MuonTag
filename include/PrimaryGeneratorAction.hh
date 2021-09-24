/// \file  PrimaryGeneratorAction.hh
/// \brief Definition of PrimaryGeneratorAction class

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class PGActionMessenger;

///Primary generator action class with particle gun
///
/// The default energy is 2 MeV, monochromatic positron

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction{
    public:
        PrimaryGeneratorAction();
        virtual ~PrimaryGeneratorAction();

        //method of the base class
        virtual void GeneratePrimaries(G4Event*);
	const G4ParticleGun* GetParticleGun() const {return fParticleGun;}
	void SetOCT(G4bool val){fOCTrun = val;}
	void SetDivergence(G4double val){fDivergence = val;}
    
    private:
        G4ParticleGun* fParticleGun; // pointer to G4 gun class
	PGActionMessenger* fMessenger;

	G4double fDivergence;

	G4bool fOCTrun;
};

#endif


