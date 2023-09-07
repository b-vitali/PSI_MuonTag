/// \file  PrimaryGeneratorAction.hh
/// \brief Definition of the PrimaryGeneratorAction class, it is the particle gun

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"

#include "TDecompChol.h"
#include "TMatrixD.h"
#include "G4AnalysisManager.hh"			// to access root stuff

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    PrimaryGeneratorAction();
    ~PrimaryGeneratorAction();

    virtual void GeneratePrimaries(G4Event*);

private:
    G4ParticleGun *fParticleGun;
};

#endif
