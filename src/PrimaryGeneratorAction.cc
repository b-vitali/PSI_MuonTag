/// \file  PrimaryGeneraotrAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"
#include "PGActionMessenger.hh"

#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include "TMath.h"


PrimaryGeneratorAction::PrimaryGeneratorAction() : 
    G4VUserPrimaryGeneratorAction(), fParticleGun(0), fDivergence(0){
    G4int n_particle = 1;
    fParticleGun = new G4ParticleGun(n_particle);

    //default particle kinematic
    auto particle = G4ParticleTable::GetParticleTable()->FindParticle("e+");
    fParticleGun->SetParticleDefinition(particle);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., -1.));
    fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.));
    fParticleGun->SetParticleEnergy(2.2*MeV);

    //optical photon for OCT

    fOCTrun = false;
    fMessenger = new PGActionMessenger(this);
}


PrimaryGeneratorAction::~PrimaryGeneratorAction(){
    delete fParticleGun;
    delete fMessenger;
}


void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent){
    G4ParticleDefinition* particle = fParticleGun->GetParticleDefinition();
	    
    // In order to avoid dependence of PrimaryGeneratorAction
    // on DetectorConstruction class we get world volume from G4LogicalVolumeStore.

    G4double world_size = 0;
    auto worldLV = G4LogicalVolumeStore::GetInstance()->GetVolume("World");
    auto scintLV = G4LogicalVolumeStore::GetInstance()->GetVolume("Element");
    auto scintPV = G4PhysicalVolumeStore::GetInstance()->GetVolume("Element");

    G4Box* worldBox = dynamic_cast<G4Box*>(worldLV->GetSolid());
    G4Box* scintBox = dynamic_cast<G4Box*>(scintLV->GetSolid());

    world_size = worldBox->GetXHalfLength();

    if(fDivergence > 0){
	world_size = scintPV->GetTranslation().z() + scintBox->GetZHalfLength();
	G4double phi = G4UniformRand() * 2 * TMath::Pi();
	G4double theta = (G4UniformRand() - 1./2) * TMath::Pi();
//	while (fabs(theta) > TMath::Pi()/2) theta = G4RandGauss::shoot(0, fDivergence);
	fParticleGun->SetParticleMomentumDirection(G4ThreeVector(sin(theta)*cos(phi), sin(theta)*sin(phi), -cos(theta)));
    }
    else if (fDivergence == 0){
	fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,-1));
    }
    // Set gun position
    if(!fOCTrun){
	    fParticleGun->SetParticleDefinition(particle);
	    fParticleGun->SetParticlePosition(G4ThreeVector(0, 0, world_size)); //0.5*scintBox->GetYHalfLength()
	    fParticleGun->GeneratePrimaryVertex(anEvent);
    }
    else{
	    particle = G4ParticleTable::GetParticleTable()->FindParticle("opticalphoton");
	    fParticleGun->SetParticleDefinition(particle);
	    fParticleGun->SetParticlePosition(G4ThreeVector(G4UniformRand() * 1.3 - 1.3*0.5,G4UniformRand() * 1.3 - 1.3*0.5, world_size));
	    fParticleGun->SetParticleEnergy(2.8*eV);
	    fParticleGun->GeneratePrimaryVertex(anEvent);
    }
}


