/// \file  PrimaryGeneratorAction.hh
/// \brief Implementation of the PrimaryGeneratorAction class, it is the particle gun

#include "PrimaryGeneratorAction.hh"

#include "Randomize.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
    fParticleGun = new G4ParticleGun(1);

    // Default in the construction so that it can be overwritten later
    G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition * particle = particleTable->FindParticle("mu+");

    G4ThreeVector pos(0.,0.,-10.*mm);
    G4ThreeVector mom(0.,0.,1.);

    fParticleGun->SetParticlePosition(pos);
    fParticleGun->SetParticleMomentumDirection(mom);
    fParticleGun->SetParticleMomentum(100.*MeV);
    fParticleGun->SetParticleDefinition(particle);

}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete fParticleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent)
{
    G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition * particle = particleTable->FindParticle("e+");

    //? Random direction of the momentum
    G4double phi = G4UniformRand();
    G4double theta = G4UniformRand();

    G4ThreeVector mom(0,0,1);
    mom.setMag(1);
    mom.setPhi(phi*2*CLHEP::pi);
    mom.setTheta(phi*CLHEP::pi); //sqrt(-2*log(theta))*cos(2*CLHEP::pi*theta)*CLHEP::pi
   
    //? Position uniform in on the 3cm circle
    G4double phi_pos = G4UniformRand();
    G4ThreeVector pos(2.5*cm*cos(phi_pos*2*CLHEP::pi),0,2.5*cm*sin(phi_pos*2*CLHEP::pi));

    fParticleGun->SetParticlePosition(pos);
    fParticleGun->SetParticleMomentumDirection(mom);
    fParticleGun->SetParticleMomentum(100.*GeV);
    fParticleGun->SetParticleDefinition(particle);

    fParticleGun->GeneratePrimaryVertex(anEvent);
}
