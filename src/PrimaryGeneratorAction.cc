/// \file  PrimaryGeneratorAction.hh
/// \brief Implementation of the PrimaryGeneratorAction class, it is the particle gun

#include "PrimaryGeneratorAction.hh"

G4bool flat_angle = false;

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
    fParticleGun = new G4ParticleGun(1);

    // Default in the construction so that it can be overwritten later
    G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition * particle = particleTable->FindParticle("mu+");

    G4ThreeVector pos(0.,0.,-15*cm);
    G4ThreeVector mom(0.,0.,1.);

    if(flat_angle) pos = G4ThreeVector(0.,0.,5.*cm);


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
    if(flat_angle){
        G4double cosTheta = - ( 1*G4UniformRand() - 1. ), phi = twopi*G4UniformRand();
        G4double sinTheta = std::sqrt(1. - cosTheta*cosTheta);
        G4double    ux = sinTheta*std::cos(phi),
                    uy = sinTheta*std::sin(phi),
                    uz = cosTheta;

    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(ux,uy,uz));
    }

    fParticleGun->GeneratePrimaryVertex(anEvent);
}
