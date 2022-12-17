/// \file  PrimaryGeneratorAction.hh
/// \brief Implementation of the PrimaryGeneratorAction class, it is the particle gun

#include "PrimaryGeneratorAction.hh"

int wich = 2;

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
    fParticleGun = new G4ParticleGun(1);

    // Default in the construction so that it can be overwritten later
    G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition * particle = particleTable->FindParticle("mu+");

    G4ThreeVector pos(0.,0.,-15*cm);
    G4ThreeVector mom(0.,0.,1.);

    fParticleGun->SetParticlePosition(pos);
    fParticleGun->SetParticleMomentumDirection(mom);
    fParticleGun->SetParticleMomentum(28.*MeV);
    fParticleGun->SetParticleDefinition(particle);

}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete fParticleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent)
{
    G4ThreeVector pos(0.,0.,-15*cm);
    G4ThreeVector mom(0.,0.,1.);
    if(wich == 0){
        bool ok=false;
        while(!ok){
            G4double cosTheta = - ( 1*G4UniformRand() - 1. ), phi = twopi*G4UniformRand();
            G4double sinTheta = std::sqrt(1. - cosTheta*cosTheta);
            G4double    ux = sinTheta*std::cos(phi),
                        uy = sinTheta*std::sin(phi),
                        uz = cosTheta;
            if(true) {
                ok = true;
                pos = G4ThreeVector(0.,0.,-10.4*cm);
                pos = G4ThreeVector(0.,0.,0*cm);
                mom = G4ThreeVector(ux,uy,uz);
            }
        }
        fParticleGun->SetParticlePosition(pos);
        fParticleGun->SetParticleMomentumDirection(mom);
    }

    else if(wich == 1){
        bool ok=false;
        while(!ok){
            G4double x = G4RandGauss::shoot(0,3);
            G4double y = G4RandGauss::shoot(0,3);

            G4double cosTheta = G4RandGauss::shoot(0,0.2), phi = twopi*G4UniformRand();
            G4double sinTheta = std::sqrt(1. - cosTheta*cosTheta);
            G4double    ux = sinTheta*std::cos(phi),
                        uy = sinTheta*std::sin(phi),
                        uz = cosTheta;
            if((x*x+y*y< 10*10)){
                ok = true;
                pos = G4ThreeVector(x,y,-15*cm);
                mom = G4ThreeVector(ux,uy,uz);
            }   
        }
        fParticleGun->SetParticlePosition(pos);
        fParticleGun->SetParticleMomentumDirection(mom);
    }

    else if(wich == 2){
        pos = G4ThreeVector(0, 0, 0*cm);
        mom = G4ThreeVector(0, 1, 0);
        fParticleGun->SetParticlePosition(pos);
        fParticleGun->SetParticleMomentumDirection(mom);
    }

    fParticleGun->GeneratePrimaryVertex(anEvent);
}
