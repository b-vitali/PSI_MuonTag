/// \file  PrimaryGeneratorAction.hh
/// \brief Implementation of the PrimaryGeneratorAction class, it is the particle gun

#include "PrimaryGeneratorAction.hh"

#include "Randomize.hh"
#include "G4AnalysisManager.hh"			// to access root stuff
#include "G4RandomDirection.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
    fParticleGun = new G4ParticleGun(1);

    // Default in the construction so that it can be overwritten later
    G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition * particle = particleTable->FindParticle("e+");

    G4ThreeVector pos(-0.3*mm,0.,-10.*mm);
    G4ThreeVector mom(0.,1.,1.);
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
/*    
    G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition * particle = particleTable->FindParticle("e+");

    // //? Position uniform in on the 3cm circle
    G4double phi_pos = G4UniformRand();
    G4ThreeVector pos(sin(phi_pos*2*CLHEP::pi),0,cos(phi_pos*2*CLHEP::pi));
    pos = 3*cm*pos;

    // //? Random direction of the momentum
    G4double phi = G4UniformRand();
    G4double theta = G4UniformRand();

    // double mom_phi = phi_pos*2*CLHEP::pi-0.5*CLHEP::pi;
    // double mom_theta = 0*CLHEP::pi/180;
    // G4ThreeVector mom(sin(mom_phi)*cos(mom_theta), sin(mom_theta), cos(mom_phi)*cos(mom_theta));

    G4ThreeVector mom = G4RandomDirection();
    G4cout<<mom.mag()<<G4endl;
    fParticleGun->SetParticlePosition(pos);

    G4double p = G4UniformRand();
    p= (p*5+30)*MeV;
    mom = p*mom;
    fParticleGun->SetParticleMomentum(mom); //28.*MeV
    fParticleGun->SetParticleDefinition(particle);
    
    G4AnalysisManager *man = G4AnalysisManager::Instance();
    man->FillNtupleDColumn(1, 0, pos.getX());
    man->FillNtupleDColumn(1, 1, pos.getY());
    man->FillNtupleDColumn(1, 2, pos.getZ());
    man->FillNtupleDColumn(1, 3, mom.getX());
    man->FillNtupleDColumn(1, 4, mom.getY());
    man->FillNtupleDColumn(1, 5, mom.getZ());
    man->AddNtupleRow(1);
*/

    fParticleGun->GeneratePrimaryVertex(anEvent);

    // fParticleGun->SetParticlePosition(G4ThreeVector(-0.3*mm,0.,-10.*mm));
    // fParticleGun->SetParticleTime(0*ns);
    // fParticleGun->GeneratePrimaryVertex(anEvent);
    // fParticleGun->SetParticlePosition(G4ThreeVector(-0.3*mm+5*mm,0.,-10.*mm));
    // fParticleGun->SetParticleTime(0*ns);
    // fParticleGun->GeneratePrimaryVertex(anEvent);
}
