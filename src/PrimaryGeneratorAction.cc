/// \file  PrimaryGeneratorAction.hh
/// \brief Implementation of the PrimaryGeneratorAction class, it is the particle gun

#include "PrimaryGeneratorAction.hh"

/*
    1   
    2   move along the z axis and vertical p
    3   Cholesky: A covariance, A = LL*, u uncorrelated variables -> Lu are correlated variables
    4  Gaus Gaus
*/
int which = 4;
int i = 0;

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
    if(which == 0){
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

    else if(which == 1){
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
                pos = G4ThreeVector(x,y,-10.5*cm);
                mom = G4ThreeVector(ux,uy,uz);
            }   
        }
        fParticleGun->SetParticlePosition(pos);
        fParticleGun->SetParticleMomentumDirection(mom);
    }

    else if(which == 2){
        pos = G4ThreeVector(0, 0, 0*cm);
        mom = G4ThreeVector(0, 1, 0);
        fParticleGun->SetParticlePosition(pos);
        fParticleGun->SetParticleMomentumDirection(mom);
    }
    else if(which == 3){
        G4double momentum = 28.0;               // MeV/c, adjust as needed
        G4double momentumSpread = 0.05;          // Spread %
        G4double z = -15*cm;
        
        // Create the covariance matrix
        TMatrixD covariance(6, 6);
        covariance(0,0) = 352.31290;   //mm2
        covariance(1,1) = 0.0031696900;    //rad2
        covariance(2,2) = 504.90090;   //mm2
        covariance(3,3) = 0.0022752900; //rad2
        covariance(4,4) = 1;
        covariance(5,5) = 28*0.05;
        covariance(0,1) = -0.25996075;
        covariance(1,0) = covariance(0,1);
        covariance(2,3) = -0.19721470;
        covariance(3,2) = covariance(2,3);
        covariance(4,5) = 0;
        covariance(5,4) = covariance(4,5);

        // Perform Cholesky decomposition
        TDecompChol chol(covariance);
        chol.Decompose();
        TMatrixD L = chol.GetU();

        // Generate random uncorrelated variables
        TMatrixD v(1, 6);
        for (int i = 0; i < 6; i++) {
            v[0][i] = G4RandGauss::shoot(0.0, 1);
        }

        v(0, 4) = 0;

        for (int i = 0; i < 6; i++) {
            G4cout << "hop "<< v[0][i] << " ";
        }
        // Apply Cholesky transformation to get correlated variables
        v = L*v;
        for (int i = 0; i < 6; i++) {
            G4cout << v[0][i] << " ";
        }

        G4double x   = v(0, 0);
        G4double xp  = v(0, 1);
        G4double y   = v(0, 2);
        G4double yp  = v(0, 3);
        G4double pz  = v(5,0)/sqrt(1+xp*xp + yp*yp);
        G4double px = xp*pz;
        G4double py = yp*pz;

        pos = G4ThreeVector(x, y, z);
        mom = G4ThreeVector(px, py, pz);

        for (int i = 0; i < 3; i++) {
            G4cout << pos[i] << " "<<mom[i]<<G4endl;
        }

        fParticleGun->SetParticlePosition(pos);    
        fParticleGun->SetParticleMomentum(mom*MeV);

        G4AnalysisManager *man = G4AnalysisManager::Instance();
        man->FillNtupleIColumn(1, 0, i);
        man->FillNtupleDColumn(1, 1, mom.mag());
        man->FillNtupleDColumn(1, 2, pos.x());
        man->FillNtupleDColumn(1, 3, pos.y());
        man->FillNtupleDColumn(1, 4, pos.z());
        man->FillNtupleDColumn(1, 5, xp);
        man->FillNtupleDColumn(1, 6, yp);
        man->FillNtupleDColumn(1, 7, mom.x());
        man->FillNtupleDColumn(1, 8, mom.y());
        man->FillNtupleDColumn(1, 9, mom.z());
        man->AddNtupleRow(1);

    }

    else if(which == 4){
        G4double momentum = 28.0;               // MeV/c, adjust as needed
        G4double momentumSpread = 0.05;          // Spread %
        G4double dx = 1;
        G4double dy = 2;
        G4double dxp = 0.005;
        G4double dyp = 0.001;

        G4double x   = G4RandGauss::shoot(0.0, dx);
        G4double y   = G4RandGauss::shoot(0.0, dy);
        G4double z = -15*cm;

        G4double xp   = G4RandGauss::shoot(0.0, dxp);
        G4double yp   = G4RandGauss::shoot(0.0, dyp);

        momentum = G4RandGauss::shoot(momentum, momentumSpread);
        G4double pz  = momentum/sqrt(1+xp*xp + yp*yp);
        G4double px = xp*pz;
        G4double py = yp*pz;

        pos = G4ThreeVector(x, y, z);
        mom = G4ThreeVector(px, py, pz);

        fParticleGun->SetParticlePosition(pos);    
        fParticleGun->SetParticleMomentum(mom*MeV);

        G4AnalysisManager *man = G4AnalysisManager::Instance();
        man->FillNtupleIColumn(1, 0, i);
        man->FillNtupleDColumn(1, 1, mom.mag());
        man->FillNtupleDColumn(1, 2, pos.x());
        man->FillNtupleDColumn(1, 3, pos.y());
        man->FillNtupleDColumn(1, 4, pos.z());
        man->FillNtupleDColumn(1, 5, xp);
        man->FillNtupleDColumn(1, 6, yp);
        man->FillNtupleDColumn(1, 7, mom.x());
        man->FillNtupleDColumn(1, 8, mom.y());
        man->FillNtupleDColumn(1, 9, mom.z());
        man->AddNtupleRow(1);

    }
    G4cout<<i<<G4endl; i+=1;

    fParticleGun->GeneratePrimaryVertex(anEvent);
}
