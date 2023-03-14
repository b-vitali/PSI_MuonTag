/// \file  PrimaryGeneratorAction.hh
/// \brief Implementation of the PrimaryGeneratorAction class, it is the particle gun

#include "PrimaryGeneratorAction.hh"

/*
    1   
    2   move along the z axis and vertical p
    3   Cholesky: A covariance, A = LL*, u uncorrelated variables -> Lu are correlated variables
*/
int which = 3;
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
        G4double mu_x  = 0.;
        G4double mu_xp = 0.;
        G4double mu_y  = 0.;
        G4double mu_yp = 0;

        G4double data[4][4] = { {1.,  0.,  0.,   0.}, {0.,  0.005,  0.,   0.}, {0.,  0.,  1.,   0.}, {0.,  0.,  0.,   0.003} };
        TMatrixD cov(4, 4, *data); // Creates a matrix with specified data

        G4double x  =  G4RandGauss::shoot(mu_x  ,sqrt(data[0][0]));
        G4double xp =  G4RandGauss::shoot(mu_xp ,sqrt(data[1][1]));
        G4double y  =  G4RandGauss::shoot(mu_y  ,sqrt(data[2][2]));
        G4double yp =  G4RandGauss::shoot(mu_yp ,sqrt(data[3][3]));


        x = mu_x; y = mu_y;
        G4double v_uncor[4] = {x,xp,y,xp};

        TDecompChol chol(cov);
        chol.Decompose();
        TMatrixD U = chol.GetU();

        TVectorD v_cor = U * TVectorD(4, v_uncor);

        // assuming z is uncoupled to xy. NB approx for now
        G4double mu_p = 28;
        G4double sigma_p = 0;
        G4double p =  G4RandGauss::shoot(mu_p ,sigma_p);
        G4ThreeVector pos(v_cor[0], v_cor[2], -10.5*cm);
        G4ThreeVector mom(sin(xp), sin(yp), 1.);

        mom = mom * p;

        fParticleGun->SetParticlePosition(pos);    
        fParticleGun->SetParticleMomentum(mom*MeV);
        // fParticleGun->SetParticleMomentum(p*MeV);
        // fParticleGun->SetParticleMomentumDirection(mom);

        G4AnalysisManager *man = G4AnalysisManager::Instance();
        man->FillNtupleIColumn(1, 0, i);
        man->FillNtupleDColumn(1, 1, p);
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
