/// \file  VirtualDetectorSD.hh
/// \brief Definition of the SensitiveDetector for the VirtualDetector

#ifndef VirtualDetectorSD_h
#define VirtualDetectorSD_h 1

// Basic class
#include "G4VSensitiveDetector.hh"

// Needed to get the infos
#include "G4ThreeVector.hh"

#include "G4RunManager.hh"

#include "G4AnalysisManager.hh"

class G4Step;

class VirtualDetectorSD : public G4VSensitiveDetector
{
    public:
        VirtualDetectorSD(G4String name);
        virtual ~VirtualDetectorSD();

    public:
        virtual G4bool ProcessHits(G4Step *, G4TouchableHistory *);
        G4bool FillHit(G4Step *aStep);

    private:
        G4int       fVDNo;      // which VirtualDetector?
        G4int       fParticleID;// pdg particle ID
        G4int      fInOut;// pdg particle ID
        G4double    fVDTime;    // time of hit
        G4double    fMom;    // entering momentum
        G4double    fMomX, fMomY, fMomZ;    // entering momentum
        G4double    fPosX, fPosY, fPosZ;    // entering position
};

#endif