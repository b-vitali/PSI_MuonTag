/// \file  VirtualDetectorSD.hh
/// \brief Definition of the SensitiveDetector for the VirtualDetector

#ifndef VirtualDetectorSD_h
#define VirtualDetectorSD_h 1

// Basic class
#include "G4VSensitiveDetector.hh"

// Needed to get the infos
#include "G4ThreeVector.hh"

#include "G4RunManager.hh"

#include "g4root.hh"

class G4Step;

class VirtualDetectorSD : public G4VSensitiveDetector
{
    public:
        VirtualDetectorSD(G4String name);
        virtual ~VirtualDetectorSD();

    public:
        virtual G4bool ProcessHits(G4Step *, G4TouchableHistory *);

    private:
        G4int       fVDNo;      // which VirtualDetector?
        G4int       fParticleID;// pdg particle ID
        G4double    fVDTime;    // time of hit
        G4double    fMom;    // entering momentum
};

#endif