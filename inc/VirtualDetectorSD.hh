/// \file  VirtualDetectorSD.hh
/// \brief Definition the SensitiveDetector for the VirtualDetector

#ifndef VirtualDetectorSD_h
#define VirtualDetectorSD_h 1

// Basic class
#include "G4VSensitiveDetector.hh"

// Needed to get the infos
#include "G4ThreeVector.hh"

class G4Step;

class VirtualDetectorSD : public G4VSensitiveDetector
{
    public:
        VirtualDetectorSD(G4String name);
        virtual ~VirtualDetectorSD();

    public:
        virtual G4bool ProcessHits(G4Step *, G4TouchableHistory *);

    private:

};

#endif