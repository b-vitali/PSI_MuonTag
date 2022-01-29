/// \file  DetectorConstruction.hh
/// \brief Class to define the experimental setup

#ifndef DetectorConstruction_h
#define DecectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

class G4Box;
class G4Material;
class G4Element;

class DetectorConstruction : public G4VUserDetectorConstruction
{
    
    public:
        DetectorConstruction();
        virtual ~DetectorConstruction();

    public:
        virtual G4VPhysicalVolume* Construct();
	
	void DefineMaterials();	
	G4VPhysicalVolume* DefineVolumes();

	// Geometry	
	G4Box* fSolidWorld;
	G4double fWorldSizeX, fWorldSizeY, fWorldSizeZ;
	G4LogicalVolume* fLogicWorld;
	G4PVPlacement* fPhysWorld;

	G4bool fCheckOverlaps;

	// Materials & Elements
	G4Material* fVacuum;
	G4Material* fBC400;
	G4Material* fLYSO;

	G4Element* fH;
	G4Element* fC;
	G4Element* fN;
	G4Element* fO;
	G4Element* fSie;
	G4Element* fY;
	G4Element* fCe;
	G4Element* fLu;

};

#endif
