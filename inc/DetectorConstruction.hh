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

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "VirtualDetectorSD.hh"

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
	
		// This is needed for SensitiveDetectors and Fields
		virtual void ConstructSDandField();

	private:
		void DefineMaterials();	
		void DefineOpticalProperties();	

		G4VPhysicalVolume* DefineVolumes();

		//? Geometry	
		// World
		G4Box* fSolidWorld;
		G4LogicalVolume* fLogicWorld;
		G4PVPlacement* fPhysWorld;

		// Scint
		G4Box* fSolidScint;
		G4LogicalVolume* fLogicScint;
		G4PVPlacement* fPhysScint;

		// VD
		G4Box* fSolidVD;
		G4LogicalVolume* fLogicVD;
		G4PVPlacement* fPhysVD;

		G4bool fCheckOverlaps;

		//? Materials & Elements
		G4Material* fAir;
		G4Material* fVacuum;
		G4Material* fBC400;
		G4Material* fLYSO;

		G4MaterialPropertiesTable* fBC400_mt;
		G4MaterialPropertiesTable* fLYSO_mt;
  		G4bool fPhotonWorldPropagation;

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
