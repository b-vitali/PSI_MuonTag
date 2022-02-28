/// \file  DetectorConstruction.hh
/// \brief Definition of the class to define the experimental setup

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

#include "DetectorMessenger.hh"

#include "G4GlobalMagFieldMessenger.hh"

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

	public:
		// Functions needed to use the Messenger
		void SetScintSize(G4double);
    	void SetScintSize(G4ThreeVector);
		void SetScintMaterial(G4String);

	private:
	
		// Material and Optical properties separated for clarity
		void DefineMaterials();	
		void DefineOpticalProperties();	

		G4VPhysicalVolume* DefineVolumes();

		DetectorMessenger* fDetectorMessenger;

    static G4ThreadLocal G4GlobalMagFieldMessenger*  fMagFieldMessenger;

		//? Geometry	
		//-----------------------------------
		// World dimentions
    	G4double fWorldSizeX;
    	G4double fWorldSizeY;
    	G4double fWorldSizeZ;

		// Scintillator dimensions
    	G4double fScintSizeX;
    	G4double fScintSizeY;
    	G4double fScintSizeZ;

		// VirtualDetector dimensions
    	G4double fVDSizeX;
    	G4double fVDSizeY;
    	G4double fVDSizeZ;

		//-----------------------------------
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

		// second VD
		G4Box* fSolidVD_2;
		G4LogicalVolume* fLogicVD_2;
		G4PVPlacement* fPhysVD_2;

		G4bool fVDOn;
		G4bool fmuEDM;
		G4bool fCheckOverlaps;

		//? Materials & Elements
		G4Material* fScintMaterial;

		//-----------------------------------
		G4Material* fAir;
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

		//-----------------------------------
		G4MaterialPropertiesTable* fBC400_mt;
		G4MaterialPropertiesTable* fLYSO_mt;
  		G4bool fPhotonWorldPropagation;

};

#endif
