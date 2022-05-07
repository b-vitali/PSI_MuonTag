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

#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "VirtualDetectorSD.hh"
#include "ScintSD.hh"

#include "DetectorMessenger.hh"

#include "G4GlobalMagFieldMessenger.hh"

#include "G4Cache.hh"
#include "G4SDManager.hh"


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

		// Scintillator Gate dimensions
    	G4double fScintSizeX_gate;
    	G4double fScintSizeY_gate;
    	G4double fScintSizeZ_gate;

		// Scintillator Telescope dimensions
    	G4double fScintSizeX_telescope;
    	G4double fScintSizeY_telescope;
    	G4double fScintSizeZ_telescope;

		G4double fScintSizeX_telescope2;
    	G4double fScintSizeY_telescope2;
    	G4double fScintSizeZ_telescope2;

		// Reads dimensions
    	G4double fReadSizeX_gate;
    	G4double fReadSizeY_gate;
    	G4double fReadSizeZ_gate;
    	
		G4double fReadSizeX_telescope;
    	G4double fReadSizeY_telescope;
    	G4double fReadSizeZ_telescope;
    	
		G4double fReadSizeX_telescope2;
    	G4double fReadSizeY_telescope2;
    	G4double fReadSizeZ_telescope2;

		// Elements dimensions
    	G4double fElementSizeX_gate;
    	G4double fElementSizeY_gate;
    	G4double fElementSizeZ_gate;
    	
		G4double fElementSizeX_telescope;
    	G4double fElementSizeY_telescope;
    	G4double fElementSizeZ_telescope;
    	
		G4double fElementSizeX_telescope2;
    	G4double fElementSizeY_telescope2;
    	G4double fElementSizeZ_telescope2;

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
		G4Box* fSolidScint_gate;
		G4LogicalVolume* fLogicScint_gate;
		G4PVPlacement* fPhysScint_gate;

		G4Box* fSolidScint_telescope;
		G4LogicalVolume* fLogicScint_telescope;
		G4PVPlacement* fPhysScint_telescope;

		G4Box* fSolidScint_telescope2;
		G4LogicalVolume* fLogicScint_telescope2;
		G4PVPlacement* fPhysScint_telescope2;

		// Readouts
		G4Box* fSolidRead_gate;
		G4LogicalVolume* fLogicRead_gate;
		G4PVPlacement* fPhysRead_gate;

		G4Box* fSolidRead_telescope;
		G4LogicalVolume* fLogicRead_telescope;
		G4PVPlacement* fPhysRead_telescope;
		
		G4Box* fSolidRead_telescope2;
		G4LogicalVolume* fLogicRead_telescope2;
		G4PVPlacement* fPhysRead_telescope2;
		
		// Elements
		G4Box* fSolidElement_gate;
		G4LogicalVolume* fLogicElement_gate;
		G4PVPlacement* fPhysElement_gate;

		G4Box* fSolidElement_telescope;
		G4LogicalVolume* fLogicElement_telescope;
		G4PVPlacement* fPhysElement_telescope;

		G4Box* fSolidElement_telescope2;
		G4LogicalVolume* fLogicElement_telescope2;
		G4PVPlacement* fPhysElement_telescope2;

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
		G4Material* fOG; // optical grease BC 631 index saint gobain

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
