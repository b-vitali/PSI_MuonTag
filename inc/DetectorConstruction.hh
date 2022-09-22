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
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "VirtualDetectorSD.hh"
#include "ScintSD.hh"
#include "SiPMSD.hh"

#include "DetectorMessenger.hh"

#include "G4GlobalMagFieldMessenger.hh"

#include "G4Cache.hh"
#include "G4SDManager.hh"


class G4Box;
class G4VSolid;
class G4SubtractionSolid;
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
		G4VPhysicalVolume* DefineVolumes_guide();
		G4VPhysicalVolume* DefineVolumes_doubleguide();

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

		// Read dimensions
    	G4double fReadSizeX;
    	G4double fReadSizeY;
    	G4double fReadSizeZ;

		// Element dimensions
    	G4double fGreaseSizeX;
    	G4double fGreaseSizeY;
    	G4double fGreaseSizeZ;

		// Element dimensions
    	G4double fElementSizeX;
    	G4double fElementSizeY;
    	G4double fElementSizeZ;

		// Element dimensions
    	G4double fElementSizeX_guide;
    	G4double fElementSizeY_guide;
    	G4double fElementSizeZ_guide;

		// Guide dimensions (Excess to Scint)
    	G4double fGuideBorderX;
    	G4double fGuideBorderY;
    	G4double fGuideBorderZ;
		G4double fGuideHoleSizeXY;
		
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

		G4Box* fSolidScint2;
		G4LogicalVolume* fLogicScint2;
		G4PVPlacement* fPhysScint2;

		G4Box* fSolidScint_guide;
		G4LogicalVolume* fLogicScint_guide;
		G4PVPlacement* fPhysScint_guide;

		// Guide
		G4Box* fSolidGuide_tmp;
		G4VSolid* fSolidGuide;
		G4LogicalVolume* fLogicGuide;
		G4PVPlacement* fPhysGuide;

		// Guide-hole
		G4Box* fSolidGuide_hole;
		G4LogicalVolume* fLogicGuide_hole;
		G4PVPlacement* fPhysGuide_hole;

		// Readout
		G4Box* fSolidRead;
		G4LogicalVolume* fLogicRead;
		G4PVPlacement* fPhysRead;

		G4Box* fSolidGrease;
		G4LogicalVolume* fLogicGrease;
		G4PVPlacement* fPhysGrease;

		// SiPM
		G4Box* fSolidSiPM;
		G4LogicalVolume* fLogicSiPM;
		G4PVPlacement* fPhysSiPM;

		G4Box* fSolidRead_guide;
		G4LogicalVolume* fLogicRead_guide;
		G4PVPlacement* fPhysRead_guide;

		// Element
		G4Box* fSolidElement;
		G4LogicalVolume* fLogicElement;
		G4PVPlacement* fPhysElement;

		G4Box* fSolidElement_guide;
		G4LogicalVolume* fLogicElement_guide;
		G4PVPlacement* fPhysElement_guide;

		// VD
		G4Box* fSolidVD;
		G4LogicalVolume* fLogicVD;
		G4PVPlacement* fPhysVD;

		// second VD
		G4Box* fSolidVD_2;
		G4LogicalVolume* fLogicVD_2;
		G4PVPlacement* fPhysVD_2;

		G4bool fVDOn;
		G4bool fCheckOverlaps;

		//? Materials & Elements
		G4Material* fScintMaterial;
		G4Material* fSiPMMaterial;

		//-----------------------------------
		G4Material* fAir;
		G4Material* fVacuum;
		G4Material* fVacuum_nogamma;
		G4Material* fBC400;
		G4Material* fLYSO;
		G4Material* fOG; // optical grease BC 631 index saint gobain
		G4Material* fSi;
		G4Material* fSiResin;

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
