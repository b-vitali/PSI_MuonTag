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
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4TwistedBox.hh"
#include "G4IntersectionSolid.hh"

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

#include "G4SubtractionSolid.hh"

#include "CreateCyFi.hh"

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
		G4VPhysicalVolume* DefineVolumes_SciFi();
		G4VPhysicalVolume* DefineVolumes_MuEDM();
		G4VPhysicalVolume* DefineVolumes_CyFi();
		G4VPhysicalVolume* DefineVolumes_G4CyFi();
		
		DetectorMessenger* fDetectorMessenger;

    static G4ThreadLocal G4GlobalMagFieldMessenger*  fMagFieldMessenger;

		//? Geometry	
		//-----------------------------------
		// World dimentions
    	G4double fWorldSizeX;
    	G4double fWorldSizeY;
    	G4double fWorldSizeZ;

		// Scintillator dimensions out
    	G4double fScintSizeX_out;
    	G4double fScintSizeY_out;
    	G4double fScintSizeZ_out;

		// Scintillator dimensions in
    	G4double fScintSizeX_in;
    	G4double fScintSizeY_in;
    	G4double fScintSizeZ_in;

		// Read dimensions
    	G4double fReadSizeX;
    	G4double fReadSizeY;
    	G4double fReadSizeZ;

		// Read dimensions
    	G4double fSiPMSizeX;
    	G4double fSiPMSizeY;
    	G4double fSiPMSizeZ;

		// Element dimensions out
    	G4double fElementSizeX_out;
    	G4double fElementSizeY_out;
    	G4double fElementSizeZ_out;

		// Element dimensions in
    	G4double fElementSizeX_in;
    	G4double fElementSizeY_in;
    	G4double fElementSizeZ_in;

		// VirtualDetector dimensions
    	G4double fVDSizeX;
    	G4double fVDSizeY;
    	G4double fVDSizeZ;

		//-----------------------------------
		// World
		G4Box* fSolidWorld;
		G4LogicalVolume* fLogicWorld;
		G4PVPlacement* fPhysWorld;

		// Scint Out
		G4Box* fSolidScint_out;
		G4LogicalVolume* fLogicScint_out;
		G4PVPlacement* fPhysScint_out;
		
		// Scint In
		G4Box* fSolidScint_in;
		G4LogicalVolume* fLogicScint_in;
		G4PVPlacement* fPhysScint_in;

		// Readout
		G4Box* fSolidRead_in;
		G4LogicalVolume* fLogicRead_in;
		G4PVPlacement* fPhysRead_in;

		// Readout
		G4Box* fSolidRead_out;
		G4LogicalVolume* fLogicRead_out;
		G4PVPlacement* fPhysRead_out;

		// Grease
		G4Box* fSolidGrease;
		G4LogicalVolume* fLogicGrease;
		G4PVPlacement* fPhysGrease;
		
		// SiPM
		G4Box* fSolidSiPM_in;
		G4LogicalVolume* fLogicSiPM_in;
		G4PVPlacement* fPhysSiPM_in;		

		// SiPM
		G4Box* fSolidSiPM_out;
		G4LogicalVolume* fLogicSiPM_out;
		G4PVPlacement* fPhysSiPM_out;

		// Element out
		G4Box* fSolidElement_out;
		G4LogicalVolume* fLogicElement_out;
		G4PVPlacement* fPhysElement_out;

		// Element in
		G4Box* fSolidElement_in;
		G4LogicalVolume* fLogicElement_in;
		G4PVPlacement* fPhysElement_in;

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
		G4Material* fBC400_noscint;
		G4Material* fLYSO;
		G4Material* fOG; // optical grease BC 631 index saint gobain
		G4Material* fSi;
		G4Material* fSiResin;

	G4Material* fBCF10;
	G4Material* fBCF12;
	G4Material* fBCF20;
	G4Material* fFClad;
	G4Material* fSClad;

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

	G4MaterialPropertiesTable* fBCF10_mt;
	G4MaterialPropertiesTable* fBCF12_mt;
	G4MaterialPropertiesTable* fBCF20_mt;
	
	G4LogicalVolume* fLogicFiber;
	G4LogicalVolume* fLogicFiberIn;

};

#endif
