/// \file  DetectorConstruction.hh
/// \brief Definition of the class to define the experimental setup

#ifndef CreateCyFi_h
#define CreateCyFi_h 1

#include "G4TessellatedSolid.hh"
#include "G4TriangularFacet.hh"
// My additional file to create the G4TessellatedSolid
#include "CreateHelix.hh"
using namespace HelixMaker;

//!
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
//!

class G4Box;
class G4Material;
class G4Element;

class CreateCyFi
{
    
    public:
		// Constructor
        CreateCyFi();
        virtual ~CreateCyFi();

	private:
		// Main function
        void Create(G4LogicalVolume * fLogicWorld);
	
		// Material and Optical properties separated for clarity
		void Materials();	
		void OpticalProperties();	

		void Volumes(G4LogicalVolume * fLogicWorld);
		
		// This is needed for SensitiveDetectors and Fields
		void SD(G4SDManager * sdManager);
		
		//? Geometry	
		//-----------------------------------
		// Read dimensions
    	G4double fReadSizeX;
    	G4double fReadSizeY;
    	G4double fReadSizeZ;

		// Read dimensions
    	G4double fSiPMSizeX;
    	G4double fSiPMSizeY;
    	G4double fSiPMSizeZ;

		//-----------------------------------
		// Readout(s)
		G4Box* fSolidRead_in;
		G4LogicalVolume* fLogicRead_in;
		G4PVPlacement* fPhysRead_in;

		G4Box* fSolidRead_out;
		G4LogicalVolume* fLogicRead_out;
		G4PVPlacement* fPhysRead_out;

		// Grease
		G4Box* fSolidGrease;
		G4LogicalVolume* fLogicGrease;
		G4PVPlacement* fPhysGrease;
		
		// SiPM(s)
		G4Box* fSolidSiPM_in;
		G4LogicalVolume* fLogicSiPM_in;
		G4PVPlacement* fPhysSiPM_in;		

		G4Box* fSolidSiPM_out;
		G4LogicalVolume* fLogicSiPM_out;
		G4PVPlacement* fPhysSiPM_out;

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

		G4MaterialPropertiesTable* fBCF10_mt;
		G4MaterialPropertiesTable* fBCF12_mt;
		G4MaterialPropertiesTable* fBCF20_mt;
};

#endif
