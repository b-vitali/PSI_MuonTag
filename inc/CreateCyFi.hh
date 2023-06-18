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
#include "G4IntersectionSolid.hh"

#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "VirtualDetectorSD.hh"
#include "ScintSD.hh"
#include "SiPMSD.hh"

#include "G4GlobalMagFieldMessenger.hh"

#include "G4Cache.hh"
#include "G4SDManager.hh"

#include "G4SubtractionSolid.hh"

class CreateCyFi
{
    public:
		//? Constructor
		CreateCyFi(double length, double radius, double thickness, double angle_in, double angle_out);
        virtual ~CreateCyFi();

		//? Main function
        void Create(G4LogicalVolume * hLogicWorld);
	
	private:
		//? All sections are separated for clarity
		void Materials(); 								// Required materials
		void OpticalProperties();						// Optical properties of the materials
		void Volumes(G4LogicalVolume * hLogicWorld);	// Geometry (needs the World Logic volume to place the CYFi)
		void SD();										// SensitiveDetectors
		
		//? Geometry	
		//-----------------------------------
		// CyFI variables
		G4double hCyFi_length;
		G4double hCyFi_radius;
		G4double hFiberThickness;
		G4double hangle_in;
		G4double hangle_out;

		//-----------------------------------
		// Read dimensions
    	G4double hReadSizeX;
    	G4double hReadSizeY;
    	G4double hReadSizeZ;

		// Read dimensions
    	G4double hSiPMSizeX;
    	G4double hSiPMSizeY;
    	G4double hSiPMSizeZ;

		//? Volume, Logical and Placement	
		//-----------------------------------
		// Readout(s)
		G4Box* hSolidRead_in;
		G4LogicalVolume* hLogicRead_in;
		G4PVPlacement* hPhysRead_in;

		G4Box* hSolidRead_out;
		G4LogicalVolume* hLogicRead_out;
		G4PVPlacement* hPhysRead_out;

		// Grease
		G4Box* hSolidGrease;
		G4LogicalVolume* hLogicGrease;
		G4PVPlacement* hPhysGrease;
		
		// SiPM(s)
		G4Box* hSolidSiPM_in;
		G4LogicalVolume* hLogicSiPM_in;
		G4PVPlacement* hPhysSiPM_in;		

		G4Box* hSolidSiPM_out;
		G4LogicalVolume* hLogicSiPM_out;
		G4PVPlacement* hPhysSiPM_out;

		G4bool hCheckOverlaps;

		//? Materials & Elements
		//-----------------------------------
		// Materials
		G4Material* hScintMaterial;
		G4Material* hSiPMMaterial;

		//-----------------------------------
		G4Material* hAir;
		G4Material* hVacuum;
		G4Material* hVacuum_nogamma;
		G4Material* hBC400;
		G4Material* hBC400_noscint;
		G4Material* hLYSO;
		G4Material* hOG; // optical grease BC 631 index saint gobain
		G4Material* hSi;
		G4Material* hSiResin;

		G4Material* hBCF10;
		G4Material* hBCF12;
		G4Material* hBCF20;
		G4Material* hFClad;
		G4Material* hSClad;

		//-----------------------------------
		// Elements
		G4Element* hH;
		G4Element* hC;
		G4Element* hN;
		G4Element* hO;
		G4Element* hSie;
		G4Element* hY;
		G4Element* hCe;
		G4Element* hLu;

		//-----------------------------------
		// Materials' Properties
		G4MaterialPropertiesTable* hBC400_mt;
		G4MaterialPropertiesTable* hLYSO_mt;

		G4MaterialPropertiesTable* hBCF10_mt;
		G4MaterialPropertiesTable* hBCF12_mt;
		G4MaterialPropertiesTable* hBCF20_mt;
};

#endif
