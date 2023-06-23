/// \file  DetectorConstruction.hh
/// \brief Definition of the class to define the experimental setup

#ifndef CreateCyFi_h
#define CreateCyFi_h 1

/*
HOW TO INCLUDE IT
You need 8 files (4 .hh & 4 .cc):
	- SiPMSD.hh/cc
	- SiPMHit.hh/cc
	- CreateHelix.hh/cc
	- CreateCyFi.hh/cc

In EventAction.hh:
	#include "CreateCyFi.hh"

In EventAction::EndOfEventAction:
	CreateCyFi * tmp_CyFi = new CreateCyFi();
	tmp_CyFi->EndOfEvent(HCE, event);

In DetectorConstruction_hh:
	#include "CreateCyFi.hh"

In DetectorConstruction::DefineVolumes:
	CreateCyFi * CyFi = new CreateCyFi(15*cm, 3.5*cm, 1*mm, 55*deg, 0*deg);
	CyFi->Create(fLogicWorld);
*/

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

#include "G4HCofThisEvent.hh"
#include "G4Event.hh"

class CreateCyFi
{
    public:
		//? Constructor
		CreateCyFi(double length, double radius, double thickness, double angle_in, double angle_out);
		CreateCyFi();
        virtual ~CreateCyFi();

		//? Main function
        void Create(G4LogicalVolume * hLogicWorld);					// This calls the private functions
	
		//? Function to write the data out
		void CreateNTuples();										// Create NTuple (to be called in EventAction)
		void FillNTuples(G4HCofThisEvent*HCE, const G4Event* event);	// Write to file (to be called in EndOfEventAction)

	private:
		//? All sections are separated for clarity
		void Materials(); 										// Required materials
		void OpticalProperties();								// Optical properties of the materials
		void Volumes();											// Geometry 
		void SD();												// SensitiveDetectors

		//? Pointer to the universe
		G4LogicalVolume * hLogicWorld;

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
		G4Material* hMaterial;

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

		//? For the EndOfEvent
		//-----------------------------------
		// Hit collection number
		G4int fCollIDSiPM_out;
		G4int fCollIDSiPM_in;

		//-----------------------------------
		// temporary SD to create the ntuple
		SiPMSD * tmp_sipm_out;
		SiPMSD * tmp_sipm_in;

		// To track if the ntuple were created
 		bool hNTUpleCreated;
		
		// to register each event just once
		G4int fEvID; 
};

#endif
