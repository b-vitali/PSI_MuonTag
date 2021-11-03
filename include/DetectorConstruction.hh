/// \file  DetectorConstruction.hh
/// \brief Definition or DetectorConstruction class

#ifndef DetectorConstruction_h
#define DecectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "ScintSD.hh"
#include "PixelSD.hh"
#include "G4Cache.hh"

class G4VPhysicalVolume;
class G4VLogicallVolume;
class G4Box;
class G4GenericMessenger;
class G4Material;
class G4Element;
class G4MaterialPropertiesTable;
class DetectorMessenger;

/// Detector construction class to define geometry

class DetectorConstruction : public G4VUserDetectorConstruction{
    
    public:
        DetectorConstruction();
        virtual ~DetectorConstruction();

    public:
        virtual G4VPhysicalVolume* Construct();
        virtual void ConstructSDandField();

    public:
    	void SetCrystalSize(G4double);
    	void SetCrystalSize(G4ThreeVector);
	void SetCrystalMaterial(G4String);
	void SetSiPMmodel(G4String);
	void SetGround(G4double);
	void SetAngle(G4double);
	void SetAngleWithOpticalGrease(G4double);
	void SetTilt(G4double);

	G4int GetNbOfPixels(){return fNbOfPixelsX * fNbOfPixelsY;}

	G4String GetSiPMmodel(){return fmodel;}
	
	
    private:
        // methods	
	
	DetectorMessenger* fDetectorMessenger;
	G4String fmodel;
	void DefineMaterials();

	G4VPhysicalVolume* DefineVolumes();
	G4VPhysicalVolume* fCrysVolume;
	G4VPhysicalVolume* fSiPMVolume;
	G4VPhysicalVolume* fElementVolume;

	G4Box* fSolidCrys;
	G4Box* fSolidWorld;
	G4Box* fSolidElement;
	G4bool fCheckOverlaps;
	G4int fNbOfPixelsX, fNbOfPixelsY;

	G4LogicalVolume* fLogicPixel;
	G4LogicalVolume* fLogicCrys;
	
	G4double fCrysSizeX, fCrysSizeY, fCrysSizeZ;
	G4double fSiPM_sizeXY, fSiPM_sizeZ, fSiPM_windowZ;
	
	// Ruin values
	G4double fGround, fAngle, fAngleWithOpticalGrease, fTilt;
	
	// Materials & Elements
	G4Material* fBC400;
	G4Material* fLYSO;
	G4Material* fMaterial;
	G4Material* fMaterialWindow;
	G4Material* fAir;
	G4Material* fVacuum;
	G4Element* fN;
	G4Element* fO;
	G4Element* fC;
	G4Element* fH;
	G4Element* fSie;
	G4Element* fY;
	G4Element* fCe;
	G4Element* fLu;
	G4Material* fSi;
	G4Material* fSiResin;
	G4Material* fOG; // optical grease BC 631 index saint gobain
	G4Material* fEpResin;

		
	G4MaterialPropertiesTable* fBC400_mt;
	G4MaterialPropertiesTable* fLYSO_mt;
	G4Cache<ScintSD*> fScint_SD;
	G4Cache<PixelSD*> fPixel_SD;

};

#endif
