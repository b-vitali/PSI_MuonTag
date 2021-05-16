#ifndef YOURDETECTORCONTRUCTION_HH
#define YOURDETECTORCONTRUCTION_HH

#include "G4VUserDetectorConstruction.hh"

class YourDetectorMessenger;
class G4Material;
class G4String;
class G4LogicalVolume;

class YourDetectorConstruction : public G4VUserDetectorConstruction {

  public:
    // CTR & DTR 
    YourDetectorConstruction();
    virtual ~YourDetectorConstruction();

    virtual G4VPhysicalVolume* Construct();

    void     SetTargetThickness(const G4double thick) { fTargetThickness = thick; }
    G4double GetTargetThickness() const               { return fTargetThickness;  }

    // Public methods to set/get the target material: G4 NIST materials
    void              SetTargetMaterial(const G4String& matName);
    const G4Material* GetTargetMaterial() const  { return fTargetMaterial; }

    // Public method to obtain the proper gun-position depending on the detector
    G4double GetGunXPosition() { return fGunXPosition; }

    // Public method to get the target logical volume pointer (used for scoring)
    const G4VPhysicalVolume* GetTargetPhysicalVolume() const { 
        return fTargetPhysicalVolume; 
    } 
    
  // Data member declaration
  private:

    // The detector messenger pointer: to set the target thickness
    YourDetectorMessenger* fDetMessenger;

    // Target material 
    G4Material*            fTargetMaterial;

    // Target logical volume pointer
    G4VPhysicalVolume*     fTargetPhysicalVolume;

    // The target thickness i.e. its (full) x-size (YZ size will be proportional)
    G4double               fTargetThickness;
    
    // The midpoint between the target and the world on the negative x-side
    G4double               fGunXPosition;
};

#endif   