#include "YourDetectorConstruction.hh"

#include "YourDetectorMessenger.hh"

// for geometry definitions 
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

// for material definitions
#include "G4Material.hh"
#include "G4NistManager.hh"

// for having units and constants
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"


YourDetectorConstruction::YourDetectorConstruction()
:   G4VUserDetectorConstruction(), 
    fTargetMaterial(nullptr), 
    fTargetPhysicalVolume(nullptr) {
    // set default target material to be Silicon
    SetTargetMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    // set default thickness
    fTargetThickness = 0.01*CLHEP::mm; 
    fGunXPosition    = 0.0; // will be set properly automaticaly 
    // create detector messenger object (for your own detector UI commands)
    fDetMessenger    = new YourDetectorMessenger(this); 
}

YourDetectorConstruction::~YourDetectorConstruction() {
    delete fDetMessenger;
}


void YourDetectorConstruction::SetTargetMaterial(const G4String& matName) {
    // try to find the material in the NIST DB
    fTargetMaterial = G4NistManager::Instance()->FindOrBuildMaterial(matName);
    if (!fTargetMaterial) {
      G4cerr<< "\n **** ERROR in YourDetectorConstruction::SetTargetMaterial() \n"
            << "        Material with the given name of < " << matName << " >  \n"
            << "        was not found in the G4 NIST material DB               \n"
            << G4endl;
      exit(-1);      
    }
}


G4VPhysicalVolume* YourDetectorConstruction::Construct() {

    G4Material* materialWorld  = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");

    // Material for the target: material pointer stored in fTargetMaterial 
    G4Material* materialTarget = fTargetMaterial;
    //you could (and you will...) define your own material

    // Define target and world sizes
    G4double targetXSize  = fTargetThickness;
    G4double targetYZSize = 50*CLHEP::mm;
    G4double worldXSize   = 1.2*targetXSize;
    G4double worldYZSize  = 1.2*targetYZSize;
    fGunXPosition         = -0.25*( worldXSize + targetXSize );

    // Create the world and the target:
    // world
    G4Box*              worldSolid   = new G4Box("solid-World",  // name
    	                                       0.5*worldXSize,   // half x-size
    	                                       0.5*worldYZSize,  // half y-size 
    	                                       0.5*worldYZSize); // half z-size
    G4LogicalVolume*    worldLogical = new G4LogicalVolume(worldSolid,     // solid
                                                           materialWorld,  // material
                                                           "logic-World"); // name
    G4VPhysicalVolume*  worldPhyscal = new G4PVPlacement(nullptr,                 // (no) rotation
                                                         G4ThreeVector(0.,0.,0.), // translation
                                                         worldLogical,            // its logical volume
                                                         "World",                 // its name
                                                         nullptr,                 // its mother volume
                                                         false,                   // not used
                                                         0);                      // cpy number
    // target
    G4Box*              targetSolid   = new G4Box("solid-Target",    // name
    	                                          0.5*targetXSize,   // half x-size
    	                                          0.5*targetYZSize,  // half y-size 
    	                                          0.5*targetYZSize); // half z-size
    G4LogicalVolume*    targetLogical = new G4LogicalVolume(targetSolid,    // solid
                                                           materialTarget,  // material
                                                           "logic-Target"); // name
    G4VPhysicalVolume*  targetPhyscal = new G4PVPlacement(nullptr,                 // (no) rotation
                                                          G4ThreeVector(0.,0.,0.), // translation
                                                          targetLogical,           // its logical volume
                                                          "Target",                // its name
                                                          worldLogical,            // its mother volume
                                                          false,                   // not used
                                                          0);                      // cpy number
    fTargetPhysicalVolume = targetPhyscal;

    // RETURN WITH THE World PHYSICAL-VOLUME POINTER:
    return worldPhyscal;
}

