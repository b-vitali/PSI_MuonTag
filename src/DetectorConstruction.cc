/// \file  DetectorConstruction.hh
/// \brief Implementation of the class to define the experimental setup

#include "DetectorConstruction.hh"
#include <math.h>  
//#include "G4MagneticField.hh"

#include "G4AutoDelete.hh"


G4ThreadLocal 
G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = 0; 

//defined here in order to have the bool muEDM
G4double 	r;
G4int 		skip;		
G4double	theta; 	
G4int 		N;		


DetectorConstruction::DetectorConstruction()
{
	fVDOn = true;

	// Scintillator dimensions
    fScintTelescopeSizeX = 2*cm;
    fScintTelescopeSizeY = 1*cm;
    fScintTelescopeSizeZ = 20*cm;

    fScintGateSizeX = 5*cm;
    fScintGateSizeY = 5*cm;
    fScintGateSizeZ = 1*mm;

	// VirtualDetector dimensions
    fVDSizeX = 20*cm;
    fVDSizeY = 20*cm;
    fVDSizeZ = 10*mm;

	// World dimentions
    fWorldSizeX = 3*std::max({fScintTelescopeSizeX,fScintTelescopeSizeX,fVDSizeX});
    fWorldSizeY = 3*std::max({fScintTelescopeSizeY,fScintTelescopeSizeY,fVDSizeY});
    fWorldSizeZ = 3*std::max({fScintTelescopeSizeZ,fScintTelescopeSizeZ,fVDSizeZ});

	// At creation it calls for the function creating the materials 
	fDetectorMessenger = new DetectorMessenger(this);
	DefineMaterials();
	DefineOpticalProperties();
}

DetectorConstruction :: ~DetectorConstruction()
{
	delete fDetectorMessenger;
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
	// At construction the DefineVolumes describes the geometry
	return DefineVolumes();
}

void DetectorConstruction::DefineMaterials()
{
	G4double a; // atomic mass
	G4double z; // atomic number
	G4double density;

	//? You can take materials already defined 
	//G4NistManager* nist = G4NistManager::Instance();
	
	//? Or you can define your elements and materials
	//! Elements
	fH  = new G4Element( "H",  "H", z =  1., a =   1.01*g/mole);
	fC  = new G4Element( "C",  "C", z =  6., a =  12.01*g/mole);
	fN  = new G4Element( "N",  "N", z =  7., a =  14.01*g/mole);
	fO  = new G4Element( "O",  "O", z =  8., a =  16.00*g/mole);
	fSie= new G4Element("Si", "Si", z = 14., a = 28.0855*g/mole);
	fY  = new G4Element( "Y",  "Y", z = 39., a = 88.90585*g/mole);
	fCe = new G4Element("Ce", "Ce", z = 58., a = 140.116*g/mole);
	fLu = new G4Element("Lu", "Lu", z = 71., a = 174.967*g/mole);

	//! Materials
	// BC400
	fBC400 = new G4Material("BC400", density = 1.023*g/cm3, 2);
	fBC400->AddElement(fC, 1000);
	fBC400->AddElement(fH, 1103);

	// LYSO
	fLYSO = new G4Material("LYSO", density = 7.1*g/cm3, 5);
	fLYSO->AddElement(fLu,  9);
	fLYSO->AddElement( fY, 10);
	fLYSO->AddElement(fSie, 5);
	fLYSO->AddElement( fO, 25);
	fLYSO->AddElement(fCe,  5);

	// Vacuuum
	fVacuum = new G4Material("Vacuum",z=1.,a=1.01*g/mole, 
		density = universe_mean_density, 
		kStateGas, 0.1 * kelvin, 1.e-19 * pascal);

	// Air
	fAir = new G4Material("Air", density = 0.0010*g/cm3, 2);
	fAir->AddElement(fN, 70 * perCent);
	fAir->AddElement(fO, 30 * perCent);

	// Assign default materials
	fScintMaterial = fBC400;
}

void DetectorConstruction::DefineOpticalProperties()
{
	//? Material properties tables
	// ----------------------------------------------------------	
	//  BC400 optics
	// ----------------------------------------------------------

	std::vector<G4double> energy, scint;
	G4double tempe, tempscint;
	std::ifstream myfile;
	myfile.open("../tables/BC400_light_out.txt");
	while(true){
		myfile >> tempe >> tempscint;
		energy.push_back(1239.84197/tempe);
		scint.push_back(tempscint);
		
		if(myfile.eof()) break;
	}
	myfile.close();

	assert(energy.size() == scint.size());
	const G4int bc400 = int(energy.size());

	G4double* BC400_Energy = new G4double[bc400];
	G4double* BC400_SCINT = new G4double[bc400];

	G4double* BC400_RIND = new G4double[bc400];
	G4double* BC400_ABSL = new G4double[bc400];
	
	for(int i = 0; i < bc400; i++){
		BC400_Energy[i] = energy.at(i)*eV;
		BC400_SCINT[i] = scint.at(i);
		BC400_RIND[i] = 1.58;
		BC400_ABSL[i] = 160*cm;
	}
	
	energy.clear();
	scint.clear();
	
	
	assert(sizeof(BC400_SCINT) == sizeof(BC400_Energy));
	
	assert(sizeof(BC400_RIND) == sizeof(BC400_Energy));
	
	assert(sizeof(BC400_ABSL) == sizeof(BC400_Energy));

	fBC400_mt = new G4MaterialPropertiesTable();
	fBC400_mt->AddProperty(       "RINDEX", BC400_Energy,  BC400_RIND, bc400);
	fBC400_mt->AddProperty(    "ABSLENGTH", BC400_Energy,  BC400_ABSL, bc400);
	fBC400_mt->AddProperty("FASTCOMPONENT", BC400_Energy, BC400_SCINT, bc400);
	
	fBC400_mt->AddConstProperty("SCINTILLATIONYIELD",        11050./MeV);
	fBC400_mt->AddConstProperty(   "RESOLUTIONSCALE",                 1);
	fBC400_mt->AddConstProperty(  "SLOWTIMECONSTANT",            2.4*ns);
	fBC400_mt->AddConstProperty(  "SLOWSCINTILLATIONRISETIME",   0.9*ns);
	fBC400_mt->AddConstProperty(  "FASTTIMECONSTANT",            2.4*ns);
	fBC400_mt->AddConstProperty(  "FASTSCINTILLATIONRISETIME",   0.9*ns);
	fBC400_mt->AddConstProperty(        "YIELDRATIO",                 0.);
	
	fBC400->SetMaterialPropertiesTable(fBC400_mt);

	//  Set Birks Constant
	fBC400->GetIonisation()->SetBirksConstant(0.117645*mm/MeV);

	// ----------------------------------------------------------
	//  LYSO optics
	// ----------------------------------------------------------

	myfile.open("../tables/LYSO_light_out.txt");
	while(true){
		myfile >> tempe >> tempscint;
		energy.push_back(1239.84197/tempe);
		scint.push_back(tempscint);
		
		if(myfile.eof()) break;
	}
	myfile.close();

	assert(energy.size() == scint.size());
	const G4int lyso = int(energy.size());

	G4double* LYSO_Energy = new G4double[lyso];
	G4double* LYSO_SCINT  = new G4double[lyso];

	G4double* LYSO_RIND = new G4double[lyso];
	G4double* LYSO_ABSL = new G4double[lyso];
	
	for(int i = 0; i < lyso; i++){
		LYSO_Energy[i] = energy.at(i)*eV;
		LYSO_SCINT[i] = scint.at(i);
		LYSO_RIND[i] = 1.81;
		LYSO_ABSL[i] = 20*cm;
	}
	
	energy.clear();
	scint.clear();
	
	
	assert(sizeof(LYSO_SCINT) == sizeof(LYSO_Energy));
	
	assert(sizeof(LYSO_RIND) == sizeof(LYSO_Energy));
	
	assert(sizeof(LYSO_ABSL) == sizeof(LYSO_Energy));

	fLYSO_mt = new G4MaterialPropertiesTable();
	fLYSO_mt->AddProperty(       "RINDEX", LYSO_Energy,  LYSO_RIND, lyso);
	fLYSO_mt->AddProperty(    "ABSLENGTH", LYSO_Energy,  LYSO_ABSL, lyso);
	fLYSO_mt->AddProperty("FASTCOMPONENT", LYSO_Energy, LYSO_SCINT, lyso);
	
	fLYSO_mt->AddConstProperty("SCINTILLATIONYIELD",        33200./MeV);
	fLYSO_mt->AddConstProperty(   "RESOLUTIONSCALE",                 1);
	fLYSO_mt->AddConstProperty(  "SLOWTIMECONSTANT",            36*ns);
	fLYSO_mt->AddConstProperty(  "FASTTIMECONSTANT",            36*ns);
	fLYSO_mt->AddConstProperty(        "YIELDRATIO",                 0.);
	
	fLYSO->SetMaterialPropertiesTable(fLYSO_mt);

	//  Set Birks Constant
	//! fLYSO->GetIonisation()->SetBirksConstant(0.117645*mm/MeV);

	// ----------------------------------------------------------	
	// Vacuum & Air
	// ----------------------------------------------------------	
	fPhotonWorldPropagation = false;
	if(fPhotonWorldPropagation){
		G4double vacuum_Energy[] = {1.5*eV, 4.*eV};
		const G4int vacnum = sizeof(vacuum_Energy) / sizeof(G4double);

		G4double vRIND = 1.;
		G4double vacuum_RIND[] = {vRIND, vRIND};
		assert(sizeof(vacuum_RIND) == sizeof(vacuum_Energy));

		G4MaterialPropertiesTable* vacuum_mt = new G4MaterialPropertiesTable();
		vacuum_mt->AddProperty("RINDEX", vacuum_Energy, vacuum_RIND, vacnum);
		fVacuum->SetMaterialPropertiesTable(vacuum_mt);
		fAir   ->SetMaterialPropertiesTable(vacuum_mt);
	}
}


G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
	
	// World Solid (size) -> Logical (material) -> PVPLacement (posiz, rotaz, to interact)
	fSolidWorld	= new G4Box("World", 0.5*fWorldSizeX, 0.5*fWorldSizeY, 0.5*fWorldSizeZ);
    fLogicWorld = new G4LogicalVolume(fSolidWorld, fVacuum, "World");
    fPhysWorld	= new G4PVPlacement(0, G4ThreeVector(0,0,0), fLogicWorld, "World", 0, false, 0, fCheckOverlaps);

	// Scintillator Solid (size) -> Logical (material) -> PVPLacement (posiz, rotaz, to interact)
	fSolidScint_telescope	= new G4Box("fSolidScint_telescope", 0.5*fScintTelescopeSizeX, 0.5*fScintTelescopeSizeY, 0.5*fScintTelescopeSizeZ);
    fLogicScint_telescope = new G4LogicalVolume(fSolidScint_telescope, fScintMaterial, "fLogicScint_telescope");

	G4double shift=5*cm;

   	for(int i=0; i<4; i=i+1)
	{
		G4Rotate3D 	rotation =  G4Rotate3D(i*90*deg, G4ThreeVector(0, 0, 1)); //i*theta*deg std::cos(theta*i)
		G4Translate3D 	translate =  G4Translate3D(G4ThreeVector(fScintTelescopeSizeY/2, fScintTelescopeSizeX/2, fScintTelescopeSizeZ/2));
		G4Transform3D transform = G4Translate3D(0,0,shift)*rotation*translate;
		fPhysScint_telescope	= new G4PVPlacement(transform, fLogicScint_telescope, "Scint_telescope", fLogicWorld, false, i, fCheckOverlaps);
	}

	fSolidScint_gate 	= new G4Box("fSolidScint_gate", 0.5*fScintGateSizeX, 0.5*fScintGateSizeY, 0.5*fScintGateSizeZ);
    fLogicScint_gate 	= new G4LogicalVolume(fSolidScint_gate, fScintMaterial, "fLogicScint_gate");
	fPhysScint_gate	= new G4PVPlacement(0, G4ThreeVector(0,0,shift-2*cm), fLogicScint_gate, "Scint_gate", fLogicWorld, false, 4, fCheckOverlaps);


	// Set how the volumes are visualized
    fLogicWorld	->SetVisAttributes(G4Colour(1, 1, 1, 0.1));
    fLogicScint_telescope	->SetVisAttributes(G4Colour(0.34, 0.57, 0.8, 0.5));
	fLogicScint_gate->SetVisAttributes(G4Colour(0.34, 0.57, 0.8, 0.5));
	
	if(fVDOn){
		// VirtualDetector Solid (size) -> Logical (material) -> PVPLacement (posiz, rotaz, to interact)
		fSolidVD	= new G4Box("VD", 0.5*fVDSizeX, 0.5*fVDSizeY, 0.5*fVDSizeZ);
    	fLogicVD 	= new G4LogicalVolume(fSolidVD, fVacuum, "VD");
    	fPhysVD		= new G4PVPlacement(0, G4ThreeVector(0., 0., shift-fVDSizeZ*0.5), fLogicVD, "VD", fLogicWorld, false, 0, fCheckOverlaps);

		// 2_VirtualDetector Solid (size) -> Logical (material) -> PVPLacement (posiz, rotaz, to interact)
		fSolidVD_2	= new G4Box("VD", 0.5*fVDSizeX, 0.5*fVDSizeY, 0.5*fVDSizeZ);
    	fLogicVD_2 	= new G4LogicalVolume(fSolidVD_2, fVacuum, "VD");
    	fPhysVD_2	= new G4PVPlacement(0, G4ThreeVector(0., 0., shift+fVDSizeZ*0.5+fScintTelescopeSizeZ), fLogicVD_2, "VD", fLogicWorld, false, 1, fCheckOverlaps);
    	fLogicVD	->SetVisAttributes(G4Colour(0.8, 0.34, 0.68, 0.2));
    	fLogicVD_2	->SetVisAttributes(G4Colour(0.8, 0.34, 0.68, 0.2));
	}
	
    return fPhysWorld;
}

void DetectorConstruction::ConstructSDandField()
{

	auto sdManager = G4SDManager::GetSDMpointer();
   	G4String SDname;

	ScintSD* scint_SD_gate = new ScintSD(SDname="Scint_gate");
  	sdManager->AddNewDetector(scint_SD_gate);
	fLogicScint_gate->SetSensitiveDetector(scint_SD_gate);

	ScintSD* scint_SD_telescope = new ScintSD(SDname="Scint_telescope");
  	sdManager->AddNewDetector(scint_SD_telescope);
  	fLogicScint_telescope->SetSensitiveDetector(scint_SD_telescope);

	if(fVDOn){
		// Create the Sensitive Detector defined in VirtualDetectorSD 
		VirtualDetectorSD * VD_SD = new VirtualDetectorSD("VirtualDetector");

		// Assign the SD to the logial volume
		fLogicVD->SetSensitiveDetector(VD_SD);

		VirtualDetectorSD * VD_SD_2 = new VirtualDetectorSD("VirtualDetector2");
		fLogicVD_2->SetSensitiveDetector(VD_SD_2);
	}
}

/*
	From here functions used through the Messenger to modify the detector
*/

void DetectorConstruction :: SetScintSize(G4double size){
	fScintGateSizeX = size;
	fScintGateSizeY = size;
	fScintGateSizeZ = size;

	fSolidScint_gate->SetXHalfLength(size*0.5);
	fSolidScint_gate->SetYHalfLength(size*0.5);
	fSolidScint_gate->SetZHalfLength(size*0.5);

	//fPhysScint->SetTranslation(G4ThreeVector(0, 0, 20*cm));
	G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

void DetectorConstruction :: SetScintSize(G4ThreeVector size){
	fScintGateSizeX = size.getX();
	fScintGateSizeY = size.getY();
	fScintGateSizeZ = size.getZ();

	fSolidScint_gate->SetXHalfLength(size.getX()*0.5);
	fSolidScint_gate->SetYHalfLength(size.getY()*0.5);
	fSolidScint_gate->SetZHalfLength(size.getZ()*0.5);

	//fPhysScint->SetTranslation(G4ThreeVector(0, 0, 20*cm));
	G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

void DetectorConstruction::SetScintMaterial(G4String name){
	if(name == "BC400") fScintMaterial = fBC400;
	else if(name == "LYSO") fScintMaterial = fLYSO;
	G4RunManager::GetRunManager()->ReinitializeGeometry();
}