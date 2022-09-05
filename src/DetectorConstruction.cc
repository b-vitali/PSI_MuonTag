/// \file  DetectorConstruction.hh
/// \brief Implementation of the class to define the experimental setup

#include "DetectorConstruction.hh"
#include <math.h>  
//#include "G4MagneticField.hh"

#include "G4AutoDelete.hh"


G4ThreadLocal 
G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = 0; 

DetectorConstruction::DetectorConstruction()
{
	//fVDOn = false;

	// Scintillator dimensions
    fScintSizeX = 20*mm;
    fScintSizeY = 20*mm;
    fScintSizeZ = 0.2*mm;

	// VirtualDetector dimensions
    fVDSizeX = 20*cm;
    fVDSizeY = 20*cm;
    fVDSizeZ = 10*mm;

	// World dimentions
    fWorldSizeX = 3*std::max(fScintSizeX,fVDSizeX);
    fWorldSizeY = 3*std::max(fScintSizeY,fVDSizeY);
    fWorldSizeZ = 2*20*std::max(fScintSizeZ,fVDSizeZ);
	
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
	fVacuum_nogamma = new G4Material("Vacuum_nogamma",z=1.,a=1.01*g/mole, 
		density = universe_mean_density, 
		kStateGas, 0.1 * kelvin, 1.e-19 * pascal);
	// Air
	fAir = new G4Material("Air", density = 0.0010*g/cm3, 2);
	fAir->AddElement(fN, 70 * perCent);
	fAir->AddElement(fO, 30 * perCent);

	// Optical grease
	fOG = new G4Material("OpticalGrease",z=1.,a=1.01*g/mole, 
		     		 density = universe_mean_density, kStateGas,
				 0.1 * kelvin, 1.e-19 * pascal);

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
	fPhotonWorldPropagation = true;
	// if(fPhotonWorldPropagation){
		G4double vacuum_Energy[] = {1.5*eV, 4.*eV};
		const G4int vacnum = sizeof(vacuum_Energy) / sizeof(G4double);

		G4double vRIND = 1.;
		G4double vacuum_RIND[] = {vRIND, vRIND};
		assert(sizeof(vacuum_RIND) == sizeof(vacuum_Energy));

		G4MaterialPropertiesTable* vacuum_mt = new G4MaterialPropertiesTable();
		vacuum_mt->AddProperty("RINDEX", vacuum_Energy, vacuum_RIND, vacnum);
		fVacuum->SetMaterialPropertiesTable(vacuum_mt);
		fAir   ->SetMaterialPropertiesTable(vacuum_mt);
	// }

	// Silicium
	G4double Si_Energy[] = {.5*eV, 9.*eV};
	const G4int Sinum = sizeof(vacuum_Energy) / sizeof(G4double);

	// ----------------------------------------------------------	
	// Optical grease
	// ----------------------------------------------------------	
	G4double OG_RIND[] = {1.465, 1.465};
	assert(sizeof(OG_RIND) == sizeof(Si_Energy));
		
	G4MaterialPropertiesTable* OG_mt = new G4MaterialPropertiesTable();
	OG_mt->AddProperty("RINDEX", Si_Energy, OG_RIND, Sinum);
	fOG->SetMaterialPropertiesTable(OG_mt);
}


G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
	fCheckOverlaps = true;

	// World Solid (size) -> Logical (material) -> PVPLacement (posiz, rotaz, to interact)
	fSolidWorld	= new G4Box("World", 0.5*fWorldSizeX, 0.5*fWorldSizeY, 0.5*fWorldSizeZ);
    fLogicWorld = new G4LogicalVolume(fSolidWorld, fVacuum_nogamma, "World");
    fPhysWorld	= new G4PVPlacement(0, G4ThreeVector(), fLogicWorld, "World", 0, false, 0, fCheckOverlaps);

	// Scintillator Solid (size) -> Logical (material) -> PVPLacement (posiz, rotaz, to interact)
	fSolidScint	= new G4Box("Scint", 0.5*fScintSizeX, 0.5*fScintSizeY, 0.5*fScintSizeZ);
    fLogicScint = new G4LogicalVolume(fSolidScint, fScintMaterial, "Scint");

	// Element
	fElementSizeX = fScintSizeX + 1*mm;
	fElementSizeY = fScintSizeY + 1*mm;
	fElementSizeZ = fScintSizeZ + 1*mm;

	fSolidElement = new G4Box("Element", 0.5*(fElementSizeX),0.5*(fElementSizeY), 0.5*(fElementSizeZ));
    fLogicElement = new G4LogicalVolume(fSolidElement, fVacuum, "Element");
	fLogicElement->SetVisAttributes(G4Colour(0, 1, 0, 0.1));

	G4Rotate3D 	rotation =  G4Rotate3D(30*deg, G4ThreeVector(0, 1, 0)); //i*theta*deg std::cos(theta*i)
	G4Translate3D 	translate =  G4Translate3D(G4ThreeVector(0., 0., 3*cm));
	G4Transform3D transform = translate*rotation;

	// Read Solid and Phys. 
	fReadSizeX = 0.4*mm;
	fReadSizeY = 3*mm;
	fReadSizeZ = fScintSizeZ;
	
	fSolidRead	= new G4Box("Read", 0.5*fReadSizeX,0.5*fReadSizeY, 0.5*fReadSizeZ);
    fLogicRead = new G4LogicalVolume(fSolidRead, fVacuum, "Read");
	fLogicRead->SetVisAttributes(G4Colour(0.0, 1.0, 1.0, 0.3));

	// grease Solid and Phys. 
	fSolidGrease = new G4Box("Read", 0.5*fReadSizeX*0.5,0.5*fReadSizeY, 0.5*fReadSizeZ);
    fLogicGrease = new G4LogicalVolume(fSolidGrease, fOG, "Grease");
	fLogicGrease->SetVisAttributes(G4Colour(1,0,0, 0.5));

	// VirtualDetector/SiPM Solid and Phys. 
	fSolidSiPM	= new G4Box("SiPM", 0.5*fReadSizeX*0.5, 0.5*fReadSizeY, 0.5*fReadSizeZ);
    fLogicSiPM 	= new G4LogicalVolume(fSolidSiPM, fOG, "SiPM");
    fLogicSiPM	->SetVisAttributes(G4Colour(0.8, 0.34, 0.68, 0.5));

	// Put Grease and SiPM in Read
    G4ThreeVector Grease_pos = G4ThreeVector(-(0.5*fReadSizeX*0.5), 0, 0); 
    G4ThreeVector SiPM_pos = G4ThreeVector(0.5*fReadSizeX*0.5, 0,0); //0.5*fReadSizeX
		
	//fPhysRead 	= new G4PVPlacement(0, G4ThreeVector(0, 0, 0), fLogicRead, "Read", fLogicWorld, true, 0, fCheckOverlaps);
	fPhysGrease	= new G4PVPlacement(0, Grease_pos, fLogicGrease, "Grease", fLogicRead, false, fCheckOverlaps);
	fPhysSiPM		= new G4PVPlacement(0, SiPM_pos, fLogicSiPM, "SiPM", fLogicRead, false, fCheckOverlaps);

	// Put Scint and Read in element
	G4ThreeVector Scint_pos = G4ThreeVector(0, 0, 0); 
   	G4ThreeVector Read_pos1 = G4ThreeVector(0.5*fReadSizeX+0.5*fScintSizeX, 0.2*fScintSizeY,0); //0.5*fReadSizeX
   	G4ThreeVector Read_pos2 = G4ThreeVector(0.5*fReadSizeX+0.5*fScintSizeX, -0*fScintSizeY,0); //0.5*fReadSizeX
   	G4ThreeVector Read_pos3 = G4ThreeVector(0.5*fReadSizeX+0.5*fScintSizeX, -0.2*fScintSizeY,0); //0.5*fReadSizeX

		fSolidScint2	= new G4Box("Scint2", 0.7*fScintSizeX, 0.7*fScintSizeY, 0.7*fScintSizeZ);
    	fLogicScint2 = new G4LogicalVolume(fSolidScint2, fScintMaterial, "Scint2");

		// Put Scint and Read in element
		fPhysElement 	= new G4PVPlacement(0, G4ThreeVector(0, 0, 0), fLogicElement, "Element", fLogicWorld, true, 0, fCheckOverlaps);
		fPhysRead		= new G4PVPlacement(0, Read_pos1, fLogicRead, "Read0", fLogicElement, true, 0, fCheckOverlaps);
		fPhysRead		= new G4PVPlacement(0, Read_pos2, fLogicRead, "Read1", fLogicElement, true, 1, fCheckOverlaps);
		fPhysRead		= new G4PVPlacement(0, Read_pos3, fLogicRead, "Read2", fLogicElement, true, 2, fCheckOverlaps);
		fPhysScint		= new G4PVPlacement(0, Scint_pos, fLogicScint, "Scint", fLogicElement, false, fCheckOverlaps);

	
		//transform = translate; //*rotation
		//fPhysElement 	= new G4PVPlacement(transform , fLogicElement, "Element", fLogicWorld, true, 1, fCheckOverlaps);
		//transform = translate*translate; //rotation*rotation
		//fPhysElement 	= new G4PVPlacement(transform , fLogicElement, "Element", fLogicWorld, true, 2, fCheckOverlaps);
		//transform = translate*translate*translate; //rotation*rotation
		//fPhysElement 	= new G4PVPlacement(transform , fLogicScint2, "Scint2", fLogicWorld, false, fCheckOverlaps);

	if(fVDOn){
		// 2_VirtualDetector Solid (size) -> Logical (material) -> PVPLacement (posiz, rotaz, to interact)
		fSolidVD_2	= new G4Box("VD", 0.5*fVDSizeX, 0.5*fVDSizeY, 0.5*fVDSizeZ);
    	fLogicVD_2 	= new G4LogicalVolume(fSolidVD_2, fVacuum, "VD");
    	fPhysVD_2	= new G4PVPlacement(0, G4ThreeVector(0., 0., 1*cm), fLogicVD_2, "VD", fLogicWorld, false, 1, fCheckOverlaps);
    	fLogicVD	->SetVisAttributes(G4Colour(0.8, 0.34, 0.68, 0.05));
    	fLogicVD_2	->SetVisAttributes(G4Colour(0.8, 0.34, 0.68, 0.05));
	}

		fPhysElement 	= new G4PVPlacement(0, G4ThreeVector(0, 0, 1*cm), fLogicElement, "Element1", fLogicWorld, true, 1, fCheckOverlaps);


/*

	// Scintillator with light guide
	fSolidScint_guide	= new G4Box("Scint_guide", 0.5*fScintSizeX, 0.5*fScintSizeY, 0.5*fScintSizeZ);
    fLogicScint_guide = new G4LogicalVolume(fSolidScint, fScintMaterial, "Scint");
	
	// Element for scintillator with light guide
	// fElementSizeX = fScintSizeX + 5*mm;
	// fElementSizeY = fScintSizeY + 5*mm;
	// fElementSizeZ = fScintSizeZ + 5*mm;

	// Guide
	fGuideBorderX=5*mm;
    fGuideBorderY=5*mm;
    fGuideBorderZ=2.5*mm;
	fGuideHoleSizeXY=15*mm;

	G4double fGuideSizeX = fGuideBorderX + fScintSizeX;
	G4double fGuideSizeY = fGuideBorderY + fScintSizeY;
	G4double fGuideSizeZ = fGuideBorderZ + fScintSizeZ;

	fSolidGuide_tmp	= new G4Box("Guide_tmp", 0.5*fGuideSizeX, 0.5*fGuideSizeY, 0.5*fGuideSizeZ);
	fSolidGuide_hole =  new G4Box("Guide_hole",0.5*fGuideHoleSizeXY, 0.5*fGuideHoleSizeXY, 0.6*fGuideSizeZ);
	fSolidGuide = new G4SubtractionSolid("Guide",fSolidGuide_tmp,fSolidScint_guide,0,G4ThreeVector(0,0,0.5*-(fGuideSizeZ-fScintSizeZ)));
	fSolidGuide = new G4SubtractionSolid("Guide",fSolidGuide,fSolidGuide_hole,0,G4ThreeVector(0,0,0));
    fLogicGuide = new G4LogicalVolume(fSolidGuide, fScintMaterial, "OpticalGrease");

	fLogicGuide->SetVisAttributes(G4Colour(0, 0, 1, 0.2));
	fPhysGuide 	= new G4PVPlacement(0, G4ThreeVector(0, 10*cm, 0), fLogicGuide, "Guide", fLogicWorld, true, 0, fCheckOverlaps);


	// G4double fGround;
	// fGround =  0.999;
	// if(fGround < 1){
		G4OpticalSurface* OpScintSurface = new G4OpticalSurface("OpScintSurface");
		OpScintSurface->SetModel(glisur);
		OpScintSurface->SetType(dielectric_dielectric);
		OpScintSurface->SetFinish(polished); //ground polished
	// 	OpScintSurface->SetPolish(fGround);
*/
	// 	new G4LogicalSkinSurface("ScintSurface", fLogicScint, OpScintSurface);	
	// }

	// G4OpticalSurface* OpScintSurface = new G4OpticalSurface("OpScintSurface");
	// OpScintSurface->SetModel(glisur);
	// OpScintSurface->SetType(dielectric_dielectric);
	// OpScintSurface->SetFinish(polished); //ground polished

	// G4OpticalSurface* OpScintSurface_Read = new G4OpticalSurface("OpScintSurface_Read");
	// OpScintSurface_Read->SetModel(glisur);
	// OpScintSurface_Read->SetType(dielectric_dielectric);
	// OpScintSurface_Read->SetFinish(ground);
	// OpScintSurface_Read->SetPolish(0.02);

	// new G4LogicalBorderSurface("OpScintSurface",  fPhysScint, fPhysWorld, OpScintSurface);
	// new G4LogicalBorderSurface("OpScintSurface_Read", fPhysScint, fPhysRead, OpScintSurface_Read);

	// G4double fGround;
	// fGround =  0.8;
	// if(fGround < 1){
	// 	G4OpticalSurface* OpScintSurface = new G4OpticalSurface("OpScintSurface");
	// 	OpScintSurface->SetModel(glisur);
	// 	OpScintSurface->SetType(dielectric_dielectric);
	// 	OpScintSurface->SetFinish(groundair); //ground polished
	// 	OpScintSurface->SetPolish(fGround);

	// 	new G4LogicalSkinSurface("ScintSurface", fLogicScint, OpScintSurface);	
	// }

	// Set how the volumes are visualized
    fLogicWorld	->SetVisAttributes(G4Colour(1, 1, 1, 0.1));
    fLogicScint	->SetVisAttributes(G4Colour(0.34, 0.57, 0.8, 0.5));
	fLogicScint2	->SetVisAttributes(G4Colour(0.8, 0.8, 0.34, 0.5));
	
    return fPhysWorld;
}

void DetectorConstruction::ConstructSDandField()
{
	auto sdManager = G4SDManager::GetSDMpointer();
   	G4String SDname;
	
	ScintSD* scint_SD = new ScintSD(SDname="Scint");
  	sdManager->AddNewDetector(scint_SD);
	fLogicScint->SetSensitiveDetector(scint_SD);


	ScintSD* scint_SD2 = new ScintSD(SDname="Scint2");
  	sdManager->AddNewDetector(scint_SD2);
	fLogicScint2->SetSensitiveDetector(scint_SD2);
		

	SiPMSD * SiPM_SD = new SiPMSD("SiPM");
  	sdManager->AddNewDetector(SiPM_SD);
	fLogicSiPM->SetSensitiveDetector(SiPM_SD);

	if(fVDOn){
		VirtualDetectorSD * VD_SD_2 = new VirtualDetectorSD("VirtualDetector2");
		fLogicVD_2->SetSensitiveDetector(VD_SD_2);
	}
}

/*
	From here functions used through the Messenger to modify the detector
*/

void DetectorConstruction :: SetScintSize(G4double size){
	fScintSizeX = size;
	fScintSizeY = size;
	fScintSizeZ = size;

	fWorldSizeX = std::max(fScintSizeX, fWorldSizeX);
	fWorldSizeY = std::max(fScintSizeY, fWorldSizeY);
	fWorldSizeZ = std::max(fScintSizeZ, fWorldSizeZ);

	fSolidScint->SetXHalfLength(size*0.5);
	fSolidScint->SetYHalfLength(size*0.5);
	fSolidScint->SetZHalfLength(size*0.5);

	fSolidWorld->SetXHalfLength(0.5*fWorldSizeX);
	fSolidWorld->SetYHalfLength(0.5*fWorldSizeY);
	fSolidWorld->SetZHalfLength(0.5*fWorldSizeZ);

	//fPhysScint->SetTranslation(G4ThreeVector(0, 0, 20*cm));
	G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

void DetectorConstruction :: SetScintSize(G4ThreeVector size){
	fScintSizeX = size.getX();
	fScintSizeY = size.getY();
	fScintSizeZ = size.getZ();

	fWorldSizeX = std::max(fScintSizeX, fWorldSizeX);
	fWorldSizeY = std::max(fScintSizeY, fWorldSizeY);
	fWorldSizeZ = std::max(fScintSizeZ, fWorldSizeZ);

	fSolidScint->SetXHalfLength(size.getX()*0.5);
	fSolidScint->SetYHalfLength(size.getY()*0.5);
	fSolidScint->SetZHalfLength(size.getZ()*0.5);

	fSolidWorld->SetXHalfLength(0.5*fWorldSizeX);
	fSolidWorld->SetYHalfLength(0.5*fWorldSizeY);
	fSolidWorld->SetZHalfLength(0.5*fWorldSizeZ);

	fSolidWorld	= new G4Box("World", 0.5*std::max(fScintSizeX, fWorldSizeX), 
		0.5*std::max(fScintSizeY, fWorldSizeY), 
		0.5*std::max(fScintSizeZ, fWorldSizeZ));

	//fPhysScint->SetTranslation(G4ThreeVector(0, 0, 20*cm));
	G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

void DetectorConstruction::SetScintMaterial(G4String name){
	if(name == "BC400") fScintMaterial = fBC400;
	else if(name == "LYSO") fScintMaterial = fLYSO;
	G4RunManager::GetRunManager()->ReinitializeGeometry();
}