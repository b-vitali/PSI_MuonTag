/// \file  DetectorConstruction.hh
/// \brief Implementation of the class to define the experimental setup

#include "DetectorConstruction.hh"
#include <math.h>  
//#include "G4MagneticField.hh"

#include "G4AutoDelete.hh"

#include "G4Trd.hh"

G4ThreadLocal 
G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = 0; 

DetectorConstruction::DetectorConstruction()
{
	//fVDOn = false;

	// Scintillator dimensions
    fScintSizeX = 50*mm;
    fScintSizeY = 50*mm;
    fScintSizeZ = 0.05*mm;

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
	G4int WHICH = 3;
	switch (WHICH)
	{
		case 0:
			return DefineVolumes();
			break;
		case 1:
			return DefineVolumes_guide();
			break;
		case 2:
			return DefineVolumes_doubleguide();
			break;
		case 3:
			return DefineVolumes_singlemuedm();
			break;
	}
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

	// Silicium
	G4NistManager* NISTman = G4NistManager::Instance();
	fSi = NISTman->FindOrBuildMaterial("G4_Si");

	// Silicon resin
	fSiResin = new G4Material("SiResin",z=1.,a=1.01*g/mole, 
		     		 density = universe_mean_density, kStateGas,
				 0.1 * kelvin, 1.e-19 * pascal);

	// Assign default materials
	fScintMaterial = fBC400;
	fSiPMMaterial  = fOG;
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

	// ----------------------------------------------------------	
	// Silicium
	// ----------------------------------------------------------	
	G4double Si_Energy[] = {.5*eV, 9.*eV};
	const G4int Sinum = sizeof(vacuum_Energy) / sizeof(G4double);

	G4double Si_RIND[] = {3.4, 3.4};
	assert(sizeof(Si_RIND) == sizeof(Si_Energy));

	G4MaterialPropertiesTable* Si_mt = new G4MaterialPropertiesTable();
	Si_mt->AddProperty("RINDEX", Si_Energy, Si_RIND, Sinum);
	fSi->SetMaterialPropertiesTable(Si_mt);

	// ----------------------------------------------------------	
	// Silicon resin
	// ----------------------------------------------------------	
	G4double SiRes_RIND[] = {1.41, 1.41};
	assert(sizeof(SiRes_RIND) == sizeof(Si_Energy));
	
	G4MaterialPropertiesTable* SiRes_mt = new G4MaterialPropertiesTable();
	SiRes_mt->AddProperty("RINDEX", Si_Energy, SiRes_RIND, Sinum);
	fSiResin->SetMaterialPropertiesTable(SiRes_mt);

	// ----------------------------------------------------------	
	// Optical grease
	// ----------------------------------------------------------	
	//? better if it was higher?
	G4double OG_RIND[] = {1.465, 1.465};
	assert(sizeof(OG_RIND) == sizeof(Si_Energy));
		
	G4MaterialPropertiesTable* OG_mt = new G4MaterialPropertiesTable();
	OG_mt->AddProperty("RINDEX", Si_Energy, OG_RIND, Sinum);
	fOG->SetMaterialPropertiesTable(OG_mt);
}

G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
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
	fElementSizeZ = 3*mm; //fScintSizeZ + 1*mm

	fSolidElement = new G4Box("Element", 0.5*(fElementSizeX),0.5*(fElementSizeY), 0.5*(fElementSizeZ));
    fLogicElement = new G4LogicalVolume(fSolidElement, fVacuum, "Element");
	fLogicElement->SetVisAttributes(G4Colour(0, 1, 0, 0.1));

	G4Rotate3D 	rotation =  G4Rotate3D(30*deg, G4ThreeVector(0, 1, 0)); //i*theta*deg std::cos(theta*i)
	G4Translate3D 	translate =  G4Translate3D(G4ThreeVector(0., 0., 3*cm));
	G4Transform3D transform = translate*rotation;

	// Read Solid and Phys. 
	fReadSizeX = 0.4*mm;
	fReadSizeY = 3*mm;
	fReadSizeZ = 3*mm;
	
	fSolidRead	= new G4Box("Read", 0.5*fReadSizeX,0.5*fReadSizeY, 0.5*fReadSizeZ);
    fLogicRead = new G4LogicalVolume(fSolidRead, fVacuum, "Read");
	fLogicRead->SetVisAttributes(G4Colour(0.0, 1.0, 1.0, 0.3));

	// grease Solid and Phys. 
	fSolidGrease = new G4Box("Read", 0.5*fReadSizeX*0.5,0.5*fReadSizeY, 0.5*fReadSizeZ);
    fLogicGrease = new G4LogicalVolume(fSolidGrease, fOG, "Grease");
	fLogicGrease->SetVisAttributes(G4Colour(1,0,0, 0.5));

	// VirtualDetector/SiPM Solid and Phys. 
	fSolidSiPM	= new G4Box("SiPM", 0.5*fReadSizeX*0.5, 0.5*fReadSizeY, 0.5*fReadSizeZ);
    fLogicSiPM 	= new G4LogicalVolume(fSolidSiPM, fSiPMMaterial, "SiPM");
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

	//fSolidScint2	= new G4Box("Scint2", 0.7*fScintSizeX, 0.7*fScintSizeY, 0.7*fScintSizeZ);
    //fLogicScint2 = new G4LogicalVolume(fSolidScint2, fScintMaterial, "Scint2");

	//? Put Scint and Read in element
	fPhysElement 	= new G4PVPlacement(0, G4ThreeVector(0, 0, 0), fLogicElement, "Element", fLogicWorld, true, 0, fCheckOverlaps);
	fPhysRead		= new G4PVPlacement(0, Read_pos1, fLogicRead, "Read0", fLogicElement, true, 0, fCheckOverlaps);
	fPhysRead		= new G4PVPlacement(0, Read_pos2, fLogicRead, "Read1", fLogicElement, true, 1, fCheckOverlaps);
	fPhysRead		= new G4PVPlacement(0, Read_pos3, fLogicRead, "Read2", fLogicElement, true, 2, fCheckOverlaps);
	fPhysScint		= new G4PVPlacement(0, Scint_pos, fLogicScint, "Scint", fLogicElement, false, fCheckOverlaps);

	fPhysElement 	= new G4PVPlacement(0, G4ThreeVector(0, 0, 2*cm), fLogicElement, "Element", fLogicWorld, true, 1, fCheckOverlaps);

	// transform = translate; //*rotation
	// fPhysElement 	= new G4PVPlacement(transform , fLogicElement, "Element", fLogicWorld, true, 1, fCheckOverlaps);
	// transform = translate*translate; //rotation*rotation
	// fPhysElement 	= new G4PVPlacement(transform , fLogicElement, "Element", fLogicWorld, true, 2, fCheckOverlaps);
	// transform = translate*translate*translate; //rotation*rotation
	//fPhysElement 	= new G4PVPlacement(transform , fLogicScint2, "Scint2", fLogicWorld, false, fCheckOverlaps);


	/*
	// G4double fGround;
	// fGround =  0.999;
	// if(fGround < 1){
		G4OpticalSurface* OpScintSurface = new G4OpticalSurface("OpScintSurface");
		OpScintSurface->SetModel(glisur);
		OpScintSurface->SetType(dielectric_dielectric);
		OpScintSurface->SetFinish(polished); //ground polished
	// 	OpScintSurface->SetPolish(fGround);

		new G4LogicalSkinSurface("ScintSurface", fLogicScint, OpScintSurface);	
	// }
	*/
	/*
	G4OpticalSurface* OpScintSurface = new G4OpticalSurface("OpScintSurface");
	OpScintSurface->SetModel(glisur);
	OpScintSurface->SetType(dielectric_dielectric);
	OpScintSurface->SetFinish(polished); //ground polished

	G4OpticalSurface* OpScintSurface_Read = new G4OpticalSurface("OpScintSurface_Read");
	OpScintSurface_Read->SetModel(glisur);
	OpScintSurface_Read->SetType(dielectric_dielectric);
	OpScintSurface_Read->SetFinish(ground);
	OpScintSurface_Read->SetPolish(0.02);

	new G4LogicalBorderSurface("OpScintSurface",  fPhysScint, fPhysWorld, OpScintSurface);
	new G4LogicalBorderSurface("OpScintSurface_Read", fPhysScint, fPhysRead, OpScintSurface_Read);

	*/
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

    fLogicWorld	->SetVisAttributes(G4Colour(1, 1, 1, 0.1));
    fLogicScint	->SetVisAttributes(G4Colour(0.34, 0.57, 0.8, 0.5));

    return fPhysWorld;
}


//----------------------------------------------------------------------
// Scintillator with light guide
//----------------------------------------------------------------------
G4VPhysicalVolume* DetectorConstruction::DefineVolumes_guide()
{
	fCheckOverlaps = true;
	G4int WHICH = 0; // various geometries

	// World Solid (size) -> Logical (material) -> PVPLacement (posiz, rotaz, to interact)
	fSolidWorld	= new G4Box("World", 0.5*fWorldSizeX, 0.5*fWorldSizeY, 0.5*fWorldSizeZ);
    fLogicWorld = new G4LogicalVolume(fSolidWorld, fVacuum_nogamma, "World");
    fPhysWorld	= new G4PVPlacement(0, G4ThreeVector(), fLogicWorld, "World", 0, false, 0, fCheckOverlaps);

	// Scintillator Solid (size) -> Logical (material) -> PVPLacement (posiz, rotaz, to interact)
	fSolidScint	= new G4Box("Scint", 0.5*fScintSizeX, 0.5*fScintSizeY, 0.5*fScintSizeZ);
    fLogicScint = new G4LogicalVolume(fSolidScint, fScintMaterial, "Scint");

	G4Rotate3D 	rotation =  G4Rotate3D(30*deg, G4ThreeVector(0, 1, 0)); //i*theta*deg std::cos(theta*i)
	G4Translate3D 	translate =  G4Translate3D(G4ThreeVector(0., 0., 3*cm));
	G4Transform3D transform = translate*rotation;

	// Read Solid and Phys. 
	fReadSizeX = 0.4*mm;
	fReadSizeY = 3*mm;
	fReadSizeZ = 3*mm;
	
	fSolidRead	= new G4Box("Read", 0.5*fReadSizeX,0.5*fReadSizeY, 0.5*fReadSizeZ);
    fLogicRead = new G4LogicalVolume(fSolidRead, fVacuum, "Read");
	fLogicRead->SetVisAttributes(G4Colour(0.0, 1.0, 1.0, 0.3));

	// grease Solid and Phys. 
	fSolidGrease = new G4Box("Read", 0.5*fReadSizeX*0.5,0.5*fReadSizeY, 0.5*fReadSizeZ);
    fLogicGrease = new G4LogicalVolume(fSolidGrease, fOG, "Grease");
	fLogicGrease->SetVisAttributes(G4Colour(1,0,0, 0.5));

	// VirtualDetector/SiPM Solid and Phys. 
	fSolidSiPM	= new G4Box("SiPM", 0.5*fReadSizeX*0.5, 0.5*fReadSizeY, 0.5*fReadSizeZ);
    fLogicSiPM 	= new G4LogicalVolume(fSolidSiPM, fSiPMMaterial, "SiPM");
    fLogicSiPM	->SetVisAttributes(G4Colour(0.8, 0.34, 0.68, 0.5));

	// Put Grease and SiPM in Read
    G4ThreeVector Grease_pos = G4ThreeVector(-(0.5*fReadSizeX*0.5), 0, 0); 
    G4ThreeVector SiPM_pos = G4ThreeVector(0.5*fReadSizeX*0.5, 0,0); //0.5*fReadSizeX
		
	fPhysRead 	= new G4PVPlacement(0, G4ThreeVector(0, 0, 0), fLogicRead, "Read", fLogicWorld, true, 0, fCheckOverlaps);
	fPhysGrease	= new G4PVPlacement(0, Grease_pos, fLogicGrease, "Grease", fLogicRead, false, fCheckOverlaps);
	fPhysSiPM		= new G4PVPlacement(0, SiPM_pos, fLogicSiPM, "SiPM", fLogicRead, false, fCheckOverlaps);

	// Scintillator
	fSolidScint_guide	= new G4Box("Scint_guide", 0.5*fScintSizeX, 0.5*fScintSizeY, 0.5*fScintSizeZ);
    fLogicScint_guide = new G4LogicalVolume(fSolidScint, fScintMaterial, "Scint");

	// Guide
	fGuideBorderX=1*mm;
    fGuideBorderY=1*mm;
    fGuideBorderZ=3*mm-fScintSizeZ;
	fGuideHoleSizeXY=10*mm;

	//? + 0.4*mm for grease between guide and scint?
	G4double fGuideSizeX = fGuideBorderX*2 + fScintSizeX;
	G4double fGuideSizeY = fGuideBorderY*2 + fScintSizeY;
	G4double fGuideSizeZ = fGuideBorderZ + fScintSizeZ;

	fSolidGuide_tmp	= new G4Box("Guide_tmp", 0.5*fGuideSizeX, 0.5*fGuideSizeY, 0.5*fGuideSizeZ);
	fSolidGuide_hole =  new G4Box("Guide_hole",0.5*fGuideHoleSizeXY, 0.5*fGuideHoleSizeXY, 0.6*fGuideSizeZ);
	G4Trd *SolidWedge = new G4Trd("wedge", fGuideHoleSizeXY*0.5, fGuideSizeX*0.5, fGuideHoleSizeXY*0.5, fGuideSizeY*0.5, 0.5*(fGuideSizeZ-fScintSizeZ));
	fSolidGuide = new G4SubtractionSolid("Guide",fSolidGuide_tmp,fSolidScint_guide,0,G4ThreeVector(0,0,-0.5*(fGuideSizeZ)));
	fSolidGuide = new G4SubtractionSolid("Guide",fSolidGuide,fSolidGuide_hole,0,G4ThreeVector(0,0,0));
	fSolidGuide = new G4SubtractionSolid("Guide",fSolidGuide,SolidWedge,0,G4ThreeVector(0,0,fScintSizeZ*0.5));

	G4LogicalVolume * LogicWedge = new G4LogicalVolume(SolidWedge, fOG, "LogicWedge");
	//G4PVPlacement* PhysWedge	= new G4PVPlacement(0, G4ThreeVector(0, 0, 1*cm), LogicWedge, "PhysWedge", fLogicWorld, false, fCheckOverlaps);

	//? material??
    fLogicGuide = new G4LogicalVolume(fSolidGuide, fOG, "fLogicGuide"); 
	fLogicGuide->SetVisAttributes(G4Colour(0, 0, 1, 0.2));

	//Element for scintillator with light guide
	fElementSizeX_guide = fGuideSizeX + 1*mm;
	fElementSizeY_guide = fGuideSizeY + 1*mm;
	fElementSizeZ_guide = fGuideSizeZ + 1*mm;

	fSolidElement_guide = new G4Box("Element_guide", 0.5*(fElementSizeX_guide),0.5*(fElementSizeY_guide), 0.5*(fElementSizeZ_guide));
    fLogicElement_guide = new G4LogicalVolume(fSolidElement_guide, fVacuum, "Element_guide");
	fLogicElement_guide->SetVisAttributes(G4Colour(0, 1, 0, 0.1));


   	G4ThreeVector Read_guide_pos1 = G4ThreeVector(0.5*fReadSizeX+0.5*fGuideSizeX, 0.2*fGuideSizeY,0); //0.5*fReadSizeX
   	G4ThreeVector Read_guide_pos2 = G4ThreeVector(0.5*fReadSizeX+0.5*fGuideSizeX, -0*fGuideSizeY,0); //0.5*fReadSizeX
   	G4ThreeVector Read_guide_pos3 = G4ThreeVector(0.5*fReadSizeX+0.5*fGuideSizeX, -0.2*fGuideSizeY,0); //0.5*fReadSizeX


	// put guide and scint in element
	fPhysElement_guide 	= new G4PVPlacement(0, G4ThreeVector(0, 0, 0), fLogicElement_guide, "Element_guide", fLogicWorld, false, fCheckOverlaps);
	fPhysGuide 	= new G4PVPlacement(0, G4ThreeVector(0, 0, 0), fLogicGuide, "Guide", fLogicElement_guide, false, fCheckOverlaps);
	fPhysScint_guide 	= new G4PVPlacement(0, G4ThreeVector(0, 0, (fScintSizeZ-fGuideSizeZ)*0.5), fLogicScint_guide, "LogicScint_guide", fLogicElement_guide, false, fCheckOverlaps);
	fPhysRead		= new G4PVPlacement(0, Read_guide_pos1, fLogicRead, "Read0", fLogicElement_guide, true, 0, fCheckOverlaps);
	fPhysRead		= new G4PVPlacement(0, Read_guide_pos2, fLogicRead, "Read1", fLogicElement_guide, true, 1, fCheckOverlaps);
	fPhysRead		= new G4PVPlacement(0, Read_guide_pos3, fLogicRead, "Read2", fLogicElement_guide, true, 2, fCheckOverlaps);

	/*
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
    fLogicScint_guide	->SetVisAttributes(G4Colour(0.34, 0.57, 0.8, 0.5));
	
    return fPhysWorld;
}

//----------------------------------------------------------------------
// Scintillator with double-light guide
//----------------------------------------------------------------------
G4VPhysicalVolume* DetectorConstruction::DefineVolumes_doubleguide()
{
	fCheckOverlaps = true;
	G4int WHICH = 0; // various geometries

	// World Solid (size) -> Logical (material) -> PVPLacement (posiz, rotaz, to interact)
	fSolidWorld	= new G4Box("World", 0.5*fWorldSizeX, 0.5*fWorldSizeY, 0.5*fWorldSizeZ);
    fLogicWorld = new G4LogicalVolume(fSolidWorld, fVacuum_nogamma, "World");
    fPhysWorld	= new G4PVPlacement(0, G4ThreeVector(), fLogicWorld, "World", 0, false, 0, fCheckOverlaps);

	// Scintillator Solid (size) -> Logical (material) -> PVPLacement (posiz, rotaz, to interact)
	fSolidScint	= new G4Box("Scint", 0.5*fScintSizeX, 0.5*fScintSizeY, 0.5*fScintSizeZ);
    fLogicScint = new G4LogicalVolume(fSolidScint, fScintMaterial, "Scint");

	G4Rotate3D 	rotation =  G4Rotate3D(30*deg, G4ThreeVector(0, 1, 0)); //i*theta*deg std::cos(theta*i)
	G4Translate3D 	translate =  G4Translate3D(G4ThreeVector(0., 0., 3*cm));
	G4Transform3D transform = translate*rotation;

	// Read Solid and Phys. 
	fReadSizeX = 0.4*mm;
	fReadSizeY = 3*mm;
	fReadSizeZ = 3*mm;
	
	fSolidRead	= new G4Box("Read", 0.5*fReadSizeX,0.5*fReadSizeY, 0.5*fReadSizeZ);
    fLogicRead = new G4LogicalVolume(fSolidRead, fVacuum, "Read");
	fLogicRead->SetVisAttributes(G4Colour(0.0, 1.0, 1.0, 0.3));

	// grease Solid and Phys. 
	fSolidGrease = new G4Box("Read", 0.5*fReadSizeX*0.5,0.5*fReadSizeY, 0.5*fReadSizeZ);
    fLogicGrease = new G4LogicalVolume(fSolidGrease, fOG, "Grease");
	fLogicGrease->SetVisAttributes(G4Colour(1,0,0, 0.5));

	// VirtualDetector/SiPM Solid and Phys. 
	fSolidSiPM	= new G4Box("SiPM", 0.5*fReadSizeX*0.5, 0.5*fReadSizeY, 0.5*fReadSizeZ);
    fLogicSiPM 	= new G4LogicalVolume(fSolidSiPM, fSiPMMaterial, "SiPM");
    fLogicSiPM	->SetVisAttributes(G4Colour(0.8, 0.34, 0.68, 0.5));

	// Put Grease and SiPM in Read
    G4ThreeVector Grease_pos = G4ThreeVector(-(0.5*fReadSizeX*0.5), 0, 0); 
    G4ThreeVector SiPM_pos = G4ThreeVector(0.5*fReadSizeX*0.5, 0,0); //0.5*fReadSizeX
		
	//fPhysRead 	= new G4PVPlacement(0, G4ThreeVector(0, 0, 0), fLogicRead, "Read", fLogicWorld, true, 0, fCheckOverlaps);
	fPhysGrease	= new G4PVPlacement(0, Grease_pos, fLogicGrease, "Grease", fLogicRead, false, fCheckOverlaps);
	fPhysSiPM		= new G4PVPlacement(0, SiPM_pos, fLogicSiPM, "SiPM", fLogicRead, false, fCheckOverlaps);

	// Scintillator
	fSolidScint_guide	= new G4Box("Scint_guide", 0.5*fScintSizeX, 0.5*fScintSizeY, 0.5*fScintSizeZ);
    fLogicScint_guide = new G4LogicalVolume(fSolidScint, fScintMaterial, "Scint");

	// Guide
	fGuideBorderX=1*mm;
    fGuideBorderY=1*mm;
    fGuideBorderZ=(3*mm)*0.5;
	fGuideHoleSizeXY=10*mm;

	//? + 0.4*mm for grease between guide and scint?
	G4double fGuideSizeX = fGuideBorderX*2 + fScintSizeX;
	G4double fGuideSizeY = fGuideBorderY*2 + fScintSizeY;
	G4double fGuideSizeZ = fGuideBorderZ;

	fSolidGuide_tmp	= new G4Box("Guide_tmp", 0.5*fGuideSizeX, 0.5*fGuideSizeY, 0.5*fGuideSizeZ);
	fSolidGuide_hole =  new G4Box("Guide_hole",0.5*fGuideHoleSizeXY, 0.5*fGuideHoleSizeXY, 0.6*fGuideSizeZ);
	G4Trd *SolidWedge = new G4Trd("wedge", fGuideHoleSizeXY*0.5, fGuideSizeX*0.5, fGuideHoleSizeXY*0.5, fGuideSizeY*0.5, 0.5*(fGuideSizeZ-fScintSizeZ));
	fSolidGuide = new G4SubtractionSolid("Guide",fSolidGuide_tmp,fSolidScint_guide,0,G4ThreeVector(0,0,-0.5*(fGuideSizeZ)));
	fSolidGuide = new G4SubtractionSolid("Guide",fSolidGuide,fSolidGuide_hole,0,G4ThreeVector(0,0,0));
	// fSolidGuide = new G4SubtractionSolid("Guide",fSolidGuide,SolidWedge,0,G4ThreeVector(0,0,fScintSizeZ*0.5));

	G4LogicalVolume * LogicWedge = new G4LogicalVolume(SolidWedge, fOG, "LogicWedge");
	//G4PVPlacement* PhysWedge	= new G4PVPlacement(0, G4ThreeVector(0, 0, 1*cm), LogicWedge, "PhysWedge", fLogicWorld, false, fCheckOverlaps);

	//? material??
    fLogicGuide = new G4LogicalVolume(fSolidGuide, fOG, "fLogicGuide"); 
	fLogicGuide->SetVisAttributes(G4Colour(0, 0, 1, 0.2));

	//Element for scintillator with light guide
	fElementSizeX_guide = fGuideSizeX + 1*mm;
	fElementSizeY_guide = fGuideSizeY + 1*mm;
	fElementSizeZ_guide = fGuideSizeZ*2 + 1*mm;

	fSolidElement_guide = new G4Box("Element_guide", 0.5*(fElementSizeX_guide),0.5*(fElementSizeY_guide), 0.5*(fElementSizeZ_guide));
    fLogicElement_guide = new G4LogicalVolume(fSolidElement_guide, fVacuum, "Element_guide");
	fLogicElement_guide->SetVisAttributes(G4Colour(0, 1, 0, 0.1));


   	G4ThreeVector Read_guide_pos1 = G4ThreeVector(0.5*fReadSizeX+0.5*fGuideSizeX, 0.2*fGuideSizeY,0); //0.5*fReadSizeX
   	G4ThreeVector Read_guide_pos2 = G4ThreeVector(0.5*fReadSizeX+0.5*fGuideSizeX, -0*fGuideSizeY,0); //0.5*fReadSizeX
   	G4ThreeVector Read_guide_pos3 = G4ThreeVector(0.5*fReadSizeX+0.5*fGuideSizeX, -0.2*fGuideSizeY,0); //0.5*fReadSizeX

	rotation =  G4Rotate3D(180*deg, G4ThreeVector(0, 1, 0)); //i*theta*deg std::cos(theta*i)
	translate =  G4Translate3D(G4ThreeVector(0., 0., -(fGuideSizeZ)*0.5));
	transform = translate*rotation;

	// put guide and scint in element
	fPhysElement_guide 	= new G4PVPlacement(0, G4ThreeVector(0, 0, 0), fLogicElement_guide, "Element_guide", fLogicWorld, false, fCheckOverlaps);
	fPhysGuide 	= new G4PVPlacement(0, G4ThreeVector(0, 0, (fGuideSizeZ)*0.5), fLogicGuide, "Guide", fLogicElement_guide, false, fCheckOverlaps);
	fPhysGuide 	= new G4PVPlacement(transform, fLogicGuide, "Guide", fLogicElement_guide, false, fCheckOverlaps);
	fPhysScint_guide 	= new G4PVPlacement(0, G4ThreeVector(0, 0, 0), fLogicScint_guide, "LogicScint_guide", fLogicElement_guide, false, fCheckOverlaps);
	fPhysRead		= new G4PVPlacement(0, Read_guide_pos1, fLogicRead, "Read0", fLogicElement_guide, true, 0, fCheckOverlaps);
	fPhysRead		= new G4PVPlacement(0, Read_guide_pos2, fLogicRead, "Read1", fLogicElement_guide, true, 1, fCheckOverlaps);
	fPhysRead		= new G4PVPlacement(0, Read_guide_pos3, fLogicRead, "Read2", fLogicElement_guide, true, 2, fCheckOverlaps);

	/*
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
    fLogicScint_guide	->SetVisAttributes(G4Colour(0.34, 0.57, 0.8, 0.5));
	
    return fPhysWorld;
}


//----------------------------------------------------------------------
// Single MuEDM cylinder's scint
//----------------------------------------------------------------------
G4VPhysicalVolume* DetectorConstruction::DefineVolumes_singlemuedm()
{
	// Scintillator dimensions
    fScintSizeX = 200*mm;
    fScintSizeY = 10*mm;
    fScintSizeZ = 3*mm;

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
	fElementSizeZ = 3*mm; //fScintSizeZ + 1*mm

	fSolidElement = new G4Box("Element", 0.5*(fElementSizeX),0.5*(fElementSizeY), 0.5*(fElementSizeZ));
    fLogicElement = new G4LogicalVolume(fSolidElement, fVacuum, "Element");
	fLogicElement->SetVisAttributes(G4Colour(0, 1, 0, 0.1));

	G4Rotate3D 	rotation =  G4Rotate3D(30*deg, G4ThreeVector(0, 1, 0)); //i*theta*deg std::cos(theta*i)
	G4Translate3D 	translate =  G4Translate3D(G4ThreeVector(0., 0., 3*cm));
	G4Transform3D transform = translate*rotation;

	// Read Solid and Phys. 
	fReadSizeX = 0.4*mm;
	fReadSizeY = 3*mm;
	fReadSizeZ = 3*mm;
	
	fSolidRead	= new G4Box("Read", 0.5*fReadSizeX,0.5*fReadSizeY, 0.5*fReadSizeZ);
    fLogicRead = new G4LogicalVolume(fSolidRead, fVacuum, "Read");
	fLogicRead->SetVisAttributes(G4Colour(0.0, 1.0, 1.0, 0.3));

	// grease Solid and Phys. 
	fSolidGrease = new G4Box("Read", 0.5*fReadSizeX*0.5,0.5*fReadSizeY, 0.5*fReadSizeZ);
    fLogicGrease = new G4LogicalVolume(fSolidGrease, fOG, "Grease");
	fLogicGrease->SetVisAttributes(G4Colour(1,0,0, 0.5));

	// VirtualDetector/SiPM Solid and Phys. 
	fSolidSiPM	= new G4Box("SiPM", 0.5*fReadSizeX*0.5, 0.5*fReadSizeY, 0.5*fReadSizeZ);
    fLogicSiPM 	= new G4LogicalVolume(fSolidSiPM, fSiPMMaterial, "SiPM");
    fLogicSiPM	->SetVisAttributes(G4Colour(0.8, 0.34, 0.68, 0.5));

	// Put Grease and SiPM in Read
    G4ThreeVector Grease_pos = G4ThreeVector(-(0.5*fReadSizeX*0.5), 0, 0); 
    G4ThreeVector SiPM_pos = G4ThreeVector(0.5*fReadSizeX*0.5, 0,0); //0.5*fReadSizeX
		
	//fPhysRead 	= new G4PVPlacement(0, G4ThreeVector(0, 0, 0), fLogicRead, "Read", fLogicWorld, true, 0, fCheckOverlaps);
	fPhysGrease	= new G4PVPlacement(0, Grease_pos, fLogicGrease, "Grease", fLogicRead, false, fCheckOverlaps);
	fPhysSiPM		= new G4PVPlacement(0, SiPM_pos, fLogicSiPM, "SiPM", fLogicRead, false, fCheckOverlaps);

	// Put Scint and Read in element
	G4ThreeVector Scint_pos = G4ThreeVector(0, 0, 0); 
   	G4ThreeVector Read_pos1 = G4ThreeVector(0.5*fReadSizeX+0.5*fScintSizeX, (3+0.25)*mm,0); //0.5*fReadSizeX
   	G4ThreeVector Read_pos2 = G4ThreeVector(0.5*fReadSizeX+0.5*fScintSizeX, 0, 0); //0.5*fReadSizeX
   	G4ThreeVector Read_pos3 = G4ThreeVector(0.5*fReadSizeX+0.5*fScintSizeX, -(3+0.25)*mm,0); //0.5*fReadSizeX

	//fSolidScint2	= new G4Box("Scint2", 0.7*fScintSizeX, 0.7*fScintSizeY, 0.7*fScintSizeZ);
    //fLogicScint2 = new G4LogicalVolume(fSolidScint2, fScintMaterial, "Scint2");

	//? Put Scint and Read in element
	fPhysElement 	= new G4PVPlacement(0, G4ThreeVector(0, 0, 0), fLogicElement, "Element", fLogicWorld, false, fCheckOverlaps);
	fPhysRead		= new G4PVPlacement(0, Read_pos1, fLogicRead, "Read0", fLogicElement, true, 0, fCheckOverlaps);
	fPhysRead		= new G4PVPlacement(0, Read_pos2, fLogicRead, "Read1", fLogicElement, true, 1, fCheckOverlaps);
	fPhysRead		= new G4PVPlacement(0, Read_pos3, fLogicRead, "Read2", fLogicElement, true, 2, fCheckOverlaps);
	
	rotation =  G4Rotate3D(180*deg, G4ThreeVector(0, 1, 0)); //i*theta*deg std::cos(theta*i)
	translate =  G4Translate3D(Read_pos1);
	transform = rotation*translate;
	fPhysRead		= new G4PVPlacement(transform, fLogicRead, "Read3", fLogicElement, true, 3, fCheckOverlaps);
	translate =  G4Translate3D(Read_pos2);
	transform = rotation*translate;
	fPhysRead		= new G4PVPlacement(transform, fLogicRead, "Read4", fLogicElement, true, 4, fCheckOverlaps);
	translate =  G4Translate3D(Read_pos3);
	transform = rotation*translate;
	fPhysRead		= new G4PVPlacement(transform, fLogicRead, "Read5", fLogicElement, true, 5, fCheckOverlaps);

	fPhysScint		= new G4PVPlacement(0, Scint_pos, fLogicScint, "Scint", fLogicElement, false, fCheckOverlaps);

	// transform = translate; //*rotation
	// fPhysElement 	= new G4PVPlacement(transform , fLogicElement, "Element", fLogicWorld, true, 1, fCheckOverlaps);
	// transform = translate*translate; //rotation*rotation
	// fPhysElement 	= new G4PVPlacement(transform , fLogicElement, "Element", fLogicWorld, true, 2, fCheckOverlaps);
	// transform = translate*translate*translate; //rotation*rotation
	//fPhysElement 	= new G4PVPlacement(transform , fLogicScint2, "Scint2", fLogicWorld, false, fCheckOverlaps);


	/*
	// G4double fGround;
	// fGround =  0.999;
	// if(fGround < 1){
		G4OpticalSurface* OpScintSurface = new G4OpticalSurface("OpScintSurface");
		OpScintSurface->SetModel(glisur);
		OpScintSurface->SetType(dielectric_dielectric);
		OpScintSurface->SetFinish(polished); //ground polished
	// 	OpScintSurface->SetPolish(fGround);

		new G4LogicalSkinSurface("ScintSurface", fLogicScint, OpScintSurface);	
	// }
	*/
	/*
	G4OpticalSurface* OpScintSurface = new G4OpticalSurface("OpScintSurface");
	OpScintSurface->SetModel(glisur);
	OpScintSurface->SetType(dielectric_dielectric);
	OpScintSurface->SetFinish(polished); //ground polished

	G4OpticalSurface* OpScintSurface_Read = new G4OpticalSurface("OpScintSurface_Read");
	OpScintSurface_Read->SetModel(glisur);
	OpScintSurface_Read->SetType(dielectric_dielectric);
	OpScintSurface_Read->SetFinish(ground);
	OpScintSurface_Read->SetPolish(0.02);

	new G4LogicalBorderSurface("OpScintSurface",  fPhysScint, fPhysWorld, OpScintSurface);
	new G4LogicalBorderSurface("OpScintSurface_Read", fPhysScint, fPhysRead, OpScintSurface_Read);

	*/
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

    fLogicWorld	->SetVisAttributes(G4Colour(1, 1, 1, 0.1));
    fLogicScint	->SetVisAttributes(G4Colour(0.34, 0.57, 0.8, 0.5));

    return fPhysWorld;
}

void DetectorConstruction::ConstructSDandField()
{
	auto sdManager = G4SDManager::GetSDMpointer();
   	G4String SDname;
	
	if(fLogicScint){
		ScintSD* scint_SD = new ScintSD(SDname="Scint");
  		sdManager->AddNewDetector(scint_SD);
		fLogicScint->SetSensitiveDetector(scint_SD);
	}

	if(fLogicScint2){
		ScintSD* scint_SD2 = new ScintSD(SDname="Scint2");
  		sdManager->AddNewDetector(scint_SD2);
		fLogicScint2->SetSensitiveDetector(scint_SD2);
	}

	if(fLogicSiPM){
		SiPMSD * SiPM_SD = new SiPMSD("SiPM");
  		sdManager->AddNewDetector(SiPM_SD);
		fLogicSiPM->SetSensitiveDetector(SiPM_SD);
	}
	if(fLogicVD_2){
		VirtualDetectorSD * VD_SD_2 = new VirtualDetectorSD("VirtualDetector2");
		fLogicVD_2->SetSensitiveDetector(VD_SD_2);
	}

	G4ThreeVector fieldValue(0.,-3*tesla,0.);
	fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
	//fMagFieldMessenger->SetVerboseLevel(1);
	
	// Register the field messenger for deleting
	G4AutoDelete::Register(fMagFieldMessenger);

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