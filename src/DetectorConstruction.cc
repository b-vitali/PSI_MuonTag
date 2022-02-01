/// \file  DetectorConstruction.cc
/// \brief Class to define the experimental setup

#include "DetectorConstruction.hh"


DetectorConstruction::DetectorConstruction()
{
	// At creation it calls for the function creating the materials 
	DefineMaterials();
	DefineOpticalProperties();
}

DetectorConstruction :: ~DetectorConstruction()
{}

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
	// World dimentions
    G4double world_sizeX = 30*cm;
    G4double world_sizeY = 30*cm;
    G4double world_sizeZ = 50*cm;
	
	// Scintillator dimensions
    G4double scint_sizeX = 20*cm;
    G4double scint_sizeY = 20*cm;
    G4double scint_sizeZ = 5*mm;

	// VirtualDetector dimensions
    G4double vd_sizeX = 20*cm;
    G4double vd_sizeY = 20*cm;
    G4double vd_sizeZ = 10*mm;

	// World Solid (size) -> Logical (material) -> PVPLacement (posiz, rotaz, to interact)
	fSolidWorld	= new G4Box("World", 0.5*world_sizeX, 0.5*world_sizeY, 0.5*world_sizeZ);
    fLogicWorld = new G4LogicalVolume(fSolidWorld, fVacuum, "World");
    fPhysWorld	= new G4PVPlacement(0, G4ThreeVector(), fLogicWorld, "World", 0, false, 0, fCheckOverlaps);

	// Scintillator Solid (size) -> Logical (material) -> PVPLacement (posiz, rotaz, to interact)
	fSolidScint	= new G4Box("Scint", 0.5*scint_sizeX, 0.5*scint_sizeY, 0.5*scint_sizeZ);
    fLogicScint = new G4LogicalVolume(fSolidScint, fBC400, "Scint");
    fPhysScint	= new G4PVPlacement(0, G4ThreeVector(0., 0., 20*cm), fLogicScint, "Scint", fLogicWorld, false, 0, fCheckOverlaps);

	// VirtualDetector Solid (size) -> Logical (material) -> PVPLacement (posiz, rotaz, to interact)
	fSolidVD	= new G4Box("VD", 0.5*vd_sizeX, 0.5*vd_sizeY, 0.5*vd_sizeZ);
    fLogicVD 	= new G4LogicalVolume(fSolidVD, fVacuum, "VD");
    fPhysVD		= new G4PVPlacement(0, G4ThreeVector(0., 0., 5*cm), fLogicVD, "VD", fLogicWorld, false, 0, fCheckOverlaps);

    fLogicWorld	->SetVisAttributes(G4Colour(1, 1, 1, 0.1));
    fLogicScint	->SetVisAttributes(G4Colour(0.34, 0.57, 0.8, 0.3));
    fLogicVD	->SetVisAttributes(G4Colour(0.8, 0.34, 0.68, 0.2));
	
    return fPhysWorld;
}

void DetectorConstruction::ConstructSDandField()
{
	VirtualDetectorSD * VD_SD = new VirtualDetectorSD("VirtualDetector");
	fLogicVD->SetSensitiveDetector(VD_SD);
}