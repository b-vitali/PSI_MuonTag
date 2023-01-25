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
G4double	theta_out; 	
G4double	theta_scint_out; 	
G4int 		N_out;		

G4double 	r_in;
G4double	theta_in; 	
G4double	theta_scint_in; 	
G4int 		N_in;		

DetectorConstruction::DetectorConstruction()
{
	fVDOn = false;

	// Scintillator out dimensions
	fScintSizeX_out = 0.2*mm;
    fScintSizeY_out = 5*cm;
    fScintSizeZ_out = 2*cm;

	// Scintillator in dimensions
	fScintSizeX_in = 0.2*mm;
    fScintSizeY_in = 5*cm;
    fScintSizeZ_in = 1*cm;

	// VirtualDetector dimensions
    fVDSizeX = 20*cm;
    fVDSizeY = 20*cm;
    fVDSizeZ = 10*mm;

	// World dimentions
    fWorldSizeX = 50*cm;
    fWorldSizeY = 50*cm;
    fWorldSizeZ = 50*cm;
	
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
	// return DefineVolumes();
	return DefineVolumes_SciFi();
	// return DefineVolumes_MuEDM();

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

	fBC400_noscint = new G4Material("BC400", density = 1.023*g/cm3, 2);
	fBC400_noscint->AddElement(fC, 1000);
	fBC400_noscint->AddElement(fH, 1103);

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
	fSiPMMaterial  = fSiResin;


	//! SciFi
	/// Materials
	// BCF10
	fBCF10 = new G4Material("BCF10", density = 1.05*g/cm3, 2);
	fBCF10->AddElement(fC, 485);
	fBCF10->AddElement(fH, 482);

	// BCF12
	fBCF12 = new G4Material("BCF12", density = 1.05*g/cm3, 2);
	fBCF12->AddElement(fC, 485);
	fBCF12->AddElement(fH, 482);

	// BCF20
	fBCF20 = new G4Material("BCF20", density = 1.05*g/cm3, 2);
	fBCF20->AddElement(fC, 485);
	fBCF20->AddElement(fH, 482);

	// First Cladding: PMMA
	fFClad = new G4Material("FClad", density = 1.2*g/cm3, 3);
	fFClad->AddElement(fC, 5);
	fFClad->AddElement(fH, 8);
	fFClad->AddElement(fO, 2);

	// Second Cladding: PMMA EMA
	fSClad = new G4Material("SClad", density = 1.2*g/cm3, 3);
	fSClad->AddElement(fC, 5);
	fSClad->AddElement(fH, 8);
	fSClad->AddElement(fO, 2);
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

	//! SciFi

	//  BCF10 optics
	myfile.open("../tables/BCF10_light_out.txt");
	while(true){
		myfile >> tempe >> tempscint;
		energy.push_back(1239.84197/tempe);
		scint.push_back(tempscint);
		
		if(myfile.eof()) break;
	}
	myfile.close();

	assert(energy.size() == scint.size());
	const G4int bcf10 = int(energy.size());

	G4double* BCF10_Energy = new G4double[bcf10];
	G4double* BCF10_SCINT = new G4double[bcf10];

	G4double* BCF10_RIND = new G4double[bcf10];
	G4double* BCF10_ABSL = new G4double[bcf10];
	
	for(int i = 0; i < bcf10; i++){
		BCF10_Energy[i] = energy.at(i)*eV;
		BCF10_SCINT[i] = scint.at(i);
		BCF10_RIND[i] = 1.6;
		BCF10_ABSL[i] = 220*cm;
	}
	
	energy.clear();
	scint.clear();
	
	
	assert(sizeof(BCF10_SCINT) == sizeof(BCF10_Energy));
	
	assert(sizeof(BCF10_RIND) == sizeof(BCF10_Energy));
	
	assert(sizeof(BCF10_ABSL) == sizeof(BCF10_Energy));

	fBCF10_mt = new G4MaterialPropertiesTable();
	fBCF10_mt->AddProperty(       "RINDEX", BCF10_Energy,  BCF10_RIND, bcf10);
	fBCF10_mt->AddProperty(    "ABSLENGTH", BCF10_Energy,  BCF10_ABSL, bcf10);
	fBCF10_mt->AddProperty("FASTCOMPONENT", BCF10_Energy, BCF10_SCINT, bcf10);
	
	fBCF10_mt->AddConstProperty("SCINTILLATIONYIELD",        8000./MeV);
	fBCF10_mt->AddConstProperty(   "RESOLUTIONSCALE",                 1);
	fBCF10_mt->AddConstProperty(  "SLOWTIMECONSTANT",            2.7*ns);
	fBCF10_mt->AddConstProperty(  "FASTTIMECONSTANT",            2.7*ns);
	fBCF10_mt->AddConstProperty(        "YIELDRATIO",                 0.);
	
	fBCF10->SetMaterialPropertiesTable(fBCF10_mt);

	//  Set Birks Constant
	fBCF10->GetIonisation()->SetBirksConstant(0.117645*mm/MeV);

	//  BCF12 optics
	myfile.open("../tables/BCF12_light_out.txt");
	while(true){
		myfile >> tempe >> tempscint;
		energy.push_back(1239.84197/tempe);
		scint.push_back(tempscint);
		
		if(myfile.eof()) break;
	}
	myfile.close();

	assert(energy.size() == scint.size());
	const G4int bcf12 = int(energy.size());

	G4double* BCF12_Energy = new G4double[bcf12];
	G4double* BCF12_SCINT = new G4double[bcf12];

	G4double* BCF12_RIND = new G4double[bcf12];
	G4double* BCF12_ABSL = new G4double[bcf12];
	
	for(int i = 0; i < bcf12; i++){
		BCF12_Energy[i] = energy.at(i)*eV;
		BCF12_SCINT[i] = scint.at(i);
		BCF12_RIND[i] = 1.6;
		BCF12_ABSL[i] = 270*cm;
	}
	
	energy.clear();
	scint.clear();
	
	
	assert(sizeof(BCF12_SCINT) == sizeof(BCF12_Energy));
	
	assert(sizeof(BCF12_RIND) == sizeof(BCF12_Energy));
	
	assert(sizeof(BCF12_ABSL) == sizeof(BCF12_Energy));

	fBCF12_mt = new G4MaterialPropertiesTable();
	fBCF12_mt->AddProperty(       "RINDEX", BCF12_Energy,  BCF12_RIND, bcf12);
	fBCF12_mt->AddProperty(    "ABSLENGTH", BCF12_Energy,  BCF12_ABSL, bcf12);
	fBCF12_mt->AddProperty("FASTCOMPONENT", BCF12_Energy, BCF12_SCINT, bcf12);
	
	fBCF12_mt->AddConstProperty("SCINTILLATIONYIELD",        8000./MeV);
	fBCF12_mt->AddConstProperty(   "RESOLUTIONSCALE",                 1);
	fBCF12_mt->AddConstProperty(  "SLOWTIMECONSTANT",            3.2*ns);
	fBCF12_mt->AddConstProperty(  "FASTTIMECONSTANT",            3.2*ns);
	fBCF12_mt->AddConstProperty(        "YIELDRATIO",                 0.);
	
	fBCF12->SetMaterialPropertiesTable(fBCF12_mt);

	//  Set Birks Constant
	fBCF12->GetIonisation()->SetBirksConstant(0.117645*mm/MeV);

	//  BCF20 optics
	myfile.open("../tables/BCF20_light_out.txt");
	while(true){
		myfile >> tempe >> tempscint;
		energy.push_back(1239.84197/tempe);
		scint.push_back(tempscint);
		
		if(myfile.eof()) break;
	}
	myfile.close();

	assert(energy.size() == scint.size());
	const G4int bcf20 = int(energy.size());

	G4double* BCF20_Energy = new G4double[bcf20];
	G4double* BCF20_SCINT = new G4double[bcf20];

	G4double* BCF20_RIND = new G4double[bcf20];
	G4double* BCF20_ABSL = new G4double[bcf20];
	
	for(int i = 0; i < bcf20; i++){
		BCF20_Energy[i] = energy.at(i)*eV;
		BCF20_SCINT[i] = scint.at(i);
		BCF20_RIND[i] = 1.6;
		BCF20_ABSL[i] = 350*cm;
	}
	
	energy.clear();
	scint.clear();
	
	
	assert(sizeof(BCF20_SCINT) == sizeof(BCF20_Energy));
	
	assert(sizeof(BCF20_RIND) == sizeof(BCF20_Energy));

	assert(sizeof(BCF20_ABSL) == sizeof(BCF20_Energy));

	fBCF20_mt = new G4MaterialPropertiesTable();
	fBCF20_mt->AddProperty(       "RINDEX", BCF20_Energy,  BCF20_RIND, bcf20);
	fBCF20_mt->AddProperty(    "ABSLENGTH", BCF20_Energy,  BCF20_ABSL, bcf20);
	fBCF20_mt->AddProperty("FASTCOMPONENT", BCF20_Energy, BCF20_SCINT, bcf20);
	
	fBCF20_mt->AddConstProperty("SCINTILLATIONYIELD",        8000./MeV);
	fBCF20_mt->AddConstProperty(   "RESOLUTIONSCALE",                 1);
	fBCF20_mt->AddConstProperty(  "SLOWTIMECONSTANT",            2.7*ns);
	fBCF20_mt->AddConstProperty(  "FASTTIMECONSTANT",            2.7*ns);
	fBCF20_mt->AddConstProperty(        "YIELDRATIO",                 0.);
	
	fBCF20->SetMaterialPropertiesTable(fBCF20_mt);

	//  Set Birks Constant
	fBCF20->GetIonisation()->SetBirksConstant(0.117645*mm/MeV);

	// Cladding
	G4double fCladRIND = 1.49;
	G4double fClad_RIND[] = {fCladRIND, fCladRIND};
	assert(sizeof(fClad_RIND) == sizeof(vacuum_Energy));

	G4MaterialPropertiesTable* fClad_mt = new G4MaterialPropertiesTable();
	fClad_mt->AddProperty("RINDEX", vacuum_Energy, fClad_RIND, vacnum);
	fFClad->SetMaterialPropertiesTable(fClad_mt);


	G4double sCladRIND = 1.42;
	G4double sClad_RIND[] = {sCladRIND, sCladRIND};
	assert(sizeof(sClad_RIND) == sizeof(vacuum_Energy));

	G4MaterialPropertiesTable* sClad_mt = new G4MaterialPropertiesTable();
	sClad_mt->AddProperty("RINDEX", vacuum_Energy, sClad_RIND, vacnum);
	fSClad->SetMaterialPropertiesTable(sClad_mt);
}


G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
	fCheckOverlaps = true;
	
	/*
		Out cylinder
	*/
	r 	= 3*cm; 		// actual path of the particle
	N_out	= 15;
	theta_out 	= 2*M_PI / N_out; 
	theta_scint_out = std::atan((fScintSizeZ_out/2) / (r - fScintSizeX_out/2)) *2;
	G4cout<<theta_out<<G4endl;
	G4cout<<N_out<<G4endl;
	
	if(N_out * theta_out >  2*M_PI) {N_out = N_out-1; theta_out = 2*M_PI / N_out;}
	while(theta_scint_out > theta_out) {
		G4cout<<"Too tight: "<<"N = "<<N_out<<" theta [deg] = "<<theta_out*180/M_PI<<" theta scint [deg] = "<<theta_scint_out*180/M_PI<<G4endl;
		N_out = N_out-1; theta_out = 2*M_PI / N_out;
	}
	G4cout<<"N = "<<N_out<<" theta [deg] = "<<theta_out*180/M_PI<<" theta scint [deg] = "<<theta_scint_out*180/M_PI<<G4endl;

	/*
		In cylinder
	*/
	
	r_in 	= 2*cm; 		// actual path of the particle
	N_in	= 15;
	theta_in 	= 2*M_PI / N_in; 
	theta_scint_in = std::atan((fScintSizeZ_in/2) / (r_in - fScintSizeX_in/2)) *2;
	G4cout<<theta_in<<G4endl;
	G4cout<<N_in<<G4endl;

	if(N_in * theta_in >  2*M_PI) {N_in = N_in-1; theta_in = 2*M_PI / N_in;}
	while(theta_scint_in > theta_in) {
		G4cout<<"Too tight: "<<"N = "<<N_in<<" theta [deg] = "<<theta_in*180/M_PI<<" theta scint [deg] = "<<theta_scint_in*180/M_PI<<G4endl;
		N_in = N_in-1; theta_in = 2*M_PI / N_in;
	}
	G4cout<<"N = "<<N_in<<" theta [deg] = "<<theta_in*180/M_PI<<" theta scint [deg] = "<<theta_scint_in*180/M_PI<<G4endl;

	// World Solid (size) -> Logical (material) -> PVPLacement (posiz, rotaz, to interact)
	fSolidWorld	= new G4Box("World", 0.5*fWorldSizeX, 0.5*fWorldSizeY, 0.5*fWorldSizeZ);
    fLogicWorld = new G4LogicalVolume(fSolidWorld, fVacuum_nogamma, "World");
    fPhysWorld	= new G4PVPlacement(0, G4ThreeVector(), fLogicWorld, "World", 0, false, 0, fCheckOverlaps);

	/*
		I will implement the two cylinder at the same time	
	*/
	//? Scintillator Solid (size) -> Logical (material) -> PVPLacement (posiz, rotaz, to interact)
	// Scintillator out
	fSolidScint_out	= new G4Box("Scint_out", 0.5*fScintSizeX_out, 0.5*fScintSizeY_out, 0.5*fScintSizeZ_out);
    fLogicScint_out = new G4LogicalVolume(fSolidScint_out, fScintMaterial, "Scint_out");
	fLogicScint_out->SetVisAttributes(G4Colour(0, 1, 0, 0.2));

	// Scintillator in
	fSolidScint_in	= new G4Box("Scint_in", 0.5*fScintSizeX_in, 0.5*fScintSizeY_in, 0.5*fScintSizeZ_in);
   	fLogicScint_in = new G4LogicalVolume(fSolidScint_in, fScintMaterial, "Scint_in");
	fLogicScint_in->SetVisAttributes(G4Colour(0.0, 1.0, 1.0, 0.3));

	//? Read Solid and Phys. (same size)
	fReadSizeX = 1.3*mm;
	fReadSizeY = 0.4*mm;
	fReadSizeZ = 1.3*mm;
	
	// out
	fSolidRead_out	= new G4Box("Read_out", 0.5*fReadSizeX,0.5*fReadSizeY, 0.5*fReadSizeZ);
    fLogicRead_out = new G4LogicalVolume(fSolidRead_out, fVacuum, "Read_out");
	fLogicRead_out->SetVisAttributes(G4Colour(0.0, 1.0, 1.0, 0.3));

	// in
	fSolidRead_in	= new G4Box("Read_in", 0.5*fReadSizeX,0.5*fReadSizeY, 0.5*fReadSizeZ);
    fLogicRead_in = new G4LogicalVolume(fSolidRead_in, fVacuum, "Read_in");
	fLogicRead_in->SetVisAttributes(G4Colour(1.0, 0.0, 1.0, 0.3));

	// grease Solid and Phys.  (same for both)
	fSolidGrease = new G4Box("Read", 0.5*fReadSizeX,0.5*0.5*fReadSizeY, 0.5*fReadSizeZ);
    fLogicGrease = new G4LogicalVolume(fSolidGrease, fOG, "Grease");
	fLogicGrease->SetVisAttributes(G4Colour(1,0,0, 0.5));

	//? VirtualDetector/SiPM Solid and Phys. 
	// out
	fSolidSiPM_out	= new G4Box("SiPM_out", 0.5*fReadSizeX, 0.5*0.5*fReadSizeY, 0.5*fReadSizeZ);
    fLogicSiPM_out 	= new G4LogicalVolume(fSolidSiPM_out, fSiPMMaterial, "SiPM_out");
    fLogicSiPM_out	->SetVisAttributes(G4Colour(0,0,1, 0.5));

	// in
	fSolidSiPM_in	= new G4Box("SiPM_in", 0.5*fReadSizeX, 0.5*0.5*fReadSizeY, 0.5*fReadSizeZ);
    fLogicSiPM_in 	= new G4LogicalVolume(fSolidSiPM_in, fSiPMMaterial, "SiPM_in");
    fLogicSiPM_in	->SetVisAttributes(G4Colour(0,1,0, 0.5));

	//? Put Grease and SiPM in Read
    G4ThreeVector Grease_pos = G4ThreeVector(0, -(0.5*fReadSizeY*0.5), 0); 
    G4ThreeVector SiPM_pos = G4ThreeVector(0, 0.5*fReadSizeY*0.5,0);
	
	// fPhysRead_in 	= new G4PVPlacement(0, G4ThreeVector(0, 0, 0), fLogicRead, "Read", fLogicWorld, true, 0, fCheckOverlaps);
	fPhysGrease	= new G4PVPlacement(0, Grease_pos, fLogicGrease, "Grease", fLogicRead_in, false, fCheckOverlaps);
	fPhysSiPM_in		= new G4PVPlacement(0, SiPM_pos, fLogicSiPM_in, "SiPM_in", fLogicRead_in, false, fCheckOverlaps);

	// fPhysRead_out 	= new G4PVPlacement(0, G4ThreeVector(0, 0, 0), fLogicRead, "Read", fLogicWorld, true, 0, fCheckOverlaps);
	fPhysGrease	= new G4PVPlacement(0, Grease_pos, fLogicGrease, "Grease", fLogicRead_out, false, fCheckOverlaps);
	fPhysSiPM_out		= new G4PVPlacement(0, SiPM_pos, fLogicSiPM_out, "SiPM_out", fLogicRead_out, false, fCheckOverlaps);

	G4Rotate3D 	flip_sipm =  G4Rotate3D(180*deg, G4ThreeVector(1, 0, 0));

	//? Elements to contain in and out system
	// Element out
	fElementSizeX_out = std::max(fScintSizeX_out, fReadSizeX);
	fElementSizeY_out = fScintSizeY_out + 1*mm;
	fElementSizeZ_out = std::max(fScintSizeZ_out, fReadSizeZ);

	fSolidElement_out = new G4Box("Element_out", 0.5*(fElementSizeX_out),0.5*(fElementSizeY_out), 0.5*(fElementSizeZ_out));
    fLogicElement_out = new G4LogicalVolume(fSolidElement_out, fVacuum, "Element_out");
	fLogicElement_out->SetVisAttributes(G4Colour(0, 1, 0, 0.2));

	// Element in
	fElementSizeX_in = std::max(fScintSizeX_in, fReadSizeX);
	fElementSizeY_in = fScintSizeY_in + 1*mm;
	fElementSizeZ_in = std::max(fScintSizeZ_in, fReadSizeZ);

	fSolidElement_in = new G4Box("Element_in", 0.5*(fElementSizeX_in),0.5*(fElementSizeY_in), 0.5*(fElementSizeZ_in));
    fLogicElement_in = new G4LogicalVolume(fSolidElement_in, fVacuum, "Element_in");
	fLogicElement_in->SetVisAttributes(G4Colour(0, 1, 0, 0.2));

	//? Put Scint and readout in the Element
	// Out
	G4int howmanySiPM_out = 5;
	G4ThreeVector Scint_pos_out = G4ThreeVector(0, 0, 0); 
	G4ThreeVector Read_pos_out = G4ThreeVector(0, 0.5*fReadSizeY + 0.5*fScintSizeY_out, 0); 
    G4Transform3D Read_transform_out;
	G4ThreeVector Read_translate_out = G4ThreeVector(0, 0, 0.5*fScintSizeZ_out/((howmanySiPM_out+1)*0.5));

	// fPhysElement_out 	= new G4PVPlacement(0, G4ThreeVector(-r, 0, 0), fLogicElement_out, "Element_out", fLogicWorld, true, 0, fCheckOverlaps);			
	fPhysScint_out		= new G4PVPlacement(0, Scint_pos_out, fLogicScint_out, "Scint_out", fLogicElement_out, false, fCheckOverlaps);
	
	// In
	G4int howmanySiPM_in = 3;
	G4ThreeVector Scint_pos_in = G4ThreeVector(0, 0, 0); 
	G4ThreeVector Read_pos_in = G4ThreeVector(0,0.5*fReadSizeY + 0.5*fScintSizeY_in, 0); 
    G4Transform3D Read_transform_in;
	G4ThreeVector Read_translate_in = G4ThreeVector(0, 0, 0.5*fScintSizeZ_in/((howmanySiPM_in+1)*0.5));

	// fPhysElement_in 	= new G4PVPlacement(0, G4ThreeVector(-r, 0, 0), fLogicElement_in, "Element_in", fLogicWorld, true, 0, fCheckOverlaps);			
	fPhysScint_in		= new G4PVPlacement(0, Scint_pos_in, fLogicScint_in, "Scint_in", fLogicElement_in, false, fCheckOverlaps);
	
	//? Loop on the number of SiPM you want. put them equally spaced on top and bottm
	for(int j=0; j<howmanySiPM_out; j += 1){
		if(j == 0) Read_transform_out =  (G4Translate3D)Read_pos_out;
		else if(j<(howmanySiPM_out+1)/2) Read_transform_out = (G4Translate3D)(j*Read_translate_out) * (G4Translate3D)Read_pos_out;
		else if(j>=(howmanySiPM_out+1)/2) Read_transform_out = (G4Translate3D)(-(howmanySiPM_out-j)*Read_translate_out) * (G4Translate3D)Read_pos_out;
		fPhysRead_out 			= new G4PVPlacement(Read_transform_out, fLogicRead_out, "Read_out", fLogicElement_out, true, j, fCheckOverlaps);
		Read_transform_out  = flip_sipm * Read_transform_out;
		fPhysRead_out 			= new G4PVPlacement(Read_transform_out, fLogicRead_out, "Read_out", fLogicElement_out, true, howmanySiPM_out+j, fCheckOverlaps);
	}

	for(int j=0; j<howmanySiPM_in; j += 1){
		if(j == 0) Read_transform_in =  (G4Translate3D)Read_pos_in;
		else if(j<(howmanySiPM_in+1)/2) Read_transform_in = (G4Translate3D)(j*Read_translate_in) * (G4Translate3D)Read_pos_in;
		else if(j>=(howmanySiPM_in+1)/2) Read_transform_in = (G4Translate3D)(-(howmanySiPM_in-j)*Read_translate_in) * (G4Translate3D)Read_pos_in;
		fPhysRead_in 			= new G4PVPlacement(Read_transform_in, fLogicRead_in, "Read_in", fLogicElement_in, true, j, fCheckOverlaps);
		Read_transform_in  = flip_sipm * Read_transform_in;
		fPhysRead_in 			= new G4PVPlacement(Read_transform_in, fLogicRead_in, "Read_in", fLogicElement_in, true, howmanySiPM_in+j, fCheckOverlaps);
	}

	// fPhysElement_out = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), fLogicElement_out, "Element_out", fLogicWorld, false, fCheckOverlaps);			
	// fPhysElement_in = new G4PVPlacement(0, G4ThreeVector(1*cm, 0, 0), fLogicElement_in, "Element_in", fLogicWorld, false, fCheckOverlaps);			

	for(int j=0; j<N_out; j += 1)
	{
		G4Rotate3D 	  rotation_out =  G4Rotate3D(j*theta_out*rad, G4ThreeVector(0, 1, 0)); //i*theta*deg std::cos(theta*i)
		G4Translate3D translate_out =  G4Translate3D(G4ThreeVector(-r, 0, 0));
		G4Transform3D transform_out= rotation_out*translate_out;
		fPhysElement_out = new G4PVPlacement(transform_out, fLogicElement_out, "Element_out", fLogicWorld, true, j, fCheckOverlaps);			
	}

	for(int j=0; j<N_in; j += 1)
	{
		G4Rotate3D 	  rotation_in =  G4Rotate3D(j*theta_in*rad, G4ThreeVector(0, 1, 0)); //i*theta*deg std::cos(theta*i)
		G4Translate3D translate_in =  G4Translate3D(G4ThreeVector(-r_in, 0, 0));
		G4Transform3D transform_in = rotation_in*translate_in;
		fPhysElement_in = new G4PVPlacement(transform_in, fLogicElement_in, "Element_in", fLogicWorld, true, j, fCheckOverlaps);			
	}
		
	// /*
		// In
	// */


	// // Element in
	// fElementSizeX_in = fScintSizeX_in + 1*mm;
	// fElementSizeY_in = fScintSizeY_in + 1*mm;
	// fElementSizeZ_in = fScintSizeZ_in + 1*mm;
	
	// fSolidElement_in = new G4Box("Element_in", 0.5*(fElementSizeX_in),0.5*(fElementSizeY_in), 0.5*(fElementSizeZ_in));
   	// fLogicElement_in = new G4LogicalVolume(fSolidElement_in, fVacuum, "Element_in");
	// fLogicElement_in->SetVisAttributes(G4Colour(0, 1, 0, 0.2));

	// fReadSizeX_in = fScintSizeX_in;
	// fReadSizeY_in = 0.5*mm;
	// fReadSizeZ_in = fScintSizeZ_in;

	// G4ThreeVector Scint_pos_in = G4ThreeVector(0, 0, 0); 
   	// G4ThreeVector Read_pos_in = G4ThreeVector(0, 0.5*fReadSizeY_in+0.5*fScintSizeY_in, 0);
    	
	// for(int j=1; j<N_in; j += 1)
	// {
		// G4Rotate3D 	  rotation_in =  G4Rotate3D(j*theta_in*rad, G4ThreeVector(0, 1, 0)); //i*theta*deg std::cos(theta*i)
		// G4Translate3D translate_in =  G4Translate3D(G4ThreeVector(-r_in, 0, 0));
		// G4Transform3D transform_in = rotation_in*translate_in;
		// fPhysElement_in = new G4PVPlacement(transform_in, fLogicElement_in, "Element_in", fLogicWorld, true, j, fCheckOverlaps);			
	// }
	

	// G4double fGround;
	// fGround =  0.999;
	// if(fGround < 1){
		G4OpticalSurface* OpScintSurface = new G4OpticalSurface("OpScintSurface");
		OpScintSurface->SetModel(glisur);
		OpScintSurface->SetType(dielectric_dielectric);
		OpScintSurface->SetFinish(polished); //ground polished
	// 	OpScintSurface->SetPolish(fGround);

	// 	new G4LogicalSkinSurface("ScintSurface", fLogicScint, OpScintSurface);	
	// 	if(!fmuEDM) new G4LogicalSkinSurface("ScintSurface", fLogicScint2, OpScintSurface);	
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
	// 	if(!fmuEDM) new G4LogicalSkinSurface("ScintSurface", fLogicScint2, OpScintSurface);	
	// }

	// Set how the volumes are visualized
    fLogicWorld	->SetVisAttributes(G4Colour(1, 1, 1, 0.1));
    if(fLogicScint_out) fLogicScint_out	->SetVisAttributes(G4Colour(0.34, 0.57, 0.8, 0.5));
    if(fLogicScint_in) fLogicScint_in	->SetVisAttributes(G4Colour(0.34, 0.57, 0.8, 0.5));


	if(fVDOn){
		// VirtualDetector Solid (size) -> Logical (material) -> PVPLacement (posiz, rotaz, to interact)
		fSolidVD	= new G4Box("VD", 0.5*fVDSizeX, 0.5*fVDSizeY, 0.5*fVDSizeZ);
    	fLogicVD 	= new G4LogicalVolume(fSolidVD, fVacuum, "VD");
    	fPhysVD		= new G4PVPlacement(0, G4ThreeVector(0., 0., 5*cm), fLogicVD, "VD", fLogicWorld, false, 0, fCheckOverlaps);

		// 2_VirtualDetector Solid (size) -> Logical (material) -> PVPLacement (posiz, rotaz, to interact)
		fSolidVD_2	= new G4Box("VD", 0.5*fVDSizeX, 0.5*fVDSizeY, 0.5*fVDSizeZ);
    	fLogicVD_2 	= new G4LogicalVolume(fSolidVD_2, fVacuum, "VD");
    	fPhysVD_2	= new G4PVPlacement(0, G4ThreeVector(0., 0., 1*cm), fLogicVD_2, "VD", fLogicWorld, false, 1, fCheckOverlaps);
    	fLogicVD	->SetVisAttributes(G4Colour(0.8, 0.34, 0.68, 0.05));
    	fLogicVD_2	->SetVisAttributes(G4Colour(0.8, 0.34, 0.68, 0.05));
	}
	
    return fPhysWorld;
}

G4VPhysicalVolume* DetectorConstruction::DefineVolumes_MuEDM()
{

	fCheckOverlaps = true;
	
	/*
		Out = Logitudinal SciFi
		In  = Transverse SciFi 
	*/
	
	G4double fFiberWidth  = 0.25*mm;
    G4double fFiberLength = 20*cm;
    G4double fFiberInLength = 3.2*cm;

	r 	= 3.14*cm; 		// actual path of the particle

	r = r+1.5*cm;
	N_out	= 8;
	theta_out 	= 2*M_PI / N_out; 
	G4cout<<theta_out<<G4endl;
	G4cout<<N_out<<G4endl;

	r_in = r-5*mm;
	N_in	= N_out; 	// Or 1 transvers SciFi every X longitudinal
	theta_in 	= 2*M_PI / N_in; 
	G4cout<<theta_in<<G4endl;
	G4cout<<N_in<<G4endl;

	// World Solid (size) -> Logical (material) -> PVPLacement (posiz, rotaz, to interact)
	fSolidWorld	= new G4Box("World", 0.5*fWorldSizeX, 0.5*fWorldSizeY, 0.5*fWorldSizeZ);
    fLogicWorld = new G4LogicalVolume(fSolidWorld, fVacuum_nogamma, "World");
    fPhysWorld	= new G4PVPlacement(0, G4ThreeVector(), fLogicWorld, "World", 0, false, 0, fCheckOverlaps);

	//? Read Solid and Phys. (same size)
	fReadSizeX = 1.3*mm;
	fReadSizeY = 0.4*mm;
	fReadSizeZ = 3.2*cm;
	
	// out
	fSolidRead_out	= new G4Box("Read_out", 0.5*fReadSizeX,0.5*fReadSizeY, 0.5*fReadSizeZ);
    fLogicRead_out = new G4LogicalVolume(fSolidRead_out, fVacuum, "Read_out");
	fLogicRead_out->SetVisAttributes(G4Colour(0.0, 1.0, 1.0, 0.3));

	// in
	fSolidRead_in	= new G4Box("Read_in", 0.5*fReadSizeX,0.5*fReadSizeY, 0.5*fReadSizeZ);
    fLogicRead_in = new G4LogicalVolume(fSolidRead_in, fVacuum, "Read_in");
	fLogicRead_in->SetVisAttributes(G4Colour(1.0, 0.0, 1.0, 0.3));

	// grease Solid and Phys.  (same for both)
	fSolidGrease = new G4Box("Read", 0.5*fReadSizeX,0.5*0.5*fReadSizeY, 0.5*fReadSizeZ);
    fLogicGrease = new G4LogicalVolume(fSolidGrease, fOG, "Grease");
	fLogicGrease->SetVisAttributes(G4Colour(1,0,0, 0.5));

	//? VirtualDetector/SiPM Solid and Phys. 
	// out
	fSolidSiPM_out	= new G4Box("SiPM_out", 0.5*fReadSizeX, 0.5*0.5*fReadSizeY, 0.5*fReadSizeZ);
    fLogicSiPM_out 	= new G4LogicalVolume(fSolidSiPM_out, fSiPMMaterial, "SiPM_out");
    fLogicSiPM_out	->SetVisAttributes(G4Colour(0,0,1, 0.5));

	// in
	fSolidSiPM_in	= new G4Box("SiPM_in", 0.5*fReadSizeX, 0.5*0.5*fReadSizeY, 0.5*fReadSizeZ);
    fLogicSiPM_in 	= new G4LogicalVolume(fSolidSiPM_in, fSiPMMaterial, "SiPM_in");
    fLogicSiPM_in	->SetVisAttributes(G4Colour(0,1,0, 0.5));

	//? Put Grease and SiPM in Read
    G4ThreeVector Grease_pos = G4ThreeVector(0, -(0.5*fReadSizeY*0.5), 0); 
    G4ThreeVector SiPM_pos = G4ThreeVector(0, 0.5*fReadSizeY*0.5,0);
	
	// fPhysRead_in 	= new G4PVPlacement(0, G4ThreeVector(0, 0, 0), fLogicRead, "Read", fLogicWorld, true, 0, fCheckOverlaps);
	fPhysGrease	= new G4PVPlacement(0, Grease_pos, fLogicGrease, "Grease", fLogicRead_in, false, fCheckOverlaps);
	fPhysSiPM_in		= new G4PVPlacement(0, SiPM_pos, fLogicSiPM_in, "SiPM_in", fLogicRead_in, false, fCheckOverlaps);

	// fPhysRead_out 	= new G4PVPlacement(0, G4ThreeVector(0, 0, 0), fLogicRead, "Read", fLogicWorld, true, 0, fCheckOverlaps);
	fPhysGrease	= new G4PVPlacement(0, Grease_pos, fLogicGrease, "Grease", fLogicRead_out, false, fCheckOverlaps);
	fPhysSiPM_out		= new G4PVPlacement(0, SiPM_pos, fLogicSiPM_out, "SiPM_out", fLogicRead_out, false, fCheckOverlaps);

	G4Rotate3D 	flip_sipm =  G4Rotate3D(180*deg, G4ThreeVector(1, 0, 0));

	//? Elements to contain in and out system
	// Element out
	fElementSizeX_out = fReadSizeX;
	fElementSizeY_out = fFiberLength + 1*mm;
	fElementSizeZ_out = fReadSizeZ;
	
	fSolidElement_out = new G4Box("Element_out", 0.5*(fElementSizeX_out),0.5*(fElementSizeY_out), 0.5*(fElementSizeZ_out));
    fLogicElement_out = new G4LogicalVolume(fSolidElement_out, fVacuum, "Element_out");
	fLogicElement_out->SetVisAttributes(G4Colour(0, 1, 0, 0.05));

	// Element in
	fElementSizeX_in = fReadSizeX;
	fElementSizeY_in = fFiberInLength + 1*mm;
	fElementSizeZ_in = fReadSizeZ;

	fSolidElement_in = new G4Box("Element_in", 0.5*(fElementSizeX_in),0.5*(fElementSizeY_in), 0.5*(fElementSizeZ_in));
    fLogicElement_in = new G4LogicalVolume(fSolidElement_in, fVacuum, "Element_in");
	fLogicElement_in->SetVisAttributes(G4Colour(0, 1, 0, 0.05));


//! OUT
	G4Material* fMaterial = fBCF20; 
	G4Box* solidElement = new G4Box("Element", 0.5 * fFiberWidth, 0.5 * fFiberLength, 0.5 * fFiberWidth);

	// Fiber
	G4Box* solidFiber = new G4Box("Fiber", 0.5 * fFiberWidth,
				               0.5 * fFiberLength,
					       0.5 * fFiberWidth);
	G4LogicalVolume* logicFiber = new G4LogicalVolume(solidFiber, fVacuum, "Fiber");
	
	// Core
	G4Box* solidCore = new G4Box("Core", 0.5 * fFiberWidth * 0.94,
				              0.5 * fFiberLength,
					      0.5 * fFiberWidth * 0.94);
	G4Box* solidCoreLong = new G4Box("CoreLong", 0.5 * fFiberWidth * 0.94,
				              0.5 * fFiberLength + 0.01*mm,
					      0.5 * fFiberWidth * 0.94);
	G4LogicalVolume* logicCore = new G4LogicalVolume(solidCore, fMaterial, "Core");
	fLogicFiber = logicCore;

	// First Cladding
	G4Box* temp = new G4Box("Temp",  0.5 * fFiberWidth * 0.98,
			                 0.5 * fFiberLength,
				         0.5 * fFiberWidth * 0.98);
	G4Box* tempLong = new G4Box("TempLong",  0.5 * fFiberWidth * 0.98,
			                 0.5 * fFiberLength + 0.01*mm,
				         0.5 * fFiberWidth * 0.98);
	G4SubtractionSolid* solidClad = new G4SubtractionSolid("fClad", temp, solidCoreLong, 0, G4ThreeVector());
	G4LogicalVolume* logicFClad = new G4LogicalVolume(solidClad, fFClad, "fClad");

	// Second Cladding
	G4SubtractionSolid* solidSClad = new G4SubtractionSolid("sClad", solidFiber, tempLong, 0, G4ThreeVector());
	G4LogicalVolume* logicSClad = new G4LogicalVolume(solidSClad, fSClad, "sClad");

	G4VPhysicalVolume* physCore = new G4PVPlacement(0,G4ThreeVector(), logicCore, "Core", logicFiber, false, 0, fCheckOverlaps);
	G4VPhysicalVolume* physFClad = new G4PVPlacement(0,G4ThreeVector(), logicFClad, "fClad", logicFiber, false, 0, fCheckOverlaps);
	G4VPhysicalVolume* physSClad = new G4PVPlacement(0,G4ThreeVector(), logicSClad, "sClad", logicFiber, false, 0, fCheckOverlaps);

	logicCore->SetVisAttributes(G4Colour(1, 1, 1, 0.2));
	logicSClad->SetVisAttributes(G4Colour(1, 0, 1, 0.05));
	logicFClad->SetVisAttributes(G4Colour(0, 1, 0, 0.2));
	logicFiber->SetVisAttributes(G4Colour(0, 0, 1, 0.0));

	// Surface
	// Core Surface
	
	G4OpticalSurface* OpCoreSurface = new G4OpticalSurface("CoreSurface");
	OpCoreSurface->SetModel(glisur);
	OpCoreSurface->SetType(dielectric_dielectric);
	OpCoreSurface->SetFinish(ground);
	OpCoreSurface->SetPolish(0.985);

	new G4LogicalBorderSurface("CoreSurface", physCore, physFClad, OpCoreSurface);

	// First Cladding Surface
	G4OpticalSurface* OpFCladSurface = new G4OpticalSurface("FCladSurface");
	OpFCladSurface->SetModel(glisur);
	OpFCladSurface->SetType(dielectric_dielectric);
	OpFCladSurface->SetFinish(ground);
	OpFCladSurface->SetPolish(0.98);

	new G4LogicalBorderSurface("FCladSurface", physFClad, physSClad, OpFCladSurface);

	// Second Cladding Surface
	G4OpticalSurface* OpSCladSurface = new G4OpticalSurface("SCladSurface");
	OpSCladSurface->SetModel(glisur);
	OpSCladSurface->SetType(dielectric_dielectric);
	OpSCladSurface->SetFinish(ground);
	OpSCladSurface->SetPolish(0.5);

	new G4LogicalBorderSurface("SCladSurface", physSClad, fPhysWorld, OpSCladSurface);
	
//!

//! IN
	G4Box* solidElementIn = new G4Box("Element", 0.5 * fFiberWidth, 0.5 * fFiberInLength, 0.5 * fFiberWidth);

	// Fiber
	G4Box* solidFiberIn = new G4Box("FiberIn", 0.5 * fFiberWidth,
				               0.5 * fFiberInLength,
					       0.5 * fFiberWidth);
	G4LogicalVolume* logicFiberIn = new G4LogicalVolume(solidFiberIn, fVacuum, "Fiber");
	
	// Core
	G4Box* solidCoreIn = new G4Box("CoreIn", 0.5 * fFiberWidth * 0.94,
				              0.5 * fFiberInLength,
					      0.5 * fFiberWidth * 0.94);
	G4Box* solidCoreLongIn = new G4Box("CoreLongIn", 0.5 * fFiberWidth * 0.94,
				              0.5 * fFiberInLength + 0.01*mm,
					      0.5 * fFiberWidth * 0.94);
	G4LogicalVolume* logicCoreIn = new G4LogicalVolume(solidCoreIn, fMaterial, "CoreIn");
	fLogicFiberIn = logicCoreIn;

	// First Cladding
	G4Box* tempIn = new G4Box("TempIn",  0.5 * fFiberWidth * 0.98,
			                 0.5 * fFiberInLength,
				         0.5 * fFiberWidth * 0.98);
	G4Box* tempLongIn = new G4Box("TempLong",  0.5 * fFiberWidth * 0.98,
			                 0.5 * fFiberInLength + 0.01*mm,
				         0.5 * fFiberWidth * 0.98);
	G4SubtractionSolid* solidCladIn = new G4SubtractionSolid("fCladIn", tempIn, solidCoreLongIn, 0, G4ThreeVector());
	G4LogicalVolume* logicFCladIn = new G4LogicalVolume(solidCladIn, fFClad, "fCladIn");

	// Second Cladding
	G4SubtractionSolid* solidSCladIn = new G4SubtractionSolid("sCladIn", solidFiberIn, tempLongIn, 0, G4ThreeVector());
	G4LogicalVolume* logicSCladIn = new G4LogicalVolume(solidSCladIn, fSClad, "sCladIn");

	G4VPhysicalVolume* physCoreIn = new G4PVPlacement(0,G4ThreeVector(), logicCoreIn, "CoreIn", logicFiberIn, false, 0, fCheckOverlaps);
	G4VPhysicalVolume* physFCladIn = new G4PVPlacement(0,G4ThreeVector(), logicFCladIn, "fCladIn", logicFiberIn, false, 0, fCheckOverlaps);
	G4VPhysicalVolume* physSCladIn = new G4PVPlacement(0,G4ThreeVector(), logicSCladIn, "sCladIn", logicFiberIn, false, 0, fCheckOverlaps);

	logicCoreIn->SetVisAttributes(G4Colour(1, 1, 1, 0.2));
	logicSCladIn->SetVisAttributes(G4Colour(1, 0, 1, 0.05));
	logicFCladIn->SetVisAttributes(G4Colour(0, 0, 1, 0.2));
	logicFiberIn->SetVisAttributes(G4Colour(0, 0, 1, 0.0));

	// Surface
	// Core Surface
	new G4LogicalBorderSurface("CoreSurfaceIn", physCoreIn, physFCladIn, OpCoreSurface);
	
	// First Cladding Surface
	new G4LogicalBorderSurface("FCladSurfaceIn", physFCladIn, physSCladIn, OpFCladSurface);

	// Second Cladding Surface
	new G4LogicalBorderSurface("SCladSurfaceIn", physSCladIn, fPhysWorld, OpSCladSurface);
	
//!

	//? Put Fiber and readout in the Element
	// Out
	G4int howmanySciFi = fElementSizeZ_out/fFiberWidth;
	G4cout<<"howmanySciFi"<<fElementSizeZ_out<<"/"<<fFiberWidth<<"="<<howmanySciFi<<G4endl;
	G4ThreeVector Read_pos_out = G4ThreeVector(0, 0.5*fReadSizeY + 0.5*fFiberLength, 0); 
    G4Transform3D Read_transform_out;
	G4ThreeVector Read_translate_out = G4ThreeVector(0, 0, fFiberWidth);
	G4ThreeVector offset = G4ThreeVector(0, 0, 0.5*fFiberWidth);
	G4ThreeVector offset_layer = G4ThreeVector(fFiberWidth, 0, 0);

	// fPhysElement_out 	= new G4PVPlacement(0, G4ThreeVector(-r, 0, 0), fLogicElement_out, "Element_out", fLogicWorld, true, 0, fCheckOverlaps);			
	for(int j=0; j<howmanySciFi-1; j += 1){
		Read_transform_out = (G4Translate3D)(j*Read_translate_out - 0.5*howmanySciFi*Read_translate_out+offset);
		new G4PVPlacement(Read_transform_out, logicFiber, "Fiber", fLogicElement_out, true , 0*howmanySciFi+j, fCheckOverlaps);
		Read_transform_out = (G4Translate3D)(offset_layer+offset)*Read_transform_out;
		new G4PVPlacement(Read_transform_out, logicFiber, "Fiber", fLogicElement_out, true , 1*howmanySciFi+j, fCheckOverlaps);
		Read_transform_out = (G4Translate3D)(-2*offset_layer)*Read_transform_out;
		new G4PVPlacement(Read_transform_out, logicFiber, "Fiber", fLogicElement_out, true , 2*howmanySciFi+j, fCheckOverlaps);
	}
	

	//?Put SiPM you want. on top and bottm
	Read_transform_out =  (G4Translate3D)Read_pos_out;
	fPhysRead_out 			= new G4PVPlacement(Read_transform_out, fLogicRead_out, "Read_out", fLogicElement_out, true, 0, fCheckOverlaps);
	Read_transform_out  = flip_sipm * Read_transform_out;
	fPhysRead_out 			= new G4PVPlacement(Read_transform_out, fLogicRead_out, "Read_out", fLogicElement_out, true, 1, fCheckOverlaps);

	for(int j=0; j<N_out; j += 1)
	{
		G4Rotate3D 	  rotation_out =  G4Rotate3D(j*theta_out*rad, G4ThreeVector(0, 1, 0)); //i*theta*deg std::cos(theta*i)
		G4Translate3D translate_out =  G4Translate3D(G4ThreeVector(-r, 0, 0));
		G4Transform3D transform_out= rotation_out*translate_out;
		fPhysElement_out = new G4PVPlacement(transform_out, fLogicElement_out, "Element_out", fLogicWorld, true, j, fCheckOverlaps);			
	}

	//? Put Fiber and readout in the Element
	// In
	howmanySciFi = fElementSizeZ_in/fFiberWidth;
	G4cout<<"howmanySciFi"<<fElementSizeZ_in<<"/"<<fFiberWidth<<"="<<howmanySciFi<<G4endl;
	G4ThreeVector Read_pos_in = G4ThreeVector(0, 0.5*fReadSizeY + 0.5*fFiberInLength, 0); 
    G4Transform3D Read_transform_in;
	G4ThreeVector Read_translate_in = G4ThreeVector(0, 0, fFiberWidth);

	// fPhysElement_out 	= new G4PVPlacement(0, G4ThreeVector(-r, 0, 0), fLogicElement_out, "Element_out", fLogicWorld, true, 0, fCheckOverlaps);			
	for(int j=0; j<howmanySciFi-1; j += 1){
		Read_transform_in = (G4Translate3D)(j*Read_translate_in - 0.5*howmanySciFi*Read_translate_in+offset);
		new G4PVPlacement(Read_transform_in, logicFiberIn, "FiberIn", fLogicElement_in, true , 0*howmanySciFi+j, fCheckOverlaps);
		Read_transform_in = (G4Translate3D)(offset_layer+offset)*Read_transform_in;
		new G4PVPlacement(Read_transform_in, logicFiberIn, "FiberIn", fLogicElement_in, true , 1*howmanySciFi+j, fCheckOverlaps);
		Read_transform_in = (G4Translate3D)(-2*offset_layer)*Read_transform_in;
		new G4PVPlacement(Read_transform_in, logicFiberIn, "FiberIn", fLogicElement_in, true , 2*howmanySciFi+j, fCheckOverlaps);
	}
	

	//?Put SiPM you want. on top and bottm
	Read_transform_in =  (G4Translate3D)Read_pos_in;
	fPhysRead_in 			= new G4PVPlacement(Read_transform_in, fLogicRead_in, "Read_in", fLogicElement_in, true, 0, fCheckOverlaps);
	Read_transform_in  = flip_sipm * Read_transform_in;
	fPhysRead_in 			= new G4PVPlacement(Read_transform_in, fLogicRead_in, "Read_in", fLogicElement_in, true, 1, fCheckOverlaps);


	//? How many stacks?
	for(int k=0; k<3; k += 1){
		for(int j=0; j<N_in; j += 1)
		{
			G4Rotate3D 	  sideway =  G4Rotate3D(90*deg, G4ThreeVector(1, 0, 0)); //i*theta*deg std::cos(theta*i)
			G4Rotate3D 	  rotation_in =  G4Rotate3D(j*theta_in*rad, G4ThreeVector(0, 1, 0)); //i*theta*deg std::cos(theta*i)
			G4Translate3D translate_in =  G4Translate3D(G4ThreeVector(-r_in, k*fElementSizeZ_in-1*fElementSizeZ_in, 0));
			if(k==0) translate_in = G4Translate3D(G4ThreeVector(-r_in, k*fElementSizeZ_in-1*fElementSizeZ_in -1*mm, 0));;
			if(k==2) translate_in = G4Translate3D(G4ThreeVector(-r_in, k*fElementSizeZ_in-1*fElementSizeZ_in +1*mm, 0));;
			G4Transform3D transform_in = rotation_in*translate_in*sideway;
			fPhysElement_in = new G4PVPlacement(transform_in, fLogicElement_in, "Element_in", fLogicWorld, true, j, fCheckOverlaps);			
		}
	}
	

	// Set how the volumes are visualized
    fLogicWorld	->SetVisAttributes(G4Colour(1, 1, 1, 0.1));
	
    return fPhysWorld;
}

G4VPhysicalVolume* DetectorConstruction::DefineVolumes_SciFi()
{

	fCheckOverlaps = true;
	
	G4double fFiberWidth  = 0.25*mm;
	// G4double fFiberWidth  = 20*mm;
    G4double fFiberLength = 2*cm;

	// World Solid (size) -> Logical (material) -> PVPLacement (posiz, rotaz, to interact)
	fSolidWorld	= new G4Box("World", 0.5*fWorldSizeX, 0.5*fWorldSizeY, 0.5*fWorldSizeZ);
    fLogicWorld = new G4LogicalVolume(fSolidWorld, fVacuum_nogamma, "World");
    fPhysWorld	= new G4PVPlacement(0, G4ThreeVector(), fLogicWorld, "World", 0, false, 0, fCheckOverlaps);

	//? Read Solid and Phys. (same size)
	fReadSizeX = 1.3*mm;
	fReadSizeY = 0.4*mm;
	fReadSizeZ = 3*cm;
	
	G4double fSiPMFrameSizeX = 1.3*mm;
	G4double fSiPMFrameSizeZ = 0.25*mm;

	// fSiPMSizeX = fSiPMFrameSizeX;
	// fSiPMSizeY = fReadSizeY;
	// fSiPMSizeZ = fSiPMFrameSizeZ;

	// G4double window_positionX = 0*mm;
	// G4double window_positionZ = 0*mm;

	fSiPMSizeX = 1.3*mm;
	fSiPMSizeY = fReadSizeY;
	fSiPMSizeZ = 0.23*mm;

	G4double window_positionX = 0;//-0.5*fSiPMFrameSizeX+0.5*fSiPMSizeX+0.1*mm;
	G4double window_positionZ = 0;//-0.5*fSiPMFrameSizeZ+0.5*fSiPMSizeZ+0.1*mm;


	// out
	fSolidRead_out	= new G4Box("Read_out", 0.5*fReadSizeX,0.5*fReadSizeY, 0.5*fReadSizeZ);
    fLogicRead_out = new G4LogicalVolume(fSolidRead_out, fSi, "Read_out"); //vacuum?
	fLogicRead_out->SetVisAttributes(G4Colour(0,0,1, 0.4));

	// grease Solid and Phys.  (same for both)
	fSolidGrease = new G4Box("Read", 0.5*fReadSizeX,0.5*0.5*fReadSizeY, 0.5*fReadSizeZ);
    fLogicGrease = new G4LogicalVolume(fSolidGrease, fOG, "Grease");
	fLogicGrease->SetVisAttributes(G4Colour(1,0,0, 0.5));

	//? VirtualDetector/SiPM Solid and Phys. 
	// out
	fSolidSiPM_out	= new G4Box("SiPM_out", 0.5*fSiPMSizeX, 0.5*0.5*fSiPMSizeY, 0.5*fSiPMSizeZ);
    fLogicSiPM_out 	= new G4LogicalVolume(fSolidSiPM_out, fSi, "SiPM_out");
    fLogicSiPM_out	->SetVisAttributes(G4Colour(0,0,1, 0.5));

	//? Put Grease and SiPM in Read
    G4ThreeVector Grease_pos = G4ThreeVector(0, -(0.5*fReadSizeY*0.5), 0); 
    G4ThreeVector SiPM_pos = G4ThreeVector(0, 0.5*fReadSizeY*0.5,0);
	
	fPhysGrease	= new G4PVPlacement(0, Grease_pos, fLogicGrease, "Grease", fLogicRead_out, false, fCheckOverlaps);
	
	// how many SiPM 
	G4int howmanySiPM = fReadSizeZ/fSiPMFrameSizeZ;
	G4cout<<"howmanySiPM"<<fElementSizeZ_out<<"/"<<fReadSizeZ<<"="<<howmanySiPM<<G4endl;
    G4Transform3D SiPM_transform_out;
	G4ThreeVector SiPM_translate_out = G4ThreeVector(0, 0, fSiPMFrameSizeZ);
	G4ThreeVector SiPM_Window_translate_out = G4ThreeVector(window_positionX, 0, window_positionZ);
	// G4ThreeVector SiPM_offset = G4ThreeVector(0, 0, fSiPMSizeX);
	// G4ThreeVector offset_layer = G4ThreeVector(fFiberWidth, 0, 0);

	for(int j=0; j<howmanySiPM; j += 1){
		SiPM_transform_out = (G4Translate3D)(SiPM_pos+j*SiPM_translate_out-(howmanySiPM-1)*0.5*SiPM_translate_out + SiPM_Window_translate_out);
		new G4PVPlacement(SiPM_transform_out, fLogicSiPM_out, "SiPM_out", fLogicRead_out, true, j, fCheckOverlaps);
	}


	G4Rotate3D 	flip_sipm =  G4Rotate3D(180*deg, G4ThreeVector(0, 0, 1));

	//? Elements to contain in and out system
	// Element out
	fElementSizeX_out = fReadSizeX;
	fElementSizeY_out = fFiberLength + 1*mm;
	fElementSizeZ_out = fReadSizeZ;

	fSolidElement_out = new G4Box("Element_out", 0.5*(fElementSizeX_out),0.5*(fElementSizeY_out), 0.5*(fElementSizeZ_out));
    fLogicElement_out = new G4LogicalVolume(fSolidElement_out, fVacuum, "Element_out");
	fLogicElement_out->SetVisAttributes(G4Colour(0, 1, 0, 0.05));


	
	G4Material* fMaterial = fBCF20; 
	G4Box* solidElement = new G4Box("Element", 0.5 * fFiberWidth, 0.5 * fFiberLength, 0.5 * fFiberWidth);

	// Fiber
	G4Box* solidFiber = new G4Box("Fiber", 0.5 * fFiberWidth,
				               0.5 * fFiberLength,
					       0.5 * fFiberWidth);
	G4LogicalVolume* logicFiber = new G4LogicalVolume(solidFiber, fVacuum, "Fiber");
	
	// Core
	G4Box* solidCore = new G4Box("Core", 0.5 * fFiberWidth * 0.94,
				              0.5 * fFiberLength,
					      0.5 * fFiberWidth * 0.94);
	G4Box* solidCoreLong = new G4Box("CoreLong", 0.5 * fFiberWidth * 0.94,
				              0.5 * fFiberLength + 0.01*mm,
					      0.5 * fFiberWidth * 0.94);
	G4LogicalVolume* logicCore = new G4LogicalVolume(solidCore, fMaterial, "Core");
	fLogicFiber = logicCore;

	// First Cladding
	G4Box* temp = new G4Box("Temp",  0.5 * fFiberWidth * 0.98,
			                 0.5 * fFiberLength,
				         0.5 * fFiberWidth * 0.98);
	G4Box* tempLong = new G4Box("TempLong",  0.5 * fFiberWidth * 0.98,
			                 0.5 * fFiberLength + 0.01*mm,
				         0.5 * fFiberWidth * 0.98);
	G4SubtractionSolid* solidClad = new G4SubtractionSolid("fClad", temp, solidCoreLong, 0, G4ThreeVector());
	G4LogicalVolume* logicFClad = new G4LogicalVolume(solidClad, fFClad, "fClad");

	// Second Cladding
	G4SubtractionSolid* solidSClad = new G4SubtractionSolid("sClad", solidFiber, tempLong, 0, G4ThreeVector());
	G4LogicalVolume* logicSClad = new G4LogicalVolume(solidSClad, fSClad, "sClad");

	G4VPhysicalVolume* physCore = new G4PVPlacement(0,G4ThreeVector(), logicCore, "Core", logicFiber, false, 0, fCheckOverlaps);
	G4VPhysicalVolume* physFClad = new G4PVPlacement(0,G4ThreeVector(), logicFClad, "fClad", logicFiber, false, 0, fCheckOverlaps);
	G4VPhysicalVolume* physSClad = new G4PVPlacement(0,G4ThreeVector(), logicSClad, "sClad", logicFiber, false, 0, fCheckOverlaps);

	logicCore->SetVisAttributes(G4Colour(0, 0.5, 1, 1));
	logicFClad->SetVisAttributes(G4Colour(0.5, 0.1, 0.1, 1));
	logicSClad->SetVisAttributes(G4Colour(0.5, 0.5, 0.2, 1));
	logicFiber->SetVisAttributes(G4Colour(0, 0, 0, 0));

	// Surface
	// Core Surface
	
	G4OpticalSurface* OpCoreSurface = new G4OpticalSurface("CoreSurface");
	OpCoreSurface->SetModel(glisur);
	OpCoreSurface->SetType(dielectric_dielectric);
	OpCoreSurface->SetFinish(ground);
	OpCoreSurface->SetPolish(0.985);

	new G4LogicalBorderSurface("CoreSurface", physCore, physFClad, OpCoreSurface);

	// First Cladding Surface
	G4OpticalSurface* OpFCladSurface = new G4OpticalSurface("FCladSurface");
	OpFCladSurface->SetModel(glisur);
	OpFCladSurface->SetType(dielectric_dielectric);
	OpFCladSurface->SetFinish(ground);
	OpFCladSurface->SetPolish(0.98);

	new G4LogicalBorderSurface("FCladSurface", physFClad, physSClad, OpFCladSurface);

	// Second Cladding Surface
	G4OpticalSurface* OpSCladSurface = new G4OpticalSurface("SCladSurface");
	OpSCladSurface->SetModel(glisur);
	OpSCladSurface->SetType(dielectric_dielectric);
	OpSCladSurface->SetFinish(ground);
	OpSCladSurface->SetPolish(0.5);

	new G4LogicalBorderSurface("SCladSurface", physSClad, fPhysWorld, OpSCladSurface);
	

	new G4PVPlacement(0, G4ThreeVector() , logicFiber, "Fiber", fLogicWorld, false, fCheckOverlaps);

	//? Put Fiber and readout in the Element
	// Out
	G4int howmanySciFi = fElementSizeZ_out/fFiberWidth;
	G4cout<<"howmanySciFi"<<fElementSizeZ_out<<"/"<<fFiberWidth<<"="<<howmanySciFi<<G4endl;
	G4ThreeVector Read_pos_out = G4ThreeVector(0, 0.5*fReadSizeY + 0.5*fFiberLength, 0); 
    G4Transform3D Read_transform_out;
	G4ThreeVector Read_translate_out = G4ThreeVector(0, 0, fFiberWidth);
	G4ThreeVector offset = G4ThreeVector(0, 0, 0.5*fFiberWidth);
	G4ThreeVector offset_layer = G4ThreeVector(fFiberWidth, 0, 0);

	// fPhysElement_out 	= new G4PVPlacement(0, G4ThreeVector(-r, 0, 0), fLogicElement_out, "Element_out", fLogicWorld, true, 0, fCheckOverlaps);			
	for(int j=0; j<howmanySciFi-1; j += 1){
		Read_transform_out = (G4Translate3D)(j*Read_translate_out - 0.5*howmanySciFi*Read_translate_out+offset);
		new G4PVPlacement(Read_transform_out, logicFiber, "Fiber", fLogicElement_out, true , 0*howmanySciFi+j, fCheckOverlaps);
		Read_transform_out = (G4Translate3D)(offset_layer+offset)*Read_transform_out;
		new G4PVPlacement(Read_transform_out, logicFiber, "Fiber", fLogicElement_out, true , 1*howmanySciFi+j, fCheckOverlaps);
		Read_transform_out = (G4Translate3D)(-2*offset_layer)*Read_transform_out;
		new G4PVPlacement(Read_transform_out, logicFiber, "Fiber", fLogicElement_out, true , 2*howmanySciFi+j, fCheckOverlaps);
	}
	
	//?Put SiPM you want. on top and bottm
	Read_transform_out =  (G4Translate3D)Read_pos_out;
	fPhysRead_out 			= new G4PVPlacement(Read_transform_out, fLogicRead_out, "Read_out", fLogicElement_out, true, 0, fCheckOverlaps);
	Read_transform_out  = flip_sipm * Read_transform_out;
	fPhysRead_out 			= new G4PVPlacement(Read_transform_out, fLogicRead_out, "Read_out", fLogicElement_out, true, 1, fCheckOverlaps);
	
	G4Rotate3D 	  rotation_out =  G4Rotate3D(90*deg, G4ThreeVector(0, 1, 0)); //i*theta*deg std::cos(theta*i)
	G4Translate3D translate_out =  G4Translate3D(G4ThreeVector(0, 0, 0));
	G4Transform3D transform_out= rotation_out*translate_out;
	fPhysElement_out = new G4PVPlacement(transform_out, fLogicElement_out, "Element_out", fLogicWorld, false, fCheckOverlaps);			

	// Set how the volumes are visualized
    fLogicWorld	->SetVisAttributes(G4Colour(1, 1, 1, 0.1));

    return fPhysWorld;
}

void DetectorConstruction::ConstructSDandField()
{
	auto sdManager = G4SDManager::GetSDMpointer();
   	G4String SDname;
	
	// if(fLogicScint_out){
		// ScintSD* scint_SD_out = new ScintSD(SDname="Scint_out");
  		// sdManager->AddNewDetector(scint_SD_out);
		// fLogicScint_out->SetSensitiveDetector(scint_SD_out);
	// }

	if(fLogicSiPM_out){
		SiPMSD * SiPM_SD_out = new SiPMSD("SiPM_out");
  		sdManager->AddNewDetector(SiPM_SD_out);
		fLogicSiPM_out->SetSensitiveDetector(SiPM_SD_out);
	}

	// if(fLogicScint_in){
		// ScintSD* scint_SD_in = new ScintSD(SDname="Scint_in");
  		// sdManager->AddNewDetector(scint_SD_in);
		// fLogicScint_in->SetSensitiveDetector(scint_SD_in);
	// }

	if(fLogicSiPM_in){
		SiPMSD * SiPM_SD_in = new SiPMSD("SiPM_in");
  		sdManager->AddNewDetector(SiPM_SD_in);
		fLogicSiPM_in->SetSensitiveDetector(SiPM_SD_in);
	}

	// G4ThreeVector fieldValue(0.,-3*tesla,0.);
	// fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
	//fMagFieldMessenger->SetVerboseLevel(1);
	
	// Register the field messenger for deleting
	G4AutoDelete::Register(fMagFieldMessenger);
	
	if(fVDOn){
		// Create the Sensitive Detector defined in VirtualDetectorSD 
		VirtualDetectorSD * VD_SD = new VirtualDetectorSD("VirtualDetector");

		// Assign the SD to the logial volume
		fLogicVD->SetSensitiveDetector(VD_SD);
	}
}

/*
	From here functions used through the Messenger to modify the detector
*/

void DetectorConstruction :: SetScintSize(G4double size){
	// fScintSizeX = size;
	// fScintSizeY = size;
	// fScintSizeZ = size;

	// fWorldSizeX = std::max(fScintSizeX, fWorldSizeX);
	// fWorldSizeY = std::max(fScintSizeY, fWorldSizeY);
	// fWorldSizeZ = std::max(fScintSizeZ, fWorldSizeZ);

	// fSolidScint->SetXHalfLength(size*0.5);
	// fSolidScint->SetYHalfLength(size*0.5);
	// fSolidScint->SetZHalfLength(size*0.5);

	// fSolidWorld->SetXHalfLength(0.5*fWorldSizeX);
	// fSolidWorld->SetYHalfLength(0.5*fWorldSizeY);
	// fSolidWorld->SetZHalfLength(0.5*fWorldSizeZ);

	// //fPhysScint->SetTranslation(G4ThreeVector(0, 0, 20*cm));
	// G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

void DetectorConstruction :: SetScintSize(G4ThreeVector size){
	// fScintSizeX = size.getX();
	// fScintSizeY = size.getY();
	// fScintSizeZ = size.getZ();

	// fWorldSizeX = std::max(fScintSizeX, fWorldSizeX);
	// fWorldSizeY = std::max(fScintSizeY, fWorldSizeY);
	// fWorldSizeZ = std::max(fScintSizeZ, fWorldSizeZ);

	// fSolidScint->SetXHalfLength(size.getX()*0.5);
	// fSolidScint->SetYHalfLength(size.getY()*0.5);
	// fSolidScint->SetZHalfLength(size.getZ()*0.5);

	// fSolidWorld->SetXHalfLength(0.5*fWorldSizeX);
	// fSolidWorld->SetYHalfLength(0.5*fWorldSizeY);
	// fSolidWorld->SetZHalfLength(0.5*fWorldSizeZ);

	// fSolidWorld	= new G4Box("World", 0.5*std::max(fScintSizeX, fWorldSizeX), 
		// 0.5*std::max(fScintSizeY, fWorldSizeY), 
		// 0.5*std::max(fScintSizeZ, fWorldSizeZ));

	// //fPhysScint->SetTranslation(G4ThreeVector(0, 0, 20*cm));
	// G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

void DetectorConstruction::SetScintMaterial(G4String name){
	// if(name == "BC400") fScintMaterial = fBC400;
	// else if(name == "LYSO") fScintMaterial = fLYSO;
	// G4RunManager::GetRunManager()->ReinitializeGeometry();
}