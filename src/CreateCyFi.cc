/// \file  CyFi.hh
/// \brief Implementation of the class to create CyFi

#include "CreateCyFi.hh"

CreateCyFi::CreateCyFi()
{}

CreateCyFi :: ~CreateCyFi()
{}

void CreateCyFi::Create(G4LogicalVolume * hLogicWorld)
{
    Materials();
	OpticalProperties();
	Volumes(hLogicWorld);
	SD();
}

void CreateCyFi::Materials()
{
	G4cout<<G4endl;
	G4cout<<"-----------------------------------------"<<G4endl;
	G4cout<<"-------- CreateCyFi::Materials() --------"<<G4endl;

	G4double a; // atomic mass
	G4double z; // atomic number
	G4double density;

	//? You can take materials already defined 
	//G4NistManager* nist = G4NistManager::Instance();
	
	//? Or you can define your elements and materials
	//! Elements
	hH  = new G4Element( "hH",  "hH", z =  1., a =   1.01*g/mole);
	hC  = new G4Element( "hC",  "hC", z =  6., a =  12.01*g/mole);
	hN  = new G4Element( "hN",  "hN", z =  7., a =  14.01*g/mole);
	hO  = new G4Element( "hO",  "hO", z =  8., a =  16.00*g/mole);
	hSie= new G4Element("hSi", "hSi", z = 14., a = 28.0855*g/mole);
	hY  = new G4Element( "hY",  "hY", z = 39., a = 88.90585*g/mole);
	hCe = new G4Element("hCe", "hCe", z = 58., a = 140.116*g/mole);
	hLu = new G4Element("hLu", "hLu", z = 71., a = 174.967*g/mole);

	//! Materials
	// BC400
	hBC400 = new G4Material("hBC400", density = 1.023*g/cm3, 2);
	hBC400->AddElement(hC, 1000);
	hBC400->AddElement(hH, 1103);

	hBC400_noscint = new G4Material("hBC400_noscint", density = 1.023*g/cm3, 2);
	hBC400_noscint->AddElement(hC, 1000);
	hBC400_noscint->AddElement(hH, 1103);

	// LYSO
	hLYSO = new G4Material("hYSO", density = 7.1*g/cm3, 5);
	hLYSO->AddElement(hLu,  9);
	hLYSO->AddElement( hY, 10);
	hLYSO->AddElement(hSie, 5);
	hLYSO->AddElement( hO, 25);
	hLYSO->AddElement(hCe,  5);

	// Vacuuum
	hVacuum = new G4Material("hVacuum",z=1.,a=1.01*g/mole, 
		density = universe_mean_density, 
		kStateGas, 0.1 * kelvin, 1.e-19 * pascal);
	hVacuum_nogamma = new G4Material("hVacuum_nogamma",z=1.,a=1.01*g/mole, 
		density = universe_mean_density, 
		kStateGas, 0.1 * kelvin, 1.e-19 * pascal);

	// Air
	hAir = new G4Material("hAir", density = 0.0010*g/cm3, 2);
	hAir->AddElement(hN, 70 * perCent);
	hAir->AddElement(hO, 30 * perCent);

	// Optical grease
	hOG = new G4Material("hOpticalGrease",z=1.,a=1.01*g/mole, 
		     		 density = universe_mean_density, kStateGas,
				 0.1 * kelvin, 1.e-19 * pascal);

	// Silicium
	G4NistManager* NISTman = G4NistManager::Instance();
	hSi = NISTman->FindOrBuildMaterial("G4_Si");

	// Silicon resin
	hSiResin = new G4Material("hSiResin",z=1.,a=1.01*g/mole, 
		     		 density = universe_mean_density, kStateGas,
				 0.1 * kelvin, 1.e-19 * pascal);

	// Assign default materials
	hScintMaterial = hBC400;
	hSiPMMaterial  = hSiResin;

	//! SciFi
	/// Materials
	// BCF10
	hBCF10 = new G4Material("hBCF10", density = 1.05*g/cm3, 2);
	hBCF10->AddElement(hC, 485);
	hBCF10->AddElement(hH, 482);

	// BCF12
	hBCF12 = new G4Material("hBCF12", density = 1.05*g/cm3, 2);
	hBCF12->AddElement(hC, 485);
	hBCF12->AddElement(hH, 482);

	// BCF20
	hBCF20 = new G4Material("hBCF20", density = 1.05*g/cm3, 2);
	hBCF20->AddElement(hC, 485);
	hBCF20->AddElement(hH, 482);

	// First Cladding: PMMA
	hFClad = new G4Material("hFClad", density = 1.2*g/cm3, 3);
	hFClad->AddElement(hC, 5);
	hFClad->AddElement(hH, 8);
	hFClad->AddElement(hO, 2);

	// Second Cladding: PMMA EMA
	hSClad = new G4Material("hSClad", density = 1.2*g/cm3, 3);
	hSClad->AddElement(hC, 5);
	hSClad->AddElement(hH, 8);
	hSClad->AddElement(hO, 2);


	G4cout<<"CreateCyFi::Materials() : Done !!"<<G4endl;
	G4cout<<"-----------------------------------------"<<G4endl;
}

void CreateCyFi::OpticalProperties()
{
	G4cout<<G4endl;
	G4cout<<"-----------------------------------------"<<G4endl;
	G4cout<<"---- CreateCyFi::OpticalProperties() ----"<<G4endl;

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

	std::reverse(energy.begin(), energy.end());
  	std::reverse(scint.begin(), scint.end());

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

	hBC400_mt = new G4MaterialPropertiesTable();
	hBC400_mt->AddProperty(       "RINDEX", BC400_Energy,  BC400_RIND, bc400);
	hBC400_mt->AddProperty(    "ABSLENGTH", BC400_Energy,  BC400_ABSL, bc400);
	hBC400_mt->AddProperty("SCINTILLATIONCOMPONENT1", BC400_Energy, BC400_SCINT, bc400);
	
	hBC400_mt->AddConstProperty("SCINTILLATIONYIELD",        11050./MeV);
	hBC400_mt->AddConstProperty(   "RESOLUTIONSCALE",                 1);
	hBC400_mt->AddConstProperty(  "SCINTILLATIONTIMECONSTANT1",            2.4*ns);
	hBC400_mt->AddConstProperty(  "SCINTILLATIONRISETIME1",   0.9*ns);
	hBC400_mt->AddConstProperty(  "SCINTILLATIONTIMECONSTANT2",            2.4*ns);
	hBC400_mt->AddConstProperty(  "SCINTILLATIONRISETIME2",   0.9*ns);
	
	hBC400->SetMaterialPropertiesTable(hBC400_mt);

	//  Set Birks Constant
	hBC400->GetIonisation()->SetBirksConstant(0.117645*mm/MeV);

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

	std::reverse(energy.begin(), energy.end());
  	std::reverse(scint.begin(), scint.end());

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

	hLYSO_mt = new G4MaterialPropertiesTable();
	hLYSO_mt->AddProperty(       "RINDEX", LYSO_Energy,  LYSO_RIND, lyso);
	hLYSO_mt->AddProperty(    "ABSLENGTH", LYSO_Energy,  LYSO_ABSL, lyso);
	hLYSO_mt->AddProperty("SCINTILLATIONCOMPONENT1", LYSO_Energy, LYSO_SCINT, lyso);
	
	hLYSO_mt->AddConstProperty("SCINTILLATIONYIELD",        33200./MeV);
	hLYSO_mt->AddConstProperty(   "RESOLUTIONSCALE",                 1);
	hLYSO_mt->AddConstProperty(  "SCINTILLATIONTIMECONSTANT1",            36*ns);
	hLYSO_mt->AddConstProperty(  "SCINTILLATIONTIMECONSTANT2",            36*ns);
	
	hLYSO->SetMaterialPropertiesTable(hLYSO_mt);

	//  Set Birks Constant
	//! fLYSO->GetIonisation()->SetBirksConstant(0.117645*mm/MeV);

	// ----------------------------------------------------------	
	// Vacuum & Air
	// ----------------------------------------------------------	
	G4double vacuum_Energy[] = {1.5*eV, 4.*eV};
	const G4int vacnum = sizeof(vacuum_Energy) / sizeof(G4double);

	G4double vRIND = 1.;
	G4double vacuum_RIND[] = {vRIND, vRIND};
	assert(sizeof(vacuum_RIND) == sizeof(vacuum_Energy));

	G4MaterialPropertiesTable* cvacuum_mt = new G4MaterialPropertiesTable();
	cvacuum_mt->AddProperty("RINDEX", vacuum_Energy, vacuum_RIND, vacnum);
	hVacuum->SetMaterialPropertiesTable(cvacuum_mt);
	hAir   ->SetMaterialPropertiesTable(cvacuum_mt);

	// ----------------------------------------------------------	
	// Silicium
	// ----------------------------------------------------------	
	G4double Si_Energy[] = {.5*eV, 9.*eV};
	const G4int Sinum = sizeof(vacuum_Energy) / sizeof(G4double);

	G4double Si_RIND[] = {3.4, 3.4};
	assert(sizeof(Si_RIND) == sizeof(Si_Energy));

	G4MaterialPropertiesTable* cSi_mt = new G4MaterialPropertiesTable();
	cSi_mt->AddProperty("RINDEX", Si_Energy, Si_RIND, Sinum);
	hSi->SetMaterialPropertiesTable(cSi_mt);

	// ----------------------------------------------------------	
	// Silicon resin
	// ----------------------------------------------------------	
	G4double SiRes_RIND[] = {1.41, 1.41};
	assert(sizeof(SiRes_RIND) == sizeof(Si_Energy));
	
	G4MaterialPropertiesTable* cSiRes_mt = new G4MaterialPropertiesTable();
	cSiRes_mt->AddProperty("RINDEX", Si_Energy, SiRes_RIND, Sinum);
	hSiResin->SetMaterialPropertiesTable(cSiRes_mt);

	// ----------------------------------------------------------	
	// Optical grease
	// ----------------------------------------------------------	
	//? better if it was higher?
	G4double OG_RIND[] = {1.465, 1.465};
	assert(sizeof(OG_RIND) == sizeof(Si_Energy));
		
	G4MaterialPropertiesTable* cOG_mt = new G4MaterialPropertiesTable();
	cOG_mt->AddProperty("RINDEX", Si_Energy, OG_RIND, Sinum);
	hOG->SetMaterialPropertiesTable(cOG_mt);

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

	std::reverse(energy.begin(), energy.end());
  	std::reverse(scint.begin(), scint.end());

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

	hBCF10_mt = new G4MaterialPropertiesTable();
	hBCF10_mt->AddProperty(       "RINDEX", BCF10_Energy,  BCF10_RIND, bcf10);
	hBCF10_mt->AddProperty(    "ABSLENGTH", BCF10_Energy,  BCF10_ABSL, bcf10);
	hBCF10_mt->AddProperty("SCINTILLATIONCOMPONENT1", BCF10_Energy, BCF10_SCINT, bcf10);
	
	hBCF10_mt->AddConstProperty("SCINTILLATIONYIELD",        8000./MeV);
	hBCF10_mt->AddConstProperty(   "RESOLUTIONSCALE",                 1);
	hBCF10_mt->AddConstProperty(  "SCINTILLATIONTIMECONSTANT1",            2.7*ns);
	hBCF10_mt->AddConstProperty(  "SCINTILLATIONTIMECONSTANT2",            2.7*ns);
	
	hBCF10->SetMaterialPropertiesTable(hBCF10_mt);

	//  Set Birks Constant
	hBCF10->GetIonisation()->SetBirksConstant(0.117645*mm/MeV);

	//  BCF12 optics
	myfile.open("../tables/BCF12_light_out.txt");
	while(true){
		myfile >> tempe >> tempscint;
		energy.push_back(1239.84197/tempe);
		scint.push_back(tempscint);
		
		if(myfile.eof()) break;
	}
	myfile.close();

	std::reverse(energy.begin(), energy.end());
  	std::reverse(scint.begin(), scint.end());

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

	hBCF12_mt = new G4MaterialPropertiesTable();
	hBCF12_mt->AddProperty(       "RINDEX", BCF12_Energy,  BCF12_RIND, bcf12);
	hBCF12_mt->AddProperty(    "ABSLENGTH", BCF12_Energy,  BCF12_ABSL, bcf12);
	hBCF12_mt->AddProperty("SCINTILLATIONCOMPONENT1", BCF12_Energy, BCF12_SCINT, bcf12);
	
	hBCF12_mt->AddConstProperty("SCINTILLATIONYIELD",        8000./MeV);
	hBCF12_mt->AddConstProperty(   "RESOLUTIONSCALE",                 1);
	hBCF12_mt->AddConstProperty(  "SCINTILLATIONTIMECONSTANT1",            3.2*ns);
	hBCF12_mt->AddConstProperty(  "SCINTILLATIONTIMECONSTANT2",            3.2*ns);
	
	hBCF12->SetMaterialPropertiesTable(hBCF12_mt);

	//  Set Birks Constant
	hBCF12->GetIonisation()->SetBirksConstant(0.117645*mm/MeV);

	//  BCF20 optics
	myfile.open("../tables/BCF20_light_out.txt");
	while(true){
		myfile >> tempe >> tempscint;
		energy.push_back(1239.84197/tempe);
		scint.push_back(tempscint);
		
		if(myfile.eof()) break;
	}
	myfile.close();

	std::reverse(energy.begin(), energy.end());
  	std::reverse(scint.begin(), scint.end());

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

	hBCF20_mt = new G4MaterialPropertiesTable();
	hBCF20_mt->AddProperty(       "RINDEX", BCF20_Energy,  BCF20_RIND, bcf20);
	hBCF20_mt->AddProperty(    "ABSLENGTH", BCF20_Energy,  BCF20_ABSL, bcf20);
	hBCF20_mt->AddProperty("SCINTILLATIONCOMPONENT1", BCF20_Energy, BCF20_SCINT, bcf20);
	
	hBCF20_mt->AddConstProperty("SCINTILLATIONYIELD",        8000./MeV);
	hBCF20_mt->AddConstProperty(   "RESOLUTIONSCALE",                 1);
	hBCF20_mt->AddConstProperty(  "SCINTILLATIONTIMECONSTANT1",            2.7*ns);
	hBCF20_mt->AddConstProperty(  "SCINTILLATIONTIMECONSTANT2",            2.7*ns);
	
	hBCF20->SetMaterialPropertiesTable(hBCF20_mt);

	//  Set Birks Constant
	hBCF20->GetIonisation()->SetBirksConstant(0.117645*mm/MeV);

	// Cladding
	G4double fCladRIND = 1.49;
	G4double fClad_RIND[] = {fCladRIND, fCladRIND};
	assert(sizeof(fClad_RIND) == sizeof(vacuum_Energy));

	G4MaterialPropertiesTable* hClad_mt = new G4MaterialPropertiesTable();
	hClad_mt->AddProperty("RINDEX", vacuum_Energy, fClad_RIND, vacnum);
	hFClad->SetMaterialPropertiesTable(hClad_mt);


	G4double sCladRIND = 1.42;
	G4double sClad_RIND[] = {sCladRIND, sCladRIND};
	assert(sizeof(sClad_RIND) == sizeof(vacuum_Energy));

	G4MaterialPropertiesTable* hsClad_mt = new G4MaterialPropertiesTable();
	hsClad_mt->AddProperty("RINDEX", vacuum_Energy, sClad_RIND, vacnum);
	hSClad->SetMaterialPropertiesTable(hsClad_mt);

	G4cout<<"CreateCyFi::OpticalProperties() : Done !!"<<G4endl;
	G4cout<<"-----------------------------------------"<<G4endl;
}

void CreateCyFi::Volumes(G4LogicalVolume * hLogicWorld)
{
	G4cout<<G4endl;
	G4cout<<"-----------------------------------------"<<G4endl;
	G4cout<<"--------- CreateCyFi::Volumes() ---------"<<G4endl;

	hCyFi_length 	= 15*cm;
	hCyFi_radius 	= 3.5*cm;
	hFiberThickness	= 5*mm;

	bool bool_CyFiOpticalGrease = true; 
	hCheckOverlaps = false;

	//? CyFi
	G4VSolid* CyFiSolid = new G4Tubs("CyFiSolid",
                               hCyFi_radius-0.5*hFiberThickness - 1*mm,			// inner radius
                               hCyFi_radius + 1.5*hFiberThickness + 1*mm + 1*mm,	// outer radius
                               0.5*hCyFi_length+3*mm,							// height
                               0.0 * deg,  360.0 * deg);  						// segment angles  
    G4LogicalVolume* CyFiLogic;
    if(bool_CyFiOpticalGrease) CyFiLogic = new G4LogicalVolume(CyFiSolid, hOG, "CyFi"); //fVacuum_nogamma
    else CyFiLogic = new G4LogicalVolume(CyFiSolid, hVacuum, "CyFi"); //fVacuum_nogamma
	G4PVPlacement* CyFiPlacement = new G4PVPlacement(0, G4ThreeVector(0,0,0), CyFiLogic, "CyFiLogic", hLogicWorld, false, 0, hCheckOverlaps);			
	CyFiLogic->SetVisAttributes(G4Colour(1, 0.45, 0, 0.1));

	//? Read Solid and Phys.
	hReadSizeX = hFiberThickness;
	hReadSizeY = hFiberThickness;
	hReadSizeZ = 0.4*mm;

	hSiPMSizeX = hFiberThickness;
	hSiPMSizeY = hFiberThickness;
	hSiPMSizeZ = 0.5*hReadSizeZ;

	// in
	hSolidRead_in	= new G4Box("Read_in", 0.5*hReadSizeX, 0.5*hReadSizeY, 0.5*hReadSizeZ);
    hLogicRead_in = new G4LogicalVolume(hSolidRead_in, hSi, "Read_in");
	hLogicRead_in->SetVisAttributes(G4Colour(1.0, 0.0, 1.0, 0.3));

	// out
	hSolidRead_out	= new G4Box("Read_out", 0.5*hReadSizeX, 0.5*hReadSizeY, 0.5*hReadSizeZ);
    hLogicRead_out = new G4LogicalVolume(hSolidRead_out, hSi, "Read_out");
	hLogicRead_out->SetVisAttributes(G4Colour(0.0, 1.0, 1.0, 0.3));

	// grease Solid and Phys.  (same for both)
	hSolidGrease = new G4Box("Read", 0.5*hReadSizeX, 0.5*hReadSizeY,0.5 * 0.5*hReadSizeZ);
    hLogicGrease = new G4LogicalVolume(hSolidGrease, hOG, "Grease");
	hLogicGrease->SetVisAttributes(G4Colour(1,0,0, 0.5));

	// SiPM out
	hSolidSiPM_out	= new G4Box("SiPM_out", 0.5*hSiPMSizeX, 0.5*hSiPMSizeY, 0.5 *hSiPMSizeZ);
    hLogicSiPM_out 	= new G4LogicalVolume(hSolidSiPM_out, hSiPMMaterial, "SiPM_out");
    hLogicSiPM_out	->SetVisAttributes(G4Colour(0,0,1, 0.5));

	// SiPM in
	hSolidSiPM_in	= new G4Box("SiPM_in", 0.5*hSiPMSizeX, 0.5*hSiPMSizeY, 0.5 *hSiPMSizeZ);
    hLogicSiPM_in 	= new G4LogicalVolume(hSolidSiPM_in, hSiPMMaterial, "SiPM_in");
    hLogicSiPM_in	->SetVisAttributes(G4Colour(0,0,1, 0.5));

	//? Put Grease and SiPM in Read
    G4ThreeVector Grease_pos = G4ThreeVector(0, 0, (0.5*0.5*hReadSizeZ)); 
    G4ThreeVector SiPM_pos = G4ThreeVector(0, 0, -0.5*0.5*hReadSizeZ);
	
	hPhysRead_in 	= new G4PVPlacement(0, G4ThreeVector(0, 0, 0), hLogicRead_in, "Read", hLogicWorld, true, 0, hCheckOverlaps);
	hPhysGrease 	= new G4PVPlacement(0, Grease_pos , hLogicGrease, "Grease", hLogicRead_in, false, hCheckOverlaps);
	hPhysSiPM_in 	= new G4PVPlacement(0, SiPM_pos, hLogicSiPM_in, "SiPM", hLogicRead_in, false, hCheckOverlaps);

	//hPhysRead_out 	= new G4PVPlacement(0, G4ThreeVector(10, 0, 0), hLogicRead_in, "Read", hLogicWorld, true, 0, hCheckOverlaps);
	hPhysGrease 	= new G4PVPlacement(0, Grease_pos , hLogicGrease, "Grease", hLogicRead_out, false, hCheckOverlaps);
	hPhysSiPM_out 	= new G4PVPlacement(0, SiPM_pos, hLogicSiPM_out, "SiPM", hLogicRead_out, false, hCheckOverlaps);


	G4Rotate3D 	flip_sipm =  G4Rotate3D(180*deg, G4ThreeVector(0, 1, 0));

	G4Material* fMaterial = hBCF20; 

	//? Elements to contain in and out system

	/*
		IN
	*/
	G4cout<<"\nInner fibers"<<G4endl;
	double r_in = hCyFi_radius;
	TVector3 center_in = TVector3(r_in,0,0);
	double angle_in = 60*deg;
	double turns_in = AngleToTurns(angle_in, hCyFi_length, r_in);
	double lenght_in = sqrt( pow(turns_in * 2*M_PI*r_in, 2) + pow(hCyFi_length, 2) ); 
	G4cout<<"Angle_in : "<<angle_in<<"; Turns_in : "<<turns_in<<"; Lenght_in : "<<lenght_in<<G4endl;
	double steps_in = 200;
	// Fiber made of fVacuum
	// G4TessellatedSolid* CreateHelix(G4String name, TVector3 center, double size, double runningangle, doule length, int steps, double extrusion)
	G4TessellatedSolid* FiberSolid_in = CreateHelix("TriangularHelix_in", center_in, hFiberThickness, angle_in, hCyFi_length, steps_in, hReadSizeZ);	
    G4LogicalVolume* FiberLogical_in = new G4LogicalVolume(FiberSolid_in, hVacuum, "HelixLogical_in");
    
	// Core made of fBCF20 and long version for cut-out
	G4TessellatedSolid* Core_long_in = CreateHelix("Core_in", center_in, 0.94*hFiberThickness, angle_in, hCyFi_length, steps_in, 5);
    G4LogicalVolume* CoreLogical_long_in = new G4LogicalVolume(Core_long_in, fMaterial, "HelixLogical_in");
	G4TessellatedSolid* Core_in = CreateHelix("Core_in", center_in, 0.94*hFiberThickness, angle_in, hCyFi_length, steps_in, 0);
    G4LogicalVolume* CoreLogical_in = new G4LogicalVolume(Core_in, fMaterial, "HelixLogical_in");

	// FirstCladding made of fFClad and long version for cut-out
	G4TessellatedSolid* FirstCladding_long_in = CreateHelix("FirstCladding_long_in", center_in, 0.98*hFiberThickness, angle_in, hCyFi_length, steps_in, 5);
    G4LogicalVolume* FirstCladdingLogical_long_in = new G4LogicalVolume(FirstCladding_long_in, hFClad, "HelixLogical_in");
	G4TessellatedSolid* tmpFirstCladding_in = CreateHelix("FirstCladding_in", center_in, 0.98*hFiberThickness, angle_in, hCyFi_length, steps_in, 0);
    G4LogicalVolume* tmpFirstCladdingLogical_in = new G4LogicalVolume(tmpFirstCladding_in, hFClad, "HelixLogical_in");

	// SecondCladding made of hSClad and long version for cut-out
	G4TessellatedSolid* tmpSecondCladding_in = CreateHelix("SecondCladding_in", center_in, hFiberThickness, angle_in, hCyFi_length, steps_in, 0);	
    G4LogicalVolume* tmpSecondCladdingLogical_in = new G4LogicalVolume(tmpSecondCladding_in, hSClad, "HelixLogical_in");

	// First Cladding: bool subtraction of temp FirstCladding and Core (long version)
	G4SubtractionSolid* FirstCladding_in = new G4SubtractionSolid("fClad", tmpFirstCladding_in, Core_long_in, 0, G4ThreeVector());
	G4LogicalVolume* FirstCladdingLogical_in = new G4LogicalVolume(FirstCladding_in, hFClad, "fClad");
	
	// Second Cladding: bool subtraction of temp SecondCladding and FirstCladding (long version)
	G4SubtractionSolid* SecondClad_in = new G4SubtractionSolid("sClad_in", tmpSecondCladding_in, FirstCladding_long_in, 0, G4ThreeVector());
	G4LogicalVolume* SecondCladdingLogical_in = new G4LogicalVolume(SecondClad_in, hFClad, "sClad_in");

	G4PVPlacement* fiberPlacement_in;
	// G4PVPlacement* fiberPlacement_in = new G4PVPlacement(0, G4ThreeVector(), FiberLogical_in, "HelixLogical_in", hLogicWorld, false, 0, hCheckOverlaps);			
	G4PVPlacement* FirstCladdingPlacement_in = new G4PVPlacement(0, G4ThreeVector(), FirstCladdingLogical_in, "fCladding_in", FiberLogical_in, false, 0, hCheckOverlaps);			
	G4PVPlacement* SecondCladdingPlacement_in = new G4PVPlacement(0, G4ThreeVector(), SecondCladdingLogical_in, "sCladding_in", FiberLogical_in, false, 0, hCheckOverlaps);			
	G4PVPlacement* CorePlacement_in = new G4PVPlacement(0, G4ThreeVector(), CoreLogical_in, "CoreLogical_in", FiberLogical_in, false, 0, hCheckOverlaps);			

	//? Surfaces
	// Core Surface
	G4OpticalSurface* OpCoreSurface_in = new G4OpticalSurface("CoreSurface_in");
	OpCoreSurface_in->SetModel(glisur);
	OpCoreSurface_in->SetType(dielectric_dielectric);
	OpCoreSurface_in->SetFinish(ground);
	OpCoreSurface_in->SetPolish(0.985);
	new G4LogicalBorderSurface("CoreSurface", CorePlacement_in, FirstCladdingPlacement_in, OpCoreSurface_in);

	// First Cladding Surface
	G4OpticalSurface* OpFCladSurface_in = new G4OpticalSurface("FCladSurface_in");
	OpFCladSurface_in->SetModel(glisur);
	OpFCladSurface_in->SetType(dielectric_dielectric);
	OpFCladSurface_in->SetFinish(ground);
	OpFCladSurface_in->SetPolish(0.98);
	new G4LogicalBorderSurface("FCladSurface", FirstCladdingPlacement_in, SecondCladdingPlacement_in, OpFCladSurface_in);

	// Second Cladding Surface
	G4OpticalSurface* OpSCladSurface_in = new G4OpticalSurface("SCladSurface_in");
	OpSCladSurface_in->SetModel(glisur);
	OpSCladSurface_in->SetType(dielectric_dielectric);
	OpSCladSurface_in->SetFinish(ground);
	OpSCladSurface_in->SetPolish(0.5);
	new G4LogicalBorderSurface("SCladSurface", SecondCladdingPlacement_in, CyFiPlacement, OpSCladSurface_in); //CyFiPlacement OR NOT?

	//? Add readout
	G4Transform3D Read_transform_in;
	TVector3 path_in = Path(0, turns_in, center_in, hCyFi_length);
	G4ThreeVector Read_translate_in = G4ThreeVector(path_in.x(), path_in.y(), path_in.z());
	G4Rotate3D 	  Read_adapt_in =  G4Rotate3D(-angle_in, G4ThreeVector(1, 0, 0));
	// Additional vector to traslat the Read along the orthogonal to the fiber  
	G4ThreeVector Read_adjust_in = G4ThreeVector(0, 0, 1);
	Read_adjust_in.rotateX(-angle_in);
	Read_adjust_in.setMag(-0.5*hReadSizeZ);
	Read_transform_in =  ((G4Translate3D)Read_adjust_in)*(G4Translate3D)Read_translate_in*Read_adapt_in;
	hPhysRead_in 	= new G4PVPlacement(Read_transform_in, hLogicRead_in, "Read", FiberLogical_in, true, 0, hCheckOverlaps);
	// Other side
	Read_adjust_in.rotateZ(2*M_PI * turns_in);
	Read_adjust_in.setMag(-0.5*hReadSizeZ);
	G4Rotate3D 	  Read_rotate_in =  G4Rotate3D(2*M_PI * turns_in, G4ThreeVector(0, 0, 1));
	// Additional vector to traslat the Read along the orthogonal to the fiber  
	Read_transform_in = ((G4Translate3D)Read_adjust_in)*(G4Translate3D)G4ThreeVector(0, 0, hCyFi_length)*Read_rotate_in*(G4Translate3D)Read_translate_in*Read_adapt_in*flip_sipm;
	hPhysRead_in 	= new G4PVPlacement(Read_transform_in, hLogicRead_in, "Read", FiberLogical_in, true, 1, hCheckOverlaps);


	int n_in = r_in * 2*M_PI / (1.5*hFiberThickness/cos(angle_in));
	n_in = n_in-1;
	G4cout<<"Fiber's arc : "<<r_in * 2*M_PI / hFiberThickness<<"; N of fibers : "<<n_in<<G4endl;
	double theta_in = 2*M_PI/n_in;
	G4cout<<"Rotation angle : "<<theta_in<<G4endl;
	for(int j=0; j<n_in; j += 1){
		G4Rotate3D 	  rotation =  G4Rotate3D(j*theta_in, G4ThreeVector(0, 0, 1));
		G4Transform3D transform = rotation;
		fiberPlacement_in = new G4PVPlacement((G4Translate3D)G4ThreeVector(0,0,-0.5*hCyFi_length)*transform, FiberLogical_in, "HelixLogical_in", CyFiLogic, true, j, hCheckOverlaps);
		if(hCheckOverlaps) print_progress_bar((double)j/(n_in-1), "HelixLogical_in");
	}

	/*
		OUT
	*/
	G4cout<<"\nOuter fibers"<<G4endl;
	double r_out = hCyFi_radius + hFiberThickness + 1*mm;
	TVector3 center_out = TVector3(r_out,0,0);
	double angle_out = 15*deg;
	double turns_out = AngleToTurns(angle_out, hCyFi_length, r_out);
	double lenght_out = sqrt( pow(turns_out * 2*M_PI*r_out, 2) + pow(hCyFi_length, 2) ); 
	G4cout<<"Angle_out : "<<angle_out<<"; Turns_out : "<<turns_out<<"; Lenght_out : "<<lenght_out<<G4endl;
	double steps_out = 200;
	// Fiber made of fVacuum
	// G4TessellatedSolid* CreateHelix(G4String name, TVector3 center, double size, double runningangle, doule length, int steps, double extrusion)
	G4TessellatedSolid* FiberSolid_out = CreateHelix("TriangularHelix_out", center_out, hFiberThickness, angle_out, hCyFi_length, steps_out, hReadSizeZ);	
    G4LogicalVolume* FiberLogical_out = new G4LogicalVolume(FiberSolid_out, hVacuum, "HelixLogical_out");
    
	// Core made of fBCF20 and long version for cut-out
	G4TessellatedSolid* Core_long_out = CreateHelix("Core_out", center_out, 0.94*hFiberThickness, angle_out, hCyFi_length, steps_out, 5);
    G4LogicalVolume* CoreLogical_long_out = new G4LogicalVolume(Core_long_out, fMaterial, "HelixLogical_out");
	G4TessellatedSolid* Core_out = CreateHelix("Core_out", center_out, 0.94*hFiberThickness, angle_out, hCyFi_length, steps_out, 0);
    G4LogicalVolume* CoreLogical_out = new G4LogicalVolume(Core_out, fMaterial, "HelixLogical_out");

	// FirstCladding made of fFClad and long version for cut-out
	G4TessellatedSolid* FirstCladding_long_out = CreateHelix("FirstCladding_long_out", center_out, 0.98*hFiberThickness, angle_out, hCyFi_length, steps_out, 5);
    G4LogicalVolume* FirstCladdingLogical_long_out = new G4LogicalVolume(FirstCladding_long_out, hFClad, "HelixLogical_out");
	G4TessellatedSolid* tmpFirstCladding_out = CreateHelix("FirstCladding_out", center_out, 0.98*hFiberThickness, angle_out, hCyFi_length, steps_out, 0);
    G4LogicalVolume* tmpFirstCladdingLogical_out = new G4LogicalVolume(tmpFirstCladding_out, hFClad, "HelixLogical_out");

	// SecondCladding made of fSClad and long version for cut-out
	G4TessellatedSolid* tmpSecondCladding_out = CreateHelix("SecondCladding_out", center_out, hFiberThickness, angle_out, hCyFi_length, steps_out, 0);	
    G4LogicalVolume* tmpSecondCladdingLogical_out = new G4LogicalVolume(tmpSecondCladding_out, hSClad, "HelixLogical_out");

	// First Cladding: bool subtraction of temp FirstCladding and Core (long version)
	G4SubtractionSolid* FirstCladding_out = new G4SubtractionSolid("fClad", tmpFirstCladding_out, Core_long_out, 0, G4ThreeVector());
	G4LogicalVolume* FirstCladdingLogical_out = new G4LogicalVolume(FirstCladding_out, hFClad, "fClad");
	
	// Second Cladding: bool subtraction of temp SecondCladding and FirstCladding (long version)
	G4SubtractionSolid* SecondClad_out = new G4SubtractionSolid("sClad_out", tmpSecondCladding_out, FirstCladding_long_out, 0, G4ThreeVector());
	G4LogicalVolume* SecondCladdingLogical_out = new G4LogicalVolume(SecondClad_out, hSClad, "sClad_out");
	
	G4PVPlacement* fiberPlacement_out;
	// fiberPlacement_out = new G4PVPlacement(0, G4ThreeVector(), FiberLogical_out, "HelixLogical_out", hLogicWorld, false, 0, hCheckOverlaps);			
	G4PVPlacement* FirstCladdingPlacement_out = new G4PVPlacement(0, G4ThreeVector(), FirstCladdingLogical_out, "fCladding_out", FiberLogical_out, false, 0, hCheckOverlaps);			
	G4PVPlacement* SecondCladdingPlacement_out = new G4PVPlacement(0, G4ThreeVector(), SecondCladdingLogical_out, "sCladding_out", FiberLogical_out, false, 0, hCheckOverlaps);			
	G4PVPlacement* CorePlacement_out = new G4PVPlacement(0, G4ThreeVector(), CoreLogical_out, "CoreLogical_out", FiberLogical_out, false, 0, hCheckOverlaps);			

	//? Surfaces
	// Core Surface
	G4OpticalSurface* OpCoreSurface_out = new G4OpticalSurface("CoreSurface_out");
	OpCoreSurface_out->SetModel(glisur);
	OpCoreSurface_out->SetType(dielectric_dielectric);
	OpCoreSurface_out->SetFinish(ground);
	OpCoreSurface_out->SetPolish(0.985);
	new G4LogicalBorderSurface("CoreSurface", CorePlacement_out, FirstCladdingPlacement_out, OpCoreSurface_out);

	// First Cladding Surface
	G4OpticalSurface* OpFCladSurface_out = new G4OpticalSurface("FCladSurface_out");
	OpFCladSurface_out->SetModel(glisur);
	OpFCladSurface_out->SetType(dielectric_dielectric);
	OpFCladSurface_out->SetFinish(ground);
	OpFCladSurface_out->SetPolish(0.98);
	new G4LogicalBorderSurface("FCladSurface", FirstCladdingPlacement_out, SecondCladdingPlacement_out, OpFCladSurface_out);

	// Second Cladding Surface
	G4OpticalSurface* OpSCladSurface_out = new G4OpticalSurface("SCladSurface_out");
	OpSCladSurface_out->SetModel(glisur);
	OpSCladSurface_out->SetType(dielectric_dielectric);
	OpSCladSurface_out->SetFinish(ground);
	OpSCladSurface_out->SetPolish(0.5);
	new G4LogicalBorderSurface("SCladSurface", SecondCladdingPlacement_out, CyFiPlacement, OpSCladSurface_out); //CyFiPlacement OR NOT?

	//? Add readout
	G4Transform3D Read_transform_out;
	TVector3 path_out = Path(0, turns_out, center_out, hCyFi_length);
	G4ThreeVector Read_translate_out = G4ThreeVector(path_out.x(), path_out.y(), path_out.z());
	G4Rotate3D 	  Read_adapt_out =  G4Rotate3D(-angle_out, G4ThreeVector(1, 0, 0));
	// Additional vector to traslat the Read along the orthogonal to the fiber  
	G4ThreeVector Read_adjust_out = G4ThreeVector(0, 0, 1);
	Read_adjust_out.rotateX(-angle_out);
	Read_adjust_out.setMag(-0.5*hReadSizeZ);
	Read_transform_out =  ((G4Translate3D)Read_adjust_out)*(G4Translate3D)Read_translate_out*Read_adapt_out;
	hPhysRead_out 	= new G4PVPlacement(Read_transform_out, hLogicRead_out, "Read", FiberLogical_out, true, 0, hCheckOverlaps);
	// Other side
	Read_adjust_out.rotateZ(2*M_PI * turns_out);
	Read_adjust_out.setMag(-0.5*hReadSizeZ);
	G4Rotate3D 	  Read_rotate_out =  G4Rotate3D(2*M_PI * turns_out, G4ThreeVector(0, 0, 1));
	// Additional vector to traslat the Read along the orthogonal to the fiber  
	Read_transform_out = ((G4Translate3D)Read_adjust_out)*(G4Translate3D)G4ThreeVector(0, 0, hCyFi_length)*Read_rotate_out*(G4Translate3D)Read_translate_out*Read_adapt_out*flip_sipm;
	hPhysRead_out 	= new G4PVPlacement(Read_transform_out, hLogicRead_out, "Read", FiberLogical_out, true, 1, hCheckOverlaps);

	int n_out = r_out * 2*M_PI / (1.5*hFiberThickness/cos(angle_out));
	n_out = n_out-1;
	G4cout<<"Fiber's arc : "<<r_out * 2*M_PI / hFiberThickness<<"; N of fibers : "<<n_out<<G4endl;
	double theta_out = 2*M_PI/n_out;
	G4cout<<"Rotation angle : "<<theta_out<<G4endl;
	for(int j=0; j<n_out; j += 1){
		G4Rotate3D 	  rotation =  G4Rotate3D(j*theta_out, G4ThreeVector(0, 0, 1));
		G4Transform3D transform = rotation;
		fiberPlacement_out = new G4PVPlacement((G4Translate3D)G4ThreeVector(0,0,-0.5*hCyFi_length)*transform, FiberLogical_out, "HelixLogical_out", CyFiLogic, true, j, hCheckOverlaps);			
    	if(hCheckOverlaps) print_progress_bar((double)(j+1)/n_out, "HelixLogical_out");
	}

	//? Set VisualAttributes
	// IN - Blue
	FiberLogical_in->SetVisAttributes(G4Colour(1, 1, 1, 0.2));
	CoreLogical_in->SetVisAttributes(G4Colour(0, 0, 1, 0.8));
	FirstCladdingLogical_in->SetVisAttributes(G4Colour(0, 1, 0, 0.6));
	SecondCladdingLogical_in->SetVisAttributes(G4Colour(1, 0, 1, 0.4));
	// OUT - Green
    FiberLogical_out->SetVisAttributes(G4Colour(1, 1, 1, 0.4));
	CoreLogical_out->SetVisAttributes(G4Colour(0, 1, 0, 0.6));
	FirstCladdingLogical_out->SetVisAttributes(G4Colour(0, 1, 0, 0.6));
	SecondCladdingLogical_out->SetVisAttributes(G4Colour(1, 0, 1, 0.4));


	G4cout<<"CreateCyFi::Volumes() : Done !!"<<G4endl;
	G4cout<<"-----------------------------------------"<<G4endl;
}

void CreateCyFi::SD()
{
	G4cout<<G4endl;
	G4cout<<"-----------------------------------------"<<G4endl;
	G4cout<<"------------ CreateCyFi::SD() -----------"<<G4endl;

	auto sdManager = G4SDManager::GetSDMpointer();
   	G4String SDname;

	if(hLogicSiPM_out){
		SiPMSD * SiPM_SD_out = new SiPMSD("SiPM_out");
  		sdManager->AddNewDetector(SiPM_SD_out);
		hLogicSiPM_out->SetSensitiveDetector(SiPM_SD_out);
	}

	if(hLogicSiPM_in){
		SiPMSD * SiPM_SD_in = new SiPMSD("SiPM_in");
  		sdManager->AddNewDetector(SiPM_SD_in);
		hLogicSiPM_in->SetSensitiveDetector(SiPM_SD_in);
	}

	G4cout<<"CreateCyFi::SD() : Done !!"<<G4endl;
	G4cout<<"-----------------------------------------"<<G4endl;
}