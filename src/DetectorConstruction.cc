/// \file  DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "RunAction.hh"
#include "SiPMModel.hh"


#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4RunManager.hh"

#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"

#include "G4Box.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4SDManager.hh"
#include "G4SDChargedFilter.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSFlatSurfaceCurrent.hh"
#include "G4PSNofSecondary.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSDoseDeposit.hh"
#include "G4PSTrackLength.hh"
#include "G4SDParticleFilter.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4GenericMessenger.hh"

/// Constructor
DetectorConstruction::DetectorConstruction() : 
	G4VUserDetectorConstruction(), fmodel("75PE"), fCrysVolume(nullptr), 
	fSiPMVolume(nullptr), fElementVolume(nullptr), fSolidCrys(nullptr), 
	fSolidWorld(nullptr), fSolidElement(nullptr), fCheckOverlaps(true), 
	fNbOfPixelsX(0), fNbOfPixelsY(0), fLogicPixel(nullptr), fLogicCrys(nullptr), 
	fCrysSizeX(2*mm), fCrysSizeY(2*mm), fCrysSizeZ(2*mm), fSiPM_sizeXY(1.3*mm), 
	fSiPM_sizeZ(0*mm), fGround(1), fAngle(0), fAngleWithOpticalGrease(0), fTilt(0)
{
	fDetectorMessenger = new DetectorMessenger(this);
	SetSiPMmodel("75PE");
	DefineMaterials();
	fMaterial = fBC400;
	fMaterialWindow = fEpResin;
}

/// Deconstructor
DetectorConstruction :: ~DetectorConstruction(){
	delete fDetectorMessenger;
}


/// Implementation

G4VPhysicalVolume* DetectorConstruction::Construct(){
	// Define volumes
	return DefineVolumes();
}

void DetectorConstruction::DefineMaterials(){
	G4NistManager* nist = G4NistManager::Instance();
	G4double a; // atomic mass
	G4double z; // atomic number
//	G4double R = 1.103; // hydrogen to carbon ratio in BC400
	G4double density;
	
	/// Elements
	fH  = new G4Element( "H",  "H", z =  1., a =   1.01*g/mole);
	fC  = new G4Element( "C",  "C", z =  6., a =  12.01*g/mole);
	fN  = new G4Element( "N",  "N", z =  7., a =  14.01*g/mole);
	fO  = new G4Element( "O",  "O", z =  8., a =  16.00*g/mole);
	fSie= new G4Element("Si", "Si", z = 14., a = 28.0855*g/mole);
	fY  = new G4Element( "Y",  "Y", z = 39., a = 88.90585*g/mole);
	fCe = new G4Element("Ce", "Ce", z = 58., a = 140.116*g/mole);
	fLu = new G4Element("Lu", "Lu", z = 71., a = 174.967*g/mole);
	
	/// Materials
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
	
	// Air
	fAir = new G4Material("Air", density = 0.0010*g/cm3, 2);
	fAir->AddElement(fN, 70 * perCent);
	fAir->AddElement(fO, 30 * perCent);

	// Silicium
	fSi = nist->FindOrBuildMaterial("G4_Si");
	
	// Silicon resin
	fSiResin = new G4Material("SiResin",z=1.,a=1.01*g/mole, 
		     		 density = universe_mean_density, kStateGas,
				 0.1 * kelvin, 1.e-19 * pascal);

	// Optical grease
	fOG = new G4Material("OpticalGrease",z=1.,a=1.01*g/mole, 
		     		 density = universe_mean_density, kStateGas,
				 0.1 * kelvin, 1.e-19 * pascal);

	// Epoxy resin
	fEpResin = new G4Material("EpoxyResin",z=1.,a=1.01*g/mole, 
		     		 density = universe_mean_density, kStateGas,
				 0.1 * kelvin, 1.e-19 * pascal);

	// Vacuuum
	fVacuum = new G4Material("Vacuum",z=1.,a=1.01*g/mole, 
		     		 density = universe_mean_density, kStateGas,
				 0.1 * kelvin, 1.e-19 * pascal);
	
	/// Material properties tables
	//  BC400 optics
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
	
	//  LYSO optics

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
	//fLYSO->GetIonisation()->SetBirksConstant(0.117645*mm/MeV);
	
	// Vacuum and air
	G4double vacuum_Energy[] = {1.5*eV, 4.*eV};
	const G4int vacnum = sizeof(vacuum_Energy) / sizeof(G4double);
	
	G4double vRIND = 1.;
	G4double vacuum_RIND[] = {vRIND, vRIND};
	assert(sizeof(vacuum_RIN) == sizeof(vacuum_Energy));
	
	G4MaterialPropertiesTable* vacuum_mt = new G4MaterialPropertiesTable();
	vacuum_mt->AddProperty("RINDEX", vacuum_Energy, vacuum_RIND, vacnum);
	fVacuum->SetMaterialPropertiesTable(vacuum_mt);
	fAir   ->SetMaterialPropertiesTable(vacuum_mt);
	
	// Silicium
	G4double Si_Energy[] = {.5*eV, 9.*eV};
	const G4int Sinum = sizeof(vacuum_Energy) / sizeof(G4double);

	G4double Si_RIND[] = {3.4, 3.4};
	assert(sizeof(Si_RIND) == sizeof(Si_Energy));
	
	G4MaterialPropertiesTable* Si_mt = new G4MaterialPropertiesTable();
	Si_mt->AddProperty("RINDEX", Si_Energy, Si_RIND, Sinum);
	fSi->SetMaterialPropertiesTable(Si_mt);

	// Silicon resin
	G4double SiRes_RIND[] = {1.41, 1.41};
	assert(sizeof(SiRes_RIND) == sizeof(Si_Energy));
	
	G4MaterialPropertiesTable* SiRes_mt = new G4MaterialPropertiesTable();
	SiRes_mt->AddProperty("RINDEX", Si_Energy, SiRes_RIND, Sinum);
	fSiResin->SetMaterialPropertiesTable(SiRes_mt);

	// Optical grease
	G4double OG_RIND[] = {1.465, 1.465};
	assert(sizeof(OG_RIND) == sizeof(Si_Energy));
	
	G4MaterialPropertiesTable* OG_mt = new G4MaterialPropertiesTable();
	OG_mt->AddProperty("RINDEX", Si_Energy, OG_RIND, Sinum);
	fOG->SetMaterialPropertiesTable(OG_mt);

	// Silicium
	G4double Ep_RIND[] = {1.55, 1.55};
	assert(sizeof(Ep_RIND) == sizeof(Si_Energy));
	
	G4MaterialPropertiesTable* Ep_mt = new G4MaterialPropertiesTable();
	Ep_mt->AddProperty("RINDEX", Si_Energy, Ep_RIND, Sinum);
	fEpResin->SetMaterialPropertiesTable(Ep_mt);
}

void DetectorConstruction::SetCrystalMaterial(G4String name){
	if(name == "BC400") fMaterial = fBC400;
	else if(name == "LYSO") fMaterial = fLYSO;
	G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void DetectorConstruction::SetSiPMmodel(G4String name){
	G4int i = 0, j = 0; //i is for the pitch, j is for the window material;
	if (name == "75PE") {i = 0; j = 1;}
	else if (name == "50PE") {i = 1; j = 1;}
	else if (name == "25PE") {i = 2; j = 1;}
	else if (name == "75CS") {i = 0; j = 0;}
	else if (name == "50CS") {i = 1; j = 0;}
	else if (name == "25CS") {i = 2; j = 0;};
	
	fmodel = name;
	if(j == 1) fMaterialWindow = fEpResin;
	else if(j == 0) fMaterialWindow = fSiResin;

	fNbOfPixelsX = Model::NbPixelsX[i];
	fNbOfPixelsY = Model::NbPixelsY[i];
	
	fSiPM_sizeZ = Model::SiPM_size_Z[j];
	fSiPM_windowZ = Model::window_size_Z[j];

	if(G4RunManager::GetRunManager()->GetUserRunAction()){
		((RunAction*) G4RunManager::GetRunManager()->GetUserRunAction())->SetDNTimeMean(1 / Model::dark_noise_rate[i]);
	}
	G4RunManager::GetRunManager()->ReinitializeGeometry();
}

G4VPhysicalVolume* DetectorConstruction::DefineVolumes(){

    /// MATERIALS AND PARAMETERS
	
    // Crystal parameters


    //SiPM parameters
    G4int nofPixels = fNbOfPixelsX * fNbOfPixelsY; 
    G4double SiPM_sizeXY = fSiPM_sizeXY;
    G4double SiPM_sizeZ = fSiPM_sizeZ;


    // World
    G4double world_sizeX = 2*std::max(fCrysSizeX, fSiPM_sizeXY);
    G4double world_sizeY = 2*std::max(fCrysSizeY, fSiPM_sizeXY);
    G4double world_sizeZ = 2*(fCrysSizeZ + fSiPM_sizeZ + fSiPM_windowZ);

    // World

    fSolidWorld = new G4Box("World", 0.5*world_sizeX, 0.5*world_sizeY, 0.5*world_sizeZ);
    G4LogicalVolume* logicWorld = new G4LogicalVolume(fSolidWorld, fVacuum, "World");
    G4PVPlacement* physWorld = new G4PVPlacement(0, G4ThreeVector(), logicWorld, "World", 0, false, 0, fCheckOverlaps);


    // Scintillator
    fSolidCrys = new G4Box("Crystal", 0.5*fCrysSizeX, 0.5*fCrysSizeY, 0.5*fCrysSizeZ);
    G4LogicalVolume* logicCrys = new G4LogicalVolume(fSolidCrys, fMaterial, "CrysLV");
    fLogicCrys = logicCrys;    

    // SiPM
    G4Box* solidSiPM = new G4Box("SiPM", 0.5*SiPM_sizeXY, 0.5*SiPM_sizeXY, 0.5*(SiPM_sizeZ + fSiPM_windowZ));
    G4LogicalVolume* logicSiPM = new G4LogicalVolume(solidSiPM, fVacuum, "SiPM");

    // SiPM window
    G4Box* solidSiPMwindow = new G4Box("SiPMwindow", 0.5*SiPM_sizeXY, 0.5*SiPM_sizeXY, 0.5*fSiPM_windowZ);
    G4LogicalVolume* logicSiPMwindow = new G4LogicalVolume(solidSiPMwindow, fMaterialWindow, "SiPMwindow");
    
    // SiPM pixel
    G4Box* solidPixel = new G4Box("Pixel", 0.5*SiPM_sizeXY / fNbOfPixelsX, 0.5*SiPM_sizeXY / fNbOfPixelsY, 0.5*SiPM_sizeZ);
    G4LogicalVolume* logicPixel = new G4LogicalVolume(solidPixel, fSi, "PixelLV");
    G4int i = 0, j = 0;
    fLogicPixel = logicPixel;

    for(G4int ipixel = 0; ipixel < nofPixels; ++ipixel){
        i = ipixel%int(fNbOfPixelsX);
        j = ipixel/int(fNbOfPixelsX);
        G4ThreeVector pixel_pos = G4ThreeVector(SiPM_sizeXY /fNbOfPixelsX * (-0.5 * fNbOfPixelsX + i + 0.5), SiPM_sizeXY /fNbOfPixelsY * (-0.5 * fNbOfPixelsY + j + 0.5), -0.5 * (fSiPM_windowZ));
        new G4PVPlacement(0, pixel_pos, logicPixel, "Pixel", logicSiPM, false, ipixel, fCheckOverlaps);
    }

    new G4PVPlacement(0, G4ThreeVector(0, 0, 0.5 * SiPM_sizeZ), logicSiPMwindow, "Window", logicSiPM, false, 0, fCheckOverlaps);

    // Full element
    fSolidElement = new G4Box("Element", 0.5 * (fCrysSizeX + SiPM_sizeXY), 0.5 * (fCrysSizeY + SiPM_sizeXY), 0.5 * (fCrysSizeZ + SiPM_sizeZ + fSiPM_windowZ) + 1*mm);
    
    G4LogicalVolume* logicElement = new G4LogicalVolume(fSolidElement, fVacuum, "Element");

    // Place SiPM in Element
    G4ThreeVector SiPM_pos = G4ThreeVector(0, 0, -0.5*fCrysSizeZ); 
    fSiPMVolume = new G4PVPlacement(0, SiPM_pos, logicSiPM, "SiPM", logicElement, false, 0, fCheckOverlaps);

    // Place Crystal in Element
    if(fAngle == 0 && fAngleWithOpticalGrease == 0){
	    G4ThreeVector crys_pos = G4ThreeVector(fTilt, 0, 0.5 * (SiPM_sizeZ + fSiPM_windowZ));
	    fCrysVolume = new G4PVPlacement(0, crys_pos, logicCrys, "Crystal", logicElement, false, 0, fCheckOverlaps);
    }

    else if(fAngle > 0){
	    G4ThreeVector crys_pos = G4ThreeVector(fTilt, 0, 0.5 * (SiPM_sizeZ + fSiPM_windowZ) + (SiPM_sizeXY*0.5 + fTilt)*tan(fAngle) + fSolidCrys->GetZHalfLength()*(1/cos(fAngle) - 1));
	    G4RotationMatrix* Rotation = new G4RotationMatrix();
	    Rotation->rotateX(0.*deg);
	    Rotation->rotateY(fAngle);
	    Rotation->rotateZ(0.*deg);

	    fCrysVolume = new G4PVPlacement(Rotation, crys_pos, logicCrys, "Crystal", logicElement, false, 0, fCheckOverlaps);
    }

    else if(fAngleWithOpticalGrease > 0){
	    G4ThreeVector crys_pos = G4ThreeVector(fTilt, 0, 0.5 * (SiPM_sizeZ + fSiPM_windowZ) + (SiPM_sizeXY*0.5 + fTilt)*tan(fAngleWithOpticalGrease) + fSolidCrys->GetZHalfLength()*(1/cos(fAngleWithOpticalGrease) - 1));
	    G4Box* solidGreaseC = new G4Box("GreaseC", 0.5*SiPM_sizeXY, 0.5*SiPM_sizeXY, 0.5*(fSiPM_sizeXY*tan(fAngleWithOpticalGrease)));
	    
	    G4RotationMatrix* Rotation = new G4RotationMatrix();
	    Rotation->rotateX(0.*deg);
	    Rotation->rotateY(fAngleWithOpticalGrease);
	    Rotation->rotateZ(0.*deg);

	    G4SubtractionSolid* solidGrease = new G4SubtractionSolid("Grease", solidGreaseC, fSolidCrys, Rotation, G4ThreeVector(fTilt, 0, -solidGreaseC->GetZHalfLength() + fSolidCrys->GetZHalfLength() + (SiPM_sizeXY*0.5 + fTilt)*tan(fAngleWithOpticalGrease) + fSolidCrys->GetZHalfLength()*(1/cos(fAngleWithOpticalGrease) - 1)));
	    
	    G4LogicalVolume* logicGrease = new G4LogicalVolume(solidGrease, fOG, "OGLV");
	    new G4PVPlacement(0, G4ThreeVector(0,0, 0.5 * (SiPM_sizeZ + fSiPM_windowZ) - fSolidCrys->GetZHalfLength() + solidGreaseC->GetZHalfLength()), logicGrease, "Grease", logicElement, false, 0, fCheckOverlaps);
	    
	    fCrysVolume = new G4PVPlacement(Rotation, crys_pos, logicCrys, "Crystal", logicElement, false, 0, fCheckOverlaps);
	    logicGrease->SetVisAttributes(G4Colour(0.,0.1, 0.5, 0.3));
    }

    //Place Element in World
    G4ThreeVector element_pos = G4ThreeVector(0, 0, - (fSolidWorld->GetZHalfLength() - fSolidElement->GetZHalfLength()) + 0.1 * mm);
    fElementVolume = new G4PVPlacement(0, element_pos, logicElement, "Element", logicWorld, false, 0, fCheckOverlaps);
    
    // Scintillator glisur
    if(fGround < 1){
    	G4OpticalSurface* OpScintSurface = new G4OpticalSurface("OpScintSurface");
	OpScintSurface->SetModel(glisur);
	OpScintSurface->SetType(dielectric_dielectric);
	OpScintSurface->SetFinish(ground);
	OpScintSurface->SetPolish(fGround);

	new G4LogicalSkinSurface("ScintSurface", logicCrys, OpScintSurface);
    }


    logicElement->SetVisAttributes(G4VisAttributes::GetInvisible());
    logicSiPM->SetVisAttributes(G4VisAttributes::GetInvisible());
    logicSiPMwindow->SetVisAttributes(G4Colour(0.,0., 0.5, 0.3));
    logicPixel->SetVisAttributes(G4Colour(0.,0.,1., 0.8));
    logicWorld->SetVisAttributes(G4Colour(1, 1, 1, 0.1));
    logicCrys->SetVisAttributes(G4Colour(1, 1, 1, 0.3));
    
    return physWorld;
}

void DetectorConstruction::ConstructSDandField(){
	if(!fScint_SD.Get()){
		G4cout << "Construction /Det/ScintSD" << G4endl;
		ScintSD* scint_SD = new ScintSD("/Det/ScintSD");
		fScint_SD.Put(scint_SD);
	};
	G4SDManager::GetSDMpointer()->AddNewDetector(fScint_SD.Get());
	SetSensitiveDetector(fLogicCrys, fScint_SD.Get());
	
	if(!fPixel_SD.Get()){
		G4cout << "Contruction /Det/PixelSD" << G4endl;
		PixelSD* pixel_SD = new PixelSD("Det/PixelSD", fmodel);
		pixel_SD->SetModel(fmodel);
		pixel_SD->DefineProperties();
		fPixel_SD.Put(pixel_SD);
	}
	else{
		PixelSD* pixel_SD = fPixel_SD.Get();
		pixel_SD->SetModel(fmodel);
		pixel_SD->DefineProperties();
	}
	G4SDManager::GetSDMpointer()->AddNewDetector(fPixel_SD.Get());
	SetSensitiveDetector(fLogicPixel, fPixel_SD.Get());
}

void DetectorConstruction :: SetCrystalSize(G4double size){
	fCrysSizeX = size;
	fCrysSizeY = size;
	fCrysSizeZ = size;
	fSolidCrys->SetXHalfLength(size*0.5);
	fSolidCrys->SetYHalfLength(size*0.5);
	fSolidCrys->SetZHalfLength(size*0.5);
	
	fSolidWorld->SetXHalfLength(std::max(fCrysSizeX, fSiPM_sizeXY));
	fSolidWorld->SetYHalfLength(std::max(fCrysSizeY, fSiPM_sizeXY));
	fSolidWorld->SetZHalfLength((fCrysSizeZ + fSiPM_sizeZ + fSiPM_windowZ));
		
	G4double SiPM_sizeZ = fSiPM_sizeZ + fSiPM_windowZ;

	fSolidElement->SetZHalfLength(0.5 * (size + SiPM_sizeZ));
	fSolidElement->SetYHalfLength(0.5 * (size + fSiPM_sizeXY));
	fSolidElement->SetXHalfLength(0.5 * (size + fSiPM_sizeXY));

	fCrysVolume->SetTranslation(G4ThreeVector(0, 0, 0.5 * SiPM_sizeZ));
	fSiPMVolume->SetTranslation(G4ThreeVector(0, 0, -0.5 * (size)));
	fElementVolume->SetTranslation(G4ThreeVector(0, 0, -0.5 * (fSolidWorld->GetZHalfLength() - fSolidElement->GetZHalfLength())));
	G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

void DetectorConstruction :: SetCrystalSize(G4ThreeVector size){
	fCrysSizeX = size.getX();
	fCrysSizeY = size.getY();
	fCrysSizeZ = size.getZ();
	fSolidCrys->SetXHalfLength(size.getX()*0.5);
	fSolidCrys->SetYHalfLength(size.getY()*0.5);
	fSolidCrys->SetZHalfLength(size.getZ()*0.5);
	
	fSolidWorld->SetXHalfLength(std::max(fCrysSizeX, fSiPM_sizeXY));
	fSolidWorld->SetYHalfLength(std::max(fCrysSizeY, fSiPM_sizeXY));
	fSolidWorld->SetZHalfLength((fCrysSizeZ + fSiPM_sizeZ + fSiPM_windowZ));
		
	G4double SiPM_sizeZ = fSiPM_sizeZ + fSiPM_windowZ;

	fSolidElement->SetZHalfLength(0.5 * (size.getZ() + SiPM_sizeZ));
	fSolidElement->SetYHalfLength(0.5 * (size.getY() + fSiPM_sizeXY));
	fSolidElement->SetXHalfLength(0.5 * (size.getX() + fSiPM_sizeXY));

	fCrysVolume->SetTranslation(G4ThreeVector(0, 0, 0.5 * SiPM_sizeZ));
	fSiPMVolume->SetTranslation(G4ThreeVector(0, 0, -0.5 * (size.getZ())));
	fElementVolume->SetTranslation(G4ThreeVector(0, 0, -(fSolidWorld->GetZHalfLength() - fSolidElement->GetZHalfLength())));
	G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

void DetectorConstruction  :: SetGround(G4double val){
	fGround = val;
	G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void DetectorConstruction :: SetAngle(G4double val){
	fAngle = val;
	fAngleWithOpticalGrease = 0;
	G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void DetectorConstruction :: SetAngleWithOpticalGrease(G4double val){
	fAngle = 0;
	fAngleWithOpticalGrease = val;
	G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void DetectorConstruction :: SetTilt(G4double val){
	fTilt = val;
	G4RunManager::GetRunManager()->ReinitializeGeometry();
}


