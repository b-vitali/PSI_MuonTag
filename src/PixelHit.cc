/// \file  PixelHit.cc
/// \brief Implementation of the PixelHit class

#include "PixelHit.hh"
#include "G4ios.hh"
#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

G4ThreadLocal G4Allocator<PixelHit>* PixelHitAllocator = nullptr;

PixelHit::PixelHit() : 
	fNCells(0), fNPhotoElectrons(0), fPhysVol(nullptr), 
	fPhysVolMother(nullptr), fDrawit(false), fPixelNumber(-1){}

PixelHit::~PixelHit(){}

PixelHit::PixelHit(const PixelHit &right) : G4VHit(){
	fNCells = right.fNCells; // number of active cells
	fNPhotoElectrons = right.fNPhotoElectrons; // number of photoelectrons generated
	fCells = right.fCells; // vector of active cells
	fCellTime = right.fCellTime; // vector of arrive time in cell
	fOCTflag = right.fOCTflag;
	fDNflag = right.fDNflag;
	fPhysVolMother = right.fPhysVolMother;
	fPhysVolGMother = right.fPhysVolGMother;
	fDrawit = right.fDrawit;
	fPixelNumber = right.fPixelNumber;
}

const PixelHit& PixelHit::operator=(const PixelHit &right){
	fNCells = right.fNCells; // number of active cells
	fNPhotoElectrons = right.fNPhotoElectrons; // number of photoelectrons generated
	fCells = right.fCells; // vector of active cells
	fCellTime = right.fCellTime; // vector of arrive time in cell
	fOCTflag = right.fOCTflag;
	fDNflag = right.fDNflag;
	fPhysVol = right.fPhysVol;
	fPhysVolMother = right.fPhysVolMother;
	fPhysVolGMother = right.fPhysVolGMother;
	fDrawit = right.fDrawit;
	fPixelNumber = right.fPixelNumber;
	return* this;
}


G4bool PixelHit::operator==(const PixelHit &) const{
	return false;
}

void PixelHit::Draw(){
	if(fDrawit && fPhysVol){
		G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
		if(pVVisManager){
			G4VisAttributes attribs(G4Colour(1.,0.,0., 0.8));
			attribs.SetForceSolid(true);
			fPhysVol->GetLogicalVolume()->SetVisAttributes(G4VisAttributes::GetInvisible());

			G4RotationMatrix rot;
			G4ThreeVector pos = fPhysVol->GetTranslation();
			if(fPhysVolMother != nullptr) pos += fPhysVolMother->GetTranslation();
			if(fPhysVolGMother != nullptr) pos += fPhysVolGMother->GetTranslation();
			G4Transform3D trans(rot, pos);
			pVVisManager->Draw(*fPhysVol, attribs, trans);

		}
	}
}

void PixelHit::Print(){}
