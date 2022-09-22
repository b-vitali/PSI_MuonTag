/// \file  SiPMHit.cc
/// \brief Implementation of the SiPMHit class

#include "SiPMHit.hh"
#include "G4ios.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"

G4ThreadLocal G4Allocator<SiPMHit>* SiPMHitAllocator = nullptr;

SiPMHit::SiPMHit() : 
	fPhysVol(nullptr){}

SiPMHit::SiPMHit(G4VPhysicalVolume* pVol) : fPhysVol(pVol){}

SiPMHit::~SiPMHit(){}

SiPMHit::SiPMHit(const SiPMHit &hit) : G4VHit(){	
	fEvent = hit.fEvent;
	fSiPMNo = hit.fSiPMNo;
	fParticleID = hit.fParticleID;

	fEin  = hit.fEin;
	fThetaIn = hit.fThetaIn;
	fPosInX = hit.fPosInX;
	fPosInY = hit.fPosInY;
	fPosInZ = hit.fPosInZ;
	fMomInX = hit.fMomInX;
	fMomInY = hit.fMomInY;
	fMomInZ = hit.fMomInZ;
	fTimeIn = hit.fTimeIn;

	fNgamma = hit.fNgamma;
	fNgammaSec = hit.fNgammaSec;

	fRight = hit.fRight;
	fLeft = hit.fLeft;
	fDown = hit.fDown;
	fUp = hit.fUp;
	fBack = hit.fBack;
	fFront = hit.fFront;
	fPhysVol = hit.fPhysVol;
}

const SiPMHit& SiPMHit::operator=(const SiPMHit &hit){
	fEvent = hit.fEvent;
	fSiPMNo = hit.fSiPMNo;
	fParticleID = hit.fParticleID;
	
	fEin  = hit.fEin;
	fThetaIn = hit.fThetaIn;
	fPosInX = hit.fPosInX;
	fPosInY = hit.fPosInY;
	fPosInZ = hit.fPosInZ;
	fMomInX = hit.fMomInX;
	fMomInY = hit.fMomInY;
	fMomInZ = hit.fMomInZ;
	fTimeIn = hit.fTimeIn;

	fNgamma = hit.fNgamma;
	fNgammaSec = hit.fNgammaSec;
	
	fRight = hit.fRight;
	fLeft = hit.fLeft;
	fDown = hit.fDown;
	fUp = hit.fUp;
	fBack = hit.fBack;
	fFront = hit.fFront;
	fPhysVol = hit.fPhysVol;
	return* this;
}

G4bool SiPMHit::operator==(const SiPMHit &) const{
	return false;
}

void SiPMHit::Draw(){}

void SiPMHit::Print(){}

