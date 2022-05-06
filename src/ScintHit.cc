/// \file  ScintHit.cc
/// \brief Implementation of the ScintHit class

#include "ScintHit.hh"
#include "G4ios.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"

G4ThreadLocal G4Allocator<ScintHit>* ScintHitAllocator = nullptr;

ScintHit::ScintHit() : 
	fPhysVol(nullptr){}

ScintHit::ScintHit(G4VPhysicalVolume* pVol) : fPhysVol(pVol){}

ScintHit::~ScintHit(){}

ScintHit::ScintHit(const ScintHit &hit) : G4VHit(){	
	fEvent = hit.fEvent;
	fScintNo = hit.fScintNo;
	fParticleID = hit.fParticleID;

	fEin  = hit.fEin;
	fEdep = hit.fEdep;
	fEout = hit.fEout;
	fDelta = hit.fDelta;
	fThetaIn = hit.fThetaIn;
	fTrackLength = hit.fTrackLength;
	fThetaOut = hit.fThetaOut;
	// fBounce = hit.fBounce;
	fPosInX = hit.fPosInX;
	fPosInY = hit.fPosInY;
	fPosInZ = hit.fPosInZ;
	fMomInX = hit.fMomInX;
	fMomInY = hit.fMomInY;
	fMomInZ = hit.fMomInZ;
	fTimeIn = hit.fTimeIn;
	fPosOutX = hit.fPosOutX;
	fPosOutY = hit.fPosOutY;
	fPosOutZ = hit.fPosOutZ;
	fMomOutX = hit.fMomOutX;
	fMomOutY = hit.fMomOutY;
	fMomOutZ = hit.fMomOutZ;
	fTimeOut = hit.fTimeOut;

	fNgamma = hit.fNgamma;
	fNgammaSec = hit.fNgammaSec;
	fNCer = hit.fNCer;
	fReflection = hit.fReflection;
	fAbsorption = hit.fAbsorption;

	fRight = hit.fRight;
	fLeft = hit.fLeft;
	fDown = hit.fDown;
	fUp = hit.fUp;
	fBack = hit.fBack;
	fFront = hit.fFront;
	fDecayTime = hit.fDecayTime;
	fPhysVol = hit.fPhysVol;
}

const ScintHit& ScintHit::operator=(const ScintHit &hit){
	fEvent = hit.fEvent;
	fScintNo = hit.fScintNo;
	fParticleID = hit.fParticleID;
	
	fEin  = hit.fEin;
	fEdep = hit.fEdep;
	fEout = hit.fEout;
	fDelta = hit.fDelta;
	fThetaIn = hit.fThetaIn;
	fTrackLength = hit.fTrackLength;
	fThetaOut = hit.fThetaOut;
	// fBounce = hit.fBounce;
	fPosInX = hit.fPosInX;
	fPosInY = hit.fPosInY;
	fPosInZ = hit.fPosInZ;
	fMomInX = hit.fMomInX;
	fMomInY = hit.fMomInY;
	fMomInZ = hit.fMomInZ;
	fTimeIn = hit.fTimeIn;
	fPosOutX = hit.fPosOutX;
	fPosOutY = hit.fPosOutY;
	fPosOutZ = hit.fPosOutZ;
	fMomOutX = hit.fMomOutX;
	fMomOutY = hit.fMomOutY;
	fMomOutZ = hit.fMomOutZ;
	fTimeOut = hit.fTimeOut;

	fNgamma = hit.fNgamma;
	fNgammaSec = hit.fNgammaSec;
	fNCer = hit.fNCer;
	fReflection = hit.fReflection;
	fAbsorption = hit.fAbsorption;
	
	fRight = hit.fRight;
	fLeft = hit.fLeft;
	fDown = hit.fDown;
	fUp = hit.fUp;
	fBack = hit.fBack;
	fFront = hit.fFront;
	fDecayTime = hit.fDecayTime;
	fPhysVol = hit.fPhysVol;
	return* this;
}

G4bool ScintHit::operator==(const ScintHit &) const{
	return false;
}

void ScintHit::Draw(){}

void ScintHit::Print(){}


