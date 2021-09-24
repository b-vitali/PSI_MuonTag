/// \file  ScintHit.cc
/// \brief Implementation of the ScintHit class

#include "ScintHit.hh"
#include "G4ios.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"

G4ThreadLocal G4Allocator<ScintHit>* ScintHitAllocator = nullptr;

ScintHit::ScintHit() : 
	fEin(0.), fEdep(0.), fEout(0.), fDelta(0), fThetaPositron(0), fBounce(0), 
	fNgamma(0), fNgammaSec(0), fNCer(0), fRight(0), fLeft(0), 
	fDown(0), fUp(0), fBack(0), fFront(0), fDecayTime(-1), fPhysVol(nullptr){}

ScintHit::ScintHit(G4VPhysicalVolume* pVol) : fPhysVol(pVol){}

ScintHit::~ScintHit(){}

ScintHit::ScintHit(const ScintHit &right) : G4VHit(){
	fEin  = right.fEin;
	fEdep = right.fEdep;
	fEout = right.fEout;
	fDelta = right.fDelta;
	fThetaIn = right.fThetaIn;
	fTrackLength = right.fTrackLength;
	fThetaPositron = right.fThetaPositron;
	fBounce = right.fBounce;
	fPosX = right.fPosX;
	fPosY = right.fPosY;
	fPosZ = right.fPosZ;
	fTime = right.fTime;
	fNgamma = right.fNgamma;
	fNgammaSec = right.fNgammaSec;
	fCer = right.fCer;
	fThetaGamma = right.fThetaGamma;
	fTimeGamma = right.fTimeGamma;
	fEGamma = right.fEGamma;
	fNCer = right.fNCer;
	fRight = right.fRight;
	fLeft = right.fLeft;
	fDown = right.fDown;
	fUp = right.fUp;
	fBack = right.fBack;
	fFront = right.fFront;
	fSiPM = right.fSiPM;
	fDecayTime = right.fDecayTime;
	fPhysVol = right.fPhysVol;
}

const ScintHit& ScintHit::operator=(const ScintHit &right){
	fEin  = right.fEin;
	fEdep = right.fEdep;
	fEout = right.fEout;
	fDelta = right.fDelta;
	fThetaIn = right.fThetaIn;
	fTrackLength = right.fTrackLength;
	fThetaPositron = right.fThetaPositron;
	fBounce = right.fBounce;
	fPosX = right.fPosX;
	fPosY = right.fPosY;
	fPosZ = right.fPosZ;
	fTime = right.fTime;
	fNgamma = right.fNgamma;
	fNgammaSec = right.fNgammaSec;
	fCer = right.fCer;
	fThetaGamma = right.fThetaGamma;
	fTimeGamma = right.fTimeGamma;
	fEGamma = right.fEGamma;
	fNCer = right.fNCer;
	fRight = right.fRight;
	fLeft = right.fLeft;
	fDown = right.fDown;
	fUp = right.fUp;
	fBack = right.fBack;
	fFront = right.fFront;
	fSiPM = right.fSiPM;
	fDecayTime = right.fDecayTime;
	fPhysVol = right.fPhysVol;
	return* this;
}

G4bool ScintHit::operator==(const ScintHit &) const{
	return false;
}

void ScintHit::Draw(){}

void ScintHit::Print(){}


