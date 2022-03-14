/// \file  ScintHit.cc
/// \brief Implementation of the ScintHit class

#include "ScintHit.hh"
#include "G4ios.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"

G4ThreadLocal G4Allocator<ScintHit>* ScintHitAllocator = nullptr;

ScintHit::ScintHit() : 
	fEvent(-1), fScintNo(-1), fParticleID(0), fEin(0.), fEdep(0.), fEout(0.), fDelta(0), fThetaOut(0), fBounce(0), 
	fNgamma(0), fNgammaSec(0), fNCer(0), fRight(0), fLeft(0), 
	fDown(0), fUp(0), fBack(0), fFront(0), fDecayTime(-1), fPhysVol(nullptr){}

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
	fBounce = hit.fBounce;
	fPosIn = hit.fPosIn;
	fTimeIn = hit.fTimeIn;
	fPosOut = hit.fPosOut;
	fTimeOut = hit.fTimeOut;

	fNgamma = hit.fNgamma;
	fNgammaSec = hit.fNgammaSec;
	fCer = hit.fCer;
	fThetaGamma = hit.fThetaGamma;
	fTimeGamma = hit.fTimeGamma;
	fEGamma = hit.fEGamma;
	fNCer = hit.fNCer;
	fRight = hit.fRight;
	fLeft = hit.fLeft;
	fDown = hit.fDown;
	fUp = hit.fUp;
	fBack = hit.fBack;
	fFront = hit.fFront;
	fSiPM = hit.fSiPM;
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
	fBounce = hit.fBounce;
	fPosIn = hit.fPosIn;
	fTimeIn = hit.fTimeIn;
	fPosOut = hit.fPosOut;
	fTimeOut = hit.fTimeOut;
	fNgamma = hit.fNgamma;
	fNgammaSec = hit.fNgammaSec;
	fCer = hit.fCer;
	fThetaGamma = hit.fThetaGamma;
	fTimeGamma = hit.fTimeGamma;
	fEGamma = hit.fEGamma;
	fNCer = hit.fNCer;
	fRight = hit.fRight;
	fLeft = hit.fLeft;
	fDown = hit.fDown;
	fUp = hit.fUp;
	fBack = hit.fBack;
	fFront = hit.fFront;
	fSiPM = hit.fSiPM;
	fDecayTime = hit.fDecayTime;
	fPhysVol = hit.fPhysVol;
	return* this;
}

G4bool ScintHit::operator==(const ScintHit &) const{
	return false;
}

void ScintHit::Draw(){}

void ScintHit::Print(){}


