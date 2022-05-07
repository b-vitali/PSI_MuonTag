/// \file  ScintHit.hh
/// \brief Definition of the ScintHit class

#ifndef ScintHit_h
#define ScintHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"

#include "TVector3.h"

#include "tls.hh"

#include <vector>

class ScintHit : public G4VHit{
	public:
		ScintHit();
		ScintHit(G4VPhysicalVolume* pVol);
		virtual ~ScintHit();

		ScintHit(const ScintHit &right);
		const ScintHit& operator=(const ScintHit &right);
		G4bool operator==(const ScintHit &right) const;

		inline void *operator new(size_t);
		inline void operator delete(void *aHit);

		virtual void Draw();
		virtual void Print();

		inline void SetEvent(std::vector<G4int> val){fEvent = val;}
		inline void SetScintNo(std::vector<G4int> val){fScintNo = val;}
		inline void SetParticleID(std::vector<G4int> val){fParticleID = val;}

		inline void SetEin (std::vector<G4double> val){fEin  = val;}
		inline void SetEdep(std::vector<G4double> val){fEdep = val;}
		inline void SetEout(std::vector<G4double> val){fEout = val;}
		inline void SetEdelta(std::vector<G4double> val){fDelta = val;}
		inline void SetThetaIn(std::vector<G4double> val){fThetaIn = val;}
		inline void SetTrackLength(std::vector<G4double> val){fTrackLength = val;}
		inline void SetThetaOut(std::vector<G4double> val){fThetaOut = val;}
		// inline void SetBounce(G4int val){fBounce = val;}
		inline void SetDecayTime(std::vector<G4double> val){fDecayTime = val;}

		inline std::vector<G4int> GetEvent(){return fEvent;}
		inline std::vector<G4int> GetScintNo(){return fScintNo;}
		inline std::vector<G4int> GetParticleID(){return fParticleID;}

		inline std::vector<G4double> GetEin (){return fEin;}
		inline std::vector<G4double> GetEdep(){return fEdep;}
		inline std::vector<G4double> GetEout(){return fEout;}
		inline std::vector<G4double> GetEdelta(){return fDelta;}
		inline std::vector<G4double> GetThetaIn(){return fThetaIn;}
		inline std::vector<G4double> GetTrackLength(){return fTrackLength;}
		inline std::vector<G4double> GetThetaOut(){return fThetaOut;}
		// inline G4int GetBounce(){return fBounce;}
		inline std::vector<G4double> GetDecayTime(){return fDecayTime;}

		inline void SetPosInX(std::vector<G4double> posX){fPosInX = posX;}
		inline void SetPosInY(std::vector<G4double> posY){fPosInY = posY;}
		inline void SetPosInZ(std::vector<G4double> posZ){fPosInZ = posZ;}
		inline void SetMomInX(std::vector<G4double> momX){fMomInX = momX;}
		inline void SetMomInY(std::vector<G4double> momY){fMomInY = momY;}
		inline void SetMomInZ(std::vector<G4double> momZ){fMomInZ = momZ;}
		inline void SetTimeIn(std::vector<G4double> time){fTimeIn = time;}
		inline std::vector<G4double> GetPosInX(){return fPosInX;}
		inline std::vector<G4double> GetPosInY(){return fPosInY;}
		inline std::vector<G4double> GetPosInZ(){return fPosInZ;}
		inline std::vector<G4double> GetMomInX(){return fMomInX;}
		inline std::vector<G4double> GetMomInY(){return fMomInY;}
		inline std::vector<G4double> GetMomInZ(){return fMomInZ;}
		inline std::vector<G4double> GetTimeIn(){return fTimeIn;}

		inline void SetPosOutX(std::vector<G4double> posX){fPosOutX = posX;}
		inline void SetPosOutY(std::vector<G4double> posY){fPosOutY = posY;}
		inline void SetPosOutZ(std::vector<G4double> posZ){fPosOutZ = posZ;}
		inline void SetMomOutX(std::vector<G4double> momX){fMomOutX = momX;}
		inline void SetMomOutY(std::vector<G4double> momY){fMomOutY = momY;}
		inline void SetMomOutZ(std::vector<G4double> momZ){fMomOutZ = momZ;}
		inline void SetTimeOut(std::vector<G4double> time){fTimeOut = time;}
		inline std::vector<G4double> GetPosOutX(){return fPosOutX;}
		inline std::vector<G4double> GetPosOutY(){return fPosOutY;}
		inline std::vector<G4double> GetPosOutZ(){return fPosOutZ;}
		inline std::vector<G4double> GetMomOutX(){return fMomOutX;}
		inline std::vector<G4double> GetMomOutY(){return fMomOutY;}
		inline std::vector<G4double> GetMomOutZ(){return fMomOutZ;}
		inline std::vector<G4double> GetTimeOut(){return fTimeOut;}

		inline void SetNgamma(std::vector<G4int> ngamma){fNgamma = ngamma;}
		inline std::vector<G4int> GetNgamma(){return fNgamma;}

		inline void SetNgammaSec(std::vector<G4int> ngammasec){fNgammaSec = ngammasec;}
		inline std::vector<G4int> GetNgammaSec(){return fNgammaSec;}

		inline void SetNCer(std::vector<G4int> ncer){fNCer = ncer;}
		inline std::vector<G4int> GetNCer(){return fNCer;}		
		
		inline void SetNAbsorption(std::vector<G4int> fabsorption){fAbsorption = fabsorption;}
		inline std::vector<G4int> GetNAbsorption(){return fAbsorption;}
			
		inline void SetNReflection(std::vector<G4int> freflection){fReflection = freflection;}
		inline std::vector<G4int> GetNReflection(){return fReflection;}

		inline void SetCurrentRight(std::vector<G4int> right){fRight = right;}
		inline std::vector<G4int> GetCurrentRight(){return fRight;}

		inline void SetCurrentLeft(std::vector<G4int> left){fLeft = left;}
		inline std::vector<G4int> GetCurrentLeft(){return fLeft;}
		
		inline void SetCurrentDown(std::vector<G4int>  down){fDown = down;}
		inline std::vector<G4int> GetCurrentDown(){return fDown;}

		inline void SetCurrentUp(std::vector<G4int> up){fUp = up;}
		inline std::vector<G4int> GetCurrentUp(){return fUp;}

		inline void SetCurrentBack(std::vector<G4int> back){fBack = back;}
		inline std::vector<G4int> GetCurrentBack(){return fBack;}

		inline void SetCurrentFront(std::vector<G4int> front){fFront = front;}
		inline std::vector<G4int> GetCurrentFront(){return fFront;}

		inline void Clear(){}
		inline const G4VPhysicalVolume* GetPhysV(){return fPhysVol;}

	private:
		// G4int fBounce; 
		std::vector<G4int> fEvent, fScintNo, fParticleID, fNgamma, fNgammaSec, fNCer;
 		std::vector<G4int> fRight, fLeft, fDown, fUp, fBack, fFront;
		std::vector<G4double> fEin, fEdep, fEout, fDelta, fThetaIn, fTrackLength, fThetaOut, fDecayTime;

		std::vector<G4int> fAbsorption, fReflection;

		std::vector<G4double> fPosInX, fPosInY, fPosInZ; 
		std::vector<G4double> fPosOutX, fPosOutY, fPosOutZ; 
		std::vector<G4double> fMomInX, fMomInY, fMomInZ; 
		std::vector<G4double> fMomOutX, fMomOutY, fMomOutZ; 
		std::vector<G4double> fTimeIn, fTimeOut;

		const G4VPhysicalVolume* fPhysVol;
};

typedef G4THitsCollection<ScintHit> ScintHitsCollection;

extern G4ThreadLocal G4Allocator<ScintHit>* ScintHitAllocator;

inline void* ScintHit::operator new(size_t){
	if(!ScintHitAllocator) ScintHitAllocator = new G4Allocator<ScintHit>;
	return (void*) ScintHitAllocator->MallocSingle();
}

inline void ScintHit::operator delete(void* aHit){
	ScintHitAllocator->FreeSingle((ScintHit*) aHit);
}

#endif


