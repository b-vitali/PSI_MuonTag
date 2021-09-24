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

		inline void SetEin (G4double val){fEin  = val;}
		inline void SetEdep(G4double val){fEdep = val;}
		inline void SetEout(G4double val){fEout = val;}
		inline void SetEdelta(G4double val){fDelta = val;}
		inline void SetThetaIn(G4double val){fThetaIn = val;}
		inline void SetTrackLength(G4double val){fTrackLength = val;}
		inline void SetThetaPositron(G4double val){fThetaPositron = val;}
		inline void SetBounce(G4int val){fBounce = val;}
		inline void SetDecayTime(G4double val){fDecayTime = val;}

		inline G4double GetEin (){return fEin;}
		inline G4double GetEdep(){return fEdep;}
		inline G4double GetEout(){return fEout;}
		inline G4double GetEdelta(){return fDelta;}
		inline G4double GetThetaIn(){return fThetaIn;}
		inline G4double GetTrackLength(){return fTrackLength;}
		inline G4double GetThetaPositron(){return fThetaPositron;}
		inline G4int GetBounce(){return fBounce;}
		inline G4double GetDecayTime(){return fDecayTime;}

		inline void SetPosX(std::vector<G4double> pos){fPosX = pos;}
		inline std::vector<G4double> GetPosX(){return fPosX;}

		inline void SetPosY(std::vector<G4double> pos){fPosY = pos;}
		inline std::vector<G4double> GetPosY(){return fPosY;}

		inline void SetPosZ(std::vector<G4double> pos){fPosZ = pos;}
		inline std::vector<G4double> GetPosZ(){return fPosZ;}

		inline void SetTime(std::vector<G4double> time){fTime = time;}
		inline std::vector<G4double> GetTime(){return fTime;}

		inline void SetNgamma(G4int ngamma){fNgamma = ngamma;}
		inline G4int GetNgamma(){return fNgamma;}

		inline void SetNgammaSec(G4int ngammasec){fNgammaSec = ngammasec;}
		inline G4int GetNgammaSec(){return fNgammaSec;}

		inline void SetCer(std::vector<G4int> cer){fCer = cer;}
		inline std::vector<G4int> GetCer(){return fCer;}

		inline void SetThetaGamma(std::vector<G4double> thetagamma){fThetaGamma = thetagamma;}
		inline std::vector<G4double> GetThetaGamma(){return fThetaGamma;}

		inline void SetTimeGamma(std::vector<G4double> timegamma){fTimeGamma = timegamma;}
		inline std::vector<G4double> GetTimeGamma(){return fTimeGamma;}

		inline void SetEGamma(std::vector<G4double> egamma){fEGamma = egamma;}
		inline std::vector<G4double> GetEGamma(){return fEGamma;}

		inline void SetNCer(G4int ncer){fNCer = ncer;}
		inline G4int GetNCer(){return fNCer;}

		inline void SetCurrentRight(G4int right){fRight = right;}
		inline G4int GetCurrentRight(){return fRight;}

		inline void SetCurrentLeft(G4int left){fLeft = left;}
		inline G4int GetCurrentLeft(){return fLeft;}
		
		inline void SetCurrentDown(G4int down){fDown = down;}
		inline G4int GetCurrentDown(){return fDown;}

		inline void SetCurrentUp(G4int up){fUp = up;}
		inline G4int GetCurrentUp(){return fUp;}

		inline void SetCurrentBack(G4int back){fBack = back;}
		inline G4int GetCurrentBack(){return fBack;}

		inline void SetCurrentFront(G4int front){fFront = front;}
		inline G4int GetCurrentFront(){return fFront;}

		inline void SetSiPM(G4int sipm){fSiPM = sipm;}
		inline G4int GetSiPM(){return fSiPM;}

		inline void Clear(){fEin = 0; fEdep = 0; fEout = 0; fDelta = 0; fThetaIn = 0; fTrackLength = 0; fThetaPositron = 0; fBounce = 0; fPosX.clear(); fPosY.clear(); fPosZ.clear(); fTime.clear(); fNgamma = 0; fNgammaSec = 0; fNCer = 0; fCer.clear(); fThetaGamma.clear(); fTimeGamma.clear(); fEGamma.clear(); fRight = 0; fLeft = 0; fDown = 0; fUp = 0; fBack = 0; fFront = 0; fSiPM = 0; fDecayTime = -1;}
		inline const G4VPhysicalVolume* GetPhysV(){return fPhysVol;}

	private:
		G4double fEin, fEdep, fEout, fDelta, fThetaIn, fTrackLength, fThetaPositron;
		G4int fBounce;
		std::vector<G4double> fPosX;
		std::vector<G4double> fPosY;
		std::vector<G4double> fPosZ;
		std::vector<G4double> fTime;
		G4int fNgamma;
		G4int fNgammaSec;
		std::vector<G4int> fCer;
		std::vector<G4double> fThetaGamma;
		std::vector<G4double> fTimeGamma;
		std::vector<G4double> fEGamma;
		G4int fNCer;
		G4int fRight;
		G4int fLeft;
		G4int fDown;
		G4int fUp;
		G4int fBack;
		G4int fFront;
		G4int fSiPM;
		G4double fDecayTime;
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


