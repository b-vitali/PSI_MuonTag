/// \file  ScintSD.hh
/// \brief Definition of the ScintSD class

#ifndef ScintSD_h
#define ScintSD_h 1

#include "ScintHit.hh"

#include "G4ThreeVector.hh"
#include "G4VSensitiveDetector.hh"
#include "G4OpticalPhoton.hh"

#include "TVector3.h"

class G4Step;
class G4HCofThisEvent;
class G4VLogicalVolume;

class ScintSD : public G4VSensitiveDetector{
	public:
		ScintSD(G4String name);
		virtual ~ScintSD();

		virtual void Initialize(G4HCofThisEvent*);
		virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory*);
		
		virtual void FillHit();

		virtual void EndOfEvent(G4HCofThisEvent*);
		virtual void clear();
		virtual void DrawAll();
		virtual void PrintAll();
		
	private:
		ScintHitsCollection* fScintCollection;
		G4int fEvent, fScintNo, fParticleID;
		G4double fEin, fEdep, fEout, fDelta, fThetaIn, fTrackLength, fThetaOut;
		G4int fBounce;
		G4ThreeVector fDirIN, fDirOUT;
		
		G4ThreeVector fPosIn;
		G4double fTimeIn;
		G4ThreeVector fPosOut;
		G4double fTimeOut;

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
		
		G4int fPhotonsCmd, fTracksCmd;
};

#endif


