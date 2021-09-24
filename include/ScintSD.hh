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
		virtual void EndOfEvent(G4HCofThisEvent*);
		virtual void clear();
		virtual void DrawAll();
		virtual void PrintAll();
		
	private:
		ScintHitsCollection* fScintCollection;
		G4double fEin, fEdep, fEout, fDelta, fThetaIn, fTrackLength, fThetaPositron;
		G4int fBounce;
		G4ThreeVector fDirIN, fDirOUT;
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
		
		G4int fPhotonsCmd, fTracksCmd;
};

#endif


