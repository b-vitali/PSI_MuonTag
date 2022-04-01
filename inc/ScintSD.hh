/// \file  ScintSD.hh
/// \brief Definition of the ScintSD class

#ifndef ScintSD_h
#define ScintSD_h 1

#include "ScintHit.hh"

#include "G4ThreeVector.hh"
#include "G4VSensitiveDetector.hh"
#include "G4OpticalPhoton.hh"

#include "TVector3.h"

#include "g4root.hh"			// to access root stuff

#include "G4SDManager.hh"

class G4Step;
class G4HCofThisEvent;
class G4VLogicalVolume;

class ScintSD : public G4VSensitiveDetector{
	public:
		ScintSD(G4String name);
		ScintSD(G4String name, G4int ntuple);
		virtual ~ScintSD();

		virtual void Initialize(G4HCofThisEvent* hitsCE);
		virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory*ROhist);
		
		virtual void FillHit();
		virtual void FillNtupla(G4AnalysisManager *man, ScintHit* scintHit, G4int ntupla);

		virtual void EndOfEvent(G4HCofThisEvent* hitsCE);
		virtual void clear();
		virtual void DrawAll();
		virtual void PrintAll();
		
	private:
		ScintHitsCollection* fScintCollection;
		G4int fEvent, fScintNo, fParticleID;
		G4int fEvent_ID, fScintNo_ID, fParticleID_ID;

		G4double fEin, fEdep, fEout, fDelta, fThetaIn, fTrackLength, fThetaOut;
		G4double fEin_ID, fEdep_ID, fEout_ID, fDelta_ID, fThetaIn_ID, fTrackLength_ID, fThetaOut_ID;

		G4int fRight, fLeft, fDown, fUp, fBack, fFront;
		G4int fRight_ID, fLeft_ID, fDown_ID, fUp_ID, fBack_ID, fFront_ID;

		G4int fBounce;

		G4ThreeVector fDirIn, fDirOut;
		G4ThreeVector fDirIn_trans, fDirOut_trans;
		
		G4ThreeVector fMomIn, fMomOut;
		G4int fMomInX_ID, fMomInY_ID, fMomInZ_ID, fMomOutX_ID, fMomOutY_ID, fMomOutZ_ID;

		G4ThreeVector fPosIn, fPosOut;
		G4int fPosInX_ID, fPosInY_ID, fPosInZ_ID, fPosOutX_ID, fPosOutY_ID, fPosOutZ_ID;
		G4double fTimeIn, fTimeOut;
		G4int fTimeIn_ID, fTimeOut_ID;

		G4int fNgamma, fNgammaSec, fNCer;
		G4int fNgamma_ID, fNgammaSec_ID, fNCer_ID;

		G4int gammaTime_ID, gammaFront_ID, gammaSide_ID, gammaTop_ID;

		G4double fDecayTime;
		G4int fDecayTime_ID;

		std::vector<G4int> fCer;
		std::vector<G4double> fThetaGamma;
		std::vector<G4double> fTimeGamma;
		std::vector<G4double> fEGamma;
		G4int fSiPM;
		
		G4int fPhotonsCmd, fTracksCmd;
		G4int hitsCID;

		std::vector<G4double> fTest;


};

#endif


