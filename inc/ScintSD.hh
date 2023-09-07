/// \file  ScintSD.hh
/// \brief Definition of the ScintSD class

#ifndef ScintSD_h
#define ScintSD_h 1

#include "ScintHit.hh"

#include "G4ThreeVector.hh"
#include "G4VSensitiveDetector.hh"
#include "G4OpticalPhoton.hh"

#include "TVector3.h"

#include "G4AnalysisManager.hh"			// to access root stuff

#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TVectorD.h"

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

		virtual void CreateEntry(G4Step *aStep);

		virtual void FillHit();
		virtual void FillNtupla(G4AnalysisManager *man, ScintHit* scintHit, G4int ntupla);

		virtual void EndOfEvent(G4HCofThisEvent* hitsCE);
		virtual void clear();
		virtual void DrawAll();
		virtual void PrintAll();
		
	private:
		ScintHitsCollection* fScintCollection;

		G4String ScintName;

		G4int Trk=0;
		G4int TrkParent=0;
		G4int EntryTrk=0;
		G4int TrkDecay=0;
		G4bool TrackOneIn = false;
		G4bool EntryCreated = false;
		G4String debug	= "";

		//G4int fBounce; 
		std::vector<G4int> fEvent, fScintNo, fParticleID, fNgamma, fNgammaSec, fNCer;
 		std::vector<G4int> fRight, fLeft, fDown, fUp, fBack, fFront;
		std::vector<G4double> fEin, fEdep, fEout, fDelta, fThetaIn, fTrackLength, fThetaOut, fDecayTime;

		std::vector<G4int> fAbsorption, fReflection;

		std::vector<G4double> fPosInX, fPosInY, fPosInZ; 
		std::vector<G4double> fPosOutX, fPosOutY, fPosOutZ; 
		std::vector<G4double> fMomInX, fMomInY, fMomInZ; 
		std::vector<G4double> fMomOutX, fMomOutY, fMomOutZ; 
		std::vector<G4double> fTimeIn, fTimeOut;
		
		G4ThreeVector fDirIn_trans, fDirOut_trans;
		G4ThreeVector fDirIn, fDirOut;

		G4int PhotonTime_ID, PhotonPositionFront_ID, PhotonPositionBack_ID, PhotonPositionLeft_ID;
		G4int PhotonPositionRight_ID, PhotonPositionUp_ID, PhotonPositionDown_ID;

		G4int hitsCID;
};

#endif


