/// \file RunAction.hh
/// \brief Definition of RunAction class

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"

#include <vector>

#include "TVector3.h"

class G4Run;
class TFile;
class TTree;
class RunActionMessenger;

/// Run action class
///
/// It prints out the data Tree


class RunAction : public G4UserRunAction {
	public:
		RunAction();
		virtual ~RunAction();

        //virtual G4Run* GenerateRun();
		virtual void BeginOfRunAction(const G4Run*);
		virtual void   EndOfRunAction(const G4Run*);
		
		void SetFileName(G4String name){fName = name;}

		inline TFile* GetFilePtr(){return fData;}
		inline TTree* GetTreePtr(){return fTree;}

		void SetEin (G4double val){fEin  = val;}
		void SetEdep(G4double val){fEdep = val;}
		void SetEout(G4double val){fEout = val;}
		void SetEdelta(G4double val){fDelta = val;}
		void SetThetaIn(G4double val){fThetaIn = val;}
		void SetTrackLength(G4double val){fTrackLength = val;}
		void SetThetaPositron(G4double val){fThetaPositron = val;}
		void SetID(G4int val){fID = val;}
		void SetBounce(G4int val){fBounce = val;}
		void SetPosX(std::vector<G4double> pos){fPosX = pos;}
		void SetPosY(std::vector<G4double> pos){fPosY = pos;}
		void SetPosZ(std::vector<G4double> pos){fPosZ = pos;}
		void SetTime(std::vector<G4double> time){fTime = time;}
		void SetNgamma(int ngamma){fNgamma = ngamma;}
		void SetNgammaSec(int ngammasec){fNgammaSec = ngammasec;}
		void SetCer(std::vector<G4int> cer){fCer = cer;}
		void SetThetaGamma(std::vector<G4double> thetagamma){fThetaGamma = thetagamma;}
		void SetTimeGamma(std::vector<G4double> timegamma){fTimeGamma = timegamma;}
		void SetEGamma(std::vector<G4double> egamma){fEGamma = egamma;}
		void SetNCer(G4int ncer){fNCer = ncer;}
		void SetDecayTime(G4double val){fDecayTime = val;}
		
		void SetCmdOCT(G4bool cmd){fCmdOCT = cmd;}
		G4bool GetCmdOCT(){return fCmdOCT;}

		void SetCmdDN(G4bool cmd){fCmdDN = cmd;}
		G4bool GetCmdDN(){return fCmdDN;}
		
		void SetCmdPhotons(G4int cmd){fCmdPhotons = cmd;}
		G4int GetCmdPhotons(){return fCmdPhotons;}

		void SetCmdTracks(G4int cmd){fCmdTracks = cmd;}
		G4int GetCmdTracks(){return fCmdTracks;}

		void SetCurrentRight(G4int right){fRight = right;}
		void SetCurrentLeft(G4int left){fLeft = left;}
		void SetCurrentDown(G4int down){fDown = down;}
		void SetCurrentUp(G4int up){fUp = up;}
		void SetCurrentBack(G4int back){fBack = back;}
		void SetCurrentFront(G4int front){fFront = front;}
		void SetSiPM(G4int sipm){fSiPM = sipm;}

		void SetNCells(G4int val){fNCells = val;}
		void SetNPhotoElectrons(G4double val){fNPhotoElectrons = val;}
		void SetCells(std::vector<G4int> val){fCells = val;}
		void SetCellTime(std::vector<G4double> val){fCellTime = val;}
		void SetOCTFlag(std::vector<G4int> val){fOCTflag = val;}
		void SetDNFlag(std::vector<G4int> val){fDNflag = val;}

		// Gun time
		inline void SetGunTimeMean(G4double val){fGunTimeMean = val;}
		inline G4double GetGunTime(){return fGunTime;}
		inline void AdvanceGunTime(){fGunTime += G4RandExponential::shoot(fGunTimeMean);}

		// Datk noise time
		inline void SetDNTimeMean(G4double val){fDNTimeMean = val;}
		inline G4double GetDNTime(){return fDNTime;}
		inline void AdvanceDNTime(){fDNTime += G4RandExponential::shoot(fDNTimeMean);}

	private:
		TFile* fData;
		TTree* fTree;
		// BC400 scorers
		G4double fEin;
		G4double fEdep;
		G4double fEout;
		G4double fDelta;
		G4double fThetaIn, fTrackLength, fThetaPositron;
		G4int fID;
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

		G4bool fCmdOCT, fCmdDN;
		G4int fCmdPhotons, fCmdTracks, fNCer;
		G4int fRight;
		G4int fLeft;
		G4int fDown;
		G4int fUp;
		G4int fBack;
		G4int fFront;
		G4int fSiPM;

		// SiPM scorers
		G4int fNCells;
		G4double fNPhotoElectrons;
		std::vector<G4int> fCells;
		std::vector<G4double> fCellTime;
		std::vector<G4int> fOCTflag;
		std::vector<G4int> fDNflag;

		//SiPM time counters
		G4double fGunTime, fDNTime;
		G4double fGunTimeMean, fDNTimeMean;
		
		G4double fDecayTime;


		G4String fName;
		RunActionMessenger* fMessenger;
};

#endif


